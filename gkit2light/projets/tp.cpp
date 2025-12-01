#include <vector>
#include <cfloat>
#include <chrono>
#include <random>
#include <cmath>
#include <cstdio>
#include <cassert>
#include <algorithm>

#include "vec.h"
#include "mat.h"
#include "color.h"
#include "image.h"
#include "image_io.h"
#include "image_hdr.h"
#include "orbiter.h"
#include "mesh.h"
#include "wavefront.h"

struct Ray
{
    Point  o;
    Vector d; 
    float  tmax;

    Ray(const Point& origine, const Point& extremite)
        : o(origine), d(Vector(origine, extremite)), tmax(1) {}

    Ray(const Point& origine, const Vector& direction)
        : o(origine), d(direction), tmax(FLT_MAX) {}

    Point point(const float t) const { return o + t * d; }
};

struct Hit
{
    float t;
    float u, v;
    int   triangle_id;

    Hit() : t(FLT_MAX), u(), v(), triangle_id(-1) {}
    Hit(const float _t, const float _u, const float _v, const int _id)
        : t(_t), u(_u), v(_v), triangle_id(_id) {}

    operator bool() const { return (triangle_id != -1); }
};

struct Triangle
{
    Point  p;     
    Vector e1,e2; 
    int    id;    

    Triangle(const TriangleData& data, const int _id)
        : p(data.a), e1(Vector(data.a, data.b)), e2(Vector(data.a, data.c)), id(_id) {}

    // Möller–Trumbore
    Hit intersect(const Ray& ray, const float tmax) const
    {
        const Vector pvec = cross(ray.d, e2);
        const float  det  = dot(e1, pvec);
        const float  eps  = 1e-8f;
        if (std::abs(det) < eps) return Hit();

        const float inv_det = 1.0f / det;
        const Vector tvec(p, ray.o);

        const float u = dot(tvec, pvec) * inv_det;
        if (u < 0.f || u > 1.f) return Hit();

        const Vector qvec = cross(tvec, e1);
        const float  v    = dot(ray.d, qvec) * inv_det;
        if (v < 0.f || u + v > 1.f) return Hit();

        const float t = dot(e2, qvec) * inv_det;
        if (t < 0.f || t > tmax) return Hit();

        return Hit(t, u, v, id);
    }
};

// Pour le rapport et pour reproduire c'est plus simple de faire avec des modes différents comme ça, j'ai mis a la fin du rapport comment l'utiliser !
enum Stage {
    STAGE_NORMALS = 0, // visu normales interp.
    STAGE_BARYCENTRIC = 1, // visu barycentres
    STAGE_DIRECT_SKY = 2, // éclairage direct ciel uniforme (directions)
    STAGE_DIRECT_EMITTER_AREA = 3, // éclairage direct sources (aire sampling)
    STAGE_INDIRECT_1_BOUNCE = 4 // direct + 1 rebond indirect
};

static inline const char* stage_name(int s) {
    switch (s) {
        case STAGE_NORMALS: return "normals";
        case STAGE_BARYCENTRIC: return "barycentric";
        case STAGE_DIRECT_SKY: return "direct_sky";
        case STAGE_DIRECT_EMITTER_AREA: return "direct_emit_area";
        case STAGE_INDIRECT_1_BOUNCE: return "indirect1";
        default: return "unknown";
    }
}

static inline bool is_emissive(const Material& m)
{
    const float eps = 1e-6f;
    return (m.emission.r > eps || m.emission.g > eps || m.emission.b > eps);
}

static inline Vector shading_normal(const Mesh& mesh, const Hit& hit)
{
    const TriangleData& data = mesh.triangle(hit.triangle_id);
    const float w = 1.f - hit.u - hit.v;
    Vector n = w * Vector(data.na) + hit.u * Vector(data.nb) + hit.v * Vector(data.nc);
    return normalize(n);
}

static inline Color albedo_at_hit(const Mesh& mesh, const Hit& hit)
{
    const Material& m = mesh.triangle_material(hit.triangle_id);
    return m.diffuse;
}

static inline Hit intersect_scene(const std::vector<Triangle>& tris, const Point& o, const Vector& d, float tmax=FLT_MAX)
{
    Ray ray(o, d);
    Hit   best;
    float best_t = tmax;
    for (const Triangle& T : tris)
        if (Hit h = T.intersect(ray, best_t))
            best = h, best_t = h.t;
    return best;
}

static inline bool visible_to_sky(const std::vector<Triangle>& tris, const Point& p, const Vector& l)
{
    return !intersect_scene(tris, p, l);
}

static inline Vector sample_uniform_hemisphere(float u1, float u2, float& cos_theta, float& pdf)
{
    cos_theta = u1;
    const float sin_theta = std::sqrt(std::max(0.f, 1.f - cos_theta*cos_theta));
    const float phi = 2.f * float(M_PI) * u2;

    pdf = 1.f / (2.f * float(M_PI));
    return Vector(std::cos(phi) * sin_theta, cos_theta, std::sin(phi) * sin_theta);
}

static inline void orthonormal_basis(const Vector& n, Vector& t, Vector& b)
{
    if (std::abs(n.y) < 0.999f) t = normalize(cross(Vector(0,1,0), n));
    else t = normalize(cross(Vector(1,0,0), n));
    b = cross(n, t);
}

static inline Vector sample_cosine_hemisphere(float u1, float u2, float& pdf)
{
    const float r   = std::sqrt(u1);
    const float phi = 2.f * float(M_PI) * u2;
    const float x = r * std::cos(phi);
    const float z = r * std::sin(phi);
    const float y = std::sqrt(std::max(0.f, 1.f - u1));

    pdf = y / float(M_PI);
    return Vector(x, y, z);
}

static inline Vector local_to_world(const Vector& l, const Vector& t, const Vector& b, const Vector& n)
{
    return l.x * t + l.y * n + l.z * b;
}

struct EmitterTri { int tri_id; float area; };

static inline float tri_area_world(const TriangleData& t, const Transform& model)
{
    const Point a = model(Point(t.a));
    const Point b = model(Point(t.b));
    const Point c = model(Point(t.c));

    return 0.5f * length(cross(Vector(a, b), Vector(a, c)));
}

static inline float build_emitters(const Mesh& mesh, const Transform& model,
                                   std::vector<EmitterTri>& emitters)
{
    emitters.clear();
    emitters.reserve(mesh.triangle_count());
    float area_sum = 0.f;

    for (int i = 0; i < mesh.triangle_count(); ++i)
    {
        const Material& mat = mesh.triangle_material(i);
        if (!is_emissive(mat)) continue;

        const TriangleData& td = mesh.triangle(i);
        const float A = tri_area_world(td, model);
        if (A <= 0.f) continue;

        emitters.push_back({ i, A });
        area_sum += A;
    }
    return area_sum;
}

struct SampledPoint {
    Point  q;
    Vector nq;
    Color  Le;
    float  pdfA;
    int    tri_id;
};

static inline Point bary_to_point(const TriangleData& t, float u, float v)
{
    const float w = 1.f - u - v;
    return Point(w*t.a.x + u*t.b.x + v*t.c.x,
                 w*t.a.y + u*t.b.y + v*t.c.y,
                 w*t.a.z + u*t.b.z + v*t.c.z);
}

static inline SampledPoint sample_emitter_point(const Mesh& mesh,
                                                const std::vector<EmitterTri>& emitters,
                                                float emitters_area_sum,
                                                std::default_random_engine& rng,
                                                std::uniform_real_distribution<float>& uni)
{
    float pick = uni(rng) * emitters_area_sum;
    int tri_id = -1;

    for (const auto& e : emitters) 
    {
        if (pick <= e.area) { tri_id = e.tri_id; break; }
        pick -= e.area;
    }

    if (tri_id < 0) tri_id = emitters.back().tri_id;

    const TriangleData& td = mesh.triangle(tri_id);
    const Material&     m  = mesh.triangle_material(tri_id);
    const float A = 0.5f * length(cross(Vector(td.a, td.b), Vector(td.a, td.c)));

    const float r1 = uni(rng), r2 = uni(rng);
    const float su = std::sqrt(r1);
    const float u = 1.f - su;
    const float v = su * r2;
    const Point  q = bary_to_point(td, u, v);

    const Vector nq = normalize(cross(Vector(td.a, td.b), Vector(td.a, td.c)));

    const float pdfA = 1.f / emitters_area_sum;

    return { q, nq, m.emission, pdfA, tri_id };
}

static inline bool visible_segment(const std::vector<Triangle>& tris, const Point& p, const Point& q, int light_tri_to_ignore)
{
    Vector d(p, q);
    float  dist = length(d);
    if (dist <= 0.f) return false;
    d = d / dist;

    const float eps = 1e-3f;
    const Point  p0 = p + eps * d;

    Ray  ray(p0, d);
    float tmax = dist - 2.f*eps;

    for (const Triangle& T : tris) 
    {
        if (T.id == light_tri_to_ignore) continue;
        if (Hit h = T.intersect(ray, tmax)) return false;
    }
    
    return true;
}


// Direct
static inline Color compute_direct_sky(const std::vector<Triangle>& tris,
                                       const Point& p, const Vector& n, const Color& albedo,
                                       std::default_random_engine& rng, std::uniform_real_distribution<float>& uni,
                                       int N)
{
    if (N <= 0) return Color(0,0,0);
    Color L(0,0,0);
    for (int k=0; k<N; ++k) {
        float u1 = uni(rng), u2 = uni(rng), cos_theta, pdf;
        const Vector l = sample_uniform_hemisphere(u1, u2, cos_theta, pdf);
        if (l.y <= 0.f) continue;
        const float cos_nl = std::max(0.f, dot(n, l));
        if (cos_nl <= 0.f) continue;
        if (!visible_to_sky(tris, p + 1e-3f*n, l)) continue;
        L = L + (albedo / float(M_PI)) * (cos_nl) * (1.f / pdf);
    }
    return L / float(N);
}

// émetteurs par échantillonnage d’aire
static inline Color compute_direct(const Mesh& mesh, const std::vector<Triangle>& tris,
                                   const std::vector<EmitterTri>& emitters, float emitters_area_sum,
                                   const Point& p, const Vector& n, const Color& albedo,
                                   std::default_random_engine& rng, std::uniform_real_distribution<float>& uni,
                                   int N)
{
    if (emitters.empty() || emitters_area_sum <= 0.f || N <= 0) return Color(0,0,0);

    Color L(0,0,0);
    for (int k=0; k<N; ++k)
    {
        const SampledPoint s = sample_emitter_point(mesh, emitters, emitters_area_sum, rng, uni);

        Vector l(p, s.q);
        float  r2 = dot(l,l);
        if (r2 <= 0.f) continue;
        const float r = std::sqrt(r2);
        l = l / r;

        const float cos_p = std::max(0.f, dot(n,   l));
        const float cos_q = std::max(0.f, dot(s.nq, -l));
        if (cos_p <= 0.f || cos_q <= 0.f) continue;

        if (!visible_segment(tris, p, s.q, s.tri_id)) continue;

        L = L + (albedo / float(M_PI)) * s.Le * (cos_p * cos_q / r2) * (1.f / s.pdfA);
    }
    return L / float(N);
}

// Indirect 1 rebond
static inline bool hit_is_emitter(const Mesh& mesh, const Hit& hit)
{
    return is_emissive(mesh.triangle_material(hit.triangle_id));
}

static inline Color compute_indirect_1bounce(const Mesh& mesh, const std::vector<Triangle>& tris,
                                             const std::vector<EmitterTri>& emitters, float emitters_area_sum,
                                             const Point& p, const Vector& n, const Color& albedo,
                                             std::default_random_engine& rng, std::uniform_real_distribution<float>& uni,
                                             int N_ind)
{
    if (N_ind <= 0) return Color(0,0,0);

    Color L(0,0,0);
    Vector t, b; orthonormal_basis(n, t, b);

    for (int j=0; j<N_ind; ++j)
    {
        float u1 = uni(rng), u2 = uni(rng), pdf_dir;
        const Vector l_local = sample_cosine_hemisphere(u1, u2, pdf_dir);
        if (pdf_dir <= 0.f) continue;

        const Vector l = local_to_world(l_local, t, b, n);

        const float eps = 1e-3f;
        const Point  p0 = p + eps * n;

        const Hit hy = intersect_scene(tris, p0, l);
        if (!hy) continue;

        if (hit_is_emitter(mesh, hy)) continue;

        const Vector ny = shading_normal(mesh, hy);
        const Color  cy = albedo_at_hit(mesh, hy);
        Point        y  = p0 + hy.t * l;
        y = y + eps * ny;

        const Color direct_y = compute_direct(mesh, tris, emitters, emitters_area_sum, y, ny, cy, rng, uni, 1);

        L = L + albedo * direct_y;
    }

    return L / float(N_ind);
}

int main(const int argc, const char** argv)
{
    const char* mesh_filename = "data/cornell.obj";
    const char* orbiter_filename = "data/cornell_orbiter.txt";
    if (argc > 1) mesh_filename = argv[1];
    if (argc > 2) orbiter_filename = argv[2];

    int stage = STAGE_INDIRECT_1_BOUNCE; // je laisse ça par défaut
    if (argc > 3) stage = std::atoi(argv[3]);

    Orbiter camera;
    if (camera.read_orbiter(orbiter_filename) < 0) return 1;

    Mesh mesh = read_mesh(mesh_filename);

    std::vector<Triangle> triangles;
    triangles.reserve(mesh.triangle_count());
    for (int i=0; i<mesh.triangle_count(); ++i)
        triangles.emplace_back(mesh.triangle(i), i);

    Image image(1024, 768);

    camera.projection(image.width(), image.height(), 45);
    const Transform model = Identity();
    const Transform view = camera.view();
    const Transform projection = camera.projection();
    const Transform viewport = camera.viewport();
    const Transform inv = Inverse(viewport * projection * view * model);

    std::vector<EmitterTri> emitters;
    const float emitters_area_sum = build_emitters(mesh, model, emitters);
    if (stage >= STAGE_DIRECT_EMITTER_AREA && emitters.empty())
        std::fprintf(stderr, "Warning: aucune source émissive dans la scène.\n");

    static thread_local std::default_random_engine rng(1337);
    static thread_local std::uniform_real_distribution<float> uni(0.f, 1.f);

    const int N_sky = 64;     // stage 2
    const int N_dir = 256;    // stage 3 & 4 (direct)
    const int N_ind = 16;     // stage 4 (indirect 1 rebond)

    auto start = std::chrono::high_resolution_clock::now();

    for (int y=0; y<image.height(); ++y)
    for (int x=0; x<image.width();  ++x)
    {
        const Point  o = inv(Point(x + 0.5f, y + 0.5f, 0));
        const Point  e = inv(Point(x + 0.5f, y + 0.5f, 1));
        const Ray    ray(o, e);

        Hit   hit;
        float tmax = ray.tmax;
        for (const Triangle& T : triangles)
            if (Hit h = T.intersect(ray, tmax))
                { assert(h.t > 0); hit = h; tmax = h.t; }

        if (!hit) { image(x,y) = Black(); continue; }

        const Vector n   = shading_normal(mesh, hit);
        const Color  alb = albedo_at_hit(mesh, hit);
        const Point  p   = ray.point(hit.t) + 1e-3f * n;

        Color L(0,0,0);

        switch (stage)
        {
            case STAGE_NORMALS:
                L = Color(std::abs(n.x), std::abs(n.y), std::abs(n.z));
                break;

            case STAGE_BARYCENTRIC:
                L = Color(1.f - hit.u - hit.v, hit.u, hit.v);
                break;

            case STAGE_DIRECT_SKY:
                L = compute_direct_sky(triangles, p, n, alb, rng, uni, N_sky);
                break;

            case STAGE_DIRECT_EMITTER_AREA:
                L = compute_direct(mesh, triangles, emitters, emitters_area_sum,
                                p, n, alb, rng, uni, N_dir);
                break;

            case STAGE_INDIRECT_1_BOUNCE:
            default:
            {
                Color Le_local(0,0,0);
                if (hit_is_emitter(mesh, hit) && dot(n, -ray.d) > 0.f)
                    Le_local = mesh.triangle_material(hit.triangle_id).emission;

                const Color Ld = compute_direct(mesh, triangles, emitters, emitters_area_sum, p, n, alb, rng, uni, N_dir);
                const Color Li = compute_indirect_1bounce(mesh, triangles, emitters, emitters_area_sum, p, n, alb, rng, uni, N_ind);
                L = Le_local + Ld + Li;
                break;
            }
        }

        image(x,y) = Color(std::min(1.f, L.r), std::min(1.f, L.g), std::min(1.f, L.b));
    }

    auto stop = std::chrono::high_resolution_clock::now();
    int cpu_ms = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
    std::printf("%dms\n", cpu_ms);

    char out_png[256]; std::snprintf(out_png, sizeof(out_png), "%s.png", stage_name(stage));
    char out_hdr[256]; std::snprintf(out_hdr, sizeof(out_hdr), "%s.hdr", stage_name(stage));

    write_image(image, out_png);
    write_image_hdr(image, out_hdr);
    std::printf("Saved: %s / %s\n", out_png, out_hdr);
    return 0;
}
