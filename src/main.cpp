//immediate TODO:
//- use Vert lerp() to do per-line interpolation?
//- split the single loop with per-line Vert pair selection into two dumb loops
    //- refactor the scanline inner loop into its own function?
//- early-out in the case that there are no rows (or no columns?) to be drawn
//- look into whether the by-value swaps are generating unnecessary code
    //- if they are, swap pointers around instead
//- write a SIMD-optimized lerp() that treats the Vert as a bag of floats?

//- figure out how MSAA is going to work
//- implement MSAA



//TODO: can we optimize metalness by treating it as a factor in the fresnel approximation?



//- sort out what code needs to get put into headers to be shared between rasterizers
//- move common.hpp, math.hpp, List.hpp, and FileBuffer.hpp into their own directory
//- label and update all those separator comments
//- figure out how to do

//global TODO:
//- sorting the render queue
//- point diffuse lighting

//- better test scene!

//- improved shadows? (resolution-dependent shading, take gradient/facing into account)
//- cascaded shadows?
//- mipmaps?
//- MSAA?
//- cubemap reflections?
//- cubemap shadows?

//- floating point framebuffer
//- bloom
//- tone mapping


//- finish debug stats
//- more specific performance counters
//- overdraw visualization
//- stats about wasted SIMD lanes

//- vsync?

#include "common.hpp"

static_assert(sizeof(size_t) == 8, "must be compiled in 64-bit mode");

#include "blit.hpp"
#include "math.hpp"
#include "List.hpp"
#include "FileBuffer.hpp"

#include "stb_image.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <SDL.h>

#ifdef _MSC_VER
    #include "intrin.h"
#else
    // #include "x86intrin.h"
    #define _bit_scan_reverse __builtin_ia32_bsrsi
    #define __rdtscp __builtin_ia32_rdtscp
#endif

////////////////////////////////////////////////////////////////////////////////
/// DYNAMIC LIBRARY LOADING                                                  ///
////////////////////////////////////////////////////////////////////////////////

#ifdef _WIN32

Uint64 (* GetPerformanceCounter) ();
Uint64 (* GetPerformanceFrequency) ();
int (* Init) (Uint32 flags);
const char * (* GetError) ();
SDL_Window * (* CreateWindow) (const char* title, int x, int y, int w, int h, Uint32 flags);
int (* PollEvent) (SDL_Event * event);
SDL_Surface * (* GetWindowSurface) (SDL_Window * window);
int (* UpdateWindowSurface) (SDL_Window * window);
int (* SetRelativeMouseMode) (SDL_bool enabled);

#include <windows.h>

template <typename TYPE>
bool load_func(HMODULE handle, TYPE * func, const char * name) {
    *func = (TYPE) GetProcAddress(handle, name);
    if (*func == nullptr) {
        printf("FAILED TO LOAD FUNCTION %s\n", name);
        printf("Error code: %d\n", GetLastError());
        return false;
    }
    return true;
}

//we don't care to unload the DLL before application exit,
//so it's fine that we just throw away the handle inside these functions
//instead of returning it like we're "supposed" to
bool load_sdl_functions(const char * filepath) {
    HMODULE handle = LoadLibraryA(filepath);
    if (!handle) {
        printf("FAILED TO LOAD DYNAMIC LIBRARY: %s\n", filepath);
        printf("Error code: %d\n", GetLastError());
        return false;
    }

    if (!load_func(handle, &GetPerformanceCounter, "SDL_GetPerformanceCounter"))     return false;
    if (!load_func(handle, &GetPerformanceFrequency, "SDL_GetPerformanceFrequency")) return false;
    if (!load_func(handle, &Init, "SDL_Init"))                                       return false;
    if (!load_func(handle, &GetError, "SDL_GetError"))                               return false;
    if (!load_func(handle, &CreateWindow, "SDL_CreateWindow"))                       return false;
    if (!load_func(handle, &PollEvent, "SDL_PollEvent"))                             return false;
    if (!load_func(handle, &GetWindowSurface, "SDL_GetWindowSurface"))               return false;
    if (!load_func(handle, &UpdateWindowSurface, "SDL_UpdateWindowSurface"))         return false;
    if (!load_func(handle, &SetRelativeMouseMode, "SDL_SetRelativeMouseMode"))       return false;

    return true;
}

#define SDL_GetPerformanceCounter GetPerformanceCounter
#define SDL_GetPerformanceFrequency GetPerformanceFrequency
#define SDL_Init Init
#define SDL_GetError GetError
#define SDL_CreateWindow CreateWindow
#define SDL_PollEvent PollEvent
#define SDL_GetWindowSurface GetWindowSurface
#define SDL_UpdateWindowSurface UpdateWindowSurface
#define SDL_SetRelativeMouseMode SetRelativeMouseMode

#endif // _WIN32

////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

u64 applicationStartupTimeValue;

double get_time() {
    u64 currentTimeValue = SDL_GetPerformanceCounter();
    u64 diffTimeValue = currentTimeValue - applicationStartupTimeValue;
    double elapsedSeconds = (double)diffTimeValue / (double)SDL_GetPerformanceFrequency();
    return elapsedSeconds;
}

//DEBUG GLOBALS
int totalTris;
int drawnTris;
int clippedTris;
int shadowTotalTris;
int shadowDrawnTris;

u64 perfFill, perfText, perfShadow, perfDraw, perfTransform, perfRasterize, perfInner, perfBlit;
uint perfDummy;

u64 perf() {
    // return __rdtsc();
    return __rdtscp(&perfDummy);
    // return __builtin_ia32_rdtscp(&perfDummy);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Image {
    Pixel * pixels;
    int w, h;
};

Image load_image(const char * filename) {
    int w, h, c;
    Pixel * pixels = (Pixel *) stbi_load(filename, &w, &h, &c, 4);
    if (!pixels) { printf("ERROR: %s missing\n", filename); }
    assert(pixels);
    return { pixels, w, h };
}

struct Tex {
    Pixel * pixels;
    Pixel * normals;
    int width, height;
    u32 wmask, hmask;
};

Tex load_texture(const char * diffuse, const char * normals) {
    Image diff = load_image(diffuse);
    Image norm = load_image(normals);
    assert(diff.w == norm.w && diff.h == norm.h);
    Tex ret = { diff.pixels, norm.pixels, diff.w, diff.h, (u32)(diff.w - 1), (u32)(diff.h - 1) };
    assert(!(ret.width & ret.wmask) && !(ret.height & ret.hmask));
    return ret;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// MODEL DATA FORMAT                                                                            ///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct Vertex {
    Vec3 p, n, t, b;
    float u, v;
};

struct Vert {
    Vec4 p; //position
    Vec3 n, l, c, s; //normal, light (directional), camera, shadow
    float u, v; //texture coords
};

typedef u16 idx_t;
struct Triangle {
    idx_t v1, v2, v3;
};

struct Model {
    size_t vertCount;
    Vertex * vertices;
    Vert * verts; //allocated within the game

    size_t triCount;
    Triangle * triangles;

    Tex tex;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// MODEL LOADING                                                                                ///
////////////////////////////////////////////////////////////////////////////////////////////////////

Model * load_model(const char * path) {
    Model * model = (Model *) read_entire_file(path);

    //patch relative pointers
    size_t base = (size_t) model;
    model->vertices = (Vertex *) ((size_t) model->vertices + base);
    model->triangles = (Triangle *) ((size_t) model->triangles + base);

    //allocate post-transform arrays
    model->verts = (Vert *) malloc(model->vertCount * sizeof(Vert));

    return model;
}

void free_model(Model * model) {
    free(model->verts);
    free(model);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// RASTERIZATION                                                                                ///
////////////////////////////////////////////////////////////////////////////////////////////////////

inline Vec4 color_sq(Vec4 v) {
    return { v.x * v.x, v.y * v.y, v.z * v.z, v.w };
}

void draw_triangle(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]) {
    ++drawnTris;

    //convert triangles to edges
    struct Edge {
        Vert v1, v2;
        float yfactor;
    } edges[3];
    float miny = canvas->height, maxy = -1;
    for (int i : range(3)) {
        Vert v1 = tri[i];
        Vert v2 = tri[(i + 1) % 3];

        //update the triangle's vertical extent
        miny = v1.p.y < miny? v1.p.y : miny;
        maxy = v1.p.y > maxy? v1.p.y : maxy;

        //sort vertices by y
        edges[i].v1 = v1.p.y < v2.p.y? v1 : v2;
        edges[i].v2 = v1.p.y < v2.p.y? v2 : v1;

        edges[i].yfactor = 1.0f/(edges[i].v2.p.y - edges[i].v1.p.y);
    }

    //convert the triangle's vertical extent to pixels
    int firstLine = miny + 1;
    int lastLine = maxy;
    //clamp vertical extent of triangle to within the screen for rasterization
    if (firstLine < 0) firstLine = 0;
    if (lastLine > canvas->height - 1) lastLine = canvas->height - 1;

    for (int y = firstLine; y <= lastLine; ++y) {
        Pixel * row = canvas->pixels + y * canvas->pitch;
        float * zrow = canvas->depth + y * canvas->zpitch;

        //the current pixel row will be within the vertical extend of only two
        //of the three edges at any time, so find those two and discard the third
        //SLOW: selecting the right edge for every line is unnecessary
        //TODO: refactor this to use two for loops, like that one guy did
        Edge e1, e2;
        if (y < edges[0].v1.p.y || y > edges[0].v2.p.y) {
            e1 = edges[1];
            e2 = edges[2];
        } else if (y < edges[1].v1.p.y || y > edges[1].v2.p.y) {
            e1 = edges[0];
            e2 = edges[2];
        } else {
            e1 = edges[0];
            e2 = edges[1];
        }

        //calculate vertical blend amounts for this scanline
        float f1a = (e1.v2.p.y - y) * e1.yfactor;
        float f2a = (e2.v2.p.y - y) * e2.yfactor;
        float f1b = 1 - f1a;
        float f2b = 1 - f2a;

        //find intersection with each edge by interpolating x along the edge
        float x1 = f1a * e1.v1.p.x + f1b * e1.v2.p.x;
        float x2 = f2a * e2.v1.p.x + f2b * e2.v2.p.x;

        //sort edges based on intersections
        float minx = x1 < x2? x1 : x2;
        float maxx = x1 > x2? x1 : x2;
        if (x1 > x2) {
            swap(e1, e2);
            swap(f1a, f2a);
            swap(f1b, f2b);
        }

        //interpolate vertex attributes at intersection points
        float z1 = f1a * e1.v1.p.z + f1b * e1.v2.p.z;
        float z2 = f2a * e2.v1.p.z + f2b * e2.v2.p.z;
        float w1 = f1a * e1.v1.p.w + f1b * e1.v2.p.w;
        float w2 = f2a * e2.v1.p.w + f2b * e2.v2.p.w;
        Vec3 l1 = f1a * e1.v1.l + f1b * e1.v2.l;
        Vec3 l2 = f2a * e2.v1.l + f2b * e2.v2.l;
        Vec3 n1 = f1a * e1.v1.n + f1b * e1.v2.n;
        Vec3 n2 = f2a * e2.v1.n + f2b * e2.v2.n;
        Vec3 c1 = f1a * e1.v1.c + f1b * e1.v2.c;
        Vec3 c2 = f2a * e2.v1.c + f2b * e2.v2.c;
        Vec3 s1 = f1a * e1.v1.s + f1b * e1.v2.s;
        Vec3 s2 = f2a * e2.v1.s + f2b * e2.v2.s;
        float u1 = f1a * e1.v1.u + f1b * e1.v2.u;
        float u2 = f2a * e2.v1.u + f2b * e2.v2.u;
        float v1 = f1a * e1.v1.v + f1b * e1.v2.v;
        float v2 = f2a * e2.v1.v + f2b * e2.v2.v;

        //factor for interpolating vertex attributes horizontally
        float xfactor = 1.0f/(maxx - minx);

        //convert horizontal extent to pixels
        int first = minx + 1;
        int last = maxx;
        //clamp horizontal extent of scanline to pixels
        if (first < 0) first = 0;
        if (last > canvas->width - 1) last = canvas->width - 1;

        u64 preInner = perf();
        for (int x = first; x <= last; ++x) {
            //calculate horizontal interpolation factor for this pixel
            float fa = (maxx - x) * xfactor;
            float fb = 1 - fa;

            //interpolate vertex attributes for this pixel
            float z = fa * z1 + fb * z2;

            //depth test early-out
            if (z > 1.0f || z < -1.0f || z >= zrow[x]) {
                continue;
            }

            float w = 1 / (fa * w1 + fb * w2);
            Vec3 l = w * (fa * l1 + fb * l2);
            Vec3 n = w * (fa * n1 + fb * n2);
            Vec3 c = w * (fa * c1 + fb * c2);
            Vec3 sh = w * (fa * s1 + fb * s2);
            float u = w * (fa * u1 + fb * u2);
            float v = w * (fa * v1 + fb * v2);

            //texture sample calculations
            float uf = u - (int)u;
            float vf = v - (int)v;
            int iu1 = (int) u      & tex->wmask;
            int iv1 = (int) v      & tex->hmask;
            int iu2 = (int)(u + 1) & tex->wmask;
            int iv2 = (int)(v + 1) & tex->hmask;

            //sample diffuse texture
            Pixel p11 = tex->pixels[tex->width * iv1 + iu1];
            Pixel p12 = tex->pixels[tex->width * iv1 + iu2];
            Pixel p21 = tex->pixels[tex->width * iv2 + iu1];
            Pixel p22 = tex->pixels[tex->width * iv2 + iu2];
            Vec4 c11 = color_sq(vec4(p11.b, p11.g, p11.r, p11.a) * (1.0f / 255));
            Vec4 c12 = color_sq(vec4(p12.b, p12.g, p12.r, p12.a) * (1.0f / 255));
            Vec4 c21 = color_sq(vec4(p21.b, p21.g, p21.r, p21.a) * (1.0f / 255));
            Vec4 c22 = color_sq(vec4(p22.b, p22.g, p22.r, p22.a) * (1.0f / 255));
            Vec4 c1 = c11 * (1 - uf) + c12 * uf;
            Vec4 c2 = c21 * (1 - uf) + c22 * uf;
            Vec4 diffRough = c1 * (1 - vf) + c2 * vf;
            Vec3 color = vec3(diffRough);
            float roughness = diffRough.w;
            float rough = roughness * roughness * roughness * roughness * roughness;

            //sample normal texture
            Pixel m11 = tex->normals[tex->width * iv1 + iu1];
            Pixel m12 = tex->normals[tex->width * iv1 + iu2];
            Pixel m21 = tex->normals[tex->width * iv2 + iu1];
            Pixel m22 = tex->normals[tex->width * iv2 + iu2];
            Vec4 n11 = vec4(m11.r, m11.g, m11.b, m11.a);
            Vec4 n12 = vec4(m12.r, m12.g, m12.b, m12.a);
            Vec4 n21 = vec4(m21.r, m21.g, m21.b, m21.a);
            Vec4 n22 = vec4(m22.r, m22.g, m22.b, m22.a);
            Vec4 n1 = n11 * (1 - uf) + n12 * uf;
            Vec4 n2 = n21 * (1 - uf) + n22 * uf;
            Vec4 norMetal = n1 * (1 - vf) + n2 * vf;
            Vec3 normal = vec3(norMetal);
            float metalness = norMetal.w * (1.0f / 255);
            // metalness = 0;

            normal = nor(normal - vec3(127.5f));
            l = nor(l);
            c = nor(c);

            float shad = 1;
            if (sh.x > 0 && sh.y > 0 && sh.x < shadow->width - 1 && sh.y < shadow->height - 1) {
                //shadow sample calculations
                float sf = sh.x - (int)sh.x;
                float tf = sh.y - (int)sh.y;
                int is1 = (int) sh.x      & shadow->wmask;
                int it1 = (int) sh.y      & shadow->hmask;
                int is2 = (int)(sh.x + 1) & shadow->wmask;
                int it2 = (int)(sh.y + 1) & shadow->hmask;

                //sample shadow texture
                float s11 = shadow->depth[shadow->width * it1 + is1];
                float s12 = shadow->depth[shadow->width * it1 + is2];
                float s21 = shadow->depth[shadow->width * it2 + is1];
                float s22 = shadow->depth[shadow->width * it2 + is2];
                float d11 = fmin(1, fmax(0, (sh.z - s11) * shadow->scale * 4));
                float d12 = fmin(1, fmax(0, (sh.z - s12) * shadow->scale * 4));
                float d21 = fmin(1, fmax(0, (sh.z - s21) * shadow->scale * 4));
                float d22 = fmin(1, fmax(0, (sh.z - s22) * shadow->scale * 4));
                float d1 = d11 * (1 - sf) + d12 * sf;
                float d2 = d21 * (1 - sf) + d22 * sf;
                shad = d1 * (1 - tf) + d2 * tf;
            }

            //directional diffuse
            float light = fmax(0, dot(l, normal));
            color *= light * shad + 0.04f;

            //directional specular
            float e = 1 / (rough + 0.0001f);
            float m = 1 / (rough + 0.02f); //multiplier
            Vec3 reflected = 2 * dot(l, normal) * normal - l;
            float spec = fmax(0, dot(reflected, c));
            //specular = M * S / (E - E * S + S)
            float specular = m * spec / (e - e * spec + spec) + 0.02f;

            //calculate fresnel factors
            float ior = 1.45; //kinda arbitrary?
            float ior2 = 20;//13.33;
            float f = sq((1 - ior) / (1 + ior));
            float f2 = sq((1 - ior2) / (1 + ior2));
            float headon = fmax(0, dot(c, normal));
            //fresnel = F + (1 - R) * (1 - F) * sq(sq(1 - C)) * (1 - C)
            float fresnel = f + sq(1 - roughness) * (1 - f) * sq(sq(1 - headon)) * (1 - headon);
            float fresnel2 = f2 + sq(1 - roughness) * (1 - f2) * sq(sq(1 - headon)) * (1 - headon);

            //calculate dielectric output
            Vec3 highlight = vec3(191, 127, 0) * (1.0f / 255) * 0.2f * fmax(0, n.y);
            highlight += vec3(specular * shad);
            highlight *= fresnel;
            Vec3 diffColor = color * (1 - fresnel);
            Vec3 dielectric = diffColor + highlight;

            //calculate metal output
            Vec3 metalHighlight = vec3(191, 127, 0) * (1.0f / 255) * 0.2f * fmax(0, n.y);
            metalHighlight += vec3(specular * shad);
            metalHighlight *= fresnel2 * vec3(diffRough);
            Vec3 metalColor = color * (1 - fresnel2);
            Vec3 metal = metalColor + metalHighlight;

            //blend using metalness
            Vec3 out = metalness * metal + (1 - metalness) * dielectric;

            row[x] = { (u8)(fmin(1, sqrtf(diffColor.x + highlight.x)) * 255),
                       (u8)(fmin(1, sqrtf(diffColor.y + highlight.y)) * 255),
                       (u8)(fmin(1, sqrtf(diffColor.z + highlight.z)) * 255), 255 };

            row[x] = { (u8)(fmin(1, sqrtf(out.x)) * 255),
                       (u8)(fmin(1, sqrtf(out.y)) * 255),
                       (u8)(fmin(1, sqrtf(out.z)) * 255), 255 };


            zrow[x] = z;
        }
        perfInner += perf() - preInner;
    }
}

//PRECOND: a != b
float inv_lerp(float a, float x, float b) {
    return (x - a)/(b - a);
}

template <typename TYPE>
TYPE lerp(TYPE a, float f, TYPE b) {
    return (1 - f) * a + f * b;
}

Vert lerp(Vert a, float f, Vert b) {
    return {
        lerp(a.p, f, b.p),
        lerp(a.n, f, b.n),
        lerp(a.l, f, b.l),
        lerp(a.c, f, b.c),
        lerp(a.s, f, b.s),
        lerp(a.u, f, b.u),
        lerp(a.v, f, b.v),
    };
}

void maybe_draw_triangle(Canvas * canv, ZBuffer * shadow, Tex * tex, Vert tri[3]) {
    //cull back faces
    //alternative method: https://www.geeksforgeeks.org/orientation-3-ordered-points/
    if ((tri[1].p.y - tri[0].p.y) * (tri[2].p.x - tri[1].p.x) -
        (tri[1].p.x - tri[0].p.x) * (tri[2].p.y - tri[1].p.y) <= 0) {
        return;
    }

    if (tri[0].p.x < 0 && tri[1].p.x < 0 && tri[2].p.x < 0) return;
    if (tri[0].p.x > canv->width && tri[1].p.x > canv->width && tri[2].p.x > canv->width) return;

    if (tri[0].p.y < 0 && tri[1].p.y < 0 && tri[2].p.y < 0) return;
    if (tri[0].p.y > canv->height && tri[1].p.y > canv->height && tri[2].p.y > canv->height) return;

    if (tri[0].p.z < -1 && tri[1].p.z < -1 && tri[2].p.z < -1) return;
    if (tri[0].p.z >  1 && tri[1].p.z >  1 && tri[2].p.z >  1) return;

    u64 preRasterize = perf();
    draw_triangle(canv, shadow, tex, tri);
    // draw_triangle_sse2_plain(canv, shadow, tex, tri);
    perfRasterize += perf() - preRasterize;
}

void draw_model(Canvas * canvas, ZBuffer * shadow, Model * model,
                Mat4 modelMat, Mat4 viewMat, Mat4 projMat, Mat4 shadowMat,
                Vec3 camPos, Vec3 lightDir)
{
    Mat4 viewproj = projMat * viewMat;
    Mat4 modelviewproj = viewproj * modelMat;
    Mat4 modelshadow = shadowMat * modelMat;
    Mat3 normalMat = inverse_transpose(mat3(modelMat));

    u64 preTransform = perf();
    //vertex transform
    for (int i : range(model->vertCount)) {
        Vertex v = model->vertices[i];
        Vec4 p = modelviewproj * vec4(v.p, 1);

        //perspective divide
        p.w = 1 / p.w;
        p.x *= p.w;
        p.y *= p.w;
        p.z *= p.w;

        //transform from clip space to screen space
        p.x = (p.x + 1) *  0.5f * canvas->width;
        p.y = (p.y - 1) * -0.5f * canvas->height;

        Mat3 tbn = mat3(noz(normalMat * v.t),
                        noz(normalMat * v.b),
                        noz(normalMat * v.n));
        //XXX: why does this work? shouldn't I have to multiply by the inverse of the TBN matrix?
        //     am I calculating the TBN wrong in such a way that it cancels out?
        Vec3 l = tbn * lightDir;
        Vec3 c = tbn * (camPos - vec3(modelMat * vec4(v.p, 1)));

        Vec3 s = vec3(modelshadow * vec4(v.p, 1));
        //transform from clip space to screen space
        s.x = (s.x + 1) *  0.5f * shadow->width;
        s.y = (s.y - 1) * -0.5f * shadow->height;

        //NOTE: assumes uv will be > -2 to make rounding faster for texture fetch
        //      (floor is slower than round-toward-zero, even in SIMD (I think))
        v.u += 2 - 0.5f / model->tex.width;
        v.v += 2 - 0.5f / model->tex.height;

        //premultiply uv by texture size
        v.u *= model->tex.width;
        v.v *= model->tex.height;

        Vec3 n = noz(normalMat * v.n);

        model->verts[i] = { p, n * p.w, l * p.w, c * p.w, s * p.w, v.u * p.w, v.v * p.w };
    }
    perfTransform += perf() - preTransform;

    //draw triangles
    for (int i : range(model->triCount)) {
        Triangle t = model->triangles[i];
        Vert tri[3] = {
            model->verts[t.v1],
            model->verts[t.v2],
            model->verts[t.v3],
        };

        ++totalTris;

        //precompute these because it makes the rest of the code more readable
        bool w0 = tri[0].p.w < 0;
        bool w1 = tri[1].p.w < 0;
        bool w2 = tri[2].p.w < 0;

        //cull triangles entirely behind the camera
        if (w0 && w1 && w2) {
            continue;
        }

        if (w0 || w1 || w2) {
            ++clippedTris;

            //rotate vert order so that tri[0] is the odd one out, while preserving winding order
            if (w1 != w0 && w1 != w2) {
                Vert tmp = tri[0];
                tri[0] = tri[1];
                tri[1] = tri[2];
                tri[2] = tmp;
            } else if (w2 != w0 && w2 != w1) {
                Vert tmp = tri[0];
                tri[0] = tri[2];
                tri[2] = tri[1];
                tri[1] = tmp;
            }

            //determine which clip scheme we need to use
            //based on whether tri[0] is behind or in front of the camera
            //NOTE: w0, w1, w2 are no longer valid since vert order may have been rotated
            if (tri[0].p.w < 0) {
                Vert a = lerp(tri[1], inv_lerp(tri[1].p.w, 100, tri[0].p.w), tri[0]);
                Vert b = lerp(tri[2], inv_lerp(tri[2].p.w, 100, tri[0].p.w), tri[0]);

                //NOTE: there are several possible orderings that would preserve winding order
                //      the one we use here was chosen arbitrarily
                Vert tri1[3] = { tri[1], tri[2], b };
                Vert tri2[3] = { tri[1], b, a };

                maybe_draw_triangle(canvas, shadow, &model->tex, tri1);
                maybe_draw_triangle(canvas, shadow, &model->tex, tri2);
            } else {
                tri[1] = lerp(tri[1], inv_lerp(tri[1].p.w, 100, tri[0].p.w), tri[0]);
                tri[2] = lerp(tri[2], inv_lerp(tri[2].p.w, 100, tri[0].p.w), tri[0]);

                maybe_draw_triangle(canvas, shadow, &model->tex, tri);
            }
        } else {
            //no need to clip
            maybe_draw_triangle(canvas, shadow, &model->tex, tri);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

void shadow_triangle(ZBuffer * canvas, Vec3 tri[3]) {
    ++shadowDrawnTris;

    //convert triangles to edges
    struct Edge {
        Vec3 v1, v2;
        float yfactor;
    } edges[3];
    float miny = canvas->height, maxy = -1;
    for (int i : range(3)) {
        Vec3 v1 = tri[i];
        Vec3 v2 = tri[(i + 1) % 3];

        //update the triangle's vertical extent
        miny = v1.y < miny? v1.y : miny;
        maxy = v1.y > maxy? v1.y : maxy;

        //sort vertices by y
        edges[i].v1 = v1.y < v2.y? v1 : v2;
        edges[i].v2 = v1.y < v2.y? v2 : v1;

        edges[i].yfactor = 1.0f/(edges[i].v2.y - edges[i].v1.y);
    }

    //convert the triangle's vertical extent to pixels
    int firstLine = miny + 1;
    int lastLine = maxy;
    //clamp vertical extent of triangle to within the screen for rasterization
    if (firstLine < 0) firstLine = 0;
    if (lastLine > canvas->height - 1) lastLine = canvas->height - 1;

    for (int y = firstLine; y <= lastLine; ++y) {
        float * zrow = canvas->depth + y * canvas->width;

        //the current pixel row will be within the vertical extend of only two
        //of the three edges at any time, so find those two and discard the third
        //SLOW: selecting the right edge for every line is unnecessary
        Edge e1, e2;
        if (y < edges[0].v1.y || y > edges[0].v2.y) {
            e1 = edges[1];
            e2 = edges[2];
        } else if (y < edges[1].v1.y || y > edges[1].v2.y) {
            e1 = edges[0];
            e2 = edges[2];
        } else {
            e1 = edges[0];
            e2 = edges[1];
        }

        //calculate vertical blend amounts for this scanline
        float f1a = (e1.v2.y - y) * e1.yfactor;
        float f2a = (e2.v2.y - y) * e2.yfactor;
        float f1b = 1 - f1a;
        float f2b = 1 - f2a;

        //find intersection with each edge by interpolating x along the edge
        float x1 = f1a * e1.v1.x + f1b * e1.v2.x;
        float x2 = f2a * e2.v1.x + f2b * e2.v2.x;

        //sort edges based on intersections
        float minx = x1 < x2? x1 : x2;
        float maxx = x1 > x2? x1 : x2;
        if (x1 > x2) {
            swap(e1, e2);
            swap(f1a, f2a);
            swap(f1b, f2b);
        }

        //interpolate vertex attributes at intersection points
        float z1 = f1a * e1.v1.z + f1b * e1.v2.z;
        float z2 = f2a * e2.v1.z + f2b * e2.v2.z;

        //factor for interpolating vertex attributes horizontally
        float xfactor = 1.0f/(maxx - minx);

        //convert horizontal extent to pixels
        int first = minx + 1;
        int last = maxx;
        //clamp horizontal extent of scanline to pixels
        if (first < 0) first = 0;
        if (last > canvas->width - 1) last = canvas->width - 1;

        for (int x = first; x <= last; ++x) {
            //calculate horizontal interpolation factor for this pixel
            float fa = (maxx - x) * xfactor;
            float fb = 1 - fa;

            //interpolate vertex attributes for this pixel
            float z = fa * z1 + fb * z2;

            //depth test early-out
            if (z <= zrow[x]) {
                continue;
            }

            zrow[x] = z;
        }
    }
}

void draw_shadow(ZBuffer * canv, Model * model, Mat4 modelMat, Mat4 viewMat) {
    Mat4 modelview = viewMat * modelMat;

    //TODO: have this pre-allocated?
    Vec3 * verts = (Vec3 *) malloc(model->vertCount * sizeof(Vec3));

    for (int i : range(model->vertCount)) {
        Vertex v = model->vertices[i];
        Vec4 p = modelview * vec4(v.p, 1);

        //transform from clip space to screen space
        p.x = (p.x + 1) *  0.5f * canv->width;
        p.y = (p.y - 1) * -0.5f * canv->height;

        verts[i] = vec3(p);
    }

    for (int i : range(model->triCount)) {
        Triangle t = model->triangles[i];
        Vec3 tri[3] = {
            verts[t.v1],
            verts[t.v2],
            verts[t.v3],
        };

        ++shadowTotalTris;

        //cull back faces
        //alternative method: https://www.geeksforgeeks.org/orientation-3-ordered-points/
        if ((tri[1].y - tri[0].y) * (tri[2].x - tri[1].x) -
            (tri[1].x - tri[0].x) * (tri[2].y - tri[1].y) <= 0) {
            continue;
        }

        if (tri[0].x < 0 && tri[1].x < 0 && tri[2].x < 0) continue;
        if (tri[0].x > canv->width && tri[1].x > canv->width && tri[2].x > canv->width) continue;

        if (tri[0].y < 0 && tri[1].y < 0 && tri[2].y < 0) continue;
        if (tri[0].y > canv->height && tri[1].y > canv->height && tri[2].y > canv->height) continue;

        shadow_triangle(canv, tri);
    }

    free(verts);
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct MonoFont {
    Pixel * pixels;
    int texWidth, texHeight;
    int glyphWidth, glyphHeight;
    int cols, rows;
};

MonoFont load_mono_font(const char * filepath, int cols, int rows) {
    Image image = load_image(filepath);
    assert(image.w % cols == 0 && image.h % rows == 0);
    return { image.pixels, image.w, image.h, image.w / cols, image.h / rows, cols, rows };
}

Pixel operator|(Pixel l, Pixel r) {
    return { (u8)(l.b | r.b), (u8)(l.g | r.g), (u8)(l.r | r.r), (u8)(l.a | r.a) };
}

Pixel & operator|=(Pixel & l, Pixel r) {
    l = l | r;
    return l;
}

//UNSAFE: does not do any bounds checking!
void draw_glyph(Canvas * canvas, MonoFont * font, int cx, int cy, int glyph) {
    int row = glyph / font->cols;
    int col = glyph % font->cols;
    int srcx = col * font->glyphWidth;
    int srcy = row * font->glyphHeight;

    for (int y : range(font->glyphHeight)) {
        Pixel * dest = canvas->pixels + (cy + y) * canvas->pitch + cx;
        float * zdest = canvas->depth + (cy + y) * canvas->width + cx;
        Pixel * src = font->pixels + (srcy + y) * font->texWidth + srcx;
        for (int x : range(font->glyphWidth)) {
            dest[x] |= src[x];
            if (src[x].b | src[x].g | src[x].r) {
                //TODO: track how many pixels covered by text?
                zdest[x] = -1;
            }
        }
    }
}

void draw_text(Canvas * canvas, MonoFont * font, int x, int y, const char * text) {
    for (int i = 0; text[i] != '\0'; ++i) {
        draw_glyph(canvas, font, x + i * font->glyphWidth, y, text[i]);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

//TODO: re-use `used` buffers calculated during meta palette choice for palette generation
//TODO: SIMD-ize "pre-amble" code
//TODO: improve memory usage of trie by only allocating as many colors per node
//      as are used globally, and only bumping the pointer by as many as are used locally

struct RawFrame {
    Pixel * base;
    Pixel * pixels;
    int pitch;
};

//TODO: find a way to write the data blocks without storing up a buffer of 255 bytes
struct BlockBuffer {
    u16 bits;
    u8 bytes[257]; //up to 12 bits can be written at once, so we need 2 extra "overflow" bytes
};

//XXX: this is very slow because of how it uses the file buffer!
void put_code(FileBuffer * buf, BlockBuffer * block, int bits, u16 code) {
    //insert new code into block buffer
    int idx = block->bits / 8;
    int bit = block->bits % 8;
    block->bytes[idx + 0] |= code <<       bit      ;
    block->bytes[idx + 1] |= code >>  (8 - bit)     ;
    block->bytes[idx + 2] |= code >> ((8 - bit) + 8);
    block->bits += bits;

    //flush the block buffer if it's full
    if (block->bits >= 255 * 8) {
        buf->check(256);
        buf->write_unsafe<u8>(255);
        buf->write_block_unsafe(block->bytes, 255);

        block->bits -= 255 * 8;
        block->bytes[0] = block->bytes[255];
        block->bytes[1] = block->bytes[256];
        memset(block->bytes + 2, 0, 255);
    }
}

//TODO: replace this array list implementation with a trie
struct Code {
    u16 len; //number of bytes represented by this code
    u16 idx; //index into byte list
};

struct TrieNode {
    i16 next[256];
};

void reset(List<TrieNode> * lzw, int tableSize) {
    memset(lzw->data, 0xFF, 4096 * sizeof(TrieNode));
    lzw->len = tableSize + 2;
}

int choose_meta_palette(List<u16 *> cookedFrames, int width, int height) {
    bool used[4096];
    bool used2[4096];
    int maxUsed[6] = {};
    int cvtMasks[6] = { 0xFFF, 0xEFF, 0xEFE, 0xEEE, 0xCEE, 0xCEC }; //favor green > red > blue
    // int cvtMasks[6] = { 0xFFF, 0xEFF, 0xEEF, 0xEEE, 0xCEE, 0xCCE }; //favor red > green > blue
    //count how many colors are used out of each of a number of possible meta-palettes
    for (int i : range(1, cookedFrames.len)) {
        u16 * frame = cookedFrames[i];
        memset(used, 0, 4096 * sizeof(bool));
        for (int i : range(width * height))
            used[frame[i]] = true;

        for (int m : range(6)) {
            memset(used2, 0, 4096 * sizeof(u8));
            for (int i : range(4096))
                used2[i & cvtMasks[m]] |= used[i];
            int count = 0;
            for (int i : range(4096))
                if (used2[i])
                    ++count;
            printf("used %3d colors out of %4d\t\t", count, 1 << (12 - m));
            maxUsed[m] = imax(maxUsed[m], count);
        }
        printf("\n");
    }

    for (int m : range(6)) {
        printf("used at most %3d colors out of %4d\n", maxUsed[m], 1 << (12 - m));
    }

    for (int m : range(6)) {
        if (maxUsed[m] < 256) {
            return cvtMasks[m];
        }
    }

    return 0; //this should never happen
}

void save_gif(int width, int height, List<RawFrame> rawFrames, int centiSeconds) {
    float preAmble = get_time();

    //cook frames (downsample to 12-bit color)
    List<u16 *> cookedFrames = create_list<u16 *>(rawFrames.len);
    cookedFrames.add((u16 *) malloc(width * height * sizeof(u16))); //dummy frame for diff base
    memset(cookedFrames[0], 0, width * height * sizeof(u16)); //set dummy frame to background color
    for (RawFrame frame : rawFrames) {
        u16 * data = (u16 *) malloc(width * height * sizeof(u16));
        for (int y : range(height)) {
            for (int x : range(width)) {
                Pixel p = frame.pixels[y * frame.pitch + x];
                data[y * width + x] = (p.b & 0xF0) << 4 | (p.g & 0xF0) | (p.r & 0xF0) >> 4;
            }
        }
        cookedFrames.add(data);
        free(frame.base);
    }

    //season the frames (apply mask)
    float preChoice = get_time();
    int cvtMask = choose_meta_palette(cookedFrames, width, height);
    printf("choice: %fs\n", get_time() - preChoice);
    printf("conversion mask: %X\n", cvtMask);
    for (u16 * frame : cookedFrames) {
        for (int i : range(width * height)) {
            frame[i] &= cvtMask;
        }
    }

    printf("preAmble: %fs\n", get_time() - preAmble);

    //header
    FileBuffer buf = create_file_buffer(2048);
    for (char c : range("GIF89a")) {
        buf.write_unsafe(c);
    }

    //logical screen descriptor
    buf.write_unsafe<u16>(width);
    buf.write_unsafe<u16>(height);
    //global color table flag, color resolution (???), sort flag, global color table size
    buf.write_unsafe<u8>(0b0'001'0'000);
    buf.write_unsafe<u8>(0); //background color index
    buf.write_unsafe<u8>(0); //pixel aspect ratio

    //application extension
    buf.write_unsafe<u8>(0x21); //extension introducer
    buf.write_unsafe<u8>(0xFF); //extension identifier
    buf.write_unsafe<u8>(11); //fixed length data size
    for (char c : range("NETSCAPE2.0")) {
        buf.write_unsafe(c);
    }
    buf.write_unsafe<u8>(3); //data block size
    buf.write_unsafe<u8>(1); //???
    buf.write_unsafe<u16>(0); //loop forever
    buf.write_unsafe<u8>(0); //block terminator

    List<TrieNode> lzw = create_list<TrieNode>(4096);
    List<u8> idxBuffer = create_list<u8>(200);
    uint largestIdxBuffer = 0; //DEBUG

    float paletteTotal = 0;
    for (int i : range(1, cookedFrames.len)) {
        u16 * pframe = cookedFrames[i - 1];
        u16 * cframe = cookedFrames[i];
        // printf("\n\n\n\nnew frame\n\n\n\n");

        float prePalette = get_time();
        //generate palette
        u8 tlb[4096] = {};
        struct Color3 { u8 r, g, b; };
        Color3 table[256] = {};
        int tableIdx = 1; //we start counting at 1 because 0 is the transparent color
        for (int i : range(width * height)) {
            if (!tlb[cframe[i]]) {
                tlb[cframe[i]] = tableIdx;
                table[tableIdx] = {
                    (u8)((cframe[i] & 0x00F) << 4 | (cframe[i] & 0x00F)     ),
                    (u8)((cframe[i] & 0x0F0)      | (cframe[i] & 0x0F0) >> 4),
                    (u8)((cframe[i] & 0xF00) >> 4 | (cframe[i] & 0xF00) >> 8),
                };
                ++tableIdx;
            }
        }
        paletteTotal += get_time() - prePalette;

        int tableBits = _bit_scan_reverse(tableIdx - 1) + 1;
        int tableSize = 1 << tableBits;
        // printf("idx: %d bits: %d size: %d\n\n\n\n", tableIdx, tableBits, tableSize);

        buf.check(8 + 10);
        //graphics control extension
        buf.write_unsafe<u8>(0x21); //extension introducer
        buf.write_unsafe<u8>(0xF9); //extension identifier
        buf.write_unsafe<u8>(4); //block size (always 4)
        //reserved, disposal method:keep, input flag, transparency flag
        buf.write_unsafe<u8>(0b000'001'0'0 | (i != 1));
        buf.write_unsafe<u16>(centiSeconds); //x/100 seconds per frame
        buf.write_unsafe<u8>(0); //transparent color index
        buf.write_unsafe<u8>(0); //block terminator

        //image descriptor
        buf.write_unsafe<u8>(0x2C); //image separator
        buf.write_unsafe<u16>(0); //image left
        buf.write_unsafe<u16>(0); //image top
        buf.write_unsafe<u16>(width);
        buf.write_unsafe<u16>(height);
        //local color table flag, interlace flag, sort flag, reserved, local color table size
        buf.write_unsafe<u8>(0b1'0'0'00'000 | (tableBits - 1));

        //local color table
        buf.write_block(table, tableSize);

        //image data
        BlockBuffer block = {};
        buf.write<u8>(tableBits); //lzw minimum code size
        reset(&lzw, tableSize);
        //XXX: do we actually need to write this?
        put_code(&buf, &block, _bit_scan_reverse(lzw.len - 1) + 1, tableSize); //clear code

        int lastCode = cframe[0] == pframe[0]? 0 : tlb[cframe[0]];
        for (int i : range(1, width * height)) {
            idxBuffer.add(cframe[i] == pframe[i]? 0 : tlb[cframe[i]]);
            int code = lzw[lastCode].next[idxBuffer[idxBuffer.len - 1]];
            if (code < 0) {
                //write to code stream
                int codeBits = _bit_scan_reverse(lzw.len - 1) + 1;
                put_code(&buf, &block, codeBits, lastCode);
                // printf("%d-%d-%d  ", lastCode, codeBits, (int) lzw.len);

                //NOTE: [I THINK] we need to leave room for 2 more codes (leftover and end code)
                //      because we don't ever reset the table after writing the leftover bits
                //XXX: is my thinking correct on this one?
                if (lzw.len > 4094) {
                    //reset buffer code table
                    put_code(&buf, &block, codeBits, tableSize);
                    reset(&lzw, tableSize);
                } else {
                    lzw[lastCode].next[idxBuffer[idxBuffer.len - 1]] = lzw.len;
                    ++lzw.len;
                }

                //reset index buffer
                idxBuffer[0] = idxBuffer[idxBuffer.len - 1];
                idxBuffer.len = 1;

                lastCode = idxBuffer[0];
            } else {
                lastCode = code;
            }

            //DEBUG
            if (idxBuffer.len > largestIdxBuffer)
                largestIdxBuffer = idxBuffer.len;
        }

        //write code for leftover index buffer contents, then the end code
        put_code(&buf, &block, _bit_scan_reverse(lzw.len - 1) + 1, lastCode);
        put_code(&buf, &block, _bit_scan_reverse(lzw.len) + 1, tableSize + 1); //end code

        //flush remaining data
        if (block.bits) {
            int bytes = (block.bits + 7) / 8; //round up
            buf.check(bytes + 1);
            buf.write_unsafe<u8>(bytes);
            buf.write_block_unsafe(block.bytes, bytes);
        }

        buf.write<u8>(0); //terminating block
        idxBuffer.len = 0; //reset encoding state
    }

    printf("palette time: %fs\n", paletteTotal); //DEBUG
    printf("largest idx buffer: %d\n", largestIdxBuffer); //DEBUG

    buf.write<u8>(0x3B); //trailing marker

    //write data to file
    FILE * fp = fopen("out.gif", "wb");
    assert(fp);
    fwrite(buf.block, buf.size(), 1, fp);
    fclose(fp);

    //cleanup
    buf.finalize();
    lzw.finalize();
    idxBuffer.finalize();
    for (u16 * frame : cookedFrames)
        free(frame);
    cookedFrames.finalize();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

struct DrawCall {
    Model * model;
    Mat4 mat;
    bool shadow;
};

//SDL defines main to SDL_Main, which causes errors when linking manually on Windows
#ifdef _WIN32
# undef main
#endif

// #include <mach-o/dyld.h>
#include <unistd.h>

int main() {
#ifdef _WIN32
    bool success = load_sdl_functions("link/SDL2.dll");
    if (!success) {
        printf("exiting application because we couldn't load SDL dynamically\n");
        exit(1);
    }
#else
    // if (chdir("/Users/stuntddude/Dropbox/c/DIWide"))
    //     printf("        !!!chdir failed!!!\n");
#endif

    //initialize timer
    applicationStartupTimeValue = SDL_GetPerformanceCounter();

    if (SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL FAILED TO INIT: %s\n", SDL_GetError());
        return 1;
    }
    printf("SDL init: %f seconds\n", get_time());

    const int gameWidth = 192*2;
    const int gameHeight = 120*2;
    const int gameScale = 4;
    assert(gameWidth % 8 == 0);

    SDL_Window * window = SDL_CreateWindow("Test Window",
        SDL_WINDOWPOS_CENTERED_DISPLAY(1), SDL_WINDOWPOS_CENTERED_DISPLAY(1),
        gameWidth * gameScale, gameHeight * gameScale, SDL_WINDOW_SHOWN);
    if (!window) {
        printf("SDL FAILED TO CREATE WINDOW: %s\n", SDL_GetError());
        return 2;
    }
    printf("SDL create window: %f seconds\n", get_time());

    Canvas canv = create_canvas(gameWidth, gameHeight, 0);
    ZBuffer depth = create_depth_buffer(1024, 1024);
    Canvas * canvas = &canv;
    ZBuffer * shadow = &depth;
    Model * sword = load_model("res/sword7.diw");
    sword->tex = load_texture("res/sword7_c.png", "res/sword7_n.png");
    Model * sword2 = load_model("res/sword9.diw");
    sword2->tex = load_texture("res/sword9_c.png", "res/sword9_n.png");
    Model * tomato = load_model("res/tomato.diw");
    tomato->tex = load_texture("res/tomato_c.png", "res/tomato_n.png");

    Vertex _vertices[4] = {
        { { 0, 0, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },   -2, -2 },
        { { 0, 1, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },   -2,  2 },
        { { 1, 0, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },    2, -2 },
        { { 1, 1, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },    2,  2 },
    };
    Triangle _triangles[2] = { { 0, 3, 1 }, { 0, 2, 3 } };
    Model test = { 4, _vertices, (Vert *) malloc(4 * sizeof(Vert)), 2, _triangles, // sword->tex };
            load_texture("res/metal1_c.png", "res/metal1_n.png") };
    Model floor = { 4, _vertices, (Vert *) malloc(4 * sizeof(Vert)), 2, _triangles,
            load_texture("res/cobblestone1_c.png", "res/cobblestone1_n.png") };
            // load_texture("res/tiles1_c.png", "res/tiles1_n.png") };

    MonoFont mono = load_mono_font("3x6-bw.png", 16, 8);
    printf("SDL full init: %f seconds\n", get_time());



    //keymap
    bool isDown[256] = {};
    bool captured = false;
    SDL_SetRelativeMouseMode((SDL_bool)captured); //yo this cast sucks shit. fuck you, SDL

    //camera
    Vec3 camPos = vec3(-2, 1.5, 7);
    float cameraYaw = 0.6;
    float cameraPitch = 0.6;

    //timestep and framerate info
    float frameTimes[100] = {};
    float time = get_time();
    float lastTime = 0;
    int frame = 0;

    List<DrawCall> drawList = {};

    bool recordingGif = false;
    float gifTimer = 0;
    int gifFps = 20;
    List<RawFrame> gifFrames = {};

    u64 perfRasterLowest = 10000000000;
    bool shouldExit = false;
    while (!shouldExit) {
        float dmx = 0, dmy = 0;
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                shouldExit = true;
            } else if (event.type == SDL_KEYDOWN || event.type == SDL_KEYUP) {
                isDown[event.key.keysym.scancode] = event.type == SDL_KEYDOWN;
                if (event.type == SDL_KEYDOWN && event.key.keysym.scancode == SDL_SCANCODE_G) {
                    int centiSeconds = 100/gifFps;
                    if (recordingGif) {
                        printf("Saving GIF...\n");
                        float preGif = get_time();
                        save_gif(gameWidth, gameHeight, gifFrames, centiSeconds);
                        printf("%fs\n", get_time() - preGif);
                        gifFrames.len = 0;
                    }
                    gifTimer = 0;
                    recordingGif = !recordingGif;
                }
            } else if (event.type == SDL_MOUSEMOTION) {
                dmx += event.motion.xrel;
                dmy += event.motion.yrel;
            }
        }

        //update timestep
        lastTime = time;
        time = get_time();
        float dt = time - lastTime;

        //camera controls
        if (isDown[SDL_SCANCODE_ESCAPE] || isDown[SDL_SCANCODE_GRAVE]) {
            captured = isDown[SDL_SCANCODE_GRAVE];
            SDL_SetRelativeMouseMode((SDL_bool)captured);
        }

        if (captured) {
            cameraYaw += dmx * 0.01f;
            cameraPitch += dmy * 0.01f;

            //lock camera pitch to within reasonable bounds
            if (cameraPitch > radians(80)) {
                cameraPitch = radians(80);
            } else if (cameraPitch < radians(-80)) {
                cameraPitch = radians(-80);
            }

            //generate 3D unit vector representing camera facing direction
            Vec3 facing = normalize({ sinf(cameraYaw), -sinf(cameraPitch), -cosf(cameraYaw) });

            float sp = dt * 10;
            if (isDown[SDL_SCANCODE_LSHIFT]) {
                sp *= 0.25f;
            }

            if (isDown[SDL_SCANCODE_A]) {
                camPos += vec3(facing.z, 0, -facing.x) * sp;
            } if (isDown[SDL_SCANCODE_D]) {
                camPos += vec3(-facing.z, 0, facing.x) * sp;
            } if (isDown[SDL_SCANCODE_W]) {
                camPos += facing * sp;
            } if (isDown[SDL_SCANCODE_S]) {
                camPos += facing * -sp;
            } if (isDown[SDL_SCANCODE_E]) {
                camPos.y += sp;
            } if (isDown[SDL_SCANCODE_Q]) {
                camPos.y -= sp;
            }
        }



        u64 preFill = perf();
        //clear the screen
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * row = canvas->pixels + y * canvas->pitch;
            float * zrow = canvas->depth + y * canvas->zpitch;
            for (int x = 0; x < canvas->width; ++x) {
                row[x] = { 191, 127, 0, 255 };
                zrow[x] = 1.0f;
            }
        }
        perfFill = perf() - preFill;



        //update sliding window filter for framerate
        float timeSum = dt;
        for (int i : range(1, ARR_SIZE(frameTimes))) {
            frameTimes[i - 1] = frameTimes[i];
            timeSum += frameTimes[i - 1];
        }
        frameTimes[ARR_SIZE(frameTimes) - 1] = dt;

        float framerate = ARR_SIZE(frameTimes) / timeSum;

#if 1
        u64 preText = perf();
        char buffer[250];
        sprintf(buffer, "total:        %5d   drawn:      %5d", totalTris, drawnTris);
        draw_text(canvas, &mono, 4, 4, buffer);
        sprintf(buffer, "shadow total: %5d   shadow:     %5d", shadowTotalTris, shadowDrawnTris);
        draw_text(canvas, &mono, 4, 12, buffer);
        sprintf(buffer, "clipped:      %5d   draw calls: %5d", clippedTris, (int) drawList.len);
        draw_text(canvas, &mono, 4, 20, buffer);
        sprintf(buffer, "fps: %3d", (int)(framerate + 0.5f));
        draw_text(canvas, &mono, canvas->width - strlen(buffer) * mono.glyphWidth - 4, 4, buffer);

        sprintf(buffer, "fill:      %6lldk", perfFill      / 1024);
        draw_text(canvas, &mono, 4, 36, buffer);
        sprintf(buffer, "text:      %6lldk", perfText      / 1024);
        draw_text(canvas, &mono, 4, 44, buffer);
        sprintf(buffer, "shadow:    %6lldk", perfShadow    / 1024);
        draw_text(canvas, &mono, 4, 52, buffer);
        sprintf(buffer, "draw:      %6lldk", perfDraw      / 1024);
        draw_text(canvas, &mono, 4, 60, buffer);
        sprintf(buffer, "blit:      %6lldk", perfBlit      / 1024);
        draw_text(canvas, &mono, 4, 68, buffer);

        sprintf(buffer, "transform: %6lldk", perfTransform / 1024);
        draw_text(canvas, &mono, 4, 84, buffer);
        sprintf(buffer, "rasterize: %6lldk", perfRasterize / 1024);
        draw_text(canvas, &mono, 4, 92, buffer);
        sprintf(buffer, "inner:     %6lldk", perfInner     / 1024);
        draw_text(canvas, &mono, 4, 100, buffer);

        sprintf(buffer, "lowest:    %6lldk", perfRasterLowest / 1024);
        draw_text(canvas, &mono, 4, 116, buffer);

        if (recordingGif) {
            draw_text(canvas, &mono, canvas->width - strlen("GIF") * mono.glyphWidth - 4,
                                     canvas->height - mono.glyphHeight - 4, "GIF");
        }
        perfText = perf() - preText;
#endif

        totalTris = drawnTris = clippedTris = shadowTotalTris = shadowDrawnTris = 0;
        perfTransform = perfRasterize = perfInner = 0;



        drawList.len = 0;
        Vec3 lightPos = vec3(-1, sinf(get_time() * 1.4f), 0);
        Vec3 lightDir = noz(vec3(sin(get_time() * 0.5f), 0.5f, cos(get_time() * 0.5f)));

        {
            Mat4 mat = translate(scale(IDENTITY_4, 16), vec3(0, -2, 4));
            drawList.add({ tomato, mat, true });
        } {
            Mat4 mat = rotateZ(scale(IDENTITY_4, 2), -HALF_PI);
            mat = translate(mat, vec3(-4, 2, -10));
            drawList.add({ sword, mat, true });
        } {
            Mat4 mat = translate(rotateY(scale(IDENTITY_4, 8), -get_time()*0.0f), vec3(4, 0, 0));
            drawList.add({ &test, mat, false });
        } {
            Mat4 mat = translate(IDENTITY_4, vec3(-0.5f, -0.5f, 0));
            mat = translate(rotateX(scale(mat, 64), -HALF_PI), vec3(0, -2, 0));
            drawList.add({ &floor, mat, false });
        } {
            Mat4 mat = translate(scale(IDENTITY_4, 2), lightPos);
            drawList.add({ tomato, mat, true });
        }



        shadow->scale = 12;
        shadow->inv = 1.0f / shadow->scale;
        Mat4 shadowMat = scale(look_at(vec3(), -lightDir, vec3(0, 1, 0)),
                               shadow->inv, shadow->inv, shadow->inv);

        u64 preShadow = perf();
        //clear the shadow buffer
        for (int y = 0; y < shadow->height; ++y) {
            float * zrow = shadow->depth + y * shadow->width;
            for (int x = 0; x < shadow->width; ++x) {
                zrow[x] = -1000;
            }
        }

        Mat4 newShadowMat = shadowMat;
        newShadowMat = translate(shadowMat, vec3(0, 0, -shadow->inv * 1.0f));
        for (DrawCall call : drawList) {
            if (call.shadow) {
                draw_shadow(shadow, call.model, call.mat, newShadowMat);
            }
        }
        perfShadow = perf() - preShadow;

        Mat4 view = translate(IDENTITY_4, vec3(-camPos.x, -camPos.y, -camPos.z));
        view = rotateY(view, cameraYaw);
        view = rotateX(view, cameraPitch);
        Mat4 proj = perspective(radians(60), (float) canvas->width / canvas->height, 1, 100);

        u64 preDraw = perf();
        for (DrawCall call : drawList) {
            draw_model(canvas, shadow,
                       call.model, call.mat, view, proj, shadowMat,
                       camPos, lightDir);
        }
        perfDraw = perf() - preDraw;



        //update lowest
        if (perfRasterize < perfRasterLowest) {
            perfRasterLowest = perfRasterize;
        }



        //upscale into the window's frame buffer
        u64 preBlit = perf();
        fast_scaled_blit(SDL_GetWindowSurface(window), canvas, gameScale);
        // debug_depth_blit(SDL_GetWindowSurface(window), shadow);
        perfBlit = perf() - preBlit;
        SDL_UpdateWindowSurface(window);



        //set aside frames for GIF
        if (recordingGif) {
            float gifFrameTime = 1.0f / gifFps;
            gifTimer += dt;
            if (gifTimer > gifFrameTime) {
                gifTimer -= gifFrameTime;
                printf("Setting aside GIF frame %d...\n", (int)(gifFrames.len + 1));
                //we steal the canvas's memory block because it's much faster than copying the data
                gifFrames.add({ canvas->base, canvas->pixels, canvas->pitch });
                size_t offset = canvas->pixels - canvas->base;
                canvas->base = (Pixel *) malloc(canvas->pixelBytes);
                canvas->pixels = canvas->base + offset;

                int centiSeconds = 100/gifFps;
                if ((int) gifFrames.len == 125 * gifFps/10) {
                    printf("Saving GIF...\n");
                    float preGif = get_time();
                    save_gif(gameWidth, gameHeight, gifFrames, centiSeconds);
                    printf("%fs\n", get_time() - preGif);
                    gifFrames.len = 0;
                    gifTimer = 0;
                    recordingGif = false;
                }
            }
        }



        fflush(stdout);
        fflush(stderr);
        frame += 1;
        SDL_Delay(12);

        //uncomment this to make the game exit immediately (good for testing compile+load times)
        // shouldExit = true;
    }

    printf("time alive: %f seconds\n", get_time());

    return 0;
}
