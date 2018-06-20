//global TODO:
//- sorting the render queue
//- point diffuse lighting
//- metalness map

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
//- text overlay for stats
//- more specific performance counters
//- overdraw visualization
//- stats about wasted SIMD lanes

//- vsync?
//- GIF exporter!

#include "common.hpp"
#include "blit.hpp"
#include "math.hpp"
#include "List.hpp"

#include "stb_image.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <SDL.h>

u64 applicationStartupTimeValue;

double get_time() {
    u64 currentTimeValue = SDL_GetPerformanceCounter();
    u64 diffTimeValue = currentTimeValue - applicationStartupTimeValue;
    double elapsedSeconds = (double)diffTimeValue / (double)SDL_GetPerformanceFrequency();
    return elapsedSeconds;
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
///                                                                                              ///
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
///                                                                                              ///
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
///                                                                                              ///
////////////////////////////////////////////////////////////////////////////////////////////////////

//DEBUG GLOBALS
int totalTris;
int drawnTris;
int clippedTris;
int shadowTotalTris;
int shadowDrawnTris;

Vec4 color_sq(Vec4 v) {
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
            u *= tex->width;
            v *= tex->height;
            float uf = u - (int)u;
            float vf = v - (int)v;
            int iu1 = (int)u & tex->wmask;
            int iv1 = (int)v & tex->hmask;
            int iu2 = (int)u + 1 & tex->wmask;
            int iv2 = (int)v + 1 & tex->hmask;

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
            Vec3 n11 = vec3(m11.r, m11.g, m11.b);
            Vec3 n12 = vec3(m12.r, m12.g, m12.b);
            Vec3 n21 = vec3(m21.r, m21.g, m21.b);
            Vec3 n22 = vec3(m22.r, m22.g, m22.b);
            Vec3 n1 = n11 * (1 - uf) + n12 * uf;
            Vec3 n2 = n21 * (1 - uf) + n22 * uf;
            Vec3 normal = n1 * (1 - vf) + n2 * vf;

            normal = nor(normal - vec3(127.5f));
            l = nor(l);
            c = nor(c);

            float shad = 1;
            if (sh.x > 0 && sh.y > 0 && sh.x < shadow->width - 1 && sh.y < shadow->height - 1) {
                //shadow sample calculations
                float sf = sh.x - (int)sh.x;
                float tf = sh.y - (int)sh.y;
                int is1 = (int)sh.x & shadow->wmask;
                int it1 = (int)sh.y & shadow->hmask;
                int is2 = (int)sh.x + 1 & shadow->wmask;
                int it2 = (int)sh.y + 1 & shadow->hmask;

                //sample diffuse texture
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
            color *= light * shad + 0.01f;

            //directional specular
            float e = 1 / (rough + 0.0001f);
            float m = 1 / (rough + 0.02f); //multiplier
            Vec3 reflected = 2 * dot(l, normal) * normal - l;
            float spec = fmax(0, dot(reflected, c));
            //specular = M * S / (E - E * S + S)
            float specular = m * spec / (e - e * spec + spec) + 0.02f;

            float ior = 1.333; //ior of water, chosen for no particular reason
            float f = sq((1 - ior) / (1 + ior));
            float headon = fmax(0, dot(c, normal));
            //fresnel = F + (1 - R) * (1 - F) * sq(sq(1 - C)) * (1 - C)
            float fresnel = f + sq(1 - roughness) * (1 - f) * sq(sq(1 - headon)) * (1 - headon);

            Vec3 highlight = vec3(191, 127, 0) * (1.0f / 255) * 0.1f * fmax(0, n.y);
            highlight += vec3(specular * shad);
            highlight *= fresnel;
            color *= 1 - fresnel;

            row[x] = { (u8)(fmin(1, sqrtf(color.x + highlight.x)) * 255),
                       (u8)(fmin(1, sqrtf(color.y + highlight.y)) * 255),
                       (u8)(fmin(1, sqrtf(color.z + highlight.z)) * 255), 255 };

            zrow[x] = z;
        }
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

    draw_triangle(canv, shadow, tex, tri);
}

void draw_model(Canvas * canvas, ZBuffer * shadow, Model * model,
                Mat4 modelMat, Mat4 viewMat, Mat4 projMat, Mat4 shadowMat,
                Vec3 camPos, Vec3 lightDir)
{
    Mat4 viewproj = projMat * viewMat;
    Mat4 modelviewproj = viewproj * modelMat;
    Mat4 modelshadow = shadowMat * modelMat;
    Mat3 normalMat = inverse_transpose(mat3(modelMat));

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

        //NOTE: assumes uv will be > -8 to make rounding faster for texture fetch
        //      (floor is slower than round-toward-zero, even in SIMD)
        v.u += 8 - 0.5f / model->tex.width;
        v.v += 8 - 0.5f / model->tex.width;

        Vec3 n = noz(normalMat * v.n);

        model->verts[i] = { p, n * p.w, l * p.w, c * p.w, s * p.w, v.u * p.w, v.v * p.w };
    }

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

struct DrawCall {
    Model * model;
    Mat4 mat;
    bool shadow;
};

int main() {
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
        { { 0, 0, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },   -1, -1 },
        { { 0, 1, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },   -1,  1 },
        { { 1, 0, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },    1, -1 },
        { { 1, 1, 0 },   { 0, 0, 1 },   { 1, 0, 0 },   { 0, 1, 0 },    1,  1 },
    };
    Triangle _triangles[2] = { { 0, 3, 1 }, { 0, 2, 3 } };
    Model test = { 4, _vertices, (Vert *) malloc(4 * sizeof(Vert)), 2, _triangles, sword2->tex };
    Model floor = { 4, _vertices, (Vert *) malloc(4 * sizeof(Vert)), 2, _triangles,
                     // load_texture("res/cobblestone1_c.png", "res/cobblestone1_n.png") };
                     load_texture("res/tiles1_c.png", "res/tiles1_n.png") };
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
    float frameTimes[10] = {};
    float time = get_time();
    float lastTime = 0;
    int frame = 0;

    bool shouldExit = false;
    while (!shouldExit) {
        float dmx = 0, dmy = 0;
        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                shouldExit = true;
            } else if (event.type == SDL_KEYDOWN || event.type == SDL_KEYUP) {
                isDown[event.key.keysym.scancode] = event.type == SDL_KEYDOWN;
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



        //clear the screen
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * row = canvas->pixels + y * canvas->pitch;
            float * zrow = canvas->depth + y * canvas->zpitch;
            for (int x = 0; x < canvas->width; ++x) {
                row[x] = { 191, 127, 0, 255 };
                zrow[x] = 1.0f;
            }
        }

        List<DrawCall> drawList = {};
        Vec3 lightPos = vec3(-1, sinf(get_time() * 1.4f), 0);
        Vec3 lightDir = noz(vec3(sin(get_time() * 0.5f), 0.5f, cos(get_time() * 0.5f)));

        {
            Mat4 mat = translate(scale(IDENTITY_4, 16), vec3(0, -2, 4));
            drawList.add({ tomato, mat, true });
        } {
            Mat4 mat = rotateZ(scale(IDENTITY_4, 2), -HALF_PI);
            mat = translate(mat, vec3(-4, 2, -10));
            drawList.add({ sword2, mat, true });
        } {
            Mat4 mat = translate(rotateY(scale(IDENTITY_4, 8), -get_time()*0.0f), vec3(4, 0, 0));
            drawList.add({ &test, mat, false });
        } {
            Mat4 mat = translate(IDENTITY_4, vec3(-0.5f, -0.5f, 0));
            mat = translate(rotateX(scale(mat, 64), -HALF_PI), vec3(0, -2, 0));
            drawList.add({ &floor, mat, false });
        } {
            Mat4 mat = translate(scale(IDENTITY_4, 2), lightPos);
            drawList.add({ tomato, mat, false });
        }



        shadow->scale = 16;
        shadow->inv = 1.0f / shadow->scale;
        Mat4 shadowMat = scale(look_at(vec3(), -lightDir, vec3(0, 1, 0)),
                               shadow->inv, shadow->inv, shadow->inv);

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

        Mat4 view = translate(IDENTITY_4, vec3(-camPos.x, -camPos.y, -camPos.z));
        view = rotateY(view, cameraYaw);
        view = rotateX(view, cameraPitch);
        Mat4 proj = perspective(radians(60), (float)canvas->width / (float)canvas->height, 1, 100);

        for (DrawCall call : drawList) {
            draw_model(canvas, shadow,
                       call.model, call.mat, view, proj, shadowMat,
                       camPos, lightDir);
        }

        drawList.finalize();



        //update sliding window filter for framerate
        float timeSum = dt;
        for (int i : range(1, ARR_SIZE(frameTimes))) {
            frameTimes[i - 1] = frameTimes[i];
            timeSum += frameTimes[i - 1];
        }
        frameTimes[ARR_SIZE(frameTimes) - 1] = dt;

        float framerate = ARR_SIZE(frameTimes) / timeSum;
        if (frame % 100 == 9) {
            printf("total: %6d   drawn: %6d   shadow total: %6d   shadow: %6d   clipped: %5d"
                   "   draw calls: %6d   fps: %3d\n",
                   totalTris, drawnTris, shadowTotalTris, shadowDrawnTris, clippedTris,
                   (int) /*drawList.len*/0, (int)(framerate + 0.5f));
        }
        totalTris = drawnTris = clippedTris = shadowTotalTris = shadowDrawnTris = 0;



        //upscale into the window's frame buffer
        fast_scaled_blit(SDL_GetWindowSurface(window), canvas, gameScale);
        // debug_depth_blit(SDL_GetWindowSurface(window), shadow);
        SDL_UpdateWindowSurface(window);

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
