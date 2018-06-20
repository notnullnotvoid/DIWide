//task list
//
//TODO: operator overloads for vector math
//TODO: better/more precise performance measurements
//TODO: build options
//TODO: put library files into external directory and symlink into projects
//
//TODO: perspective correctness
//TODO: global material database
//TODO: switch to view space for lighting calculations
//TODO: more physically accurate fresnel
//TODO: per-material reflectivity
//TODO: fresnel calculation based on IOR?
//TODO: deferred rendering?
//
//TODO: Windows build
//TODO: packaged Mac build
//
//TODO: basic player control
//TODO: fixed update tick
//TODO: fixed camera
//TODO: basic terrain pieces
//TODO: random terrain gen?
//
//TODO: SIMD rendering
//TODO: multithreaded rendering?

//NOTE(miles): I'm not sure if this is right! <SDL.h> alone doesn't compile for me on OSX,
//but I think that's what you're supposed to use, at least on other platforms.
//let me know if you have problems with this include
#ifdef __APPLE__
# include <SDL2/SDL.h>
#else
# include <SDL.h>
#endif

#include "system.hpp"
#include "std.hpp"
#include "math.hpp"
#include "ArrayList.hpp"

#include "lib/stb_image.h"

//funky stuffs
#include "x86intrin.h"
#include "immintrin.h"

template <typename TYPE>
inline void swap(TYPE & a, TYPE & b) {
    TYPE temp = a;
    a = b;
    b = temp;
}

//constants
const float PI = 3.1415926535;
const float TWO_PI = PI * 2;
const float HALF_PI = PI / 2;
const float QUARTER_PI = PI / 4;

inline float degrees(float radians) {
    return radians * (1.0f / PI * 180.0f);
}

inline float radians(float degrees) {
    return degrees * (1.0f / 180.0f * PI);
}

inline float sq(float f) {
    return f * f;
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

u64 applicationStartupTimeValue;

double get_time() {
    u64 currentTimeValue = SDL_GetPerformanceCounter();
    u64 diffTimeValue = currentTimeValue - applicationStartupTimeValue;
    double elapsedSeconds = (double)diffTimeValue / (double)SDL_GetPerformanceFrequency();
    return elapsedSeconds;
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

//debug
inline void print_matrix(Mat4 m) {
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m00, m.m01, m.m02, m.m03);
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m10, m.m11, m.m12, m.m13);
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m20, m.m21, m.m22, m.m23);
    printf("%8.3f%8.3f%8.3f%8.3f\n", m.m30, m.m31, m.m32, m.m33);
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

struct Vert {
    float x, y, z, w; //fully transformed position
    float u, v;
};

struct Triangle {
    u16 idx[3];
};

struct Model {
    size_t vertexCount;
    Vec3 * vertices;
    Vert * verts; //post-transform coordinates

    size_t triangleCount;
    Triangle * triangles;

    ArrayList<Triangle> clippedTriangles;
};

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

typedef struct Pixel {
    u8 r, g, b, a;
} Color;

// const int WIDTH = 2040/3; //85*8
// const int WIDTH = 256*4;
// const int HEIGHT = 192*4;
const int WIDTH = 80*12;
const int HEIGHT = 50*12;
const int SCALE = 1;

//DEBUG GLOBALS
int totalTris;
int drawnTris;
int globalFrame;

void transform_model(Model * model, Mat4 modelMatrix, Mat4 viewMatrix, Mat4 projectionMatrix) {
    Mat4 viewproj = mul(projectionMatrix, viewMatrix);
    Mat4 modelviewproj = mul(viewproj, modelMatrix);

    //transform positions
    for (int i = 0; i < model->vertexCount; ++i) {
        Vec4 a = vec4(model->vertices[i], 1);

        //matrix transforms ("vertex shader")
        Vec4 b = mul(modelviewproj, a);

        //perspective divide
        b.w = 1 / b.w;
        b.x *= b.w;
        b.y *= b.w;
        b.z *= b.w;

        //transform from clip space to screen space
        //TODO: incorporate this into the modelviewproj matrix
        //      (or should we calculate this during the rasterization phase?)
        b.x = (b.x + 1) *  0.5f * WIDTH;
        b.y = (b.y - 1) * -0.5f * HEIGHT;

        model->verts[i] = { b.x, b.y, b.z, b.w, a.x * 4 * b.w, a.y * 4 * b.w };
    }
}

void draw_static_model(Model * model,
    SDL_Surface * canvas, float * depthBuffer,
    Color * texData, int texWidth, int texHeight)
{
    //draw triangles
    for (int t = 0; t < model->triangleCount; ++t) {
        Triangle tri = model->triangles[t];
        Vert triangle[3] = {
            model->verts[tri.idx[0]],
            model->verts[tri.idx[1]],
            model->verts[tri.idx[2]],
        };

        ++totalTris;

        //alternative method: https://www.geeksforgeeks.org/orientation-3-ordered-points/
        // if ((triangle[1].y - triangle[0].y) * (triangle[2].x - triangle[1].x) -
        //     (triangle[1].x - triangle[0].x) * (triangle[2].y - triangle[1].y) <= 0) {
        //     continue;
        // }

        if (triangle[0].w < 0 || triangle[1].w < 0 || triangle[2].w < 0) {
            continue;
        }

        if (triangle[0].x < 0 && triangle[1].x < 0 && triangle[2].x < 0) {
            continue;
        }
        if (triangle[0].x > WIDTH && triangle[1].x > WIDTH && triangle[2].x > WIDTH) {
            continue;
        }

        if (triangle[0].y < 0 && triangle[1].y < 0 && triangle[2].y < 0) {
            continue;
        }
        if (triangle[0].y > HEIGHT && triangle[1].y > HEIGHT && triangle[2].y > HEIGHT) {
            continue;
        }

        if (triangle[0].z < -1 && triangle[1].z < -1 && triangle[2].z < -1) {
            continue;
        }
        if (triangle[0].z > 1 && triangle[1].z > 1 && triangle[2].z > 1) {
            continue;
        }

        ++drawnTris;

        struct Edge {
            Vert v1, v2; //sorted by y
            float yfactor; //factor for interpolating vertex attributes vertically
        };

        //URGENT TODO: keep normals attached to vertices!

        Edge edges[3];

        float miny = HEIGHT;
        float maxy = 0;

        //convert triangles to edges
        for (int i = 0; i < 3; ++i) {
            Vert v1 = triangle[i];
            Vert v2 = triangle[(i + 1) % 3];

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
        if (lastLine > HEIGHT - 1) lastLine = HEIGHT - 1;

        for (int y = firstLine; y <= lastLine; ++y) {
            Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
            float * zrow = depthBuffer + y * WIDTH;

            //the current pixel row will be within the vertical extend of only two
            //of the three edges at any time, so find those two and discard the third
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
            float w1 = f1a * e1.v1.w + f1b * e1.v2.w;
            float w2 = f2a * e2.v1.w + f2b * e2.v2.w;

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
            if (last > WIDTH - 1) last = WIDTH - 1;

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
                float u = w * (fa * u1 + fb * u2);
                float v = w * (fa * v1 + fb * v2);



                //"fragment shader"
                float r, g, b;

#if 1
                int iu = u * texWidth;
                int iv = v * texHeight;
                if (iu >= 0 && iv >= 0 && iu < texWidth && iv < texHeight) {
                    int idx = iv * texWidth + iu;
                    Color color = texData[idx];
                    r = sq(color.r / 256.0f);
                    g = sq(color.g / 256.0f);
                    b = sq(color.b / 256.0f);
                } else {
                    Color color = texData[0];
                    r = sq(color.r / 256.0f);
                    g = sq(color.g / 256.0f);
                    b = sq(color.b / 256.0f);
                }
#else
                int iu = u * texWidth;
                int iv = v * texHeight;
                int valid = iu >= 0 & iv >= 0 & iu < texWidth & iv < texHeight? 0xFFFFFFFF : 0;
                int candidate = iv * texWidth + iu;
                int idx = candidate & valid;
                Color color = texData[idx];
                r = sq(color.r / 256.0f);
                g = sq(color.g / 256.0f);
                b = sq(color.b / 256.0f);
#endif

                // r = g = b = z;
                // r = u;
                // g = v;
                // b = z;

                //apply gamma correction
                r = sqrtf(r);
                g = sqrtf(g);
                b = sqrtf(b);

                //clamp
                r = fmin(r, 1);
                g = fmin(g, 1);
                b = fmin(b, 1);

                row[x] = { (u8)(r * 255), (u8)(g * 255), (u8)(b * 255), row[x].a };
                zrow[x] = z;
            }
        }
    }
}

void draw_static_model_avx(Model * model,
    SDL_Surface * canvas, float * depthBuffer,
    Color * texData, int texWidth, int texHeight)
{
    //draw triangles
    for (int t = 0; t < model->triangleCount; ++t) {
        Triangle tri = model->triangles[t];
        Vert triangle[3] = {
            model->verts[tri.idx[0]],
            model->verts[tri.idx[1]],
            model->verts[tri.idx[2]],
        };

        ++totalTris;

        //alternative method: https://www.geeksforgeeks.org/orientation-3-ordered-points/
        // if ((triangle[1].y - triangle[0].y) * (triangle[2].x - triangle[1].x) -
        //     (triangle[1].x - triangle[0].x) * (triangle[2].y - triangle[1].y) <= 0) {
        //     continue;
        // }

        if (triangle[0].w < 0 || triangle[1].w < 0 || triangle[2].w < 0) {
            continue;
        }

        if (triangle[0].x < 0 && triangle[1].x < 0 && triangle[2].x < 0) {
            continue;
        }
        if (triangle[0].x > WIDTH && triangle[1].x > WIDTH && triangle[2].x > WIDTH) {
            continue;
        }

        if (triangle[0].y < 0 && triangle[1].y < 0 && triangle[2].y < 0) {
            continue;
        }
        if (triangle[0].y > HEIGHT && triangle[1].y > HEIGHT && triangle[2].y > HEIGHT) {
            continue;
        }

        if (triangle[0].z < -1 && triangle[1].z < -1 && triangle[2].z < -1) {
            continue;
        }
        if (triangle[0].z > 1 && triangle[1].z > 1 && triangle[2].z > 1) {
            continue;
        }

        ++drawnTris;

        struct Edge {
            Vert v1, v2; //sorted by y
            float yfactor; //factor for interpolating vertex attributes vertically
        };

        Edge edges[3];

        float miny = HEIGHT;
        float maxy = 0;

        //convert triangles to edges
        for (int i = 0; i < 3; ++i) {
            Vert v1 = triangle[i];
            Vert v2 = triangle[(i + 1) % 3];

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
        if (lastLine > HEIGHT - 1) lastLine = HEIGHT - 1;

        for (int y = firstLine; y <= lastLine; ++y) {
            Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
            float * zrow = depthBuffer + y * WIDTH;

            //the current pixel row will be within the vertical extend of only two
            //of the three edges at any time, so find those two and discard the third
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
            float w1 = f1a * e1.v1.w + f1b * e1.v2.w;
            float w2 = f2a * e2.v1.w + f2b * e2.v2.w;

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
            if (last > WIDTH - 1) last = WIDTH - 1;

            for (int x = first; x <= last; x += 8) {
                //generate inclusion mask for pixels
                uint inclusionMaskSource[8] = {};
                for (int i = 0; i < 8 && x + i <= last; ++i) {
                    inclusionMaskSource[i] = 0xFFFFFFFF;
                }
                __m256i inclusionMask = _mm256_lddqu_si256((__m256i *) inclusionMaskSource);

                //calculate horizontal interpolation factor for each pixel
                __m256 xx = _mm256_set_ps(7, 6, 5, 4, 3, 2, 1, 0);
                __m256 fx = _mm256_set1_ps(x);
                xx = _mm256_add_ps(xx, fx);
                __m256 fmaxx = _mm256_set1_ps(maxx);
                __m256 fa = _mm256_sub_ps(fmaxx, xx);
                __m256 fxfactor = _mm256_set1_ps(xfactor);
                fa = _mm256_mul_ps(fa, fxfactor);
                __m256 f1 = _mm256_set1_ps(1.0f);
                __m256 fb = _mm256_sub_ps(f1, fa);

                //interpolate z for each pixel
                __m256 fz1 = _mm256_mul_ps(_mm256_set1_ps(z1), fa);
                __m256 fz2 = _mm256_mul_ps(_mm256_set1_ps(z2), fb);
                __m256 z = _mm256_add_ps(fz1, fz2);

                //TODO: depth test masking & early-out

                //calculate reciprocal w component for each pixel
                __m256 fw1 = _mm256_mul_ps(_mm256_set1_ps(w1), fa);
                __m256 fw2 = _mm256_mul_ps(_mm256_set1_ps(w2), fb);
                __m256 w = _mm256_rcp_ps(_mm256_add_ps(fw1, fw2));
                // __m256 w = _mm256_div_ps(_mm256_set1_ps(1), _mm256_add_ps(fw1, fw2));

                //interpolate u for each pixel
                __m256 fu1 = _mm256_mul_ps(_mm256_set1_ps(u1), fa);
                __m256 fu2 = _mm256_mul_ps(_mm256_set1_ps(u2), fb);
                __m256 u = _mm256_mul_ps(_mm256_add_ps(fu1, fu2), w);

                //interpolate u for each pixel
                __m256 fv1 = _mm256_mul_ps(_mm256_set1_ps(v1), fa);
                __m256 fv2 = _mm256_mul_ps(_mm256_set1_ps(v2), fb);
                __m256 v = _mm256_mul_ps(_mm256_add_ps(fv1, fv2), w);

#if 0
                //texture fetch
                //convert uv coords to pixel coords
                __m256i iu = _mm256_cvtps_epi32(_mm256_mul_ps(u, _mm256_set1_ps(texWidth)));
                __m256i iv = _mm256_cvtps_epi32(_mm256_mul_ps(v, _mm256_set1_ps(texHeight)));
                //test whether pixel coords are inside texture
                __m256i t1 = _mm256_cmpgt_epi32(iu, _mm256_set1_epi32(-1));
                __m256i t2 = _mm256_cmpgt_epi32(iv, _mm256_set1_epi32(-1));
                __m256i t3 = _mm256_cmpgt_epi32(_mm256_set1_epi32(texWidth), iu);
                __m256i t4 = _mm256_cmpgt_epi32(_mm256_set1_epi32(texHeight), iv);
                //combine tests together with inclusion mask into a fetch mask
                //TODO: is ANDing with the inclusion mask actually faster, or slower?
                __m256i fetchMask = _mm256_and_si256(inclusionMask, t1);
                fetchMask = _mm256_and_si256(fetchMask, t2);
                fetchMask = _mm256_and_si256(fetchMask, t3);
                fetchMask = _mm256_and_si256(fetchMask, t4);
                //generate pixel indices
                __m256i idx = _mm256_add_epi32(_mm256_mullo_epi32(iv, _mm256_set1_epi32(texWidth)), iu);
                //masked texture fetch
                __m256i tex = _mm256_mask_i32gather_epi32(_mm256_set1_epi32(0xFFFF00FF), texData, idx, fetchMask, 4);

                //unpack texture
                __m256i mask8 = _mm256_set1_epi32(0x000000FF);
                __m256 tr = _mm256_cvtepi32_ps(_mm256_and_si256(tex, mask8));
                __m256 tg = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex, 8), mask8));
                __m256 tb = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex, 16), mask8));
#endif

#if 0
                //clip UV coords to within [0, 1]
                u = _mm256_sub_ps(u, _mm256_floor_ps(u));
                v = _mm256_sub_ps(v, _mm256_floor_ps(v));
                //calculate texture coords in pixel space
                __m256i iu = _mm256_cvtps_epi32(_mm256_mul_ps(u, _mm256_set1_ps(texWidth - 1)));
                __m256i iv = _mm256_cvtps_epi32(_mm256_mul_ps(v, _mm256_set1_ps(texHeight - 1)));
                //texture fetch
                __m256i idx = _mm256_add_epi32(_mm256_mullo_epi32(iv, _mm256_set1_epi32(texWidth)), iu);
                __m256i tex = _mm256_i32gather_epi32(texData, idx, 4);

                //unpack texture
                __m256i mask8 = _mm256_set1_epi32(0x000000FF);
                __m256 tr = _mm256_cvtepi32_ps(_mm256_and_si256(tex, mask8));
                __m256 tg = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex, 8), mask8));
                __m256 tb = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex, 16), mask8));
#endif

#if 1
                //clip UV coords to within [0, 1]
                u = _mm256_sub_ps(u, _mm256_floor_ps(u));
                v = _mm256_sub_ps(v, _mm256_floor_ps(v));
                //generate bilinear blend factors
                __m256 tu = _mm256_mul_ps(u, _mm256_set1_ps(texWidth - 1));
                __m256 tv = _mm256_mul_ps(v, _mm256_set1_ps(texHeight - 1));
                __m256 blendx2 = _mm256_sub_ps(tu, _mm256_floor_ps(tu));
                __m256 blendy2 = _mm256_sub_ps(tv, _mm256_floor_ps(tv));
                __m256 blendx1 = _mm256_sub_ps(_mm256_set1_ps(1), blendx2);
                __m256 blendy1 = _mm256_sub_ps(_mm256_set1_ps(1), blendy2);
                //calculate texture coords in pixel space
                __m256i iu1 = _mm256_cvtps_epi32(_mm256_floor_ps(tu));
                __m256i iv1 = _mm256_cvtps_epi32(_mm256_floor_ps(tv));
                __m256i iu2 = _mm256_add_epi32(iu1, _mm256_set1_epi32(1));
                __m256i iv2 = _mm256_add_epi32(iv1, _mm256_set1_epi32(1));
                //generate pixel indices
                __m256i idx11 = _mm256_add_epi32(_mm256_mullo_epi32(iv1, _mm256_set1_epi32(texWidth)), iu1);
                __m256i idx12 = _mm256_add_epi32(_mm256_mullo_epi32(iv1, _mm256_set1_epi32(texWidth)), iu2);
                __m256i idx21 = _mm256_add_epi32(_mm256_mullo_epi32(iv2, _mm256_set1_epi32(texWidth)), iu1);
                __m256i idx22 = _mm256_add_epi32(_mm256_mullo_epi32(iv2, _mm256_set1_epi32(texWidth)), iu2);
                //do texture fetches
                __m256i mask8 = _mm256_set1_epi32(0x000000FF);

                __m256i tex11 = _mm256_i32gather_epi32(texData, idx11, 4);
                __m256 tr11 = _mm256_cvtepi32_ps(_mm256_and_si256(tex11, mask8));
                __m256 tg11 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex11, 8), mask8));
                __m256 tb11 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex11, 16), mask8));

                __m256i tex12 = _mm256_i32gather_epi32(texData, idx12, 4);
                __m256 tr12 = _mm256_cvtepi32_ps(_mm256_and_si256(tex12, mask8));
                __m256 tg12 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex12, 8), mask8));
                __m256 tb12 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex12, 16), mask8));

                __m256i tex21 = _mm256_i32gather_epi32(texData, idx21, 4);
                __m256 tr21 = _mm256_cvtepi32_ps(_mm256_and_si256(tex21, mask8));
                __m256 tg21 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex21, 8), mask8));
                __m256 tb21 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex21, 16), mask8));

                __m256i tex22 = _mm256_i32gather_epi32(texData, idx22, 4);
                __m256 tr22 = _mm256_cvtepi32_ps(_mm256_and_si256(tex22, mask8));
                __m256 tg22 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex22, 8), mask8));
                __m256 tb22 = _mm256_cvtepi32_ps(_mm256_and_si256(_mm256_srli_epi32(tex22, 16), mask8));

                //blend texture samples together
                __m256 tr1 = _mm256_add_ps(_mm256_mul_ps(tr11, blendx1), _mm256_mul_ps(tr12, blendx2));
                __m256 tr2 = _mm256_add_ps(_mm256_mul_ps(tr21, blendx1), _mm256_mul_ps(tr22, blendx2));
                __m256 tr = _mm256_add_ps(_mm256_mul_ps(tr1, blendy1), _mm256_mul_ps(tr2, blendy2));

                __m256 tg1 = _mm256_add_ps(_mm256_mul_ps(tg11, blendx1), _mm256_mul_ps(tg12, blendx2));
                __m256 tg2 = _mm256_add_ps(_mm256_mul_ps(tg21, blendx1), _mm256_mul_ps(tg22, blendx2));
                __m256 tg = _mm256_add_ps(_mm256_mul_ps(tg1, blendy1), _mm256_mul_ps(tg2, blendy2));

                __m256 tb1 = _mm256_add_ps(_mm256_mul_ps(tb11, blendx1), _mm256_mul_ps(tb12, blendx2));
                __m256 tb2 = _mm256_add_ps(_mm256_mul_ps(tb21, blendx1), _mm256_mul_ps(tb22, blendx2));
                __m256 tb = _mm256_add_ps(_mm256_mul_ps(tb1, blendy1), _mm256_mul_ps(tb2, blendy2));

                // __m256 tr = blendx1;
                // __m256 tg = blendy1;
                // __m256 tb = blendx1;
#endif

                //map [0, 255] to [0, 1]
                __m256 fr255 = _mm256_set1_ps(1.0f / 255.0f);
                tr = _mm256_mul_ps(tr, fr255);
                tg = _mm256_mul_ps(tg, fr255);
                tb = _mm256_mul_ps(tb, fr255);

                //TODO: convert texture gamma to light-linear space

                //set color
                // __m256 r = _mm256_set1_ps(1.0f);
                // __m256 g = _mm256_set1_ps(0.0f);
                // __m256 b = _mm256_set1_ps(1.0f);
                __m256 r = tr;
                __m256 g = tg;
                __m256 b = tb;
                __m256 a = _mm256_set1_ps(1.0f);

                //TODO: convert light-linear to gamma space

                //TODO: clamp RGB values

                //convert to packed RGBA8
                //TODO: there may be a faster way to do this using swizzle stuff, look into it?
                __m256 f255 = _mm256_set1_ps(255.0f);
                r = _mm256_mul_ps(r, f255);
                g = _mm256_mul_ps(g, f255);
                b = _mm256_mul_ps(b, f255);
                a = _mm256_mul_ps(a, f255);
                __m256i ir = _mm256_cvtps_epi32(r);
                __m256i ig = _mm256_cvtps_epi32(g);
                __m256i ib = _mm256_cvtps_epi32(b);
                __m256i ia = _mm256_cvtps_epi32(a);
                ig = _mm256_slli_epi32(ig, 8);
                ib = _mm256_slli_epi32(ib, 16);
                ia = _mm256_slli_epi32(ia, 24);
                __m256i out = _mm256_or_si256(ir, ig);
                out = _mm256_or_si256(out, ib);
                out = _mm256_or_si256(out, ia);

                _mm256_maskstore_epi32((int *) (&row[x]), inclusionMask, out);

                //TODO: write to the depth buffer



                // //"fragment shader"
                // float r, g, b;

                // int iu = u * texWidth;
                // int iv = v * texHeight;
                // int valid = iu >= 0 & iv >= 0 & iu < texWidth & iv < texHeight? 0xFFFFFFFF : 0;
                // int candidate = iv * texWidth + iu;
                // int idx = candidate & valid;
                // Color color = texData[idx];
                // r = sq(color.r / 256.0f);
                // g = sq(color.g / 256.0f);
                // b = sq(color.b / 256.0f);

                // //apply gamma correction
                // r = sqrtf(r);
                // g = sqrtf(g);
                // b = sqrtf(b);

                // //clamp
                // r = fmin(r, 1);
                // g = fmin(g, 1);
                // b = fmin(b, 1);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////////////

struct DrawCall {
    Model * model;
    Mat4 matrix; //model-to-world transform
    bool shadow; //we will probably replace this with a flags mask later
};

int main() {
    //initialize timer
    applicationStartupTimeValue = SDL_GetPerformanceCounter();

    if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_JOYSTICK | SDL_INIT_GAMECONTROLLER)) {
        printf("SDL FAILED TO INIT: %s\n", SDL_GetError());
        return 1;
    }

    printf("SDL init: %f seconds\n", get_time());

    SDL_Window * window = SDL_CreateWindow("Test Window",
        SDL_WINDOWPOS_CENTERED_DISPLAY(1), SDL_WINDOWPOS_CENTERED_DISPLAY(1),
        WIDTH * SCALE, HEIGHT * SCALE,
        SDL_WINDOW_SHOWN);
    if (window == nullptr) {
        printf("SDL FAILED TO CREATE WINDOW: %s\n", SDL_GetError());
        return 1;
    }

    printf("SDL create window: %f seconds\n", get_time());

    //TODO: match window surface's byte order for faster blit?
    SDL_Surface * canvas = SDL_CreateRGBSurface(0, WIDTH, HEIGHT, 32,
        0x000000FF,
        0x0000FF00,
        0x00FF0000,
        0xFF000000);
    if (canvas == nullptr) {
        printf("SDL FAILED TO CREATE SURFACE: %s\n", SDL_GetError());
        return 1;
    }

    printf("SDL full init: %f seconds\n", get_time());

    Vec3 vertices[] = {
        { 0, 0, 0 },
        { 0, 1, 0 },
        { 1, 0, 0 },
        { 1, 1, 0 },
    };

    Triangle triangles[] = {
        { { 0, 1, 2 } },
        { { 1, 3, 2 } },
    };

    Model testModel = {};
    testModel.vertexCount = 4;
    testModel.vertices = vertices;
    testModel.verts = (Vert *) malloc(4 * sizeof(Vert));
    testModel.triangleCount = 2;
    testModel.triangles = triangles;
    testModel.clippedTriangles.init();

    //allocate buffers
    float * depthBuffer = (float *) malloc(WIDTH * HEIGHT * sizeof(float));



    //timestep and framerate info
    float frameTimes[10] = {};
    float time = get_time();
    float lastTime = 0;
    int frame = 0;



    //initialize render list
    ArrayList<DrawCall> drawList = create_array_list<DrawCall>(100);



    //keymap
    //TODO: make maps for tracking key events which are reset every frame
    //      (so that we can check for keydown and keyup events outside the event loop)
    bool isDown[256] = {};

    //input mode state
    bool captured = false;
    if (captured) {
        SDL_SetRelativeMouseMode(SDL_TRUE);
    }

    Vec3 camPos = vec3(0, 0, 5);
    float cameraYaw = 0;
    float cameraPitch = 0;



    //load image
    int x, y, n;
    u8 * data = stbi_load("res/biyori.png", &x, &y, &n, 4);

    printf("done initializing: %f seconds\n", get_time());

    const int PERF_WINDOW_SAMPLES = 50;
    u32 fillTimes[PERF_WINDOW_SAMPLES] = {};
    u32 transformTimes[PERF_WINDOW_SAMPLES] = {};
    u32 drawTimes[PERF_WINDOW_SAMPLES] = {};
    u32 blitTimes[PERF_WINDOW_SAMPLES] = {};







    bool shouldExit = false;
    while (!shouldExit) {
        float dmx = 0;
        float dmy = 0;

        SDL_Event event;
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_QUIT) {
                shouldExit = true;
            } else if (event.type == SDL_KEYDOWN) {
                isDown[event.key.keysym.scancode] = true;
            } else if (event.type == SDL_KEYUP) {
                isDown[event.key.keysym.scancode] = false;
            } else if (event.type == SDL_MOUSEMOTION) {
                dmx += event.motion.xrel;
                dmy += event.motion.yrel;
            }
        }

        //update timestep
        lastTime = time;
        time = get_time();
        float dt = time - lastTime;



        if (isDown[SDL_SCANCODE_ESCAPE]) {
            captured = false;
            SDL_SetRelativeMouseMode(SDL_FALSE);
        } else if (isDown[SDL_SCANCODE_GRAVE]) {
            captured = true;
            SDL_SetRelativeMouseMode(SDL_TRUE);
        }



        //camera controls
        if (captured) {
            //first-person debug cam
            cameraYaw += dmx * 0.01f;
            cameraPitch += dmy * 0.01f;

            //generate 3D unit vector representing camera facing direction
            Vec3 facing = normalize({ sinf(cameraYaw), -sinf(cameraPitch), -cosf(cameraYaw) });

            float sp = dt * 10;
            if (isDown[SDL_SCANCODE_LSHIFT]) {
                sp *= 0.25f;
            }

            if (isDown[SDL_SCANCODE_A]) {
                camPos.x += facing.z * sp;
                camPos.z -= facing.x * sp;
            }
            if (isDown[SDL_SCANCODE_D]) {
                camPos.x -= facing.z * sp;
                camPos.z += facing.x * sp;
            }
            if (isDown[SDL_SCANCODE_W]) {
                camPos.x += facing.x * sp;
                camPos.y += facing.y * sp;
                camPos.z += facing.z * sp;
            }
            if (isDown[SDL_SCANCODE_S]) {
                camPos.x -= facing.x * sp;
                camPos.y -= facing.y * sp;
                camPos.z -= facing.z * sp;
            }
            if (isDown[SDL_SCANCODE_UP]) { //if up is down, something is very wrong... or not
                camPos.y += sp;
            }
            if (isDown[SDL_SCANCODE_DOWN]) { //BUG: this should always return true, but doesn't
                camPos.y -= sp;
            }
        }

        //lock camera pitch to within reasonable bounds
        if (cameraPitch > radians(80)) {
            cameraPitch = radians(80);
        } else if (cameraPitch < radians(-80)) {
            cameraPitch = radians(-80);
        }


        float beforeFrame = get_time();

        SDL_LockSurface(canvas);

        //clear the screen
        u64 beforeFill = __rdtsc();
#if 0
        for (int y = 0; y < HEIGHT; ++y) {
            Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
            float * zrow = depthBuffer + y * WIDTH;
            for (int x = 0; x < WIDTH; ++x) {
                row[x] = { 63, 191, 191, 255 };
                // row[x] = { 0, 0, 0, 255 };
                zrow[x] = 1.0f;
                // row[x] = { (u8)x, (u8)y, (u8)(x * y), 255 }; //test pattern
            }
        }
#else
        assert(WIDTH % 8 == 0);
        __m256i colorFill = _mm256_set1_epi32(0xFFFFFF00);
        __m256 depthFill = _mm256_set1_ps(1.0f);
        for (int y = 0; y < HEIGHT; ++y) {
            Pixel * row = (Pixel *) ((u8 *)canvas->pixels + y * canvas->pitch);
            float * zrow = depthBuffer + y * WIDTH;
            for (int x = 0; x < WIDTH; x += 8) {
                _mm256_store_si256((__m256i *) (&row[x]), colorFill);
                _mm256_store_si256((__m256i *) (&zrow[x]), depthFill);
            }
        }
#endif
        u64 afterFill = __rdtsc();



        //clear the render queue
        drawList.len = 0;

        //place things in the render queue
        {
            Mat4 modelMatrix = IDENTITY_4;
            modelMatrix = scale(modelMatrix, 5);
            // modelMatrix = translate(modelMatrix, Vec4{ sinf(get_time() * 0.1f) * 10, 0, 0 });
            // modelMatrix = rotateX(modelMatrix, get_time());
            // modelMatrix = rotateY(modelMatrix, get_time() * 0.5f);
            // modelMatrix = rotateZ(modelMatrix, get_time() * 0.25f);

            drawList.add({ &testModel, modelMatrix, true });
        }



        //apply view transforms
        Mat4 view = IDENTITY_4;
        Mat4 projection = IDENTITY_4;

        view = translate(view, vec3(-camPos.x, -camPos.y, -camPos.z));
        view = rotateY(view, cameraYaw);
        view = rotateX(view, cameraPitch);

        projection = perspective(radians(60), (float)WIDTH / (float)HEIGHT, 1, 100);


        //render to frame buffer
        // float beforeTime = get_time();
        i64 transformCycles = 0;
        i64 drawCycles = 0;
        for (int i = 0; i < drawList.len; ++i) {
            i64 beforeCycles = __rdtsc();
            transform_model(drawList[i].model, drawList[i].matrix, view, projection);
            i64 afterCycles = __rdtsc();
            transformCycles += afterCycles - beforeCycles;

            beforeCycles = __rdtsc();
            // draw_static_model(drawList[i].model, canvas, depthBuffer, (Color *) data, x, y);
            draw_static_model_avx(drawList[i].model, canvas, depthBuffer, (Color *) data, x, y);
            afterCycles = __rdtsc();
            drawCycles += afterCycles - beforeCycles;
        }
        // float afterTime = get_time();

        // printf("%8.6f seconds\n", afterTime - beforeTime);
        // printf("%5lldk cycles\n", (afterCycles - beforeCycles)/1000);
        // printf("%8lld cycles\n", (afterCycles - beforeCycles)/1);



        SDL_UnlockSurface(canvas);



        //update sliding window filter for framerate
        float timeSum = 0;
        for (int i = 1; i < ARR_SIZE(frameTimes); ++i) {
            frameTimes[i - 1] = frameTimes[i];
            timeSum += frameTimes[i - 1];
        }
        frameTimes[ARR_SIZE(frameTimes) - 1] = dt;
        timeSum += dt;

        float framerate = ARR_SIZE(frameTimes) / timeSum;
        if (frame % 100 == 99) {
            // printf("%d fps\n", (int)(framerate + 0.5f));
            printf("total: %6d   drawn: %6d   draw calls: %6d   fps: %3d\n",
                   totalTris, drawnTris,
                   (int) drawList.len, (int)(framerate + 0.5f));
        }

        totalTris = drawnTris = 0;



        SDL_Surface * surface = SDL_GetWindowSurface(window);

        u64 beforeBlit = __rdtsc();
#if 0
        if (SCALE > 1) {
            SDL_BlitScaled(canvas, nullptr, surface, nullptr);
        } else {
            SDL_BlitSurface(canvas, nullptr, surface, nullptr);
        }
#else
        assert(WIDTH % 8 == 0);
        u8 shuffleMaskSource[32] = {
            2 +  0, 1 +  0, 0 +  0, 3 +  0,
            2 +  4, 1 +  4, 0 +  4, 3 +  4,
            2 +  8, 1 +  8, 0 +  8, 3 +  8,
            2 + 12, 1 + 12, 0 + 12, 3 + 12,
            2 +  0, 1 +  0, 0 +  0, 3 +  0,
            2 +  4, 1 +  4, 0 +  4, 3 +  4,
            2 +  8, 1 +  8, 0 +  8, 3 +  8,
            2 + 12, 1 + 12, 0 + 12, 3 + 12,
        };
        __m256i shuffleMask = _mm256_lddqu_si256((__m256i *) shuffleMaskSource);

        if (SCALE == 4) {
            for (int y = 0; y < HEIGHT; ++y) {
                Pixel * src = (Pixel *) ((u8 *) canvas->pixels + y * canvas->pitch);

                for (int x = 0; x < WIDTH; x += 8) {
                    //swizzle from ABGR to ARGB big endian byte order
                    __m256i chunk = _mm256_lddqu_si256((__m256i *) &src[x]);
                    chunk = _mm256_shuffle_epi8(chunk, shuffleMask);
                    //4x dupe pixels (notes follow, written in big endian order)
                    // 7777|6666   5555|4444   3333|2222   1111|0000
                    // 7777|3333   6666|2222   5555|1111   4444|0000
                    // 4.hi,3.hi   2.hi,1.hi   4.lo,3.lo   2.lo,1.lo
                    // src3,src4   src1,src2   src3,src4   src1,src2
                    // 0001 0011   0001 0011   0000 0010   0000 0010
                    __m256i shuf1 = _mm256_shuffle_epi32(chunk, 0b00000000);
                    __m256i shuf2 = _mm256_shuffle_epi32(chunk, 0b01010101);
                    __m256i shuf3 = _mm256_shuffle_epi32(chunk, 0b10101010);
                    __m256i shuf4 = _mm256_shuffle_epi32(chunk, 0b11111111);
                    __m256i perm1 = _mm256_permute2x128_si256(shuf1, shuf2, 0b00100000);
                    __m256i perm2 = _mm256_permute2x128_si256(shuf3, shuf4, 0b00100000);
                    __m256i perm3 = _mm256_permute2x128_si256(shuf1, shuf2, 0b00110001);
                    __m256i perm4 = _mm256_permute2x128_si256(shuf3, shuf4, 0b00110001);
                    //write pixels to memory
                    for (int i = 0; i < 4; ++i) {
                        u8 * dest = ((u8 *) surface->pixels + (4 * y + i) * surface->pitch);
                        _mm256_store_si256((__m256i *) &dest[x * 4 * 4 +  0 * 4], perm1);
                        _mm256_store_si256((__m256i *) &dest[x * 4 * 4 +  8 * 4], perm2);
                        _mm256_store_si256((__m256i *) &dest[x * 4 * 4 + 16 * 4], perm3);
                        _mm256_store_si256((__m256i *) &dest[x * 4 * 4 + 24 * 4], perm4);
                    }
                }
            }
        } else if (SCALE == 3) {
            for (int y = 0; y < HEIGHT; ++y) {
                Pixel * src = (Pixel *) ((u8 *) canvas->pixels + y * canvas->pitch);

                for (int x = 0; x < WIDTH; x += 8) {
                    //swizzle from ABGR to ARGB big endian byte order
                    __m256i chunk = _mm256_lddqu_si256((__m256i *) &src[x]);
                    chunk = _mm256_shuffle_epi8(chunk, shuffleMask);
                    //3x dupe pixels (notes follow, written in big endian order)
                    // desired             3:7776|6655   2:5444|3332   1:2211|1000
                    // uncombined          3:6655,2211   2:7776,3332   1:5444,1000
                    // shuffle mask (0b)   3:1010 0101   2:1111 1110   1:0100 0000
                    // combination mask    3:2.hi,3.hi   2:1.hi,2.lo   1:3.lo,1.lo
                    // a.k.a (arg)         3:src2,src3   2:src1,src2   1:src1,src3
                    // a.k.a (0b)          3:0001 0011   2:0001 0010   1:0010 0000
                    __m256i shuf1 = _mm256_shuffle_epi32(chunk, 0b01000000);
                    __m256i shuf2 = _mm256_shuffle_epi32(chunk, 0b11111110);
                    __m256i shuf3 = _mm256_shuffle_epi32(chunk, 0b10100101);
                    __m256i perm1 = _mm256_permute2x128_si256(shuf1, shuf3, 0b00100000);
                    __m256i perm2 = _mm256_permute2x128_si256(shuf1, shuf2, 0b00010010);
                    __m256i perm3 = _mm256_permute2x128_si256(shuf2, shuf3, 0b00010011);
                    //write pixels to memory
                    for (int i = 0; i < 3; ++i) {
                        u8 * dest = ((u8 *) surface->pixels + (3 * y + i) * surface->pitch);
                        _mm256_store_si256((__m256i *) &dest[x * 3 * 4 +  0 * 4], perm1);
                        _mm256_store_si256((__m256i *) &dest[x * 3 * 4 +  8 * 4], perm2);
                        _mm256_store_si256((__m256i *) &dest[x * 3 * 4 + 16 * 4], perm3);
                    }
                }
            }
        } else if (SCALE == 2) {
            for (int y = 0; y < HEIGHT; ++y) {
                Pixel * src = (Pixel *) ((u8 *) canvas->pixels + y * canvas->pitch);
                Pixel * dest1 = (Pixel *) ((u8 *) surface->pixels + (2 * y + 0) * surface->pitch);
                Pixel * dest2 = (Pixel *) ((u8 *) surface->pixels + (2 * y + 1) * surface->pitch);

                for (int x = 0; x < WIDTH; x += 8) {
                    __m256i chunk = _mm256_lddqu_si256((__m256i *) &src[x]);
                    chunk = _mm256_shuffle_epi8(chunk, shuffleMask);
                    __m256i srcLow = _mm256_unpacklo_epi32(chunk, chunk);
                    __m256i srcHigh = _mm256_unpackhi_epi32(chunk, chunk);
                    __m256i destLow = _mm256_permute2x128_si256(srcLow, srcHigh, 0x20);
                    __m256i destHigh = _mm256_permute2x128_si256(srcLow, srcHigh, 0x31);
                    _mm256_store_si256((__m256i *) &dest1[x * 2 + 0], destLow);
                    _mm256_store_si256((__m256i *) &dest1[x * 2 + 8], destHigh);
                    _mm256_store_si256((__m256i *) &dest2[x * 2 + 0], destLow);
                    _mm256_store_si256((__m256i *) &dest2[x * 2 + 8], destHigh);
                }
            }
        } else {
            for (int y = 0; y < HEIGHT; ++y) {
                Pixel * src = (Pixel *) ((u8 *) canvas->pixels + y * canvas->pitch);
                Pixel * dest = (Pixel *) ((u8 *) surface->pixels + y * surface->pitch);

                for (int x = 0; x < WIDTH; x += 8) {
                    __m256i chunk = _mm256_lddqu_si256((__m256i *) &src[x]);
                    chunk = _mm256_shuffle_epi8(chunk, shuffleMask);
                    _mm256_store_si256((__m256i *) &dest[x], chunk);
                }
            }
        }
#endif
        u64 afterBlit = __rdtsc();

        float afterFrame = get_time();

        //track performance measurements
        for (int i = PERF_WINDOW_SAMPLES - 1; i > 0; --i) {
            fillTimes[i] = fillTimes[i - 1];
            transformTimes[i] = transformTimes[i - 1];
            drawTimes[i] = drawTimes[i - 1];
            blitTimes[i] = blitTimes[i - 1];
        }
        fillTimes[0] = afterFill - beforeFill;
        transformTimes[0] = transformCycles;
        drawTimes[0] = drawCycles;
        blitTimes[0] = afterBlit - beforeBlit;
        //average performance measurements
        u64 fillSum = 0;
        u64 transformSum = 0;
        u64 drawSum = 0;
        u64 blitSum = 0;
        for (int i = 0; i < PERF_WINDOW_SAMPLES; ++i) {
            fillSum += fillTimes[i];
            transformSum += transformTimes[i];
            drawSum += drawTimes[i];
            blitSum += blitTimes[i];
        }
        printf("%5lldk fill   %5lldk transform   %5lldk draw   %5lldk blit   frame time %fs\n",
            fillSum / PERF_WINDOW_SAMPLES / 1000, transformSum / PERF_WINDOW_SAMPLES / 1000,
            drawSum / PERF_WINDOW_SAMPLES / 1000, blitSum / PERF_WINDOW_SAMPLES / 1000,
            afterFrame - beforeFrame);



        SDL_UpdateWindowSurface(window);

        fflush(stdout);
        fflush(stderr);

        frame += 1;

        globalFrame = frame;

        //uncomment this to make the game exit immediately (good for testing compile+load times)
        // shouldExit = true;
    }







    printf("time alive: %f seconds\n", get_time());

    //free resources (not strictly necessary)
    drawList.finalize();

    SDL_FreeSurface(canvas);
    SDL_DestroyWindow(window);
    SDL_Quit();

    printf("shutdown: %f seconds\n", get_time());

    return 0;
}
