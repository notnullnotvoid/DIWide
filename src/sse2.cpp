#include "common.hpp"
#include "math.hpp"
#include "blit.hpp"


#ifdef _MSC_VER
    #include "intrin.h"
#else
    #include "x86intrin.h"
#endif

struct Tex {
    Pixel * pixels;
    Pixel * normals;
    int width, height;
    u32 wmask, hmask;
};

extern int totalTris;
extern int drawnTris;
extern int clippedTris;
extern int shadowTotalTris;
extern int shadowDrawnTris;

extern u64 perfInner;
extern uint perfDummy;

inline u64 perf() {
    // return __rdtsc();
    return __rdtscp(&perfDummy);
}



struct Vert {
    Vec4 p; //position
    Vec3 n, l, c, s; //normal, light (directional), camera, shadow
    float u, v; //texture coords
};

#include "sse2.hpp"

inline __m128i gather_u32(u32 * addr, __m128i idx) {
    u32 idx11[4]; _mm_store_si128((__m128i *) idx11, idx);
    u32 ic11[4]; //masking this load only seems to make this slower... so we don't bother!
    for (int i = 0; i < 4; ++i) {
        ic11[i] = addr[idx11[i]];
    }
    return _mm_load_si128((__m128i *) ic11);
}

inline __m128 gather_f32(float * addr, __m128i idx) {
    u32 idx11[4]; _mm_store_si128((__m128i *) idx11, idx);
    float ic11[4]; //masking this load only seems to make this slower... so we don't bother!
    for (int i = 0; i < 4; ++i) {
        ic11[i] = addr[idx11[i]];
    }
    return _mm_load_ps(ic11);
}

void draw_triangle_sse2_plain(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]) {
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
        for (int x = first; x <= last; x += 4) {
            //generate inclusion mask to avoid overwriting pixels
            __m128i ix = _mm_add_epi32(_mm_set_epi32(3, 2, 1, 0), _mm_set1_epi32(x));
            __m128i inclusionMask = _mm_cmpgt_epi32(_mm_set1_epi32(last + 1), ix);

            //calculate horizontal interpolation factor for this pixel
            __m128 fx = _mm_add_ps(_mm_set_ps(3, 2, 1, 0), _mm_set1_ps(x));
            __m128 fa = _mm_mul_ps(_mm_sub_ps(_mm_set1_ps(maxx), fx), _mm_set1_ps(xfactor));
            __m128 fb = _mm_sub_ps(_mm_set1_ps(1), fa);
            __m128 z = _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(z1)),
                                  _mm_mul_ps(fb, _mm_set1_ps(z2)));

            //depth test
            __m128 zmask = _mm_and_ps(_mm_cmplt_ps(z, _mm_set1_ps(1)),
                           _mm_and_ps(_mm_cmpgt_ps(z, _mm_set1_ps(-1)),
                                      _mm_cmplt_ps(z, _mm_loadu_ps(zrow + x))));
            inclusionMask = _mm_and_si128(_mm_castps_si128(zmask), inclusionMask);

            //depth test early-out
            if(!_mm_movemask_epi8(inclusionMask)) {
                continue;
            }

            //interpolate vertex attributes
            // __m128 w = _mm_rcp_ps(_mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(w1)),
            //                                  _mm_mul_ps(fb, _mm_set1_ps(w2))));
            __m128 w = _mm_div_ps(_mm_set1_ps(1), _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(w1)),
                                                             _mm_mul_ps(fb, _mm_set1_ps(w2))));

            __m128 lx = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(l1.x)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(l2.x))));
            __m128 ly = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(l1.y)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(l2.y))));
            __m128 lz = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(l1.z)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(l2.z))));

            __m128 cx = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(c1.x)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(c2.x))));
            __m128 cy = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(c1.y)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(c2.y))));
            __m128 cz = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(c1.z)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(c2.z))));

            __m128 sx = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(s1.x)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(s2.x))));
            __m128 sy = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(s1.y)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(s2.y))));
            __m128 sz = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(s1.z)),
                                                 _mm_mul_ps(fb, _mm_set1_ps(s2.z))));

            __m128 u = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(u1)),
                                                _mm_mul_ps(fb, _mm_set1_ps(u2))));
            __m128 v = _mm_mul_ps(w, _mm_add_ps(_mm_mul_ps(fa, _mm_set1_ps(v1)),
                                                _mm_mul_ps(fb, _mm_set1_ps(v2))));



            //texture coord calculations
            __m128i wmask = _mm_set1_epi32(tex->wmask);
            __m128i hmask = _mm_set1_epi32(tex->hmask);
            __m128i iu1 = _mm_and_si128(_mm_cvttps_epi32(u), wmask);
            __m128i iv1 = _mm_and_si128(_mm_cvttps_epi32(v), hmask);
            __m128i iu2 =
                _mm_and_si128(_mm_add_epi32(_mm_cvttps_epi32(u), _mm_set1_epi32(1)), wmask);
            __m128i iv2 =
                _mm_and_si128(_mm_add_epi32(_mm_cvttps_epi32(v), _mm_set1_epi32(1)), hmask);

            __m128 uf1 = _mm_sub_ps(u, _mm_cvtepi32_ps(_mm_cvttps_epi32(u)));
            __m128 vf1 = _mm_sub_ps(v, _mm_cvtepi32_ps(_mm_cvttps_epi32(v)));
            __m128 uf2 = _mm_sub_ps(_mm_set1_ps(1), uf1);
            __m128 vf2 = _mm_sub_ps(_mm_set1_ps(1), vf1);

            //compute texture indices
            __m128i idx11 = _mm_add_epi32(_mm_madd_epi16(_mm_set1_epi32(tex->width), iv1), iu1);
            __m128i idx12 = _mm_add_epi32(_mm_madd_epi16(_mm_set1_epi32(tex->width), iv1), iu2);
            __m128i idx21 = _mm_add_epi32(_mm_madd_epi16(_mm_set1_epi32(tex->width), iv2), iu1);
            __m128i idx22 = _mm_add_epi32(_mm_madd_epi16(_mm_set1_epi32(tex->width), iv2), iu2);



            //fetch texture data
            __m128i p11 = gather_u32((u32 *) tex->pixels, idx11);
            __m128i p12 = gather_u32((u32 *) tex->pixels, idx12);
            __m128i p21 = gather_u32((u32 *) tex->pixels, idx21);
            __m128i p22 = gather_u32((u32 *) tex->pixels, idx22);

            //extract color channels as floats
            __m128i mask8 = _mm_set1_epi32(0x000000FF);
            __m128 inv255 = _mm_set1_ps(1.0f / 255);

            __m128 b11 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(p11,      mask8)), inv255);
            __m128 g11 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p11,  8), mask8)), inv255);
            __m128 r11 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p11, 16), mask8)), inv255);
            __m128 a11 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p11, 24), mask8)), inv255);

            __m128 b12 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(p12,      mask8)), inv255);
            __m128 g12 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p12,  8), mask8)), inv255);
            __m128 r12 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p12, 16), mask8)), inv255);
            __m128 a12 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p12, 24), mask8)), inv255);

            __m128 b21 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(p21,      mask8)), inv255);
            __m128 g21 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p21,  8), mask8)), inv255);
            __m128 r21 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p21, 16), mask8)), inv255);
            __m128 a21 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p21, 24), mask8)), inv255);

            __m128 b22 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(p22,      mask8)), inv255);
            __m128 g22 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p22,  8), mask8)), inv255);
            __m128 r22 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p22, 16), mask8)), inv255);
            __m128 a22 = _mm_mul_ps(_mm_cvtepi32_ps(_mm_and_si128(
                                                   _mm_srli_epi32(p22, 24), mask8)), inv255);

            //convert to linear space
            b11 = _mm_mul_ps(b11, b11);
            g11 = _mm_mul_ps(g11, g11);
            r11 = _mm_mul_ps(r11, r11);

            b12 = _mm_mul_ps(b12, b12);
            g12 = _mm_mul_ps(g12, g12);
            r12 = _mm_mul_ps(r12, r12);

            b21 = _mm_mul_ps(b21, b21);
            g21 = _mm_mul_ps(g21, g21);
            r21 = _mm_mul_ps(r21, r21);

            b22 = _mm_mul_ps(b22, b22);
            g22 = _mm_mul_ps(g22, g22);
            r22 = _mm_mul_ps(r22, r22);

            //blend texture samples
            __m128 b1 = _mm_add_ps(_mm_mul_ps(b11, uf2), _mm_mul_ps(b12, uf1));
            __m128 g1 = _mm_add_ps(_mm_mul_ps(g11, uf2), _mm_mul_ps(g12, uf1));
            __m128 r1 = _mm_add_ps(_mm_mul_ps(r11, uf2), _mm_mul_ps(r12, uf1));
            __m128 a1 = _mm_add_ps(_mm_mul_ps(a11, uf2), _mm_mul_ps(a12, uf1));

            __m128 b2 = _mm_add_ps(_mm_mul_ps(b21, uf2), _mm_mul_ps(b22, uf1));
            __m128 g2 = _mm_add_ps(_mm_mul_ps(g21, uf2), _mm_mul_ps(g22, uf1));
            __m128 r2 = _mm_add_ps(_mm_mul_ps(r21, uf2), _mm_mul_ps(r22, uf1));
            __m128 a2 = _mm_add_ps(_mm_mul_ps(a21, uf2), _mm_mul_ps(a22, uf1));

            __m128 cb = _mm_add_ps(_mm_mul_ps(b1, vf2), _mm_mul_ps(b2, vf1));
            __m128 cg = _mm_add_ps(_mm_mul_ps(g1, vf2), _mm_mul_ps(g2, vf1));
            __m128 cr = _mm_add_ps(_mm_mul_ps(r1, vf2), _mm_mul_ps(r2, vf1));
            __m128 roughness = _mm_add_ps(_mm_mul_ps(a1, vf2), _mm_mul_ps(a2, vf1));



            //fetch normal data
            __m128i n11 = gather_u32((u32 *) tex->normals, idx11);
            __m128i n12 = gather_u32((u32 *) tex->normals, idx12);
            __m128i n21 = gather_u32((u32 *) tex->normals, idx21);
            __m128i n22 = gather_u32((u32 *) tex->normals, idx22);

            //extract normal vectors as floats
            __m128 nz11 = _mm_cvtepi32_ps(_mm_and_si128(               n11,      mask8));
            __m128 ny11 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n11,  8), mask8));
            __m128 nx11 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n11, 16), mask8));

            __m128 nz12 = _mm_cvtepi32_ps(_mm_and_si128(               n12,      mask8));
            __m128 ny12 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n12,  8), mask8));
            __m128 nx12 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n12, 16), mask8));

            __m128 nz21 = _mm_cvtepi32_ps(_mm_and_si128(               n21,      mask8));
            __m128 ny21 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n21,  8), mask8));
            __m128 nx21 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n21, 16), mask8));

            __m128 nz22 = _mm_cvtepi32_ps(_mm_and_si128(               n22,      mask8));
            __m128 ny22 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n22,  8), mask8));
            __m128 nx22 = _mm_cvtepi32_ps(_mm_and_si128(_mm_srli_epi32(n22, 16), mask8));

            //blend normal samples
            __m128 nx1 = _mm_add_ps(_mm_mul_ps(nx11, uf2), _mm_mul_ps(nx12, uf1));
            __m128 ny1 = _mm_add_ps(_mm_mul_ps(ny11, uf2), _mm_mul_ps(ny12, uf1));
            __m128 nz1 = _mm_add_ps(_mm_mul_ps(nz11, uf2), _mm_mul_ps(nz12, uf1));
            __m128 nx2 = _mm_add_ps(_mm_mul_ps(nx21, uf2), _mm_mul_ps(nx22, uf1));
            __m128 ny2 = _mm_add_ps(_mm_mul_ps(ny21, uf2), _mm_mul_ps(ny22, uf1));
            __m128 nz2 = _mm_add_ps(_mm_mul_ps(nz21, uf2), _mm_mul_ps(nz22, uf1));
            __m128 nx = _mm_add_ps(_mm_mul_ps(nx1, vf2), _mm_mul_ps(nx2, vf1));
            __m128 ny = _mm_add_ps(_mm_mul_ps(ny1, vf2), _mm_mul_ps(ny2, vf1));
            __m128 nz = _mm_add_ps(_mm_mul_ps(nz1, vf2), _mm_mul_ps(nz2, vf1));



            //generate shadow mask
            __m128 shMask = _mm_castsi128_ps(inclusionMask);
            shMask = _mm_and_ps(shMask, _mm_cmpgt_ps(sx, _mm_set1_ps(0)));
            shMask = _mm_and_ps(shMask, _mm_cmpgt_ps(sy, _mm_set1_ps(0)));
            shMask = _mm_and_ps(shMask, _mm_cmplt_ps(sx, _mm_set1_ps(shadow->width - 1)));
            shMask = _mm_and_ps(shMask, _mm_cmplt_ps(sy, _mm_set1_ps(shadow->height - 1)));
            __m128 shad = _mm_set1_ps(1);
            if (_mm_movemask_epi8(_mm_castps_si128(shMask))) {
                //shadow sample coords
                __m128i ix = _mm_cvttps_epi32(sx);
                __m128i iy = _mm_cvttps_epi32(sy);
                __m128i wmask = _mm_set1_epi32(shadow->wmask);
                __m128i hmask = _mm_set1_epi32(shadow->hmask);
                __m128i is1 = _mm_and_si128(ix, wmask);
                __m128i it1 = _mm_and_si128(iy, hmask);
                __m128i is2 = _mm_and_si128(_mm_add_epi32(ix, _mm_set1_epi32(1)), wmask);
                __m128i it2 = _mm_and_si128(_mm_add_epi32(iy, _mm_set1_epi32(1)), hmask);

                //calculate sample indices
                __m128i width = _mm_set1_epi32(shadow->width);
                __m128i idx11 = _mm_add_epi32(_mm_madd_epi16(width, it1), is1);
                __m128i idx12 = _mm_add_epi32(_mm_madd_epi16(width, it1), is2);
                __m128i idx21 = _mm_add_epi32(_mm_madd_epi16(width, it2), is1);
                __m128i idx22 = _mm_add_epi32(_mm_madd_epi16(width, it2), is2);

                //sample shadow texture
                __m128 s11 = gather_f32(shadow->depth, idx11);
                __m128 s12 = gather_f32(shadow->depth, idx12);
                __m128 s21 = gather_f32(shadow->depth, idx21);
                __m128 s22 = gather_f32(shadow->depth, idx22);

                //adjust samples
                s11 = _mm_mul_ps(_mm_sub_ps(sz, s11), _mm_set1_ps(shadow->scale * 4));
                s12 = _mm_mul_ps(_mm_sub_ps(sz, s12), _mm_set1_ps(shadow->scale * 4));
                s21 = _mm_mul_ps(_mm_sub_ps(sz, s21), _mm_set1_ps(shadow->scale * 4));
                s22 = _mm_mul_ps(_mm_sub_ps(sz, s22), _mm_set1_ps(shadow->scale * 4));

                s11 = _mm_min_ps(_mm_set1_ps(1), _mm_max_ps(_mm_set1_ps(0), s11));
                s12 = _mm_min_ps(_mm_set1_ps(1), _mm_max_ps(_mm_set1_ps(0), s12));
                s21 = _mm_min_ps(_mm_set1_ps(1), _mm_max_ps(_mm_set1_ps(0), s21));
                s22 = _mm_min_ps(_mm_set1_ps(1), _mm_max_ps(_mm_set1_ps(0), s22));

                //blend samples
                __m128 sf1 = _mm_sub_ps(sx, _mm_cvtepi32_ps(ix));
                __m128 tf1 = _mm_sub_ps(sy, _mm_cvtepi32_ps(iy));
                __m128 sf2 = _mm_sub_ps(_mm_set1_ps(1), sf1);
                __m128 tf2 = _mm_sub_ps(_mm_set1_ps(1), tf1);

                __m128 s1 = _mm_add_ps(_mm_mul_ps(s11, sf2), _mm_mul_ps(s12, sf1));
                __m128 s2 = _mm_add_ps(_mm_mul_ps(s21, sf2), _mm_mul_ps(s22, sf1));
                __m128 s  = _mm_add_ps(_mm_mul_ps(s1, tf2), _mm_mul_ps(s2, tf1));

                //fold back into shadow mask
                __m128 invMask = _mm_xor_ps(shMask, _mm_castsi128_ps(_mm_set1_epi32(0xFFFFFFFF)));
                shad = _mm_or_ps(_mm_and_ps(shMask, s), _mm_and_ps(_mm_set1_ps(1), invMask));
            }



            //normalize normal vector
            nx = _mm_sub_ps(nx, _mm_set1_ps(127.5f));
            ny = _mm_sub_ps(ny, _mm_set1_ps(127.5f));
            nz = _mm_sub_ps(nz, _mm_set1_ps(127.5f));
            __m128 nmag = _mm_rsqrt_ps(_mm_add_ps(_mm_mul_ps(nx, nx),
                                       _mm_add_ps(_mm_mul_ps(ny, ny), _mm_mul_ps(nz, nz))));
            nx = _mm_mul_ps(nx, nmag);
            ny = _mm_mul_ps(ny, nmag);
            nz = _mm_mul_ps(nz, nmag);

            //normalize directional light vector
            __m128 lmag = _mm_rsqrt_ps(_mm_add_ps(_mm_mul_ps(lx, lx),
                                       _mm_add_ps(_mm_mul_ps(ly, ly), _mm_mul_ps(lz, lz))));
            lx = _mm_mul_ps(lx, lmag);
            ly = _mm_mul_ps(ly, lmag);
            lz = _mm_mul_ps(lz, lmag);

            //normalize camera direction vector
            __m128 cmag = _mm_rsqrt_ps(_mm_add_ps(_mm_mul_ps(cx, cx),
                                       _mm_add_ps(_mm_mul_ps(cy, cy), _mm_mul_ps(cz, cz))));
            cx = _mm_mul_ps(cx, cmag);
            cy = _mm_mul_ps(cy, cmag);
            cz = _mm_mul_ps(cz, cmag);

            //calculate direction diffuse light (n . l)
            __m128 lightDot = _mm_add_ps(_mm_mul_ps(nx, lx),
                              _mm_add_ps(_mm_mul_ps(ny, ly), _mm_mul_ps(nz, lz)));
            __m128 light = _mm_mul_ps(shad, _mm_max_ps(lightDot, _mm_set1_ps(0)));

            //apply light to color
            cb = _mm_mul_ps(cb, _mm_add_ps(light, _mm_set1_ps(0.04f)));
            cg = _mm_mul_ps(cg, _mm_add_ps(light, _mm_set1_ps(0.04f)));
            cr = _mm_mul_ps(cr, _mm_add_ps(light, _mm_set1_ps(0.04f)));

            //specular aspects
            __m128 rough = _mm_mul_ps(roughness, roughness);
            rough = _mm_mul_ps(roughness, _mm_mul_ps(rough, rough));
            __m128 e = _mm_rcp_ps(_mm_add_ps(rough, _mm_set1_ps(0.0001f)));
            __m128 m = _mm_rcp_ps(_mm_add_ps(rough, _mm_set1_ps(0.02f)));

            //reflection ray
            __m128 rx = _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(2), _mm_mul_ps(lightDot, nx)), lx);
            __m128 ry = _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(2), _mm_mul_ps(lightDot, ny)), ly);
            __m128 rz = _mm_sub_ps(_mm_mul_ps(_mm_set1_ps(2), _mm_mul_ps(lightDot, nz)), lz);

            //specular calculations
            __m128 refl = _mm_add_ps(_mm_mul_ps(rx, cx),
                          _mm_add_ps(_mm_mul_ps(ry, cy), _mm_mul_ps(rz, cz)));
            __m128 spec = _mm_max_ps(_mm_set1_ps(0), refl);
            __m128 spec_subexpr = _mm_rcp_ps(_mm_add_ps(_mm_sub_ps(e, _mm_mul_ps(e, spec)), spec));
            __m128 specular = _mm_add_ps(_mm_mul_ps(_mm_mul_ps(m, spec), spec_subexpr),
                                         _mm_set1_ps(0.02f));

            __m128 ior = _mm_set1_ps(1.333f); //ior of water, chosen for no particular reason
            __m128 f = _mm_mul_ps(_mm_sub_ps(_mm_set1_ps(1), ior),
                                  _mm_rcp_ps(_mm_add_ps(_mm_set1_ps(1), ior)));
            f = _mm_mul_ps(f, f);

            __m128 headon = _mm_add_ps(_mm_mul_ps(nx, cx),
                            _mm_add_ps(_mm_mul_ps(ny, cy), _mm_mul_ps(nz, cz)));
            headon = _mm_sub_ps(_mm_set1_ps(1), _mm_max_ps(_mm_set1_ps(0), headon));
            headon = _mm_mul_ps(_mm_mul_ps(_mm_mul_ps(headon, headon),
                                           _mm_mul_ps(headon, headon)), headon);

            __m128 invRough = _mm_sub_ps(_mm_set1_ps(1), roughness);
            __m128 fresnel = _mm_add_ps(f, _mm_mul_ps(_mm_sub_ps(_mm_set1_ps(1), f),
                                                      _mm_mul_ps(_mm_mul_ps(invRough, invRough),
                                                                 headon)));

            __m128 highlight = _mm_mul_ps(_mm_mul_ps(specular, shad), fresnel);
            cb = _mm_mul_ps(cb, _mm_sub_ps(_mm_set1_ps(1), fresnel));
            cg = _mm_mul_ps(cg, _mm_sub_ps(_mm_set1_ps(1), fresnel));
            cr = _mm_mul_ps(cr, _mm_sub_ps(_mm_set1_ps(1), fresnel));



            //convert to gamma space
            __m128 b = _mm_sqrt_ps(_mm_add_ps(cb, highlight));
            __m128 g = _mm_sqrt_ps(_mm_add_ps(cg, highlight));
            __m128 r = _mm_sqrt_ps(_mm_add_ps(cr, highlight));

            //write b, g, r to color buffer
            __m128 f255 = _mm_set1_ps(255);
            __m128i ib = _mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(b, _mm_set1_ps(1)), f255));
            __m128i ig = _mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(g, _mm_set1_ps(1)), f255));
            __m128i ir = _mm_cvtps_epi32(_mm_mul_ps(_mm_min_ps(r, _mm_set1_ps(1)), f255));
            // __m128i ib = _mm_and_si128(mask8, _mm_cvtps_epi32(_mm_mul_ps(b, f255)));
            // __m128i ig = _mm_and_si128(mask8, _mm_cvtps_epi32(_mm_mul_ps(g, f255)));
            // __m128i ir = _mm_and_si128(mask8, _mm_cvtps_epi32(_mm_mul_ps(r, f255)));
            __m128i out = _mm_or_si128(_mm_slli_epi32(ir, 16),
                          _mm_or_si128(_mm_slli_epi32(ig,  8), ib));
            _mm_maskmoveu_si128(out, inclusionMask, (char *)(row + x));

            //write to z buffer
            _mm_maskmoveu_si128(_mm_castps_si128(z), inclusionMask, (char *)(zrow + x));
        }
        perfInner += perf() - preInner;
    }
}

inline __m128i add(__m128i l, __m128i r) {
    return _mm_add_epi32(l, r);
}

inline __m128 add(__m128 l, __m128 r) {
    return _mm_add_ps(l, r);
}

inline __m128 sub(__m128 l, __m128 r) {
    return _mm_sub_ps(l, r);
}

inline __m128 mul(__m128 l, __m128 r) {
    return _mm_mul_ps(l, r);
}

inline __m128 _and(__m128 l, __m128 r) {
    return _mm_and_ps(l, r);
}

inline __m128i _and(__m128i l, __m128i r) {
    return _mm_and_si128(l, r);
}

inline __m128i epi32(int i) {
    return _mm_set1_epi32(i);
}

inline __m128i epi32(__m128 f) {
    return _mm_cvttps_epi32(f);
}

inline __m128 ps(float f) {
    return _mm_set1_ps(f);
}

inline __m128 ps(__m128i i) {
    return _mm_cvtepi32_ps(i);
}

inline __m128 lt(__m128 l, __m128 r) {
    return _mm_cmplt_ps(l, r);
}

inline __m128 gt(__m128 l, __m128 r) {
    return _mm_cmpgt_ps(l, r);
}

inline __m128i srli(__m128i l, int imm8) {
    return _mm_srli_epi32(l, imm8);
}

inline __m128 min(__m128 l, __m128 r) {
    return _mm_min_ps(l, r);
}

inline __m128 max(__m128 l, __m128 r) {
    return _mm_max_ps(l, r);
}

inline __m128 _xor(__m128 l, __m128 r) {
    return _mm_xor_ps(l, r);
}

inline __m128 rsqrt(__m128 f) {
    return _mm_rsqrt_ps(f);
}

inline __m128 rcp(__m128 f) {
    return _mm_rcp_ps(f);
}

inline __m128 sqrt(__m128 f) {
    return _mm_sqrt_ps(f);
}

void draw_triangle_sse2(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]) {
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
        for (int x = first; x <= last; x += 4) {
            //generate inclusion mask to avoid overwriting pixels
            __m128i ix = add(_mm_set_epi32(3, 2, 1, 0), epi32(x));
            __m128i inclusionMask = _mm_cmpgt_epi32(epi32(last + 1), ix);

            //calculate horizontal interpolation factor for this pixel
            __m128 fx = add(_mm_set_ps(3, 2, 1, 0), ps(x));
            __m128 fa = mul(sub(ps(maxx), fx), ps(xfactor));
            __m128 fb = sub(ps(1), fa);
            __m128 z = add(mul(fa, ps(z1)), mul(fb, ps(z2)));

            //depth test
            __m128 zmask = _and(lt(z, ps(1)), _and(gt(z, ps(-1)), lt(z, _mm_loadu_ps(zrow + x))));
            inclusionMask = _and(_mm_castps_si128(zmask), inclusionMask);

            //depth test early-out
            if(!_mm_movemask_epi8(inclusionMask)) {
                continue;
            }

            //interpolate vertex attributes
            // __m128 w = rcp(add(mul(fa, ps(w1)), mul(fb, ps(w2))));
            __m128 w = _mm_div_ps(ps(1), add(mul(fa, ps(w1)), mul(fb, ps(w2))));

            __m128 lx = mul(w, add(mul(fa, ps(l1.x)), mul(fb, ps(l2.x))));
            __m128 ly = mul(w, add(mul(fa, ps(l1.y)), mul(fb, ps(l2.y))));
            __m128 lz = mul(w, add(mul(fa, ps(l1.z)), mul(fb, ps(l2.z))));

            __m128 cx = mul(w, add(mul(fa, ps(c1.x)), mul(fb, ps(c2.x))));
            __m128 cy = mul(w, add(mul(fa, ps(c1.y)), mul(fb, ps(c2.y))));
            __m128 cz = mul(w, add(mul(fa, ps(c1.z)), mul(fb, ps(c2.z))));

            __m128 sx = mul(w, add(mul(fa, ps(s1.x)), mul(fb, ps(s2.x))));
            __m128 sy = mul(w, add(mul(fa, ps(s1.y)), mul(fb, ps(s2.y))));
            __m128 sz = mul(w, add(mul(fa, ps(s1.z)), mul(fb, ps(s2.z))));

            __m128 u = mul(w, add(mul(fa, ps(u1)), mul(fb, ps(u2))));
            __m128 v = mul(w, add(mul(fa, ps(v1)), mul(fb, ps(v2))));



            //texture coord calculations
            __m128i iu1 = _and(epi32(u), epi32(tex->wmask));
            __m128i iv1 = _and(epi32(v), epi32(tex->hmask));
            __m128i iu2 = _and(add(epi32(u), epi32(1)), epi32(tex->wmask));
            __m128i iv2 = _and(add(epi32(v), epi32(1)), epi32(tex->hmask));

            __m128 uf1 = sub(u, ps(epi32(u)));
            __m128 vf1 = sub(v, ps(epi32(v)));
            __m128 uf2 = sub(ps(1), uf1);
            __m128 vf2 = sub(ps(1), vf1);

            //compute texture indices
            __m128i idx11 = add(_mm_madd_epi16(epi32(tex->width), iv1), iu1);
            __m128i idx12 = add(_mm_madd_epi16(epi32(tex->width), iv1), iu2);
            __m128i idx21 = add(_mm_madd_epi16(epi32(tex->width), iv2), iu1);
            __m128i idx22 = add(_mm_madd_epi16(epi32(tex->width), iv2), iu2);



            //fetch texture data
            __m128i p11 = gather_u32((u32 *) tex->pixels, idx11);
            __m128i p12 = gather_u32((u32 *) tex->pixels, idx12);
            __m128i p21 = gather_u32((u32 *) tex->pixels, idx21);
            __m128i p22 = gather_u32((u32 *) tex->pixels, idx22);

            //extract color channels as floats
            __m128i mask8 = epi32(0x000000FF);
            __m128 inv255 = ps(1.0f / 255);

            __m128 b11 = mul(ps(_and(p11,      mask8)), inv255);
            __m128 g11 = mul(ps(_and(srli(p11,  8), mask8)), inv255);
            __m128 r11 = mul(ps(_and(srli(p11, 16), mask8)), inv255);
            __m128 a11 = mul(ps(_and(srli(p11, 24), mask8)), inv255);

            __m128 b12 = mul(ps(_and(p12,      mask8)), inv255);
            __m128 g12 = mul(ps(_and(srli(p12,  8), mask8)), inv255);
            __m128 r12 = mul(ps(_and(srli(p12, 16), mask8)), inv255);
            __m128 a12 = mul(ps(_and(srli(p12, 24), mask8)), inv255);

            __m128 b21 = mul(ps(_and(p21,      mask8)), inv255);
            __m128 g21 = mul(ps(_and(srli(p21,  8), mask8)), inv255);
            __m128 r21 = mul(ps(_and(srli(p21, 16), mask8)), inv255);
            __m128 a21 = mul(ps(_and(srli(p21, 24), mask8)), inv255);

            __m128 b22 = mul(ps(_and(p22,      mask8)), inv255);
            __m128 g22 = mul(ps(_and(srli(p22,  8), mask8)), inv255);
            __m128 r22 = mul(ps(_and(srli(p22, 16), mask8)), inv255);
            __m128 a22 = mul(ps(_and(srli(p22, 24), mask8)), inv255);

            //convert to linear space
            b11 = mul(b11, b11);
            g11 = mul(g11, g11);
            r11 = mul(r11, r11);

            b12 = mul(b12, b12);
            g12 = mul(g12, g12);
            r12 = mul(r12, r12);

            b21 = mul(b21, b21);
            g21 = mul(g21, g21);
            r21 = mul(r21, r21);

            b22 = mul(b22, b22);
            g22 = mul(g22, g22);
            r22 = mul(r22, r22);

            //blend texture samples
            __m128 b1 = add(mul(b11, uf2), mul(b12, uf1));
            __m128 g1 = add(mul(g11, uf2), mul(g12, uf1));
            __m128 r1 = add(mul(r11, uf2), mul(r12, uf1));
            __m128 a1 = add(mul(a11, uf2), mul(a12, uf1));

            __m128 b2 = add(mul(b21, uf2), mul(b22, uf1));
            __m128 g2 = add(mul(g21, uf2), mul(g22, uf1));
            __m128 r2 = add(mul(r21, uf2), mul(r22, uf1));
            __m128 a2 = add(mul(a21, uf2), mul(a22, uf1));

            __m128 cb = add(mul(b1, vf2), mul(b2, vf1));
            __m128 cg = add(mul(g1, vf2), mul(g2, vf1));
            __m128 cr = add(mul(r1, vf2), mul(r2, vf1));
            __m128 roughness = add(mul(a1, vf2), mul(a2, vf1));



            //fetch normal data
            __m128i n11 = gather_u32((u32 *) tex->normals, idx11);
            __m128i n12 = gather_u32((u32 *) tex->normals, idx12);
            __m128i n21 = gather_u32((u32 *) tex->normals, idx21);
            __m128i n22 = gather_u32((u32 *) tex->normals, idx22);

            //extract normal vectors as floats
            __m128 nz11 = ps(_and(     n11,      mask8));
            __m128 ny11 = ps(_and(srli(n11,  8), mask8));
            __m128 nx11 = ps(_and(srli(n11, 16), mask8));

            __m128 nz12 = ps(_and(     n12,      mask8));
            __m128 ny12 = ps(_and(srli(n12,  8), mask8));
            __m128 nx12 = ps(_and(srli(n12, 16), mask8));

            __m128 nz21 = ps(_and(     n21,      mask8));
            __m128 ny21 = ps(_and(srli(n21,  8), mask8));
            __m128 nx21 = ps(_and(srli(n21, 16), mask8));

            __m128 nz22 = ps(_and(     n22,      mask8));
            __m128 ny22 = ps(_and(srli(n22,  8), mask8));
            __m128 nx22 = ps(_and(srli(n22, 16), mask8));

            //blend normal samples
            __m128 nx1 = add(mul(nx11, uf2), mul(nx12, uf1));
            __m128 ny1 = add(mul(ny11, uf2), mul(ny12, uf1));
            __m128 nz1 = add(mul(nz11, uf2), mul(nz12, uf1));
            __m128 nx2 = add(mul(nx21, uf2), mul(nx22, uf1));
            __m128 ny2 = add(mul(ny21, uf2), mul(ny22, uf1));
            __m128 nz2 = add(mul(nz21, uf2), mul(nz22, uf1));
            __m128 nx = add(mul(nx1, vf2), mul(nx2, vf1));
            __m128 ny = add(mul(ny1, vf2), mul(ny2, vf1));
            __m128 nz = add(mul(nz1, vf2), mul(nz2, vf1));



            //generate shadow mask
            __m128 shMask = _mm_castsi128_ps(inclusionMask);
            shMask = _and(shMask, gt(sx, ps(0)));
            shMask = _and(shMask, gt(sy, ps(0)));
            shMask = _and(shMask, lt(sx, ps(shadow->width - 1)));
            shMask = _and(shMask, lt(sy, ps(shadow->height - 1)));
            __m128 shad = ps(1);
            if (_mm_movemask_epi8(_mm_castps_si128(shMask))) {
                //shadow sample coords
                __m128i ix = epi32(sx);
                __m128i iy = epi32(sy);
                __m128i wmask = epi32(shadow->wmask);
                __m128i hmask = epi32(shadow->hmask);
                __m128i is1 = _and(ix, wmask);
                __m128i it1 = _and(iy, hmask);
                __m128i is2 = _and(add(ix, epi32(1)), wmask);
                __m128i it2 = _and(add(iy, epi32(1)), hmask);

                //calculate sample indices
                __m128i width = epi32(shadow->width);
                __m128i idx11 = add(_mm_madd_epi16(width, it1), is1);
                __m128i idx12 = add(_mm_madd_epi16(width, it1), is2);
                __m128i idx21 = add(_mm_madd_epi16(width, it2), is1);
                __m128i idx22 = add(_mm_madd_epi16(width, it2), is2);

                //sample shadow texture
                __m128 s11 = gather_f32(shadow->depth, idx11);
                __m128 s12 = gather_f32(shadow->depth, idx12);
                __m128 s21 = gather_f32(shadow->depth, idx21);
                __m128 s22 = gather_f32(shadow->depth, idx22);

                //adjust samples
                s11 = min(ps(1), max(ps(0), mul(sub(sz, s11), ps(shadow->scale * 4))));
                s12 = min(ps(1), max(ps(0), mul(sub(sz, s12), ps(shadow->scale * 4))));
                s21 = min(ps(1), max(ps(0), mul(sub(sz, s21), ps(shadow->scale * 4))));
                s22 = min(ps(1), max(ps(0), mul(sub(sz, s22), ps(shadow->scale * 4))));

                //blend samples
                __m128 sf1 = sub(sx, ps(ix));
                __m128 tf1 = sub(sy, ps(iy));
                __m128 sf2 = sub(ps(1), sf1);
                __m128 tf2 = sub(ps(1), tf1);

                __m128 s1 = add(mul(s11, sf2), mul(s12, sf1));
                __m128 s2 = add(mul(s21, sf2), mul(s22, sf1));
                __m128 s  = add(mul(s1, tf2), mul(s2, tf1));

                //fold back into shadow mask
                __m128 invMask = _xor(shMask, _mm_castsi128_ps(epi32(0xFFFFFFFF)));
                shad = _xor(_and(shMask, s), _and(ps(1), invMask));
            }



            //normalize normal vector
            nx = sub(nx, ps(127.5f));
            ny = sub(ny, ps(127.5f));
            nz = sub(nz, ps(127.5f));
            __m128 nmag = rsqrt(add(mul(nx, nx), add(mul(ny, ny), mul(nz, nz))));
            nx = mul(nx, nmag);
            ny = mul(ny, nmag);
            nz = mul(nz, nmag);

            //normalize directional light vector
            __m128 lmag = rsqrt(add(mul(lx, lx), add(mul(ly, ly), mul(lz, lz))));
            lx = mul(lx, lmag);
            ly = mul(ly, lmag);
            lz = mul(lz, lmag);

            //normalize camera direction vector
            __m128 cmag = rsqrt(add(mul(cx, cx), add(mul(cy, cy), mul(cz, cz))));
            cx = mul(cx, cmag);
            cy = mul(cy, cmag);
            cz = mul(cz, cmag);

            //calculate direction diffuse light (n . l)
            __m128 lightDot = add(mul(nx, lx), add(mul(ny, ly), mul(nz, lz)));
            __m128 light = mul(shad, max(lightDot, ps(0)));

            //apply light to color
            cb = mul(cb, add(light, ps(0.04f)));
            cg = mul(cg, add(light, ps(0.04f)));
            cr = mul(cr, add(light, ps(0.04f)));

            //specular aspects
            __m128 rough = mul(roughness, roughness);
            rough = mul(roughness, mul(rough, rough));
            __m128 e = rcp(add(rough, ps(0.0001f)));
            __m128 m = rcp(add(rough, ps(0.02f)));

            //reflection ray
            __m128 rx = sub(mul(ps(2), mul(lightDot, nx)), lx);
            __m128 ry = sub(mul(ps(2), mul(lightDot, ny)), ly);
            __m128 rz = sub(mul(ps(2), mul(lightDot, nz)), lz);

            //specular calculations
            __m128 refl = add(mul(rx, cx), add(mul(ry, cy), mul(rz, cz)));
            __m128 spec = max(ps(0), refl);
            __m128 spec_subexpr = rcp(add(sub(e, mul(e, spec)), spec));
            __m128 specular = add(mul(mul(m, spec), spec_subexpr), ps(0.02f));

            __m128 ior = ps(1.333f); //ior of water, chosen for no particular reason
            __m128 f = mul(sub(ps(1), ior), rcp(add(ps(1), ior)));
            f = mul(f, f);

            __m128 headon = add(mul(nx, cx), add(mul(ny, cy), mul(nz, cz)));
            headon = sub(ps(1), max(ps(0), headon));
            headon = mul(mul(mul(headon, headon), mul(headon, headon)), headon);

            __m128 invRough = sub(ps(1), roughness);
            __m128 fresnel = add(f, mul(sub(ps(1), f), mul(mul(invRough, invRough), headon)));

            __m128 highlight = mul(mul(specular, shad), fresnel);
            cb = mul(cb, sub(ps(1), fresnel));
            cg = mul(cg, sub(ps(1), fresnel));
            cr = mul(cr, sub(ps(1), fresnel));



            //convert to gamma space
            __m128 b = sqrt(add(cb, highlight));
            __m128 g = sqrt(add(cg, highlight));
            __m128 r = sqrt(add(cr, highlight));

            //write b, g, r to color buffer
            __m128 f255 = ps(255);
            __m128i ib = epi32(mul(min(b, ps(1)), f255));
            __m128i ig = epi32(mul(min(g, ps(1)), f255));
            __m128i ir = epi32(mul(min(r, ps(1)), f255));
            // __m128i ib = _and(mask8, epi32(mul(b, f255)));
            // __m128i ig = _and(mask8, epi32(mul(g, f255)));
            // __m128i ir = _and(mask8, epi32(mul(r, f255)));
            __m128i out = _mm_or_si128(_mm_slli_epi32(ir, 16),
                          _mm_or_si128(_mm_slli_epi32(ig,  8), ib));
            _mm_maskmoveu_si128(out, inclusionMask, (char *)(row + x));

            //write to z buffer
            _mm_maskmoveu_si128(_mm_castps_si128(z), inclusionMask, (char *)(zrow + x));
        }
        perfInner += perf() - preInner;
    }
}

struct f128;

struct i128 {
    __m128i m;
    operator f128();
};

struct f128 {
    __m128 m;
    operator i128() { return { _mm_cvtps_epi32(m) }; }
};

i128::operator f128() { return { _mm_cvtepi32_ps(m) }; }

i128 cast_i128(f128 f) { return { _mm_castps_si128(f.m) }; }
f128 cast_f128(i128 i) { return { _mm_castsi128_ps(i.m) }; }

i128 ii128(int i) { return { _mm_set1_epi32(i) }; }
i128 ii128(int a, int b, int c, int d) { return { _mm_set_epi32(a, b, c, d) }; }
f128 ff128(float f) { return { _mm_set1_ps(f) }; }
f128 ff128(float a, float b, float c, float d) { return { _mm_set_ps(a, b, c, d) }; }

i128 operator<<(i128 l, int i8) { return { _mm_slli_epi32 (l.m, i8 ) }; }
i128 operator>>(i128 l, int i8) { return { _mm_srli_epi32 (l.m, i8 ) }; }
i128 operator+ (i128 l, i128 r) { return { _mm_add_epi32  (l.m, r.m) }; }
i128 operator- (i128 l, i128 r) { return { _mm_sub_epi32  (l.m, r.m) }; }
i128 operator& (i128 l, i128 r) { return { _mm_and_si128  (l.m, r.m) }; }
i128 operator| (i128 l, i128 r) { return { _mm_or_si128   (l.m, r.m) }; }
i128 operator^ (i128 l, i128 r) { return { _mm_xor_si128  (l.m, r.m) }; }
i128 operator> (i128 l, i128 r) { return { _mm_cmpgt_epi32(l.m, r.m) }; }
i128 operator< (i128 l, i128 r) { return { _mm_cmplt_epi32(l.m, r.m) }; }

i128 operator+ (i128 l, int r) { return { _mm_add_epi32  (l.m, _mm_set1_epi32(r)) }; }
i128 operator- (i128 l, int r) { return { _mm_sub_epi32  (l.m, _mm_set1_epi32(r)) }; }
i128 operator& (i128 l, int r) { return { _mm_and_si128  (l.m, _mm_set1_epi32(r)) }; }
i128 operator| (i128 l, int r) { return { _mm_or_si128   (l.m, _mm_set1_epi32(r)) }; }
i128 operator^ (i128 l, int r) { return { _mm_xor_si128  (l.m, _mm_set1_epi32(r)) }; }
i128 operator> (i128 l, int r) { return { _mm_cmpgt_epi32(l.m, _mm_set1_epi32(r)) }; }
i128 operator< (i128 l, int r) { return { _mm_cmplt_epi32(l.m, _mm_set1_epi32(r)) }; }

i128 operator+ (int l, i128 r) { return { _mm_add_epi32  (_mm_set1_epi32(l), r.m) }; }
i128 operator- (int l, i128 r) { return { _mm_sub_epi32  (_mm_set1_epi32(l), r.m) }; }
i128 operator& (int l, i128 r) { return { _mm_and_si128  (_mm_set1_epi32(l), r.m) }; }
i128 operator| (int l, i128 r) { return { _mm_or_si128   (_mm_set1_epi32(l), r.m) }; }
i128 operator^ (int l, i128 r) { return { _mm_xor_si128  (_mm_set1_epi32(l), r.m) }; }
i128 operator> (int l, i128 r) { return { _mm_cmpgt_epi32(_mm_set1_epi32(l), r.m) }; }
i128 operator< (int l, i128 r) { return { _mm_cmplt_epi32(_mm_set1_epi32(l), r.m) }; }

f128 operator* (f128 l, f128 r) { return { _mm_mul_ps    (l.m, r.m) }; }
f128 operator/ (f128 l, f128 r) { return { _mm_div_ps    (l.m, r.m) }; }
f128 operator+ (f128 l, f128 r) { return { _mm_add_ps    (l.m, r.m) }; }
f128 operator- (f128 l, f128 r) { return { _mm_sub_ps    (l.m, r.m) }; }
f128 operator& (f128 l, f128 r) { return { _mm_and_ps    (l.m, r.m) }; }
f128 operator| (f128 l, f128 r) { return { _mm_or_ps     (l.m, r.m) }; }
f128 operator^ (f128 l, f128 r) { return { _mm_xor_ps    (l.m, r.m) }; }
f128 operator> (f128 l, f128 r) { return { _mm_cmpgt_ps  (l.m, r.m) }; }
f128 operator< (f128 l, f128 r) { return { _mm_cmplt_ps  (l.m, r.m) }; }

f128 operator* (f128 l, float r) { return { _mm_mul_ps  (l.m, _mm_set1_ps(r)) }; }
f128 operator/ (f128 l, float r) { return { _mm_div_ps  (l.m, _mm_set1_ps(r)) }; }
f128 operator+ (f128 l, float r) { return { _mm_add_ps  (l.m, _mm_set1_ps(r)) }; }
f128 operator- (f128 l, float r) { return { _mm_sub_ps  (l.m, _mm_set1_ps(r)) }; }
f128 operator& (f128 l, float r) { return { _mm_and_ps  (l.m, _mm_set1_ps(r)) }; }
f128 operator| (f128 l, float r) { return { _mm_or_ps   (l.m, _mm_set1_ps(r)) }; }
f128 operator^ (f128 l, float r) { return { _mm_xor_ps  (l.m, _mm_set1_ps(r)) }; }
f128 operator> (f128 l, float r) { return { _mm_cmpgt_ps(l.m, _mm_set1_ps(r)) }; }
f128 operator< (f128 l, float r) { return { _mm_cmplt_ps(l.m, _mm_set1_ps(r)) }; }

f128 operator* (float l, f128 r) { return { _mm_mul_ps  (_mm_set1_ps(l), r.m) }; }
f128 operator/ (float l, f128 r) { return { _mm_div_ps  (_mm_set1_ps(l), r.m) }; }
f128 operator+ (float l, f128 r) { return { _mm_add_ps  (_mm_set1_ps(l), r.m) }; }
f128 operator- (float l, f128 r) { return { _mm_sub_ps  (_mm_set1_ps(l), r.m) }; }
f128 operator& (float l, f128 r) { return { _mm_and_ps  (_mm_set1_ps(l), r.m) }; }
f128 operator| (float l, f128 r) { return { _mm_or_ps   (_mm_set1_ps(l), r.m) }; }
f128 operator^ (float l, f128 r) { return { _mm_xor_ps  (_mm_set1_ps(l), r.m) }; }
f128 operator> (float l, f128 r) { return { _mm_cmpgt_ps(_mm_set1_ps(l), r.m) }; }
f128 operator< (float l, f128 r) { return { _mm_cmplt_ps(_mm_set1_ps(l), r.m) }; }

f128 min(f128 l, f128 r) { return { _mm_min_ps(l.m, r.m) }; }
f128 max(f128 l, f128 r) { return { _mm_max_ps(l.m, r.m) }; }

f128 min(f128 l, float r) { return { _mm_min_ps(l.m, _mm_set1_ps(r)) }; }
f128 max(f128 l, float r) { return { _mm_max_ps(l.m, _mm_set1_ps(r)) }; }

f128 min(float l, f128 r) { return { _mm_min_ps(_mm_set1_ps(l), r.m) }; }
f128 max(float l, f128 r) { return { _mm_max_ps(_mm_set1_ps(l), r.m) }; }

i128 floori(f128 f) { return { _mm_cvttps_epi32(f.m) }; }
f128 floor(f128 f) { return { _mm_cvtepi32_ps(_mm_cvttps_epi32(f.m)) }; }
f128 rsqrt(f128 f) { return { _mm_rsqrt_ps(f.m) }; }
f128 sqrt(f128 f) { return { _mm_sqrt_ps(f.m) }; }
f128 rcp(f128 f) { return { _mm_rcp_ps(f.m) }; }

void draw_triangle_sse2_struct(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]) {
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
        for (int x = first; x <= last; x += 4) {
            //generate inclusion mask to avoid overwriting pixels
            i128 ix = ii128(3, 2, 1, 0) + ii128(x);
            i128 inclusionMask = ii128(last + 1) > ix;

            //calculate horizontal interpolation factor for this pixel
            f128 fx = ix;
            f128 fa = (maxx - fx) * xfactor;
            f128 fb = 1.0f - fa;
            f128 z = fa * z1 + fb * z2;

            //depth test
            f128 zmask = (z < 1.0f) & (z > -1.0f) & (z < f128{_mm_loadu_ps(zrow + x)});
            inclusionMask = cast_i128(zmask) & inclusionMask;

            //depth test early-out
            if(!_mm_movemask_epi8(inclusionMask.m)) {
                continue;
            }

            //interpolate vertex attributes
            // f128 w = rcp(fa * w1 + fb * w2);
            f128 w = 1.0f / (fa * w1 + fb * w2);

            f128 lx = w * (fa * l1.x + fb * l2.x);
            f128 ly = w * (fa * l1.y + fb * l2.y);
            f128 lz = w * (fa * l1.z + fb * l2.z);
            f128 cx = w * (fa * c1.x + fb * c2.x);
            f128 cy = w * (fa * c1.y + fb * c2.y);
            f128 cz = w * (fa * c1.z + fb * c2.z);
            f128 sx = w * (fa * s1.x + fb * s2.x);
            f128 sy = w * (fa * s1.y + fb * s2.y);
            f128 sz = w * (fa * s1.z + fb * s2.z);
            f128 u  = w * (fa * u1   + fb * u2  );
            f128 v  = w * (fa * v1   + fb * v2  );



            //texture coord calculations
            i128 iu1 =  floori(u)      & tex->wmask;
            i128 iv1 =  floori(v)      & tex->hmask;
            i128 iu2 = (floori(u) + 1) & tex->wmask;
            i128 iv2 = (floori(v) + 1) & tex->hmask;

            f128 uf1 = u - floor(u);
            f128 vf1 = v - floor(v);
            f128 uf2 = 1.0f - uf1;
            f128 vf2 = 1.0f - vf1;

            //compute texture indices
            i128 idx11 = i128{_mm_madd_epi16(_mm_set1_epi32(tex->width), iv1.m)} + iu1;
            i128 idx12 = i128{_mm_madd_epi16(_mm_set1_epi32(tex->width), iv1.m)} + iu2;
            i128 idx21 = i128{_mm_madd_epi16(_mm_set1_epi32(tex->width), iv2.m)} + iu1;
            i128 idx22 = i128{_mm_madd_epi16(_mm_set1_epi32(tex->width), iv2.m)} + iu2;



            //fetch texture data
            i128 p11 = { gather_u32((u32 *) tex->pixels, idx11.m) };
            i128 p12 = { gather_u32((u32 *) tex->pixels, idx12.m) };
            i128 p21 = { gather_u32((u32 *) tex->pixels, idx21.m) };
            i128 p22 = { gather_u32((u32 *) tex->pixels, idx22.m) };

            //extract color channels as floats in linear space
            f128 b11 = sq(( p11        & 0xFF) * (1.0f / 255));
            f128 g11 = sq(((p11 >>  8) & 0xFF) * (1.0f / 255));
            f128 r11 = sq(((p11 >> 16) & 0xFF) * (1.0f / 255));
            f128 a11 =    ((p11 >> 24) & 0xFF) * (1.0f / 255) ;
            f128 b12 = sq(( p12        & 0xFF) * (1.0f / 255));
            f128 g12 = sq(((p12 >>  8) & 0xFF) * (1.0f / 255));
            f128 r12 = sq(((p12 >> 16) & 0xFF) * (1.0f / 255));
            f128 a12 =    ((p12 >> 24) & 0xFF) * (1.0f / 255) ;
            f128 b21 = sq(( p21        & 0xFF) * (1.0f / 255));
            f128 g21 = sq(((p21 >>  8) & 0xFF) * (1.0f / 255));
            f128 r21 = sq(((p21 >> 16) & 0xFF) * (1.0f / 255));
            f128 a21 =    ((p21 >> 24) & 0xFF) * (1.0f / 255) ;
            f128 b22 = sq(( p22        & 0xFF) * (1.0f / 255));
            f128 g22 = sq(((p22 >>  8) & 0xFF) * (1.0f / 255));
            f128 r22 = sq(((p22 >> 16) & 0xFF) * (1.0f / 255));
            f128 a22 =    ((p22 >> 24) & 0xFF) * (1.0f / 255) ;

            //blend texture samples
            f128 b1 = b11 * uf2 + b12 * uf1;
            f128 g1 = g11 * uf2 + g12 * uf1;
            f128 r1 = r11 * uf2 + r12 * uf1;
            f128 a1 = a11 * uf2 + a12 * uf1;
            f128 b2 = b21 * uf2 + b22 * uf1;
            f128 g2 = g21 * uf2 + g22 * uf1;
            f128 r2 = r21 * uf2 + r22 * uf1;
            f128 a2 = a21 * uf2 + a22 * uf1;
            f128 cb = b1  * vf2 + b2  * vf1;
            f128 cg = g1  * vf2 + g2  * vf1;
            f128 cr = r1  * vf2 + r2  * vf1;
            f128 roughness = a1 * vf2 + a2 * vf1;



            //fetch normal data
            i128 n11 = { gather_u32((u32 *) tex->normals, idx11.m) };
            i128 n12 = { gather_u32((u32 *) tex->normals, idx12.m) };
            i128 n21 = { gather_u32((u32 *) tex->normals, idx21.m) };
            i128 n22 = { gather_u32((u32 *) tex->normals, idx22.m) };

            //extract normal vectors as floats
            f128 nz11 =  n11        & 0xFF;
            f128 ny11 = (n11 >>  8) & 0xFF;
            f128 nx11 = (n11 >> 16) & 0xFF;
            f128 nz12 =  n12        & 0xFF;
            f128 ny12 = (n12 >>  8) & 0xFF;
            f128 nx12 = (n12 >> 16) & 0xFF;
            f128 nz21 =  n21        & 0xFF;
            f128 ny21 = (n21 >>  8) & 0xFF;
            f128 nx21 = (n21 >> 16) & 0xFF;
            f128 nz22 =  n22        & 0xFF;
            f128 ny22 = (n22 >>  8) & 0xFF;
            f128 nx22 = (n22 >> 16) & 0xFF;

            //blend normal samples
            f128 nx1 = nx11 * uf2 + nx12 * uf1;
            f128 ny1 = ny11 * uf2 + ny12 * uf1;
            f128 nz1 = nz11 * uf2 + nz12 * uf1;
            f128 nx2 = nx21 * uf2 + nx22 * uf1;
            f128 ny2 = ny21 * uf2 + ny22 * uf1;
            f128 nz2 = nz21 * uf2 + nz22 * uf1;
            f128 nx  = nx1  * vf2 + nx2  * vf1;
            f128 ny  = ny1  * vf2 + ny2  * vf1;
            f128 nz  = nz1  * vf2 + nz2  * vf1;



            //generate shadow mask
            f128 shMask = cast_f128(inclusionMask) &
                    (sx > 0.0f) & (sy > 0.0f) &
                    (sx < shadow->width - 1.0f) & (sy < shadow->height - 1.0f);
            f128 shad = ff128(1);
            if (_mm_movemask_epi8(_mm_castps_si128(shMask.m))) {
                //shadow sample coords
                i128 is1 =  floori(sx)      & shadow->wmask;
                i128 it1 =  floori(sy)      & shadow->hmask;
                i128 is2 = (floori(sx) + 1) & shadow->wmask;
                i128 it2 = (floori(sy) + 1) & shadow->hmask;

                //calculate sample indices
                i128 idx11 = i128{_mm_madd_epi16(ii128(shadow->width).m, it1.m)} + is1;
                i128 idx12 = i128{_mm_madd_epi16(ii128(shadow->width).m, it1.m)} + is2;
                i128 idx21 = i128{_mm_madd_epi16(ii128(shadow->width).m, it2.m)} + is1;
                i128 idx22 = i128{_mm_madd_epi16(ii128(shadow->width).m, it2.m)} + is2;

                //sample shadow texture
                f128 s11 = { gather_f32(shadow->depth, idx11.m) };
                f128 s12 = { gather_f32(shadow->depth, idx12.m) };
                f128 s21 = { gather_f32(shadow->depth, idx21.m) };
                f128 s22 = { gather_f32(shadow->depth, idx22.m) };

                //adjust samples
                s11 = min(1, max(0, (sz - s11) * (shadow->scale * 4)));
                s12 = min(1, max(0, (sz - s12) * (shadow->scale * 4)));
                s21 = min(1, max(0, (sz - s21) * (shadow->scale * 4)));
                s22 = min(1, max(0, (sz - s22) * (shadow->scale * 4)));

                //blend samples
                f128 sf = sx - floor(sx);
                f128 tf = sy - floor(sy);
                f128 s1 = s11 * (1.0f - sf) + s12 * sf;
                f128 s2 = s21 * (1.0f - sf) + s22 * sf;
                f128 s  = s1  * (1.0f - tf) + s2  * tf;

                //fold back into shadow mask
                f128 invMask = shMask ^ cast_f128(ii128(0xFFFFFFFF));
                shad = (shMask & s) | (1.0f & invMask);
            }



            //normalize normal vector
            nx = nx - 127.5f;
            ny = ny - 127.5f;
            nz = nz - 127.5f;
            f128 nmag = rsqrt(nx * nx + ny * ny + nz * nz);
            nx = nx * nmag;
            ny = ny * nmag;
            nz = nz * nmag;

            //normalize directional light vector
            // f128 lmag = 1 / sqrt(lx * lx + ly * ly + lz * lz);
            f128 lmag = rsqrt(lx * lx + ly * ly + lz * lz);
            lx = lx * lmag;
            ly = ly * lmag;
            lz = lz * lmag;

            //normalize camera direction vector
            // f128 cmag = 1 / sqrt(cx * cx + cy * cy + cz * cz);
            f128 cmag = rsqrt(cx * cx + cy * cy + cz * cz);
            cx = cx * cmag;
            cy = cy * cmag;
            cz = cz * cmag;

            //calculate direction diffuse light (n . l)
            f128 lightDot = nx * lx + ny * ly + nz * lz;
            f128 light = shad * max(lightDot, 0);

            //apply light to color
            cb = cb * (light + 0.04f);
            cg = cg * (light + 0.04f);
            cr = cr * (light + 0.04f);

            //specular aspects
            f128 rough = sq(sq(roughness)) * roughness;
            f128 e = rcp(rough + 0.0001f);
            f128 m = rcp(rough + 0.03f);

            //reflection ray
            f128 rx = 2 * lightDot * nx - lx;
            f128 ry = 2 * lightDot * ny - ly;
            f128 rz = 2 * lightDot * nz - lz;

            //specular calculations//
            f128 spec = max(0, rx * cx + ry * cy + rz * cz);
            f128 specular = m * spec * rcp(e - e * spec + spec) + 0.02f;

            f128 ior = ff128(1.333f); //ior of water, chosen for no particular reason
            f128 f = sq((1.0f - ior) / (1.0f + ior));

            f128 headon = nx * cx + ny * cy + nz * cz;
            headon = 1.0f - max(0, headon);
            headon = headon * headon * headon * headon * headon;

            f128 fresnel = f + (1.0f - f) * sq(1.0f - roughness) * headon;

            f128 highlight = specular * shad * fresnel;
            cb = cb * (1.0f - fresnel);
            cg = cg * (1.0f - fresnel);
            cr = cr * (1.0f - fresnel);



            //convert to gamma space
            f128 b = rcp(rsqrt(cb + highlight));
            f128 g = rcp(rsqrt(cg + highlight));
            f128 r = rcp(rsqrt(cr + highlight));

            //write b, g, r to color buffer
            i128 ib = min(b, 1) * 255;
            i128 ig = min(g, 1) * 255;
            i128 ir = min(r, 1) * 255;
            i128 out = (ir << 16) | (ig <<  8) | ib;
            _mm_maskmoveu_si128(out.m, inclusionMask.m, (char *)(row + x));

            //write to z buffer
            _mm_maskmoveu_si128(cast_i128(z).m, inclusionMask.m, (char *)(zrow + x));
        }
        perfInner += perf() - preInner;
    }
}

