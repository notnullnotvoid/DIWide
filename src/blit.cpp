#include "blit.hpp"

#include <stdlib.h>

void * aligned_alloc(size_t align, size_t size) {
    size_t mask = align - 1;
    assert(!(align & mask)); //ensure alignment is power-of-two
    void * ret = malloc(size);
    size_t addr = (size_t) ret;
    assert(!(addr & mask)); //ensure allocation is aligned
    return ret;
}

Canvas create_canvas(int width, int height, int margin) {
    //NOTE: we allocate extra buffer margin around the explicitly used canvas area
    //      so that operations can safely read from slightly outside the used area
    //      of the canvas. this simplifies and speeds up some operations because it
    //      eliminates the need for bounds checks and explicit handling of edge cases
    int canvasBytes = (width + 2 * margin) * (height + 2 * margin) * sizeof(Pixel);
    int bufferBytes = width * height * sizeof(float);
    assert(canvasBytes % 32 == 0); //aligned alloc will fail if this is not true!
    assert(bufferBytes % 32 == 0); //aligned alloc will fail if this is not true!
    //NOTE: we allocate on 32-byte boundaries for AVX speed,
    //      and pad extra bytes on the end so that unmasked SIMD reads won't go out of bounds
    Pixel * canvasData = (Pixel *) aligned_alloc(32, canvasBytes + 32);
    float * depthData = (float *) aligned_alloc(32, bufferBytes + 32);
    Canvas canv = {
        canvasData + margin * (width + 2 * margin) + margin, depthData,
        width, height, width + 2 * margin, width
    };
    //fill whole canvas (including margins) with solid black
    for (int y = 0; y < canv.height; ++y) {
        Pixel * row = canv.pixels + y * canv.pitch;
        for (int x = 0; x < canv.width; ++x) {
            row[x] = { 255, 0, 255, 255 };
        }
    }
    return canv;
}

//same as above, but doesn't allocate a color buffer for the canvas
ZBuffer create_depth_buffer(int width, int height) {
    int bufferBytes = width * height * sizeof(Pixel);
    assert(bufferBytes % 32 == 0); //aligned alloc will fail if this is not true!
    ZBuffer ret = { (float *) aligned_alloc(32, bufferBytes + 32), 0, 0,
                    width, height, (u32)(width - 1), (u32)(height - 1) };
    assert(!(ret.width & ret.wmask) && !(ret.height & ret.hmask));
    return ret;
}

void debug_depth_blit(SDL_Surface * surface, ZBuffer * depth) {
    int maxy = depth->height < surface->h? depth->height : surface->h;
    int maxx = depth->width < surface->w? depth->width : surface->w;
    for (int y = 0; y < maxy; ++y) {
        Pixel * row = (Pixel *) ((u8 *) surface->pixels + y * surface->pitch);
        float * zrow = depth->depth + y * depth->width;
        for (int x = 0; x < maxx; ++x) {
            u8 z = zrow[x] * 255;
            row[x] = { z, z, z, 255 };
        }
    }
}

void fast_scaled_blit(SDL_Surface * surface, Canvas * canvas, int scale) {
    if (scale == 1) {
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * src = canvas->pixels + y * canvas->pitch;
            Pixel * dest = (Pixel *) ((u8 *) surface->pixels + y * surface->pitch);
            for (int x = 0; x < canvas->width; ++x) {
                dest[x] = src[x];
            }
        }
    } else if (scale == 2) {
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * src = canvas->pixels + y * canvas->pitch;
            Pixel * dest1 = (Pixel *) ((u8 *) surface->pixels + (y * 2 + 0) * surface->pitch);
            Pixel * dest2 = (Pixel *) ((u8 *) surface->pixels + (y * 2 + 1) * surface->pitch);
            for (int x = 0; x < canvas->width; ++x) {
                Pixel p = src[x];
                // p = { (u8)(p.b & 0xE0), (u8)(p.g & 0xE0), (u8)(p.r & 0xE0), p.a };
                // p = { (u8)(p.b | p.b >> 3), (u8)(p.g | p.g >> 3), (u8)(p.r | p.r >> 3), p.a };
                dest1[x * 2 + 0] = p;
                dest1[x * 2 + 1] = p;
                dest2[x * 2 + 0] = p;
                dest2[x * 2 + 1] = p;
            }
        }
    } else if (scale == 3) {
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * src = canvas->pixels + y * canvas->pitch;
            Pixel * dest1 = (Pixel *) ((u8 *) surface->pixels + (y * 3 + 0) * surface->pitch);
            Pixel * dest2 = (Pixel *) ((u8 *) surface->pixels + (y * 3 + 1) * surface->pitch);
            Pixel * dest3 = (Pixel *) ((u8 *) surface->pixels + (y * 3 + 2) * surface->pitch);
            for (int x = 0; x < canvas->width; ++x) {
                Pixel p = src[x];
                // p = { (u8)(p.b & 0xC0), (u8)(p.g & 0xD0), (u8)(p.r & 0xD0), p.a };
                // p = { (u8)(p.b | p.b >> 3), (u8)(p.g | p.g >> 3), (u8)(p.r | p.r >> 3), p.a };
                dest1[x * 3 + 0] = p;
                dest1[x * 3 + 1] = p;
                dest1[x * 3 + 2] = p;
                dest2[x * 3 + 0] = p;
                dest2[x * 3 + 1] = p;
                dest2[x * 3 + 2] = p;
                dest3[x * 3 + 0] = p;
                dest3[x * 3 + 1] = p;
                dest3[x * 3 + 2] = p;
            }
        }
    } else if (scale == 4) {
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * src = canvas->pixels + y * canvas->pitch;
            Pixel * dest1 = (Pixel *) ((u8 *) surface->pixels + (y * 4 + 0) * surface->pitch);
            Pixel * dest2 = (Pixel *) ((u8 *) surface->pixels + (y * 4 + 1) * surface->pitch);
            Pixel * dest3 = (Pixel *) ((u8 *) surface->pixels + (y * 4 + 2) * surface->pitch);
            Pixel * dest4 = (Pixel *) ((u8 *) surface->pixels + (y * 4 + 3) * surface->pitch);
            for (int x = 0; x < canvas->width; ++x) {
                Pixel p = src[x];
                // p = { (u8)(p.b & 0xF0), (u8)(p.g & 0xF0), (u8)(p.r & 0xF0), p.a };
                // p = { (u8)(p.b | p.b >> 4), (u8)(p.g | p.g >> 4), (u8)(p.r | p.r >> 4), p.a };

                // p = { (u8)(p.b & 0xC0), (u8)(p.g & 0xC0), (u8)(p.r & 0xC0), p.a };
                // p = { (u8)(p.b | p.b >> 2), (u8)(p.g | p.g >> 2), (u8)(p.r | p.r >> 2), p.a };
                // p = { (u8)(p.b | p.b >> 4), (u8)(p.g | p.g >> 4), (u8)(p.r | p.r >> 4), p.a };

                // p = { (u8)(p.b & 0x80), (u8)(p.g & 0x80), (u8)(p.r & 0x80), p.a };
                // p = { (u8)(p.b | p.b >> 1), (u8)(p.g | p.g >> 1), (u8)(p.r | p.r >> 1), p.a };
                // p = { (u8)(p.b | p.b >> 2), (u8)(p.g | p.g >> 2), (u8)(p.r | p.r >> 2), p.a };
                // p = { (u8)(p.b | p.b >> 4), (u8)(p.g | p.g >> 4), (u8)(p.r | p.r >> 4), p.a };
                dest1[x * 4 + 0] = p;
                dest1[x * 4 + 1] = p;
                dest1[x * 4 + 2] = p;
                dest1[x * 4 + 3] = p;
                dest2[x * 4 + 0] = p;
                dest2[x * 4 + 1] = p;
                dest2[x * 4 + 2] = p;
                dest2[x * 4 + 3] = p;
                dest3[x * 4 + 0] = p;
                dest3[x * 4 + 1] = p;
                dest3[x * 4 + 2] = p;
                dest3[x * 4 + 3] = p;
                dest4[x * 4 + 0] = p;
                dest4[x * 4 + 1] = p;
                dest4[x * 4 + 2] = p;
                dest4[x * 4 + 3] = p;
            }
        }
    } else {
        for (int y = 0; y < canvas->height; ++y) {
            Pixel * row = canvas->pixels + y * canvas->pitch;
            for (int x = 0; x < canvas->width; ++x) {
                //write into dest pixels
                for (int yy = 0; yy < scale; ++yy) {
                    Pixel * dest =
                        (Pixel *)((u8 *)surface->pixels + (y * scale + yy) * surface->pitch);
                    for (int xx = 0; xx < scale; ++xx) {
                        dest[x * scale + xx] = row[x];
                    }
                }
            }
        }
    }
}
