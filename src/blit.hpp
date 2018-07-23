// HEAVYWEIGHT CANVAS OPERATIONS DEFINED HERE

#ifndef BLIT_HPP
#define BLIT_HPP

#include "common.hpp"

#include <SDL_surface.h>

//match window surface's byte order for faster blit
typedef struct Pixel {
    u8 b, g, r, a;
} Color;

struct Canvas {
	Pixel * base; //without the base address, we can't free the memory
    Pixel * pixels;
    size_t pixelBytes;
    float * depth;
    int width, height;
    int pitch, zpitch; //number of pixels, NOT number of bytes!
};

struct ZBuffer {
	float * depth;
	float scale, inv;
	int width, height;
	u32 wmask, hmask;
};

Canvas create_canvas(int width, int height, int margin);
ZBuffer create_depth_buffer(int width, int height);
void debug_depth_blit(SDL_Surface * canvas, ZBuffer * depth);
void fast_scaled_blit(SDL_Surface * surface, Canvas * canvas, int scale);

#endif // BLIT_HPP
