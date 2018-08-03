#ifndef GIF_HPP
#define GIF_HPP

#include "List.hpp"
#include "blit.hpp" //TODO: get rid of this dependency so it will be easy to drop into new projects

struct RawFrame {
    Pixel * base;
    Pixel * pixels;
    int pitch;
};

void save_gif(int width, int height, List<RawFrame> rawFrames, int centiSeconds);

#endif
