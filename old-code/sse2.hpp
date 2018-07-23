#ifndef SSE2_HPP
#define SSE2_HPP

void draw_triangle_sse2_plain(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]);
void draw_triangle_sse2(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]);
void draw_triangle_sse2_struct(Canvas * canvas, ZBuffer * shadow, Tex * tex, Vert tri[3]);

#endif // SSE2_HPP
