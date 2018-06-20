#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "common.hpp"
#include "List.hpp"
#include "FileBuffer.hpp"
#include "math.hpp"

//match window surface's byte order for faster blit
typedef struct Pixel {
    u8 b, g, r, a;
} Color;

struct Tex {
    Pixel * pixels;
    Pixel * normals;
    int width, height;
    u32 wmask, hmask;
};

struct Vertex {
    Vec3 p, n, t, b;
    float u, v;
};

typedef u16 idx_t;
struct Triangle {
    idx_t v0, v1, v2;
};

struct Model {
    size_t vertCount;
    Vertex * vertices;
    void * verts; //allocated within the game

    size_t triangleCount;
    Triangle * triangles;

    Tex tex;
};

int main(int argc, char ** argv) {
    const char * infile = nullptr;
    const char * outfile = nullptr;

    enum ArgType {
        ARG_NONE, ARG_INFILE, ARG_OUTFILE,
    };

    //parse command line arguments
    ArgType type = ARG_NONE;
    for (int i : range(1, argc)) {
        if (strcmp(argv[i], "-i") == 0) {
            type = ARG_INFILE;
        } else if (strcmp(argv[i], "-o") == 0) {
            type = ARG_OUTFILE;
        } else if (type == ARG_INFILE) {
            infile = argv[i];
            type = ARG_NONE;
        } else if (type == ARG_OUTFILE) {
            outfile = argv[i];
            type = ARG_NONE;
        }
    }



    //argument checking
    if (!infile) {
        printf("ERROR: No input file specified!\n");
        return 1;
    }
    if (!outfile) {
        printf("DEBUG: No output file specified, continuing in debug mode\n");
    }



    char * file = read_entire_file(infile);

    //split file into lines
    //TOOD: migrate this into a helper split() in common.hpp?
    List<char *> lines = {};
    char * line = strtok(file, "\r\n");
    while (line != nullptr) {
        lines.add(line);
        line = strtok(nullptr, "\r\n");
    }



    enum ParseState {
        PARSE_NONE, PARSE_HEADER, PARSE_VERT, PARSE_TRI,
    };

    List<Vertex> verts = {};
    List<Triangle> tris = {};
    int vcount = 0;
    int fcount = 0;

    ParseState state = PARSE_HEADER;
    for (char * line : lines) {
        if (state == PARSE_HEADER) {
            char * tok = strtok(line, " \t");

            if (!strcmp(tok, "format")) {
                assert(!strcmp(strtok(nullptr, " \t"), "ascii"));
            } else if (!strcmp(tok, "comment")) {
                //nop (ignore comment lines)
            } else if (!strcmp(tok, "element")) {
                char * etype = strtok(nullptr, " \t");
                char * count = strtok(nullptr, " \t");
                assert(etype && count);

                if (!strcmp(etype, "vertex")) {
                    vcount = atoi(count);
                    verts = create_list<Vertex>(vcount);
                } else if (!strcmp(etype, "face")) {
                    fcount = atoi(count);
                    tris = create_list<Triangle>(fcount);
                } else {
                    assert(0 && "element type must be vertex or face");
                }
            //TODO: configurable properties?
            } else if (!strcmp(tok, "end_header")) {
                state = PARSE_VERT;
            }
        } else if (state == PARSE_VERT) {
            assert(vcount);
            verts.add({ vec3(atof(strtok(line   , " \t")),
                             atof(strtok(nullptr, " \t")),
                             atof(strtok(nullptr, " \t"))),
                        vec3(atof(strtok(nullptr, " \t")),
                             atof(strtok(nullptr, " \t")),
                             atof(strtok(nullptr, " \t"))),
                        vec3(), vec3(),
                        (float) atof(strtok(nullptr, " \t")),
                        (float) atof(strtok(nullptr, " \t")) });
            --vcount;
            if (!vcount) {
                state = PARSE_TRI;
            }
        } else if (state == PARSE_TRI) {
            assert(fcount);

            int n = atoi(strtok(line, " \t"));
            assert(n >= 3);

            List<int> indices = create_list<int>(n);
                for (int _ UNUSED : range(n)) {
                    indices.add(atoi(strtok(nullptr, " \t")));
                }

                for (int i : range(2, n)) {
                    tris.add({ (idx_t) indices[0], (idx_t) indices[i - 1], (idx_t) indices[i] });
                }
            indices.finalize();

            --fcount;
            if (!fcount) {
                state = PARSE_NONE;
            }
        }
    }

    //calculate tangents and bitangents
    for (Triangle t : tris) {
        Vertex & v0 = verts[t.v0];
        Vertex & v1 = verts[t.v1];
        Vertex & v2 = verts[t.v2];

        Vec3 e1 = v1.p - v0.p;
        Vec3 e2 = v2.p - v0.p;

        float du1 = v1.u - v0.u;
        float dv1 = v1.v - v0.v;
        float du2 = v2.u - v0.u;
        float dv2 = v2.v - v0.v;

        float f = 1.0f / (du1 * dv2 - du2 * dv1);
        Vec3 tangent = vec3(f * (dv2 * e1.x - dv1 * e2.x),
                            f * (dv2 * e1.y - dv1 * e2.y),
                            f * (dv2 * e1.z - dv1 * e2.z));
        Vec3 bitangent = vec3(f * (du1 * e2.x - du2 * e1.x),
                              f * (du1 * e2.y - du2 * e1.y),
                              f * (du1 * e2.z - du2 * e1.z));
        v0.t += tangent;
        v1.t += tangent;
        v2.t += tangent;
        v0.b += bitangent;
        v1.b += bitangent;
        v2.b += bitangent;
    }

    //normalize, orthogonalize, and re-normalize the tangent and bitangent
    for (Vertex & v : verts) {
        v.t = noz(v.t);
        v.b = noz(v.b);
        // v.t = noz(v.t - dot(v.n, v.t) * v.n); //orthogonalize tangent wrt normal
        // v.b = noz(v.b - dot(v.n, v.b) * v.n); //orthogonalize bitangent wrt normal
        // v.b = noz(v.b - dot(v.t, v.b) * v.t); //orthogonalize bitangent wrt tangent
    }

    //cleanup
    lines.finalize();
    free(file);



    FileBuffer buf = {};
    auto model = buf.write<Model>({ verts.len, nullptr, nullptr, tris.len, nullptr });

    buf.update(model + offsetof(Model, vertices), buf.size());
    for (Vertex v : verts) {
        buf.write(v);
    }

    buf.update(model + offsetof(Model, triangles), buf.size());
    for (Triangle t : tris) {
        buf.write(t);
    }

    //cleanup
    verts.finalize();
    tris.finalize();



    if (outfile) {
        FILE * fp = fopen(outfile, "wb");
        assert(fp);
        fwrite(buf.block, buf.size(), 1, fp);
        fclose(fp);
    } else {
        //TODO: debug print
    }
}
