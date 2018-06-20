//this program walks the root and lib/ directories of the project to generate a .ninja file
//yes this is a pun on "bob the builder" don't @ me

//bob/build.sh builds bob.
//you only need to build bob once before you can build the rest of the project

//TODO: don't overwrite ninja file if it's identical to the generated output

//TODO: walk folder structures recursively
//TODO: generate .app bundles for macOS
//TODO: test and make sure this works with symlinks

//wishful thinking section
//TODO: build for Android
//TODO: build for iOS

#include <stdlib.h>
#include <stdio.h>
#include <dirent.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <time.h>
#include <stdarg.h>

#include "common.hpp"
#include "List.hpp"

bool verbose;

enum FileType {
    FILE_NONE,
    FILE_TEXTURE,
    FILE_TEX_LINEAR,
    FILE_MODEL,
    FILE_CPPLIB,
    FILE_CLIB,
    FILE_SRC,
};

struct File {
    FileType type;
    char * dir;
    char * name;
    char * ext;
};

List<File> add_files(List<File> files, FileType type, const char * dir, const char * ext) {
    //open directory
    DIR * dp = opendir(dir);
    if (!dp) {
        if (verbose) {
            printf("WARNING: Could not open directory %s, assuming empty\n", dir);
        }
        return files;
    }

    //scan for source files
    dirent * ep = readdir(dp);
    while (ep) {
        char * dot = strrchr(ep->d_name, '.');
        if (dot && !strcmp(dot + 1, ext)) {
            files.add({ type, dup(dir), dup(ep->d_name, dot), dup(ext) });
        }
        ep = readdir(dp);
    }

    closedir(dp);

    return files;
}

enum ArgType {
    ARG_NONE,
    ARG_ROOT,
    ARG_TEMP,
    ARG_OUTPUT,
};

//NOTE: this is no longer necessary, but may be useful in a future project...
//      should I maybe include it in common.hpp?
char * dsprintf(char * buf, const char * fmt, ...) {
    size_t len = buf? strlen(buf) : 0;
    va_list args1, args2;
    va_start(args1, fmt);
    va_copy(args2, args1);
    buf = (char *) realloc(buf, len + vsnprintf(nullptr, 0, fmt, args1) + 1);
    vsprintf(buf + len, fmt, args2);
    va_end(args1);
    va_end(args2);
    return buf;
}

int main(int argc, char ** argv) {
    const char * root = ".";
    const char * temp = "template.ninja";
    const char * output = "build.ninja";

    //parse command line arguments
    ArgType type = ARG_NONE;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-v")) {
            verbose = true;
        } else if (!strcmp(argv[i], "-r")) {
            type = ARG_ROOT;
        } else if (!strcmp(argv[i], "-t")) {
            type = ARG_TEMP;
        } else if (!strcmp(argv[i], "-o")) {
            type = ARG_OUTPUT;
        } else if (type == ARG_ROOT) {
            root = argv[i];
            type = ARG_NONE;
        } else if (type == ARG_TEMP) {
            temp = argv[i];
            type = ARG_NONE;
        } else if (type == ARG_OUTPUT) {
            output = argv[i];
            type = ARG_NONE;
        }
    }

    //change directory
    chdir(root);

    //read template ninja file
    char * part = read_entire_file(temp);
    if (!part) {
        printf("Error while reading file %s\n", temp);
        return 1;
    }

    //gather and classify files
    List<File> files = create_list<File>(1);
    files = add_files(files, FILE_TEXTURE, "asset", "png");
    files = add_files(files, FILE_TEXTURE, "asset", "tga");
    files = add_files(files, FILE_TEX_LINEAR, "asset/linear", "png");
    files = add_files(files, FILE_TEX_LINEAR, "asset/linear", "tga");
    files = add_files(files, FILE_MODEL, "asset", "ply");
    files = add_files(files, FILE_CPPLIB, "lib", "cpp");
    files = add_files(files, FILE_CLIB, "lib", "c");
    files = add_files(files, FILE_SRC, ".", "cpp");
    files = add_files(files, FILE_SRC, "src", "cpp");

    //write ninja build file
    //paste template into output file
    FILE * out = fopen(output, "w");
    fprintf(out, "%s\n", part);
    free(part);

    //write asset conversion commands
    for (File f : files) {
        if (f.type == FILE_TEXTURE) {
            fprintf(out, "build res/src/%s.png: img %s/%s.%s\n", f.name, f.dir, f.name, f.ext);
        } else if (f.type == FILE_TEX_LINEAR) {
            fprintf(out, "build res/src/%s.png: img_lin %s/%s.%s\n", f.name, f.dir, f.name, f.ext);
        } else if (f.type == FILE_MODEL) {
            fprintf(out, "build res/%s.diw: ply %s/%s.%s\n", f.name, f.dir, f.name, f.ext);
        }
    }
    fprintf(out, "\n");

    //write compile commands
    for (File f : files) {
        auto rule = match_pair<const char *>({
            { "lib", FILE_CPPLIB },
            { "clib", FILE_CLIB },
            { "src", FILE_SRC },
        }, f.type, nullptr);

        if (rule) {
            fprintf(out, "build $builddir/%s.o: %s %s/%s.%s\n",
                f.name, rule, f.dir, f.name, f.ext);
        }
    }

    //write link command
    fprintf(out, "\nbuild game: link");
    for (File f : files) {
        if (one_of({ FILE_CPPLIB, FILE_CLIB, FILE_SRC }, f.type)) {
            fprintf(out, " $builddir/%s.o", f.name);
        }
    }
    fprintf(out, "\n");

    fclose(out);

/*
    //windows build
    //TODO: specify windows build script location (and whether to generate it at all?)
    FILE * win = fopen("build.bat", "w");
    fprintf(win, "cl /Ox /ISDL2 /Isoloud /DWITH_SDL2");
    for (File f : files) {
        if (one_of({ FILE_CPPLIB, FILE_CLIB, FILE_SRC }, f.type)) {
            fprintf(win, " %s/%s.%s", f.dir, f.name, f.ext);
        }
    }
    fprintf(win, " /link /out:game.exe /SUBSYSTEM:CONSOLE\r\ndel *.obj\r\ngame.exe\r\n");

    fclose(win);
*/

    return 0;
}
