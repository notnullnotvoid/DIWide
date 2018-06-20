#ifndef COMMON_HPP
#define COMMON_HPP

#include <cstdint>
#include <cstddef>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cassert>

#include <initializer_list>

////////////////////////////////////////////////////////////////////////////////
/// TYPES                                                                    ///
////////////////////////////////////////////////////////////////////////////////

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

typedef int8_t i8;
typedef int16_t i16;
typedef int32_t i32;
typedef int64_t i64;

typedef uint8_t b8;
typedef uint16_t b16;
typedef uint32_t b32;
typedef uint64_t b64;

typedef float f32;
typedef double f64;
typedef long double f80;

typedef unsigned int uint;

template <typename T1, typename T2>
struct Pair {
    T1 first;
    T2 second;
};

////////////////////////////////////////////////////////////////////////////////
/// RANGES                                                                   ///
////////////////////////////////////////////////////////////////////////////////

struct _range {
    int first;
    int last;
    int step;
    struct iterator {
        int value;
        int step;
        int operator*() { return value; }
        bool operator!=(iterator & other) { return value != other.value; }
        iterator& operator++() { value += step; return *this; }
    };
    iterator begin() { return { first, step }; }
    iterator end() { return { last, step }; }
};
inline _range range(int last) { return { 0, last, 1 }; }
inline _range range(int first, int last) { return { first, last, first > last? -1 : 1 }; }

//this uses a bogus comparison operator to iterate a null-terminated string
//in a range for loop without needing to call strlen() first
struct _strange {
    char * str;
    struct iterator {
        char * ch;
        char & operator*() { return *ch; }
        bool operator!=(iterator & other) { return *ch != '\0' && *other.ch != '\0'; }
        iterator& operator++() { ++ch; return *this; }
    };
    iterator begin() { return { str }; }
    iterator end() { return { str }; }
};
struct _const_strange {
    const char * str;
    struct iterator {
        const char * ch;
        const char & operator*() { return *ch; }
        bool operator!=(iterator & other) { return *ch != '\0' && *other.ch != '\0'; }
        iterator& operator++() { ++ch; return *this; }
    };
    iterator begin() { return { str }; }
    iterator end() { return { str }; }
};
inline _strange range(char * str) { return { str }; }
inline _const_strange range(const char * str) { return { str }; }

////////////////////////////////////////////////////////////////////////////////
/// STRING UTILS                                                             ///
////////////////////////////////////////////////////////////////////////////////

inline char * dup(const char * src, int len) {
    char * ret = (char *) malloc(len + 1);
    strncpy(ret, src, len);
    ret[len] = '\0';
    return ret;
}

inline char * dup(const char * src, char * end) {
    return dup(src, end - src);
}

inline char * dup(const char * src) {
    return dup(src, strlen(src));
}

inline bool same(const char * first, const char * second) {
    return !strcmp(first, second);
}

////////////////////////////////////////////////////////////////////////////////
/// DEFER                                                                    ///
////////////////////////////////////////////////////////////////////////////////

//a simplified version of https://gist.github.com/p2004a/045726d70a490d12ad62
template <typename F>
struct _defer {
    F f;
    _defer(F f) : f(f) {}
    ~_defer() { f(); }
};

struct {
    template <typename F>
    _defer<F> operator<<(F f) { return _defer<F>(f); };
} _deferrer;

#define TOKENPASTE2(x, y) x ## y
#define TOKENPASTE(x, y) TOKENPASTE2(x, y)
#define defer auto TOKENPASTE(__deferred_lambda_call, __COUNTER__) = _deferrer << [&]

////////////////////////////////////////////////////////////////////////////////
/// ???                                                                      ///
////////////////////////////////////////////////////////////////////////////////

//NOTE: the correct behavior of this function is unfortunately not guaranteed by the standard
inline char * read_entire_file(const char * filepath) {
    FILE * f = fopen(filepath, "rb");
    assert(f);
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);  //same as rewind(f);

    char * string = (char *) malloc(fsize + 1);
    fread(string, fsize, 1, f);
    fclose(f);

    string[fsize] = 0;

    return string;
}

template <typename VAL, typename KEY>
VAL match_pair(std::initializer_list<Pair<VAL, KEY>> pairs, KEY key, VAL defaultVal) {
    for (auto it = pairs.begin(); it != pairs.end(); ++it) {
        if (it->second == key) {
            return it->first;
        }
    }
    return defaultVal;
}

template <typename TYPE>
bool one_of(std::initializer_list<TYPE> list, TYPE key) {
    for (auto it = list.begin(); it != list.end(); ++it) {
        if (*it == key) {
            return true;
        }
    }
    return false;
}

template <typename TYPE>
inline void swap(TYPE & a, TYPE & b) {
    TYPE temp = a;
    a = b;
    b = temp;
}

#define ARR_SIZE(x) (sizeof(x) / sizeof((x)[0]))

#if defined(__GNUC__)
#  define UNUSED __attribute__ ((unused))
#elif defined(_MSC_VER)
#  define UNUSED __pragma(warning(suppress:4100))
#else
#  define UNUSED
#endif

#endif
