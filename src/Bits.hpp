#ifndef BITS_HPP
#define BITS_HPP

struct Bits {
    using u8 = uint8_t;

    u8 * data;
    size_t len;
    size_t max;

    void init(size_t reserve = 512) {
        assert(reserve > 0);
        data = (TYPE *) malloc(reserve / 8);
        max = reserve;
    }
};

#endif
