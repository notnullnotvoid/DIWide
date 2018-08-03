gcc-8 -std=c++14 -Ilib -ISDL2 -Ofast -o game -lSDL2 lib/stb_image.cpp src/blit.cpp src/gif.cpp src/main.cpp || exit
./game
