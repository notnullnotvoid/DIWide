clang -std=c++11 -stdlib=libc++ -I../lib -Wall -Ofast -o bob bob.cpp -lc++ || exit
./bob -r .. -t base.ninja -o bob/out.ninja -v
