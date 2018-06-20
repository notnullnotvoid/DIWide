rm ./bob
clang -std=c++1z -stdlib=libc++ -Wall -Ofast -o bob bob.cpp -lc++
./bob -r .. -t base.ninja -o bob/out.ninja -v