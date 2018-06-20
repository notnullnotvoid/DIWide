rm ./ply
clang -std=c++1z -stdlib=libc++ -Wall -Ofast -o ply ply.cpp -lc++
./ply -i sword7.ply -o sword7.diw