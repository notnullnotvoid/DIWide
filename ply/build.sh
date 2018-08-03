clang -std=c++1z -stdlib=libc++ -I../lib -Wall -Ofast -o ply ply.cpp -lc++ || exit
./ply -i sword7.ply -o sword7.diw