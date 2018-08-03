cl /std:c++14 /MP /Ox /ISDL2 /Ilib lib/stb_image.cpp src/blit.cpp src/gif.cpp src/main.cpp /link /out:game.exe /SUBSYSTEM:CONSOLE || exit /b
del *.obj
game.exe
