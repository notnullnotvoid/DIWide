cl /Ox /ISDL2 /Ilib lib/stb_image.cpp src/blit.cpp src/main.cpp /link /out:game.exe /SUBSYSTEM:CONSOLE
del *.obj
game.exe
