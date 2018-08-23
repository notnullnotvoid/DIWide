# DIWide
## Features
- perspective-correct vertex attribute interpolation
- correct near-plane clipping
- viewport and backface culling
- fully gamma-correct pipeline
- bilinear texture filtering
- diffuse, normal, and roughness textures
- high-resolution shadow mapping
- simplified physically-based lighting model
- post-processing bit-depth reduction
- built-in custom GIF exporter

## How to build
- **macOS**: First install GCC 8.x and alias it to `gcc-8`. This alias will be created automatically if you install GCC via [Homebrew](https://brew.sh/) (`brew install gcc`). SDL2 must also be installed (`brew install sdl2`). After that, run `mac-build.sh`.
- **Windows**: Visual Studio command line tools must be installed. If you don't have Visual Studio installed already, you can get an installer for the standalone build tools [here](https://visualstudio.microsoft.com/downloads/#build-tools-for-visual-studio-2017). After that, run `win-build.bat` from the 64-bit developer command prompt.
- **Advanced build (macOS only for now)**: In addition to GCC and SDL2, you will need to install ninja (`brew install ninja`). First, build bob via `bob/buid.sh`. This only needs to happen once. Then, run `build.sh` in the root of the repository for subsequent builds.

Currently, only x86-64 CPUs are supported. GCC is used instead of clang on macOS because it produces slightly (around 10-20%) faster code than clang. Right now the renderer runs significantly (around 2x) slower on Windows, because MSVC does a poorer job of optimizing our particular code. I'll add an option for building with GCC on Windows soon. I will also try to get the advanced build working on Windows.

## Screenshots
![](https://i.imgur.com/lATL7sO.gif)
![](https://i.imgur.com/aFhFUMK.gif)
![](https://i.imgur.com/9QAmmLm.png)
![](https://i.imgur.com/p627VsA.png)
