cd "$(dirname "$0")"

bob/bob -t base.ninja -o build.ninja
export NINJA_STATUS='[%t/%f] %es '
ninja || exit
rm build.ninja

./game
