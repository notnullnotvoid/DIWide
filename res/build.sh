cd "$(dirname "$0")"

convert "src/CobbleStone_01_BC.png" -separate -delete 3 -swap 0,2 \( "src/CobbleStone_01_R.png" -separate \) -combine -flip "cobblestone1_c.png"
convert "src/CobbleStone_01_N.png" -separate -delete 3 -swap 0,2 \( "src/CobbleStone_01_M.png" -separate \) -combine -flip "cobblestone1_n.png"

convert "src/BrokenTiles_01_BC.png" -separate -delete 3 -swap 0,2 \( "src/BrokenTiles_01_R.png" -separate \) -combine -flip "tiles1_c.png"
convert "src/BrokenTiles_01_N.png" -separate -delete 3 -swap 0,2 \( "src/BrokenTiles_01_M.png" -separate \) -combine -flip "tiles1_n.png"

convert "src/Sword7_diffuse.png" -separate -delete 3 -swap 0,2 \( "src/Sword7_roughness.png" -separate \) -combine -flip "sword7_c.png"
convert "src/Sword7_normal.png" -separate -delete 3 -swap 0,2 \( "src/Sword7_metallic.png" -separate \) -combine -flip "sword7_n.png"

convert "src/Sword9_diffuse.png" -separate -delete 3 -swap 0,2 \( "src/Sword9_roughness.png" -separate \) -combine -flip "sword9_c.png"
convert "src/Sword9_normal.png" -separate -delete 3 -swap 0,2 \( "src/Sword9_metallic.png" -separate \) -combine -flip "sword9_n.png"

convert "src/Tomato_diffuse.png" -separate -delete 3 -swap 0,2 \( "src/Tomato_MetalSmooth.png" -separate -delete 0 -negate \) -combine -flip "tomato_c.png"
convert "src/Tomato_normal.png" -separate -delete 3 -swap 0,2 \( "src/Tomato_MetalSmooth.png" -separate \) -combine -flip "tomato_n.png"

convert "src/MetalSheets_01_BC.png" -separate -delete 3 -swap 0,2 \( "src/MetalSheets_01_R.png" -separate \) -combine -flip "metal1_c.png"
convert "src/MetalSheets_01_N.png" -separate -delete 3 -swap 0,2 \( "src/MetalSheets_01_M.png" -separate \) -combine -flip "metal1_n.png"
