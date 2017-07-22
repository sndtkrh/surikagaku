# Surikagaku-tokuron (数理科学特論)

## Overview
This is a software to compute simplicial homology and cubical homology.

単体複体のホモロジーと方体集合のホモロジーを計算するソフトウェアです．

## Building
To compile C++ code using this library, you need to use C++ compiler compatible with C++11.
GCC G++ 6.3 or later is recomended.
For example, you can compile `examples/simplicial.cpp` as follows:

```
$ g++ -std=c++11 ./examples/simplicial.cpp -o simplicial
```

このライブラリを用いたC++コードをコンパイルするためにはC++11に対応したコンパイラが必要です．
G++6.3以上を推奨します．
例えば`examples/simplicial.cpp`をコンパイルするためには次のようにします．

```
$ g++ -std=c++11 ./examples/simplicial.cpp -o simplicial
```

## Homology Computation

### Simplicial Homology
To compute simplicial homology,
you include `src/simplicial_homology.hpp`.
You can see examples in `examples/simplicial.cpp`.

単体複体のホモロジーを計算するには
`src/simplicial_homology.hpp`をインクルードします．
具体的な例は`examples/simplicial.cpp`にあります．

### Cubical Homology
To compute cubical homology,
you include `src/cubical_homology.hpp`.
You can see examples in `examples/cubical.cpp`.

方体集合のホモロジーを計算するには
`src/cubical_homology.hpp`をインクルードします．
具体的な例は`examples/cubical.cpp`にあります．