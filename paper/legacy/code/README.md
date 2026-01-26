# RefleX triangulation (paper package)

This folder contains a **self-contained** C++ reference implementation used in the paper.

## Build

From `paper/code/`:

```bash
g++ -O3 -std=c++17 test.cpp -o test
```

## Run

```bash
./test ../../polygons/benchmark/random_10000.poly
```

The program prints:

- `n`: number of vertices
- `triangles`: number of output triangles
- `expected`: `n-2`
- `reflex`: number of reflex vertices

## Input format

`.poly` text file:

```
n
x0 y0
x1 y1
...
```

Vertices must describe a **simple polygon** in CCW or CW order (we normalize internally).


