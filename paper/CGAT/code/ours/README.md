# Our implementation (CGAT artifact)

## Build

```bash
g++ -O3 -std=c++17 reflex_cli.cpp -o reflex_cli
```

## Run

```bash
./reflex_cli --input polygon.poly --output out.tri --algo chain
```

The program prints a CSV-like line including:
- `k_count`: number of local maxima (paper’s \(k\))
- `reflex_count`: number of reflex vertices (paper’s \(r\))
- `time_ms`


