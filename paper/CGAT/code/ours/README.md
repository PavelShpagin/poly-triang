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

Algorithm selection:
- `--algo chain` (default): chain-based sweep with hybrid fallback when \(k \ge \max(8, n/8)\)
- `--algo chain_only`: chain-based sweep only (no fallback)
- `--algo linked`: edge-based sweep with linked representation (fallback path)


