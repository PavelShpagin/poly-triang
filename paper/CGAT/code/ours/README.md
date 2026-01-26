# Our implementation

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
- `--algo chain` (default): chain-based sweep (no fallback)
- `--algo chain_only`: alias of `chain` (kept for backwards compatibility)
- `--algo linked`: edge-based sweep with linked representation
- `--algo seidel`: Seidel baseline (C code; for debugging/benchmarking only)


