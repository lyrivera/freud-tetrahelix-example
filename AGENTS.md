# AGENTS.md

## Project: Tetrahelix GSD Generator

### Purpose
Generate GSD v5 files with Boerdijk-Coxeter (tetrahelix) sphere packings for molecular dynamics simulations.

### Key Script
- `tetrahelix_dual_gsd.py`: Generates dual helices (left + right handed) in a cubic periodic box

### Dependencies
```bash
pip install numpy gsd  # gsd>=3.0
```

### Usage
Edit variables at top of `tetrahelix_dual_gsd.py`:
- `PARTICLE_RADIUS`: Sphere radius (default: 1.0)
- `NUM_PER_CHAIN`: Particles per helix (default: 10, total = 2×N)
- `BUFFER_SURFACE`: Surface-to-surface separation (default: 0.0 = touching)

Run: `python3 tetrahelix_dual_gsd.py`

Output: `tetrahelix_dual.gsd` (20 particles, cubic periodic box, L≈44.05 for defaults)

### Geometry Notes
- Uses verified Boerdijk-Coxeter helix formulas from Wikipedia
- Helix radius: r = 3√3/10 ≈ 0.5196 (for unit edge length)
- Twist angle: θ = arccos(-2/3) ≈ 131.81° per particle
- Z-increment: h = 1/√10 ≈ 0.3162 (for unit edge length)
- Helices separated by 10 particle diameters (20×radius between axes)
- Periodic boundaries set in simulation when loading (not in GSD file for v5 format)
