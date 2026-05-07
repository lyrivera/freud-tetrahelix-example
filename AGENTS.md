# AGENTS.md

## Project: Tetrahelix GSD Generator & Analysis

### Purpose
Generate and analyze GSD v5 files with Boerdijk-Coxeter (tetrahelix) sphere packings for molecular dynamics simulations.

### Key Scripts
- `tetrahelix_dual_gsd.py`: Generates dual helices (left + right handed) in a cubic periodic box
- `analyze_tetrahelix.py`: Analyzes GSD files using freud's EnvironmentCluster to identify tetrahedral environments

### Dependencies
```bash
pip install numpy gsd freud-analysis  # gsd>=3.0, freud-analysis>=3.0
```

### Usage: Generate
Edit variables at top of `tetrahelix_dual_gsd.py`:
- `PARTICLE_RADIUS`: Sphere radius (default: 1.0)
- `NUM_PER_CHAIN`: Particles per helix (default: 10, total = 2×N)
- `BUFFER_SURFACE`: Surface-to-surface separation (default: 0.0 = touching)

Run: `python3 tetrahelix_dual_gsd.py`

Output: `tetrahelix_dual.gsd` (40 particles for defaults, cubic periodic box, L≈36.12)

### Usage: Analyze
Edit variables at top of `analyze_tetrahelix.py`:
- `PARTICLE_RADIUS`, `BUFFER_SURFACE`: Must match GSD generation
- `R_MAX_FACTOR`: Multiplier for neighbor cutoff (default: 1.5)
- `THRESHOLD`: Matching threshold (default: 0.3)

Run: `python3 analyze_tetrahelix.py`

Output:
- `tetrahelix_analyzed.gsd`: New GSD with Type A (non-helix) and Type B (helix) for visualization
- `tetrahelix_matching_particles.txt`: Indices of matching particles

### Geometry Notes
- Uses verified Boerdijk-Coxeter helix formulas from Wikipedia
- Helix radius: r = 3√3/10 ≈ 0.5196 (for unit edge length)
- Twist angle: θ = arccos(-2/3) ≈ 131.81° per particle
- Z-increment: h = 1/√10 ≈ 0.3162 (for unit edge length)
- Helices separated by 10 particle diameters (20×radius between axes)
- Particles form tetrahedral environments with edge length = 2R + buffer
- Analysis uses freud EnvironmentCluster to find particles with similar tetrahedral environments
