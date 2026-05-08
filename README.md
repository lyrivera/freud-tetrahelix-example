# Tetrahelix GSD Generator & Analysis

Generate and analyze GSD v5 files with Boerdijk-Coxeter (tetrahelix) sphere packings for molecular dynamics simulations.

## Scripts

### 1. Generate Tetrahelix (`tetrahelix_dual_gsd.py`)
Creates a GSD file with two parallel tetrahelix chains (left + right handed) in a cubic periodic box.

**Config (top of script):**
- `PARTICLE_RADIUS`: Sphere radius (default: 1.0)
- `NUM_PER_CHAIN`: Particles per helix (default: 20, total = 2×N)
- `BUFFER_SURFACE`: Surface-to-surface separation (default: 0.0 = touching)
- `NOISE_SIGMA`: Gaussian noise sigma (default: 0.05)
- `NUM_RANDOM_PARTICLES`: Extra random particles (default: 0)
- `RANDOM_SEED`: RNG seed for reproducibility (default: None)

**Run:** `python tetrahelix_dual_gsd.py`

**Output:** `tetrahelix_dual.gsd` (40 particles for defaults, cubic periodic box, L≈36.12)

---

### 2. Analyze Tetrahelix (`analyze_tetrahelix.py`)
Identifies tetrahedral environments using freud's `EnvironmentMotifMatch` and writes a new GSD with type differentiation.

**Config (top of script):**
- `PARTICLE_RADIUS`, `BUFFER_SURFACE`: Must match GSD generation
- `R_MAX_FACTOR`: Multiplier for neighbor cutoff (default: 1.5)
- `THRESHOLD`: Matching threshold (default: 0.5)

**Run:** `python analyze_tetrahelix.py`

**Output:** `tetrahelix_analyzed.gsd` with types:
- Type A = other particles
- Type L = left-handed helix
- Type R = right-handed helix

---

## Geometry Notes
- Helix radius: r = 3√3/10 ≈ 0.5196 (unit edge length)
- Twist angle: θ = arccos(-2/3) ≈ 131.81° per particle
- Z-increment: h = 1/√10 ≈ 0.3162 (unit edge length)
- Helices separated by 10 particle diameters (20×radius between axes)
- Particles form tetrahedral environments with edge length = 2R + buffer

---

## Visualization
Open `tetrahelix_analyzed.gsd` in OVITO and color particles by type:
- Type A: gray (other)
- Type L: blue (left-handed)
- Type R: red (right-handed)
