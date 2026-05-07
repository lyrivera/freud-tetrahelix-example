"""
Generate a GSD file with two parallel Boerdijk-Coxeter (tetrahelix) chains:
one left-handed, one right-handed, separated by 10 particle diameters.
Periodic cubic simulation box.
"""

import math
import numpy as np
import gsd.hoomd

# ==================== User-Configurable Variables ====================
PARTICLE_RADIUS = 1.0          # Sphere radius R
NUM_PER_CHAIN = 20             # Particles per helix (total = 2 * NUM_PER_CHAIN)
BUFFER_SURFACE = 0.0           # Surface-to-surface separation between adjacent particles (0 = touching)
OUTPUT_FILENAME = "tetrahelix_dual.gsd"
# ==============================================================

# ==================== Boerdijk-Coxeter Helix Geometry ====================
# From Wikipedia: https://en.wikipedia.org/wiki/Boerdijk%E2%80%93Coxeter_helix
# Unit-edge tetrahedra (center-to-center distance between adjacent vertices = 1)
r_unit = 3 * math.sqrt(3) / 10    # Helix radius for unit edge length, ~0.5196
theta_right = math.acos(-2/3)     # Right-handed twist angle, ~2.300 rad (131.81°)
theta_left = -theta_right         # Left-handed twist angle
h_unit = 1 / math.sqrt(10)        # Z-increment per particle for unit edge length, ~0.3162
# ========================================================================

# ==================== Scaling ====================
# Center-to-center distance between adjacent particles
d_center = 2 * PARTICLE_RADIUS + BUFFER_SURFACE
# Scale factor to map unit helix to desired spacing
scale = d_center  # Unit helix has center-to-center distance 1, so scale by d_center
# =================================================

# ==================== Generate Positions ====================
positions = []
axis_offset = 10 * PARTICLE_RADIUS  # 10 diameters between axes: each axis is ±10R from origin (20R total)
# Generate helix symmetrically around z=0
z_start = -scale * (NUM_PER_CHAIN - 1) * h_unit / 2
for n in range(NUM_PER_CHAIN):
    z = z_start + scale * n * h_unit
    # Right-handed helix (x = +axis_offset)
    x_r = scale * r_unit * math.cos(n * theta_right) + axis_offset
    y_r = scale * r_unit * math.sin(n * theta_right)
    positions.append([x_r, y_r, z])
    # Left-handed helix (x = -axis_offset)
    x_l = scale * r_unit * math.cos(n * theta_left) - axis_offset
    y_l = scale * r_unit * math.sin(n * theta_left)
    positions.append([x_l, y_l, z])
positions = np.array(positions, dtype=np.float32)
# ============================================================

# ==================== Cubic Periodic Box ====================
# Account for particle radii in bounds
particle_diameter = 2 * PARTICLE_RADIUS
min_vals = positions.min(axis=0) - PARTICLE_RADIUS
max_vals = positions.max(axis=0) + PARTICLE_RADIUS

# Calculate required spans
Lx_needed = max_vals[0] - min_vals[0]
Ly_needed = max_vals[1] - min_vals[1]
Lz_needed = max_vals[2] - min_vals[2]

# Cubic box length: maximum of all spans, then multiply by 1.5
L = max(Lx_needed, Ly_needed, Lz_needed) * 1.5

# HOOMD boxes are centered at origin: [-L/2, L/2] in each dimension
# Particles are already generated around origin, so no shift needed
positions_shifted = positions
# ===========================================================

# ==================== Write GSD File ====================
with gsd.hoomd.open(OUTPUT_FILENAME, "w") as traj:
    frame = gsd.hoomd.Frame()
    frame.particles.N = len(positions_shifted)
    frame.particles.position = positions_shifted
    frame.particles.diameter = np.full(frame.particles.N, particle_diameter, dtype=np.float32)
    frame.particles.typeid = np.zeros(frame.particles.N, dtype=np.uint32)
    frame.particles.types = ["A"]
    # Cubic orthorhombic box: Lx=Ly=Lz=L, no tilt (0,0,0)
    frame.configuration.box = [L, L, L, 0.0, 0.0, 0.0]
    traj.append(frame)
    print("Note: Set periodic boundaries in your MD simulation when loading this GSD file.")

    # Store values for printing
    num_particles = frame.particles.N
    box_length = L

print(f"Generated {OUTPUT_FILENAME}")
print(f"Total particles: {num_particles}")
print(f"Cubic box length: {box_length:.4f}")
print(f"Helix axis separation: {2 * axis_offset:.1f} (10 diameters = {10 * particle_diameter:.1f})")
print()
print("Chain positions (verification):")
print(f"  Left chain X: [{positions_shifted[1::2][:,0].min():.2f}, {positions_shifted[1::2][:,0].max():.2f}]")
print(f"  Right chain X: [{positions_shifted[0::2][:,0].min():.2f}, {positions_shifted[0::2][:,0].max():.2f}]")
print(f"  Box X: [-{box_length/2:.2f}, {box_length/2:.2f}] (centered at origin)")
# Verify particles are inside box (HOOMD uses [-L/2, L/2])
box_half = box_length / 2
inside = np.all((positions_shifted >= -box_half) & (positions_shifted <= box_half), axis=1)
if np.all(inside):
    print(f"  All {num_particles} particles inside box: ✓")
else:
    print(f"  WARNING: {np.sum(~inside)} particles outside box!")
# ========================================================
