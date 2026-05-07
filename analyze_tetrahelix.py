#!/usr/bin/env python3
"""
Analyze tetrahelix GSD file using freud's EnvironmentMotifMatch.
Identifies particles in tetrahedral environments and writes a new GSD
with type differentiation for visualization.
"""

import numpy as np
import gsd.hoomd
import freud.box
import freud.environment
import freud.locality

# ==================== Configuration ====================
PARTICLE_RADIUS = 1.0
BUFFER_SURFACE = 0.0
GSD_FILENAME = "tetrahelix_dual.gsd"
OUTPUT_FILENAME = "tetrahelix_analyzed.gsd"
R_MAX_FACTOR = 1.5
THRESHOLD = 0.3  # 30% of d_center (2.0) = 0.6
# =====================================================

# ==================== Load GSD File ====================
with gsd.hoomd.open(GSD_FILENAME, 'r') as traj:
    frame = traj[0]
    positions = frame.particles.position
    box_L = frame.configuration.box[0]

bx = freud.box.Box(Lx=box_L, Ly=box_L, Lz=box_L)
# =====================================================

# ==================== Define Tetrahedron Motif ====================
d_center = 2 * PARTICLE_RADIUS + BUFFER_SURFACE

# Motif: 6 neighbor vectors for tetrahelix environment
# User-provided motif (already has correct distances for d_center=2.0)
motif = np.array([
    [-1.0,        -np.sqrt(3) / 3,        -2 * np.sqrt(6) / 3],
    [ 1.0,        -np.sqrt(3) / 3,        -2 * np.sqrt(6) / 3],
    [ 0.0,         2 * np.sqrt(3) / 3,    -2 * np.sqrt(6) / 3],
    [ 5 / 3,       5 * np.sqrt(3) / 9,    -2 * np.sqrt(6) / 9],
    [ 1 / 9,       31 * np.sqrt(3) / 27,   2 * np.sqrt(6) / 27],
    [ 32 / 27,     38 * np.sqrt(3) / 81,   46 * np.sqrt(6) / 81],
], dtype=np.float32)

print(f"Motif (6 neighbor vectors, edge lengths should be ~{d_center:.1f}):")
for i, v in enumerate(motif):
    print(f"  v{i}: dist={np.linalg.norm(v):.3f}")
# =====================================================

# ==================== Run EnvironmentMotifMatch ====================
r_max = d_center * R_MAX_FACTOR  # e.g., 1.5 * 2.0 = 3.0

emmatch = freud.environment.EnvironmentMotifMatch()
emmatch.compute(
    system=frame,
    motif=motif,
    threshold=THRESHOLD,
    env_neighbors={'r_max': r_max, 'num_neighbors': 10},
    registration=True  # User specified True
)

matches = emmatch.matches
# =====================================================

# ==================== Write NEW GSD with Type Differentiation ====================
# Type 0 = "A" (non-matching), Type 1 = "B" (matching tetrahelix)
typeids = np.zeros(len(positions), dtype=np.uint32)
typeids[matches] = 1

with gsd.hoomd.open(OUTPUT_FILENAME, "w") as f:
    new_frame = gsd.hoomd.Frame()
    new_frame.particles.N = len(positions)
    new_frame.particles.position = positions
    new_frame.particles.diameter = np.full(len(positions), 2*PARTICLE_RADIUS, dtype=np.float32)
    new_frame.particles.typeid = typeids
    new_frame.particles.types = ["A", "B"]
    new_frame.configuration.box = [box_L, box_L, box_L, 0.0, 0.0, 0.0]
    f.append(new_frame)

print(f"\nWritten {OUTPUT_FILENAME}")
print(f"  Type A (other): {np.sum(~matches)} particles")
print(f"  Type B (tetrahelix): {np.sum(matches)} particles")
# =====================================================

# ==================== Console Output ====================
print(f"\nAnalysis Results:")
print(f"  Total particles: {len(positions)}")
print(f"  Tetrahelix particles: {np.sum(matches)}")
if np.sum(matches) > 0:
    print(f"  Tetrahelix indices: {np.where(matches)[0]}")

# Save indices to file
np.savetxt('tetrahelix_matching_particles.txt', np.where(matches)[0], fmt='%d')
print(f"  Saved indices to: tetrahelix_matching_particles.txt")
# =====================================================
