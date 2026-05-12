#!/usr/bin/env python3
"""
Analyze tetrahelix GSD file using freud's EnvironmentMotifMatch.
Identifies particles in tetrahedral environments and writes a new GSD
with type differentiation for visualization.
Types: A=other, L=left-handed, R=right-handed
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
THRESHOLD = 0.5  # Increased to catch particles with fewer neighbors
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
# Unit-scaled, then multiplied by d_center for generality
motif = d_center * np.array([
    [-0.5,        -np.sqrt(3) / 6,        -np.sqrt(6) / 3],
    [ 0.5,        -np.sqrt(3) / 6,        -np.sqrt(6) / 3],
    [ 0.0,         np.sqrt(3) / 3,        -np.sqrt(6) / 3],
    [ 5 / 6,       5 * np.sqrt(3) / 18,   -np.sqrt(6) / 9],
    [ 1 / 18,      31 * np.sqrt(3) / 54,   np.sqrt(6) / 27],
    [ 16 / 27,     19 * np.sqrt(3) / 81,   23 * np.sqrt(6) / 81],
], dtype=np.float32)

print(f"Motif (6 neighbor vectors, edge lengths should be ~{d_center:.1f}):")
for i, v in enumerate(motif):
    print(f"  v{i}: dist={np.linalg.norm(v):.3f}")
# =====================================================

# ==================== Run EnvironmentMotifMatch (Right-Handed) ====================
r_max = d_center * R_MAX_FACTOR

# Create neighbor list for all particles
aq = freud.locality.AABBQuery(bx, positions)
nlist = aq.query(positions, {'r_max': r_max, 'num_neighbors': 6}).toNeighborList()

emmatch = freud.environment.EnvironmentMotifMatch()
emmatch.compute(
    system=frame,
    motif=motif,
    threshold=THRESHOLD,
    env_neighbors={'r_max': r_max, 'num_neighbors': 6},
    registration=True
)

matches_R = emmatch.matches  # Right-handed matches
# =====================================================

# ==================== Match Left-Handed Chain ====================
# Convert positions to right-handed coordinates (mirror x-axis)
# Left-handed helix becomes right-handed when x -> -x
positions_mirrored = positions.copy()
positions_mirrored[:, 0] = -positions_mirrored[:, 0]

# Create new frame with mirrored positions
frame_mirrored = gsd.hoomd.Frame()
frame_mirrored.particles.N = len(positions_mirrored)
frame_mirrored.particles.position = positions_mirrored
frame_mirrored.configuration.box = [box_L, box_L, box_L, 0.0, 0.0, 0.0]

# Create neighbor list for mirrored positions
aq_mirrored = freud.locality.AABBQuery(bx, positions_mirrored)
nlist_mirrored = aq_mirrored.query(positions_mirrored, {'r_max': r_max, 'num_neighbors': 6}).toNeighborList()

# Run EnvironmentMotifMatch on mirrored positions
emmatch_left = freud.environment.EnvironmentMotifMatch()
emmatch_left.compute(
    system=frame_mirrored,
    motif=motif,
    threshold=THRESHOLD,
    env_neighbors={'r_max': r_max, 'num_neighbors': 6},
    registration=True
)

matches_L = emmatch_left.matches  # Left-handed matches (from mirrored coords)
# =====================================================

# ==================== Assign Types: A=other, L=left, R=right ====================
typeids = np.zeros(len(positions), dtype=np.uint32)

# Right-handed matches get type 2 ("R")
typeids[matches_R] = 2

# Left-handed matches get type 1 ("L")
# If a particle matches both, prioritize R (user specified)
if np.any(matches_R & matches_L):
    print("Note: Some particles match both L and R, prioritizing R")
    typeids[matches_L] = 1  # Will be overwritten by R if both match
    typeids[matches_R] = 2  # R takes precedence
else:
    typeids[matches_L] = 1
# =====================================================

# ==================== Write NEW GSD with Type Differentiation ====================
with gsd.hoomd.open(OUTPUT_FILENAME, "w") as f:
    new_frame = gsd.hoomd.Frame()
    new_frame.particles.N = len(positions)
    new_frame.particles.position = positions
    new_frame.particles.diameter = np.full(len(positions), 2*PARTICLE_RADIUS, dtype=np.float32)
    new_frame.particles.typeid = typeids
    new_frame.particles.types = ["A", "L", "R"]  # A=other, L=left, R=right
    new_frame.configuration.box = [box_L, box_L, box_L, 0.0, 0.0, 0.0]
    f.append(new_frame)

print(f"\nWritten {OUTPUT_FILENAME}")
print(f"  Type A (other): {np.sum(typeids == 0)} particles")
print(f"  Type L (left-handed): {np.sum(typeids == 1)} particles")
print(f"  Type R (right-handed): {np.sum(typeids == 2)} particles")
# =====================================================

# ==================== Console Output ====================
print(f"\nAnalysis Results:")
print(f"  Total particles: {len(positions)}")
print(f"  Left-handed matches: {np.sum(matches_L)}")
print(f"  Right-handed matches: {np.sum(matches_R)}")

if np.sum(matches_R) > 0:
    print(f"  R indices: {np.where(matches_R)[0]}")
if np.sum(matches_L) > 0:
    print(f"  L indices: {np.where(matches_L)[0]}")
# =====================================================
