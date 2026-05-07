#!/usr/bin/env python3
"""
Analyze tetrahelix GSD file using freud's EnvironmentCluster.
Identifies particles with tetrahedral environments and writes a new GSD
with type differentiation for visualization.
"""

import numpy as np
import gsd.hoomd
from freud import box as freud_box
from freud import environment as freud_env
from freud import locality

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

bx = freud_box.Box(Lx=box_L, Ly=box_L, Lz=box_L)
# =====================================================

# ==================== Run EnvironmentCluster ====================
d_center = 2 * PARTICLE_RADIUS + BUFFER_SURFACE
r_max = d_center * R_MAX_FACTOR

env_cluster = freud_env.EnvironmentCluster()
env_cluster.compute(
    system=frame,
    threshold=THRESHOLD,
    cluster_neighbors={'r_max': r_max},
    env_neighbors={'r_max': r_max},
    registration=True
)

cluster_idx = env_cluster.cluster_idx
num_clusters = env_cluster.num_clusters
print(f"Found {num_clusters} clusters")
for i in range(num_clusters):
    count = np.sum(cluster_idx == i)
    print(f"  Cluster {i}: {count} particles")
# =====================================================

# ==================== Identify Tetrahelix Cluster ====================
# The tetrahelix particles should form the largest cluster
# (or clusters with 6 neighbors at distance d_center)
unique, counts = np.unique(cluster_idx, return_counts=True)
largest_cluster = unique[np.argmax(counts)]
print(f"\nLargest cluster: {largest_cluster} with {counts.max()} particles")

# Mark particles in largest cluster as "tetrahelix"
matches = (cluster_idx == largest_cluster)
# =====================================================

# ==================== Write NEW GSD with Type Differentiation ====================
# Type 0 = "A" (other), Type 1 = "B" (tetrahelix)
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
