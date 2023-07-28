import numpy as np
from numpy import loadtxt
import time
from collections import defaultdict
import csv

## read data
with open(f"coordinate.txt", "r") as f:
    info = f.readline().strip('\n').split(' ')
pt_data = loadtxt(f"coordinates.txt", comments="#", delimiter=" ", unpack=False)

coordinates = np.vstack([pt_data[:,0]*float(info[0]), 
                         pt_data[:,1]*float(info[1]), 
                         pt_data[:,2]*float(info[2])]).T
resid = pt_data[:,-1]
atomid = pt_data[:,3]

## clustering
def periodic_distance(p1, p2, bx, by, bz):
    dx = np.abs(p1[:, 0, None] - p2[:, 0])
    dy = np.abs(p1[:, 1, None] - p2[:, 1])
    dz = np.abs(p1[:, 2, None] - p2[:, 2])

    dx = np.minimum(dx, bx - dx)
    dy = np.minimum(dy, by - dy)
    dz = np.minimum(dz, bz - dz)

    return np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

def find_root(roots, node):
    while roots[node] != node:
        node = roots[node]
    return node

def union(roots, node1, node2):
    root1 = find_root(roots, node1)
    root2 = find_root(roots, node2)
    if root1 != root2:
        roots[root2] = root1

def group_points(points, threshold, bx, by, bz, batch_size=1000):
    groups = defaultdict(list)
    roots = {i: i for i in range(len(points))}
    num_points = len(points)

    # Calculate distance matrix in batches
    for i in range(0, num_points, batch_size):
        j_end = min(i + batch_size, num_points)
        dist_matrix_batch = periodic_distance(points[i:j_end], points, bx, by, bz)
        dist_matrix_batch[np.tril_indices_from(dist_matrix_batch, k=0)] = None
        close_pairs = np.where(dist_matrix_batch <= threshold)

        # Perform unions based on the mask
        for close_i, close_j in zip(close_pairs[0], close_pairs[1]):
            union(roots, i + close_i, close_j)

    # Compress roots to obtain final groups
    for i in range(num_points):
        root = find_root(roots, i)
        groups[root].append(i)

    # Sort groups based on the first member and sort members within each group
    sorted_groups = sorted(groups.items(), key=lambda x: x[1][0])   
    
    final_groups = []
    group_size = []
    for group_id, members in sorted_groups:
        group_size.append(len(members))
        final_groups.append((group_id, sorted(members)))

    return final_groups, group_size

start = time.time()
result, cluster_size = group_points(coordinates, 0.94, float(info[0]), float(info[1]), float(info[2]))
print(f'Time: {(time.time()-start):.3f} ses.')

## write files
with open(f"cluster_sizes.txt", "w") as f:
    for group_id, members in result:
        f.write(str(len(members)))
        f.write('\n')

clusters = np.zeros([len(result), max(cluster_size)])
clusters_resid = np.zeros([len(result), max(cluster_size)])

i = 0
for group_id, members in result:
    clusters[i][:len(members)] = np.array(members)+1
    clusters_resid[i][:len(members)] = resid[np.array(members)]
    i += 1
clusters = clusters.astype(int).tolist()
clusters_resid = clusters_resid.astype(int).tolist()

with open(f"clusters.txt", "w") as f:
    csv.writer(f, delimiter=' ').writerows(clusters)

with open(f"clusters_resid.txt", "w") as f:
    csv.writer(f, delimiter=' ').writerows(clusters_resid)