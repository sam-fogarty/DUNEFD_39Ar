import numpy as np
import matplotlib.pyplot as plt
import random

data = np.load('all_hit_data.npz')
pds_t = data['pds_t'][0:100]
print('loaded data')
y = np.ones_like(pds_t)

clusters = []
eps = 0.5
points_sorted = pds_t
curr_point = points_sorted[0]
curr_cluster = [curr_point]
for point in points_sorted[1:]:
    if point <= curr_point + eps:
        curr_cluster.append(point)
    else:
        clusters.append(curr_cluster)
        curr_cluster = [point]
    curr_point = point
clusters.append(curr_cluster)
#print(clusters)
min_t = min(pds_t)
fig, axes = plt.subplots(figsize=(10,6))
for cluster in clusters:
    y = np.ones_like(cluster)
    axes.plot(cluster-min_t, y, 'o', markersize=random.randint(3, 8))
#axes.plot(pds_t-min(pds_t), y, 'bo', markersize=3)
#axes.set_xlim(,2000)
plt.show()
