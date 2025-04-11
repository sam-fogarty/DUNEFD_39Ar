import uproot
import numpy as np
import matplotlib.pyplot as plt

npz_data = np.load('matched_hit_data_smallerProxCut_TwoOrThreePlanes.npz')
Z_npz = npz_data['Z']
Y_npz = npz_data['Y']
dt_npz = npz_data['DelT']

f = uproot.open('pdhd_028850_LowECLMatching_test.root')
tree = f["results"]
dt_root = tree["dt"].array()
Z_root = tree["charge_z"].array()
Y_root = tree["charge_y"].array()

bins=25
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(10,6))
axes[0].hist(dt_root, label='root', bins=bins, alpha=0.6)
axes[0].hist(dt_npz, label='npz', bins=bins, alpha=0.6)

axes[1].hist(Z_root, label='root', bins=bins, alpha=0.6)
axes[1].hist(Z_npz, label='npz', bins=bins, alpha=0.6)

axes[2].hist(Y_root, label='root', bins=bins, alpha=0.6)
axes[2].hist(Y_npz, label='npz', bins=bins, alpha=0.6)
axes[0].legend()
axes[1].legend()
axes[2].legend()
plt.show()

