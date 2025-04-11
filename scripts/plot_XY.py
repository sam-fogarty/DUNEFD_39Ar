import matplotlib.pyplot as plt
import numpy as np

data = np.load('matched_hit_data.npz')
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12,10))
#axes[0].hist(data['cluster_t'], bins=100, label='charge')
#axes[0].hist(data['pds_t'], bins=100, label='pds')
#print(len(np.unique(data['pds_t'])))
#axes[0].legend()
print(len(data['Z']))
#mask = (data['Z'] > 390) & (data['Z'] < 460) & (data['Y'] > 0) & (data['Y'] < 80)
mask = np.ones(len(data['X']), dtype=bool)
X = data['X'][mask]
Y = data['Y'][mask]
APA = data['APA'][mask]
axes[0].hist2d(X[APA == 2], Y[APA==2], bins=75)
axes[1].hist2d(X[APA == 3], Y[APA==3], bins=75)
axes[2].hist2d(X[APA == 4], Y[APA==4], bins=75)

#axes[0].set_xlim(230, 460)
#axes[1].set_xlim(0, 230)
#axes[2].set_xlim(230, 460)
#axes[0].set_ylim(10, 605)
#axes[1].set_ylim(10, 605)
#axes[2].set_ylim(10, 605)
#axes[0].hist(data['X'], bins=50)
#print(len(data['X']))
#axes[0].hist2d(data['Z'], data['Y'], bins=100)
#axes[1].hist2d(data['cluster_z'][data['cluster_apa'] > 2], data['cluster_y'][data['cluster_apa'] > 2], bins=50)
#axes[0].set_ylabel('Y [cm]')
#axes[1].set_ylabel('Y [cm]')
#axes[0].set_xlabel('Z [cm]')
#axes[1].set_xlabel('Z [cm]')


plt.show()
