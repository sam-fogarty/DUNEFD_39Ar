import matplotlib.pyplot as plt
import numpy as np

data = np.load('matched_hit_data_smallerProxCut_ThreePlanesOnly.npz')
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(6,12))
#axes[0].hist(data['cluster_t'], bins=100, label='charge')
#axes[0].hist(data['pds_t'], bins=100, label='pds')
#print(len(np.unique(data['pds_t'])))
#axes[0].legend()
print(len(data['Z']))
bins=263
#mask = (data['Z'] > 400) & (data['Z'] < 460) & (data['APA'] != 3) & (data['Y'] > 0) & (data['Y'] < 50)
mask = np.ones(len(data['Z']), dtype=bool)
DelT = data['DelT'][mask]
APA = data['APA'][mask]
#mask = (data['Z'] > 390) & (data['Z'] < 460) \
#        & (data['Y'] > 0) & (data['Y'] < 80)
axes[0].hist(DelT[APA == 2], bins=bins, label='APA 2')
#axes[0].hist(data['X'][data['APA'] == 2], bins=200)
#axes[0].set_xlim(-360, -150)
#axes[0].set_xlim(150, 360)
#axes[0].set_xlabel('Z [cm]')
axes[0].set_ylabel('Counts')
#axes[0].set_title('APA 2')
axes[0].legend()

axes[1].hist(DelT[APA == 3], bins=bins, label='APA 3')
#axes[1].hist(data['X'][data['APA'] == 3], bins=200)
#axes[1].set_xlim(-360, -150)
#axes[1].set_xlim(150, 360)
#axes[1].set_xlabel('Z [cm]')
axes[1].set_ylabel('Counts')
#axes[1].set_title('APA 3')
axes[1].legend()

axes[2].hist(DelT[APA == 4], bins=bins, label='APA 4')
#axes[2].hist(data['X'][data['APA'] == 4], bins=200)
#axes[2].set_xlim(-360, -150)
#axes[2].set_xlim(150, 360)
axes[2].set_xlabel(r'Charge t - PDS t [$\mu$s]')
axes[2].set_ylabel('Counts')
axes[2].legend()
#axes[2].set_title('APA 4')

fig.suptitle('Charge Light Matched Low Energy Clusters (Three Plane Matches)')

S_Range = (250, 500)
B_Range = (1000, 2000)

T = np.sum((DelT < S_Range[1]) & (DelT > S_Range[0]))
B = np.sum((DelT < B_Range[1]) & (DelT > B_Range[0]))/4
S = T - B
print('Purity = ', S/T)

#print(len(data['X']))
#axes[0].hist2d(data['Z'], data['Y'], bins=100)
#axes[1].hist2d(data['cluster_z'][data['cluster_apa'] > 2], data['cluster_y'][data['cluster_apa'] > 2], bins=50)
#axes[0].set_ylabel('Y [cm]')
#axes[1].set_ylabel('Y [cm]')
#axes[0].set_xlabel('Z [cm]')
#axes[1].set_xlabel('Z [cm]')


plt.show()
