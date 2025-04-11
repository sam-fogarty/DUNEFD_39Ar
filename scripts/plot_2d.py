import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import uproot

#data = np.load('matched_hit_data_smallerProxCut_ThreePlanesOnly.npz')
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(13,7), sharey=True)
#axes[0].hist(data['cluster_t'], bins=100, label='charge')
#axes[0].hist(data['pds_t'], bins=100, label='pds')
#print(len(np.unique(data['pds_t'])))
#axes[0].legend()
#f = uproot.open('pdhd_028086_LowECLMatching_r2.root')
f = uproot.open('pdhd_028850_LowECLMatching_test2.root')
tree = f["results"]
dt = tree["dt"].array() - 256
charge_z = tree["charge_z"].array()
charge_y = tree["charge_y"].array()
APA = tree["APA"].array()
charge_nplanes = tree["charge_nplanes"].array()
pds_amplitude = tree["pds_amplitude"].array()
charge_energy = tree["charge_energy"].array()
f.close()
#mask = (data['Z'] > 390) & (data['Z'] < 460) & (data['Y'] > 0) & (data['Y'] < 80)
#mask = np.ones(len(data['Z']), dtype=bool)

#mask = (data['Z'] > 400) & (data['Z'] < 460) & (data['Y'] > 0) & (data['Y'] < 50)
#mask = (data['DelT'] > 250) & (data['DelT'] < 500)
Z = np.array(charge_z)#data['Z']#[~mask]
Y = np.array(charge_y)#data['Y']#[~mask]
APA = np.array(APA)#data['APA']#[~mask]
bins=[30,50]
Range_apa2 = [[230, 460], [0, 605]]
Range_apa3 = [[0, 230], [0, 605]]
Range_apa4 = [[230, 460], [0, 605]]
mask = (charge_nplanes == 3) #& (pds_amplitude > 120) & (charge_energy > 15) & (dt > 0) & (dt < 60)
# Plot for APA == 2
h2d_apa2 = axes[0].hist2d(Z[(APA == 2) & mask], Y[(APA == 2) & mask], bins=bins, range=Range_apa2, norm=LogNorm())
axes[0].set_xlabel('Z [cm]', fontsize=8)
axes[0].set_ylabel('Y [cm]', fontsize=8)
plt.colorbar(h2d_apa2[3], ax=axes[0])

# Plot for APA == 3
h2d_apa3 = axes[1].hist2d(Z[(APA == 3) & mask], Y[(APA == 3) & mask], bins=bins, range=Range_apa3, norm=LogNorm())
axes[1].set_xlabel('Z [cm]', fontsize=8)
axes[1].set_ylabel('Y [cm]', fontsize=8)
plt.colorbar(h2d_apa3[3], ax=axes[1])

# Plot for APA == 4
h2d_apa4 = axes[2].hist2d(Z[(APA == 4) & mask], Y[(APA == 4) & mask], bins=bins, range=Range_apa4, norm=LogNorm())
axes[2].set_xlabel('Z [cm]', fontsize=8)
axes[2].set_ylabel('Y [cm]', fontsize=8)
plt.colorbar(h2d_apa4[3], ax=axes[2])

#plt.tight_layout()

#axes[0].hist2d(Z[APA == 2], Y[APA==2], bins=75)
#axes[1].hist2d(Z[APA == 3], Y[APA==3], bins=75)
#axes[2].hist2d(Z[APA == 4], Y[APA==4], bins=75)
#axes[0].set_xlabel('Z [cm]')
#axes[1].set_xlabel('Z [cm]')
#axes[2].set_xlabel('Z [cm]')
#axes[0].set_ylabel('Y [cm]')

axes[0].set_title('APA 2')
axes[1].set_title('APA 3')
axes[2].set_title('APA 4')

axes[0].set_xlim(230, 460)
axes[1].set_xlim(0, 230)
axes[2].set_xlim(230, 460)
axes[0].set_ylim(10, 605)
axes[1].set_ylim(10, 605)
axes[2].set_ylim(10, 605)

fig.suptitle('Charge Light Matched Low Energy Clusters (Three Plane Matches)')
#axes[0].hist(data['X'], bins=50)
#print(len(data['X']))
#axes[0].hist2d(data['Z'], data['Y'], bins=100)
#axes[1].hist2d(data['cluster_z'][data['cluster_apa'] > 2], data['cluster_y'][data['cluster_apa'] > 2], bins=50)
#axes[0].set_ylabel('Y [cm]')
#axes[1].set_ylabel('Y [cm]')
#axes[0].set_xlabel('Z [cm]')
#axes[1].set_xlabel('Z [cm]')


plt.show()
