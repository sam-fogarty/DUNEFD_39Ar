import matplotlib.pyplot as plt
import numpy as np
import uproot

#data = np.load('matched_hit_data_smallerProxCut_ThreePlanesOnly.npz')
fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(12,12))

f = uproot.open('pdhd_028113_LowECLMatching_r2.root')
tree = f["results"]
charge_x = tree["charge_x"].array()
dt = tree["dt"].array() - 256
apa = tree["APA"].array()
charge_nplanes = tree["charge_nplanes"].array()
pds_amplitude = tree["pds_amplitude"].array()
charge_energy = tree["charge_energy"].array()
pds_channel = tree["pds_ch"].array()
charge_peaktime = tree["charge_peaktime"].array()
print(charge_x)
#all_charge_peaktime = f["isolated_charge"]["peaktime"].array()
#all_pds_ch = f["isolated_pds"]["ch"].array()
#all_pds_amplitude = f["isolated_pds"]["amplitude"].array()
#tree_2 = f["results_twoPDShit"]
#dt_twoHit = tree_2['dt'].array()
#pds_amplitude_twoHit_1 = tree_2['pds_amplitude_1'].array()
#pds_amplitude_twoHit_2 = tree_2['pds_amplitude_2'].array()
#charge_nplanes_twoHit = tree_2['charge_nplanes'].array()
#print('# of two hit events = ', len(dt_twoHit))
#print('# of two hit events after applying cuts = ', np.sum((pds_amplitude_twoHit_1 > 120) & (pds_amplitude_twoHit_2 > 120) & (charge_nplanes_twoHit == 3)))
f.close()
#from tqdm import tqdm
#x = []
#for i in tqdm(range(len(charge_x))):
#    x.append(charge_x[i])
S_Range = (0, 150)
B_Range = (-200, -100)
amplitude_cut = 120
energy_cut=15
print(charge_peaktime)
mask = (charge_nplanes == 3) & (pds_amplitude > amplitude_cut) & (charge_energy > energy_cut)#& (charge_peaktime > 500)#& (pds_channel == 89) & (pds_amplitude > amplitude_cut) & (charge_energy > energy_cut) 

T = np.sum((dt[mask] < S_Range[1]) & (dt[mask] > S_Range[0]))
B = np.sum((dt[mask] < B_Range[1]) & (dt[mask] > B_Range[0]))/((B_Range[1] - B_Range[0])/(S_Range[1]-S_Range[0]))
S = T - B
print('3 view:')
print('Total T = ', T)
print('B = ', B)

purity_1 = S/T
print('purity = ', purity_1)
mask2 = (charge_nplanes > 1) #& (charge_peaktime > 500) #& (pds_amplitude > amplitude_cut) & (charge_energy > energy_cut)

T = np.sum((dt[mask2] < S_Range[1]) & (dt[mask2] > S_Range[0]))
B = np.sum((dt[mask2] < B_Range[1]) & (dt[mask2] > B_Range[0]))/((B_Range[1] - B_Range[0])/(S_Range[1]-S_Range[0]))
S = T - B
print('2 and 3 view:')
print('Total T = ', T)
print('B = ', B)
print('total signal = ', S)
print('Purity, all views = ', S/T)
purity_2 = S/T

axes[0][0].hist(dt[mask], bins=60)
axes[0][1].hist(dt[mask2], bins=60)
#axes[0][0].fill_between([S_Range[0], S_Range[1]], 0, 150, color='gray', alpha=0.4, label=f'Est. Purity = {purity_1:.2f}')
axes[0][0].set_xlabel(r'$\Delta t$ [cm]')
#axes[0][1].set_xlabel(r'$\Delta t$ [cm]')

axes[0][0].legend()
axes[0][1].legend()
factor = 5E-3
R = 0.67
Wions = 23.6*10E-6
calib = Wions/(factor*R)
#axes[1][0].hist(charge_energy[(charge_nplanes == 2) & (dt > 255) & (dt < 450)]*calib, range=(0,2),bins=50, color='b', alpha=0.6, label='2-view matches')
#axes[1][0].hist(charge_energy[mask]*calib, range=(0,3),bins=50, color='r', alpha=0.6, label='3-view matches')
axes[1][0].hist(pds_amplitude, bins=50)
axes[1][0].set_xlabel('Cluster Energy (~MeV)')
axes[0][0].set_title(f'3-view clusters; pds ampl. > {amplitude_cut}; energy > {(energy_cut*calib):.3f} MeV')
#axes[0][1].set_title(f'2 and 3-view clusters; pds ampl. > {amplitude_cut}; energy > {(energy_cut*calib):.3f} MeV')
#axes[1][0].set_xlim(0, 3)
axes[1][0].legend()
#print('len(all_pds_ch) = ', len(all_pds_ch))
axes[1][1].set_title('Channel counts without the amplitude and energy cuts')
#unique_ch, ch_counts = np.unique(all_pds_ch, return_counts=True)
#axes[1][1].step(unique_ch, ch_counts)
#axes[1][1].set_xlabel('PDS Offline Channel')
#print('pds channel with most counts = ', unique_ch[np.argmax(ch_counts)])
#axes[1][0].hist(all_charge_peaktime, bins=50)
#axes[2].hist(charge_energy, bins=1500)
#axes[0].set_xlim(0, 750)
#axes[1].set_xlim(0, 750)
#axes[2].set_xlim(0, 100)
#axes[0].hist(data['X'][data['APA'] == 2], bins=200)
#axes[0].set_xlim(-360, -150)
#axes.set_ylabel('Counts')
#axes[0].set_xlim(150, 360)
#axes[0].set_xlabel('Z [cm]')
#axes[0].set_ylabel('X reco [cm]')
#axes[0].set_title('APA 2')
#axes[0].legend()

#axes[1].hist(data['X'][data['APA'] == 3]-5, bins=200, label='APA 3')
#axes[1].hist(data['X'][data['APA'] == 3], bins=200)
#axes[1].set_xlim(-360, -150)
#axes[1].set_xlim(150, 360)
#axes[1].set_xlabel('Z [cm]')
#axes[1].set_ylabel('X reco [cm]')
#axes[1].set_ylabel('Counts')
#axes[1].set_title('APA 3')
#axes[1].legend()

#axes[2].hist(data['X'][data['APA'] == 4]-5, bins=200, label='APA 4')
#axes[2].hist(data['X'][data['APA'] == 4], bins=200)
#axes[2].set_xlim(-360, -150)
#axes[2].set_xlim(150, 360)
#axes[2].set_xlabel('X reco [cm]')
#axes[2].set_ylabel('Counts')
#axes[2].legend()
#axes[2].set_title('APA 4')
#fig.suptitle('Charge Light Matched Low Energy Clusters (Three Plane Matches)')
#print(len(data['X']))
#axes[0].hist2d(data['Z'], data['Y'], bins=100)
#axes[1].hist2d(data['cluster_z'][data['cluster_apa'] > 2], data['cluster_y'][data['cluster_apa'] > 2], bins=50)
#axes[0].set_ylabel('Y [cm]')
#axes[1].set_ylabel('Y [cm]')
#axes[0].set_xlabel('Z [cm]')
#axes[1].set_xlabel('Z [cm]')


plt.show()
