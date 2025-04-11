import ROOT
import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm
import h5py

copyFiles = True

run = 29218
folder_descr = 'test2'
if copyFiles:
    folder = f'/pnfs/dune/scratch/users/sfogarty/JOBSUB/CLMatching/{run}_{folder_descr}/'
else:
    folder=f'{run}_{folder_descr}/'

filepaths_temp = os.listdir(folder)
filepaths = []
for file in filepaths_temp:
    if file.endswith('.root'):
        if copyFiles:
            file = 'root://fndca1.fnal.gov:1094//pnfs/fnal.gov/usr/' + (folder + file).strip('/pnfs/')
            filepaths.append(file)
        else:
            file = os.path.join(folder, file)
            filepaths.append(file)

firstSave = True
charge_z, charge_y, charge_apa, charge_energy, charge_nplanes, charge_dt = [],[],[],[],[],[]
light_amplitude, light_z, light_y, light_apa, light_ch, light_nhits = [],[],[],[],[],[]
for i,filepath in tqdm(enumerate(filepaths)):
    if copyFiles:
        print(f'Copying file {i}/{len(filepaths)}')
        os.system(f'mkdir -p 29218_{folder_descr}')
        if not os.path.exists(f'./29218_{folder_descr}/'+os.path.basename(filepath)):
            os.system(f'xrdcp {filepath} ./29218_{folder_descr}/')
    else:
        if not os.path.exists(filepath):
            continue
        try:
            f = ROOT.TFile.Open(filepath, "READ")
        except:
            print(f"Error: Could not open file {filepath}")
            continue
        if not f or f.IsZombie():
            print(f"Error: Could not open file {filepath}")
            continue
        LightOffsetsTree = f.Get("Offsets_Light_CLMatched")
        #LightOffsetsTree.Print()
        ChargeOffsetsTree = f.Get("Offsets_Charge_CLMatched")
        #ChargeOffsetsTree.Print()

        light_indices = []
        for entry in LightOffsetsTree:
            start_index = entry.start_index
            stop_index = entry.stop_index
            groupsize = entry.groupsize
            light_indices.append((start_index, stop_index))
        charge_indices = []
        for entry in ChargeOffsetsTree:
            start_index = entry.start_index
            stop_index = entry.stop_index
            groupsize = entry.groupsize
            charge_indices.append((start_index, stop_index))
            #print(f"start_index: {start_index}, stop_index: {stop_index}")
        
        PDSHitsTree = f.Get("IsolatedPDSGroups_Matched")
        ChargeClustersTree = f.Get("IsolatedChargeClusters_Matched")

        for j in range(len(charge_indices)):
            t_avg = 0
            light_amplitude_temp = []
            light_z_temp = []
            light_y_temp = []
            light_apa_temp = []
            light_ch_temp = []
            for index in range(light_indices[j][0], light_indices[j][1]):
                PDSHitsTree.GetEntry(index)
                t_avg += PDSHitsTree.t
                light_amplitude_temp.append(PDSHitsTree.amplitude)
                light_z_temp.append(PDSHitsTree.z)
                light_y_temp.append(PDSHitsTree.y)
                light_apa_temp.append(PDSHitsTree.apa)
                light_ch_temp.append(PDSHitsTree.ch)
            if light_indices[j][1] - light_indices[j][0] == 1:
                light_amplitude_temp += [-1,-1]
                light_z_temp += [-1,-1]
                light_y_temp += [-1,-1]
                light_apa_temp += [-1,-1]
                light_ch_temp += [-1,-1]
            elif light_indices[j][1] - light_indices[j][0] == 2:
                light_amplitude_temp.append(-1)
                light_z_temp.append(-1)
                light_y_temp.append(-1)
                light_apa_temp.append(-1)
                light_ch_temp.append(-1)
            
            t_avg /= (light_indices[j][1] - light_indices[j][0])

            if charge_indices[j][1] - charge_indices[j][0] == 1:
                light_nhits.append(light_indices[j][1] - light_indices[j][0])
                light_amplitude.append(light_amplitude_temp)
                light_z.append(light_z_temp)
                light_y.append(light_y_temp)
                light_apa.append(light_apa_temp)
                light_ch.append(light_ch_temp)
                for index in range(charge_indices[j][0], charge_indices[j][1]):
                    ChargeClustersTree.GetEntry(index)
                    charge_dt.append(ChargeClustersTree.t - t_avg)
                    charge_z.append(ChargeClustersTree.z)
                    charge_y.append(ChargeClustersTree.y)
                    charge_apa.append(ChargeClustersTree.apa)
                    charge_energy.append(ChargeClustersTree.E)
                    charge_nplanes.append(ChargeClustersTree.nplanes)

        f.Close()

    if not copyFiles:
        if firstSave and i != 0:
            light_amplitude = np.array(light_amplitude)
            light_z = np.array(light_z)
            light_y = np.array(light_y)
            light_apa = np.array(light_apa)
            light_ch = np.array(light_ch)
            light_nhits = np.array(light_nhits)
            charge_z = np.array(charge_z)
            charge_y = np.array(charge_y)
            charge_apa = np.array(charge_apa)
            charge_energy = np.array(charge_energy)
            charge_nplanes = np.array(charge_nplanes)
            charge_dt = np.array(charge_dt)
            firstSave = False
            with h5py.File(f'run_{run}_light_charge_{folder_descr}.h5', 'w') as f:
                f.create_dataset('light_amplitude', data=light_amplitude, maxshape=(None, light_amplitude.shape[1]))
                f.create_dataset('light_z', data=light_z, maxshape=(None, light_z.shape[1]))
                f.create_dataset('light_y', data=light_y, maxshape=(None, light_y.shape[1]))
                f.create_dataset('light_apa', data=light_apa, maxshape=(None, light_apa.shape[1]))
                f.create_dataset('light_ch', data=light_ch, maxshape=(None, light_ch.shape[1]))
                f.create_dataset('light_nhits', data=light_nhits, maxshape=(None,))
                f.create_dataset('charge_z', data=charge_z, maxshape=(None,))
                f.create_dataset('charge_y', data=charge_y, maxshape=(None,))
                f.create_dataset('charge_apa', data=charge_apa, maxshape=(None,))
                f.create_dataset('charge_energy', data=charge_energy, maxshape=(None,))
                f.create_dataset('charge_nplanes', data=charge_nplanes, maxshape=(None,))
                f.create_dataset('charge_dt', data=charge_dt, maxshape=(None,))
            charge_z, charge_y, charge_apa, charge_energy, charge_nplanes, charge_dt = [],[],[],[],[],[]
            light_amplitude, light_z, light_y, light_apa, light_ch, light_nhits = [],[],[],[],[],[]
        elif not firstSave and i != 0:
            light_amplitude = np.array(light_amplitude)
            light_z = np.array(light_z)
            light_y = np.array(light_y)
            light_apa = np.array(light_apa)
            light_ch = np.array(light_ch)
            light_nhits = np.array(light_nhits)
            charge_z = np.array(charge_z)
            charge_y = np.array(charge_y)
            charge_apa = np.array(charge_apa)
            charge_energy = np.array(charge_energy)
            charge_nplanes = np.array(charge_nplanes)
            charge_dt = np.array(charge_dt)
            with h5py.File(f'run_{run}_light_charge_{folder_descr}.h5', 'a') as f:
                f['light_amplitude'].resize((f['light_amplitude'].shape[0] + light_amplitude.shape[0]), axis=0)
                f['light_amplitude'][-len(light_amplitude):] = light_amplitude
                f['light_z'].resize((f['light_z'].shape[0] + light_z.shape[0]), axis=0)
                f['light_z'][-len(light_z):] = light_z
                f['light_y'].resize((f['light_y'].shape[0] + light_y.shape[0]), axis=0)
                f['light_y'][-len(light_y):] = light_y
                f['light_apa'].resize((f['light_apa'].shape[0] + light_apa.shape[0]), axis=0)
                f['light_apa'][-len(light_apa):] = light_apa
                f['light_ch'].resize((f['light_ch'].shape[0] + light_ch.shape[0]), axis=0)
                f['light_ch'][-len(light_ch):] = light_ch
                f['light_nhits'].resize((f['light_nhits'].shape[0] + light_nhits.shape[0]), axis=0)
                f['light_nhits'][-len(light_nhits):] = light_nhits
                f['charge_z'].resize((f['charge_z'].shape[0] + charge_z.shape[0]), axis=0)
                f['charge_z'][-len(charge_z):] = charge_z
                f['charge_y'].resize((f['charge_y'].shape[0] + charge_y.shape[0]), axis=0)
                f['charge_y'][-len(charge_y):] = charge_y
                f['charge_apa'].resize((f['charge_apa'].shape[0] + charge_apa.shape[0]), axis=0)
                f['charge_apa'][-len(charge_apa):] = charge_apa
                f['charge_energy'].resize((f['charge_energy'].shape[0] + charge_energy.shape[0]), axis=0)
                f['charge_energy'][-len(charge_energy):] = charge_energy
                f['charge_nplanes'].resize((f['charge_nplanes'].shape[0] + charge_nplanes.shape[0]), axis=0)
                f['charge_nplanes'][-len(charge_nplanes):] = charge_nplanes
                f['charge_dt'].resize((f['charge_dt'].shape[0] + charge_dt.shape[0]), axis=0)
                f['charge_dt'][-len(charge_dt):] = charge_dt
            charge_z, charge_y, charge_apa, charge_energy, charge_nplanes, charge_dt = [],[],[],[],[],[]
            light_amplitude, light_z, light_y, light_apa, light_ch, light_nhits = [],[],[],[],[],[]
            
     
"""
if not copyFiles:
    plt.figure(figsize=(10, 6)) 

    bins = 100
    binsize = 10
    Range = np.array([0, binsize * bins])-100
    hist, bin_edges = np.histogram(all_charge_DeltaT, bins=bins, range=Range)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    errors = np.sqrt(hist)
    plt.step(bin_centers, hist, where='mid', color='blue', label='Charge Delta T')
    plt.errorbar(bin_centers, hist, yerr=errors, fmt='o', color='blue', markersize=2, label='Error')
    plt.xlabel('Delta T [us]')
    plt.ylabel('Counts')
    plt.title('Histogram of Charge Delta T')
    plt.legend()
    plt.show()
"""
"""
plt.figure(figsize=(10, 6))
plt.hist(light_groupsize_all, bins=5, range=(0, 5), color='blue', alpha=0.7)
plt.xlabel('Group Size')
plt.ylabel('Counts')
plt.title('Histogram of Charge Group Size')
plt.grid()
plt.show()
"""

