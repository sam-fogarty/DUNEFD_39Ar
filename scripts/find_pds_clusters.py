import numpy as np
from waffles.input import raw_hdf5_reader as rhr
import waffles.np04_data.ProtoDUNE_HD_APA_maps as p
from sklearn.cluster import DBSCAN
from tqdm import tqdm
import os
import h5py

# setup DBSCAN parameters
eps_T = 2*16 # nsec
min_samples_T = 1
eps_cm = 120 # cm
min_samples_cm = 1

data_dtype = np.dtype([('X', '<f4'), ('Y', '<f4'), ('Z', '<f4'), ('T', '<i8'), ('APA', '<i2'), \
                       ('channel', 'i4'), ('adcs', 'i4', (1024)), ('record_number', 'i4')])
header_dtype = np.dtype([('run_number', 'i4'), ('eps_T', 'i4'), ('min_samples_T', 'i4'), ('eps_cm', 'i4'), ('min_samples_cm', 'i4')])

run = 28850
filepaths_file=f'{run}_filepaths.txt'
if not os.path.exists(filepaths_file):
    raise Exception(f"Filepaths file '{filepaths_file}' does not exist")
else:
    print(f'Looking for isolated PDS hits from the file in {filepaths_file}')

header_data = np.zeros(1, dtype=header_dtype)
header_data['run_number'] = run
header_data['eps_T'] = eps_T
header_data['eps_cm'] = eps_cm
header_data['min_samples_T'] = min_samples_T
header_data['min_samples_cm'] = min_samples_cm

output_filename = f'run{run}_PDS_hitdata_test.hdf5'
with h5py.File(output_filename, 'w') as output_file:
    output_file.create_dataset('header', data=header_data)
    output_file.create_dataset('hits', data=np.zeros((0,), dtype=data_dtype), maxshape=(None,))
print(f'Output being saved to {output_filename}')

limit = 2
with open(filepaths_file, 'r') as f:
    filenum = 0
    for file in tqdm(f):
        if filenum >= limit:
            break
        #file='/eos' + file.split('/eos')[2].strip()
        #if not os.path.exists(file):
        #    raise Exception(f"File '{file} not found'")
        waveforms = rhr.WaveformSet_from_hdf5_file(file, False, 0.0, 1.0, 1)

        # look up positions for waveforms
        X, Y, Z, T, APA, R = [],[],[],[],[],[]
        channels=[]
        wvfms_list = []
        wvfms = waveforms.waveforms
        for i in range(len(wvfms)):
            channel = wvfms[i].channel + wvfms[i].endpoint*100
            try:
                channelInfo = p.position_map[channel]
                channels.append(channel)
            except:
                continue
            T.append(wvfms[i].timestamp*16*1e-3)
            X.append(channelInfo.X)
            Y.append(channelInfo.Y)
            Z.append(channelInfo.Z)
            APA.append(channelInfo.APA)
            R.append(wvfms[i].record_number)
            wvfms_list.append(wvfms[i].adcs)
        
        T = np.array(T)
        X = np.array(X)
        Y = np.array(Y)
        Z = np.array(Z)
        APA = np.array(APA)
        R = np.array(R)
        channels = np.array(channels)
        wvfms_list = np.array(wvfms_list)

        # time clustering
        db_T = DBSCAN(eps=eps_T, min_samples=min_samples_T).fit(T.reshape(-1, 1))

        # space clustering
        fake_dim = db_T.labels_*10000
        positions = np.hstack((X.reshape(-1, 1),Y.reshape(-1, 1),Z.reshape(-1, 1), fake_dim.reshape(-1, 1)))
        db_cm = DBSCAN(eps=eps_cm, min_samples=min_samples_cm).fit(positions)

        # select single hit clusters
        labels = db_cm.labels_
        unique_labels, counts = np.unique(labels, return_counts=True)
        max_cluster_size = 1
        unique_labels_keep = unique_labels[counts == max_cluster_size]
        unique_labels_keep_mask = np.isin(labels, unique_labels_keep)

        # put data in numpy array
        data = np.zeros(np.sum(unique_labels_keep_mask), dtype=data_dtype)
        data['T'] = T[unique_labels_keep_mask]
        data['X'] = X[unique_labels_keep_mask]
        data['Y'] = Y[unique_labels_keep_mask]
        data['Z'] = Z[unique_labels_keep_mask]
        data['APA'] = APA[unique_labels_keep_mask]
        data['record_number'] = R[unique_labels_keep_mask]
        data['channel'] = channels[unique_labels_keep_mask]
        data['adcs'] = wvfms_list[unique_labels_keep_mask]

        # store data in hdf5
        with h5py.File(output_filename, 'a') as f:
            f['hits'].resize((f['hits'].shape[0] + data.shape[0]), axis=0)
            f['hits'][-data.shape[0]:] = data

        filenum+=1
        
                

