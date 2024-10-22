import numpy as np
from waffles.input import raw_hdf5_reader as rhr
from tqdm import tqdm
import os
import h5py
import pickle
import time

run = 28850
filepaths_file=f'{run}_filepaths.txt'
if not os.path.exists(filepaths_file):
    raise Exception(f"Filepaths file '{filepaths_file}' does not exist")
os.system(f'mkdir -p {run}_pickle')

limit = 1000
with open(filepaths_file, 'r') as f:
    filenum = 0
    for file in tqdm(f):
        if filenum == limit:
            break
        file = file.strip()
        filename = os.path.basename(file)
        pickle_filename = f'{run}_pickle/'+filename.split('.')[0]+'.pkl'
        hdf5_filename = f'{run}_pickle/'+filename
        if os.path.exists(pickle_filename):
            print(f'pickled file already exists, skipping: {pickle_filename}')
            filenum+=1
            continue
        else:
            time.sleep(2)
        print(f'xrdcp {file} {run}_pickle/')
        os.system(f'xrdcp {file} {run}_pickle/')
        #file='/eos' + file.split('/eos')[2].strip()
        #if not os.path.exists(file):
        #    raise Exception(f"File '{file} not found'")
        print('filename = ', filename)
        waveforms = rhr.WaveformSet_from_hdf5_file(hdf5_filename, False, 0.0, 1.0, 1)
        #print(filename.split('.')[0]+'.pkl')
        with open(pickle_filename, "wb") as f2:
            pickle.dump(waveforms, f2)
        os.system(f'rm {hdf5_filename}')
        filenum+=1
