import h5py
import uproot
import numpy as np
import os
import sys

def hdf5_to_root(input_filepath):
    if not os.path.isfile(input_filepath) or not input_filepath.endswith('.hdf5'):
        print("Please provide a valid HDF5 file.")
        return
    
    output_filepath = input_filepath.replace('.hdf5', '.root')
    
    with h5py.File(input_filepath, 'r') as hdf5_file:
        dataset = np.array(hdf5_file['hits'])
        
        z_data = dataset['Z'][:]
        y_data = dataset['Y'][:]
        t_data = dataset['T'][:]
        apa_data = dataset['APA'][:]
        data = {
            "Z": z_data,
            "Y": y_data,
            "T": t_data,
            "APA": apa_data
        }
        
        with uproot.recreate(output_filepath) as root_file:
            root_file["hits"] = {
                "Z": z_data,
                "Y": y_data,
                "T": t_data,
                "APA": apa_data
            }
            #root_file["hits"].extend(data)
        
    print(f"Data successfully written to {output_filepath}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python hdf5_to_root.py <input_filepath>")
    else:
        hdf5_to_root(sys.argv[1])
