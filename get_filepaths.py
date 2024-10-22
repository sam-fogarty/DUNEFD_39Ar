# script to find reco files for a given run using metacat and rucio
# make sure to setup metacat and rucio first

import argparse
import os
#from tqdm import tqdm

def main(run):#, nfile_start, nfile_stop):

    #nfile_start = int(nfile_start)
    #nfile_stop = int(nfile_stop)    

    #file_keyword='hd-protodune-det-reco'	
    #data_tier = 'full-reconstructed'
    file_keyword='hd-protodune-det-reco'
    #data_tier='raw'
    #Format='hdf5'
    Format='root'
    data_tier='full-reconstructed'
    metacat_query = f'metacat query "files where {run} in core.runs and core.data_tier={data_tier}"'
    filenames_file = f'{file_keyword}_{run}_filenames_{data_tier}.txt'
    os.system(f'rm {filenames_file}')
    os.system(metacat_query+f' > {filenames_file}')
        
    file_list = []
    i = 0
    totalLines=0
    with open(filenames_file, 'r') as f:
        for line in f:
            totalLines+=1
    #location='FNAL_DCACHE'
    location='DUNE_US_FNAL_DISK_STAGE'
    with open(filenames_file, 'r') as f:
        for k, line in enumerate(f):
            if k % 10 == 0:
                print(f'{k}/{totalLines}')
            if line.strip().split(':')[0] == file_keyword:
                #if i >= nfile_start and i <= nfile_stop:
                os.system('rm -rf rucio_output.txt')
                #os.system(f'rucio list-file-replicas {line.strip()}')
                os.system(f'rucio list-file-replicas {line.strip()} > rucio_output.txt')
                with open('rucio_output.txt', 'r') as f2:
                    fulltxt=f2.read()
                filepath = fulltxt.split(f'{location}:')[-1].strip().split('|')[0].strip()
                file_list.append(filepath)	
                i += 1
    if len(file_list) == 0:
        print('No files found!')
    else:
        with open(f'{run}_filepaths.txt', 'w') as f:
            f.writelines(f'{line}\n' for line in file_list)
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to find reco file paths")
    parser.add_argument('run', help='Run number XXXXX to look for.')
    #parser.add_argument('nfile_start', help='Start number of filepath to get')
    #parser.add_argument('nfile_stop', help='Stop number of filepath to get')
    args = parser.parse_args()
    main(args.run)#, args.nfile_start, args.nfile_stop)
