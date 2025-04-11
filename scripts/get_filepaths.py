# script to find reco files for a given run using metacat and rucio
# make sure to setup metacat and rucio first

import argparse
import os
import time

def main(run):
    file_keyword='hd-protodune-det-reco'
    Format='root'
    data_tier='full-reconstructed'
    metacat_query = f'metacat query "files where {run} in core.runs and core.data_tier={data_tier} and core.run_type=hd-protodune"'
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
            time.sleep(0.2)
            if k % 10 == 0:
                print(f'{k}/{totalLines}')

            if line.strip().split(':')[0] == file_keyword:
                os.system('rm -rf rucio_output.txt')
                os.system(f'rucio list-file-replicas {line.strip()} > rucio_output.txt')
                with open('rucio_output.txt', 'r') as f2:
                    fulltxt=f2.read()
                filepath = fulltxt.split(f'{location}:')[-1].strip().split('|')[0].strip()
                file_list.append(filepath)	
                i += 1
            if k % 100 and len(file_list) != 0:
                with open(f'{run}_filepaths.txt', 'a') as f:
                    f.writelines(f'{line}\n' for line in file_list)
                file_list = []
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Script to find reco file paths")
    parser.add_argument('run', help='Run number XXXXX to look for.')
    args = parser.parse_args()
    main(args.run)
