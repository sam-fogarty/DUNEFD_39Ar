import numpy as np
import uproot
import h5py
import os
from tqdm import tqdm

charge_dir = '/pnfs/dune/scratch/users/lavaut/03684/1'
light_dir = '/exp/dune/data/users/sfogarty/DUNEFD_39Ar/28850_pickle'

run = 28850
T_lower_window = 0
T_upper_window = 2500
z_window = 60
y_window = 15
charge_dir_files = os.listdir(charge_dir)
charge_files = []
for charge_file in charge_dir_files:
   if charge_file.split('.')[-1] == 'root':
       charge_files.append(charge_file)

light_files = os.listdir(light_dir)

NearOrFar = []
Z_position = []
Y_position = []
X_position = []
cluster_peaktime = []
cluster_match_type = []
energy = []
abs_time = []
DelT_array = []
cluster_APA = []

light_file_dict = {}
for file in light_files:
    if not file.split('.')[-1] == 'pkl':
        continue
    words = file.split(f'np04hd_raw_run0{run}')[1].split('_')
    subrun, dataflow = int(words[1]), int(words[2].split('dataflow')[1])
    light_file_dict[(subrun, dataflow)] = file.split('.')[0]+'_isolatedPDSHits.hdf5'
total_matches = 0

all_cluster_t = []
all_cluster_z = []
all_cluster_y = []
all_cluster_apa = []

all_pds_t = np.zeros((0,))
all_pds_y = np.zeros((0,))
all_pds_z = np.zeros((0,))
all_pds_apa= np.zeros((0,))
nfiles = 1
for F in range(nfiles):
    charge_file = charge_files[F]
    print(charge_file)
    if not charge_file.split('.')[-1] == 'root':
        continue
    words = charge_file.split(f'np04hd_raw_run0{run}')[1].split('_')
    subrun, dataflow = int(words[1]), int(words[2].split('dataflow')[1])
    try:
        light_file = light_file_dict[(subrun, dataflow)]
    except:
        continue
    if not os.path.exists(light_dir+'/'+light_file):
        continue
    try:
        with h5py.File(light_dir+'/'+light_file, 'r') as f:
            pds_hits = np.array(f['hits'])
    except:
        print('Dataset hits does not exists, skipping file.')
        continue
    pds_t = pds_hits['T']
    pds_y = pds_hits['Y']
    pds_z = pds_hits['Z']
    pds_x = pds_hits['X']
    pds_apa = pds_hits['APA']
    all_pds_t = np.concatenate((all_pds_t, pds_t))
    all_pds_y = np.concatenate((all_pds_y, pds_y))
    all_pds_z = np.concatenate((all_pds_z, pds_z))
    all_pds_apa = np.concatenate((all_pds_apa, pds_apa))
    
    # TEST
    #pds_t = np.array([1e16*0.016])
    #pds_y = np.array([275.1595])
    #pds_z = np.array([35.33625])
    #pds_x = np.array([357.295655])
    #pds_apa = np.array([3])

    f_charge = uproot.open(charge_dir+'/'+charge_file)
    tree = f_charge["ana"]["ClusterTree"]

    NumberOfCollection = tree["NumberOfCollection"].array()
    NumberOfPlane0 = tree["NumberOfPlane0"].array()
    NumberOfPlane1 = tree["NumberOfPlane1"].array()
    NearOrFarToTheBeam = tree["NearOrFarToTheBeam"].array()
    EnergyCollection = tree["EnergyCollection"].array()
    Z = tree["Z"].array()
    Y = tree["Y"].array()
    
    # TEST
    #NumberOfCollection = np.array([[1]])
    #NumberOfPlane0 = np.array([[1]])
    #NumberOfPlane1 = np.array([[1]])
    #NearOrFarToTheBeam = np.array([[-1]])
    #EnergyCollection = np.array([[1]])
    #Z = np.array([[35.0]])
    #Y = np.array([[275]])
    
    #print(NearOrFarToTheBeam)
    #lastStart = 0
    #for k in range(len(NearOrFarToTheBeam)):
    #    print(len(NearOrFarToTheBeam[k]), '                       ', len(Z[k]))
    #    print(NearOrFarToTheBeam[k][lastStart:lastStart+10])
    #    print(len(NearOrFarToTheBeam[k]) - lastStart)
    #    lastStart = len(NearOrFarToTheBeam[k])
    #    #print(len(Z[k]))

    PeakTime = tree["PeakTime"].array()
    CRP_T0 = tree["CRP_T0"].array()
    eventID = tree["eventID"].array()
    clusterNumber = tree["NumberOfCluster"].array()
    
    # TEST
    #PeakTime = np.array([[195]])
    #CRP_T0 = np.array([1e16])
    #eventID = np.array([0])    

    withCollection = 0
    withCandP0 = 0
    withCandP1 = 0
    withAll = 0
    totalClusters = 0
    event_t0 = []
    event_ids = []
    for i in range(len(CRP_T0)):
        event_t0.append(CRP_T0[i])
        event_ids.append(eventID[i])

    # loop through events
    lastLength = 0
    for i in range(len(NumberOfCollection)):
        # loop through clusters in event
        for j in range(len(NumberOfCollection[i])):
            match_type = ''
            if NumberOfCollection[i][j]:
                withCollection += 1
            if NumberOfCollection[i][j] and NumberOfPlane0[i][j] and not NumberOfPlane1[i][j]:
                withCandP0 += 1
                match_type = 'C-P0'
            elif NumberOfCollection[i][j] and NumberOfPlane1[i][j] and not NumberOfPlane0[i][j]:
                withCandP1 += 1
                match_type = 'C-P1'
            elif NumberOfCollection[i][j] and NumberOfPlane0[i][j] and NumberOfPlane1[i][j]:
                withAll += 1
                match_type = 'C-P0-P1'
            totalClusters += 1
            cluster_time = PeakTime[i][j]*0.512 + event_t0[i]*0.016

            if Z[i][j] < 230 and NearOrFarToTheBeam[i][j] == 1:
                clusterAPA = 1
            elif Z[i][j] < 230 and NearOrFarToTheBeam[i][j] == -1:
                clusterAPA = 3
            elif Z[i][j] > 230 and NearOrFarToTheBeam[i][j] == 1:
                clusterAPA = 2
            elif Z[i][j] > 230 and NearOrFarToTheBeam[i][j] == -1:
                clusterAPA = 4
            else:
                print('Cannot find APA for cluster, skipping')
                continue
            #print('clusterAPA = ', clusterAPA)
            #print('cluster_time = ', cluster_time)
            all_cluster_t.append(cluster_time)
            all_cluster_z.append(Z[i][j])
            all_cluster_y.append(Y[i][j])
            all_cluster_apa.append(clusterAPA)
            
            #if np.sum(clusterAPA == pds_apa):
            #    print('APA match')
            #if np.sum(pds_t > cluster_time - T_lower_window):
            #    print('lower window match')
            #if np.sum(pds_t < cluster_time + T_upper_window):
            #    print('upper window match')
            #print('pds_t = ', pds_t)
            #print('pds_z = ', pds_z)
            #print('pds_y = ', pds_y)
            #print('q Z = ', Z[i][j])
            #print('q Y = ', Y[i][j]) 
            matched_mask = (clusterAPA == pds_apa) & \
                    (pds_t - T_lower_window < cluster_time) & \
                    (pds_t + T_upper_window > cluster_time) & \
                    (pds_z - z_window < Z[i][j]) & \
                    (pds_z + z_window > Z[i][j]) & \
                    (pds_y - y_window < Y[i][j]) & \
                    (pds_y + y_window > Y[i][j])
            #print('number of matched clusters = ', np.sum(matched_mask))
            nmatches = np.sum(matched_mask)
            total_matches += nmatches
            #print(np.sum(matched_mask))
            if not nmatches:
                continue
            NearOrFar.append(NearOrFarToTheBeam[i][lastLength+j])
            for index in np.where(matched_mask)[0]:
                Z_position.append(Z[i][j])
                Y_position.append(Y[i][j])
                DelT = cluster_time-pds_t[index]
                DelT_array.append(DelT)
                if pds_x[index] < 0:
                    X_position.append(pds_x[index]+DelT*0.16)
                else:
                    X_position.append(pds_x[index]-DelT*0.16)
                cluster_APA.append(clusterAPA) 
            cluster_peaktime.append(PeakTime[i][j])
            cluster_match_type.append(match_type)
            energy.append(EnergyCollection[i][j])
            #cluster_APA.append(clusterAPA)
            abs_time.append(cluster_time)
        lastLength = len(NearOrFarToTheBeam[i])

np.savez('all_hit_data_smallerProxCut_TwoOrThreePlanes.npz', cluster_t=all_cluster_t, cluster_y=all_cluster_y, cluster_z=all_cluster_z, cluster_apa=all_cluster_apa, pds_t=all_pds_t, pds_z=all_pds_z, pds_y=all_pds_y, pds_apa=all_pds_apa)
print(f'Del T = {DelT_array}')
np.savez('matched_hit_data_smallerProxCut_TwoOrThreePlanes.npz', APA=cluster_APA, DelT=DelT_array, X=X_position, Y=Y_position, Z=Z_position)
#print(f'Total clusters = {totalClusters}')
#print(f'Total with any collection = {withCollection}')
#print(f'Total with collection and P0 = {withCandP0}')
#print(f'Total with collection and P1 = {withCandP1}')
#print(f'Total with all three planes = {withAll}')
#print(f'Sum = {withCandP0+withCandP1+withAll}')
print(f'total matched clusters = {total_matches}')
