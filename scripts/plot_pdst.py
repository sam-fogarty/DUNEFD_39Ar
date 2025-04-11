import matplotlib.pyplot as plt
import numpy as np
import uproot
import h5py

f = uproot.open('all_timestamps.root')
t_pds = f["light"]['t'].array()
eventID_pds = f["light"]['eventID'].array()
t_charge = f["light"]['t'].array()
eventID_charge = f["charge"]['eventID'].array()
f.close()

print(t_pds[0:10])
print(t_charge[0:10])
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(10,6))
eventID = 10
print(eventID_pds)
print(len(t_pds[eventID_pds == eventID+1]))
print(len(t_charge[eventID_charge == eventID]))
axes[0].hist(t_pds[eventID_pds == eventID+1], bins=20, label='light')
axes[0].hist(t_charge[eventID_charge == eventID], bins=20, label='charge')
print(min(t_pds[eventID_pds == eventID+1]) - min(t_charge[eventID_charge == eventID]))
#plt.hist(t_charge[eventID_charge == eventID], bins=20, label='charge')
#axes[0].yscale('log')
#axes[1].yscale('log')
axes[0].legend()
axes[1].legend()
plt.show()
