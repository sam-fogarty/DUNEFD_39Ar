import uproot
import numpy as np
import matplotlib.pyplot as plt

filepath = '/pnfs/dune/scratch/users/lavaut/03684/1/np04hd_raw_run028850_0074_dataflow3_datawriter_0_20240822T131102_reco_stage1_reco_stage2_20240822T180717_keepup_singleHit_2024-10-14T_152233Z.root'
f_charge = uproot.open(filepath)
tree = f_charge["ana"]["ClusterTree"]
t0 = tree["CRP_T0"].array()
print(t0)
plt.hist(t0, bins=50)
plt.show()

