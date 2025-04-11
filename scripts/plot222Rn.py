import uproot
import numpy as np
import matplotlib.pyplot as plt

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,12))

f = uproot.open('Rn222_PDHD_test_028850.root')
tree = f["results"]
DelT = tree["DelT"].array()
E_1 = tree["E_1"].array()
E_2 = tree["E_2"].array()
print('total matches = ', len(DelT))
axes[0].hist(-1*DelT, bins=25, range=(0,25))
factor = 5E-3
R = 0.67
Wions = 23.6*10E-6
calib = Wions/(factor*R)

axes[1].hist(E_1*calib, range=(0, 5), bins=50)
axes[1].hist(E_2*calib, range=(0,5), bins=50)
plt.show()
