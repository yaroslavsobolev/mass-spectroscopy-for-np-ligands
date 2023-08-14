import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit
from mass_spec import *

ref_peak_mz = 121.050873
peak_width = 10e-02
ms_filename = 'data/diana/Standard/ST1/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')

ref_list = np.loadtxt('data/reference_spectra/APCI_calibrant.txt', delimiter='\t', skiprows=1)
pos_refs = ref_list[:,1]
pos_refs = pos_refs[np.logical_and(pos_refs < np.max(data[:,0]), pos_refs > np.min(data[:,0]))]
diffs = []
for ref_peak_mz in pos_refs:
    fit_mean, _, _, _ = fit_peak_twice(data, ref_peak_mz, peak_width)
    diffs.append(fit_mean-ref_peak_mz)

#2-acetylanthracene
target = 'C16H12O'
ref_peak_mz = 221.09665
ms_filename = 'data/diana/Standard/ST2/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
fit_mean, _, _, _ = fit_peak_twice(data, ref_peak_mz, peak_width)
print(fit_mean - ref_peak_mz)

diffs = np.array(diffs)
pos_refs = np.insert(pos_refs, 1, ref_peak_mz)
diffs = np.insert(diffs, 1, fit_mean-ref_peak_mz)
plt.plot(pos_refs, diffs, 'o-')
plt.show()
# plt.plot(data[:,0], data[:,1])
# plt.show()

