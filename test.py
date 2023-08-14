import matplotlib.pyplot as plt

from mass_spec import *

ref_peak_mz = 246.22554
peak_width = 0.5
ms_filename = 'data/diana/TMA DART/TMA#1/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
fit_mean, _, _, a = fit_peak_twice(data, ref_peak_mz, peak_width)

fig, ax = plt.subplots()
plt.plot(data[:,0], data[:,1])
# ymax = 3e5
# deltay = 1e4
# plt.plot([fit_mean, fit_mean], [a+deltay, ymax], color='black', alpha=0.3)

an2 = ax.annotate("{0:.3f}".format(fit_mean), xy=(fit_mean, a),
                  xytext=(0, 30), textcoords="offset points",
                  arrowprops=dict(arrowstyle="->"),
                  ha='center')

plt.show()
print(fit_mean - ref_peak_mz)