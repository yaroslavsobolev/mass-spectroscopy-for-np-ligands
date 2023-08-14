import matplotlib.pyplot as plt

from mass_spec import *

ms_filename = 'data/diana/TMA ESI/TMA/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# # apply recalibration
# data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
# make figure
scalefactor = 0.9
fig, ax = plt.subplots(figsize=(25.5*scalefactor,2.5*scalefactor))
# plt.plot(data[:,0], data[:,1])
# plt.ylim(0, np.max(data[:,1])*1.5)
# ms_filename2 = 'data/diana/for Yaroslav/ESI MS/TMA_MUP_2 17.06.21/+MS.txt'
ms_filename2 = 'data/diana/for Yaroslav/ESI MS/TMA_MUS_2 17.06.21/+MS.txt'
# ms_filename = master_folder + '-MS.txt'
data2 = np.loadtxt(ms_filename2, delimiter='\t')
data2[:,1] = data2[:,1]/np.max(data2[:,1])
plt.plot(data2[:,0], data2[:,1], label='MUP/TMA', color='C4')

matching = find_maching_peaks(data2, data, thresh=0.08)
for p in matching:
    plt.axvline(x=p, linestyle='--', color='grey', alpha=0.5, zorder=-10)
plt.ylim(0, np.max(data2[:,1])*1.1)
plt.xlim(100, 520)
# fig.savefig('data/diana/for Yaroslav/ESI MS/TMA 03.15.19/TMA_MUP_comb_.png', dpi=300)
fig.savefig('data/diana/for Yaroslav/ESI MS/TMA 03.15.19/TMA_MUS_comb_.png', dpi=300)
plt.show()