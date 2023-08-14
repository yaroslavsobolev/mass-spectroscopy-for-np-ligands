import matplotlib.pyplot as plt
from mass_spec import *

# min_window_for_fitting = 0.01

ms_filename = 'data/diana/for Yaroslav/ESI MS/MUS 11.06.20/-MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# # apply recalibration
# data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
# make figure
scalefactor = 1
fig, ax = plt.subplots(figsize=(8*scalefactor,8*scalefactor))
plt.plot(data[:,0], data[:,1], label='MUP')
plt.ylim(0, np.max(data[:,1])*1.5)
plt.xlim(242, 305)
smileses = ['[O-]S(CCCCCCCCCCCSSCCCCCCCCCCCS(=O)([O-])=O)(=O)=O'
            ]
offsets = [ [370, 30],
            [370, 30],
            [165, 30],
            [155, 30],
[173, 30]
]
xoffsets = [0, 0, 0]
# xoffsets[1] = -30
# xoffsets[3] = 40
exactmz = []
measuredmz = []
for i,smiles in enumerate(smileses):
    z = 2
    ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i])
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)

# smiles = 'S=CCCCCCCCCCCP(O)[O-]'
# z = 1
# i = 1
# ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
#                                            xoffset=xoffsets[i], peak_width=0.05, initial_mz_guess=249.05)
# exactmz.append(ref_peak_mz)
# measuredmz.append(fit_mean)

ms_filename2 = 'data/diana/for Yaroslav/ESI MS/TMA_MUS_2 17.06.21/-MS.txt'
# ms_filename = master_folder + '-MS.txt'
data2 = np.loadtxt(ms_filename2, delimiter='\t')
data2[:,1] = data2[:,1]/np.max(data2[:,1])
plt.plot(data2[:,0], -data2[:,1], label='MUP/TMA')

smileses = ['[O-]S(CCCCCCCCCCCS)(=O)=O',
            '[O-]S(CCCCCCCCCC=CSS)(=O)=O'
            ]
offsets = [ [-365, -50],
            [-335, -60],
            [-185, -30]
]
xoffsets = [0]*len(smileses)
# xoffsets[1] = 20
# xoffsets[2] = -20
for i,smiles in enumerate(smileses):
    z=1
    ref_peak_mz, fit_mean, intensity = add_mol(data2, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i], flip=True)
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)

# matching = find_maching_peaks(data, data2)
# for p in matching:
#     plt.axvline(x=p, linestyle='--', color='grey', alpha=0.5, zorder=-10)
plt.ylim(np.min(-1*data2[:,1])*1.5, np.max(data[:,1])*1.5)

# plt.legend(loc='lower left')
plt.xlabel('m/z')
plt.ylabel('Intensity, normalized')
# Use absolute value for y-ticks
ticks =  ax.get_yticks()
ax.set_yticklabels([abs(float(tick)) for tick in ticks])
plt.savefig('data/diana/for Yaroslav/ESI MS/TMA_MUS_2 17.06.21/MUS_TMA_comb_ESI_.png', dpi=300)
plt.show()

# # make recalibration curve
xs = np.array(measuredmz)
ys = np.array(exactmz)-np.array(measuredmz)
plt.plot(xs, ys, 'o', alpha=0.5, label='Data')
z = np.polyfit(xs, ys, 1)
# np.save('data/diana/TMA ESI/TMA/recalibration.npy', z)
# print(z)
plt.plot(xs, np.poly1d(z)(xs), alpha=0.5, label='Linear fit')
plt.ylabel('m/z error (expected minus measured), Da')
plt.xlabel('m/z')
plt.legend()
plt.show()

smileses = ['[O]S(CCCCCCCCCCCSSCCCCCCCCCCCS(=O)([O])=O)(=O)=O',
            '[O]S(CCCCCCCCCCCS)(=O)=O']
zs = [2,1]
for i,smiles in enumerate(smileses):
    fig, ax = plt.subplots(figsize=(3.2, 2))
    analyze_isotopic_pattern(data, smiles, z=zs[i], dm=2.5, maincolor='C0', theor_color='C2')
    fig.savefig('data/diana/for Yaroslav/ESI MS/TMA_MUS_2 17.06.21/MUS_TMA_comb_ESI_isot_{0}.png'.format(smiles), dpi=300)
    plt.show()