import matplotlib.pyplot as plt
from mass_spec import *

# min_window_for_fitting = 0.01
ms_filename = 'data/diana/for Yaroslav/ESI MS/MUP 11.13.20/-MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# apply recalibration
data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/for Yaroslav/ESI MS/MUP 11.13.20/recalibration.npy'))(data[:,0])
# make figure
scalefactor = 1
# fig, ax = plt.subplots(figsize=(8*scalefactor,8*scalefactor))
fig, ax = plt.subplots(figsize=(8*scalefactor,12*scalefactor)) # high mz
plt.plot(data[:,0], data[:,1], label='MUP')
plt.ylim(0, np.max(data[:,1])*1.5)
plt.xlim(150, 1000)
smileses = ['[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O',
            'OP(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O',
            '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C[N+](C)(C)C',
            # '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C[N+](C)(C)CCCCCCCCCCCS[S]',
            'O=P([O-])(O)CCCCCCCCCCCSSCCCCCCCCCCCP([O-])([O-])=O.C[N+](C)(C)C.C[N+](C)(C)C'
            # 'C=CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O.C[N+](C)(C)C',
            # 'C=CCCCCCSSC=CCCCCCCCCC=P([O-])([O-])[O-].C[N+](C)(C)C.C[N+](C)(C)C'
            # 'C[N+](C)(C)C.C=CCCCCCCCCCP(O)(OOP(CCCCCCCCCC=C)([O-])=O)=O',
            # '[S-]CCCCCCCCCC',
            # '[S-]CCCCCCCCC=C',
            # '[O-]P(C)(OP(O)(C)=O)=O'
            # '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C=CCCCC[N+](C)(C)C'
            ]
offsets = [ [165, 30],
            [135, 30],
            [175, 30],
            [135, 30],
            [173, 30],
            [173, 30],
            [173, 30],
            [173, 30],
            [173, 30],
            [173, 30],
            [173, 30]
            ]
xoffsets = [0]*len(smileses)
# xoffsets[5] = 70
# xoffsets[2] = -20
exactmz = []
measuredmz = []
for i,smiles in enumerate(smileses):
    if i == 0:
        z = 2
    else:
        z = 1
    ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i], ppm=True)
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)

# ##add the 249 organics
# smiles = 'S=CCCCCCCCCCCP(O)[O-]'
# z = 1
# i = 1
# offsets_here = [130, 30]
# ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets_here[0], text_offset=offsets_here[1],
#                                            xoffset=-30, peak_width=0.05, initial_mz_guess=249.05)
# exactmz.append(ref_peak_mz)
# measuredmz.append(fit_mean)

ms_filename2 = 'data/diana/20211012 ESI-MS/TMA_MUP_4 with 70 percent MUP/+MS.txt'
# ms_filename2 = 'data/diana/20211012 ESI-MS/TMA_MUP_3 with 5 percent MUP/-MS.txt'
# ms_filename2 = 'data/diana/20211012 ESI-MS/TMA_MUP_2 with 4 percent MUP/-MS.txt'
# ms_filename2 = 'data/diana/20211012 ESI-MS/TMA_MUP_1 with 4 percent MUP/-MS.txt'


data2 = np.loadtxt(ms_filename2, delimiter='\t')
data2[:,1] = data2[:,1]/np.max(data2[:,1])
# apply recalibration
data2[:,0] = data2[:,0] + np.poly1d(np.load('data/diana/20211012 ESI-MS/recalibration.npy'))(data2[:,0])
plt.plot(data2[:,0], -data2[:,1], label='MUP/TMA')

smileses = ['[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O',
            'OP(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O',
            '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C[N+](C)(C)C',
            # '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C[N+](C)(C)CCCCCCCCCCCS[S]',
            'O=P([O-])(O)CCCCCCCCCCCSSCCCCCCCCCCCP([O-])([O-])=O.C[N+](C)(C)C.C[N+](C)(C)C']
for x in offsets:
    x[1] = -x[1]
xoffsets = [0]*len(smileses)
# xoffsets[1] = 20
# # xoffsets[2] = -20

# for i,smiles in enumerate(smileses):
#     z=1
#     ref_peak_mz, fit_mean, intensity = add_mol(data2, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
#             xoffset=xoffsets[i], flip=True)

# matching = find_maching_peaks(data, data2) # this is for low-mz part
matching = find_maching_peaks(data, data2, prominence=0.005, thresh=0.01) # high-mz part
for p in matching:
    plt.axvline(x=p, linestyle='--', color='grey', alpha=0.5, zorder=-10)
plt.ylim(np.min(-1*data2[:,1])*1.1, np.max(data[:,1])*1.5)

# plt.legend(loc='lower left')
plt.xlabel('m/z')
plt.ylabel('Intensity, normalized')
# # Use absolute value for y-ticks
# ticks =  ax.get_yticks()
# ax.set_yticklabels([abs(float(tick)) for tick in ticks])

# plt.xlim(225, 275)
# plt.savefig('data/diana/20211012 ESI-MS/TMA_MUP_4 with 70 percent MUP/MUP_TMA_70percent_comb_ESI_.png', dpi=300)

# plt.ylim(-0.07, 0.2)
# plt.xlim(495, 720)
# plt.savefig('data/diana/20211012 ESI-MS/TMA_MUP_4 with 70 percent MUP/MUP_TMA_70percent_comb_ESI_HighMZ.png', dpi=300)

plt.show()

# # # make recalibration curve
# xs = np.array(measuredmz)
# ys = np.array(exactmz)-np.array(measuredmz)
# plt.plot(xs, ys, 'o', alpha=0.5, label='Data')
# z = np.polyfit(xs, ys, 1)
# # np.save('data/diana/for Yaroslav/ESI MS/MUP 11.13.20/recalibration.npy', z)
# print(z)
# plt.plot(xs, np.poly1d(z)(xs), alpha=0.5, label='Linear fit')
# plt.ylabel('m/z error (expected minus measured), Da')
# plt.xlabel('m/z')
# plt.legend()
# plt.show()

smileses = ['[O]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O])(O)=O)(O)=O',
            'OP(CCCCCCCCCCCSSCCCCCCCCCCCP([O])(O)=O)(O)=O',
            'C=CCCCCCSSC=CCCCCCCCCC=P([O])([O])[O].C[N+](C)(C)C.C[N+](C)(C)C',
            'C[N+](C)(C)C.C=CCCCCCCCCCP(O)(OOP(CCCCCCCCCC=C)([O-])=O)=O',
            ]
zs = [2, 1, 1, 1]
for i,smiles in enumerate(smileses):
    fig, ax = plt.subplots(figsize=(3.2, 2))
    plt.title(i)
    analyze_isotopic_pattern(data, smiles, z=zs[i], dm=4, maincolor='C0', theor_color='C2')
    # fig.savefig('data/diana/for Yaroslav/ESI MS/MUP 11.13.20/MUP_TMA_comb_ESI_isot_{0}.png'.format(smiles), dpi=300)
    plt.show()