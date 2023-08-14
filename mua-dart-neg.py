import matplotlib.pyplot as plt
from mass_spec import *

# min_window_for_fitting = 0.01

ms_filename = 'data/diana/for Yaroslav/DART MS/MUA 01.14.21/-MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# # apply recalibration
# data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
# make figure
scalefactor = 1
fig, ax = plt.subplots(figsize=(8*scalefactor,10*scalefactor))
plt.plot(data[:,0], data[:,1])
plt.ylim(0, np.max(data[:,1])*3.95)
plt.xlim(170, 280)
smileses = [#'C#CCCCCCCCCC(O)[O-]'
            '[O-]C(CCCCCCCCC=C)=O',
            # 'OC#CCCCCCCCCCS[S-]',

            'S=CCCCCCCCCCC([O-])=O',
            '[O-]C#CCCCCCCCCCSS',
            'CCCCCCCCCCCC([O-])=O.O=O',
            # '[S-]SCCCCCCCCCC=C=O',
            # '[O-]C(O)CCCCCCCCC=CSSC',
            # '[O-]C(CCCCCCCCC=CSSC)=O',
            '[O-]C(CCCCCCCCCCSSC)=O',
            # 'OC#CCCCCCCCCCS[S-].O=O'
            ]
offsets = [ [150, 30],
            [160, 30],
            [865, 60],
            [855, 60],
[173, 30]
]
# xoffsets = [0]*len(smileses)
# xoffsets[2] = -30
# xoffsets[3] = 40
xoffsets = [-10,
            -40,
            -40,
            40,
            0]
exactmz = []
measuredmz = []
for i,smiles in enumerate(smileses):
    z = 1
    ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i])
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)



ms_filename2 = 'data/diana/for Yaroslav/DART MS/TMA_MUA_2 17.06.21/-MS.txt'
# ms_filename = master_folder + '-MS.txt'
data2 = np.loadtxt(ms_filename2, delimiter='\t')
data2[:,1] = data2[:,1]/np.max(data2[:,1])
plt.plot(data2[:,0], -data2[:,1])

matching = find_maching_peaks(data, data2)
for p in matching:
    plt.vlines(x=p, linestyle='--', color='grey', alpha=0.5, zorder=-10, ymin=-1.1, ymax=0)#data2[np.where(data2[:,0]==p)[0], 1]*2)
plt.ylim(np.min(-1*data2[:,1])*1.1, np.max(data[:,1])*1.45) #3.23

plt.xlabel('m/z')
plt.ylabel('Intensity, normalized')
ticks =  ax.get_yticks()
ax.set_yticklabels([abs(float(tick)) for tick in ticks])
plt.savefig('data/diana/for Yaroslav/DART MS/MUA 01.14.21/MUA_DART_comb_.png', dpi=300)
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

smileses = ['O=O.CCCCCCCCCCCC([O])=O',
            '[O]C#CCCCCCCCCCSS']
zs = [1]*len(smileses)
for i,smiles in enumerate(smileses):
    fig, ax = plt.subplots(figsize=(2.7, 2))
    analyze_isotopic_pattern(data, smiles, z=zs[i], dm=3, maincolor='C0', theor_color='C2')
    fig.savefig('data/diana/for Yaroslav/DART MS/MUA 01.14.21/MUA_DART_comb_isot_{0}.png'.format(smiles), dpi=300)
    plt.show()