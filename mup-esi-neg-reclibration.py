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
fig, ax = plt.subplots(figsize=(8*scalefactor,8*scalefactor))
plt.plot(data[:,0], data[:,1], label='MUP')
plt.ylim(0, np.max(data[:,1])*1.5)
plt.xlim(150, 1000)
smileses = ['[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O',
            'OP(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O',
            '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C[N+](C)(C)C'
            # '[O-]P(CCCCCCCCCCCSSCCCCCCCCCCCP([O-])(O)=O)(O)=O.C[N+](C)(C)CCCCCCCCCCCS[S]',
            # 'O=P([O-])(O)CCCCCCCCCCCSSCCCCCCCCCCCP([O-])([O-])=O.C[N+](C)(C)C.C[N+](C)(C)C'
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
            [185, 30],
            [185, 30],
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
            xoffset=xoffsets[i])
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)

ms_filename2 = 'data/diana/20211012 ESI-MS/TMA_MUP_4 with 70 percent MUP/-MS.txt'
# ms_filename = master_folder + '-MS.txt'
data = np.loadtxt(ms_filename2, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
plt.plot(data[:,0], -data[:,1], label='MUP/TMA')
# apply recalibration
# data2[:,0] = data2[:,0] + np.poly1d(np.load('data/diana/for Yaroslav/ESI MS/MUP 11.13.20/recalibration.npy'))(data2[:,0])

measuredmz2 = []
for i,smiles in enumerate(smileses):
    if i == 0:
        z = 2
    else:
        z = 1
    ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i], flip=True)
    exactmz.append(ref_peak_mz)
    measuredmz2.append(fit_mean)
plt.show()

# # make recalibration curve
xs = np.array(measuredmz2)
ys = np.array(measuredmz) - np.array(measuredmz2)
plt.plot(xs, ys, 'o', alpha=0.5, label='Data')
z = np.polyfit(xs, ys, 1)
np.save('data/diana/20211012 ESI-MS/recalibration.npy', z)
print(z)
plt.plot(xs, np.poly1d(z)(xs), alpha=0.5, label='Linear fit')
plt.ylabel('m/z error (expected minus measured), Da')
plt.xlabel('m/z')
plt.legend()
plt.show()