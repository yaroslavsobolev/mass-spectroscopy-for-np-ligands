import matplotlib.pyplot as plt
from mass_spec import *

ms_filename = 'data/diana/TMA ESI/TMA/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
print('Sampling rate: {0}'.format(data[10,0] - data[9,0]))
data[:,1] = data[:,1]/np.max(data[:,1])
# apply recalibration
data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
scalefactor = 0.9
fig, ax = plt.subplots(figsize=(25.5*scalefactor,11*scalefactor)) # was 20 by 11
plt.plot(data[:,0], data[:,1])
plt.ylim(0, np.max(data[:,1])*1.95)
plt.xlim(100, 520)
smileses = ['[H+].C[N](C)(C)CCCCCCCCCCCSSC=CCCCCCCCCC[N](C)(C)C',
            'CN(C)CCCCCCCCCCCSSCCCCCCCCCCC[N+](C)(C)C',
            'C[N+](C)(C)CCCCCCCCCCCSSCCCCCCCCCC=C',
            '[S]SCCCCCCCCCCC[N+](C)(C)C',
            'C[N+](C)(C)CCCCCCCCCCCSSCCCCCCCCCCC[N+](C)(C)C',
            'C=CCCCCCCCCC[N+](C)(C)C',
            'C=CCCCCCCCC[N+](C)(C)C',
            'C=CCCCCCCC[N+](C)(C)C',
            'C=CCCCCCC[N+](C)(C)C',
            'C=CCCCCC[N+](C)(C)C',
            'C=CCCCC[N+](C)(C)C',
            'C=CCCC[N+](C)(C)C',
            'C=CCC[N+](C)(C)C']
offsets = [ [165, 30],
            [165, 30],
            [165, 30],
            [175, 30],
            [175, 30],
            [148, 22],
            [175+30+25, 80+30],
            [140, 30],
            [130, 30],
            [125, 30],
            [102+60+15, 30+60],
            [95+12, 30],
            [90+40+7, 30+40]
]
xoffsets = [0]*len(smileses)
xoffsets[0] = 20
xoffsets[1] = -20
exactmz = []
measuredmz = []
for i,smiles in enumerate(smileses):
    if i == 4:
        z = 2
    else:
        z = 1
    ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i], ppm=True)
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)
plt.xlabel('m/z')
plt.ylabel('Intensity, normalized')
# plt.savefig('TMA_ESI_comb_.png', dpi=300)
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
plt.xlabel('m/z, Da')
plt.legend()
plt.savefig('TMA_ESI_error.png', dpi=300)
plt.show()