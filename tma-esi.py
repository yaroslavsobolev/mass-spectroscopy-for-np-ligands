import matplotlib.pyplot as plt
from mass_spec import *

ms_filename = 'data/diana/TMA ESI/TMA/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# # apply recalibration
# data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
# make figure
fig, ax = plt.subplots(figsize=(20,9))
plt.plot(data[:,0], data[:,1])
plt.ylim(0, np.max(data[:,1])*1.85)
plt.xlim(100, 520)
smileses = ['[H+].C[N](C)(C)CCCCCCCCCCCSSC=CCCCCCCCCC[N](C)(C)C',
            'CN(C)CCCCCCCCCCCSSCCCCCCCCCCC[N+](C)(C)C',
            'C[N+](C)(C)CCCCCCCCCCCSSCCCCCCCCCC=C',
            'C[N+](C)(C)CCCCCCCCCCCS[S-]',
            'C[N+](C)(C)CCCCCCCCCCCSSCCCCCCCCCCC[N+](C)(C)C',
            'C[N+](C)(C)CCCCCCCCCC=C',
            'C[N+](C)(C)CCCCCCCCC=C',
            'C[N+](C)(C)CCCCCCCC=C',
            'C[N+](C)(C)CCCCCCC=C',
            'C[N+](C)(C)CCCCCC=C',
            'C[N+](C)(C)CCCCC=C',
            'C[N+](C)(C)CCCC=C',
            'C[N+](C)(C)CCC=C']
offsets = [ [140, 30],
            [140, 30],
            [140, 30],
            [150, 30],
            [150, 30],
            [135, 30],
            [175+30, 80+30],
            [120, 30],
            [115, 30],
            [110, 30],
            [102+60, 30+60],
            [95, 30],
            [90+40, 30+40]
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
            xoffset=xoffsets[i])
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)
plt.xlabel('m/z')
plt.ylabel('Intensity, normalized')
plt.savefig('temp.png', dpi=300)
plt.show()

# make recalibration curve
xs = np.array(measuredmz)
ys = np.array(exactmz)-np.array(measuredmz)
plt.plot(xs, ys, 'o', alpha=0.5, label='Data')
z = np.polyfit(xs, ys, 1)
np.save('data/diana/TMA ESI/TMA/recalibration.npy', z)
print(z)
plt.plot(xs, np.poly1d(z)(xs), alpha=0.5, label='Linear fit')
plt.ylabel('m/z error (expected minus measured), Da')
plt.xlabel('m/z')
plt.legend()
plt.show()