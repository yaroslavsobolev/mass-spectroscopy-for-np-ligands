import matplotlib.pyplot as plt
from mass_spec import *

min_window_for_fitting = 0.01

master_folder = 'data/diana/for Yaroslav/ESI MS/TMA_MUP_2 17.06.21/'
ms_filename = master_folder + '-MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# # apply recalibration
# data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
# make figure
scalefactor = 1
fig, ax = plt.subplots(figsize=(8*scalefactor,5*scalefactor))
plt.plot(data[:,0], data[:,1])
plt.ylim(0, np.max(data[:,1])*1.5)
plt.xlim(150, 280)
smileses = ['S=CCCCCCCCCCCP(O)[O-]'
            ]
offsets = [ [370, 30],
            [160, 30],
            [165, 30],
            [155, 30],
            [173, 30]
]
xoffsets = [0]*len(smileses)
# xoffsets[2] = -30
# xoffsets[3] = 40
exactmz = []
measuredmz = []
for i,smiles in enumerate(smileses):
    z=1
    ref_peak_mz, fit_mean, intensity = add_mol(data, ax, smiles, z=z, img_offset=offsets[i][0], text_offset=offsets[i][1],
            xoffset=xoffsets[i])
    exactmz.append(ref_peak_mz)
    measuredmz.append(fit_mean)
plt.xlabel('m/z')
plt.ylabel('Intensity, normalized')
plt.savefig(master_folder + 'MUP_TMA_ESI_.png', dpi=300)
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

smileses = ['S=CCCCCCCCCCCP(O)[O]']
zs = [1]
for i,smiles in enumerate(smileses):
    fig, ax = plt.subplots(figsize=(3.2, 2))
    analyze_isotopic_pattern(data, smiles, z=zs[i], dm=2.5, maincolor='C1', theor_color='C2')
    fig.savefig(master_folder + 'MUP_TMA_ESI_isot_{0}.png'.format(smiles), dpi=300)
    plt.show()