# import openbabel
import matplotlib.pyplot as plt
import pybel
from mass_spec import *


ms_filename = 'data/diana/TMA ESI/TMA/+MS.txt'
data = np.loadtxt(ms_filename, delimiter='\t')
data[:,1] = data[:,1]/np.max(data[:,1])
# apply recalibration
data[:,0] = data[:,0] + np.poly1d(np.load('data/diana/TMA ESI/TMA/recalibration.npy'))(data[:,0])
zs = [1,2,1]
smileses = ['C[N+](C)(C)CCCCCCCCCCCS[S-]',
            'C[N+](C)(C)CCCCCCCCCCCSSCCCCCCCCCCC[N+](C)(C)C',
            'C[N+](C)(C)CCCCCCCCCC=C']
for i,smiles in enumerate(smileses):
    fig, ax = plt.subplots(figsize=(4, 2))
    analyze_isotopic_pattern(data, smiles, z=zs[i], maincolor='C0', theor_color='C2')
    fig.savefig('data/diana/TMA ESI/TMA/{0}.png'.format(smiles), dpi=300)
    plt.show()


# mymol.draw(show=False, filename='data/molecule_images/{0}.svg'.format(smiles))