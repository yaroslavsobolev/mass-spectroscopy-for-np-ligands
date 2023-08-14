import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox
import matplotlib.image as mpimg
from scipy.stats import norm
from scipy.optimize import curve_fit
import pybel
from molmass import Formula
import pickle
from scipy.signal import find_peaks

## VALUE FOR TMA
min_window_for_fitting = 0.05
number_of_sigmas_for_second_window = 1.5

# # VALUE FOR MUA
# min_window_for_fitting = 0.006
# number_of_sigmas_for_second_window = 1.2



def asymm_gauss(x, mean, s1, s2):
    '''
    Asymmetric gaussian function. It has two sigmas instead of one.

    Parameters
    ----------
    x : ndarray
        input agrument
    mean : float
        not actually the mean, but the value of x where the function reaches maximum
    s1 : float
        sigma on the left side
    s2 : float
        sigma on the right side

    Returns
    -------
    ndarray
        The function evaluated at input x
    '''
    return np.heaviside(x-mean, 0.5)*np.exp(-(x-mean)**2/(2*s2**2)) + np.heaviside(mean-x, 0.5)*np.exp(-(x-mean)**2/(2*s1**2))
#
#
def fit_peak(data, ref_peak_mz, peak_width, do_plot=False):
    xdata = data[:,0]
    ydata = data[np.logical_and(xdata > ref_peak_mz-peak_width, xdata < ref_peak_mz+peak_width),1]
    xdata = xdata[np.logical_and(xdata > ref_peak_mz-peak_width, xdata < ref_peak_mz+peak_width)]
    def func(x, mean, s1, s2, a):
        return a*asymm_gauss(x, mean, s1, s2)
    popt, pcov = curve_fit(func, xdata, ydata, p0=(ref_peak_mz, peak_width, peak_width, np.max(ydata)),
                           bounds=[(0, 0, 0, 0), (np.inf, np.inf, np.inf, np.inf)])
    _, s1, s2, a = popt
    xs = np.linspace(ref_peak_mz-peak_width, ref_peak_mz+peak_width, 100)
    print(popt)
    if do_plot:
        plt.plot(xdata, ydata, 'o')
        plt.plot(xs, func(xs, *popt), 'r-')
        plt.title('{0}'.format(peak_width))
        plt.show()
    return popt

def fit_peak_twice(data, ref_peak_mz, peak_width, do_plot=False):
    mean, s1, s2, a = fit_peak(data, ref_peak_mz, peak_width, do_plot=do_plot)
    popt = fit_peak(data, mean, peak_width=max(min_window_for_fitting, number_of_sigmas_for_second_window*max([s1, s2])), do_plot=do_plot)
    return popt

def monomass(smiles):
    return pybel.readstring('smi', smiles).exactmass

def fit_mol(data, smiles, z=1, peak_width = 0.4, initial_mz_guess = False):
    ref_peak_mz = monomass(smiles)/z
    if initial_mz_guess == False:
        initial_mz_guess = ref_peak_mz
    # peak_width = 0.1
    # if z==2:
    #     do_plot = True
    #     fig2 = plt.figure(2)
    #     fig3 = plt.figure(3)
    # else:
    #     do_plot = False
    fit_mean, _, _, a = fit_peak_twice(data, initial_mz_guess, peak_width)
    print('z={0}, dm={1}'.format(z, fit_mean-ref_peak_mz))
    return ref_peak_mz, fit_mean, a

# ## Version where I still used Marvin JS for images
# def annotate_mol(smiles, ref_peak_mz, fit_mean, intensity, ax, img_offset=130, text_offset=30,
#                  xoffset=0, z=1, ppm=False):
#     if not ppm:
#         if z == 1:
#             an2 = ax.annotate("{0:.4f}\n{1:.2f} mDa".format(fit_mean, 1e3 * (fit_mean - ref_peak_mz)), xy=(fit_mean, intensity),
#                               xytext=(xoffset, text_offset), textcoords="offset points",
#                               arrowprops=dict(arrowstyle="->", alpha=0.5),
#                               ha='center')
#         else:
#             an2 = ax.annotate("{0:.4f}\nz={2}\n{1:.2f} mDa".format(fit_mean, 1e3 * (fit_mean - ref_peak_mz), z), xy=(fit_mean, intensity),
#                               xytext=(xoffset, text_offset), textcoords="offset points",
#                               arrowprops=dict(arrowstyle="->", alpha=0.5),
#                               ha='center')
#     else:
#         if z == 1:
#             an2 = ax.annotate("{0:.4f}\n{1:.1f} ppm".format(fit_mean, 1e6 * (fit_mean - ref_peak_mz)/ref_peak_mz),
#                               xy=(fit_mean, intensity),
#                               xytext=(xoffset, text_offset), textcoords="offset points",
#                               arrowprops=dict(arrowstyle="->", alpha=0.5),
#                               ha='center')
#         else:
#             an2 = ax.annotate("{0:.4f}\nz={2}\n{1:.1f} ppm".format(fit_mean, 1e6 * (fit_mean - ref_peak_mz)/ref_peak_mz, z),
#                               xy=(fit_mean, intensity),
#                               xytext=(xoffset, text_offset), textcoords="offset points",
#                               arrowprops=dict(arrowstyle="->", alpha=0.5),
#                               ha='center')
#     arr_lena = mpimg.imread('data/molecule_images/{0}.png'.format(smiles))
#     imagebox = OffsetImage(arr_lena, zoom=0.4)
#     ab = AnnotationBbox(imagebox, xy=(fit_mean, intensity),
#                         xybox=(xoffset, img_offset),
#                         xycoords='data',
#                         boxcoords="offset points",
#                         pad=0.5, frameon=False)
#     ab.zorder = -100
#     ax.add_artist(ab)

def annotate_mol(smiles, ref_peak_mz, fit_mean, intensity, ax, img_offset=130, text_offset=30,
                 xoffset=0, z=1, ppm=False):
    if not ppm:
        if z == 1:
            an2 = ax.annotate("{0:.4f}\n{1:.0f} mDa".format(fit_mean, 1e3 * (fit_mean - ref_peak_mz)), xy=(fit_mean, intensity),
                              xytext=(xoffset, text_offset), textcoords="offset points",
                              arrowprops=dict(arrowstyle="->", alpha=0.5),
                              ha='center')
        else:
            an2 = ax.annotate("{0:.4f}\nz={2}\n{1:.0f} mDa".format(fit_mean, 1e3 * (fit_mean - ref_peak_mz), z), xy=(fit_mean, intensity),
                              xytext=(xoffset, text_offset), textcoords="offset points",
                              arrowprops=dict(arrowstyle="->", alpha=0.5),
                              ha='center')
    else:
        if z == 1:
            an2 = ax.annotate("{0:.4f}\n{1:.1f} ppm".format(fit_mean, 1e6 * (fit_mean - ref_peak_mz)/ref_peak_mz),
                              xy=(fit_mean, intensity),
                              xytext=(xoffset, text_offset), textcoords="offset points",
                              arrowprops=dict(arrowstyle="->", alpha=0.5),
                              ha='center')
        else:
            an2 = ax.annotate("{0:.4f}\nz={2}\n{1:.1f} ppm".format(fit_mean, 1e6 * (fit_mean - ref_peak_mz)/ref_peak_mz, z),
                              xy=(fit_mean, intensity),
                              xytext=(xoffset, text_offset), textcoords="offset points",
                              arrowprops=dict(arrowstyle="->", alpha=0.5),
                              ha='center')
    arr_lena = mpimg.imread('data/molecule_images/{0}.png'.format(smiles))
    imagebox = OffsetImage(arr_lena, zoom=0.15)
    ab = AnnotationBbox(imagebox, xy=(fit_mean, intensity),
                        xybox=(xoffset, img_offset),
                        xycoords='data',
                        boxcoords="offset points",
                        pad=0.5, frameon=False)
    ab.zorder = -100
    ax.add_artist(ab)

def add_mol(data, ax, smiles, z=1, img_offset=130, text_offset=130, xoffset=0, peak_width = 0.4, initial_mz_guess = False,
            flip=False, ppm=False):
    ref_peak_mz, fit_mean, intensity = fit_mol(data, smiles, z=z, peak_width=peak_width, initial_mz_guess=initial_mz_guess)
    if not flip:
        annotate_mol(smiles, ref_peak_mz, fit_mean, intensity, ax, z=z, img_offset=img_offset,
                     text_offset=text_offset, xoffset=xoffset, ppm=ppm)
    else:
        annotate_mol(smiles, ref_peak_mz, fit_mean, -intensity, ax, z=z, img_offset=img_offset,
                     text_offset=text_offset, xoffset=xoffset, ppm=ppm)
    return ref_peak_mz, fit_mean, intensity


def scaled_asymm_gauss(x, mean, s1, s2, a):
    return a*asymm_gauss(x, mean, s1, s2)

def sim_spec(x, ref, s1, s2, a):
    """Function that simulates the mass spectrum

    It takes in the txt file with m/z values and abundanced generated by
    sisweb.com isotope pattern simulator. Then it convolutes
    the Dirac delta functions at these m/z with a peak function.
    In this case the peak function is an asymmetric gaussian.

    Parameters
    ----------
    x : 1d numpy array
        The m/z at which the simulated mass spectrum is to be evaluated.
    ref : 2d numpy array
        list of central m/z and abundances for peaks

    Returns
    -------
    spectrum
        Array of count values corresponding to input x.
    """
    res = 0
    for i in range(ref.shape[0]):
        res += scaled_asymm_gauss(x, ref[i,0], s1, s2, ref[i,1]/100)
        # res += voigt(x, 1, ref[i,0], lw, shape)*ref[i,1]/100
        # res += norm.pdf(x, ref[i,0], lw)*np.sqrt(2*np.pi)*lw*ref[i,1]/100
    return res

def analyze_isotopic_pattern(data, smiles, dm = 5, peak_width = 0.5, z=1,
                             maincolor='C0', theor_color='C1'):
    # smiles = "CN(C)CCCCCCCCCCCSSCCCCCCCCCCC[N+](C)(C)C"
    mymol = pybel.readstring('smi', smiles)
    print(mymol.formula)
    print(mymol.exactmass)
    f = Formula(mymol.formula)
    print(f.isotope.mass)
    spec = f.spectrum(minfract=1e-15)
    ref = [spec[x] for x in spec.keys()]
    print(spec)
    # ms_filename = 'data/diana/TMA ESI/TMA/+MS.txt'
    # data = np.loadtxt(ms_filename, delimiter='\t')
    # data[:,1] = data[:,1]/np.max(data[:,1])

    # ref_peak_mz, fit_mean, intensity = fit_mol(data, smiles, z=z)
    # mzerror = fit_mean-ref_peak_mz

    xdata = data[:,0]
    ydata = data[:,1]
    # ydata = data[np.logical_and(xdata > ref_peak_mz-peak_width, xdata < ref_peak_mz+dm),1]
    # xdata = xdata[np.logical_and(xdata > ref_peak_mz-peak_width, xdata < ref_peak_mz+dm)]
    mzmin = np.min(xdata)
    mzmax = np.max(xdata)
    ymax = np.max(ydata)

    intensity=1

    plt.plot(xdata, ydata/ymax, '-',alpha=1, label='Measured', color=maincolor)

    for i,r in enumerate(ref):
        mz, abun = r
        mz = mz/z
        intensity_here = intensity*abun/ref[0][1]/ymax
        if i==0:
            # label='Expected (Δm={0:.2f}mDa)'.format(mzerror*1000)
            label = 'Expected'
        else:
            label=None
        plt.plot([mz+mzerror, mz+mzerror], [0, intensity_here], color=theor_color, alpha=0.4, linewidth=10,
                 solid_capstyle="butt", label=label)
    plt.xlabel('m/z')
    plt.ylabel('Normalized intensity')
    plt.xlim(mzmin, mzmax)
    # plt.title(smiles)
    plt.legend()
    plt.tight_layout()
    plt.show()

def make_rainbow_table(maxHCNOS = [50, 20, 3, 0, 2]):
    masses = []
    formulas = []
    elements = ['H', 'C', 'N', 'O', 'S']
    for h in range(maxHCNOS[0] + 1):
        print('H={0}'.format(h))
        for c in range(maxHCNOS[1] + 1):
            for n in range(maxHCNOS[2] + 1):
                for o in range(maxHCNOS[3] + 1):
                    for s in range(maxHCNOS[4] + 1):
                        if h+c+n+o+s == 0:
                            continue
                        molecule = ''
                        counts = [h,c,n,o,s]
                        for k, element in enumerate(elements):
                            if counts[k] > 0:
                                molecule += element + str(counts[k])
                        f = Formula(molecule)
                        mass = f.isotope.mass
                        masses.append(mass)
                        formulas.append(molecule)

    masses = np.array(masses)
    np.save('data/rainbow/table1.npy', masses)
    pickle.dump(formulas, open('data/rainbow/table1_names.p', 'wb'))

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_values_close_to_other_list(source, target, thresh=0.04):
    res = []
    for v in source:
        nearest = find_nearest(target, v)
        if abs(nearest-v) < thresh:
            res.append(v)
    return res

def find_peaks_mz(data1, prominence=0.02):
    peaks, _ = find_peaks(data1[:,1], prominence=prominence)
    return data1[peaks,0]

def find_maching_peaks(data1, data2, prominence=0.02, thresh=0.01):
    peaks1 = find_peaks_mz(data1, prominence=prominence)
    peaks2 = find_peaks_mz(data2, prominence=prominence)
    matching = find_values_close_to_other_list(peaks1, peaks2, thresh=thresh)
    return matching

# make_rainbow_table(maxHCNOS = [50, 20, 3, 0, 2])
if __name__ == '__main__':
    pass
    # rainbow_masses = np.load('data/rainbow/table1.npy')
    # rainbow_name = pickle.load(open('data/rainbow/table1_names.p', 'rb'))

    # def find_mf(mass):
    #     index = np.argmin(np.abs(rainbow_masses-mass))
    #     return rainbow_name[index],rainbow_masses[index]

    # targmass = 200.2727
    # formula, mass = find_mf(targmass)
    # print(formula)
    # print(mass)
    # print((targmass-mass)*1000)

# # Import the package into python
# import pyisopach
# # Create Molecule object for Sucrose
# mol = pyisopach.Molecule("C56H53N2O18S2")
# # Return molecular weight
# mol.molecular_weight
# # Calculate isotopic distribution/pattern
# xxx = np.array(mol.isotopic_distribution(electrons=0)).T
# #
#

#
#
# xdata = data[:,0]
# ydata = data[np.logical_and(xdata > ref_peak_mz-peak_width, xdata < ref_peak_mz+peak_width),1]
# xdata = xdata[np.logical_and(xdata > ref_peak_mz-peak_width, xdata < ref_peak_mz+peak_width)]
# xdata = data[:,0]
# ydata = data[np.logical_and(xdata > m_range[0], xdata < m_range[1]),1]
# xdata = xdata[np.logical_and(xdata > m_range[0], xdata < m_range[1])]
# plt.plot(xdata, ydata, 'o',alpha=0.5, label='Measured')
# xs = np.linspace(m_range[0], m_range[1], 10000)
# plt.plot(xs, sim_spec(xs, ref), label='Expected')
# plt.xlabel('m/z')
# plt.ylabel('Counts')
# plt.title(targ_mol)
# plt.legend()
# plt.show()
# fit LW
# from scipy.optimize import curve_fit
# targ_ms = np.array([[x, data[i,1]] for i,x in enumerate(data[:,0]) if x>m_range[0] and x<m_range[1]])
#
# def func(x, a, b, c, d, e):
#     return c*sim_spec(x+b, refs[targ_mol], lw=a, shape=e) + c*d*sim_spec(x+b+2*1.00783, refs['C35H58N2OSH'], lw=a, shape=e)
#
# popt, pcov = curve_fit(func, targ_ms[:,0], targ_ms[:,1], p0=[0.2, 0, 1, 0.12, 0.2], maxfev=1000,
#                        bounds=[[-np.inf, 0, -np.inf, 0, 0], [np.inf, 1e-5, np.inf, np.inf, np.inf]])
# print(popt)
# xs = np.linspace(m_range[0], m_range[1], 500)
# ys = func(xs, popt[0], popt[1], popt[2], popt[3], popt[4])
# # ys = ys/np.max(ys)
# plt.stem(refs[targ_mol][:,0], refs[targ_mol][:,1]/100, markerfmt=' ', basefmt='none', linefmt='grey',
#          label='Expected, mass & abundance')
# plt.stem(refs[targ_mol][:,0]-2*1.00783, refs[targ_mol][:,1]/100*popt[3], markerfmt=' ', basefmt='none', linefmt='grey')
# plt.plot(xs, ys, alpha=0.5, label='Expected')
# plt.xlim(m_range[0], m_range[1])
# plt.legend()
# plt.xlabel('m/z')
# plt.show()