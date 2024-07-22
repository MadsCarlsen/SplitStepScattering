import numpy as np 
import matplotlib.pyplot as plt 
from scipy.io import FortranFile
from numpy.fft import fftshift, fftfreq
from matplotlib.colors import LogNorm

# PLOT STUFF 
major = 7.5
minor = 3
width = 1.25
plt.rcParams["figure.figsize"] = (6,6)
plt.rc('text', usetex=True)
plt.rc("axes", labelsize=16) # 18
plt.rc("xtick", labelsize=16, top=True, direction="in")
plt.rc("ytick", labelsize=16, right=True, direction="in")
plt.rc("axes", titlesize=16)
plt.rc("legend", fontsize=14)
plt.rcParams['font.family'] = "serif"
plt.rcParams['axes.linewidth'] = width
plt.rcParams['xtick.minor.width'] = width
plt.rcParams['xtick.major.width'] = width
plt.rcParams['ytick.minor.width'] = width
plt.rcParams['ytick.major.width'] = width
plt.rcParams['xtick.major.size'] = major
plt.rcParams['xtick.minor.size'] = minor
plt.rcParams['ytick.major.size'] = major
plt.rcParams['ytick.minor.size'] = minor

def trapz_2D(f, x_list, y_list): 
    return np.trapz(np.trapz(f, y_list,), x_list)

# Loads a binary fortran data file 
def load_fortran_binary(path, dim=1): 
    f = FortranFile(path, 'r')
    data = f.read_reals(float)
    dim1 = int(len(data)/dim) 
    data = data.reshape((dim1, dim), order="F")  # Reshaping for 2D arrays 
    f.close()
    return data

def load_dims(path='data/dim_imag.txt'): 
    data = np.loadtxt(path)
    Nx = int(data[0])
    Ny = int(data[1])
    x_low = data[2]
    x_high = data[3]
    y_low = data[4]
    y_high = data[5]
    dkx = data[6]
    dky = data[7]
    kx_max = data[8]
    ky_max = data[9]
    return Nx, Ny, x_low, x_high, y_low, y_high, dkx, dky, kx_max, ky_max

def V(x, y):
    a = 0.63  # 0.63 is H like 
    r = x**2 + y**2
    return -1 / np.sqrt(r + a) 

Nx, Ny, x_low, x_high, y_low, y_high, dkx, dky, kx_max, ky_max = load_dims()

x_list = np.linspace(x_low, x_high, Nx)
y_list = np.linspace(y_low, y_high, Ny)
dx = x_list[1] - x_list[0]
dy = y_list[1] - y_list[0]

kx_list = np.linspace(-kx_max, kx_max, Nx)
ky_list = np.linspace(-ky_max, ky_max, Ny)

# Real space wave functions
psi0_r = load_fortran_binary('data/imag_psi0_r.dat', dim=Ny)
psi0_i = load_fortran_binary('data/imag_psi0_i.dat', dim=Ny)
psi1_r = load_fortran_binary('data/imag_psi1_r.dat', dim=Ny)
psi1_i = load_fortran_binary('data/imag_psi1_i.dat', dim=Ny)
psi0 = psi0_r + 1j*psi0_i
psi1 = psi1_r + 1j*psi1_i

# Momentum space wave functions 
#psi0_r = load_fortran_binary('data/imag_psi0k_r.dat', dim=Ny)
#psi0_i = load_fortran_binary('data/imag_psi0k_i.dat', dim=Ny)
#psi1_r = load_fortran_binary('data/imag_psi1k_r.dat', dim=Ny)
#psi1_i = load_fortran_binary('data/imag_psi1k_i.dat', dim=Ny)
#psi0_k = psi0_r + 1j*psi0_i
#psi1_k = psi1_r + 1j*psi1_i

del(psi0_r, psi0_i, psi1_r, psi1_i)

# Test normalization
print('psi0 norm: ', trapz_2D(np.abs(psi0)**2, x_list, y_list))
print('psi1 norm: ', trapz_2D(np.abs(psi1)**2, x_list, y_list))
#print('psi0_k norm: ', trapz_2D(np.abs(psi0_k)**2, kx_list, ky_list))
#print('psi1_k norm: ', trapz_2D(np.abs(psi1_k)**2, kx_list, ky_list))

vmax, vmin = 1, 1e-8

XX, YY = np.meshgrid(x_list, y_list, indexing='ij')
#KX, KY = np.meshgrid(kx_list, ky_list, indexing='ij')

V_grid = V(XX, YY)

fig, axs = plt.subplots(1,2, figsize=(10,4))
axs[0].contour(XX, YY, V_grid, 10, origin='lower', extent=[x_list[0], x_list[-1], y_list[0], y_list[-1]], colors=['w'], zorder=10, alpha=0.3)
im = axs[0].pcolormesh(XX, YY, np.abs(psi0)**2, cmap='inferno', zorder=0, shading='gouraud', rasterized=True)#, vmin=vmin, vmax=vmax)
axs[0].axis('equal')
axs[0].set_xlim(-10,10)
axs[0].set_ylim(-10,10)
axs[0].set_xlabel('$x$ (a.u.)')
axs[0].set_ylabel('$y$ (a.u.)', labelpad=-10)
cbar = fig.colorbar(im, ax=axs[0])

axs[0].set_xticks([-10, -5, 0, 5, 10])
axs[0].tick_params(axis='both', which='both', color='w', width=0.75)


axs[1].contour(XX, YY, V_grid, 10, origin='lower', extent=[x_list[0], x_list[-1], y_list[0], y_list[-1]], colors=['w'], zorder=10, alpha=0.3)
im1 = axs[1].pcolormesh(XX, YY, np.abs(psi1)**2, cmap='inferno', zorder=0, shading='gouraud', rasterized=True)#, vmin=vmin, vmax=vmax)

axs[1].axis('equal')
axs[1].set_xlim(-10,10)
axs[1].set_ylim(-10,10)
axs[1].set_xlabel('$x$ (a.u.)')
axs[1].set_ylabel('$y$ (a.u.)', labelpad=-10)
cbar1 = fig.colorbar(im1, ax=axs[1])
axs[1].set_xticks([-10, -5, 0, 5, 10])
axs[1].tick_params(axis='both', which='both', color='w', width=0.9)


axs[0].text(0.1, 0.85, '(a)', color='w', horizontalalignment='left', verticalalignment='center', transform=axs[0].transAxes, fontsize=20)
axs[0].text(0.95, 0.85, r'$| \langle \mathbf{r} | 0 \rangle |^2$', color='w', horizontalalignment='right', verticalalignment='center', transform=axs[0].transAxes, fontsize=20)

axs[1].text(0.1, 0.85, '(b)', color='w', horizontalalignment='left', verticalalignment='center', transform=axs[1].transAxes, fontsize=20)
axs[1].text(0.95, 0.85, r'$| \langle \mathbf{r} | 1 \rangle |^2$', color='w', horizontalalignment='right', verticalalignment='center', transform=axs[1].transAxes, fontsize=20)

#axs[0].set_title(r'$| \langle \mathbf{r} | i \rangle |^2$', color='k', rotation='vertical',x=-0.3,y=0.3, fontsize=20)

plt.tight_layout()
#plt.savefig('ss_eigenstates.pdf', bbox_inches='tight')
plt.show()
plt.show()