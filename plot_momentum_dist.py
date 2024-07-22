import numpy as np 
import matplotlib.pyplot as plt
from scipy.io import FortranFile
import matplotlib.colors as colors

# PLOT STUFF 
major = 7.5
minor = 3
width = 1.25
plt.rcParams["figure.figsize"] = (6,6)
plt.rc('text', usetex=True)
plt.rc("axes", labelsize=18) # 18
plt.rc("xtick", labelsize=16, top=True, direction="in")
plt.rc("ytick", labelsize=16, right=True, direction="in")
plt.rc("axes", titlesize=18)
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


# 2D trapz integral
def trapz_2D(f, x_list, y_list):
    return np.trapz(np.trapz(f, y_list), x_list)

# Loads a binary fortran data file 
def load_fortran_binary(path, dim=1): 
    f = FortranFile(path, 'r')
    data = f.read_reals(float)
    dim1 = int(len(data)/dim) 
    data = data.reshape((dim1, dim), order="F")  # Reshaping for 2D arrays 
    f.close()
    return data

def load_dim(path): 
    dim_data = np.loadtxt(path)
    kx_lower = dim_data[0]
    kx_upper = dim_data[1]
    Nx = int(dim_data[2])
    ky_lower = dim_data[3]
    ky_upper = dim_data[4]
    Ny = int(dim_data[5])
    return kx_lower, kx_upper, Nx, ky_lower, ky_upper, Ny 

def load_WP(wp_name, dim, path='simulation_data/'): 
    data_r = load_fortran_binary(path + wp_name + '_r.dat', dim=dim)
    data_i = load_fortran_binary(path + wp_name + '_i.dat', dim=dim)
    return data_r + 1j * data_i

path = 'data/' 

kx_lower, kx_upper, Nx, ky_lower, ky_upper, Ny = load_dim(path + 'dim.txt')

wp0_name = 'psi0_k'
wp1_name = 'psi1_k'
wp0_free_name = 'psi0_k_free'

print(Nx, Ny)


data0 = load_WP(wp0_name, dim=Ny, path=path)
data0_free = load_WP(wp0_free_name, dim=Ny, path=path)
data1 = load_WP(wp1_name, dim=Ny, path=path)

kx_list = np.linspace(kx_lower, kx_upper, Nx)
ky_list = np.linspace(ky_lower, ky_upper, Ny)
KX, KY = np.meshgrid(kx_list, ky_list, indexing='ij')

# Try to integrate and see if norm is conserved 
P0 = trapz_2D(np.abs(data0)**2, kx_list, ky_list)
P1 = trapz_2D(np.abs(data1)**2, kx_list, ky_list)
print(P0, P1, P0 + P1)

fig, axs = plt.subplots(1,2, figsize=(10,4))#, sharex=True)
#im = ax.pcolormesh(KX, KY, np.abs(data), cmap='inferno_r', shading='gouraud') 
im = axs[0].pcolormesh(KX, KY, np.abs(data0)**2, cmap='turbo', shading='gouraud', rasterized=True) 
#axs[0].axis('equal')
#ax.set_aspect('equal', adjustable='box')
axs[0].set_xlim(3.5, 4.5)
axs[0].set_ylim(-1., 1.)
cbar = fig.colorbar(im, ax=axs[0])
axs[0].set_xlabel('$k_x$ (a.u.)')
axs[0].set_ylabel('$k_y$ (a.u.)')
#axs[0].minorticks_on()
#axs[0].tick_params(axis='y', which='minor', left=False, right=False)
axs[0].tick_params(axis='both', which='both', color='w', width=0.9)


#im = ax.pcolormesh(KX, KY, np.abs(data), cmap='inferno_r', shading='gouraud') 
im = axs[1].pcolormesh(KX, KY, np.abs(data1)**2, cmap='turbo', shading='gouraud', rasterized=True)
#im = axs[1].pcolormesh(KX, KY, np.abs(data1)**2 + np.abs(data0)**2 - np.abs(data1_nofield)**2 - np.abs(data0_nofield)**2, shading='gouraud', rasterized=True, cmap='RdBu', norm=colors.CenteredNorm())
#im = axs[1].pcolormesh(KX, KY, np.abs(data1)**2 - np.abs(data1_nofield)**2, cmap='turbo', shading='gouraud', rasterized=True)
#im = axs[1].pcolormesh(KX, KY, np.abs(data0)**2 - np.abs(data0_nofield)**2, cmap='turbo', shading='gouraud', rasterized=True)

#axs[1].axis('equal') 
#ax.set_aspect('equal', adjustable='box')
axs[0].set_xlim(3.5, 4.5)
axs[0].set_ylim(-1., 1.)
cbar = fig.colorbar(im, ax=axs[1])
axs[1].set_xlabel('$k_x$ (a.u.)')
axs[1].set_ylabel('$k_y$ (a.u.)')
axs[1].minorticks_on()
axs[1].tick_params(axis='y', which='minor', left=False, right=False)
axs[1].tick_params(axis='both', which='both', color='w', width=0.9)

axs[0].text(0.08, 0.875, '(a)', color='w', horizontalalignment='left', verticalalignment='center', transform=axs[0].transAxes, fontsize=20)
axs[0].text(0.92, 0.875, '$|\psi_{0}|^2$', color='w', horizontalalignment='right', verticalalignment='center', transform=axs[0].transAxes, fontsize=20)

axs[1].text(0.08, 0.875, '(b)', color='w', horizontalalignment='left', verticalalignment='center', transform=axs[1].transAxes, fontsize=20)
axs[1].text(0.92, 0.875, '$|\psi_{1}|^2$', color='w', horizontalalignment='right', verticalalignment='center', transform=axs[1].transAxes, fontsize=20)

plt.tight_layout()
#plt.savefig('momentum_dist.pdf', bbox_inches='tight')
plt.show()