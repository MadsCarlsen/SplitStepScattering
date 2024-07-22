import numpy as np 
import matplotlib.pyplot as plt 
from scipy.io import FortranFile
from numpy.fft import fftshift, fftfreq
from matplotlib.colors import LogNorm

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

def load_dims(path='scatter_pots/dim.txt'): 
    data = np.loadtxt(path)
    Nx = int(data[0])
    Ny = int(data[1])
    x_low = data[6]
    x_high = data[7]
    y_low = data[8]
    y_high = data[9]
    return Nx, Ny, x_low, x_high, y_low, y_high, 

Nx, Ny, x_low, x_high, y_low, y_high = load_dims()

x_list = np.linspace(x_low, x_high, Nx)
y_list = np.linspace(y_low, y_high, Ny)
dx = x_list[1] - x_list[0]
dy = y_list[1] - y_list[0]

# Real space wave functions
V00_r = load_fortran_binary('scatter_pots/V00_r.dat', dim=Ny)
V00_i = load_fortran_binary('scatter_pots/V00_i.dat', dim=Ny)
V11_r = load_fortran_binary('scatter_pots/V11_r.dat', dim=Ny)
V11_i = load_fortran_binary('scatter_pots/V11_i.dat', dim=Ny)
V01_r = load_fortran_binary('scatter_pots/V01_r.dat', dim=Ny)
V01_i = load_fortran_binary('scatter_pots/V01_i.dat', dim=Ny)

V00 = V00_r + 1j*V00_i
V11 = V11_r + 1j*V11_i
V01 = V01_r + 1j*V01_i

del(V00_r, V00_i, V11_r, V11_i, V01_r, V01_i)

vmax, vmin = 1, 1e-8

XX, YY = np.meshgrid(x_list, y_list, indexing='ij')

fig, ax = plt.subplots()
ax.plot(y_list, np.real(V00[Nx//2,:]), c='C0')
ax.plot(x_list, np.real(V00[:,Ny//2]), c='C0', ls='--')
ax.plot(y_list, np.real(V11[Nx//2,:]), c='C1')
ax.plot(x_list, np.real(V11[:,Ny//2]), c='C1', ls='--')
ax.plot(y_list, np.real(V01[Nx//2,:]), c='C2')
ax.plot(x_list, np.real(V01[:,Ny//2]), c='C2', ls='--')


fig, axs = plt.subplots(1,2, figsize=(10, 5))
im0 = axs[0].pcolormesh(XX, YY, np.real(V00), cmap='turbo')# norm=LogNorm(vmin=vmin, vmax=vmax))
im1 = axs[1].pcolormesh(XX, YY, np.imag(V00), cmap='turbo')
#im0 = axs[0].pcolormesh(XX, YY, np.real(psi1), cmap='turbo')# norm=LogNorm(vmin=vmin, vmax=vmax))
#im1 = axs[1].pcolormesh(XX, YY, np.imag(psi1), cmap='turbo')
fig.colorbar(im0, ax=axs[0])
fig.colorbar(im1, ax=axs[1])

fig, axs = plt.subplots(1,2, figsize=(10, 5))
im0 = axs[0].pcolormesh(XX, YY, np.real(V11), cmap='turbo')# norm=LogNorm(vmin=vmin, vmax=vmax))
im1 = axs[1].pcolormesh(XX, YY, np.imag(V11), cmap='turbo')
#im0 = axs[0].pcolormesh(XX, YY, np.real(psi1), cmap='turbo')# norm=LogNorm(vmin=vmin, vmax=vmax))
#im1 = axs[1].pcolormesh(XX, YY, np.imag(psi1), cmap='turbo')
fig.colorbar(im0, ax=axs[0])
fig.colorbar(im1, ax=axs[1])

fig, axs = plt.subplots(1,2, figsize=(10, 5))
im0 = axs[0].pcolormesh(XX, YY, np.real(V01), cmap='turbo')# norm=LogNorm(vmin=vmin, vmax=vmax))
im1 = axs[1].pcolormesh(XX, YY, np.imag(V01), cmap='turbo')
#im0 = axs[0].pcolormesh(XX, YY, np.real(psi1), cmap='turbo')# norm=LogNorm(vmin=vmin, vmax=vmax))
#im1 = axs[1].pcolormesh(XX, YY, np.imag(psi1), cmap='turbo')
fig.colorbar(im0, ax=axs[0])
fig.colorbar(im1, ax=axs[1])


plt.show()