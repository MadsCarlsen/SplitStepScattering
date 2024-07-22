# SplitStepScattering
Code for simulating scattering of an electron wave packet (WP) on a two level system (2LS) in the presence of a laser field. The WP is represented on a 2D Cartesian grid and propagation is performed using split-step FFT.

## Breif description
The code contains three programs, found in the *soruce* directory. All three are needed to perfrom a simulation, and be used in the following order of presentation.

### Imaginary propagation
***imag_prop.f90*** is a program for determining the ground state and first excited state of a potential through imaginary time propagation. The standard potential is a soft Coulomb potential, but other single active electron potentials can be implemented. Settings for imaginary time propagation is set in the *settings_imag.nml* file, which also contains a (very) short description of the various parameters.

The output wave functions (WF) are saved in *data* with the prefix *imag*. The WFs can be visualized using the *plot_eigenstates.py* script.

### Calculation of scattering potentials 
***calculate_scatt_pots.f90*** is a program for calculating the potentials an external projectile electron will feel due to the precense of the 2LS. The code assumes the target is a core with charge +1 and an electron with charge -1 (i.e. hydrogen like), but this can easily be modified. Note, the code does not account for exchange effects.

The program uses the eigenstates found thorugh imaginary time propagation as input, and outputs the scattering potentials for the different channels in the directory *scatter_pots*. The potentials are calculated on the grid used in the real time propagation. The scattering potentials can be visualized using the *plot_scatt_pots.py* script.

### Real time propagation
***split_step_A_field.f90*** is a program for real time propagation of a projectile electron WP in the presence of a 2LS and an electric field. The 2LS interaction happens through the scattering potentials found above. The laser field is represented by its vector potential in the dipole approximation. Settings for real time propagation are set in *settings.nml*, which again contains a short description of the various parameters.

Currently the code only works when the entire laser pules is contained within the simulation time. Furthermore the field is assumed to be linearly polarized along the y-direction, but this can easily be changed.

The code contains two standard choices for the projectile WP: A standard Gaussian WP or a WP laser modulated in the longitudinal (x) direction. 

As output the momentum space WFs in the two channels are saved at the end of propagation in *data*. Furthermore the corresponding 'free' solutions (when no 2LS is present) are saved. The momentum distributions can be visualized using the *plot_momentum_dist.py* script.

NB: Running this script requires providing the setting file as a command line input, i.e. `./split_step.exe settings.nml`. This allows one to save several different settings files. Note, however that if the grid is changed, new scattering potentials must be calculated for each configuration.


## Dependencies and compilation
- GNU Fortran (tested on Ubuntu gfortran 11.4.0)
- OpenMP
- FFTW3

Compilation should be as simple as running the bash script *compile*. 


