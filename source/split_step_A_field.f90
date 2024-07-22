program split_step_scattering
use, intrinsic :: iso_c_binding
use utils
use omp_lib
implicit none
include "fftw3.f03" 

integer, parameter :: wp = dp 

! Define parameters to be loaded 
integer :: Nx, Ny
real(wp) :: x_low, x_high, y_low, y_high, r_max, r_boundary 
real(wp) :: dt, t_interract
real(wp) :: mask_bound, mask_max
integer :: cos_power
real(wp) :: x0, y0, k0x, k0y, sigma_time, sigma_y, g_factor, omega
integer :: max_N 
logical :: PINEM_packet 
real(wp) :: theta, phi 
real(wp) :: I0, cep, omega_l
integer :: N_cycles

real(wp) :: E0, E1  ! Energies of the 2LS

! Other variables 
real(wp) :: period, sigma_x, sigma_kx, sigma_ky, t_final, t_bunch, k_2LS_shift
real(wp), allocatable :: x_list(:), y_list(:), kx_list(:), ky_list(:)
real(wp) :: dx, dy, dkx, dky, kx_max, ky_max 
complex(wp), allocatable :: psi0(:,:), psi1(:,:), psi0_new(:,:), psi1_new(:,:), psi0_save(:,:), psi1_save(:,:)  ! State matrices
real(wp), allocatable :: V00(:,:), V01(:,:), V11(:,:)  ! Matrices to store the potential 
real(wp) :: r_pot_max 
complex(wp), allocatable :: UT(:,:), sin_fac(:,:), cos_fac(:,:)  ! Matrices to store propagation operators
real(wp), allocatable :: r_mat(:,:), D_sqrt(:,:), V_diff(:,:), mask(:)
integer :: i,j, Nt 
real(wp) :: r

! Variables for field propagation 
real(wp) :: t_i, D_sqrt_T, sin_T, cos_T, A_ti 
complex(wp) :: dip01, dip10
real(wp) :: dip_r, dip_i  ! Used for loading 
real(wp) :: rtUp, Up, E_amp, laser_duration, t_shift 

! Variables for saved data 
real(wp), allocatable :: kx_save(:), ky_save(:), diff_arr(:)
real(wp) :: kx_lower, kx_upper, ky_lower, ky_upper 
integer :: kx_low_i, kx_high_i, ky_low_i, ky_high_i, Nx_save, Ny_save  

! FFT variables 
type(c_ptr) :: plan0_f, plan0_b, plan1_f, plan1_b
integer :: void 

! LOAD SETTINGS  
character(len=100) :: configFileName
integer :: ierr
namelist /SPLIT_STEP_GRID/ Nx, Ny, x_low, x_high, y_low, y_high, dt, r_max, r_boundary, mask_bound, mask_max, cos_power

namelist /SPLIT_STEP_OTHER/ t_interract, omega, PINEM_packet, theta, phi, x0, y0, k0x, k0y, &
                            sigma_time, sigma_y, g_factor, max_N, omega_l, I0, cep, N_cycles, &
                            kx_lower, kx_upper, ky_lower, ky_upper 

! Check if the command line argument is provided
if (command_argument_count() /= 1) then
    write(*,*) 'Provide config file name!'
else
    ! Get the command line argument
    call get_command_argument(1, configFileName, ierr)
end if 

open(file=configFileName, unit=1)
read(nml=SPLIT_STEP_GRID, unit=1)
read(nml=SPLIT_STEP_OTHER, unit=1)
close(1)

write(*,*) 'Settings loaded from file: ', trim(configFileName)

! SETUP MULTITHREADING FFT 
call omp_set_num_threads(2)  ! <---- Change this for other nr. of cores! 
void = fftw_init_threads()
call fftw_plan_with_nthreads(omp_get_max_threads())
write(*,*) 'Using threads: ', omp_get_max_threads()


! SETUP GRIDS 
allocate(x_list(Nx), y_list(Ny), kx_list(Nx), ky_list(Ny))
call linspace(x_list, x_low, x_high, Nx)
call linspace(y_list, y_low, y_high, Ny)

dx = x_list(2) - x_list(1)
dy = y_list(2) - y_list(1)
dkx = 2._wp*pi / (Nx*dx)
dky = 2._wp*pi / (Ny*dy)
kx_max = (Nx-1._wp)/2._wp * dkx
ky_max = (Ny-1._wp)/2._wp * dky

do i=1, Nx 
    kx_list(i) = -kx_max + (i-1._wp) * dkx
end do
do j=1, Ny 
    ky_list(j) = -ky_max + (j-1._wp) * dky
end do


! LOAD AND CALCULATE POTENTIAL MATRICES 
allocate(V00(Nx,Ny))
allocate(V01(Nx,Ny))
allocate(V11(Nx,Ny))

! Build matrix containing radial values (used later as mask anyway)
allocate(r_mat(Nx,Ny))
do j=1, Ny 
    do i=1, Nx
        r_mat(i,j) = sqrt(x_list(i)**2 + (y_list(j)-y0)**2)  ! Introduce y0 shift in potential (impact parameter of WP)
    end do 
end do

call load_scatt_pots(V00, V01, V11, E0, E1, dip10, r_pot_max, Nx, Ny)
dip01 = conjg(dip10)

write(*,*) 'k0x : ', k0x 
write(*,*) 'Energies :', E0, E1
write(*,*) 'Dipole is ', dip10 
write(*,*) 'r_pot_max : ', r_pot_max

! Save the potential grids to check they are ok? 
!open(file='data/V00.dat', unit=1, form='unformatted')
!open(file='data/V01.dat', unit=2, form='unformatted')
!open(file='data/V11.dat', unit=3, form='unformatted')
!write(1) V00 
!write(2) V01
!write(3) V11
!close(1)
!close(2)
!close(3)


! CALCULATE DIFFERENT VARIABLES
! WP parameters 
period = 2._wp*pi / (E1-E0)
sigma_x = sigma_time*k0x
sigma_kx = 1._wp/(2._wp * sigma_x)
sigma_ky = 1._wp/(2._wp * sigma_y)
k_2LS_shift = sqrt(k0x**2 + 2._wp*(E1-E0)) - k0x  

! Field settings 
cep = cep * pi 
!omega_l = 2._wp*pi * 137.036_wp / (lambda * 1.0e-9_wp / 5.29177e-11_wp)
E_amp = sqrt(I0 / 3.50945e16_wp)  ! Maximum electric field amplitude in a.u.
Up = E_amp**2 / (4._wp * omega_l**2)  ! Ponderomotive energy in a.u.
rtUp = sqrt(Up)
laser_duration = 2._wp*N_cycles*pi/omega_l 
write(*,*) 'Up :', Up

! Time parameters 
t_interract = t_interract * sigma_x/k0x 
t_final = 2._wp * t_interract 
t_bunch = k0x**2 / (g_factor * omega**2)  
Nt = int(t_final / dt) 

! Check that laser fits within simulation time 
if (2._wp * t_interract > laser_duration) then 
    t_shift = t_interract - laser_duration/2._wp  ! Center laser pulse on WP collision time 
else 
    write(*,*) 'Laser field does not fit within simulation time!'
    stop 
end if 
!t_shift = 0._wp  ! No shift for now 

! 2LS parameters
theta = theta * pi 
phi = phi * pi 

write(*,*) 't_final : ', t_final


! CALCULATE THE PROPAGATION OPERATORS 
allocate(UT(Nx,Ny), D_sqrt(Nx,Ny), sin_fac(Nx,Ny), cos_fac(Nx,Ny), V_diff(Nx,Ny))

! The kinetic operator (split the exp calculation to avoid expensive calls)
do i=1, Nx 
    UT(i,:) = exp(cmplx(0._wp, -0.5_wp*kx_list(i)**2*dt, wp))
end do 
do j=1, Ny 
    UT(:,j) = UT(:,j) * exp(cmplx(0._wp, -0.5_wp*ky_list(j)**2*dt, wp))
end do
UT = UT * exp(cmplx(0._wp, -(E0+E1)*dt/2._wp, wp))

! Potential operator 
dt = 0.5_wp * dt   ! Should really just use dt/2 in the following, but eh... 
D_sqrt = sqrt(4._wp * abs(V01)**2 + (V00-V11)**2)
cos_fac = exp(cmplx(0._wp, -0.5_wp*(V00+V11)*dt, wp)) * cos(0.5_wp*D_sqrt*dt)

! Handle divergence in sin_fac when potential is zero 
sin_fac = exp(cmplx(0._wp, -0.5_wp*(V00+V11)*dt, wp))
where (r_mat <= r_pot_max)
    sin_fac = sin_fac * cmplx(0._wp, sin(0.5_wp*D_sqrt*dt) / D_sqrt, wp)
elsewhere 
    sin_fac = sin_fac * cmplx(0._wp, 0.5_wp * dt, wp)
end where 

V_diff = V11 - V00
dt = 2._wp * dt  

! Now free up memory from the arrays we don't need anymore 
deallocate(V00, V11, D_sqrt, r_mat)


! CALCULATE MASK TO REMOVE WAVE FUNCTION IN Y-DIRECTION
allocate(mask(Ny))
mask = 1._wp  ! 1 is no absorption  
do j=1, Ny 
    if (abs(y_list(j)) < mask_bound) then 
        cycle 
    else if (y_list(j) >= mask_bound) then 
        mask(j) = 1._wp - cos(0.5_wp*pi * (mask_max - y_list(j)) / (mask_max-mask_bound))**cos_power
    else if (y_list(j) <= -mask_bound) then 
        mask(j) = 1._wp - cos(0.5_wp*pi * (mask_max + y_list(j)) / (mask_max-mask_bound))**cos_power
    end if 
end do


! CALCULATE SLICES TO SAVE 
allocate(diff_arr(Nx))
kx_low_i = minloc(abs(kx_list - kx_lower), 1)
kx_high_i = minloc(abs(kx_list - kx_upper), 1)
ky_low_i = minloc(abs(ky_list - ky_lower), 1)
ky_high_i = minloc(abs(ky_list - ky_upper), 1)
kx_save = kx_list(kx_low_i:kx_high_i)
ky_save = ky_list(ky_low_i:ky_high_i)
Nx_save = size(kx_save)
Ny_save = size(ky_save)


! PREPARE THE FFT. Must happen BEFORE initialization of initial states. 
write(*,*) 'Planning FFTs...'
allocate(psi0(Nx,Ny), psi1(Nx,Ny), psi0_new(Nx,Ny), psi1_new(Nx,Ny))
plan0_f = fftw_plan_dft_2d(Ny, Nx, psi0_new, psi0, FFTW_FORWARD, FFTW_MEASURE)  !PATIENT is faster than MEASURE overall (but takes longer to prepare) 
plan1_f = fftw_plan_dft_2d(Ny, Nx, psi1_new, psi1, FFTW_FORWARD, FFTW_MEASURE)

plan0_b = fftw_plan_dft_2d(Ny, Nx, psi0_new, psi0_new, FFTW_BACKWARD, FFTW_MEASURE)
plan1_b = fftw_plan_dft_2d(Ny, Nx, psi1_new, psi1_new, FFTW_BACKWARD, FFTW_MEASURE)


! CALCULATE THE INTIAL WAVEFUNCTIONS 
if (PINEM_packet) then 
    write(*,*) 'Performing calculations for a PINEM WP'
else 
    write(*,*) 'Performing calculations for a Gauss WP'
end if

! This part is slow, might as well do it in parallel...
!$OMP PARALLEL DO PRIVATE(i)
do j=1, Ny 
    do i=1, Nx 
        if (PINEM_packet) then 
            ! PINEM WP 
            psi0(i,j) = PINEM_WP(x_list(i), t_bunch-t_interract, t_bunch*k0x, k0x, sigma_kx, omega, max_N, g_factor) *& 
                        Gauss_WP(y_list(j), -t_interract, 0._wp, 0._wp, sigma_ky)
            psi1(i,j) = psi0(i,j) * exp(cmplx(0._wp, phi, wp)) * cos(theta/2._wp)  ! Adding coefficients to populate 2LS 
            psi0(i,j) = psi0(i,j) * sin(theta/2._wp)
        else  
            ! GAUSS WP (NB: longitudinal waist happens at target, unlike PINEM WP above!) 
            psi0(i,j) = Gauss_WP(x_list(i), -t_interract, 0._wp, k0x, sigma_kx) * &
                        Gauss_WP(y_list(j), -t_interract, 0._wp, k0y, sigma_ky)
            psi1(i,j) = psi0(i,j) * exp(cmplx(0._wp, phi, wp)) * cos(theta/2._wp)  ! Adding coefficients to populate 2LS 
            psi0(i,j) = psi0(i,j) * sin(theta/2._wp)
        end if 
    end do 
end do 
!$OMP END PARALLEL DO 

! Now send the WF into the FFT-ready form 
call convert_space_to_FFT(psi0, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
call convert_space_to_FFT(psi1, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)


! FREE PROPAGATION OF INTIAL STATE 
! Go to momentum space 
psi0_new = psi0
psi1_new = psi1
call fftw_execute_dft(plan0_f, psi0_new, psi0)
call fftw_execute_dft(plan1_f, psi1_new, psi1)
call convert_FFT_to_momentum(psi0, x_list(1), y_list(1), dx, dy, dkx, dky, Nx, Ny)
call convert_FFT_to_momentum(psi1, x_list(1), y_list(1), dx, dy, dkx, dky, Nx, Ny)

! First phases of 2LS 
psi0 = psi0 * exp(cmplx(0._wp, -E0*Nt*dt, wp))
psi1 = psi1 * exp(cmplx(0._wp, -E1*Nt*dt, wp))

! Then time phases of the WP 
do i=1, Nx 
    psi0(i,:) = psi0(i,:) * exp(cmplx(0._wp, -kx_list(i)**2/2._wp * Nt*dt, wp))
    psi1(i,:) = psi1(i,:) * exp(cmplx(0._wp, -kx_list(i)**2/2._wp * Nt*dt, wp))
end do
do j=1, Ny 
    psi0(:,j) = psi0(:,j) * exp(cmplx(0._wp, -ky_list(j)**2/2._wp * Nt*dt, wp))
    psi1(:,j) = psi1(:,j) * exp(cmplx(0._wp, -ky_list(j)**2/2._wp * Nt*dt, wp))
end do

! Save the initial state in k-space 
allocate(psi0_save(Nx_save, Ny_save))
allocate(psi1_save(Nx_save, Ny_save))

psi0_save = psi0(kx_low_i:kx_high_i, ky_low_i:ky_high_i)
psi1_save = psi1(kx_low_i:kx_high_i, ky_low_i:ky_high_i)

open(file='data/psi0_k_free_r.dat', unit=1, form='unformatted')
open(file='data/psi0_k_free_i.dat', unit=2, form='unformatted')
open(file='data/psi1_k_free_r.dat', unit=3, form='unformatted')
open(file='data/psi1_k_free_i.dat', unit=4, form='unformatted')
write(1) real(psi0_save)
write(2) aimag(psi0_save) 
write(3) real(psi1_save)
write(4) aimag(psi1_save)
close(1)
close(2)
close(3)
close(4)
deallocate(psi0_save, psi1_save)
psi0 = psi0_new  ! Reset to be ready for full propagation 
psi1 = psi1_new


! FULL PROPAGATION USING SPLIT-STEP FFT 
do i=1, Nt 
    write(*,*) i, '/', Nt
    t_i = (i-0.5_wp)*dt  ! Time value at half time step (where time dependent operators are evaluated)

    ! Half step in real space 
    psi0_new = cos_fac*psi0 + sin_fac*(V_diff*psi0 - 2._wp*V01*psi1)
    psi1_new = cos_fac*psi1 + sin_fac*(-V_diff*psi1 - 2._wp*V01*psi0)

    ! Go to momentum space 
    call fftw_execute_dft(plan0_f, psi0_new, psi0)
    call fftw_execute_dft(plan1_f, psi1_new, psi1)
    
    ! First calculate A-field operators
    A_ti = A_field(t_i - t_shift, omega_l, N_cycles, rtUp, cep)  
    D_sqrt_T = sqrt(4._wp*A_ti**2*abs(dip10)**2 + (E0 - E1)**2) 
    sin_T = sin(D_sqrt_T * dt/2._wp) / D_sqrt_T
    cos_T = cos(D_sqrt_T * dt/2._wp)

    ! Now take the step regarding UT and 2LS 
    psi0_new = UT * (cos_T*psi0 + cmplx(0._wp,1._wp,wp)*sin_T*((E1-E0)*psi0 - 2._wp*A_ti*dip01*psi1)) * & 
               exp(cmplx(0._wp, -A_ti**2*dt, wp))

    psi1_new = UT * (cos_T*psi1 + cmplx(0._wp,1._wp,wp)*sin_T*((E0-E1)*psi1 - 2._wp*A_ti*dip10*psi0)) * &
               exp(cmplx(0._wp, -A_ti**2*dt, wp))

    ! To complete step add exp with A-field. Again splitting the exp to minimize expensive calls
    do j=1, Ny  ! Assuming field along y here 
        psi0_new(:,j) = psi0_new(:,j) * exp(cmplx(0._wp, -ky_list(j)*A_ti*dt,wp))
        psi1_new(:,j) = psi1_new(:,j) * exp(cmplx(0._wp, -ky_list(j)*A_ti*dt,wp))
    end do 

    ! Go back to real space 
    call fftw_execute_dft(plan0_b, psi0_new, psi0_new)
    call fftw_execute_dft(plan1_b, psi1_new, psi1_new)
    psi0_new = psi0_new / (Nx*Ny)  ! Take care of normalization (FFTW does not do this)
    psi1_new = psi1_new / (Nx*Ny)

    ! Last half step in real space again 
    psi0 = cos_fac*psi0_new + sin_fac*(V_diff*psi0_new - 2._wp*V01*psi1_new)
    psi1 = cos_fac*psi1_new + sin_fac*(-V_diff*psi1_new - 2._wp*V01*psi0_new)
    
    ! Apply mask to remove small parts going in direct y-direction
    do j=1, Nx 
        psi0(j,:) = psi0(j,:) * mask
        psi1(j,:) = psi1(j,:) * mask
    end do 

    ! Break loop if t > t_interract 
    !if (i*dt > 16.9_wp) exit
end do 

! Transform to momentum space to save momentum dist 
psi0_new = psi0 
psi1_new = psi1
call fftw_execute_dft(plan0_f, psi0_new, psi0)
call fftw_execute_dft(plan1_f, psi1_new, psi1)

! DESTROY THE FFT PLANS 
call fftw_destroy_plan(plan0_f)
call fftw_destroy_plan(plan1_f)
call fftw_destroy_plan(plan0_b)
call fftw_destroy_plan(plan1_b)

! Free up memory before saving 
deallocate(psi0_new, psi1_new, UT, sin_fac, cos_fac, V_diff, V01)


! SAVE DATA 
! First convert WFs to real momentum space 
call convert_FFT_to_momentum(psi0, x_list(1), y_list(1), dx, dy, dkx, dky, Nx, Ny)
call convert_FFT_to_momentum(psi1, x_list(1), y_list(1), dx, dy, dkx, dky, Nx, Ny)

allocate(psi0_save(Nx_save, Ny_save))
allocate(psi1_save(Nx_save, Ny_save))
psi0_save = psi0(kx_low_i:kx_high_i, ky_low_i:ky_high_i)
psi1_save = psi1(kx_low_i:kx_high_i, ky_low_i:ky_high_i)

! Save the data
open(file='data/psi0_k_r.dat', unit=1, form='unformatted')
open(file='data/psi0_k_i.dat', unit=2, form='unformatted')
open(file='data/psi1_k_r.dat', unit=3, form='unformatted')
open(file='data/psi1_k_i.dat', unit=4, form='unformatted')
write(1) real(psi0_save)
write(2) aimag(psi0_save) 
write(3) real(psi1_save)
write(4) aimag(psi1_save)

! Save file with dimensions of grid 
open(8, file='data/dim.txt')
write(8,*) kx_list(kx_low_i)
write(8,*) kx_list(kx_high_i) 
write(8,*) Nx_save
write(8,*) ky_list(ky_low_i)
write(8,*) ky_list(ky_high_i)
write(8,*) Ny_save
close(8)

contains 

! Sin2 vector potential 
function A_field(t, omega_l, N_cycles, rtUp, cep) result(res)
    real(wp) :: res 
    real(wp) :: t, omega_l, rtUp, cep
    integer :: N_cycles 
    
    if (t >= 0 .and. t <= 2*N_cycles*pi/omega_l) then 
        res = 2*rtUp * sin(omega_l * t / (2._wp*N_cycles))**2 * cos(omega_l*t + cep)
    else 
        res = 0._wp
    end if
end function A_field 

! A simple Gaussian WP in real space, waist at x0, propagated a time t, central momentum k0 and width in k-space sigma_k 
function Gauss_WP(x, t, x0, k0, sigma_k) result(res)
    real(wp) :: x, t, x0, k0, sigma_k
    complex(wp) :: res 
    complex(wp) :: gamma, N
    gamma = 1._wp / cmplx(1._wp, 2._wp*t*sigma_k**2, wp)
    N = sigma_k * sqrt(2._wp*gamma) / (2._wp*pi*sigma_k**2)**(0.25_wp)
    res = N * exp(gamma * (-sigma_k**2 * (x+x0)**2 + cmplx(0._wp, k0*(x+x0), wp)) + k0**2 / (4._wp*sigma_k**2)*(gamma - 1._wp))
end function Gauss_WP 

! Laser modulated WP in real space. Modulated in x-direction, Gaussian in y-direction.
function PINEM_WP(x, t, x0, k0, sigma_k, omega, max_N, g_factor) result(res)
    real(wp) :: x, t, x0, k0, sigma_k, omega, g_factor
    integer :: max_N 
    complex(wp) :: res 
    real(wp) :: phi_N, delta
    integer :: i 

    phi_N = 0._wp 
    delta = omega / k0 
    res = bessel_jn(0, g_factor) * exp(phi_N) * Gauss_WP(x,t,x0,k0,sigma_k)
    do i=1, max_N
        res = res + bessel_jn(i, g_factor) * exp(phi_N) * Gauss_WP(x,t,x0,k0+i*delta,sigma_k)
        res = res + bessel_jn(-i, g_factor) * exp(phi_N) * Gauss_WP(x,t,x0,k0-i*delta,sigma_k)
    end do
end function PINEM_WP

! Routine to load the scatter potential and build them into the matrices used in the calculation
subroutine load_scatt_pots(V00, V01, V11, E0, E1, dipole10, r_pot_max, Nx, Ny) 
    integer, intent(in) :: Nx, Ny 
    real(wp), intent(out) :: V00(Nx, Ny), V01(Nx, Ny), V11(Nx, Ny)
    real(wp), intent(out) :: E0, E1, r_pot_max
    complex(wp), intent(out) :: dipole10
    real(wp), allocatable :: V00_grid(:,:), V01_grid(:,:), V11_grid(:,:) 
    integer :: Nx_grid, Ny_grid, Nx_start, Ny_start, Nx_end, Ny_end
    real(wp) :: x_grid_low, x_grid_high, y_grid_low, y_grid_high, dip_real, dip_imag

    ! Load dimensions of the scatter potential grid 
    open(1, file='scatter_pots/dim.txt')
    read(1,*) Nx_grid, Ny_grid, Nx_start, Nx_end, Ny_start, Ny_end, x_grid_low, x_grid_high, y_grid_low, y_grid_high, r_pot_max
    close(1)

    ! Perform checks to see that scatter potential grid fits with the simulation grid
    if (abs(x_list(Nx_start) - x_grid_low) > 1e-3*dx) then 
        write(*,*) 'Mismatch in x grid start!'
        stop 
    else if (abs(x_list(Nx_end) - x_grid_high) > 1e-3*dx) then 
        write(*,*) 'Mismatch in x grid end!'
        stop 
    else if (abs(y_list(Ny_start) - y_grid_low) > 1e-3*dy) then 
        write(*,*) 'Mismatch in y grid start!'
        stop 
    else if (abs(y_list(Ny_end) - y_grid_high) > 1e-3*dy) then
        write(*,*) 'Mismatch in y grid end!'
        stop
    end if
    !write(*,*) x_list(Nx_start), x_grid_low, x_list(Nx_end), x_grid_high

    ! Load the scatter potentials 
    allocate(V00_grid(Nx_grid, Ny_grid), V01_grid(Nx_grid, Ny_grid), V11_grid(Nx_grid, Ny_grid))
    open(1, file='scatter_pots/V00_r.dat', form='unformatted')
    open(2, file='scatter_pots/V01_r.dat', form='unformatted')
    open(3, file='scatter_pots/V11_r.dat', form='unformatted')
    read(1) V00_grid
    read(2) V01_grid
    read(3) V11_grid
    close(1)
    close(2)
    close(3)

    ! Load energies and dipole element 
    open(1, file='scatter_pots/energies.txt')
    open(2, file='scatter_pots/dip10_element.txt')
    read(1,*) E0, E1
    read(2,*) dip_real, dip_imag
    close(1)
    close(2)
    dipole10 = cmplx(dip_real, dip_imag, wp)

    ! TODO: Should take care of potential shift in y0?... 

    ! Set the potential matrices used in calculation to 0, and build in the non-zero part 
    V00 = 0._wp
    V01 = 0._wp
    V11 = 0._wp
    V00(Nx_start:Nx_end, Ny_start:Ny_end) = V00_grid
    V01(Nx_start:Nx_end, Ny_start:Ny_end) = V01_grid
    V11(Nx_start:Nx_end, Ny_start:Ny_end) = V11_grid
end subroutine load_scatt_pots

! Routine to convert function on real grid to funciton ready for repeated FFT evaluation 
subroutine convert_space_to_FFT(psi, x_list, y_list, kx_left, ky_left, Nx, Ny)
    complex(wp), intent(inout) :: psi(Nx, Ny)
    real(wp), intent(in) :: x_list(Nx), y_list(Ny), kx_left, ky_left
    integer, intent(in) :: Nx, Ny 
    integer :: i
    do i=1, Nx
        psi(i,:) = psi(i,:) * exp(cmplx(0._wp, -kx_left*x_list(i), wp))
    end do 
    do j=1, Ny 
        psi(:,j) = psi(:,j) * exp(cmplx(0._wp, -ky_left*y_list(j), wp))
    end do
end subroutine convert_space_to_FFT

! Routine to convert FFT function in momentum space to real momentum function 
subroutine convert_FFT_to_momentum(psi, x_left, y_left, dx, dy, dkx, dky, Nx, Ny)
    complex(wp), intent(inout) :: psi(Nx, Ny)
    real(wp), intent(in) :: x_left, y_left, dx, dy, dkx, dky 
    integer, intent(in) :: Nx, Ny 
    integer :: i
    do i=1, Nx 
        psi(i,:) = psi(i,:) * exp(cmplx(0._wp, -x_left*dkx * (i-1), wp))
    end do 
    do i=1, Ny 
        psi(:,i) = psi(:,i) * exp(cmplx(0._wp, -y_left*dky * (i-1), wp))
    end do
    psi = psi * dx*dy/(2._wp*pi)
end subroutine convert_FFT_to_momentum

! A subroutine to calculate the frequencies used by FFT. For now only even number of points! (Currently not used)
subroutine fftfreq(arr, N, delta)
    integer, intent(in) :: N
    real(wp), intent(in) :: delta 
    real(wp), intent(out) :: arr(N)
    integer :: i 
    arr(1) = 0._wp 
    do i=1, N/2-1
        arr(i+1) = i 
    end do
    do i = N/2, N-1 
        arr(i+1) = -N/2 + (i-N/2)
    end do 
    arr = arr / (delta*N)
end subroutine fftfreq

end program split_step_scattering 
