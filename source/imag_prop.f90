program imag_prop
use, intrinsic :: iso_c_binding
use utils  
use omp_lib
implicit none 
include "fftw3.f03" 

integer, parameter :: wp = dp

real(wp) :: x_low, x_high, y_low, y_high, dt, t_final, soft_param
real(wp) :: dx, dy, dkx, dky, kx_max, ky_max 
real(wp) :: dt_imag
integer :: Nx, Ny, Nt, Nt_imag, Nt_check, N_threads 
complex(wp), allocatable :: psi_guess(:,:), psi0(:,:), psi1(:,:), psi_project(:,:,:)
real(wp), allocatable :: UV(:,:), UT(:,:), E_T(:,:), E_V(:,:)  ! Lists to hold the energy, needed to check convergence
real(wp), allocatable :: x_list(:), y_list(:), kx_list(:), ky_list(:)
integer :: i,j, void 
complex(wp) :: E0, E1
real(wp) :: convergence_0, convergence_1, R 
type(c_ptr) :: plan_forward


! LOAD VARIABLES
namelist /IMAG_SETTINGS/ Nx, Ny, x_low, x_high, y_low, y_high, t_final, dt, Nt, soft_param, N_threads 
open(file='settings_imag.nml', unit=1)
read(nml=IMAG_SETTINGS, unit=1)
close(1)

dt_imag = dt 
Nt_imag = Nt
Nt_check = int(Nt_imag / 10)


! PREPARE IMAGINARY TIME PROPAGATION
! Setup multithread FFT 
call omp_set_num_threads(N_threads)
void = fftw_init_threads()
call fftw_plan_with_nthreads(omp_get_max_threads())

! Define the coordinate/momentum lists 
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

! Build the potential and kinetic operators 
R = 3._wp ! This is equillibrium : 2.424_wp 
allocate(UT(Nx,Ny), UV(Nx,Ny), E_T(Nx,Ny), E_V(Nx,Ny))
do j=1, Ny
    do i=1, Nx 
        UT(i,j) = exp(-(kx_list(i)**2 + ky_list(j)**2)/2._wp*dt_imag)
        E_T(i,j) = (kx_list(i)**2 + ky_list(j)**2)/(2._wp)

        UV(i,j) = exp(-soft_Coulomb(x_list(i), y_list(j), soft_param) * dt_imag * 0.5_wp) 
        E_V(i,j) = soft_Coulomb(x_list(i), y_list(j), soft_param)
        !UV(i,j) = exp(-H2plus(x_list(i), y_list(j), 0._wp, R, soft_param) * dt_imag * 0.5_wp) 
        !E_V(i,j) = H2plus(x_list(i), y_list(j), 0._wp, R, soft_param)
    end do 
end do 

! Initialize psi guess - just take a Guassian for the ground state? 
allocate(psi_guess(Nx, Ny)) 
do j=1, Ny 
    do i=1, Nx 
        !psi_guess(i,j) = exp(-((x_list(i)+ 0.1_wp)**2  + (y_list(j) - 0.05_wp)**2)/4._wp)
        psi_guess(i,j) = exp(-((x_list(i))**2  + (y_list(j))**2)/4._wp) 
    end do 
end do
psi_guess = psi_guess / sqrt(trapz_2D(abs(psi_guess)**2, dx, dy, Nx, Ny))

write(*,*) 'Initial norm : ', trapz_2D(abs(psi_guess)**2, dx, dy, Nx, Ny)


! IMAGINARY TIME PROPAGATION 
allocate(psi0(Nx, Ny))
plan_forward = fftw_plan_dft_2d(Ny, Nx, psi0, psi0, FFTW_FORWARD, FFTW_MEASURE)  ! Make this here, before we put anything in psi0

! Obtain the ground state 
call imag_prop_groundstate(psi_guess, psi0, E0, convergence_0)

write(*,*)
write(*,*) 'Starting calculation for the excited state...'
! Prepare guess for excited state 
do j=1, Ny 
    do i=1, Nx 
        psi_guess(i,j) = exp(-((x_list(i))**2 + (y_list(j)-0.3_wp)**2)) 
    end do 
end do

allocate(psi1(Nx, Ny), psi_project(Nx, Ny, 1)) 
psi_project(:,:,1) = psi0
call convert_space_to_FFT(psi_project(:,:,1), x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)

! Obtain the excited state by projecting out ground state 
call imag_prop_excited(psi_guess, psi1, E1, convergence_1, psi_project, 1)

write(*,*)

! SAVE RESULTS 
write(*,*) 'Saving final states to files...'
call save_arr('data/imag_psi0_r.dat', real(psi0), Nx, Ny)
call save_arr('data/imag_psi0_i.dat', aimag(psi0), Nx, Ny)
call save_arr('data/imag_psi1_r.dat', real(psi1), Nx, Ny)
call save_arr('data/imag_psi1_i.dat', aimag(psi1), Nx, Ny)

! We actually also want to save the momentum space states?...
write(*,*) 'Saving momentum space states to files...'
call convert_space_to_FFT(psi0, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
call fftw_execute_dft(plan_forward, psi0, psi0)
call convert_FFT_to_momentum(psi0, x_low, y_low, dx, dy, dkx, dky, Nx, Ny)
call save_arr('data/imag_psi0k_r.dat', real(psi0), Nx, Ny)
call save_arr('data/imag_psi0k_i.dat', aimag(psi0), Nx, Ny)

! Do the same for the exicted state 
psi0 = psi1
call convert_space_to_FFT(psi0, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
call fftw_execute_dft(plan_forward, psi0, psi0)
call convert_FFT_to_momentum(psi0, x_low, y_low, dx, dy, dkx, dky, Nx, Ny)
call save_arr('data/imag_psi1k_r.dat', real(psi0), Nx, Ny)
call save_arr('data/imag_psi1k_i.dat', aimag(psi0), Nx, Ny)

! Save energies and convergence! 
write(*,*) 'Saving energies and convergence to file...'
open(7, file='data/energies.txt')
write(7,*) real(E0) 
write(7,*) aimag(E0)
write(7,*) convergence_0
write(7,*) real(E1)
write(7,*) aimag(E1)
write(7,*) convergence_1
close(7)

! Also save energies to scatter_pots folder
open(7, file='scatter_pots/energies.txt')
write(7,*) real(E0)
write(7,*) real(E1) 
close(7)

! Save grid dimensions 
write(*,*) 'Saving grid dimensions to file...'
open(8, file='data/dim_imag.txt')
write(8,*) Nx 
write(8,*) Ny 
write(8,*) x_low
write(8,*) x_high
write(8,*) y_low 
write(8,*) y_high
write(8,*) dkx
write(8,*) dky
write(8,*) kx_max
write(8,*) ky_max
write(8,*) R 
write(8,*) soft_param
close(8)

contains 

function soft_Coulomb(x, y, a) result(V)
    real(wp) :: x, y, a, V
    real(wp) :: r
    
    r = x**2 + y**2
    V = -1._wp / sqrt(r + a)
end function soft_Coulomb

function H2plus(x, y, Rx, Ry, a) result(V)
    real(wp) :: x, y, Rx, Ry, a, V
    real(wp) :: r1, r2
    
    r1 = sqrt((x - Rx/2)**2 + (y - Ry/2)**2 + a**2)
    r2 = sqrt((x + Rx/2)**2 + (y + Ry/2)**2 + a**2)
    V = -1._wp / r1 - 1._wp / r2
end function H2plus

subroutine save_arr(save_path, arr, Nx, Ny)
    integer, intent(in) :: Nx, Ny
    character(len=*), intent(in) :: save_path
    real(wp), intent(in) :: arr(Nx,Ny)
    open(file=save_path, unit=1, form='unformatted')
    write(1) arr
    close(1)
end subroutine save_arr

subroutine imag_prop_excited(psi_in, psi_out, E, conv, psi_project, N_project)
    integer, intent(in) :: N_project
    complex(wp), intent(in) :: psi_in(Nx, Ny)
    complex(wp), intent(out) :: psi_out(Nx, Ny)
    complex(wp), intent(in) :: psi_project(Nx, Ny, N_project)
    complex(wp), intent(out) :: E 
    real(wp), intent(out) :: conv
    complex(wp) :: psi_H(Nx, Ny)
    type(c_ptr) :: plan_forward, plan_backward, plan_forward_H, plan_backward_H 
    integer :: i, j

    ! Make the FFTW plans (Nx and Ny must be swapped in plans compared to arrays)
    ! NB: psi/psi_k is destroyed by plans, so must initialize state AFTER plans are made 
    plan_forward = fftw_plan_dft_2d(Ny, Nx, psi_out, psi_out, FFTW_FORWARD, FFTW_MEASURE)  ! FFTW_MEASURE to get real speed
    plan_backward = fftw_plan_dft_2d(Ny, Nx, psi_out, psi_out, FFTW_BACKWARD, FFTW_MEASURE)
    plan_forward_H = fftw_plan_dft_2d(Ny, Nx, psi_H, psi_H, FFTW_FORWARD, FFTW_MEASURE)  
    plan_backward_H = fftw_plan_dft_2d(Ny, Nx, psi_H, psi_H, FFTW_BACKWARD, FFTW_MEASURE)
    
    psi_out = psi_in

    !Send the initial state to FFT ready form 
    call convert_space_to_FFT(psi_out, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
    
    !PROPAGATION (split the steps!)
    do i=1, Nt_imag   
        ! First half step in real space 
        psi_out = psi_out * UV 
        
        ! Then we go to momentum space and take full step 
        call fftw_execute_dft(plan_forward, psi_out, psi_out)  
        psi_out = psi_out * UT 
    
        ! Then go back to real space again 
        call fftw_execute_dft(plan_backward, psi_out, psi_out)
        psi_out = psi_out / (Nx*Ny) * UV

        ! Normalize after each time step! 
        psi_out = psi_out / sqrt(trapz_2D(abs(psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny))
    
        ! Project out states 
        psi_H = 0._wp  ! Use psi_H as temp storage
        do j=1, N_project 
            psi_H = psi_H + psi_project(:,:,j) * trapz_2D(conjg(psi_project(:,:,j)) * psi_out, x_list(2)-x_list(1), & 
                                                          y_list(2)-y_list(1), Nx, Ny)
        end do 
        !write(*,*) 'SIZE OF psi_H = ', trapz_2D(conjg(psi_project(:,:,1)) * psi_out, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
        psi_out = psi_out - psi_H
        psi_out = psi_out / sqrt(trapz_2D(abs(psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny))
    
        ! Check convergence and energy ~ 100 times during the propagation
        if (mod(i, Nt_check) == 0) then 
            ! Calculate H on psi
            psi_H = psi_out 
            call fftw_execute_dft(plan_forward_H, psi_H, psi_H) 
            psi_H = psi_H * E_T 
            call fftw_execute_dft(plan_backward_H, psi_H, psi_H)
            psi_H = psi_H / (Nx*Ny)  ! Fix FFT normliazation 
            psi_H = psi_H + psi_out * E_V  ! Add the potential energy 
    
            ! Now get E and check convergence 
            E = trapz_2D(conjg(psi_out) * psi_H, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
            conv = trapz_2D(abs(psi_H - E*psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
            write(*,*) i
            write(*,*) '    Energy = ', E
            write(*,*) '    Conv. = ', conv
        end if 
    end do 
    
    write(*,*) 'Propagation done!'
    write(*,*) 'Final norm : ', trapz_2D(abs(psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    
    ! Delete the plans
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
    
    ! Final calculation of energy and convergence! 
    psi_H = psi_out 
    call fftw_execute_dft(plan_forward_H, psi_H, psi_H) 
    psi_H = psi_H * E_T 
    call fftw_execute_dft(plan_backward_H, psi_H, psi_H)
    psi_H = psi_H / (Nx*Ny)  ! Fix FFT normliazation 
    psi_H = psi_H + psi_out * E_V  ! Add the potential energy 
    E = trapz_2D(conjg(psi_out) * psi_H, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    conv = trapz_2D(abs(psi_H - E*psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    
    call convert_FFT_to_space(psi_out, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
    
    ! Delete the plans
    call fftw_destroy_plan(plan_forward_H)
    call fftw_destroy_plan(plan_backward_H)
end subroutine imag_prop_excited

subroutine imag_prop_groundstate(psi_in, psi_out, E, conv) 
    complex(wp), intent(in) :: psi_in(Nx, Ny)
    complex(wp), intent(out) :: psi_out(Nx, Ny)
    complex(wp), intent(out) :: E 
    real(wp), intent(out) :: conv
    complex(wp) :: psi_H(Nx, Ny)
    type(c_ptr) :: plan_forward, plan_backward, plan_forward_H, plan_backward_H 
    integer :: i 

    ! Make the FFTW plans (Nx and Ny must be swapped in plans compared to arrays)
    ! NB: psi/psi_k is destroyed by plans, so must initialize state AFTER plans are made 
    plan_forward = fftw_plan_dft_2d(Ny, Nx, psi_out, psi_out, FFTW_FORWARD, FFTW_MEASURE)  ! FFTW_MEASURE to get real speed
    plan_backward = fftw_plan_dft_2d(Ny, Nx, psi_out, psi_out, FFTW_BACKWARD, FFTW_MEASURE)
    plan_forward_H = fftw_plan_dft_2d(Ny, Nx, psi_H, psi_H, FFTW_FORWARD, FFTW_MEASURE)  
    plan_backward_H = fftw_plan_dft_2d(Ny, Nx, psi_H, psi_H, FFTW_BACKWARD, FFTW_MEASURE)
    
    psi_out = psi_in

    !Send the initial state to FFT ready form 
    call convert_space_to_FFT(psi_out, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
    

    !PROPAGATION (split the steps!)
    do i=1, Nt_imag   
        ! First half step in real space 
        psi_out = psi_out * UV 
        
        ! Then we go to momentum space and take full step 
        call fftw_execute_dft(plan_forward, psi_out, psi_out)  
        psi_out = psi_out * UT 
    
        ! Then go back to real space again 
        call fftw_execute_dft(plan_backward, psi_out, psi_out)
        psi_out = psi_out / (Nx*Ny) * UV
    
        !write(*,*) trapz_2D(abs(psi)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    
        ! Normalize after each time step! 
        psi_out = psi_out / sqrt(trapz_2D(abs(psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny))
    
        !write(*,*) trapz_2D(abs(psi)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
        ! Check convergence and energy ~ 100 times during the propagation
        if (mod(i, Nt_check) == 0) then 
            ! Calculate H on psi
            psi_H = psi_out 
            call fftw_execute_dft(plan_forward_H, psi_H, psi_H) 
            psi_H = psi_H * E_T 
            call fftw_execute_dft(plan_backward_H, psi_H, psi_H)
            psi_H = psi_H / (Nx*Ny)  ! Fix FFT normliazation 
            psi_H = psi_H + psi_out * E_V  ! Add the potential energy 
    
            ! Now get E and check convergence 
            E = trapz_2D(conjg(psi_out) * psi_H, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
            conv = trapz_2D(abs(psi_H - E*psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
            write(*,*) i
            write(*,*) '    Energy = ', E
            write(*,*) '    Conv. = ', conv
        end if 
    end do 
    
    write(*,*) 'Propagation done!'
    write(*,*) 'Final norm : ', trapz_2D(abs(psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    
    ! Delete the plans
    call fftw_destroy_plan(plan_forward)
    call fftw_destroy_plan(plan_backward)
    
    ! Final calculation of energy and convergence! 
    psi_H = psi_out 
    call fftw_execute_dft(plan_forward_H, psi_H, psi_H) 
    psi_H = psi_H * E_T 
    call fftw_execute_dft(plan_backward_H, psi_H, psi_H)
    psi_H = psi_H / (Nx*Ny)  ! Fix FFT normliazation 
    psi_H = psi_H + psi_out * E_V  ! Add the potential energy 
    E = trapz_2D(conjg(psi_out) * psi_H, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    conv = trapz_2D(abs(psi_H - E*psi_out)**2, x_list(2)-x_list(1), y_list(2)-y_list(1), Nx, Ny)
    
    call convert_FFT_to_space(psi_out, x_list, y_list, kx_list(1), ky_list(1), Nx, Ny)
    
    ! Delete the plans
    call fftw_destroy_plan(plan_forward_H)
    call fftw_destroy_plan(plan_backward_H)
end subroutine imag_prop_groundstate

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

! Inverse of the function above 
subroutine convert_FFT_to_space(psi, x_list, y_list, kx_left, ky_left, Nx, Ny)
    complex(wp), intent(inout) :: psi(Nx, Ny)
    real(wp), intent(in) :: x_list(Nx), y_list(Ny), kx_left, ky_left
    integer, intent(in) :: Nx, Ny 
    integer :: i
    do i=1, Nx
        psi(i,:) = psi(i,:) * exp(cmplx(0._wp, kx_left*x_list(i), wp))
    end do 
    do j=1, Ny 
        psi(:,j) = psi(:,j) * exp(cmplx(0._wp, ky_left*y_list(j), wp))
    end do
end subroutine convert_FFT_to_space

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


end program imag_prop 