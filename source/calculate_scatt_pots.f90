program calculate_scatt_pots 
use utils 
implicit none

integer, parameter :: wp = dp 

integer :: Nx, Ny 
real(wp) :: x_low, x_high, y_low, y_high, dx, dy
real(wp) :: dkx, dky, kx_max, ky_max, R, soft_param
real(wp), allocatable :: x_list(:), y_list(:), kx_list(:), ky_list(:)
real(wp), allocatable :: psi_r(:,:), psi_i(:,:)
complex(wp), allocatable :: psi0(:,:), psi1(:,:)
complex(wp), allocatable :: V00(:,:), V11(:,:), V01(:,:)  
complex(wp) :: dip10
integer :: i,j 

! Params for calculation grid 
integer :: Nx_grid, Ny_grid, Nx_cal_low, Nx_cal_high, Ny_cal_low, Ny_cal_high, Nx_cal, Ny_cal
real(wp) :: r_cal, x_low_grid, x_high_grid, y_low_grid, y_high_grid, dx_cal, dy_cal, r_boundary, ri, mask_i 
real(wp), allocatable :: x_list_grid(:), y_list_grid(:), x_list_cal(:), y_list_cal(:)


! Load the parameters form the imaginary propagation dim file 
open(1, file='data/dim_imag.txt')
read(1, *) Nx, Ny, x_low, x_high, y_low, y_high, dkx, dky, kx_max, ky_max, R, soft_param
close(1)

! Create spatial and momentum grids
allocate(x_list(Nx), y_list(Ny), kx_list(Nx), ky_list(Ny))
call linspace(x_list, x_low, x_high, Nx)
call linspace(y_list, y_low, y_high, Ny)
call linspace(kx_list, -kx_max, kx_max, Nx)
call linspace(ky_list, -ky_max, ky_max, Ny)

dx = x_list(2) - x_list(1)
dy = y_list(2) - y_list(1)


! LOAD MOMENTUM WF AND DETERMINE DIPOLE 
allocate(psi_r(Nx, Ny), psi_i(Nx, Ny), psi0(Nx, Ny), psi1(Nx, Ny)) 
open(1, file='data/imag_psi0k_r.dat', form='unformatted') 
open(2, file='data/imag_psi0k_i.dat', form='unformatted')
read(1) psi_r
read(2) psi_i
close(1)
close(2)
psi0 = cmplx(psi_r, psi_i, wp)

open(1, file='data/imag_psi1k_r.dat', form='unformatted')
open(2, file='data/imag_psi1k_i.dat', form='unformatted')
read(1) psi_r
read(2) psi_i
close(1)
close(2)
psi1 = cmplx(psi_r, psi_i, wp)

! Dipole element 
dip10 = 0.0_wp
do i=1, Ny
    dip10 = dip10 + sum(conjg(psi1(:,i)) * psi0(:,i) * ky_list(i)) 
end do
dip10 = dip10 * dkx * dky

write(*,*) 'Dipole element: ', dip10

open(1, file='scatter_pots/dip10_element.txt')
write(1,*) real(dip10)
write(1,*) aimag(dip10)
close(1)


! LOAD WF AND DETERMINE CALCULATION GRID 
open(1, file='data/imag_psi0_r.dat', form='unformatted') 
open(2, file='data/imag_psi0_i.dat', form='unformatted')
read(1) psi_r
read(2) psi_i
close(1)
close(2)
psi0 = cmplx(psi_r, psi_i, wp)

open(1, file='data/imag_psi1_r.dat', form='unformatted')
open(2, file='data/imag_psi1_i.dat', form='unformatted')
read(1) psi_r
read(2) psi_i
close(1)
close(2)
psi1 = cmplx(psi_r, psi_i, wp)

deallocate(psi_r, psi_i)

call load_grid_settings(Nx_grid, Ny_grid, x_low_grid, x_high_grid, y_low_grid, y_high_grid, r_cal, r_boundary)

allocate(x_list_grid(Nx_grid), y_list_grid(Ny_grid))
call linspace(x_list_grid, x_low_grid, x_high_grid, Nx_grid)
call linspace(y_list_grid, y_low_grid, y_high_grid, Ny_grid)

dx_cal = x_list_grid(2) - x_list_grid(1)
dy_cal = y_list_grid(2) - y_list_grid(1)

! Find the x-index where x is greater than +- r_cal
do i=1, Nx_grid 
    if (x_list_grid(i) >= -r_cal) then 
        Nx_cal_low = i 
        exit 
    end if
end do 
do i=Nx_cal_low, Nx_grid
    if (x_list_grid(i) > r_cal) then 
        Nx_cal_high = i-1
        exit 
    end if
end do
do i=1, Ny_grid
    if (y_list_grid(i) >= -r_cal) then 
        Ny_cal_low = i 
        exit 
    end if
end do
do i=Ny_cal_low, Ny_grid
    if (y_list_grid(i) > r_cal) then 
        Ny_cal_high = i-1
        exit 
    end if
end do
!write(*,*) x_list_grid(Nx_cal_low), x_list_grid(Nx_cal_high), y_list_grid(Ny_cal_low), y_list_grid(Ny_cal_high) 
Nx_cal = Nx_cal_high - Nx_cal_low + 1
Ny_cal = Ny_cal_high - Ny_cal_low + 1

write(*,*) 'Nx_cal : ', Nx_cal, ' Ny_cal : ', Ny_cal


! CALCULATE SCATTERING POTENTIALS 
allocate(V00(Nx_cal, Ny_cal), V11(Nx_cal, Ny_cal), V01(Nx_cal, Ny_cal)) 
allocate(x_list_cal(Nx_cal), y_list_cal(Ny_cal))
x_list_cal = x_list_grid(Nx_cal_low:Nx_cal_high)
y_list_cal = y_list_grid(Ny_cal_low:Ny_cal_high)

write(*,*) 'Caculating scattering potentials!'
call determine_scatt_potentials(V00, V11, V01, psi0, psi1, soft_param, R, &
                                 x_list_cal, y_list_cal, Nx_cal, Ny_cal)

! Let's apply smooth cutoff to the potentials 
do i=1, Nx_cal 
    do j=1, Ny_cal 
        ri = sqrt(x_list_cal(i)**2 + y_list_cal(j)**2)
        if (ri >= r_boundary) then 
            mask_i = 1._wp - cos(0.5_wp*pi * (r_cal - ri) / (r_cal-r_boundary))**8
            V00(i, j) = V00(i, j) * mask_i
            V11(i, j) = V11(i, j) * mask_i
            V01(i, j) = V01(i, j) * mask_i
        end if 
    end do
end do

! Now save 
write(*,*) 'Saving potentials...'
call save_arr('scatter_pots/V00_r.dat', real(V00), Nx_cal, Ny_cal)
call save_arr('scatter_pots/V00_i.dat', aimag(V00), Nx_cal, Ny_cal)
call save_arr('scatter_pots/V11_r.dat', real(V11), Nx_cal, Ny_cal)
call save_arr('scatter_pots/V11_i.dat', aimag(V11), Nx_cal, Ny_cal)
call save_arr('scatter_pots/V01_r.dat', real(V01), Nx_cal, Ny_cal)
call save_arr('scatter_pots/V01_i.dat', aimag(V01), Nx_cal, Ny_cal)

! Also save dim file with the cal indeices 
open(1, file='scatter_pots/dim.txt')
write(1,*) Nx_cal 
write(1,*) Ny_cal
write(1,*) Nx_cal_low 
write(1,*) Nx_cal_high
write(1,*) Ny_cal_low
write(1,*) Ny_cal_high
write(1,*) x_list_cal(1)
write(1,*) x_list_cal(Nx_cal)
write(1,*) y_list_cal(1)
write(1,*) y_list_cal(Ny_cal)
write(1,*) r_cal 
close(1)

contains 

! Stupid wrapper function to load the grid settings from the split-step settings file 
subroutine load_grid_settings(Nx_grid, Ny_grid, x_low_grid, x_high_grid, y_low_grid, y_high_grid, & 
                              r_cal, r_boundary)
    integer, intent(out) :: Nx_grid, Ny_grid
    real(wp), intent(out) :: x_low_grid, x_high_grid, y_low_grid, y_high_grid, r_cal, r_boundary
    integer :: Nx, Ny
    real(wp) :: x_low, x_high, y_low, y_high, r_max, dt, mask_bound, mask_max, cos_power

    namelist /SPLIT_STEP_GRID/ Nx, Ny, x_low, x_high, y_low, y_high, r_max, r_boundary, dt, mask_bound, mask_max, cos_power
    open(file='settings.nml', unit=1)
    read(nml=SPLIT_STEP_GRID, unit=1)
    close(1)

    Nx_grid = Nx 
    Ny_grid = Ny 
    x_low_grid = x_low 
    x_high_grid = x_high
    y_low_grid = y_low
    y_high_grid = y_high
    r_cal = r_max
end subroutine load_grid_settings

subroutine save_arr(save_path, arr, Nx, Ny)
    integer, intent(in) :: Nx, Ny
    character(len=*), intent(in) :: save_path
    real(wp), intent(in) :: arr(Nx,Ny)
    open(file=save_path, unit=1, form='unformatted')
    write(1) arr
    close(1)
end subroutine save_arr

subroutine determine_scatt_potentials(V00, V11, V01, psi0, psi1, soft_param, R, &
                                      x_list_cal, y_list_cal, Nx_cal, Ny_cal)
    integer, intent(in) :: Nx_cal, Ny_cal
    real(wp), intent(in) :: x_list_cal(Nx_cal), y_list_cal(Ny_cal)
    complex(wp), intent(out) :: V00(Nx_cal, Ny_cal), V11(Nx_cal, Ny_cal), V01(Nx_cal, Ny_cal)
    complex(wp), intent(in) :: psi0(Nx, Ny), psi1(Nx, Ny)
    real(wp), intent(in) :: soft_param, R
    integer :: i, j 
    complex(wp) :: V00_i, V11_i, V01_i

    !$OMP PARALLEL DO PRIVATE(j, i, V00_i, V11_i, V01_i)
    do j=1, Ny_cal
        write(*,*) 'Calculating row ', j, '/', Ny_cal
        do i=1, Nx_cal
            ! First test we are not outside r_cal radius 
            if (sqrt(x_list_cal(i)**2 + y_list_cal(j)**2) > r_cal) then 
                V00(i, j) = 0.0_wp
                V11(i, j) = 0.0_wp
                V01(i, j) = 0.0_wp
                cycle
            end if

            ! If not then actually calculate 
            call integrate_scatt_potential(V00_i, V11_i, V01_i, x_list_cal(i), y_list_cal(j), psi0, psi1, soft_param, R)
            V00(i, j) = V00_i
            V11(i, j) = V11_i
            V01(i, j) = V01_i
        end do 
    end do
    !$OMP END PARALLEL DO
end subroutine determine_scatt_potentials

subroutine integrate_scatt_potential(V00_i, V11_i, V01_i, x, y, psi0, psi1, soft_param, R)
    complex(wp), intent(out) :: V00_i, V11_i, V01_i
    real(wp), intent(in) :: x, y 
    complex(wp), intent(in) :: psi0(Nx, Ny), psi1(Nx, Ny)
    real(wp), intent(in) :: soft_param, R
    integer :: i, j
    real(wp) :: V_ee, V_a1, V_a2 
    complex(wp) :: overlap_int00, overlap_int11, overlap_int01, overlap00, overlap11, overlap01

    V00_i = 0.0_wp
    V11_i = 0.0_wp
    V01_i = 0.0_wp
    overlap_int00 = 0.0_wp
    overlap_int11 = 0.0_wp
    overlap_int01 = 0.0_wp

    ! First calculate the electron-electron interaction and determine the overlap 
    do j=1, Ny 
        do i=1, Nx
            ! Determine scattering interaction at each point 
            overlap00 = conjg(psi0(i, j)) * psi0(i, j)
            overlap11 = conjg(psi1(i, j)) * psi1(i, j)
            overlap01 = conjg(psi0(i, j)) * psi1(i, j)
            overlap_int00 = overlap_int00 + overlap00
            overlap_int11 = overlap_int11 + overlap11
            overlap_int01 = overlap_int01 + overlap01

            V_ee = 1.0_wp / sqrt((x_list(i) - x)**2 + (y_list(j) - y)**2 + soft_param)
            
            V00_i = V00_i + V_ee * overlap00 
            V11_i = V11_i + V_ee * overlap11
            V01_i = V01_i + V_ee * overlap01
        end do 
    end do
    
    ! Scale the numerical integration 
    V00_i = V00_i * dx*dy 
    V11_i = V11_i * dx*dy
    V01_i = V01_i * dx*dy
    overlap_int00 = overlap_int00 * dx*dy
    overlap_int11 = overlap_int11 * dx*dy
    overlap_int01 = overlap_int01 * dx*dy

    ! Now add the core integrals using the overlap integral
    V_a1 = 1._wp / sqrt(x**2 + y**2 + soft_param) 
    V_a2 = 0._wp 
    !V_a1 = 1._wp / sqrt(x**2 + (y-R/2._wp)**2 + soft_param) 
    !V_a2 = 1._wp / sqrt(x**2 + (y+R/2._wp)**2 + soft_param) 

    V00_i = V00_i - overlap_int00 * (V_a1 + V_a2)
    V11_i = V11_i - overlap_int11 * (V_a1 + V_a2)
    V01_i = V01_i - overlap_int01 * (V_a1 + V_a2) 
end subroutine integrate_scatt_potential

end program calculate_scatt_pots