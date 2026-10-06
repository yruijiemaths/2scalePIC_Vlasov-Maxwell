!File: Module Particules
module particles

use mesh, only: field, deltat => dt,  &
                dimx, &
                nx, dx, pi!, eext

!use pusher, only: output_particles!, calpov

use interpolations, only: particle
use mpi
use pic_config, only: cfg_initial_condition, cfg_particles_per_cell, &
                      cfg_density_perturbation, cfg_two_stream_sigma_vx, &
                      cfg_two_stream_sigma_vy, cfg_weibel_sigma_vx, &
                      cfg_weibel_sigma_vy, cfg_stream_velocity, cfg_wave_number, &
                      cfg_bump_tail_amplitude, cfg_bump_drift_velocity, &
                      cfg_bump_thermal_velocity, cfg_bump_bulk_sigma_vx, &
                      cfg_bump_sigma_vy

implicit none

integer(8)    :: plocal
real(8)    :: ptotal

! Shared physical initial-condition selector.  Main programs use the same
! value for their field initialization.
integer, parameter :: initial_original_two_stream = 1
integer, parameter :: initial_anisotropic_weibel = 2
integer, parameter :: initial_bump_on_tail = 3
integer            :: initial_condition = initial_anisotropic_weibel

contains

subroutine initialize(p)

integer         :: code, psize, prank
type (particle) :: p
integer :: nbpm                                  ! particles / cell
integer(8)       :: kini
integer(4)       :: nn
real(8)    :: xi, yi, zi, temm, sigma1, sigma2, theta, alpha, kpara, vshift
integer, allocatable :: seed(:)


nbpm = cfg_particles_per_cell

call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)  
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code) 

initial_condition = cfg_initial_condition
if (initial_condition /= initial_original_two_stream .and. &
    initial_condition /= initial_anisotropic_weibel .and. &
    initial_condition /= initial_bump_on_tail) then
    if (prank == 0) print *, 'Invalid cfg_initial_condition; expected 1, 2 or 3.'
    call MPI_ABORT(MPI_COMM_WORLD, 1, code)
endif

ptotal = real(nbpm*nx,8)  ! number of particles on all processors
! Distribute the remainder instead of silently dropping particles when the
! global particle count is not divisible by the MPI process count.
plocal = nbpm*nx/int(psize,8)
if (int(prank,8) < modulo(nbpm*nx,int(psize,8))) plocal = plocal+1_8
p%nbpa = plocal

allocate(p%pos(plocal,1)); p%pos = 0.0_8
allocate(p%vit(plocal,2)); p%vit = 0.0_8
allocate(p%w(plocal))  ; p%w = dimx / ptotal! 
allocate(p%ele_x(plocal)); p%ele_x = 0.0_8
allocate(p%ele_y(plocal)); p%ele_y = 0.0_8
allocate(p%mag_z(plocal)); p%mag_z = 0.0_8

nn=1
call random_seed(size = nn)
allocate(seed(nn))
seed = 1
seed = seed + prank
call random_seed(put=seed)
!print*, 'Prank =', prank, ', Final seed =', seed, 'Particles =', plocal


kpara = cfg_wave_number
select case (initial_condition)
case (initial_original_two_stream)
    alpha = 0.d0
    sigma1 = cfg_two_stream_sigma_vx
    sigma2 = cfg_two_stream_sigma_vy
    vshift = cfg_stream_velocity
    if (prank == 0) print *, 'Initial condition 1: original two-stream'
case (initial_anisotropic_weibel)
    alpha = cfg_density_perturbation
    sigma1 = cfg_weibel_sigma_vx
    sigma2 = cfg_weibel_sigma_vy
    vshift = 0.d0
    if (prank == 0) write(*,'(a,3es12.4)') &
      'Initial condition 2: anisotropic Maxwellian sigma_vx/sigma_vy/alpha = ', &
      sigma1, sigma2, alpha
case (initial_bump_on_tail)
    alpha = cfg_density_perturbation
    sigma1 = cfg_bump_bulk_sigma_vx
    sigma2 = cfg_bump_sigma_vy
    vshift = cfg_bump_drift_velocity
    if (prank == 0) print *, 'Initial condition 3: bump-on-tail'
end select

kini=1
do while (kini<=plocal)
    call random_number(xi)
    xi=dimx*xi
    call random_number(zi)
    zi=(1+alpha)*zi
    temm=(1.0d0 + alpha*dcos(kpara*xi))
    if (temm>=zi) then
    p%pos(kini,1)=xi
    kini=kini+1
    endif
enddo

select case (initial_condition)
case (initial_original_two_stream)
    kini=1
    do while (kini<=plocal)
        call random_number(xi)
        xi=(xi-0.5d0)*12*sigma1
        call random_number(yi)
        yi=(yi-0.5d0)*12*sigma2
        call random_number(zi)
        zi=zi/(2*pi*sigma1*sigma2)
        temm = dexp(-0.5d0*(yi/sigma2)**2) * &
               (dexp(-0.5d0*((xi-vshift)/sigma1)**2) + &
                dexp(-0.5d0*((xi+vshift)/sigma1)**2)) / &
               (2.d0*pi*sigma1*sigma2)
        if (temm>=zi) then
            p%vit(kini,1)=xi
            p%vit(kini,2)=yi
            kini=kini+1
        endif
    enddo
case (initial_anisotropic_weibel)
    kini=1
    do while (kini<=plocal)
        call random_number(xi)
        call random_number(yi)
        if (xi <= 1.d-14) cycle
        theta = 2.d0*pi*yi
        p%vit(kini,1) = sigma1*sqrt(-2.d0*log(xi))*cos(theta)
        p%vit(kini,2) = sigma2*sqrt(-2.d0*log(xi))*sin(theta)
        kini=kini+1
    enddo
case (initial_bump_on_tail)
    kini=1
    do while (kini<=plocal)
        call random_number(zi)
        call random_number(xi)
        call random_number(yi)
        if (xi <= 1.d-14) cycle
        theta = 2.d0*pi*yi
        if (zi < cfg_bump_tail_amplitude*cfg_bump_thermal_velocity / &
                 (sigma1 + cfg_bump_tail_amplitude*cfg_bump_thermal_velocity)) then
            p%vit(kini,1) = vshift + cfg_bump_thermal_velocity * &
                             sqrt(-2.d0*log(xi))*cos(theta)
        else
            p%vit(kini,1) = sigma1*sqrt(-2.d0*log(xi))*cos(theta)
        endif
        call random_number(xi)
        call random_number(yi)
        if (xi <= 1.d-14) cycle
        theta = 2.d0*pi*yi
        p%vit(kini,2) = sigma2*sqrt(-2.d0*log(xi))*sin(theta)
        kini=kini+1
    enddo
end select

end subroutine initialize

end module particles
