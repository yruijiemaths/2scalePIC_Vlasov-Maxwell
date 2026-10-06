module poisson
use mesh, only: dx, dimx
use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f'

private
integer(8)                                   :: i, k, forward, backward
real(8)                                   :: pi
complex(8), dimension(:), allocatable :: rhs_t, psi_t

public init_poisson, solve_poisson, free_poisson, calcul_dx, solve_poisson_face, calcul_dx_face

contains


subroutine init_poisson(nx)

integer(8), intent(in) :: nx

pi = 4d0 * atan(1d0)

allocate(rhs_t(nx))
allocate(psi_t(nx))
call dfftw_plan_dft_1d(forward, nx, rhs_t, psi_t, fftw_forward, fftw_estimate)
call dfftw_plan_dft_1d(backward, nx, psi_t, rhs_t, fftw_backward, fftw_estimate)

end subroutine init_poisson


subroutine solve_poisson(rho, nx, dimx, ex)

integer(8), intent(in)    :: nx
real(8), intent(in)    :: rho(0:nx), dimx
real(8), intent(inout) :: ex(0:nx)
integer(8)                :: ind_x
real(8)                :: kx

pi = 4d0 * atan(1d0)
psi_t = cmplx(rho(0:nx-1), 0d0, kind=8)

call dfftw_execute_dft(forward, psi_t, rhs_t)
do i = 1, nx
    if (i <= nx/2) then
      ind_x = i - 1
    else
      ind_x = -nx + (i - 1)
    end if

    kx = 2d0 * pi * real(ind_x,8) / dimx

    if (ind_x == 0) then
      rhs_t(i) = cmplx(0d0,0d0,kind=8)     ! zero mode = 0
    else
      rhs_t(i) = -cmplx(0d0, kx,kind=8) * rhs_t(i) / cmplx(kx*kx, 0d0,kind=8)
    end if
end do
call dfftw_execute_dft(backward, rhs_t, psi_t)

ex(0:nx-1) = real(psi_t)
ex(nx) = ex(0)
ex = ex / real(nx,kind=8)

end subroutine solve_poisson


subroutine solve_poisson_face(rho, nx, dx, ex)

integer(8), intent(in)  :: nx
real(8),    intent(in)  :: dx
real(8),    intent(in)  :: rho(0:nx)
real(8),    intent(out) :: ex(0:nx)

integer(8) :: i
real(8)    :: rho_mean, ex_mean, accum
real(8), allocatable :: g(:)

allocate(g(0:nx-1))

rho_mean = sum(rho(0:nx-1)) / real(nx,8)

do i = 0, nx-1
    g(i) = rho(i) - rho_mean
end do

ex(0) = 0.0d0

do i = 1, nx-1
    ex(i) = ex(i-1) + dx * g(i)
end do

ex_mean = sum(ex(0:nx-1)) / real(nx,8)

do i = 0, nx-1
    ex(i) = ex(i) - ex_mean
end do

ex(nx) = ex(0)

deallocate(g)

end subroutine solve_poisson_face


subroutine calcul_dx(nx, dimx, f, dfdx)

integer(8), intent(in)    :: nx
real(8), intent(in)       :: dimx
real(8), intent(in)       :: f(0:nx)
real(8), intent(out)      :: dfdx(0:nx)
integer(8) :: i, ind_k
real(8) :: kx, pi
complex(8), allocatable :: fhat(:), dfhat(:)

pi = 4d0 * atan(1d0)

allocate(fhat(0:nx-1))
allocate(dfhat(0:nx-1))

do i = 0, nx-1
    fhat(i) = cmplx(f(i), 0d0, kind=8)
end do


call dfftw_execute_dft(forward, fhat, fhat)
do i = 0, nx-1
    if (i <= nx/2) then
        ind_k = i
    else
        ind_k = i - nx
    end if
    kx = 2d0 * pi * real(ind_k, kind=8) / dimx
    dfhat(i) = cmplx(0d0, kx, kind=8) * fhat(i)
end do
call dfftw_execute_dft(backward, dfhat, dfhat)
do i = 0, nx-1
    dfdx(i) = real(dfhat(i)) / real(nx, kind=8)
end do
dfdx(nx) = dfdx(0)  

deallocate(fhat)
deallocate(dfhat)

end subroutine calcul_dx


subroutine calcul_dx_face(nx, dimx, f, dfdx)

integer(8), intent(in)    :: nx
real(8),    intent(in)    :: dimx
real(8),    intent(in)    :: f(0:nx)
real(8),    intent(out)   :: dfdx(0:nx)

integer(8) :: i, im1
real(8)    :: dx


dx = dimx / real(nx, kind=8)

do i = 0, nx-1
    im1 = modulo(i-1, nx)
    dfdx(i) = ( f(i) - f(im1) ) / dx
end do

dfdx(nx) = dfdx(0)

end subroutine calcul_dx_face


subroutine free_poisson()

call dfftw_destroy_plan(forward)
call dfftw_destroy_plan(backward)
deallocate(psi_t)
deallocate(rhs_t)

end subroutine free_poisson

end module poisson
