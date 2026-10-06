module poisson
use mesh, only: dx, dimx
use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f'

private
integer(8)                            :: i, k, forward, backward
real(8)                               :: pi
complex(8), dimension(:), allocatable :: rhs_t, psi_t

public init_poisson, solve_poisson, free_poisson, calcul_dx, calcul_inv_dx, remove_nyquist_from_rho

contains


subroutine init_poisson(nx)

integer(8), intent(in) :: nx

pi = 4d0 * atan(1d0)

allocate(rhs_t(0:nx-1))
allocate(psi_t(0:nx-1))

call dfftw_plan_dft_1d(forward, nx, rhs_t, rhs_t, fftw_forward, fftw_estimate)
call dfftw_plan_dft_1d(backward, nx, psi_t, psi_t, fftw_backward, fftw_estimate)

end subroutine init_poisson


subroutine solve_poisson(rho, nx, dimx, ex)
implicit none
integer(8), intent(in) :: nx
real(8), intent(in) :: dimx
real(8), intent(in) :: rho(0:nx)
real(8), intent(inout) :: ex(0:nx)
integer(8) :: i, ind_x
real(8) :: kx, pi
complex(8), allocatable :: rho_t(:), ex_t(:)

pi = 4d0 * atan(1d0)

allocate(rho_t(0:nx-1))
allocate(ex_t(0:nx-1))

rho_t(0:nx-1) = cmplx(rho(0:nx-1), 0d0, kind=8)

call dfftw_execute_dft(forward, rho_t, rho_t)

do i = 0, nx-1

    if (i <= nx/2) then
        ind_x = i
    else
        ind_x = i - nx
    end if

    kx = 2d0 * pi * real(ind_x,8) / dimx

    if (ind_x == 0) then
        ex_t(i) = cmplx(0d0, 0d0, kind=8)
    else if (mod(nx,2_8) == 0 .and. i == nx/2) then
        ex_t(i) = cmplx(0d0, 0d0, kind=8)
    else
        ex_t(i) = rho_t(i) / cmplx(0d0, kx, kind=8)
    end if

end do

call dfftw_execute_dft(backward, ex_t, ex_t)

do i = 0, nx-1
    ex(i) = real(ex_t(i)) / real(nx,8)
end do

ex(nx) = ex(0)

deallocate(rho_t)
deallocate(ex_t)

end subroutine solve_poisson


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


subroutine calcul_inv_dx(nx, dimx, dfdx, mean_f, f)
implicit none
integer(8), intent(in)    :: nx
real(8), intent(in)       :: dimx
real(8), intent(in)       :: dfdx(0:nx)
real(8), intent(in)       :: mean_f
real(8), intent(out)      :: f(0:nx)

integer(8) :: i, ind_k
real(8) :: kx, pi, mean_dfdx
complex(8), allocatable :: dfhat(:), fhat(:)

pi = 4d0 * atan(1d0)

allocate(dfhat(0:nx-1))
allocate(fhat(0:nx-1))

mean_dfdx = sum(dfdx(0:nx-1)) / real(nx,8)

do i = 0, nx-1
    dfhat(i) = cmplx(dfdx(i) - mean_dfdx, 0d0, kind=8)
end do

call dfftw_execute_dft(forward, dfhat, dfhat)

do i = 0, nx-1
    if (i == 0) then
        fhat(i) = cmplx(real(nx,8) * mean_f, 0d0, kind=8)

    else if (mod(nx,2_8) == 0 .and. i == nx/2) then
        ! Nyquist 模不能稳定地除以 ik，直接置零
        fhat(i) = cmplx(0d0, 0d0, kind=8)

    else
        if (i <= nx/2) then
            ind_k = i
        else
            ind_k = i - nx
        end if

        kx = 2d0 * pi * real(ind_k, kind=8) / dimx
        fhat(i) = dfhat(i) / cmplx(0d0, kx, kind=8)
    end if
end do

call dfftw_execute_dft(backward, fhat, fhat)

do i = 0, nx-1
    f(i) = real(fhat(i)) / real(nx,8)
end do

f(nx) = f(0)

deallocate(dfhat)
deallocate(fhat)

end subroutine calcul_inv_dx


subroutine remove_nyquist_from_rho(rho, nx)

implicit none

integer(8), intent(in) :: nx
real(8), intent(inout) :: rho(0:nx)

integer(8) :: i
complex(8), allocatable :: rho_hat(:)

allocate(rho_hat(0:nx-1))

do i = 0, nx-1
    rho_hat(i) = cmplx(rho(i), 0.d0, kind=8)
enddo

call dfftw_execute_dft(forward, rho_hat, rho_hat)

if (mod(nx,2_8) == 0) then
    rho_hat(nx/2) = cmplx(0.d0, 0.d0, kind=8)
endif
call dfftw_execute_dft(backward, rho_hat, rho_hat)

do i = 0, nx-1
    rho(i) = real(rho_hat(i)) / real(nx,8)
enddo

rho(nx) = rho(0)

deallocate(rho_hat)

end subroutine remove_nyquist_from_rho


subroutine free_poisson()

call dfftw_destroy_plan(forward)
call dfftw_destroy_plan(backward)
deallocate(psi_t)
deallocate(rhs_t)

end subroutine free_poisson

end module poisson
