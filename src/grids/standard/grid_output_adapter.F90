module grid_output_adapter
  use mesh, only: field, nx, dx, dimx
  use poisson, only: calcul_dx
  implicit none
  private
  public :: write_field_snapshot, compute_field_energy, compute_field_l2_norms, &
            compute_poisson_residual
contains
  subroutine write_field_snapshot(f, unit)
    type(field), intent(in) :: f
    integer, intent(in) :: unit
    integer :: i
    do i = 0, nx
      write(unit,*) f%ex(i), f%ey(i), f%bz(i)
    end do
  end subroutine write_field_snapshot

  subroutine compute_field_energy(f, value)
    type(field), intent(in) :: f
    real(8), intent(out) :: value
    value = 0.5d0*sum(f%ex(0:nx-1)**2 + f%ey(0:nx-1)**2 + &
                     f%bz(0:nx-1)**2)*dx
  end subroutine compute_field_energy

  subroutine compute_field_l2_norms(f, values)
    type(field), intent(in) :: f
    real(8), intent(out) :: values(3)
    values(1) = sqrt(sum(f%ex(0:nx-1)**2)*dx)
    values(2) = sqrt(sum(f%ey(0:nx-1)**2)*dx)
    values(3) = sqrt(sum(f%bz(0:nx-1)**2)*dx)
  end subroutine compute_field_l2_norms

  subroutine compute_poisson_residual(f, rho, value)
    type(field), intent(in) :: f
    real(8), intent(in) :: rho(0:nx)
    real(8), intent(out) :: value
    real(8) :: derivative(0:nx), rho_mean
    call calcul_dx(nx, dimx, f%ex, derivative)
    rho_mean = sum(rho(0:nx-1))*dx/dimx
    value = sqrt(sum((derivative(0:nx-1)-rho(0:nx-1)+rho_mean)**2)*dx)
  end subroutine compute_poisson_residual
end module grid_output_adapter
