module grid_output_adapter
  use mesh, only: field, nx, dx, dimx
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
      write(unit,*) f%exhalf(i), f%ey(i), f%bz(i)
    end do
  end subroutine write_field_snapshot

  subroutine compute_field_energy(f, value)
    type(field), intent(in) :: f
    real(8), intent(out) :: value
    value = 0.5d0*sum(f%exhalf(0:nx-1)**2 + f%ey(0:nx-1)**2 + &
                     f%bz(0:nx-1)**2)*dx
  end subroutine compute_field_energy

  subroutine compute_field_l2_norms(f, values)
    type(field), intent(in) :: f
    real(8), intent(out) :: values(3)
    values(1) = sqrt(sum(f%exhalf(0:nx-1)**2)*dx)
    values(2) = sqrt(sum(f%ey(0:nx-1)**2)*dx)
    values(3) = sqrt(sum(f%bz(0:nx-1)**2)*dx)
  end subroutine compute_field_l2_norms

  subroutine compute_poisson_residual(f, rho, value)
    type(field), intent(in) :: f
    real(8), intent(in) :: rho(0:nx)
    real(8), intent(out) :: value
    real(8) :: rho_mean
    integer :: i, im1
    value = 0.d0
    rho_mean = sum(rho(0:nx-1))*dx/dimx
    do i = 0, nx-1
      im1 = modulo(i-1, int(nx))
      value = value + ((f%exhalf(i)-f%exhalf(im1))/dx-rho(i)+rho_mean)**2
    end do
    value = sqrt(value*dx)
  end subroutine compute_poisson_residual
end module grid_output_adapter
