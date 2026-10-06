module interpolations

use mesh, only: field, dx, nx, dimx, dvx, dvy, nvx, nvy

implicit none

integer(8), private :: ipart

type particle
   integer(8)           :: nbpa
   real(8), allocatable :: pos(:,:), vit(:,:)
   real(8), allocatable :: ele_x(:), ele_y(:), mag_z(:)
   real(8), allocatable :: w(:)
end type particle

contains


function f_m6( q )
real(8), intent(in) :: q
real(8)             :: f_m6

if ( q < 1.0d0 ) then
f_m6 = (3.0d0-q)**5-6.0d0*(2.0d0-q)**5+15.0d0*(1.0d0-q)**5
else if ( q >= 1.0d0 .and. q < 2.0d0 ) then
f_m6 = (3.0d0-q)**5-6.0d0*(2.0d0-q)**5
else if ( q >= 2.0d0 .and. q < 3.0d0 ) then
f_m6 = (3.0d0-q)**5
else
f_m6 = 0.0d0
end if

f_m6 = f_m6 / 120.0d0

return
end function f_m6


function df_m6(q)
real(8), intent(in) :: q
real(8)             :: df_m6

if ( q < 1.0d0 ) then
   df_m6 = -5.d0*(3.d0-q)**4+30.d0*(2.d0-q)**4-75.d0*(1.d0-q)**4
else if ( q < 2.0d0 ) then
   df_m6 = -5.d0*(3.d0-q)**4+30.d0*(2.d0-q)**4
else if ( q < 3.0d0 ) then
   df_m6 = -5.d0*(3.d0-q)**4
else
   df_m6 = 0.d0
end if

df_m6 = df_m6 / 120.d0
return
end function df_m6


subroutine interpol_m6(f, p)

type(particle), intent(inout)      ::p
type(field),    intent(in)    :: f
integer(8)            :: k
integer(8)            :: i
integer(8)            :: im1, im2, im3, ip1, ip2, ip3
real(8)           :: dpx,xp
real(8)          :: cm3x, cm2x, cm1x, cx, cp1x, cp2x, cp3x
real(8)           :: e(-2:2)

real(8), allocatable :: ex(:), ey(:), mz(:)

allocate(ex(0:nx-1))
allocate(ey(0:nx-1))
allocate(mz(0:nx-1))

e(-2) =   1.0_8/120.0_8
e(-1) =  26.0_8/120.0_8
e( 0) =  66.0_8/120.0_8
e(+1) =  26.0_8/120.0_8
e(+2) =   1.0_8/120.0_8

do i = 0,nx-1

im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)

ex(i) = e(-2)*f%ex(im2) + e(-1)*f%ex(im1) + e(0)*f%ex(i) + e(1)*f%ex(ip1) + e(2)*f%ex(ip2)
ey(i) = e(-2)*f%ey(im2) + e(-1)*f%ey(im1) + e(0)*f%ey(i) + e(1)*f%ey(ip1) + e(2)*f%ey(ip2)
mz(i) = e(-2)*f%bz(im2) + e(-1)*f%bz(im1) + e(0)*f%bz(i) + e(1)*f%bz(ip1) + e(2)*f%bz(ip2)

end do


do k=1,p%nbpa

xp = p%pos(k,1)/dx
i = floor(xp)
dpx = xp - real(i, kind=8)
i=modulo(i,nx)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.0d0+dpx)
cp3x = f_m6(3.0d0-dpx)
cm2x = f_m6(2.0d0+dpx)
cp2x = f_m6(2.0d0-dpx)
cm1x = f_m6(1.0d0+dpx)
cp1x = f_m6(1.0d0-dpx)
cx   = f_m6(dpx)

p%ele_x(k) = cm3x*ex(im3) + cm2x*ex(im2) + cm1x*ex(im1) + cx*ex(i) + cp1x*ex(ip1) + cp2x*ex(ip2) + cp3x*ex(ip3)
p%ele_y(k) = cm3x*ey(im3) + cm2x*ey(im2) + cm1x*ey(im1) + cx*ey(i) + cp1x*ey(ip1) + cp2x*ey(ip2) + cp3x*ey(ip3)
p%mag_z(k) = cm3x*mz(im3) + cm2x*mz(im2) + cm1x*mz(im1) + cx*mz(i) + cp1x*mz(ip1) + cp2x*mz(ip2) + cp3x*mz(ip3)

end do

deallocate(ex)
deallocate(ey)
deallocate(mz)

end subroutine interpol_m6


subroutine calcul_rho_m6( p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8)  :: i, k
integer(8)  :: im1, im2, im3, ip1, ip2, ip3
real(8)  :: dpx
real(8) :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)  :: weight, xp

f%rho = 0.0_8

do k = 1, p%nbpa

xp = p%pos(k,1)/dx
i = floor(xp)
dpx = xp - real(i, kind=8)
i = modulo(i, nx)
weight = p%w(k)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.0d0+dpx)
cp3x = f_m6(3.0d0-dpx)
cm2x = f_m6(2.0d0+dpx)
cp2x = f_m6(2.0d0-dpx)
cm1x = f_m6(1.0d0+dpx)
cp1x = f_m6(1.0d0-dpx)
cx   = f_m6(dpx)

f%rho(im3) = f%rho(im3) + cm3x * weight
f%rho(im2) = f%rho(im2) + cm2x * weight
f%rho(im1) = f%rho(im1) + cm1x * weight
f%rho(i  ) = f%rho(i  ) + cx   * weight
f%rho(ip1) = f%rho(ip1) + cp1x * weight
f%rho(ip2) = f%rho(ip2) + cp2x * weight
f%rho(ip3) = f%rho(ip3) + cp3x * weight

end do

f%rho(nx) = f%rho(0)
f%rho = f%rho / (dx)

end subroutine calcul_rho_m6


subroutine calcul_rho_n_m6( p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8)  :: i, k
integer(8)  :: im1, im2, im3, ip1, ip2, ip3
real(8)  :: dpx
real(8) :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)  :: weight, rho_total, xp

f%rho = 0.0_8

do k = 1, p%nbpa

xp = p%pos(k,1)/dx
i = floor(xp)
dpx = xp - real(i, kind=8)
i = modulo(i, nx)
weight = p%w(k)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.0d0+dpx)
cp3x = f_m6(3.0d0-dpx)
cm2x = f_m6(2.0d0+dpx)
cp2x = f_m6(2.0d0-dpx)
cm1x = f_m6(1.0d0+dpx)
cp1x = f_m6(1.0d0-dpx)
cx   = f_m6(dpx)

f%rho(im3) = f%rho(im3) + cm3x * weight
f%rho(im2) = f%rho(im2) + cm2x * weight
f%rho(im1) = f%rho(im1) + cm1x * weight
f%rho(i  ) = f%rho(i  ) + cx   * weight
f%rho(ip1) = f%rho(ip1) + cp1x * weight
f%rho(ip2) = f%rho(ip2) + cp2x * weight
f%rho(ip3) = f%rho(ip3) + cp3x * weight

end do

f%rho(nx) = f%rho(0)
f%rho = f%rho / (dx)
rho_total = sum(f%rho(0:nx-1))*dx
f%rho = f%rho - rho_total/dimx

end subroutine calcul_rho_n_m6


subroutine calcul_current_m6( p, f )

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8)  :: i, k
integer(8)  :: im1, im2, im3, ip1, ip2, ip3
real(8)  :: dpx
real(8) :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)  :: weight, xp, vx, vy

f%current = 0.0_8

do k = 1, p%nbpa

xp = p%pos(k,1)/dx
i = floor(xp)
dpx = xp - real(i, kind=8)
i = modulo(i, nx)
weight = p%w(k)
vx = p%vit(k,1)
vy = p%vit(k,2)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.0d0+dpx)
cp3x = f_m6(3.0d0-dpx)
cm2x = f_m6(2.0d0+dpx)
cp2x = f_m6(2.0d0-dpx)
cm1x = f_m6(1.0d0+dpx)
cp1x = f_m6(1.0d0-dpx)
cx   = f_m6(dpx)

f%current(1, im3) = f%current(1, im3) + cm3x * weight * vx
f%current(1, im2) = f%current(1, im2) + cm2x * weight * vx
f%current(1, im1) = f%current(1, im1) + cm1x * weight * vx
f%current(1, i  ) = f%current(1, i  ) + cx   * weight * vx
f%current(1, ip1) = f%current(1, ip1) + cp1x * weight * vx
f%current(1, ip2) = f%current(1, ip2) + cp2x * weight * vx
f%current(1, ip3) = f%current(1, ip3) + cp3x * weight * vx

f%current(2, im3) = f%current(2, im3) + cm3x * weight * vy
f%current(2, im2) = f%current(2, im2) + cm2x * weight * vy
f%current(2, im1) = f%current(2, im1) + cm1x * weight * vy
f%current(2, i  ) = f%current(2, i  ) + cx   * weight * vy
f%current(2, ip1) = f%current(2, ip1) + cp1x * weight * vy
f%current(2, ip2) = f%current(2, ip2) + cp2x * weight * vy
f%current(2, ip3) = f%current(2, ip3) + cp3x * weight * vy

end do

f%current(:, nx) = f%current(:, 0)
f%current = f%current / (dx)

end subroutine calcul_current_m6


subroutine calcul_dtcurrent_m6( p, dXdt, dVdt, f )

type(particle), intent(in)    :: p
real(8),        intent(in)    :: dXdt(:)
real(8),        intent(in)    :: dVdt(:, :)
type(field),    intent(inout) :: f

integer(8) :: k, i
integer(8) :: im1, im2, im3, ip1, ip2, ip3
real(8)    :: xp, dpx, weight
real(8)    :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)    :: dcx, dcm1x, dcm2x, dcm3x, dcp1x, dcp2x, dcp3x
real(8)    :: vx, vy, dvx, dvy, dXt

f%dtcurrent = 0.d0   

do k = 1, p%nbpa

xp   = p%pos(k,1) / dx
i    = floor(xp)
dpx  = xp - real(i,8)
i    = modulo(i, nx)

weight = p%w(k)
vx = p%vit(k,1)
vy = p%vit(k,2)
dvx = dVdt(k,1)
dvy = dVdt(k,2)
dXt  = dXdt(k)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.d0 + dpx)
cm2x = f_m6(2.d0 + dpx)
cm1x = f_m6(1.d0 + dpx)
cx   = f_m6(dpx)
cp1x = f_m6(1.d0 - dpx)
cp2x = f_m6(2.d0 - dpx)
cp3x = f_m6(3.d0 - dpx)

dcm3x = df_m6(3.d0 + dpx)
dcm2x = df_m6(2.d0 + dpx)
dcm1x = df_m6(1.d0 + dpx)
dcx   = df_m6(dpx)
dcp1x = -df_m6(1.d0 - dpx)
dcp2x = -df_m6(2.d0 - dpx)
dcp3x = -df_m6(3.d0 - dpx)

f%dtcurrent(1, im3) = f%dtcurrent(1, im3) + weight * ( dvx * cm3x - vx * dXt * dcm3x / dx )
f%dtcurrent(1, im2) = f%dtcurrent(1, im2) + weight * ( dvx * cm2x - vx * dXt * dcm2x / dx )
f%dtcurrent(1, im1) = f%dtcurrent(1, im1) + weight * ( dvx * cm1x - vx * dXt * dcm1x / dx )
f%dtcurrent(1, i  ) = f%dtcurrent(1, i  ) + weight * ( dvx * cx   - vx * dXt * dcx   / dx )
f%dtcurrent(1, ip1) = f%dtcurrent(1, ip1) + weight * ( dvx * cp1x - vx * dXt * dcp1x / dx )
f%dtcurrent(1, ip2) = f%dtcurrent(1, ip2) + weight * ( dvx * cp2x - vx * dXt * dcp2x / dx )
f%dtcurrent(1, ip3) = f%dtcurrent(1, ip3) + weight * ( dvx * cp3x - vx * dXt * dcp3x / dx )

f%dtcurrent(2, im3) = f%dtcurrent(2, im3) + weight * ( dvy * cm3x - vy * dXt * dcm3x / dx )
f%dtcurrent(2, im2) = f%dtcurrent(2, im2) + weight * ( dvy * cm2x - vy * dXt * dcm2x / dx )
f%dtcurrent(2, im1) = f%dtcurrent(2, im1) + weight * ( dvy * cm1x - vy * dXt * dcm1x / dx )
f%dtcurrent(2, i  ) = f%dtcurrent(2, i  ) + weight * ( dvy * cx   - vy * dXt * dcx   / dx )
f%dtcurrent(2, ip1) = f%dtcurrent(2, ip1) + weight * ( dvy * cp1x - vy * dXt * dcp1x / dx )
f%dtcurrent(2, ip2) = f%dtcurrent(2, ip2) + weight * ( dvy * cp2x - vy * dXt * dcp2x / dx )
f%dtcurrent(2, ip3) = f%dtcurrent(2, ip3) + weight * ( dvy * cp3x - vy * dXt * dcp3x / dx )

end do

f%dtcurrent(:, nx) = f%dtcurrent(:, 0)
f%dtcurrent = f%dtcurrent / dx

end subroutine calcul_dtcurrent_m6


subroutine calcul_dxcurrent_m6(p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8) :: k, i
integer(8) :: im1, im2, im3, ip1, ip2, ip3
real(8)    :: xp, dpx, weight, vx, vy
real(8)    :: dcm3x, dcm2x, dcm1x, dcx, dcp1x, dcp2x, dcp3x

f%dxcurrent = 0.d0

do k = 1, p%nbpa

xp  = p%pos(k,1) / dx
i   = floor(xp)
dpx = xp - real(i, kind=8)
i   = modulo(i, nx)

weight = p%w(k)
vx     = p%vit(k,1)
vy     = p%vit(k,2)

im3 = modulo(i-3, nx)
im2 = modulo(i-2, nx)
im1 = modulo(i-1, nx)
ip1 = modulo(i+1, nx)
ip2 = modulo(i+2, nx)
ip3 = modulo(i+3, nx)

dcm3x =  df_m6(3.d0 + dpx)
dcm2x =  df_m6(2.d0 + dpx)
dcm1x =  df_m6(1.d0 + dpx)
dcx   =  df_m6(dpx)
dcp1x = -df_m6(1.d0 - dpx)
dcp2x = -df_m6(2.d0 - dpx)
dcp3x = -df_m6(3.d0 - dpx)

f%dxcurrent(1, im3) = f%dxcurrent(1, im3) + weight * vx * dcm3x
f%dxcurrent(1, im2) = f%dxcurrent(1, im2) + weight * vx * dcm2x
f%dxcurrent(1, im1) = f%dxcurrent(1, im1) + weight * vx * dcm1x
f%dxcurrent(1, i  ) = f%dxcurrent(1, i  ) + weight * vx * dcx
f%dxcurrent(1, ip1) = f%dxcurrent(1, ip1) + weight * vx * dcp1x
f%dxcurrent(1, ip2) = f%dxcurrent(1, ip2) + weight * vx * dcp2x
f%dxcurrent(1, ip3) = f%dxcurrent(1, ip3) + weight * vx * dcp3x

f%dxcurrent(2, im3) = f%dxcurrent(2, im3) + weight * vy * dcm3x
f%dxcurrent(2, im2) = f%dxcurrent(2, im2) + weight * vy * dcm2x
f%dxcurrent(2, im1) = f%dxcurrent(2, im1) + weight * vy * dcm1x
f%dxcurrent(2, i  ) = f%dxcurrent(2, i  ) + weight * vy * dcx
f%dxcurrent(2, ip1) = f%dxcurrent(2, ip1) + weight * vy * dcp1x
f%dxcurrent(2, ip2) = f%dxcurrent(2, ip2) + weight * vy * dcp2x
f%dxcurrent(2, ip3) = f%dxcurrent(2, ip3) + weight * vy * dcp3x

end do

f%dxcurrent(:, nx) = f%dxcurrent(:, 0)
f%dxcurrent = f%dxcurrent / (dx * dx)

end subroutine calcul_dxcurrent_m6


subroutine calcul_fvxvy_m6(p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8) :: k, i, j
integer(8) :: im1, im2, im3, ip1, ip2, ip3
integer(8) :: jm1, jm2, jm3, jp1, jp2, jp3
real(8)    :: xv, yv, dpx, dpy, weight, vxmin, vymin, sigma1, sigma2
real(8)    :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)    :: cy, cm1y, cm2y, cm3y, cp1y, cp2y, cp3y

vxmin = -6d0
vymin = -6d0

f%fvxvy = 0.0_8

do k = 1, p%nbpa
xv = (p%vit(k,1) - vxmin) / (dvx)
yv = (p%vit(k,2) - vymin) / (dvy)
i = floor(xv)
j = floor(yv)
dpx = xv - real(i,8)
dpy = yv - real(j,8)
weight = p%w(k)

im3 = i-3; im2 = i-2; im1 = i-1
ip1 = i+1; ip2 = i+2; ip3 = i+3
jm3 = j-3; jm2 = j-2; jm1 = j-1
jp1 = j+1; jp2 = j+2; jp3 = j+3

if (im3 < 0 .or. ip3 > nvx .or. jm3 < 0 .or. jp3 > nvy) cycle

cm3x = f_m6(3.d0 + dpx)
cm2x = f_m6(2.d0 + dpx)
cm1x = f_m6(1.d0 + dpx)
cx   = f_m6(dpx)
cp1x = f_m6(1.d0 - dpx)
cp2x = f_m6(2.d0 - dpx)
cp3x = f_m6(3.d0 - dpx)
cm3y = f_m6(3.d0 + dpy)
cm2y = f_m6(2.d0 + dpy)
cm1y = f_m6(1.d0 + dpy)
cy   = f_m6(dpy)
cp1y = f_m6(1.d0 - dpy)
cp2y = f_m6(2.d0 - dpy)
cp3y = f_m6(3.d0 - dpy)

f%fvxvy(im3,jm3) = f%fvxvy(im3,jm3) + cm3x * cm3y * weight
f%fvxvy(im3,jm2) = f%fvxvy(im3,jm2) + cm3x * cm2y * weight
f%fvxvy(im3,jm1) = f%fvxvy(im3,jm1) + cm3x * cm1y * weight
f%fvxvy(im3,j )  = f%fvxvy(im3,j )  + cm3x * cy   * weight
f%fvxvy(im3,jp1) = f%fvxvy(im3,jp1) + cm3x * cp1y * weight
f%fvxvy(im3,jp2) = f%fvxvy(im3,jp2) + cm3x * cp2y * weight
f%fvxvy(im3,jp3) = f%fvxvy(im3,jp3) + cm3x * cp3y * weight

f%fvxvy(im2,jm3) = f%fvxvy(im2,jm3) + cm2x * cm3y * weight
f%fvxvy(im2,jm2) = f%fvxvy(im2,jm2) + cm2x * cm2y * weight
f%fvxvy(im2,jm1) = f%fvxvy(im2,jm1) + cm2x * cm1y * weight
f%fvxvy(im2,j  ) = f%fvxvy(im2,j )  + cm2x * cy   * weight
f%fvxvy(im2,jp1) = f%fvxvy(im2,jp1) + cm2x * cp1y * weight
f%fvxvy(im2,jp2) = f%fvxvy(im2,jp2) + cm2x * cp2y * weight
f%fvxvy(im2,jp3) = f%fvxvy(im2,jp3) + cm2x * cp3y * weight

f%fvxvy(im1,jm3) = f%fvxvy(im1,jm3) + cm1x * cm3y * weight
f%fvxvy(im1,jm2) = f%fvxvy(im1,jm2) + cm1x * cm2y * weight
f%fvxvy(im1,jm1) = f%fvxvy(im1,jm1) + cm1x * cm1y * weight
f%fvxvy(im1,j ) = f%fvxvy(im1,j  ) + cm1x * cy   * weight
f%fvxvy(im1,jp1) = f%fvxvy(im1,jp1) + cm1x * cp1y * weight
f%fvxvy(im1,jp2) = f%fvxvy(im1,jp2) + cm1x * cp2y * weight
f%fvxvy(im1,jp3) = f%fvxvy(im1,jp3) + cm1x * cp3y * weight

f%fvxvy(i  ,jm3) = f%fvxvy(i  ,jm3) + cx   * cm3y * weight
f%fvxvy(i  ,jm2) = f%fvxvy(i  ,jm2) + cx   * cm2y * weight
f%fvxvy(i  ,jm1) = f%fvxvy(i  ,jm1) + cx   * cm1y * weight
f%fvxvy(i  ,j  ) = f%fvxvy(i  ,j  ) + cx   * cy   * weight
f%fvxvy(i  ,jp1) = f%fvxvy(i  ,jp1) + cx   * cp1y * weight
f%fvxvy(i  ,jp2) = f%fvxvy(i  ,jp2) + cx   * cp2y * weight
f%fvxvy(i  ,jp3) = f%fvxvy(i  ,jp3) + cx   * cp3y * weight

f%fvxvy(ip1,jm3) = f%fvxvy(ip1,jm3) + cp1x * cm3y * weight
f%fvxvy(ip1,jm2) = f%fvxvy(ip1,jm2) + cp1x * cm2y * weight
f%fvxvy(ip1,jm1) = f%fvxvy(ip1,jm1) + cp1x * cm1y * weight
f%fvxvy(ip1,j ) = f%fvxvy(ip1,j  ) + cp1x * cy   * weight
f%fvxvy(ip1,jp1) = f%fvxvy(ip1,jp1) + cp1x * cp1y * weight
f%fvxvy(ip1,jp2) = f%fvxvy(ip1,jp2) + cp1x * cp2y * weight
f%fvxvy(ip1,jp3) = f%fvxvy(ip1,jp3) + cp1x * cp3y * weight

f%fvxvy(ip2,jm3) = f%fvxvy(ip2,jm3) + cp2x * cm3y * weight
f%fvxvy(ip2,jm2) = f%fvxvy(ip2,jm2) + cp2x * cm2y * weight
f%fvxvy(ip2,jm1) = f%fvxvy(ip2,jm1) + cp2x * cm1y * weight
f%fvxvy(ip2,j ) = f%fvxvy(ip2,j ) + cp2x * cy   * weight
f%fvxvy(ip2,jp1) = f%fvxvy(ip2,jp1) + cp2x * cp1y * weight
f%fvxvy(ip2,jp2) = f%fvxvy(ip2,jp2) + cp2x * cp2y * weight
f%fvxvy(ip2,jp3) = f%fvxvy(ip2,jp3) + cp2x * cp3y * weight

f%fvxvy(ip3,jm3) = f%fvxvy(ip3,jm3) + cp3x * cm3y * weight
f%fvxvy(ip3,jm2) = f%fvxvy(ip3,jm2) + cp3x * cm2y * weight
f%fvxvy(ip3,jm1) = f%fvxvy(ip3,jm1) + cp3x * cm1y * weight
f%fvxvy(ip3,j) = f%fvxvy(ip3,j) + cp3x * cy   * weight
f%fvxvy(ip3,jp1) = f%fvxvy(ip3,jp1) + cp3x * cp1y * weight
f%fvxvy(ip3,jp2) = f%fvxvy(ip3,jp2) + cp3x * cp2y * weight
f%fvxvy(ip3,jp3) = f%fvxvy(ip3,jp3) + cp3x * cp3y * weight

enddo

f%fvxvy = f%fvxvy / (dvx*dvy)

end subroutine calcul_fvxvy_m6


subroutine calcul_fxvx_m6(p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8) :: k
integer(8) :: i, iv
integer(8) :: im1, im2, im3, ip1, ip2, ip3
integer(8) :: ivm1, ivm2, ivm3, ivp1, ivp2, ivp3
real(8)    :: xp, dpx
real(8)    :: vp, dpv
real(8)    :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)    :: cv, cm1v, cm2v, cm3v, cp1v, cp2v, cp3v
real(8)    :: weight
real(8)    :: vxmin

vxmin = -6.0d0

f%fxvx = 0.0_8

do k = 1, p%nbpa
xp  = p%pos(k,1) / dx
i   = floor(xp)
dpx = xp - real(i,8)
i   = modulo(i, nx)

im3 = modulo(i-3, nx)
im2 = modulo(i-2, nx)
im1 = modulo(i-1, nx)
ip1 = modulo(i+1, nx)
ip2 = modulo(i+2, nx)
ip3 = modulo(i+3, nx)

cm3x = f_m6(3.d0 + dpx)
cm2x = f_m6(2.d0 + dpx)
cm1x = f_m6(1.d0 + dpx)
cx   = f_m6(dpx)
cp1x = f_m6(1.d0 - dpx)
cp2x = f_m6(2.d0 - dpx)
cp3x = f_m6(3.d0 - dpx)

vp  = (p%vit(k,1) - vxmin) / dvx
iv  = floor(vp)
dpv = vp - real(iv,8)

ivm3 = iv-3; ivm2 = iv-2; ivm1 = iv-1
ivp1 = iv+1; ivp2 = iv+2; ivp3 = iv+3

if (ivm3 < 0 .or. ivp3 > nvx) cycle

cm3v = f_m6(3.d0 + dpv)
cm2v = f_m6(2.d0 + dpv)
cm1v = f_m6(1.d0 + dpv)
cv   = f_m6(dpv)
cp1v = f_m6(1.d0 - dpv)
cp2v = f_m6(2.d0 - dpv)
cp3v = f_m6(3.d0 - dpv)

weight = p%w(k)

f%fxvx(im3,ivm3) = f%fxvx(im3,ivm3) + cm3x * cm3v * weight
f%fxvx(im3,ivm2) = f%fxvx(im3,ivm2) + cm3x * cm2v * weight
f%fxvx(im3,ivm1) = f%fxvx(im3,ivm1) + cm3x * cm1v * weight
f%fxvx(im3,iv  ) = f%fxvx(im3,iv  ) + cm3x * cv   * weight
f%fxvx(im3,ivp1) = f%fxvx(im3,ivp1) + cm3x * cp1v * weight
f%fxvx(im3,ivp2) = f%fxvx(im3,ivp2) + cm3x * cp2v * weight
f%fxvx(im3,ivp3) = f%fxvx(im3,ivp3) + cm3x * cp3v * weight

f%fxvx(im2,ivm3) = f%fxvx(im2,ivm3) + cm2x * cm3v * weight
f%fxvx(im2,ivm2) = f%fxvx(im2,ivm2) + cm2x * cm2v * weight
f%fxvx(im2,ivm1) = f%fxvx(im2,ivm1) + cm2x * cm1v * weight
f%fxvx(im2,iv  ) = f%fxvx(im2,iv  ) + cm2x * cv   * weight
f%fxvx(im2,ivp1) = f%fxvx(im2,ivp1) + cm2x * cp1v * weight
f%fxvx(im2,ivp2) = f%fxvx(im2,ivp2) + cm2x * cp2v * weight
f%fxvx(im2,ivp3) = f%fxvx(im2,ivp3) + cm2x * cp3v * weight

f%fxvx(im1,ivm3) = f%fxvx(im1,ivm3) + cm1x * cm3v * weight
f%fxvx(im1,ivm2) = f%fxvx(im1,ivm2) + cm1x * cm2v * weight
f%fxvx(im1,ivm1) = f%fxvx(im1,ivm1) + cm1x * cm1v * weight
f%fxvx(im1,iv  ) = f%fxvx(im1,iv  ) + cm1x * cv   * weight
f%fxvx(im1,ivp1) = f%fxvx(im1,ivp1) + cm1x * cp1v * weight
f%fxvx(im1,ivp2) = f%fxvx(im1,ivp2) + cm1x * cp2v * weight
f%fxvx(im1,ivp3) = f%fxvx(im1,ivp3) + cm1x * cp3v * weight

f%fxvx(i  ,ivm3) = f%fxvx(i  ,ivm3) + cx * cm3v * weight
f%fxvx(i  ,ivm2) = f%fxvx(i  ,ivm2) + cx * cm2v * weight
f%fxvx(i  ,ivm1) = f%fxvx(i  ,ivm1) + cx * cm1v * weight
f%fxvx(i  ,iv  ) = f%fxvx(i  ,iv  ) + cx * cv   * weight
f%fxvx(i  ,ivp1) = f%fxvx(i  ,ivp1) + cx * cp1v * weight
f%fxvx(i  ,ivp2) = f%fxvx(i  ,ivp2) + cx * cp2v * weight
f%fxvx(i  ,ivp3) = f%fxvx(i  ,ivp3) + cx * cp3v * weight

f%fxvx(ip2,ivm3) = f%fxvx(ip2,ivm3) + cp2x * cm3v * weight
f%fxvx(ip2,ivm2) = f%fxvx(ip2,ivm2) + cp2x * cm2v * weight
f%fxvx(ip2,ivm1) = f%fxvx(ip2,ivm1) + cp2x * cm1v * weight
f%fxvx(ip2,iv  ) = f%fxvx(ip2,iv  ) + cp2x * cv   * weight
f%fxvx(ip2,ivp1) = f%fxvx(ip2,ivp1) + cp2x * cp1v * weight
f%fxvx(ip2,ivp2) = f%fxvx(ip2,ivp2) + cp2x * cp2v * weight
f%fxvx(ip2,ivp3) = f%fxvx(ip2,ivp3) + cp2x * cp3v * weight

f%fxvx(ip1,ivm3) = f%fxvx(ip1,ivm3) + cp1x * cm3v * weight
f%fxvx(ip1,ivm2) = f%fxvx(ip1,ivm2) + cp1x * cm2v * weight
f%fxvx(ip1,ivm1) = f%fxvx(ip1,ivm1) + cp1x * cm1v * weight
f%fxvx(ip1,iv  ) = f%fxvx(ip1,iv  ) + cp1x * cv   * weight
f%fxvx(ip1,ivp1) = f%fxvx(ip1,ivp1) + cp1x * cp1v * weight
f%fxvx(ip1,ivp2) = f%fxvx(ip1,ivp2) + cp1x * cp2v * weight
f%fxvx(ip1,ivp3) = f%fxvx(ip1,ivp3) + cp1x * cp3v * weight


f%fxvx(ip3,ivm3) = f%fxvx(ip3,ivm3) + cp3x * cm3v * weight
f%fxvx(ip3,ivm2) = f%fxvx(ip3,ivm2) + cp3x * cm2v * weight
f%fxvx(ip3,ivm1) = f%fxvx(ip3,ivm1) + cp3x * cm1v * weight
f%fxvx(ip3,iv  ) = f%fxvx(ip3,iv  ) + cp3x * cv   * weight
f%fxvx(ip3,ivp1) = f%fxvx(ip3,ivp1) + cp3x * cp1v * weight
f%fxvx(ip3,ivp2) = f%fxvx(ip3,ivp2) + cp3x * cp2v * weight
f%fxvx(ip3,ivp3) = f%fxvx(ip3,ivp3) + cp3x * cp3v * weight

end do

f%fxvx = f%fxvx / (dx * dvx)

end subroutine calcul_fxvx_m6



subroutine calcul_energy(p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8)  :: i, k
integer(8)  :: im1, im2, im3, ip1, ip2, ip3
real(8)  :: dpx
real(8) :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8)  :: weight, rho_total, xp, velocity

f%rho = 0.0_8
do k = 1, p%nbpa

velocity=p%vit(k,1)**2+p%vit(k,2)**2

xp = p%pos(k,1)/dx
i = floor(xp)
dpx = xp - real(i, kind=8)
weight = p%w(k)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.0d0+dpx)
cp3x = f_m6(3.0d0-dpx)
cm2x = f_m6(2.0d0+dpx)
cp2x = f_m6(2.0d0-dpx)
cm1x = f_m6(1.0d0+dpx)
cp1x = f_m6(1.0d0-dpx)
cx   = f_m6(dpx)

f%rho(im3) = f%rho(im3) + cm3x * weight * velocity
f%rho(im2) = f%rho(im2) + cm2x * weight * velocity
f%rho(im1) = f%rho(im1) + cm1x * weight * velocity
f%rho(i   ) = f%rho(i   ) + cx   * weight * velocity
f%rho(ip1) = f%rho(ip1) + cp1x * weight * velocity
f%rho(ip2) = f%rho(ip2) + cp2x * weight * velocity
f%rho(ip3) = f%rho(ip3) + cp3x * weight * velocity

end do

f%rho(nx) = f%rho(0)
f%rho = f%rho / dx

end subroutine calcul_energy


subroutine calcul_DPRPIR(p, f, tau_now, H)

type(particle), intent(in)    :: p
type(field),    intent(in)    :: f
real(8),        intent(in)    :: tau_now
real(8),        intent(inout) :: H(3,0:nx)

integer(8) :: i, k
integer(8) :: im1, im2, im3, ip1, ip2, ip3
real(8) :: xp, dpx, weight
real(8) :: cx, cm1x, cm2x, cm3x, cp1x, cp2x, cp3x
real(8) :: W1, W2, Bp
real(8) :: coef1, coef2

H = 0.0_8

do k = 1, p%nbpa

xp = p%pos(k,1)/dx
i = floor(xp)
dpx = xp - real(i, kind=8)
i = modulo(i, nx)

weight = p%w(k)

W1 = p%vit(k,1)
W2 = p%vit(k,2)
Bp = p%mag_z(k)

coef1 = -weight * Bp * ( -dsin(tau_now)*W1 + dcos(tau_now)*W2 )
coef2 = -weight * Bp * ( -dcos(tau_now)*W1 - dsin(tau_now)*W2 )

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)

cm3x = f_m6(3.0d0+dpx)
cp3x = f_m6(3.0d0-dpx)
cm2x = f_m6(2.0d0+dpx)
cp2x = f_m6(2.0d0-dpx)
cm1x = f_m6(1.0d0+dpx)
cp1x = f_m6(1.0d0-dpx)
cx   = f_m6(dpx)

H(1,im3) = H(1,im3) + coef1 * cm3x
H(1,im2) = H(1,im2) + coef1 * cm2x
H(1,im1) = H(1,im1) + coef1 * cm1x
H(1,i  ) = H(1,i  ) + coef1 * cx
H(1,ip1) = H(1,ip1) + coef1 * cp1x
H(1,ip2) = H(1,ip2) + coef1 * cp2x
H(1,ip3) = H(1,ip3) + coef1 * cp3x

H(2,im3) = H(2,im3) + coef2 * cm3x
H(2,im2) = H(2,im2) + coef2 * cm2x
H(2,im1) = H(2,im1) + coef2 * cm1x
H(2,i  ) = H(2,i  ) + coef2 * cx
H(2,ip1) = H(2,ip1) + coef2 * cp1x
H(2,ip2) = H(2,ip2) + coef2 * cp2x
H(2,ip3) = H(2,ip3) + coef2 * cp3x

end do

H(3,:) = 0.0_8
H(:,nx) = H(:,0)
H = H / dx

end subroutine calcul_DPRPIR

end module interpolations
