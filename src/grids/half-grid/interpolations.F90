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


function ddf_m6(q)
real(8), intent(in) :: q
real(8)             :: ddf_m6

if ( q < 1.0d0 ) then
   ddf_m6 = 20.d0*(3.d0-q)**3-120.d0*(2.d0-q)**3+300.d0*(1.d0-q)**3
else if ( q < 2.0d0 ) then
   ddf_m6 = 20.d0*(3.d0-q)**3-120.d0*(2.d0-q)**3
else if ( q < 3.0d0 ) then
   ddf_m6 = 20.d0*(3.d0-q)**3
else
   ddf_m6 = 0.d0
end if

ddf_m6 = ddf_m6 / 120.d0
return
end function ddf_m6


function f_m5(q)
real(8), intent(in) :: q
real(8)             :: f_m5

if ( q < 0.5d0 ) then
   f_m5 = (2.5d0-q)**4 - 5.0d0*(1.5d0-q)**4 + 10.0d0*(0.5d0-q)**4
else if ( q >= 0.5d0 .and. q < 1.5d0 ) then
   f_m5 = (2.5d0-q)**4 - 5.0d0*(1.5d0-q)**4
else if ( q >= 1.5d0 .and. q < 2.5d0 ) then
   f_m5 = (2.5d0-q)**4
else
   f_m5 = 0.0d0
end if

f_m5 = f_m5 / 24.0d0
return
end function f_m5


function df_m5(q)
real(8), intent(in) :: q
real(8)             :: df_m5

if ( q < 0.5d0 ) then
   df_m5 = -4.0d0*(2.5d0-q)**3 + 20.0d0*(1.5d0-q)**3 - 40.0d0*(0.5d0-q)**3
else if ( q >= 0.5d0 .and. q < 1.5d0 ) then
   df_m5 = -4.0d0*(2.5d0-q)**3 + 20.0d0*(1.5d0-q)**3
else if ( q >= 1.5d0 .and. q < 2.5d0 ) then
   df_m5 = -4.0d0*(2.5d0-q)**3
else
   df_m5 = 0.0d0
end if

df_m5 = df_m5 / 24.0d0

return
end function df_m5


subroutine interpol_m6(f, p)

type(particle), intent(inout) :: p
type(field),    intent(in)    :: f

integer(8) :: k
integer(8) :: i, ix
integer(8) :: im1, im2, im3, ip1, ip2, ip3
real(8)    :: dpx, xp
real(8)    :: dpx_ex, xp_ex
real(8)    :: cm3x, cm2x, cm1x, cx, cp1x, cp2x, cp3x
real(8)    :: cm3x_ex, cm2x_ex, cm1x_ex, cx_ex, cp1x_ex, cp2x_ex, cp3x_ex
real(8)    :: e(-2:2)

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

ex(i) = e(-2)*f%exhalf(im2) + e(-1)*f%exhalf(im1) + e(0)*f%exhalf(i) + e(1)*f%exhalf(ip1) + e(2)*f%exhalf(ip2)
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
p%ele_y(k) = cm3x*ey(im3) + cm2x*ey(im2) + cm1x*ey(im1) + cx*ey(i) + cp1x*ey(ip1) + cp2x*ey(ip2) + cp3x*ey(ip3)
p%mag_z(k) = cm3x*mz(im3) + cm2x*mz(im2) + cm1x*mz(im1) + cx*mz(i) + cp1x*mz(ip1) + cp2x*mz(ip2) + cp3x*mz(ip3)

xp_ex  = p%pos(k,1) / dx - 0.5d0
ix = floor(xp_ex)
dpx_ex = xp_ex - real(ix, kind=8)
ix = modulo(ix, nx)
im3 = modulo(ix-3, nx)
im2 = modulo(ix-2, nx)
im1 = modulo(ix-1, nx)
ip1 = modulo(ix+1, nx)
ip2 = modulo(ix+2, nx)
ip3 = modulo(ix+3, nx)
cm3x_ex = f_m6(3.0d0 + dpx_ex)
cm2x_ex = f_m6(2.0d0 + dpx_ex)
cm1x_ex = f_m6(1.0d0 + dpx_ex)
cx_ex   = f_m6(dpx_ex)
cp1x_ex = f_m6(1.0d0 - dpx_ex)
cp2x_ex = f_m6(2.0d0 - dpx_ex)
cp3x_ex = f_m6(3.0d0 - dpx_ex)
p%ele_x(k) = cm3x_ex*ex(im3) + cm2x_ex*ex(im2) + cm1x_ex*ex(im1) + cx_ex*ex(ix) + cp1x_ex*ex(ip1) + cp2x_ex*ex(ip2) + cp3x_ex*ex(ip3)

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


subroutine calcul_current(p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8)  :: i, ix, k
integer(8)  :: im1, im2, im3, ip1, ip2, ip3
integer(8)  :: im1x, im2x, ip1x, ip2x
real(8)     :: dpx, dpx_x
real(8)     :: cx, cm1x, cm2x, cp1x, cp2x
real(8)     :: cx_y, cm1y, cm2y, cm3y, cp1y, cp2y, cp3y
real(8)     :: weight, xp, xp_x, vx, vy

f%currenty = 0.0_8
f%currentxhalf = 0.0_8

do k = 1, p%nbpa

weight = p%w(k)
vx     = p%vit(k,1)
vy     = p%vit(k,2)

xp  = p%pos(k,1)/dx
i   = floor(xp)
dpx = xp - real(i, kind=8)
i   = modulo(i, nx)

im3 = modulo(i-3, nx)
im2 = modulo(i-2, nx)
im1 = modulo(i-1, nx)
ip1 = modulo(i+1, nx)
ip2 = modulo(i+2, nx)
ip3 = modulo(i+3, nx)
cm3y = f_m6(3.0d0 + dpx)
cm2y = f_m6(2.0d0 + dpx)
cm1y = f_m6(1.0d0 + dpx)
cx_y = f_m6(dpx)
cp1y = f_m6(1.0d0 - dpx)
cp2y = f_m6(2.0d0 - dpx)
cp3y = f_m6(3.0d0 - dpx)

f%currenty(im3) = f%currenty(im3) + cm3y * weight * vy
f%currenty(im2) = f%currenty(im2) + cm2y * weight * vy
f%currenty(im1) = f%currenty(im1) + cm1y * weight * vy
f%currenty(i  ) = f%currenty(i  ) + cx_y * weight * vy
f%currenty(ip1) = f%currenty(ip1) + cp1y * weight * vy
f%currenty(ip2) = f%currenty(ip2) + cp2y * weight * vy
f%currenty(ip3) = f%currenty(ip3) + cp3y * weight * vy


xp_x  = p%pos(k,1)/dx - 0.5d0
ix    = floor(xp_x)
dpx_x = xp_x - real(ix, kind=8)
ix    = modulo(ix, nx)

im2x = modulo(ix-2, nx)
im1x = modulo(ix-1, nx)
ip1x = modulo(ix+1, nx)
ip2x = modulo(ix+2, nx)
cm2x = f_m5(2.0d0 + dpx_x)
cm1x = f_m5(1.0d0 + dpx_x)
cx   = f_m5(dpx_x)
cp1x = f_m5(1.0d0 - dpx_x)
cp2x = f_m5(2.0d0 - dpx_x)

f%currentxhalf(im2x) = f%currentxhalf(im2x) + cm2x * weight * vx
f%currentxhalf(im1x) = f%currentxhalf(im1x) + cm1x * weight * vx
f%currentxhalf(ix  ) = f%currentxhalf(ix  ) + cx   * weight * vx
f%currentxhalf(ip1x) = f%currentxhalf(ip1x) + cp1x * weight * vx
f%currentxhalf(ip2x) = f%currentxhalf(ip2x) + cp2x * weight * vx

end do

f%currentxhalf(nx) = f%currentxhalf(0)
f%currentxhalf = f%currentxhalf / dx

f%currenty(nx) = f%currenty(0)
f%currenty = f%currenty / dx

end subroutine calcul_current


subroutine calcul_dtcurrent(p, dXdt, dVdt, f)

type(particle), intent(in)    :: p
real(8),        intent(in)    :: dXdt(:)
real(8),        intent(in)    :: dVdt(:, :)
type(field),    intent(inout) :: f

integer(8) :: k, i, ix
integer(8) :: im1, im2, im3, ip1, ip2, ip3
integer(8) :: im1x, im2x, ip1x, ip2x
real(8)    :: xp, xp_x, dpx, dpx_x, weight
real(8)    :: vx, vy, dvx, dvy, dXt

real(8)    :: cx_y, cm1y, cm2y, cm3y, cp1y, cp2y, cp3y
real(8)    :: dcx_y, dcm1y, dcm2y, dcm3y, dcp1y, dcp2y, dcp3y

real(8)    :: cx_x, cm1x, cm2x, cp1x, cp2x
real(8)    :: dcx_x, dcm1x, dcm2x, dcp1x, dcp2x

f%dtcurrenty     = 0.d0
f%dtcurrentxhalf = 0.d0

do k = 1, p%nbpa

weight = p%w(k)
vx     = p%vit(k,1)
vy     = p%vit(k,2)
dvx    = dVdt(k,1)
dvy    = dVdt(k,2)
dXt    = dXdt(k)

xp   = p%pos(k,1) / dx
i    = floor(xp)
dpx  = xp - real(i,8)
i    = modulo(i, nx)

im3 = modulo(i-3,nx)
im2 = modulo(i-2,nx)
im1 = modulo(i-1,nx)
ip1 = modulo(i+1,nx)
ip2 = modulo(i+2,nx)
ip3 = modulo(i+3,nx)
cm3y = f_m6(3.d0 + dpx)
cm2y = f_m6(2.d0 + dpx)
cm1y = f_m6(1.d0 + dpx)
cx_y = f_m6(dpx)
cp1y = f_m6(1.d0 - dpx)
cp2y = f_m6(2.d0 - dpx)
cp3y = f_m6(3.d0 - dpx)
dcm3y = df_m6(3.d0 + dpx)
dcm2y = df_m6(2.d0 + dpx)
dcm1y = df_m6(1.d0 + dpx)
dcx_y = df_m6(dpx)
dcp1y = -df_m6(1.d0 - dpx)
dcp2y = -df_m6(2.d0 - dpx)
dcp3y = -df_m6(3.d0 - dpx)

f%dtcurrenty(im3) = f%dtcurrenty(im3) + weight * ( dvy * cm3y - vy * dXt * dcm3y / dx )
f%dtcurrenty(im2) = f%dtcurrenty(im2) + weight * ( dvy * cm2y - vy * dXt * dcm2y / dx )
f%dtcurrenty(im1) = f%dtcurrenty(im1) + weight * ( dvy * cm1y - vy * dXt * dcm1y / dx )
f%dtcurrenty(i  ) = f%dtcurrenty(i  ) + weight * ( dvy * cx_y - vy * dXt * dcx_y / dx )
f%dtcurrenty(ip1) = f%dtcurrenty(ip1) + weight * ( dvy * cp1y - vy * dXt * dcp1y / dx )
f%dtcurrenty(ip2) = f%dtcurrenty(ip2) + weight * ( dvy * cp2y - vy * dXt * dcp2y / dx )
f%dtcurrenty(ip3) = f%dtcurrenty(ip3) + weight * ( dvy * cp3y - vy * dXt * dcp3y / dx )

xp_x  = p%pos(k,1) / dx - 0.5d0
ix    = floor(xp_x)
dpx_x = xp_x - real(ix,8)
ix    = modulo(ix, nx)

im2x = modulo(ix-2, nx)
im1x = modulo(ix-1, nx)
ip1x = modulo(ix+1, nx)
ip2x = modulo(ix+2, nx)
cm2x = f_m5(2.d0 + dpx_x)
cm1x = f_m5(1.d0 + dpx_x)
cx_x = f_m5(dpx_x)
cp1x = f_m5(1.d0 - dpx_x)
cp2x = f_m5(2.d0 - dpx_x)
dcm2x = df_m5(2.d0 + dpx_x)
dcm1x = df_m5(1.d0 + dpx_x)
dcx_x = df_m5(dpx_x)
dcp1x = -df_m5(1.d0 - dpx_x)
dcp2x = -df_m5(2.d0 - dpx_x)

f%dtcurrentxhalf(im2x) = f%dtcurrentxhalf(im2x) + weight * ( dvx * cm2x - vx * dXt * dcm2x / dx )
f%dtcurrentxhalf(im1x) = f%dtcurrentxhalf(im1x) + weight * ( dvx * cm1x - vx * dXt * dcm1x / dx )
f%dtcurrentxhalf(ix  ) = f%dtcurrentxhalf(ix  ) + weight * ( dvx * cx_x - vx * dXt * dcx_x / dx )
f%dtcurrentxhalf(ip1x) = f%dtcurrentxhalf(ip1x) + weight * ( dvx * cp1x - vx * dXt * dcp1x / dx )
f%dtcurrentxhalf(ip2x) = f%dtcurrentxhalf(ip2x) + weight * ( dvx * cp2x - vx * dXt * dcp2x / dx )

end do

f%dtcurrenty(nx)     = f%dtcurrenty(0)
f%dtcurrentxhalf(nx) = f%dtcurrentxhalf(0)

f%dtcurrenty     = f%dtcurrenty / dx
f%dtcurrentxhalf = f%dtcurrentxhalf / dx

end subroutine calcul_dtcurrent


subroutine calcul_dxcurrent(p, f)

type(particle), intent(in)    :: p
type(field),    intent(inout) :: f

integer(8) :: k, i
integer(8) :: im1, im2, im3, ip1, ip2, ip3
real(8)    :: xp, dpx, weight, vx, vy
real(8)    :: dcm3x, dcm2x, dcm1x, dcx, dcp1x, dcp2x, dcp3x

f%dxcurrentxhalf = 0.d0
f%dxcurrenty = 0.d0

do i = 0, nx-1
    f%dxcurrentxhalf(i) = ( f%currentxhalf(i) - f%currentxhalf(modulo(i-1,nx)) ) / dx
end do
f%dxcurrentxhalf(nx) = f%dxcurrentxhalf(0)


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

f%dxcurrenty(im3) = f%dxcurrenty(im3) + weight * vy * dcm3x
f%dxcurrenty(im2) = f%dxcurrenty(im2) + weight * vy * dcm2x
f%dxcurrenty(im1) = f%dxcurrenty(im1) + weight * vy * dcm1x
f%dxcurrenty(i  ) = f%dxcurrenty(i  ) + weight * vy * dcx
f%dxcurrenty(ip1) = f%dxcurrenty(ip1) + weight * vy * dcp1x
f%dxcurrenty(ip2) = f%dxcurrenty(ip2) + weight * vy * dcp2x
f%dxcurrenty(ip3) = f%dxcurrenty(ip3) + weight * vy * dcp3x

end do

f%dxcurrenty(nx) = f%dxcurrenty(0)
f%dxcurrenty = f%dxcurrenty / (dx * dx)

end subroutine calcul_dxcurrent


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
i = modulo(i, nx)
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


end module interpolations
