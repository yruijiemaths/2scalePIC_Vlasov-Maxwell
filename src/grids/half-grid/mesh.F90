module mesh

type field
  real(8), dimension(:), allocatable :: exhalf, ey, deydx
  real(8), dimension(:), allocatable :: bz, dbzdx
  real(8), dimension(:), allocatable :: rho
  real(8), dimension(:,:), allocatable :: fvxvy, fxvx
  real(8), dimension(:), allocatable :: currentxhalf, dtcurrentxhalf, dxcurrentxhalf, currenty, dtcurrenty, dxcurrenty
end type field

 
integer(8) :: nx = 32
integer(8) :: nvx = 40, nvy = 40
real(8) :: ep
real(8) :: dx
real(8) :: dvx, dvy
real(8) :: dimx
real(8) :: tfinal,dt
real(8) :: pi


!real(8), dimension(3) :: eext

!contains !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!subroutine readin( donfil )
!
!implicit none
!
!character(len=*) :: donfil
!
!namelist/nlchem/ dimx,      &   ! domain dimensions
!                 dimy,      & 
!                 dimz,      &
!                 nx,        &   ! points along x
!                 ny,        &   ! points along y
!                 nz,        &   ! points along z
!                 dt,        & 
!                 tfinal,    &   ! final time
!                 nstepmax,  &   
!                 ndiag,     &   ! diagnostics
!                 c,         &   ! speed of light
!                 e0,        &   ! permittivy
!                 ldiach,    &   ! enable field plots
!                 !eext,      &   ! external electric field
!                 wp,        &   ! Plasma frequency
!                 v0,        &   ! Stream velocity
!                 vt             ! Thermal speed
!
!ldiach  = .false.
!!eext(:) = 0.0_8
!
!!rapbol  = pcharg / pmass ! default value, could be overwrite later
!
!!open(10,file=donfil,status='old')
!!read(10,nlchem) 
!!close(10)
!
!!write(6,donnees)
!dimx=4d0*pi
!dimy=2d0*pi
!dimz=1d0
!nx=64
!ny=32
!nz=4
!tfinal=pi/2.d0
!dt      = tfinal/2**3
!
!csq = c * c
!
!write(*,*)" size dimx       = ", dimx
!write(*,*)" size dimy       = ", dimy
!write(*,*)" size dimz       = ", dimz
!write(*,*)" time            = ", tfinal
!write(*,*)" nx              = ", nx
!write(*,*)" ny              = ", ny
!write(*,*)" nz              = ", nz
!write(*,*)" time step       = ", dt
!write(*,*)" speed of light  = ", c
!write(*,*)" tfinal          = ", tfinal
!
!dx = dimx / real(nx, kind=8)
!dy = dimy / real(ny, kind=8)
!dz = dimz / real(nz, kind=8)
!
!nstep = floor(tfinal/dt)
!!idiag = nstep/ndiag
!
!write(*,*)
!write(*,*) " dx = ", dx
!write(*,*) " dy = ", dy
!write(*,*) " dz = ", dz
!!write(*,*) " external electriques field "
!!write(*,*) " eext = ", eext(1:3)
!write(*,*)
!
!!nstep = merge( nstep , nstepmax, nstep < nstepmax ) 
!!write(*,*) " nstep = ", nstep
!
!end subroutine readin




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module mesh
