program vlasov_maxwell_pic_1dx2dv

use mesh
use particles
use interpolations
use poisson
use mpi
use pic_config, only: load_pic_config, apply_pic_config, print_pic_config, cfg_ntau, &
                      cfg_bz_two_stream, cfg_bz_weibel, cfg_wave_number
use pic_output, only: snapshot_output_due, diagnostic_output_due, &
                      poisson_tau_output_due, write_snapshot_outputs, &
                      write_diagnostic_outputs, write_poisson_tau
use, intrinsic :: iso_c_binding

implicit none
include 'fftw3.f'

type (field)         :: f
type (particle)      :: p
real(8)              :: cc, ss, ni, res_poi
integer              :: nnx
integer              :: prank
integer(8)           :: kpar, nstep, istep
integer              :: psize
integer              :: code
real(8)              :: t_start, dimvx, dimvy
real(8)              :: kpara=1d0
integer              :: i, j, k, l
real(8)              :: htau
integer(8)           :: fwtau, bwtau
integer              :: Ntau
real(8), allocatable :: rho_total(:), rho_total_new(:), fvxvy_total(:,:), fxvx_total(:,:), current_total(:,:), x_pos(:,:), v_vec(:,:), E_ini(:,:), B_ini(:), dExdx(:)
real(8), allocatable :: dxdt(:), dvdt(:,:), D_FRPIR(:,:), D_PRPIR(:,:,:), d2Eydx2(:,:,:), d2Bzdx2(:,:,:), d2Eydxdtau(:,:,:), d2Bzdxdtau(:,:,:), tmp_vector(:)
real(8), allocatable :: dExdx_tau(:,:,:), d2Exdxdtau_tau(:,:,:), dcurrentxdx_tau(:,:,:), rho_tau(:,:,:), drhodtau_tau(:,:,:)
complex(8), allocatable  :: dExdx_tau_temp(:,:,:), rho_tau_temp(:,:,:)
complex(8), allocatable  :: Up(:,:,:), Utempp(:,:,:), Upn_1(:,:,:), Utemppn_1(:,:,:), Uc(:,:,:), Utempc(:,:,:), Ucn_1(:,:,:), dUpdt(:,:,:), Uc_dExdx(:,:), Uc_dExdx_temp(:,:), Uc_Ex(:,:), Uc_Ex_n_1(:,:)
complex(8), allocatable  :: Sp(:,:,:), Sfp(:,:,:), Sc(:,:,:), Sfc(:,:,:), Sc1st(:,:,:), Sc2nd(:,:,:), PiSc(:,:)
complex(8), allocatable  :: Ele(:,:,:), Mag(:,:,:)!Particle
complex(8), allocatable  :: Cur(:,:,:), dEydx(:,:,:), dBzdx(:,:,:), dCurdt(:,:,:), dCurdx(:,:,:)!Cell
complex(8), allocatable  :: rho_total_tau(:,:), rho_total_new_tau(:,:), rho_total_old_tau(:,:)
complex(8), allocatable :: ptau(:), ptauh(:), fftin(:), fftout(:)
complex(8), allocatable :: Stau(:,:), Stautemp(:,:), Stau1(:,:)
complex(8), allocatable :: dUcdt(:,:,:)
real(8) , allocatable    ::  ktau(:), tau(:)
real(8) , allocatable    ::  energy(:), momentum(:,:), mass(:)
real(8) :: cpu_start, cpu_end, elapsed
integer :: file_unit
character(len=20) :: time_str
character(len=50) :: filename


call MPI_INIT(code)                             
t_start = MPI_WTIME()
call MPI_COMM_RANK(MPI_COMM_WORLD,prank,code)  
call MPI_COMM_SIZE(MPI_COMM_WORLD,psize,code)

!---parameters----
pi = 4.d0 * datan(1.d0)
dimx = 2.d0*pi/kpara
tfinal = 1000d0
!ep = 1.0d0/(2d0**10)
ep = 1.0d-2/4
!dt = 1.0d0/(2d0**4)
dt = 1.0d-1/2
call load_pic_config(prank)
Ntau=cfg_ntau
call apply_pic_config(nx,nvx,nvy,dimx,dimvx,dimvy,ep,dt,tfinal)
call print_pic_config(prank)
allocate(ptau(0:Ntau-1),ptauh(0:Ntau-1),fftin(0:Ntau-1),fftout(0:Ntau-1))
allocate(Stau(6,0:Ntau-1),Stautemp(6,0:Ntau-1),Stau1(6,0:Ntau-1))
nstep=nint(tfinal/dt)
if (nstep < 1_8) error stop 'tfinal/dt must produce at least one time step.'
dt=tfinal/real(nstep,8)
dx = dimx / real(nx, kind=8)

dvx = dimvx / real(nvx, kind=8)
dvy = dimvy / real(nvy, kind=8)


call MPI_BCAST(nx, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, code)
call initialize(p)
call init_poisson(nx)
nnx=nx

allocate(rho_total(0:nx))
allocate(rho_total_new(0:nx))
allocate(tmp_vector(0:nx))
allocate(fvxvy_total(0:nvx,0:nvy))
allocate(fxvx_total(0:nx,0:nvx))
allocate(current_total(2,0:nx))
allocate(dExdx(0:nx))
allocate(f%ex(0:nx)); f%ex = 0.0_8
allocate(f%ey(0:nx)); f%ey = 0.0_8
allocate(f%deydx(0:nx)); f%deydx = 0.0_8
allocate(f%bz(0:nx)); f%bz = 0.0_8
allocate(f%dbzdx(0:nx)); f%dbzdx = 0.0_8
allocate(f%rho(0:nx)); f%rho = 0.0_8
allocate(f%fvxvy(0:nvx,0:nvy)); f%fvxvy = 0.0_8
allocate(f%fxvx(0:nx,0:nvx)); f%fxvx = 0.0_8
allocate(f%current(2,0:nx)); f%current = 0.0_8
allocate(f%dtcurrent(2,0:nx)); f%dtcurrent = 0.0_8
allocate(f%dxcurrent(2,0:nx)); f%dxcurrent = 0.0_8

allocate(E_ini(2,0:nx)); E_ini = 0.0_8
allocate(B_ini(0:nx)); B_ini = 0.0_8

allocate(Up(3,0:Ntau-1,plocal))
allocate(dUpdt(3,0:Ntau-1,plocal))
allocate(dxdt(plocal))
allocate(dvdt(plocal,2))
allocate(Upn_1(3,0:Ntau-1,plocal))
allocate(Utemppn_1(3,0:Ntau-1,plocal))
allocate(Utempp(3,0:Ntau-1,plocal))
allocate(Sp(3,0:Ntau-1,plocal))
allocate(Sfp(3,0:Ntau-1,plocal))
allocate(Ele(2,0:Ntau-1,plocal))
allocate(Mag(1,0:Ntau-1,plocal))

allocate(Uc(3,0:Ntau-1,0:nx))
allocate(Uc_dExdx(0:Ntau-1,0:nx))
allocate(Uc_dExdx_temp(0:Ntau-1,0:nx))
allocate(Uc_Ex(0:Ntau-1,0:nx))
allocate(Uc_Ex_n_1(0:Ntau-1,0:nx))
allocate(Ucn_1(3,0:Ntau-1,0:nx))
allocate(Utempc(3,0:Ntau-1,0:nx))
allocate(Sc(3,0:Ntau-1,0:nx))
allocate(Sc1st(3,0:Ntau-1,0:nx))
allocate(Sc2nd(3,0:Ntau-1,0:nx))
allocate(PiSc(3,0:nx)); PiSc = 0.0_8
allocate(D_FRPIR(3,0:nx)); D_FRPIR = 0.0_8
allocate(D_PRPIR(3,0:Ntau-1,0:nx)); D_PRPIR = 0.0_8
allocate(Sfc(3,0:Ntau-1,0:nx))
allocate(Cur(2,0:Ntau-1,0:nx))
allocate(dCurdt(2,0:Ntau-1,0:nx))
allocate(dCurdx(2,0:Ntau-1,0:nx))
allocate(dEydx(1,0:Ntau-1,0:nx))
allocate(dBzdx(1,0:Ntau-1,0:nx))
allocate(d2Eydx2(1,0:Ntau-1,0:nx))
allocate(d2Bzdx2(1,0:Ntau-1,0:nx))
allocate(d2Eydxdtau(1,0:Ntau-1,0:nx))
allocate(d2Bzdxdtau(1,0:Ntau-1,0:nx))
allocate(dExdx_tau(0:nstep,0:Ntau-1,0:nx))
allocate(d2Exdxdtau_tau(0:nstep,0:Ntau-1,0:nx))
allocate(dcurrentxdx_tau(0:nstep,0:Ntau-1,0:nx))
allocate(rho_tau(0:nstep,0:Ntau-1,0:nx))
allocate(drhodtau_tau(0:nstep,0:Ntau-1,0:nx))
allocate(dExdx_tau_temp(0:nstep,0:Ntau-1,0:nx))
allocate(rho_tau_temp(0:nstep,0:Ntau-1,0:nx))
allocate(rho_total_tau(0:Ntau-1,0:nx))
allocate(rho_total_new_tau(0:Ntau-1,0:nx))
allocate(rho_total_old_tau(0:Ntau-1,0:nx))


allocate(x_pos(plocal,2))
allocate(v_vec(plocal,2))
allocate(energy(0:nstep))
allocate(mass(0:nstep))

allocate(dUcdt(3, 0:Ntau-1, 0:nx))

if (prank == 0) then
    print"(a,9g10.3)", ' dims = ', dimx, nx, dx
    print"(a,9g10.3)", ' dt= ',  dt, nstep
    print"(a,9g10.3)", ' ep= ',  ep
endif


htau=2*pi/Ntau
allocate(ktau(0:Ntau-1))
allocate(tau(0:Ntau-1))
call dfftw_plan_dft_1d(fwtau,Ntau,fftin,fftout,fftw_forward,fftw_estimate)
call dfftw_plan_dft_1d(bwtau,Ntau,fftin,fftout,fftw_backward,fftw_estimate)
do j=1,ntau
    if (j<=ntau/2) then
        i = j-1
    else
        i = -ntau+(j-1)
    endif
    ktau(j-1) = real(i,8)
    tau(j-1)=(j-1)*htau
enddo
do j=1,ntau-1
    ptauh(j)=(cdexp(-dt*cmplx(0d0,ktau(j),kind=8)/ep)-1.d0)/(-dt*cmplx(0d0,ktau(j),kind=8)/ep)!q_l
    ptau(j)=(cdexp(dt*cmplx(0d0,ktau(j),kind=8)/ep)-1.d0)/(dt*cmplx(0d0,ktau(j),kind=8)/ep)!q_l
enddo
ptau(0)=1d0
ptauh(0)=1d0
ptau=ptau*dt
ptauh=ptauh*dt

!-- Initialize fields and diagnostics in memory; no diagnostic files are written.
if (.true.) then
call calcul_ey(nx, dimx, f%ey)
call calcul_bz(nx, dimx, f%bz)
call calcul_rho_m6(p, f)
call local_to_global(f%rho, rho_total)
call remove_nyquist_from_rho(rho_total, nx)
call solve_poisson(rho_total, nx, dimx, f%ex)
ni=sum(rho_total(0:nx-1))*dx/dimx
mass(0)=sum(rho_total(0:nx-1))*dx

call calcul_dx(nx, dimx, f%ex, dExdx)

call calcul_energy(p, f)
call local_to_global(f%rho, rho_total)
energy(0)=sum(rho_total(0:nx-1))*dx/2.d0
do i=0,nx-1
    energy(0)=energy(0)+(f%ex(i)**2+f%ey(i)**2)*dx/2.d0+(f%bz(i)**2)*dx/2.d0
end do
endif

!---output solution---
if (.false.) then
call calcul_fvxvy_m6(p, f)
call local_to_global_fvxvy(f%fvxvy, fvxvy_total)
if (prank == 0 ) then
    open(unit=849,file='fvxvy_T0_ep0.05_dt0.05.dat')
    do i = 0, nvx
    do j = 0, nvy
        write(849,*) fvxvy_total(i,j)
    end do
    end do
    close(849)
endif

call calcul_fxvx_m6(p, f)
call local_to_global_fxvx(f%fxvx, fxvx_total)
if (prank == 0 ) then
    open(unit=850,file='fxvx_T0_ep0.05_dt0.05.dat')
    do i = 0, nx
    do j = 0, nvx
        write(850,*) fxvx_total(i,j)
    end do
    end do
    close(850)
endif
endif

call calcul_rho_m6(p, f)
call local_to_global(f%rho, rho_total)
call remove_nyquist_from_rho(rho_total, nx)
if (snapshot_output_due(0_8, nstep, dt)) call write_snapshot_outputs(p, f, 0.d0, prank, rho_total)
if (diagnostic_output_due(0_8, nstep)) call write_diagnostic_outputs(p, f, 0.d0, prank, rho_total)

!---prepare initial data for two-scale eq.---
x_pos(:,1)=p%pos(:,1)
v_vec(:,1)=p%vit(:,1)
v_vec(:,2)=p%vit(:,2)
call output_particles(p)
call calcul_ey(nx, dimx, f%ey)
call calcul_bz(nx, dimx, f%bz)
call calcul_rho_m6(p, f)
call local_to_global(f%rho, rho_total)
call solve_poisson(rho_total, nx, dimx, f%ex)
E_ini(1,:)=f%ex(:)
E_ini(2,:)=f%ey(:)
B_ini(:)=f%bz(:)

!---compute dJdt, dJdx at t=0---
do i=0,Ntau-1
    call interpol_m6(f, p)
    Ele(1,i,:)=p%ele_x
    Ele(2,i,:)=p%ele_y
    Mag(1,i,:)=p%mag_z
enddo
do kpar=1,plocal
    do i=0,Ntau-1
        cc=dcos(tau(i))
        ss=dsin(tau(i))
        dUpdt(1,i,kpar)=cc*v_vec(kpar,1)+ss*v_vec(kpar,2)
        dUpdt(2,i,kpar)=cc*Ele(1,i,kpar)-ss*Ele(2,i,kpar)+Mag(1,i,kpar)*v_vec(kpar,2)
        dUpdt(3,i,kpar)=ss*Ele(1,i,kpar)+cc*Ele(2,i,kpar)-Mag(1,i,kpar)*v_vec(kpar,1)
    enddo
enddo
do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    do kpar=1,plocal
        p%pos(kpar,1)=x_pos(kpar,1)
        p%vit(kpar,1)=cc*v_vec(kpar,1)+ss*v_vec(kpar,2)
        p%vit(kpar,2)=-ss*v_vec(kpar,1)+cc*v_vec(kpar,2)

        dxdt(kpar)=dreal(dUpdt(1,i,kpar))
        dvdt(kpar,1)=dreal(cc*dUpdt(2,i,kpar)-ss*dUpdt(3,i,kpar))
        dvdt(kpar,2)=dreal(ss*dUpdt(2,i,kpar)+cc*dUpdt(3,i,kpar))
    enddo
    call calcul_dtcurrent_m6(p, dxdt, dvdt, f)
    call local_to_global(f%dtcurrent(1,:), current_total(1,:))
    call local_to_global(f%dtcurrent(2,:), current_total(2,:))
    dCurdt(1,i,:)=current_total(1,:)
    dCurdt(2,i,:)=current_total(2,:)

    call calcul_dxcurrent_m6(p, f)
    call local_to_global(f%dxcurrent(1,:), current_total(1,:))
    call local_to_global(f%dxcurrent(2,:), current_total(2,:))
    dCurdx(1,i,:)=current_total(1,:)
    dCurdx(2,i,:)=current_total(2,:)
enddo
!---end compute dJdt, dJdx at t=0---

!for particle
do i=0,Ntau-1
    p%pos(:,1)=x_pos(:,1)
    p%vit(:,1)=dcos(tau(i))*v_vec(:,1)+dsin(tau(i))*v_vec(:,2)
    p%vit(:,2)=-dsin(tau(i))*v_vec(:,1)+dcos(tau(i))*v_vec(:,2)
    !!particle to cell
    call output_particles(p)
    call calcul_ey(nx, dimx, f%ey)
    call calcul_bz(nx, dimx, f%bz)
    call calcul_dx(nx, dimx, f%ey, f%deydx)
    call calcul_dx(nx, dimx, f%bz, f%dbzdx)
    call calcul_current_m6(p, f)
    call local_to_global(f%current(1,:), current_total(1,:))
    call local_to_global(f%current(2,:), current_total(2,:))
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total)
    call solve_poisson(rho_total, nx, dimx, f%ex)
    Cur(1,i,:)=current_total(1,:)
    Cur(2,i,:)=current_total(2,:)
    dEydx(1,i,:)=f%deydx
    dBzdx(1,i,:)=f%dbzdx
    call calcul_dx(nx, dimx, f%deydx, d2Eydx2(1,i,:))
    call calcul_dx(nx, dimx, f%dbzdx, d2Bzdx2(1,i,:))
    !!cell to particle
    call interpol_m6(f, p)!compute p%ele_xy, p%cur_xy, p%mag_z, p%dele_ydx, p%dmag_zdx by interpolation
    Ele(1,i,:)=p%ele_x
    Ele(2,i,:)=p%ele_y
    Mag(1,i,:)=p%mag_z
    call calcul_DPRPIR(p, f, tau(i), D_PRPIR(:,i,0:nx))
    call calcul_dx(nx, dimx, f%ex, dExdx_tau(0,i,0:nx))
    call calcul_dx(nx, dimx, current_total(1,:), dcurrentxdx_tau(0,i,0:nx))
    rho_tau(0,i,:)=rho_total
enddo

do k=0,nx
    fftin=dExdx_tau(0,:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    d2Exdxdtau_tau(0,:,k)=fftout
    d2Exdxdtau_tau(0,0,k) = 0.d0
    do i=1,Ntau-1
        d2Exdxdtau_tau(0,i,k)=cmplx(0d0,1d0,kind=8) * ktau(i) * d2Exdxdtau_tau(0,i,k) / Ntau
    enddo
    fftin=d2Exdxdtau_tau(0,:,k)
    call dfftw_execute_dft(bwtau, fftin, fftout)
    d2Exdxdtau_tau(0,:,k)=fftout

    fftin=rho_tau(0,:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    drhodtau_tau(0,:,k)=fftout
    drhodtau_tau(0,0,k) = 0.d0
    do i=1,Ntau-1
        drhodtau_tau(0,i,k)=cmplx(0d0,1d0,kind=8) * ktau(i) * drhodtau_tau(0,i,k) / Ntau
    enddo
    fftin=drhodtau_tau(0,:,k)
    call dfftw_execute_dft(bwtau, fftin, fftout)
    drhodtau_tau(0,:,k)=fftout
enddo

do kpar=1,plocal
    do i=0,Ntau-1
        cc=dcos(tau(i))
        ss=dsin(tau(i))
        Sp(1,i,kpar)=cc*v_vec(kpar,1)+ss*v_vec(kpar,2)!S
        Sp(2,i,kpar)=cc*Ele(1,i,kpar)-ss*Ele(2,i,kpar)+Mag(1,i,kpar)*v_vec(kpar,2)
        Sp(3,i,kpar)=ss*Ele(1,i,kpar)+cc*Ele(2,i,kpar)-Mag(1,i,kpar)*v_vec(kpar,1)
    enddo
    do l=1,3
        fftin=Sp(l,:,kpar)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Sfp(l,:,kpar)=fftout
    enddo
    do i=1,Ntau-1
        Sfp(:,i,kpar) = -cmplx(0d0,1d0,kind=8) * Sfp(:,i,kpar) / ktau(i)/ Ntau
    enddo
    Sfp(:,0,kpar)=0.d0
    do l=1,3
        fftin=Sfp(l,:,kpar)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Sp(l,:,kpar)=fftout
    enddo
    Up(1,:,kpar)=x_pos(kpar,1)+(Sp(1,:,kpar)-Sp(1,0,kpar))*ep
    Up(2,:,kpar)=v_vec(kpar,1)+(Sp(2,:,kpar)-Sp(2,0,kpar))*ep
    Up(3,:,kpar)=v_vec(kpar,2)+(Sp(3,:,kpar)-Sp(3,0,kpar))*ep
enddo

!for field
do k=0,nx
    do i=0,Ntau-1
        Sc(1,i,k)=-Cur(1,i,k)
        Sc(2,i,k)=-dBzdx(1,i,k)-Cur(2,i,k)
        Sc(3,i,k)=-dEydx(1,i,k)
    enddo

    !---compute \Pi Sc(0,tau,F(0))---
    do l = 1, 3
        PiSc(l,k) = 0.d0
        do i = 0, Ntau-1
            PiSc(l,k) = PiSc(l,k) + Sc(l,i,k)
        enddo
        PiSc(l,k) = PiSc(l,k) / Ntau
    enddo
    !---end compute \Pi Sc(0,tau,F(0))---

    !---compute \partial_{\tau,x}Ey, \partial_{\tau,x}Bz---
    fftin=dEydx(1,:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    d2Eydxdtau(1,:,k)=fftout
    d2Eydxdtau(1,0,k) = 0.d0
    fftin=dBzdx(1,:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    d2Bzdxdtau(1,:,k)=fftout
    d2Bzdxdtau(1,0,k) = 0.d0
    do i=1,Ntau-1
        d2Eydxdtau(1,i,k)=cmplx(0d0,1d0,kind=8) * ktau(i) * d2Eydxdtau(1,i,k) / Ntau
        d2Bzdxdtau(1,i,k)=cmplx(0d0,1d0,kind=8) * ktau(i) * d2Bzdxdtau(1,i,k) / Ntau
    enddo
    fftin=d2Eydxdtau(1,:,k)
    call dfftw_execute_dft(bwtau, fftin, fftout)
    d2Eydxdtau(1,:,k)=fftout
    fftin=d2Bzdxdtau(1,:,k)
    call dfftw_execute_dft(bwtau, fftin, fftout)
    d2Bzdxdtau(1,:,k)=fftout
    !---end compute \partial_{\tau,x}Ey, \partial_{\tau,x}Bz---

    do l=1,3
        fftin=Sc(l,:,k)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Sfc(l,:,k)=fftout
    enddo
    do i=1,Ntau-1
        Sfc(:,i,k)=-cmplx(0d0,1d0,kind=8) * Sfc(:,i,k) / ktau(i) / Ntau
    enddo
    Sfc(:,0,k) = 0.d0
    do l=1,3
        fftin=Sfc(l,:,k)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Sc(l,:,k)=fftout
    enddo
    Uc(1,:,k) = E_ini(1,k) + (Sc(1,:,k) - Sc(1,0,k)) * ep
    Uc(2,:,k) = E_ini(2,k) + (Sc(2,:,k) - Sc(2,0,k)) * ep
    Uc(3,:,k) = B_ini(k) + (Sc(3,:,k) - Sc(3,0,k)) * ep
enddo

!---second order---
if (.true.) then
do k=0,nx
    do i=0,Ntau-1
        cc=dcos(tau(i))
        ss=dsin(tau(i))
        p%pos(:,1)=dreal(Up(1,i,:))
        p%vit(:,1)=dreal(cc*Up(2,i,:)+ss*Up(3,i,:))
        p%vit(:,2)=dreal(-ss*Up(2,i,:)+cc*Up(3,i,:))
        f%ey(:)=dreal(Uc(2,i,:))
        f%bz(:)=dreal(Uc(3,i,:))
        call output_particles(p)
        call calcul_dx(nx, dimx, f%ey, f%deydx)
        call calcul_dx(nx, dimx, f%bz, f%dbzdx)
        call calcul_current_m6(p, f)
        call local_to_global(f%current(1,:), current_total(1,:))
        call local_to_global(f%current(2,:), current_total(2,:))
        Cur(1,i,:)=current_total(1,:)
        Cur(2,i,:)=current_total(2,:)
        dEydx(1,i,:)=f%deydx
        dBzdx(1,i,:)=f%dbzdx
        Sc1st(1,i,k)=-Cur(1,i,k)
        Sc1st(2,i,k)=-dBzdx(1,i,k)-Cur(2,i,k)
        Sc1st(3,i,k)=-dEydx(1,i,k)
    enddo
    do l=1,3
        fftin=Sc1st(l,:,k)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Sfc(l,:,k)=fftout
    enddo
    do i=1,Ntau-1
        Sfc(:,i,k) = -cmplx(0d0,1d0,kind=8) * Sfc(:,i,k) / ktau(i) / Ntau
    enddo
    Sfc(:,0,k) = 0.d0
    do l=1,3
        fftin=Sfc(l,:,k)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Sc1st(l,:,k)=fftout
    enddo
enddo
!---calcul \partial_F\Pi \Pi Sc(0,tau,F(0))---
D_FRPIR(1,:) = 0.d0
call calcul_dx(nx, dimx, -dreal(PiSc(3,:)), D_FRPIR(2,:))
call calcul_dx(nx, dimx, -dreal(PiSc(2,:)), D_FRPIR(3,:))

do k=0,nx
    do i=0,Ntau-1
        Sc(1,i,k) = D_FRPIR(1,k) + D_PRPIR(1,i,k) 
        Sc(2,i,k) = D_FRPIR(2,k) + D_PRPIR(2,i,k) 
        Sc(3,i,k) = D_FRPIR(3,k) + D_PRPIR(3,i,k) 
    enddo
    do l=1,3
        fftin=Sc(l,:,k)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Sfc(l,:,k)=fftout
    enddo
    do i = 1, Ntau-1
        Sfc(:,i,k) = - Sfc(:,i,k) / (ktau(i)*ktau(i)) / Ntau
    enddo
    Sfc(:,0,k) = 0.d0
    do l=1,3
        fftin=Sfc(l,:,k)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Sc2nd(l,:,k)=fftout
    enddo
    !Uc(1,:,k) = E_ini(1,k) + (Sc1st(1,:,k) - Sc1st(1,0,k)) * ep - (Sc2nd(1,:,k) - Sc2nd(1,0,k)) * ep * ep
    Uc(2,:,k) = E_ini(2,k) + (Sc1st(2,:,k) - Sc1st(2,0,k)) * ep - (Sc2nd(2,:,k) - Sc2nd(2,0,k)) * ep * ep
    Uc(3,:,k) = B_ini(k) + (Sc1st(3,:,k) - Sc1st(3,0,k)) * ep - (Sc2nd(3,:,k) - Sc2nd(3,0,k)) * ep * ep
enddo

do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    p%pos(:,1)=dreal(Up(1,i,:))
    p%vit(:,1)=dreal(cc*Up(2,i,:)+ss*Up(3,i,:))
    p%vit(:,2)=dreal(-ss*Up(2,i,:)+cc*Up(3,i,:))
    call output_particles(p)
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total)
    call remove_nyquist_from_rho(rho_total, nx)
    rho_tau(0,i,:)=rho_total
    call solve_poisson(rho_total, nx, dimx, f%ex)
    Uc(1,i,:) = dcmplx(f%ex(:), 0.d0)
enddo

endif

if (prank == 0) then
    print"(a,9g10.3)", 'Initialization done, start time loop'
endif

if (.true.) then
    do i=0,Ntau-1
        call calcul_dx(nx, dimx, dreal(Uc(1,i,0:nx)), dExdx_tau(0,i,0:nx))
    enddo

    res_poi = 0.d0
    do i=0,Ntau-1
        res_poi = max(sqrt(sum((dExdx_tau(0,i,0:nx-1) - rho_tau(0,i,0:nx-1) &
        + sum(rho_tau(0,i,0:nx-1))*dx/dimx)**2)*dx), res_poi )
        !res_poi=max(abs(dExdx_tau(0,i,k) - rho_tau(0,i,k) + sum(rho_tau(0,i,0:nx-1)) / real(nx,8)), res_poi)
    enddo

endif

if (poisson_tau_output_due(0_8, nstep, dt)) call write_poisson_tau(0.d0, res_poi, prank)

!---Loop over time---
call cpu_time(cpu_start)

tfinal=0.d0
do istep = 1, nstep
    tfinal=tfinal+dt
    if (istep == 1) then
    Utempp=Up
    Utempc=Uc
    do i=0,Ntau-1
        call calcul_dx(nx, dimx, dreal(Utempc(1,i,0:nx)), tmp_vector(0:nx))
        Uc_dExdx_temp(i,0:nx) = dcmplx(tmp_vector, 0.d0)
    enddo
    do kpar=1,plocal
        do l=1,3
            fftin=Up(l,:,kpar)
            call dfftw_execute_dft(fwtau, fftin, fftout)
            Stau(l,:)=fftout
        enddo
        Up(:,:,kpar)=Stau(1:3,:)/Ntau
    enddo
    do i=0,nx
        do l=1,3
            fftin=Uc(l,:,i)
            call dfftw_execute_dft(fwtau, fftin, fftout)
            Stau(l+3,:)=fftout
        enddo
        Uc(:,:,i)=Stau(4:6,:)/Ntau
    enddo
    do i=0,nx
        fftin=Uc_dExdx_temp(:,i)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Uc_dExdx_temp(:,i)=fftout/Ntau
    enddo
    call prediction()
    Up=Utempp
    Uc=Utempc

    else

    Utempp=Up
    Utemppn_1=Upn_1
    Utempc=Uc
    do i=0,Ntau-1
        call calcul_dx(nx, dimx, dreal(Ucn_1(1,i,0:nx)), tmp_vector(0:nx))
        Uc_dExdx_temp(i,0:nx) = dcmplx(tmp_vector, 0.d0)
    enddo
    do kpar=1,plocal
        do l=1,3
            fftin=Up(l,:,kpar)
            call dfftw_execute_dft(fwtau, fftin, fftout)
            Stau(l,:)=fftout
        enddo
        Up(:,:,kpar)=Stau(1:3,:)/Ntau
        do l=1,3
            fftin=Upn_1(l,:,kpar)
            call dfftw_execute_dft(fwtau, fftin, fftout)
            Stau(l,:)=fftout
        enddo
        Upn_1(:,:,kpar)=Stau(1:3,:)/Ntau
    enddo
    do i=0,nx
        do l=1,3
            fftin=Uc(l,:,i)
            call dfftw_execute_dft(fwtau, fftin, fftout)
            Stau(l+3,:)=fftout
        enddo
        Uc(:,:,i)=Stau(4:6,:)/Ntau
        do l=1,3
            fftin=Ucn_1(l,:,i)
            call dfftw_execute_dft(fwtau, fftin, fftout)
            Stau(l+3,:)=fftout
        enddo
        Ucn_1(:,:,i)=Stau(4:6,:)/Ntau
    enddo
    do i=0,nx
        fftin=Uc_dExdx_temp(:,i)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Uc_dExdx_temp(:,i)=fftout/Ntau
    enddo
    call correction()

    endif
    

    if (.false.) then
    do i=0,Ntau-1
        call calcul_dx(nx, dimx, dreal(Uc(1,i,0:nx)), dExdx_tau(istep,i,0:nx))
    enddo

    do i=0,Ntau-1
        cc=dcos(tau(i))
        ss=dsin(tau(i))
        p%pos(:,1)=dreal(Up(1,i,:))
        p%vit(:,1)=dreal(cc*Up(2,i,:)+ss*Up(3,i,:))
        p%vit(:,2)=dreal(-ss*Up(2,i,:)+cc*Up(3,i,:))
        call output_particles(p)
        call calcul_rho_m6(p, f)
        call local_to_global(f%rho, rho_total)
        call remove_nyquist_from_rho(rho_total, nx)
        rho_tau(istep,i,:)=rho_total
    enddo

    res_poi = 0.d0
    do i=0,Ntau-1
        res_poi = max(sqrt(sum((dExdx_tau(istep,i,0:nx-1) - rho_tau(istep,i,0:nx-1) &
        + sum(rho_tau(istep,i,0:nx-1))*dx/dimx)**2)*dx), res_poi )
        !res_poi=max(abs(dExdx_tau(istep,i,k) - rho_tau(istep,i,k) + sum(rho_tau(istep,i,0:nx-1))*dx/dimx), res_poi)
    enddo

    if ( prank == 0 ) then
        open(unit=85,file='T1000_ep0.0025_dt0.05_Poissonres_tau.dat', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        write(85,*) res_poi
        close(85)
    endif

    call energyuse()
    mass(istep)=sum(rho_total(0:nx-1))*dx
    if (prank == 0 ) then
        open(unit=85,file='T1000_ep0.0025_dt0.05_mass.dat', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        write(85,*) mass(istep)
        close(85)
    endif

    call calcul_dx(nx, dimx, f%ex, dExdx)
    if (prank == 0 ) then
        open(unit=85,file='T1000_ep0.0025_dt0.05_Poissonres.dat', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        write(85,*) sqrt(sum((dExdx(0:nx-1)-rho_total(0:nx-1)+ni)**2d0)*dx)
        close(85)
    endif
    
    call calcul_energy(p, f)
    call local_to_global(f%rho, rho_total)
    energy(istep)=sum(rho_total(0:nx-1))*dx/2.d0
    do i=0,nx-1
        energy(istep)=energy(istep)+(f%ex(i)**2+f%ey(i)**2)*dx/2.d0+(f%bz(i)**2)*dx/2.d0
    end do
    if (prank == 0 ) then
        open(unit=852,file='T1000_ep0.0025_dt0.05_energy.dat', STATUS='OLD', ACTION='WRITE', POSITION='APPEND')
        write(852,*) energy(istep)
        close(852)
    endif


    endif

    if (.false.) then
    if (mod(istep, 2000_8) == 0) then
        call energyuse()
        call calcul_fvxvy_m6(p, f)
        call local_to_global_fvxvy(f%fvxvy, fvxvy_total)
        if (prank == 0) then
        write(filename, '("fvxvy_T", I0, "_ep0.01_dt0.05.dat")') nint(tfinal)
        open(newunit=file_unit, file=trim(filename), status='replace', action='write')
        do i = 0, nvx
        do j = 0, nvy
            write(file_unit,*) fvxvy_total(i,j)
        end do
        end do
        close(file_unit)
        end if

        call calcul_fxvx_m6(p, f)
        call local_to_global_fxvx(f%fxvx, fxvx_total)
        if (prank == 0) then
        write(filename, '("fxvx_T", I0, "_ep0.01_dt0.05.dat")') nint(tfinal)
        open(newunit=file_unit, file=trim(filename), status='replace', action='write')
        do i = 0, nx
        do j = 0, nvx
            write(file_unit,*) fxvx_total(i,j)
        end do
        end do
        close(file_unit)
        end if

    endif
    endif

    if (poisson_tau_output_due(istep, nstep, dt)) then
        do i=0,Ntau-1
            call calcul_dx(nx, dimx, dreal(Uc(1,i,0:nx)), dExdx_tau(istep,i,0:nx))
        enddo
        do i=0,Ntau-1
            cc=dcos(tau(i)); ss=dsin(tau(i))
            p%pos(:,1)=dreal(Up(1,i,:))
            p%vit(:,1)=dreal(cc*Up(2,i,:)+ss*Up(3,i,:))
            p%vit(:,2)=dreal(-ss*Up(2,i,:)+cc*Up(3,i,:))
            call output_particles(p)
            call calcul_rho_m6(p, f)
            call local_to_global(f%rho, rho_total)
            call remove_nyquist_from_rho(rho_total, nx)
            rho_tau(istep,i,:)=rho_total
        enddo
        res_poi=0.d0
        do i=0,Ntau-1
            res_poi=max(res_poi,sqrt(sum((dExdx_tau(istep,i,0:nx-1)- &
                rho_tau(istep,i,0:nx-1)+sum(rho_tau(istep,i,0:nx-1))*dx/dimx)**2)*dx))
        enddo
        call write_poisson_tau(tfinal, res_poi, prank)
    endif

    if (snapshot_output_due(istep, nstep, dt) .or. &
        diagnostic_output_due(istep, nstep)) then
        call energyuse()
        if (snapshot_output_due(istep, nstep, dt)) &
            call write_snapshot_outputs(p, f, real(istep,8)*dt, prank, rho_total)
        if (diagnostic_output_due(istep, nstep)) &
            call write_diagnostic_outputs(p, f, real(istep,8)*dt, prank, rho_total)
    endif

    if (prank == 0 .and. mod(istep, 100_8) == 0) then
        print"(a,9g10.3)", ' process, iter= ',  istep, nstep
        !print"(a,9g10.3)", ' energy = ', dabs(energy(istep)-energy(0))/dabs(energy(0))
    endif

end do

call cpu_time(cpu_end)
elapsed = cpu_end - cpu_start

if (prank == 0) then
    print '(A, F10.3)', 'Total CPU time (s): ', elapsed
end if

if (.false.) then
call energyuse()
if (prank == 0) then
    open(unit=85,file='rho_T1_ep2-10_dt0.0005.dat')
    do i = 0, nx
        write(85,*) rho_total(i)
    end do
    close(85)

    open(unit=86,file='ExEy_T1_ep2-10_dt0.0005.dat')
    do i = 0, nx
        write(86,*) f%ex(i), f%ey(i)
    end do
    close(86)

    open(unit=86,file='Bz_T1_ep2-10_dt0.0005.dat')
    do i = 0, nx
        write(86,*) f%bz(i)
    end do
    close(86)
end if
call calcul_energy(p, f)
call local_to_global(f%rho, rho_total)
if (prank == 0 ) then
    open(unit=87,file='rhov_T1_ep2-10_dt0.0005.dat')
    do i = 0, nx
        write(87,*) rho_total(i)
    end do
    close(87)
endif
endif

call free_poisson()
call dfftw_destroy_plan(fwtau)
call dfftw_destroy_plan(bwtau)
deallocate(rho_total)
deallocate(rho_total_new)
deallocate(tmp_vector)
deallocate(fvxvy_total)
deallocate(current_total)
deallocate(dExdx)
deallocate(tau)
deallocate(ktau)
deallocate(energy)
deallocate(mass)
deallocate(Up)
deallocate(dUpdt)
deallocate(dxdt)
deallocate(dvdt)
deallocate(Upn_1)
deallocate(Utemppn_1)
deallocate(Utempp)
deallocate(Sp)
deallocate(Sfp)
deallocate(Ele)
deallocate(Mag)
deallocate(Uc)
deallocate(Ucn_1)
deallocate(Uc_dExdx)
deallocate(Uc_dExdx_temp)
deallocate(Uc_Ex)
deallocate(Uc_Ex_n_1)
deallocate(Utempc)
deallocate(Sc)
deallocate(Sfc)
deallocate(Sc1st)
deallocate(Sc2nd)
deallocate(PiSc)
deallocate(D_FRPIR)
deallocate(D_PRPIR)
deallocate(Cur)
deallocate(dCurdt)
deallocate(dCurdx)
deallocate(dEydx)
deallocate(dBzdx)
deallocate(d2Eydx2)
deallocate(d2Bzdx2)
deallocate(d2Eydxdtau)
deallocate(d2Bzdxdtau)
deallocate(dExdx_tau)
deallocate(d2Exdxdtau_tau)
deallocate(dcurrentxdx_tau)
deallocate(rho_tau)
deallocate(drhodtau_tau)
deallocate(dExdx_tau_temp)
deallocate(rho_tau_temp)
deallocate(x_pos)
deallocate(v_vec)
deallocate(dUcdt)
deallocate(E_ini)
deallocate(B_ini)
deallocate(rho_total_tau)
deallocate(rho_total_new_tau)
deallocate(rho_total_old_tau)
call MPI_BARRIER(MPI_COMM_WORLD,code)
call MPI_FINALIZE(code)

stop

!-----------------------------------------------------------------------

contains

!
! Sum values across all processors and put
! the result on the processor 0
!
subroutine local_to_global( local_array2d, global_array2d )

real(8), dimension(:), intent(in)  :: local_array2d
real(8), dimension(:), intent(out) :: global_array2d

call MPI_ALLREDUCE(local_array2d,        &
global_array2d,       &
nnx+1,        &
MPI_REAL8,            &
MPI_SUM,              &
MPI_COMM_WORLD,       &
code)

end subroutine local_to_global


subroutine local_to_global_fvxvy(fvxvy_local, fvxvy_global)

use mpi
implicit none

real(8), dimension(:,:), intent(in)  :: fvxvy_local
real(8), dimension(:,:), intent(out) :: fvxvy_global

integer :: ierr
integer :: nsend

nsend = size(fvxvy_local)

call MPI_ALLREDUCE( fvxvy_local,   &
                    fvxvy_global,  &
                    nsend,        &
                    MPI_REAL8,    &
                    MPI_SUM,      &
                    MPI_COMM_WORLD, &
                    ierr )

end subroutine local_to_global_fvxvy


subroutine local_to_global_fxvx(fxvx_local, fxvx_global)

use mpi
implicit none

real(8), dimension(:,:), intent(in)  :: fxvx_local
real(8), dimension(:,:), intent(out) :: fxvx_global

integer :: ierr
integer :: nsend

nsend = size(fxvx_local)

call MPI_ALLREDUCE( fxvx_local,   &
                    fxvx_global,  &
                    nsend,        &
                    MPI_REAL8,    &
                    MPI_SUM,      &
                    MPI_COMM_WORLD, &
                    ierr )

end subroutine local_to_global_fxvx


subroutine output_particles( p )

type(particle), intent(inout) :: p
integer :: ipart

do ipart = 1, p%nbpa !periodic BC

do while (p%pos(ipart,1)>=dimx)
p%pos(ipart,1) = p%pos(ipart,1) - dimx
end do

do while (p%pos(ipart,1)<0.0_8)
p%pos(ipart,1) = p%pos(ipart,1) + dimx
end do

end do

end subroutine output_particles


subroutine prediction()
do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    p%pos(:,1)=dreal(Utempp(1,i,:))
    p%vit(:,1)=dreal(cc*Utempp(2,i,:)+ss*Utempp(3,i,:))
    p%vit(:,2)=dreal(-ss*Utempp(2,i,:)+cc*Utempp(3,i,:))
    f%ex(:)=dreal(Utempc(1,i,:))
    f%ey(:)=dreal(Utempc(2,i,:))
    f%bz(:)=dreal(Utempc(3,i,:))
    call output_particles(p)
    call calcul_dx(nx, dimx, f%ey, f%deydx)
    call calcul_dx(nx, dimx, f%bz, f%dbzdx)
    call calcul_current_m6(p, f)
    call local_to_global(f%current(1,:), current_total(1,:))
    call local_to_global(f%current(2,:), current_total(2,:))
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total)
    call remove_nyquist_from_rho(rho_total, nx)
    Cur(1,i,:)=current_total(1,:)
    Cur(2,i,:)=current_total(2,:)
    dEydx(1,i,:)=f%deydx
    dBzdx(1,i,:)=f%dbzdx
    call interpol_m6(f, p)
    Ele(1,i,:)=p%ele_x
    Ele(2,i,:)=p%ele_y
    Mag(1,i,:)=p%mag_z
    rho_total_tau(i,:)=rho_total(:)
enddo
do kpar=1,plocal
    do i=0,Ntau-1
        cc=dcos(tau(i))
        ss=dsin(tau(i))
        Stautemp(1,i)=cc*Utempp(2,i,kpar)+ss*Utempp(3,i,kpar)
        Stautemp(2,i)=cc*Ele(1,i,kpar)-ss*Ele(2,i,kpar)+Mag(1,i,kpar)*Utempp(3,i,kpar)
        Stautemp(3,i)=ss*Ele(1,i,kpar)+cc*Ele(2,i,kpar)-Mag(1,i,kpar)*Utempp(2,i,kpar)
    enddo
    do l=1,3
        fftin=Stautemp(l,:)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Stau(l,:)=fftout
    enddo
    Stau=Stau/Ntau 
    do j=0,Ntau-1
        Stau(1:3,j)=cdexp(-dt*cmplx(0d0,ktau(j),kind=8)/ep)*Up(:,j,kpar)+ptauh(j)*Stau(1:3,j)
    enddo
    do l=1,3
        fftin=Stau(l,:)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Stau1(l,:)=fftout
    enddo
    Upn_1(:,:,kpar)=Utempp(:,:,kpar)
    Utempp(:,:,kpar)=Stau1(1:3,:)
enddo

do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    p%pos(:,1)=dreal(Utempp(1,i,:))
    p%vit(:,1)=dreal(cc*Utempp(2,i,:)+ss*Utempp(3,i,:))
    p%vit(:,2)=dreal(-ss*Utempp(2,i,:)+cc*Utempp(3,i,:))
    call output_particles(p)
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total_new)
    call remove_nyquist_from_rho(rho_total_new, nx)
    rho_total_new_tau(i,:)=rho_total_new(:)
enddo
do k=0,nx
    fftin=rho_total_tau(:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    rho_total_tau(:,k)=fftout/Ntau
    fftin=rho_total_new_tau(:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    rho_total_new_tau(:,k)=fftout/Ntau
    do i=0,Ntau-1
        Uc_dExdx(i,k)=cdexp(-dt*cmplx(0d0,ktau(i),kind=8)/ep)*(Uc_dExdx_temp(i,k)-rho_total_tau(i,k))+rho_total_new_tau(i,k)
    enddo
    fftin=Uc_dExdx(:,k)
    call dfftw_execute_dft(bwtau, fftin, fftout)
    Uc_dExdx(:,k)=fftout
enddo
do i=0,Ntau-1
    call calcul_inv_dx(nx, dimx, dreal(Uc_dExdx(i,:)), 0.d0, tmp_vector)
    Uc_Ex(i,:) = dcmplx(tmp_vector, 0.d0)
enddo

do k=0,nx
    do i=0,Ntau-1
        Stautemp(5,i)=-dBzdx(1,i,k)-Cur(2,i,k)
        Stautemp(6,i)=-dEydx(1,i,k)
    enddo
    do l=5,6
        fftin=Stautemp(l,:)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Stau(l,:)=fftout
    enddo
    Stau=Stau/Ntau
    do j=0,Ntau-1
        Stau(5:6,j)=cdexp(-dt*cmplx(0d0,ktau(j),kind=8)/ep)*Uc(2:3,j,k)+ptauh(j)*Stau(5:6,j)
    enddo
    do l=5,6
        fftin=Stau(l,:)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Stau1(l,:)=fftout
    enddo
    Ucn_1(:,:,k)=Utempc(:,:,k)
    Utempc(1,:,k)=Uc_Ex(:,k)
    Utempc(2:3,:,k)=Stau1(5:6,:)
enddo
end subroutine prediction


subroutine correction()
do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    p%pos(:,1)=dreal(Utempp(1,i,:))
    p%vit(:,1)=dreal(cc*Utempp(2,i,:)+ss*Utempp(3,i,:))
    p%vit(:,2)=dreal(-ss*Utempp(2,i,:)+cc*Utempp(3,i,:))
    f%ex(:)=dreal(Utempc(1,i,:))
    f%ey(:)=dreal(Utempc(2,i,:))
    f%bz(:)=dreal(Utempc(3,i,:))
    call output_particles(p)
    call calcul_dx(nx, dimx, f%ey, f%deydx)
    call calcul_dx(nx, dimx, f%bz, f%dbzdx)
    call calcul_current_m6(p, f)
    call local_to_global(f%current(1,:), current_total(1,:))
    call local_to_global(f%current(2,:), current_total(2,:))
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total)
    call remove_nyquist_from_rho(rho_total, nx)
    Cur(1,i,:)=current_total(1,:)
    Cur(2,i,:)=current_total(2,:)
    dEydx(1,i,:)=f%deydx
    dBzdx(1,i,:)=f%dbzdx
    call interpol_m6( f, p)
    Ele(1,i,:)=p%ele_x
    Ele(2,i,:)=p%ele_y
    Mag(1,i,:)=p%mag_z
enddo

do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    p%pos(:,1)=dreal(Utemppn_1(1,i,:))
    p%vit(:,1)=dreal(cc*Utemppn_1(2,i,:)+ss*Utemppn_1(3,i,:))
    p%vit(:,2)=dreal(-ss*Utemppn_1(2,i,:)+cc*Utemppn_1(3,i,:))
    call output_particles(p)
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total)
    call remove_nyquist_from_rho(rho_total, nx)
    rho_total_tau(i,:)=rho_total(:)
enddo

do kpar=1,plocal
    do i=0,Ntau-1
        cc=dcos(tau(i))
        ss=dsin(tau(i))
        Stautemp(1,i)=cc*Utempp(2,i,kpar)+ss*Utempp(3,i,kpar)
        Stautemp(2,i)=cc*Ele(1,i,kpar)-ss*Ele(2,i,kpar)+Mag(1,i,kpar)*Utempp(3,i,kpar)
        Stautemp(3,i)=ss*Ele(1,i,kpar)+cc*Ele(2,i,kpar)-Mag(1,i,kpar)*Utempp(2,i,kpar)
    enddo
    do l=1,3
        fftin=Stautemp(l,:)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Stau(l,:)=fftout
    enddo
    Stau=Stau/Ntau 
    do j=0,Ntau-1
        Stau(1:3,j)=cdexp(-2*dt*cmplx(0d0,ktau(j),kind=8)/ep)*Upn_1(:,j,kpar)+ptauh(j)*Stau(1:3,j)+cdexp(-2*dt*cmplx(0d0,ktau(j),kind=8)/ep)*ptau(j)*Stau(1:3,j)
    enddo
    do l=1,3
        fftin=Stau(l,:)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Stau1(l,:)=fftout
    enddo  
    Upn_1(:,:,kpar)=Utempp(:,:,kpar)
    Up(:,:,kpar)=Stau1(1:3,:)
enddo

do i=0,Ntau-1
    cc=dcos(tau(i))
    ss=dsin(tau(i))
    p%pos(:,1)=dreal(Up(1,i,:))
    p%vit(:,1)=dreal(cc*Up(2,i,:)+ss*Up(3,i,:))
    p%vit(:,2)=dreal(-ss*Up(2,i,:)+cc*Up(3,i,:))
    call output_particles(p)
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total_new)
    call remove_nyquist_from_rho(rho_total_new, nx)
    rho_total_new_tau(i,:)=rho_total_new(:)
enddo
do k=0,nx
    fftin=rho_total_tau(:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    rho_total_tau(:,k)=fftout/Ntau
    fftin=rho_total_new_tau(:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    rho_total_new_tau(:,k)=fftout/Ntau
    do i=0,Ntau-1
        Uc_dExdx(i,k)=cdexp(-2*dt*cmplx(0d0,ktau(i),kind=8)/ep)*(Uc_dExdx_temp(i,k)-rho_total_tau(i,k))+rho_total_new_tau(i,k)
    enddo
    fftin=Uc_dExdx(:,k)
    call dfftw_execute_dft(bwtau, fftin, fftout)
    Uc_dExdx(:,k)=fftout
enddo
do i=0,Ntau-1
    call calcul_inv_dx(nx, dimx, dreal(Uc_dExdx(i,:)), 0.d0, tmp_vector)
    Uc_Ex(i,:) = dcmplx(tmp_vector, 0.d0)
enddo

do k=0,nx
    do i=0,Ntau-1
        Stautemp(5,i)=-dBzdx(1,i,k)-Cur(2,i,k)
        Stautemp(6,i)=-dEydx(1,i,k)
    enddo
    do l=5,6
        fftin=Stautemp(l,:)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Stau(l,:)=fftout
    enddo
    Stau=Stau/Ntau
    do j=0,Ntau-1
        Stau(5:6,j)=cdexp(-2*dt*cmplx(0d0,ktau(j),kind=8)/ep)*Ucn_1(2:3,j,k)+ptauh(j)*Stau(5:6,j)+cdexp(-2*dt*cmplx(0d0,ktau(j),kind=8)/ep)*ptau(j)*Stau(5:6,j)
    enddo
    do l=5,6
        fftin=Stau(l,:)
        call dfftw_execute_dft(bwtau, fftin, fftout)
        Stau1(l,:)=fftout
    enddo
    Ucn_1(:,:,k)=Utempc(:,:,k)
    Uc(1,:,k)=Uc_Ex(:,k)
    Uc(2:3,:,k)=Stau1(5:6,:)
enddo
end subroutine correction


subroutine energyuse()
! Reconstruct charge from the two-scale density at tau=t/epsilon.  This is
! the same spectral trace used for the fields and therefore preserves the
! two-scale Gauss relation.  Depositing after tracing particle positions is
! not equivalent because the particle shape function is nonlinear in x.
do i=0,Ntau-1
    p%pos(:,1)=dreal(Up(1,i,:))
    call calcul_rho_m6(p, f)
    call local_to_global(f%rho, rho_total)
    call remove_nyquist_from_rho(rho_total, nx)
    rho_total_tau(i,:)=dcmplx(rho_total(:),0.d0)
enddo
do k=0,nx
    fftin=rho_total_tau(:,k)
    call dfftw_execute_dft(fwtau, fftin, fftout)
    rho_total(k)=0.d0
    do i=0,Ntau-1
        rho_total(k)=rho_total(k)+dreal(fftout(i)/Ntau* &
            cdexp(cmplx(0d0,ktau(i),kind=8)*tfinal/ep))
    enddo
enddo
do kpar=1,plocal
    do l=1,3
        fftin=Up(l,:,kpar)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Stautemp(l,:)=fftout
    enddo
    Stau(1:3,0)=0.d0
    do i=0,Ntau-1
        Stau(1:3,0)=Stautemp(1:3,i)/Ntau*cdexp(cmplx(0d0,ktau(i),kind=8)*tfinal/ep)+Stau(1:3,0)
    enddo
    cc=dcos(tfinal/ep)
    ss=dsin(tfinal/ep)
    p%pos(kpar,1)=dreal(Stau(1,0))
    p%vit(kpar,1)=dreal(cc*Stau(2,0)+ss*Stau(3,0))
    p%vit(kpar,2)=dreal(-ss*Stau(2,0)+cc*Stau(3,0))
enddo
do k=0,nx
    do l=1,3
        fftin=Uc(l,:,k)
        call dfftw_execute_dft(fwtau, fftin, fftout)
        Stautemp(l+3,:)=fftout
    enddo
    Stau(4:6,0)=0.d0
    do i=0,Ntau-1
        Stau(4:6,0)=Stautemp(4:6,i)/Ntau*cdexp(cmplx(0d0,ktau(i),kind=8)*tfinal/ep)+Stau(4:6,0)
    enddo
    f%ex(k)=dreal(Stau(4,0))
    f%ey(k)=dreal(Stau(5,0))
    f%bz(k)=dreal(Stau(6,0))
enddo
call output_particles(p)
call calcul_current_m6(p, f)
end subroutine energyuse


subroutine calcul_ey(nx, dimx, ey)
!---calculate electronic field on the cell---
integer(8), intent(in)    :: nx
real(8), intent(in)    :: dimx
real(8), intent(inout) :: ey(0:nx)
integer     :: i
real(8)     :: x, dx

dx = dimx / real(nx, kind=8)
select case (initial_condition)
case (initial_original_two_stream)
    ! Original two-stream field initial condition.
    do i=1,nx
        x=i*dx
        ey(i)=0.d0
    enddo
case (initial_anisotropic_weibel)
    ! Anisotropic Maxwellian field initial condition.
    do i=1,nx
        x=i*dx
        ey(i)=0.d0
    enddo
case (initial_bump_on_tail)
    ey(1:nx)=0.d0
end select
ey(0) = ey(nx)

end subroutine calcul_ey


subroutine calcul_bz(nx, dimx, bz)
!---calculate magnetic field on the cell---
integer(8), intent(in)    :: nx
real(8), intent(in)    :: dimx
real(8), intent(inout) :: bz(0:nx)
integer     :: i
real(8)     :: x, dx, beta, k

dx = dimx / real(nx, kind=8)
k=cfg_wave_number
select case (initial_condition)
case (initial_original_two_stream)
    beta=cfg_bz_two_stream
case (initial_anisotropic_weibel)
    beta=cfg_bz_weibel
case (initial_bump_on_tail)
    beta=0.d0
end select
do i=1,nx
    x=i*dx
    !bz(i)=beta*dcos(k*x)
    bz(i)=beta*dsin(k*x)
enddo
bz(0) = bz(nx)

end subroutine calcul_bz

end program vlasov_maxwell_pic_1dx2dv
