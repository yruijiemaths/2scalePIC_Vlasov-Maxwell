program vlasov_maxwell_pic_1dx2dv

use mesh
use particles
use interpolations
use poisson, only: init_poisson, free_poisson, calcul_dx, solve_poisson_face
use mpi
use pic_config, only: load_pic_config, apply_pic_config, print_pic_config, cfg_ntau, &
                      cfg_bz_two_stream, cfg_bz_weibel, cfg_wave_number
use pic_output, only: snapshot_output_due, diagnostic_output_due, &
                      poisson_tau_output_due, write_snapshot_outputs, &
                      write_diagnostic_outputs, write_poisson_tau
use, intrinsic :: iso_fortran_env, only: output_unit

implicit none
include 'fftw3.f'

integer :: Ntau
type(field) :: f
type(particle) :: p

integer :: ierr, prank, psize
integer :: itau, ix, icomp, istep
integer(8) :: ip, nstep
integer(8) :: plan_tau_forward, plan_tau_backward
integer(8) :: plan_x_forward, plan_x_backward
real(8) :: htau, kpara, dimvx, dimvy
real(8) :: continuity_linf, gauss_linf
real(8) :: particle_solve_residual
real(8) :: ni
real(8) :: cpu_start, cpu_end, elapsed
real(8) :: stage_time_local(6), stage_time_global(6), wall_mark
logical, parameter :: output_final_solution = .false.

real(8), allocatable :: tau(:), ktau(:)
real(8), allocatable :: rho_old(:,:), rho_new(:,:)
real(8), allocatable :: rho_base(:,:)
real(8), allocatable :: jx_face(:,:)
real(8), allocatable :: work_face(:)
real(8), allocatable :: currenty_global(:)
real(8), allocatable :: jy_old(:,:)
real(8), allocatable :: rho_total(:), x_pos(:,:), v_vec(:,:)
real(8), allocatable :: E_ini(:,:), B_ini(:)

complex(8), allocatable :: Up(:,:,:), Up_new(:,:,:), Fp(:,:,:)
complex(8), allocatable :: Ufield(:,:,:), Ufield_new(:,:,:)
complex(8), allocatable :: Fey(:,:), Fbz(:,:)





complex(8), allocatable :: field_second_source(:,:,:)
complex(8), allocatable :: xfftin(:), xfftout(:)
complex(8), allocatable :: fftin(:), fftout(:)
complex(8), allocatable :: anti_tau(:), anti_tau2(:), anti_tau_twice(:)
complex(8), allocatable :: dtau_matrix(:,:)
character(len=256) :: arg, state_file, ep_label, dt_label, diag_prefix
character(len=256) :: poisson_tau_file, poisson_file, energy_file
integer :: arg_status

call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, prank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, psize, ierr)

!-----------------------------------------------------------------------
! Parameters.  They intentionally mirror the existing TSI drivers.
!-----------------------------------------------------------------------
pi = 4.d0*datan(1.d0)
kpara = 1.d0
dimx = 2.d0*pi/kpara
tfinal = 1000.d0
ep = 1.d-3
dt = 1.d-4
state_file = ''
call load_pic_config(prank)
Ntau=cfg_ntau
dimvx = 12.d0
dimvy = 12.d0
call apply_pic_config(nx,nvx,nvy,dimx,dimvx,dimvy,ep,dt,tfinal)
call print_pic_config(prank)
allocate(fftin(0:Ntau-1),fftout(0:Ntau-1))
allocate(anti_tau(0:Ntau-1),anti_tau2(0:Ntau-1),anti_tau_twice(0:Ntau-1))
call get_command_argument(1,arg,status=arg_status)
if (arg_status == 0 .and. len_trim(arg) > 0) then
    read(arg,*) ep
endif
call get_command_argument(2,arg,status=arg_status)
if (arg_status == 0 .and. len_trim(arg) > 0) then
    read(arg,*) dt
endif
call get_command_argument(3,arg,status=arg_status)
if (arg_status == 0 .and. len_trim(arg) > 0) read(arg,*) tfinal
call get_command_argument(4,state_file,status=arg_status)
if (ep <= 0.d0 .or. dt <= 0.d0 .or. tfinal <= 0.d0) &
    error stop 'epsilon, dt and tfinal must be positive.'
nstep = nint(tfinal/dt)
if (nstep < 1_8) error stop 'tfinal/dt must produce at least one time step.'
dt = tfinal/real(nstep,8)
call compact_real_label(ep,ep_label)
call compact_real_label(dt,dt_label)
write(diag_prefix,'("T",i0,"_ep",a,"_dt",a)') nint(tfinal),trim(ep_label),trim(dt_label)
poisson_tau_file = trim(diag_prefix)//'_Poissonres_tau.dat'
poisson_file = trim(diag_prefix)//'_Poissonres.dat'
energy_file = trim(diag_prefix)//'_energy.dat'
dx = dimx/real(nx,8)
dvx = dimvx/real(nvx,8)
dvy = dimvy/real(nvy,8)

call MPI_BCAST(nx, 1, MPI_INTEGER8, 0, MPI_COMM_WORLD, ierr)
call initialize(p)
call allocate_field_arrays()

allocate(tau(0:Ntau-1), ktau(0:Ntau-1))
allocate(rho_old(0:Ntau-1,0:nx), rho_new(0:Ntau-1,0:nx))
allocate(jx_face(0:Ntau-1,0:nx)); jx_face = 0.d0
allocate(work_face(0:nx))
allocate(currenty_global(0:nx))
allocate(jy_old(0:Ntau-1,0:nx))
allocate(rho_total(0:nx))
allocate(x_pos(plocal,1), v_vec(plocal,2))
allocate(E_ini(2,0:nx), B_ini(0:nx))

allocate(Up(3,0:Ntau-1,1:plocal))
allocate(Up_new(3,0:Ntau-1,1:plocal))
allocate(Fp(3,0:Ntau-1,1:plocal))
allocate(Ufield(3,0:Ntau-1,0:nx))
allocate(Ufield_new(3,0:Ntau-1,0:nx))
allocate(Fey(0:Ntau-1,0:nx), Fbz(0:Ntau-1,0:nx))





allocate(field_second_source(2,0:Ntau-1,0:nx))
allocate(xfftin(0:nx-1),xfftout(0:nx-1))

htau = 2.d0*pi/real(Ntau,8)
do itau = 0, Ntau-1
    tau(itau) = real(itau,8)*htau
    if (itau <= Ntau/2) then
        ktau(itau) = real(itau,8)
    else
        ktau(itau) = real(itau-Ntau,8)
    endif
enddo
ktau(Ntau/2) = 0.d0

call dfftw_plan_dft_1d(plan_tau_forward, Ntau, fftin, fftout, fftw_forward, fftw_estimate)
call dfftw_plan_dft_1d(plan_tau_backward, Ntau, fftin, fftout, fftw_backward, fftw_estimate)
call dfftw_plan_dft_1d(plan_x_forward, int(nx), xfftin, xfftout, fftw_forward, fftw_estimate)
call dfftw_plan_dft_1d(plan_x_backward, int(nx), xfftin, xfftout, fftw_backward, fftw_estimate)
allocate(dtau_matrix(0:Ntau-1,0:Ntau-1))
call build_dtau_matrix(dtau_matrix)
call init_poisson(nx)

!-----------------------------------------------------------------------
! Physical initial data.  This block follows the layout used by the
! original TSI main programs.
!-----------------------------------------------------------------------
x_pos(:,1) = p%pos(:,1)
v_vec(:,1) = p%vit(:,1)
v_vec(:,2) = p%vit(:,2)

call calcul_ey(nx, dimx, f%ey)
call calcul_bz(nx, dimx, f%bz)
call calcul_rho_m6(p, f)
call allreduce_real(f%rho, rho_total)
call solve_poisson_face(rho_total,nx,dx,E_ini(1,:))
E_ini(2,:) = f%ey(:)
B_ini(:) = f%bz(:)

!-----------------------------------------------------------------------
! Prepare initial data for the two-scale equation.
!-----------------------------------------------------------------------
allocate(rho_base(0:Ntau-1,0:nx))

! Zeroth-order two-scale profile.
do itau = 0, Ntau-1
    Up(1,itau,:) = cmplx(x_pos(:,1),0.d0,kind=8)
    Up(2,itau,:) = cmplx(v_vec(:,1),0.d0,kind=8)
    Up(3,itau,:) = cmplx(v_vec(:,2),0.d0,kind=8)
    Ufield(1,itau,:) = cmplx(E_ini(1,:),0.d0,kind=8)
    Ufield(2,itau,:) = cmplx(E_ini(2,:),0.d0,kind=8)
    Ufield(3,itau,:) = cmplx(B_ini(:),0.d0,kind=8)
enddo

! Evaluate the zeroth-order right-hand side.
call evaluate_explicit_rhs(Up, Ufield, Fp, Fey, Fbz)

! Terms D_F Pi S + D_P Pi S required by the second-order
! well-prepared transverse-field profile.  They must be evaluated at
! the uncorrected physical initial data.
call compute_second_order_field_source(Up, Ufield, Fey, Fbz, &
                                       field_second_source)

! First-order prepared particle profile:
! Up(0,tau) = Up(0,0) + ep*(A*Fp(tau)-A*Fp(0)).
do ip = 1, plocal
    do icomp = 1, 3
        call tau_zero_mean_antiderivative(Fp(icomp,:,ip),anti_tau)
        Up(icomp,:,ip) = Up(icomp,:,ip) + ep* &
                         (anti_tau-anti_tau(0))
    enddo
enddo

! First-order predictor for the transverse-field profile.  This
! predictor is used to evaluate S^[1] in the second-order correction.
do ix = 0, nx-1
    call tau_zero_mean_antiderivative(Fey(:,ix),anti_tau)
    Ufield(2,:,ix) = Ufield(2,:,ix) + ep*(anti_tau-anti_tau(0))
    call tau_zero_mean_antiderivative(Fbz(:,ix),anti_tau)
    Ufield(3,:,ix) = Ufield(3,:,ix) + ep*(anti_tau-anti_tau(0))
enddo
Ufield(2:3,:,nx) = Ufield(2:3,:,0)

! Second-order prepared transverse fields, following the construction
! used by the TSI drivers:
!   Uc = uc0 + ep*(A*S^[1](tau)-A*S^[1](0))
!             - ep^2*(A^2*G(tau)-A^2*G(0)),
! where G = D_F Pi S + D_P Pi S.  The subtraction at tau=0 preserves
! exactly the prescribed physical initial Ey and Bz.
call evaluate_explicit_rhs(Up, Ufield, Fp, Fey, Fbz)
do ix = 0, nx-1
    call tau_zero_mean_antiderivative(Fey(:,ix),anti_tau)
    call tau_zero_mean_antiderivative(field_second_source(1,:,ix),anti_tau2)
    call tau_zero_mean_antiderivative(anti_tau2,anti_tau_twice)
    Ufield(2,:,ix) = cmplx(E_ini(2,ix),0.d0,kind=8) + ep*(anti_tau-anti_tau(0)) - ep*ep*(anti_tau_twice-anti_tau_twice(0))

    call tau_zero_mean_antiderivative(Fbz(:,ix),anti_tau)
    call tau_zero_mean_antiderivative(field_second_source(2,:,ix),anti_tau2)
    call tau_zero_mean_antiderivative(anti_tau2,anti_tau_twice)
    Ufield(3,:,ix) = cmplx(B_ini(ix),0.d0,kind=8) + ep*(anti_tau-anti_tau(0)) - ep*ep*(anti_tau_twice-anti_tau_twice(0))
enddo
Ufield(2:3,:,nx) = Ufield(2:3,:,0)

! Recompute Ex from the corrected particle positions so that the
! prepared profile satisfies the discrete Gauss law at every tau.
call deposit_density_all_tau(Up, rho_base)
call solve_gauss_all_tau(rho_base, Ufield)
deallocate(rho_base)

call deposit_density_all_tau(Up, rho_old)

ni = sum(rho_old(0,0:nx-1))*dx/dimx
call solve_gauss_all_tau(rho_old, Ufield)

call compute_diagnostics(rho_old, rho_old, jx_face, Ufield, continuity_linf, gauss_linf, .true.)
if (poisson_tau_output_due(0_8, nstep, dt)) call write_poisson_tau(0.d0, gauss_linf, prank)
! Initialization routines reuse p and f while sweeping over tau.  Reconstruct
! the physical t=0 slice before the initial snapshot and diagnostics are
! emitted so their contents do not depend on the last tau node visited.
call extract_physical_state(0.d0)
if (snapshot_output_due(0_8, nstep, dt)) call write_snapshot_outputs(p, f, 0.d0, prank, rho_total)
if (diagnostic_output_due(0_8, nstep)) call write_diagnostic_outputs(p, f, 0.d0, prank, rho_total)
call cpu_time(cpu_start)
if (prank == 0) then



    write(*,'(a,es12.4,a,i0,a,i0)') 'TSFV-CPimex: ep=',ep, ', Ntau=',Ntau,', MPI ranks=',psize

    write(*,'(a,f10.4,a,es12.4,a,i0)') 'tfinal=',tfinal,', dt=',dt,', nstep=',nstep
    write(*,'(a)') 'phase-space output is controlled by config/pic.nml'
    write(*,'(a)') 'particle solver : direct dense complex solve'
    flush(output_unit)
endif

!-----------------------------------------------------------------------
! Time loop.
!-----------------------------------------------------------------------
do istep = 1, nstep
    stage_time_local = 0.d0
    wall_mark = MPI_WTIME()
    call evaluate_explicit_rhs(Up, Ufield, Fp, Fey, Fbz)
    stage_time_local(1) = MPI_WTIME()-wall_mark

    ! First-order IMEX Euler.  The frozen-coefficient magnetic rotation
    ! and the linear Maxwell wave operator are on the implicit side;
    ! interpolation coefficients and Jy are evaluated at time n.
    wall_mark = MPI_WTIME()
    call particle_imex_update(Up, Ufield, Up_new)
    stage_time_local(2) = MPI_WTIME()-wall_mark
    wall_mark = MPI_WTIME()
    call maxwell_imex_update(Ufield, jy_old, Ufield_new)
    stage_time_local(3) = MPI_WTIME()-wall_mark

    wall_mark = MPI_WTIME()
    call deposit_density_all_tau(Up_new, rho_new)
    stage_time_local(4) = MPI_WTIME()-wall_mark
    wall_mark = MPI_WTIME()
    call deposit_charge_preserving_jx(Up(1,:,:), Up_new(1,:,:), jx_face)
    stage_time_local(5) = MPI_WTIME()-wall_mark
    wall_mark = MPI_WTIME()
    call update_exhalf(Ufield(1,:,:), jx_face, Ufield_new(1,:,:))
    stage_time_local(6) = MPI_WTIME()-wall_mark

    call MPI_ALLREDUCE(stage_time_local,stage_time_global,6,MPI_REAL8,MPI_MAX, &
                       MPI_COMM_WORLD,ierr)

    call compute_diagnostics(rho_old, rho_new, jx_face, Ufield_new, continuity_linf, gauss_linf, .false.)
    if (poisson_tau_output_due(int(istep,8), nstep, dt)) &
        call write_poisson_tau(real(istep,8)*dt, gauss_linf, prank)

    if (prank == 0) then
        if (mod(istep,1000) == 0 .or. istep == 1 .or. istep == nstep) then
            write(*,'(a,i0,a,es12.4,a,2es12.4)') 'step ',istep, &
                ': particle_residual=',particle_solve_residual, &
                ', cont/gauss=',continuity_linf,gauss_linf
            write(*,'(a,6f10.4)') '  seconds rhs/particle/maxwell/rho/jx/ex =', &
                stage_time_global
            flush(output_unit)
        endif
    endif

    Up = Up_new
    Ufield = Ufield_new
    rho_old = rho_new
    if (snapshot_output_due(int(istep,8), nstep, dt) .or. &
        diagnostic_output_due(int(istep,8), nstep)) then
        call extract_physical_state(real(istep,8)*dt)
        if (snapshot_output_due(int(istep,8), nstep, dt)) &
            call write_snapshot_outputs(p, f, real(istep,8)*dt, prank, rho_total)
        if (diagnostic_output_due(int(istep,8), nstep)) &
            call write_diagnostic_outputs(p, f, real(istep,8)*dt, prank, rho_total)
    end if
enddo

if (output_final_solution) call write_final_solution()
call cpu_time(cpu_end)
elapsed = cpu_end-cpu_start
if (prank == 0) write(*,'(a,f12.3)') 'Total CPU time (s): ',elapsed

if (prank == 0 .and. len_trim(state_file) > 0) then
    open(unit=71,file=trim(state_file),status='replace',action='write')
    do itau = 0, Ntau-1
        do ix = 0, nx-1
            write(71,'(2i6,6es25.16)') itau,ix, &
                real(Ufield(1,itau,ix),8),aimag(Ufield(1,itau,ix)), &
                real(Ufield(2,itau,ix),8),aimag(Ufield(2,itau,ix)), &
                real(Ufield(3,itau,ix),8),aimag(Ufield(3,itau,ix))
        enddo
    enddo
    close(71)
endif

call dfftw_destroy_plan(plan_tau_forward)
call dfftw_destroy_plan(plan_tau_backward)
call dfftw_destroy_plan(plan_x_forward)
call dfftw_destroy_plan(plan_x_backward)
call free_poisson()
call MPI_FINALIZE(ierr)

contains

!-----------------------------------------------------------------------
subroutine compact_real_label(value,label)

real(8), intent(in) :: value
character(len=*), intent(out) :: label
character(len=64) :: buffer
integer :: last

write(buffer,'(f24.12)') value
buffer = adjustl(buffer)
last = len_trim(buffer)
do while (last > 1 .and. buffer(last:last) == '0')
    last = last-1
enddo
if (buffer(last:last) == '.') last = last-1
label = buffer(:last)

end subroutine compact_real_label

!-----------------------------------------------------------------------
subroutine allocate_field_arrays()

allocate(f%exhalf(0:nx));          f%exhalf = 0.d0
allocate(f%ey(0:nx));              f%ey = 0.d0
allocate(f%deydx(0:nx));           f%deydx = 0.d0
allocate(f%bz(0:nx));              f%bz = 0.d0
allocate(f%dbzdx(0:nx));           f%dbzdx = 0.d0
allocate(f%rho(0:nx));             f%rho = 0.d0
allocate(f%currentxhalf(0:nx));     f%currentxhalf = 0.d0
allocate(f%currenty(0:nx));         f%currenty = 0.d0
allocate(f%fvxvy(0:nvx,0:nvy));     f%fvxvy = 0.d0
allocate(f%fxvx(0:nx,0:nvx));       f%fxvx = 0.d0

end subroutine allocate_field_arrays

subroutine solve_gauss_all_tau(rho_tau, field_tau)

real(8), intent(in) :: rho_tau(0:Ntau-1,0:nx)
complex(8), intent(inout) :: field_tau(3,0:Ntau-1,0:nx)

do itau = 0, Ntau-1
    call solve_poisson_face(rho_tau(itau,:),nx,dx,work_face)
    field_tau(1,itau,:) = cmplx(work_face,0.d0,kind=8)
enddo

end subroutine solve_gauss_all_tau

!-----------------------------------------------------------------------
subroutine compute_second_order_field_source(particle_u, field_u, &
                                             ey_rhs_zero, bz_rhs_zero, source)

complex(8), intent(in) :: particle_u(3,0:Ntau-1,1:plocal)
complex(8), intent(in) :: field_u(3,0:Ntau-1,0:nx)
complex(8), intent(in) :: ey_rhs_zero(0:Ntau-1,0:nx)
complex(8), intent(in) :: bz_rhs_zero(0:Ntau-1,0:nx)
complex(8), intent(out) :: source(2,0:Ntau-1,0:nx)

real(8) :: mean_ey_rhs(0:nx), mean_bz_rhs(0:nx)
real(8) :: dfr_ey(0:nx), dfr_bz(0:nx)
real(8) :: dpr_ey_local(0:nx), dpr_ey_global(0:nx)
real(8) :: c, s

! D_F Pi S: the Maxwell part couples the averaged Ey and Bz right
! hand sides through their spatial derivatives.
mean_ey_rhs = sum(real(ey_rhs_zero,kind=8),dim=1)/real(Ntau,8)
mean_bz_rhs = sum(real(bz_rhs_zero,kind=8),dim=1)/real(Ntau,8)
call calcul_dx(nx,dimx,-mean_bz_rhs,dfr_ey)
call calcul_dx(nx,dimx,-mean_ey_rhs,dfr_bz)

do itau = 0, Ntau-1
    c = dcos(tau(itau))
    s = dsin(tau(itau))
    p%pos(:,1) = real(particle_u(1,itau,:),kind=8)
    p%vit(:,1) = real(c*particle_u(2,itau,:)+ &
                      s*particle_u(3,itau,:),kind=8)
    p%vit(:,2) = real(-s*particle_u(2,itau,:)+ &
                       c*particle_u(3,itau,:),kind=8)
    call wrap_particle_copy()

    f%exhalf = real(field_u(1,itau,:),kind=8)
    f%ey = real(field_u(2,itau,:),kind=8)
    f%bz = real(field_u(3,itau,:),kind=8)
    call interpol_m6(f,p)
    call deposit_dprpir_ey(tau(itau),dpr_ey_local)
    call allreduce_real(dpr_ey_local,dpr_ey_global)

    source(1,itau,:) = cmplx(dfr_ey+dpr_ey_global,0.d0,kind=8)
    source(2,itau,:) = cmplx(dfr_bz,0.d0,kind=8)
enddo

source(:,:,nx) = source(:,:,0)

end subroutine compute_second_order_field_source

!-----------------------------------------------------------------------
subroutine deposit_dprpir_ey(tau_now,h_ey)

real(8), intent(in) :: tau_now
real(8), intent(out) :: h_ey(0:nx)

integer(8) :: k, inode, inode_unwrapped
real(8) :: xp, q, coef, c, s

h_ey = 0.d0
c = dcos(tau_now)
s = dsin(tau_now)

do k = 1, p%nbpa
    ! This is the Ey component of D_P Pi S used by the existing TSI
    ! second-order initialization.  p%vit contains the physical,
    ! tau-rotated velocity and p%mag_z the interpolated Bz.
    coef = -p%w(k)*p%mag_z(k)* &
           (-c*p%vit(k,1)-s*p%vit(k,2))
    xp = p%pos(k,1)/dx
    do inode_unwrapped = floor(xp)-3, floor(xp)+3
        q = abs(real(inode_unwrapped,8)-xp)
        if (q < 3.d0) then
            inode = modulo(inode_unwrapped,nx)
            h_ey(inode) = h_ey(inode)+coef*f_m6(q)
        endif
    enddo
enddo

h_ey(0:nx-1) = h_ey(0:nx-1)/dx
h_ey(nx) = h_ey(0)

end subroutine deposit_dprpir_ey

!-----------------------------------------------------------------------
subroutine evaluate_explicit_rhs(particle_u, field_u, particle_rhs, ey_rhs, bz_rhs)

complex(8), intent(in) :: particle_u(3,0:Ntau-1,1:plocal)
complex(8), intent(in) :: field_u(3,0:Ntau-1,0:nx)
complex(8), intent(out) :: particle_rhs(3,0:Ntau-1,1:plocal)
complex(8), intent(out) :: ey_rhs(0:Ntau-1,0:nx)
complex(8), intent(out) :: bz_rhs(0:Ntau-1,0:nx)

real(8) :: c, s

do itau = 0, Ntau-1
    c = dcos(tau(itau))
    s = dsin(tau(itau))

    p%pos(:,1) = real(particle_u(1,itau,:),kind=8)
    p%vit(:,1) = real(c*particle_u(2,itau,:)+s*particle_u(3,itau,:),kind=8)
    p%vit(:,2) = real(-s*particle_u(2,itau,:)+c*particle_u(3,itau,:),kind=8)
    call wrap_particle_copy()

    f%exhalf = real(field_u(1,itau,:),kind=8)
    f%ey = real(field_u(2,itau,:),kind=8)
    f%bz = real(field_u(3,itau,:),kind=8)
    call interpol_m6(f,p)
    call calcul_current(p,f)
    call allreduce_real(f%currenty,currenty_global)
    jy_old(itau,:) = currenty_global

    call calcul_dx(nx,dimx,f%ey,f%deydx)
    call calcul_dx(nx,dimx,f%bz,f%dbzdx)

    do ip = 1, plocal
        particle_rhs(1,itau,ip) = c*particle_u(2,itau,ip) + s*particle_u(3,itau,ip)
        particle_rhs(2,itau,ip) = c*p%ele_x(ip) - s*p%ele_y(ip) + p%mag_z(ip)*particle_u(3,itau,ip)
        particle_rhs(3,itau,ip) = s*p%ele_x(ip) + c*p%ele_y(ip) - p%mag_z(ip)*particle_u(2,itau,ip)
    enddo

    ey_rhs(itau,:) = cmplx(-f%dbzdx-currenty_global,0.d0,kind=8)
    bz_rhs(itau,:) = cmplx(-f%deydx,0.d0,kind=8)
enddo

end subroutine evaluate_explicit_rhs

!-----------------------------------------------------------------------
subroutine build_dtau_matrix(dmat)

complex(8), intent(out) :: dmat(0:Ntau-1,0:Ntau-1)
integer :: row, col

do col = 0, Ntau-1
    fftin = cmplx(0.d0,0.d0,kind=8)
    fftin(col) = cmplx(1.d0,0.d0,kind=8)
    call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
    do row = 0, Ntau-1
        fftout(row) = cmplx(0.d0,ktau(row),kind=8)*fftout(row)
    enddo
    call dfftw_execute_dft(plan_tau_backward,fftout,fftin)
    dmat(:,col) = fftin/real(Ntau,8)
enddo

end subroutine build_dtau_matrix

!-----------------------------------------------------------------------
subroutine particle_imex_update(particle_old, field_old, particle_new)

complex(8), intent(in) :: particle_old(3,0:Ntau-1,1:plocal)
complex(8), intent(in) :: field_old(3,0:Ntau-1,0:nx)
complex(8), intent(out) :: particle_new(3,0:Ntau-1,1:plocal)

real(8), allocatable :: e1(:,:), e2(:,:), bcoef(:,:)
complex(8) :: amat(0:2*Ntau-1,0:2*Ntau-1), amat0(0:2*Ntau-1,0:2*Ntau-1)
complex(8) :: vrhs(0:2*Ntau-1), vrhs0(0:2*Ntau-1), vnew(0:2*Ntau-1)
complex(8) :: rhsx(0:Ntau-1)
real(8) :: c, s, local_error, global_error
integer :: info

allocate(e1(0:Ntau-1,1:plocal),e2(0:Ntau-1,1:plocal), &
         bcoef(0:Ntau-1,1:plocal))

! Freeze the interpolation coefficients at time level n.
do itau = 0, Ntau-1
    c = dcos(tau(itau)); s = dsin(tau(itau))
    p%pos(:,1) = real(particle_old(1,itau,:),kind=8)
    call wrap_particle_copy()
    f%exhalf = real(field_old(1,itau,:),kind=8)
    f%ey = real(field_old(2,itau,:),kind=8)
    f%bz = real(field_old(3,itau,:),kind=8)
    call interpol_m6(f,p)
    e1(itau,:) = c*p%ele_x(:)-s*p%ele_y(:)
    e2(itau,:) = s*p%ele_x(:)+c*p%ele_y(:)
    bcoef(itau,:) = p%mag_z(:)
enddo

local_error = 0.d0
do ip = 1, plocal
    ! Solve the full two-component complex system.  Keeping both blocks
    ! also preserves any tiny imaginary parts introduced by FFT roundoff.
    amat = cmplx(0.d0,0.d0,kind=8)
    amat(0:Ntau-1,0:Ntau-1) = (dt/ep)*dtau_matrix
    amat(Ntau:2*Ntau-1,Ntau:2*Ntau-1) = (dt/ep)*dtau_matrix
    do itau = 0, Ntau-1
        amat(itau,itau) = amat(itau,itau)+cmplx(1.d0,0.d0,kind=8)
        amat(itau+Ntau,itau+Ntau) = amat(itau+Ntau,itau+Ntau)+cmplx(1.d0,0.d0,kind=8)
        amat(itau,itau+Ntau) = cmplx(-dt*bcoef(itau,ip),0.d0,kind=8)
        amat(itau+Ntau,itau) = cmplx(dt*bcoef(itau,ip),0.d0,kind=8)
    enddo
    vrhs(0:Ntau-1) = particle_old(2,:,ip)+dt*cmplx(e1(:,ip),0.d0,kind=8)
    vrhs(Ntau:2*Ntau-1) = particle_old(3,:,ip)+dt*cmplx(e2(:,ip),0.d0,kind=8)
    amat0 = amat
    vrhs0 = vrhs
    call solve_complex_dense(amat,vrhs,vnew,info)
    if (info /= 0) then
        write(*,'(a,i0,a,i0)') 'ERROR: singular direct particle matrix on rank ',prank,', particle ',ip
        call MPI_ABORT(MPI_COMM_WORLD,info,ierr)
    endif
    local_error = max(local_error,maxval(abs(matmul(amat0,vnew)-vrhs0)))
    particle_new(2,:,ip) = vnew(0:Ntau-1)
    particle_new(3,:,ip) = vnew(Ntau:2*Ntau-1)
    do itau = 0, Ntau-1
        c = dcos(tau(itau)); s = dsin(tau(itau))
        rhsx(itau) = c*particle_new(2,itau,ip)+s*particle_new(3,itau,ip)
    enddo
    call be_tau_update(particle_old(1,:,ip),rhsx,particle_new(1,:,ip))
enddo
call MPI_ALLREDUCE(local_error,global_error,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
particle_solve_residual = global_error
deallocate(e1,e2,bcoef)

end subroutine particle_imex_update

!-----------------------------------------------------------------------
subroutine solve_complex_dense(a,b,x,info)

complex(8), intent(inout) :: a(0:,0:)
complex(8), intent(inout) :: b(0:)
complex(8), intent(out) :: x(0:)
integer, intent(out) :: info
complex(8) :: rowtmp(0:size(b)-1), factor, btmp
real(8) :: pivot_size
integer :: i, j, k, pivot, nsys

info = 0
nsys = size(b)
do k = 0, nsys-2
    pivot = k
    pivot_size = abs(a(k,k))
    do i = k+1, nsys-1
        if (abs(a(i,k)) > pivot_size) then
            pivot = i
            pivot_size = abs(a(i,k))
        endif
    enddo
    if (pivot_size <= tiny(1.d0)) then
        info = k+1
        return
    endif
    if (pivot /= k) then
        rowtmp = a(k,:); a(k,:) = a(pivot,:); a(pivot,:) = rowtmp
        btmp = b(k); b(k) = b(pivot); b(pivot) = btmp
    endif
    do i = k+1, nsys-1
        factor = a(i,k)/a(k,k)
        a(i,k) = cmplx(0.d0,0.d0,kind=8)
        do j = k+1, nsys-1
            a(i,j) = a(i,j)-factor*a(k,j)
        enddo
        b(i) = b(i)-factor*b(k)
    enddo
enddo
if (abs(a(nsys-1,nsys-1)) <= tiny(1.d0)) then
    info = nsys
    return
endif
x = cmplx(0.d0,0.d0,kind=8)
do i = nsys-1, 0, -1
    x(i) = b(i)
    do j = i+1, nsys-1
        x(i) = x(i)-a(i,j)*x(j)
    enddo
    x(i) = x(i)/a(i,i)
enddo

end subroutine solve_complex_dense

!-----------------------------------------------------------------------
subroutine maxwell_imex_update(field_old, current_old, field_new)

complex(8), intent(in) :: field_old(3,0:Ntau-1,0:nx)
real(8), intent(in) :: current_old(0:Ntau-1,0:nx)
complex(8), intent(inout) :: field_new(3,0:Ntau-1,0:nx)

complex(8), allocatable :: eyhat(:,:), bzhat(:,:), rhsey(:,:), rhsbz(:,:)
complex(8) :: a, offdiag, determinant, rhs_e, rhs_b
real(8) :: kx
integer :: ik

allocate(eyhat(0:Ntau-1,0:nx-1),bzhat(0:Ntau-1,0:nx-1))
allocate(rhsey(0:Ntau-1,0:nx-1),rhsbz(0:Ntau-1,0:nx-1))

! Transform the right-hand sides first in tau and then in x.
do ix = 0, nx-1
    fftin = field_old(2,:,ix)-dt*cmplx(current_old(:,ix),0.d0,kind=8)
    call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
    rhsey(:,ix) = fftout
    fftin = field_old(3,:,ix)
    call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
    rhsbz(:,ix) = fftout
enddo
do itau = 0, Ntau-1
    xfftin = rhsey(itau,:)
    call dfftw_execute_dft(plan_x_forward,xfftin,xfftout)
    rhsey(itau,:) = xfftout
    xfftin = rhsbz(itau,:)
    call dfftw_execute_dft(plan_x_forward,xfftin,xfftout)
    rhsbz(itau,:) = xfftout
enddo

! [a, i*dt*k; i*dt*k, a] [Ey;Bz] = [Ey^n-dt*Jy^n;Bz^n].
! This is backward Euler for both tau/eps and the Maxwell wave operator.
do itau = 0, Ntau-1
    a = 1.d0+dt*cmplx(0.d0,ktau(itau),kind=8)/ep
    do ik = 0, nx-1
        if (mod(nx,2_8) == 0 .and. ik == nx/2) then
            ! Match calcul_dx: the real Nyquist derivative is discarded.
            kx = 0.d0
        else if (ik <= nx/2) then
            kx = 2.d0*pi*real(ik,8)/dimx
        else
            kx = 2.d0*pi*real(ik-nx,8)/dimx
        endif
        offdiag = cmplx(0.d0,dt*kx,kind=8)
        determinant = a*a-offdiag*offdiag
        rhs_e = rhsey(itau,ik); rhs_b = rhsbz(itau,ik)
        eyhat(itau,ik) = (a*rhs_e-offdiag*rhs_b)/determinant
        bzhat(itau,ik) = (a*rhs_b-offdiag*rhs_e)/determinant
    enddo
enddo

! Inverse x and tau transforms.
do itau = 0, Ntau-1
    xfftin = eyhat(itau,:)
    call dfftw_execute_dft(plan_x_backward,xfftin,xfftout)
    eyhat(itau,:) = xfftout/real(nx,8)
    xfftin = bzhat(itau,:)
    call dfftw_execute_dft(plan_x_backward,xfftin,xfftout)
    bzhat(itau,:) = xfftout/real(nx,8)
enddo
do ix = 0, nx-1
    fftin = eyhat(:,ix)
    call dfftw_execute_dft(plan_tau_backward,fftin,fftout)
    field_new(2,:,ix) = fftout/real(Ntau,8)
    fftin = bzhat(:,ix)
    call dfftw_execute_dft(plan_tau_backward,fftin,fftout)
    field_new(3,:,ix) = fftout/real(Ntau,8)
enddo
field_new(2:3,:,nx) = field_new(2:3,:,0)

deallocate(eyhat,bzhat,rhsey,rhsbz)

end subroutine maxwell_imex_update

!-----------------------------------------------------------------------
subroutine be_tau_update(old_value, explicit_rhs, new_value)

complex(8), intent(in) :: old_value(0:Ntau-1)
complex(8), intent(in) :: explicit_rhs(0:Ntau-1)
complex(8), intent(out) :: new_value(0:Ntau-1)

do itau = 0, Ntau-1
    fftin(itau) = old_value(itau) + dt*explicit_rhs(itau)
enddo
call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
do itau = 0, Ntau-1
    fftout(itau) = fftout(itau)/(1.d0+dt*cmplx(0.d0,ktau(itau),kind=8)/ep)
enddo
call dfftw_execute_dft(plan_tau_backward,fftout,fftin)
new_value = fftin/real(Ntau,8)

end subroutine be_tau_update


!-----------------------------------------------------------------------
subroutine tau_zero_mean_antiderivative(value, anti)

complex(8), intent(in) :: value(0:Ntau-1)
complex(8), intent(out) :: anti(0:Ntau-1)

fftin = value
call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
fftout(0) = cmplx(0.d0,0.d0,kind=8)
do itau = 1, Ntau-1
    if (ktau(itau) == 0.d0) then
        fftout(itau) = cmplx(0.d0,0.d0,kind=8)
    else
        fftout(itau) = fftout(itau)/ &
                       cmplx(0.d0,ktau(itau),kind=8)
    endif
enddo
call dfftw_execute_dft(plan_tau_backward,fftout,fftin)
anti = fftin/real(Ntau,8)

end subroutine tau_zero_mean_antiderivative

!-----------------------------------------------------------------------
subroutine deposit_density_all_tau(particle_u, rho_tau)

complex(8), intent(in) :: particle_u(3,0:Ntau-1,1:plocal)
real(8), intent(out) :: rho_tau(0:Ntau-1,0:nx)

do itau = 0, Ntau-1
    p%pos(:,1) = real(particle_u(1,itau,:),kind=8)
    call wrap_particle_copy()
    call calcul_rho_m6(p,f)
    call allreduce_real(f%rho,rho_tau(itau,:))
enddo

end subroutine deposit_density_all_tau

!-----------------------------------------------------------------------
subroutine deposit_charge_preserving_jx(x_old, x_new, jx_global)

complex(8), intent(in) :: x_old(0:Ntau-1,1:plocal)
complex(8), intent(in) :: x_new(0:Ntau-1,1:plocal)
real(8), intent(out) :: jx_global(0:Ntau-1,0:nx)

real(8), allocatable :: jx_local(:,:)
real(8) :: primitive_tau(0:Ntau-1), dprimitive_tau(0:Ntau-1)
real(8) :: xa, xb, xface, qold, qnew
real(8) :: weight, jt
integer(8) :: iface_unwrapped, iface_min, iface_max, iface

allocate(jx_local(0:Ntau-1,0:nx))
jx_local = 0.d0

do ip = 1, plocal
    weight = p%w(ip)
    do itau = 0, Ntau-1
        xa = real(x_old(itau,ip),kind=8)
        xb = real(x_new(itau,ip),kind=8)

        ! Only faces intersecting the swept M5 support can contribute.
        iface_min = floor(min(xa,xb)/dx-3.d0)
        iface_max = ceiling(max(xa,xb)/dx+3.d0)
        do iface_unwrapped = iface_min, iface_max
            xface = (real(iface_unwrapped,8)+0.5d0)*dx
            qold = (xface-xa)/dx
            qnew = (xface-xb)/dx
            jt = weight/dt*(primitive_m5_stable(qold)-primitive_m5_stable(qnew))
            iface = modulo(iface_unwrapped,nx)
            jx_local(itau,iface) = jx_local(itau,iface) + jt
        enddo

    enddo
enddo

! Construct the tau-current particle by particle.  For a fixed face set
!
!   P_p(tau) = primitive_m5((x_face-X_p^{n+1}(tau))/dx),
!   J_tau,p  = -w_p/ep * D_tau P_p.
!
! D_tau is exactly the FFT differentiation used for rho in diagnostics.
! Since the M6/M5 identity gives D_x(-w_p P_p)=-rho_p, the two linear
! discrete derivatives commute and therefore
!
!   D_x J_tau,p = -(D_tau rho_p^{n+1})/ep
!
! without invoking a continuous chain rule at the tau collocation nodes.
do ip = 1, plocal
    weight = p%w(ip)
    iface_min = floor(minval(real(x_new(:,ip),kind=8))/dx-3.d0)
    iface_max = ceiling(maxval(real(x_new(:,ip),kind=8))/dx+3.d0)
    do iface_unwrapped = iface_min, iface_max
        xface = (real(iface_unwrapped,8)+0.5d0)*dx
        do itau = 0, Ntau-1
            qnew = (xface-real(x_new(itau,ip),kind=8))/dx
            primitive_tau(itau) = primitive_m5_stable(qnew)
        enddo
        call spectral_tau_derivative_real(primitive_tau,dprimitive_tau)
        iface = modulo(iface_unwrapped,nx)
        jx_local(:,iface) = jx_local(:,iface) - &
            weight*dprimitive_tau(:)/ep
    enddo
enddo

! Reduce J_x=J_t+J_tau only after both particlewise contributions have
! been assembled, so the current is not multiplied by the MPI rank count.
do itau = 0, Ntau-1
    jx_local(itau,nx) = jx_local(itau,0)
    call allreduce_real(jx_local(itau,:),jx_global(itau,:))
enddo

deallocate(jx_local)

end subroutine deposit_charge_preserving_jx


!-----------------------------------------------------------------------
real(8) function primitive_m5(q)

real(8), intent(in) :: q
real(8), parameter :: coef(0:5) = &
    (/1.d0,-5.d0,10.d0,-10.d0,5.d0,-1.d0/)
real(8) :: z
integer :: m

primitive_m5 = 0.d0
do m = 0, 5
    z = q+2.5d0-real(m,8)
    if (z > 0.d0) primitive_m5 = primitive_m5+coef(m)*z**5
enddo
primitive_m5 = primitive_m5/120.d0

end function primitive_m5

!-----------------------------------------------------------------------
real(8) function primitive_m5_stable(q)

real(8), intent(in) :: q

! Outside the compact M5 support the cumulative shape is exactly 0 or 1.
! Returning these constants avoids cancellation in the alternating-power
! representation and makes their tau derivatives exactly zero.
if (q <= -2.5d0) then
    primitive_m5_stable = 0.d0
else if (q >= 2.5d0) then
    primitive_m5_stable = 1.d0
else
    primitive_m5_stable = primitive_m5(q)
endif

end function primitive_m5_stable

!-----------------------------------------------------------------------
subroutine update_exhalf(ex_old, jx, ex_new)

complex(8), intent(in) :: ex_old(0:Ntau-1,0:nx)
real(8), intent(in) :: jx(0:Ntau-1,0:nx)
complex(8), intent(out) :: ex_new(0:Ntau-1,0:nx)

complex(8) :: rhs_tau(0:Ntau-1)
real(8) :: mean_current

do ix = 0, nx-1
    do itau = 0, Ntau-1
        mean_current = sum(jx(itau,0:nx-1))/real(nx,8)
        rhs_tau(itau) = cmplx(-(jx(itau,ix)-mean_current),0.d0,kind=8)
    enddo
    call be_tau_update(ex_old(:,ix),rhs_tau,ex_new(:,ix))
enddo
ex_new(:,nx) = ex_new(:,0)

end subroutine update_exhalf


!-----------------------------------------------------------------------
subroutine compute_diagnostics(rho_n, rho_np1, jx, field_np1, cont_norm, gauss_norm, initial)

real(8), intent(in) :: rho_n(0:Ntau-1,0:nx)
real(8), intent(in) :: rho_np1(0:Ntau-1,0:nx)
real(8), intent(in) :: jx(0:Ntau-1,0:nx)
complex(8), intent(in) :: field_np1(3,0:Ntau-1,0:nx)
real(8), intent(out) :: cont_norm, gauss_norm
logical, intent(in) :: initial

real(8) :: local_cont, local_gauss
real(8) :: drhodtau(0:Ntau-1)
integer(8) :: im1

local_cont = 0.d0
local_gauss = 0.d0

do ix = 0, nx-1
    call spectral_tau_derivative_real(rho_np1(:,ix),drhodtau)
    do itau = 0, Ntau-1
        im1 = modulo(int(ix,8)-1_8,nx)
        if (.not.initial) then
            local_cont = max(local_cont,abs( &
                (rho_np1(itau,ix)-rho_n(itau,ix))/dt + &
                drhodtau(itau)/ep + &
                (jx(itau,ix)-jx(itau,im1))/dx))
        endif
        local_gauss = max(local_gauss,abs( &
            (real(field_np1(1,itau,ix),kind=8)- &
             real(field_np1(1,itau,im1),kind=8))/dx - &
            rho_np1(itau,ix)+ni))
    enddo
enddo

call MPI_ALLREDUCE(local_cont,cont_norm,1,MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,ierr)
call MPI_ALLREDUCE(local_gauss,gauss_norm,1,MPI_REAL8,MPI_MAX, MPI_COMM_WORLD,ierr)

end subroutine compute_diagnostics

!-----------------------------------------------------------------------
subroutine spectral_tau_derivative_real(value, derivative)

real(8), intent(in) :: value(0:Ntau-1)
real(8), intent(out) :: derivative(0:Ntau-1)

fftin = cmplx(value,0.d0,kind=8)
call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
do itau = 0, Ntau-1
    fftout(itau) = cmplx(0.d0,ktau(itau),kind=8)*fftout(itau)
enddo
call dfftw_execute_dft(plan_tau_backward,fftout,fftin)
derivative = real(fftin,kind=8)/real(Ntau,8)

end subroutine spectral_tau_derivative_real

!-----------------------------------------------------------------------
subroutine extract_physical_state(time_now)

real(8), intent(in) :: time_now
complex(8) :: particle_value(3), field_value(3), phase_factor
real(8) :: phase, c, s

phase = modulo(time_now/ep,2.d0*pi)
c = dcos(phase)
s = dsin(phase)

do ip = 1, plocal
    do icomp = 1, 3
        fftin = Up(icomp,:,ip)
        call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
        particle_value(icomp) = cmplx(0.d0,0.d0,kind=8)
        do itau = 0, Ntau-1
            phase_factor = exp(cmplx(0.d0,ktau(itau)*phase,kind=8))
            particle_value(icomp) = particle_value(icomp)+fftout(itau)*phase_factor/real(Ntau,8)
        enddo
    enddo
    p%pos(ip,1) = modulo(real(particle_value(1),kind=8),dimx)
    p%vit(ip,1) = real(c*particle_value(2)+s*particle_value(3),kind=8)
    p%vit(ip,2) = real(-s*particle_value(2)+c*particle_value(3),kind=8)
enddo

do ix = 0, nx
    do icomp = 1, 3
        fftin = Ufield(icomp,:,ix)
        call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
        field_value(icomp) = cmplx(0.d0,0.d0,kind=8)
        do itau = 0, Ntau-1
            phase_factor = exp(cmplx(0.d0,ktau(itau)*phase,kind=8))
            field_value(icomp) = field_value(icomp)+fftout(itau)*phase_factor/real(Ntau,8)
        enddo
    enddo
    f%exhalf(ix) = real(field_value(1),kind=8)
    f%ey(ix) = real(field_value(2),kind=8)
    f%bz(ix) = real(field_value(3),kind=8)
enddo

do ix = 0, nx
    fftin = cmplx(rho_old(:,ix),0.d0,kind=8)
    call dfftw_execute_dft(plan_tau_forward,fftin,fftout)
    rho_total(ix) = 0.d0
    do itau = 0, Ntau-1
        phase_factor = exp(cmplx(0.d0,ktau(itau)*phase,kind=8))
        rho_total(ix) = rho_total(ix)+real(fftout(itau)*phase_factor,kind=8)/real(Ntau,8)
    enddo
enddo

end subroutine extract_physical_state

!-----------------------------------------------------------------------
subroutine write_observables(time_now,replace_files)

real(8), intent(in) :: time_now
logical, intent(in) :: replace_files
real(8) :: poisson_tau, poisson_slice, energy_now, rho_mean
real(8) :: kinetic_global(0:nx), residual_sum
integer(8) :: im1

poisson_tau = 0.d0
do itau = 0, Ntau-1
    rho_mean = sum(rho_old(itau,0:nx-1))*dx/dimx
    residual_sum = 0.d0
    do ix = 0, nx-1
        im1 = modulo(int(ix,8)-1_8,nx)
        residual_sum = residual_sum+ &
            ((real(Ufield(1,itau,ix),kind=8)- &
              real(Ufield(1,itau,im1),kind=8))/dx- &
              rho_old(itau,ix)+rho_mean)**2
    enddo
    poisson_tau = max(poisson_tau,sqrt(residual_sum*dx))
enddo

call extract_physical_state(time_now)
rho_mean = sum(rho_total(0:nx-1))*dx/dimx
residual_sum = 0.d0
do ix = 0, nx-1
    im1 = modulo(int(ix,8)-1_8,nx)
    residual_sum = residual_sum+((f%exhalf(ix)-f%exhalf(im1))/dx-rho_total(ix)+rho_mean)**2
enddo
poisson_slice = sqrt(residual_sum*dx)

call calcul_energy(p,f)
call allreduce_real(f%rho,kinetic_global)
energy_now = sum(kinetic_global(0:nx-1))*dx/2.d0
do ix = 0, nx-1
    energy_now = energy_now+(f%exhalf(ix)**2+f%ey(ix)**2+f%bz(ix)**2)*dx/2.d0
enddo

if (prank == 0) then
    if (replace_files) then
        open(unit=81,file=trim(poisson_tau_file),status='replace',action='write')
        open(unit=83,file=trim(poisson_file),status='replace',action='write')
        open(unit=84,file=trim(energy_file),status='replace',action='write')
    else
        open(unit=81,file=trim(poisson_tau_file),status='old',action='write',position='append')
        open(unit=83,file=trim(poisson_file),status='old',action='write',position='append')
        open(unit=84,file=trim(energy_file),status='old',action='write',position='append')
    endif
    write(81,*) poisson_tau
    write(83,*) poisson_slice
    write(84,*) energy_now
    close(81)
    close(83)
    close(84)
endif

end subroutine write_observables

!-----------------------------------------------------------------------
subroutine write_final_solution()

real(8) :: rhov_global(0:nx)
character(len=256) :: filename
integer :: file_unit

call extract_physical_state(tfinal)

if (prank == 0) then
    filename = 'rho_'//trim(diag_prefix)//'.dat'
    open(newunit=file_unit,file=trim(filename),status='replace',action='write')
    do ix = 0, nx
        write(file_unit,*) rho_total(ix)
    enddo
    close(file_unit)

    filename = 'ExEy_'//trim(diag_prefix)//'.dat'
    open(newunit=file_unit,file=trim(filename),status='replace',action='write')
    do ix = 0, nx
        write(file_unit,*) f%exhalf(ix),f%ey(ix)
    enddo
    close(file_unit)

    filename = 'Bz_'//trim(diag_prefix)//'.dat'
    open(newunit=file_unit,file=trim(filename),status='replace',action='write')
    do ix = 0, nx
        write(file_unit,*) f%bz(ix)
    enddo
    close(file_unit)
endif

call calcul_energy(p,f)
call allreduce_real(f%rho,rhov_global)
if (prank == 0) then
    filename = 'rhov_'//trim(diag_prefix)//'.dat'
    open(newunit=file_unit,file=trim(filename),status='replace',action='write')
    do ix = 0, nx
        write(file_unit,*) rhov_global(ix)
    enddo
    close(file_unit)
endif

end subroutine write_final_solution

!-----------------------------------------------------------------------
subroutine allreduce_real(local_value,global_value)

real(8), intent(in) :: local_value(0:nx)
real(8), intent(out) :: global_value(0:nx)

call MPI_ALLREDUCE(local_value,global_value,int(nx+1),MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,ierr)

end subroutine allreduce_real

!-----------------------------------------------------------------------
subroutine wrap_particle_copy()

do ip = 1, p%nbpa
    p%pos(ip,1) = modulo(p%pos(ip,1),dimx)
enddo

end subroutine wrap_particle_copy

!-----------------------------------------------------------------------
subroutine calcul_ey(nx_in, dimx_in, ey)

integer(8), intent(in) :: nx_in
real(8), intent(in) :: dimx_in
real(8), intent(out) :: ey(0:nx_in)
integer(8) :: i

do i = 1, nx_in
    ey(i) = 0.d0
enddo
ey(0) = ey(nx_in)

end subroutine calcul_ey

!-----------------------------------------------------------------------
subroutine calcul_bz(nx_in, dimx_in, bz)

integer(8), intent(in) :: nx_in
real(8), intent(in) :: dimx_in
real(8), intent(out) :: bz(0:nx_in)
integer(8) :: i
real(8) :: x, dx_local, beta, k

dx_local = dimx_in/real(nx_in,8)
k = cfg_wave_number
select case (initial_condition)
case (initial_original_two_stream)
    beta = cfg_bz_two_stream
case (initial_anisotropic_weibel)
    beta = cfg_bz_weibel
case (initial_bump_on_tail)
    beta = 0.d0
end select
do i = 1, nx_in
    x = real(i,8)*dx_local
    bz(i) = beta*dsin(k*x)
enddo
bz(0) = bz(nx_in)

end subroutine calcul_bz

end program vlasov_maxwell_pic_1dx2dv
