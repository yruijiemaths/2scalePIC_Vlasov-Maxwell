module pic_output
  use mesh, only: field, nx, nvx, nvy, dx, mesh_ep => ep, mesh_dt => dt
  use interpolations, only: particle, calcul_fvxvy_m6, calcul_fxvx_m6, &
                            calcul_rho_m6, calcul_energy
  use grid_output_adapter, only: write_field_snapshot, compute_field_energy, &
                                 compute_field_l2_norms, compute_poisson_residual
  use pic_config, only: cfg_output_fxvx, cfg_output_fvxvy, &
                        cfg_output_rho, cfg_output_rhov, cfg_output_fields, &
                        cfg_output_mass, cfg_output_charge, cfg_output_energy, cfg_output_poisson, &
                        cfg_output_kinetic_components, cfg_output_field_l2, &
                        cfg_output_poisson_tau, cfg_output_momentum, &
                        cfg_output_initial, cfg_output_final, &
                        cfg_snapshot_interval, cfg_diagnostic_interval_steps, &
                        cfg_output_prefix
  use mpi
  implicit none
  private

  logical, save :: scalar_files_started = .false.
  logical, save :: momentum_directory_ready = .false.

  public :: snapshot_output_due, diagnostic_output_due, poisson_tau_output_due, &
            write_snapshot_outputs, write_diagnostic_outputs, write_poisson_tau

contains

  logical function snapshot_output_due(step, nstep, dt)
    integer(8), intent(in) :: step, nstep
    real(8), intent(in) :: dt
    integer(8) :: stride

    snapshot_output_due = .false.
    if (.not. snapshot_output_enabled()) return
    if (step == 0_8) then
      snapshot_output_due = cfg_output_initial
      return
    end if

    stride = max(1_8, nint(cfg_snapshot_interval/dt, kind=8))
    snapshot_output_due = mod(step, stride) == 0_8
    if (cfg_output_final .and. step == nstep) snapshot_output_due = .true.
  end function snapshot_output_due

  logical function diagnostic_output_due(step, nstep)
    integer(8), intent(in) :: step, nstep
    diagnostic_output_due = .false.
    if (.not. diagnostic_output_enabled()) return
    if (step == 0_8) then
      diagnostic_output_due = cfg_output_initial
      return
    end if
    diagnostic_output_due = mod(step, int(cfg_diagnostic_interval_steps,8)) == 0_8
    if (cfg_output_final .and. step == nstep) diagnostic_output_due = .true.
  end function diagnostic_output_due

  logical function poisson_tau_output_due(step, nstep, dt)
    integer(8), intent(in) :: step, nstep
    real(8), intent(in) :: dt
    poisson_tau_output_due = cfg_output_poisson_tau .and. diagnostic_output_due(step, nstep)
  end function poisson_tau_output_due

  subroutine write_poisson_tau(time_now, value, rank)
    real(8), intent(in) :: time_now, value
    integer, intent(in) :: rank
    if (rank == 0 .and. cfg_output_poisson_tau) &
      call append_scalar('Poissonres_tau', time_now, value)
  end subroutine write_poisson_tau

  subroutine write_diagnostic_outputs(p, f, time_now, rank, rho_override)
    type(particle), intent(inout) :: p
    type(field), intent(inout) :: f
    real(8), intent(in) :: time_now
    integer, intent(in) :: rank
    real(8), intent(in), optional :: rho_override(0:nx)
    real(8), allocatable :: rho_global(:), rhov_global(:)
    integer :: ierr
    real(8) :: mass_value, energy_value, field_energy, poisson_value
    real(8) :: local_particle_values(2), particle_values(2)
    real(8) :: local_momentum(2), momentum(2)
    real(8) :: kinetic_1, kinetic_2
    real(8) :: field_l2(3)
    logical :: need_rho

    need_rho = cfg_output_mass .or. cfg_output_charge .or. &
      cfg_output_poisson
    if (need_rho) then
      allocate(rho_global(0:nx))
      if (present(rho_override)) then
        rho_global = rho_override
      else
        call calcul_rho_m6(p, f)
        call MPI_REDUCE(f%rho, rho_global, size(rho_global), MPI_REAL8, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      end if
    end if

    if (cfg_output_energy) then
      allocate(rhov_global(0:nx))
      call calcul_energy(p, f)
      call MPI_REDUCE(f%rho, rhov_global, size(rhov_global), MPI_REAL8, &
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    end if

    if (rank == 0 .and. cfg_output_mass) then
      mass_value = sum(rho_global(0:nx-1))*dx
      call append_scalar('mass', time_now, mass_value)
    end if

    if (rank == 0 .and. cfg_output_charge) then
      mass_value = sum(rho_global(0:nx-1))*dx
      call append_scalar('charge', time_now, mass_value)
    end if

    if (cfg_output_momentum) then
      local_momentum(1) = sum(p%w*p%vit(:,1))
      local_momentum(2) = sum(p%w*p%vit(:,2))
      call MPI_REDUCE(local_momentum, momentum, 2, MPI_REAL8, MPI_SUM, 0, &
                      MPI_COMM_WORLD, ierr)
      if (rank == 0) call append_momentum(momentum)
    end if

    if (rank == 0 .and. cfg_output_energy) then
      call compute_field_energy(f, field_energy)
      energy_value = 0.5d0*sum(rhov_global(0:nx-1))*dx + field_energy
      call append_scalar('energy', time_now, energy_value)
    end if

    if (cfg_output_kinetic_components) then
      local_particle_values(1) = sum(p%w*p%vit(:,1)**2)
      local_particle_values(2) = sum(p%w*p%vit(:,2)**2)
      call MPI_REDUCE(local_particle_values, particle_values, 2, MPI_REAL8, &
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (rank == 0) then
        kinetic_1 = 0.5d0*particle_values(1)
        kinetic_2 = 0.5d0*particle_values(2)
        call append_kinetic_components(kinetic_1, kinetic_2)
      end if
    end if

    if (rank == 0 .and. cfg_output_field_l2) then
      call compute_field_l2_norms(f, field_l2)
      call append_field_l2(field_l2)
    end if

    if (rank == 0 .and. cfg_output_poisson) then
      call compute_poisson_residual(f, rho_global, poisson_value)
      call append_scalar('Poissonres', time_now, poisson_value)
    end if

    if (allocated(rho_global)) deallocate(rho_global)
    if (allocated(rhov_global)) deallocate(rhov_global)
    scalar_files_started = .true.
  end subroutine write_diagnostic_outputs

  subroutine write_snapshot_outputs(p, f, time_now, rank, rho_override)
    type(particle), intent(inout) :: p
    type(field), intent(inout) :: f
    real(8), intent(in) :: time_now
    integer, intent(in) :: rank
    real(8), intent(in), optional :: rho_override(0:nx)
    real(8), allocatable :: global_values(:,:), line_global(:)
    character(len=512) :: filename
    integer :: ierr, unit, i, j

    if (cfg_output_rho) then
      allocate(line_global(0:nx))
      if (present(rho_override)) then
        line_global = rho_override
      else
        call calcul_rho_m6(p, f)
        call MPI_REDUCE(f%rho, line_global, size(line_global), MPI_REAL8, &
                        MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      end if
      if (rank == 0) then
        call output_filename('rho', time_now, filename)
        open(newunit=unit, file=trim(filename), status='replace', action='write')
        do i = 0, nx
          write(unit,*) line_global(i)
        end do
        close(unit)
      end if
      deallocate(line_global)
    end if

    if (cfg_output_rhov) then
      allocate(line_global(0:nx))
      call calcul_energy(p, f)
      call MPI_REDUCE(f%rho, line_global, size(line_global), MPI_REAL8, &
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (rank == 0) then
        call output_filename('rhov', time_now, filename)
        open(newunit=unit, file=trim(filename), status='replace', action='write')
        do i = 0, nx
          write(unit,*) line_global(i)
        end do
        close(unit)
      end if
      deallocate(line_global)
    end if

    if (rank == 0 .and. cfg_output_fields) then
      call output_filename('ExEyBz', time_now, filename)
      open(newunit=unit, file=trim(filename), status='replace', action='write')
      call write_field_snapshot(f, unit)
      close(unit)
    end if

    if (cfg_output_fvxvy) then
      allocate(global_values(0:nvx,0:nvy))
      call calcul_fvxvy_m6(p, f)
      call MPI_REDUCE(f%fvxvy, global_values, size(global_values), MPI_REAL8, &
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (rank == 0) then
        call output_filename('fvxvy', time_now, filename)
        open(newunit=unit, file=trim(filename), status='replace', action='write')
        do i = 0, nvx
          do j = 0, nvy
            write(unit,*) global_values(i,j)
          end do
        end do
        close(unit)
      end if
      deallocate(global_values)
    end if

    if (cfg_output_fxvx) then
      allocate(global_values(0:nx,0:nvx))
      call calcul_fxvx_m6(p, f)
      call MPI_REDUCE(f%fxvx, global_values, size(global_values), MPI_REAL8, &
                      MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      if (rank == 0) then
        call output_filename('fxvx', time_now, filename)
        open(newunit=unit, file=trim(filename), status='replace', action='write')
        do i = 0, nx
          do j = 0, nvx
            write(unit,*) global_values(i,j)
          end do
        end do
        close(unit)
      end if
      deallocate(global_values)
    end if
  end subroutine write_snapshot_outputs

  logical function snapshot_output_enabled()
    snapshot_output_enabled = cfg_output_fxvx .or. cfg_output_fvxvy .or. &
      cfg_output_rho .or. cfg_output_rhov .or. cfg_output_fields
  end function snapshot_output_enabled

  logical function diagnostic_output_enabled()
    diagnostic_output_enabled = cfg_output_mass .or. cfg_output_charge .or. cfg_output_energy .or. &
      cfg_output_kinetic_components .or. cfg_output_field_l2 .or. &
      cfg_output_poisson .or. cfg_output_poisson_tau .or. cfg_output_momentum
  end function diagnostic_output_enabled

  subroutine append_momentum(momentum)
    real(8), intent(in) :: momentum(2)
    character(len=512) :: basename, filename
    integer :: unit, exit_status

    if (.not. momentum_directory_ready) then
      call execute_command_line('mkdir -p momentum', wait=.true., &
                                exitstat=exit_status)
      if (exit_status /= 0) error stop 'Cannot create momentum output directory.'
      momentum_directory_ready = .true.
    end if
    call scalar_filename('momentum', basename)
    filename = 'momentum/'//trim(basename)
    call open_diagnostic_file(filename, unit)
    write(unit,*) momentum(1), momentum(2)
    close(unit)
  end subroutine append_momentum

  subroutine append_kinetic_components(kinetic_1, kinetic_2)
    real(8), intent(in) :: kinetic_1, kinetic_2
    character(len=512) :: filename
    integer :: unit
    call scalar_filename('kinetic_components', filename)
    call open_diagnostic_file(filename, unit)
    write(unit,*) kinetic_1, kinetic_2
    close(unit)
  end subroutine append_kinetic_components

  subroutine append_field_l2(field_l2)
    real(8), intent(in) :: field_l2(3)
    character(len=512) :: filename
    integer :: unit
    call scalar_filename('field_l2', filename)
    call open_diagnostic_file(filename, unit)
    write(unit,*) field_l2
    close(unit)
  end subroutine append_field_l2

  subroutine append_scalar(quantity, time_now, value)
    character(len=*), intent(in) :: quantity
    real(8), intent(in) :: time_now, value
    character(len=512) :: filename
    integer :: unit

    call scalar_filename(quantity, filename)
    call open_diagnostic_file(filename, unit)
    write(unit,*) time_now, value
    close(unit)
  end subroutine append_scalar

  subroutine open_diagnostic_file(filename, unit)
    character(len=*), intent(in) :: filename
    integer, intent(out) :: unit
    if (scalar_files_started) then
      open(newunit=unit, file=trim(filename), status='unknown', &
           action='write', position='append')
    else
      open(newunit=unit, file=trim(filename), status='replace', action='write')
    end if
  end subroutine open_diagnostic_file

  subroutine scalar_filename(quantity, filename)
    character(len=*), intent(in) :: quantity
    character(len=*), intent(out) :: filename
    character(len=64) :: ep_label, dt_label
    call compact_real_label(mesh_ep, ep_label)
    call compact_real_label(mesh_dt, dt_label)
    if (len_trim(cfg_output_prefix) > 0) then
      write(filename,'(a,"_",a,"_ep",a,"_dt",a,".dat")') &
        trim(cfg_output_prefix), trim(quantity), trim(ep_label), trim(dt_label)
    else
      write(filename,'(a,"_ep",a,"_dt",a,".dat")') &
        trim(quantity), trim(ep_label), trim(dt_label)
    end if
  end subroutine scalar_filename

  subroutine output_filename(quantity, time_now, filename)
    character(len=*), intent(in) :: quantity
    real(8), intent(in) :: time_now
    character(len=*), intent(out) :: filename
    character(len=64) :: time_label, ep_label, dt_label

    call compact_real_label(time_now, time_label)
    call compact_real_label(mesh_ep, ep_label)
    call compact_real_label(mesh_dt, dt_label)
    if (len_trim(cfg_output_prefix) > 0) then
      write(filename,'(a,"_",a,"_T",a,"_ep",a,"_dt",a,".dat")') &
        trim(cfg_output_prefix), trim(quantity), trim(time_label), &
        trim(ep_label), trim(dt_label)
    else
      write(filename,'(a,"_T",a,"_ep",a,"_dt",a,".dat")') &
        trim(quantity), trim(time_label), trim(ep_label), trim(dt_label)
    end if
  end subroutine output_filename

  subroutine compact_real_label(value, label)
    real(8), intent(in) :: value
    character(len=*), intent(out) :: label
    character(len=64) :: buffer
    integer :: last

    write(buffer,'(f24.12)') value
    buffer = adjustl(buffer)
    last = len_trim(buffer)
    do while (last > 1 .and. buffer(last:last) == '0')
      last = last-1
    end do
    if (buffer(last:last) == '.') last = last-1
    label = buffer(:last)
  end subroutine compact_real_label

end module pic_output
