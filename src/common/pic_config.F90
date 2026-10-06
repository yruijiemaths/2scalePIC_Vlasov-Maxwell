module pic_config
  implicit none
  private

  integer, public :: cfg_initial_condition = 2
  integer, public :: cfg_particles_per_cell = 10000
  integer, public :: cfg_nx = 32
  integer, public :: cfg_nvx = 40
  integer, public :: cfg_nvy = 40
  integer, public :: cfg_ntau = 16
  real(8), public :: cfg_domain_length = -1.d0
  real(8), public :: cfg_velocity_x_length = 12.d0
  real(8), public :: cfg_velocity_y_length = 12.d0
  real(8), public :: cfg_epsilon = -1.d0
  real(8), public :: cfg_dt = -1.d0
  real(8), public :: cfg_tfinal = -1.d0
  real(8), public :: cfg_density_perturbation = 5.d-2
  real(8), public :: cfg_two_stream_sigma_vx = 1.d0
  real(8), public :: cfg_two_stream_sigma_vy = 1.4142135623730951d0
  real(8), public :: cfg_weibel_sigma_vx = 1.d0
  real(8), public :: cfg_weibel_sigma_vy = 2.d0
  real(8), public :: cfg_stream_velocity = 2.d0
  real(8), public :: cfg_wave_number = 1.d0
  real(8), public :: cfg_bz_two_stream = 1.d-6
  real(8), public :: cfg_bz_weibel = 1.d-2
  real(8), public :: cfg_bump_tail_amplitude = 2.d0/9.d0
  real(8), public :: cfg_bump_drift_velocity = 4.5d0
  real(8), public :: cfg_bump_thermal_velocity = 0.5d0
  real(8), public :: cfg_bump_bulk_sigma_vx = 1.d0
  real(8), public :: cfg_bump_sigma_vy = 1.d0
  logical, public :: cfg_output_fxvx = .true.
  logical, public :: cfg_output_fvxvy = .true.
  logical, public :: cfg_output_rho = .false.
  logical, public :: cfg_output_rhov = .false.
  logical, public :: cfg_output_fields = .false.
  logical, public :: cfg_output_mass = .false.
  logical, public :: cfg_output_charge = .false.
  logical, public :: cfg_output_energy = .false.
  logical, public :: cfg_output_kinetic_components = .false.
  logical, public :: cfg_output_field_l2 = .false.
  logical, public :: cfg_output_poisson = .false.
  logical, public :: cfg_output_poisson_tau = .false.
  logical, public :: cfg_output_momentum = .false.
  logical, public :: cfg_output_initial = .false.
  logical, public :: cfg_output_final = .true.
  real(8), public :: cfg_snapshot_interval = 100.d0
  integer, public :: cfg_diagnostic_interval_steps = 1
  character(len=128), public :: cfg_output_prefix = ''
  character(len=16), public :: cfg_implicit_solver = 'iterative'
  real(8), public :: cfg_iterative_tolerance = 1.d-12
  integer, public :: cfg_iterative_restart = 30
  integer, public :: cfg_iterative_max_iterations = 300

  public :: load_pic_config, apply_pic_config, print_pic_config

contains

  subroutine load_pic_config(rank)
    integer, intent(in) :: rank
    integer :: unit, ios, env_status, env_length
    logical :: exists
    character(len=512) :: filename, message
    namelist /pic/ cfg_initial_condition, cfg_particles_per_cell, &
      cfg_nx, cfg_nvx, cfg_nvy, cfg_ntau, cfg_domain_length, &
      cfg_velocity_x_length, cfg_velocity_y_length, cfg_epsilon, cfg_dt, &
      cfg_tfinal, cfg_density_perturbation, cfg_two_stream_sigma_vx, &
      cfg_two_stream_sigma_vy, cfg_weibel_sigma_vx, cfg_weibel_sigma_vy, &
      cfg_stream_velocity, cfg_wave_number, cfg_bz_two_stream, cfg_bz_weibel, &
      cfg_bump_tail_amplitude, cfg_bump_drift_velocity, &
      cfg_bump_thermal_velocity, cfg_bump_bulk_sigma_vx, cfg_bump_sigma_vy, &
      cfg_implicit_solver, cfg_iterative_tolerance, cfg_iterative_restart, &
      cfg_iterative_max_iterations
    namelist /pic_output/ cfg_output_fxvx, cfg_output_fvxvy, cfg_output_rho, &
      cfg_output_rhov, cfg_output_fields, cfg_output_mass, cfg_output_charge, cfg_output_energy, &
      cfg_output_kinetic_components, cfg_output_field_l2, &
      cfg_output_poisson, cfg_output_poisson_tau, cfg_output_momentum, &
      cfg_output_initial, cfg_output_final, cfg_snapshot_interval, &
      cfg_diagnostic_interval_steps, cfg_output_prefix

    filename = 'config/pic.nml'
    call get_environment_variable('PIC_CONFIG', filename, length=env_length, &
                                  status=env_status)
    if (env_status /= 0 .or. env_length == 0) filename = 'config/pic.nml'

    inquire(file=trim(filename), exist=exists)
    if (.not. exists) then
      if (rank == 0) write(*,'(a)') 'PIC config not found; using method defaults: '//trim(filename)
      return
    end if

    open(newunit=unit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'Cannot open PIC configuration file.'
    read(unit, nml=pic, iostat=ios, iomsg=message)
    if (ios /= 0) then
      if (rank == 0) write(*,'(a)') trim(message)
      error stop 'Invalid PIC configuration file.'
    end if

    rewind(unit)
    read(unit, nml=pic_output, iostat=ios, iomsg=message)
    close(unit)
    if (ios /= 0 .and. ios /= -1) then
      if (rank == 0) write(*,'(a)') trim(message)
      error stop 'Invalid PIC output configuration.'
    end if

    call validate_pic_config()
    if (rank == 0) write(*,'(a)') 'Loaded PIC configuration: '//trim(filename)
  end subroutine load_pic_config

  subroutine apply_pic_config(nx, nvx, nvy, dimx, dimvx, dimvy, ep, dt, tfinal)
    integer(8), intent(inout) :: nx, nvx, nvy
    real(8), intent(inout) :: dimx, dimvx, dimvy, ep, dt, tfinal

    nx = int(cfg_nx,8)
    nvx = int(cfg_nvx,8)
    nvy = int(cfg_nvy,8)
    if (cfg_domain_length > 0.d0) then
      dimx = cfg_domain_length
    else
      ! One wavelength is the natural default periodic domain.
      dimx = 8.d0*datan(1.d0)/cfg_wave_number
    end if
    dimvx = cfg_velocity_x_length
    dimvy = cfg_velocity_y_length
    if (cfg_epsilon > 0.d0) ep = cfg_epsilon
    if (cfg_dt > 0.d0) dt = cfg_dt
    if (cfg_tfinal > 0.d0) tfinal = cfg_tfinal
  end subroutine apply_pic_config

  subroutine validate_pic_config()
    if (cfg_initial_condition < 1 .or. cfg_initial_condition > 3) &
      error stop 'cfg_initial_condition must be 1, 2 or 3.'
    if (cfg_particles_per_cell <= 0) error stop 'cfg_particles_per_cell must be positive.'
    if (cfg_nx <= 0 .or. cfg_nvx <= 0 .or. cfg_nvy <= 0 .or. cfg_ntau <= 0) &
      error stop 'cfg_nx, cfg_nvx, cfg_nvy and cfg_ntau must be positive.'
    if (mod(cfg_ntau,2) /= 0) error stop 'cfg_ntau must be even.'
    if (cfg_velocity_x_length <= 0.d0 .or. cfg_velocity_y_length <= 0.d0) &
      error stop 'Velocity-domain lengths must be positive.'
    if (cfg_wave_number <= 0.d0) error stop 'cfg_wave_number must be positive.'
    if (cfg_epsilon == 0.d0 .or. cfg_dt == 0.d0 .or. cfg_tfinal == 0.d0) &
      error stop 'epsilon, dt and tfinal must be positive when specified.'
    if (cfg_epsilon < 0.d0 .and. cfg_epsilon /= -1.d0) &
      error stop 'cfg_epsilon must be positive or -1 to use the method default.'
    if (cfg_dt < 0.d0 .and. cfg_dt /= -1.d0) &
      error stop 'cfg_dt must be positive or -1 to use the method default.'
    if (cfg_tfinal < 0.d0 .and. cfg_tfinal /= -1.d0) &
      error stop 'cfg_tfinal must be positive or -1 to use the method default.'
    if (cfg_dt > 0.d0 .and. cfg_tfinal > 0.d0) then
      if (nint(cfg_tfinal/cfg_dt,kind=8) < 1_8) &
        error stop 'cfg_tfinal/cfg_dt must produce at least one time step.'
    end if
    if (cfg_two_stream_sigma_vx <= 0.d0 .or. cfg_two_stream_sigma_vy <= 0.d0 .or. &
        cfg_weibel_sigma_vx <= 0.d0 .or. cfg_weibel_sigma_vy <= 0.d0) &
      error stop 'Velocity standard deviations must be positive.'
    if (cfg_bump_tail_amplitude < 0.d0 .or. cfg_bump_thermal_velocity <= 0.d0 .or. &
        cfg_bump_bulk_sigma_vx <= 0.d0 .or. cfg_bump_sigma_vy <= 0.d0) &
      error stop 'Bump-on-tail parameters must have nonnegative tail fraction and positive widths.'
    if (cfg_snapshot_interval <= 0.d0) error stop 'cfg_snapshot_interval must be positive.'
    if (cfg_diagnostic_interval_steps <= 0) &
      error stop 'cfg_diagnostic_interval_steps must be positive.'
    if (trim(cfg_implicit_solver) /= 'iterative' .and. &
        trim(cfg_implicit_solver) /= 'direct') &
      error stop "cfg_implicit_solver must be 'iterative' or 'direct'."
    if (cfg_iterative_tolerance <= 0.d0) &
      error stop 'cfg_iterative_tolerance must be positive.'
    if (cfg_iterative_restart <= 0 .or. cfg_iterative_max_iterations <= 0) &
      error stop 'GMRES restart and maximum iterations must be positive.'
  end subroutine validate_pic_config

  subroutine print_pic_config(rank)
    integer, intent(in) :: rank
    if (rank /= 0) return
    write(*,'(a,3(i0,1x))') 'grid nx/nvx/nvy = ', cfg_nx, cfg_nvx, cfg_nvy
    write(*,'(a,i0)') 'tau grid points = ', cfg_ntau
    write(*,'(a,i0,a,i0)') 'initial condition = ', cfg_initial_condition, &
      ', particles/cell = ', cfg_particles_per_cell
    write(*,'(a,l1,a,l1)') 'output fxvx/fvxvy = ', cfg_output_fxvx, '/', cfg_output_fvxvy
    write(*,'(a,11(l1,1x))') &
      'output rho/rhov/fields/mass/Q/energy/kinetic_components/field_l2/Poisson/Poisson_tau/momentum = ', &
      cfg_output_rho, cfg_output_rhov, cfg_output_fields, cfg_output_mass, &
      cfg_output_charge, cfg_output_energy, cfg_output_kinetic_components, cfg_output_field_l2, &
      cfg_output_poisson, cfg_output_poisson_tau, cfg_output_momentum
    write(*,'(a,f12.4)') 'snapshot time interval = ', cfg_snapshot_interval
    write(*,'(a,i0)') 'diagnostic step interval = ', cfg_diagnostic_interval_steps
    write(*,'(a,a)') 'implicit particle solver = ', trim(cfg_implicit_solver)
    if (trim(cfg_implicit_solver) == 'iterative') &
      write(*,'(a,es10.3,a,i0,a,i0)') 'GMRES tolerance/restart/maxiter = ', &
        cfg_iterative_tolerance, '/', cfg_iterative_restart, '/', &
        cfg_iterative_max_iterations
  end subroutine print_pic_config

end module pic_config
