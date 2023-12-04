module solar_data
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use shr_kind_mod,   only: shr_kind_cl
  use spmd_utils,     only: masterproc
  use cam_abortutils, only: endrun
  use pio
  use cam_pio_utils, only : cam_pio_openfile
  use cam_logfile,    only: iulog
  use infnan,         only: nan, assignment(=)
! use phys_control,   only: use_simple_phys

  implicit none

  save
  private
  public :: solar_data_readnl
  public :: solar_data_init
  public :: solar_data_advance

  character(len=shr_kind_cl) :: solar_irrad_data_file = 'NONE'
  character(len=shr_kind_cl) :: solar_parms_data_file = 'NONE'
  character(len=shr_kind_cl) :: solar_euv_data_file   = 'NONE'
  character(len=shr_kind_cl) :: solar_wind_data_file  = 'NONE'

  character(len=8)   :: solar_data_type = 'SERIAL'      ! "FIXED" or "SERIAL"
  integer            :: solar_data_ymd = -99999999      ! YYYYMMDD for "FIXED" type
  integer            :: solar_data_tod = 0              ! seconds of day for "FIXED" type
  real(r8)           :: solar_const = -9999._r8         ! constant TSI (W/m2)
  logical            :: solar_htng_spctrl_scl = .false. ! do rad heating spectral scaling


!#ifdef WACCM_TSMLT
!  public :: nbins, we  ! number of wavelength samples of spectrum, wavelength endpoints
!  public :: sol_etf
!  public :: sol_irrad
!  public :: sol_tsi
!  public :: do_spctrl_scaling
!  public :: has_spectrum
!  public :: has_ref_spectrum
!  public :: ssi_ref
!  public :: ref_tsi
!
!  integer :: nbins
!  integer :: ntimes
!  real(r8), allocatable :: sol_etf(:)
!  real(r8), allocatable :: irradi(:,:)
!  real(r8)              :: itsi(2)
!  real(r8), allocatable :: ssi_ref(:)  ! a reference spectrum constructed from 3 solar cycles of data
!
!  real(r8), allocatable :: we(:)
!  real(r8), allocatable :: data_times(:)
!
!  integer :: last_index = 1
!  type(file_desc_t) :: file_id
!  integer :: ssi_vid
!  integer :: tsi_vid
!  integer :: ref_vid
!  integer :: tsi_ref_vid
!
!  logical :: initialized = .false.
!  logical :: has_spectrum = .false.
!  logical :: has_ref_spectrum = .false.
!  logical :: has_tsi = .false.
!
!  real(r8), allocatable :: sol_irrad(:)
!  real(r8)              :: sol_tsi = -1.0_r8
!  real(r8)              :: ref_tsi
!
!  real(r8), allocatable :: irrad_fac(:)
!  real(r8), allocatable :: etf_fac(:)
!  real(r8), allocatable :: dellam(:)
!
!  logical  :: fixed_scon
!  logical  :: fixed_solar
!  real(r8) :: offset_time
!
!  logical            :: do_spctrl_scaling = .false.
!
!  character(len=256) :: solar_data_file = ''
!#endif

 contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine solar_data_readnl( nlfile )
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use spmd_utils,      only: mpicom, masterprocid, mpi_character, mpi_integer, mpi_logical, mpi_real8
    use solar_parms_data,only: solar_parms_on
    use solar_wind_data, only: solar_wind_on
    
    ! arguments
    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! local vars
    integer :: unitn, ierr

    namelist /solar_data_opts/ &
         solar_irrad_data_file, solar_parms_data_file, solar_euv_data_file, solar_wind_data_file, &
         solar_data_type, solar_data_ymd, solar_data_tod, solar_const, solar_htng_spctrl_scl
    
!   if (use_simple_phys) return

    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'solar_data_opts', status=ierr)
       if (ierr == 0) then
          read(unitn, solar_data_opts, iostat=ierr)
          if (ierr /= 0) then
             call endrun('solar_data_readnl: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

!#ifdef WACCM_TSMLT
!    if ( (len_trim(solar_data_file) > 0) .and. (solar_const>0._r8) ) then
!       call endrun('solar_data_readnl: ERROR cannot specify both solar_data_file and solar_const')
!    endif
!
!    if ( len_trim(solar_data_file) > 0 ) then
!       fixed_scon = .false.
!    else
!       fixed_scon = .true.
!    endif
!
!    if ( solar_const>0._r8 ) then
!       sol_tsi = solar_const
!    endif
!    ref_tsi = nan
!#endif

    ! broadcast the options to all MPI tasks
    call mpi_bcast(solar_irrad_data_file, len(solar_irrad_data_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(solar_parms_data_file, len(solar_parms_data_file), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(solar_euv_data_file,   len(solar_euv_data_file),   mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(solar_wind_data_file,  len(solar_wind_data_file),  mpi_character, masterprocid, mpicom, ierr)

    call mpi_bcast(solar_data_type, len(solar_data_type), mpi_character, masterprocid, mpicom, ierr)
    call mpi_bcast(solar_data_ymd,  1,                    mpi_integer,   masterprocid, mpicom, ierr)
    call mpi_bcast(solar_data_tod,  1,                    mpi_integer,   masterprocid, mpicom, ierr)
    call mpi_bcast(solar_const,     1,                    mpi_real8 ,    masterprocid, mpicom, ierr)
    call mpi_bcast(solar_htng_spctrl_scl,1,               mpi_logical,   masterprocid, mpicom, ierr)

    if ( (solar_irrad_data_file.ne.'NONE') .and. (solar_const>0._r8) ) then
       call endrun('solar_data_readnl: ERROR cannot specify both solar_irrad_data_file and solar_const')      
    endif

    if ( (solar_data_ymd>0 .or. solar_data_tod>0) .and. trim(solar_data_type)=='SERIAL' ) then
       call endrun('solar_data_readnl: ERROR cannot set solar_data_ymd or solar_data_tod with solar_data_type=SERIAL')
    endif

    if (masterproc) then
       write(iulog,*) 'solar_data_readnl: solar_const (W/m2) = ', solar_const
       write(iulog,*) 'solar_data_readnl: solar_irrad_data_file = ',trim(solar_irrad_data_file)
       write(iulog,*) 'solar_data_readnl: solar_parms_data_file = ',trim(solar_parms_data_file)
       write(iulog,*) 'solar_data_readnl: solar_euv_data_file = ',trim(solar_euv_data_file)
       write(iulog,*) 'solar_data_readnl: solar_wind_data_file = ',trim(solar_wind_data_file)
       write(iulog,*) 'solar_data_readnl: solar_data_type = ',trim(solar_data_type)
       write(iulog,*) 'solar_data_readnl: solar_data_ymd  = ',solar_data_ymd
       write(iulog,*) 'solar_data_readnl: solar_data_tod  = ',solar_data_tod
    endif

    solar_parms_on = solar_parms_data_file.ne.'NONE'
    solar_wind_on = solar_wind_data_file.ne.'NONE'

  end subroutine solar_data_readnl

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
  subroutine solar_data_init()
    use solar_irrad_data, only: solar_irrad_init
    use solar_parms_data, only: solar_parms_init
    use solar_wind_data,  only: solar_wind_init
    use solar_euv_data,   only: solar_euv_init

    logical :: fixed_solar
    fixed_solar = trim(solar_data_type) == 'FIXED'

    call solar_irrad_init( solar_irrad_data_file, fixed_solar, solar_data_ymd, solar_data_tod, &
                           solar_const, solar_htng_spctrl_scl )
    call solar_parms_init( solar_parms_data_file, fixed_solar, solar_data_ymd, solar_data_tod )
    call solar_wind_init( solar_wind_data_file, fixed_solar, solar_data_ymd, solar_data_tod )
    call solar_euv_init( solar_euv_data_file, fixed_solar, solar_data_ymd, solar_data_tod )

   end subroutine solar_data_init

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
   subroutine solar_data_advance()

     use solar_irrad_data, only: solar_irrad_advance
     use solar_parms_data, only: solar_parms_advance
     use solar_wind_data,  only: solar_wind_advance
     use solar_euv_data,   only: solar_euv_advance

     call solar_irrad_advance()
     call solar_parms_advance()
     call solar_wind_advance()
     call solar_euv_advance()

   end subroutine solar_data_advance

end module solar_data
