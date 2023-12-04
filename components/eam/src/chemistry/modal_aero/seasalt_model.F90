!===============================================================================
! Seasalt for Modal Aerosol Model
!===============================================================================
module seasalt_model
  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use ppgrid,         only: pcols, pver
  use modal_aero_data,only: ntot_amode, nslt=>nSeaSalt
  use cam_abortutils, only: endrun
  use tracer_data,    only: trfld, trfile
  use cam_logfile,    only: iulog

  implicit none
  private

  public :: seasalt_nbin
  public :: seasalt_nnum
  public :: seasalt_names
  public :: seasalt_indices
  public :: seasalt_init
  public :: seasalt_emis
  public :: seasalt_active

#ifdef WACCM_TSMLT
  public :: n_ocean_data
  public :: nslt_om
  public :: F_eff_out                 ! Output effective enrichment ratio?
                                      !  (logical, currently set to FALSE)
  public :: has_mam_mom               ! run with marine organics?
                                      !  (logical, set to TRUE if user supplies file)
  public :: advance_ocean_data        ! advance ocean data in time
  public :: init_ocean_data           ! initialize ocean data variables
! public :: ocean_data_readnl         ! read ocean data namelist
#endif

  integer, protected :: seasalt_nbin ! = nslt
  integer, protected :: seasalt_nnum ! = nnum

  character(len=6), protected, allocatable :: seasalt_names(:)
  integer, protected, allocatable :: seasalt_indices(:)

  logical :: seasalt_active = .false.

  real(r8):: emis_scale

! TODO SMB: Implement better mechanism for setting this switch.
#if (defined MODAL_AERO_9MODE || defined MODAL_AERO_4MODE_MOM || defined WACCM_TSMLT)
   logical :: has_mam_mom = .true.
#else
   logical :: has_mam_mom = .false.
#endif

#ifdef WACCM_TSMLT
  integer, parameter :: & ! number of ocean data fields
       n_ocean_data = 4
  logical, parameter   :: F_eff_out = .false.
  type(trfld), pointer :: fields(:)
  type(trfile)         :: file
  real(r8), parameter :: small_oceanorg = 1.0e-30
! Namelist variables related to dataset specification
  character(len=32)   :: specifier(n_ocean_data) = ''
  character(len=256)  :: filename = ' '
  character(len=256)  :: filelist = ' '
  character(len=256)  :: datapath = ' '
  character(len=32)   :: datatype = 'CYCLICAL'
  integer             :: data_cycle_yr = 0
  logical             :: rmv_file = .false.
  integer             :: fixed_ymd = 0
  integer             :: fixed_tod = 0
  logical :: debug_mam_mom = .false.

#if  ( defined MODAL_AERO_7MODE )
  integer, parameter :: nslt_om = 0
  integer, parameter :: nnum_om = 0
  integer, parameter :: om_num_modes = 0
! character(len=6),parameter :: seasalt_names(nslt+nslt_om+nnum+nnum_om) = &
!      (/ 'ncl_a1', 'ncl_a2', 'ncl_a4', 'ncl_a6', 'num_a1', 'num_a2', 'num_a4', 'num_a6' /)
  integer, parameter :: om_num_ind = 0
#elif( defined MODAL_AERO_3MODE || defined MODAL_AERO_4MODE )
  integer, parameter :: nslt_om = 0
  integer, parameter :: nnum_om = 0
  integer, parameter :: om_num_modes = 0
! character(len=6),parameter :: seasalt_names(nslt+nslt_om+nnum+nnum_om) = &
!      (/ 'ncl_a1', 'ncl_a2', 'ncl_a3', &
!         'num_a1', 'num_a2', 'num_a3'/)
  integer, parameter :: om_num_ind = 0
#elif( defined MODAL_AERO_4MODE_MOM || defined WACCM_TSMLT )
  integer, parameter :: nslt_om = 3
  integer, parameter :: nnum_om = 1
  integer, parameter :: om_num_modes = 3
! character(len=6),parameter :: seasalt_names(nslt+nslt_om+nnum+nnum_om) = &
!      (/ 'ncl_a1', 'ncl_a2', 'ncl_a3', &
!      'mom_a1', 'mom_a2', 'mom_a4', &
!      'num_a1', 'num_a2', 'num_a3', 'num_a4'/)
  integer, dimension(om_num_modes), parameter :: om_num_ind =  (/ 1, 2, 4 /)
#elif (defined MODAL_AERO_9MODE)
  integer, parameter :: nslt_om = 12
  integer, parameter :: nnum_om = 2
  integer, parameter :: om_num_modes = 4
! character(len=8),parameter :: seasalt_names(nslt+nslt_om+nnum+nnum_om) = &
!      (/'ncl_a1  ', 'ncl_a2  ', 'ncl_a4  ', 'ncl_a6  ', &
!        'mpoly_a1', 'mpoly_a2', 'mpoly_a8', 'mpoly_a9', &
!        'mprot_a1', 'mprot_a2', 'mprot_a8', 'mprot_a9', &
!        'mlip_a1 ', 'mlip_a2 ', 'mlip_a8 ', 'mlip_a9 ', &
!        'num_a1  ', 'num_a2  ', 'num_a4  ', 'num_a6  ', &
!        'num_a8  ', 'num_a9  ' &
!        /)
  integer, dimension(om_num_modes), parameter :: om_num_ind =  (/ 1, 2, 5, 6 /)
#endif

#endif

contains
  
  !=============================================================================
  !=============================================================================
  subroutine seasalt_init(seasalt_emis_scale)
    use sslt_sections, only: sslt_sections_init
    use constituents,  only: cnst_get_ind
    use rad_constituents, only: rad_cnst_get_info

    real(r8), intent(in) :: seasalt_emis_scale
    integer :: m, l, nspec, ndx
    character(len=32) :: spec_name
    
    seasalt_nbin = nslt
    seasalt_nnum = nslt
    allocate(seasalt_names(2*nslt))
    allocate(seasalt_indices(2*nslt))

    ndx=0
    do m = 1, ntot_amode
       call rad_cnst_get_info(0, m, nspec=nspec)
       do l = 1, nspec
          call rad_cnst_get_info(0, m, l, spec_name=spec_name )
          if (spec_name(:3) == 'ncl') then
             ndx=ndx+1
             seasalt_names(ndx) = spec_name
             seasalt_names(nslt+ndx) = 'num_'//spec_name(5:)
             call cnst_get_ind(seasalt_names(     ndx), seasalt_indices(     ndx))
             call cnst_get_ind(seasalt_names(nslt+ndx), seasalt_indices(nslt+ndx))
          endif
       enddo
    enddo

    seasalt_active = any(seasalt_indices(:) > 0)
    if (.not.seasalt_active) return

    call sslt_sections_init()

    emis_scale = seasalt_emis_scale

  end subroutine seasalt_init

  !=============================================================================
  !=============================================================================
  subroutine seasalt_emis( u10cubed,  srf_temp, ocnfrc, ncol, cflx )

    use sslt_sections, only: nsections, fluxes, Dg, rdry
    use mo_constants,  only: dns_aer_sst=>seasalt_density, pi

    ! dummy arguments
    real(r8), intent(in) :: u10cubed(:)
    real(r8), intent(in) :: srf_temp(:)
    real(r8), intent(in) :: ocnfrc(:)
    integer,  intent(in) :: ncol
    real(r8), intent(inout) :: cflx(:,:)

    ! local vars
    integer  :: mn, mm, ibin, isec, i
    real(r8) :: fi(ncol,nsections)

    real(r8) :: sst_sz_range_lo (nslt)
    real(r8) :: sst_sz_range_hi (nslt)

    if (nslt==4) then
       sst_sz_range_lo (:) = (/ 0.08e-6_r8, 0.02e-6_r8, 0.3e-6_r8,  1.0e-6_r8 /) ! accu, aitken, fine, coarse
       sst_sz_range_hi (:) = (/ 0.3e-6_r8,  0.08e-6_r8, 1.0e-6_r8, 10.0e-6_r8 /)
    else if (nslt==3) then
       sst_sz_range_lo (:) =  (/ 0.08e-6_r8,  0.02e-6_r8,  1.0e-6_r8 /)  ! accu, aitken, coarse
       sst_sz_range_hi (:) =  (/ 1.0e-6_r8,   0.08e-6_r8, 10.0e-6_r8 /)
    endif

    fi(:ncol,:nsections) = fluxes( srf_temp, u10cubed, ncol )

    do ibin = 1,nslt
       mm = seasalt_indices(ibin)
       mn = seasalt_indices(nslt+ibin)
       
       if (mn>0) then
          do i=1, nsections
             if (Dg(i).ge.sst_sz_range_lo(ibin) .and. Dg(i).lt.sst_sz_range_hi(ibin)) then
                cflx(:ncol,mn)=cflx(:ncol,mn)+fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  !++ ag: scale sea-salt
             endif
          enddo
       endif

       cflx(:ncol,mm)=0.0_r8
       do i=1, nsections
          if (Dg(i).ge.sst_sz_range_lo(ibin) .and. Dg(i).lt.sst_sz_range_hi(ibin)) then
             cflx(:ncol,mm)=cflx(:ncol,mm)+fi(:ncol,i)*ocnfrc(:ncol)*emis_scale  &   !++ ag: scale sea-salt
                  *4._r8/3._r8*pi*rdry(i)**3*dns_aer_sst  ! should use dry size, convert from number to mass flux (kg/m2/s)
          endif
       enddo

    enddo

  end subroutine seasalt_emis

#ifdef WACCM_TSMLT
!-------------------------------------------------------------------
!! READ INPUT FILES, CREATE FIELDS, and horizontally interpolate ocean data
!-------------------------------------------------------------------
subroutine init_ocean_data()
    !-----------------------------------------------------------------------
    !
    ! Purpose: register advected constituents for all aerosols
    !
    ! Method:
    ! <Describe the algorithm(s) used in the routine.>
    ! <Also include any applicable external references.>
    !
    ! Author: S. M. Burrows, adapted from dust_initialize
    !
    !-----------------------------------------------------------------------

    use tracer_data,      only : trcdata_init
    use cam_history,      only : addfld, horiz_only, add_default
    use spmd_utils,       only : masterproc
    use sslt_sections,    only : nsections

    !-----------------------------------------------------------------------
    !    ... local variables
    !-----------------------------------------------------------------------

!    type(interp_type)     :: lon_wgts, lat_wgts
!    real(r8), parameter   :: zero=0._r8, twopi=2._r8*pi

    integer :: i, m, m_om
    integer :: number_flds

    if ( masterproc ) then
       write(iulog,*) 'ocean organics are prescribed in :'//trim(filename)
    endif

!    allocate (file%in_pbuf(size(specifier)))
     allocate (file%in_pbuf(n_ocean_data))
     file%in_pbuf(:) = .false.
!    file%in_pbuf(:) = .true.
!
!       fields(i)%pbuf_ndx = pbuf_get_index(fields(i)%fldnam,errcode)

    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, data_cycle_yr, fixed_ymd, fixed_tod, datatype)

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if ( number_flds .eq. n_ocean_data ) then
       if ( masterproc ) then
          write(iulog,"(A21,I3,A20)") 'Successfully read in ',number_flds,' ocean data fields'
       endif
    else if( number_flds < 1 ) then
       if ( masterproc ) then
          write(iulog,*) 'Failed to read in any ocean data'
          write(iulog,*) ' '
       endif
       return
    else if ( number_flds .ne. n_ocean_data ) then
       if ( masterproc ) then
          write(iulog,"(A8,I3,A20)") 'Read in ',number_flds,' ocean data fields'
          write(iulog,"(A9,I3,A20)") 'Expected ',n_ocean_data,' ocean data fields'
          write(iulog,*) ' '
          return
       endif
    end if

    ! Following loop adds fields for output.
    !   Note that the units are given in fields(i)%units, avgflag='A' indicates output mean
    fldloop:do i = 1,n_ocean_data

       if ( masterproc ) then
          write(iulog,*) 'adding field '//fields(i)%fldnam//' ...'
       endif

!!$       ! Set names of variable tendencies and declare them as history variables
!!$       !    addfld(fname,                 unite,              numlev, avgflag, long_name, decomp_type, ...)
       if ( trim(fields(i)%fldnam) == "chla" ) then
          call addfld(trim(fields(i)%fldnam), horiz_only, 'A', 'mg L-1 ', 'ocean input data: '//fields(i)%fldnam )
          call add_default (fields(i)%fldnam, 1, ' ')
       else
          call addfld(trim(fields(i)%fldnam), horiz_only, 'A', 'uM C ', 'ocean input data: '//fields(i)%fldnam )
          call add_default (fields(i)%fldnam, 1, ' ')
       endif

    enddo fldloop

! FOR DEBUGGING
    debug: if (debug_mam_mom) then
       call addfld('mpoly_debug', horiz_only, 'A', ' ', 'mpoly_debug' )
       call add_default ('mpoly_debug', 1, ' ')

       call addfld('mass_frac_bub_tot', horiz_only, 'A', ' ', 'total organic mass fraction of bubble' )
       call add_default ('mass_frac_bub_tot', 1, ' ')

       call addfld('mass_frac_bub_poly', horiz_only, 'A', ' ', 'total organic mass fraction (poly)' )
       call add_default ('mass_frac_bub_poly', 1, ' ')
       call addfld('mass_frac_bub_prot', horiz_only, 'A', ' ', 'total organic mass fraction (prot)' )
       call add_default ('mass_frac_bub_prot', 1, ' ')

       call addfld('mass_frac_bub_lip', horiz_only, 'A', ' ', 'total organic mass fraction (lip)' )
       call add_default ('mass_frac_bub_lip', 1, ' ')

       om_mode_loop: do m_om=1,nslt_om
#if ( defined MODAL_AERO_9MODE )
          m = nslt+(n-1)*om_num_modes+m_om
#elif ( defined MODAL_AERO_4MODE_MOM || defined WACCM_TSMLT)
          m = nslt+m_om
#endif
          call addfld('cflx_'//trim(seasalt_names(m))//'_debug', horiz_only, 'A', ' ', 'accumulation organic mass emissions' )
          call add_default ('cflx_'//trim(seasalt_names(m))//'_debug', 1, ' ')
       enddo om_mode_loop

       call addfld('omf_bub_section_mpoly', horiz_only, 'A',' ', 'omf poly' )
       call add_default ('omf_bub_section_mpoly', 1, ' ')

       call addfld('omf_bub_section_mprot', horiz_only, 'A',' ', 'omf prot' )
       call add_default ('omf_bub_section_mprot', 1, ' ')

       call addfld('omf_bub_section_mlip', horiz_only, 'A',' ', 'omf lip' )
       call add_default ('omf_bub_section_mlip', 1, ' ')

    endif debug

    if ( masterproc ) then
       write(iulog,*) 'Done initializing marine organics data'
    endif

  end subroutine init_ocean_data

!-------------------------------------------------------------------
! Advance ocean data fields to the current time step
!
! Adapted from prescribed_aero_adv
!
! Author: Susannah M. Burrows
! Date: 13 Jan 2015
!-------------------------------------------------------------------
subroutine advance_ocean_data(state, pbuf2d)
    use physics_types,  only : physics_state
    use tracer_data,    only : advance_trcdata, get_fld_data, put_fld_data
    use ppgrid,         only : begchunk, endchunk
    use physics_buffer, only : physics_buffer_desc, pbuf_get_chunk
    use cam_history,    only : outfld
    use spmd_utils,     only : masterproc

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf_chnk(:)

    integer :: i,c,ncol
!    real(r8),pointer :: outdata(:,:)
!    real(r8) :: outdata(pcols,begchunk:endchunk)
    real(r8) :: outdata(pcols,1)
    integer lchnk

!    write(iulog,*) 'Advancing ocean data ...' ! for debugging
!   if ( masterproc ) then
!      write(iulog,*) 'before advance_trcdata in ocean_data'
!   endif

    call advance_trcdata( fields, file, state, pbuf2d )

!   if ( masterproc ) then
!      write(iulog,*) 'after  advance_trcdata in ocean_data'
!   endif
!    write(iulog,*) 'Done advancing ocean data ...' ! for debugging

! Add new values to history files
    fldloop:do i = 1,n_ocean_data

       chnkloop: do c = begchunk,endchunk
          ncol = state(c)%ncol
          pbuf_chnk => pbuf_get_chunk(pbuf2d, c)
          lchnk = state(c)%lchnk

          call get_fld_data( fields, fields(i)%fldnam, outdata(:ncol,:), ncol, lchnk, pbuf_chnk)

          ! work-around for interpolation errors that introduce negative values
          ! near coasts: reset negative values to zero.
          where (outdata(:ncol,1) < small_oceanorg)
             outdata(:ncol,1) = 0.0_r8
          end where

          call put_fld_data( fields, fields(i)%fldnam, outdata(:ncol,:), ncol, lchnk, pbuf_chnk)

          ! The following line is probably redundant but is included for safety
          call get_fld_data( fields, fields(i)%fldnam, outdata(:ncol,:), ncol, lchnk, pbuf_chnk)

          call outfld( trim(fields(i)%fldnam), outdata(:ncol,1), ncol, lchnk )
       enddo chnkloop

    enddo fldloop

end subroutine advance_ocean_data

#endif

end module seasalt_model
