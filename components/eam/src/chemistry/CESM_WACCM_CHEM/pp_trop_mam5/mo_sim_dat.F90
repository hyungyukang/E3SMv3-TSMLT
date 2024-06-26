
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,     only : diag_map
      use chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,     only : pht_alias_lst, pht_alias_mult
      use chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid, num_rnts, rxntot
      use cam_abortutils,only : endrun
      use mo_tracname,   only : solsym
      use chem_mods,     only : frc_from_dataset
      use chem_mods,     only : is_scalar, is_vector
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use cam_logfile,   only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      is_scalar = .true.
      is_vector = .false.

      clscnt(:) = (/      0,     0,     0,    27,     0 /)

      cls_rxt_cnt(:,4) = (/      1,     6,     0,    27 /)

      solsym(: 27) = (/ 'H2O2            ','H2SO4           ','SO2             ','DMS             ','SOAG            ', &
                        'so4_a1          ','pom_a1          ','soa_a1          ','bc_a1           ','dst_a1          ', &
                        'ncl_a1          ','num_a1          ','so4_a2          ','dst_a2          ','soa_a2          ', &
                        'ncl_a2          ','num_a2          ','dst_a3          ','ncl_a3          ','so4_a3          ', &
                        'num_a3          ','pom_a4          ','bc_a4           ','num_a4          ','so4_a5          ', &
                        'num_a5          ','H2O             ' /)

      adv_mass(: 27) = (/    34.013600_r8,    98.078400_r8,    64.064800_r8,    62.132400_r8,    12.011000_r8, &
                            115.107340_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,   135.064039_r8, &
                             58.442468_r8,     1.007400_r8,   115.107340_r8,   135.064039_r8,    12.011000_r8, &
                             58.442468_r8,     1.007400_r8,   135.064039_r8,    58.442468_r8,   115.107340_r8, &
                              1.007400_r8,    12.011000_r8,    12.011000_r8,     1.007400_r8,   115.107340_r8, &
                              1.007400_r8,    18.014200_r8 /)

      crb_mass(: 27) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,    24.022000_r8,    12.011000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,    12.011000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8 /)

      fix_mass(:  7) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 47.9982000_r8, 17.0068000_r8, &
                          62.0049400_r8, 33.0062000_r8 /)

      clsmap(: 27,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  15,  16,  14,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26,  27 /)

      permute(: 27,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                             11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                             21,  22,  23,  24,  25,  26,  27 /)

      diag_map(: 27) = (/    1,   3,   5,   7,   8,   9,  10,  11,  12,  13, &
                            14,  15,  16,  17,  18,  19,  20,  21,  22,  23, &
                            24,  25,  26,  27,  28,  29,  30 /)

      extfrc_lst(:  9) = (/ 'SO2             ','so4_a1          ','so4_a2          ','pom_a4          ','bc_a4           ', &
                            'num_a1          ','num_a2          ','num_a4          ','H2O             ' /)

      frc_from_dataset(:  9) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true., .true. /)

      inv_lst(:  7) = (/ 'M               ', 'N2              ', 'O2              ', 'O3              ', 'OH              ', &
                         'NO3             ', 'HO2             ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(     1:     4) = (/ 'jh2o2                           ', 'usr_HO2_HO2                     ', &
                                      'usr_SO2_OH                      ', 'usr_DMS_OH                      ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   4,   6 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ '                ' /)
      pht_alias_lst(:,2) = (/ '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8 /)
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate num_rnts; error = ',ios
         call endrun
      end if
      num_rnts(:) = (/      2,     2,     2,     2,     2,     2 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
