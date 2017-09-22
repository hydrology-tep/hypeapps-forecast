!> \file hype_wbout.f90
!> Contains module hype_waterbalance.

!>Handles specific output of waterbalance for HYPE model
MODULE HYPE_WATERBALANCE

  !Copyright 2014,2016-2107 SMHI
  !
  !This file is part of HYPE.

  !HYPE is free software: you can redistribute it and/or modify it under
  !the terms of the Lesser GNU General Public License as published by
  !the Free Software Foundation, either version 3 of the License, or (at
  !your option) any later version. HYPE is distributed in the hope that
  !it will be useful, but WITHOUT ANY WARRANTY; without even the implied
  !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
  !the Lesser GNU General Public License for more details. You should
  !have received a copy of the Lesser GNU General Public License along
  !with HYPE. If not, see <http://www.gnu.org/licenses/>.
  !-----------------------------------------------------------------------------------------

  USE LibDate
  USE READWRITE_ROUTINES, ONLY : write_mathsep
  USE WORLDVAR, ONLY : fileunit_free, &
                       fileunit_get, &
                       get_seq_filename, &
                       writematlab, &   !be used here?
                       timeformat       !be used here?
  !Also uses MODVAR                       

  IMPLICIT NONE
  PRIVATE :: num_wbflows,num_wbstores,num_wbirrflows,num_wbfpflows,   &
             fname_wbflows,fname_wbstores,fname_wbirrflows,fname_wbfpflows,   &
             funit_wbflows,funit_wbstores,funit_wbirrflows,funit_wbfpflows,   &
             nwbirrbasin,printindex_wbirrflows, &
             nwbfpbasin,printindex_wbfpflows

!> \name Water balance flow variable indices
!> \{
  INTEGER,PARAMETER :: w_sfallsnow = 1   !<Index for snowfall on snow (on land)
  INTEGER,PARAMETER :: w_smeltsl1  = 2   !<Index for snow melt to infiltration
  INTEGER,PARAMETER :: w_smeltsr   = 3   !<Index for snow melt to surface runoff
  INTEGER,PARAMETER :: w_smeltmp1  = 4   !<Index for snow melt to macropore to soil layer 1
  INTEGER,PARAMETER :: w_smeltmp2  = 5   !<Index for snow melt to macropore to soil layer 2
  INTEGER,PARAMETER :: w_smeltmp3  = 6   !<Index for snow melt to macropore to soil layer 3
  INTEGER,PARAMETER :: w_stoice    = 7   !<Index for snow incorporated in glacier
  INTEGER,PARAMETER :: w_precglac  = 8   !<Index for precipitation on glacier
  INTEGER,PARAMETER :: w_gmeltsl1  = 9   !<Index for glacier melt to infiltration
  INTEGER,PARAMETER :: w_gmeltsr   = 10  !<Index for glacier melt to surface runoff
  INTEGER,PARAMETER :: w_gmeltmp1  = 11  !<Index for glacier melt to macropore to soil layer 1
  INTEGER,PARAMETER :: w_gmeltmp2  = 12  !<Index for glacier melt to macropore to soil layer 2
  INTEGER,PARAMETER :: w_gmeltmp3  = 13  !<Index for glacier melt to macropore to soil layer 3
  INTEGER,PARAMETER :: w_infrain   = 14  !<Index for infiltration of rain
  INTEGER,PARAMETER :: w_rainsr    = 55  !<Index for rain to surface runoff
  INTEGER,PARAMETER :: w_rainmp1   = 15  !<Index for rain to macropore to soil layer 1
  INTEGER,PARAMETER :: w_rainmp2   = 16  !<Index for rain to macropore to soil layer 2
  INTEGER,PARAMETER :: w_rainmp3   = 17  !<Index for rain to macropore to soil layer 3
  INTEGER,PARAMETER :: w_evap1     = 18  !<Index for evapotranspiration from soil layer 1
  INTEGER,PARAMETER :: w_evap2     = 19  !<Index for evapotranspiration from soil layer 2
  INTEGER,PARAMETER :: w_surfrf    = 20  !<Index for surface runoff (saturated overland flow)
  INTEGER,PARAMETER :: w_tile1     = 21  !<Index for tile drainage pipes runoff from soil layer 1
  INTEGER,PARAMETER :: w_tile2     = 22  !<Index for tile drainage pipes runoff from soil layer 2
  INTEGER,PARAMETER :: w_tile3     = 23  !<Index for tile drainage pipes runoff from soil layer 3
  INTEGER,PARAMETER :: w_gwrunf1   = 24  !<Index for groundwater runoff from soil layer 1
  INTEGER,PARAMETER :: w_gwrunf2   = 25  !<Index for groundwater runoff from soil layer 2
  INTEGER,PARAMETER :: w_gwrunf3   = 26  !<Index for groundwater runoff from soil layer 3
  INTEGER,PARAMETER :: w_perc1     = 27  !<Index for percolation from soil layer 1
  INTEGER,PARAMETER :: w_perc2     = 28  !<Index for percolation from soil layer 2
  INTEGER,PARAMETER :: w_upwell1   = 29  !<Index for upwelling to soil layer 1
  INTEGER,PARAMETER :: w_upwell2   = 30  !<Index for upwelling to soil layer 2
  INTEGER,PARAMETER :: w_rural1    = 31  !<Index for flow from rural load to soil layer 1
  INTEGER,PARAMETER :: w_rural2    = 32  !<Index for flow from rural load to soil layer 2
  INTEGER,PARAMETER :: w_rural3    = 33  !<Index for flow from rural load to soil layer 3
  INTEGER,PARAMETER :: w_rgrwto1   = 34  !<Index for regional groundwater flow to soil layer 1
  INTEGER,PARAMETER :: w_rgrwto2   = 35  !<Index for regional groundwater flow to soil layer 2
  INTEGER,PARAMETER :: w_rgrwto3   = 36  !<Index for regional groundwater flow to soil layer 3
  INTEGER,PARAMETER :: w_piriver   = 37  !<Index for precipitation on local stream
  INTEGER,PARAMETER :: w_pilake    = 38  !<Index for precipitation on internal lake
  INTEGER,PARAMETER :: w_pmriver   = 39  !<Index for precipitation on main river
  INTEGER,PARAMETER :: w_polake    = 40  !<Index for precipitation on outlet lake
  INTEGER,PARAMETER :: w_eiriver   = 41  !<Index for evaporation of local stream
  INTEGER,PARAMETER :: w_eilake    = 42  !<Index for evaporation of internal lake
  INTEGER,PARAMETER :: w_emriver   = 43  !<Index for evaporation of main river
  INTEGER,PARAMETER :: w_eolake    = 44  !<Index for evaporation of outlet lake
  INTEGER,PARAMETER :: w_irtoil    = 45  !<Index for flow from local stream to internal lake
  INTEGER,PARAMETER :: w_irtomr    = 46  !<Index for flow from local stream to main river
  INTEGER,PARAMETER :: w_iltomr    = 47  !<Index for flow from internal lake to main river
  INTEGER,PARAMETER :: w_mrtool    = 48  !<Index for flow from main river to outlet lake
  INTEGER,PARAMETER :: w_oltomb    = 49  !<Index for flow from outlet lake to main downstream subbasin
  INTEGER,PARAMETER :: w_oltob     = 50  !<Index for flow from outlet lake to branch downstream
  INTEGER,PARAMETER :: w_rural4    = 51  !<Index for rural flow to local river
  INTEGER,PARAMETER :: w_pstomr    = 52  !<Index for point source flow to main river
  INTEGER,PARAMETER :: w_abstmr    = 53  !<Index for abstraction from main river
  INTEGER,PARAMETER :: w_abstol    = 54  !<Index for abstraction from outlet lake
  INTEGER,PARAMETER :: w_rgrwtool  = 56  !<Index for regional groundwater flow to outlet lake
  INTEGER,PARAMETER :: w_rgrwof1   = 57  !<Index for regional groundwater flow from soil layer 1
  INTEGER,PARAMETER :: w_rgrwof2   = 58  !<Index for regional groundwater flow from soil layer 2
  INTEGER,PARAMETER :: w_rgrwof3   = 59  !<Index for regional groundwater flow from soil layer 3
  INTEGER,PARAMETER :: w_rgrwtoos  = 60  !<Index for regional groundwater flow to outside model domain (grwdown=0)
  INTEGER,PARAMETER :: w_rgrwofmr  = 61  !<Index for regional groundwater flow from main river
  INTEGER,PARAMETER :: w_rgrwtomr  = 62  !<Index for regional groundwater flow to main river
  INTEGER,PARAMETER :: w_evap3     = 63  !<Index for evaporation from snow
  INTEGER,PARAMETER :: w_evap4     = 64  !<Index for evaporation from glacier
  INTEGER,PARAMETER :: num_wbflows = 64  !<Number of water balance flows
!>\}

!> \name Water balance store variable indices
!> \{
  INTEGER,PARAMETER :: w_snow     = 1   !<Index for snow
  INTEGER,PARAMETER :: w_glacier  = 2   !<Index for glacier
  INTEGER,PARAMETER :: w_soil1    = 3   !<Index for soil layer 1 (uppermost)
  INTEGER,PARAMETER :: w_soil2    = 4   !<Index for soil layer 2
  INTEGER,PARAMETER :: w_soil3    = 5   !<Index for soil layer 3 (deepest)
  INTEGER,PARAMETER :: w_iriver   = 6   !<Index for local stream
  INTEGER,PARAMETER :: w_ilake    = 7   !<Index for internal lake
  INTEGER,PARAMETER :: w_mriver   = 8   !<Index for main river
  INTEGER,PARAMETER :: w_olake    = 9   !<Index for outlet lake
  INTEGER,PARAMETER :: w_irrcanal = 10  !<Index for irrigation canal
  INTEGER,PARAMETER :: w_aquifer  = 11  !<Index for aquifer
  INTEGER,PARAMETER :: w_riverplain = 12  !<Index for main river floodplain
  INTEGER,PARAMETER :: w_lakeplain  = 13  !<Index for outlet lake floodplain
  INTEGER,PARAMETER :: num_wbstores = 13  !<Number of water balance stores
!>\}

!> \name Water balance selected subbasin flow variable indices (irrigation)
!> \{
  INTEGER,PARAMETER :: w_apply1    = 1   !<Index for irrigation applied to soil layer 1
  INTEGER,PARAMETER :: w_apply2    = 2   !<Index for irrigation applied to soil layer 2
  INTEGER,PARAMETER :: w_wdfromil  = 3   !<Index for withdrawal of water from internal lake
  INTEGER,PARAMETER :: w_wdfromdg  = 4   !<Index for withdrawal of water from deep groundwater
  INTEGER,PARAMETER :: w_wdfromol  = 5   !<Index for withdrawal of water from outlet lake as local source
  INTEGER,PARAMETER :: w_wdfrommr  = 6   !<Index for withdrawal of water from main river volume and inflow as local source
  INTEGER,PARAMETER :: w_evapirrc  = 7   !<Index for evaporation losses of (local) irrigation canal
  INTEGER,PARAMETER :: w_wdregol   = 8   !<Index for withdrawal of water from regional source (outlet lake) of this subbasin
  INTEGER,PARAMETER :: w_wdregmr   = 9   !<Index for withdrawal of water from regional source (main river volume and inflow) of this subbasin
  INTEGER,PARAMETER :: w_evapregol = 10  !<Index for evaporation losses at withdrawal of water from outlet lake (regional source)
  INTEGER,PARAMETER :: w_evapregmr = 11  !<Index for evaporation losses at withdrawal of water from river volume and inflow (regional source)
  INTEGER,PARAMETER :: w_wdoutside = 12  !<Index for withdrawal of water from outside model domain to cover all demands
  INTEGER,PARAMETER :: w_rgrwtoir  = 13  !<Index for withdrawal of water from aquifer
  INTEGER,PARAMETER :: num_wbirrflows = 13  !<Number of water balance flows regarding irrigation
!>\}

!> \name Water balance selected subbasin flow variable indices (floodplain)
!> \{
  INTEGER,PARAMETER :: w_pmrfp     = 1   !<Index for precipitation on main river floodplain
  INTEGER,PARAMETER :: w_polfp     = 2   !<Index for precipitation on outlet lake floodplain
  INTEGER,PARAMETER :: w_emrfp     = 3   !<Index for evaporation of main river floodplain
  INTEGER,PARAMETER :: w_eolfp     = 4   !<Index for evaporation of outlet lake floodplain
  INTEGER,PARAMETER :: w_infmrfp   = 5   !<Index for infiltration from main river floodplain to soil layer 1
  INTEGER,PARAMETER :: w_infolfp   = 6   !<Index for infiltration from outlet lake floodplain to soil layer 1
  INTEGER,PARAMETER :: w_smtsrmrfp = 7   !<Index for snow melt to surface runoff to main river floodplain
  INTEGER,PARAMETER :: w_smtsrolfp = 8   !<Index for snow melt to surface runoff to outlet lake floodplain
  INTEGER,PARAMETER :: w_rtsrmrfp  = 9   !<Index for rain to surface runoff to main river floodplain
  INTEGER,PARAMETER :: w_rtsrolfp  = 10  !<Index for rain to surface runoff to outlet lake floodplain
  INTEGER,PARAMETER :: w_srftmrfp  = 11  !<Index for surface runoff (saturated overland flow) to main river floodplain
  INTEGER,PARAMETER :: w_srftolfp  = 12  !<Index for surface runoff (saturated overland flow) to outlet lake floodplain
  INTEGER,PARAMETER :: w_grf1mrfp  = 13  !<Index for groundwater runoff from soil layer 1 to main river floodplain
  INTEGER,PARAMETER :: w_grf2mrfp  = 14  !<Index for groundwater runoff from soil layer 2 to main river floodplain
  INTEGER,PARAMETER :: w_grf3mrfp  = 15  !<Index for groundwater runoff from soil layer 3 to main river floodplain
  INTEGER,PARAMETER :: w_grf1olfp  = 16  !<Index for groundwater runoff from soil layer 1 to outlet lake floodplain
  INTEGER,PARAMETER :: w_grf2olfp  = 17  !<Index for groundwater runoff from soil layer 2 to outlet lake floodplain
  INTEGER,PARAMETER :: w_grf3olfp  = 18  !<Index for groundwater runoff from soil layer 3 to outlet lake floodplain
  INTEGER,PARAMETER :: w_stomrfp   = 19  !<Index for snow incorporated in main river floodplain
  INTEGER,PARAMETER :: w_stoolfp   = 20  !<Index for snow incorporated in outlet lake floodplain
  INTEGER,PARAMETER :: w_fptomr    = 21  !<Index for flow from floodplain to main river 
  INTEGER,PARAMETER :: w_mrtofp    = 22  !<Index for flow from main river to floodplain 
  INTEGER,PARAMETER :: w_fptool    = 23  !<Index for flow from floodplain to outlet lake 
  INTEGER,PARAMETER :: w_oltofp    = 24  !<Index for flow from outlet lake to floodplain 
  INTEGER,PARAMETER :: w_trf1mrfp  = 25  !<Index for tile runoff to main river floodplain
  INTEGER,PARAMETER :: w_trf2mrfp  = 26  !<Index for tile runoff to main river floodplain
  INTEGER,PARAMETER :: w_trf3mrfp  = 27  !<Index for tile runoff to main river floodplain
  INTEGER,PARAMETER :: w_trf1olfp  = 28  !<Index for tile runoff to outlet lake floodplain
  INTEGER,PARAMETER :: w_trf2olfp  = 29  !<Index for tile runoff to outlet lake floodplain
  INTEGER,PARAMETER :: w_trf3olfp  = 30  !<Index for tile runoff to outlet lake floodplain
  INTEGER,PARAMETER :: w_rgrwofmrfp  = 31  !<Index for regional groundwater flow leaving main river floodplain
  INTEGER,PARAMETER :: num_wbfpflows = 31  !<Number of water balance flows regarding floodplains
!>\}

!>\name Variables for waterbalance output
!>\{
! !Variables to hold all waterbalance values. 
! !Variables for file names and file units, and for which subbasins has irrigation.
  INTEGER :: funit_wbstores(num_wbstores)          !<fileunits waterbalance stores
  INTEGER :: funit_wbflows(num_wbflows)            !<fileunits waterbalance flows
  INTEGER :: funit_wbirrflows(num_wbirrflows)      !<fileunits waterbalance irrigation flows
  INTEGER :: funit_wbfpflows(num_wbfpflows)        !<fileunits waterbalance floodplain flows
  CHARACTER(LEN=20) :: fname_wbstores(num_wbstores) = (/ &
      'WBs_snow.txt      ', &
      'WBs_glacier.txt   ','WBs_soillayer1.txt', &
      'WBs_soillayer2.txt','WBs_soillayer3.txt', &
      'WBs_lstream.txt   ','WBs_ilake.txt     ', &
      'WBs_mriver.txt    ','WBs_olake.txt     ', &
      'WBs_irrcanal.txt  ','WBs_aquifer.txt   ', &
      'WBs_riverplain.txt','WBs_lakeplain.txt '/)   !<filenames waterbalance stores
  CHARACTER(LEN=50) :: fname_wbflows(num_wbflows) = (/ &
      'WBf_snowfall__snow.txt                          ', &
      'WBf_snowmelt_snow_soillayer1.txt                ', &
      'WBf_snowmelt_surfacerunoff_snow_lstream.txt     ', &
      'WBf_snowmelt_via_macropore_snow_soillayer1.txt  ', &
      'WBf_snowmelt_via_macropore_snow_soillayer2.txt  ', &
      'WBf_snowmelt_via_macropore_snow_soillayer3.txt  ', &
      'WBf_growingice_snow_glacier.txt                 ', &
      'WBf_precipitation__glacier.txt                  ', &
      'WBf_melt_glacier_soillayer1.txt                 ', &
      'WBf_melt_surfacerunoff_glacier_lstream.txt      ', &
      'WBf_melt_via_macropore_glacier_sollayer1.txt    ', &
      'WBf_melt_via_macropore_glacier_sollayer2.txt    ', &
      'WBf_melt_via_macropore_glacier_sollayer3.txt    ', &
      'WBf_rain__soillayer1.txt                        ', &
      'WBf_rain_via_macropore__soillayer1.txt          ', &
      'WBf_rain_via_macropore__soillayer2.txt          ', &
      'WBf_rain_via_macropore__soillayer3.txt          ', &
      'WBf_evaporation_soillayer1_.txt                 ', &
      'WBf_evaporation_soillayer2_.txt                 ', &
      'WBf_satsurfaceflow_soillayer1_lstream.txt       ', &
      'WBf_tilerunoff_soillayer1_lstream.txt           ', &
      'WBf_tilerunoff_soillayer2_lstream.txt           ', &
      'WBf_tilerunoff_soillayer3_lstream.txt           ', &
      'WBf_soilrunoff_soillayer1_lstream.txt           ', &
      'WBf_soilrunoff_soillayer2_lstream.txt           ', &
      'WBf_soilrunoff_soillayer3_lstream.txt           ', &
      'WBf_percolation_soillayer1_soillayer2.txt       ', &
      'WBf_percolation_soillayer2_soillayer3.txt       ', &
      'WBf_upwell_soillayer2_soillayer1.txt            ', &
      'WBf_upwell_soillayer3_soillayer2.txt            ', &
      'WBf_ruralflow__soillayer1.txt                   ', &
      'WBf_ruralflow__soillayer2.txt                   ', &
      'WBf_ruralflow__soillayer3.txt                   ', &
      'WBf_regionalgroundwater_reservoir_soillayer1.txt', &
      'WBf_regionalgroundwater_reservoir_soillayer2.txt', &
      'WBf_regionalgroundwater_reservoir_soillayer3.txt', &
      'WBf_precipitation__lstream.txt                  ', &
      'WBf_precipitation__ilake.txt                    ', &
      'WBf_precipitation__mriver.txt                   ', &
      'WBf_precipitation__olake.txt                    ', &
      'WBf_evaporation_lstream_.txt                    ', &
      'WBf_evaporation_ilake_.txt                      ', &
      'WBf_evaporation_mriver_.txt                     ', &
      'WBf_evaporation_olake_.txt                      ', &
      'WBf_flow_lstream_ilake.txt                      ', &
      'WBf_flow_lstream_mriver.txt                     ', &
      'WBf_flow_ilake_mriver.txt                       ', &
      'WBf_flow_mriver_olake.txt                       ', &
      'WBf_flow_olake_mriver_maindownstream.txt        ', &
      'WBf_flow_olake_mriver_branchdownstream.txt      ', &
      'WBf_ruralflow__lstream.txt                      ', &
      'WBf_pointsource__mriver.txt                     ', &
      'WBf_abstraction_mriver_.txt                     ', &
      'WBf_abstraction_olake_.txt                      ', &
      'WBf_rain_surfacerunoff__lstream.txt             ', &
      'WBf_regionalgroundwater_reservoir_olake.txt     ', &
      'WBf_regionalgroundwater_soillayer1_reservoir.txt', &
      'WBf_regionalgroundwater_soillayer2_reservoir.txt', &
      'WBf_regionalgroundwater_soillayer3_reservoir.txt', &
      'WBf_regionalgroundwater_reservoir_.txt          ', &
      'WBf_regionalgroundwater_mriver_reservoir.txt    ', &
      'WBf_regionalgroundwater_reservoir_mriver.txt    ', &
      'WBf_evaporation_snow_.txt                       ', &
      'WBf_evaporation_glacier_.txt                    '/) !<filenames waterbalance flows
  CHARACTER(LEN=50) :: fname_wbirrflows(num_wbirrflows) = (/ &
      'WBfs_apply_irrcanal_soillayer1.txt             ', &
      'WBfs_apply_irrcanal_soillayer2.txt             ', &
      'WBfs_local_withdraw_ilake_irrcanal.txt         ', &
      'WBfs_local_withdraw_deepgrw_irrcanal.txt       ', &
      'WBfs_local_withdraw_olake_irrcanal.txt         ', &
      'WBfs_local_withdraw_mriver_irrcanal.txt        ', &
      'WBfs_evaporation_irrcanal_.txt                 ', &
      'WBfs_regsrc_withdraw_olake_irrcanal.txt        ', &
      'WBfs_regsrc_withdraw_mriver_irrcanal.txt       ', &
      'WBfs_regsrc_evaporation_olake_.txt             ', &
      'WBfs_regsrc_evaporation_mriver_.txt            ', &
      'WBfs_unlimited_withdraw__irrcanal.txt          ', &
      'WBfs_regionalgroundwater_reservoir_irrcanal.txt'/) !<filenames waterbalance irrigation flows
  CHARACTER(LEN=50) :: fname_wbfpflows(num_wbfpflows) = (/ &
      'WBff_precipitation__mriver_floodplain.txt      ', &
      'WBff_precipitation__olake_floodplain.txt       ', &
      'WBff_evaporation_mriver_floodplain_.txt        ', &
      'WBff_evaporation_olake_floodplain_.txt         ', &
      'WBff_infilt_mriver_floodplain_soillayer1.txt   ', &
      'WBff_infilt_olake_floodplain_soillayer1.txt    ', &
      'WBff_melt_surfrunoff_snow_mriver_floodplain.txt', &
      'WBff_melt_surfrunoff_snow_olake_floodplain.txt ', &
      'WBff_rain_surfrunoff_mriver_floodplain.txt     ', &
      'WBff_rain_surfrunoff_olake_floodplain.txt      ', &
      'WBff_satsfflow_soillayer1_mriver_floodplain.txt', &
      'WBff_satsfflow_soillayer1_olake_floodplain.txt ', &
      'WBff_soilrunf_soillayer1_mriver_floodplain.txt ', &
      'WBff_soilrunf_soillayer2_mriver_floodplain.txt ', &
      'WBff_soilrunf_soillayer3_mriver_floodplain.txt ', &
      'WBff_soilrunf_soillayer1_olake_floodplain.txt  ', &
      'WBff_soilrunf_soillayer2_olake_floodplain.txt  ', &
      'WBff_soilrunf_soillayer3_olake_floodplain.txt  ', &
      'WBff_flooded_snow_mriver_floodplain.txt        ', &
      'WBff_flooded_snow_olake_floodplain.txt         ', &
      'WBff_flow_mriver_floodplain_mriver.txt         ', &
      'WBff_flow_mriver_mriver_floodplain.txt         ', &
      'WBff_flow_olake_floodplain_olake.txt           ', &
      'WBff_flow_olake_olake_floodplain.txt           ', &
      'WBff_tilerunf_soillayer1_mriver_floodplain.txt ', &
      'WBff_tilerunf_soillayer2_mriver_floodplain.txt ', &
      'WBff_tilerunf_soillayer3_mriver_floodplain.txt ', &
      'WBff_tilerunf_soillayer1_olake_floodplain.txt  ', &
      'WBff_tilerunf_soillayer2_olake_floodplain.txt  ', &
      'WBff_tilerunf_soillayer3_olake_floodplain.txt  ', &
      'WBff_reggroundwater_mr_floodplain_reservoir.txt'/) !<filenames waterbalance floodplain flows
  REAL, ALLOCATABLE, PUBLIC :: wbflows(:,:)     !<Current time step water balance flows (num_wbflows,nsub)
  REAL, ALLOCATABLE, PUBLIC :: wbirrflows(:,:)  !<Current time step water balance irrigation flows
  REAL, ALLOCATABLE, PUBLIC :: wbfpflows(:,:)   !<Current time step water balance floodplain flows
  REAL, ALLOCATABLE, PUBLIC :: wbstores(:,:)    !<Current time step water balance stores (at the end of the time step)
  INTEGER :: nwbirrbasin                        !<Number of subbasins with irrigation flows
  INTEGER :: nwbfpbasin                         !<Number of subbasins with floodplains
  INTEGER, ALLOCATABLE :: printindex_wbirrflows(:)              !<index of subbasins with irrigation, to be printed
  INTEGER, ALLOCATABLE :: printindex_wbfpflows(:)     !<index of subbasins with floodplains, to be printed
!>\}

CONTAINS

  !>Open and write heading for files of waterbalance flows and stores
  !>Calculate which subbasins need to print irrigation flows and floodplain flows
  !----------------------------------------------------------------------
  SUBROUTINE prepare_waterbalance_files(dir,n,na,subid)

    USE MODVAR, ONLY : irrigationsystem,floodindex

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: subid(n)     !<subid of subbasins

    !Local variables
    INTEGER i,j   !loop index
    CHARACTER(LEN=50) filename
    INTEGER lt,lout,lout2,lout3,lout4
    INTEGER dimirr    !size of irrigation data
    INTEGER lirrindex(n)
    CHARACTER (LEN=16)   t
    CHARACTER (LEN=3400000) outtxt,outtxt2,outtxt3,outtxt4

    !Calculate subbasins with irrigation flows (both irrigated areas and regional sources should be in MgmtData)
    nwbirrbasin = 0
    IF(ALLOCATED(irrigationsystem))THEN
      lirrindex = 0
      dimirr = SIZE(irrigationsystem)
      DO i = 1,n
        DO j = 1, dimirr
          IF(subid(i)==irrigationsystem(j)%subid)THEN
            lirrindex(i) = j           !this basin (i) has irrigation found on row j in the irrigationsystem object
            nwbirrbasin = nwbirrbasin + 1   !count number of subbasin with irrigation applied
          ENDIF
        ENDDO
      ENDDO
      ALLOCATE(printindex_wbirrflows(nwbirrbasin))
      j = 1
      DO i = 1,n
        IF(lirrindex(i)>0)THEN
          printindex_wbirrflows(j) = i
          j = j + 1
        ENDIF
      ENDDO
      IF(j-1/=nwbirrbasin) WRITE(6,*) 'ERROR: number of irrbasin for wb output'   !test
    ENDIF
    
    !Calculate number of subbasins with floodplains
    nwbfpbasin = 0
    IF(ALLOCATED(floodindex))THEN
      DO i = 1,n
        IF(floodindex(i)>0)THEN
          nwbfpbasin = nwbfpbasin + 1   !count number of subbasin with floodplains
        ENDIF
      ENDDO
      ALLOCATE(printindex_wbfpflows(nwbfpbasin))
      j = 1
      DO i = 1,n
        IF(floodindex(i)>0)THEN
          printindex_wbfpflows(j) = i
          j = j + 1
        ENDIF
      ENDDO
      IF(j-1/=nwbfpbasin) WRITE(6,*) 'ERROR: number of floodplain basin for wb output'   !test
    ENDIF

    !Prepare headings (all subbasins)
    outtxt(1:5) = 'DATE'//CHAR(9)
    lout = 5
    DO j = 1,n
      WRITE(t,'(i16)') subid(j)
      t = ADJUSTL(t)
      lt = LEN_TRIM(t)
      outtxt(lout+1:lout+lt) = t(1:lt)
      lout = lout+lt
      IF(j < n) THEN
        outtxt(lout+1:lout+1) = CHAR(9)    !Also the last one necessary
        lout = lout+1                      !for reading in "free format"
      ENDIF
    ENDDO

    !Prepare headings (irrigation subbasins)
    IF(nwbirrbasin>0)THEN
      outtxt2(1:5) = 'DATE'//CHAR(9)
      lout2 = 5
      DO j = 1,n
        IF(lirrindex(j)==0) CYCLE    !skip subbasins without irrigation
        WRITE(t,'(i16)') subid(j)
        t = ADJUSTL(t)
        lt = LEN_TRIM(t)
        outtxt2(lout2+1:lout2+lt) = t(1:lt)
        lout2 = lout2+lt
        IF(j < n) THEN
          outtxt2(lout2+1:lout2+1) = CHAR(9)    !Also the last one necessary
          lout2 = lout2+1                      !for reading in "free format"
        ENDIF
      ENDDO
    ENDIF

    !Prepare headings (aquifers)
    IF(na>0)THEN
      outtxt3(1:5) = 'DATE'//CHAR(9)
      lout3 = 5
      DO j = 1,na
        WRITE(t,'(i16)') j
        t = ADJUSTL(t)
        lt = LEN_TRIM(t)
        outtxt3(lout3+1:lout3+lt) = t(1:lt)
        lout3 = lout3+lt
        IF(j < na) THEN
          outtxt3(lout3+1:lout3+1) = CHAR(9)    !Also the last one necessary
          lout3 = lout3+1                      !for reading in "free format"
        ENDIF
      ENDDO
    ENDIF

    !Prepare headings (floodplain subbasins)
    IF(nwbfpbasin>0)THEN
      outtxt4(1:5) = 'DATE'//CHAR(9)
      lout4 = 5
      DO j = 1,n
        IF(floodindex(j)==0) CYCLE    !skip subbasins without floodplain
        WRITE(t,'(i16)') subid(j)
        t = ADJUSTL(t)
        lt = LEN_TRIM(t)
        outtxt4(lout4+1:lout4+lt) = t(1:lt)
        lout4 = lout4+lt
        IF(j < n) THEN
          outtxt4(lout4+1:lout4+1) = CHAR(9)    !Also the last one necessary
          lout4 = lout4+1                      !for reading in "free format"
        ENDIF
      ENDDO
    ENDIF

    !Open water balance flow files for writing and write heading
    DO i = 1,num_wbflows
      funit_wbflows(i) = fileunit_get()
      filename = fname_wbflows(i)
  !    CALL get_seq_filename(filename)    
      OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=funit_wbflows(i),STATUS='unknown',FORM='formatted')

      !Write headings
      IF(i==w_rgrwtoos.AND.na>0)THEN
        WRITE(funit_wbflows(i),'(a)') outtxt3(1:lout3)
      ELSE
        WRITE(funit_wbflows(i),'(a)') outtxt(1:lout)
      ENDIF
    ENDDO
  
    !Open water balance storage files for writing and write heading
    DO i = 1,num_wbstores
      filename = fname_wbstores(i)
      IF(i==w_snow.OR.i==w_glacier.OR.i==w_soil1.OR.i==w_soil2.OR.i==w_soil3.OR.  &
         i==w_iriver.OR.i==w_ilake.OR.i==w_mriver.OR.i==w_olake.OR. &
         (i==w_irrcanal.AND.nwbirrbasin>0).OR.  &
         (i==w_riverplain.AND.nwbfpbasin>0).OR.  &
         (i==w_lakeplain.AND.nwbfpbasin>0).OR.  &
         (i==w_aquifer.AND.na>0))THEN
        funit_wbstores(i) = fileunit_get()
        OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=funit_wbstores(i),STATUS='unknown',FORM='formatted')
      ENDIF

      !Write headings
      IF(i==w_irrcanal)THEN
        IF(nwbirrbasin>0) WRITE(funit_wbstores(i),'(a)') outtxt2(1:lout2)
      ELSEIF(i==w_aquifer)THEN
        IF(na>0) WRITE(funit_wbstores(i),'(a)') outtxt3(1:lout3)
      ELSEIF(i==w_riverplain.OR.i==w_lakeplain)THEN
        IF(nwbfpbasin>0) WRITE(funit_wbstores(i),'(a)') outtxt4(1:lout4)
      ELSE
        WRITE(funit_wbstores(i),'(a)') outtxt(1:lout)    
      ENDIF
    ENDDO

    !Open water balance selected flow files for writing and write heading
    IF(nwbirrbasin>0)THEN
      DO i = 1,num_wbirrflows
        funit_wbirrflows(i) = fileunit_get()
        filename = fname_wbirrflows(i)
    !    CALL get_seq_filename(filename)
        OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=funit_wbirrflows(i),STATUS='unknown',FORM='formatted')

        !Write headings
        WRITE(funit_wbirrflows(i),'(a)') outtxt2(1:lout2)
      ENDDO
    ENDIF
    IF(nwbfpbasin>0)THEN
      DO i = 1,num_wbfpflows
        funit_wbfpflows(i) = fileunit_get()
        filename = fname_wbfpflows(i)
        OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=funit_wbfpflows(i),STATUS='unknown',FORM='formatted')

        !Write headings
        WRITE(funit_wbfpflows(i),'(a)') outtxt4(1:lout4)
      ENDDO
    ENDIF
    
  END SUBROUTINE prepare_waterbalance_files

  !>Allocate variables for water balance calculations
  !----------------------------------------------------------------------
  SUBROUTINE initiate_waterbalance_output(n)

    !Argument declaration
    INTEGER, INTENT(IN) :: n            !<Number of subbasins

    !Local variables

    IF(.NOT.ALLOCATED(wbflows)) ALLOCATE(wbflows(num_wbflows,n))
    IF(.NOT.ALLOCATED(wbstores)) ALLOCATE(wbstores(num_wbstores,n))
    IF(.NOT.ALLOCATED(wbirrflows)) ALLOCATE(wbirrflows(num_wbirrflows,n))
    IF(.NOT.ALLOCATED(wbfpflows)) ALLOCATE(wbfpflows(num_wbfpflows,n))
  
  END SUBROUTINE initiate_waterbalance_output

  !>Write one timestep of waterbalance flows and stores
  !----------------------------------------------------------------------
  SUBROUTINE print_waterbalance_timestep(n,na,date)

    !Argument declaration
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    TYPE(DateType), INTENT(IN) :: date  !<Current time in DateType format

    !Local variables
    INTEGER, PARAMETER :: ndec = 9  !? 7 siffror signifikanta för REAL, skriv ut 8?
    INTEGER i                     !loop-variable
    CHARACTER(LEN=16) textin

    !Format timestamp
    IF(writematlab)THEN
      IF(timeformat==0)THEN
        CALL format_date(date,'yyyymmdd',textin)
      ELSE
        CALL format_date(date,'yyyymmddHHMM',textin)
      ENDIF
    ELSE
      IF(timeformat==0)THEN
        CALL format_date(date,'yyyy-mm-dd',textin)
      ELSE
        CALL format_date(date,'yyyy-mm-dd HH:MM',textin)
      ENDIF
    ENDIF

    !Write water balance flows to files
    DO i = 1,num_wbflows
      IF(i==w_rgrwtoos.AND.na>0)THEN
        CALL write_mathsep(funit_wbflows(i),na,wbflows(i,1:na),ndec,textin,CHAR(9))
      ELSE
        CALL write_mathsep(funit_wbflows(i),n,wbflows(i,:),ndec,textin,CHAR(9))
      ENDIF
    ENDDO

    !Write water balance stores to files
    DO i = 1,num_wbstores
      IF(i==w_irrcanal)THEN
        IF(nwbirrbasin>0) CALL write_mathsep(funit_wbstores(i),nwbirrbasin,wbstores(i,printindex_wbirrflows(:)),ndec,textin,CHAR(9))
      ELSEIF(i==w_riverplain)THEN
        IF(nwbfpbasin>0) CALL write_mathsep(funit_wbstores(i),nwbfpbasin,wbstores(i,printindex_wbfpflows(:)),ndec,textin,CHAR(9))
      ELSEIF(i==w_lakeplain)THEN
        IF(nwbfpbasin>0) CALL write_mathsep(funit_wbstores(i),nwbfpbasin,wbstores(i,printindex_wbfpflows(:)),ndec,textin,CHAR(9))
      ELSEIF(i==w_aquifer)THEN
        IF(na>0) CALL write_mathsep(funit_wbstores(i),na,wbstores(i,1:na),ndec,textin,CHAR(9))
      ELSE
        CALL write_mathsep(funit_wbstores(i),n,wbstores(i,:),ndec,textin,CHAR(9))
      ENDIF
    ENDDO

    !Write water balance irrigation flows to files
    IF(nwbirrbasin>0)THEN
      DO i = 1,num_wbirrflows
        CALL write_mathsep(funit_wbirrflows(i),nwbirrbasin,wbirrflows(i,printindex_wbirrflows(:)),ndec,textin,CHAR(9))
      ENDDO
    ENDIF

    !Write water balance floodplain flows to files
    IF(nwbfpbasin>0)THEN
      DO i = 1,num_wbfpflows
        CALL write_mathsep(funit_wbfpflows(i),nwbfpbasin,wbfpflows(i,printindex_wbfpflows(:)),ndec,textin,CHAR(9))
      ENDDO
    ENDIF

  END SUBROUTINE print_waterbalance_timestep

  !>Write one timestep of waterbalance stores
  !----------------------------------------------------------------------
  SUBROUTINE print_initial_waterbalance_stores(n,na)

    !Argument declaration
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
!    TYPE(DateType), INTENT(IN) :: bdate  !<Begin time in DateType format

    !Local constants and variables
    INTEGER, PARAMETER :: ndec = 9  !?
    INTEGER i                     !loop-variable
    CHARACTER(LEN=16) textin

    !Format timestamp
    textin = 'initial'

    !Write water balance stores to files
    DO i = 1,num_wbstores
      IF(i==w_irrcanal)THEN
        IF(nwbirrbasin>0) CALL write_mathsep(funit_wbstores(i),nwbirrbasin,wbstores(i,printindex_wbirrflows(:)),ndec,textin,CHAR(9))
      ELSEIF(i==w_riverplain)THEN
        IF(nwbfpbasin>0) CALL write_mathsep(funit_wbstores(i),nwbfpbasin,wbstores(i,printindex_wbfpflows(:)),ndec,textin,CHAR(9))
      ELSEIF(i==w_lakeplain)THEN
        IF(nwbfpbasin>0) CALL write_mathsep(funit_wbstores(i),nwbfpbasin,wbstores(i,printindex_wbfpflows(:)),ndec,textin,CHAR(9))
      ELSEIF(i==w_aquifer)THEN
        IF(na>0) CALL write_mathsep(funit_wbstores(i),na,wbstores(i,1:na),ndec,textin,CHAR(9))
      ELSE
        CALL write_mathsep(funit_wbstores(i),n,wbstores(i,:),ndec,textin,CHAR(9))
      ENDIF
    ENDDO


  END SUBROUTINE print_initial_waterbalance_stores

  !>Close files of waterbalance flows and stores
  !----------------------------------------------------------------------
  SUBROUTINE close_waterbalance_files(na)

    !Argument declarations
    INTEGER, INTENT(IN) :: na           !<Number of aquifers

    !Local variables
    INTEGER i   !loop index

    !Close water balance flow files
    DO i = 1,num_wbflows
      CLOSE(funit_wbflows(i))
      CALL fileunit_free(funit_wbflows(i))
    ENDDO
  
    !Close water balance storage files
    DO i = 1,num_wbstores
      IF(i==w_snow.OR.i==w_glacier.OR.i==w_soil1.OR.i==w_soil2.OR.i==w_soil3.OR.  &
         i==w_iriver.OR.i==w_ilake.OR.i==w_mriver.OR.i==w_olake.OR. &
         (i==w_irrcanal.AND.nwbirrbasin>0).OR.  &
         (i==w_riverplain.AND.nwbfpbasin>0).OR.  &
         (i==w_lakeplain.AND.nwbfpbasin>0).OR.  &
         (i==w_aquifer.AND.na>0))THEN
        CLOSE(funit_wbstores(i))
        CALL fileunit_free(funit_wbstores(i))
      ENDIF
    ENDDO

    !Close water balance selected flow files
    IF(nwbirrbasin>0)THEN
      DO i = 1,num_wbirrflows
        CLOSE(funit_wbirrflows(i))
        CALL fileunit_free(funit_wbirrflows(i))
      ENDDO
    ENDIF
    IF(nwbfpbasin>0)THEN
      DO i = 1,num_wbfpflows
        CLOSE(funit_wbfpflows(i))
        CALL fileunit_free(funit_wbfpflows(i))
      ENDDO
    ENDIF

  END SUBROUTINE close_waterbalance_files



END MODULE HYPE_WATERBALANCE
