!> \file soil_proc.f90
!> Contains module soil_processes.

!>Water processes in soil in HYPE and some more
MODULE SOIL_PROCESSES

  !Copyright 2012-2017 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !-----------------------------------------------------------------------------------------

  USE STATETYPE_MODULE, ONLY :soilstatetype,snowicestatetype
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water, &
                                          error_remove_water, &
                                          add_water
  USE NPC_SOIL_PROCESSES, ONLY : atmdep_in_loss, &
                                 doc_percolation_reduction, &
                                 onpp_percolation_reduction
  USE ATMOSPHERIC_PROCESSES, ONLY : deltasaturationpressure_function
  !Uses also modvar, hypevariables, hype_indata

  IMPLICIT NONE
  PRIVATE 
!Private procedures
!-------------------
! calculate_fractional_snowcover
! airpressure_elevationfunction
! psychrometric_constant
! get_current_petmodel
! calculate_weighted_temperature
  PUBLIC initiate_soil_water_state, &
         reinitiate_soil_depth, &
         initiate_soil_water, &
         calculate_snow, &
         latentheat_tempfunction, &
         calculate_potential_evaporation, &
         calculate_actual_soil_evapotranspiration, &
         calculate_tile_drainage, &
         calculate_soil_runoff, &
         infiltration, &
         add_infiltration, &
         flood_infiltration, &
         percolation, &
         add_macropore_flow, &
         calculate_groundwater_table, &
         calculate_snowmelt, &
         calculate_snowdepth, &
         calculate_soiltemp , &
         calculate_weighted_temperature, &
         calculate_frostdepth, &
         calculate_soil_moisture_deficit, &
         snowalbedo_function
         
  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(21)  !error message for location of remove_water call, 1-6 not used
  PARAMETER (errstring = (/'regional groundwater flow, soillayer 3        ',    & !1
                           'regional groundwater flow, soillayer 2        ',    &
                           'regional groundwater flow, soillayer 1        ',    &
                           'evapotranspiration lake, less than lake volume',    &
                           'evapotranspiration lake, more than lake volume',    &
                           'evapotranspiration lake, slowlake part used   ',       &
                           'evapotranspiration soil, soillayer 1          ',      & !7
                           'evapotranspiration soil, soillayer 2          ',      & !8
                           'runoff from soillayer 1, stream in layer 1    ',      &
                           'runoff from soillayer 2, stream in layer 2    ',      &
                           'runoff from soillayer 1, stream in layer 2    ',      &
                           'runoff from soillayer 3, stream in layer 3    ',      &
                           'runoff from soillayer 1, stream in layer 3    ',      &
                           'runoff from soillayer 2, stream in layer 3    ',      &   !14
                           'tile runoff, drainage pipe in soillayer 1     ',      &
                           'tile runoff, drainage pipe in soillayer 2     ',      &
                           'tile runoff, drainage pipe in soillayer 3     ',      &  !17
                           'percolation from soillayer 1                  ',      &  
                           'percolation from soillayer 2                  ',      &  !19
                           'snow melt                                     ',      &  !20
                           'snow sublimation                              '/))

CONTAINS

  !>\brief Initiate soil water state variables when no saved state exist.
  !!
  !>\b Reference ModelDescription Chapter Land routines (Basic assumptions)
  !---------------------------------------------------------------------------
  SUBROUTINE initiate_soil_water_state(soilstate)

    USE HYPEVARIABLES, ONLY : m_wcfc,m_wcwp,m_wcep,  &
                              m_wcfc1,m_wcfc2,m_wcfc3, &
                              m_wcwp1,m_wcwp2,m_wcwp3, &
                              m_wcep1,m_wcep2,m_wcep3, &
                              m_soilstretch
    USE MODVAR, ONLY : classdata,     &
                       nsub,          &
                       nclass,        &
                       maxsoillayers, &
                       soiliniwet,    &
                       soildepthstretch, &
                       soilthick,     &
                       soilpar,landpar
    
    !Argument declaration
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate   !<Soil states
    
    !Local variables
    INTEGER i,j           !loop-variables
    REAL :: iniwater(maxsoillayers,nclass)
    REAL :: tempsoilthick(maxsoillayers,nclass)

    !>\b Algoritm \n
    !>Calculate temporary soil layer thickness (mm)
    tempsoilthick = soilthick * 1000.
    IF(soildepthstretch)THEN
      DO j = 1,nclass
        IF(landpar(m_soilstretch,classdata(j)%luse)>0.)THEN
          tempsoilthick(2,j) = tempsoilthick(2,j) * landpar(m_soilstretch,classdata(j)%luse)
          tempsoilthick(3,j) = tempsoilthick(3,j) * landpar(m_soilstretch,classdata(j)%luse)
        ENDIF
      ENDDO
    ENDIF
    
    !>Calculate size of water storage in soil (wp,fc and ep) in mm
    iniwater = 0.
    DO j = 1,nclass
      IF(soilpar(m_wcfc1,classdata(j)%soil) > 0.)THEN                               !Field capacity (mm)
        iniwater(1,j)=soilpar(m_wcfc1,classdata(j)%soil) * tempsoilthick(1,j)       !First layer
        iniwater(2,j)=soilpar(m_wcfc2,classdata(j)%soil) * tempsoilthick(2,j)       !Second layer
        iniwater(3,j)=soilpar(m_wcfc3,classdata(j)%soil) * tempsoilthick(3,j)       !Third layer
      ELSE
        iniwater(:,j)=soilpar(m_wcfc,classdata(j)%soil) * tempsoilthick(:,j)        !All layers
      ENDIF
      IF(soilpar(m_wcwp1,classdata(j)%soil) > 0.)THEN                               !Wilting point (mm)
        iniwater(1,j)=iniwater(1,j)+soilpar(m_wcwp1,classdata(j)%soil) * tempsoilthick(1,j)        !First layer
        iniwater(2,j)=iniwater(2,j)+soilpar(m_wcwp2,classdata(j)%soil) * tempsoilthick(2,j)        !Second layer
        iniwater(3,j)=iniwater(3,j)+soilpar(m_wcwp3,classdata(j)%soil) * tempsoilthick(3,j)        !Third layer
      ELSE
        iniwater(:,j)=iniwater(:,j)+soilpar(m_wcwp,classdata(j)%soil) * tempsoilthick(:,j)         !All layers
      ENDIF
      IF(soiliniwet)THEN
        IF(soilpar(m_wcep1,classdata(j)%soil) > 0.)THEN                             !Effective porosity (mm)
          iniwater(1,j)=iniwater(1,j)+soilpar(m_wcep1,classdata(j)%soil) * tempsoilthick(1,j)     !First layer
          iniwater(2,j)=iniwater(2,j)+soilpar(m_wcep2,classdata(j)%soil) * tempsoilthick(2,j)     !Second layer
          iniwater(3,j)=iniwater(3,j)+soilpar(m_wcep3,classdata(j)%soil) * tempsoilthick(3,j)     !Third layer
        ELSE
          iniwater(:,j)=iniwater(:,j)+soilpar(m_wcep,classdata(j)%soil) * tempsoilthick(:,j)       !All layers
        ENDIF
      ENDIF
    ENDDO
    !>Initiate soil state to saturation or plant available water
    DO i = 1,nsub
      soilstate%water(:,:,i) = iniwater(:,:)
    ENDDO

  END SUBROUTINE initiate_soil_water_state

  !>\brief Reinitiate soil layer measures with soil stretch parameter.
  !!
  !>\b Reference ModelDescription Chapter Land routines (Basic assumptions)
  !>\b Consequences Module modvar variables soildepth and soilthick may change.
  !---------------------------------------------------------------------------
  SUBROUTINE reinitiate_soil_depth(nsl,nc,soillayerdepth,soillayerthick)
  
    USE HYPEVARIABLES, ONLY : m_soilstretch
    USE MODVAR, ONLY : classdata, &
                       soildepthstretch, &
                       landpar
    
    !Argument declarations
    INTEGER,INTENT(IN) :: nsl  !<first dimension of soil layer variables (maxsoillayers)
    INTEGER,INTENT(IN) :: nc   !<second dimension of soil layer variables (nclass)
    REAL,INTENT(INOUT) :: soillayerdepth(nsl,nc)   !<Lower border of soil layers (m) (layer,class)
    REAL,INTENT(INOUT) :: soillayerthick(nsl,nc)   !<Thickness of soil layers (m) (layer,class)
    
    !Local variables
    INTEGER j
  
    !>\b Algoritm \n
    IF(.NOT.soildepthstretch) RETURN
    
    !>Calculate new soil layer thickness (m)
    DO j = 1,nc
      IF(landpar(m_soilstretch,classdata(j)%luse)>0.)THEN
        soillayerthick(2,j) = soillayerthick(2,j) * landpar(m_soilstretch,classdata(j)%luse)
        soillayerthick(3,j) = soillayerthick(3,j) * landpar(m_soilstretch,classdata(j)%luse)
        soillayerdepth(2,j) = soillayerdepth(1,j) + soillayerthick(2,j)
        soillayerdepth(3,j) = soillayerdepth(2,j) + soillayerthick(3,j)
      ENDIF
    ENDDO
  
  END SUBROUTINE reinitiate_soil_depth

  !>\brief Initiate soil water help parameters.
  !!
  !>\b Consequences Module hypevariables epotdist, soilrc, basinrrcscorr, basincevpam, 
  !> basincevpph, basinlp, pwmm, wpmm, fcmm and epmm is set.
  !!
  !>\b Reference ModelDescription Chapter Land routines (Basic assumptions, Soil water -
  !!Groundwater runoff) and Processes above ground (Evaporation)
  !---------------------------------------------------------------------------
  SUBROUTINE initiate_soil_water()

    USE HYPEVARIABLES, ONLY : epotdist, &         !OUT
                              soilrc,   &         !OUT
                              wpmm,fcmm,epmm, &   !OUT
                              pwmm,  & !OUT
                              basinrrcscorr,   & !OUT
                              basincevpam,   & !OUT
                              basincevpph,   & !OUT
                              basinlp,   & !OUT
                              m_epotdist,    &
                              m_cevpam,m_cevpph,m_lp,  &
                              m_wcfc,m_wcwp,m_wcep,  &
                              m_wcfc1,m_wcfc2,m_wcfc3, &
                              m_wcwp1,m_wcwp2,m_wcwp3, &
                              m_wcep1,m_wcep2,m_wcep3, &
                              m_rrcs1,m_rrcs2,m_rrcs3, &
                              m_rrcscorr,  &
                              n_rrcsc,n_rrcs3,n_cevpa,n_cevpp,n_lp
    USE MODVAR, ONLY : classdata,         &
                       basin,             &
                       nsub, nclass,      &
                       maxsoillayers,     &
                       soildepth,         &
                       soilthick,         &
                       regiondivision, &
                       genpar,soilpar,regpar, &
                       conductregest
    USE HYPE_INDATA, ONLY : set_regest_parameter
         
    !Local variables
    INTEGER i,j,isb           !loop-variables (subbasin,class)
    REAL    coeff
    REAL    sums
    REAL    rc0,rc1,rc2,b   !help variables for recession coefficient calculation

    !>\b Algoritm \n
    !>Calculate distribution of potential evaporation between soil layers
    IF(.NOT.ALLOCATED(epotdist)) ALLOCATE(epotdist(2,nclass))
    coeff=genpar(m_epotdist)
    DO j = 1, nclass
      sums = soilthick(1,j)*EXP(-coeff*soildepth(1,j)/2.) + soilthick(2,j)*EXP(-coeff*(soildepth(1,j)+(soildepth(2,j)-soildepth(1,j))/2.))
      epotdist(1,j) = soilthick(1,j)*EXP(-coeff*soildepth(1,j)/2.) / sums
    ENDDO
    epotdist(2,:) = 1 - epotdist(1,:)

    !>Initiate soil water content parameters
    IF(.NOT.ALLOCATED(wpmm)) ALLOCATE(wpmm(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(fcmm)) ALLOCATE(fcmm(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(epmm)) ALLOCATE(epmm(maxsoillayers,nclass))
    IF(.NOT.ALLOCATED(pwmm)) ALLOCATE(pwmm(maxsoillayers,nclass))
    DO j = 1,nclass
      IF(soilpar(m_wcfc1,classdata(j)%soil) > 0)THEN                               !Field capacity (mm)
        fcmm(1,j)=soilpar(m_wcfc1,classdata(j)%soil) * soilthick(1,j) * 1000.       !First layer
        fcmm(2,j)=soilpar(m_wcfc2,classdata(j)%soil) * soilthick(2,j) * 1000.       !Second layer
        fcmm(3,j)=soilpar(m_wcfc3,classdata(j)%soil) * soilthick(3,j) * 1000.       !Third layer
      ELSE
        fcmm(:,j)=soilpar(m_wcfc,classdata(j)%soil) * soilthick(:,j) * 1000.        !All layers
      ENDIF
      IF(soilpar(m_wcwp1,classdata(j)%soil) > 0)THEN                               !Wilting point (mm)
        wpmm(1,j)=soilpar(m_wcwp1,classdata(j)%soil) * soilthick(1,j) * 1000.       !First layer
        wpmm(2,j)=soilpar(m_wcwp2,classdata(j)%soil) * soilthick(2,j) * 1000.       !Second layer
        wpmm(3,j)=soilpar(m_wcwp3,classdata(j)%soil) * soilthick(3,j) * 1000.       !Third layer
      ELSE
        wpmm(:,j)=soilpar(m_wcwp,classdata(j)%soil) * soilthick(:,j) * 1000.        !All layers
      ENDIF
      IF(soilpar(m_wcep1,classdata(j)%soil) > 0)THEN                               !Effectiv porosity (mm)
        epmm(1,j)=soilpar(m_wcep1,classdata(j)%soil) * soilthick(1,j) * 1000.       !First layer
        epmm(2,j)=soilpar(m_wcep2,classdata(j)%soil) * soilthick(2,j) * 1000.       !Second layer
        epmm(3,j)=soilpar(m_wcep3,classdata(j)%soil) * soilthick(3,j) * 1000.       !Third layer
      ELSE
        epmm(:,j)=soilpar(m_wcep,classdata(j)%soil) * soilthick(:,j) * 1000.        !All layers
      ENDIF
    ENDDO
    pwmm = wpmm + fcmm + epmm

    !Set soil runoff recession correction
    IF(.NOT.ALLOCATED(basinrrcscorr)) ALLOCATE(basinrrcscorr(nsub))
    DO i = 1,nsub
      IF(basin(i)%parregion(regiondivision(m_rrcscorr))>0)THEN
        basinrrcscorr(i) = 1. + regpar(m_rrcscorr,basin(i)%parregion(regiondivision(m_rrcscorr)))   !Correction of recession coefficients
      ELSE
        basinrrcscorr(i)  = 1.
      ENDIF
        
      !Replace parameter values with regional parameter estimates
      IF(conductregest) CALL set_regest_parameter(i,n_rrcsc,basinrrcscorr(i),1.)
    ENDDO

    !Initiate soil runoff recession coeffcients
    IF(.NOT.ALLOCATED(soilrc)) ALLOCATE(soilrc(maxsoillayers,nclass,nsub))
    !>Calculate adjustment factors
    rc2 = genpar(m_rrcs3)                         !Runoff coefficient slope dependence   
    DO i = 1,nsub
      !>Replace parameter values with regional parameter estimates
      IF(conductregest) CALL set_regest_parameter(i,n_rrcs3,rc2)

      !>Calculate soil runoff recession coeffcients for each soil layer, class and subbasin
      DO j=1,nclass
        rc0 = soilpar(m_rrcs1,classdata(j)%soil)*basinrrcscorr(i)       !Runoff coefficient in surface layer 
        IF(rc0>1.) rc0 = 1.
        rc0 = rc0+rc2*basin(i)%slope    !runoff coefficient in upper soil layer (slope dependent)
        IF(rc0>1.) rc0 = 1.
        rc1 = soilpar(m_rrcs2,classdata(j)%soil)*basinrrcscorr(i)       !Runoff coefficient in bottom layer 
        IF(rc1>1.) rc1 = 1.
        IF(rc1==0) rc1 = rc0
        b = LOG(rc0/rc1)/((soildepth(3,j)-soilthick(3,j)/ 2.) - soilthick(1,j)/ 2.)
        soilrc(1,j,i) = rc0
        soilrc(3,j,i) = rc1
        soilrc(2,j,i) = rc0 * EXP(-b* (soildepth(2,j) - soilthick(2,j)/2. - soilthick(1,j)/ 2.))
      ENDDO
    ENDDO
    
    !>Set evaporation seasonal corrections
    IF(.NOT.ALLOCATED(basincevpam))THEN
      ALLOCATE(basincevpam(nsub))
      ALLOCATE(basincevpph(nsub))
    ENDIF
    basincevpam = genpar(m_cevpam)
    basincevpph = genpar(m_cevpph)
    !Replace parameter values with regional parameter estimates
    IF(conductregest)THEN
      DO isb = 1,nsub
        CALL set_regest_parameter(isb,n_cevpa,basincevpam(isb))
        CALL set_regest_parameter(isb,n_cevpp,basincevpph(isb))
      ENDDO
    ENDIF

    !>Set evaporation subbasin parameter
    IF(.NOT.ALLOCATED(basinlp)) ALLOCATE(basinlp(nsub))
    basinlp = genpar(m_lp)
    IF(conductregest)THEN   !Replace parameter value with regional parameter estimates
      DO isb = 1,nsub
        CALL set_regest_parameter(isb,n_lp,basinlp(isb))
      ENDDO
    ENDIF      

  END SUBROUTINE initiate_soil_water

  !>Subroutine for calculation of changes in snow pack; snowfall addition, snow 
  !>pack melting and snow age
  !!
  !>\b Reference ModelDescription Chapter Land routines (Snow routines)
  !------------------------------------------------------------------------
  SUBROUTINE calculate_snow(i,j,subid,iluse,snowfall,csnowfall,snow,csnow,temp,melt, &
                            cmelt,swrad,snowage,snowmax,snowdepth,snowcover,epot,evap, &
                            cevap,effcov)
  
    USE MODVAR, ONLY : genpar,          &
                       numsubstances,   &
                       missing_value,   &
                       modeloption,     &
                       p_snowmelt,      &
                       i_t1,i_t2,       &
                       p_snowevap
    USE HYPEVARIABLES, ONLY : m_cmrad,m_fsceff,m_cmrefr, &
                              m_dsndens,m_T1evap
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: j          !<index of current class
    INTEGER, INTENT(IN) :: subid      !<subbasin id
    INTEGER, INTENT(IN) :: iluse      !<index of landuse
    REAL, INTENT(IN)    :: snowfall   !<precipitation as snow (mm/timestep)
    REAL, INTENT(IN)    :: csnowfall(numsubstances) !<concentration of precipitation as snow 
    REAL, INTENT(INOUT) :: snow       !<snow pack (mm)
    REAL, INTENT(INOUT) :: csnow(numsubstances) !<concentration of snow 
    REAL, INTENT(IN)    :: temp       !<air temperature (C)
    REAL, INTENT(OUT)   :: melt       !<snow melt (mm/timestep)
    REAL, INTENT(OUT)   :: cmelt(numsubstances)     !<substances of snow melt
    REAL, INTENT(IN)    :: swrad      !<shortwave radiation (MJ/m2/day?)
    REAL, INTENT(INOUT) :: snowage    !<age of snow (timesteps)
    REAL, INTENT(INOUT) :: snowmax    !<maximum snow pack during winter (mm)
    REAL, INTENT(INOUT) :: snowdepth  !<current depth of snow (cm)
    REAL, INTENT(INOUT) :: snowcover  !<snowcover fraction
    REAL, INTENT(IN)    :: epot       !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap       !<snow sublimation (mm/timestep)
    REAL, INTENT(OUT)   :: effcov     !<effective snowcover used for scaling snow and soil evaporation (0 if snowevap is switched off)
    REAL, INTENT(OUT)   :: cevap(numsubstances)   !<concentrations in snow sublimation
    
    !Local variables
    INTEGER status
    REAL fsceff,abla,abla0 
    REAL oldsnow    !snow pack at beginning of calculations (cm)

    !>\b Algorithm \n
    !>Set parameter values and default output
    fsceff = genpar(m_fsceff)   !efficiency of fractional snow cover to reduce melt and evap
    effcov = 1.-fsceff*(1.-snowcover) !effective snow cover used to scale melt and evap
    oldsnow = snow
    cmelt = 0.
    cevap = 0.
    
    !>Calculate potential snow melt
    CALL calculate_snowmelt(iluse,temp,swrad,snow,snowage,effcov,melt)
    IF(melt>snow) melt = snow
    
    !>Calculate potential evaporation from snow (sublimation)
    IF(modeloption(p_snowevap)>1)THEN
      IF(snow>0)THEN
        evap = epot * effcov
      ELSE
        evap = 0.
        effcov = 0. !make sure effcov = 0 if there is no snow, otherwise there will be no soil evaporation
      ENDIF
    ELSE
      evap = 0.
      effcov = 0. !make sure effcov = 0 if snowevap is switched off, otherwise there will be no soil evaporation
    ENDIF
    IF(evap>snow) evap = snow

    !>Calculate ablation and check against snow pack
    abla0 = melt + evap      !potential ablation (melt+evap)
    abla = MIN(abla0, snow)  !limit ablation to available snow
    IF(abla0>0.)THEN
      IF(abla<abla0)THEN
        melt = melt * abla/abla0   !distribute ablation on melt and evap
        evap = evap * abla/abla0
        IF(melt+evap>snow) evap = snow - melt !extra safe?
      ENDIF
    ELSE
      melt = 0.
      evap = 0.
    ENDIF
    
    !>Update snow pack with snowfall, melting and sublimation
    IF(snowfall>0.) CALL add_water(numsubstances,snow,csnow,snowfall,csnowfall)
    IF(melt>0.)THEN
      cmelt = csnow
      IF(melt<snow)THEN
        CALL remove_water(snow,numsubstances,csnow,melt,cmelt,status)
        IF(status/=0) CALL error_remove_water(errstring(20),subid,i,j)
      ELSE
        snow = 0.
        csnow = 0.  !remove last traces
      ENDIF
      IF(i_t2>0) cmelt(i_t2)=0. !temp.conc. in meltwater = 0
    ENDIF
    IF(evap>0.)THEN
      cevap=0.
      IF(i_t1>0) cevap(i_t1)= genpar(m_T1evap) * csnow(i_t1)
      IF(evap<snow)THEN
        CALL remove_water(snow,numsubstances,csnow,evap,cevap,status)
        IF(status/=0) CALL error_remove_water(errstring(21),subid,i,j)
      ELSE
        snow = 0.
        csnow = 0.  !remove last traces
      ENDIF
    ENDIF  

    !>Calculate degree of snowcover and update maximum snow pach during winter
    CALL calculate_fractional_snowcover(iluse,0.,snow,snowmax,snowcover)

    !>Calculate snow age and snow depth
    CALL calculate_snowdepth(iluse,snow,oldsnow,snowfall,temp,genpar(m_dsndens),snowage,snowdepth)

  END SUBROUTINE calculate_snow
  
  !>Subroutine for calculation of snow melt by different methods
  !!
  !> \b Reference ModelDescription Chapter Land routines (Snow routines)
  !------------------------------------------------------------------------
  SUBROUTINE calculate_snowmelt(iluse,temp,swrad,snow,snowage,effcov,melt)
  
    USE MODVAR, ONLY : landpar,   &
                       genpar,          &
                       modeloption,     &
                       p_snowmelt
    USE HYPEVARIABLES, ONLY : m_ttmp,m_cmlt,  &
                              m_snalbmin,m_snalbmax,m_snalbkexp,  &
                              m_cmrad,m_cmrefr
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iluse      !<index of landuse
    REAL, INTENT(IN)    :: temp       !<air temperature (C)
    REAL, INTENT(IN)    :: swrad      !<shortwave radiation (MJ/m2/day?)
    REAL, INTENT(IN)    :: snow       !<snow pack (mm)
    REAL, INTENT(IN)    :: snowage    !<age of snow (timesteps)
    REAL, INTENT(IN)    :: effcov     !<effective snowcover used for scaling of snow melt
    REAL, INTENT(OUT)   :: melt       !<snow melt (mm/timestep)
    
    !Local variables
    REAL tt       !threshold temperature for snow melt (and evaporation) (C)
    REAL cm       !coefficient for snow melt (mm/C/timestep)
    REAL snowalbedo
    REAL snalbmax, snalbmin, snalbkexp
    REAL cmrad     !radiation index snow melt factor (mm/MJ/timestep)
    REAL cmrefr  

    !>\b Algorithm \n
    !>Set parameter values and default output
    tt = landpar(m_ttmp,iluse)  !threshold temperature for snow melt
    cm = landpar(m_cmlt,iluse)  !Coefficient for snow melt
    melt = 0.
    
    !Melting calculation sequence:
    ! -> To reduce the number of snowmelt options, snow cover melt scaling 
    !    is now included in all snowmelt models. Instead of options for each combination
    !    of melt and sublimation, the fraction of snow cover reduction is controlled
    !    by a new parameter ffscred (fraction of fsc reduction)
    ! -> For backward compitability, the previous snowmelt options values are still used 
    !    (0,1 temp index, 2 temp+rad index) - but note that p_snowmelt = 0 may now also 
    !    imply snowcover scaling and sublimation, depending on parameter ffscred and fepotsnow.
    ! -> Sublimation is calculated separately after the melt subroutine, and is also controlled by ffscred.
    ! -> Ablation = melt + sublimation is introduced
    ! -> Minimization of ablation to current snow is made after calculation of (potential) 
    !    melt and sublimation. The reduction from potential to actual ablation is finally
    !    distributed on melt and sublimation.

    !>Select snow melt model
    SELECT CASE(modeloption(p_snowmelt))
    CASE(0,1) 
      !>\li Case 0 and 1: Original temperature index model, calculated with or without snowcover scaling
      IF(snow>0 .AND. temp >= tt) THEN
        melt = cm   * (temp - tt)  !potential melt
        melt = melt * effcov       !snowcover melt scaling (no reduction of snowcover=1 and/or fsceff=0)
      ENDIF
    CASE(2)
      !>\li Temperature and radiation index model, with/without snowcover scaling and refreezing
      !Set parameter values
      snalbmin  = landpar(m_snalbmin,iluse)
      snalbmax  = landpar(m_snalbmax,iluse)
      snalbkexp = landpar(m_snalbkexp,iluse)
      cmrad     = landpar(m_cmrad,iluse)
      cmrefr    = genpar(m_cmrefr)
    
      !Radiation melt component
      snowalbedo = snowalbedo_function(snowage,snalbmin,snalbmax,snalbkexp)
      melt = cmrad * swrad * (1.-snowalbedo)
      
      !Add Temperature component
      IF(snow>0. .AND. temp >= tt)THEN
        melt = melt + cm * (temp - tt)
      ENDIF
      
      !Refreezing component when temperatures below tt, as a fraction cmrefr of cm 
      IF(snow>0. .AND. temp < tt .AND. melt > 0.)THEN
        melt = melt - cmrefr * cm * (tt - temp)
        IF(melt<0.) melt = 0.
      ENDIF
      
      !Scale melt with fractional snow cover
      melt = melt * effcov
    CASE DEFAULT
      !>\li Case default (0): Original temperature index model, calculated with or without snowcover scaling
      IF(snow>0 .AND. temp >= tt) THEN
        melt = cm   * (temp - tt)  !potential melt
        melt = melt * effcov       !snowcover melt scaling
      ENDIF
    END SELECT
    
  END SUBROUTINE calculate_snowmelt
  
!>Function to calculate snow albedo depending on the snow age        
!>
!> \b Reference ModelDescription Chapter Land routines (Snow routines)
!------------------------------------------------------------------------
  FUNCTION snowalbedo_function(snowage,albmin,albmax,kexp) RESULT(albedo)
  
  !Argument declarations
  REAL,INTENT(IN)  :: snowage !<snow age (timesteps)
  REAL,INTENT(IN)  :: albmin  !<minimum albedo (typical value 0.4)
  REAL,INTENT(IN)  :: albmax  !<maximum albedo (typical value 0.9)
  REAL,INTENT(IN)  :: kexp    !<exponential factor (1/time step) (typical value 0.1 for daily timesteps)
  REAL             :: albedo  ! albedo, fractional reflection of shortwave radiation (-)
  !David Gustafsson, 2013-02-05
  
  !Calculate albedo with a simplified exponential function
  albedo = albmin+(albmax-albmin)*EXP(-kexp*snowage)
  
  END FUNCTION snowalbedo_function

  !>\brief Subroutine for calculation of fractional snow cover area
  !>The maximum snow of the winter season is also calculated
  !>
  !Based on Lindström&Gardelin(1999;2000) following the implementation in 
  !the Rossby centre RCA-model (Samuelsson et al, 2006)
  !
  !> \b Reference ModelDescription Chapter Land routines (Snow routines - Snow cover)
  !------------------------------------------------------------------------
  SUBROUTINE calculate_fractional_snowcover(iluse,elevstd,snow,snowmax,fsc)
  
    USE MODVAR, ONLY : landpar, genpar, seconds_per_timestep
    USE HYPEVARIABLES, ONLY :  m_fscmax,m_fscmin,m_fsclim, & 
                               m_fscdistmax, m_fscdist0,m_fscdist1, &
                               m_fsck1, m_fsckexp
    !Argument declarations
    INTEGER, INTENT(IN) :: iluse       !<index of landuse
    REAL, INTENT(IN)    :: elevstd     !<standard deviation of elevation (m)
    REAL, INTENT(IN)    :: snow        !<snow pack (mm)
    REAL, INTENT(INOUT) :: snowmax     !<maximum snow pack during winter (mm)
    REAL, INTENT(OUT)   :: fsc         !<fractional snowcover area (-)

    !Local variables
    REAL timestep_seconds,fscdist
    timestep_seconds = REAL(seconds_per_timestep)
    
    !>\b Algorithm \n
    !>Check snow pack status
    IF(snow.gt.0.)THEN   !Snow present
      !>Check if snowcover model is used
      IF(genpar(m_fscmax)==0)THEN
        fsc = 1.
      ELSE
        !>Check snowpack development phase, and select corresponding FSC function
        IF(snowmax.le.0.)THEN
          !1) Accumulation phase, snowmax = 0
          !1.1) fsc = tangens-hyperbolic function, Eq 28 (Samuelsson 2006)
          fsc = MAX(genpar(m_fscmin),genpar(m_fscmax) * TANH(0.1 * snow))
          !1.2) Set snowmax = snow, if fsc >= fscmax - fsclim
          IF(fsc.ge.(genpar(m_fscmax)-genpar(m_fsclim)))THEN
            snowmax = snow
          ENDIF
        ELSE
          !2) Melting phase, snowmax>0 (onset in previous timesteps)
          !2.1) update snowmax
          IF(snow.GT.snowmax)THEN
            !update snowmax to new maximum snow value
            snowmax = snow
          ELSE
            !decrease snowmax towards end of melt season, eq. 31 (Samuelsson 2006)
            IF(snow.LT.genpar(m_fsck1)*snowmax)THEN
              snowmax = snowmax - (genpar(m_fsck1) * snowmax - snow)*(1.-EXP(-genpar(m_fsckexp) * timestep_seconds)) / genpar(m_fsck1)
            ENDIF 
          ENDIF
          !2.2) calculate snow distribution factor, Eq 30 (Samuelsson 2006)
          fscdist = MIN(landpar(m_fscdistmax,iluse),landpar(m_fscdist0,iluse) + landpar(m_fscdist1,iluse) * elevstd)
          !2.3) fsc=linear function, Eq 29 (Samuelsson 2006)
          fsc = MAX(genpar(m_fscmin),MIN(genpar(m_fscmax),snow / (snowmax * fscdist)))
        ENDIF
      ENDIF 
    ELSE  !No snow
      snowmax = 0.
      fsc = 0.
    ENDIF
    
  END SUBROUTINE calculate_fractional_snowcover
  
  !>Calculation of snowdepth and age of snow. Snow depth depends on age of snow.
  !>
  !>\b Reference ModelDescription Chapter Land routines (Snow routines - Soil temperature and snow depth)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_snowdepth(iluse,snow,oldsnow,snowfall,temp,snowdensdt,snowage,snowdepth)

    USE MODVAR, ONLY : timesteps_per_day, &
                       modeloption,   &
                       p_snowdensity,  &
                       genpar,landpar
    USE HYPEVARIABLES, ONLY : m_ttmp, &
                              m_sndens0,m_sdnsmax, &
                              m_sdnsrate,m_sdnsradd
    
    !Argument declarations
    INTEGER, INTENT(IN) :: iluse      !<index of landuse
    REAL, INTENT(IN)    :: snow       !<snow water equivalent (mm)
    REAL, INTENT(IN)    :: oldsnow    !<snow water equivalent before snowfall/melt this timestep (mm)
    REAL, INTENT(IN)    :: snowfall   !<precipitation as snow (mm/timestep)
    REAL, INTENT(IN)    :: temp       !<air temperature (C)
    REAL, INTENT(IN)    :: snowdensdt !<snow density increase due to ageing (g/cm3.timestep)
    REAL, INTENT(INOUT) :: snowage    !<help variable for snow; age of snow pack (timesteps)
    REAL, INTENT(INOUT) :: snowdepth  !<current depth of snow (cm)

    !Local variables
    REAL snowdens    !snow density (g/cm3)
    REAL ddens       !change in snow density (g/cm3/timestep)
    REAL snowdens0   !snow density at snowfall, model parameter
    REAL tt          !threshold temperature for snow melt (and evaporation) (C)

    !>\b Algorithm \n
    !>Set parameter values

    !> \b Algorithm \n
    !>Set model parameters
    tt = landpar(m_ttmp,iluse)
    snowdens0 = genpar(m_sndens0)
    
    !>Update snow age with time step and new snow
    IF(snow==0)THEN
      snowage = 0.
    ELSE
      snowage = snowage + 1.
      IF(oldsnow < snow) snowage = snowage * oldsnow / snow     !Assume that melt is drawn from snowfall in first hand
    ENDIF
    
    !>Calculate snow depth, depends of choice of snow density model
    IF(snow<=0.)THEN
      snowdepth  = 0.
    ELSE  
      SELECT CASE(modeloption(p_snowdensity))
      !>\li Case of age depending snow density model (0,default):
      CASE(0)
        snowdens = snowdens0 + snowdensdt * snowage / REAL(timesteps_per_day)
        IF(snowdens>0.) snowdepth  = 0.1 * snow / snowdens     !0.1 cm/mm
      !>\li Case of compacting factor snow density model (1):
      CASE(1)
        snowdepth = snowdepth + snowfall/10./snowdens0  !add snowfall to old snow
        snowdens = (oldsnow+snowfall)/10./snowdepth
        IF(temp>tt)THEN
          ddens = (genpar(m_sdnsrate)+genpar(m_sdnsradd)) * (genpar(m_sdnsmax) - snowdens)  !snow compactation, warm days
        ELSE
          ddens = genpar(m_sdnsrate) * (genpar(m_sdnsmax) - snowdens)  !snow compactation, cold days
        ENDIF
        snowdens = snowdens + ddens
        snowdepth = snow/10./snowdens
      ENDSELECT
    ENDIF

  END SUBROUTINE calculate_snowdepth

  !>Calculate air pressure (kPa) as a function of elevation, FAO(7)
  !-------------------------------------------------------------------------------
  REAL FUNCTION airpressure_elevationfunction(elev)
     
    !Argument decalaration
    REAL, INTENT(IN) :: elev   !<elevation
     
    airpressure_elevationfunction = 101.3 * ((293. - 0.0065 * elev)/293. ) ** 5.26
     
  END FUNCTION airpressure_elevationfunction
  
  !>Calculate latent heat of vaporization (MJ kg-1) as a function of temperature
  !-------------------------------------------------------------------------------
  REAL FUNCTION latentheat_tempfunction(temp)

    !Argument decalaration
    REAL, INTENT(IN) :: temp !<temperature (C)
     
    latentheat_tempfunction = 2.501 - 0.002361 * temp  ![MJ/kg]
     
  END FUNCTION latentheat_tempfunction
  
  !> Calculate psychrometric constant (kPa C^-1) as a function of
  !! air pressure and latent heat of vaporization, FAO
  !-------------------------------------------------------------------------------
  REAL FUNCTION psychrometric_constant(airpressure,lambda)
  
    !Argument declarations
    REAL, INTENT(IN) :: airpressure !<air pressure [kPa]
    REAL, INTENT(IN) :: lambda      !<latent heat of vaporization [MJ/kg]
    
    !Parameter declaration
    REAL, PARAMETER :: cp = 0.001013  !specific heat of moist air at constant pressure (MJ kg^-1 C^-1)
 
    psychrometric_constant = cp * airpressure / (0.622 * lambda)
 
  END FUNCTION psychrometric_constant

  !>Calculates potential evaporation or uses a value supplied as input
  !TODO: Why is this in soil_proc. Move to atm_proc?
  !
  !> \b Reference ModelDescription Processes above ground (Evaporation)
  !--------------------------------------------------------------
  SUBROUTINE calculate_potential_evaporation(i,j,temp,epot,radext,swrad,netrad,actvap,satvap,wind,epotsnow)
  
    USE MODVAR, ONLY : basin,classdata, &
                       landpar, &
                       genpar, &
                       dayno, &
                       pi, &
                       xobsi,xobsindex, &
                       classbasin, &
                       tsofday,timesteps_per_day
    USE HYPEVARIABLES, ONLY : o_reepot, &
                              m_ttmp,m_cevp, &
                              basincevpam, &
                              basincevpph, &
                              m_kc, &
                              m_krs,m_jhtadd,m_jhtscale,m_alfapt,m_fepotsnow

    !Argument declarations
    INTEGER, INTENT(IN) :: i      !<index of current subbasin
    INTEGER, INTENT(IN) :: j      !<index of current class 
    REAL, INTENT(IN)    :: temp   !<air temperature
    REAL, INTENT(OUT)   :: epot   !<potential evapotranspiration [mm/timestep]
    REAL, INTENT(IN)    :: radext !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: swrad  !<downward shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: netrad !<net downward radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: actvap !<actual vapor pressure [kPa]
    REAL, INTENT(IN)    :: satvap !<saturated vapour pressure [kPa]
    REAL, INTENT(IN)    :: wind   !<wind speed [m/s]
    REAL, INTENT(OUT)   :: epotsnow !<potential evapotranspiration for snow mm/timestep
    
    !Local variables
    REAL tt       !threshold temperature for melting (C)
    REAL ce       !coefficient for evaporation (mm/C/timestep)
    REAL dsatvap  !Slope of saturation pressure curve [kPa/C]
    REAL gamma    !psychrometric constant
    REAL lambda   !latent heat of evaporation [MJ/kg]
    REAL airpressure  !atmospheric pressure [kPa]
    REAL kc       !crop coefficient used for the new PET functions
    REAL elev     !elevation
    REAL turbidity ! atmospheric turbidity
    REAL fepotsnow ! fraction of potential evaporation used for snow
    INTEGER current_petmodel
    
    !>\b Algorithm \n
    !>Set local parameters and corrections
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation
    ce = landpar(m_cevp,classdata(j)%luse)       !Coefficient for potential evaporation
    fepotsnow = landpar(m_fepotsnow,classdata(j)%luse)       !Coefficient for potential evaporation for snow
!    ce = ce * (1 + basincevpam(i)*SIN(2.*pi*(dayno-basincevpph(i))/365.))
    ce = ce * (1 + basincevpam(i)*SIN(2.*pi*(dayno-1+REAL(tsofday)/REAL(timesteps_per_day)-basincevpph(i))/365.))
    
    !>Calculate additional input variables for the alternative PET functions
    current_petmodel = get_current_petmodel(i)
    IF(current_petmodel.GT.1)THEN
      dsatvap = deltasaturationpressure_function(temp)  !Slope of saturated vapour pressure curve, using mean temperature
      lambda = latentheat_tempfunction(temp)  !Latent heat of vaporization
      elev = basin(i)%elev+classbasin(i,j)%deltah
      airpressure = airpressure_elevationfunction(elev)  !Air pressure, assuming normal pressure at sea level
      gamma = psychrometric_constant(airpressure,lambda) !Psychrometric constant
      kc = landpar(m_kc(current_petmodel),classdata(j)%luse)  !PET-model specific Landuse scaling parameter, "crop coefficient"
      IF(kc==0) kc = landpar(m_kc(1),classdata(j)%luse)  !Default Landuse scaling parameter, "crop coefficient"
      turbidity = swrad / radext
    ENDIF
      
    !>Calculate potential evaporation with the selected petmodel
    epot = 0.
    SELECT CASE(current_petmodel)
      CASE(0) !HYPE original model (with Xobs replacement, if available)
        IF(xobsindex(o_reepot,i)>0)THEN       
          epot = xobsi(xobsindex(o_reepot,i))
        ELSEIF(temp>tt)THEN
          epot = ce*(temp-tt)
        ELSE
          epot = 0.  
        ENDIF
      CASE(1) !HYPE original model (without Xobs replacement)
        IF(temp>tt)THEN
          epot = ce*(temp-tt)
        ELSE
          epot = 0.  
        ENDIF
      CASE(2) !Modified Jensen-Haise/McGuinness following Oudin et al (2005)
        !parameters suggested by Oudin et al, jhtadd = 5, jhtscale = 100
        epot = kc * MAX(0.,radext / (lambda) * (temp + genpar(m_jhtadd)) / genpar(m_jhtscale))
      CASE(3) !Hargreaves-Samani (known to overpredict in humid areas)
        ! The function is modified by DG to limit the "turbidity-factor" with the Ångström formula:
        ! 
        !   The original Hargreaves function is:
        !     epot = 0.0023 * radext / (lambda*rho) * (Tmax-Tmin)^0.5 * (temp + 17.8)
        !   and the Hargreaves turbidity for estimating swrad = krs * (Tmax-Tmin)^0.5
        !
        !   Thus, by replacing (Tmax-Tmin)^2 with turbidity/krs, we get a reasonable limitation of the Hargreaves (tmax-Tmin) impact
        !  (furthermore, if Tmax-min was missing, we actually use the clearsky turbidity at this point)
        !
        ! also note that rho = 1 and excluded in equations below...
        epot = max(0.,kc * 0.0023 * radext /(lambda) * turbidity / genpar(m_krs) * (temp+17.8))
      CASE(4) ! Priestly Taylor (known to underpredict in arid and semi-arid areas)
        epot = max(0.,kc * genpar(m_alfapt) * dsatvap * netrad / (lambda * (dsatvap+gamma)))
      CASE(5) ! FAO Penman Monteith reference crop evapotranspiration
        epot = max(0., kc * ((0.408 * dsatvap * netrad + gamma*900./(temp+273.)*wind*(satvap-actvap))/(dsatvap+gamma*(1.+0.34*wind))))
      END SELECT
      !>Calculate potential evaporation for snow evaporation (sublimation)
      epotsnow = fepotsnow * epot

  END SUBROUTINE calculate_potential_evaporation

  !>\brief Get the pet model for the current subbasin
  !--------------------------------------------------------------------
  INTEGER FUNCTION get_current_petmodel(i)

    USE MODVAR, ONLY : conductbasinpetmodel, &
                       petmodel, &
                       modeloption, &
                       p_petmodel
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i  !<current subbasin
    
    IF(conductbasinpetmodel)THEN
      get_current_petmodel = petmodel(i)
    ELSE
      get_current_petmodel = modeloption(p_petmodel)
    ENDIF
    
  END FUNCTION get_current_petmodel
  
  !>\brief Calculate and remove evapotranspiration from the soil upper
  !>two layers
  !
  !> \b Reference ModelDescription Chapter Processes above ground (Evaporation)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_actual_soil_evapotranspiration(i,j,temp,epot,wp,fc,  &
                                      epotfrac,soilstate,evap,evapflows,cevap,barefrac)

    USE MODVAR, ONLY : basin,classdata, &
                       genpar,   &
                       landpar,  &
                       numsubstances,   &
                       maxsoillayers,   &
                       i_t1,            &
                       soilthick, &
                       realzero
    USE HYPEVARIABLES, ONLY : basinlp, &
                              m_ttmp,  &
                              m_T1evap, &
                              m_ttrig,m_tredA,m_tredB

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<index of current subbasin
    INTEGER, INTENT(IN) :: j                        !<index of current class
    REAL, INTENT(IN)    :: temp                     !<air temperature
    REAL, INTENT(IN)    :: epot                     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<field capacity (mm)
    REAL, INTENT(IN)    :: epotfrac(2)              !<relative distribution of potential evaporation between upper two soil layers (-)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    REAL, INTENT(OUT)   :: evap                     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evapflows(2)             !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubstances)     !<concentration in evapotranspiration (?)
    REAL, INTENT(IN)    :: barefrac                 !<fraction of soil that has evapotranspiration (-)

    !Local variables
    INTEGER k   !loop-variable
    INTEGER status  !error status of subroutine
    REAL evap1,evap2  !evapotranspiration of soillayer 1 and 2 (mm/timestep)
    REAL cevap1(numsubstances),cevap2(numsubstances)  !concentration of evapotranspiration
    REAL soiltemp_reduction(maxsoillayers),soiltemp

    !>\b Algorithm \n
    !Default values output variables
    evap = 0.
    cevap = 0.
    evapflows = 0.

    !>Calculate soil temperature reduction
    soiltemp_reduction=1.
    IF(landpar(m_tredA,classdata(j)%luse).GT.0.)THEN
      DO k=1,2
        IF(soilthick(k,j)>0)THEN
          soiltemp = soilstate%temp(k,j,i)
          IF(soiltemp.GT.landpar(m_ttrig,classdata(j)%luse))THEN
            soiltemp_reduction(k) = 1.-EXP(-landpar(m_tredA,classdata(j)%luse) * (soiltemp-landpar(m_ttrig,classdata(j)%luse))**landpar(m_tredB,classdata(j)%luse))
          ELSE
            soiltemp_reduction(k) = 0.
          ENDIF
        ENDIF
      ENDDO
    ENDIF

    !>If temperature above threshold:
    IF(temp>landpar(m_ttmp,classdata(j)%luse))THEN

      !>\li calculate actual evapotranspiration in the uppermost layer 
      IF(soilstate%water(1,j,i) - wp(1)> basinlp(i) * fc(1)) THEN
        evap1 = epot*epotfrac(1)*soiltemp_reduction(1)
      ELSEIF(soilstate%water(1,j,i)-wp(1) <= 0.0) THEN
        evap1 = 0.0
      ELSE
        evap1 = epot*epotfrac(1)*((soilstate%water(1,j,i)-wp(1))/(basinlp(i) * fc(1)))*soiltemp_reduction(1)
      ENDIF
      IF(evap1+realzero>soilstate%water(1,j,i)-wp(1)) evap1 = soilstate%water(1,j,i)-wp(1)
      DO k=1,numsubstances
        cevap1(k) = 0.
        IF(k==i_t1) cevap1(k) = genpar(m_T1evap) * soilstate%conc(k,1,j,i)
      ENDDO
      evap1 = evap1*MIN(1.,barefrac)  !Scale evapotranspiration with fraction of bare soil

      !>\li Remove evapotranspiration of soillayer 1
      IF(evap1+realzero<soilstate%water(1,j,i))THEN
        CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),evap1,cevap1,status)
        IF(status.NE.0) CALL error_remove_water(errstring(7),basin(i)%subid,i,j)
      ELSE
        soilstate%water(1,j,i) = 0.
        IF(numsubstances>0) soilstate%conc(:,1,j,i) = 0.    !remove last traces, safe for wp=0
      ENDIF
      
      !Second soillayer:
      IF(soilthick(2,j)>0)THEN

        !>\li Calculate actual evapotranspiration in the second soillayer
        IF(soilstate%water(2,j,i)-wp(2) > basinlp(i) * fc(2)) THEN
          evap2 = epot*epotfrac(2)*soiltemp_reduction(2)
        ELSEIF(soilstate%water(2,j,i)-wp(2) <= 0.0) THEN
          evap2 = 0.0
        ELSE
          evap2 = epot*epotfrac(2)*((soilstate%water(2,j,i)-wp(2))/(basinlp(i) * fc(2)))*soiltemp_reduction(2)
        ENDIF
        IF(evap2+realzero>soilstate%water(2,j,i)-wp(2)) evap2 = soilstate%water(2,j,i) - wp(2)
        DO k=1,numsubstances
          cevap2(k) = 0.
          IF(k==i_t1) cevap2(k) = genpar(m_T1evap) * soilstate%conc(k,2,j,i)
        ENDDO
        evap2 = evap2*MIN(1.,barefrac)  !Scale evapotranspiration with fraction of bare soil

        !>\li Remove evapotranspiration of soillayer 2
        IF(evap2+realzero<soilstate%water(2,j,i))THEN
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),evap2,cevap2,status)
          IF(status.NE.0) CALL error_remove_water(errstring(8),basin(i)%subid,i,j)
        ELSE
          soilstate%water(2,j,i) = 0.
          IF(numsubstances>0) soilstate%conc(:,2,j,i) = 0.  !remove last traces, safe for wp=0
        ENDIF
      ELSE
        evap2 = 0.
        IF(numsubstances>0) cevap2 = 0.
      ENDIF

      !>Set output variables
      evap = evap1 + evap2
      evapflows(1) = evap1
      evapflows(2) = evap2
      IF(i_t1>0 .AND. evap>0) cevap(i_t1) = (cevap1(i_t1)*evap1 + cevap2(i_t1)*evap2)/evap
    ENDIF

  END SUBROUTINE calculate_actual_soil_evapotranspiration

  !>\brief Drainage level runoff: tile or drainage pipe
  !!
  !> \b Reference ModelDescription Chapter Land routines (Soil water - Runoff through drainage pipes)
  !------------------------------------------------------------------
  SUBROUTINE calculate_tile_drainage(i,j,isoil,subid,wp,fc,ep,&
       sdepth,sthick,tdepth,rrcscorr,soilstate,runoffd,crunoffd,cweights)

    USE MODVAR, ONLY : soilpar,  &
                       numsubstances,   &
                       maxsoillayers
    USE HYPEVARIABLES, ONLY : m_trrcs

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of current subbasin
    INTEGER, INTENT(IN) :: j        !<index of current class
    INTEGER, INTENT(IN) :: isoil    !<soil type index
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    REAL, INTENT(IN)    :: wp(maxsoillayers) !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers) !<"field capacity" volume (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers) !<effective porosity volume (mm)
    REAL, INTENT(IN)    :: sdepth(maxsoillayers) !<Lower border of soil layers (m)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<Thickness of soil layers (m)
    REAL, INTENT(IN)    :: tdepth                !<Tiledepth (m)
    REAL, INTENT(IN)    :: rrcscorr    !<correction of recession coefficients
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    REAL, INTENT(OUT)   :: runoffd                    !<runoff
    REAL, INTENT(OUT)   :: crunoffd(numsubstances)    !<concentration of runoff 
    REAL, INTENT(OUT)   :: cweights(maxsoillayers)    !<weights for calc. drain.conc from layer.conc (zero or one)
    
    !Local variables 
    INTEGER status    !error status of subroutine call
    REAL    deltah    !groundwater level above tile drainage pipe level (m)
    REAL    trc       !coefficient for runoff recession tile flow (fraction per timestep)

    !>\b Algorithm \n
    !>Set default output values
    runoffd = 0.
    crunoffd = 0.
    cweights(:) = 0.

    IF(tdepth<=0) RETURN   !no tile drainage

    !>Set local parameters
    trc = soilpar(m_trrcs,isoil)*rrcscorr       !Runoff coefficient for tile runoff
    IF(trc==0)RETURN
    IF(trc>1.) trc = 1.

    !> epending on depth of tile drainage pip calculate:
    IF(tdepth>0 .AND. trc>0.)THEN
      IF(tdepth<=sdepth(1))THEN       !Drainage pipe in uppermost layer

        !>\li drainage in uppermost soil layer
        deltah = (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1) - (sdepth(1) - tdepth)  !m
        IF(deltah>0.)THEN
          runoffd = trc * deltah / sthick(1) * ep(1)
          IF(numsubstances>0) crunoffd(:)=soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoffd,crunoffd,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(15),subid,i,j)
          cweights(1) = 1.
        ENDIF
      ELSEIF(tdepth<=sdepth(2))THEN   
        !>\li drainage in middle soil layer
        deltah = (soilstate%water(2,j,i)-wp(2)-fc(2))/ep(2) * sthick(2) - (sdepth(2) - tdepth)
        IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN
          deltah = deltah + (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1)
        ENDIF
        IF(deltah>0.)THEN
          runoffd = trc * deltah / sthick(2) * ep(2)
          IF(runoffd > soilstate%water(2,j,i)-wp(2)-fc(2))  runoffd = soilstate%water(2,j,i)-wp(2)-fc(2)
          IF(numsubstances>0) crunoffd(:)=soilstate%conc(:,2,j,i)
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),runoffd,crunoffd,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(16),subid,i,j)
          cweights(2) = 1.
        ENDIF
      ELSE                                  
        !>\li drainage in deepest soil layer
        deltah = (soilstate%water(3,j,i)-wp(3)-fc(3))/ep(3) * sthick(3) - (sdepth(3) - tdepth)
        IF(soilstate%water(3,j,i)-wp(3)-fc(3)-ep(3)>=0.)THEN
          deltah = deltah + (soilstate%water(2,j,i)-wp(2)-fc(2))/ep(2) * sthick(2)
          IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN
            deltah = deltah + (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1)
          ENDIF
        ENDIF
        IF(deltah>0.)THEN
          runoffd = trc * deltah / sthick(3) * ep(3)
          IF(runoffd > soilstate%water(3,j,i)-wp(3)-fc(3))  runoffd = soilstate%water(3,j,i)-wp(3)-fc(3)
          IF(numsubstances>0) crunoffd(:)=soilstate%conc(:,3,j,i)
          CALL remove_water(soilstate%water(3,j,i),numsubstances,soilstate%conc(:,3,j,i),runoffd,crunoffd,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(17),subid,i,j)
          cweights(3) = 1.
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE calculate_tile_drainage

  !>\brief Soil runoff from all soil layers
  !!
  !> \b Reference ModelDescription Chapter Land routines (Soil water - Groundwater runoff)
  !--------------------------------------------------------------------------------
  SUBROUTINE calculate_soil_runoff(i,j,subid,wp,fc,ep,sdepth,sthick,ddepth,soilstate,soilrunoff,csoilrunoff)

    USE MODVAR, ONLY : numsubstances,   &
                       maxsoillayers
    USE HYPEVARIABLES, ONLY : soilrc

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of current subbasin
    INTEGER, INTENT(IN) :: j        !<index of current class
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    REAL, INTENT(IN)    :: wp(maxsoillayers) !<volume below wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers) !<"field capacity" volume (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers) !<effective porosity volume (mm)
    REAL, INTENT(IN)    :: sdepth(maxsoillayers) !<Lower border of soil layers (m)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<Thickness of soil layers (m)
    REAL, INTENT(IN)    :: ddepth         !<Depth of stream, drainagedepth (m)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    REAL, INTENT(OUT)   :: soilrunoff(maxsoillayers)  !<runoff
    REAL, INTENT(OUT)   :: csoilrunoff(numsubstances,maxsoillayers)    !<concentration of runoff
    
    !Local variables
    INTEGER status    !error status of subroutine call
    REAL    deltah,deltah1,deltah2    !groundwater level above drainage level (streamdepth) (m)
    REAL    runoff1,runoff2,runoff3
    REAL    crunoff1(numsubstances),crunoff2(numsubstances),crunoff3(numsubstances)

    !Soil runoff: ditch or local stream
    runoff1 = 0.
    crunoff1 = 0.
    runoff2 = 0.
    crunoff2 = 0.
    runoff3 = 0.
    crunoff3 = 0.
    IF(ddepth>0)THEN
      IF(ddepth<=sdepth(1))THEN       !Drainage level in uppermost layer
        deltah = (soilstate%water(1,j,i)-wp(1)-fc(1))/ep(1) * sthick(1) - (sdepth(1) - ddepth)  !m
        IF(deltah>0.)THEN
          runoff1 = soilrc(1,j,i) * deltah / sthick(1) * ep(1)
          IF(numsubstances>0) crunoff1 = soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoff1,crunoff1,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(9),subid,i,j)
        ENDIF
      ELSEIF(ddepth<=sdepth(2))THEN   !Drainage level in middle layer
        deltah = (soilstate%water(2,j,i)-wp(2)-fc(2))/ep(2) * sthick(2) - (sdepth(2) - ddepth)
        deltah1 = soilstate%water(1,j,i)-wp(1)-fc(1)           
        IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN  
          IF(deltah1>0.) deltah = deltah + deltah1/ep(1) * sthick(1)
        ENDIF
        IF(deltah>0.)THEN
          runoff2 = MIN(soilstate%water(2,j,i)-wp(2)-fc(2), soilrc(2,j,i) * deltah / sthick(2) * ep(2))
          IF(numsubstances>0) crunoff2 = soilstate%conc(:,2,j,i)
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),runoff2,crunoff2,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(10),subid,i,j)
        ENDIF
        IF(deltah1>0.)THEN          
          runoff1 = soilrc(1,j,i) * deltah1   
          IF(numsubstances>0) crunoff1=soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoff1,crunoff1,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(11),subid,i,j)
        ENDIF
      ELSE                                      !Drainage level in deepest layer
        deltah = (soilstate%water(3,j,i)-wp(3)-fc(3))/ep(3) * sthick(3) - (sdepth(3) - ddepth)
        deltah2 = soilstate%water(2,j,i)-wp(2)-fc(2)        
        deltah1 = soilstate%water(1,j,i)-wp(1)-fc(1)         
        IF(soilstate%water(3,j,i)-wp(3)-fc(3)-ep(3)>=0.)THEN
          IF(deltah2>0.) deltah = deltah + deltah2/ep(2) * sthick(2)
          IF(soilstate%water(2,j,i)-wp(2)-fc(2)-ep(2)>=0.)THEN
            IF(deltah1>0.) deltah = deltah + deltah1/ep(1) * sthick(1)
          ENDIF
        ENDIF
        IF(deltah>0.)THEN
          runoff3 = MIN(soilstate%water(3,j,i)-wp(3)-fc(3), soilrc(3,j,i) * deltah / sthick(3) * ep(3))
          IF(numsubstances>0) crunoff3 = soilstate%conc(:,3,j,i)
          CALL remove_water(soilstate%water(3,j,i),numsubstances,soilstate%conc(:,3,j,i),runoff3,crunoff3,status) 
          IF(status.NE.0) CALL error_remove_water(errstring(12),subid,i,j)
        ENDIF
        IF(deltah1>0.)THEN           
          runoff1 = soilrc(1,j,i) * deltah1 
          IF(numsubstances>0) crunoff1=soilstate%conc(:,1,j,i)
          CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),runoff1,crunoff1,status)
          IF(status.NE.0) CALL error_remove_water(errstring(13),subid,i,j)
        ENDIF
        IF(deltah2>0.)THEN           
          deltah = deltah+sdepth(2)-ddepth
          IF(deltah>0.)THEN
            runoff2 = MIN(deltah2,soilrc(2,j,i) * deltah / sthick(2) * ep(2))
            IF(numsubstances>0) crunoff2=soilstate%conc(:,2,j,i)
          ELSE  !(perched watertable)
            runoff2 = soilrc(2,j,i) * deltah2
            IF(numsubstances>0) crunoff2=soilstate%conc(:,2,j,i)
          ENDIF
        ENDIF
        IF(runoff2>0.)THEN
          CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),runoff2,crunoff2,status)
          IF(status.NE.0) CALL error_remove_water(errstring(14),subid,i,j)
        ENDIF
      ENDIF
    ENDIF

    soilrunoff(1) = runoff1
    soilrunoff(2) = runoff2
    soilrunoff(3) = runoff3
    csoilrunoff(:,1) = crunoff1
    csoilrunoff(:,2) = crunoff2
    csoilrunoff(:,3) = crunoff3

  END SUBROUTINE calculate_soil_runoff

  !>\brief Calculate infiltration to soil and surface flow and macropore flow 
  !>due to limited infiltration capacity. 
  !>
  !>\b Reference ModelDescription Chapter Land routines (Soil water - 
  !> Diversion of surface runoff and macropore flow, Infiltration)
  !-----------------------------------------------------------------------------------
  SUBROUTINE infiltration(i,j,isoil,wp,fc,ep,ginfilt,cginfilt,temp,tmin,tmax,  &
       infilt,cinfilt,surfaceflow,csurfaceflow,macroflow,cmacroflow,frozenstate,soilstate)

    USE MODVAR, ONLY : numsubstances, &
                       maxsoillayers, &
                       modeloption, &
                       p_infiltration, &
                       soilpar, &
                       timesteps_per_day, &
                       missing_value
    USE HYPEVARIABLES, ONLY : m_macrate,m_mactrinf,m_mactrsm, &
                              m_srrate,m_bfroznsoil

    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: isoil              !<index of soil type
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<effective porosity volume (mm) (water available for runoff)
    REAL, INTENT(IN)    :: ginfilt            !<gross infiltration; rain+snowmelt (mm/timestep)
    REAL, INTENT(IN)    :: cginfilt(numsubstances)      !<concentration of gross infiltration
    REAL, INTENT(IN)    :: temp            !<current air temperature (degree Celsius)
    REAL, INTENT(IN)    :: tmin            !<current daily minimum air temperature (degree Celsius)
    REAL, INTENT(IN)    :: tmax            !<current daily maximum air temperature (degree Celsius)
    REAL, INTENT(OUT)   :: infilt             !<infiltration (mm/timestep)
    REAL, INTENT(OUT)   :: cinfilt(numsubstances)       !<concentration of infiltration
    REAL, INTENT(OUT)   :: surfaceflow        !<surface runoff due to limited infiltration capacity (mm/timestep)
    REAL, INTENT(OUT)   :: csurfaceflow(numsubstances)  !<concentration of surface flow
    REAL, INTENT(OUT)   :: macroflow          !<macropore flow (mm/timestep)
    REAL, INTENT(OUT)   :: cmacroflow(numsubstances)    !<concentration of macropore flow
    TYPE(snowicestatetype),INTENT(IN) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate   !<Soil states
    
    !Local variables
    REAL help,avail
    REAL macrate,inflowtres,srrate   !infiltration excess runoff and macropore flow parameters
    REAL smtresmm
    REAL tmin_check, tmax_check   !temperature for checking for ice lens
    REAL potinfilt
    REAL cuminfilt  !cumulative infiltration from snow melt during rest of season according to Zhao and Gray [mm]
    REAL t0         !time for this infiltration, opportunity time [hours]

    !>\b Algorithm \n
    !>Set default output values
    infilt       = ginfilt
    cinfilt      = cginfilt
    surfaceflow  = 0.
    csurfaceflow = 0.
    macroflow    = 0.
    cmacroflow   = 0.

    !>If no incoming water; return
    IF(ginfilt == 0.)  RETURN

    !>Set parameter values
    macrate = soilpar(m_macrate,isoil)         !macropore flow fraction (-)
    inflowtres = soilpar(m_mactrinf,isoil)     !threshold for macropore (and surface runoff) flow (mm/timestep)
    smtresmm = (wp(1)+ fc(1))*soilpar(m_mactrsm,isoil)       !soil moisture threshold for macropore (and surface runoff) flow (mm)
    srrate = soilpar(m_srrate,isoil)           !surface runoff flow fraction (-)
    IF(macrate+srrate>1.)THEN
      help = macrate + srrate
      macrate = macrate/help
      srrate = srrate/help
    ENDIF

    !>Calculate surface flow and macropore flow due to limited infiltration capacity 
    avail = ginfilt - inflowtres 
    IF(avail>0. .AND. soilstate%water(1,j,i) > smtresmm) THEN
      macroflow = macrate * avail
      surfaceflow = srrate * avail
      cmacroflow = cginfilt
      csurfaceflow = cginfilt
    ENDIF

    !>Calculate net infiltration
    infilt = ginfilt - macroflow - surfaceflow
    cinfilt = cginfilt
    
    !>Calculate effect of frozen soil on infiltration and surface runoff
    !>Based on Zhao & Gray 1999 Estimate snowmelt infiltration into frozen soils
    !>coded by M.K. MacDonald (27 October 2015)
    IF(modeloption(p_infiltration)==1)THEN
      
      !> Presence of icelens depends on daily maximum and minimum air temperature. 
      tmin_check = tmin
      tmax_check = tmax
      IF(tmin_check == missing_value) tmin_check = temp - 5. !based on winter climate normal data for Winnipeg Richardson Int'l A (M.K.MacDonald)
      IF(tmax_check == missing_value) tmax_check = temp + 5.
      IF((tmin_check<-10. .AND. infilt>=5.) .OR. (tmax_check<0. .AND. soilstate%icelens(j,i)==1)) THEN  
        !> Ice lens restricted infiltration, flow redirected to macroflow & surfaceflow, no infiltration
        soilstate%icelens(j,i) = 1
        surfaceflow = surfaceflow + infilt*srrate/(macrate+srrate)
        macroflow = macroflow + infilt*macrate/(macrate+srrate)
        infilt = 0.
      ELSEIF(soilstate%temp(1,j,i)<=0.) THEN
        !> Frozen soil limited infiltration, no ice lens
        soilstate%icelens(j,i) = 0
        IF(soilstate%water(1,j,i)>=ep(1)+wp(1)+fc(1))THEN
          potinfilt = 0.
        ELSE
          t0 = MAX(1., 0.65*frozenstate%snow(j,i)-5) !opportunity time [hours]
          cuminfilt = soilpar(m_bfroznsoil,isoil) * (0.99**2.92) * ((1.-soilstate%water(1,j,i)/(ep(1)+wp(1)+fc(1)))**1.64) * ((0.-soilstate%temp(1,j,i))/273.15)**(-0.45) * t0**0.44 !infiltration [mm]
          potinfilt = cuminfilt/(t0/(24./timesteps_per_day)) ! potential infiltration [mm/timestep]
          potinfilt = MAX(potinfilt, 0.)
        ENDIF
        IF(potinfilt < infilt)THEN ! reduce infilt to potential and redirected water to macroflow & surfaceflow
         surfaceflow = surfaceflow + (infilt - potinfilt)*srrate/(macrate+srrate)
         macroflow = macroflow + (infilt - potinfilt)*macrate/(macrate+srrate)
         infilt = potinfilt
        ENDIF
      ELSE
        !Unfrozen soils (regular HYPE infiltration)
        soilstate%icelens(j,i) = 0
      ENDIF

      !Concentration is the same for all flows, but may not have been set for all flows above
      cmacroflow = cginfilt
      csurfaceflow = cginfilt
    
    ENDIF

  END SUBROUTINE infiltration

  !>Add infiltration to the upper soillayer soil, including
  !>transfering of IN in infiltration to solid ON in soil
  !------------------------------------------------------------------       
  SUBROUTINE add_infiltration(i,j,iluse,infilt,cinfilt,soilstate)
       
    USE MODVAR, ONLY : i_in, numsubstances
       
    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: iluse              !<index of landuse
    REAL, INTENT(IN)    :: infilt             !<infiltration (mm/timestep)
    REAL, INTENT(INOUT)   :: cinfilt(numsubstances)  !<concentration of infiltration
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate   !<soil states

    !>Add infiltration to the upper soillayer soil, including
    !>transfering of IN in infiltration to solid ON in soil
    IF(infilt>0)THEN
      IF(i_in>0) CALL atmdep_in_loss(iluse,infilt,soilstate%fastN(:,j,i),cinfilt(i_in))   !(Atmospheric) IN moved to fastN 
      CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),infilt,cinfilt)
    ENDIF

  END SUBROUTINE add_infiltration

  !>Calculate infiltration from flooded floodplain to floodplain soil layer 1
  SUBROUTINE flood_infiltration(i,j,pw,infilt,outfilt,cinfilt,soilstate)
       
    USE MODVAR, ONLY : numsubstances
       
    !Argument declaration
    INTEGER, INTENT(IN) :: i            !<index of current subbasin
    INTEGER, INTENT(IN) :: j            !<index of current class 
    REAL, INTENT(IN)    :: pw           !<pore volume (mm)
    REAL, INTENT(IN)    :: infilt       !<available for infiltration (mm/timestep)
    REAL, INTENT(OUT)   :: outfilt      !<actual infiltration (mm/timestep)
    REAL, INTENT(INOUT) :: cinfilt(numsubstances)       !<concentration of infiltration
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states

    !Local variables
    REAL space
    
    outfilt = infilt
    !>Add infiltration to the upper soillayer soil if there is room
    IF(outfilt>0.)THEN
      space = pw - soilstate%water(1,j,i) 
      IF(space>0.)THEN
        IF(space<outfilt) outfilt = space
        CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),outfilt,cinfilt)
      ELSE
        outfilt = 0.
      ENDIF  
    ENDIF

  END SUBROUTINE flood_infiltration

  !>Calculate percolation down through the soil layers considering previous percolation 
  !>and the maximum percolation of the whole timestep
  !>Includes change in concentration of percolating water. 
  !>
  !>\b Reference ModelDescription Chapter Land routines (Soil water - Percolation)
  !----------------------------------------------------------------------
  SUBROUTINE percolation(i,j,isoil,subid,wp,fc,ep,sthick,percflow,cpercflow,soilstate)

    USE MODVAR, ONLY : numsubstances,   &
                       maxsoillayers,   &
                       genpar,soilpar,landpar, &            
                       classdata
    USE HYPEVARIABLES, ONLY : m_perc1,m_perc2,    &
                              m_crate5,   &
                              m_onpercred, m_pppercred

    !Argument declaration
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    INTEGER, INTENT(IN) :: isoil              !<index of soil type
    INTEGER, INTENT(IN) :: subid              !<subbasin id
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<effective porosity volume (mm) (water avaliable for runoff)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<thickness of soil layers (m)
    REAL, INTENT(INOUT) :: percflow(2)        !<percolation (mm/time step)
    REAL, INTENT(INOUT) :: cpercflow(2,numsubstances) !<concentration of percolation (mm/time step)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    
    !Local variables
    INTEGER status            !error status of subroutine
    REAL perc1,perc2          !percolation soillayer 1 and 2 (mm/timestep)
    REAL maxperc1,maxperc2    !maximum percolation from soillayer 1 and maximum percolation to soillayer 2 (mm/timestep)
    REAL cperc(numsubstances) !concentration of percolation
    REAL firstpercflow(2),cfirstpercflow(2,numsubstances)

    !>\b Algorithm \n
    !Start calculate percolation
    firstpercflow = percflow
    cfirstpercflow = cpercflow
    IF(sthick(2)>0)THEN

      !>Calculate limitations for second percolation              
      IF(soilstate%water(1,j,i)-wp(1)>fc(1))THEN
        maxperc1 = MIN((soilstate%water(1,j,i)-wp(1)-fc(1)),soilpar(m_perc1,isoil)-firstpercflow(1))      !Maximum percolation from uppermost soil layer
      ELSE
        maxperc1 = 0.
      ENDIF
      maxperc2 = MIN(wp(3)+fc(3)+ep(3)-soilstate%water(3,j,i),soilpar(m_perc2,isoil)-firstpercflow(2))   !Maximum percolation deep soil layer can accept

      !>Calculate percolation amount of water
      IF(soilstate%water(2,j,i)+maxperc1<=wp(2)+fc(2))THEN
        perc1=maxperc1
        perc2=0.
      ELSEIF(soilstate%water(2,j,i)+maxperc1<=wp(2)+fc(2)+ep(2))THEN
        perc1=maxperc1
        perc2=MIN(soilstate%water(2,j,i)+perc1-wp(2)-fc(2),maxperc2)
      ELSE
        perc2=MIN(soilstate%water(2,j,i)+maxperc1-wp(2)-fc(2),maxperc2)
        IF(maxperc1+soilstate%water(2,j,i)-perc2<=wp(2)+fc(2)+ep(2))THEN
          perc1=maxperc1
        ELSE  !perc2=maxperc2
          perc1=wp(2)+fc(2)+ep(2)-soilstate%water(2,j,i)+perc2
        ENDIF
      ENDIF

      !>Move percolation water to underlaying soillayer and reduce the concentrations:
      IF(perc1>0)THEN
        IF(numsubstances>0) cperc=soilstate%conc(:,1,j,i)
        !>\li Reduce OC, ON and PP concentration of water percolating from upper soillayer
        CALL doc_percolation_reduction(numsubstances,cperc,genpar(m_crate5),soilstate%temp(1,j,i),   &
             soilstate%water(1,j,i),wp(1),wp(1)+fc(1)+ep(1),sthick(1))
        CALL onpp_percolation_reduction(numsubstances,cperc,landpar(m_onpercred,classdata(j)%luse),   &
             landpar(m_pppercred,classdata(j)%luse))  
        !>\li Remove water from upper soillayer and add to second soillayer
        CALL remove_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),perc1,cperc,status) 
        IF(status.NE.0) CALL error_remove_water(errstring(18),subid,i,j)
        CALL add_water(numsubstances,soilstate%water(2,j,i),soilstate%conc(:,2,j,i),perc1,cperc)
        percflow(1) = perc1 + firstpercflow(1)
        IF(numsubstances>0.) cpercflow(1,:) = (cperc(:)*perc1+cfirstpercflow(1,:)*firstpercflow(1))/percflow(1)
      ENDIF
      IF(perc2>0)THEN
        IF(numsubstances>0) cperc=soilstate%conc(:,2,j,i)
        !>\li Reduce OC, ON and PP concentration of water percolating from middle soillayer
        CALL doc_percolation_reduction(numsubstances,cperc,genpar(m_crate5),soilstate%temp(2,j,i),   &
             soilstate%water(2,j,i),wp(2),wp(2)+fc(2)+ep(2),sthick(2))
        CALL onpp_percolation_reduction(numsubstances,cperc,landpar(m_onpercred,classdata(j)%luse),   &
             landpar(m_pppercred,classdata(j)%luse))
        !>\li Remove water from middle soillayer and add to third soillayer
        CALL remove_water(soilstate%water(2,j,i),numsubstances,soilstate%conc(:,2,j,i),perc2,cperc,status) 
        IF(status.NE.0) CALL error_remove_water(errstring(19),subid,i,j)
        CALL add_water(numsubstances,soilstate%water(3,j,i),soilstate%conc(:,3,j,i),perc2,cperc)
        percflow(2) = perc2 + firstpercflow(2)
        IF(numsubstances>0.) cpercflow(2,:) = (cperc(:)*perc2+cfirstpercflow(2,:)*firstpercflow(2))/percflow(2)
      ENDIF

    ENDIF

  END SUBROUTINE percolation

  !>\brief Add macropore water flow to soil layer with groundwater level. 
  !>
  !!If this soillayer can't take all water the rest is added to the
  !!soillayer(s) above. 
  !>
  !>\b Reference ModelDescription Chapter Land routines (Soil water - Macropore flow)
  !---------------------------------------------------------------------
  SUBROUTINE add_macropore_flow(i,j,macroflow,cmacroflow,wp,fc,ep,pw,sdepth,sthick,slmacroflows,soilstate)

    USE MODVAR, ONLY : numsubstances,   &
                       maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: i                  !<index of current subbasin
    INTEGER, INTENT(IN) :: j                  !<index of current class 
    REAL, INTENT(IN)    :: macroflow          !<macropore flow to be added (mm/timestep)
    REAL, INTENT(IN)    :: cmacroflow(numsubstances) !<concentration of macropore flow 
    REAL, INTENT(IN)    :: wp(maxsoillayers)  !<wilting point volume (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)  !<"field capacity" volume (mm) (water available for evaporation but not for runoff)
    REAL, INTENT(IN)    :: ep(maxsoillayers)  !<"effective porosity" volume (mm) (water avaliable for runoff)
    REAL, INTENT(IN)    :: pw(maxsoillayers)  !<total pore volume (mm)
    REAL, INTENT(IN)    :: sdepth(maxsoillayers) !<lower border of soil layers (m)
    REAL, INTENT(IN)    :: sthick(maxsoillayers) !<thickness of soil layers (m)
    REAL, INTENT(OUT)   :: slmacroflows(maxsoillayers) !<macropore flow to each soil layer (mm/timestep)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    
    !Local variables
    REAL gwat                   !groundwater table (m) (negative)
    REAL newsoil                
    REAL fill,fill2

    !Start subroutine
    !>\b Algorithm \n
    !>If no macropore flow: return
    slmacroflows = 0
    IF(macroflow<=0) RETURN   

    !>Find soillayer of groundwater table
    CALL calculate_groundwater_table(soilstate%water(:,j,i),wp,fc,ep,pw,sdepth(:),sthick(:),gwat)
    !>If groundwater table in soillayer three:
    IF(-gwat>sdepth(2)) THEN
      !>\li Check if soillayer three has room and add the water there is room for
      newsoil = soilstate%water(3,j,i) + macroflow
      IF(newsoil > pw(3)) THEN
        fill = newsoil - pw(3)
        IF(numsubstances>0) soilstate%conc(:,3,j,i) = (soilstate%water(3,j,i)*soilstate%conc(:,3,j,i) + (macroflow - fill)*cmacroflow)/pw(3)
        soilstate%water(3,j,i) = pw(3)
        slmacroflows(3) = macroflow - fill
      ELSE
        fill = 0.
        IF(numsubstances>0) soilstate%conc(:,3,j,i) = (soilstate%water(3,j,i)*soilstate%conc(:,3,j,i) + macroflow*cmacroflow)/newsoil
        soilstate%water(3,j,i) = newsoil
        slmacroflows(3) = macroflow
      ENDIF
      !>\li If too much water, check if soillayer 2 has room and add the water there is room for
      IF(fill > 0.) THEN
        newsoil = soilstate%water(2,j,i) + fill
        IF(newsoil > pw(2)) THEN
          fill2 = newsoil - pw(2)
          IF(numsubstances>0) soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + (fill-fill2)*cmacroflow)/pw(2)
          soilstate%water(2,j,i) = pw(2)
          slmacroflows(2) = fill - fill2
        ELSE
          fill2 = 0.
          IF(numsubstances>0) soilstate%conc(:,2,j,i) = ((newsoil-fill)*soilstate%conc(:,2,j,i) + fill*cmacroflow)/newsoil
          soilstate%water(2,j,i) = newsoil
          slmacroflows(2) = fill
        ENDIF
        !>\li If still too much water add the rest to soillayer 1
        IF(fill2 > 0.) THEN
          newsoil = soilstate%water(1,j,i) + fill2
          IF(numsubstances>0) soilstate%conc(:,1,j,i) = (soilstate%water(1,j,i)*soilstate%conc(:,1,j,i) + fill2*cmacroflow)/newsoil
          soilstate%water(1,j,i) = newsoil
          slmacroflows(1) = fill2
        ENDIF
      ENDIF

    !>Elseif groundwater table in soillayer two:
    ELSEIF(-gwat>sdepth(1)) THEN
      newsoil = soilstate%water(2,j,i) + macroflow
      !>\li Check if soillayer 2 has room and add the water there is room for
      IF(newsoil > pw(2)) THEN
        fill = newsoil - pw(2)
        IF(numsubstances>0) soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + (macroflow - fill)*cmacroflow)/pw(2)
        soilstate%water(2,j,i) = pw(2)
        slmacroflows(2) = macroflow - fill
      ELSE
        fill = 0.
        IF(numsubstances>0) soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + macroflow*cmacroflow)/newsoil
        soilstate%water(2,j,i) = newsoil
        slmacroflows(2) = macroflow
      ENDIF
      !>\li If too much water add the rest to soillayer 1
      IF(fill > 0.) THEN
        CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),fill,cmacroflow)
        slmacroflows(1) = fill
      ENDIF

    !>Elseif groundwater table in soillayer one:
    ELSE
      !>\li Add macropore flow to soillayer 1
      CALL add_water(numsubstances,soilstate%water(1,j,i),soilstate%conc(:,1,j,i),macroflow,cmacroflow)
      slmacroflows(1) = macroflow
    ENDIF

  END SUBROUTINE add_macropore_flow

  !>Subroutine for calculation of ground water table level (metres
  !>above land surface)
  !!
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions - Diagnostic variables)
  !--------------------------------------------------------------------
  SUBROUTINE calculate_groundwater_table(soil,wpvol,fcvol,&
       epvol,totvol,soildep,thickness,gwat)
    USE MODVAR, ONLY : maxsoillayers

    !Argument declarations
    REAL, INTENT(IN)     :: soil(maxsoillayers)     !<soil moisture (mm)
    REAL, INTENT(IN)     :: wpvol(maxsoillayers)    !<wilting point volume in all layers (mm) 
    REAL, INTENT(IN)     :: fcvol(maxsoillayers)    !<field capacity volume in all layers (mm)
    REAL, INTENT(IN)     :: epvol(maxsoillayers)    !<effective porosity volume in all layers (mm) 
    REAL, INTENT(IN)     :: totvol(maxsoillayers)   !<maximum volume of water in soil layer (pore volume) (mm) 
    REAL, INTENT(IN)     :: soildep(maxsoillayers)  !<depth of soil layers (m)
    REAL, INTENT(IN)     :: thickness(maxsoillayers)  !<thickness of soil layers (m)
    REAL, INTENT(OUT)    :: gwat !<ground water table (m), measured from land surface and up
    
    !Local parameters
    REAL, PARAMETER :: mindiff = 0.000005 !safeguard in choosing groundwater table from soil layer tables
    
    !Local variables
    REAL gwat1,gwat2,gwat3      !groundwater table of each soil layer

    gwat = 0.

    IF(thickness(2)>0)THEN
      IF(soil(1)-wpvol(1)-fcvol(1)>0.0)THEN
        gwat1 = (soil(1) - totvol(1))/epvol(1) * thickness(1)     !negative (m)
        IF(gwat1 > 0) gwat1 = (soil(1) - totvol(1)) * 0.001     !100% porositet above land surface
      ELSE
        gwat1 = -soildep(1)
      ENDIF
      IF(soil(2)-wpvol(2)-fcvol(2)>0.0)THEN
        gwat2 = (soil(2)-totvol(2))/epvol(2) * thickness(2) - soildep(1)     !negative (m)
      ELSE
        gwat2 = -soildep(2)
      ENDIF
      IF(soil(3)-wpvol(3)-fcvol(3)>0.0)THEN
        gwat3 = (soil(3)-totvol(3))/epvol(3) * thickness(3) - soildep(2)     !negative (m)
      ELSE
        gwat3 = -soildep(3)
      ENDIF
      IF(-gwat3>soildep(2)+mindiff)THEN     !Find ground water table as lowest level with filled pores below
        gwat = gwat3
      ELSEIF(-gwat2>soildep(1)+mindiff)THEN
        gwat = gwat2
      ELSE
        gwat = gwat1
      ENDIF
    ELSE
      IF(soil(1)-wpvol(1)-fcvol(1)>0.0)THEN
        gwat1 = (soil(1) - totvol(1))/epvol(1) * thickness(1)     !negative (m)
        IF(gwat1 > 0) gwat1 = (soil(1) - totvol(1)) * 0.001     !100% porositet above land surface
      ELSE
        gwat1 = -soildep(1)
      ENDIF
      gwat = gwat1
    ENDIF

  END SUBROUTINE calculate_groundwater_table

  !>Calculation of soil temperature in soil layers and deep soil
  !>
  !> \b Reference ModelDescription Chapter Land routines (Snow routines - Soil temperature and snow depth)
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_soiltemp(n,airtemp,snowdepth,soilmemdeep,soilmemlayer,deeptemp,soiltemp)

    USE MODVAR, ONLY : timesteps_per_day
    
    !Argument declarations
    INTEGER, INTENT(IN)  :: n               !<number of soil layers
    REAL, INTENT(IN)     :: airtemp         !<air temperature (degree Celcius) 
    REAL, INTENT(IN)     :: snowdepth       !<snow depth (cm)
    REAL, INTENT(IN)     :: soilmemdeep     !<parameter, temperature memory of deep soil (days)
    REAL, INTENT(IN)     :: soilmemlayer(n) !<parameter, temperature memory of soil layer (timesteps)
    REAL, INTENT(INOUT)  :: deeptemp        !<deep soil temperature (degree Celcius)
    REAL, INTENT(INOUT)  :: soiltemp(n)     !<soil temperature (degree Celcius)
    
    !Local parameters
    REAL, PARAMETER :: spfrost = 10.      !coefficient of equation for weight of air temperature in soil temperature calculation (days/cm)
    REAL, PARAMETER :: weightdeep = 0.001 !weight of deep soil temperature for soil temperature calculation (dimensionless)
    
    !Local variables
    INTEGER k        !layer index
    REAL weightair   !weight of air temperature for soil temperature calculation (dimensionless)

    !> \b Algorithm \n
    !>Calculate deep soil temperature
    weightair = 1./ ((soilmemdeep + spfrost * snowdepth)*timesteps_per_day)
    CALL calculate_weighted_temperature(airtemp,weightair,0.0,0.0,deeptemp)
    !>Calculate soil layer temperature for each soil layer
    DO k = 1,n
      weightair = 1./ (soilmemlayer(k) + spfrost * snowdepth*timesteps_per_day)
      CALL calculate_weighted_temperature(airtemp,weightair,deeptemp,weightdeep,soiltemp(k))
    ENDDO

   END SUBROUTINE calculate_soiltemp

  !>Calculation of soil temperature as an average of three temperatures: 
  !>air temperature, deep soil temperature and previous soil temperature.
  !>
  !> \b Reference ModelDescription Chapter Land routines (Snow routines - Soil temperature and snow depth)
  !--------------------------------------------------------------------------------
  SUBROUTINE calculate_weighted_temperature(temp1,weight1,temp2,weight2,soiltemp)

    !Argument declarations
    REAL, INTENT(IN)     :: temp1     !<air temperature (degree Celcius) 
    REAL, INTENT(IN)     :: weight1   !<weight of temp1 (dimensionless)
    REAL, INTENT(IN)     :: temp2     !<temperature of deep soil (ca 1m) (degree Celcius)
    REAL, INTENT(IN)     :: weight2   !<weight of temp2 (dimensionless)
    REAL, INTENT(INOUT)  :: soiltemp  !<soil layer temperature (degree Celcius)

    !> \b Algorithm \n
    !>Calculate weighted temperature
    soiltemp = soiltemp * (1. - weight1 - weight2) + &
         temp1    * weight1  +   &
         temp2    * weight2

  END SUBROUTINE calculate_weighted_temperature

  !>Calculation of soil frost depth depending on temperature of soil
  !>
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions - Diagnostic variables)
  !-------------------------------------------------------------------------------------
  SUBROUTINE calculate_frostdepth(fc,cfrost,sfrost,soil,frostdepth,soiltemp,thickness)
  
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    REAL, INTENT(IN)  :: fc           !<water content at field capacity (mm)
    REAL, INTENT(IN)  :: cfrost       !<soil frost coefficient, land use dependent (cm/degree)
    REAL, INTENT(IN)  :: sfrost       !<soil frost coefficient, soil type dependent (cm/degree)
    REAL, INTENT(IN)  :: soil         !<soil water (mm) 
    REAL, INTENT(OUT) :: frostdepth   !<depth of soil frost (cm)
    REAL, INTENT(IN)  :: soiltemp(:)  !<soil temperature (degree Celcius)
    REAL, INTENT(IN)  :: thickness(:) !<soil layer thickness (m)
    
    !Local variables 
    INTEGER dim
    REAL uppertemp  !average temperature upper soil (ca 50 cm)

    !> \b Algorithm \n
    !>If soil frost parameters are set:
    IF(cfrost>0 .AND. sfrost>0)THEN
      dim = SIZE(soiltemp)
      !> \li Calculate average temperature of upper two soil layers
      IF(dim==1)THEN
        uppertemp = soiltemp(1)
      ELSE
        uppertemp = (soiltemp(1)*thickness(1)+soiltemp(2)*thickness(2))/(thickness(1)+thickness(2))
      ENDIF
      !> \li If temperature is negative, calculate soil frost depth
      IF(uppertemp<0)THEN
        frostdepth = cfrost * sfrost * uppertemp * fc / soil
      ELSE
        frostdepth = 0.
      ENDIF
    ELSE
    !Else soil frost is set to missing
      frostdepth = missing_value
    ENDIF

  END SUBROUTINE calculate_frostdepth

  !> \brief Calculation of soil moisture deficit (mm left to field capacity)
  !>in top two layers
  !>
  !> \b Reference ModelDescription Chapter Land routines (Basic assumptions - Diagnostic variables)
  !---------------------------------------------------------------------
  SUBROUTINE calculate_soil_moisture_deficit(soil,wpvol,fcvol,thickness,smdef)
    USE MODVAR, ONLY : maxsoillayers

    !Argument declaration
    REAL, INTENT(IN)     :: soil(maxsoillayers)       !<soil moisture  (mm)
    REAL, INTENT(IN)     :: wpvol(maxsoillayers)      !<wilting point volume in all layers (mm) 
    REAL, INTENT(IN)     :: fcvol(maxsoillayers)      !<"field capacity" volume in all layers (mm)
    REAL, INTENT(IN)     :: thickness(maxsoillayers)  !<thickness of soil layers (m)
    REAL, INTENT(OUT)    :: smdef                     !<soil moisture deficit (mm)
    
    !Local variables
    REAL smdef1,smdef2      !soil moisture deficit of each soil layer

    !> \b Algorithm \n
    !>Initate soil moisture deficit to zero
    smdef = 0.

    !>Calculate soil moisture deficit in top soil layer, add to total
    smdef1 = fcvol(1)+wpvol(1)-soil(1)
    IF(smdef1>0.) smdef = smdef + smdef1

    !>Calculate soil moisture deficit in second soil layer, add to total
    IF(thickness(2)>0)THEN
      smdef2 = fcvol(2)+wpvol(2)-soil(2)
      IF(smdef2>0.) smdef = smdef + smdef2
    ENDIF

  END SUBROUTINE calculate_soil_moisture_deficit

END MODULE SOIL_PROCESSES
