!> \file sw_proc.f90
!> Contains module surfacewater_processes.

!>Lake and river water related subroutines in HYPE
MODULE SURFACEWATER_PROCESSES

  !Copyright 2012-2017 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
  !------------------------------------------------------------------------
  USE STATETYPE_MODULE
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water,       &
                                          error_remove_water, &
                                          add_water
  USE SOIL_PROCESSES, ONLY : calculate_snowmelt,  &
                             calculate_snowdepth,  &
                             snowalbedo_function,   &
                             latentheat_tempfunction
       
  !Also uses hypevariables, modvar,general_functions
  IMPLICIT NONE
  PRIVATE
  !--------------------------------------
  !Private procedures 
  !--------------------------------------
  ! set_rivertemp
  ! calc_qbank 
  ! update_qbank
  ! get_current_lake_outflow_parameters
  ! get_current_production_flow
  ! apply_seasonal_factor_on_production_flow
  ! adjust_threshold_for_seasonal_variation
  ! get_current_rating_parameters
  ! calculate_lake_outlet_outflow
  ! calculate_maxprod_outflow
  ! average_flow_rating_curve
  ! recalculate_branched_flow
  ! riverice_riverwater_interaction
  ! calculate_snow_on_ice
  ! calculate_lakeice_lakewater_interaction
  ! calculate_icedepth
  ! calculate_T2_transfer
  ! calculate_T2_transfer_upper2lower
  ! calculate_floodplain_volume, &
  ! calculate_floodplain_equilibriumlevel, &
  !-------------------------------------
  PUBLIC :: calculate_landarea_riverlength,  &
            add_precipitation_to_river, &
            add_precipitation_to_floodplain, &
            calculate_river_evaporation, &
            calculate_floodplain_evaporation, &
            calculate_actual_lake_evaporation, &
            sum_upstream_area, &
            set_general_rating_k,  &
            calculate_water_temperature, &
            set_water_temperature,  &
            calculate_river_characteristics, &
            translation_in_river, &
            point_abstraction_from_main_river_inflow, &
            point_abstraction_from_main_river, &
            point_abstraction_from_outlet_lake, &
            calculate_ilake_outflow, &
            calculate_outflow_from_outlet_lake, &
            calculate_flow_from_outlet_lake_waterstage, &
            remove_outflow_from_lake, &
            calculate_flow_within_lake, &
            calculate_olake_waterstage, &
            calculate_regamp_adjusted_waterstage, &
            calculate_branched_flow,  &
            calculate_branched_flow_new, &
            set_lake_outlets, &
            calculate_lake_volume, &
            calculate_lake_epilimnion_depth, &
            T2_processes_in_river, &
            T2_processes_in_lake, &
            ice_processes_in_lake, &
            ice_processes_in_river, &
            add_T2_concentration_in_precipitation_on_water, &
            get_rivertempvol, &
            calculate_floodplain_waterlevel, &
            calculate_waterbody_floodplain_interflow

  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(11)  !error message for location of remove_water call
  PARAMETER (errstring = (/'evapotranspiration lake, less than lake volume',    &   !1
                           'evapotranspiration lake, more than lake volume',    &   !2
                           'evapotranspiration lake, slowlake part used   ',    &   !3
                           'lake outflow, no NPC simulation               ',    &   !4 
                           'lake outflow, no division in parts (NPC sim)  ',    &   !5
                           'lake outflow, from fastlake part              ',    &   !6
                           'lake outflow, from slowlake part              ',    &   !7
                           'flow between fast- and slowlake parts         ',    &   !8
                           'flow between waterbody and floodplain         ',    &   !9
                           'flow between floodplain and waterbody         ',    &   !10
                           'evapotranspiration river                      ' /))     !11

  CONTAINS

  !>\brief Calculate land area of subbasins and determine riverlength 
  !!for local streams and main rivers.
  !!The landare does not include floodplains (dry or flooded).
  !>
  !\b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions)
  !----------------------------------------------------------------------------
  SUBROUTINE calculate_landarea_riverlength(nsub,landarea,rivlength)

    USE MODVAR, ONLY : basin,classbasin, &
                       slc_ilake,slc_olake, &
                       slc_lriver,slc_mriver, &
                       conductflood, &
                       floodindex,flooding
    
    !Argument declarations
    INTEGER, INTENT(IN) :: nsub   !<Number of subbasins
    REAL, INTENT(OUT)   :: landarea(nsub)       !<land area [m2]
    REAL, INTENT(OUT)   :: rivlength(2,nsub)    !<river length [m]
   
    !Local variables
    INTEGER i
    REAL default_rivlen(nsub)

    !>\b Algorithm \n
    !>Calculate land area of subbasin
    landarea = basin%area
    IF(slc_ilake>0)  landarea = landarea - basin(:)%area * classbasin(:,slc_ilake)%part
    IF(slc_lriver>0) landarea = landarea - basin%area * classbasin(:,slc_lriver)%part
    IF(.NOT.conductflood)THEN
      IF(slc_olake>0)  landarea = landarea - basin%area * classbasin(:,slc_olake)%part
      IF(slc_mriver>0) landarea = landarea - basin%area * classbasin(:,slc_mriver)%part
    ELSE
      DO i = 1,nsub
        IF(floodindex(i)>0)THEN
          IF(flooding(floodindex(i))%fpfmr>0.)THEN
            landarea(i) = landarea(i) - basin(i)%area * classbasin(i,slc_mriver)%part * flooding(floodindex(i))%fpfmr
          ELSE
            IF(slc_mriver>0) landarea(i) = landarea(i) - basin(i)%area * classbasin(i,slc_mriver)%part
          ENDIF
          IF(flooding(floodindex(i))%fpfol>0.)THEN
            landarea(i) = landarea(i) - basin(i)%area * classbasin(i,slc_olake)%part * flooding(floodindex(i))%fpfol
          ELSE
            IF(slc_olake>0) landarea(i) = landarea(i) - basin(i)%area * classbasin(i,slc_olake)%part
          ENDIF
        ELSE
          IF(slc_mriver>0) landarea(i) = landarea(i) - basin(i)%area * classbasin(i,slc_mriver)%part
          IF(slc_olake>0) landarea(i) = landarea(i) - basin(i)%area * classbasin(i,slc_olake)%part
        ENDIF
      ENDDO
    ENDIF
    !Calculate square root of landarea
    DO i = 1,nsub
      IF(landarea(i)<0) landarea(i)=0.   !Safe for all lake subbasin (1-1<0)
      default_rivlen(i) = SQRT(landarea(i))
    ENDDO
    !>Set river length from GeoData, or if zero use default value, i.e. square root of land area
    rivlength(1,:) = basin(:)%rivlen(1)  !local river length
    WHERE(rivlength(1,:)==0) rivlength(1,:) = default_rivlen(:)
    rivlength(2,:) = basin(:)%rivlen(2)  !main river length
    WHERE(rivlength(2,:)==0) rivlength(2,:) = default_rivlen(:)

  END SUBROUTINE calculate_landarea_riverlength

  !>\brief Add precipitation to river, according to volume of watercourse elements.
  !>
  !\b Reference ModelDescription Chapter Rivers and lakes (Rivers - Common river processes)
  !----------------------------------------------------------------------------
  SUBROUTINE add_precipitation_to_river(i,pooltype,area,prec,cprec,riverstate)

    USE MODVAR, ONLY : numsubstances
    USE HYPEVARIABLES, ONLY : ttpart,ttstep

    !Argument declarations
    INTEGER, INTENT(IN) :: i                          !<index of subbasin
    INTEGER, INTENT(IN) :: pooltype                   !<rivertype: 1=lriver, 2=mriver
    REAL, INTENT(IN)    :: area                       !<river area (m2)
    REAL, INTENT(IN)    :: prec                       !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: cprec(numsubstances)       !<concentration of precipitation
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state

    !Local variables
    INTEGER l
    REAL precm3
    REAL totvol
    REAL waterfrac

    !>\b Algorithm \n
    precm3 = prec * 1.E-3 * area
    !>Calculate total volume of river to use fractions of river water i different compartment
    totvol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))

    !>Add precipitation to river watercourse for each compartment in relation to its volume fraction
    IF(totvol>0)THEN
      IF(riverstate%water(pooltype,i)>0)THEN
        waterfrac = riverstate%water(pooltype,i)/totvol
        CALL add_water(numsubstances,riverstate%water(pooltype,i),riverstate%conc(:,pooltype,i),waterfrac*precm3,cprec)
      ENDIF
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
          CALL add_water(numsubstances,riverstate%qqueue(l,pooltype,i),riverstate%cqueue(:,l,pooltype,i),waterfrac*precm3,cprec)
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that remaining outflow will be correct
          CALL add_water(numsubstances,riverstate%qqueue(l,pooltype,i),riverstate%cqueue(:,l,pooltype,i),waterfrac*precm3,cprec)
        ENDIF
      ENDIF
    ELSE
      !>If no river volume add all precipitation to river water compartment
      riverstate%water(pooltype,i) = precm3
    ENDIF  

  END SUBROUTINE add_precipitation_to_river

  !>\brief Add precipitation to river floodplain
  !>
  !\b Reference ModelDescription Chapter Rivers and lakes (Floodplains)
  !----------------------------------------------------------------------------
  SUBROUTINE add_precipitation_to_floodplain(i,pooltype,area,prec,cprec,miscstate,load)

    USE MODVAR, ONLY : numsubstances

    !Argument declarations
    INTEGER, INTENT(IN) :: i                          !<index of subbasin
    INTEGER, INTENT(IN) :: pooltype                   !<type: 1=mriver, 2=olake
    REAL, INTENT(IN)    :: area                       !<flooded area (m2)
    REAL, INTENT(IN)    :: prec                       !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: cprec(numsubstances)       !<concentration of precipitation
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate    !<Floodplain state
    REAL, INTENT(OUT)   :: load(numsubstances)        !<load of precipitation

    !Local variables
    REAL precm3   ![m3]

    precm3 = prec * 1.E-3 * area
    load = 0.

    !Add precipitation to river watercourse
    CALL add_water(numsubstances,miscstate%floodwater(pooltype,i),miscstate%cfloodwater(:,pooltype,i),precm3,cprec)
      
    !Calculate load
    IF(numsubstances>0) load = cprec*prec*area

  END SUBROUTINE add_precipitation_to_floodplain

  !>\brief Calculate and remove evaporation from river
  !>
  !>\b Reference ModelDescription Chapters Rivers and lakes (Rivers - Common river processes)
  !> and Processes above ground (Evaporation)
  !----------------------------------------------------------------------------------------
  SUBROUTINE calculate_river_evaporation(i,j,pooltype,numsubst,area,temp,epot,evap,cevap,riverstate)

    USE MODVAR, ONLY : basin,classdata, &
                       genpar,    &
                       landpar,   &
                       cwater,    &
                       i_t1,i_t2
    USE HYPEVARIABLES, ONLY : m_ttmp,ttpart,ttstep, &
                              m_T1evap

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: j        !<class index
    INTEGER, INTENT(IN) :: pooltype !<river type (local or main)
    INTEGER, INTENT(IN) :: numsubst !<number of substances modelled
    REAL, INTENT(IN)    :: area     !<river area (m2)
    REAL, INTENT(IN)    :: temp     !<air temperature
    REAL, INTENT(IN)    :: epot     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubst) !<concentration in evapotranspiration (eg. mg/L)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<Lake state

    !Local variables
    INTEGER k,l     !loop-variable, substance/queue
    INTEGER status  !error status of subroutine
    REAL tt         !threshold temperature for evaporation (C)
    REAL evapm3     !actual evaporation in m3
    REAL totvol     !total river watercourse volume (m3)
    REAL waterfrac  !fraction of water to be removed
    
    !>\b Algorithm \n
    !>Set default values output variables (zero evaporation)
    evap = 0.
    cevap = 0.

    !Set local parameter
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation

    IF(temp > tt) THEN
      !>If temperature is above threshold river evaporation is potential
      evapm3 = epot*area*1.E-3
      !>Set concentration of evaporation
      DO k=1,numsubst
        cevap(k) = 0.
        IF(k==i_t1) cevap(k) = genpar(m_T1evap) * riverstate%conc(k,pooltype,i)
        !For T2==temperature, set the evaporation concentration = L/c = latent heat of vaporization divided by heat capacity
        IF(k==i_t2) cevap(k) = 1000. * latentheat_tempfunction(riverstate%conc(k,pooltype,i)) / cwater ! 1000 * [MJ/kg] / [KJ/kg/C] = C
      ENDDO

      !>Calculate total river volume to use fractions of river water i different compartments
      totvol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
      IF(totvol<=0.) RETURN

      !>Check if enough water is available for evaporation in each compartment
      IF(evapm3<totvol)THEN
        !>Remove evaporation from river watercourse compartments
        waterfrac = riverstate%water(pooltype,i)/totvol
        CALL remove_water(riverstate%water(pooltype,i),numsubst,riverstate%conc(:,pooltype,i),waterfrac*evapm3,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(11),basin(i)%subid,i,j)
        DO l = 1,ttstep(pooltype,i)
          IF(riverstate%qqueue(l,pooltype,i)>0)THEN
            waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
            CALL remove_water(riverstate%qqueue(l,pooltype,i),numsubst,riverstate%cqueue(:,l,pooltype,i),waterfrac*evapm3,cevap,status)
            IF(status.NE.0) CALL error_remove_water(errstring(11),basin(i)%subid,i,j)
          ENDIF
        ENDDO
        IF(ttpart(pooltype,i)>0)THEN
          l = ttstep(pooltype,i) + 1
          IF(riverstate%qqueue(l,pooltype,i)>0)THEN
            waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that pool get correct concentration change
            CALL remove_water(riverstate%qqueue(l,pooltype,i),numsubst,riverstate%cqueue(:,l,pooltype,i),waterfrac*evapm3,cevap,status)
            IF(status.NE.0) CALL error_remove_water(errstring(11),basin(i)%subid,i,j)
          ENDIF
        ENDIF
        evap = epot
      ELSE
        !>If less water than wanted, remove last traces of substance with the evaporation
        evapm3 = totvol
        riverstate%water(pooltype,i) = 0.
        IF(numsubst>0.) riverstate%conc(:,pooltype,i) = 0.
        DO l = 1,ttstep(pooltype,i)
          riverstate%qqueue(l,pooltype,i) = 0.
          IF(numsubst>0.) riverstate%cqueue(:,l,pooltype,i) = 0.
        ENDDO
        IF(ttpart(pooltype,i)>0)THEN
          l = ttstep(pooltype,i) + 1
          riverstate%qqueue(l,pooltype,i) = 0.
          IF(numsubst>0.) riverstate%cqueue(:,l,pooltype,i) = 0.
        ENDIF
        evap = evapm3/area*1000.
      ENDIF   
    ENDIF

  END SUBROUTINE calculate_river_evaporation

  !>\brief Calculate and remove evaporation from floodplain
  !>
  !\b Reference ModelDescription Chapter Rivers and Lakes (Rivers - Common river processes)
  !----------------------------------------------------------------------------------------
  SUBROUTINE calculate_floodplain_evaporation(i,j,pooltype,numsubst,area,temp,epot,evap,cevap,miscstate)

    USE MODVAR, ONLY : basin,classdata, &
                       genpar,   &
                       landpar,  &
                       i_t1
    USE HYPEVARIABLES, ONLY : m_ttmp, &
                              m_T1evap

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: j        !<class index
    INTEGER, INTENT(IN) :: pooltype !<type 1= main river, 2=olake
    INTEGER, INTENT(IN) :: numsubst !<number of substances modelled
    REAL, INTENT(IN)    :: area    !<floodplain area (m2)
    REAL, INTENT(IN)    :: temp     !<air temperature
    REAL, INTENT(IN)    :: epot     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubst) !<concentration in evapotranspiration (eg. mg/L)
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<Floodplain state

    !Local variables
    INTEGER k       !loop-variable, substance/queue
    INTEGER status  !error status of subroutine
    REAL tt         !threshold temperature for evaporation (C)
    REAL evapm3     !actual evaporation in m3
    REAL totvol     !total river watercourse volume (m3)
    
    !Default values output variables
    evap = 0.
    cevap = 0.

    !Available water
    totvol = miscstate%floodwater(pooltype,i)
    IF(totvol<=0.) RETURN

    !Set local parameter
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation

    IF(temp > tt) THEN
      !Potential evaporation is default for temperature above threshold
      evapm3 = epot*area*1.E-3
      !Set concentration of evaporation (T1)
      DO k=1,numsubst
        cevap(k) = 0.
        IF(k==i_t1) cevap(k) = genpar(m_T1evap) * miscstate%cfloodwater(k,pooltype,i)
!        IF(k==i_t2) cevap(k) = 1000. * latentheat_tempfunction(riverstate%conc(k,pooltype,i)) / cwater ! 1000 * [MJ/kg] / [KJ/kg/C] = C
      ENDDO

      !Check if enough water is available for evaporation
      IF(evapm3<totvol)THEN
        !Remove evaporation from river watercourse
        CALL remove_water(miscstate%floodwater(pooltype,i),numsubst,miscstate%cfloodwater(:,pooltype,i),evapm3,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,j)
        evap = epot
      ELSE
        !Remove last traces of substance
        evapm3 = totvol
        miscstate%floodwater(pooltype,i) = 0.
        IF(numsubst>0) miscstate%cfloodwater(:,pooltype,i) = 0.
        evap = evapm3/area*1000.
      ENDIF   
    ENDIF

  END SUBROUTINE calculate_floodplain_evaporation

  !>\brief Calculate total volume and mean T2 temperature concentration in river
  !----------------------------------------------------------------------------------------
  SUBROUTINE get_rivertempvol(i,pooltype,riverstate,meanrivertemp,totrivervol)

    USE MODVAR, ONLY : i_t2
    USE HYPEVARIABLES, ONLY : ttpart,ttstep
  
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: pooltype !<river type (local or main)
    TYPE(riverstatetype),INTENT(IN) :: riverstate  !<River state
    REAL, INTENT(OUT)   :: meanrivertemp  !<temperature of river water
    REAL, INTENT(OUT)   :: totrivervol    !<volume of river water
    
    INTEGER l
    
    !Total volume in all river elements (translation boxes and river volume)
    totrivervol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
    
    IF(totrivervol.GT.0.)THEN
      !Weighted average T2 concentration
      meanrivertemp = riverstate%conc(i_t2,pooltype,i) * riverstate%water(pooltype,i)/totrivervol
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          meanrivertemp = meanrivertemp + riverstate%cqueue(i_t2,l,pooltype,i) * riverstate%qqueue(l,pooltype,i)/totrivervol
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          meanrivertemp = meanrivertemp + ttpart(pooltype,i) * riverstate%cqueue(i_t2,l,pooltype,i) * riverstate%qqueue(l,pooltype,i)/totrivervol
        ENDIF
      ENDIF
    ELSE
      meanrivertemp = 0.
    ENDIF

  END SUBROUTINE get_rivertempvol
  
  !>\brief Set a T2 temperature concentration to all river elements
  !----------------------------------------------------------------------------------------
  SUBROUTINE set_rivertemp(i,pooltype,riverstate,meanrivertemp)

    USE MODVAR, ONLY : i_t2
    USE HYPEVARIABLES, ONLY : ttpart,ttstep
  
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: pooltype !<river type (local or main)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state
    REAL, INTENT(IN)   :: meanrivertemp !<temperature of river
    
    INTEGER l
    
    !Riverbox
    IF(riverstate%water(pooltype,i).GT.0.)THEN
      riverstate%conc(i_t2,pooltype,i) = meanrivertemp
    ELSE
      riverstate%conc(i_t2,pooltype,i) = 0.
    ENDIF
      
    !Translation boxes
    DO l = 1,ttstep(pooltype,i)
      IF(riverstate%qqueue(l,pooltype,i)>0)THEN
        riverstate%cqueue(i_t2,l,pooltype,i) = meanrivertemp
      ELSE
        riverstate%cqueue(i_t2,l,pooltype,i) = 0.
      ENDIF
    ENDDO

    IF(ttpart(pooltype,i)>0)THEN
      l = ttstep(pooltype,i) + 1
      IF(riverstate%qqueue(l,pooltype,i)>0)THEN
        riverstate%cqueue(i_t2,l,pooltype,i) = meanrivertemp
      ELSE
        riverstate%cqueue(i_t2,l,pooltype,i) = 0.
      ENDIF
    ENDIF
    
  END SUBROUTINE set_rivertemp
  
  !>\brief Calculate and remove evaporation from lake
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Evaporation)
  !----------------------------------------------------------------------------------------
  SUBROUTINE calculate_actual_lake_evaporation(i,j,itype,numsubst,temp,epot,evap,cevap,lakestate)

    USE MODVAR, ONLY : basin,classdata, &
                       genpar,  &
                       landpar, &
                       cwater,  &
                       i_t1,i_t2
    USE HYPEVARIABLES, ONLY : m_ttmp,m_T1evap

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<subbasin index
    INTEGER, INTENT(IN) :: j        !<class index
    INTEGER, INTENT(IN) :: itype    !<lake type (ilake or olake)
    INTEGER, INTENT(IN) :: numsubst !<number of substances modelled
    REAL, INTENT(IN)    :: temp     !<air temperature
    REAL, INTENT(IN)    :: epot     !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: evap     !<actual evapotranspiration (mm/timestep)
    REAL, INTENT(OUT)   :: cevap(numsubst) !<concentration in evapotranspiration (eg. mg/L)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state

    !Local variables
    INTEGER k       !loop-variable, substance
    INTEGER status  !error status of subroutine
    REAL tt         !threshold temperature for evaporation (C)

    !Default values output variables
    evap = 0.
    cevap = 0.

    !Set local parameter
    tt = landpar(m_ttmp,classdata(j)%luse)       !Threshold temperature for snow melt and evaporation

    IF(temp>tt .AND. epot>0) THEN

      !Calculate actual evaporation, potential evaporation is default for temperature above threshold
      evap = epot
      DO k=1,numsubst
        cevap(k) = 0.
        IF(k==i_t1) cevap(k) = genpar(m_T1evap) * lakestate%conc(k,itype,i)
        !For T2==temperature, set the evaporation concentration = L/c = latent heat of vaporization divided by heat capacity
        IF(k==i_t2) cevap(k) = 1000. * latentheat_tempfunction(lakestate%conc(k,itype,i)) / cwater ! 1000 * [MJ/kg] / [KJ/kg/C] = C
      ENDDO

      !Remove evaporation from lake, check if enough water is available              
      IF(evap<lakestate%water(itype,i))THEN
        CALL remove_water(lakestate%water(itype,i),numsubst,lakestate%conc(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),basin(i)%subid,i,j)
      ELSEIF(numsubst==0)THEN
        evap = lakestate%water(itype,i)
        !cevap = lakestate%conc(:,itype,i)    !remove last traces of substances when lake dries out
        CALL remove_water(lakestate%water(itype,i),numsubst,lakestate%conc(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(2),basin(i)%subid,i,j)
      ELSEIF(lakestate%slowwater(itype,i)>0.)THEN
        !lake divided in two parts for N-, P- and C-simulations
        !mix lake parts and move to slowlake-part, evaporation from mixed lake volume
        CALL add_water(numsubst,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
        lakestate%water(itype,i)=0.
        IF(numsubst>0) lakestate%conc(:,itype,i)=0.
        IF(evap>=lakestate%slowwater(itype,i))THEN
          evap = lakestate%slowwater(itype,i)
          IF(numsubst>0) cevap = lakestate%concslow(:,itype,i)   !remove last traces of substances when lake dries out
        ENDIF
        CALL remove_water(lakestate%slowwater(itype,i),numsubst,lakestate%concslow(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(3),basin(i)%subid,i,j)
      ELSE
        !no water in slowlake-part, empty lake
        evap = lakestate%water(itype,i)
        IF(numsubst>0) cevap = lakestate%conc(:,itype,i)    !remove last traces of substances when lake dries out
        CALL remove_water(lakestate%water(itype,i),numsubst,lakestate%conc(:,itype,i),evap,cevap,status)
        IF(status.NE.0) CALL error_remove_water(errstring(3),basin(i)%subid,i,j)
      ENDIF
    ENDIF

  END SUBROUTINE calculate_actual_lake_evaporation

  !>Subroutine for summation of the area upstream of the outlet of all
  !>subbasins of the catchment
  !-------------------------------------------------------------------
  SUBROUTINE sum_upstream_area(n,areasum)
  
    USE MODVAR, ONLY : path,        &
                       basin,       &
                       branchdata,  &
                       branchindex

    !Argument declarations
    INTEGER, INTENT(IN)  :: n             !<number of subbasins
    REAL, INTENT(OUT)    :: areasum(n)    !<upstream area (m2)
    
    !Local variables
    INTEGER i,j,k,m         !loop variables
    REAL usarea             !summation variable for upstream area
    INTEGER, DIMENSION(n) :: A
    LOGICAL branchexists    !flag for branchdata available

    branchexists = .FALSE.
    IF(ALLOCATED(branchdata)) branchexists = .TRUE.
    A = 0
    areasum = 0.

    DO i = 1,n
      k = 0
      m = 0 
      USarea = 0. 
      DO j = 1,n
        IF(branchexists)THEN
          IF(branchindex(j)>0)THEN  !branch for this subbasin
            IF(i == path(j)%main)THEN
              m = m + A(j)
              USarea = USarea + areasum(j) * branchdata(branchindex(j))%mainpart
              k = k + 1      
            ENDIF
            IF(i == branchdata(branchindex(j))%branch)THEN
              m = m + a(j)
              usarea = usarea + areasum(j)*(1.-branchdata(branchindex(j))%mainpart)
              k = k + 1      
            ENDIF
          ELSE  !no branch for this subbasin
            IF(i == path(j)%main)THEN
              m = m + A(j)
              usarea = usarea + areasum(j)
              k = k + 1      
            ENDIF
          ENDIF
        ELSE    !no branches at all in the model set-up
           IF(i == path(j)%main)THEN
              m = m + a(j)
              usarea = usarea + areasum(j)
              k = k + 1      
           ENDIF
        ENDIF
      ENDDO
      IF(k==0) THEN                 !no inflows
        A(i) = 1
        areasum(i) = basin(i)%area
      ELSEIF(k==m) THEN             !k inflow, m (all) have their upstream area ready
        A(i) = 1
        areasum(i) = USarea + basin(i)%area
      ELSE                          !not all inflow have their upstream area ready (m<k)
        A(i) = 0                    !this indicates an error in coupling
        areasum(i) = 0.

        WRITE(6,*) 'ERROR in coupling of subbasins, some downstream basin before upstream basin'
        WRITE(6,*) 'i= ',i,' subid= ',basin(i)%subid
        WRITE(6,*) '(k= ',k,' m= ',m,')'
        STOP 1
      ENDIF
    ENDDO

  END SUBROUTINE sum_upstream_area

  !>Subroutine for calculation of general rating curve k-value for each lake
  !!
  !\b Reference ModelDescription Chapter Rivers and lakes (Lakes - Common lake processes)
  !-------------------------------------------------------------------
  SUBROUTINE set_general_rating_k(nl,n,locarea,areasum,rating)
    
    USE HYPEVARIABLES, ONLY : m_grat1,  &
                              m_grat3,   &
                              m_ratcorr, &
                              m_ilrrat1,m_olrrat1
    USE MODVAR, ONLY : genpar,  &
                       regpar,  &
                       basin, &
                       regiondivision

    !Argument declarations
    INTEGER, INTENT(IN)  :: nl            !<number of lake types
    INTEGER, INTENT(IN)  :: n             !<number of subbasins
    REAL, INTENT(IN)     :: locarea(n)    !<landarea of subbasin [m2]
    REAL, INTENT(IN)     :: areasum(n)    !<upstream area [m2]
    REAL, INTENT(OUT)    :: rating(nl,n)  !<k-value of rating curve
    
    !Local variables
    INTEGER i         !loop variables
    REAL ratcorr

    !Initial value
    rating = 0

    DO i = 1,n
      IF(basin(i)%parregion(regiondivision(m_ratcorr))>0)THEN
        ratcorr = 1. + regpar(m_ratcorr,basin(i)%parregion(regiondivision(m_ratcorr)))   !Correction of general rating curve grat1-parameter
      ELSE
        ratcorr = 1.
      ENDIF
      IF(basin(i)%parregion(regiondivision(m_ilrrat1))>0) rating(1,i) = regpar(m_ilrrat1,basin(i)%parregion(regiondivision(m_ilrrat1))) * ratcorr
      IF(basin(i)%parregion(regiondivision(m_olrrat1))>0) rating(2,i) = regpar(m_olrrat1,basin(i)%parregion(regiondivision(m_olrrat1))) * ratcorr
      IF(rating(1,i)<=0.) rating(1,i) = genpar(m_grat1) * ratcorr
      IF(rating(2,i)<=0.) rating(2,i) = genpar(m_grat1) * ratcorr
    ENDDO

    IF(genpar(m_grat3)>0.)THEN
      DO i = 1,n
        IF(locarea(i)*basin(i)%ilakecatch>0.) rating(1,i) = rating(1,i)*(locarea(i)*basin(i)%ilakecatch)**genpar(m_grat3)
        IF(areasum(i)>0.) rating(2,i) = rating(2,i)*(areasum(i))**genpar(m_grat3)
      ENDDO
    ENDIF
      
  END SUBROUTINE set_general_rating_k

  !>\brief Calculates temperature of rivers and lakes and other temperature variables
  !> rivertemp: 20-day moving average of air temperature
  !> laketemp: 5-day moving average of air temperature
  !> Also 10- and 20-day moving average of lake- and river temperature is calculated.
  !-----------------------------------------------------------------------
  SUBROUTINE calculate_water_temperature(i,airtemp,riverstate,lakestate)

    USE MODVAR, ONLY : genpar,        &
                       basin,         &
                       classbasin,    &
                       slc_ilake,     &
                       slc_olake,     &
                       i_in,i_sp,i_oc,i_ae,  &
                       timesteps_per_day
    USE HYPEVARIABLES, ONLY : m_laketemp

    !Argument declarations
    INTEGER, INTENT(IN) :: i                 !<index of current subbasin
    REAL, INTENT(IN)    :: airtemp           !<air temperature for subbasin
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
    
    !Local parameters
    REAL, PARAMETER :: rivertemp_days = 20.     !Number of days for river temperature calculation
    REAL, PARAMETER :: laketemp_days  = 5.      !Number of days for lake temperature calculation
    REAL, PARAMETER :: T10day_parameter = 10.
    REAL, PARAMETER :: T20day_parameter = 20.
    
    !Local variables
    INTEGER watertype               !Internal or main/outlet
    REAL    mtimesteps,mtimesteps2  !Number of timesteps temperature is averaged over

    !>\b Algorithm \n
    !>Calculate river temperature, same for local and main river
    IF(timesteps_per_day==1)THEN
      riverstate%temp(:,i) = riverstate%temp(:,i) + ((airtemp - riverstate%temp(:,i)) / rivertemp_days)
    ELSE
      mtimesteps = timesteps_per_day*rivertemp_days
      riverstate%temp(:,i) = riverstate%temp(:,i) + ((airtemp - riverstate%temp(:,i)) / mtimesteps)
    ENDIF

    !>Calculate lake temperature, same for internal and outlet lake (if exist)
    IF(genpar(m_laketemp)==0)THEN
      !>\li If parameter not set: as 5-day moving average
      IF(slc_ilake>0)THEN
        IF(classbasin(i,slc_ilake)%part>0)THEN
          IF(timesteps_per_day==1)THEN
            lakestate%temp(1,i) = lakestate%temp(1,i) + ((airtemp - lakestate%temp(1,i)) / laketemp_days)
          ELSE
            mtimesteps = timesteps_per_day*laketemp_days
            lakestate%temp(1,i) = lakestate%temp(1,i) + ((airtemp - lakestate%temp(1,i)) / mtimesteps)
          ENDIF
        ENDIF
      ENDIF
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0)THEN
          IF(timesteps_per_day==1)THEN
            lakestate%temp(2,i) = lakestate%temp(2,i) + ((airtemp - lakestate%temp(2,i)) / laketemp_days)
          ELSE
            mtimesteps = timesteps_per_day*laketemp_days
            lakestate%temp(2,i) = lakestate%temp(2,i) + ((airtemp - lakestate%temp(2,i)) / mtimesteps)
          ENDIF
        ENDIF
      ENDIF
    ELSE
      !>\li Elseif parameter set: as a moving average of a period determined by lake depth
      IF(slc_ilake>0)THEN
        IF(classbasin(i,slc_ilake)%part>0)THEN
          mtimesteps = timesteps_per_day*MIN(MAX(basin(i)%lakedepth(1),5.),5.+genpar(m_laketemp))
          lakestate%temp(1,i) = lakestate%temp(1,i) + ((airtemp - lakestate%temp(1,i)) / mtimesteps)
        ENDIF
      ENDIF
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0)THEN
          mtimesteps = timesteps_per_day*MIN(MAX(basin(i)%lakedepth(2),5.),5.+genpar(m_laketemp))
          lakestate%temp(2,i) = lakestate%temp(2,i) + ((airtemp - lakestate%temp(2,i)) / mtimesteps)
        ENDIF
      ENDIF
    ENDIF

    !>Calculate 10- and 20-day mean of water temperature for N,P,C or S processes
    IF(i_in>0 .OR. i_sp>0 .OR. i_oc>0 .OR. i_ae>0)THEN
      mtimesteps = timesteps_per_day*t10day_parameter
      mtimesteps2 = timesteps_per_day*t20day_parameter
      DO watertype = 1,2                   !(1=local/internal, 2=main/outlet)
        lakestate%temp10(watertype,i) = lakestate%temp10(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp10(watertype,i)) / mtimesteps)
        lakestate%temp20(watertype,i) = lakestate%temp20(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp20(watertype,i)) / mtimesteps2)
        riverstate%temp10(watertype,i) = riverstate%temp10(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp10(watertype,i)) / mtimesteps)
        riverstate%temp20(watertype,i) = riverstate%temp20(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp20(watertype,i)) / mtimesteps2)
      ENDDO
    ENDIF

  END SUBROUTINE calculate_water_temperature

  !>\brief Set temperature of rivers and lakes from T2 and calculate other temperature variables
  !Note: This is for today temperature
  !-----------------------------------------------------------------------
  SUBROUTINE set_water_temperature(waterbody,i,riverstate,lakestate)

    USE MODVAR, ONLY : classbasin,  &
                       slc_ilake,     &
                       slc_olake,     &
                       i_in,i_sp,i_oc,i_t2,  &
                       timesteps_per_day

    !Argument declarations
    INTEGER, INTENT(IN) :: waterbody         !<flag for waterbody, 1=lstream,2=ilake,3=main river,4=olake
    INTEGER, INTENT(IN) :: i                 !<index of current subbasin
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
    
    !Local parameters
    REAL, PARAMETER :: rivertemp_days = 20.     !Number of days for river temperature calculation
    REAL, PARAMETER :: laketemp_days  = 5.      !Number of days for lake temperature calculation
    REAL, PARAMETER :: T10day_parameter = 10.
    REAL, PARAMETER :: T20day_parameter = 20.
    
    !Local variables
    INTEGER watertype               !Internal or main/outlet
    REAL    mtimesteps,mtimesteps2  !Number of timesteps temperature is averaged over

    !>\b Algorithm \n
    !>Set river temperature to T2 temperature
    IF(waterbody==1) riverstate%temp(1,i) = riverstate%conc(i_t2,1,i)
    IF(waterbody==3) riverstate%temp(2,i) = riverstate%conc(i_t2,2,i)

    !>Set lake temperature (if exist)
    IF(waterbody==2)THEN
      IF(slc_ilake>0)THEN
        IF(classbasin(i,slc_ilake)%part>0)THEN
          lakestate%temp(1,i) = lakestate%conc(i_t2,1,i)
        ENDIF
      ENDIF
    ENDIF
    IF(waterbody==4)THEN
      IF(slc_olake>0)THEN
        IF(classbasin(i,slc_olake)%part>0)THEN
          lakestate%temp(2,i) = lakestate%conc(i_t2,2,i)
        ENDIF
      ENDIF
    ENDIF

    !>Calculate 10- and 20-day mean of water temperature for N,P or C processes
    IF(i_in>0 .OR. i_sp>0 .OR. i_oc>0)THEN
      mtimesteps = timesteps_per_day*t10day_parameter
      mtimesteps2 = timesteps_per_day*t20day_parameter
      watertype = 1                   !local/internal
      IF(waterbody==2) lakestate%temp10(watertype,i) = lakestate%temp10(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp10(watertype,i)) / mtimesteps)
      IF(waterbody==2) lakestate%temp20(watertype,i) = lakestate%temp20(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp20(watertype,i)) / mtimesteps2)
      IF(waterbody==1) riverstate%temp10(watertype,i) = riverstate%temp10(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp10(watertype,i)) / mtimesteps)
      IF(waterbody==1) riverstate%temp20(watertype,i) = riverstate%temp20(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp20(watertype,i)) / mtimesteps2)
      watertype = 2                   !main/outlet)
      IF(waterbody==4) lakestate%temp10(watertype,i) = lakestate%temp10(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp10(watertype,i)) / mtimesteps)
      IF(waterbody==4) lakestate%temp20(watertype,i) = lakestate%temp20(watertype,i) + ((lakestate%temp(watertype,i) - lakestate%temp20(watertype,i)) / mtimesteps2)
      IF(waterbody==3) riverstate%temp10(watertype,i) = riverstate%temp10(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp10(watertype,i)) / mtimesteps)
      IF(waterbody==3) riverstate%temp20(watertype,i) = riverstate%temp20(watertype,i) + ((riverstate%temp(watertype,i) - riverstate%temp20(watertype,i)) / mtimesteps2)
    ENDIF

  END SUBROUTINE set_water_temperature

  !>\brief Calculate river characteristics
  !>River characteristics include depth, area, bankful flow and 365-day-average-Q
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_river_characteristics(i,itype,flow,calcNPT,riverstate,depth,riverarea,qbank)
  
    USE MODVAR, ONLY : basin, &
                       regpar, &
                       genpar, &
                       regiondivision
    USE HYPEVARIABLES, ONLY : riverlength,  &
                              deadwidth,  &
                              m_velpar1,  &
                              m_velpar2,  &
                              m_velpar3,  &
                              m_widpar1,  &
                              m_widpar2,  &
                              m_widpar3,  &
                              m_maxwidth

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index of current subbasin
    INTEGER, INTENT(IN) :: itype      !<lake type (ilake or olake)
    REAL, INTENT(IN)    :: flow       !<river flow (m3/s) 
    LOGICAL, INTENT(IN) :: calcNPT    !<status of N,P,T1 simulation (to calculate bankful flow)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    REAL, INTENT(OUT)   :: depth      !<river depth (m)
    REAL, INTENT(OUT)   :: riverarea  !<river surface area (m2)
    REAL, INTENT(OUT)   :: qbank      !<flow at bank-ful river channel (m3/s)
    
    !Local variables
    REAL par(6)      !model parameters for velocity and width of river
    REAL rlength     !river length (m)
    REAL velocity    !river velocity (m/s)
    REAL width       !river width (m)

    !>\b Algorithm \n
    !>Set parameter values
    IF(basin(i)%parregion(regiondivision(m_velpar1))>0)THEN
      par = (/regpar(m_velpar1,basin(i)%parregion(regiondivision(m_velpar1))),regpar(m_velpar2,basin(i)%parregion(regiondivision(m_velpar2))),regpar(m_velpar3,basin(i)%parregion(regiondivision(m_velpar3))),regpar(m_widpar1,basin(i)%parregion(regiondivision(m_widpar1))),regpar(m_widpar2,basin(i)%parregion(regiondivision(m_widpar2))),regpar(m_widpar3,basin(i)%parregion(regiondivision(m_widpar3))) /)
    ELSE
      par = (/0,0,0,0,0,0/)   !OK? gives width=1 and depth=flow
    ENDIF

    !>Update state variable 365-day mean river discharge (m3/s)
    riverstate%Qmean(itype,i) = riverstate%Qmean(itype,i) + (flow-riverstate%Qmean(itype,i))/365.

    !>River length,depth and width, depend on velocity
    rlength = riverlength(itype,i)
    depth = 0.020          !low flow default value
    width = depth * 5.     !low flow default value
    IF(riverstate%Qmean(itype,i)>0.01.AND.flow>0.) THEN        
      velocity = (10**par(1)) * (riverstate%Qmean(itype,i)**par(2)) * ((flow/riverstate%Qmean(itype,i))**par(3))
      IF(velocity>0.2) THEN   
        width = (10**par(4)) * (flow/velocity)**(par(5)+par(6)*LOG10(flow/velocity))
        depth = (flow / velocity) / width             
      ENDIF
    ENDIF

    !>River (surface/bottom) area
    IF(genpar(m_maxwidth)>0)THEN
      riverarea = min(max(width,deadwidth(itype,i)),genpar(m_maxwidth)) * rlength
    ELSE  
      riverarea = max(width,deadwidth(itype,i)) * rlength
    ENDIF

    !TODO: separate into two subroutines: Qmean,depth,area resp Qbank. Depth,area only needed if a=0? Qmean more? depth more?
    !>Calculate new bankfull flow, stored in Q2max. 
    IF(calcNPT) CALL calc_qbank(flow,i,itype,riverstate%Q365(:,itype,i),riverstate%Qdayacc(:,itype,i),qbank)  !Subroutine also updates Qmax, Qdayacc and riverQ365.
    
  END SUBROUTINE calculate_river_characteristics

  !>Estimates the bank full flow by the second highest q from the
  !>daily values of last year
  !>
  !>\b Consequences Module hypevariables variables qmax, q2mqx, iqmax, and iq2max 
  !> may change.
  !>
  !>\b Reference ModelDescription Chapter Rivers and lakes (Rivers - Common river processes)
  !---------------------------------------------------------------------
  SUBROUTINE calc_qbank(flow,i,itype,riverq365,Qdayacc,Qbank)
    
    USE HYPEVARIABLES, ONLY : qmax,q2max,   &  !OUT
                              iqmax,iq2max     !OUT
    USE MODVAR, ONLY : dayno,     &
                       timesteps_per_day, &
                       endofday,  &
                       tsofday

    !Argument declarations
    REAL, INTENT(IN)     :: flow    !<flow current time step (m3/s)
    INTEGER, INTENT(IN)  :: i       !<index of current subbasin
    INTEGER, INTENT(IN)  :: itype   !<river type 1=local, 2=main
    REAL,INTENT(INOUT)   :: riverq365(366)  !<river flow last 365 days (m3/s)
    REAL,INTENT(INOUT)   :: Qdayacc(timesteps_per_day)  !<river flow last day (m3/s)
    REAL, INTENT(OUT)    :: qbank   !<bankfull flow
    
    !local variables
    REAL q        !average flow for day (m3/s)

    !Accumulate flow values for calculation of daily mean
    Qdayacc(tsofday) = flow

    IF(endofday)THEN
      q = SUM(Qdayacc(:))/REAL(timesteps_per_day)
      riverq365(dayno) = q !First year: initial assignment, following years: overwrite oldest value
      !Estimate river bankful flow with second highest flow
      IF(dayno==iqmax(itype,i) .OR. dayno==iq2max(itype,i))THEN !too old values, search whole array for new
        CALL update_qbank(riverq365(:),qmax(itype,i),q2max(itype,i),iqmax(itype,i),iq2max(itype,i))
      ELSEIF(q > qmax(itype,i))THEN
        q2max(itype,i) = qmax(itype,i)     !new highest flow
        iq2max(itype,i) = iqmax(itype,i)
        qmax(itype,i) = q
        iqmax(itype,i) = dayno
      ELSEIF(q > q2max(itype,i))THEN    !new second highest flow
        q2max(itype,i) = q
        iq2max(itype,i) = dayno
      ENDIF
    ENDIF
    qbank = q2max(itype,i)

  END SUBROUTINE calc_qbank

  !>Update highest and second highest flow when one of them reach
  !>retiring age
  !---------------------------------------------------------------------
  SUBROUTINE update_qbank(q_array,qmax,q2,imax,i2)

    !Argument declarations
    REAL, INTENT(IN)     :: q_array(366)  !<flow all days last year
    REAL, INTENT(OUT)    :: qmax          !<highest flow all days last year
    REAL, INTENT(OUT)    :: q2            !<second highest flow all days last year
    INTEGER, INTENT(OUT) :: imax          !<index of highest flow all days last year
    INTEGER, INTENT(OUT) :: i2            !<index of second highest flow all days last year
    
    !Local variables
    INTEGER i

    qmax = 0.
    q2 = 0. 

    DO i = 1, 366
      IF(q_array(i) >= qmax)THEN 
        q2 = qmax
        i2 = imax
        qmax = q_array(i)
        imax = i
      ELSEIF(q_array(i) > q2)THEN
        q2 = q_array(i)
        i2 = i
      ENDIF
    ENDDO

  END SUBROUTINE update_qbank

  !>Translation (delay) in river       
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Rivers - Common river processes)
  !-------------------------------------------------------------------
  SUBROUTINE translation_in_river(i,itype,qin,cin,qout,cout,riverstate)
  
    USE MODVAR, ONLY : numsubstances,     &
                       realzero,  &
                       seconds_per_timestep
    USE HYPEVARIABLES, ONLY : transtime,  &
                              ttstep,     &
                              ttpart

    !Argument declaration
    INTEGER, INTENT(IN) :: i                    !<index of current subbasin
    INTEGER, INTENT(IN) :: itype                !<river type (local or main)
    REAL, INTENT(IN)    :: qin                  !<inflow to river train (m3/s)
    REAL, INTENT(IN)    :: cin(numsubstances)   !<concentration of inflow to river train
    REAL, INTENT(OUT)   :: qout                 !<outflow of river train (m3/s)
    REAL, INTENT(OUT)   :: cout(numsubstances)  !<concentration of outflow of river train
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    !Local variables
    INTEGER y     !translation, whole time steps
    REAL    x     !translation, additional part of time step

    !>\b Algoritm \n
    !>Add new inflow to translation variable (river train)
    riverstate%qqueue(0,itype,i) = qin * seconds_per_timestep
    IF(numsubstances>0) riverstate%cqueue(:,0,itype,i) = cin

    IF(transtime(itype,i)>0)THEN
      !>Calculate outflow from river train
      y = ttstep(itype,i)
      x = ttpart(itype,i)
      qout = (1.-x)*riverstate%qqueue(y,itype,i) + x*riverstate%qqueue(y+1,itype,i) !Calculate flow (m3) from river after translation
      IF(qout>realzero .AND. numsubstances>0)THEN
        cout = ((1.-x)*riverstate%qqueue(y,itype,i)*riverstate%cqueue(:,y,itype,i) + &
                 x*riverstate%qqueue(y+1,itype,i)*riverstate%cqueue(:,y+1,itype,i))/qout
      ELSE
        cout = 0.
      ENDIF
      qout = qout / seconds_per_timestep  !flow (m3/s)

      !>Translate the flows in the river train
      riverstate%qqueue(1:y+1,itype,i) = riverstate%qqueue(0:y,itype,i)
      IF(numsubstances>0) riverstate%cqueue(:,1:y+1,itype,i) = riverstate%cqueue(:,0:y,itype,i)
    ELSE
      !Elseif no delay, outflow = inflow
      qout = qin
      IF(numsubstances>0) cout = cin
    ENDIF

  END SUBROUTINE translation_in_river


  !>\brief Abstraction of water from main river inflow and river
  !>Abstraction is taken from main river (before adding current inflow) and 
  !>local and upstream river inflow
  !>
  !> \b Reference ModelDescription Chapter Water management (Point sources - Negative point source)
  !-------------------------------------------------------------------
  SUBROUTINE point_abstraction_from_main_river_inflow(i,pooltype,q,riverstate,removedflow)
  
    USE MODVAR, ONLY : basin, load,  &
                       seconds_per_timestep
    USE HYPEVARIABLES, ONLY : ttstep,     &
                              ttpart

    !Argument declaration
    INTEGER, INTENT(IN) :: i                          !<index of current subbasin
    INTEGER, INTENT(IN) :: pooltype                   !<river type (local or main)
    REAL, INTENT(INOUT) :: q                          !<inflow (m3/s)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<river states
    REAL, INTENT(OUT)   :: removedflow                !<removed flow (m3/timestep)

    !Local variables
    INTEGER l     
    REAL    totvol      !volume in river (m3)
    REAL    absvol      !abstraction volume to be removed (m3)
    REAL    waterfrac   !fraction of water in current pool

    !>\b Algoritm \n
    removedflow = 0.
    IF(load(i)%abstrind/=3) RETURN  !not main river inflow abstraction
    IF(load(i)%abstrvol==0) RETURN
    
    !>If abstraction of water: Calculate amount
    totvol = q*seconds_per_timestep + riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
    absvol = load(i)%abstrvol*seconds_per_timestep  !m3
    !>Remove abstraction water proportionally from river and queue
    IF(absvol<totvol)THEN
      IF(q>0.)THEN
        waterfrac = q/totvol  !*seconds_per_timestep
        q = q - waterfrac*absvol  !/seconds_per_timestep
      ENDIF
      IF(riverstate%water(pooltype,i)>0.)THEN
        waterfrac = riverstate%water(pooltype,i)/totvol
        riverstate%water(pooltype,i) = riverstate%water(pooltype,i) - waterfrac*absvol
      ENDIF
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0.)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
          riverstate%qqueue(l,pooltype,i) = riverstate%qqueue(l,pooltype,i) - waterfrac*absvol
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that remaining outflow will be correct
          riverstate%qqueue(l,pooltype,i) = riverstate%qqueue(l,pooltype,i) - waterfrac*absvol
        ENDIF
      ENDIF
      removedflow = absvol
    ELSE
      riverstate%water(pooltype,i) = 0.
      riverstate%qqueue(1:ttstep(pooltype,i)+1,pooltype,i) = 0.
      IF(absvol>totvol)THEN
        WRITE(6,*) 'Warning: Point source abstraction from river could not be fulfilled, not enough water in river.'
        WRITE(6,*) 'Warning: subbasin ',basin(i)%subid, 'abstracted volume: ',totvol, 'of wished volume: ',absvol
      ENDIF
      removedflow = totvol
    ENDIF

  END SUBROUTINE point_abstraction_from_main_river_inflow

  !>\brief Abstraction of water from main river
  !>Abstraction is taken from main river volume and queue after inflow from local and upstream rivers
  !>
  !> \b Reference ModelDescription Chapter Water management (Point sources - Negative point source)
  !-------------------------------------------------------------------
  SUBROUTINE point_abstraction_from_main_river(i,pooltype,riverstate,removedflow)
  
    USE MODVAR, ONLY : basin, load,  &
                       seconds_per_timestep
    USE HYPEVARIABLES, ONLY : ttstep,     &
                              ttpart

    !Argument declaration
    INTEGER, INTENT(IN) :: i                          !<index of current subbasin
    INTEGER, INTENT(IN) :: pooltype                   !<river type (local or main)
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<river states
    REAL, INTENT(OUT)   :: removedflow                !<removed flow (m3/timestep)

    !Local variables
    INTEGER l     
    REAL    totvol      !volume in river (m3)
    REAL    absvol      !abstraction volume to be removed (m3)
    REAL    waterfrac   !fraction of water in current pool

    !>\b Algoritm \n
    removedflow = 0.
    IF(load(i)%abstrind/=1) RETURN  !not main river abstraction
    IF(load(i)%abstrvol==0) RETURN
    
    !>If abstraction of water: Calculate amount
    totvol = riverstate%water(pooltype,i) + (SUM(riverstate%qqueue(1:ttstep(pooltype,i),pooltype,i)) + riverstate%qqueue(ttstep(pooltype,i)+1,pooltype,i) * ttpart(pooltype,i))
    absvol = load(i)%abstrvol*seconds_per_timestep  !m3
    !>Remove abstraction water proportionally from river and queue
    IF(absvol<totvol)THEN
      IF(riverstate%water(pooltype,i)>0.)THEN
        waterfrac = riverstate%water(pooltype,i)/totvol
        riverstate%water(pooltype,i) = riverstate%water(pooltype,i) - waterfrac*absvol
      ENDIF
      DO l = 1,ttstep(pooltype,i)
        IF(riverstate%qqueue(l,pooltype,i)>0.)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol
          riverstate%qqueue(l,pooltype,i) = riverstate%qqueue(l,pooltype,i) - waterfrac*absvol
        ENDIF
      ENDDO
      IF(ttpart(pooltype,i)>0)THEN
        l = ttstep(pooltype,i) + 1
        IF(riverstate%qqueue(l,pooltype,i)>0)THEN
          waterfrac = riverstate%qqueue(l,pooltype,i)/totvol    !Note whole volume so that remaining outflow will be correct
          riverstate%qqueue(l,pooltype,i) = riverstate%qqueue(l,pooltype,i) - waterfrac*absvol
        ENDIF
      ENDIF
      removedflow = absvol
    ELSE
      riverstate%water(pooltype,i) = 0.
      riverstate%qqueue(1:ttstep(pooltype,i)+1,pooltype,i) = 0.
      IF(absvol>totvol)THEN
        WRITE(6,*) 'Warning: Point source abstraction from river could not be fulfilled, not enough water in river.'
        WRITE(6,*) 'Warning: subbasin ',basin(i)%subid, 'abstracted volume: ',totvol, 'of wished volume: ',absvol
      ENDIF
      removedflow = totvol
    ENDIF

  END SUBROUTINE point_abstraction_from_main_river

  !>Abstraction of water from outlet lake
  !>
  !> \b Reference ModelDescription Chapter Water management (Point sources - Negative point source)
  !-------------------------------------------------------------------
  SUBROUTINE point_abstraction_from_outlet_lake(i,pooltype,qunitfactor,lakestate,removedflow)
  
    USE MODVAR, ONLY : load,  &
                       seconds_per_timestep

    !Argument declaration
    INTEGER, INTENT(IN) :: i                    !<index of current subbasin
    INTEGER, INTENT(IN) :: pooltype             !<lake type (local or outlet)
    REAL, INTENT(IN)    :: qunitfactor          !<transformation factor m3/s->mm/timestep
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake states
    REAL, INTENT(OUT)   :: removedflow          !<removed flow (m3/timestep)

    !Local variables
    REAL    lakevol     !volume in lake (mm)
    REAL    absvol      !abstraction volume to be removed (mm)
    REAL    waterfrac   !fraction of water in current pool

    !>\b Algoritm \n
    removedflow = 0.
    IF(load(i)%abstrind/=2) RETURN  !not lake abstraction
    IF(load(i)%abstrvol==0) RETURN

    !>If abstraction of water: Calculate amount
    lakevol = lakestate%water(pooltype,i)       !mm
    IF(ALLOCATED(lakestate%slowwater)) lakevol = lakevol + lakestate%slowwater(pooltype,i)       !mm
    absvol = load(i)%abstrvol * qunitfactor  !mm/ts
    
    !>Remove abstraction water proportionally from fast and slow lake part
    IF(absvol<lakevol)THEN
      IF(ALLOCATED(lakestate%slowwater))THEN
        waterfrac = lakestate%water(pooltype,i)/lakevol
        lakestate%water(pooltype,i) = lakestate%water(pooltype,i) - absvol*waterfrac
        lakestate%slowwater(pooltype,i) = lakestate%slowwater(pooltype,i) - absvol*(1.-waterfrac)
      ELSE
        lakestate%water(pooltype,i) = lakestate%water(pooltype,i) - absvol
      ENDIF
      removedflow = load(i)%abstrvol * seconds_per_timestep
    ELSE
      lakestate%water(pooltype,i) = 0.
      IF(ALLOCATED(lakestate%slowwater)) lakestate%slowwater(pooltype,i) = 0.
      IF(absvol>lakevol) WRITE(6,*) 'Warning: Wanted abstraction from lake could not be fulfilled, not enough water in lake.'
      removedflow = lakevol / qunitfactor * seconds_per_timestep
    ENDIF

  END SUBROUTINE point_abstraction_from_outlet_lake

  !>\brief Subroutine for finding current lake outflow parameters. 
  !------------------------------------------------------------------------------
  SUBROUTINE get_current_lake_outflow_parameters(i,itype,lakeareain,olakewst,   &
                                         upstreamlakebasin,have2outlets,ratck,ratcexp, &
                                         w0Today,wmin,damProd,maxProd,minProd,out2ratck,out2ratcexp,  &
                                         out2w0Today,out2wmin,out2damProd,out2maxProd,out2minProd,qin)
       
    USE HYPEVARIABLES, ONLY : ratingk,      &  
                              m_grat2,      &  
                              m_ilrrat2, &
                              m_krelflood,  &
                              m_kthrflood,  &
                              m_klowflood
    USE MODVAR, ONLY : missing_value, &
                       lake, &
                       dam, &
                       basin, &
                       lakebasin, &
                       genpar, &
                       regpar, &
                       regiondivision, &
                       lakeindex, &
                       lakeout2index, &
                       damindex, &
                       lakebasinindex

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: itype         !<lake type (local or main)
    REAL, INTENT(IN)    :: lakeareain    !<lakearea (m2) (from GeoData)
    REAL, INTENT(IN)    :: olakewst      !<outlet lake water stage (m)
    LOGICAL, INTENT(OUT):: upstreamlakebasin   !<Upstream lake basin?
    LOGICAL, INTENT(OUT):: have2outlets   !<Lake with two outlets?
    REAL, INTENT(OUT)   :: ratck        !<current rating curve parameter rate
    REAL, INTENT(OUT)   :: ratcexp      !<current rating curve parameter exponent
    REAL, INTENT(OUT)   :: w0Today      !<current water level threshold in w-reference system (m)
    REAL, INTENT(OUT)   :: wmin         !<minimum water level threshold (sänkningsgräns) (m) 
    REAL, INTENT(OUT)   :: damProd      !<current dam production flow (m3/s)
    REAL, INTENT(OUT)   :: maxProd      !<(current) maximum production flow (m3/s)
    REAL, INTENT(OUT)   :: minProd      !<(current) minimum flow (m3/s)
    REAL, INTENT(OUT)   :: out2ratck    !<current rating curve parameter rate of outlet 2
    REAL, INTENT(OUT)   :: out2ratcexp  !<current rating curve parameter exponent of outlet 2
    REAL, INTENT(OUT)   :: out2w0Today  !<current water level threshold in w-reference system (m) of outlet 2
    REAL, INTENT(OUT)   :: out2wmin     !<minimum water level threshold (sänkningsgräns) (m) of outlet 2
    REAL, INTENT(OUT)   :: out2damProd  !<current dam production flow (m3/s) of outlet 2
    REAL, INTENT(OUT)   :: out2maxProd  !<(current) maximum production flow (m3/s) of outlet 2
    REAL, INTENT(OUT)   :: out2minProd  !<(current) minimum flow (m3/s) of outlet 2
    REAL,OPTIONAL,INTENT(IN) :: qin     !<current inflow to lake (m3/s)
    
    !Local variables
    REAL wlmr                   !water level lake (m)
    REAL qamp,qpha              !parameters for regulation of lake
    REAL qprodToday   !Production flow (m3/s)
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 
    REAL rating2        !general rating curve parameters outlet lake
    REAL out2w0rel  !-"-
    INTEGER dampurpose          ! Purpose of dam, 1=irrigation, 2=water supply, 3=flood control, 4=hydropower
    INTEGER current_lake  !Index of current lake in lake
    REAL qinfmax                ! Max mean monthly inflow
    REAL snowfrac               ! Fraction of precipitaiton falling as snow in catchment upstream of dam
    REAL regvol 
    REAL qinftoday              ! Current inflow for 'today'
    REAL qthresh                ! Threshold inflow over which a flood control dam save water
    REAL lthresh                ! Threshold reservoir level over which flood control dam releases extra flow
    REAL qinfmin                ! Min mean monthly inflow
    LOGICAL minflow,out2minflow ! Flag for minimum flow
    
    !Initial values
    wlmr=0.
    qprodToday = 0.
    qamp = 0.
    qpha = 0.
    regvol=0.
    qinfmax = 0.
    qinfmin=0.
    snowfrac = 0.
    dampurpose = 0
    lakearea = lakeareain    
    qinftoday=0.  
    lthresh = 0.  
    qthresh = 0.  
    current_lake = 0
    minflow = .FALSE.
    out2w0rel = 0.
    out2minflow = .FALSE.

    !Default output
    damProd = 0. 
    maxProd = 0. 
    minProd = 0. 
    ratck = 0.
    ratcexp = 0.
    wmin = missing_value
    out2damProd = 0. 
    out2maxProd = 0. 
    out2minProd = 0. 
    out2ratck = 0.
    out2ratcexp = 0.
    out2wmin = missing_value
    upstreamlakebasin = .FALSE.
    have2outlets = .FALSE.

    !Local lake parameter values
    IF(itype==1)THEN
      IF(basin(i)%parregion(regiondivision(m_ilrrat2))>0) rating2=regpar(m_ilrrat2,basin(i)%parregion(regiondivision(m_ilrrat2)))  !TODO: check subroutine for only olake, probably not used
      IF(rating2<=0.) rating2 = genpar(m_grat2)
      w0Today = basin(i)%lakedepth(1)      !ilake depth = threshold (m)
      ratck   = ratingk(itype,i)
      ratcexp = rating2
      RETURN
    ENDIF
    
    !Outlet lake or dam
    
    IF(PRESENT(qin))THEN
      qinftoday=qin
    ENDIF
     
    !Current lake parameter values
      IF(ALLOCATED(lakeindex))THEN
        IF(lakeindex(i)>0)THEN
          current_lake = lakeindex(i)
          maxProd = lake(lakeindex(i))%mqprod         !Maximum production flow
          minflow = lake(lakeindex(i))%minflow        !Minimum flow
          wmin = lake(lakeindex(i))%wmin              !lake threshold/"sänkningsgräns"
        ENDIF
      ENDIF
      IF(ALLOCATED(lakeout2index))THEN  !Second outlet of lake/dam
        IF(lakeout2index(i)>0)THEN
          have2outlets = .TRUE.
          out2maxProd = lake(lakeout2index(i))%mqprod         !Maximum production flow
          out2minflow = lake(lakeout2index(i))%minflow        !Minimum flow
          out2wmin = lake(lakeout2index(i))%wmin              !lake threshold/"sänkningsgräns"
          out2w0rel = lake(lakeout2index(i))%w0ref            !w0 of outlet 2 relative to outlet 1 (w0ref)
        ENDIF
      ENDIF
      IF(ALLOCATED(lakebasinindex))THEN
        IF(lakebasinindex(i)>0)THEN  
          IF(lakebasin(lakebasinindex(i))%last)THEN          !Recalculate water stage for last lake basin
            current_lake = lakebasin(lakebasinindex(i))%ilk
            wmin = lake(lakebasin(lakebasinindex(i))%ilk)%wmin              !lake threshold/"sänkningsgräns"
          ELSE
            upstreamlakebasin = .TRUE.
            current_lake = lakebasin(lakebasinindex(i))%ilk
            wmin = lake(lakebasin(lakebasinindex(i))%ilk)%wmin              !lake threshold/"sänkningsgräns"
          ENDIF
        ENDIF
      ENDIF
      IF(ALLOCATED(damindex))THEN
        IF(damindex(i)>0)THEN
          regvol = dam(damindex(i))%regvol          !Regvol
          snowfrac = dam(damindex(i))%snowfrac      ! Fraction of prec falling as snow upstream of dam
          qamp = dam(damindex(i))%qamp              !Amplitude of sin-adjustment of qprod
          qpha = dam(damindex(i))%qpha              !Phase of sin-adjustment of qprod
          wmin = dam(damindex(i))%wmin              !lake threshold/"sänkningsgräns"
          qinfmin = dam(damindex(i))%qinfmin
          qinfmax = dam(damindex(i))%qinfmax
          dampurpose = dam(damindex(i))%purpose              !dam purpose
          lthresh = 0.-genpar(m_klowflood)*regvol*1000000./lakearea       !threshold level for extra flood control releases (typical 1/3 of regvol)
          qthresh = genpar(m_kthrflood)*qinfmax 
        ENDIF
      ENDIF

      !Water level for outlet lake
      wlmr = olakewst
     
      !Dam in DamData
      IF(damindex(i)>0)THEN    !The below code is repeated first for dams, then for lakes
              
        !Current production flow for dam with regulation volume (calibrated dams)
        IF(wmin.NE.missing_value)THEN       ! i.e. RegVol>0
          CALL get_current_production_flow(0,damindex(i),wlmr,qprodToday)

          !Calculate dam outflow depending on dam purpose and current production flow
          IF(dampurpose==4)THEN    !hydroelectric dam
            damProd=qprodToday              ! If not a snow/seasonal redist dam and qamp not given, damProd=constant (Qinf or Qprod1)
            IF(qamp>0)THEN
              damProd = apply_seasonal_factor_on_production_flow(0,damindex(i),qprodToday)
            ELSEIF(snowfrac > 0.35)THEN     !CD2014 This determines if dams seasonally change flow by comparing regulation capacity to inflows (i.e. for dams wtih more than 35 % of precip that is snow
              IF(qpha>0) THEN
                qamp = 0.71                 !CD2014 based on regression from data, used if Qamp not given.                               
                damProd = apply_seasonal_factor_on_production_flow(0,damindex(i),qprodToday,qamp)
              ENDIF
            ENDIF
          ELSEIF(dampurpose==1)THEN   !irrigation dam
            damProd=qprodToday
          ELSEIF(dampurpose==2)THEN   !water supply dam
            damProd=qprodToday
          ELSEIF(dampurpose==3)THEN   !flood control dam (aim is to maintain dam as empty as possible)
            IF(qinftoday < qthresh)THEN        ! IF inflow today < threshold inflow
              IF(wlmr<lthresh)THEN             ! IF water level today < threshold level
                damProd=Qinftoday              ! Release the inflow
              ELSE      
                IF(qinftoday < qinfmin)THEN
                  damProd=MIN(Qinftoday*genpar(m_krelflood),qthresh)     ! If water level above threshold, release more than inflows (i.e. try empty the dam)
                ENDIF
              ENDIF
            ELSE                                 ! If inflow today >= threshold inflow
              damProd=qthresh                    ! Release maximum allowable flow
            ENDIF
          ELSE
            damProd=qprodToday                                                         
          ENDIF
          damProd=MIN(damProd,regvol*1000000./86400.)          ! Test: Limit damProd to the Regvol for one day
        ENDIF   !Close wmin.ne.missing

        !Ordinary outflow threshold
        w0Today = 0.
        
        !Set rating curve parameters
        IF(wlmr>w0Today)THEN
          CALL get_current_rating_parameters(i,0,damindex(i),ratck,ratcexp)
        ENDIF     

      !Lake/Dam in LakeData or GeoData
      ELSE
   
        !Current production flow for dam with regulation volume
        IF(wmin.NE.missing_value)THEN
          CALL get_current_production_flow(current_lake,0,wlmr,qprodToday)
          damProd = apply_seasonal_factor_on_production_flow(current_lake,0,qprodToday)
          IF(minflow) minProd = damProd     !Set minimum flow to current production flow
        ENDIF
 
        !Outflow threshold
        w0Today = 0.                  !Threshold w0 applies in general
        CALL adjust_threshold_for_seasonal_variation(current_lake,w0Today)
        IF(upstreamlakebasin)THEN  !Upstream lake basin; all water above threshold transported to last lakebasin
          IF(wmin==missing_value)THEN
            w0Today = w0Today     
          ELSE
            w0Today = wmin      
          ENDIF
        ENDIF
        
        IF(.NOT.upstreamlakebasin .AND. wlmr>w0Today)THEN
          CALL get_current_rating_parameters(i,current_lake,0,ratck,ratcexp)
        ENDIF
        
      ENDIF   ! ENDIF for lakes (i.e. damindex =/0)
      
      IF(have2outlets)THEN    !Calculate the second outlet
        !Current production flow for dam with regulation volume
        IF(out2wmin.NE.missing_value)THEN
          CALL get_current_production_flow(lakeout2index(i),0,wlmr,qprodToday)
          out2damProd = apply_seasonal_factor_on_production_flow(lakeout2index(i),0,qprodToday)
          IF(out2minflow) out2minProd = out2damProd     !Set minimum flow to current production flow
        ENDIF
 
        !Outflow threshold
        out2w0Today = 0. + out2w0rel        !Ordinary outflow threshold
        CALL adjust_threshold_for_seasonal_variation(lakeout2index(i),out2w0Today)
 
        IF(wlmr>out2w0Today)THEN
          CALL get_current_rating_parameters(i,lakeout2index(i),0,out2ratck,out2ratcexp)
        ENDIF
        
      ENDIF

  END SUBROUTINE get_current_lake_outflow_parameters
 
  !>\brief Subroutine for finding current production flow
  !------------------------------------------------------------------------------
  SUBROUTINE get_current_production_flow(current_lake,current_dam,wlmr,prodflow)
       
    USE HYPEVARIABLES, ONLY : m_limprod

    USE MODVAR, ONLY : missing_value, &
                       dayno, &
                       lake, &
                       dam, &
                       genpar

    !Argument declarations
    INTEGER, INTENT(IN) :: current_lake !<index of lake for current lake
    INTEGER, INTENT(IN) :: current_dam  !<index in dam for current dam
    REAL, INTENT(IN)    :: wlmr         !<outlet lake water stage (m)
    REAL, INTENT(OUT)   :: prodflow     !<current production flow (m3/s)
    
    !Local variables
    REAL qprod0            !Production flow which is equal to mean inflow
    REAL qprod1, qprod2    !Production flow (m3/s) for production periods 1 and 2
    REAL fracLevel         !Actual reservoir situation, in fraction of (w0 - wmin)
    REAL wmin              !minimum water level threshold (sänkningsgräns) (m) 
    REAL fillDamThreshold  !Percentage of reservoir capacity bellow which economy regime starts
    INTEGER dayno1, dayno2 !Starting day nr. for production periods 1 and 2
    
    !Default output
    prodflow = 0. 
    IF(current_lake==0 .AND. current_dam==0) RETURN
   
    !Current lake parameter values
    IF(current_lake>0)THEN
      qprod1 = lake(current_lake)%qprod1
      qprod2 = lake(current_lake)%qprod2
      dayno1 = lake(current_lake)%datum1
      dayno2 = lake(current_lake)%datum2
      wmin = lake(current_lake)%wmin    
      fillDamThreshold = lake(current_lake)%limprod
    ENDIF
    IF(current_dam>0)THEN   !Dam has priority (LakeData can hold NP-parameters)
      qprod0 = dam(current_dam)%qinfmed
      qprod1 = dam(current_dam)%qprod1
      qprod2 = dam(current_dam)%qprod2
      dayno1 = dam(current_dam)%datum1
      dayno2 = dam(current_dam)%datum2
      wmin = dam(current_dam)%wmin    
      fillDamThreshold = dam(current_dam)%limprod
      IF(qprod1<=0.)THEN
        qprod1=qprod0
        qprod2=qprod0
      ENDIF
    ENDIF
    IF(fillDamThreshold<0.) fillDamThreshold=genpar(m_limprod) !in case missing value, use general value

    !Calculate current production flow
    prodflow = qprod1                                              !Production rate 1 applies in general
    IF(dayno1*dayno2 > 0) THEN                                       !If both dates for different production regimes are non-zero, ...
      IF(dayno < dayno1 .OR. dayno >= dayno2)  prodflow = qprod2   !... and today is not within the periode of production rate 1, then production rate 2 applies
    ENDIF
    IF(fillDamThreshold>0)THEN
      fracLevel = (wlmr - wmin)/(0. - wmin) 
      IF(fracLevel > 0. .AND. fracLevel < fillDamThreshold) THEN
        prodflow = fracLevel/fillDamThreshold * prodflow        !Economy regime, if reservoir fractional filling is lower than threshold
      ENDIF
    ENDIF

  END SUBROUTINE get_current_production_flow
 
  !>\brief Function for applying seasonal variation on production flow
  !------------------------------------------------------------------------------
  REAL FUNCTION apply_seasonal_factor_on_production_flow(current_lake,current_dam,prodflow,qampin)
       
    USE MODVAR, ONLY : dayno,pi, &
                       lake,dam, &
                       tsofday,timesteps_per_day

    !Argument declarations
    INTEGER, INTENT(IN) :: current_lake !<index of lake for current lake
    INTEGER, INTENT(IN) :: current_dam  !<index in dam for current dam
    REAL, INTENT(IN)    :: prodflow     !<current production flow (m3/s)
    REAL, INTENT(IN), OPTIONAL :: qampin !<amplitude of seasonal variation to be used instead of parameter value
    
    !Local variables
    REAL qamp    !Amplitude of sin-adjustment of qprod
    REAL qpha    !Phase of sin-adjustment of qprod
    REAL seasonalfactor
    
    !Initial values
    qamp = 0.
    qpha = 0.
    seasonalfactor = 1.
    apply_seasonal_factor_on_production_flow = 0. 
    
    !Current lake parameter values
    IF(current_lake>0)THEN
      qamp = lake(current_lake)%qamp
      qpha = lake(current_lake)%qpha
    ENDIF
    IF(current_dam>0)THEN   !Dam has priority (LakeData then hold NP-parameters)
      qamp = dam(current_dam)%qamp
      qpha = dam(current_dam)%qpha
    ENDIF
    IF(PRESENT(qampin)) qamp = qampin

    !Calculate seasonal factor
!    IF(qamp>0.) seasonalfactor = 1. + qamp * SIN(2.*pi*(dayno+qpha)/365.)
    IF(qamp>0.) seasonalfactor = 1. + qamp * SIN(2.*pi*(dayno-1+REAL(tsofday)/REAL(timesteps_per_day)+qpha)/365.)
    apply_seasonal_factor_on_production_flow = seasonalfactor * prodflow
   
  END FUNCTION apply_seasonal_factor_on_production_flow
 
  !>\brief Subroutine for calculation current threshold 
  !------------------------------------------------------------------------------
  SUBROUTINE adjust_threshold_for_seasonal_variation(current_lake,w0)
       
    USE MODVAR, ONLY : dayno, &
                       lake

    !Argument declarations
    INTEGER, INTENT(IN) :: current_lake !<index of lake for current lake
    REAL, INTENT(INOUT) :: w0           !<current threshold (m)
    
    !Local variables
    INTEGER dayno1, dayno2      !Starting day nr. for production periods 1 and 2
    REAL deltaw0                !difference in water level threshold for period 2 (m)
    
    !Initial values
    deltaw0 = 0.
    
    !Current lake parameter values
    IF(current_lake>0)THEN
      dayno1 = lake(current_lake)%datum1
      dayno2 = lake(current_lake)%datum2
      deltaw0 = lake(current_lake)%deltaw0
    ENDIF

    !Calculate current threshold
    IF(deltaw0/=0.)THEN
      IF(dayno1*dayno2 > 0) THEN                             
        IF(dayno<dayno1 .OR. dayno>=dayno2)  w0 = w0 + deltaw0   !threshold for period 2 applies
      ENDIF
    ENDIF
   
  END SUBROUTINE adjust_threshold_for_seasonal_variation
 
  !>\brief Subroutine for finding current rating curve parameters for outlet lake
  !------------------------------------------------------------------------------
  SUBROUTINE get_current_rating_parameters(i,current_lake,current_dam,ratck,ratcexp)
       
    USE HYPEVARIABLES, ONLY : ratingk,      &  
                              m_grat2,      &  
                              m_olrrat2
    USE MODVAR, ONLY : missing_value, &
                       lake, &
                       dam, &
                       basin, &
                       regiondivision, &
                       genpar, &
                       regpar

    !Argument declarations
    INTEGER, INTENT(IN) :: i            !<index of current subbasin
    INTEGER, INTENT(IN) :: current_lake !<index of lake for current lake
    INTEGER, INTENT(IN) :: current_dam  !<index in dam for current dam
    REAL, INTENT(OUT)   :: ratck        !<current rating curve parameter rate
    REAL, INTENT(OUT)   :: ratcexp      !<current rating curve parameter exponent
    
    !Local variables
    REAL rating2        !general rating curve parameters outlet lake
    REAL regrate,regexp  !current parameters for specific rating curve
    REAL wmin
    
    INTEGER, PARAMETER :: itype = 2   !olake

    !Default output
    ratck = 0.
    ratcexp = 0.
    
    !Initial values current parameters
    regrate = 0.
    regexp = 0.
    wmin = missing_value

    !General rating curve exponent
    rating2 = 0.
    IF(basin(i)%parregion(regiondivision(m_olrrat2))>0) rating2=regpar(m_olrrat2,basin(i)%parregion(regiondivision(m_olrrat2))) 
    IF(rating2<=0.) rating2 = genpar(m_grat2)
     
    !Current lake or dam parameter values
    IF(current_lake>0)THEN
      regrate = lake(current_lake)%rate
      regexp = lake(current_lake)%exp
      wmin = lake(current_lake)%wmin
    ENDIF
    IF(current_dam>0)THEN   !Dam has priority (LakeData can hold NP-parameters)
      regrate = dam(current_dam)%rate
      regexp = dam(current_dam)%exp
      wmin = dam(current_dam)%wmin
    ENDIF
       
    !Set rating curve parameters
    IF(regrate>0.)THEN                 !Specific rating curve for lake or dam spill
      ratck = regrate
      ratcexp = regexp
    ELSEIF(wmin.NE.missing_value)THEN  !Dam without rating curve for spill
    ELSE                               !General rating curve for lake
      ratck = ratingk(itype,i)
      ratcexp = rating2
    ENDIF
        
  END SUBROUTINE get_current_rating_parameters
 
  !>\brief Subroutine for calculation and removal of outflow from local lake. 
  !!General rating curve is used for ilakes. 
  !----------------------------------------------------------------------------
  SUBROUTINE calculate_ilake_outflow(i,subid,ns,dividedlake,qin,lakearea,qunitfactor, &
                                     outflowm3s,coutflow,lakestate)

    USE HYPEVARIABLES, ONLY : ratingk, &  
                              m_grat2, &  
                              m_ilrrat2
    USE MODVAR, ONLY : basin, &
                       genpar, &
                       regpar, &
                       regiondivision

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: subid         !<subid of current subbasin
    INTEGER, INTENT(IN) :: ns            !<number of substances
    LOGICAL, INTENT(IN) :: dividedlake   !<status of lake divided in fast and slow part
    REAL, INTENT(IN)    :: qin           !<inflow of lake (m3/s) 
    REAL, INTENT(IN)    :: lakearea      !<lakearea (m2)
    REAL, INTENT(IN)    :: qunitfactor   !<factor for transforming flow for lake from m3/s to mm/timestep and back
    REAL, INTENT(OUT)   :: outflowm3s    !<outflow of lake (m3/s)
    REAL, INTENT(OUT)   :: coutflow(ns)  !<concentration of outflow of lake
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state

    !Local parameters
    INTEGER, PARAMETER :: itype = 1        !lake type (local)
    REAL, PARAMETER :: lakedepth_not_used = 1.
    REAL, PARAMETER :: epidepth_not_used = 1.

    !Local variables
    REAL lakewstmm              !lake water stage (mm)
    REAL wlmr,wlmr0             !water level lake  (m)
    REAL w0Today                !water level threshold  (m)
    REAL ratingc,ratinge        !current rating curve parameters
    REAL outflowmm              !outflow of lake (mm)

    !>\b Algorithm \n
    !Initial values
    outflowm3s = 0.
    lakewstmm=lakestate%water(itype,i)
    IF(dividedlake) lakewstmm=lakewstmm+lakestate%slowwater(itype,i)
    wlmr = lakewstmm*0.001              !Water in lake [m]
    
    !Current parameter values
    ratingc = ratingk(itype,i)
    ratinge = 0.
    IF(basin(i)%parregion(regiondivision(m_ilrrat2))>0) ratinge=regpar(m_ilrrat2,basin(i)%parregion(regiondivision(m_ilrrat2)))  !TODO: check subroutine for only olake, probably not used
    IF(ratinge<=0.) ratinge = genpar(m_grat2)
    w0Today = basin(i)%lakedepth(itype)
   
    !>Calculate outflow from general rating curve
    wlmr0 = wlmr - w0Today
    IF(wlmr0>0.)  outflowm3s = average_flow_rating_curve(qin,lakearea,wlmr0,ratingc,ratinge)  ![m3/s]
    outflowmm = outflowm3s * qunitfactor        ![mm/timestep]

    !>Check outflow against lake volume
    IF(outflowmm*0.001>wlmr0)THEN   !Check for enough water in lake (bad rating curve or numerical problems)
      IF(wlmr0>0.)THEN
        outflowmm = wlmr0*1000.
      ELSE
        outflowmm = 0.
      ENDIF
      IF(outflowmm>lakewstmm) outflowmm = lakewstmm   !Safety for rounded wlmr and ldepth = 0
      outflowm3s = outflowmm/qunitfactor
    ENDIF
    
    !>Remove outflow from lake
    CALL remove_outflow_from_lake(i,itype,ns,outflowmm,subid,lakedepth_not_used,epidepth_not_used,lakewstmm,coutflow,lakestate)

  END SUBROUTINE calculate_ilake_outflow

  !>\brief Subroutine for calculation outflow from outlet lake. 
  !>
  !!For outlet lakes several options exist: 
  !!Specific rating curve, general rating curve, all water above threshold for 
  !!upstream lake basin, regulation with spill by rating curve, 
  !!constant production flow depending on date or two separate rating curves for olake.
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_outflow_from_outlet_lake(i,qin,lakeareain,lakewstmm,   &
                                         qunitfactorin,outflowm3s,outflowmm, &
                                         outflow1,outflow2,maxQprodOUT,minFlowOUT,lakestate)
       
    USE HYPEVARIABLES, ONLY : lakeoutlet
    USE MODVAR, ONLY : missing_value, &
                       branchindex, &
                       basin,   &
                       lakebasin,  &
                       lakebasinindex

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: qin           !<inflow of lake (m3/s) 
    REAL, INTENT(IN)    :: lakeareain    !<lakearea (m2)
    REAL, INTENT(IN)    :: lakewstmm     !<lake water stage (mm)
    REAL, INTENT(IN)    :: qunitfactorin !<factor for transforming flow for lake from m3/s to mm/timestep and back
    REAL, INTENT(OUT)   :: outflowm3s    !<outflow of lake (m3/s)
    REAL, INTENT(OUT)   :: outflowmm     !<outflow of lake (mm)
    REAL, INTENT(OUT)   :: outflow1      !<outflow of first outlet (m3/s) or all outflow if not lake with 2 outlets in LD
    REAL, INTENT(OUT)   :: outflow2      !<outflow of second outlet (m3/s)
    REAL, INTENT(OUT)   :: maxQprodOUT   !<temporary maximum production for updated flow (m3/s)
    REAL, INTENT(OUT)   :: minFlowOUT    !<current minimum flow for updated flow (m3/s)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state

    !Local variables
    LOGICAL upstreamlakebasin   !Upstream lake basin
    LOGICAL have2outlets        !Lake with two outlets
    REAL wlmr,wlmr0             !water level lake  (m)
    REAL wmin                   !water levels threshold for production (m)
    REAL w0ref                  !water stage reference level (not used)
    REAL w0Today                !water level threshold  (m)
    REAL wthresh                !water level threshold for checking volume (m)
    REAL ratingc,ratinge        !current rating curve parameters outlet lake
    REAL damProd,maxQprod       !Current and maximum dam production flow (m3/s)
    REAL minflow,out2minflow    !Minimum flows (m3/s)
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 
    REAL qunitfactor            !factor for transforming flow for lake from m3/s to mm/timestep and back (adjusted for last lakebasin) 
    REAL out2ratingc,out2ratinge,out2w0Today,out2wmin,out2damProd,out2maxQprod !parameters for outlet 2
    REAL wcheck1, wcheck2       !water level thresholds for check of flow against volume
    REAL outflowmmnew,outflownew

    !Local parameters
    INTEGER, PARAMETER :: itype = 2  !lake type (outlet lake)

    !>\b Algorithm \n
    !> Set initial values
    outflowm3s = 0.
    outflowmm = 0.
    outflow1 = 0.
    outflow2 = 0.
    maxQprodOUT = 0.
    minFlowOUT = 0.
    wlmr = lakewstmm*0.001              !Water in lake (m) (default=absolute)
    lakearea = lakeareain      
    qunitfactor = qunitfactorin
    
    !>Calculate water level for outlet lake
    CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)

    !>Get current parameter values
    CALL get_current_lake_outflow_parameters(i,itype,lakeareain,wlmr,   &
                 upstreamlakebasin,have2outlets,ratingc,ratinge,w0Today,  &
                 wmin,damProd,maxQprod,minflow,out2ratingc,out2ratinge,out2w0Today,out2wmin, &
                 out2damProd,out2maxQprod,out2minflow,qin)
   
    !>Outflow determination (and check) for two outlet lake
    IF(have2outlets)THEN    !Outlet lake with two outlets
      CALL calculate_lake_outlet_outflow(lakeoutlet(branchindex(i))%otype(1), &
                                    qin,lakeareain,wlmr,ratingc,  &
                                    ratinge,w0Today,wmin,damProd,wcheck1,outflow1)
      CALL calculate_lake_outlet_outflow(lakeoutlet(branchindex(i))%otype(2), &
                                    qin,lakeareain,wlmr,out2ratingc,  &
                                    out2ratinge,out2w0Today,out2wmin,out2damProd,wcheck2,outflow2)
      IF(lakeoutlet(branchindex(i))%otype(1)==2 .OR. lakeoutlet(branchindex(i))%otype(1)==7) &  !Check against max production
        CALL calculate_maxprod_outflow(outflow1,outflow2,maxQprod,out2minflow)
      IF(lakeoutlet(branchindex(i))%otype(2)==2 .OR. lakeoutlet(branchindex(i))%otype(2)==7) &
        CALL calculate_maxprod_outflow(outflow2,outflow1,out2maxQprod,minflow)
      IF(lakeoutlet(branchindex(i))%change==4 .OR. &
          lakeoutlet(branchindex(i))%change==6) maxQprodOUT = damProd   !Set temporary maxProd for recalculating branched flow (step 1)
      IF(lakeoutlet(branchindex(i))%change==5 .OR. &
          lakeoutlet(branchindex(i))%change==7) maxQprodOUT = out2damProd
      IF(lakeoutlet(branchindex(i))%change==7 .OR. &
          lakeoutlet(branchindex(i))%change==9) minflowOUT = minflow   !Set current minflow for recalculating branched flow
      IF(lakeoutlet(branchindex(i))%change==6 .OR. &
          lakeoutlet(branchindex(i))%change==8) minflowOUT = out2minflow
      outflowm3s = outflow1 + outflow2
      outflowmm = outflowm3s * qunitfactor        !to mm/ts
      IF(outflowmm>(0.-MIN(wcheck1,wcheck2))*1.E3 .OR. outflowmm>lakewstmm)THEN !Check against available water
        outflowmmnew = MIN((0.-MIN(wcheck1,wcheck2))*1.E3,lakewstmm)
        outflownew = outflowmmnew / qunitfactor
        IF(outflownew/=outflowm3s)THEN
          IF(lakeoutlet(branchindex(i))%change==4 .OR. &
              lakeoutlet(branchindex(i))%change==6) maxQprod = maxQprodOUT   !Set temporary maxProd for recalculating branched flow (step 2)
          IF(lakeoutlet(branchindex(i))%change==5 .OR. &
              lakeoutlet(branchindex(i))%change==7) out2maxQprod = maxQprodOUT
          CALL recalculate_branched_flow(lakeoutlet(branchindex(i))%change, &
                      outflownew,maxQprod,out2maxQprod,minflow,out2minflow,outflow1,outflow2)
          outflowm3s = outflownew
          outflowmm = outflowmmnew
        ENDIF
      ENDIF
      RETURN
    ENDIF

    !>Outflow determination for one outlet lake
    wlmr0 = wlmr -w0Today
    IF(upstreamlakebasin)THEN     !Upstream lake basin; all water above threshold
      IF(wlmr0>0.) outflowm3s = wlmr0 * 1000. / qunitfactor
      outflow1 = outflowm3s
    ELSEIF(wmin==missing_value)THEN   !Not regulated lake
      IF(wlmr0>0.)THEN            !Rating curve used for water above threshold
        IF(ratingc>0.)THEN
          outflowm3s = MAX(average_flow_rating_curve(qin,lakearea,wlmr0,ratingc,ratinge), damProd) !damProd=0 here
        ELSE
          WRITE(6,*) 'Error: Ended in else that is not acceptable. Outflow of outlet lake'
          WRITE(6,*) 'Check input data for this lake.'
          WRITE(6,*) 'i',i,'itype',itype
          WRITE(6,*) 'More info: wlmr',wlmr,'w0Today',w0Today,'wmin',wmin,'ratingc',ratingc
        ENDIF
      ENDIF
      outflow1 = outflowm3s
    ELSE          !Regulated lake
      IF(wlmr0>0.)THEN
        IF(ratingc>0)THEN     !Specific rating curve for production or dam spill
          outflowm3s = MAX(average_flow_rating_curve(qin,lakearea,wlmr0,ratingc,ratinge), damProd)
        ELSE      !Dam without rating curve for spill; all water above threshold but at least production
          outflowm3s = MAX(wlmr0*1000./qunitfactor, damProd)
        ENDIF
      ELSEIF(wlmr>wmin)THEN             !Production flow to lower threshold
        outflowm3s = MIN((wlmr-wmin) * 1000. / qunitfactor, damProd)
      ENDIF
      outflow1 = outflowm3s
    ENDIF
    outflowmm = outflowm3s * qunitfactor        !to mm/ts
    
    !>Check outflow against lake volume (bad rating curve parameters or numerical problems)
    
    !>\li Calculate threshold for checking
    IF(wmin==missing_value .OR. upstreamlakebasin)THEN
      wthresh = w0Today
    ELSE  
      wthresh = wmin
    ENDIF

    !>\li Recalculate current water stage for last lakebasin of lakebasin lake, because this is the volume available
    IF(ALLOCATED(lakebasinindex))THEN
      IF(lakebasinindex(i)>0)THEN  
        IF(lakebasin(lakebasinindex(i))%last)THEN
          wlmr = lakewstmm*0.001 - basin(i)%lakedepth(2) 
        ENDIF
      ENDIF
    ENDIF
      
    !>\li Check against lowest water level allowed and lake volume
    IF(outflowmm*0.001>wlmr-wthresh)THEN
      IF(wlmr>wthresh)THEN
        outflowmm = (wlmr-wthresh)*1000.
      ELSE
        outflowmm = 0.
      ENDIF
      IF(outflowmm>lakewstmm) outflowmm = lakewstmm   !Safety for rounded wlmr used. 
      outflowm3s = outflowmm/qunitfactor
    ENDIF

  END SUBROUTINE calculate_outflow_from_outlet_lake

  !>\brief Subroutine for calculating current outflow of one lake outlet
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_lake_outlet_outflow(otype,qin,lakearea,wlmr,ratc,ratexp,w0Today,  &
                                      wmin,damProd,wcheck,outflow)

    !Argument declarations
    INTEGER, INTENT(IN) :: otype    !<outlet type as defined in lakeoutlet (1-7)
    REAL, INTENT(IN)    :: qin      !<inflow of lake (m3/s) 
    REAL, INTENT(IN)    :: lakearea !<lakearea (m2)
    REAL, INTENT(IN)    :: wlmr     !<lake water stage (m)
    REAL, INTENT(IN)    :: ratc     !<rating curve coefficient
    REAL, INTENT(IN)    :: ratexp   !<rating curve exponent
    REAL, INTENT(IN)    :: w0Today  !<upper threshold (m)
    REAL, INTENT(IN)    :: wmin     !<lower threshold (m)
    REAL, INTENT(IN)    :: damProd  !<current production flow (m3/s)
    REAL, INTENT(OUT)   :: wcheck   !<current outflow threshold (m)
    REAL, INTENT(OUT)   :: outflow  !<current outflow (m3/s)

    !Local variables

    outflow = 0.
    
    !>Calculate outflow for current outlet type
    IF(otype==1 .OR. otype==2 .OR. otype==8)THEN
      IF(wlmr>wmin)THEN             !Production flow to lower threshold
        outflow = damProd
        wcheck = wmin
      ENDIF
    ELSEIF(otype==3 .OR. otype==4 .OR. otype==6 .OR. otype==7)THEN
      IF(wlmr>w0Today)THEN          !Rating curve for water above threshold
        outflow = average_flow_rating_curve(qin,lakearea,wlmr-w0Today,ratc,ratexp)
        wcheck = w0Today
      ENDIF
    ELSEIF(otype==5 .OR. otype==9)THEN
      IF(wlmr>w0Today)THEN          !Production flow and rating curve above threshold
        outflow = average_flow_rating_curve(qin,lakearea,wlmr-w0Today,ratc,ratexp)
        outflow = MAX(outflow, damProd)
      ELSEIF(wlmr>wmin)THEN
        outflow = damProd
      ENDIF
      wcheck = wmin
    ENDIF

  END SUBROUTINE calculate_lake_outlet_outflow

  !>\brief Subroutine for increasing production flow to maximum before using
  !!overflow branch
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_maxprod_outflow(outflow1,outflow2,maxQprod,minflow2)

    !Argument declarations
    REAL, INTENT(INOUT) :: outflow1   !<outflow of first outlet (m3/s)
    REAL, INTENT(INOUT) :: outflow2   !<outflow of second outlet (m3/s)
    REAL, INTENT(IN)    :: maxQprod   !<maximum production (m3/s), second priority
    REAL, INTENT(IN)    :: minflow2   !<minimum flow (m3/s), first priority

    !Local variables
    REAL outflowm3s   !total flow

    !>Redistribute flow for maximum production
    outflowm3s = outflow1 + outflow2
    IF(outflowm3s<minflow2)THEN
      outflow1 = 0.
      outflow2 = outflowm3s
    ELSEIF(outflowm3s>maxQprod+minflow2)THEN
      outflow1 = maxQprod
      outflow2 = outflowm3s - outflow1
    ELSE
      outflow1 = outflowm3s - minflow2
      outflow2 = minflow2
    ENDIF

  END SUBROUTINE calculate_maxprod_outflow

  !>\brief Momentanous flow by rating curve
  !>
  !>Subroutine for calculation momentanous outflow from lake from current lake 
  !>water stage by simple lake rating curve equation. 
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_flow_from_outlet_lake_waterstage(i,lakeareain,lakewstmm,outflowm3s,lakestate)
       
    USE GENERAL_FUNCTIONS, ONLY : simple_rating_curve

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: lakeareain    !<lakearea (m2)
    REAL, INTENT(IN)    :: lakewstmm     !<lake water stage (mm)
    REAL, INTENT(OUT)   :: outflowm3s    !<outflow of lake (m3/s)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state

    !Local variables
    LOGICAL upstreamlakebasin   !Upstream lake basin; Tappa tillrinning
    LOGICAL have2outlets        !Lake with two outlets
    REAL wlmr                   !water level lake (m)
    REAL wmin                   !water level threshold in w-reference system (m)
    REAL w0Today                !water level threshold in w-reference system (m)
    REAL w0ref                  !water stage reference level (m) (not used)
    REAL ratingc,ratinge        !general rating curve parameters outlet lake
    REAL damProd                !Dam production flow
    REAL lakearea               !lakearea (m2) (adjusted for last lakebasin) 
    REAL out2ratingc,out2ratinge,out2w0Today,out2wmin,out2damProd !parameters for outlet 2
    REAL maxProd,out2maxProd    !parameters not used here
    REAL minflow,out2minflow    !parameters not used here

    !Local constant
    INTEGER, PARAMETER :: itype = 2   !lake type (outlet lake)

    !Initial values
    outflowm3s = 0.
    lakearea = lakeareain     
    wlmr = lakewstmm*0.001              !Water in lake (m) (default=absolute)
    CALL calculate_olake_waterstage(i,lakewstmm,lakeareain,lakearea,wlmr,lakestate,w0ref)
    
    !Current parameter values
    CALL get_current_lake_outflow_parameters(i,itype,lakeareain,wlmr,   &
                 upstreamlakebasin,have2outlets,ratingc,ratinge,w0Today,  &
                 wmin,damProd,maxProd,minflow,out2ratingc,out2ratinge,out2w0Today,out2wmin, &
                 out2damProd,out2maxProd,out2minflow)

    !Outflow determination
    IF(ratingc>0)THEN
      outflowm3s = simple_rating_curve(wlmr,ratingc,ratinge,w0Today)
    ELSE
      outflowm3s = 0. !Error in indata reaching this else?
    ENDIF

  END SUBROUTINE calculate_flow_from_outlet_lake_waterstage

  !>\brief Removal of outflow from lake and setting of
  !>concentration of outflow.
  !>
  !>\b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions)
  !-----------------------------------------------------------------------
  SUBROUTINE remove_outflow_from_lake(i,itype,ns,outflowmm,subid,ldepthm,epidepth,lakewstmm,coutflow,lakestate)

    USE HYPEVARIABLES, ONLY : m_lddeeplake, &
                              m_ldfastlake, &
                              m_gt2mix,     &
                              m_ldt2mix
    USE MODVAR, ONLY : lakedatapar,     &
                       lakedataparindex, &
                       genpar, &
                       i_t2

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: itype         !<lake type (local or main)
    INTEGER, INTENT(IN) :: ns            !<number of substances
    REAL, INTENT(IN)    :: outflowmm     !<outflow of lake (mm/timestep)
    INTEGER, INTENT(IN) :: subid         !<subid of current subbasin, for error output
    REAL, INTENT(IN)    :: ldepthm       !<lake depth (m)
    REAL, INTENT(IN)    :: epidepth      !<lake epilimnion depth (m)
    REAL, INTENT(IN)    :: lakewstmm     !<lake water stage (mm)
    REAL, INTENT(OUT)   :: coutflow(ns)  !<concentration of outflow of lake
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    
    !Local variables
    INTEGER status
    REAL fastoutpart  !fraction of outflow from fast lake part (if possible)
    REAL q1,q2        !outflow from lake parts (mm/timestep)
    REAL wstabovethr  !water stage above lake threshold/dam crest (m)
    REAL t2conc       !outflow T2 temperature
    REAL upperpart    !part of outflow that is taken from the upper part of the lake at T2 calculations
    REAL tempconc(ns) !help variable used in t2 calculations
    LOGICAL usestrattemp  !set water temp in outflow of lake to depend on lake temperature stratification

    !>\b Algorithm \n
    !>Preparations: default values
    coutflow = 0.
    usestrattemp = .FALSE.

    !>Preparations: T2-simulation with special outflow temperature of outlet lakes
    !Calculate the fraction of the outflow that will be taken from upper lake part
    IF(i_t2 > 0 .AND. itype == 2) THEN              
      usestrattemp = .TRUE.
      wstabovethr = lakewstmm * 0.001 - ldepthm   !waterstage above lake threshold
      IF(wstabovethr > 0) THEN
        upperpart = MAX(0.,MIN(1.,epidepth / wstabovethr))
      ELSE
        upperpart = 1.    !Dams with production flow will get all uppertemp!?
      ENDIF         
      !Decide if mean temperature is used as outflow t2 concentration    
      IF(genpar(m_gt2mix)==1)THEN
        usestrattemp = .FALSE.
      ELSEIF(lakedatapar(lakedataparindex(i,itype),m_ldt2mix)==1)THEN
        usestrattemp = .FALSE.
      ENDIF
      IF(usestrattemp) t2conc = upperpart * lakestate%uppertemp(itype,i) + (1-upperpart) * lakestate%lowertemp(itype,i)   
    ENDIF

    !>Remove outflow and set outflow concentrations:
    !>If lake is not divided:
    IF(ns==0)THEN
      IF(outflowmm>0)THEN
        !>\li Outflow is removed from lake
        !coutflow = lakestate%conc(:,itype,i)
        CALL remove_water(lakestate%water(itype,i),ns,lakestate%conc(:,itype,i),outflowmm,coutflow,status)
        IF(status.NE.0) CALL error_remove_water(errstring(4),subid,i,itype)
      ENDIF
    ELSE
      IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)==0)THEN
        !>If divided totally mixed lake (only slowlake):
        !>\li Add all water to slowlake and mix
        IF(lakestate%water(itype,i)>0)THEN   !add lakewater to slowwater and mix
          CALL add_water(ns,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
          lakestate%water(itype,i)=0.
          IF(ns>0) lakestate%conc(:,itype,i)=0.
        ENDIF
        !>\li Outflow is removed from slowlake
        IF(outflowmm>0)THEN
          IF(ns>0) coutflow(:) = lakestate%concslow(:,itype,i)
          IF(usestrattemp) coutflow(i_t2) = t2conc
          CALL remove_water(lakestate%slowwater(itype,i),ns,lakestate%concslow(:,itype,i),outflowmm,coutflow,status)
          IF(status.NE.0) CALL error_remove_water(errstring(5),subid,i,itype)
        ENDIF
      ELSE
        !>If lake is divided in fast and slow part and outflow occur from both parts:
        !>\li Calculate the outflow fraction from each part
        IF(outflowmm>0)THEN
          fastoutpart=lakedatapar(lakedataparindex(i,itype),m_ldfastlake)*lakestate%water(itype,i)/(lakestate%water(itype,i)+lakestate%slowwater(itype,i))
          q1 = fastoutpart*outflowmm
          q2 = (1.-fastoutpart)*outflowmm
          IF(q1>lakestate%water(itype,i))THEN
            q2 = q2 + (q1 - lakestate%water(itype,i))
            q1 = lakestate%water(itype,i)
          ENDIF
          IF(q2>lakestate%slowwater(itype,i))THEN
            q1 = q1 + (q2 - lakestate%slowwater(itype,i))
            q2 = lakestate%slowwater(itype,i)
          ENDIF
          !>\li Remove the outflow from both parts
          IF(q1>0.)THEN
            IF(ns>0) tempconc(:) = lakestate%conc(:,itype,i)
            IF(usestrattemp) tempconc(i_t2) = t2conc
            CALL remove_water(lakestate%water(itype,i),ns,lakestate%conc(:,itype,i),q1,tempconc(:),status)
            IF(status.NE.0) CALL error_remove_water(errstring(6),subid,i,itype)
          ENDIF
          IF(q2>0)THEN 
            IF(ns>0) tempconc(:) = lakestate%concslow(:,itype,i)
            IF(usestrattemp) tempconc(i_t2) = t2conc
            CALL remove_water(lakestate%slowwater(itype,i),ns,lakestate%concslow(:,itype,i),q2,tempconc(:),status)
            IF(status.NE.0) CALL error_remove_water(errstring(7),subid,i,itype)
          ENDIF
          !>\li Calculate the concentration of outflow
          IF(ns>0) coutflow(:) = (q1*lakestate%conc(:,itype,i)+q2*lakestate%concslow(:,itype,i))/outflowmm
          IF(usestrattemp) coutflow(i_t2) = t2conc
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE remove_outflow_from_lake

  !>\brief Flow between lake parts for divided lake
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions)
  !-------------------------------------------------------------------
  SUBROUTINE calculate_flow_within_lake(i,itype,subid,lakestate)

    USE HYPEVARIABLES, ONLY : slowlakeini,  &
                              m_lddeeplake
    USE MODVAR, ONLY : numsubstances,   &
                       lakedatapar, &
                       lakedataparindex
 
    !Argument declarations
    INTEGER, INTENT(IN) :: i     !<index of current subbasin
    INTEGER, INTENT(IN) :: itype !<lake type (local or main)
    INTEGER, INTENT(IN) :: subid !<subid of current subbasin, for error output
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    
    !Local variables
    INTEGER status
    REAL fill   !flow from lakewi to slowlake

    IF(.NOT.ALLOCATED(lakestate%slowwater)) RETURN !no lake division, could change to if ns=0
    
    !> \b Algorithm \n
    !> If the lake is divided and the slow part has room for more water:
    IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)>0.)THEN
      IF(lakestate%water(itype,i)>slowlakeini(itype,i)-lakestate%slowwater(itype,i))THEN    !fill up slowlake
        !>\li If the fast part has enough water; fill up the slow part
        fill=slowlakeini(itype,i)-lakestate%slowwater(itype,i)
        lakestate%concslow(:,itype,i)=(lakestate%conc(:,itype,i)*fill+lakestate%concslow(:,itype,i)*lakestate%slowwater(itype,i))/slowlakeini(itype,i)
        lakestate%slowwater(itype,i)=slowlakeini(itype,i)
        CALL remove_water(lakestate%water(itype,i),numsubstances,lakestate%conc(:,itype,i),fill,lakestate%conc(:,itype,i),status)
        IF(status.NE.0) CALL error_remove_water(errstring(8),subid,i,itype)
      !>\li Else empty all water from fast part in the slow part
      ELSEIF(lakestate%water(itype,i)>0.AND.slowlakeini(itype,i)-lakestate%slowwater(itype,i)>0)THEN    !move all lakewi-water to slowlake
        CALL add_water(numsubstances,lakestate%slowwater(itype,i),lakestate%concslow(:,itype,i),lakestate%water(itype,i),lakestate%conc(:,itype,i))
        lakestate%water(itype,i)=0.
        lakestate%conc(:,itype,i)=0.
      ENDIF  
    ENDIF

  END SUBROUTINE calculate_flow_within_lake

  !>\brief Flow from rating curve.
  !>Estimates average lake outflow (m3/s) during one timestep by 
  !>linearization of rating equation q = k*(w-w0)**p (further developed 
  !>from Lindström, G., 2016. Lake water levels for calibration of the 
  !> S-HYPE model. Hydrology Research 47,4, pp. 672-682. doi: 10.2166/nh.2016.019).
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Lakes - Common lake processes)
  !-----------------------------------------------------------------
  REAL FUNCTION average_flow_rating_curve(q_in,l_area,wst,k,p)

    USE MODVAR, ONLY : seconds_per_timestep, &
                       doublezero

    !Argument declarations
    REAL, INTENT(IN) :: q_in   !<inflow (m3/s)
    REAL, INTENT(IN) :: l_area !<lake area (m2)
    REAL, INTENT(IN) :: wst    !<current water level above threshold (m)
    REAL, INTENT(IN) :: k      !<rating curve coefficient
    REAL, INTENT(IN) :: p      !<rating curve exponent
    
    !Local variables
    DOUBLE PRECISION dh,h,h0,hr,qut,r,t1,t2,z

    qut = 0.D0
    dh = DBLE(q_in)*DBLE(seconds_per_timestep)/DBLE(l_area) !Inflow added in HYPE
    h0 = DBLE(wst)-dh !Initial height (m)
    IF(h0>0.D0) THEN
      t2 = DBLE(seconds_per_timestep)
      hr = h0
    ELSEIF (h0+dh>0.D0) THEN
      t1 = -DBLE(l_area)*h0/DBLE(q_in)
      t2 = DBLE(seconds_per_timestep)-t1
      hr = DBLE(q_in)*t2/DBLE(l_area)/10.D0
    ELSE
      t2 = 0.D0
    ENDIF

    IF(t2>0.D0) THEN
      r = DBLE(p)*DBLE(k)*(hr**(DBLE(p-1.)))/DBLE(l_area) !Linearized recession rate (1/sec)
      IF(r>doublezero)THEN
        z = hr+DBLE(q_in)/r/DBLE(l_area)-hr/DBLE(p)   !Auxiliary variable (m)
        h = (hr-z)*EXP(-r*t2)+z  !New height above threshold (m)
        qut = DBLE(q_in)-DBLE(l_area)*(h-h0)/DBLE(seconds_per_timestep)
        IF(qut<0.D0) qut = 0.D0
      ENDIF
    ENDIF
    average_flow_rating_curve = REAL(qut)

  END FUNCTION average_flow_rating_curve

  !>Calculate outlet lake water stage (m) in local reference system and for w-reference system
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions, Lakes - Outlet lake (olake) as a lake basin)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_olake_waterstage(i,lakewatermm,lakeareain,lakearea,lakewst,lakestate,w0ref)

    USE MODVAR, ONLY : basin,           &
                       classbasin,      &
                       dam,             &
                       damindex,        &
                       lake,            &
                       lakeindex,       &
                       lakebasin,       &
                       lakebasinindex,  &
                       missing_value,   &
                       dayno,           &
                       nsub,            &
                       slc_olake

    !Arguments declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: lakewatermm   !<outlet lake water content (mm)
    REAL, INTENT(IN)    :: lakeareain    !<outlet lake area of subbasin (m2)
    REAL, INTENT(OUT)   :: lakearea      !<outlet lake area (of subbasin or whole lake for last lakebasin) (m2)
    REAL, INTENT(OUT)   :: lakewst       !<outlet lake water stage (m)
    REAL, INTENT(OUT)   :: w0ref         !<level to be added for w-ref outlet lake water stage (m)
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<Lake state
    
    !Local variables
    INTEGER isb        !subbasin-loop index
    INTEGER dayno1, dayno2      !Starting day nr. for production periods 1 and 2    
    REAL lakewaterm    !lake water (m)
    REAL wmin          !lake threshold (sänkningsgräns)
    REAL deltaw0       !difference in lake threshold period 2 (m)
    REAL deltaw        !distance between thresholds (m)
    REAL isb_lakewst   !lake water above threshold of lake basins in lake (m)
    REAL lack          !water lacking for lake basin to reach threshold (m3)
    REAL overwmin      !water volume above wmin (to be spread over whole lake area)

    !> \b Algoritm \n
    !>Check for lake existance; return if not found    
    lakearea = lakeareain   !Default output lake area
    w0ref = 0.
    IF(lakeareain==0)THEN
      lakewst = missing_value
      RETURN
    ENDIF

    !Lake water level of current subbasin
    lakewaterm = lakewatermm * 0.001

    !Lake water reference
    IF(ALLOCATED(lakeindex))THEN
      IF(lakeindex(i)>0)THEN
        w0ref = lake(lakeindex(i))%w0ref
      ENDIF
    ENDIF
    IF(ALLOCATED(damindex))THEN
      IF(damindex(i)>0)THEN
        w0ref = dam(damindex(i))%w0ref
      ENDIF
    ENDIF

    !>Calculate lake water stage (for single olake) in local reference system
    lakewst = lakewaterm - basin(i)%lakedepth(2)

    !Lakebasin lake
    IF(ALLOCATED(lakebasinindex))THEN
      IF(lakebasinindex(i)>0)THEN  
        IF(lakebasin(lakebasinindex(i))%last)THEN
          !>If outlet of lakebasin lake; waterstage is calculated as mean for whole lake
          !>\li Calculate lake water reference at threshold
          w0ref = lake(lakebasin(lakebasinindex(i))%ilk)%w0ref
          deltaw0 = lake(lakebasin(lakebasinindex(i))%ilk)%deltaw0  
          wmin = lake(lakebasin(lakebasinindex(i))%ilk)%wmin
          dayno1 = lake(lakebasin(lakebasinindex(i))%ilk)%datum1          !Starting day nr. for production period 1
          dayno2 = lake(lakebasin(lakebasinindex(i))%ilk)%datum2          !                                period 2
          IF(wmin==missing_value)THEN
            deltaw = 0
            IF (dayno1*dayno2 > 0) THEN                                       !If both dates for different production regimes are non-zero,
              IF (dayno < dayno1 .OR. dayno >= dayno2)  deltaw = - deltaw0    !and today is not within the period 1, then threshold for period 2 applies
            ENDIF
          ELSE
            deltaw = 0. - wmin
          ENDIF
          !Check for lake basin with water stage below threshold and 
          !>\li Calculate lack of lake water in lakebasins of the lake
          lack = 0
          DO isb=1,nsub
            IF(lakebasinindex(isb)>0)THEN
              IF(lakebasin(lakebasinindex(isb))%ilk == lakebasin(lakebasinindex(i))%ilk)THEN
                isb_lakewst = lakestate%water(2,isb)
                IF(ALLOCATED(lakestate%slowwater)) isb_lakewst = isb_lakewst + lakestate%slowwater(2,isb)
                isb_lakewst = isb_lakewst * 0.001 - (basin(isb)%lakedepth(2) - deltaw)
                IF(isb_lakewst<0)THEN
                  lack = lack - isb_lakewst*classbasin(isb,slc_olake)%part*basin(isb)%area
                ENDIF
              ENDIF
            ENDIF
          ENDDO
          !>\li Calculate average water stage and area for whole lake            
          overwmin = (lakewaterm - (basin(i)%lakedepth(2)-deltaw))*lakearea - lack
          IF(overwmin>0)THEN
            lakewst = overwmin/lake(lakebasin(lakebasinindex(i))%ilk)%area - deltaw   !W>wmin
          ELSE
            overwmin = (lakewaterm - (basin(i)%lakedepth(2)-deltaw))*lakearea
            IF(overwmin>0)THEN
              lakewst = 0. - deltaw   !W=wmin
            ELSE
              lakewst = overwmin/lakearea - deltaw  !W<wmin
            ENDIF
          ENDIF
          lakearea = lake(lakebasin(lakebasinindex(i))%ilk)%area  
        ELSE
          lakewst = lakewaterm - basin(i)%lakedepth(2)
          w0ref = 0.
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE calculate_olake_waterstage

  !>Calculate outlet lake water stage (m) in local reference system adjusted for "real" regulation amplitude
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_regamp_adjusted_waterstage(i,lakeareain,lakewst,lakewstadj)

    USE MODVAR, ONLY : dam,             &
                       damindex,        &
                       lake,            &
                       lakeindex,       &
                       lakeout2index,   &
                       lakebasin,       &
                       lakebasinindex,  &
                       missing_value

    !Arguments declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: lakeareain    !<outlet lake area of subbasin (m2)
    REAL, INTENT(IN)    :: lakewst       !<outlet lake water stage (m)
    REAL, INTENT(OUT)   :: lakewstadj    !<outlet lake water stage adjusted for "real" amplitude of regulation volume (m)
    
    !Local variables
    REAL wfactor       !regulation amplitude scaling factor

    !> \b Algoritm \n
    wfactor = missing_value
    lakewstadj = missing_value
    
    !Check for lake existance; return if not found    
    IF(lakeareain==0) RETURN

    !>Get regulation amplitude adjustment factor
    !Case of single lake:
    IF(ALLOCATED(lakeindex))THEN
      IF(lakeindex(i)>0)THEN
        wfactor = lake(lakeindex(i))%wampcoeff
      ENDIF
    ENDIF
    IF(ALLOCATED(lakeout2index))THEN
      IF(lakeout2index(i)>0)THEN
        IF(wfactor==missing_value) wfactor = lake(lakeout2index(i))%wampcoeff
      ENDIF
    ENDIF
    IF(ALLOCATED(damindex))THEN
      IF(damindex(i)>0)THEN
        wfactor = dam(damindex(i))%wampcoeff
      ENDIF
    ENDIF

    !Case of lakebasin lake (last basin):
    IF(ALLOCATED(lakebasinindex))THEN
      IF(lakebasinindex(i)>0)THEN  
        IF(lakebasin(lakebasinindex(i))%last)THEN
          wfactor = lake(lakebasin(lakebasinindex(i))%ilk)%wampcoeff
        ENDIF
      ENDIF
    ENDIF

    !>Calculate adjusted lake water stage
    lakewstadj = lakewst
    IF(wfactor/=missing_value .AND. lakewst<0.) lakewstadj = lakewst*wfactor

  END SUBROUTINE calculate_regamp_adjusted_waterstage

  !>Calculate division of subbasin outlet flow into main channel and branch
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions, Rivers - Main river)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_branched_flow(i,totflow,mainflow,branchflow)

    USE MODVAR, ONLY : branchdata,    &
                       branchindex

    !Argument declaration
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: totflow       !<outflow of subbasin
    REAL, INTENT(OUT)   :: mainflow      !<flow in main channel
    REAL, INTENT(OUT)   :: branchflow    !<flow in branch

    !Local variables
    REAL part, minQmain, maxQmain, maxQbranch
    
    !> \b Algorithm \n
    !>Initialisation, default is all flow in main (single) channel
    mainflow   = totflow
    branchflow = 0.

    !>Check for branch existance and flow>0
    IF(.NOT.ALLOCATED(branchdata)) RETURN
    IF(branchindex(i)==0) RETURN
    IF(totflow == 0) RETURN
   
    !>Set current parameter values
    part = branchdata(branchindex(i))%mainpart
    maxQmain = branchdata(branchindex(i))%maxQ
    minQmain = branchdata(branchindex(i))%minQ
    maxQbranch = branchdata(branchindex(i))%maxQbranch
    
    !>Calculate flow in main channel and in branch
    mainflow = totflow
    IF(totflow>minQmain)THEN
      mainflow = part * (totflow - minQmain) + minQmain
    ENDIF
    IF(maxQmain>0 .AND. mainflow>maxQmain)THEN
      mainflow = maxQmain
    ELSEIF(maxQbranch>0 .AND. (1.-part)*(totflow-minQmain)>maxQbranch)THEN
      mainflow = totflow - maxQbranch
    ENDIF
    branchflow = totflow - mainflow
    
  END SUBROUTINE calculate_branched_flow

  !>Calculate subbasin outlet flow division into main channel and branch
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions, Rivers - Main river)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_branched_flow_new(i,totflow,simflow1,simflow2,maxProdin,minflowin,mainflow,branchflow)

    USE HYPEVARIABLES, ONLY : lakeoutlet
    USE MODVAR, ONLY : branchdata,    &
                       branchindex,   &
                       lake,   &
                       lakeindex,   &
                       lakeout2index

    !Argument declaration
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: totflow       !<outflow of subbasin (updated)
    REAL, INTENT(IN)    :: simflow1      !<outflow of main branch (not updated) if lake with 2 outlets, otherwise all flow
    REAL, INTENT(IN)    :: simflow2      !<outflow of second branch (not updated) if lake with 2 outlets
    REAL, INTENT(IN)    :: maxProdin     !<temporary maximum production flow for recalculating flow
    REAL, INTENT(IN)    :: minflowin     !<current minimum flow for recalculating flow
    REAL, INTENT(OUT)   :: mainflow      !<flow in main channel
    REAL, INTENT(OUT)   :: branchflow    !<flow in branch

    !Local variables
    REAL simtotflow
    REAL maxQprod,maxQprod2
    REAL minflow,minflow2
    LOGICAL have2outflows
    
    !> \b Algorithm \n
    !>Initialisation, default is all flow in main (single) channel
    mainflow   = totflow
    branchflow = 0.

    !>Check for branch existance and flow>0
    IF(.NOT.ALLOCATED(branchdata)) RETURN
    IF(branchindex(i)==0) RETURN
    IF(totflow == 0) RETURN
   
    !>Set current parameter values for two outlet lake
    have2outflows = .FALSE.
    IF(ALLOCATED(lakeout2index))THEN  !Second outlet of lake/dam
      IF(lakeout2index(i)>0)THEN
        have2outflows = .TRUE.
        maxQprod = lake(lakeindex(i))%mqprod
        maxQprod2 = lake(lakeout2index(i))%mqprod
        IF(lakeoutlet(branchindex(i))%change==4 .OR. &
           lakeoutlet(branchindex(i))%change==6) maxQprod = maxProdin   !Set temporary maxProd for recalculating branched flow
        IF(lakeoutlet(branchindex(i))%change==5 .OR. &
           lakeoutlet(branchindex(i))%change==7) maxQprod2 = maxProdin
        minflow = 0.
        minflow2 = 0.
        IF(lakeoutlet(branchindex(i))%change==7 .OR. &
           lakeoutlet(branchindex(i))%change==9) minflow = minflowin   !Set current minflow for recalculating branched flow
        IF(lakeoutlet(branchindex(i))%change==6 .OR. &
           lakeoutlet(branchindex(i))%change==8) minflow2 = minflowin
      ENDIF
    ENDIF
    simtotflow = simflow1 + simflow2
    
    !>Calculate main and branched flow for lake/dam with two outlets
    IF(have2outflows)THEN
      mainflow = simflow1
      branchflow = simflow2
      IF(simtotflow == totflow) RETURN
      CALL recalculate_branched_flow(lakeoutlet(branchindex(i))%change,totflow,maxQprod,maxQprod2,minflow,minflow2,mainflow,branchflow)
    ELSE
    !>Calculate flow in main channel and in branch from BranchData fractions
      CALL calculate_branched_flow(i,totflow,mainflow,branchflow)
    ENDIF
    
  END SUBROUTINE calculate_branched_flow_new

  !>Recalculate two outlet lake outflow division after total outflow has been changed
  !>
  !> \b Reference ModelDescription Chapter Rivers and lakes (Basic assumptions, Rivers - Main river)
  !------------------------------------------------------------------------------
  SUBROUTINE recalculate_branched_flow(cmethod,totflow,maxQprod,maxQprod2,minflow1,minflow2,simflow1,simflow2)

    !Argument declaration
    INTEGER, INTENT(IN) :: cmethod       !<code for change flow method
    REAL, INTENT(IN)    :: totflow       !<outflow of subbasin (updated)
    REAL, INTENT(IN)    :: maxQprod      !<maximum production flow, outlet 1
    REAL, INTENT(IN)    :: maxQprod2     !<maximum production flow, outlet 2
    REAL, INTENT(IN)    :: minflow1      !<minimum flow, outlet 1
    REAL, INTENT(IN)    :: minflow2      !<minimum flow, outlet 2
    REAL, INTENT(INOUT) :: simflow1      !<outflow of outlet 1
    REAL, INTENT(INOUT) :: simflow2      !<outflow of outlet 2

    !Local variables
    REAL simtotflow
    REAL outflow1,outflow2
    REAL add
    
    !> \b Algorithm \n
    !>Simple case; new flow is zero
    IF(totflow == 0)THEN
      simflow1 = 0.
      simflow2 = 0.
      RETURN
    ENDIF
   
    !>Initialisation
    !>Set current parameter values
    outflow1 = 0.
    outflow2 = 0.
    simtotflow = simflow1 + simflow2
    
    !>Recalculate main and branched flow for lake/dam with two outlets
    IF(cmethod==1 .OR. cmethod==4)THEN
    !>\li Hydropower plant with production flow and spill in different branches: check for maximum/current production flow
      outflow1=MIN(totflow,maxQprod)
      IF(totflow>outflow1) outflow2 = totflow - outflow1
    ELSEIF(cmethod==2 .OR. cmethod==5)THEN
      outflow2=MIN(totflow,maxQprod2)
      IF(totflow>outflow2) outflow1 = totflow - outflow2
    !>\li Minimum flow together with production flow and spill in different branches
    ELSEIF(cmethod==6 .OR. cmethod==8)THEN
      outflow2=MIN(totflow,minflow2)
      outflow1=MIN(totflow-outflow2,maxQprod)
      IF(totflow>outflow1+outflow2)THEN
        add = totflow-outflow1-outflow2
        outflow2=outflow2+add
      ENDIF
    ELSEIF(cmethod==7 .OR. cmethod==9)THEN
      outflow1=MIN(totflow,minflow1)
      outflow2=MIN(totflow-outflow1,maxQprod2)
      IF(totflow>outflow1+outflow2)THEN
        add = totflow-outflow1-outflow2
        outflow1=outflow1+add
      ENDIF
    !>\li Other outflow division: keep relation between flow in the different branches
    ELSEIF(cmethod==3)THEN
      IF(simtotflow>0.)THEN
        outflow1 = simflow1 * (totflow/simtotflow)
        outflow2 = simflow2 * (totflow/simtotflow)
      ELSE
        outflow1 = totflow   !Zero flow calculated: put flow in main channel
      ENDIF
    ELSE
      !? Unknown how to recalculate branched flow
      IF(simtotflow>0.)THEN
        outflow1 = simflow1 * (totflow/simtotflow)
        outflow2 = simflow2 * (totflow/simtotflow)
      ELSE
        outflow1 = totflow
      ENDIF
    ENDIF

    !Set changed output flows
    simflow1 = outflow1
    simflow2 = outflow2

    
  END SUBROUTINE recalculate_branched_flow

  !>Find information about how to calculate subbasin outlet flows and flow division 
  !>into main channel and branch
  !>
  !>\b Consequences Module hypevariable structure lakeoutlet will be allocated and set
  !------------------------------------------------------------------------------
  SUBROUTINE set_lake_outlets()

    USE HYPEVARIABLES, ONLY : lakeoutlet  !OUT
    USE MODVAR, ONLY : basin, &
                       branchdata,   &
                       branchindex,   &
                       missing_value, &
                       nsub, &
                       lake,    &
                       lakeindex, &
                       lakeout2index

    !Argument declaration

    !Local variables
    INTEGER i,nbranch
    REAL deltaw0,wmin             !parameter values outlet 1
    REAL regrate,maxProd          !-"-
    REAL out2maxProd,out2regrate  !parameter values outlet 2
    REAL out2deltaw0,out2wmin,out2w0rel !-"-
    LOGICAL minflow,out2minflow      !flag for minimum flow outlet 1 and 2
    

    !> \b Algorithm \n
    !>Initialisation; check and allocate
    IF(.NOT.ALLOCATED(lakeout2index)) RETURN !?
    nbranch = SIZE(branchdata)
    IF(.NOT.ALLOCATED(lakeoutlet)) ALLOCATE(lakeoutlet(nbranch))

    DO i=1,nsub
      IF(lakeout2index(i)>0)THEN
 
        !Initial values
        regrate = 0.
        maxProd = 0. 
        minflow = .FALSE.
        deltaw0 = 0.
        wmin = missing_value
        out2regrate = 0.
        out2maxProd = 0. 
        out2minflow = .FALSE.
        out2deltaw0 = 0.
        out2wmin = missing_value
        out2w0rel = 0.

        !Parameter values
        maxProd = lake(lakeindex(i))%mqprod         !Maximum production flow
        minflow = lake(lakeindex(i))%minflow        !Minimum flow
        regrate = lake(lakeindex(i))%rate           !Rating curve parameter
        deltaw0 = lake(lakeindex(i))%deltaw0        !difference in lake threshold/"dämningsgräns" period 2
        wmin = lake(lakeindex(i))%wmin              !lake threshold/"sänkningsgräns"
        out2maxProd = lake(lakeout2index(i))%mqprod         !Maximum production flow
        out2minflow = lake(lakeout2index(i))%minflow        !Minimum flow
        out2regrate = lake(lakeout2index(i))%rate           !Rating curve parameter
        out2deltaw0 = lake(lakeout2index(i))%deltaw0        !difference in lake threshold/"dämningsgräns" period 2
        out2wmin = lake(lakeout2index(i))%wmin              !lake threshold/"sänkningsgräns"
        out2w0rel = lake(lakeout2index(i))%w0ref            !w0 of outlet 2 relative to outlet 1 (w0ref)

        !Outlet types:
        !1: Production flow (defined by wmin)
        !2: Production flow with higher maximum production (defined by wmin and mqprod)
        !3: Overflow (rate)
        !4: Additional overflow in branch (rate, w0rel)
        !5: Production flow and overflow (wmin, rate)
        !6: Production flow based on two-time thresholds (defined by rate, deltaw0)
        !7: Production flow based on two-time thresholds with higher maximum production (rate, deltaw0, mqprod)
        !8: Minimum flow (defined by wmin, minflow) ("mintappning")
        !9: Minimum flow and overflow (wmin, rate, minflow) ("mintappning" and "spill")
        !>Find lake outlet type, outlet 1
        IF(wmin/=missing_value)THEN
          IF(maxProd>0.)THEN
            lakeoutlet(branchindex(i))%otype(1) = 2
          ELSEIF(minflow)THEN
            IF(regrate>0.)THEN
              lakeoutlet(branchindex(i))%otype(1) = 9
            ELSE
              lakeoutlet(branchindex(i))%otype(1) = 8
            ENDIF
          ELSEIF(regrate>0.)THEN
            lakeoutlet(branchindex(i))%otype(1) = 5
          ELSE
            lakeoutlet(branchindex(i))%otype(1) = 1
          ENDIF
        ELSEIF(regrate>0.)THEN
          IF(deltaw0>0.)THEN
            IF(maxProd>0.)THEN
              lakeoutlet(branchindex(i))%otype(1) = 7
            ELSE
              lakeoutlet(branchindex(i))%otype(1) = 6
            ENDIF
          ELSE
            lakeoutlet(branchindex(i))%otype(1) = 3
          ENDIF
        ELSE
          WRITE(6,*) 'ERROR: Not allowed combination of LakeData parameters' !?
          WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
          STOP 1
        ENDIF
       !>Find lake outlet type, outlet 2
        IF(out2wmin/=missing_value)THEN
          IF(out2maxProd>0.)THEN
            lakeoutlet(branchindex(i))%otype(2) = 2
          ELSEIF(out2minflow)THEN
            IF(out2regrate>0.)THEN
              lakeoutlet(branchindex(i))%otype(2) = 9
            ELSE
              lakeoutlet(branchindex(i))%otype(2) = 8
            ENDIF
          ELSEIF(out2regrate>0.)THEN
            lakeoutlet(branchindex(i))%otype(2) = 5
          ELSE
            lakeoutlet(branchindex(i))%otype(2) = 1
          ENDIF
        ELSEIF(out2regrate>0.)THEN
          IF(out2deltaw0>0.)THEN
            IF(out2maxProd>0.)THEN
              lakeoutlet(branchindex(i))%otype(2) = 7
            ELSE
              lakeoutlet(branchindex(i))%otype(2) = 6
            ENDIF
          ELSEIF(out2w0rel/=0)THEN
            lakeoutlet(branchindex(i))%otype(2) = 4
          ELSE
            lakeoutlet(branchindex(i))%otype(2) = 3
          ENDIF
        ELSE
          WRITE(6,*) 'ERROR: Not allowed combination of LakeData parameters' !?
          WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
          STOP 1
        ENDIF
          
        !>Set method for changing flow (after updating totalflow).
        !Method 1: main has priority and a maximum allowed value
        !Method 2: branch has priority and a maximum allowed value
        !Method 3: main and branch equal, change proportionally to old values
        !Method 4: main has priority and a wanted value
        !Method 5: branch has priority and a wanted value
        !Method 6: minflow in branch has priority, then prod in main, last outflow in branch
        !Method 7: minflow in main has priority, then prod in branch, last outflow in main
        !Method 8: minflow in branch has priority, then maxprod in main, last outflow in branch
        !Method 9: minflow in main has priority, then maxprod in branch, last outflow in main
        IF(lakeoutlet(branchindex(i))%otype(1) == 1)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 1)THEN
            lakeoutlet(branchindex(i))%change = 3   !proportionally
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 3 .OR. &
                 lakeoutlet(branchindex(i))%otype(2) == 6)THEN
            lakeoutlet(branchindex(i))%change = 4   !mainprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 5)THEN
            lakeoutlet(branchindex(i))%change = 3   !proportionally in lack of better method
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 7)THEN
            lakeoutlet(branchindex(i))%change = 2   !branchprodflag, maxProd exist
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 8)THEN
            lakeoutlet(branchindex(i))%change = 5   !minflow in branch, calculated as branchprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            lakeoutlet(branchindex(i))%change = 6   !minflow in branch, prod in main
          ELSE
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !2 and 4
            WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
            STOP 1
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 2)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 5)THEN
            lakeoutlet(branchindex(i))%change = 3   !proportionally in lack of better method
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            lakeoutlet(branchindex(i))%change = 8   !minflow in branch, prod in main
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 1 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 4 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 7 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 8)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !1,2,4,7,8
            WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
            STOP 1
          ELSE
            lakeoutlet(branchindex(i))%change = 1   !mainprodflag, maxProd exist (3,6)
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 3)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 1)THEN
            lakeoutlet(branchindex(i))%change = 5   !branchprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 8)THEN
            lakeoutlet(branchindex(i))%change = 5   !minflow in branch, calculated as branchprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
                 lakeoutlet(branchindex(i))%otype(2) == 7)THEN
            lakeoutlet(branchindex(i))%change = 2   !branchprodflag, maxProd exist
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 5 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !5,9
            WRITE(6,*) 'ERROR: Check LakeData.txt subbasin: ',basin(i)%subid
            STOP 1
          ELSE
            lakeoutlet(branchindex(i))%change = 3   !proportionally (3,4,6?)
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 5)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 2)THEN
            lakeoutlet(branchindex(i))%change = 3   !proportionally in lack of better method
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 1 .OR. &
                 lakeoutlet(branchindex(i))%otype(2) == 4)THEN
            lakeoutlet(branchindex(i))%change = 3   !proportionally in lack of better method
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 7)THEN
            lakeoutlet(branchindex(i))%change = 2   !branchprodflag, maxProd finns
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 8)THEN
            lakeoutlet(branchindex(i))%change = 5   !minflow in branch, calculated as branchprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 3 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !3,9
            WRITE(6,*) 'ERROR: Check LakeData.txt subbasin: ',basin(i)%subid
            STOP 1
          ELSE
            lakeoutlet(branchindex(i))%change = 3   !proportionally (5,6?)
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 6)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
             lakeoutlet(branchindex(i))%otype(2) == 7)THEN
            lakeoutlet(branchindex(i))%change = 2   !branchprodflag, maxProd finns
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 1)THEN
            lakeoutlet(branchindex(i))%change = 5   !branchprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 8)THEN
            lakeoutlet(branchindex(i))%change = 5   !minflow in branch, calculated as branchprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            lakeoutlet(branchindex(i))%change = 6   !minflow in branch, prod in main
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 4)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !4
            WRITE(6,*) 'ERROR: Check LakeData.txt subbasin: ',basin(i)%subid
            STOP 1
          ELSE
            lakeoutlet(branchindex(i))%change = 3   !proportionally (3?,5?,6)
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 7)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
             lakeoutlet(branchindex(i))%otype(2) == 4 .OR. &
             lakeoutlet(branchindex(i))%otype(2) == 7 .OR. &
             lakeoutlet(branchindex(i))%otype(2) == 8)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !2,4,7,8
            WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
            STOP 1
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            lakeoutlet(branchindex(i))%change = 8   !minflow in branch, maxprod in main
          ELSE
            lakeoutlet(branchindex(i))%change = 1   !mainprodflag, maxProd finns (1,3,5,6)
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 8)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 1 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 3 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 5 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 6)THEN
            lakeoutlet(branchindex(i))%change = 4   !minflow in main, calculated as mainprod without max
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 7 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 8 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !2,7,8,9
            WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
            STOP 1
          ENDIF
        ELSEIF(lakeoutlet(branchindex(i))%otype(1) == 9)THEN
          IF(lakeoutlet(branchindex(i))%otype(2) == 1 .OR. &
             lakeoutlet(branchindex(i))%otype(2) == 6)THEN
            lakeoutlet(branchindex(i))%change = 7   !minflow in main, prod in branch
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 2 .OR. &
             lakeoutlet(branchindex(i))%otype(2) == 7)THEN
            lakeoutlet(branchindex(i))%change = 9   !minflow in main, maxprod in branch
          ELSEIF(lakeoutlet(branchindex(i))%otype(2) == 3 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 4 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 5 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 8 .OR. &
            lakeoutlet(branchindex(i))%otype(2) == 9)THEN
            WRITE(6,*) 'ERROR: Outlet combination of lake not allowed.'  !3,4,5,8,9
            WRITE(6,*) 'ERROR: Check LakeData. Subbasin: ',basin(i)%subid
            STOP 1
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE set_lake_outlets

  !>Calculate different volumes of lakes for print out.
  !>Volume for ilakes, volume for olakes, and volume for whole lakes (basindivided).
  !------------------------------------------------------------------
  SUBROUTINE calculate_lake_volume(itype,i,a,lakewi,lakebasinvol,lakevol,lakevolsum)
  
    USE MODVAR, ONLY : lakebasin, &
                       lakebasinindex, &
                       nbasinlakes, &
                       missing_value                   

    !Argument declarations
    INTEGER, INTENT(IN)   :: itype            !<lake type; ilake=1, olake=2
    INTEGER, INTENT(IN)   :: i                !<index of current subbasin
    REAL, INTENT(IN)      :: a                !<lake area (m2)
    REAL, INTENT(IN)      :: lakewi           !<lake water stage (mm)
    REAL, INTENT(INOUT)   :: lakebasinvol(2)  !<volume of olake and ilake
    REAL, INTENT(INOUT)   :: lakevol          !<volume of olake/volume for lake with basins in outlet basin
    REAL, INTENT(INOUT)   :: lakevolsum(nbasinlakes)    !<to sum lakebasins to outlet basin
    
    !Local variables
    INTEGER lakeid !lake id, used for lakes with many basins (lakeid i lakedata)

    !> \b Algoritm \n
    !> Calculate lake volume for current lake
    lakebasinvol(itype) = lakewi * 0.001 * a  !volume in lake (m3)

    !> If outlet lake
    IF(itype==2)THEN
      !\li Set volume for olake
      lakevol = lakebasinvol(itype)
      !\li For basin-lakes: calculate volume of outlet, upstream sub-lakebasins volumes set to missing
      IF(ALLOCATED(lakebasinindex))THEN    
        IF(lakebasinindex(i) .NE. 0) THEN   !lakebasin
          lakeid = lakebasin(lakebasinindex(i))%ilk
          lakevol = missing_value
          lakevolsum(lakeid) = lakevolsum(lakeid) + lakebasinvol(itype)
          IF(lakebasin(lakebasinindex(i))%last) THEN  !outlet of lakebasin-lake
            lakevol = lakevolsum(lakeid)
          ENDIF
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE calculate_lake_volume

  
  !>Calculate temperature(T2) processes in rivers
  !----------------------------------------------------------
   SUBROUTINE T2_processes_in_river(i,itype,temp,swrad,riversurft,riverarea,frozenstate,riverstate,freezeupday,freezeuparea)
  
    USE MODVAR, ONLY: genpar, &
                      i_t2, &
                      modeloption, &
                      p_lakeriverice, &
                      cwater
    USE HYPEVARIABLES, ONLY: m_t2trriver, &
                             m_riceTf, &
                             m_tcfriver, &
                             m_scfriver, &
                             m_ccfriver, &
                             m_lcfriver, &
                             m_stbcorr1, &
                             m_stbcorr2, &
                             m_stbcorr3, &
                             m_limt2exch
    
    !Argument variables
    INTEGER, INTENT(IN) :: i               !<index of subbasin
    INTEGER, INTENT(IN) :: itype           !<index of river type (local = 1, main = 2)
    REAL, INTENT(IN)    :: temp            !<air temperature
    REAL, INTENT(IN)    :: swrad           !<solar radiation
    REAL, INTENT(INOUT) :: riversurft(2)   !<water surface temperature
    REAL, INTENT(IN)    :: riverarea       !<river area
    TYPE(snowicestatetype),INTENT(IN)  :: frozenstate   !<Snow and ice states
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    INTEGER,INTENT(OUT) :: freezeupday
    REAL, INTENT(INOUT) :: freezeuparea     !<fraction of riverarea with newice formation
    
    !Local variables    
    REAL t2transfcorr
    REAL watertemp, watervol,icefreefraction

    !Initiate heat deficit and freezeup flag and surface temp variables
    freezeuparea = 0.
    freezeupday = 0
    riversurft(itype) = 0.  

    !Get total river water volume and mean T2 temperature
    CALL get_rivertempvol(i,itype,riverstate,watertemp,watervol)
    watervol = watervol * 1000. / riverarea    !scale volume [m3] to depth [mm]

    IF(watervol.GT.0.)THEN    !Skip calculations if there is no water in the river

      icefreefraction = 1. - frozenstate%rivericecov(itype,i)      !Fraction of icefree river surface area
      t2transfcorr = 1.      !Seasonal correction of T2 exchange coefficient   

      !River-Atmosphere T2 exchange, only in ice-free conditions and if there is some water in the river
      IF(icefreefraction.GT.0.)THEN    
        !River-atmosphere exchange 
        ! optional models  (will be reduced to one option after some initial testing for EHYPE3.0 and SHYPE2012)
        SELECT CASE(modeloption(p_lakeriverice))
        CASE(2) ! new model based on Piccolroaz et al 2013, with modifications for fractional ice cover, and calculation of fractional freezup area
          CALL calculate_watersurface_heatbalance(temp,swrad,watertemp,watervol*riverarea*0.001,riverarea*icefreefraction, & 
                                                  genpar(m_tcfriver),genpar(m_scfriver),genpar(m_ccfriver),genpar(m_lcfriver), &
                                                  genpar(m_limt2exch),freezeuparea,genpar(m_riceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))
        CASE(1) ! the simple air-water temperature exchange model (Johan/David), with modifications for fractional ice cover, and calculation of fractional freezup area
          CALL calculate_T2_transfer(temp,watertemp,watervol*riverarea*0.001,riverarea*icefreefraction,genpar(m_t2trriver)*t2transfcorr, & 
                                     freezeuparea,genpar(m_riceTf))
        ENDSELECT
        
        !Check the freezeup conditions
        IF(freezeuparea.GT.0.)THEN
          !freezup area is the fraction of previously unfrozen area (riverarea*icefreefraction), where new ice formation is triggered
          !re-scale to a fraction of the entire river area:
          freezeuparea = freezeuparea * icefreefraction
          freezeupday = 1
        ENDIF
       
        !Assign update values to the riverstate variables
        CALL set_rivertemp(i,itype,riverstate,watertemp)
      
        !Assign river surface temperature if (partly) icefree conditions - it's later rescaled after ice calculations
        riversurft(itype) = riverstate%conc(i_t2,itype,i)
      ENDIF
    ELSE
      !Set T2 temperature to 0. if there is no water
      CALL set_rivertemp(i,itype,riverstate,0.)
    ENDIF
    
  END SUBROUTINE T2_processes_in_river

  !>Calculate ice processes in rivers
  !----------------------------------------------------------
  SUBROUTINE ice_processes_in_river(i,itype,iluse,snowfall,temp,riversurftemp,  &
                  riverarea,swrad,frozenstate,riverstate, &
                  freezeupday,breakupday,freezeuparea)
     
    USE MODVAR, ONLY: genpar
    USE HYPEVARIABLES, ONLY: m_sndens0,&
                             m_ricesndens, &
                             m_ricetf, &
                             m_ricekika, &
                             m_ricekexp, &
                             m_ricetmelt,  &
                             m_ricewme,  &
                             m_riceTf

    !Argument declaration    
    INTEGER, INTENT(IN) :: i                 !<index of subbasin
    INTEGER, INTENT(IN) :: itype             !<index of lake/river type
    INTEGER, INTENT(IN) :: iluse             !<index of landuse
    REAL,INTENT(IN)     :: snowfall          !<snowfall
    REAL,INTENT(IN)     :: temp              !<air temperature
    REAL,INTENT(INOUT)  :: riversurftemp(2)  !<water surface temperature
    REAL,INTENT(IN)     :: riverarea         !<river area, m2
    REAL,INTENT(IN)     :: swrad             !<shortwave radiation
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    INTEGER, INTENT(IN) :: freezeupday    !<status freeze-up day
    INTEGER, INTENT(OUT) :: breakupday  !<status ice break-up day
    REAL, INTENT(IN)     :: freezeuparea     !<fraction of river area with newice formation (calculated by surface heat balance function)
    
    !Local variables
    REAL driverwidt, dsnowdt
    REAL oldsnow, melt
    REAL newicesurftemp,newice,newbice,newicesnow,newicesnowdepth,oldsurftemp
    INTEGER newbreakup

    !Local parameters
    REAL, PARAMETER :: dice = 0.917    !density of ice, fraction of water
    REAL, PARAMETER :: mm2cm = 0.1
    REAL, PARAMETER :: cm2mm = 10.
    
    !RiverIceModel: Initialization of some variables
    breakupday = 0
    newicesurftemp = 0.
    newice = 0.
    newbice = 0.
    newicesnow = 0.
    newicesnowdepth = 0.
    newbreakup = 0
    
    !New ice formation on "freezeuparea" (calculated by surface heat balance function)
    IF(freezeuparea>0.)THEN
      CALL calculate_icedepth(newicesurftemp, newice, &
                              newbice,newicesnow,newicesnowdepth, & 
                              temp,driverwidt,dsnowdt,freezeupday,newbreakup, &
                              genpar(m_ricetf),genpar(m_ricekika),genpar(m_ricekexp),genpar(m_ricetmelt))
    ENDIF    

    !Calculate development of the old river ice
    IF(frozenstate%riverice(itype,i)>0.)THEN
       
      !first guess is that the old ice (or snow) is melting at 0 degrees
      oldsurftemp = 0.
      
      !Snow on riverice calculation
      oldsnow = frozenstate%riversnow(itype,i)
      CALL calculate_snow_on_ice(iluse,snowfall,temp,melt,swrad,frozenstate%riversnow(itype,i),  &
                                 frozenstate%riversnowage(itype,i))
         
      !Update snow age and calculate snow depth for snow on ice
      CALL calculate_snowdepth(iluse,frozenstate%riversnow(itype,i),oldsnow,snowfall,temp,  &
                               genpar(m_ricesndens),frozenstate%riversnowage(itype,i),frozenstate%riversnowdepth(itype,i))

      !Ice depth calculation (incl. update of skin temperature)
      CALL calculate_icedepth(oldsurftemp, frozenstate%riverice(itype,i), &
                              frozenstate%riverbice(itype,i),frozenstate%riversnow(itype,i),frozenstate%riversnowdepth(itype,i), & 
                              temp,driverwidt,dsnowdt,freezeupday,breakupday, &
                              genpar(m_ricetf),genpar(m_ricekika),genpar(m_ricekexp),genpar(m_ricetmelt))

      !If river temperature is above freezing, use the excess heat to melt some river ice from below (see further in corresponding lake routine)
      CALL riverice_riverwater_interaction(i,itype,riverstate,frozenstate,riverarea,breakupday,driverwidt)
    ENDIF

    !Add new ice to the old ice
    IF(newice>0.)THEN
      IF(frozenstate%riverice(itype,i)>0.)THEN
        frozenstate%riversnow(itype,i) = frozenstate%riversnow(itype,i)* frozenstate%rivericecov(itype,i)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        frozenstate%riversnowdepth(itype,i) = frozenstate%riversnowdepth(itype,i) * frozenstate%rivericecov(itype,i)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        frozenstate%riverice(itype,i) = (frozenstate%riverice(itype,i)*frozenstate%rivericecov(itype,i) + newice*freezeuparea)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        frozenstate%riverbice(itype,i) = (frozenstate%riverbice(itype,i)*frozenstate%rivericecov(itype,i) + newbice*freezeuparea)/(frozenstate%rivericecov(itype,i)+freezeuparea)
        riversurftemp(itype) = newicesurftemp * freezeuparea + oldsurftemp * frozenstate%rivericecov(itype,i) + riversurftemp(itype)*(1.-freezeuparea-frozenstate%rivericecov(itype,i))
        frozenstate%rivericecov(itype,i) = (frozenstate%rivericecov(itype,i)+freezeuparea)      
      ELSE
        frozenstate%riversnow(itype,i) = 0.0
        frozenstate%riversnowage(itype,i) = 0.0
        frozenstate%riversnowdepth(itype,i) = 0.0
        frozenstate%riverice(itype,i) = newice
        frozenstate%riverbice(itype,i) = newbice
        riversurftemp(itype) = newicesurftemp * freezeuparea + riversurftemp(itype)*(1.-freezeuparea)
        frozenstate%rivericecov(itype,i) = freezeuparea
        !Make sure breakupday is 0 (strange situation with complete meltout of old ice and newice formation at the same time)
        IF(breakupday==1) breakupday=0
      ENDIF
    ELSE
      !Only old ice remaining
      IF(frozenstate%riverice(itype,i).GT.0.)THEN
        !weighted surface temperature (oldice and open water surface temperature)
        riversurftemp(itype) = oldsurftemp * frozenstate%rivericecov(itype,i) + riversurftemp(itype)*(1.-frozenstate%rivericecov(itype,i))
      ELSE
        !no new snow and no old snow
        !check if there was complete meltout today, in that case make sure all variables are reset
        IF(breakupday.EQ.1)THEN
          frozenstate%riverice(itype,i) = 0.
          frozenstate%riverbice(itype,i) = 0.
          frozenstate%riversnow(itype,i) = 0.
          frozenstate%riversnowage(itype,i) = 0.
          frozenstate%riversnowdepth(itype,i) = 0.
          riversurftemp(itype) = genpar(m_riceTf) * frozenstate%rivericecov(itype,i) + riversurftemp(itype)*(1.-frozenstate%rivericecov(itype,i))
          frozenstate%rivericecov(itype,i) = 0.
        ENDIF
      ENDIF
    ENDIF

  END SUBROUTINE ice_processes_in_river
  
  !>Calculate interaction between river water and river ice 
  ! - heat from water temperature above freezing is used to melt river ice 
  !   by reducing the fractional area, rather than reducing ice depth
  ! - latent heat correspondning to ice meltwater is also added to the water
  !--------------------------------------------------------------------------
  SUBROUTINE riverice_riverwater_interaction(i, itype, riverstate, frozenstate, riverarea, breakupday, driverwidt)

    USE MODVAR, ONLY : genpar, cwater
    USE HYPEVARIABLES, ONLY : m_riceTf,m_ricewme
    
    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<current subbasin index
    INTEGER, INTENT(IN) :: itype    !<river type
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate      !<River state
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
    REAL, INTENT(IN) :: riverarea   !<river area (m2)
    INTEGER, INTENT(INOUT) :: breakupday !<status of river ice break up
    REAL, INTENT(IN) :: driverwidt
    
    
    !local variables
    REAL watertemp, watervol, icewater, meltheat, waterheat, meltwater, newwatertemp,oldicecover
    
    !parameters
    REAL, PARAMETER :: L = 3.35E5     !latent heat of freezing, J/kg
    REAL, PARAMETER :: dice = .917    !density of ice, fraction of water
    REAL, PARAMETER :: mm2cm = 0.1
    REAL, PARAMETER :: cm2mm = 10.
  
    !Get total river water volume and mean T2 temperature
    CALL get_rivertempvol(i,itype,riverstate,watertemp,watervol)
    
    !scale volume [m3] to depth [mm], volume water per unit area
    watervol = watervol * 1000. / (riverarea)
    
    oldicecover = frozenstate%rivericecov(itype,i)
    
    IF(watervol.GT.0.)THEN
    
      !available heat for melting (C * KG/M2 * 1000 * KJ/KG/C = J/M2), per unit area
      waterheat = (watertemp-genpar(m_riceTf)) * watervol * 1000. * cwater
      
      IF(waterheat.GT.0.)THEN
        !Try to melt some ice, but only if the ice did not already melt completely (breakupday==1)
        IF(breakupday.EQ.0)THEN
          ! !melt the ice, from below, in cm ice
          ! bottommelt = min(frozenstate%riverice(itype,i),waterheat/(L*dice)*mm2cm)
          ! meltheat   = bottommelt * (L*dice) * cm2mm
          ! meltwater = bottommelt * dice *cm2mm
          
          !river ice and snow mass, in mm water, per unit area of ice covered river
          icewater = frozenstate%riverice(itype,i)*dice*cm2mm + frozenstate%riversnow(itype,i)
          
          !ice melt, in mm per unit area of ice-covered river 
          ! - it is thus unly the water below the ice which is interacting with the ice
          ! - the available heat is scaled with the "meltefficiency" parameter 
          meltwater = min(icewater,genpar(m_ricewme)*waterheat/L)
          meltheat = meltwater * L
             
! 3) update the frozen states with bottom melt
          
          !frozenstate%riverice(itype,i)=max(0.,frozenstate%riverice(itype,i)-bottommelt)
          !IF(frozenstate%riverice(itype,i).GT.0.)THEN
          IF((icewater-meltwater).GT.0.)THEN
            !some ice remains, reduce ice content by reducing the fractional area
            frozenstate%rivericecov(itype,i) = min(1.,max(0.,frozenstate%rivericecov(itype,i)*(1.- meltwater/icewater)))           
!            frozenstate%riverbice(itype,i)=max(0.,frozenstate%riverbice(itype,i)-bottommelt)
          ELSE
            !complete melt of the riverice
            frozenstate%riverice(itype,i) =0.
            frozenstate%riverbice(itype,i)=0.
            
            !add heat needed to melt the riversnow to the meltheat
            !meltheat = meltheat + frozenstate%riversnow(itype,i) * L
              
            !add snow to the meltwater
            !meltwater = meltwater + frozenstate%riversnow(itype,i)
              
            !reset the snow states
            frozenstate%riversnow(itype,i)=0.
            frozenstate%riversnowage(itype,i)=0.
            
            !and the ice cover state
            frozenstate%rivericecov(itype,i)=0.
            
            !set breakup flag to 1
            breakupday = 1
          ENDIF
        ELSE
          !Ice was already melted away by the icedepth function
          meltheat   = 0.
          meltwater  = 0.
        ENDIF
      ELSE
        meltheat = 0.
        meltwater = 0.
      ENDIF
! 4) use any remaining heat and the zero degree melt water to update the river state
        
      !remove melt heat from heat content of the lake water (this is now per unit area previously ice covered river)
      waterheat = waterheat - meltheat
      
      !add any previous surface melt water to the meltwater
      IF(driverwidt.GT.0)THEN
        meltwater = meltwater + driverwidt
      ENDIF
      
      !temperature of water from remaining heat content
      newwatertemp=max(waterheat/(watervol * 1000. * cwater) + genpar(m_riceTf),genpar(m_riceTf))
        
      !dilute with the meltwater, which is at freezing point
      newwatertemp=max(genpar(m_ricetf),newwatertemp * (watervol - meltwater)/watervol)
       
      !weighted temperature between (previously) ice covered and ice free water
      watertemp = watertemp * (1.-oldicecover) + newwatertemp * oldicecover
       
      !finally, assign update values to the riverstate variables
      CALL set_rivertemp(i,itype,riverstate,watertemp)
   ENDIF
  
  END SUBROUTINE riverice_riverwater_interaction
  
  !>Calculate lake (typical or current) epilimnion depth
  !>
  !>\b Reference Hanna (1990) Evaluation of models predicting mixing
  !>depth. Can. J. Fish. Aquat. Sci. 47:940-947.
  !----------------------------------------------------------
  SUBROUTINE calculate_lake_epilimnion_depth(lakearea,prec,evapl,qinmm,epidepth)

    !Argument declarations
    REAL,INTENT(IN)   :: lakearea   !<lake area (m2)
    REAL,INTENT(IN)   :: prec       !<precipitation (mm/ts)
    REAL,INTENT(IN)   :: evapl      !<evaporation (mm/ts)
    REAL,INTENT(IN)   :: qinmm      !<inflow (mm/ts)
    REAL, INTENT(OUT) :: epidepth   !<lake epilimnion depth (m)

    !Typical depth to thermocline, function of lake area (Hanna, 1990)
    epidepth = 6.95 * (lakearea * 1.0E-6)**0.185
    
    !Adapted epilimnion depth due to precipitation, evaporation and lake inflow
    !to avoid problems with low turnover time
    epidepth = epidepth + prec * 0.001 - evapl * 0.001
    epidepth = epidepth + qinmm * 0.001
    
  END SUBROUTINE calculate_lake_epilimnion_depth

  !>Subroutine for calculation of snow on ice changes; snowfall addition and
  !>snow pack melting
  !------------------------------------------------------------------------
  SUBROUTINE calculate_snow_on_ice(iluse,snowfall,temp,melt,swrad,snow,snowage)
  
    USE MODVAR, ONLY : genpar
    USE HYPEVARIABLES, ONLY : m_licewcorr

    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    REAL, INTENT(IN)    :: snowfall !<precipitation as snow (mm/timestep) 
    REAL, INTENT(IN)    :: temp     !<air temperature (C)
    REAL, INTENT(OUT)   :: melt     !<snow melt (mm/timestep)
    REAL, INTENT(IN)    :: swrad    !<shortwave radiation (MJ/m2/day?)
    REAL, INTENT(INOUT) :: snow     !<snow pack (mm)
    REAL, INTENT(INOUT) :: snowage  !<snowage (timesteps)

    !Local variables
    REAL newsnow
    REAL snowcover

    !>\b Algorithm \n
    !>Set parameter values
    snowcover = 1.    !just set snowcover = 1., and introduce snowcover calculation on lake ice later...
    
    !>Calculate snow melt
    CALL calculate_snowmelt(iluse,temp,swrad,snow,snowage,snowcover,melt)
    melt = MAX(0.,MIN(melt, snow))  !Safeguard
    
    !>Update the snow pack with snowfall and melting
    newsnow = MAX(0.,snow + genpar(m_licewcorr)*snowfall  - melt)
    snow = newsnow

  END SUBROUTINE calculate_snow_on_ice

  !>Calculate lake ice processes
  !----------------------------------------------------------
  SUBROUTINE ice_processes_in_lake(i,itype,iluse,snowfall,temp,lakesurftemp,  &
                                   swrad,frozenstate,lakestate,freezeupday, &
                                   breakupday,epidepth,freezeuparea)
    
    USE MODVAR, ONLY: genpar
    USE HYPEVARIABLES, ONLY: m_sndens0, &
                             m_licesndens, &
                             m_licetf,   &
                             m_licekika, &
                             m_licekexp, &
                             m_licetmelt,  &
                             m_licewme,  &
                             m_liceTf
           
    !Argument declarations
    INTEGER, INTENT(IN) :: i                !<index of subbasin
    INTEGER, INTENT(IN) :: itype            !<index of lake/river type
    INTEGER, INTENT(IN) :: iluse            !<index of landuse
    REAL,INTENT(IN)     :: snowfall         !<snowfall
    REAL,INTENT(IN)     :: temp             !<air temp
    REAL,INTENT(INOUT)  :: lakesurftemp(2)  !<water surface temperature
    REAL,INTENT(IN)     :: swrad            !<shortwave radiation
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    INTEGER, INTENT(IN) :: freezeupday   !<status freeze-up day
    INTEGER, INTENT(OUT) :: breakupday   !<status ice break-up day
    REAL, INTENT(IN)    :: epidepth   !<epilimnion depth (m)
    REAL, INTENT(IN)    :: freezeuparea    !<fractional water surface area with newice formation, given by temperature routine
    
    !Local variables
    REAL dlakewidt, dsnowdt
    REAL oldsnow, melt
    REAL newicesurftemp,newice,newbice,newicesnow,newicesnowdepth,oldsurftemp
    INTEGER newbreakupday
   
    !Initialization
    breakupday  = 0
    newicesurftemp = 0.
    newice = 0.
    newbice = 0.
    newicesnow = 0.
    newicesnowdepth = 0.
    newbreakupday=0
    
    !Newice formation on "freezeuparea" (calculated by surface heat balance function)
    IF(freezeuparea.GT.0.)THEN
      CALL calculate_icedepth(newicesurftemp, newice, &
                              newbice,newicesnow,newicesnowdepth, & 
                              temp,dlakewidt,dsnowdt,freezeupday,newbreakupday, &
                              genpar(m_licetf),genpar(m_licekika),genpar(m_licekexp),genpar(m_licetmelt))
    ENDIF

    !Calculate development of the old lake ice
    IF(frozenstate%lakeice(itype,i).GT.0)THEN
       !FROZEN LAKE

       !first guess is that the ice (or snow) is melting at 0 degrees
       oldsurftemp = 0.0
       
       !snow on lakeice calculation
       oldsnow = frozenstate%lakesnow(itype,i)
       CALL calculate_snow_on_ice(iluse,snowfall,temp,melt,swrad,frozenstate%lakesnow(itype,i), &
                                  frozenstate%lakesnowage(itype,i))
                  
       !Update snow age and snow depth for snow on ice
       CALL calculate_snowdepth(iluse,frozenstate%lakesnow(itype,i),oldsnow,snowfall,temp, &
                                genpar(m_licesndens),frozenstate%lakesnowage(itype,i),frozenstate%lakesnowdepth(itype,i))
      
       !Ice depth calculation (inlc. update of skin temperature)
       CALL calculate_icedepth(oldsurftemp,frozenstate%lakeice(itype,i),frozenstate%lakebice(itype,i), &
                               frozenstate%lakesnow(itype,i),frozenstate%lakesnowdepth(itype,i), & 
                               temp,dlakewidt,dsnowdt,freezeupday,breakupday, &
                               genpar(m_licetf),genpar(m_licekika),genpar(m_licekexp),genpar(m_licetmelt))
       
       !Calculate bottom melt due to heat from lake water temperatures above freezing, as well as influence of surface melt on lake water temperature
       CALL calculate_lakeice_lakewater_interaction(itype,i,frozenstate,lakestate,dlakewidt,epidepth,breakupday)
       
    ENDIF
    
    !Add new ice to the old ice
    IF(newice.GT.0.)THEN
      IF(frozenstate%lakeice(itype,i).GT.0.)THEN
         frozenstate%lakesnow(itype,i) = frozenstate%lakesnow(itype,i)* frozenstate%lakeicecov(itype,i)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         frozenstate%lakesnowdepth(itype,i) = frozenstate%lakesnowdepth(itype,i) * frozenstate%lakeicecov(itype,i)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         frozenstate%lakeice(itype,i) = (frozenstate%lakeice(itype,i)*frozenstate%lakeicecov(itype,i) + newice*freezeuparea)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         frozenstate%lakebice(itype,i) = (frozenstate%lakebice(itype,i)*frozenstate%lakeicecov(itype,i) + newbice*freezeuparea)/(frozenstate%lakeicecov(itype,i)+freezeuparea)
         lakesurftemp(itype) = newicesurftemp * freezeuparea + oldsurftemp * frozenstate%lakeicecov(itype,i) + lakesurftemp(itype)*(1. - freezeuparea - frozenstate%lakeicecov(itype,i))
         frozenstate%lakeicecov(itype,i) = (frozenstate%lakeicecov(itype,i)+freezeuparea)  
      ELSE
         frozenstate%lakesnow(itype,i) = 0.0
         frozenstate%lakesnowage(itype,i) = 0.0
         frozenstate%lakesnowdepth(itype,i) = 0.0
         frozenstate%lakeice(itype,i) = newice
         frozenstate%lakebice(itype,i) = newbice
         lakesurftemp(itype) = newicesurftemp * freezeuparea + lakesurftemp(itype)*(1. - freezeuparea)
         frozenstate%lakeicecov(itype,i) = freezeuparea
         !Make sure breakupday is 0 (strange situation with complete meltout of old ice and newice formation at the same time)
         IF(breakupday==1) breakupday=0
      ENDIF
    ELSE
      !Or just check breakup conditions of old ice, and/or update the lakesurf temperature
      IF(frozenstate%lakeice(itype,i).GT.0.)THEN
        lakesurftemp(itype) = oldsurftemp * frozenstate%lakeicecov(itype,i) + lakesurftemp(itype)*(1. - frozenstate%lakeicecov(itype,i))
      ELSE
        !no new snow and no old snow
        !check if there was complete meltout today, in that case make sure all variables are reset
        IF(breakupday.EQ.1)THEN
          frozenstate%lakeice(itype,i) = 0.
          frozenstate%lakebice(itype,i) = 0.
          frozenstate%lakesnow(itype,i) = 0.
          frozenstate%lakesnowage(itype,i) = 0.
          frozenstate%lakesnowdepth(itype,i) = 0.0
          lakesurftemp(itype) = genpar(m_liceTf) * frozenstate%lakeicecov(itype,i) + lakesurftemp(itype)*(1.-frozenstate%lakeicecov(itype,i))
          frozenstate%lakeicecov(itype,i) = 0.
        ENDIF
      ENDIF   
    ENDIF
    
  END SUBROUTINE ice_processes_in_lake
  
  !>Calculate lake ice melt from heat from lake water, as well as influence of ice surface melt on lake water temperature
  ! - depending on lake type (fast and slow split or not, deep or shallow), a mean water temperature and water volume is
  !   calculated for the interaction with the lake ice. The resulting watertemperature is then assigned to the 
  !   various lake water components
  ! - heat from water temperature above freezing is used to melt lake ice 
  !   by reducing the fractional area, rather than reducing ice depth
  ! - latent heat correspondning to ice meltwater is also added to the water
  SUBROUTINE calculate_lakeice_lakewater_interaction(itype,i,frozenstate,lakestate,dlakewidt,epidepth,breakupday)

    USE MODVAR, ONLY: genpar, lakedataparindex,lakedatapar,i_t2,cwater
    USE HYPEVARIABLES, ONLY: m_lddeeplake, m_ldfastlake, m_licetf, m_licewme

    !Argument declarations
    INTEGER,INTENT(IN) :: i             !<index of subbasin
    INTEGER,INTENT(IN) :: itype         !<index of lake type (ilake = 1, olake = 2)
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    REAL, INTENT(IN)   :: dlakewidt
    REAL, INTENT(IN)   :: epidepth    !<epilimnion depth (m)
    INTEGER, INTENT(INOUT) :: breakupday
    
    !Local variables    
    INTEGER laketype
    REAL meantemp,meanwater,waterheat,meltheat,watertemp,watervol,icewater,meltwater,newwatertemp,oldicecov
    REAL fastoutpart

    !parameters
    real, parameter :: L = 3.35E5     ! latent heat of freezing, J/kg
    real, parameter :: dice = .917    ! density of ice, fraction of water
    real, parameter :: mm2cm = 0.1
    real, parameter :: cm2mm = 10.
    
!--------------------------------------------------------------------------------------
! lakewater-lakeice interaction:
!
! 1) find out how much water and at what temperature we have for melting ice from below
! 2) melt corresponding ice (from below: black ice, slush ice, snow)
! 3) update the frozen states
! 4) use remaining heat and heat from melt water to update the lake state
!
! the first and last step is complicated by the various lake water storage configurations
!--------------------------------------------------------------------------------------    
    oldicecov = frozenstate%lakeicecov(itype,i)
    
!1)find out how much water and at what temperature we have for melting ice from below
    
    IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)==0)THEN
       !--------------------------------------
       !Lake model WITHOUT fast and slow split
       !--------------------------------------
       
       !Calculate total water stage (WATER+SLOWWATER) and average temperature (CONC*WATER+CONCSLOW*SLOWWATER)/(WATER+SLOWWATER)
       meanwater = lakestate%water(itype,i)+lakestate%slowwater(itype,i)
       IF(meanwater.GT.0.)THEN
         meantemp  = (lakestate%water(itype,i)*lakestate%conc(i_t2,itype,i)+lakestate%slowwater(itype,i)*lakestate%concslow(i_t2,itype,i))/meanwater
       
         !Check lake depth, if thermal stratification
         IF(itype==2 .AND. epidepth < meanwater*0.001 .AND. epidepth>0.)THEN !why is this only possible for olakes?
           !Two-layer olake, waterdepth > thermocline, olake
    
           !->derive lake uppertemp(t) from meantemp(t, preliminary) and lowertemp(t-1)
           lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth

           !temperature and water volume interacting with the ice
           watertemp = lakestate%uppertemp(itype,i)
           watervol  = epidepth*1000.
           laketype  = 2
         ELSE
           !one-layer olake or ilakes
           watertemp = meantemp
           watervol  = meanwater
           laketype  = 1
         ENDIF
       ELSE
         !no water in the lake, do nothing
         laketype = 0
       ENDIF
    ELSE
       !--------------------------------------
       !Lake model WITH fast and slow split
       !--------------------------------------
       fastoutpart = 0.
       IF(lakestate%water(itype,i)+lakestate%slowwater(itype,i)>0.) fastoutpart=lakedatapar(lakedataparindex(i,itype),m_ldfastlake)*lakestate%water(itype,i)/(lakestate%water(itype,i)+lakestate%slowwater(itype,i))
       
       !areally weighted averaged water stage
       meanwater = lakestate%slowwater(itype,i)*(1-fastoutpart)+lakestate%water(itype,i)*fastoutpart
       
       IF(meanwater.GT.0.)THEN
         !calculate weighted average mean temperature
         meantemp = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / meanwater
       
         IF(itype==2 .AND. epidepth < meanwater*0.001 .AND. epidepth>0.)THEN
           !Deep lake with thermal stratification (assume they have a common lowertemp)
           
           !calculate upper temp from meantemp and the lowertemp
           lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
           
           !temperature and water volume interacting with the ice
           watertemp = lakestate%uppertemp(itype,i)
           watervol  = epidepth*1000.
           laketype  = 4
         ELSE
           !One-layer lake
           watertemp = meantemp
           watervol  = meanwater
           laketype  = 3
         ENDIF
       ELSE
         !no water in the lake - do nothing
         laketype = 0
       ENDIF        
    ENDIF
    
! 2) melt corresponding ice (from below: black ice, slush ice, snow), takin fractional ice cover into account
    IF(laketype.GT.0)THEN
      !available heat for melting (C * KG/M2 * 1000 * KJ/KG/C = J/M2)
      waterheat = (watertemp-genpar(m_liceTf)) * watervol * 1000. * cwater 
       
      IF(waterheat.GT.0.)THEN
      
        !Try bottom melt only if there was not already complete meltout (breakupday==1)
        IF(breakupday.EQ.0)THEN
          !!melt the ice from below, in cm ice
          !bottommelt = min(frozenstate%lakeice(itype,i),waterheat/(L*dice)*mm2cm)
          !meltheat   = bottommelt * (L*dice) * cm2mm
          !meltwater  = bottommelt * dice *cm2mm
          
          !lake ice and snow mass, in mm water, per unit area of ice covered lake
          icewater = frozenstate%lakeice(itype,i)*dice*cm2mm + frozenstate%lakesnow(itype,i)
          
          !ice melt, in mm per unit area of ice-covered river
          ! - it is thus only the water below the ice which is interacting with the ice
          ! - the available heat is scaled with a "Meltefficiency" parameter
          meltwater = MIN(icewater,genpar(m_licewme)*waterheat/L)
          meltheat = meltwater * L
        
! 3) update the frozen states
          IF((icewater-meltwater).GT.0.)THEN
            !some ice remains, redice icemass by reducing fractional coverage
            frozenstate%lakeicecov(itype,i) = MIN(1.,MAX(0.,frozenstate%lakeicecov(itype,i)*(1-meltwater/icewater)))
          ELSE
            !complete melt of the lakeice
            frozenstate%lakeice(itype,i)=0.
            frozenstate%lakebice(itype,i)=0.
        
            !add heat needed to melt the lakesnow to the meltheat
            !meltheat = meltheat + frozenstate%lakesnow(itype,i) * L
            
            !add snow to the meltwater
            !meltwater = meltwater + frozenstate%lakesnow(itype,i)
            
            !reset the snow states
            frozenstate%lakesnow(itype,i)=0.
            frozenstate%lakesnowage(itype,i)=0.
            
            !and ice cover area
            frozenstate%lakeicecov(itype,i) = 0.
            
            !set breakupflag to 1
            breakupday = 1
          ENDIF
        ELSE
          meltheat = 0.
          meltwater = 0.
        ENDIF
      ELSE
        meltheat = 0.
        meltwater = 0.
      ENDIF
! 4) use any remaining heat and the zero degree melt water to update the lake state
      
      !remove melt heat from heat content of the lake water
      waterheat = waterheat - meltheat
      
      !add any previous surface melt water to the meltwater
      IF(dlakewidt.GT.0)THEN
        meltwater = meltwater + dlakewidt
      ENDIF
      
      !temperature of water from remaining heat content
      newwatertemp=MAX(waterheat/(watervol * 1000. * cwater) + genpar(m_liceTf),genpar(m_liceTf))
      
      !dilute with the meltwater, which is at freezing point
      newwatertemp = MAX(genpar(m_liceTf),newwatertemp * (watervol - meltwater)/watervol)
      
      !weighted temperature, between icefree and icecovered water
      watertemp = oldicecov * newwatertemp + (1.-oldicecov)*watertemp
      
      !finally, assign update values to the real state variable
      SELECT CASE(laketype)
      
        CASE(1) !single layer without split
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = watertemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = watertemp
          lakestate%uppertemp(itype,i) = watertemp
          lakestate%lowertemp(itype,i) = watertemp
        
        CASE(2) !two-layer without split
          lakestate%uppertemp(itype,i) = watertemp
          meantemp = (lakestate%uppertemp(itype,i) * epidepth + (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i) )/ (meanwater * 0.001)
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = meantemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = meantemp
          
        CASE(3) !single layer with split
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = watertemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = watertemp
          lakestate%uppertemp(itype,i) = watertemp
          lakestate%lowertemp(itype,i) = watertemp
        
        CASE(4) !two-layer with split
          lakestate%uppertemp(itype,i) = watertemp
          meantemp = (lakestate%uppertemp(itype,i) * epidepth + (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i) )/ (meanwater * 0.001)
          IF(lakestate%water(itype,i).GT.0.)THEN
            lakestate%conc(i_t2,itype,i) = meantemp
          ELSE
            lakestate%conc(i_t2,itype,i) = 0.
          ENDIF
          lakestate%concslow(i_t2,itype,i) = meantemp
      
      END SELECT
    
    ENDIF !if laketype = 0, no water in lake -> do nothing

  END SUBROUTINE calculate_lakeice_lakewater_interaction

  !>Calculate lake T2 temperature processes
  !
  ! the concept with lakestate%WATER & lakestate%SLOWWATER and the corresponding concentrations lakestate%CONC & lakestate%CONCSLOW
  ! is partly incompatible with the conceptual model for lake temperature and it's vertical distribution (uppertemp and lowertemp)
  !
  ! it's also important to notice that even a lake without partitioning into fastpart and slowpart is still using 
  !    the state variables SLOWWATER and CONCSLOW. The difference is that for a non-splitted lake, everything left in WATER and CONC is transfered
  !    to SLOWWATER and CONCSLOW at the end of the timestep (SUBROUTINE calculate_flow_within_lake). However, this takes
  !    place AFTER the lake temperature (this routine) and ice calculations. Thus, at this point the total water stage and average temperature
  !    must take both WATER and SLOWWATER into account to be correct. If we apply the same concept on the splitted lakes, we can calculate an 
  !    average upper and lower lake temperature, but keeping individual average temperatures in the two parts of the lake.
  !  This comment is also valid for the lakeice subroutine.
  !----------------------------------------------------------
  SUBROUTINE T2_processes_in_lake(i,itype,temp,swrad,lakesurft,lakearea,epidepth,frozenstate,lakestate,freezeup,freezeuparea)

    USE MODVAR, ONLY: genpar, &
                      lakedataparindex, &
                      lakedatapar, &
                      i_t2,                &
                      modeloption,         &
                      p_lakeriverice
    USE HYPEVARIABLES, ONLY: m_lddeeplake, &
                             m_t2trlake,   &
                             m_ldfastlake, &
                             m_upper2deep, &
                             m_liceTf,     &
                             m_tcflake,    &
                             m_scflake,    &
                             m_ccflake,    &
                             m_lcflake,    &
                             m_stbcorr1,  &
                             m_stbcorr2,  &
                             m_stbcorr3,  &
                             m_limt2exch

    !Argument declarations
    INTEGER,INTENT(IN) :: i             !<index of subbasin
    INTEGER,INTENT(IN) :: itype         !<index of lake type (ilake = 1, olake = 2)
    REAL,INTENT(IN)    :: temp          !<air temp
    REAL,INTENT(IN)    :: swrad         !<shortwave radiation, MJ/m2/day
    REAL,INTENT(INOUT) :: lakesurft(2)  !<water surface temperature
    REAL,INTENT(IN)    :: lakearea      !<lake area (m2)
    REAL,INTENT(IN)    :: epidepth      !<epilimnion depth (m)
    TYPE(snowicestatetype),INTENT(IN)  :: frozenstate   !<Snow and ice states
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state
    INTEGER, INTENT(OUT) :: freezeup    !<is water cooling below freezing piont (1 yes, 0 no)?
    REAL, INTENT(OUT)  :: freezeuparea  !<fraction of lake area with newice formation

    !Local variables    
    LOGICAL epilimnion
    REAL meantemp, meanwater
    REAL t2transfcorr, fastoutpart
    REAL icefreefraction, freezeuparea2

    !0 Some initializations
    freezeup = 0
    epilimnion = .FALSE.
    freezeuparea = 0.
    freezeuparea2 = 0.

    !1 Lake-atmosphere T2 exchange
    ! 
    !1.1 Seasonal correction of T2 exchange coefficient   
    t2transfcorr = 1.  !Modify according to Johans suggestion below?
    
    !1.2 Depth to thermocline, function of lake area (REF)
!    epidepth = 6.95 * (lakearea / 1000000)**0.185
    
    !1.3 Lake-Atmosphere T2 exchange, depending on lake distribution type
    IF(lakedatapar(lakedataparindex(i,itype),m_lddeeplake)==0)THEN

      !Lake model without fast and slow split
       
      !Calculate total water stage (WATER+SLOWWATER) and average temperature (CONC*WATER+CONCSLOW*SLOWWATER)/(WATER+SLOWWATER)
      meanwater = lakestate%water(itype,i)+lakestate%slowwater(itype,i)
      IF(meanwater.GT.0.)THEN
        meantemp  = (lakestate%water(itype,i)*lakestate%conc(i_t2,itype,i)+lakestate%slowwater(itype,i)*lakestate%concslow(i_t2,itype,i))/meanwater
       
        !Check lake depth, if thermal stratification
        IF(itype==2 .AND. epidepth < meanwater*0.001 .AND. epidepth>0.)THEN !why is this only possible for olakes?
          !!Two-layer olake, waterdepth > thermocline, olake
          epilimnion = .TRUE.
          
          !->derive lake uppertemp(t) from meantemp(t, preliminary) and lowertemp(t-1)
          lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
!          olduppertemp = lakestate%uppertemp(itype,i) 

          !Introducing fractional ice cover to get smoother transition over the freezing point
          icefreefraction = 1. - frozenstate%lakeicecov(itype,i)
           
          !->exchange with atmosphere - if no ice - updating meantemp(t) and uppertemp(t)
          IF(icefreefraction.GT.0.)THEN
            !temperature flow calculated from (temp-uppertemp), updating the mean temperature (concslow(i_t2,:,:))
            ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
            SELECT CASE(modeloption(p_lakeriverice))
            CASE(2) !new model based on Piccolroaz et al 2013, modified for fractional ice cover and newice formation
              CALL calculate_watersurface_heatbalance(temp-1. * (lakestate%uppertemp(itype,i) - meantemp),swrad,meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction, & 
                                                      genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                      genpar(m_limt2exch),freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CASE(1) !original function of Johan/David, modified for fractional ice cover and newice formation
              CALL calculate_T2_transfer(temp-1. * (lakestate%uppertemp(itype,i) - meantemp),meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction, & 
                                         genpar(m_t2trlake)*t2transfcorr,freezeuparea,genpar(m_liceTf)) !JS4
            END SELECT

            !->re-calculate uppertemp
            lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth

            !Check freezeup conditions, indicated by relative freezeuparea
            IF(freezeuparea.GT.0.)THEN
              !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
              !re-scale to a fraction of the entire waterarea:
              freezeuparea = freezeuparea * icefreefraction
              freezeup = 1
            ENDIF
!            !check for freezeup conditions on the new uppertemp as well, CP160620 test
!            IF(lakestate%uppertemp(itype,i).LT.genpar(m_liceTf))THEN
!              freezeup(itype) = 1
!              !estimate a freezuparea (reduction in the open water surface area) so that the result of the surface heat balance is equal to the freezing point 
!              freezeuparea2 = MAX(0.,MIN(1.,(genpar(m_liceTf) - lakestate%uppertemp(itype,i))/(olduppertemp-lakestate%uppertemp(itype,i))))
!              IF(freezeuparea==0)THEN
!                freezeuparea = freezeuparea2 * icefreefraction
!              ELSE
!                freezeuparea = MAX(0.,MIN(1.,freezeuparea+freezeuparea2*(icefreefraction-freezeuparea)))
!              ENDIF
!              lakestate%uppertemp(itype,i) = genpar(m_liceTf)
!              lakestate%lowertemp(itype,i) = (meantemp * meanwater * 0.001 - genpar(m_liceTf)*epidepth)/(meanwater * 0.001 - epidepth)
!!              WRITE(6,*) i,itype,lakestate%lowertemp(itype,i)
!            ENDIF           
            
          ENDIF
        ELSE
          !Otherwise, single-layer, ilake

          !Introducing fractional ice cover to get smoother transition over the freezing point
          icefreefraction = 1. - frozenstate%lakeicecov(itype,i)
          !!->exchange with atmosphere, if no ice, update meantemp(t)
          IF(icefreefraction.GT.0.)THEN
            ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
            SELECT CASE(modeloption(p_lakeriverice))
            CASE(2) ! new model based on Piccolroaz et al 2013
              CALL calculate_watersurface_heatbalance(temp,swrad,meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction,genpar(m_tcflake), & 
                                                      genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                      genpar(m_limt2exch),freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CASE(1)
              CALL calculate_T2_transfer(temp,meantemp,meanwater*lakearea*0.001,lakearea*icefreefraction,genpar(m_t2trlake)*t2transfcorr, &
                                         freezeuparea,genpar(m_liceTf))
            END SELECT

            !Check freezeup conditions, indicated by relative freezeuparea
            IF(freezeuparea.GT.0.)THEN
              !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
              !re-scale to a fraction of the entire waterarea:
              freezeuparea = freezeuparea * icefreefraction
              freezeup = 1
            ENDIF
          ENDIF

          lakestate%uppertemp(itype,i) = meantemp
          lakestate%lowertemp(itype,i) = meantemp
        ENDIF
      ELSE
        !no water in the lake, set temperature to 0
        meantemp  = 0.
        icefreefraction = 0.  !CP150506
      ENDIF
 
      !Finally, assign the updated meantemp to the lakestate%conc AND lakestate%concslow 
      lakestate%concslow(i_t2,itype,i) = meantemp
      lakestate%conc(i_t2,itype,i)     = meantemp
             
    ELSE
      !
      !Lakes with split in fast and slow part
      fastoutpart = 0.
      IF(lakestate%water(itype,i)+lakestate%slowwater(itype,i)>0.) fastoutpart=lakedatapar(lakedataparindex(i,itype),m_ldfastlake)*lakestate%water(itype,i)/(lakestate%water(itype,i)+lakestate%slowwater(itype,i))
       
      !areally weighted averaged water stage
      meanwater = lakestate%slowwater(itype,i)*(1-fastoutpart)+lakestate%water(itype,i)*fastoutpart
       
      IF(itype==2 .AND. epidepth < meanwater*0.001 .AND. epidepth>0.)THEN
        !Deep lake with thermal stratification (assume they have a common lowertemp)
        epilimnion = .TRUE.
          
        !calculate weighted average mean temperature
        IF(meanwater.GT.0.)THEN
          meantemp = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / meanwater
        ELSE
          meantemp=0.   !This else is unnecsessary now
          lakestate%lowertemp(itype,i) = 0.
          lakestate%uppertemp(itype,i) = 0.
        ENDIF        
         
        !calculate upper temp from meantemp and the lowertemp
        lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
!        olduppertemp = lakestate%uppertemp(itype,i) 
        !Introducing fractional ice cover to get smoother transition over the freezing point
        icefreefraction = 1. - frozenstate%lakeicecov(itype,i)

        !atmosphere-lake T2 transfer, if there is no ice
        IF(icefreefraction.GT.0.)THEN
          !atmosphere-lake interaction calculated for fastpart and slowpart separately (the same temperature gradient is effectively used here: temp-uppertemp)
          ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
          SELECT CASE(modeloption(p_lakeriverice))
          CASE(2) ! new model based on Piccolroaz et al 2013
            CALL calculate_watersurface_heatbalance(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%conc(i_t2,itype,i)),swrad, & 
                                                    lakestate%conc(i_t2,itype,i),lakestate%water(itype,i)*lakearea*0.001*fastoutpart, & 
                                                    lakearea*fastoutpart*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    genpar(m_limt2exch),freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CALL calculate_watersurface_heatbalance(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%concslow(i_t2,itype,i)),swrad, & 
                                                    lakestate%concslow(i_t2,itype,i),lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart), & 
                                                    lakearea*(1-fastoutpart)*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    genpar(m_limt2exch),freezeuparea2,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
          CASE(1)
            CALL calculate_T2_transfer(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%conc(i_t2,itype,i)),lakestate%conc(i_t2,itype,i), & 
                                       lakestate%water(itype,i)*lakearea*0.001*fastoutpart,lakearea*fastoutpart*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                                    freezeuparea,genpar(m_liceTf))
            CALL calculate_T2_transfer(temp-1.*(lakestate%uppertemp(itype,i) - lakestate%concslow(i_t2,itype,i)),lakestate%concslow(i_t2,itype,i), & 
                                       lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart),lakearea*(1-fastoutpart)*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                                    freezeuparea2,genpar(m_liceTf))
          END SELECT
            
          !Check freezeup conditions, indicated by relative freezeuparea
          IF(freezeuparea.GT.0. .OR. freezeuparea2.GT.0.)THEN
            !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
            !re-scale to a fraction of the entire waterarea:
            freezeuparea = (freezeuparea*fastoutpart+freezeuparea2*(1-fastoutpart)) * icefreefraction
            freezeup = 1
          ENDIF

          !->re-calculate averaged meantemp and upper temp, then check for freezeup conditions on the new uppertemp as well
          meantemp = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / meanwater
          lakestate%uppertemp(itype,i) = (meantemp * meanwater * 0.001 - (meanwater * 0.001 - epidepth) * lakestate%lowertemp(itype,i)) / epidepth
          IF(lakestate%uppertemp(itype,i).LT.genpar(m_liceTf))THEN
            freezeup = 1
            !!estimate a freezuparea (reduction in the open water surface area) so that the result of the surface heat balance is equal to the freezing point , CP160620 test
            !freezeuparea2 = MAX(0.,MIN(1.,(genpar(m_liceTf) - lakestate%uppertemp(itype,i))/(olduppertemp-lakestate%uppertemp(itype,i))))
            !IF(freezeuparea==0)THEN
            !  freezeuparea = freezeuparea2 * icefreefraction
            !ELSE
            !  freezeuparea = MAX(0.,MIN(1.,freezeuparea+freezeuparea2*(icefreefraction-freezeuparea)))
            !ENDIF
            lakestate%uppertemp(itype,i) = genpar(m_liceTf) !recalculated later again
            !lakestate%lowertemp(itype,i) = (meantemp * meanwater * 0.001 - genpar(m_liceTf)*epidepth)/(meanwater * 0.001 - epidepth)
          ENDIF           
        ENDIF 
      ELSE
        !Shallow lake without thermal stratification

        !Introducing fractional ice cover to get smoother transition over the freezing point
        icefreefraction = 1. - frozenstate%lakeicecov(itype,i)
        IF(icefreefraction.GT.0.)THEN
          ! optional models  (will be reduced to one option after som initial testing for EHYPE3.0 and SHYPE2012)
          SELECT CASE(modeloption(p_lakeriverice))
          CASE(2) ! new model based on Piccolroaz et al 2013
            CALL calculate_watersurface_heatbalance(temp,swrad, & 
                                                    lakestate%conc(i_t2,itype,i),lakestate%water(itype,i)*lakearea*0.001*fastoutpart, & 
                                                    lakearea*fastoutpart*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    genpar(m_limt2exch),freezeuparea,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
            CALL calculate_watersurface_heatbalance(temp,swrad, & 
                                                    lakestate%concslow(i_t2,itype,i),lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart), & 
                                                    lakearea*(1-fastoutpart)*icefreefraction,genpar(m_tcflake),genpar(m_scflake),genpar(m_ccflake),genpar(m_lcflake), & 
                                                    genpar(m_limt2exch),freezeuparea2,genpar(m_liceTf),genpar(m_stbcorr1),genpar(m_stbcorr2),genpar(m_stbcorr3))     
          CASE(1)
            CALL calculate_T2_transfer(temp,lakestate%conc(i_t2,itype,i),lakestate%water(itype,i)*lakearea*0.001*fastoutpart, & 
                                       lakearea*fastoutpart*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                       freezeuparea,genpar(m_liceTf))
            CALL calculate_T2_transfer(temp,lakestate%concslow(i_t2,itype,i),lakestate%slowwater(itype,i)*lakearea*0.001*(1-fastoutpart), & 
                                       lakearea*(1-fastoutpart)*icefreefraction,genpar(m_t2trlake)*t2transfcorr, & 
                                       freezeuparea2,genpar(m_liceTf))
          END SELECT 
         
          !Check freezeup conditions, indicated by relative freezeuparea
          IF(freezeuparea.GT.0. .OR. freezeuparea2.GT.0.)THEN
            !freezup area is the fraction of previously unfrozen area (waterarea*icefreefraction), where new ice formation is triggered
            !re-scale to a fraction of the entire waterarea:
            freezeuparea = (freezeuparea*fastoutpart+freezeuparea2*(1-fastoutpart)) * icefreefraction
            freezeup = 1
          ENDIF
       
          !assign new weighted average temperature for the lowertemp variable
          lakestate%lowertemp(itype,i) = (fastoutpart * lakestate%conc(i_t2,itype,i) * lakestate%water(itype,i) + (1.-fastoutpart) * lakestate%concslow(i_t2,itype,i) * lakestate%slowwater(itype,i)) / &
                              (lakestate%slowwater(itype,i)*(1-fastoutpart)+lakestate%water(itype,i)*fastoutpart)
           
          !set uppertemp to the same, for the output
          lakestate%uppertemp(itype,i) = lakestate%lowertemp(itype,i) 
        ENDIF
      ENDIF
    ENDIF
          
    !2: Upper->Lower Lake T2 exchange
                    
    IF(epilimnion)THEN 
      !Autumn circulation
      IF(lakestate%uppertemp(itype,i)< lakestate%lowertemp(itype,i) .AND. lakestate%uppertemp(itype,i) > 3.95)THEN  !autumn circulation
        lakestate%lowertemp(itype,i) = meantemp
        lakestate%uppertemp(itype,i) = lakestate%lowertemp(itype,i)
      ELSE
      !Spring circulation
        IF(lakestate%uppertemp(itype,i)> lakestate%lowertemp(itype,i) .AND. lakestate%uppertemp(itype,i) < 3.95)THEN  !spring circulation  
          lakestate%lowertemp(itype,i) = meantemp
          lakestate%uppertemp(itype,i) = lakestate%lowertemp(itype,i)
        ELSE
          !Startification; Heat transfer between upper and lower (new function)
          CALL calculate_T2_transfer_upper2lower(lakestate%uppertemp(itype,i),lakestate%lowertemp(itype,i),epidepth*lakearea, & 
                                                (meanwater*0.001-epidepth)*lakearea,lakearea,genpar(m_upper2deep)) 
        ENDIF
      ENDIF
    ENDIF
    
    !Assign lake surface temperature if icefree conditions
    IF((icefreefraction-freezeuparea).GT.0.)lakesurft(itype) = lakestate%uppertemp(itype,i)
    
  END SUBROUTINE T2_processes_in_lake
  
  !>Calculate temperature(T2) "concentration" in lake/river precipitation
  !>due to ice presens
  !Changes the default T2 concentration, set for class
  !Modified to use fractional ice cover
  !----------------------------------------------------------
   SUBROUTINE add_T2_concentration_in_precipitation_on_water(prec,temp,snowfall,rainfall,watertemp,cprec,icecover)
  
    REAL, INTENT(IN)      :: prec       !<precipitation
    REAL, INTENT(IN)      :: temp       !<air temperature
    REAL, INTENT(IN)      :: snowfall   !<snow fall
    REAL, INTENT(IN)      :: rainfall   !<rain fall
    REAL, INTENT(IN)      :: watertemp  !<temperature of water
    REAL, INTENT(INOUT)   :: cprec      !<T2 concentration of precipitation
    REAL, INTENT(IN)      :: icecover   !<ice cover
    
    !This is now much more straight forward, using the fractional ice cover:
    ! Rainfall has always cprec = airtemp (but not lower than freezing point)
    ! Snowfall on the ice-free fraction has cprec = latentheat of freezing + sensible heat content
    ! Snowfall in the ice-covered fraction has cprec = laketemp
    IF(prec.GT.0)THEN
      !Rainfall temperature = max(0,air temp)
      cprec = rainfall * MAX(0.0,temp)
      
      !Snowfalltemp on ice   = watertemp (temporary), negative latent heat is added later when snow is melting in the ice routine
      !Snowfalltemp on water = airtemp + negative latent heat, taking into account diff spec.heat of ice and water
      cprec = cprec + snowfall * (watertemp * icecover + (MIN(temp,0.0)*2.1/4.2 - 335./4.2)*(1.-icecover))

      !Weighting by total precipitation
      cprec = cprec/prec  
    ELSE
      cprec = 0.
    ENDIF
    
  END SUBROUTINE add_T2_concentration_in_precipitation_on_water
  
!>Subroutine to calculate growth of ice on lakes and rivers (only after freezeup has been identified)
! Developers: David Gustafsson(SMHI)
!
! Model is largely based on the review of thermodynamical ice models by Leppäranta (1993)
!
!    Ice growth: "freezing-degree-day"
!    Ice melt:   "positive degree-days" with constant 0.1-0.5 degC day/cm
!
! Snow on ice is considered, however, snowaccumulation and snowmelt is supposed to be calcluated outside of this routine:
!    the routine is calculating freezing of slush ice in case the snowmass is large enough to submerge the ice surface.
!
! Input/output to the model is icedepths, snowmass, snowdepth, and mass rate of lakewater and snowwater transformation from ice
!
! Model was calibrated with Swedish lake ice depth and river ice depth data for the North Hydrology project (Gustafsson et al, 2012)
!---------------------------------------
  SUBROUTINE calculate_icedepth(tsurf,iced,biced,snowm,snowd,tair,dlakewidt,dsnowdt,ifreezeup,ibreakup,tf,kika,kexp,pm)

    !Argument declarations
    REAL, INTENT(INOUT) :: tsurf       !<lake surface temperature, when the lake is ice and/or snowcovered, Tsurf is back calculated from ice growth, unless its melting, then its set to 0
    REAL, INTENT(INOUT) :: iced        !<ice depth, cm (black ice + snowice)
    REAL, INTENT(INOUT) :: biced       !<black ice, cm
    REAL, INTENT(INOUT) :: snowm       !<snowmass, mm
    REAL, INTENT(INOUT) :: snowd       !<snowdepth, cm
    REAL, INTENT(IN)    :: tair        !<air temperature, C
    REAL, INTENT(OUT)   :: dlakewidt   !<transformation of lake ice to lake water (positive direction from ice to water)
    REAL, INTENT(OUT)   :: dsnowdt     !<transformation of snow to lake ice       (positive direction from ice to snow)
    INTEGER, INTENT(IN) :: ifreezeup   !<freeze-up day flag (1=yes, 0=no)
    INTEGER, INTENT(OUT):: ibreakup    !<break-up day flag  (1=yes, 0=no)
    REAL, INTENT(IN)    :: tf          !<tf (~0.)  , freezing point temperature of the lake/river water, °C
    REAL, INTENT(IN)    :: kika        !<kika(~10.), ratio between thermal conductivity of ice and heat exchange coef in air
    REAL, INTENT(IN)    :: kexp        !<kiks(~10.), as above but for snow, actually dependent on snow density, but we use a fixed value
    REAL, INTENT(IN)    :: pm          !<pm (~0.5) , degree-day melt factor for ice, cm/°C
    
    !Local variables
    REAL :: slushd    !slush depth,     cm
    REAL :: siced     !snow ice depth, cm
    REAL :: dsnow     !snow density, g/cm3
    REAL :: S         !freezing degree days, Cday
    REAL :: M         !melting degree days, Cday
    REAL :: dHsicedt  !snowice growth, cm/day, potential
    REAL :: dHsicedt2 !snowice growth, cm/day, actual
    REAL :: dHicedt   !blackice growth, cm/day, actual
    REAL :: oldiced   !icedepth at start of calculaton, cm
    !Parameters, calculated in the code from the input parameters
    REAL :: ki       ! ki, thermal conductivity of ice, W/m/C, ki ~= 1.9 (see Leppäranta-93)
    REAL :: ai       ! a, degreeday factor, see Leppäranta(1993)~=3.3 if ki=1.9 W/m/C
    REAL :: ka       ! ka = ki/kika, heat exchange in air
    REAL :: ks       ! ks = ki*(rhosnow/rhoice)^ksexp, thermal conductivity in snow
    REAL :: kiks
               
    !Physical Constants
    REAL, PARAMETER :: L = 335.      !latent heat of freezing, J/g
    REAL, PARAMETER :: dice = 0.917  !density of ice, g/cm3
    REAL, PARAMETER :: mm2cm = 0.1, cm2mm = 10. !parameters for transformation from mm to cm
    
    !Conversion of some parameters
    ki = 2.2                ! standard value 2.2 for fresh water black ice
    ki = ki * 86400. / 100. ! (W/m/oC) -> (J/d/cm/C)
    ai  = (2*ki/dice/L)**0.5      
    ka = ki/kika    ! kika and kexp can be calibration parameters if needed

    !Initialization of some variables
    dHicedt   = 0.
    dHsicedt  = 0.
    dHsicedt2 = 0.
    dlakewidt = 0.
    dsnowdt   = 0.
    oldiced = iced
    siced = iced - biced
    
    !If there is old ice or if freeze-up condition has been met, calculate ice growth and ice melt
    IF(iced.gt.0. .OR. ifreezeup.EQ.1)THEN ! ifreezeup eq. to tsurf < 0...

      !Freezing and melting degree days (actually, HYPE is using a threshold temperature for snow melt)
      S = AMAX1(0.,-Tair)   ! freezing degree days
      M = AMAX1(0.,Tair)     ! melting degree days

      !Accumulation and Melt of snow on ice is treated outside of this function
 
      !Snow density
      IF(snowm.GT.0.0 .AND. snowd.GT.0.0)THEN
        dsnow = snowm * mm2cm / snowd
        ks = ki * (dsnow/dice)**kexp  
      ELSE
        dsnow = 0.0
        ks = ki * (0.1/dice)**kexp  
      ENDIF
      kiks = ki/ks
        
      !Ice growth
      IF(snowm*mm2cm.GT.iced *(1.-dice))THEN
        !Submerged snow on ice, snowmass exceeds floating capacity of the ice
        !slush depth [cm] above ice surface (depends on snow
        !density, snow mass, and ice mass (assuming no capillary rise
        !in snow), limited by snow depth (check density if there is problem):
        slushd = (snowm*mm2cm - iced * (1.-dice))/(dsnow/dice)
        IF(slushd.gt.snowd)THEN
          WRITE(6,*) 'WARNING: slushdepth > snowdepth. slushdepth, snowdepth, dsnow:',slushd,snowd,dsnow
          slushd = snowd
        ENDIF

        ! Snow-ice growth (d(Hsi)/dt), see Leppäranta(1993), eq 21
        IF(Tair.LT.Tf)THEN
          ! height change, of the snow ice, limited by the slush depth
          dHsicedt  = ks * (Tf - Tair)/(snowd+kika)/(dice*L*(1.-dsnow/dice))
          dHsicedt2 = amin1(slushd,dHsicedt)  ! only valid for daily time steps
          
          ! update surface temperature
          tsurf = dHsicedt*(1.-dsnow/dice)*L*dice / ka + Tair
          
          !update ice depths and masses:
          snowd  = amax1(0.,snowd - dHsicedt2)          ! snow depth, cm
          snowm  = amax1(0.,snowm - dHsicedt2 * dsnow * cm2mm)  ! snowmass, mm

          siced  = amax1(0.,siced + dHsicedt2)          ! snow ice depth, cm
          slushd = amax1(0.,slushd - dHsicedt2)         ! slush depth, cm
      
          iced   = amax1(0.,biced + siced)             ! total ice depth, cm

          !how much lake water (mm) and snow mass (mm) is transformed to snow-ice?
          dlakewidt  = dlakewidt - dHsicedt2 * (1.-dsnow/dice) * cm2mm
          dsnowdt    = dsnowdt   - dHsicedt2 * dsnow * cm2mm

          ! if the potential snow-ice growth was larger than the
          ! slushdepth, it means that we have additional heat loss to
          ! freeze also the black ice, which could be used to calculate black ice
          ! growth at this point:
          
        ENDIF
      ELSE
        ! ICE SURFACE ABOVE WATER SURFACE, AND WE MAY ESTIMATE BLACK ICE GROWTH 
        slushd = 0.
        ! (black) ice growth, including insulation of snow on ice (see Leppäranta(1983), dHdt = 0.5*a^2*S/(H + ki/ka + kiks * h)
        dHicedt = 0.5*ai**2 * S /(iced + kika + snowd*kiks)
        iced = iced + dHicedt
        
        ! update surface temperature
        tsurf = dHicedt * L * dice / ka + Tair
        
        ! we do the calculation for the total ice depth (then separate
        ! snow ice from clack ice)
        biced = iced - siced
        
        !how much lake water (mm) and snow mass (mm) is transformed to snow-ice?
        dlakewidt  = dlakewidt - dHicedt * dice * cm2mm
            
      ENDIF

      ! ICE MELT, simple degree day (if there is ice, if there is no snow, and if there is positive degree days)
      IF(M.GT.0. .AND. iced.GT.0. .AND. snowd.LE.0.)THEN
        dHicedt = - AMIN1(iced,M*pm)         ! pm, degree day factor [cm/C/day]
        iced    = AMAX1(0.,iced + dHicedt)   ! total ice depth [cm]
        siced   = AMAX1(0.,siced + dHicedt); ! snow ice [cm], is melted before the black iace
        biced   = AMAX1(0.,iced-siced);      ! black ice [cm]
        
        ! how much lake water is generated?
        dlakewidt  = dlakewidt - dHicedt * dice * cm2mm

        ! set surface temperature to 0
        tsurf = 0.0
      ENDIF
        
      ! BREAK UP DAY
      IF(iced.LE.0. .AND. oldiced .GT. 0.)THEN
        ibreakup = 1
        iced = 0.
        snowd = 0.
        biced = 0.
        snowm=0.
        slushd=0.
        siced = 0.
        tsurf = 0.0
      ENDIF
    ENDIF
      
  END SUBROUTINE calculate_icedepth

!>Subroutine to calculate transfer of heat from air to water
!---------------------------------------------------------------
  SUBROUTINE calculate_T2_transfer(airtemp,watertemp,watervol,waterarea,  &
                                   T2transfer,freezeuparea,freezingpoint)

    USE MODVAR, ONLY: cwater,seconds_per_timestep

    !Argument declaration
    REAL, INTENT(IN)    :: airtemp       !<air temperature (deg Celsius)
    REAL, INTENT(INOUT) :: watertemp     !<water temperature (deg Celsius)
    REAL, INTENT(IN)    :: watervol      !<surface water volume (m3 or mm)
    REAL, INTENT(IN)    :: waterarea     !<surface water area (m2)
    REAL, INTENT(IN)    :: T2transfer    !<heat transfer parmeter from air to water (J/m2/s/deg)
    REAL, INTENT(OUT)   :: freezeuparea  !fractional area were ice formation is trigered (fraction, 0-1)
    REAL, INTENT(IN)    :: freezingpoint !freezingpoint temperature, deg C

    !Local variable declarations
    REAL t2_transf              !T2 transfer    
    REAL density
    REAL heatcapacity, thermcond
      
    density = 1000.
    heatcapacity = cwater * density * 1000.
    thermcond = T2transfer * seconds_per_timestep   !J/m2/deg/timestep
    freezeuparea = 0.
      
    IF(airtemp > watertemp)THEN
      t2_transf = MIN((airtemp - watertemp) * watervol * heatcapacity,(airtemp - watertemp)* waterarea * thermcond)
    ELSE
      t2_transf = MAX((airtemp - watertemp) * watervol * heatcapacity,(airtemp - watertemp)* waterarea * thermcond)
    ENDIF
    IF(watervol>0.)THEN
      !evaluate ice formation conditions (new temperature<freezing point)
      IF((watertemp * watervol * heatcapacity + t2_transf) / (watervol * heatcapacity).LT.freezingpoint)THEN
        !estimate a freezup area (reduction in the open water surface area) so that the result of the surface heat balance is equal to the freezing point 
        freezeuparea=  max(0.,min(1.,1. - (freezingpoint * (watervol * heatcapacity) - watertemp * watervol * heatcapacity)/t2_transf))
        watertemp = freezingpoint
      ELSE
        !calculate new temperature, water volume must be in m3!
        watertemp = (watertemp * watervol * heatcapacity + t2_transf) / (watervol * heatcapacity)
      ENDIF
    ENDIF
 
  END SUBROUTINE calculate_T2_transfer
      
!>Subroutine to calculate transfer of heat(temperature) between upper and lower layer in lakes 
!---------------------------------------
  SUBROUTINE calculate_T2_transfer_upper2lower(uppertemp,lowertemp,uppervol,lowervol,waterarea,T2transfer) 

    USE MODVAR, ONLY: cwater,seconds_per_timestep

    !Argument declaration
    REAL, INTENT(INOUT) :: uppertemp     !<upper water temperature (deg Celsius)
    REAL, INTENT(INOUT) :: lowertemp     !<lower water temperature (deg Celsius)
    REAL, INTENT(IN)    :: uppervol      !<upepr layer water volume (m3)
    REAL, INTENT(IN)    :: lowervol      !<upepr layer water volume (m3)
    REAL, INTENT(IN)    :: waterarea     !<surface water area (m2)
    REAL, INTENT(IN)    :: T2transfer    !<heat transfer parmeter from air to water (J/m2/s/deg)

    !Local variable declarations
    REAL t2_transf              !T2 transfer    
    REAL density
    REAL heatcapacity, thermcond
    REAL equiltemp
    
    density = 1000.
    heatcapacity = cwater * density * 1000.
    thermcond = T2transfer * seconds_per_timestep   !J/m2/deg/timestep
    
    !Calculate equilibrium temperature, when heat is evenly distributed
    IF((uppervol+lowervol).GT.0.)THEN
      equiltemp = (uppertemp*uppervol + lowertemp*lowervol)/(uppervol+lowervol)
    
      !calculate heatflow and update temperatures, depending on initial gradient:
      IF(uppertemp > lowertemp)THEN
        !heat flow from upper to lower
        t2_transf = (uppertemp - lowertemp)* waterarea * thermcond
        !Upper and lower temperatures, limited by equilibrium temperature
        uppertemp = MAX(equiltemp,(uppertemp * uppervol * heatcapacity - t2_transf) / (uppervol * heatcapacity))
        lowertemp = MIN(equiltemp,(lowertemp * lowervol * heatcapacity + t2_transf) / (lowervol * heatcapacity))
      ELSE
        !heat flow from lower to upper
        t2_transf = (lowertemp - uppertemp)* waterarea * thermcond
        !Upper and lower temperatures, limited by equilibrium temperature
        uppertemp = MIN(equiltemp,(uppertemp * uppervol * heatcapacity + t2_transf) / (uppervol * heatcapacity))
        lowertemp = MAX(equiltemp,(lowertemp * lowervol * heatcapacity - t2_transf) / (lowervol * heatcapacity))
      ENDIF
    ELSE
      uppertemp=0.
      lowertemp=0.
    ENDIF
 
  END SUBROUTINE calculate_T2_transfer_upper2lower

!>\brief Subroutine to calculate transfer of heat from air to water including a solar radiation term and a residual term.
!>
!>The routine is based on the model sugested by Piccolroaz et al (2013), with modifications 
!>to use real (or estimated) shortwave radiation. 
!>Partly ice covered situations can be taken into account by reducing the input waterarea
!>If the heat balance is negative enough to lower temperature below freezing, a reduction in
!>the surface area is estimated, which shows at how large area the ice is forming.
!
!TODO: make T2 subroutines work for other timestep than day
!---------------------------------------
  SUBROUTINE calculate_watersurface_heatbalance(airtemp,swrad,watertemp,watervol, & 
                  waterarea,tempcoef,radcoef,constcoef,lincoef,limt2exch, &
                  freezeuparea,freezingpoint,stabpar1,stabpar2,stabpar3)

  USE MODVAR, ONLY: cwater,seconds_per_timestep
    
    !Argument declaration
    REAL, INTENT(IN)    :: airtemp       !<air temperature (deg Celsius)
    REAL, INTENT(IN)    :: swrad         !<shortwave radiation (MJ/m2/day)
    REAL, INTENT(INOUT) :: watertemp     !<water temperature (deg Celsius)
    REAL, INTENT(IN)    :: watervol      !<water volume (m3)
    REAL, INTENT(IN)    :: waterarea     !<water surface area (m2)
    REAL, INTENT(IN)    :: tempcoef      !<heat transfer parameter from air to water (J/m2/s/deg)
    REAL, INTENT(IN)    :: radcoef       !<heat transfer parameter from radiation to water (fraction, 0-1)
    REAL, INTENT(IN)    :: constcoef     !<heat transfer parameter, constant residual term (J/m2/s)
    REAL, INTENT(IN)    :: lincoef       !<heat transfer parameter, linear residualterm (J/m2/s/deg)
    REAL, INTENT(IN)    :: limt2exch     !<heat transfer parameter, limit depth for only temperature exchange (m)
    REAL, INTENT(OUT)   :: freezeuparea  !<fractional area were ice formation is trigered (fraction, 0-1)
    REAL, INTENT(IN)    :: freezingpoint !<freezingpoint temperature, deg C
    REAL, INTENT(IN)    :: stabpar1      !<Stability parameter, affects both heating and cooling. No correction if set to zero
    REAL, INTENT(IN)    :: stabpar2      !<Stability parameter, affects cooling. No correction if set to zero
    REAL, INTENT(IN)    :: stabpar3      !<Stability parameter, affects heating. No correction if set to zero
       
    !Local variable declarations
    REAL netheat                  !Net heat flux to the water (J/timestep)    
    REAL density                  !Water density (kg/m3)
    REAL heatcapacity             !heat capacity of water (J/m3/deg)
    REAL tempdiff                 !Temperature difference
    REAL stabfunction             !Stability correction function

    REAL, PARAMETER :: seconds_per_day = 86400.  

    !> \b Algorithm \n
    density      = 1000.                    ! kg/m3, density of water
    heatcapacity = cwater * density * 1000. ! J/m3/deg  [kJ/kg/deg * kg/m3 * 1/k]
    freezeuparea = 0.
                          
    IF(watervol/waterarea>limt2exch)THEN
      IF(watervol>0.)THEN   !make calculation only if the water has a volume
        netheat = 0.   !initialize the net heat flux, J/timestep
      
        !>Calculate stability correction for heat exchange between air and water
        tempdiff = airtemp - watertemp
        IF(tempdiff>0.) THEN
          stabfunction = 1./(1. + stabpar1 * tempdiff)**stabpar3
        ELSE
          stabfunction = 1./(1. - stabpar1 * tempdiff)**(-stabpar2)
        ENDIF
        !>Add the air temperature term
        IF(airtemp > watertemp)THEN
          netheat = netheat + MIN((airtemp - watertemp) * watervol * heatcapacity, stabfunction * (airtemp - watertemp)* waterarea * tempcoef * seconds_per_timestep) !J/timestep
        ELSE
          netheat = netheat + MAX((airtemp - watertemp) * watervol * heatcapacity, stabfunction * (airtemp - watertemp)* waterarea * tempcoef * seconds_per_timestep) !J/timestep
        ENDIF
        !>Add the radiation term, MJ/m2/day => J/m2/s and then multiplied with timestep in s.
        netheat = netheat + 1.E6 * swrad /seconds_per_day * waterarea * radcoef * seconds_per_timestep
        !>Add the residual term, same units as temperature equation
        netheat = netheat + (watertemp*lincoef + constcoef) * waterarea * tempcoef * seconds_per_timestep
      
        !>Evaluate ice formation conditions (new temperature<freezing point) and calculate new water temperature
        IF((watertemp * watervol * heatcapacity + netheat) / (watervol * heatcapacity).LT.freezingpoint)THEN
          !estimate a freezup area (reduction in the open water surface area) so that the result of the surface heat balance is equal to the freezing point 
          freezeuparea = MAX(0.,MIN(1.,1. - (freezingpoint * (watervol * heatcapacity) - watertemp * watervol * heatcapacity)/netheat))
          watertemp = freezingpoint
        ELSE
          !calculate new temperature, water volume must be in m3!
          watertemp = (watertemp * watervol * heatcapacity + netheat) / (watervol * heatcapacity)
        ENDIF
      ENDIF
    ELSE
      !Use only temperature heat exchange for shallow waters
      CALL calculate_T2_transfer(airtemp,watertemp,watervol,waterarea,tempcoef,freezeuparea,freezingpoint)
    ENDIF
 
  END SUBROUTINE calculate_watersurface_heatbalance

!>\brief Calculate water level and flooded area from water volume for a sloping floodplain
!-----------------------------------------------------------------------------------------
SUBROUTINE calculate_floodplain_waterlevel(vol,amax,ymax,y,a)

  !Argument declarations
  REAL, INTENT(IN)  :: vol  !<water volume in floodplain [m3]
  REAL, INTENT(IN)  :: amax !<area at maximum areal extent [m2]
  REAL, INTENT(IN)  :: ymax !<water level at maximum areal extent [m]
  REAL, INTENT(OUT) :: y    !<calculated water level at volume volm3 [m]
  REAL, INTENT(OUT) :: a    !<calculated area [m2]
  
  !Local variables
  REAL volmax             !water volume when a=amax and y=ymax [m3]

  !volume at maximum extent
  volmax    = amax * ymax * 0.5

  !Calculate water depth and area
  IF(vol <= volmax)THEN 
    !case 1: water volume smaller than volume at maximum extent
    y = SQRT(vol * ymax * 2. / amax )
    a = y * amax / ymax
  ELSE
    !case 2: water volume larger than volume at maximum extent
    y = ymax + (vol-volmax)/amax
    a = amax
  ENDIF 

END SUBROUTINE calculate_floodplain_waterlevel

!>\brief Calculate water volume and flooded area from water level for a sloping floodplain
!-----------------------------------------------------------------------------------------
SUBROUTINE calculate_floodplain_volume(y,amax,ymax,vol,a)

  !Argument declarations
  REAL, INTENT(IN)  :: y    !<water level [m]
  REAL, INTENT(IN)  :: ymax !<water level at maximum areal extent [m]
  REAL, INTENT(IN)  :: amax !<area at maximum areal extent [m2]
  REAL, INTENT(OUT) :: vol  !<water volume in floodplain [m3]
  REAL, INTENT(OUT) :: a    !<calculated area [m2]

  !Calculate water volume and area
  IF(y <= ymax)THEN 
    !case 1: water level below level at maximum extent
    vol = y*y*amax/(ymax*2) 
    a = y * amax / ymax
  ELSE
    !case 2: water level above level at maximum extent
    vol = amax * ymax * 0.5 !volume at ymax
    vol = vol + (y-ymax)*amax !adding volume in water above ymax
    a = amax
  ENDIF 

END SUBROUTINE calculate_floodplain_volume


!>\brief Calculate equilibrium water level in river(or lake) and flooded area 
!         for a sloping floodplain with maximum extent amax and corresponding level ymax.
!          - the equilibrium level is given in the reference system for the contributing water storage
!-----------------------------------------------------------------------------------------------------
SUBROUTINE calculate_floodplain_equilibriumlevel(volp,volr,flr,flp,ar,amax,ymax,yeq,r2p)

  !Argument declarations
  REAL, INTENT(IN)  :: volp   !<water volume in floodplain [m3]
  REAL, INTENT(IN)  :: volr   !<water volume in river (or lake) [m3]
  REAL, INTENT(IN)  :: flr    !<flooding level for the river (or lake) [m]
  REAL, INTENT(IN)  :: flp    !<flooding level for the floodplain [m]
  REAL, INTENT(IN)  :: ar     !<area of the river (or lake) [m2]
  REAL, INTENT(IN)  :: amax   !<area at maximum areal extent [m2]
  REAL, INTENT(IN)  :: ymax   !<water level at maximum areal extent [m]
  REAL, INTENT(OUT) :: yeq    !<equilibrium water level [m]
  INTEGER, INTENT(OUT) :: r2p !<flow direction flag, 1=river2plain, 2=plain2river,0=no flow
  
  !Local variables
  REAL :: voltot ! [m3] total volume to distribute
  REAL :: yr
  REAL :: yp
  REAL :: ap
  
  !water level in river (or lake)
  yr = volr/ar
  
  !water level in floodplain
  CALL calculate_floodplain_waterlevel(volp,amax,ymax,yp,ap)
  
  !make sure at least one water surfac is above its flooding level, otherwise nothing to do
  IF(yr<=flr .AND. yp<= flp)THEN
    !Nothing to do
    yeq = -9999
    r2p = 0
    RETURN
  ELSE !Continue
  
    !determine which water surface is higher than the other in a common reference system
    IF((yr-flr)>=(yp-flp))THEN
      r2p = 1
    ELSE
      r2p = 2
    ENDIF
    
    !total water volume to distribute on river (or lake) and floodplain
    voltot = volp + volr
    
    !Exclude/replace river water volume below the lowest common level, in other words: 
    ! modify the system into a solvable problem, where the river and floodplain bottom are
    ! at the same level.
    voltot = voltot - (flr-flp)*ar
    
    !If the remaining volume is zero or negative, it means that the equilibriuim level is 
    ! below the lowest common level, thus it is set to the flooding level of the water 
    ! storage with the highest current surface (river or plain)
    IF(voltot<=0.)THEN
      !Select appropriate flooding level as equilibrium level, and return
      SELECT CASE(r2p)
        CASE(1) !river to plain
          yeq = flr
        CASE(2) !plain to river
          yeq = flp
        CASE DEFAULT
      ENDSELECT
      RETURN
    ELSE !continue
      
      !solve second degree equation for the equilibrium water level - with regard to the common bottom level
      yeq = ( -ar + sqrt(ar*ar + 2. * amax * voltot / ymax) ) / (amax / ymax)

      !adjust the equilibrium level to common floodlevel reference, maximize to 0 
      yeq = MAX(0.,yeq - flp)
      
      !finally, adjust equilibrium level to appropriate reference system, depending on flow direction:
      SELECT CASE(r2p)
        CASE(1) !river to plain
          yeq = flr+yeq
        CASE(2) !plain to river
          yeq = flp+yeq
        CASE DEFAULT
      ENDSELECT      
    ENDIF
  ENDIF
  RETURN
  
END SUBROUTINE calculate_floodplain_equilibriumlevel

  !>\brief Calculate interflow between water body and floodplain
  !-------------------------------------------------------------
  SUBROUTINE calculate_waterbody_floodplain_interflow(i,fpamax,warea,ifpar,volp,concp,volw,concw,fpdepth,fpdegree,interflow)

    USE MODVAR, ONLY : numsubstances, &
                       basin
    
    !Argument declarations
    INTEGER,INTENT(IN) :: i         !<index of subbasin
    REAL,INTENT(IN)    :: fpamax    !<maximum area of the flood plain [m2] 
    REAL,INTENT(IN)    :: warea     !<(maximum=constant) area of the water body [m2] 
    REAL,INTENT(IN)    :: ifpar(5)  !<current interflow parameters (flmrr/floll,flmrp/flolp,fymmr/fymol,rcr2fp/rcl2fp,rcfp2r/rcfp2l
    REAL,INTENT(INOUT) :: volp      !<water volume in floodplain [m3]
    REAL,INTENT(INOUT) :: concp(numsubstances) !<floodplain concentrations
    REAL,INTENT(INOUT) :: volw      !<water volume in water body (river or lake) [m3]
    REAL,INTENT(INOUT) :: concw(numsubstances) !<water body (river or lake) concentration
    REAL,INTENT(OUT)   :: fpdepth   !<flood plain water depth [m]
    REAL,INTENT(OUT)   :: fpdegree  !<flood plain degree of flood [%]
    REAL,INTENT(OUT)   :: interflow !<flow from waterbody to floodplain [m3/timestep] (can be negative)

    !Local variables
    REAL wlmw     !water body water depth [m] 
    REAL wlmfp    !floodplain water depth [m] 
    REAL fparea   !floodplain water surface area [m2]
    REAL wlequil  !equilibrium level [m]
    REAL qmrflood !interflow [m3 (m)]
    REAL qmrfloodc(numsubstances) !concentration of interflow
    REAL voltemp, atemp
    INTEGER fpflowdir !direction of interflow
    INTEGER status    !error status of subroutine calls
    
      interflow = 0.
      
      !water depth in river [m]
      wlmw = volw/warea !water level in main river [m]

      !Calculate current floodplain water depth [m] and water surface area [m2]:
      CALL calculate_floodplain_waterlevel(volp,fpamax,ifpar(3),wlmfp,fparea)
            
      !Flow occurs  when the river and/or the floodplain are higher than
      !their respective flooding thresholds:
      IF(wlmw>ifpar(1) .OR. wlmfp>ifpar(2))THEN
            
        !Get flow direction and equilibrium level:
        CALL calculate_floodplain_equilibriumlevel(volp,volw,     & 
                                                    ifpar(1),ifpar(2),warea,                        &
                                                    fpamax, ifpar(3), &
                                                    wlequil,fpflowdir)
        !wlequil is now defined in the reference system of the river or the floodplain depending on direction
              
        IF(fpflowdir == 1)THEN !flow from waterbody (river or lake) to floodplain
                
          !water level change in river due to flow from river to plain, calculated as
          !a fraction of the potential water level change defined by parameter 
          qmrflood = ifpar(4) * (wlmw - wlequil) * warea
                
          !remove water and substances from waterbody and add it to flood plain
          qmrfloodc(:) = concw
          CALL remove_water(volw,numsubstances,concw,qmrflood,qmrfloodc,status)
          IF(status.NE.0) CALL error_remove_water(errstring(9),basin(i)%subid,i,0)
          CALL add_water(numsubstances,volp,concp,qmrflood,qmrfloodc)
          interflow = qmrflood
              
        ELSEIF(fpflowdir == 2)THEN !flow from floodplain to river
                
          !water level change in plain due to flow from plain to river (m):
          qmrflood = ifpar(5) * (wlmfp - wlequil)
          !Floodplain volume at the new water level:
          CALL calculate_floodplain_volume(wlmfp-qmrflood, fpamax, &
                                            ifpar(3),voltemp,atemp)

          !Calculate flow as change in volume [m3]:
          qmrflood = volp - voltemp
                
          !Remove qmrflood from flood plain and add it to waterbody
          qmrfloodc(:) = concp
          CALL remove_water(volp,numsubstances,concp,qmrflood,qmrfloodc,status)
          IF(status.NE.0) CALL error_remove_water(errstring(10),basin(i)%subid,i,0)
          CALL add_water(numsubstances,volw,concw,qmrflood,qmrfloodc)
          interflow = - qmrflood

        ENDIF
      ELSE
        !No interflow this time:
        !qmrflood = 0.
      ENDIF
          
      !Calculate output variables floodplain water level and water area
      CALL calculate_floodplain_waterlevel(volp,fpamax,ifpar(3),wlmfp,fparea)
      fpdepth = wlmfp
      fpdegree = fparea/fpamax*100.
      

  END SUBROUTINE calculate_waterbody_floodplain_interflow

        
END MODULE SURFACEWATER_PROCESSES
