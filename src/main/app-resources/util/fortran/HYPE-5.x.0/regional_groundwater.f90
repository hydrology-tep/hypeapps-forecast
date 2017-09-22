!> \file regional_groundwater.f90
!> Contains module regional_groundwater_module.

!>\brief Regional groundwater calculations in HYPE
!>
!>Two options for regional groundwater: 
!>1-instant transport to outlet of same or other subbasin
!>2-aquifer with delay of groundwater flow
!!
!>Model 1: Regional groundwater flow can move water (and concentrations) from lowest soil
!>layer to the outlet lake of the subbasin or to the lowest soil layer in another
!>subbasin. It uses the surface water subbasin coupling if no separate coupling 
!>for regional groundwater flow is given.
!>
!>Model 2: The aquifer modelled is a sub-surface, hydrologically connected, unconfined storage.
!>It is treated as a single large store f water and reactor for substances.
MODULE REGIONAL_GROUNDWATER_MODULE


  !Copyright 2012-2017 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

  USE STATETYPE_MODULE
  USE GENERAL_WATER_CONCENTRATION
  USE MODVAR, ONLY : aquifer,           & !AquiferData
                     path,              & !-"-
                     naquifers,         & !-"-
                     nsub,              &
                     numsubstances,     &
                     missing_value,     &
                     timesteps_per_day, &
                     seconds_per_timestep
  USE NPC_SOIL_PROCESSES, ONLY : soil_denitrification
  !Uses also modvar, hypevariables, npc_soil_processes

  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------
  !Private procedures 
  !----------------------------------------------
  ! calculate_soillayer_groundwaterflow_removal 
  ! set_percolation_concentration
  ! calculate_percolation_delay
  ! calculate_aquifer_outflow
  !----------------------------------------------
  PUBLIC :: initiate_regional_groundwater_flow, &
            initiate_aquifer_state,   &
            initiate_aquifer_model,   &
            calculate_regional_groundwater_flow, &
            calculate_river_groundwaterflow_removal, &
            calculate_river_floodplain_groundwaterflow_removal, &
            calculate_current_flow_from_aquifer, &
            add_regional_groundwater_flow_to_soil, &
            add_regional_groundwater_flow_to_olake, &
            calculate_aquifer,        &
            calculate_delayed_water,  &
            add_aquifer_flow_to_river,  &
            calculate_aquifer_waterlevel

  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(5)  !error message for location of remove_water call
  PARAMETER(errstring = (/'regional groundwater flow from soillayer         ', &
                          'outflow from aquifer (not subid, but aquid)      ', &
                          'groundwater outflow from main river              ', &
                          'groundwater outflow from main river (queue)      ', &
                          'aquifer outflow for irrigation (not subid, aquid)'/))

  !Private variables, global in this module
  SAVE
  REAL, ALLOCATABLE :: regrwarea(:)  !Receiving area of groundwater flow (land area for now) (m2) (nsub)

  REAL, ALLOCATABLE :: grwflow(:)    !Outflow from groundwater (mm*m2/timestep) (nsub)
  REAL, ALLOCATABLE :: cgrwflow(:,:) !Concentration of outflow from groundwater (numsubstances,nsub)

  REAL, ALLOCATABLE :: dcoeff(:)      !Coefficient for recharge delay equation
  REAL, ALLOCATABLE :: addflow(:)     !Flow from aquifer to main rivers (m3/timestep)
  REAL, ALLOCATABLE :: caddflow(:,:)  !Concentration of outflow from aquifer to main rivers

CONTAINS

  !>Initiate regional groundwater flow variables
  !>
  !>\b Consequences Module variables from regional_groundwater_module may be allocated
  !--------------------------------------------------------------------------
  SUBROUTINE initiate_regional_groundwater_flow(n,ns,cb_ilake,cb_olake,cb_lriver,cb_mriver)

    USE MODVAR, ONLY : basin,      &
                       classbasin

    !Argument declarations
    INTEGER, INTENT(IN)   :: n           !<number of subbasin
    INTEGER, INTENT(IN)   :: ns          !<number of substances
    INTEGER, INTENT(IN)   :: cb_ilake    !<slc-index for internal lake
    INTEGER, INTENT(IN)   :: cb_olake    !<slc-index for outlet lake
    INTEGER, INTENT(IN)   :: cb_lriver   !<slc-index for local stream
    INTEGER, INTENT(IN)   :: cb_mriver   !<slc-index for main river

    !Local variables
    REAL lpart(n)    !lake fraction of subbasin area

    !>\b Algorithm \n
    !>Set receiving area of groundwater flow (land area for now)
    IF(.NOT.ALLOCATED(regrwarea)) ALLOCATE(regrwarea(n))
    regrwarea = basin(1:n)%area 
    IF(cb_olake>0 .OR. cb_ilake>0 .OR. cb_lriver>0 .OR. cb_mriver>0)THEN 
       lpart = 0.
       IF(cb_olake>0)  lpart = lpart + classbasin(1:n,cb_olake)%part
       IF(cb_ilake>0)  lpart = lpart + classbasin(1:n,cb_ilake)%part
       IF(cb_lriver>0) lpart = lpart + classbasin(1:n,cb_lriver)%part
       IF(cb_mriver>0) lpart = lpart + classbasin(1:n,cb_mriver)%part
       regrwarea = regrwarea * (1. - lpart)
       WHERE(regrwarea<0)regrwarea=0.   !Safe for all lake subbasin (1-1<0)
    ENDIF
    !>Allocate groundwater flow variables
    IF(.NOT.ALLOCATED(grwflow))  ALLOCATE(grwflow(n))
    IF(.NOT.ALLOCATED(cgrwflow)) ALLOCATE(cgrwflow(ns,n))

  END SUBROUTINE initiate_regional_groundwater_flow

  !>Initiate aquifer state variables
  !--------------------------------------------------------------------------
  SUBROUTINE initiate_aquifer_state(aquiferstate)

    USE MODVAR, ONLY : i_in,i_sp,i_t2
    
    !Argument declarations
    TYPE(aquiferstatetype),INTENT(INOUT):: aquiferstate !<Aquifer states

    !Local variables
    INTEGER ia

    !>\b Algorithm \n
    !>Calculate aquifer average volyme = initial volume
    DO ia=1,naquifers
      aquiferstate%water(ia) = aquifer(ia)%inivol   !m3
      IF(i_in>0) aquiferstate%conc(i_in,ia) = aquifer(ia)%conc_IN   !mg/L
      IF(i_sp>0) aquiferstate%conc(i_sp,ia) = aquifer(ia)%conc_SP   !mg/L
      IF(i_t2>0) aquiferstate%conc(i_t2,ia) = aquifer(ia)%temperature   !degree C
    ENDDO

  END SUBROUTINE initiate_aquifer_state

  !>Initiate regional groundwater flow variables
  !>
  !>\b Consequences Module variables from regional_groundwater_module may be allocated
  !--------------------------------------------------------------------------
  SUBROUTINE initiate_aquifer_model(n,ns,na)

    USE HYPEVARIABLES, ONLY : m_aqdelcorr
    USE MODVAR, ONLY : regpar, &
                       regiondivision
    
    !Argument declarations
    INTEGER, INTENT(IN)   :: n           !<number of subbasin
    INTEGER, INTENT(IN)   :: ns          !<number of substances
    INTEGER, INTENT(IN)   :: na          !<number of aquifers

    !Local variables
    INTEGER ia

    !>\b Algorithm \n
    !>Allocate aquifer flow variables
    IF(.NOT.ALLOCATED(grwflow))  ALLOCATE(grwflow(n))
    IF(.NOT.ALLOCATED(cgrwflow)) ALLOCATE(cgrwflow(ns,n))
    
    !Dessa behöver bättre namn!!
    IF(.NOT.ALLOCATED(addflow)) ALLOCATE(addflow(n))
    IF(.NOT.ALLOCATED(caddflow)) ALLOCATE(caddflow(ns,n))
    
    !>Calculate coefficient for deep percolation delay
    IF(.NOT.(ALLOCATED(dcoeff))) ALLOCATE(dcoeff(na))
    dcoeff = 0. !no delay
    DO ia = 1,na
      IF(aquifer(ia)%percdelay>0. .AND. aquifer(ia)%parregion(regiondivision(m_aqdelcorr))>0)THEN
        dcoeff(ia) = EXP(-1./(aquifer(ia)%percdelay*timesteps_per_day*(1.+regpar(m_aqdelcorr,aquifer(ia)%parregion(regiondivision(m_aqdelcorr))))))
      ENDIF
    ENDDO

  END SUBROUTINE initiate_aquifer_model

  !>\brief Calculate regional groundwater flow from all subbasins and land
  !classes. Removes the flow from the soil
  !!
  !>\b Consequences Module variables grwflow and cgrwflow is set.
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Regional groundwater flow
  !>and Aquifers - Deep percolation)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_regional_groundwater_flow(soilstate,outputflow,classload)

    USE HYPEVARIABLES, ONLY : m_rcgrwst,  &
                              m_rcgrw,  &
                              m_aqpercorr, &
                              wpmm,     &
                              fcmm
    USE MODVAR, ONLY : basin,      &
                       classbasin, &
                       classdata,  &
                       soilthick,  &
                       nsub,       &
                       nclass,     &
                       numsubstances, &
                       maxsoillayers, &
                       slc_ilake,     &
                       slc_olake,     &
                       slc_lriver,    &
                       slc_mriver,    &
                       regiondivision, &
                       genpar,        &
                       regpar,        &
                       soilpar

    !Arument declarations
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    REAL, INTENT(OUT)   :: outputflow(4,nsub)         !<Outflow from groundwater (m3/timestep) (1=all,2-4=per soillayer)
    REAL, INTENT(OUT) :: classload(nclass,numsubstances,nsub) !<Load removed from soil by regional groundwater flow (kg/timestep)

    !Local variables
    INTEGER i,j   !loop-variables i-subbasin, j-class
    INTEGER k     !soil layer
    INTEGER isoil !soil type
    REAL a        !class' area fraction
    REAL rcgrwj   !model parameters
    REAL grwflowj,cgrwflowj(numsubstances)
    REAL soillflow(3,nsub)  !groundwater flow of subbasin accumulated per soillayer

    !>\b Algorithm \n
    !>Initiations
    grwflow  = 0.
    cgrwflow = 0.
    outputflow = 0.
    classload = 0.

    !>Return if no regional groundwater is simulated
    IF(genpar(m_rcgrw)==0 .AND. SUM(soilpar(m_rcgrwst,:))==0) RETURN

    !Initiations
    soillflow = 0.

    !>For every subbasin and land class:
    DO i = 1,nsub
      IF(naquifers==0.OR.path(i)%rechargebasin)THEN
        DO j = 1,nclass
          IF(j==slc_ilake .OR. j==slc_olake .OR.  &   !only land classes
             j==slc_lriver .OR. j==slc_mriver) CYCLE
          a = classbasin(i,j)%part
          IF(a>0)THEN
            isoil = classdata(j)%soil
            IF(soilpar(m_rcgrwst,isoil)>0 .AND. basin(i)%parregion(regiondivision(m_aqpercorr))>0)THEN
              rcgrwj = soilpar(m_rcgrwst,isoil)*(1.+regpar(m_aqpercorr,basin(i)%parregion(regiondivision(m_aqpercorr))))
            ELSE
              rcgrwj = genpar(m_rcgrw)       !default coefficient for regional groundwater flow
            ENDIF
            !>\li Calculate soil groundwater removal from the lowest soil layer 
            IF(soilthick(3,j)>0)THEN
              CALL calculate_soillayer_groundwaterflow_removal(i,j,   &
                   basin(i)%subid,numsubstances,soilstate%water(3,j,i),   &
                   soilstate%conc(:,3,j,i),rcgrwj,fcmm(3,j)+wpmm(3,j),   &
                   grwflowj,cgrwflowj)
              k = 3
            ELSEIF(soilthick(2,j)>0)THEN
              CALL calculate_soillayer_groundwaterflow_removal(i,j,   &
                   basin(i)%subid,numsubstances,soilstate%water(2,j,i),   &
                   soilstate%conc(:,2,j,i),rcgrwj,fcmm(2,j)+wpmm(2,j),   &
                   grwflowj,cgrwflowj)
              k = 2
            ELSE
              CALL calculate_soillayer_groundwaterflow_removal(i,j,   &
                   basin(i)%subid,numsubstances,soilstate%water(1,j,i),   &
                   soilstate%conc(:,1,j,i),rcgrwj,fcmm(1,j)+wpmm(1,j),   &
                   grwflowj,cgrwflowj)
              k = 1
            ENDIF
            !>\li Accumulate to subbasin regional groundwater flow
            classload(j,:,i) = cgrwflowj(:) * grwflowj * a
            grwflow(i) = grwflow(i) + grwflowj * a
            cgrwflow(:,i) = cgrwflow(:,i) + classload(j,:,i)    !sum c*q
            soillflow(k,i) = soillflow(k,i) + grwflowj * a
          ENDIF   !a>0
        ENDDO

        !>Calculate load removed from each subbasin 
        IF(numsubstances>0) classload(:,:,i)=classload(:,:,i) * basin(i)%area * 1.E-6          !Load removed (kg/timestep)
        !Finalise groundwater flow for subbasin, i.e. change units           
        IF(numsubstances>0 .AND. grwflow(i)>0.) cgrwflow(:,i) = cgrwflow(:,i) / grwflow(i)       !conc, mg/L (or other)
        grwflow(i) = grwflow(i) * basin(i)%area        !mm*m2 
        soillflow(:,i) = soillflow(:,i) * basin(i)%area * 1.E-3   !m3/timestep
      ENDIF

    ENDDO
    !Set output argument for printout
    outputflow(1,:) = grwflow * 1.E-3   !m3/timestep
    outputflow(2:4,:) = soillflow

  END SUBROUTINE calculate_regional_groundwater_flow

  !>\brief Calculate regional groundwater flow from a soillayer 
  !>Reduce OC-, ON-, and PP-concentration of flow and remove the water from the soil layer.
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Regional groundwater flow
  !>and Aquifers - Deep percolation)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_soillayer_groundwaterflow_removal(i,j,subid,n,soil,  &
       csoil,rcgrwpar,plantwater,flow,cflow)

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index for current subbasin
    INTEGER, INTENT(IN) :: j          !<index for current class
    INTEGER, INTENT(IN) :: subid      !<subid for current subbasin
    INTEGER, INTENT(IN) :: n          !<array-dimension, number of substances
    REAL, INTENT(INOUT) :: soil       !<current soil water (mm)
    REAL, INTENT(INOUT) :: csoil(n)   !<concentration of soil water
    REAL, INTENT(IN)    :: rcgrwpar   !<parameter for regional groundwater removal
    REAL, INTENT(IN)    :: plantwater !<volume of field capacity (mm)
    REAL, INTENT(INOUT) :: flow       !<calculated regional groundwater flow contribution from this soillayer (mm/timestep)
    REAL, INTENT(INOUT) :: cflow(n)   !<concentration of calculated regional groundwater flow 

    !Local variables
    INTEGER status

    !>\b Algorithm \n
    !Default, no flow
    flow = 0.
    cflow(:) = 0. 

    !>If water above field capacity:
    IF(soil>plantwater)THEN
      !>Calculate regional groundwater flow
      flow = rcgrwpar * (soil-plantwater)
      CALL set_percolation_concentration(n,csoil,cflow) !Reduce the concentration of the flows
      !>Remove flow from soil
      CALL remove_water(soil,n,csoil,flow,cflow,status)
      IF(status.NE.0) CALL error_remove_water(errstring(1),subid,i,j)
    ENDIF

  END SUBROUTINE calculate_soillayer_groundwaterflow_removal

  !>\brief Calculate regional groundwater flow from river. 
  !>Reduce OC-, ON-, and PP-concentration of flow and remove the water from the river.
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_river_groundwaterflow_removal(i,j,subid,n,riverstate,flow)

    USE HYPEVARIABLES, ONLY : m_rcgrwst, &
                              m_rcgrw, &
                              m_aqpercorr, &
                              ttpart, ttstep
   USE MODVAR, ONLY : basin, &
                      classdata, &
                      path, &
                      regiondivision, &
                      genpar, &
                      regpar, &
                      soilpar

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index for current subbasin
    INTEGER, INTENT(IN) :: j          !<index for class
    INTEGER, INTENT(IN) :: subid      !<subid for current subbasin
    INTEGER, INTENT(IN) :: n          !<array-dimension, number of substances
    TYPE(riverstatetype),INTENT(INOUT)  :: riverstate  !<River states
    REAL, INTENT(OUT) :: flow         !<calculated regional groundwater flow from river (m3/timestep)

    !Local variables
    INTEGER, PARAMETER :: itype=2
    INTEGER status
    
    INTEGER isoil, l
    REAL cflow(n)
    REAL rcgrwj
    REAL rivervol,waterfrac
    REAL accsubst(n)

    !>\b Algorithm \n
    !Default
    flow = 0.

    IF(path(i)%rechargebasin)THEN
      
      !>Calculate regional groundwater flow from river
      rivervol = riverstate%water(itype,i) + (SUM(riverstate%qqueue(1:ttstep(itype,i),itype,i)) + riverstate%qqueue(ttstep(itype,i)+1,itype,i) * ttpart(itype,i))
      IF(rivervol<=0.)RETURN
      rcgrwj = genpar(m_rcgrw)       !default coefficient for regional groundwater flow
      IF(j>0)THEN 
        isoil = classdata(j)%soil
        IF(soilpar(m_rcgrwst,isoil)>0 .AND. basin(i)%parregion(regiondivision(m_aqpercorr))>0)THEN
          rcgrwj = soilpar(m_rcgrwst,isoil)*(1.+regpar(m_aqpercorr,basin(i)%parregion(regiondivision(m_aqpercorr))))
        ENDIF
      ENDIF
      flow = rcgrwj * rivervol   !m3/ts
    
      IF(flow>0.)THEN
        accsubst = cgrwflow(:,i)*grwflow(i)
        !>Remove flow from river damping box...
        IF(riverstate%water(itype,i)>0)THEN
          waterfrac = riverstate%water(itype,i)/rivervol*flow
          CALL set_percolation_concentration(n,riverstate%conc(:,itype,i),cflow)  !remove_water can't handle leaving substances behind
          IF(waterfrac<riverstate%water(itype,i))THEN
            CALL remove_water(riverstate%water(itype,i),n,riverstate%conc(:,itype,i),waterfrac,cflow,status)
            IF(status.NE.0) CALL error_remove_water(errstring(3),subid,i,itype)
          ELSE
            riverstate%water(itype,i)=0.
            riverstate%conc(:,itype,i)=0.
          ENDIF
          accsubst = accsubst + cflow * waterfrac*1.E3
          grwflow(i) = grwflow(i) + waterfrac*1.E3   !mm*m2
        ENDIF
        !>...and from the queue proportional to their volume.
        DO l = 1,ttstep(itype,i)
          IF(riverstate%qqueue(l,itype,i)>0)THEN
            waterfrac = riverstate%qqueue(l,itype,i)/rivervol*flow
            CALL set_percolation_concentration(n,riverstate%cqueue(:,l,itype,i),cflow)  !remove_water can't handle leaving substances behind
            IF(waterfrac<riverstate%qqueue(l,itype,i))THEN
              CALL remove_water(riverstate%qqueue(l,itype,i),n,riverstate%cqueue(:,l,itype,i),waterfrac,cflow,status)
              IF(status.NE.0) CALL error_remove_water(errstring(4),subid,i,itype)
            ELSE
              riverstate%qqueue(l,itype,i)=0.
              riverstate%cqueue(:,l,itype,i)=0.
            ENDIF
            accsubst = accsubst + cflow * waterfrac*1.E3
            grwflow(i) = grwflow(i) + waterfrac*1.E3   !mm*m2
          ENDIF
        ENDDO
        IF(ttpart(itype,i)>0)THEN
          l = ttstep(itype,i) + 1
          IF(riverstate%qqueue(l,itype,i)>0)THEN
            waterfrac = riverstate%qqueue(l,itype,i)/rivervol*flow
            CALL set_percolation_concentration(n,riverstate%cqueue(:,l,itype,i),cflow)  !remove_water can't handle leaving substances behind
            IF(waterfrac<riverstate%qqueue(l,itype,i))THEN
              CALL remove_water(riverstate%qqueue(l,itype,i),n,riverstate%cqueue(:,l,itype,i),waterfrac,cflow,status)
              IF(status.NE.0) CALL error_remove_water(errstring(4),subid,i,itype)
            ELSE
              riverstate%qqueue(l,itype,i)=0.
              riverstate%cqueue(:,l,itype,i)=0.
            ENDIF
            accsubst = accsubst + cflow * waterfrac*1.E3
            grwflow(i) = grwflow(i) + waterfrac*1.E3   !mm*m2
          ENDIF
        ENDIF

        !>Calculate concentration of flow
        cgrwflow(:,i) = accsubst/grwflow(i)
      ENDIF
    ENDIF

  END SUBROUTINE calculate_river_groundwaterflow_removal

  !>\brief Calculate regional groundwater flow from main river flooded floodplain. 
  !>Reduce OC-, ON-, and PP-concentration of flow and remove the water from the floodplain.
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_river_floodplain_groundwaterflow_removal(i,j,subid,n,miscstate,flow)

    USE HYPEVARIABLES, ONLY : m_rcgrwst, &
                              m_rcgrw, &
                              m_aqpercorr
   USE MODVAR, ONLY : basin, &
                      classdata, &
                      path, &
                      regiondivision, &
                      genpar, &
                      regpar, &
                      soilpar

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<index for current subbasin
    INTEGER, INTENT(IN) :: j          !<index for class
    INTEGER, INTENT(IN) :: subid      !<subid for current subbasin
    INTEGER, INTENT(IN) :: n          !<array-dimension, number of substances
    TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
    REAL, INTENT(OUT) :: flow         !<calculated regional groundwater flow from river (m3/timestep)

    !Local variables
    INTEGER, PARAMETER :: itype=1 !main river
    INTEGER status
    
    INTEGER isoil
    REAL cflow(n)
    REAL rcgrwj   !model parameter
    REAL accsubst(n)

    !>\b Algorithm \n
    !Default
    flow = 0.
    IF(path(i)%rechargebasin)THEN

      !>Calculate regional groundwater flow from river floodplain
      IF(miscstate%floodwater(itype,i)<=0.)RETURN
      rcgrwj = genpar(m_rcgrw)       !default coefficient for regional groundwater flow
      IF(j>0)THEN 
        isoil = classdata(j)%soil
        IF(soilpar(m_rcgrwst,isoil)>0)THEN
          rcgrwj = soilpar(m_rcgrwst,isoil)*(1.+regpar(m_aqpercorr,basin(i)%parregion(regiondivision(m_aqpercorr))))
        ENDIF
      ENDIF
      flow = rcgrwj * miscstate%floodwater(itype,i)   !m3/ts
    
      IF(flow>0.)THEN
        accsubst = cgrwflow(:,i)*grwflow(i)
        !>Remove flow from river
        CALL set_percolation_concentration(n,miscstate%cfloodwater(:,itype,i),cflow)
        IF(flow<miscstate%floodwater(itype,i))THEN
          CALL remove_water(miscstate%floodwater(itype,i),n,miscstate%cfloodwater(:,itype,i),flow,cflow,status)
          IF(status.NE.0) CALL error_remove_water(errstring(3),subid,i,itype)
        ELSE
          miscstate%floodwater(itype,i)=0.
          miscstate%cfloodwater(:,itype,i)=0.
        ENDIF
        accsubst = accsubst + cflow * flow*1.E3
        grwflow(i) = grwflow(i) + flow*1.E3   !mm*m2
        cgrwflow(:,i) = accsubst/grwflow(i)
      ENDIF
    ENDIF

  END SUBROUTINE calculate_river_floodplain_groundwaterflow_removal

  !>During deep percolation ON, PP and OC are assumed to stay in the water source
  !>
  !> \b Reference ModelDescription Chapter Deep processes (Aquifers)
  !--------------------------------------------------------------------------
  SUBROUTINE set_percolation_concentration(n,waterconc,percconc)

    USE MODVAR, ONLY : i_on,i_pp,i_oc
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n            !<number of substances
    REAL, INTENT(IN)    :: waterconc(n) !<concentration of water source
    REAL, INTENT(OUT)   :: percconc(n)  !<concentration of outflowing water
    
    IF(n==0)RETURN
    percconc = waterconc
    IF(i_on>0) percconc(i_on) = 0.
    IF(i_pp>0) percconc(i_pp) = 0.
    IF(i_oc>0) percconc(i_oc) = 0.
    
  END SUBROUTINE set_percolation_concentration
  
  !>\brief Calculate percolation from soil layers and rivers that flows toward each aquifer.
  !!
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers - Deep percolation)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_aquifer_percolation(perc,cperc)

    !Argument declarations
    REAL, INTENT(OUT)   :: perc(naquifers)            !<water leaving soillayer (m3/timestep) 
    REAL, INTENT(OUT)   :: cperc(numsubstances,naquifers) !<concentration of water leaving soillayer 

    !Local variables
    INTEGER i,ia   !loop-variables i-subbasin, ia-aquifer

    !>\b Algorithm \n
    !>Initiations
    perc = 0.
    cperc = 0.

    !>For every subbasin:
    DO i = 1,nsub
      IF(path(i)%rechargebasin)THEN
        !>\li Add subbasin groundwater seepage to aquifer percolation
        perc(path(i)%aquid) = perc(path(i)%aquid) + grwflow(i)*1.E-3
        cperc(:,path(i)%aquid) = cperc(:,path(i)%aquid) + cgrwflow(:,i)*grwflow(i)*1.E-3  !amount
        IF(grwflow(i)>0.) cgrwflow(:,i) = cgrwflow(:,i) / grwflow(i)       !conc, mg/L (or other)
      ENDIF
    ENDDO

    !>Calculate concentration of percolation
    DO ia=1,naquifers
      IF(perc(ia)>0) cperc(:,ia) = cperc(:,ia)/perc(ia)   !concentration
    ENDDO

  END SUBROUTINE calculate_aquifer_percolation

  !>\brief Calculate aquifer flows.
  !>
  !>Calculate aquifer flows, including removal of regional groundwater flow 
  !!from all subbasins and classes (not lakes), delay of recharge, mixing of 
  !!aquifer and outflow from aquifer.
  !!
  !>\b Consequences Module variables grwflow and cgrwflow may change.
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_aquifer(aquiferstate,aqoutflow,remflow2,aquiferirrloss)

    USE HYPEVARIABLES, ONLY : m_denitaq,m_hsatINsoil
    USE MODVAR, ONLY : genpar
    
    !Arument declarations
    TYPE(aquiferstatetype),INTENT(INOUT):: aquiferstate !<Aquifer states
    REAL, INTENT(OUT)   :: aqoutflow(naquifers)        !<Outflow from aquifer (m3/timestep)
    REAL, INTENT(OUT)   :: remflow2(nsub)              !<Flow from soillayers toward aquifer (m3/timestep)
    REAL, INTENT(IN)    :: aquiferirrloss(naquifers)   !<Flow from aquifer to irrigation (m3/timestep)

    !Local variables
    INTEGER status
    INTEGER ia   !loop-variable ia-aquifer
    REAL deepperc(naquifers)    !water leaving soillayer-model (m3/timestep) 
    REAL cdeepperc(numsubstances,naquifers) !concentration of water -"- 
    REAL recharge(naquifers)    !water recharging aquifer (m3/timestep) 
    REAL crecharge(numsubstances,naquifers) !concentration of water -"- 
    REAL chelp(numsubstances)   !concentration
    REAL caqoutflow(numsubstances) !concentration of outflow from aquifer
    REAL sink(numsubstances)    !not used

    !>\b Algorithm \n
    !>Initiations for time step
    aqoutflow = 0.

    !>Accumulate percolation from soil and river (all subbasins) to aquifer percolation
    CALL calculate_aquifer_percolation(deepperc,cdeepperc)
    
    !Recharge delay and substances
    CALL calculate_percolation_delay(aquiferstate,naquifers,numsubstances,deepperc,cdeepperc,recharge,crecharge)
    
    !>Add recharge to aquifer
    DO ia=1,naquifers
      chelp = crecharge(:,ia)
      CALL add_water(numsubstances,aquiferstate%water(ia),aquiferstate%conc(:,ia),recharge(ia),chelp)
    ENDDO
    
    !Calculate aquifer nutrient processes
    DO ia=1,naquifers
      CALL soil_denitrification(aquifer(ia)%maxvol,genpar(m_denitaq),genpar(m_hsatINsoil),aquiferstate%water(ia),aquifer(ia)%temperature,aquiferstate%conc(:,ia),sink)
    ENDDO
    
    !>Calculate aquifer outflow
    DO ia=1,naquifers
      CALL calculate_aquifer_outflow(ia,aquiferstate,aqoutflow(ia),caqoutflow)
      aquiferstate%nextoutflow(ia) = aqoutflow(ia)
      IF(numsubstances>0) aquiferstate%cnextoutflow(:,ia) = caqoutflow
      caqoutflow = aquiferstate%conc(:,ia)
      CALL remove_water(aquiferstate%water(ia),numsubstances,aquiferstate%conc(:,ia),aquiferirrloss(ia),caqoutflow,status)
      IF(status.NE.0) CALL error_remove_water(errstring(5),ia,0,0)
    ENDDO

    !>Set output variables
    remflow2 = grwflow*1.E-3
    
  END SUBROUTINE calculate_aquifer

  !>\brief Calculate deep percolation delay
  !!
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers - Deep percolation delay and aquifer recharge)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_percolation_delay(aquiferstate,na,ns,perc,cperc,recharge,crecharge)

    !Argument declarations
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    INTEGER, INTENT(IN) :: na               !<Number of aquifers
    INTEGER, INTENT(IN) :: ns               !<Number of substances
    REAL, INTENT(IN)    :: perc(na)         !<water leaving soillayer-model (m3/timestep) 
    REAL, INTENT(IN)    :: cperc(ns,na)     !<concentration of water -"- 
    REAL, INTENT(OUT)   :: recharge(na)     !<water entering aquifer (m3/timestep) 
    REAL, INTENT(OUT)   :: crecharge(ns,na) !<concentration of water -"- 

    !Local variables
    INTEGER ia    !loop-variables ia-aquifer

    !>\b Algorithm \n
    crecharge = 0.
    
    !>Calculate delay
    DO ia = 1,na
      recharge(ia) = (1.-dcoeff(ia))*perc(ia)+dcoeff(ia)*aquiferstate%lastrecharge(ia)
      IF(ns>0)THEN
        IF(recharge(ia)>0) crecharge(:,ia) = ((1.-dcoeff(ia))*perc(ia)*cperc(:,ia)+dcoeff(ia)* &
          aquiferstate%lastrecharge(ia)*aquiferstate%clastrecharge(:,ia))/recharge(ia)
        aquiferstate%clastrecharge(:,ia) = crecharge(:,ia)
      ENDIF
      aquiferstate%lastrecharge(ia) = recharge(ia)
    ENDDO

  END SUBROUTINE calculate_percolation_delay

  !>\brief Calculate the percolating water in delay of reaching aquifer
  !!
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers - Deep percolation delay and aquifer recharge)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_delayed_water(aquiferstate,na,water)

    !Argument declarations
    TYPE(aquiferstatetype),INTENT(IN) :: aquiferstate   !<Aquifer states
    INTEGER, INTENT(IN) :: na               !<Number of aquifers
    REAL, INTENT(OUT)   :: water(na)        !<water in delay (m3/timestep) 

    !Local variables
    INTEGER ia    !loop-variable, ia-aquifer

    !>\b Algorithm \n
    !>Calculate currently delayed water
    DO ia = 1,na
      water(ia) = aquiferstate%lastrecharge(ia) * dcoeff(ia)/(1.-dcoeff(ia))
    ENDDO

  END SUBROUTINE calculate_delayed_water

  !>\brief Subroutine for calculation of outflow from aquifer.
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_aquifer_outflow(ia,aquiferstate,aqoutflow,caqoutflow)

    USE MODVAR, ONLY : regpar, &
                       regiondivision
    USE HYPEVARIABLES, ONLY : m_aqretcorr

    !Argument declarations
    INTEGER, INTENT(IN) :: ia             !<index of current aquifer
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer state
    REAL, INTENT(OUT)   :: aqoutflow            !<aquifer outflow (m3/timestep)
    REAL, INTENT(OUT)   :: caqoutflow(numsubstances)        !<Concentration of outflow from aquifer

    !Local variables
    INTEGER status      !status of routine
    REAL outpar         !recession coefficient for current aquifer outflow

    !Initial values

    !Current parameter values
    IF(aquifer(ia)%parregion(regiondivision(m_aqretcorr))>0) outpar = aquifer(ia)%retrate*(1.+regpar(m_aqretcorr,aquifer(ia)%parregion(regiondivision(m_aqretcorr))))
    IF(outpar>1) outpar = 1.

    !Outflow from aquifer (recession coefficient method)
    aqoutflow = outpar * aquiferstate%water(ia)
    caqoutflow = aquiferstate%conc(:,ia)
    CALL remove_water(aquiferstate%water(ia),numsubstances,aquiferstate%conc(:,ia),aqoutflow,caqoutflow,status)
    IF(status.NE.0) CALL error_remove_water(errstring(2),ia,0,0)
    
  END SUBROUTINE calculate_aquifer_outflow

  !>\brief Subroutine for division of outflow from aquifer to main river.
  !
  !>\b Consequences Module variables addflow and caddflow will be set
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_current_flow_from_aquifer(na,ns,aquiferstate,outsideflow)

    !Argument declarations
    INTEGER, INTENT(IN) :: na             !<number of aquifers
    INTEGER, INTENT(IN) :: ns             !<number of substances
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer state
    REAL, INTENT(OUT)   :: outsideflow(na)         !<Outflow from aquifer to outside model system (m3/timestep)

    !Local variables
    INTEGER i,ia                !Subbasin and aquifer loop index
    REAL aqoutflow              !aquifer outflow (m3/timestep)
    REAL caqoutflow(ns)         !Concentration of outflow from aquifer
    REAL totoutflow             !Aquifer summed outflow (for check)

    !Initial values
    addflow = 0.    !module variable
    caddflow = 0.   !module variable
    outsideflow = 0.

    !Outflow from aquifer (recession coefficient method)
    DO ia = 1,naquifers
      
      aqoutflow = aquiferstate%nextoutflow(ia)
      caqoutflow = aquiferstate%cnextoutflow(:,ia)
      
      IF(aqoutflow>0)THEN
        !Relay outflow to recieving rivers
        totoutflow = 0.
        DO i=1,nsub
          IF(path(i)%recievefraction>0.)THEN
            IF(path(i)%aquid==ia)THEN
              caddflow(:,i) = caddflow(:,i)*addflow(i) + caqoutflow(:)*path(i)%recievefraction*aqoutflow
              addflow(i) = addflow(i) + path(i)%recievefraction*aqoutflow
              caddflow(:,i) = caddflow(:,i)/addflow(i)
              totoutflow = totoutflow + path(i)%recievefraction*aqoutflow
            ENDIF
          ENDIF
        ENDDO
        
        !Check flow outside model system
        IF(totoutflow<aqoutflow)THEN
          outsideflow(ia) = aqoutflow - totoutflow
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE calculate_current_flow_from_aquifer

  !>\brief Add regional lateral groundwater flow from other subbasins
  !>to soil
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Regional groundwater flow)
  !------------------------------------------------------------------------
  SUBROUTINE add_regional_groundwater_flow_to_soil(i,j,classarea,pw,soilstate,load,upwardflow,addedflow)

    USE MODVAR, ONLY : path,    &
                       nsub,            &
                       numsubstances,   &
                       maxsoillayers,   &
                       soilthick

    !Argument declarations
    INTEGER, INTENT(IN) :: i                   !<current subbasin index
    INTEGER, INTENT(IN) :: j                   !<current class index
    REAL, INTENT(IN)    :: classarea           !<class area (km2)
    REAL, INTENT(IN)    :: pw(maxsoillayers)   !<pore volume (mm)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    REAL, INTENT(INOUT) :: load(numsubstances) !<load into soil (kg/timestep)
    REAL, INTENT(OUT) :: upwardflow(2)         !<upwelling due to overflowing lower soil layers (mm/timestep)
    REAL, INTENT(INOUT) :: addedflow(maxsoillayers,nsub) !<added regional groundwater flow (m3/timestep)
    
    !Local variables
    INTEGER k       !loop-variable
    INTEGER lowestsoillayer !soil layer rural load was added to
    REAL grwflowj   !flow to this class 
    REAL localupward(2),localadded(maxsoillayers)

    !>\b Algorithm \n
    upwardflow = 0.

    !>For every subbasin that has regional flow to current subbasin:
    DO k = 1,nsub
      IF(path(k)%grw1==i)THEN
        !>\li Add regional groundwater flow to soil
        grwflowj = grwflow(k) * (1. - path(k)%grwtolake) / regrwarea(i)   !->mm
        CALL inflow_lowest_soillayer(numsubstances,maxsoillayers, &
                cgrwflow(:,k),grwflowj,pw,soilthick(:,j),soilstate%water(:,j,i),  &
                soilstate%conc(:,:,j,i),localupward,localadded,lowestsoillayer)
        !>\li Accumulate load of regional groundwater flow into soil (kg/timestep)
        load = load + cgrwflow(:,k)*grwflowj*classarea  
        upwardflow = upwardflow + localupward
        addedflow(:,k) = addedflow(:,k) + localadded * classarea * 1.E3   !m3
      ENDIF
    ENDDO

  END SUBROUTINE add_regional_groundwater_flow_to_soil

  !>Add regional groundwater flow to outlet lake of subbasin
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Regional groundwater flow - Regional groundwater to outlet lake)
  !--------------------------------------------------------------------------
  SUBROUTINE add_regional_groundwater_flow_to_olake(i,itype,numsubst,   &
       qfactor,lakev,lakec,qcinfl,grwload,addedwater)

    USE MODVAR, ONLY : path,      &
         seconds_per_timestep

    !Argument declarations
    INTEGER, INTENT(IN) :: i              !<index for subbasin
    INTEGER, INTENT(IN) :: itype          !<lake type (ilake or olake)
    INTEGER, INTENT(IN) :: numsubst       !<number of substances (numsubstances)
    REAL, INTENT(IN)    :: qfactor        !<flow unit transformation factor (m3/s->mm/timestep)
    REAL, INTENT(INOUT) :: lakev          !<lake volume (i.e. water stage) (mm)
    REAL, INTENT(INOUT) :: lakec(numsubst)     !<concentration of lake
    REAL, INTENT(INOUT) :: qcinfl         !<total inflow to olake (m3/s) (print out variable)
    REAL, INTENT(OUT)   :: grwload(numsubst) !<load added (kg/timestep)
    REAL, INTENT(OUT)   :: addedwater        !<water added (m3/timestep)
    
    !Local variables
    REAL grwflowj   !reg. grw flow to current olake

    !>\b Algorithm \n
    !>Return if not outlet lake
    grwload = 0.;addedwater = 0.
    IF(itype/=2) RETURN
    
    !>If regional groundwater flow to current lake:
    IF(path(i)%grwtolake>0 .AND. grwflow(i)>0)THEN

      !>\li Calculate regional groundwater flow to lake and add to variable for inflow to lake
      addedwater = grwflow(i) * path(i)%grwtolake * 1.E-3   !mm*m2 -> m3/timestep
      grwflowj = addedwater / seconds_per_timestep   !m3/s
      qcinfl = qcinfl + grwflowj
      grwflowj = grwflowj * qfactor                                     !m3/s -> mm/timestep

      !>Add regional groundwater flow to lake
      CALL add_water(numsubst,lakev,lakec,grwflowj,cgrwflow(:,i))

      !>Calculate (nutrient) load added to lake by regional groundwater flow
      grwload = cgrwflow(:,i)*grwflow(i)*path(i)%grwtolake*1.E-6      !Load of regional groundwater flow into this subbasin olake (kg/timestep)
    ENDIF 

  END SUBROUTINE add_regional_groundwater_flow_to_olake

  !>\brief Add aquifer outflow to main river
  !>
  !>\b Reference ModelDescription Chapter Deep processes (Aquifers)
  !--------------------------------------------------------------------------
  SUBROUTINE add_aquifer_flow_to_river(i,numsubst,qin,cin,addedwater,addedload)

    !Argument declarations
    INTEGER, INTENT(IN) :: i              !<index for subbasin
    INTEGER, INTENT(IN) :: numsubst       !<number of substances (numsubstances)
    REAL, INTENT(INOUT) :: qin            !<inflow to river (m3/s)
    REAL, INTENT(INOUT) :: cin(numsubst)  !<concentration of river inflow
    REAL, INTENT(OUT)   :: addedwater     !<water added (m3/timestep)
    REAL, INTENT(OUT)   :: addedload(numsubst)  !<load to river inflow
    
    !Local variables
    REAL qadd   !aquifer flow to current subbasin (m3/s)
    REAL cadd(numsubst)   !conc of aquifer flow

    !>\b Algorithm \n
    addedwater = addflow(i)
    addedload = 0.
    
    !>Return if no outflow here
    IF(addedwater==0) RETURN
    
    !>Add aquifer flow to river inflow
    qadd = addedwater / seconds_per_timestep
    cadd = caddflow(:,i)
    IF(qin>0)THEN
      cin = (qin * cin + qadd * cadd)/(qin + qadd)
      qin = qin + qadd
    ELSE
      qin = qadd
      cin = cadd
    ENDIF
    addedload = addedwater*cadd*1.E-3   !kg/ts

  END SUBROUTINE add_aquifer_flow_to_river

  !>Calculate aquifer water level below ground
  !----------------------------------------------
  FUNCTION calculate_aquifer_waterlevel(ns,na,aquiferstate) RESULT(subbwl)

    !Argument declarations  
    INTEGER, INTENT(IN) :: ns             !<Number of subbasins
    INTEGER, INTENT(IN) :: na             !<Number of aquifers
    TYPE(aquiferstatetype),INTENT(IN) :: aquiferstate  !<Aquifer state
    REAL subbwl(ns)                       !<Aquifer water depth (m)
    
    INTEGER i,ia    !loop-index, i-subbasin,ia-aquifer 
    REAL aqwl(na) !aquifer water depth (m)
    
    DO ia=1,na
      aqwl(ia) = aquifer(ia)%basedepth + aquiferstate%water(ia)/aquifer(ia)%porosity/aquifer(ia)%area
    ENDDO
    subbwl = missing_value
    DO i=1,ns
      IF(path(i)%rechargebasin .OR. path(i)%recievefraction>0) subbwl(i) = aqwl(path(i)%aquid)
    ENDDO

  END FUNCTION calculate_aquifer_waterlevel

END MODULE REGIONAL_GROUNDWATER_MODULE
