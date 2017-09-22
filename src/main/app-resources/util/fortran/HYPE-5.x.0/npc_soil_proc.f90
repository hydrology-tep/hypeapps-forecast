!> \file npc_soil_proc.f90
!> Contains module npc_soil_processes.

!>Nitrogen, phosphorus and organic carbon processes in soil in HYPE
MODULE NPC_SOIL_PROCESSES

!Copyright 2012-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.


  !Used modules
  USE GENERAL_WATER_CONCENTRATION
  USE GENERAL_FUNCTIONS
  USE STATETYPE_MODULE, ONLY : snowicestatetype,soilstatetype,miscstatetype
  !Uses also modvar, hypevariables
  IMPLICIT NONE
  PRIVATE
  !-------------------------------------------------
  ! Private procedures
  !-------------------------------------------------
  ! crop_sources
  ! calculate_erosion 
  ! calculate_hbvsed_erosion
  ! calculate_pp_transport
  ! soil_pool_transformations
  ! plant_uptake 
  ! freundlich 
  ! soil_carbon_pool_transformations
  ! riparian_moisturefactor 
  !-------------------------------------------------
  PUBLIC :: initiate_soil_npc_state, &
            set_class_precipitation_concentration_and_load, &
            add_dry_deposition_to_landclass, &
            atmdep_in_loss, &
            calculate_plant, &
            particle_processes_for_runoff, &
            soil_np_processes, &
            soil_denitrification,  &
            balance_spsoil, &
            croprotation_soilpoolaverage, &
            local_diffuse_source, &
            soil_carbon_processes, &
            doc_percolation_reduction, &
            onpp_percolation_reduction, &
            class_riparian_zone_processes
  
CONTAINS
  
  !>\brief Initiation of soil for nutrient concentrations and soil pools.
  !!Nutrient pools in kg/km2. 
  !---------------------------------------------------------------------------
  SUBROUTINE initiate_soil_npc_state(initN,initP,initC,maxlayers,soilstate)

    USE HYPEVARIABLES, ONLY : m_humusP0,m_humusN0,  &
                              m_fastN0,m_fastP0,    &
                              m_partP0,m_partP1,m_partP2,m_partP3,  &
                              m_hNhalf,m_hPhalf,m_pPhalf,  &
                              m_humusC1,m_fastC1, &
                              m_humusC2,m_fastC2, &
                              m_humusC3,m_fastC3, &
                              m_phoscorr,         &
                              m_onconc0,m_ppconc0,  &
                              m_occonc0, &
                              m_soilstretch
    USE MODVAR, ONLY : classdata,     &
                       basin,         &
                       nsub, nclass,  &
                       regiondivision, &
                       maxsoillayers, &
                       soilthick,     &
                       soildepthstretch, &
                       i_on,i_pp,i_oc, & 
                       genpar,landpar,regpar

    !Argument declarations
    LOGICAL, INTENT(IN) :: initN     !<flag for initiation of nitrogen model
    LOGICAL, INTENT(IN) :: initP     !<flag for initiation of phosphorus model
    LOGICAL, INTENT(IN) :: initC     !<flag for initiation of organic carbon model
    INTEGER, INTENT(IN) :: maxlayers !<maximum number of soil layers in the model set-up
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    
    !Local variables
    INTEGER i,j           !loop-variables (subbasin,class)
    REAL    halfHN,halfHP,halfPP   !help variables
    REAL    phoscorr        !correction of phosphorus level
    REAL :: tempsoilthick(maxsoillayers,nclass)
    REAL :: help(maxsoillayers)

    !>\b Algorithm \n
    !>Calculate temporary soil layer thickness (m)
    tempsoilthick = soilthick
    IF(soildepthstretch)THEN
      DO j = 1,nclass
        IF(landpar(m_soilstretch,classdata(j)%luse)>0.)THEN
          tempsoilthick(2,j) = tempsoilthick(2,j) * landpar(m_soilstretch,classdata(j)%luse)
          tempsoilthick(3,j) = tempsoilthick(3,j) * landpar(m_soilstretch,classdata(j)%luse)
        ENDIF
      ENDDO
    ENDIF
    !>Initialize soil nitrogen pools (kg/km2)
    IF(initN)THEN
      !Initiate first subbasin
      DO j = 1,nclass   !Initiate with top layer concentration
        soilstate%humusN(:,j,1) = landpar(m_humusN0,classdata(j)%luse) * tempsoilthick(:,j)
        soilstate%fastN(:,j,1)  = genpar(m_fastN0) * tempsoilthick(:,j)
        IF(maxlayers>1)THEN   !Adjust according to soillayer decline in concentration        
          halfHN = log(2.) / landpar(m_hNhalf,classdata(j)%luse)
          soilstate%humusN(2,j,1) = soilstate%humusN(2,j,1) * exp(-halfHN * ((tempsoilthick(1,j)/2.) + (tempsoilthick(2,j)/2.)))
          soilstate%humusN(3,j,1) = soilstate%humusN(3,j,1) * exp(-halfHN * ((tempsoilthick(1,j)/2.) + tempsoilthick(2,j) +  (tempsoilthick(3,j)/2.)))
        ENDIF
        !Initiate the rest of the subbasins
        DO i = 2,nsub
          soilstate%humusN(:,j,i) = soilstate%humusN(:,j,1)
          soilstate%fastN(:,j,i)  = soilstate%fastN(:,j,1)
        ENDDO
      ENDDO
    ENDIF

    !>Initialize soil phosphorus pools (kg/km2)
    IF(initP)THEN
      DO i = 1,nsub
        IF(basin(i)%parregion(regiondivision(m_phoscorr))>0)THEN
          phoscorr = 1. + regpar(m_phoscorr,basin(i)%parregion(regiondivision(m_phoscorr)))   !Correction of phosphorus
        ELSE
          phoscorr  = 1.
        ENDIF
        DO j = 1,nclass   
          soilstate%humusP(:,j,i) = phoscorr*landpar(m_humusP0,classdata(j)%luse) * tempsoilthick(:,j) !Initiate with top layer concentration
          halfHP = log(2.) / landpar(m_hPhalf,classdata(j)%luse)    
          soilstate%humusP(2,j,i) = soilstate%humusP(2,j,i) * exp(-halfHP * ((tempsoilthick(1,j)/2.) + (tempsoilthick(2,j)/2.)))           
          soilstate%humusP(3,j,i) = soilstate%humusP(3,j,i) * exp(-halfHP * ((tempsoilthick(1,j)/2.) + tempsoilthick(2,j) +  (tempsoilthick(3,j)/2.)))  
          soilstate%fastP(:,j,i)  = phoscorr*genpar(m_fastP0) * tempsoilthick(:,j)
          IF(landpar(m_partP0,classdata(j)%luse)>0)THEN
            soilstate%partP(:,j,i)  = phoscorr*landpar(m_partP0,classdata(j)%luse) * tempsoilthick(:,j)
            halfPP = log(2.) / landpar(m_pPhalf,classdata(j)%luse)
            soilstate%partP(2,j,i) = soilstate%partP(2,j,i) * exp(-halfPP * ((tempsoilthick(1,j)/2.) +  (tempsoilthick(2,j)/2.)))            
            soilstate%partP(3,j,i) = soilstate%partP(3,j,i) * exp(-halfPP * ((tempsoilthick(1,j)/2.) + tempsoilthick(2,j) + (tempsoilthick(3,j)/2.)))            
          ELSE
            soilstate%partP(1,j,i)  = phoscorr*landpar(m_partP1,classdata(j)%luse) * tempsoilthick(1,j)
            soilstate%partP(2,j,i)  = phoscorr*landpar(m_partP2,classdata(j)%luse) * tempsoilthick(2,j)
            soilstate%partP(3,j,i)  = phoscorr*landpar(m_partP3,classdata(j)%luse) * tempsoilthick(3,j)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    !>Initialize soil organic carbon pools (kg/km2)
    IF(initC)THEN
      DO j = 1,nclass
        help(1) = landpar(m_humusC1,classdata(j)%luse) * tempsoilthick(1,j)
        help(2) = landpar(m_humusC2,classdata(j)%luse) * tempsoilthick(2,j)
        help(3) = landpar(m_humusC3,classdata(j)%luse) * tempsoilthick(3,j)
        soilstate%humusC(:,j,1) = help(:)
      ENDDO  
      DO i = 2,nsub
        soilstate%humusC(:,:,i) = soilstate%humusC(:,:,1)
      ENDDO
      DO j = 1,nclass
        help(1) = landpar(m_fastC1,classdata(j)%luse) * tempsoilthick(1,j)
        help(2) = landpar(m_fastC2,classdata(j)%luse) * tempsoilthick(2,j)
        help(3) = landpar(m_fastC3,classdata(j)%luse) * tempsoilthick(3,j)
        soilstate%fastC(:,j,1) = help(:)
      ENDDO  
      DO i = 2,nsub
        soilstate%fastC(:,:,i) = soilstate%fastC(:,:,1)
      ENDDO
    ENDIF

    !>Initialize soil organic carbon soil water concentration (mg/L)
    IF(initC)THEN
      DO i = 1,nsub
        DO j = 1,nclass
          soilstate%conc(i_oc,:,j,i) = landpar(m_occonc0,classdata(j)%luse)
        ENDDO
      ENDDO
    ENDIF

    !>Initilize nutrient concentration of soil (mg/L)
    IF(initN)THEN
      DO i = 1,nsub
        DO j = 1,nclass
          soilstate%conc(i_on,:,j,i) = landpar(m_onconc0,classdata(j)%luse)
        ENDDO
      ENDDO
    ENDIF
    IF(initP)THEN
      DO i = 1,nsub
        DO j = 1,nclass
          soilstate%conc(i_pp,:,j,i) = landpar(m_ppconc0,classdata(j)%luse)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE initiate_soil_npc_state

  !>Calculate class concentration of precipitation and nutrient load
  !!
  !>\b Reference ModelDescription Chapter Processes above ground (Atmospheric deposition of nitrogen and phosphorus)
  !---------------------------------------------------------------------
  SUBROUTINE set_class_precipitation_concentration_and_load(ns,area,precorg,temp,prec,cprecorg,cprec,precload)
    USE MODVAR, ONLY : i_in,i_sp,i_t2,    &
                       genpar
    USE HYPEVARIABLES, ONLY : m_atmload

    !Argument declarations
    INTEGER,INTENT(IN)    :: ns           !<number of substances, array dimension
    REAL,   INTENT(IN)    :: area         !<class area (km2)
    REAL,   INTENT(IN)    :: precorg      !<precipitation from Pobs.txt (mm/timestep)
    REAL,   INTENT(IN)    :: temp         !<temperature of class (degree Celsius)
    REAL,   INTENT(IN)    :: prec         !<precipitation of class (mm/timestep)
    REAL,   INTENT(IN)    :: cprecorg(ns) !<concentration of precipitation (Pobs)
    REAL,   INTENT(OUT)   :: cprec(ns)    !<concentration of precipitation (prec)
    REAL,   INTENT(INOUT) :: precload(ns) !<nutrient load of precipitation (kg/timestep)

    cprec = 0
    !Wet deposition NP
    IF(prec>0)THEN
      IF(genpar(m_atmload)==1)THEN
        cprec = cprecorg*(precorg/prec)
      ELSE
        cprec = cprecorg
      ENDIF    
      IF(i_in>0)THEN
        precload(i_in) = cprec(i_in)*prec*area
      ENDIF
      IF(i_sp>0)THEN
        precload(i_sp) = cprec(i_sp)*prec*area
      ENDIF

      !Temp.conc. in prec = temp(air)|max=0 for landclasses [DG/JS Temp.model May-2013]
      IF(i_t2>0) cprec(i_t2)= MAX(0.,temp) 
    ENDIF

  END SUBROUTINE set_class_precipitation_concentration_and_load

  !>\brief Calculate dry atmospheric deposition of N and P and add it
  !>to snow or soil on land classes
  !!
  !> \b Reference ModelDescription Processes above ground (Atmospheric deposition of nitrogen and phosphorus)
  !---------------------------------------------------------------------
  SUBROUTINE add_dry_deposition_to_landclass(i,j,iluse,conductN,conductP, &
                                  areaij,source,veg,dryin,drypp,frozenstate,soilstate)

    USE MODVAR, ONLY : numsubstances,     &
         i_in,i_pp

    !Argument declarations
    INTEGER, INTENT(IN) :: i                          !<current subbasin
    INTEGER, INTENT(IN) :: j                          !<current class
    INTEGER, INTENT(IN) :: iluse                      !<index landuse
    LOGICAL, INTENT(IN) :: conductN                   !<status of N simulation
    LOGICAL, INTENT(IN) :: conductP                   !<status of P simulation
    REAL, INTENT(IN)    :: areaij                     !<classarea (km2)
    REAL, INTENT(INOUT) :: source(numsubstances)      !<dry deposition (kg/timestep)
    INTEGER, INTENT(IN) :: veg                        !<index vegetation
    REAL, INTENT(IN)    :: dryin(3)                   !<dry deposition IN (kg/km2/timestep)
    REAL, INTENT(IN)    :: drypp                      !<dry deposition PP (kg/km2/timestep)
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    
    !Local variables
    REAL fakeflow
    REAL sourceN
    REAL, DIMENSION(1) :: sourceP
    REAL :: sourceDD(numsubstances),sourceD2(numsubstances)

    IF(.NOT.conductN.AND..NOT.conductP) RETURN

    !Prepare for dry deposition
    sourceDD = 0.
    IF(i_in>0 .AND. veg>0) sourceDD(i_in) = dryin(veg)
    IF(i_pp>0) sourceDD(i_pp) = drypp

    !Add dry deposition of Inorg-N and PartP to snow or soil
    IF(frozenstate%snow(j,i)>0)THEN
      CALL add_source_to_water(frozenstate%snow(j,i),numsubstances,frozenstate%csnow(:,j,i),sourceDD)
    ELSE 
      IF(conductN .AND. veg>0)THEN
        sourceN = dryin(veg)
        fakeflow = 1.
        CALL atmdep_in_loss(iluse,fakeflow,soilstate%fastN(:,j,i),sourceN)   !add some of dep to fastN, changes sourceN
        sourceD2 = 0.
        sourceD2(i_in) = sourceN
        CALL add_source_to_water(soilstate%water(1,j,i),numsubstances,soilstate%conc(:,1,j,i),sourceD2)
      ENDIF
      IF(conductP)THEN
        sourceP = drypp
        CALL production_pool(1,soilstate%partP(1,j,i),sourceP)
      ENDIF
    ENDIF

    !Calculate atmospheric dry deposition loads (kg/timestep)
    source(:) = sourceDD*areaij

  END SUBROUTINE add_dry_deposition_to_landclass

  !>Transformation of inorganic N from atmospheric deposition (wet and
  !>dry through infiltration) to solid soil organic N.
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines 
  !>(Vegetation and soil surface processes - Transformation of nitrogen from atmospheric deposition)
  !----------------------------------------------------------
  SUBROUTINE atmdep_in_loss(iluse,flow,fastN,conc)

    USE MODVAR, ONLY : maxsoillayers,   & 
         landpar
    USE HYPEVARIABLES, ONLY : m_ponatm

    !Argument declarations
    INTEGER, INTENT(IN) :: iluse   !<index landuse
    REAL, INTENT(IN)    :: flow    !<infiltration (mm/timestep)
    REAL, INTENT(INOUT) :: fastN(maxsoillayers)    !<pool fastN (kg/km2)
    REAL, INTENT(INOUT) :: conc    !<concentration of IN in infiltration (mg/L)

    !Variable declarations
    REAL, DIMENSION(1) :: sourceN  !mg/L*mm=kg/km2
    REAL onpart !fraction of IN (from deposition) that is transformed to ON when reaching soil (-)

    !Set and check parameter for transformation
    onpart = landpar(m_ponatm,iluse)
    IF(onpart==0) RETURN

    !Part of infiltration IN goes to organic fastN-pool
    sourceN = conc*onpart*flow
    CALL production_pool(1,fastN(1),sourceN)
    conc = conc*(1.-onpart)

  END SUBROUTINE atmdep_in_loss

  !>Calculate the plant growth per crop and take the average of the
  !nutrient uptake for the class.
  !!
  !!TODO: Calculate plant do not work for shorter time steps
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines 
  !>(Vegetation and soil surface processes - Potential vegetation uptake of nitrogen)
  !------------------------------------------------------------
  SUBROUTINE calculate_plant(i,j,dayno,temp,daylength,common_uptake,miscstate)

    USE MODVAR, ONLY : maxsoillayers,   &
                       cropdata,        &
                       find_croppart,   &
                       modeloption,p_growthstart, &
                       i_in,i_sp

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<index of subbasin
    INTEGER, INTENT(IN) :: j                        !<index of class
    INTEGER, INTENT(IN) :: dayno                    !<pseudo day number of the year
    REAL, INTENT(IN)    :: temp                     !<class temperature
    REAL, INTENT(IN)    :: daylength                !<day length (hours)
    REAL, INTENT(OUT)   :: common_uptake(2,2)       !<plant uptake (kg/km2/day), (N:P,sl1:sl2)
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<Miscstate states
    
    !Local variables
    INTEGER k,kcrop
    LOGICAL calcN,calcP
    REAL uptakeN
    REAL part                           !part of class area that is the current crop
    REAL up1,up2,up3,bd2,bd3,bd5        !crop and cultivation parameters
    REAL sl1part,PNratio                !uptake parameters
    REAL help, tmpfcn

    calcN = .FALSE.
    calcP = .FALSE.
    IF(i_in>0) calcN = .TRUE.
    IF(i_sp>0) calcP = .TRUE.

    common_uptake = 0.
    IF(.NOT.calcN .AND. .NOT.calcP) RETURN

    DO kcrop = 1, 2                             !main and secondary crop
      CALL find_croppart(i,j,kcrop,k,part)
      IF(k==0) CYCLE                            !no crop given, no uptake
      IF(part>0)THEN

        !Find the growing start dayno
        bd2 = cropdata(k)%baredayno2          !sow date, default
        IF(modeloption(p_growthstart)==1)THEN
          IF(dayno == 1) THEN            !reset (ok for the southern hemisphere also)
            miscstate%gdd(kcrop,j,i) = 0.
            miscstate%gsbegin(kcrop,j,i) = 0
          ENDIF
          bd2 = miscstate%gsbegin(kcrop,j,i)
          IF(daylength >= cropdata(k)%daylength .AND. dayno >= cropdata(k)%firstday)THEN        !criteria met, accumulate GDD
            miscstate%gdd(kcrop,j,i) = miscstate%gdd(kcrop,j,i) + MAX(0. , temp - cropdata(k)%basetemp)   !Don't work for shorter time steps!?
          ENDIF    
          IF(miscstate%gdd(kcrop,j,i) > cropdata(k)%gddsow)THEN
            IF(miscstate%gsbegin(kcrop,j,i) == 0)THEN 
              bd2 = dayno
              miscstate%gsbegin(kcrop,j,i) = dayno
            ENDIF
          ELSE
            bd2 = 366       !no uptake yet
          ENDIF
        ENDIF
        
        !Data of this crop
        up1 = cropdata(k)%uptake1
        up2 = cropdata(k)%uptake2
        up3 = cropdata(k)%uptake3
        bd3 = cropdata(k)%baredayno3          !harvest, end of growth
        bd5 = cropdata(k)%baredayno5          !autumn sow data
        sl1part = cropdata(k)%uptakeupper  !fraction of uptake in soil layer 1
        PNratio = cropdata(k)%PNuptakeRatio   !PN ration of uptake
        help = (up1-up2)*exp(-up3*(dayno-bd2))
        uptakeN = up1*up2*up3*help/(up2+help)/(up2+help)!potential uptake of N, from sawing to autumn ploughing
        IF(dayno<bd2) uptakeN = 0.
        IF(dayno>bd3) THEN
          uptakeN = 0.
          IF(bd5>0 .AND. dayno>bd5) THEN        !autumn crops grow after autumn sawing
            IF(temp<5.) THEN       
              tmpfcn=0.
            ELSE
              tmpfcn = min(1., (1./20.)*temp-(5./20.)) 
            ENDIF
            help = (up1-up2)*exp(-up3*(dayno-(bd5+25)))
            uptakeN = tmpfcn * up1*up2*up3*help/(up2+help)/(up2+help) 
          ENDIF
        ENDIF

        !Sum to common uptake function
        common_uptake(1,1) = common_uptake(1,1) + part * uptakeN * sl1part                  !g N/m2/day
        common_uptake(1,2) = common_uptake(1,2) + part * uptakeN * (1.-sl1part)             !g N/m2/day
        common_uptake(2,1) = common_uptake(2,1) + part * uptakeN * sl1part * PNratio        !g P/m2/day
        common_uptake(2,2) = common_uptake(2,2) + part * uptakeN * (1.-sl1part) * PNratio   !g P/m2/day

      ENDIF
    ENDDO
    common_uptake = common_uptake * 1000.       !g/m2->kg/km2

  END SUBROUTINE calculate_plant

  !>\brief Calculate and add the crops sources of nutrient to the soil
  !(fertilization and crop residues).
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Nutrient sources - 
  !!Fertilizer, Plant residues)
  !--------------------------------------------------------------
  SUBROUTINE crop_sources(i,j,dayno,calcN,calcP,calcC,ndays,ldays,areaij,thickness,soilstate,sources)

    USE MODVAR, ONLY : numsubstances, &
                       maxsoillayers, &
                       find_croppart, &
                       cropdata, &
                       i_in,i_on,i_sp,i_pp,i_oc, &
                       prevdoy, &
                       tsofday

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<subbasin
    INTEGER, INTENT(IN) :: j                        !<class
    INTEGER, INTENT(IN) :: dayno                    !<pseudo day number of the year
    LOGICAL, INTENT(IN) :: calcN                    !<Status of nitrogen simulation
    LOGICAL, INTENT(IN) :: calcP                    !<Status of phosphorus simulation
    LOGICAL, INTENT(IN) :: calcC                    !<Status of organic carbon simulation
    REAL, INTENT(IN)    :: ndays                    !<number of days to spread fertilizer
    REAL, INTENT(IN)    :: ldays                    !<number of days to spread plant residuals
    REAL, INTENT(IN)    :: areaij                   !<class area (km2)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers (m)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    REAL, INTENT(OUT)   :: sources(2,numsubstances) !<load from fertilizer and plant residues (kg/timestep)
    
    !Local variables
    INTEGER k,kcrop,isl
    INTEGER fertperiod,litterperiod
    REAL common_add(maxsoillayers,numsubstances)       !fert and manure, kg/km2
    REAL common_res(maxsoillayers,numsubstances)       !residuals, kg/km2
    REAL part                           !part of class area that is the current crop
    REAL littdays
    
    !Local parameters
    REAL, PARAMETER :: inorgpart = 0.5  !For manure IN-ON
    REAL, PARAMETER :: fastppart = 0.5  !For manure fastP-humusP

    !>\b Algorithm \n
    sources(:,:) = 0.0
    IF(.NOT.(calcN.OR.calcP.OR.calcC)) RETURN   !no substance modelled
    IF(tsofday /= 1) RETURN     !apply crop sources once per day (first time step)

    common_add = 0.
    common_res = 0.
    fertperiod = INT(ndays)
    litterperiod = INT(ldays)
    littdays = ldays
    !>For each crop on the class:
    DO kcrop = 1,2                             !main and secondary crop
      CALL find_croppart(i,j,kcrop,k,part)
      IF(k==0) CYCLE                            !no crop given, no uptake
      IF(part>0)THEN

        !>Calculate current anount of fertiliser and manure application
        IF((dayno >= cropdata(k)%fertday1 .AND. dayno < cropdata(k)%fertday1 + fertperiod) .OR. &
             (dayno < cropdata(k)%fertday1 + fertperiod - prevdoy))THEN
          IF(thickness(2)>0)THEN
            IF(calcN) common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%fertnamount1 * (1. - cropdata(k)%fertdown1) / ndays
            IF(calcP) common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%fertpamount1 * (1. - cropdata(k)%fertdown1) / ndays
            IF(calcN) common_add(2,i_in) = common_add(2,i_in) + part * cropdata(k)%fertnamount1 * cropdata(k)%fertdown1 / ndays
            IF(calcP) common_add(2,i_sp) = common_add(2,i_sp) + part * cropdata(k)%fertpamount1 * cropdata(k)%fertdown1 / ndays
          ELSE
            IF(calcN) common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%fertnamount1 / ndays
            IF(calcP) common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%fertpamount1 / ndays
          ENDIF
        ENDIF
        IF((dayno >= cropdata(k)%fertday2 .AND. dayno < cropdata(k)%fertday2 + fertperiod) .OR. &
             (dayno < cropdata(k)%fertday2 + fertperiod - prevdoy))THEN
          IF(thickness(2)>0)THEN
            IF(calcN) common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%fertnamount2 * (1.-cropdata(k)%fertdown2) / ndays
            IF(calcP) common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%fertpamount2 * (1.-cropdata(k)%fertdown2) / ndays
            IF(calcN) common_add(2,i_in) = common_add(2,i_in) + part * cropdata(k)%fertnamount2 * cropdata(k)%fertdown2 / ndays
            IF(calcP) common_add(2,i_sp) = common_add(2,i_sp) + part * cropdata(k)%fertpamount2 * cropdata(k)%fertdown2 / ndays
          ELSE
            IF(calcN) common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%fertnamount2 / ndays
            IF(calcP) common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%fertpamount2 / ndays
          ENDIF
        ENDIF
        IF((dayno >= cropdata(k)%manday1 .AND. dayno < cropdata(k)%manday1 + fertperiod) .OR.  &
             (dayno < cropdata(k)%manday1 + fertperiod - prevdoy))THEN
          IF(thickness(2)>0)THEN
            IF(calcN)THEN
              common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%mannamount1 * (1.-cropdata(k)%mandown1) * inorgpart / ndays
              common_add(1,i_on) = common_add(1,i_on) + part * cropdata(k)%mannamount1 * (1.-cropdata(k)%mandown1) * (1.-inorgpart) / ndays
              common_add(2,i_in) = common_add(2,i_in) + part * cropdata(k)%mannamount1 * cropdata(k)%mandown1 * inorgpart / ndays
              common_add(2,i_on) = common_add(2,i_on) + part * cropdata(k)%mannamount1 * cropdata(k)%mandown1 * (1.-inorgpart) / ndays
            ENDIF
            IF(calcP)THEN
              common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%manpamount1 * (1.-cropdata(k)%mandown1) * fastppart / ndays
              common_add(1,i_pp) = common_add(1,i_pp) + part * cropdata(k)%manpamount1 * (1.-cropdata(k)%mandown1) * (1.-fastppart) / ndays
              common_add(2,i_sp) = common_add(2,i_sp) + part * cropdata(k)%manpamount1 * cropdata(k)%mandown1 * fastppart / ndays
              common_add(2,i_pp) = common_add(2,i_pp) + part * cropdata(k)%manpamount1 * cropdata(k)%mandown1 * (1.-fastppart) / ndays
            ENDIF
          ELSE
            IF(calcN) common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%mannamount1 * inorgpart / ndays
            IF(calcP) common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%manpamount1 * fastppart / ndays
            IF(calcN) common_add(1,i_on) = common_add(1,i_on) + part * cropdata(k)%mannamount1 * (1.-inorgpart) / ndays
            IF(calcP) common_add(1,i_pp) = common_add(1,i_pp) + part * cropdata(k)%manpamount1 * (1.-fastppart) / ndays
          ENDIF
        ENDIF
        IF((dayno >= cropdata(k)%manday2 .AND. dayno < cropdata(k)%manday2 + fertperiod) .OR.  &
             (dayno < cropdata(k)%manday2 + fertperiod - prevdoy))THEN
          IF(thickness(2)>0)THEN
            IF(calcN)THEN
              common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%mannamount2 * (1.-cropdata(k)%mandown2) * inorgpart / ndays
              common_add(1,i_on) = common_add(1,i_on) + part * cropdata(k)%mannamount2 * (1.-cropdata(k)%mandown2) * (1.-inorgpart) / ndays
              common_add(2,i_in) = common_add(2,i_in) + part * cropdata(k)%mannamount2 * cropdata(k)%mandown2 * inorgpart / ndays
              common_add(2,i_on) = common_add(2,i_on) + part * cropdata(k)%mannamount2 * cropdata(k)%mandown2 * (1.-inorgpart) / ndays
            ENDIF
            IF(calcP)THEN
              common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%manpamount2 * (1.-cropdata(k)%mandown2) * fastppart / ndays
              common_add(1,i_pp) = common_add(1,i_pp) + part * cropdata(k)%manpamount2 * (1.-cropdata(k)%mandown2) * (1.-fastppart) / ndays
              common_add(2,i_sp) = common_add(2,i_sp) + part * cropdata(k)%manpamount2 * cropdata(k)%mandown2 * fastppart / ndays
              common_add(2,i_pp) = common_add(2,i_pp) + part * cropdata(k)%manpamount2 * cropdata(k)%mandown2 * (1.-fastppart) / ndays
            ENDIF
          ELSE
            IF(calcN) common_add(1,i_in) = common_add(1,i_in) + part * cropdata(k)%mannamount2 * inorgpart / ndays
            IF(calcP) common_add(1,i_sp) = common_add(1,i_sp) + part * cropdata(k)%manpamount2 * fastppart / ndays
            IF(calcN) common_add(1,i_on) = common_add(1,i_on) + part * cropdata(k)%mannamount2 * (1.-inorgpart) / ndays
            IF(calcP) common_add(1,i_pp) = common_add(1,i_pp) + part * cropdata(k)%manpamount2 * (1.-fastppart) / ndays
          ENDIF
        ENDIF
        !>Calculate current residue application
        IF(cropdata(k)%resdayno==0) littdays = 365.   !litterfall every day
        IF(cropdata(k)%resdayno==0 .OR.   &
             ((dayno >= cropdata(k)%resdayno .AND. dayno < cropdata(k)%resdayno + litterperiod) .OR.  &
             (dayno < cropdata(k)%resdayno + litterperiod - prevdoy)))THEN
          IF(thickness(2)>0)THEN
            IF(calcN)THEN
              common_res(1,i_in) = common_res(1,i_in) + cropdata(k)%resfast * part * cropdata(k)%resnamount * (1. - cropdata(k)%resdown) / littdays
              common_res(1,i_on) = common_res(1,i_on) + (1. - cropdata(k)%resfast) * part * cropdata(k)%resnamount * (1. - cropdata(k)%resdown) / littdays
              common_res(2,i_in) = common_res(2,i_in) + cropdata(k)%resfast * part * cropdata(k)%resnamount * cropdata(k)%resdown / littdays
              common_res(2,i_on) = common_res(2,i_on) + (1. - cropdata(k)%resfast) * part * cropdata(k)%resnamount * cropdata(k)%resdown / littdays
            ENDIF
            IF(calcP)THEN
              common_res(1,i_sp) = common_res(1,i_sp) + cropdata(k)%resfast * part * cropdata(k)%respamount * (1. - cropdata(k)%resdown) / littdays
              common_res(1,i_pp) = common_res(1,i_pp) + (1. - cropdata(k)%resfast) * part * cropdata(k)%respamount * (1. - cropdata(k)%resdown) / littdays
              common_res(2,i_sp) = common_res(2,i_sp) + cropdata(k)%resfast * part * cropdata(k)%respamount * cropdata(k)%resdown / littdays  
              common_res(2,i_pp) = common_res(2,i_pp) + (1. - cropdata(k)%resfast) * part * cropdata(k)%respamount * cropdata(k)%resdown / littdays
            ENDIF
            IF(calcC)THEN
              common_res(1,i_oc) = common_res(1,i_oc) + part * cropdata(k)%rescamount * (1. - cropdata(k)%resdown) / littdays
              common_res(2,i_oc) = common_res(2,i_oc) + part * cropdata(k)%rescamount * cropdata(k)%resdown / littdays
            ENDIF
          ELSE
            IF(calcN) common_res(1,i_in) = common_res(1,i_in) + cropdata(k)%resfast * part * cropdata(k)%resnamount / littdays
            IF(calcN) common_res(1,i_on) = common_res(1,i_on) + (1. - cropdata(k)%resfast) * part * cropdata(k)%resnamount / littdays
            IF(calcP) common_res(1,i_sp) = common_res(1,i_sp) + cropdata(k)%resfast * part * cropdata(k)%respamount / littdays
            IF(calcP) common_res(1,i_pp) = common_res(1,i_pp) + (1. - cropdata(k)%resfast) * part * cropdata(k)%respamount / littdays
            IF(calcC) common_res(1,i_oc) = common_res(1,i_oc) + part * cropdata(k)%rescamount / littdays 
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    !>Add fertiliser and manure to soil pools
    IF(calcN)THEN
      IF(SUM(common_add(:,i_in))>0.)THEN
        DO isl=1,2
          IF(soilstate%water(isl,j,i)>0.)THEN
            CALL add_source_to_water(soilstate%water(isl,j,i),1,soilstate%conc(i_in,isl,j,i),common_add(isl,i_in))
            sources(1,i_in) = sources(1,i_in) + common_add(isl,i_in)*areaij  !InorgN load from fertilizer from soil layer 1 (kg)
          ELSE
            CALL production_pool(1,soilstate%fastN(isl,j,i),common_add(isl,i_in))
            sources(1,i_on) = sources(1,i_on) + common_add(isl,i_in)*areaij  !OrgN load from fertilizer from soil layer 1 (kg)
          ENDIF
        ENDDO
      ENDIF
      IF(SUM(common_add(:,i_on))>0.)THEN
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%fastN(:,j,i),common_add(:,i_on))
          sources(1,i_on) = sources(1,i_on) + SUM(common_add(:,i_on))*areaij  !OrgN load from fertilizer, both soil layers (kg)
        ELSE
          CALL production_pool(1,soilstate%fastN(1,j,i),common_add(1,i_on))
          sources(1,i_on) = sources(1,i_on) + common_add(1,i_on)*areaij       !OrgN load from fertilizer, soil layer 1 (kg)
        ENDIF
      ENDIF
    ENDIF
    IF(calcP)THEN
      IF(SUM(common_add(:,i_sp))>0.)THEN
        DO isl = 1,2
          IF(soilstate%water(isl,j,i)>0.)THEN
            CALL add_source_to_water(soilstate%water(isl,j,i),1,soilstate%conc(i_sp,isl,j,i),common_add(isl,i_sp))
            sources(1,i_sp) = sources(1,i_sp) + common_add(isl,i_sp)*areaij        !SRP load from fertilizer, soil layer 1 (kg)
          ELSE
            CALL production_pool(1,soilstate%fastP(isl,j,i),common_add(isl,i_sp))
            sources(1,i_pp) = sources(1,i_pp) + common_add(isl,i_sp)*areaij        !PartP load from fertilizer, soil layer 1 (kg)
          ENDIF
        ENDDO
      ENDIF
      IF(SUM(common_add(:,i_pp))>0.)THEN
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%fastP(:,j,i),common_add(:,i_pp))
          sources(1,i_pp) = sources(1,i_pp) + SUM(common_add(:,i_pp))*areaij  !PartP load from fertilizer, both soil layers (kg)
        ELSE
          CALL production_pool(1,soilstate%fastP(1,j,i),common_add(1,i_pp))
          sources(1,i_pp) = sources(1,i_pp) + common_add(1,i_pp)*areaij        !PartP load from fertilizer, soil layer 1 (kg)
        ENDIF
      ENDIF
    ENDIF

    !>Add residues to soil pools
    IF(calcN)THEN
      IF(SUM(common_res(:,i_in))>0.)THEN
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%fastN(:,j,i),common_res(:,i_in))
        ELSE
          CALL production_pool(1,soilstate%fastN(1,j,i),common_res(1,i_in))
        ENDIF
        sources(2,i_on) = sources(2,i_on) + SUM(common_res(:,i_in))*areaij  !OrgN load from residuals, all soil layers, fastN type (kg)
      ENDIF
      IF(SUM(common_res(:,i_on))>0.)THEN
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%humusN(:,j,i),common_res(:,i_on))
        ELSE
          CALL production_pool(1,soilstate%humusN(1,j,i),common_res(1,i_on))
        ENDIF
        sources(2,i_on) = sources(2,i_on) + SUM(common_res(:,i_on))*areaij  !OrgN load from residuals, all soil layers, humusN (kg)
      ENDIF
    ENDIF
    IF(calcP)THEN
      IF(SUM(common_res(:,i_sp))>0.)THEN                
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%fastP(:,j,i),common_res(:,i_sp))     
        ELSE
          CALL production_pool(1,soilstate%fastP(1,j,i),common_res(1,i_sp))                
        ENDIF
        sources(2,i_pp) = sources(2,i_pp) + SUM(common_res(:,i_sp))*areaij                !PartP load from residuals, all soil layers, fastP type (kg)
      ENDIF
      IF(SUM(common_res(:,i_pp))>0.)THEN                
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%humusP(:,j,i),common_res(:,i_pp))     
        ELSE
          CALL production_pool(1,soilstate%humusP(1,j,i),common_res(1,i_pp))       
        ENDIF
        sources(2,i_pp) = sources(2,i_pp) + SUM(common_res(:,i_pp))*areaij                !PartP load from residuals, all soil layers, humusP type (kg)
      ENDIF
    ENDIF
    IF(calcC)THEN
      IF(SUM(common_res(:,i_oc))>0.)THEN                
        IF(thickness(2)>0)THEN
          CALL production_pool(maxsoillayers,soilstate%fastC(:,j,i),common_res(:,i_oc))     
        ELSE
          CALL production_pool(1,soilstate%fastC(1,j,i),common_res(1,i_oc))     
        ENDIF
        sources(2,i_oc) = sources(2,i_oc) + SUM(common_res(:,i_oc))*areaij                !OrgC load from residuals, all soil layers (kg)
      ENDIF
    ENDIF

  END SUBROUTINE crop_sources

  !>\brief Erosion of soil particles and phosphorus with fast flow, surface flow and 
  !!macropore flow. Delay of particles in temporary storage. The 
  !!released P is transported as PP by all runoff flows. Both partP and humusP of the soil
  !!contribute to eroded P. 
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Vegetation 
  !! and soil surface processes - Erosion calculations)
  !---------------------------------------------------------------------------------------
  SUBROUTINE particle_processes_for_runoff(i,j,isoil,iluse,dayno,prec,surfacerunoff,macroporflow,    &
                          tilerunoff,totflow,phoscorr,csurface,cmacro,ctile,crunoff1,crunoff2,  &
                          crunoff3,snow,soilstate)

    USE MODVAR, ONLY : soilthick, &
                       basin, &
                       modeloption, &
                       p_erosion, &
                       numsubstances, &
                       i_pp,i_ss, &
                       genpar,landpar,soilpar
    USE HYPEVARIABLES, ONLY : m_soilcoh,m_soilerod,m_sreroexp,           &
                              m_filtPbuf,m_filtPinner,m_filtPother,m_macfilt, &
                              m_pprelmax,m_pprelexp,m_erodluse,m_erodsoil, &
                              m_erodslope,m_erodexp,m_erodindex

    !Argument declarations
    INTEGER, INTENT(IN) :: i         !<index of current subbasin
    INTEGER, INTENT(IN) :: j         !<index of current class
    INTEGER, INTENT(IN) :: isoil     !<index of soil type
    INTEGER, INTENT(IN) :: iluse     !<index of landuse
    INTEGER, INTENT(IN) :: dayno     !<pseudo day number of the year
    REAL, INTENT(IN)    :: prec      !<precipitation (rainfall only)    
    REAL, INTENT(IN)    :: surfacerunoff   !<surface runoff (mm/timestep)
    REAL, INTENT(IN)    :: macroporflow    !<macropore flow (mm/timestep)
    REAL, INTENT(IN)    :: tilerunoff      !<tile drainage runoff (mm/timestep)
    REAL, INTENT(IN)    :: totflow     !<total runoff (surfacerunoff, tilerunoff, soilrunoff layer 1-3)
    REAL, INTENT(IN)    :: phoscorr    !<correction of phosphorus level
    REAL, INTENT(INOUT) :: csurface(numsubstances)    !<PP concentration surface runoff (mg/L)
    REAL, INTENT(INOUT) :: cmacro(numsubstances)      !<PP concentration macropores (mg/L)
    REAL, INTENT(INOUT) :: ctile(numsubstances)       !<PP concentration tile drainage (mg/L)
    REAL, INTENT(INOUT) :: crunoff1(numsubstances)    !<PP concentration soil runoff layer 1 (mg/L)
    REAL, INTENT(INOUT) :: crunoff2(numsubstances)    !<PP concentration soil runoff layer 2 (mg/L)
    REAL, INTENT(INOUT) :: crunoff3(numsubstances)    !<PP concentration soil runoff layer 3 (mg/L)
    REAL,INTENT(IN)     :: snow        !<snow water (mm)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states

    !Local variables
    REAL PPrel,SSrel          !PP/SS released from delay pool
    REAL srfilt               !total filtration of surface runoff PP
    REAL erodingflow          !Flow eroding the surface
    REAL erodedsed            !eroded (transported) sediment (kg/km2)
    REAL erodedP              !total eroded PP (kg/km2)
    REAL fracminP             !fraction of erodedP in mineral form (from partP)     
    REAL removePP(1)          !actually lost PP (kg/km2)
    REAL removeHP(1)          !actually lost humusP (kg/km2)            
    REAL macroppconc,surfrppconc  !variables to temporary hold PP-konc of macropore flow and surfacerunoff 
    REAL macrossconc,surfrssconc  !variables to temporary hold SS-konc of macropore flow and surfacerunoff 
    REAL newconc                !concentration for released material from delay pool

    !>\b Algorithm \n
    !>Calculate erosion of soil and transport of sediment by runoff
    macroppconc = 0; surfrppconc = 0.
    macrossconc = 0; surfrssconc = 0.
    erodingflow = surfacerunoff + macroporflow
    SELECT CASE(modeloption(p_erosion))
    CASE(0)
      CALL calculate_erosion(i,j,dayno,prec,surfacerunoff, &   !calculate eroded sediment
             soilpar(m_soilcoh,isoil),phoscorr*soilpar(m_soilerod,isoil), &
             snow,genpar(m_sreroexp),erodingflow,erodedsed)
    CASE(1)
      CALL calculate_hbvsed_erosion(i,prec,landpar(m_erodluse,iluse),soilpar(m_erodsoil,isoil), &   !calculate eroded sediment basedon HBV-sed
                      genpar(m_erodslope),genpar(m_erodexp),genpar(m_erodindex),erodedsed)
    END SELECT
    !>Calculate particulate phosphorus of eroded sediment
    IF(i_pp>0) CALL calculate_pp_transport(erodingflow,erodedsed,soilstate%partP(1,j,i),soilstate%humusP(1,j,i),soilthick(1,j),fracminP,erodedP)

    !>The eroded sediment is divided between surface runoff and macropore flow (if positive flows) 
    !>and reduced by filtering/deposition on the way to the stream
    IF(surfacerunoff>0.) THEN     !erodedP goes back to soil if no surface runoff
      srfilt = landpar(m_filtPother,iluse) + basin(i)%closewater * (1. + basin(i)%buffer * (landpar(m_filtPbuf,iluse) - 1.)) + landpar(m_filtPinner,iluse) * (1. - basin(i)%closewater)
      IF(i_pp>0) surfrppconc = srfilt * erodedP / erodingflow
      IF(i_ss>0) surfrssconc = srfilt * erodedsed / erodingflow
      IF(tilerunoff>0.) THEN       !only when runoffd>0  will PP from macropor reach runoff
        IF(macroporflow>0. )THEN
          IF(i_pp>0) macroppconc =  soilpar(m_macfilt,isoil) * erodedP / erodingflow
          IF(i_ss>0) macrossconc =  soilpar(m_macfilt,isoil) * erodedsed / erodingflow
        ENDIF
      ENDIF
      IF(i_pp>0)THEN
        !>The final amount eroded P is removed from soil pools
        removePP(1) = surfacerunoff*surfrppconc*fracminP+macroporflow*macroppconc*fracminP !kg/km2                 
        removeHP(1) = surfacerunoff*surfrppconc*(1-fracminP)+macroporflow*macroppconc*(1-fracminP) !kg/km2        
        CALL retention_pool(1,soilstate%partP(1,j,i),removePP)
        CALL retention_pool(1,soilstate%humusP(1,j,i),removeHP)
      ENDIF
    ENDIF

    !>Sediment and PP of macropore flow is 100% filtered and all added to the temporary pools
    IF(macroporflow>0.) THEN
      IF(i_pp>0)THEN
        macroppconc = macroppconc + cmacro(i_pp)    
        cmacro(i_pp) = 0.
        soilstate%PPrelpool(j,i) = soilstate%PPrelpool(j,i) + macroppconc*macroporflow
      ENDIF
      IF(i_ss>0)THEN
        macrossconc = macrossconc + cmacro(i_ss)    
        cmacro(i_ss) = 0.
        soilstate%Srelpool(j,i) = soilstate%Srelpool(j,i) + macrossconc*macroporflow
      ENDIF
    ENDIF

    !Temporary pool delays surface runoff eroded PP and tile runoff PP. 
    !The phosphorus/sediment is released during high total runoff equally distributed over all runoff paths.
    IF(totflow>0.)THEN
      IF(i_pp>0)THEN
        csurface(i_pp) = csurface(i_pp) + surfrppconc   !Add eroded PP to surface runoff concentration 
        !>Add PP of tile drainage and surface flow to temporary PP pool 
        soilstate%PPrelpool(j,i) = soilstate%PPrelpool(j,i) + ctile(i_pp) * tilerunoff + csurface(i_pp) * surfacerunoff   !add sources
        !>Calculate release from temporary PP pool and new concentrations of the flows
        PPrel = MIN(soilstate%PPrelpool(j,i), soilstate%PPrelpool(j,i)* ((totflow  / genpar(m_pprelmax))**genpar(m_pprelexp)))    !export
        soilstate%PPrelpool(j,i) = soilstate%PPrelpool(j,i) - PPrel
        newconc = PPrel / totflow
        csurface(i_pp) = newconc
        ctile(i_pp) = newconc
        crunoff1(i_pp) = crunoff1(i_pp) + newconc
        crunoff2(i_pp) = crunoff2(i_pp) + newconc
        crunoff3(i_pp) = crunoff3(i_pp) + newconc
      ENDIF
      IF(i_ss>0)THEN
        csurface(i_ss) = csurface(i_ss) + surfrssconc   !Add eroded sed to surface runoff concentration 
        !>Add SS of tile drainage and surface flow to temporary pool 
        soilstate%Srelpool(j,i) = soilstate%Srelpool(j,i) + ctile(i_ss) * tilerunoff + csurface(i_ss) * surfacerunoff
        !>Calculate release from temporary pool and new concentrations of the flows
        SSrel = MIN(soilstate%Srelpool(j,i), soilstate%Srelpool(j,i)* ((totflow  / genpar(m_pprelmax))**genpar(m_pprelexp)))    !export
        soilstate%Srelpool(j,i) = soilstate%Srelpool(j,i) - SSrel
        newconc = SSrel / totflow
        csurface(i_ss) = newconc
        ctile(i_ss) = newconc
        crunoff1(i_ss) = crunoff1(i_ss) + newconc
        crunoff2(i_ss) = crunoff2(i_ss) + newconc
        crunoff3(i_ss) = crunoff3(i_ss) + newconc
      ENDIF
    ENDIF

  END SUBROUTINE particle_processes_for_runoff

  !>Calculate eroded particles from soil and attached P
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Vegetation 
  !! and soil surface processes - Erosion calculations, Crop cover and ground cover)
  !------------------------------------------------------------------
  SUBROUTINE calculate_erosion(i,j,dayno,prec,surfacerunoff,cohesion, &
                               erodibility,snow,sreroexp,flow, &
                               erodedsed)

    USE MODVAR, ONLY : find_croppart, &
                       cropdata, &
                       basin, &
                       pi                              

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    INTEGER, INTENT(IN) :: j             !<index of current class
    INTEGER, INTENT(IN) :: dayno         !<current pseudo day number of the year
    REAL, INTENT(IN)    :: prec          !<precipitation (rainfall only)    
    REAL, INTENT(IN)    :: surfacerunoff !<saturated overland flow and excess infiltration (mm)
    REAL, INTENT(IN)    :: cohesion      !<(kPa)
    REAL, INTENT(IN)    :: erodibility   !<(g/J)
    REAL, INTENT(IN)    :: snow          !<snow water (mm)
    REAL, INTENT(IN)    :: sreroexp      !<surface runoff erosion exponent 
    REAL, INTENT(IN)    :: flow          !<fast flow
    REAL, INTENT(OUT)   :: erodedsed     !<eroded (transported) sediment (kg/km2)

    !Local variables
    INTEGER k,kcrop
    REAL Rainfall_energy, cropcover, groundcover
    REAL intensity
    REAL part
    REAL bd1, bd2, bd3, bd4, bd5  !cultivation dates
    REAL common_cropcover, common_groundcover,maxday1, maxday2
    REAL transportfactor
    REAL mobilisedsed         !mobilised suspended sediment from rainfall and surface runoff (g/m2)
    
    !Local parameters
    REAL, PARAMETER :: trans1 = 4.0
    REAL, PARAMETER :: trans2 = 1.3  

    !>/b Algorithm
    mobilisedsed = 0.
    erodedsed = 0.
    common_cropcover = 1.
    common_groundcover = 1.
    IF(cohesion==0.OR.erodibility==0) RETURN      !no parameter values -> no erosion

    !>Calculate current cropcover and groundcover, will limit erosion
    DO kcrop = 1, 2
      CALL find_croppart(i,j,kcrop,k,part)
      IF(k==0) CYCLE                              !no crop given
      IF(part>0)THEN
        bd1 = cropdata(k)%baredayno1        !spring ploughing
        bd2 = cropdata(k)%baredayno2        !sow date / beginning of growing season 
        bd3 = cropdata(k)%baredayno3        !harvest
        bd4 = cropdata(k)%baredayno4        !autumn ploughing
        bd5 = cropdata(k)%baredayno5        !winter crops sowing date
        IF(bd1==0 .AND. bd4==0) THEN                   !year-round-crop
          cropcover = cropdata(k)%ccmax1
          groundcover = cropdata(k)%gcmax1
        ELSE                                           !spring, winter and row crops
          maxday1 = REAL(INT(bd2 + (bd3-bd2)/2.))   !day of maximum crop and ground cover in summer
          maxday2 = REAL(INT(bd5 + (365-bd5)/2.))   !day of maximum crop and ground cover for winter crops in autumn
          IF(bd5 > 0)THEN                              !winter crop
            IF(dayno < bd2)THEN
              cropcover   = cropdata(k)%ccmax2
              groundcover = cropdata(k)%gcmax2
            ELSEIF(dayno < maxday1)THEN
              cropcover   = cropdata(k)%ccmax2 + (cropdata(k)%ccmax1-cropdata(k)%ccmax2)*((dayno - bd2) /(maxday1 - bd2)) 
              groundcover = cropdata(k)%gcmax2 + (cropdata(k)%gcmax1-cropdata(k)%gcmax2)*((dayno - bd2) /(maxday1 - bd2))
            ELSEIF(dayno < bd3)THEN
              cropcover   = cropdata(k)%ccmax1
              groundcover = cropdata(k)%gcmax1
            ELSEIF(dayno < bd4)THEN
              cropcover   = cropdata(k)%gcmax1
              groundcover = cropdata(k)%gcmax1    
            ELSEIF(dayno < bd5)THEN
              cropcover = 0.
              groundcover = 0.
            ELSEIF(dayno < maxday2)THEN
              cropcover   = cropdata(k)%ccmax2 * ((dayno - bd5) /(maxday2 - bd5))
              groundcover = cropdata(k)%gcmax2 * ((dayno - bd5) /(maxday2 - bd5))
            ELSE
              cropcover   = cropdata(k)%ccmax2
              groundcover = cropdata(k)%gcmax2
            ENDIF
          ELSEIF(bd1 > 0)THEN                          !spring crop (or row crop) with spring ploughing
            IF(dayno < bd1)THEN
              cropcover   = cropdata(k)%gcmax1
              groundcover = cropdata(k)%gcmax1             
            ELSEIF(dayno < bd2)THEN
              cropcover   = 0.
              groundcover = 0.
            ELSEIF(dayno < maxday1) THEN
              cropcover   = cropdata(k)%ccmax1 * ((dayno - bd2) /(maxday1 - bd2)) 
              groundcover = cropdata(k)%gcmax1 * ((dayno - bd2) /(maxday1 - bd2))
            ELSEIF(dayno < bd3) THEN
              cropcover   = cropdata(k)%ccmax1
              groundcover = cropdata(k)%gcmax1
            ELSE
              cropcover   = cropdata(k)%gcmax1      
              groundcover = cropdata(k)%gcmax1    
            ENDIF
          ELSE                                         !spring crop (or row crop) with autumn ploughing
            IF(dayno < bd2)THEN
              cropcover   = 0.
              groundcover = 0.
            ELSEIF(dayno < maxday1) THEN
              cropcover   = cropdata(k)%ccmax1 * ((dayno - bd2) /(maxday1 - bd2)) 
              groundcover = cropdata(k)%gcmax1 * ((dayno - bd2) /(maxday1 - bd2))
            ELSEIF(dayno < bd3) THEN        
              cropcover   = cropdata(k)%ccmax1
              groundcover = cropdata(k)%gcmax1
            ELSEIF(dayno < bd4) THEN        
              cropcover   = cropdata(k)%gcmax1
              groundcover = cropdata(k)%gcmax1    
            ELSE        
              cropcover   = 0.
              groundcover = 0.
            ENDIF
          ENDIF
        ENDIF
        common_cropcover   = common_cropcover * (1. - part * cropcover)
        common_groundcover = common_groundcover * (1. - part * groundcover)
      ENDIF !part>0
    ENDDO
    common_cropcover   = 1. - common_cropcover
    common_groundcover = 1. - common_groundcover

    !Check for snow limiting erosion
    intensity = 1. !intenspar
    IF(snow>0.) intensity = 0.    !snow

    !>Calculate particles that is eroded by rain splash detachment and by overland flow (mobilised sediment)
    IF(prec > 0.) THEN
      IF(intensity > 0.) THEN
        IF(prec>5.0) THEN     !TODO: shorter timestep, other threshold?, holds for all over the world?, reference?
          Rainfall_energy = 8.95+8.44*LOG10(prec*(0.257+sin(2*3.14*((dayno-70.)/365.))*0.09)*2.)
        ELSE
          Rainfall_energy = 0.
        ENDIF
        Rainfall_energy = prec * Rainfall_energy        !J/m2
        mobilisedsed = Rainfall_energy * (1. - common_cropcover) * erodibility  !g/m2
      ENDIF
    ENDIF
    IF(surfacerunoff > 0.) THEN   
       mobilisedsed = mobilisedsed + (((surfacerunoff * 365.) ** sreroexp) * (1. - common_groundcover) * (1./(0.5 * cohesion)) * SIN(basin(i)%slope / 100.)) / 365. !g/m2   
    ENDIF

    !>Transport capacity of fast flowing water may limit transport of sediment
    IF(flow>0.)THEN
      transportfactor = MIN(1.,(flow / trans1)**trans2)   
    ELSE
      transportfactor = 1.
    ENDIF
      
    !>Eroded sediment calculated from mobilised sediment, possibly limited by the transport capacity
    erodedsed = 1000. * mobilisedsed * transportfactor  !kg/km2

  END SUBROUTINE calculate_erosion

  !>Calculate eroded particles from soil by HBV-sed based model
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Vegetation 
  !! and soil surface processes - Erosion calculations)
  !------------------------------------------------------------------
  SUBROUTINE calculate_hbvsed_erosion(i,prec,lusepar,soilpar,slopepar, &
                                      precexppar,eroindexpar,erodedsed)

    USE MODVAR, ONLY : basin

    !Argument declarations
    INTEGER, INTENT(IN) :: i             !<index of current subbasin
    REAL, INTENT(IN)    :: prec          !<precipitation (rainfall only)    
    REAL, INTENT(IN)    :: lusepar       !<soil erosion factor (land use dependence)
    REAL, INTENT(IN)    :: soilpar       !<soil erosion factor (soil dependence)
    REAL, INTENT(IN)    :: slopepar      !<slope erosion factor (exponent) 
    REAL, INTENT(IN)    :: precexppar    !<erosion precipitation dependence factor (exponent)
    REAL, INTENT(IN)    :: eroindexpar   !<model parameter for scaling of erosion index
    REAL, INTENT(OUT)   :: erodedsed     !<eroded (transported) sediment (kg/km2)

    !Local variables
    REAL erosionindex
    REAL mobilisedsed   !mobilised suspended sediment from rainfall (g/m2)
    
    !>/b Algorithm
    !>Calculate mobilised sediments from erosion index, slope and rainfall
    erosionindex = basin(i)%eroindex / eroindexpar    
    mobilisedsed = (basin(i)%slope / 5.)**slopepar * lusepar * soilpar * erosionindex * prec**precexppar       !tonnes/km2 = g/m2
      
    !>Eroded sediment calculated from mobilised sediment
    erodedsed = 1000. * mobilisedsed   !kg/km2

  END SUBROUTINE calculate_hbvsed_erosion

  !> Calculate how much particulate phosphorus is transported
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Vegetation 
  !! and soil surface processes - Erosion calculations)
  !------------------------------------------------------------------
  SUBROUTINE calculate_pp_transport(flow,erodedsed,soilPartP,soilHumusP,thickness,fracminP,erodedP)
    
    USE HYPEVARIABLES, ONLY : bulkdensity
    
    !Argument declarations
    REAL, INTENT(IN)  :: flow        !<fast flow
    REAL, INTENT(IN)  :: erodedsed   !<mobilised sediment (kg/km2)
    REAL, INTENT(IN)  :: soilPartP   !<partP upper soil layer (kg/m2)
    REAL, INTENT(IN)  :: soilHumusP  !<humusP upper soil layer (kg/m2)
    REAL, INTENT(IN)  :: thickness   !<upper soillayer thickness (m)    
    REAL, INTENT(OUT) :: fracminP    !<fraction of eroded P that is mineral
    REAL, INTENT(OUT) :: erodedP     !<eroded PP (kg/m2)
    
    !Local variables
    REAL enrichment

    !Local parameters
    REAL, PARAMETER :: maxenrichment = 4.  
    REAL, PARAMETER :: stab = 1.5
    REAL, PARAMETER :: flowstab = 4.

    !>\b Algorithm
    !>Calculate enrichment factor - finer soil particles are more likely to be eroded and contain more P per weight unit
    IF(flow>0.)THEN     
      IF(flow>flowstab)THEN 
        enrichment = stab
      ELSE 
        enrichment = maxenrichment - (maxenrichment-stab)* flow / flowstab
      ENDIF
    ELSE
      enrichment = 0.
    ENDIF    
     
    !>Calculate erodedP as P in eroded sediments and effect of enrichment
    erodedP = 1.0E-6 * erodedsed * ((soilPartP + soilHumusP) / (thickness * bulkdensity)) * enrichment
    !>Calculate mineral fraction of erodedP
    fracminP = soilPartP / (soilPartP + soilHumusP)
  
  END SUBROUTINE calculate_pp_transport
  
  !>Calculate nutrient processes in soil; sources, transformations,
  !and retention
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes)
  !------------------------------------------------------------------------------
  SUBROUTINE soil_np_processes(i,j,iluse,calcN,calcP,calcC,dayno,area,wp,fc,ep, &
                      common_uptake,thickness,ndays,ldays,source,sink,nitrification,  &
                      denitrification,cropuptake,cropsources,pardisfN,pardisfP, &
                      pardishN,pardishP,parminfN,parminfP,pardegrhN,denpar1,  &
                      pardegrhP,soilstate)       

    USE HYPEVARIABLES, ONLY : m_denitr3
    USE MODVAR, ONLY : numsubstances,   &
                       maxsoillayers,   &
                       i_in

    INTEGER, INTENT(IN) :: i                        !<index of subbasin 
    INTEGER, INTENT(IN) :: j                        !<index of class
    INTEGER, INTENT(IN) :: iluse                    !<index of land use
    LOGICAL, INTENT(IN) :: calcN                    !<Status of nitrogen simulation
    LOGICAL, INTENT(IN) :: calcP                    !<Status of phosphorus simulation
    LOGICAL, INTENT(IN) :: calcC                    !<Status of organic carbon simulation
    INTEGER, INTENT(IN) :: dayno                    !<pseudo day number of the year
    REAL, INTENT(IN)    :: area                     !<class area (km2)
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<water content at wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<water content at field capacity (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers)        !<water content in effective porosity (mm)
    REAL, INTENT(IN)    :: common_uptake(2,2)       !<uptake by plants (kg NP/km2/timestep)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(IN)    :: ndays                    !<model parameter, days for fertilization
    REAL, INTENT(IN)    :: ldays                    !<model parameter, days for litterfall
    REAL, INTENT(OUT)   :: source(numsubstances)    !<source of soil water nutrients through mineralization
    REAL, INTENT(OUT)   :: sink(numsubstances)      !<OBS: denitrification sink only
    REAL, INTENT(OUT)   :: nitrification            !<nitrification
    REAL, INTENT(OUT)   :: denitrification(maxsoillayers)  !<denitrification
    REAL, INTENT(OUT)   :: cropuptake               !<crop uptake of IN
    REAL, INTENT(OUT)   :: cropsources(2,numsubstances)    !<loads of fertilizer and residuals (kg/timestep)
    REAL, INTENT(IN)    :: pardisfN            !<model parameter mineralisation fastN to dissolved ON
    REAL, INTENT(IN)    :: pardisfP            !<model parameter mineralisation fastP to dissolved OP
    REAL, INTENT(IN)    :: pardishN            !<model parameter mineralisation humusN to dissolved ON
    REAL, INTENT(IN)    :: pardishP            !<model parameter mineralisation humusP to dissolved OP
    REAL, INTENT(IN)    :: parminfN            !<model parameter mineralisation fastN
    REAL, INTENT(IN)    :: parminfP            !<model parameter mineralisation fastP
    REAL, INTENT(IN)    :: pardegrhN           !<model parameter degradation humusN
    REAL, INTENT(IN)    :: denpar1             !<generally adjusted denitrification rate
    REAL, INTENT(IN)    :: pardegrhP           !<model parameter degradation humusP          
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states

    !No substances modelled
    IF (numsubstances == 0) RETURN

    !Calculate the nutrient processes
    CALL crop_sources(i,j,dayno,calcN,calcP,calcC,ndays,ldays,area,thickness,soilstate,cropsources)                                              !fertilizer and residues, also orgC
    CALL soil_pool_transformations(i,j,calcN,calcP,wp,fc,ep,thickness,source,pardisfN,  &
         pardisfP,pardishN,pardishP,parminfN,parminfP,pardegrhN,pardegrhP,  &
         soilstate)   !transformation from humusN->fastN->IN, ON->IN and fastP->SRP humusP -> fastP  
    IF(i_in>0) nitrification = source(i_in)
    CALL plant_uptake(calcN,calcP,common_uptake,wp,thickness,soilstate%water(:,j,i),soilstate%conc(:,:,j,i),sink)    
    IF(i_in>0) cropuptake = sink(i_in)
    CALL calculate_soil_profile_denitrification(calcN,i,j,iluse,wp,fc,ep,thickness,denpar1,soilstate,denitrification)

  END SUBROUTINE soil_np_processes

  !>\brief Transformation between NP-pools in soil. Degradation/delay
  !functions.
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes)
  !------------------------------------------------------------------------------
  SUBROUTINE soil_pool_transformations(i,j,calcN,calcP,wp,fc,ep,thickness,source,        &
       pardisfN,pardisfP,pardishN,pardishP,  &
       minfNpar,minfPpar,degrhNpar,degrhPpar,soilstate)        

    USE MODVAR, ONLY : &
         numsubstances,maxsoillayers,          &
         i_in,i_on,i_sp,i_pp
    USE HYPEVARIABLES, ONLY : satact,Thetapow,Thetalow,Thetaupp

    !Argument declarations
    INTEGER, INTENT(IN) :: i                        !<current subbasin
    INTEGER, INTENT(IN) :: j                        !<current slc-class
    LOGICAL, INTENT(IN) :: calcN    !<flag for nitrogen simulation
    LOGICAL, INTENT(IN) :: calcP    !<flag for phosphorus simulation
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<water content at wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<water content at field capacity (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers)        !<water content: effectiv porosity (mm)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(OUT)   :: source(numsubstances)    !<source of soil water nutrients through mineralization
    REAL, INTENT(IN)    :: pardisfN            !<model parameter mineralisation fastN to dissolved ON
    REAL, INTENT(IN)    :: pardisfP            !<model parameter mineralisation fastP to dissolved OP
    REAL, INTENT(IN)    :: pardishN            !<model parameter mineralisation humusN to dissolved ON
    REAL, INTENT(IN)    :: pardishP            !<model parameter mineralisation humusP to dissolved OP
    REAL, INTENT(IN)    :: minfNpar            !<model parameter mineralisation fastN
    REAL, INTENT(IN)    :: minfPpar            !<model parameter mineralisation fastP
    REAL, INTENT(IN)    :: degrhNpar           !<model parameter degradation humusN
    REAL, INTENT(IN)    :: degrhPpar           !<model parameter degradation humusP          
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states

    !Local variables
    INTEGER k   !soillayer
    REAL inorgNpool(maxsoillayers),SRPpool(maxsoillayers) !mg/L*mm=kg/km2
    REAL DONpool(maxsoillayers),DOPpool(maxsoillayers) !mg/L*mm=kg/km2
    REAL degradhN(maxsoillayers),transfN(maxsoillayers)
    REAL transfP(maxsoillayers),degradhP(maxsoillayers)       
    REAL dissolfN(maxsoillayers),dissolfP(maxsoillayers)       
    REAL dissolhN(maxsoillayers),dissolhP(maxsoillayers)       
    REAL tmpfcn(maxsoillayers)
    REAL sm(maxsoillayers), smfcn(maxsoillayers)

    source = 0.

    !Current pools of N and P dissolved in soil water
    IF(calcN)THEN
      inorgNpool(1) = soilstate%conc(i_in,1,j,i) * soilstate%water(1,j,i)
      DONpool(1) = soilstate%conc(i_on,1,j,i) * soilstate%water(1,j,i)
      IF(thickness(2)>0)THEN
        inorgNpool(2) = soilstate%conc(i_in,2,j,i) * soilstate%water(2,j,i)
        DONpool(2) = soilstate%conc(i_on,2,j,i) * soilstate%water(2,j,i)
      ENDIF
      IF(thickness(3)>0)THEN
        inorgNpool(3) = soilstate%conc(i_in,3,j,i) * soilstate%water(3,j,i)
        DONpool(3) = soilstate%conc(i_on,3,j,i) * soilstate%water(3,j,i)
      ENDIF
    ENDIF
    IF(calcP)THEN
      SRPpool(1) = soilstate%conc(i_sp,1,j,i) * soilstate%water(1,j,i)
      DOPpool(1) = soilstate%conc(i_pp,1,j,i) * soilstate%water(1,j,i)
      IF(thickness(2)>0)THEN
        SRPpool(2) = soilstate%conc(i_sp,2,j,i) * soilstate%water(2,j,i)
        DOPpool(2) = soilstate%conc(i_pp,2,j,i) * soilstate%water(2,j,i)
      ENDIF
      IF(thickness(3)>0)THEN
        SRPpool(3) = soilstate%conc(i_sp,3,j,i) * soilstate%water(3,j,i)
        DOPpool(3) = soilstate%conc(i_pp,3,j,i) * soilstate%water(3,j,i)
      ENDIF
    ENDIF

    !Temperature dependence factor
    DO k = 1,maxsoillayers
      IF(thickness(k)>0) tmpfcn(k) = tempfactor(soilstate%temp(k,j,i))
    ENDDO

    !Soil moisture dependence factor
    sm(:) = soilstate%water(:,j,i)
    smfcn(1) = moisturefactor(sm(1),wp(1),wp(1)+fc(1)+ep(1),thickness(1),satact,thetapow,thetalow,thetaupp)
    smfcn(2) = moisturefactor(sm(2),wp(2),wp(2)+fc(2)+ep(2),thickness(2),satact,thetapow,thetalow,thetaupp)
    smfcn(3) = moisturefactor(sm(3),wp(3),wp(3)+fc(3)+ep(3),thickness(3),satact,thetapow,thetalow,thetaupp)

    !Degradation of humusN to fastN
    degradhN = 0.
    IF(calcN)THEN
      degradhN(:) = degrhNpar * tmpfcn(:) * smfcn(:) * soilstate%humusN(:,j,i)
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%humusN(:,j,i),degradhN(:))      !may change degradhN
        CALL production_pool(maxsoillayers,soilstate%fastN(:,j,i),degradhN(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%humusN(1:2,j,i),degradhN(1:2))      !may change degradhN
        CALL production_pool(2,soilstate%fastN(1:2,j,i),degradhN(1:2))
      ELSE
        CALL retention_pool(1,soilstate%humusN(1,j,i),degradhN(1))      !may change degradhN
        CALL production_pool(1,soilstate%fastN(1,j,i),degradhN(1))
      ENDIF
    ENDIF

    !Transformation of fastN to inorganic N
    transfN = 0.
    IF(calcN)THEN
      transfN(:) = minfNpar * tmpfcn(:) * smfcn(:) * soilstate%fastN(:,j,i)     !kanske inte gödseltillskott borde bero på Temp??
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%fastN(:,j,i),transfN(:)) !transfN may change in retention_pool
        CALL production_pool(maxsoillayers,inorgNpool(:),transfN(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%fastN(1:2,j,i),transfN(1:2)) !transfN may change in retention_pool
        CALL production_pool(2,inorgNpool(1:2),transfN(1:2))
      ELSE
        CALL retention_pool(1,soilstate%fastN(1,j,i),transfN(1)) !transfN may change in retention_pool
        CALL production_pool(1,inorgNpool(1),transfN(1))
      ENDIF
    ENDIF

    !Calculate the new soil concentrations of IN
    IF(calcN)THEN
      CALL new_concentration(inorgNpool(1),soilstate%water(1,j,i),soilstate%conc(i_in,1,j,i))
      IF(thickness(2)>0) CALL new_concentration(inorgNpool(2),soilstate%water(2,j,i),soilstate%conc(i_in,2,j,i))
      IF(thickness(3)>0) CALL new_concentration(inorgNpool(3),soilstate%water(3,j,i),soilstate%conc(i_in,3,j,i))
    ENDIF

    !Degradation of humusP to fastP           
    degradhP = 0.
    IF(calcP)THEN
      degradhP(:) = degrhPpar * tmpfcn(:) * smfcn(:) * soilstate%humusP(:,j,i)
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%humusP(:,j,i),degradhP(:))      !may change degradhP
        CALL production_pool(maxsoillayers,soilstate%fastP(:,j,i),degradhP(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%humusP(1:2,j,i),degradhP(1:2))      !may change degradhP
        CALL production_pool(2,soilstate%fastP(1:2,j,i),degradhP(1:2))
      ELSE
        CALL retention_pool(1,soilstate%humusP(1,j,i),degradhP(1))      !may change degradhP
        CALL production_pool(1,soilstate%fastP(1,j,i),degradhP(1))
      ENDIF
    ENDIF

    !Transformation of fastP to SRP
    transfP = 0.
    IF(calcP)THEN
      transfP(:) = minfPpar * tmpfcn(:) * smfcn(:) * soilstate%fastP(:,j,i)
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%fastP(:,j,i),transfP(:)) !transfP may change in retention_pool
        CALL production_pool(maxsoillayers,SRPpool(:),transfP(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%fastP(1:2,j,i),transfP(1:2)) !transfP may change in retention_pool
        CALL production_pool(2,SRPpool(1:2),transfP(1:2))
      ELSE
        CALL retention_pool(1,soilstate%fastP(1,j,i),transfP(1)) !transfP may change in retention_pool
        CALL production_pool(1,SRPpool(1),transfP(1))
      ENDIF
    ENDIF

    !Calculate the new soil concentrations of SRP
    IF(calcP)THEN
      CALL new_concentration(SRPpool(1),soilstate%water(1,j,i),soilstate%conc(i_sp,1,j,i))
      IF(thickness(2)>0) CALL new_concentration(SRPpool(2),soilstate%water(2,j,i),soilstate%conc(i_sp,2,j,i))
      IF(thickness(3)>0) CALL new_concentration(SRPpool(3),soilstate%water(3,j,i),soilstate%conc(i_sp,3,j,i))
    ENDIF

    !Transformation of fastN to dissolved organic N
    dissolfN = 0.
    IF(calcN)THEN
      dissolfN(:) = pardisfN * tmpfcn(:) * smfcn(:) * soilstate%fastN(:,j,i)  
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%fastN(:,j,i),dissolfN(:)) !dissolfN may change in retention_pool
        CALL production_pool(maxsoillayers,DONpool(:),dissolfN(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%fastN(1:2,j,i),dissolfN(1:2)) !dissolfN may change in retention_pool
        CALL production_pool(2,DONpool(1:2),dissolfN(1:2))
      ELSE
        CALL retention_pool(1,soilstate%fastN(1,j,i),dissolfN(1)) !dissolfN may change in retention_pool
        CALL production_pool(1,DONpool(1),dissolfN(1))
      ENDIF
    ENDIF

    !Transformation of humusN to dissolved organic N
    dissolhN = 0.
    IF(calcN)THEN
      dissolhN(:) = pardishN * tmpfcn(:) * smfcn(:) * soilstate%humusN(:,j,i)
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%humusN(:,j,i),dissolhN(:)) !dissolhN may change in retention_pool
        CALL production_pool(maxsoillayers,DONpool(:),dissolhN(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%humusN(1:2,j,i),dissolhN(1:2)) !dissolhN may change in retention_pool
        CALL production_pool(2,DONpool(1:2),dissolhN(1:2))
      ELSE
        CALL retention_pool(1,soilstate%humusN(1,j,i),dissolhN(1)) !dissolhN may change in retention_pool
        CALL production_pool(1,DONpool(1),dissolhN(1))
      ENDIF
    ENDIF

    !Calculate the new soil concentrations of ON
    IF(calcN)THEN
      CALL new_concentration(DONpool(1),soilstate%water(1,j,i),soilstate%conc(i_on,1,j,i))
      IF(thickness(2)>0) CALL new_concentration(DONpool(2),soilstate%water(2,j,i),soilstate%conc(i_on,2,j,i))
      IF(thickness(3)>0) CALL new_concentration(DONpool(3),soilstate%water(3,j,i),soilstate%conc(i_on,3,j,i))
    ENDIF

    !Transformation of fastP to dissolved organic P (DOP/PP)
    dissolfP = 0.
    IF(calcP)THEN
      dissolfP(:) = pardisfP * tmpfcn(:) * smfcn(:) * soilstate%fastP(:,j,i)  
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%fastP(:,j,i),dissolfP(:)) !dissolfP may change in retention_pool
        CALL production_pool(maxsoillayers,DOPpool(:),dissolfP(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%fastP(1:2,j,i),dissolfP(1:2)) !dissolfP may change in retention_pool
        CALL production_pool(2,DOPpool(1:2),dissolfP(1:2))
      ELSE
        CALL retention_pool(1,soilstate%fastP(1,j,i),dissolfP(1)) !dissolfP may change in retention_pool
        CALL production_pool(1,DOPpool(1),dissolfP(1))
      ENDIF
    ENDIF

    !Transformation of humusP to dissolved organic P
    dissolhP = 0.
    IF(calcP)THEN
      dissolhP(:) = pardishP * tmpfcn(:) * smfcn(:) * soilstate%humusP(:,j,i)
      IF(thickness(3)>0)THEN
        CALL retention_pool(maxsoillayers,soilstate%humusP(:,j,i),dissolhP(:)) !dissolhP may change in retention_pool
        CALL production_pool(maxsoillayers,DOPpool(:),dissolhP(:))
      ELSEIF(thickness(2)>0)THEN
        CALL retention_pool(2,soilstate%humusP(1:2,j,i),dissolhP(1:2)) !dissolhP may change in retention_pool
        CALL production_pool(2,DOPpool(1:2),dissolhP(1:2))
      ELSE
        CALL retention_pool(1,soilstate%humusP(1,j,i),dissolhP(1)) !dissolhP may change in retention_pool
        CALL production_pool(1,DOPpool(1),dissolhP(1))
      ENDIF
    ENDIF

    !Calculate the new soil concentrations of PP (DOP)
    IF(calcP)THEN
      CALL new_concentration(DOPpool(1),soilstate%water(1,j,i),soilstate%conc(i_pp,1,j,i))
      IF(thickness(2)>0) CALL new_concentration(DOPpool(2),soilstate%water(2,j,i),soilstate%conc(i_pp,2,j,i))
      IF(thickness(3)>0) CALL new_concentration(DOPpool(3),soilstate%water(3,j,i),soilstate%conc(i_pp,3,j,i))
    ENDIF

    !Sum the sources (kg/km2)
    IF(calcN)THEN
       source(i_in)=transfN(1)+transfN(2)+transfN(3)
       source(i_on)=SUM(dissolfN)+SUM(dissolhN)
    ENDIF
    IF(calcP)THEN
       source(i_sp)=transfP(1)+transfP(2)+transfP(3)
       source(i_pp)=SUM(dissolfP)+SUM(dissolhP)
    ENDIF

  END SUBROUTINE soil_pool_transformations

  !>\brief Calculate retention in soil due to plant uptake
  !> Inorganic nitrogen and phosphorus is removed.
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Vegetation nutrient uptake)
  !-------------------------------------------------------------------
  SUBROUTINE plant_uptake(calcN,calcP,common_uptake,wp,thickness,soil,csoil,sink)

    USE MODVAR, ONLY : numsubstances, &
         maxsoillayers,       &
         i_in,i_sp
         
    !Argument declarations
    LOGICAL, INTENT(IN) :: calcN                    !<flag for nitrogen simulation
    LOGICAL, INTENT(IN) :: calcP                    !<flag for phosphorus simulation
    REAL, INTENT(IN)    :: common_uptake(2,2)       !<plant uptake (kg/km2/d), (N:P,sl1:sl2)
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<water content at wilting point (mm)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(IN)    :: soil(maxsoillayers)      !<soil moisture (mm)
    REAL, INTENT(INOUT) :: csoil(numsubstances,maxsoillayers)     !<concentration in soil moisture (mg/L etc)
    REAL, INTENT(OUT)   :: sink(numsubstances)      !<sink of nutrients in this subroutine (kg/km2)
    
    !Local variables
    REAL uptakeN(maxsoillayers),uptakeP(maxsoillayers)      !kg/km2
    REAL inorgNpool(maxsoillayers),SRPpool(maxsoillayers) !mg/L*mm=kg/km2
    REAL maxpooluptake(2)

    sink = 0.
    IF(numsubstances == 0) RETURN

    !>\b Algorithm \n
    !>Initial calculations; IN and SRP pools, maximum uptake
    inorgNpool = 0.
    SRPpool    = 0.
    IF(calcN)THEN
      inorgNpool(1) = csoil(i_in,1) * soil(1)
      IF(thickness(2)>0) inorgNpool(2) = csoil(i_in,2) * soil(2)
    ENDIF
    IF(calcP)THEN
      SRPpool(1) = csoil(i_sp,1) * soil(1)
      IF(thickness(2)>0) SRPpool(2) = csoil(i_sp,2) * soil(2)
    ENDIF
    maxpooluptake(1) = (soil(1)-wp(1))/soil(1)                      !Maximum uptake dependent on soil moisture
    IF(thickness(2)>0) maxpooluptake(2) = (soil(2)-wp(2))/soil(2)

    !>Calculate plant nutrient uptake
    uptakeN = 0.
    uptakeP = 0.
    IF(calcN)THEN
      uptakeN(1) = MIN(common_uptake(1,1), maxpooluptake(1)*inorgNpool(1))
      IF(thickness(2)>0)THEN
        uptakeN(2) = MIN(common_uptake(1,2), maxpooluptake(2)*inorgNpool(2))
      ENDIF
      IF(thickness(2)>0)THEN
        CALL retention_pool(2,inorgNpool(1:2),uptakeN(1:2)) !uptake may change in retention_pool
      ELSE
        CALL retention_pool(1,inorgNpool(1),uptakeN(1)) !uptake may change in retention_pool
      ENDIF
    ENDIF
    IF(calcP)THEN    
      uptakeP(1) = MIN(common_uptake(2,1), maxpooluptake(1)*SRPpool(1))
      IF(thickness(2)>0)THEN
        uptakeP(2) = MIN(common_uptake(2,2), maxpooluptake(2)*SRPpool(2))
      ENDIF
      IF(thickness(2)>0)THEN
        CALL retention_pool(2,SRPpool(1:2),uptakeP(1:2))
      ELSE
        CALL retention_pool(1,SRPpool(1),uptakeP(1))
      ENDIF
    ENDIF

    !>Sum up the sinks (kg/km2)
    IF(calcN) THEN
      IF(thickness(2)>0)THEN
        sink(i_in)=uptakeN(1)+uptakeN(2)
      ELSE
        sink(i_in)=uptakeN(1)
      ENDIF
    ENDIF
    IF(calcP)THEN
      IF(thickness(2)>0)THEN
        sink(i_sp)=uptakeP(1)+uptakeP(2)
      ELSE
        sink(i_sp)=uptakeP(1)
      ENDIF
    ENDIF

    !>Calculate the new soil concentrations
    IF(calcN)THEN
      CALL new_concentration(inorgNpool(1),soil(1),csoil(i_in,1))
      IF(thickness(2)>0) CALL new_concentration(inorgNpool(2),soil(2),csoil(i_in,2))
    ENDIF
    IF(calcP)THEN
      CALL new_concentration(SRPpool(1),soil(1),csoil(i_sp,1))
      IF(thickness(2)>0) CALL new_concentration(SRPpool(2),soil(2),csoil(i_sp,2))
    ENDIF

  END SUBROUTINE plant_uptake

  !>\brief Calculate retention in soil due to denitritfication
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Denitrification)
  !----------------------------------------------------------------------
  SUBROUTINE calculate_soil_profile_denitrification(calcN,i,j,iluse,wp,fc,ep,thickness,denpar_basic,soilstate,denitrification)
    
    USE HYPEVARIABLES, ONLY : m_denitr3, m_hsatINsoil, m_denit3reg
    USE MODVAR, ONLY : basin, &
                       numsubstances, &
                       maxsoillayers, &
                       regiondivision, &
                       landpar, &
                       regpar, &
                       genpar, &
                       i_in

    LOGICAL, INTENT(IN) :: calcN                    !<Status of nitrogen simulation
    INTEGER, INTENT(IN) :: i                        !<index of subbasin 
    INTEGER, INTENT(IN) :: j                        !<index of class
    INTEGER, INTENT(IN) :: iluse                    !<index of land use
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<water content at wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<water content at field capacity (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers)        !<water content in effective porosity (mm)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(IN)    :: denpar_basic             !<generally adjusted denitrification rate
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate !<Soil states
    REAL, INTENT(OUT)   :: denitrification(maxsoillayers)  !<denitrification
   
    !Local variables 
    REAL denpar, porevolume(maxsoillayers)
    REAL sink(numsubstances)
    
    IF(calcN)THEN
      porevolume = wp+fc+ep
      denitrification = 0.
      denpar = denpar_basic
      IF(thickness(1)>0.) CALL soil_denitrification(porevolume(1),denpar,genpar(m_hsatINsoil),soilstate%water(1,j,i),soilstate%temp(1,j,i),soilstate%conc(:,1,j,i),sink)
      denitrification(1) = sink(i_in)
      IF(thickness(2)>0.) CALL soil_denitrification(porevolume(2),denpar,genpar(m_hsatINsoil),soilstate%water(2,j,i),soilstate%temp(2,j,i),soilstate%conc(:,2,j,i),sink)
      denitrification(2) = sink(i_in)
      IF(thickness(3)>0.)THEN
        IF(landpar(m_denitr3,iluse)>0.) denpar=landpar(m_denitr3,iluse)
        IF(landpar(m_denitr3,iluse)<0.) denpar=0.
        IF(regpar(m_denit3reg,basin(i)%parregion(regiondivision(m_denit3reg)))>0.) denpar=regpar(m_denit3reg,basin(i)%parregion(regiondivision(m_denit3reg)))
        CALL soil_denitrification(porevolume(3),denpar,genpar(m_hsatINsoil),soilstate%water(3,j,i),soilstate%temp(3,j,i),soilstate%conc(:,3,j,i),sink)
        denitrification(3) = sink(i_in)
      ENDIF
    ELSE
      denitrification = 0.
    ENDIF
    
  END SUBROUTINE calculate_soil_profile_denitrification

  !>\brief Calculate retention in soil due to denitritfication
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Denitrification)
  !----------------------------------------------------------------------
  SUBROUTINE soil_denitrification(maxwc,denpar,halfsatIN,soil,stemp,csoil,sink)

    USE HYPEVARIABLES, ONLY : smfdenitlim,    &
                              smfdenitpow
    USE MODVAR, ONLY : numsubstances,       &
                       i_in

    !Argument declarations
    REAL, INTENT(IN)    :: maxwc        !<Maximum water content of soil (mm)
    REAL, INTENT(IN)    :: denpar       !<model parameter denitrification in soil
    REAL, INTENT(IN)    :: halfsatIN    !<model parameter half saturation IN of soil (mg/L)
    REAL, INTENT(IN)    :: soil         !<soil water (mm)
    REAL, INTENT(IN)    :: stemp        !<soil temperature (degree Celcius)
    REAL, INTENT(INOUT) :: csoil(numsubstances)     !<concentration of soil water
    REAL, INTENT(OUT)   :: sink(numsubstances)      !<sink of nutrient in this subroutine (kg/km2)

    !Local variables
    REAL denitr(1)
    REAL inorgNpool(1)                   !mg/L*mm=kg/km2
    REAL tmpfcn, smfcn, concfcn

    !Substance modelled
    sink = 0.
    IF(i_in==0) RETURN
    inorgNpool(1) = csoil(i_in) * soil    !Initial IN pool

    !>Dependence factors of denitrification
    tmpfcn = tempfactor(stemp)
    smfcn = exponential_moisturefactor(soil,maxwc,smfdenitlim,smfdenitpow)
    concfcn = halfsatconcfactor(csoil(i_in),halfsatIN)

    !>Denitrification
    denitr(1) = denpar * inorgNpool(1) * tmpfcn * smfcn * concfcn
    CALL retention_pool(1,inorgNpool(1),denitr(1)) !denitr may change in retention_pool

    !>Set the sink (kg/km2)
    sink(i_in) = denitr(1)

    !>Calculate the new soil concentrations
    CALL new_concentration(inorgNpool(1),soil,csoil(i_in))

  END SUBROUTINE soil_denitrification

  !>\brief Exchange of SRP in soil water with soil PartP pool,
  !>adsorption-desorption
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Balance SP-PartP)
  !-----------------------------------------------------------------------
  SUBROUTINE balance_spsoil(i,j,calcP,thickness,kfr,nfr,kadsdes,soilstate)

    USE MODVAR, ONLY : i_sp,      &
         maxsoillayers

    !Argument declaration
    INTEGER, INTENT(IN)  :: i                        !<current subbasin index
    INTEGER, INTENT(IN)  :: j                        !<current soiltype-landuse class index
    LOGICAL, INTENT(IN)  :: calcP                    !<flag for nitrogen simulation
    REAL, INTENT(IN)     :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(IN)     :: kfr                      !<model parameter freundlich
    REAL, INTENT(IN)     :: nfr                      !<model parameter freundlich
    REAL, INTENT(IN)     :: kadsdes                  !<model parameter freundlich, adsorption or desportion rate (1/d) 
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate  !<Soil states
    
    !Calculate exchange between SRP and PartP
    IF(calcP)THEN
      CALL freundlich(i,j,1,soilstate%partP(1,j,i), soilstate%conc(i_sp,1,j,i), soilstate%water(1,j,i), thickness(1), Kfr, Nfr, Kadsdes)  
      IF(thickness(2)>0)  CALL freundlich(i,j,2,soilstate%partP(2,j,i), soilstate%conc(i_sp,2,j,i), soilstate%water(2,j,i), thickness(2), Kfr, Nfr, Kadsdes)     
      IF(thickness(3)>0)  CALL freundlich(i,j,3,soilstate%partP(3,j,i), soilstate%conc(i_sp,3,j,i), soilstate%water(3,j,i), thickness(3), Kfr, Nfr, Kadsdes)
    ENDIF

  END SUBROUTINE balance_spsoil

  !>Simulate crop rotation by averaging soil nutrient pools once a year
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Crop rotation effect on soil pools)
  !-------------------------------------------------------------------
  SUBROUTINE croprotation_soilpoolaverage(i,ncrot,pool)

    USE MODVAR, ONLY : classbasin,    &
                       classdata,     &
                       soilthick,     &
                       nsub,          &
                       nclass,        &
                       maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<subbasin index
    INTEGER, INTENT(IN) :: ncrot      !<number of crop rotation sequences
    REAL, INTENT(INOUT) :: pool(maxsoillayers,nclass,nsub)     !<soil solid nutrient pool
    
    !Local variables
    INTEGER j,k   !loop-variables
    REAL sumpool(maxsoillayers)
    REAL sumarea(maxsoillayers), a

    DO k=1,ncrot
      sumpool = 0.
      sumarea = 0.
      DO j = 1,nclass
        a = classbasin(i,j)%part
        IF(k==classdata(j)%rotation .AND. a>0.)THEN                  
          sumpool(1) = sumpool(1) + pool(1,j,i)*a
          sumarea(1) = sumarea(1) + a
          IF(soilthick(2,j)>0)THEN
            sumpool(2) = sumpool(2) + pool(2,j,i)*a
            sumarea(2) = sumarea(2) + a
          ENDIF
          IF(soilthick(3,j)>0)THEN
            sumpool(3) = sumpool(3) + pool(3,j,i)*a
            sumarea(3) = sumarea(3) + a
          ENDIF
        ENDIF
      ENDDO !j 
      IF(sumarea(1)>0.)THEN
        sumpool(1) = sumpool(1) / sumarea(1)
      ENDIF
      IF(sumarea(2)>0.)THEN
        sumpool(2) = sumpool(2) / sumarea(2)
      ENDIF
      IF(sumarea(3)>0.)THEN
        sumpool(3) = sumpool(3) / sumarea(3)
      ENDIF
      DO j = 1,nclass
        a = classbasin(i,j)%part
        IF(k==classdata(j)%rotation .AND. a>0.)THEN
          pool(1,j,i) = sumpool(1)
          IF(soilthick(2,j)>0 .AND. sumarea(2)>0.)THEN
            pool(2,j,i) = sumpool(2)
          ENDIF
          IF(soilthick(3,j)>0 .AND. sumarea(3)>0.)THEN
            pool(3,j,i) = sumpool(3)
          ENDIF
        ENDIF
      ENDDO !j
    ENDDO !k       

  END SUBROUTINE croprotation_soilpoolaverage

  !>\brief Calculate the balance between P in solution and mineral P
  !@todo Reference to literature for freundlich is needed.
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Balance SP-PartP)
  !--------------------------------------------------------------------------
  SUBROUTINE  freundlich(iin,jin,lin,PoolPP,SRP_Conc,Vol,LayerThick,Kfr,Nfr,Kadsdes)

    USE HYPEVARIABLES, ONLY : bulkdensity

    !Argument declarations
    INTEGER, INTENT(IN) :: iin         !<index of subbasin
    INTEGER, INTENT(IN) :: jin         !<index of class
    INTEGER, INTENT(IN) :: lin         !<index of soillayer
    REAL, INTENT(INOUT) :: poolPP      !<Pool PP  (kg/km2)
    REAL, INTENT(INOUT) :: SRP_conc    !<SRP conc of soil water (mg/l)
    REAL, INTENT(IN)    :: vol         !<Soil water content in layer (mm)
    REAL, INTENT(IN)    :: layerthick  !<soil layer thickness (m)
    REAL, INTENT(IN)    :: kfr         !<freundlich adsorption isoterm (l/kg)
    REAL, INTENT(IN)    :: nfr         !<freundlich exponential coefficient
    REAL, INTENT(IN)    :: kadsdes     !<adsorption/desorption coefficient (1/d)
      
    !Local parameters
    REAL, PARAMETER :: Limit = 0.00001 !Threshold for breaking iterations
    !Local variables
    INTEGER i
    REAL totalP, PPequi_conc
    REAL conc_sol, adsdes
    REAL x0, xn, xn_1, fxn, fprimxn, dx
    REAL coeff
    REAL nfrloc
    REAL help

    IF(Vol==0) RETURN
    IF(Kfr==0 .OR. Nfr==0 .OR. Kadsdes==0)THEN
      WRITE(6,*) 'Error: Values for freundlich parameters missing'
      STOP 1
    ENDIF

    nfrloc = Nfr
    totalP = PoolPP + SRP_conc * vol  !Total amount of P (kg/km2) 
    IF(totalP==0.) RETURN
    conc_sol = (poolPP / bulkdensity) / layerthick        !mg P /kg soil
    IF(conc_sol<=0.)THEN
      nfrloc = 1.
      write(6,*) 'Warning: soil partP <=0. Freundlich will give error, take shortcut.'
    ENDIF
    coeff = kfr * bulkdensity * layerthick

    IF(nfrloc==1.)THEN
      xn_1 = totalP / (vol + coeff)
      PPequi_conc = kfr * xn_1  
    ELSE  
      !Newton-Raphson method used to calculate equilibrium concentration
      x0 = exp((log(conc_sol)-log(kfr))/ nfr)  !Initial guess of equilibrium liquid conc
      fxn = x0 * vol + coeff * (x0 ** nfr) - totalP 
      xn = x0
      xn_1 = xn
      i = 0
      DO WHILE (ABS(fxn) > limit .AND. i < 20) !iteration to calculate equilibrium concentations
        fxn = xn * vol + coeff * (xn ** nfr) - totalP
        fprimxn = vol + nfr * coeff * ( xn **(nfr-1)) 
        dx = fxn / fprimxn
        IF(ABS(dx)<0.000001*xn) EXIT
        xn_1 = xn - dx
        IF(xn_1 <= 0.) THEN 
          xn_1 = 1.0E-10
        ENDIF
        xn = xn_1
        i = i+1
      ENDDO
      PPequi_conc = kfr * (xn_1 ** nfr)    
    ENDIF

    !Calculate new pool and concentration, depends on the equilibrium concentration
    IF(ABS(PPequi_conc - conc_sol) > 10E-6)THEN
      adsdes = (PPequi_conc - conc_sol)*(1-EXP(-kadsdes)) !kinetic adsorption/desorption
      help = adsdes * bulkdensity * layerthick
      IF(-help>poolPP.OR.SRP_conc<(help/ vol))THEN
        IF(-help>poolPP) help = -poolPP
        IF(SRP_conc<(help/ vol)) help = SRP_conc*vol
        WRITE(6,*) 'Warning: freundlich flow adjusted, was larger than pool'
      ENDIF
      poolPP = poolPP + help   !new Pool PP
      SRP_conc = SRP_conc - (help/ vol)  !new liquid conc
    ENDIF

    !Liten säkerhetsåtgärd ifall negativa SRP_conc skulle dyka upp igen. Så det märks.
    !Ta bort om ett år (idag 140506) ifall inte använd innan dess.
    IF(SRP_conc<0.) THEN
      WRITE(6,*) 'ERROR: SRP_Conc in freundlich negative!',SRP_conc
      WRITE(6,*) iin,jin,lin
      STOP 1
    ENDIF

  END SUBROUTINE freundlich

  !>\brief Add load from local diffuse sources to the lowest soil
  !>layer
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land 
  !!routines (Nutrient sources - Rural household diffuse source)
  !-----------------------------------------------------------------
  SUBROUTINE local_diffuse_source(i,j,pw,classarea,soilstate,ruralaload,upwardflow,ruralflow,cruralflow)

    USE MODVAR, ONLY : genpar,  &
                       load,  &
                       numsubstances,   &
                       maxsoillayers,   &
                       i_in,i_on,i_sp,i_pp,   &
                       soilthick, &
                       landarea
    USE HYPEVARIABLES, ONLY : m_locsoil

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<current subbasin index
    INTEGER, INTENT(IN) :: j        !<current slc-class index
    REAL, INTENT(IN)    :: pw(maxsoillayers)  !<pore volume of soil (mm)
    REAL, INTENT(IN)    :: classarea  !<area of class (km2)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    REAL, INTENT(OUT) :: ruralaload(numsubstances)    !<Load from rural households (kg/timestep)
    REAL, INTENT(OUT) :: upwardflow(2)   !<upwelling due to overflowing lower soil layers (mm/timestep)
    REAL, INTENT(OUT) :: ruralflow(maxsoillayers)    !<added flow from rural (mm/timestep)
    REAL, INTENT(OUT) :: cruralflow(numsubstances)    !<concentration from rural households (mg/L)
    
    !Local variables
    INTEGER soillayer !soil layer rural load was added to
    REAL cin(numsubstances),qin    !concentration and water flow for rural source to soil
    REAL ppload(1)  !PP-load of rural soice to soil

    !Default output value
    upwardflow = 0.
    ruralflow = 0.
    ruralaload = 0.
    cruralflow = 0.
    
    !Set local parameter
    IF(genpar(m_locsoil)>0)THEN
      qin = genpar(m_locsoil) * load(i)%volloc / landarea(i) * 1000.     !->mm/d !fel vid korta tidsteg!
      IF(qin>0)THEN

        !Set concentration of rural inflow to lowest soillayer (IN,ON, and SP)
        cin = 0.  !for other substances, maybe T2 should be something else     
        IF(i_in>0)THEN 
          cin(i_in) = load(i)%locconc(i_in)
          cin(i_on) = load(i)%locconc(i_on)
          ruralaload(i_in) = cin(i_in) * qin
          ruralaload(i_on) = cin(i_on) * qin
        ENDIF
        IF(i_sp>0)THEN 
          cin(i_sp) = load(i)%locconc(i_sp)   !PP-conc = 0
          ruralaload(i_sp) = cin(i_sp) * qin
        ENDIF
        cruralflow = load(i)%locconc

        !Add rural inflow to lowest soillayer
        CALL inflow_lowest_soillayer(numsubstances,maxsoillayers,cin,qin,pw,  &
                   soilthick(:,j),soilstate%water(:,j,i),soilstate%conc(:,:,j,i),upwardflow,ruralflow,soillayer)

        !Rural inflow of PP added to fastP of lowest soillayer            
        IF(i_sp>0)THEN 
          ppload(1) = load(i)%locconc(i_pp) * qin       !kg/km2=mg/L*mm
          CALL production_pool(1,soilstate%fastP(soillayer,j,i),ppload)
          ruralaload(i_pp) = ppload(1)   !kg/km2/timestep
        ENDIF

        ruralaload(:) = ruralaload(:) * classarea         !Transform diffuse load, ruralA to kg/timestep
      ENDIF
    ENDIF

  END SUBROUTINE local_diffuse_source

  !>\brief Subroutine for organic carbon transformation processes in
  !>the soil
  !!
  !> \b Reference ModelDescription Chapter Organic carbon (Source of organic material and Soil processes)
  !-----------------------------------------------------------------
  SUBROUTINE soil_carbon_processes(i,j,wp,fc,ep,pw,thickness,klh,klo,&
       kho,kof,koflim,minc,soimf,soimr,soilstate)       

    USE MODVAR, ONLY : maxsoillayers

    INTEGER, INTENT(IN) :: i                        !<current subbasin index
    INTEGER, INTENT(IN) :: j                        !<current slc-class index
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<water content at wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<water content at field capacity (mm)
    REAL, INTENT(IN)    :: ep(maxsoillayers)        !<water content: effectiv porosity (mm)
    REAL, INTENT(IN)    :: pw(maxsoillayers)        !<water content: total porosity (mm)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(IN)    :: klh      !<transformation parameter
    REAL, INTENT(IN)    :: klo      !<rate of transformation
    REAL, INTENT(IN)    :: kho      !<rate of transformation
    REAL, INTENT(IN)    :: kof      !<rate of transformation
    REAL, INTENT(IN)    :: koflim   !<threshold for transformation
    REAL, INTENT(IN)    :: minc     !<fraction mineralisation to DIC (-)
    REAL, INTENT(IN)    :: soimf    !<saturation soilmoisture factor (-)
    REAL, INTENT(IN)    :: soimr    !<rate soilmoisture factor (-)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states

    !Calculate the nutrient processes
    CALL soil_carbon_pool_transformations(i,j,wp,fc,pw,thickness,klh,klo,kho,kof,koflim,minc,soimf,soimr,soilstate)

  END SUBROUTINE soil_carbon_processes
  
  !>Subroutine for organic carbon transformation between pools in the
  !>soil
  !!
  !> \b Reference ModelDescription Chapter Organic carbon (Soil processes)
  !--------------------------------------------------------------------
  SUBROUTINE soil_carbon_pool_transformations(i,j,wp,fc,pw, &
       thickness,klh,klo,kho,kof,koflim,minc,soimf,soimr,soilstate)  
 
    USE MODVAR, ONLY : maxsoillayers,i_oc
    USE HYPEVARIABLES, ONLY :   &
         thetapow,thetalow,thetaupp  

    !Argument declaration
    INTEGER, INTENT(IN) :: i      !<current subbasin index
    INTEGER, INTENT(IN) :: j      !<current class index
    REAL, INTENT(IN)    :: wp(maxsoillayers)        !<water content at wilting point (mm)
    REAL, INTENT(IN)    :: fc(maxsoillayers)        !<water content at field capacity (mm)
    REAL, INTENT(IN)    :: pw(maxsoillayers)        !<water content: total porosity (mm)
    REAL, INTENT(IN)    :: thickness(maxsoillayers) !<thickness of soil layers
    REAL, INTENT(IN)    :: klh      !<transformation rate of fastC to humusC (d-1)
    REAL, INTENT(IN)    :: klo      !<degradation fastC (d-1)
    REAL, INTENT(IN)    :: kho      !<degradation humusC (d-1)
    REAL, INTENT(IN)    :: kof      !<transformation to fastC (d-1)
    REAL, INTENT(IN)    :: koflim   !<threshold for transformation to fastC (-)
    REAL, INTENT(IN)    :: minc     !<fraction mineralisation to DIC (-)      
    REAL, INTENT(IN)    :: soimf    !<satuaration soilmoisture factor (-)
    REAL, INTENT(IN)    :: soimr    !<rate soilmoisture factor (-)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states

    !Local variables
    INTEGER k   !soillayer
    REAL DOCpool(maxsoillayers)
    REAL fasttohumus(maxsoillayers)
    REAL doctofast(maxsoillayers)
    REAL transhC(maxsoillayers)
    REAL transfC(maxsoillayers)
    REAL tmpfcn(maxsoillayers),smfcn(maxsoillayers)
    REAL real1temp(1)  !helpvariable, needed for gfortran
    REAL fracprod     

    !Initialisation
    DOCpool = soilstate%water(:,j,i) * soilstate%conc(i_oc,:,j,i)      !Current DOC pool
    fracprod = 1. - minc    !fraction of degradationn that is not mineralised

    !Temperature dependence factor
    DO k = 1,maxsoillayers
       IF(thickness(k)>0) tmpfcn(k) = tempfactor(soilstate%temp(k,j,i))
    ENDDO

    !Soil moisture dependence factor
    smfcn(1) = moisturefactor(soilstate%water(1,j,i),wp(1),pw(1),thickness(1),soimf,thetapow,soimr,thetaupp)
    smfcn(2) = moisturefactor(soilstate%water(2,j,i),wp(2),pw(2),thickness(2),soimf,thetapow,soimr,thetaupp)
    smfcn(3) = moisturefactor(soilstate%water(3,j,i),wp(3),pw(3),thickness(3),soimf,thetapow,soimr,thetaupp)

    !Transformation between fastC and humusC
    fasttohumus = 0.
    IF(thickness(3)>0)THEN
      fasttohumus(:) = klh * tmpfcn(:) * smfcn(:) * soilstate%fastC(:,j,i)
      CALL retention_pool(3,soilstate%fastC(:,j,i),fasttohumus(:)) !fasttohumus may change in retention_pool
      CALL production_pool(3,soilstate%humusC(:,j,i),fracprod*fasttohumus(:))
    ELSEIF(thickness(2)>0)THEN
      fasttohumus(1:2) = klh * tmpfcn(1:2) * smfcn(1:2) * soilstate%fastC(1:2,j,i)
      CALL retention_pool(2,soilstate%fastC(1:2,j,i),fasttohumus(1:2)) !fasttohumus may change in retention_pool
      CALL production_pool(2,soilstate%humusC(1:2,j,i),fracprod*fasttohumus(1:2))
    ELSE
      fasttohumus(1) = klh * tmpfcn(1) * smfcn(1) * soilstate%fastC(1,j,i)
      CALL retention_pool(1,soilstate%fastC(1,j,i),fasttohumus(1)) !fasttohumus may change in retention_pool
      real1temp(1)=fracprod*fasttohumus(1)
      CALL production_pool(1,soilstate%humusC(1,j,i),real1temp(1))
    ENDIF

    !Transformation of DOC to fastC
    doctofast = 0.
    DO k = 1,maxsoillayers
      IF(thickness(k)>0.AND.smfcn(k)<koflim.AND.soilstate%water(k,j,i)<fc(k)+wp(k).AND.soilstate%temp(k,j,i)<5.)THEN
        doctofast(k)=kof*DOCpool(k)
      ENDIF
    ENDDO
    IF(thickness(3)>0)THEN
      CALL retention_pool(3,DOCpool(:),doctofast(:)) !doctofast may change in retention_pool
      CALL production_pool(3,soilstate%fastC(:,j,i),fracprod*doctofast(:))
    ELSEIF(thickness(2)>0)THEN
      CALL retention_pool(2,DOCpool(1:2),doctofast(1:2)) !doctofast may change in retention_pool
      CALL production_pool(2,soilstate%fastC(1:2,j,i),fracprod*doctofast(1:2))
    ELSE
      CALL retention_pool(1,DOCpool(1),doctofast(1)) !doctofast may change in retention_pool
      real1temp(1)=fracprod*doctofast(1)
      CALL production_pool(1,soilstate%fastC(1,j,i),real1temp(1))
    ENDIF

    !Transformation of fastC to DOC
    transfC = 0.
    IF(thickness(3)>0)THEN
      transfC(:) = klo * tmpfcn(:) * smfcn(:) * soilstate%fastC(:,j,i)
      CALL retention_pool(maxsoillayers,soilstate%fastC(:,j,i),transfC(:)) !transfC may change in retention_pool
      CALL production_pool(maxsoillayers,DOCpool(:),fracprod*transfC(:))
    ELSEIF(thickness(2)>0)THEN
      transfC(1:2) = klo * tmpfcn(1:2) * smfcn(1:2) * soilstate%fastC(1:2,j,i)
      CALL retention_pool(2,soilstate%fastC(1:2,j,i),transfC(1:2)) !transfC may change in retention_pool
      CALL production_pool(2,DOCpool(1:2),fracprod*transfC(1:2))
    ELSE
      transfC(1) = klo * tmpfcn(1) * smfcn(1) * soilstate%fastC(1,j,i)
      CALL retention_pool(1,soilstate%fastC(1,j,i),transfC(1)) !transfC may change in retention_pool
      real1temp(1)=fracprod*transfC(1)
      CALL production_pool(1,DOCpool(1),real1temp(1))
    ENDIF

    !Transformation of humusC to DOC
    transhC = 0.
    IF(thickness(3)>0)THEN
      transhC(:) = kho * tmpfcn(:) * smfcn(:) * soilstate%humusC(:,j,i)
      CALL retention_pool(maxsoillayers,soilstate%humusC(:,j,i),transhC(:)) !transhC may change in retention_pool
      CALL production_pool(maxsoillayers,DOCpool(:),fracprod*transhC(:))
    ELSEIF(thickness(2)>0)THEN
      transhC(1:2) = kho * tmpfcn(1:2) * smfcn(1:2) * soilstate%humusC(1:2,j,i)
      CALL retention_pool(2,soilstate%humusC(1:2,j,i),transhC(1:2)) !transhC may change in retention_pool
      CALL production_pool(2,DOCpool(1:2),fracprod*transhC(1:2))
    ELSE
      transhC(1) = kho * tmpfcn(1) * smfcn(1) * soilstate%humusC(1,j,i)
      CALL retention_pool(1,soilstate%humusC(1,j,i),transhC(1))     !transhC may change in retention_pool
      real1temp(1)=fracprod*transhC(1)
      CALL production_pool(1,DOCpool(1),real1temp(1))
    ENDIF

    !Calculate the new soil concentrations of DOC
    CALL new_concentration(DOCpool(1),soilstate%water(1,j,i),soilstate%conc(i_oc,1,j,i))
    IF(thickness(2)>0) CALL new_concentration(DOCpool(2),soilstate%water(2,j,i),soilstate%conc(i_oc,2,j,i))
    IF(thickness(3)>0) CALL new_concentration(DOCpool(3),soilstate%water(3,j,i),soilstate%conc(i_oc,3,j,i))

  END SUBROUTINE soil_carbon_pool_transformations

  !------------------------------------------------------------------
  !>Subroutine for dissolved organic carbon reduction during
  !percolation down in the soil layers
  !!
  !> \b Reference ModelDescription Chapter Organic carbon (Soil processes - Percolation)
  !------------------------------------------------------------------
  SUBROUTINE doc_percolation_reduction(n,conc,kpar,temp,soilm,wiltp,porew,thickness)

    USE MODVAR, ONLY : i_oc
    USE HYPEVARIABLES, ONLY : satact,thetapow,thetalow,thetaupp  

    !Argument declarations
    INTEGER, INTENT(IN) :: n         !<number of substances, array dimension
    REAL, INTENT(INOUT) :: conc(n)   !<concentration of flow (mg/l)
    REAL, INTENT(IN)    :: kpar      !<reduction parameter
    REAL, INTENT(IN)    :: temp      !<soil layer temperature (degree Celcius)
    REAL, INTENT(IN)    :: soilm     !<soil water in soillayer (mm)
    REAL, INTENT(IN)    :: wiltp     !<wilting point of soillayer (mm)
    REAL, INTENT(IN)    :: porew     !<pore wolume of soillayer (mm)
    REAL, INTENT(IN)    :: thickness !<thickness of soillayer (m)
    
    !Local variables
    REAL newconc
    REAL tmpfcn
    REAL smfcn

    !No carbon modelled
    IF (i_oc == 0) RETURN

    !Temperature dependence factor
    tmpfcn = tempfactor(temp)

    !Soil moisture dependence factor
    smfcn = moisturefactor(soilm,wiltp,porew,thickness,satact,thetapow,thetalow,thetaupp)

    !Calculate the concentration reduction
    IF(conc(i_oc)>0)THEN
       newconc = conc(i_oc)*(1. - kpar*tmpfcn*smfcn)
    ELSE
       newconc = conc(i_oc)
    ENDIF
    conc(i_oc) = newconc

  END SUBROUTINE doc_percolation_reduction

  !>Subroutine for dissolved organic nitrogen and particulate
  !phosphorus reduction during percolation down in the soil layers
  !!
  !> \b Reference ModelDescription Chapter Nitrogen and phosphorus in land routines (Soil processes - Percolation)
  !------------------------------------------------------------
  SUBROUTINE onpp_percolation_reduction(n,conc,kparN,kparP)         

    USE MODVAR, ONLY : i_on, i_pp

    !Argument declarations
    INTEGER, INTENT(IN) :: n        !<number of substances, array dimension
    REAL, INTENT(INOUT) :: conc(n)  !<concentration of percolation water (mg/l)
    REAL, INTENT(IN)    :: kparN    !<reduction parameter on
    REAL, INTENT(IN)    :: kparP    !<reduction parameter pp

    !No nitrogen and phosphorus modelled
    IF (i_pp == 0.AND.i_on == 0) RETURN

    !Calculate the concentration reduction
    IF(i_on > 0)THEN
       conc(i_on) = conc(i_on)*(1.-kparN)
    ENDIF

    IF(i_pp > 0)THEN
       conc(i_pp) = conc(i_pp)*(1.-kparP)
    ENDIF

  END SUBROUTINE onpp_percolation_reduction

  !>Subroutine for dissolved organic carbon change in runoff from soil
  !!due to riparian zone processes
  !!
  !> \b Reference ModelDescription Chapter Organic carbon (Riparian zone)
  !-----------------------------------------------------------------------
  SUBROUTINE class_riparian_zone_processes(n,elev,q,conc1,conc2,conc3,oldgrw, &
       kpar,stemp,apar,spar,gwat,t10,t20,soilm,wiltp,porew,depthm)

    USE MODVAR, ONLY : i_oc,    &
         maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: n         !<number of substances, array dimension
    REAL, INTENT(IN)    :: elev      !<class elevation (möh)
    REAL, INTENT(IN)    :: q         !<runoff through riparian zone (mm)
    REAL, INTENT(INOUT) :: conc1(n)  !<concentration of runoff 1 (mg/l*mm)
    REAL, INTENT(INOUT) :: conc2(n)  !<concentration of runoff 2 (mg/l*mm)
    REAL, INTENT(INOUT) :: conc3(n)  !<concentration of runoff 3 (mg/l*mm)
    REAL, INTENT(INOUT) :: oldgrw    !<ground water table last time step (m)
    REAL, INTENT(IN)    :: kpar      !<model parameter; amount of riparian zone influence
    REAL, INTENT(IN)    :: stemp(maxsoillayers)  !<soil temperature of class (same all soil layers)
    REAL, INTENT(IN)    :: apar      !<model parameter
    REAL, INTENT(IN)    :: spar      !<model parameter
    REAL, INTENT(IN)    :: gwat      !<groundwater level of class (m)
    REAL, INTENT(IN)    :: t10       !<10-day mean air temperature
    REAL, INTENT(IN)    :: t20       !<20-day mean air temperature
    REAL, INTENT(IN)    :: soilm     !<soil moisture all soil layers (mm)
    REAL, INTENT(IN)    :: wiltp     !<wilting point of class (mm)
    REAL, INTENT(IN)    :: porew     !<pore wolume of class (all soil layers) (mm)
    REAL, INTENT(IN)    :: depthm    !<thickness of soil layer (all soillayers) (m)
    
    !Local variables
    INTEGER k  !soillayer
    LOGICAL grwrising
    REAL grwfcn, seafcn, smfcn
    REAL tmpfcn(maxsoillayers)
    REAL help,factor(maxsoillayers)

    !Riparian zone not active
    IF(kpar == 0) RETURN         !No riparian zone processes

    !Changes in groundwater table
    IF(gwat>=oldgrw)THEN   !Obs: normally gwat negative
      grwrising = .TRUE.
    ELSE
      grwrising = .FALSE.
    ENDIF
    oldgrw = gwat
    grwfcn = exp(apar * gwat)

    IF(q == 0) RETURN            !No soil layer runoff    

    !Seasonal dependence
    seafcn = 1.
    IF(t10<t20) seafcn = spar   !autumn

    !Temperature dependence factor based on soil temperature
    DO k = 1, maxsoillayers
      tmpfcn(k) = tempfactor(stemp(k))
    ENDDO

    !Soil moisture dependence factor 
    CALL riparian_moisturefactor(smfcn,soilm,wiltp,porew,depthm,grwrising)

    !Calculate increase in DOC (infinite source)
    help = kpar * (elev/1000.) * grwfcn * seafcn * smfcn
    factor = 1. + tmpfcn * help
    IF(conc1(i_oc)>0.) conc1(i_oc) = factor(1) * conc1(i_oc)
    IF(conc2(i_oc)>0.) conc2(i_oc) = factor(2) * conc2(i_oc)
    IF(conc3(i_oc)>0.) conc3(i_oc) = factor(3) * conc3(i_oc)

  END SUBROUTINE class_riparian_zone_processes

  !>Calculates the soil moisture dependence factor, different for
  !!rising and sinking ground water level
  !!
  !> \b Reference ModelDescription Chapter Organic carbon (Riparian zone)
  !-------------------------------------------------------------------
  SUBROUTINE riparian_moisturefactor(smfcn,sm,wp,pw,depthm,rising)
    USE HYPEVARIABLES, ONLY : satact,       &
         thetaupp,     &
         thetalow

    !Argument declarations
    REAL, INTENT(OUT)   :: smfcn  !<soil moisture dependence factor
    REAL, INTENT(IN)    :: sm     !<soil moisture
    REAL, INTENT(IN)    :: wp     !<wilting point pore wolume
    REAL, INTENT(IN)    :: pw     !<total pore wolume
    REAL, INTENT(IN)    :: depthm !<thickness of soil layers (m)
    LOGICAL, INTENT(IN) :: rising !<ground water table is rising?
    
    !Local variables
    REAL thickness

    !Start of subroutine
    smfcn = 0.
    thickness = depthm * 1000.      !mm
    IF(sm >= pw) THEN   
       smfcn = satact
    ELSEIF(sm <= wp) THEN
       smfcn=0.
    ELSE 
       IF(rising)THEN
          smfcn = MIN(1., (1-satact)*(pw-sm)/(thetaupp/100.*thickness) &
               + satact, (sm-wp)/(thetalow/100.*thickness))    
       ELSE  !sinking
          smfcn = MIN(1., (1-satact)*(pw-sm)/(thetaupp/100.*thickness) &
               + satact, (sm-wp)/(thetalow/100.*thickness)*satact)
       ENDIF
    ENDIF

  END SUBROUTINE riparian_moisturefactor

END MODULE NPC_SOIL_PROCESSES
