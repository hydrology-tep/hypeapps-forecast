!> \file irrigation.f90
!> Contains module irrigation_module.

!>Irrigation calculations in HYPE
MODULE IRRIGATION_MODULE

  !Copyright 2011-2015,2017 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

  !Used modules
  USE STATETYPE_MODULE
  !hypevariables
  !modvar
  IMPLICIT NONE
  PRIVATE
  !----------------------------------
  ! Private procedures
  !----------------------------------
  ! check_for_irrigated_classes
  ! set_regional_irrigation 
  ! set_regional_collectors 
  ! calculate_kcb 
  ! set_irrigation_season_end
  ! irrigation_season 
  ! immersion_season 
  ! irrigation_abstraction_sink
  !----------------------------------
  PUBLIC :: initiate_irrigation, &
            initiate_timestep_irrigation, &
            calculate_irrigation_water_demand, &
            apply_irrigation, &
            calculate_irrigation, &
            check_for_irrigated_classes

  !Private variable declarations
  SAVE
  INTEGER, ALLOCATABLE :: regirrwat_source(:)      !<Index i of subbasin which is the source of regional water for irrigation (nsub)
  INTEGER, ALLOCATABLE :: regirrwat_collectors(:,:)!<Index i of subbasins which has this subbasin as a source of regional water for irrigation (max number of collectors,nsub)
  INTEGER, ALLOCATABLE :: nregirrwat_collectors(:) !<Number of subbasins which has this subbasin as a source of regional water for irrigation (nsub)
  INTEGER, ALLOCATABLE :: irrindex(:)              !<Index to irrigation variable (nsub)
  INTEGER, ALLOCATABLE :: irrtype(:)               !<Irrigation type for class (0=no irrigation) (nclass)
  REAL,    ALLOCATABLE :: totaldemand_regsrc(:)    !<Total water demand at regional source (m3) (nsub)
  REAL,    ALLOCATABLE :: regionaldemand(:)        !<Water demand at regional source for the local subbasin (m3) (nsub)
  REAL,    ALLOCATABLE :: fieldneed(:,:)           !<Plant water deficit for irrigated classes (mm) (nclass,nsub)
  REAL,    ALLOCATABLE :: irrigation(:,:)          !<Water to be applied at this timestep for irrigated classes (mm) (nclass,nsub)
  REAL,    ALLOCATABLE :: cirrigation(:,:,:)       !<Concentration of irrigation object (numsubstances,nclass,nsub)

  PUBLIC irrtype
CONTAINS

  !>Initiation of irrigation. Calculated irrigation network.
  !>Called from initiate_model once per simulation
  !---------------------------------------------------------
  SUBROUTINE initiate_irrigation(n,nc)

    USE HYPEVARIABLES, ONLY : m_sswcorr

    USE MODVAR, ONLY : numsubstances

    !Argument declarations
    INTEGER, INTENT(IN) :: n      !<number of subbasin (nsub)
    INTEGER, INTENT(IN) :: nc     !<number of classes (nclass)
    
    !Local variables
    INTEGER nmax

    !>\b Algorithm \n
    !>First call, first simulation:
    IF(.NOT.ALLOCATED(fieldneed))THEN
      !>\li Allocate variables for irrigation calculation
      ALLOCATE(regirrwat_source(n))
      regirrwat_source = 0        !no regional source
      ALLOCATE(nregirrwat_collectors(n))
      nregirrwat_collectors = 0
      ALLOCATE(totaldemand_regsrc(n))
      ALLOCATE(regionaldemand(n))
      ALLOCATE(irrigation(nc,n))  
      ALLOCATE(cirrigation(numsubstances,nc,n))  
      ALLOCATE(fieldneed(nc,n))  
      ALLOCATE(irrindex(n))  
      irrindex = 0

      !>\li Determine regional irrigation coupling
      CALL set_regional_irrigation(n,nmax)

      !>\li Allocate and set variable for regional sources 
      IF(nmax>0)THEN
        IF(.NOT.ALLOCATED(regirrwat_collectors)) ALLOCATE(regirrwat_collectors(nmax,n))
        regirrwat_collectors = 0
        CALL set_regional_collectors(n)
      ENDIF
       
      !>\li Initialize season variable
      CALL set_irrigation_season_end()

    ENDIF

  END SUBROUTINE initiate_irrigation

  !>\brief Determine which classes are irrigated. 
  !!Calculate array with flag for irrigated classes. 
  !>
  !>\b Consequences Module irrigation_module variable irrtype is set
  !----------------------------------------------------------------------
  SUBROUTINE check_for_irrigated_classes(nc)

    USE MODVAR, ONLY : classdata,  &
                       cropindex,     &
                       cropirrdata

    !Argument declarations
    INTEGER, INTENT(IN) :: nc     !<number of classes (nclass)
    
    !Local variables
    INTEGER ir,j

    irrtype = 0

    DO j=1,nc
       DO ir=1,SIZE(cropindex,2)
          IF(classdata(j)%crop>0)THEN
             IF(cropindex(classdata(j)%crop,ir)>0)THEN
                IF(cropirrdata(cropindex(classdata(j)%crop,ir))%plantingdayno>0) irrtype(j)=1
             ENDIF
          ENDIF
          IF(classdata(j)%crop2>0)THEN
             IF(cropindex(classdata(j)%crop2,ir)>0)THEN
                IF(cropirrdata(cropindex(classdata(j)%crop2,ir))%plantingdayno>0) irrtype(j)=1
             ENDIF
          ENDIF
          IF(irrtype(j)==1)EXIT
       ENDDO
    ENDDO

  END SUBROUTINE check_for_irrigated_classes

  !>\brief Calculate index arrays for the regional source and the
  !>irrigation variable.
  !!Count number of collectors from each regional source.
  !!
  !>\b Consequences Module irrigation_module variables irrindex and regirrwat_source are set
  !----------------------------------------------------------------------------
  SUBROUTINE set_regional_irrigation(n,m)

    USE MODVAR, ONLY : basin,     &
                       irrigationsystem

    !Argument decalarations
    INTEGER, INTENT(IN)  :: n      !<number of subbasin (nsub)
    INTEGER, INTENT(OUT) :: m      !<maximum number of subbasins from a regional source
    
    !Local variables
    INTEGER i,j
    INTEGER dimirr    !size of irrigation indata variable
    INTEGER,ALLOCATABLE :: subindex(:)
    INTEGER numc(n)

    dimirr = SIZE(irrigationsystem)
    numc = 0
    ALLOCATE(subindex(dimirr))

    !Set index to find parameters in irrigation
    DO i = 1,n    !collector
       DO j = 1, dimirr
          IF(basin(i)%subid==irrigationsystem(j)%subid)THEN
             irrindex(i) = j           !this basin (i) has irrigation found on row j in the irrigationsystem object
             subindex(j) = i           !this row (j) belong to subbasin i
          ENDIF
       ENDDO
    ENDDO

    DO j = 1, dimirr
       IF(irrigationsystem(j)%regsourceid>0)THEN
          DO i = 1,n
             IF(basin(i)%subid==irrigationsystem(j)%regsourceid)THEN
                regirrwat_source(subindex(j)) = i      !source
                numc(i) = numc(i) + 1   !sum up number of collectors for this source
             ENDIF
          ENDDO
       ENDIF
    ENDDO

    m = MAXVAL(numc(:))
    IF(ALLOCATED(subindex)) DEALLOCATE(subindex)

  END SUBROUTINE set_regional_irrigation

  !>\brief Determine collectors for all subbasins, i.e. which
  !>subbasins has this as a regional source of irrigation
  !!
  !>\b Consequences Module irrigation_module variables regirrwat_collectors 
  !> and nregirrwat_collectors are set.
  !---------------------------------------------------------------------------
  SUBROUTINE set_regional_collectors(n)

    !Argument declarations
    INTEGER, INTENT(IN)  :: n      !<number of subbasin
    
    !Local variables
    INTEGER i

    DO i = 1,n
       IF(regirrwat_source(i)>0)THEN
          nregirrwat_collectors(regirrwat_source(i)) = nregirrwat_collectors(regirrwat_source(i)) + 1
          regirrwat_collectors(nregirrwat_collectors(regirrwat_source(i)),regirrwat_source(i)) = i
       ENDIF
    ENDDO

  END SUBROUTINE set_regional_collectors

  !\brief Initiation of irrigation variables for the time step calculation.
  !!Called from model once per timestep. 
  !!
  !>\b Consequences Module irrigation_module variables totaldemand_regsrc, 
  !!regionaldemand, fieldneed, irrigation and cirrigation are set
  !---------------------------------------------------------------------------
  SUBROUTINE initiate_timestep_irrigation(nsubst,miscstate)

    USE MODVAR, ONLY : doirrigation
    
    !Argument declarations
    INTEGER,INTENT(IN)    :: nsubst   !<Number of substances
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    !Irrigation included in model set-up?
    IF(.NOT.doirrigation) RETURN

    !Initialisations
    totaldemand_regsrc = 0.
    regionaldemand = 0.
    fieldneed = 0.

    !Irrigation calculated last time step is now to be applied    
    irrigation(:,:) = miscstate%nextirrigation(:,:) 
    miscstate%nextirrigation(:,:) = 0.
    IF(nsubst>0)THEN
      cirrigation(:,:,:) = miscstate%cnextirrigation(:,:,:)
      miscstate%cnextirrigation(:,:,:) = 0.
    ENDIF

  END SUBROUTINE initiate_timestep_irrigation

  !>\brief Calculate irrigation water demand of this class this timestep.
  !!Called from model once per timestep/subbasin/land class.
  !!
  !>\b Reference ModelDescription Water management (Irrigation - Irrigation water demand)
  !!
  !>\b Consequences Module irrigation_module variable fieldneed is set.
  !---------------------------------------------------------------------------
  SUBROUTINE calculate_irrigation_water_demand(i,j,dn,area,sswcorrpar,immdeppar,  &
                        iwdfracpar,wdpar,soil,wpl,fcl,epl,epot,epotfrac,iwdv)

    USE MODVAR, ONLY : classdata,      &
                       basin,          &
                       cropindex,      &
                       cropirrdata,    &
                       irrigationsystem,  &
                       maxsoillayers

    !Argument declarations
    INTEGER, INTENT(IN) :: i        !<index of subbasin
    INTEGER, INTENT(IN) :: j        !<index of slc-class
    INTEGER, INTENT(IN) :: dn       !<current pseudo day number
    REAL, INTENT(IN)    :: area     !<class area (km2)
    REAL, INTENT(IN)    :: sswcorrpar   !<parameter for scaling of irrigation threshold (-)
    REAL, INTENT(IN)    :: immdeppar    !<parameter for target immersion depth (mm)
    REAL, INTENT(IN)    :: iwdfracpar   !<parameter for scaling of irrigation water demand relative to threshold (-)
    REAL, INTENT(IN)    :: wdpar        !<parameter for constant water demand (mm/timestep)
    REAL, INTENT(IN)    :: soil(maxsoillayers)  !<soil moisture per soil layer (mm)
    REAL, INTENT(IN)    :: wpl(maxsoillayers)   !<wilting point per soil layer (mm)
    REAL, INTENT(IN)    :: fcl(maxsoillayers)   !<field capacity per soil layer (mm)
    REAL, INTENT(IN)    :: epl(maxsoillayers)   !<efficient porosity per soil layer (mm)
    REAL, INTENT(IN)    :: epot            !<potential evapotranspiration (mm/timestep)
    REAL, INTENT(IN)    :: epotfrac(2)     !<relative distribution of potential evaporation between upper two soil layers (-)
    REAL, INTENT(OUT)   :: iwdv            !<subbasin irrigation water demand for the soil (volume) (m3/timestep)
    
    !Local variables
    INTEGER cropiindex            !index in cropirrdata for current subbasin-class (if irrigated)
    LOGICAL immcrop, immseason    !flags for immersed crop
    REAL sm1,sm2        !soil moisture in layer 1 and layer 2 (mm)    
    REAL pottrans, pottrans1,pottrans2  !potential transpiration (mm), pottrans is for both layers whereas pottrans2 is for layer 2 and pottrans1 for layer 1
    REAL kcb       !crop basal coefficient
    REAL ssw       !original soil moisture threshold for irrigation (fraction of plant available soil water, used same for each layer)
    REAL sswcorr   !adjusted soil moisture threshold for irrigation, adjusted with sswcorrpar (fraction of plant available soil water, used same for each layer)
    REAL iwd,iwd1,iwd2       !irrigation water demand for the soil (mm/timestep). iwd1 and iwd2 are for layers 1 and 2 respectively while iwd is for the total for the rootzone

    !>\b Algorithm \n
    !Initialisations
    iwd = 0.
    iwd1 = 0.
    iwd2 = 0.
    iwdv = 0.

    !>Check if the class is irrigated today
    IF(.NOT.(irrtype(j)==1))RETURN  !Irrigated classes har irrigation type 1

    !Find irrigation data on crop
    cropiindex = cropindex(classdata(j)%crop,basin(i)%region)

    !Irrigation season?
    IF(.NOT.irrigation_season(cropiindex,dn)) RETURN

    !Immersion crop and season?
    immcrop = .FALSE.
    immseason = .FALSE.
    IF(cropirrdata(cropiindex)%imm_start>0) immcrop = .TRUE.
    IF(immcrop) immseason = immersion_season(cropiindex,dn)

    !>Depending on irrigation and crop type:
    IF(immcrop .AND. immseason)THEN
      !>Calculate irrigation water demand for immersed crop irrigation
      iwd = wpl(1) + fcl(1) + epl(1) + immdeppar - soil(1)     !mm/timestep
      IF(iwd<0.) iwd = 0.
    ELSE  
      !>Calculate irrigation water demand by soil moisture based irrigation
      !>\li Calculate current basal crop coefficient, kcb-value
      CALL calculate_kcb(cropiindex,dn,kcb)  

      !>\li Calculate soil moisture threshold for irrigation
      !Generic for both layers
      pottrans = kcb*epot                                                    !potential transpiration (all layers)
      ssw = 1. - (cropirrdata(cropiindex)%dlref + 0.04*(5. - pottrans/0.95)) !Soil water stress threshold, original
      IF(ssw<0.2) ssw=0.2
      IF(ssw>0.9) ssw=0.9
      sswcorr = ssw*sswcorrpar
      IF(sswcorr>1.) sswcorr = 1. !soil moisture threshold for irrigation, possibly modified by a user set parameter 

      !Layer 1
      pottrans1 = epot*epotfrac(1)*kcb                            ! Potential transpiration of the first layer  
      sm1 = soil(1) - wpl(1)
      IF(sm1<0) sm1 = 0.
      IF(sm1<sswcorr*fcl(1))THEN
        !>\li Irrigation water demand Layer 1, mm/timestep
        IF(irrigationsystem(irrindex(i))%demandtype==1)THEN         !Alternative: constant application, only applying in Layer 1 
          iwd1 = wdpar
        ELSEIF(irrigationsystem(irrindex(i))%demandtype==2)THEN     !Alternative: soil moisture deficit, fill the field capacity
          iwd1 = fcl(1) - sm1
        ELSEIF(irrigationsystem(irrindex(i))%demandtype==3)THEN     !Alternative: threshold dependent
          iwd1 = MIN((sswcorr*fcl(1) - sm1)*iwdfracpar,fcl(1) - sm1)    
        ENDIF
      ENDIF

      !Layer 2            
      pottrans2 = epot*epotfrac(2)*kcb
      sm2 = 0.
      IF(fcl(2)+wpl(2)>0) sm2 = soil(2) - wpl(2)        
      IF(sm2<0) sm2 = 0.        
      IF(sm2<sswcorr*fcl(2))THEN
        !>\li Irrigation water demand Layer 2, mm/timestep
        IF(irrigationsystem(irrindex(i))%demandtype==1)THEN         !Alternative: constant application. Nothing here, only specifying in Layer one
          iwd2 = 0.0
        ELSEIF(irrigationsystem(irrindex(i))%demandtype==2)THEN     !Alternative: soil moisture deficit, fill the field capacity
          iwd2 = fcl(2) - sm2
        ELSEIF(irrigationsystem(irrindex(i))%demandtype==3)THEN     !Alternative: threshold dependent
          iwd2 = MIN((sswcorr*fcl(2) - sm2)*iwdfracpar,fcl(2) - sm2)     
          !          ELSEIF(irrigationsystem(irrindex(i))%demandtype==4)THEN     !Alternative: difference between potential transpiration and current soil moisture (proxy for evap2)
          !            iwd2 = MAX(0.0,MIN(pottrans2-sm2,fcl(2)-sm2))            OBS WRONG SO FAR - SHOULD BE ASSESSED AGAINST EVAP2*KCB INSTEAD OF sm2!
        ENDIF
      ENDIF

      !>\li Sum irrigation water demand for the upper two soil layers
      iwd = iwd1 + iwd2
    ENDIF

    !>Set output variables
    iwdv = iwd*area*1000.  !Total irrigation water demand (m3/timestep)
    fieldneed(j,i) = iwd   !Classwise irrigation water demand  (mm/timestep)

  END SUBROUTINE calculate_irrigation_water_demand

  !>\brief Calculate basal crop coefficient, kcb-parameter 
  !!Equation 66 from Allen et al. (1998) FAO Irrigation and Drainage
  !!Paper No.56
  !>\b Reference ModelDescription Water management (Irrigation - Irrigation water demand)
  !------------------------------------------------------------------------------
  SUBROUTINE calculate_kcb(k,dn,kcb)
  
    USE MODVAR, ONLY : cropirrdata

    !Argument declarations
    INTEGER, INTENT(IN)  :: k        !<index in cropdata
    INTEGER, INTENT(IN)  :: dn       !<current day number
    REAL, INTENT(OUT)    :: kcb      !<basal crop coefficient
    
    !Local variables
    INTEGER growthday     !day in growing season

    growthday = dn - cropirrdata(k)%plantingdayno + 1
    IF(growthday <= 0) growthday = growthday + 365  
    IF(growthday<=cropirrdata(k)%lengthini)THEN
       kcb = cropirrdata(k)%kcbini
    ELSEIF(growthday<=cropirrdata(k)%lengthini+cropirrdata(k)%lengthdev)THEN
       kcb = cropirrdata(k)%kcbini + (growthday - cropirrdata(k)%lengthini)*(cropirrdata(k)%kcbmid-cropirrdata(k)%kcbini)/REAL(cropirrdata(k)%lengthdev)
    ELSEIF(growthday<cropirrdata(k)%lengthini+cropirrdata(k)%lengthdev+cropirrdata(k)%lengthmid)THEN
       kcb = cropirrdata(k)%kcbmid 
    ELSEIF(growthday<cropirrdata(k)%lengthini+cropirrdata(k)%lengthdev+cropirrdata(k)%lengthmid+cropirrdata(k)%lengthlate)THEN
       kcb = cropirrdata(k)%kcbmid + (growthday - (cropirrdata(k)%lengthini+cropirrdata(k)%lengthdev+cropirrdata(k)%lengthmid))/REAL(cropirrdata(k)%lengthlate)*(cropirrdata(k)%kcbend-cropirrdata(k)%kcbmid)
    ELSE
       kcb = 0.
    ENDIF

  END SUBROUTINE calculate_kcb

  !>Calculate and set end of irrigation season and immersion season
  !>
  !>\b Consequences Module modvar variable cropirrdata is changed.
  !----------------------------------------------------------------
  SUBROUTINE set_irrigation_season_end()

    USE MODVAR, ONLY : cropirrdata,   & !OUT
         ncrop         
              
    !Local variables
    INTEGER k

    DO k = 1,ncrop   
       IF(cropirrdata(k)%plantingdayno>0)THEN
          cropirrdata(k)%season_end = cropirrdata(k)%plantingdayno + cropirrdata(k)%lengthini + &
               cropirrdata(k)%lengthdev + cropirrdata(k)%lengthmid + cropirrdata(k)%lengthlate
       ELSE
          cropirrdata(k)%season_end = 0
       ENDIF
       IF(cropirrdata(k)%imm_start>0)THEN
          IF(cropirrdata(k)%imm_end < cropirrdata(k)%imm_start) cropirrdata(k)%imm_end = cropirrdata(k)%imm_end + 365
       ENDIF
    ENDDO

  END SUBROUTINE set_irrigation_season_end

  !>Function to determine if it is irrigation season
  !-------------------------------------------------------------------
  LOGICAL FUNCTION irrigation_season(k,dn)

    USE MODVAR, ONLY : cropirrdata

    !Argument declarations
    INTEGER, INTENT(IN) :: k      !<index in cropdata
    INTEGER, INTENT(IN) :: dn     !<day number
    
    !Local variables
    LOGICAL status

    status = .FALSE.
    IF(cropirrdata(k)%season_end > 365)THEN
       IF(dn>=cropirrdata(k)%plantingdayno) status = .TRUE.
       IF(dn<=cropirrdata(k)%season_end - 365) status = .TRUE.
    ELSE
       IF(dn>=cropirrdata(k)%plantingdayno .AND. dn<=cropirrdata(k)%season_end) status = .TRUE.
    ENDIF
    irrigation_season = status

  END FUNCTION irrigation_season

  !>Function to determine if it is immersion season
  !--------------------------------------------------------------------
  LOGICAL FUNCTION immersion_season(k,dn)

    USE MODVAR, ONLY : cropirrdata

    !Argument declarations
    INTEGER, INTENT(IN) :: k      !<index in cropdata
    INTEGER, INTENT(IN) :: dn     !<day number
    
    !Local variables
    LOGICAL status

    status = .FALSE.
    IF(cropirrdata(k)%imm_end > 365)THEN
       IF(dn>=cropirrdata(k)%imm_start) status = .TRUE.
       IF(dn<=cropirrdata(k)%imm_end - 365) status = .TRUE.
    ELSE
       IF(dn>=cropirrdata(k)%imm_start .AND. dn<=cropirrdata(k)%imm_end) status = .TRUE.
    ENDIF
    immersion_season = status

  END FUNCTION immersion_season

  !>\brief Add irrigation water to the soil of one class.
  !!Called from model once per timestep/subbasin/land class.
  !!
  !>\b Reference ModelDescription Water management (Irrigation - Irrigation water application)
  !--------------------------------------------------------------------
  SUBROUTINE apply_irrigation(i,j,layerfrac,porvolume,soilstate,appirr,appirrmass)

    USE MODVAR, ONLY : numsubstances,&
                       classbasin,   &
                       basin

    !Argument declarations
    INTEGER, INTENT(IN) :: i                         !<index of subbasin
    INTEGER, INTENT(IN) :: j                         !<index of slc-class
    REAL, INTENT(IN)    :: layerfrac(2)              !<division of irrigation between two top soil layers (fraction in each)
    REAL, INTENT(IN)    :: porvolume                 !<pore volume of soil layer 2 (mm)
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate  !<Soil states
    REAL, INTENT(OUT)   :: appirr(2)                 !<applied irrigation (mm)
    REAL, INTENT(OUT)   :: appirrmass(numsubstances) !<mass of applied irrigation (kg/timestep if the conc is in mg/L)    
    
    !Local variables
    REAL apply1, apply2    !water applied to soil layer 1 and 2
    REAL maxfill           !available room in soil layer 2 for irrigation water

    !>\b Algorithm \n
    !>Initialisations
    appirr = 0.
    appirrmass = 0.

    IF(irrigation(j,i)>0)THEN  

      !>Add irrigation to second soillayer
      apply2 = irrigation(j,i)*layerfrac(2)
      maxfill = porvolume - soilstate%water(2,j,i)
      IF(apply2>maxfill) apply2 = maxfill
      IF(numsubstances>0) soilstate%conc(:,2,j,i) = (soilstate%water(2,j,i)*soilstate%conc(:,2,j,i) + apply2*cirrigation(:,j,i)) / (soilstate%water(2,j,i) + apply2) ! (old mass + applied mass)/(old volume + applied volume)
      soilstate%water(2,j,i) = soilstate%water(2,j,i) + apply2

      !>Add irrigation to first layer
      apply1 = irrigation(j,i) - apply2
      IF(numsubstances>0) soilstate%conc(:,1,j,i) = (soilstate%water(1,j,i)*soilstate%conc(:,1,j,i) + apply1*cirrigation(:,j,i)) / (soilstate%water(1,j,i) + apply1)
      soilstate%water(1,j,i) = soilstate%water(1,j,i) + apply1

      !>Accumulate applied water for subbasin print out
      appirr(1) = apply1
      appirr(2) = apply2
      IF(numsubstances>0) appirrmass = cirrigation(:,j,i) * irrigation(j,i) * basin(i)%area * classbasin(i,j)%part * 1.E-6 ! mass of the the application (mg/L -> kg for e.g. nitrogen (mg/L=g/m3->g->kg)
    ENDIF

  END SUBROUTINE apply_irrigation

  !>\brief Calculate irrigation water withdrawal and losses on the way
  !>to the irrigated field for a subbasin.
  !!Called from model once per timestep/subbasin. 
  !!
  !>\b Reference ModelDescription Water management (Irrigation - Irrigation water withdrawal, Irrigation water application)
  !!
  !>\b Consequences Module irrigation_module variables
  !!totaldemand_regsrc, regionaldemand, and fieldneed may change.
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_irrigation(i,na,irrigationpar,q,lakew,riverv,aquiferv, &
       accpw,irrnetloss,gwabstr,ldabstr,lrabstr,rsabstr,cq,   &
       irrsinkmass,cilake,colake,criver,caquifer,soilstate, &
       nextirrigation,cnextirrigation,wblirrflows,wbrirrflows,wbrirrevap,aquiferloss)    

    USE HYPEVARIABLES, ONLY : epotdist                              
    USE MODVAR, ONLY : nsub,         &
                       nclass,       &
                       path,         &
                       seconds_per_timestep,  &
                       slc_ilake,    &
                       slc_olake,    &
                       classbasin,   &
                       basin,        &
                       irrigationsystem, &
                       numsubstances,&
                       irrunlimited

    !Argument declarations
    INTEGER, INTENT(IN) :: i                 !<index of current subbasin
    INTEGER, INTENT(IN) :: na                !<number of aquifers
    REAL, INTENT(IN)    :: irrigationpar(5)         !<parameters
    REAL, INTENT(INOUT) :: q                 !<flow in main river (m3/s)
    REAL, INTENT(INOUT) :: lakew(2)          !<lake water stage (ilake,olake) (mm)
    REAL, INTENT(INOUT) :: riverv            !<river volume (type 2) (m3)
    REAL, INTENT(INOUT) :: aquiferv(na)      !<aquifer volume (m3)
    REAL, INTENT(IN)    :: accpw             !<irrigation water need for this subbasin (m3/timestep)  =sum(pwn(i,:))
    REAL, INTENT(INOUT) :: irrnetloss(nsub)  !<accumulated irrigation losses (m3/timestep)
    REAL, INTENT(OUT)   :: gwabstr           !<groundwater abstraction for irrigation (m3/timestep) 
    REAL, INTENT(OUT)   :: ldabstr           !<local dam water abstraction for irrigation (m3/timestep)
    REAL, INTENT(OUT)   :: lrabstr           !<local river water abstraction for irrigation (m3/timestep)
    REAL, INTENT(OUT)   :: rsabstr           !<regional surface water abstraction for irrigation in other subbasins (m3/timestep)
    REAL, INTENT(IN)    :: cq(numsubstances) !<concentration of flow in main river (typically mg/L)
    REAL, INTENT(OUT)   :: irrsinkmass(numsubstances) !<mass of substances in the irrigation sink in this subbasin (typically in kg)
    REAL, INTENT(IN)    :: cilake(numsubstances) !<concentration laketype 1 (typically mg/L)
    REAL, INTENT(IN)    :: colake(numsubstances) !<concentration laketype 2 (typically mg/L)
    REAL, INTENT(IN)    :: criver(numsubstances) !<concentration river type 2 (typically mg/L)
    REAL, INTENT(IN)    :: caquifer(numsubstances,na) !<concentration aquifer (typically mg/L)
    TYPE(soilstatetype),INTENT(IN)  :: soilstate !<Soil states
    REAL, INTENT(INOUT) :: nextirrigation(nclass,nsub)       !<irrigation for next time step (mm/timestep)
    REAL, INTENT(INOUT) :: cnextirrigation(numsubstances,nclass,nsub) !<concentration irrigation next time step
    REAL, INTENT(OUT)   :: wblirrflows(7)       !<irrigations flows needed for water balance output (m3/timestep)
    REAL, INTENT(INOUT) :: wbrirrflows(4,nsub)  !<irrigations flows needed for water balance output, regional sources (m3/timestep)
    REAL, INTENT(OUT)   :: wbrirrevap(nsub)     !<irrigations flows needed for water balance output, regional sources evaporation (m3/timestep)
    REAL, INTENT(INOUT) :: aquiferloss(na)      !<accumulated aquifer removal (m3/timestep)

    !Local variables
    INTEGER ia,j,k,iadd
    REAL help
    REAL org_q,org_rv                     !river inflow and volume before irrigation water is removed
    REAL org_lake                         !lake waterstage before irrigation water is removed
    REAL groundwater_part, surfacewater_part    !these will be subbasin or other region dependent
    REAL local_efficiency                 !efficiency of irrigation (field and local irrigation net)
    REAL twn                              !total water need for the subbasin
    REAL twn_left                         !remaining water need for the subbasin (total)
    REAL twn_left_gw                      !remaining groundwater need for the subbasin
    REAL twn_left_gw_comp                 !The part of the remaining groundwater need for the subbasin which can be withdrawn through local source compensation
    REAL twn_left_sw                      !remaining surface water need for the subbasin
    REAL sw_left                          !sum of available surface water after initial local withdrawal and before source compensation
    REAL local_dw_avail                   !local dam water available for irrigation (olake) (m3)
    REAL local_dw2_avail                  !local dam water available for irrigation (ilake) (m3)
    REAL local_sw_avail                   !local surface water available for irrigation (main river) (m3)
    REAL local_rv_avail                   !local surface water available for irrigation (river volume) (m3)
    REAL local_gw_avail                   !local aquifer water available for irrigation (m3)
    REAL local_dw_irr                     !water for irrigation from local dam (olake) (m3)
    REAL local_dw2_irr                    !water for irrigation from local dam (ilake) (m3)
    REAL local_sw_irr                     !water for irrigation from river inflow (m3)
    REAL local_rv_irr                     !water for irrigation from river volume (m3)
    REAL local_gw_irr                     !water for irrigation from local groundwater (m3)
    REAL clocal_dw_irr(numsubstances)     !concentration of local_dw_irr (mg/L typically)
    REAL clocal_dw2_irr(numsubstances)    !concentration of local_dw2_irr (mg/L typically)
    REAL clocal_rv_irr(numsubstances)     !concentration of local_rv_irr (mg/L typically)
    REAL clocal_gw_irr(numsubstances)     !concentration of local_gw_irr (mg/L typically)
    REAL local_add_irr                    !water for irrigation that actually reaches the soil (sw+gw) (m3)
    REAL local_add_irr_factor             !fraction of what the plant need that is gets from local surface water (-)
    REAL clocal_add_irr(numsubstances)    !concentration of local_add_irr (mg/L typically)
    REAL regional_dw_avail                !dam water available for irrigation in the regional surface water (olake) (m3)
    REAL regional_sw_avail                !water available for irrigation in the regional surface water (main river inflow) (m3)
    REAL regional_rv_avail                !water available for irrigation in the regional surface water (main river volume)(m3)
    REAL regional_dw_irr                  !water for irrigation taken from regional surface water (olake) (m3)
    REAL regional_sw_irr                  !water for irrigation taken from regional surface water (main river inflow)(m3)
    REAL regional_rv_irr                  !water for irrigation taken from regional surface water (main river volume)(m3)
    REAL regional_dw_fraction             !fraction of water for irrigation taken from regional surface water (lake part)
    REAL regional_rw_fraction             !fraction of water for irrigation taken from regional surface water (river part)
    REAL cregional_dw_irr(numsubstances)
    REAL cregional_rv_irr(numsubstances)
    REAL regional_add_irr                 !water for irrigation taken from regional surface water (all)(m3)
    REAL regional_add_irr_factor          !fraction of water need that is fullfilled 
    REAL cregional_add_irr(numsubstances)
    REAL watertobasin                     !what the local subbasin gets from the regional source (m3)
    REAL sum_wtb                          !sum of what the local subbasin gets from the regional source (m3)
    REAL watertofield                     !what the field of local subbasin gets from the regional source (m3)
    REAL cwatertofield(numsubstances)                     
    REAL previous                         !help variable, value before a certain operation 
    REAL add                              !help variable, value added in a certain operation
    REAL a3,csoil3(numsubstances)         !help variables for average third soil layer concentration
    REAL revap                            !help variable, evaporation losses regional source/regional network
    REAL regirrpar         !parameter scaling the water withdrawal from regional source (-)
    REAL pirrspar          !parameter scaling the irrigation water abstractions for surface water(-)
    REAL pirrgpar          !parameter scaling the irrigation water abstractions for groundwater(-)
    REAL cirrsinkpar       !parameter: concentration reduction fraction in settlement tanks (at irrigation abstraction points)
    REAL irrcomppar        !parameter regulating the degree of comensation allowed between local groundwater/surfacewater sources (1 = allow 100% compensation, 0=do not allow compensation, default)

    REAL,PARAMETER :: safe_factor = 0.9999    !factor multiplied by available water source to avoid zero (negative) volumes

    !Initialisations output variables
    gwabstr = 0.     !OUT
    ldabstr = 0.     !OUT
    lrabstr = 0.     !OUT
    rsabstr = 0.     !OUT
    wblirrflows = 0. !OUT
    wbrirrevap = 0.  !OUT
    
    !Initiation local parameters
    regirrpar = irrigationpar(1)        !parameter scaling the water withdrawal from regional source (-)
    pirrspar = irrigationpar(2)         !parameter scaling the irrigation water abstractions for surface water(-)
    pirrgpar = irrigationpar(3)         !parameter scaling the irrigation water abstractions for groundwater(-)
    cirrsinkpar = irrigationpar(4)      !parameter: concentration reduction fraction in settlement tanks (at irrigation abstraction points)
    irrcomppar = irrigationpar(5)       !parameter regulating the degree of comensation allowed between local groundwater/surfacewater sources (1 = allow 100% compensation, 0=do not allow compensation, default)

    !Initiations for regional irrigation source
    org_lake = lakew(2)                              
    org_q = q
    org_rv = riverv                             
    local_dw_irr = 0.
    local_sw_irr = 0.
    local_rv_irr = 0.
    cregional_dw_irr = colake
    cregional_rv_irr = criver
    irrsinkmass = 0.

    !Irrigation with water from local sources
    !---------------------------------------- 
    IF(accpw>0)THEN

      !!Special case of unlimited irrigation
       IF(irrunlimited) THEN
          local_add_irr_factor = 1.0     !Unlimited irrigation always allows full fieldneed irrigation
          irrnetloss(i) = 0.0            !No irrigation losses with unlimited irrigation         

          !Save irrigation amount and concentration to classes for next time step, spread in accordance with fieldneed.
          DO j = 1,nclass
             previous = nextirrigation(j,i)                                                    !Current depth in the relevant nextirrigation object
             nextirrigation(j,i) = nextirrigation(j,i) + local_add_irr_factor * fieldneed(j,i) !Accumulated depth (mm) to be added, scaling fieldneed to local availability
             add = nextirrigation(j,i)-previous                                                !Amount to be added in this iteration of the loop 
             clocal_add_irr = soilstate%conc(:,2,j,i)*epotdist(2,j) + soilstate%conc(:,1,j,i)*epotdist(1,j)     ! Concentration of the added water = concentration of the soil layers where the water is put (i.e. no conc. change), scaled by the proportion of water ending up in each layer
             IF (nextirrigation(j,i)>0) cnextirrigation(:,j,i) = (cnextirrigation(:,j,i)*previous + clocal_add_irr*add)/ nextirrigation(j,i) ! New concentration (as below but here based on the existing concentrations)
          ENDDO
          fieldneed(:,i)=0.   !Zero the field need since you've provided everything in unlimited irrigation
          wblirrflows(6) = accpw
          RETURN  !go back to calling subroutine
       ENDIF

      !!Case of normal source-limited irrigation

       !Initialisations parameters
       groundwater_part = irrigationsystem(irrindex(i))%gw_part
       surfacewater_part = irrigationsystem(irrindex(i))%sw_part
       local_efficiency = irrigationsystem(irrindex(i))%local_eff

       !Calculate local water need for the subbasin (m3), adjusting for local inefficiencies
       twn = accpw / local_efficiency

       !Calculate local water availability from surface water (dam, river inflow, river volume) for irrigation
       local_dw_avail  = 0.
       local_dw2_avail = 0.
       local_sw_avail  = 0.
       local_rv_avail  = 0.
       local_dw_irr  = 0.
       local_dw2_irr = 0.
       local_sw_irr  = 0.
       local_rv_irr  = 0.
       clocal_dw_irr = colake  ! Initiating concentrations here since else they will have no value in the equation below, and even with a non-zero value if they are not to be used, they will not affect the calc since the volumes are 0
       clocal_dw2_irr = cilake          
       clocal_rv_irr = criver      

       IF(surfacewater_part>0)THEN
         !Calculate local water availability from irrigation dam (olake) (m3)
         IF(irrigationsystem(irrindex(i))%dam .AND. slc_olake>0)THEN
           IF(classbasin(i,slc_olake)%part>0.)THEN
             local_dw_avail = lakew(2) *1.E-3 * classbasin(i,slc_olake)%part * basin(i)%area * safe_factor
             local_dw_irr = MIN(twn*surfacewater_part,local_dw_avail)
             local_dw_avail = local_dw_avail - local_dw_irr
           ENDIF
         ENDIF

         !Calculate local water availability from irrigation dam (ilake) (m3)
         IF(irrigationsystem(irrindex(i))%dam .AND. slc_ilake>0)THEN
           IF(classbasin(i,slc_ilake)%part>0.)THEN
             local_dw2_avail = lakew(1) *1.E-3 * classbasin(i,slc_ilake)%part * basin(i)%area * safe_factor
             local_dw2_irr = MIN(twn*surfacewater_part-local_dw_irr,local_dw2_avail)
             local_dw2_avail = local_dw2_avail - local_dw2_irr
           ENDIF
         ENDIF

          !Calculate local water availability from river inflow (m3)
          local_sw_avail = q * seconds_per_timestep * safe_factor
          IF(local_dw_avail==0.)THEN  ! Jafet: This is primarily a test if there is an olake in the area, and secondly if there is any water demand remaining. Ignores ilakes, but calc. below still accounts for it.       
             local_sw_irr = MIN(twn*surfacewater_part - local_dw_irr - local_dw2_irr,local_sw_avail)
             local_sw_avail = local_sw_avail - local_sw_irr
          ENDIF

          !Calculate local water availability from river volume (m3)
          local_rv_avail = riverv * safe_factor
          IF(local_dw_avail==0. .AND. local_sw_avail==0.)THEN
             local_rv_irr = MIN(twn*surfacewater_part - local_dw_irr - local_dw2_irr - local_sw_irr,local_rv_avail)
             local_rv_avail = local_rv_avail - local_rv_irr
          ENDIF
       ENDIF

       !Calculate local water availability from groundwater (m3)
       local_gw_avail  = 0.
       local_gw_irr = 0.
       clocal_gw_irr = 0.

       !Version withdrawing from an deep aquifer source (simulated or external)
       IF(groundwater_part>0)THEN
         IF(na>0)THEN
           IF(path(i)%rechargebasin.OR.path(i)%recievefraction>0.)THEN
             ia = path(i)%aquid
             local_gw_avail = (aquiferv(ia)-aquiferloss(ia)) * safe_factor
             IF(local_gw_avail>0.)THEN
               local_gw_irr = MIN(twn*groundwater_part,local_gw_avail)
               clocal_gw_irr = caquifer(:,ia)
               local_gw_avail = local_gw_avail - local_gw_irr
             ENDIF
           ELSE
             local_gw_irr = twn*groundwater_part   ! Amount available equals the demand (i.e. unlimited source)
             csoil3 = 0.; a3 = 0.
             DO j = 1,nclass
               IF(classbasin(i,j)%part>0)THEN
                 csoil3 = csoil3 + soilstate%conc(:,3,j,i) * classbasin(i,j)%part
                 a3 = a3 + classbasin(i,j)%part
               ENDIF
             ENDDO
             IF(a3>0) clocal_gw_irr = csoil3/a3       ! Assume average concentration of third soil layer, Change to deep aquifer when available
           ENDIF
         ELSE
           local_gw_irr = twn*groundwater_part   ! Amount available equals the demand (i.e. unlimited source)
           csoil3 = 0.; a3 = 0.
           DO j = 1,nclass
             IF(classbasin(i,j)%part>0)THEN
                csoil3 = csoil3 + soilstate%conc(:,3,j,i) * classbasin(i,j)%part
                a3 = a3 + classbasin(i,j)%part
             ENDIF
           ENDDO
           IF(a3>0) clocal_gw_irr = csoil3/a3       ! Assume average concentration of third soil layer, Change to deep aquifer when available
         ENDIF
       ENDIF

       !Calculate if/how much water demand remains from surface water and groundwater. Has the water need been met by local sources?
       twn_left_gw = twn*groundwater_part - local_gw_irr
       twn_left_sw = twn*surfacewater_part - local_dw_irr - local_dw2_irr - local_sw_irr - local_rv_irr 
       sw_left = local_dw_avail + local_dw2_avail + local_sw_avail + local_rv_avail

       !Compensation from the other local source 
       IF (irrcomppar>0)THEN                           ! If compensation is allowed in par.txt
          IF(twn_left_gw>0 .AND. sw_left>0)THEN        ! If groundwater demand remains, and surface water is available. (Basically never happens if using unlimited GW source. Kept for future changes)
             twn_left_gw_comp =  twn_left_gw*irrcomppar ! Define how much can be withdrawn through compensation
             IF(local_dw_avail>0)THEN                   !Check for more water available in irrigation dam (olake)
                help = MIN(local_dw_avail,twn_left_gw_comp)
                local_dw_irr = local_dw_irr + help
                local_dw_avail = local_dw_avail - help
                twn_left_gw = twn_left_gw - help
                twn_left_gw_comp = twn_left_gw_comp - help
             ENDIF
             IF(twn_left_gw_comp>0 .AND. local_dw2_avail>0)THEN    !Check for more water available in irrigation dam (ilake)
                help = MIN(local_dw2_avail,twn_left_gw_comp)
                local_dw2_irr = local_dw2_irr + help
                local_dw2_avail = local_dw2_avail - help
                twn_left_gw = twn_left_gw - help
                twn_left_gw_comp = twn_left_gw_comp - help          
             ENDIF
             IF(twn_left_gw_comp>0 .AND. local_sw_avail>0)THEN    !Check for more water available as river inflow
                help = MIN(local_sw_avail,twn_left_gw_comp)
                local_sw_irr = local_sw_irr + help
                local_sw_avail = local_sw_avail - help
                twn_left_gw = twn_left_gw - help
                twn_left_gw_comp = twn_left_gw_comp - help          
             ENDIF
             IF(twn_left_gw_comp>0 .AND. local_rv_avail>0)THEN    !Check for more water available in river
                help = MIN(local_rv_avail,twn_left_gw_comp)
                local_rv_irr = local_rv_irr + help
                local_rv_avail = local_rv_avail - help
                twn_left_gw = twn_left_gw - help
                twn_left_gw_comp = twn_left_gw_comp - help          
             ENDIF
          ELSEIF(twn_left_sw>0)THEN      ! If surface water demand remains, compensate from GW source
             help = twn_left_sw*irrcomppar
             IF(na>0)THEN
               IF(path(i)%rechargebasin.OR.path(i)%recievefraction>0.)THEN
                 help = MIN(help,local_gw_avail)
                 local_gw_irr = local_gw_irr + help
                 local_gw_avail = local_gw_avail - help
                 twn_left_sw = twn_left_sw - help
               ELSE
                 local_gw_irr = local_gw_irr + help
                 twn_left_sw = twn_left_sw - help
               ENDIF
             ELSE  
               local_gw_irr = local_gw_irr + help
               twn_left_sw = twn_left_sw - help
             ENDIF
          ENDIF
       ENDIF

       !Total remaining local demand
       twn_left = twn_left_gw + twn_left_sw    !Normal case

       !Remove water from the sources. Apply general pirrspar and pirrgpar reduction.
       IF(local_dw_irr>0)THEN
         lakew(2) = lakew(2) - pirrspar * local_dw_irr / (classbasin(i,slc_olake)%part * basin(i)%area) * 1.E3   !mm
         ldabstr = pirrspar * local_dw_irr    !m3
         wblirrflows(3) = ldabstr
       ENDIF
       IF(local_dw2_irr>0)THEN
         lakew(1) = lakew(1) - pirrspar * local_dw2_irr / (classbasin(i,slc_ilake)%part * basin(i)%area) * 1.E3   !mm
         ldabstr = ldabstr + pirrspar * local_dw2_irr    !m3
         wblirrflows(2) = pirrspar * local_dw2_irr
       ENDIF
       IF(local_sw_irr>0.)THEN
         q = q - pirrspar * local_sw_irr/seconds_per_timestep     !m3/s
         lrabstr = pirrspar * local_sw_irr    !m3
         wblirrflows(4) =  lrabstr
       ENDIF
       IF(local_rv_irr>0)THEN
         riverv = riverv - pirrspar * local_rv_irr     !m3
         lrabstr = lrabstr + pirrspar * local_rv_irr    !m3
         wblirrflows(4) = lrabstr   !both river withdrawal included
       ENDIF

       !Groundwater - withdrawing from modelled aquifer or deep unlimited aquifer
       IF(local_gw_irr>0)THEN
         gwabstr = pirrgpar * local_gw_irr    !groundwater removal from aquifer (m3), regulated by pirrgpar
         IF(na>0)THEN
           IF(path(i)%rechargebasin.OR.path(i)%recievefraction>0.)THEN
!             aquiferv(ia) = aquiferv(ia) - gwabstr
             aquiferloss(ia) = aquiferloss(ia) + gwabstr    !save aquifer removal for after all subbasins
             wblirrflows(7) = gwabstr
           ELSE
             wblirrflows(1) = gwabstr
           ENDIF
         ELSE
           wblirrflows(1) = gwabstr
         ENDIF
       ENDIF

       !Calculate water actually reaching the fields
       local_add_irr = (ldabstr+lrabstr+gwabstr)*local_efficiency  !m3
       local_add_irr_factor = local_add_irr / accpw
       wblirrflows(5) = (ldabstr+lrabstr+gwabstr) - local_add_irr

       !Calculate losses in the irrigation network and the field      
       irrnetloss(i) = irrnetloss(i) + wblirrflows(5)     !m3/timestep

       !Calculate weighted mass and concentrations of withdrawals
       IF(local_add_irr >0.)THEN
         clocal_add_irr = ((clocal_dw_irr*local_dw_irr + clocal_dw2_irr*local_dw2_irr + cq*local_sw_irr + clocal_rv_irr*local_rv_irr)*pirrspar + clocal_gw_irr*local_gw_irr*pirrgpar)/local_add_irr ! Mass / volume
         CALL irrigation_abstraction_sink(clocal_add_irr,cirrsinkpar,local_add_irr/local_efficiency,irrsinkmass)         ! Modifies the concentrations of the withdrawal with a settlement basin 
         clocal_add_irr = clocal_add_irr / local_efficiency    ! Apply evapo-concentration of the concentrations (i.e. more concentrated due to water evaporating)
       ELSE
         clocal_add_irr = 0.
       ENDIF

       !Save irrigation amount and concentration to classes for next time step, spread in accordance with fieldneed.
       DO j = 1,nclass
         previous = nextirrigation(j,i)                                                                         ! Current depth in the relevant nextirrigation object
         nextirrigation(j,i) = nextirrigation(j,i) + local_add_irr_factor * fieldneed(j,i)                      ! Accumulated depth (mm) to be added, scaling fieldneed to local availability
         add = nextirrigation(j,i)-previous                                                                     ! Amount to be added in this iteration of the loop 
         IF(nextirrigation(j,i)>0) cnextirrigation(:,j,i) = (cnextirrigation(:,j,i)*previous + clocal_add_irr*add)/ nextirrigation(j,i)   ! New concentration (based on mass/volume but simplified to "depthconc"/depth: mm*mg/L / mm), Condition is there to avoid NaN.
       ENDDO

    ENDIF   !IF(accpw>0)

    !Remaining water need for this subbasin after irrigation with local water
    IF(accpw>0 .AND. twn_left>0 .AND. regirrwat_source(i)>0)THEN
      regionaldemand(i) = twn_left * regirrpar / irrigationsystem(irrindex(i))%reg_eff                          ! apply factor controlling the strength of the regional connection (regirrpar) and scale up to regional scale
      totaldemand_regsrc(regirrwat_source(i)) = totaldemand_regsrc(regirrwat_source(i)) + regionaldemand(i)      !m3
      IF(local_add_irr>0)THEN
        fieldneed(:,i) = fieldneed(:,i) * twn_left/twn * regirrpar                             ! calculate remaining fieldneed (need to be proportional to regionaldemand as this is used for setting nextirrigation), and apply factor controlling the strength of the regional connection
      ELSE
        fieldneed(:,i) = fieldneed(:,i) * regirrpar                             ! calculate remaining fieldneed, apply factor controlling the strength of the regional connection
      ENDIF
    ELSE
      fieldneed(:,i)=0.
    ENDIF

    !Irrigation with water from regional sources
    !-------------------------------------------
    IF(totaldemand_regsrc(i)>0 .AND. regionaldemand(i)==0.)THEN             !Jafet: Checks for regionaldemand(i) because otherwise you may take water from what was protected by pirrgpar and pirrspar above
       !No need to change for irrunlimited since totaldemand_regsrc(i) is then always 0 (unchanged since initiation) + we return before
       twn = totaldemand_regsrc(i)
       !Calculate regional water availability from irrigation dam (olake) (m3)
       IF(irrigationsystem(irrindex(i))%dam)THEN
          regional_dw_avail = org_lake *1.E-3 * classbasin(i,slc_olake)%part * basin(i)%area * safe_factor - local_dw_irr  !volume available after local withdrawals (before pirrspar)
          IF(regional_dw_avail<0) regional_dw_avail = 0.
          regional_dw_irr = MIN(twn,regional_dw_avail)
          regional_dw_avail = regional_dw_avail - regional_dw_irr
       ELSE
          regional_dw_avail = 0.
          regional_dw_irr = 0.
       ENDIF
       !Calculate regional water availability from river inflow (m3)
       IF(regional_dw_avail==0.)THEN
          regional_sw_avail = org_q * seconds_per_timestep * safe_factor - local_sw_irr    !m3
          IF(regional_sw_avail<0) regional_sw_avail = 0.
          regional_sw_irr = MIN(twn - regional_dw_irr,regional_sw_avail)
          regional_sw_avail = regional_sw_avail - regional_sw_irr
       ELSE
          regional_sw_avail = 0.  
          regional_sw_irr = 0.       
       ENDIF
       !Calculate regional water availability from river volume (m3)
       IF(regional_sw_avail==0.)THEN
          regional_rv_avail = org_rv * safe_factor - local_rv_irr      !m3
          IF(regional_rv_avail<0) regional_rv_avail = 0.
          regional_rv_irr = MIN(twn - regional_dw_irr - regional_sw_irr,regional_rv_avail)
       ELSE
          regional_rv_irr = 0.
       ENDIF

       !Calculate water actually taken for irrigation
       regional_dw_irr = regional_dw_irr * pirrspar
       regional_sw_irr = regional_sw_irr * pirrspar
       regional_rv_irr = regional_rv_irr * pirrspar
       regional_add_irr = regional_dw_irr + regional_sw_irr + regional_rv_irr  !m3
       regional_dw_fraction = regional_dw_irr/regional_add_irr
       regional_rw_fraction = (regional_sw_irr + regional_rv_irr)/regional_add_irr

       IF(regional_add_irr>0.)THEN

          !Remove water from the sources
          IF(regional_dw_irr>0)THEN
             lakew(2) = lakew(2) - regional_dw_irr / (classbasin(i,slc_olake)%part * basin(i)%area) * 1.E3   !mm
             rsabstr = regional_dw_irr
          ENDIF
          IF(regional_sw_irr>0)THEN
             q = q - regional_sw_irr/seconds_per_timestep     !m3/s
             rsabstr = rsabstr + regional_sw_irr
          ENDIF
          IF(regional_rv_irr>0.)THEN
             riverv = riverv - regional_rv_irr     !m3
             rsabstr = rsabstr + regional_rv_irr
          ENDIF

        !Calculate concentrations in the removed water
          cregional_add_irr = (cregional_dw_irr*regional_dw_irr + cq*regional_sw_irr + cregional_rv_irr*regional_rv_irr)/regional_add_irr
          CALL irrigation_abstraction_sink(cregional_add_irr,cirrsinkpar,regional_add_irr,irrsinkmass)

          !Calculate water & concentrations actually reaching the fields for each subbasin that has this subbasin as its regional source
          regional_add_irr_factor = regional_add_irr / twn
          sum_wtb = 0.
          DO k=1,nregirrwat_collectors(i)
             iadd = regirrwat_collectors(k,i)          ! the SUBBASIN index of the current "collector" subbasin.
             watertobasin = regionaldemand(iadd) * regional_add_irr_factor * irrigationsystem(irrindex(iadd))%reg_eff                      !what the local subbasin gets from the regional source (m3) after regional losses
             watertofield = watertobasin * irrigationsystem(irrindex(iadd))%local_eff                                                      !what the field of local subbasin gets (m3)
             cwatertofield = cregional_add_irr / (irrigationsystem(irrindex(iadd))%reg_eff * irrigationsystem(irrindex(iadd))%local_eff)   ! Concentration of watertofield (evaporation-> concentration both regionally and locally).
             wbrirrflows(1,iadd) = regional_dw_fraction * watertobasin
             wbrirrflows(2,iadd) = regional_rw_fraction * watertobasin
             
             !Calculate losses in the regional irrigation network, connected to source subbasin
             revap = regionaldemand(iadd) * regional_add_irr_factor - watertobasin
             wbrirrflows(3,i) = wbrirrflows(3,i) + regional_dw_fraction * revap
             wbrirrflows(4,i) = wbrirrflows(4,i) + regional_rw_fraction * revap

             !Calculate losses in the local irrigation network and the field
             wbrirrevap(iadd) = watertobasin - watertofield
             irrnetloss(iadd) = irrnetloss(iadd) + wbrirrevap(iadd)       !m3/timestep
             sum_wtb = sum_wtb + watertobasin

             !Save irrigation to classes for next time step (mm)
             DO j = 1,nclass
                previous = nextirrigation(j,iadd)                                                                               ! Current depth in the relevant nextirrigation object
                nextirrigation(j,iadd) = nextirrigation(j,iadd) + regional_add_irr_factor * fieldneed(j,iadd)                   ! Accumulated depth (mm) to be added, scaling fieldneed to local availability  ! Note fieldneed has been scaled already above with regirrpar
                add = nextirrigation(j,iadd)-previous                                                                           ! Amount to be added in this iteration of the loop             
                IF(nextirrigation(j,iadd)>0.) cnextirrigation(:,j,iadd) = (cnextirrigation(:,j,iadd)*previous + cwatertofield*add)/ nextirrigation(j,iadd)    ! New concentration (based on mass/volume but simplified to "depthconc"/depth: mm*mg/L / mm)
             ENDDO
          ENDDO

          !Calculate losses in the regional irrigation network
          irrnetloss(i) = irrnetloss(i) + regional_add_irr - sum_wtb       !m3/timestep
       ENDIF  ! regional_add_irr>0

    ENDIF

  END SUBROUTINE calculate_irrigation

  !>\brief Change concentration of the water abstracted for irrigation
  !>from all sources (rivers, lakes/dams, groundwater)
  !!This simulates a settlement basin at the point of the withdrawal
  !!Called once per timestep & subbasin & abstraction (in the
  !!CALCULATE_IRRIGATION subroutine)
  !!
  !>\b Reference ModelDescription Water management (Irrigation - Irrigation water withdrawal)  
  !--------------------------------------------------------------------------
  SUBROUTINE irrigation_abstraction_sink(conc,cirrsinkpar,abstr,sinkmass)

    USE MODVAR, ONLY :  numsubstances,  &
                       i_pp,i_on 

    !Argument declarations
    REAL, INTENT(INOUT) :: conc(numsubstances)      !<concentration of water abstracted from the source (typically in mg/L)
    REAL, INTENT(IN)    :: cirrsinkpar              !<irrsink parameter (-)
    REAL, INTENT(IN)    :: abstr                    !<amount of water abstracted from the source  (m3)
    REAL, INTENT(INOUT) :: sinkmass(numsubstances)  !<mass of the withdrawn substances in this sink (typically in kg/timestep)
    
    !Local variables
    REAL concin(numsubstances)                      !incoming concentrations    
    REAL csink(numsubstances)                       !concentrations of the sink (well, only the part that is calculated in this call)

    !Check if needed
    IF(abstr==0) RETURN

    !Initiations
    concin = conc    
    csink = 0.

    !Modify concentrations according to cirrsinkpar for PP and ON
    IF(i_pp>0) conc(i_pp) = conc(i_pp)*(1-cirrsinkpar)
    IF(i_on>0) conc(i_on) = conc(i_on)*(1-cirrsinkpar)

    !Calculate withdrawn mass, for summation elsewhere
    csink = concin - conc
    sinkmass = sinkmass + (csink * abstr / 1000)        ! mass of the the application (mg/L -> kg for e.g. nitrogen (mg/L=g/m3->g->kg)

  END SUBROUTINE irrigation_abstraction_sink

END MODULE IRRIGATION_MODULE
