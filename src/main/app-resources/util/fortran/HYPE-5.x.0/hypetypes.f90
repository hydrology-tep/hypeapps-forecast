!> \file hypetypes.f90
!> Contains module statetype_module.

!>Types for state variables for the HYPE model (HYdrological Predictions for the Environment)
!>
!>Five types are distinguished; four for different hydrological features in the 
!>landscape (snow and ice, soil, river and lake), and one for miscellaneous states.
!>Not all state variables are states of the hydrological model, but hold memory of 
!>older values or averages that is needed if the model is restarted.
!>Each statetype has seven subroutines that handle the state variables as a group; 
!>allocation, deallocation, initiation (to zero), initiation for submodel simulation, 
!>get size of states/array, write to array and read from array.
!!
MODULE STATETYPE_MODULE

!Copyright 2013-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!------------------------------------------------------------------------------
!Used by modules
!------------------------------------------------------------------------------
!"main"
!"optim"
!datamodule 
!modelmodule
!soilmodel_default
!glacier_soilmodel
!irrigation
!modelmodule
!soil_processes
!surfacewater_processes
!npc_soil_processes
!npc_surfacewater_processes
!regional_groundwater
!------------------------------------------------------------------------------

  IMPLICIT NONE
  PRIVATE

!(No) Parameter declarations

!Type declarations for model state variables
!> \brief Type for state variables related to snow and ice
  TYPE SNOWICESTATETYPE
    REAL,ALLOCATABLE :: snow(:,:)           !<snow on soil, snow water equivalent per soil-landuse combination (class,subbasin)
    REAL,ALLOCATABLE :: csnow(:,:,:)        !<concentration of snow on soil, snow water equivalent per soil-landuse combination (substance,class,subbasin)
    REAL,ALLOCATABLE :: snowage(:,:)        !<age of snow pack on soil per land-soil combination (class,subbasin)
    REAL,ALLOCATABLE :: snowdepth(:,:)      !<depth of snow pack on soil per land-soil combination (class,subbasin)
    REAL,ALLOCATABLE :: snowcov(:,:)        !<degree of snow cover on soil per land-soil combination (class,subbasin)
    REAL,ALLOCATABLE :: snowmax(:,:)        !<max snow pack on soil per land-soil combination (class,subbasin)
    REAL,ALLOCATABLE :: glacvol(:)          !<glacier ice volume, one glacier class per subbasin (m3)
    REAL,ALLOCATABLE :: riversnow(:,:)      !<snow on river, snow water equivalent per river type (type,subbasin)
    REAL,ALLOCATABLE :: riversnowage(:,:)   !<age of snow on river per river type (type,subbasin)
    REAL,ALLOCATABLE :: riversnowdepth(:,:) !<depth of snow pack on river per river type (type,subbasin)
    REAL,ALLOCATABLE :: riverice(:,:)       !<ice on river, total (cm) per river type (type,subbasin)
    REAL,ALLOCATABLE :: riverbice(:,:)      !<black ice on river (cm) per river type (type,subbasin)
    REAL,ALLOCATABLE :: rivericecov(:,:)    !<degree of ice cover on river (fraction) per river type (type,subbasin)
    REAL,ALLOCATABLE :: lakesnow(:,:)       !<snow on lake, snow water equivalent per lake type (type,subbasin) !!OR snow(slc_ilake,i)? 
    REAL,ALLOCATABLE :: lakesnowage(:,:)    !<age of snow on lake per lake type (type,subbasin) !!OR snow(slc_mriver,i) no used when a=0?
    REAL,ALLOCATABLE :: lakesnowdepth(:,:)  !<depth of snow pack on lake per lake type (type,subbasin)
    REAL,ALLOCATABLE :: lakeice(:,:)        !<ice on lake, total (cm) per lake type (type,subbasin)
    REAL,ALLOCATABLE :: lakebice(:,:)       !<black ice on lake (cm) per lake type (type,subbasin)
    REAL,ALLOCATABLE :: lakeicecov(:,:)     !<degree of ice cover on lake (fraction) per lake type (type,subbasin)
  END TYPE SNOWICESTATETYPE
!> \brief Type for state variables related to soil
  TYPE SOILSTATETYPE
    REAL,ALLOCATABLE :: water(:,:,:)   !<soil moisture per soil-landuse combination  (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: temp(:,:,:)    !<temperature of soil (degree Celcius) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: deeptemp(:,:)  !<temperature of deep soil (degree Celcius) (class,subbasin)
    REAL,ALLOCATABLE :: conc(:,:,:,:)  !<concentration of soil moisture per land-soil combination (substance,soillayer,class,subbasin)
    REAL,ALLOCATABLE :: humusN(:,:,:)  !<humusN pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: fastN(:,:,:)   !<fastN pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: partP(:,:,:)   !<partP pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: fastP(:,:,:)   !<fastP pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: humusP(:,:,:)  !<humusN pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: fastC(:,:,:)   !<fastC pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: humusC(:,:,:)  !<humusC pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin)
    REAL,ALLOCATABLE :: PPrelpool(:,:) !<pool of delayed PP in runoff (mg/m2) (class,subbasin)  
    REAL,ALLOCATABLE :: Srelpool(:,:)  !<pool of delayed SS in runoff (mg/m2) (class,subbasin)  
    REAL,ALLOCATABLE :: oldgrw(:,:)    !<class ground water table yesterday (m) (class,subbasin)
    REAL,ALLOCATABLE :: partT1(:,:,:)  !<adsorbed T1 pool in soil (kg/km2=mg/m2) (kg/km2/mm=mg/L) (soillayer,class,subbasin) 
    INTEGER,ALLOCATABLE :: icelens(:,:) !<presence of ice lens for infiltration into frozen soils (nclass,subbasin)
  END TYPE SOILSTATETYPE
!> \brief Type for state variables related to groundwater/aquifers
  TYPE AQUIFERSTATETYPE
    REAL,ALLOCATABLE :: water(:)           !<water volume of aquifer (m3 or m) (aquifer)
    REAL,ALLOCATABLE :: conc(:,:)          !<concentration of aquifer water (substance,aquifer)
    REAL,ALLOCATABLE :: lastrecharge(:)    !<recharge flow last timestep (m3 or m) (aquifer)
    REAL,ALLOCATABLE :: clastrecharge(:,:) !<concentration of recharge flow last time step (substance,aquifer)
    REAL,ALLOCATABLE :: nextoutflow(:)     !<outflow of aquifer to be applied next time step (m3) (aquifer)
    REAL,ALLOCATABLE :: cnextoutflow(:,:)  !<concentration of outflow of aquifer (substance,aquifer)
  END TYPE AQUIFERSTATETYPE
!> \brief Type for state variables related to rivers
  TYPE RIVERSTATETYPE
    REAL,ALLOCATABLE :: water(:,:)       !<water volume of river (rivertype,subbasin) (m3)
    REAL,ALLOCATABLE :: temp(:,:)        !<river temperature (rivertype,subbasin)
    REAL,ALLOCATABLE :: conc(:,:,:)      !<concentration in river box (substance,rivertype,subbasin)
    REAL,ALLOCATABLE :: TPmean(:,:)      !<365-day mean Total Phosphorus concentration of river (rivertype,subbasin)
    REAL,ALLOCATABLE :: temp10(:,:)      !<10-day mean river temperature (rivertype,subbasin)
    REAL,ALLOCATABLE :: temp20(:,:)      !<20-day mean river temperature (rivertype,subbasin)
    REAL,ALLOCATABLE :: Psed(:,:)        !<phosphorus sediment in river (kg) (rivertype,subbasin)
    REAL,ALLOCATABLE :: qqueue(:,:,:)    !<q of water translated in river (m3) (0:maxlag,rivertype,subbasin)
    REAL,ALLOCATABLE :: cqueue(:,:,:,:)  !<c of water translated in river (mg/l) (substance,0:maxlag,rivertype,subbasin)
    REAL,ALLOCATABLE :: cwetland(:,:,:)  !<concentration of river wetland volume (mg/l) (substance,rivertype,subbasin)
    REAL,ALLOCATABLE :: Qdayacc(:,:,:)   !<Q on time step<1 day, used in calculation of daily mean (m3/s) (timestep,rivertype,subbasin)
    REAL,ALLOCATABLE :: Q365(:,:,:)      !<Q river last 365 days (m3/s) (dayno,rivertype,subbasin)
    REAL,ALLOCATABLE :: Qmean(:,:)       !<average discharge in river (365-days-MA) (m3/s) (rivertype,subbasin)
    REAL,ALLOCATABLE :: T1sed(:,:)       !<T1 in sediment in river (kg) (rivertype,subbasin)
    REAL,ALLOCATABLE :: Ssed(:,:)        !<sediment in river (kg) (rivertype,subbasin)
  END TYPE RIVERSTATETYPE
!> \brief Type for state variables related to lake
  TYPE LAKESTATETYPE
    REAL,ALLOCATABLE :: water(:,:)      !<water stage of lake (mm) (laketype,subbasin)
    REAL,ALLOCATABLE :: temp(:,:)       !<lake temperature (laketype,subbasin)
    REAL,ALLOCATABLE :: conc(:,:,:)     !<concentration in lakes (mg/L etc) (substances,laketype,subbasin)
    REAL,ALLOCATABLE :: TPmean(:,:)     !<365-day mean Total Phosphorus concentration of lake (laketype,subbasin)
    REAL,ALLOCATABLE :: slowwater(:,:)  !<volume or level of long turnover lake part (mm) (laketype,subbasin)
    REAL,ALLOCATABLE :: concslow(:,:,:) !<concentrations of long turnover lake part (mg/l) (substance,laketype,subbasin)
    REAL,ALLOCATABLE :: temp10(:,:)     !<10-day mean lake temperature (laketype,subbasin)
    REAL,ALLOCATABLE :: temp20(:,:)     !<20-day mean lake temperature (laketype,subbasin)
    REAL,ALLOCATABLE :: uppertemp(:,:)  !<upper layer (epilimnion) temperature (laketype,subbasin)
    REAL,ALLOCATABLE :: lowertemp(:,:)  !<lower layer (hypolimnion) temperature (laketype,subbasin)
  END TYPE LAKESTATETYPE
!> \brief Type for miscellaneous state variables
  TYPE MISCSTATETYPE
    REAL,ALLOCATABLE :: temp5(:)      !<5-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: temp30(:)     !<30-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: temp10(:)     !<10-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: temp20(:)     !<20-day mean air temperature (subbasin)
    REAL,ALLOCATABLE :: gdd(:,:,:)    !<degree days for growing season start (croporder,class,subbasin)
    INTEGER, ALLOCATABLE :: gsbegin(:,:,:)      !<growing seasons start pseudo dayno (croporder,class,subbasin)
    REAL,ALLOCATABLE :: nextirrigation(:,:)     !<irrigation water to be applied next timestep (mm) (class,subbasin)
    REAL,ALLOCATABLE :: cnextirrigation(:,:,:)  !<concentration of the irrigation water (substance,class,subbasin)
    REAL,ALLOCATABLE :: updatestationsarcorr(:) !<AR-error of each subbasin (update function qar) (subbasin)
    REAL,ALLOCATABLE :: floodwater(:,:)         !<volume of water flooded (m3) (mainriver/olake,subbasin)
    REAL,ALLOCATABLE :: cfloodwater(:,:,:)      !<concentration of water flooded (m3) (substances,mainriver/olake,subbasin)
    REAL,ALLOCATABLE :: partT1sf(:,:) !<T1 on soil surface (kg/km2=mg/m2) (kg/km2/mm=mg/L) (class,subbasin) 
  END TYPE MISCSTATETYPE

!Type declarations for model state information variables
!> \brief Type for pointers to state variables 
    TYPE POINTERTYPE
      REAL,POINTER :: d1(:) => NULL()
      REAL,POINTER :: d2(:,:) => NULL()
      INTEGER,POINTER :: d2i(:,:) => NULL()
      REAL,POINTER :: d3(:,:,:) => NULL()
      INTEGER,POINTER :: d3i(:,:,:) => NULL()
      REAL,POINTER :: d4(:,:,:,:) => NULL()
    ENDTYPE
!> \brief Type for model state information 
    TYPE STATEINFOTYPE
      LOGICAL :: allok             !state variable allocation status
      INTEGER :: ndim              !state variable number of dimensions
      INTEGER :: dims(4)           !max dimension for state variable is 4
      CHARACTER(LEN=20) :: svname  !state variable name
      CHARACTER(LEN=20) :: stname  !state type name
      TYPE(POINTERTYPE) :: svpoint !pointer to state variable
    ENDTYPE

  PUBLIC SNOWICESTATETYPE
  PUBLIC SOILSTATETYPE
  PUBLIC AQUIFERSTATETYPE
  PUBLIC RIVERSTATETYPE
  PUBLIC LAKESTATETYPE
  PUBLIC MISCSTATETYPE
  PUBLIC STATEINFOTYPE
  PUBLIC set_stateinfo
  PUBLIC allocate_model_states                !Eventuellt ska dessa flyttas ut i HYSS, i en data-module med starttillstånd-hantering
  PUBLIC deallocate_model_states
  PUBLIC initiate_state_zero
  PUBLIC initiate_frozenstate_submodel
  PUBLIC initiate_soilstate_submodel
  PUBLIC initiate_aquiferstate_submodel
  PUBLIC initiate_riverstate_submodel
  PUBLIC initiate_lakestate_submodel
  PUBLIC initiate_miscstate_submodel
  PUBLIC get_frozenstate_variables_arraysize
  PUBLIC set_frozenstate_variables_to_array
  PUBLIC set_frozenstate_variables_from_array
  PUBLIC get_soilstate_variables_arraysize
  PUBLIC set_soilstate_variables_to_array
  PUBLIC set_soilstate_variables_from_array
  PUBLIC get_aquiferstate_variables_arraysize
  PUBLIC set_aquiferstate_variables_to_array
  PUBLIC set_aquiferstate_variables_from_array
  PUBLIC get_riverstate_variables_arraysize
  PUBLIC set_riverstate_variables_to_array
  PUBLIC set_riverstate_variables_from_array
  PUBLIC get_lakestate_variables_arraysize
  PUBLIC set_lakestate_variables_to_array
  PUBLIC set_lakestate_variables_from_array
  PUBLIC get_miscstate_variables_arraysize
  PUBLIC set_miscstate_variables_to_array
  PUBLIC set_miscstate_variables_from_array

CONTAINS

  !>Allocate state information variable and assign the information
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE set_stateinfo(n,ns,nc,na,nsl,nr,nl,ml,mt,stateinfo,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    INTEGER,INTENT(IN) :: ml      !<Maximum river lag steps
    INTEGER,INTENT(IN) :: mt      !<Number of timesteps per day
    TYPE(stateinfotype),ALLOCATABLE,INTENT(INOUT) :: stateinfo(:)
    TYPE(snowicestatetype),TARGET,INTENT(IN) :: frozenstate  !<Snow and ice states
    TYPE(soilstatetype),TARGET,INTENT(IN)    :: soilstate    !<Soil states
    TYPE(aquiferstatetype),TARGET,INTENT(IN) :: aquiferstate !<Aquifer states
    TYPE(riverstatetype),TARGET,INTENT(IN)   :: riverstate   !<River states
    TYPE(lakestatetype),TARGET,INTENT(IN)    :: lakestate    !<Lake states
    TYPE(miscstatetype),TARGET,INTENT(IN)    :: miscstate    !<Misc states
      
    !Local variables 
    INTEGER i,iadd
    INTEGER nstates    !number of states in statetypes
    
    !Local parameter (if these are changed, corresponding changes must be done to assimInfo%assim_categories!)
    INTEGER, PARAMETER :: nfrozenstates = 19
    INTEGER, PARAMETER :: nsoilstates = 16
    INTEGER, PARAMETER :: naquiferstates = 6
    INTEGER, PARAMETER :: nriverstates = 15
    INTEGER, PARAMETER :: nlakestates = 10
    INTEGER, PARAMETER :: nmiscstates = 12
    
    !Allocate and initialise state information variable
    nstates = nfrozenstates+nsoilstates+naquiferstates+nriverstates+nlakestates+nmiscstates
    ALLOCATE(stateinfo(nstates))
    DO i = 1,4    !4 is maximum number of dimensions among states at the moment
      stateinfo(:)%dims(i) = 0
    ENDDO
    stateinfo(:)%allok = .FALSE.
    
    !Set state information for snowicestatetype states (1:nfrozenstates)
    
    stateinfo(1:nfrozenstates)%stname = 'frozenstate'
    stateinfo(1)%svname = 'snow'
    stateinfo(2)%svname = 'csnow'
    stateinfo(3)%svname = 'snowage'
    stateinfo(4)%svname = 'snowdepth'
    stateinfo(5)%svname = 'snowcov'
    stateinfo(6)%svname = 'snowmax'
    stateinfo(7)%svname = 'glacvol'
    stateinfo(8)%svname = 'riversnow'
    stateinfo(9)%svname = 'riversnowage'
    stateinfo(10)%svname = 'riversnowdepth'
    stateinfo(11)%svname = 'riverice'
    stateinfo(12)%svname = 'riverbice'
    stateinfo(13)%svname = 'rivericecov'
    stateinfo(14)%svname = 'lakesnow'
    stateinfo(15)%svname = 'lakesnowage'
    stateinfo(16)%svname = 'lakesnowdepth'
    stateinfo(17)%svname = 'lakeice'
    stateinfo(18)%svname = 'lakebice'
    stateinfo(19)%svname = 'lakeicecov'
    stateinfo(1)%ndim = 2
    stateinfo(2)%ndim = 3
    stateinfo(3:6)%ndim = 2
    stateinfo(7)%ndim = 1
    stateinfo(8:19)%ndim = 2
    stateinfo(1)%dims(1:2) = (/nc,n/)
    stateinfo(2)%dims(1:3) = (/ns,nc,n/)
    DO i = 3,6
      stateinfo(i)%dims(1:2) = (/nc,n/)
    ENDDO
    stateinfo(7)%dims(1)   = n
    DO i = 8,13
      stateinfo(i)%dims(1:2) = (/nr,n/)
    ENDDO
    DO i = 14,19
      stateinfo(i)%dims(1:2) = (/nl,n/)
    ENDDO
    IF(ALLOCATED(frozenstate%snow))THEN
      stateinfo(1)%allok = .TRUE.
      ALLOCATE(stateinfo(1)%svpoint%d2(stateinfo(1)%dims(1),stateinfo(1)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%csnow))THEN
      stateinfo(2)%allok = .TRUE.
      ALLOCATE(stateinfo(2)%svpoint%d3(stateinfo(2)%dims(1),stateinfo(2)%dims(2),stateinfo(2)%dims(3)))
    ENDIF
    IF(ALLOCATED(frozenstate%snowage))THEN
      stateinfo(3)%allok = .TRUE.
      ALLOCATE(stateinfo(3)%svpoint%d2(stateinfo(3)%dims(1),stateinfo(3)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%snowdepth))THEN
      stateinfo(4)%allok = .TRUE.
      ALLOCATE(stateinfo(4)%svpoint%d2(stateinfo(4)%dims(1),stateinfo(4)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%snowcov))THEN
      stateinfo(5)%allok = .TRUE.
      ALLOCATE(stateinfo(5)%svpoint%d2(stateinfo(5)%dims(1),stateinfo(5)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%snowmax))THEN
      stateinfo(6)%allok = .TRUE.
      ALLOCATE(stateinfo(6)%svpoint%d2(stateinfo(6)%dims(1),stateinfo(6)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%glacvol))THEN
      stateinfo(7)%allok = .TRUE.
      ALLOCATE(stateinfo(7)%svpoint%d1(stateinfo(7)%dims(1)))
    ENDIF
    IF(ALLOCATED(frozenstate%riversnow))THEN
      stateinfo(8)%allok = .TRUE.
      ALLOCATE(stateinfo(8)%svpoint%d2(stateinfo(8)%dims(1),stateinfo(8)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%riversnowage))THEN
      stateinfo(9)%allok = .TRUE.
      ALLOCATE(stateinfo(9)%svpoint%d2(stateinfo(9)%dims(1),stateinfo(9)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%riversnowdepth))THEN
      stateinfo(10)%allok = .TRUE.
      ALLOCATE(stateinfo(10)%svpoint%d2(stateinfo(10)%dims(1),stateinfo(10)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%riverice))THEN
      stateinfo(11)%allok = .TRUE.
      ALLOCATE(stateinfo(11)%svpoint%d2(stateinfo(11)%dims(1),stateinfo(11)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%riverbice))THEN
      stateinfo(12)%allok = .TRUE.
      ALLOCATE(stateinfo(12)%svpoint%d2(stateinfo(12)%dims(1),stateinfo(12)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%rivericecov))THEN
      stateinfo(13)%allok = .TRUE.
      ALLOCATE(stateinfo(13)%svpoint%d2(stateinfo(13)%dims(1),stateinfo(13)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%lakesnow))THEN
      stateinfo(14)%allok = .TRUE.
      ALLOCATE(stateinfo(14)%svpoint%d2(stateinfo(14)%dims(1),stateinfo(14)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%lakesnowage))THEN
      stateinfo(15)%allok = .TRUE.
      ALLOCATE(stateinfo(15)%svpoint%d2(stateinfo(15)%dims(1),stateinfo(15)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%lakesnowdepth))THEN
      stateinfo(16)%allok = .TRUE.
      ALLOCATE(stateinfo(16)%svpoint%d2(stateinfo(16)%dims(1),stateinfo(16)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%lakeice))THEN
      stateinfo(17)%allok = .TRUE.
      ALLOCATE(stateinfo(17)%svpoint%d2(stateinfo(17)%dims(1),stateinfo(17)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%lakebice))THEN
      stateinfo(18)%allok = .TRUE.
      ALLOCATE(stateinfo(18)%svpoint%d2(stateinfo(18)%dims(1),stateinfo(18)%dims(2)))
    ENDIF
    IF(ALLOCATED(frozenstate%lakeicecov))THEN
      stateinfo(19)%allok = .TRUE.
      ALLOCATE(stateinfo(19)%svpoint%d2(stateinfo(19)%dims(1),stateinfo(19)%dims(2)))
    ENDIF
    stateinfo(1)%svpoint%d2  => frozenstate%snow
    stateinfo(2)%svpoint%d3  => frozenstate%csnow
    stateinfo(3)%svpoint%d2  => frozenstate%snowage
    stateinfo(4)%svpoint%d2  => frozenstate%snowdepth
    stateinfo(5)%svpoint%d2  => frozenstate%snowcov
    stateinfo(6)%svpoint%d2  => frozenstate%snowmax
    stateinfo(7)%svpoint%d1  => frozenstate%glacvol
    stateinfo(8)%svpoint%d2  => frozenstate%riversnow
    stateinfo(9)%svpoint%d2  => frozenstate%riversnowage
    stateinfo(10)%svpoint%d2 => frozenstate%riversnowdepth
    stateinfo(11)%svpoint%d2 => frozenstate%riverice
    stateinfo(12)%svpoint%d2 => frozenstate%riverbice
    stateinfo(13)%svpoint%d2 => frozenstate%rivericecov
    stateinfo(14)%svpoint%d2 => frozenstate%lakesnow
    stateinfo(15)%svpoint%d2 => frozenstate%lakesnowage
    stateinfo(16)%svpoint%d2 => frozenstate%lakesnowdepth
    stateinfo(17)%svpoint%d2 => frozenstate%lakeice
    stateinfo(18)%svpoint%d2 => frozenstate%lakebice
    stateinfo(19)%svpoint%d2 => frozenstate%lakeicecov

    !Set state information for soilstatetype states
    iadd=nfrozenstates
    stateinfo(iadd+1:iadd+nsoilstates)%stname = 'soilstate'
    stateinfo(iadd+1)%svname = 'water'
    stateinfo(iadd+2)%svname = 'temp'
    stateinfo(iadd+3)%svname = 'deeptemp'
    stateinfo(iadd+4)%svname = 'conc'
    stateinfo(iadd+5)%svname = 'humusN'
    stateinfo(iadd+6)%svname = 'fastN'
    stateinfo(iadd+7)%svname = 'partP'
    stateinfo(iadd+8)%svname = 'fastP'
    stateinfo(iadd+9)%svname = 'humusP'
    stateinfo(iadd+10)%svname = 'fastC'
    stateinfo(iadd+11)%svname = 'humusC'
    stateinfo(iadd+12)%svname = 'PPrelpool'
    stateinfo(iadd+13)%svname = 'Srelpool'
    stateinfo(iadd+14)%svname = 'oldgrw'
    stateinfo(iadd+15)%svname = 'icelens'
    stateinfo(iadd+16)%svname = 'partT1'
    stateinfo(iadd+1:iadd+2)%ndim = 3
    stateinfo(iadd+3)%ndim = 2
    stateinfo(iadd+4)%ndim = 4
    stateinfo(iadd+5:iadd+11)%ndim = 3
    stateinfo(iadd+12:iadd+14)%ndim = 2
    stateinfo(iadd+1)%dims(1:3) = (/nsl,nc,n/)
    stateinfo(iadd+2)%dims(1:3) = (/nsl,nc,n/)
    stateinfo(iadd+3)%dims(1:2) = (/nc,n/)
    stateinfo(iadd+4)%dims(1:4) = (/ns,nsl,nc,n/)
    DO i = iadd+5,iadd+11
      stateinfo(i)%dims(1:3) = (/nsl,nc,n/)
    ENDDO
    DO i = iadd+12,iadd+15
      stateinfo(i)%dims(1:2) = (/nc,n/)
    ENDDO
    stateinfo(iadd+16)%dims(1:3) = (/nsl,nc,n/)
    IF(ALLOCATED(soilstate%water))THEN
      stateinfo(iadd+1)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+1)%svpoint%d3(stateinfo(iadd+1)%dims(1),stateinfo(iadd+1)%dims(2),stateinfo(iadd+1)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%temp))THEN
      stateinfo(iadd+2)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+2)%svpoint%d3(stateinfo(iadd+2)%dims(1),stateinfo(iadd+2)%dims(2),stateinfo(iadd+2)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%deeptemp))THEN
      stateinfo(iadd+3)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+3)%svpoint%d2(stateinfo(iadd+3)%dims(1),stateinfo(iadd+3)%dims(2)))
    ENDIF
    IF(ALLOCATED(soilstate%conc))THEN
      stateinfo(iadd+4)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+4)%svpoint%d4(stateinfo(iadd+4)%dims(1),stateinfo(iadd+4)%dims(2),stateinfo(iadd+4)%dims(3),stateinfo(iadd+4)%dims(4)))
    ENDIF
    IF(ALLOCATED(soilstate%humusN))THEN
      stateinfo(iadd+5)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+5)%svpoint%d3(stateinfo(iadd+5)%dims(1),stateinfo(iadd+5)%dims(2),stateinfo(iadd+5)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%fastN))THEN
      stateinfo(iadd+6)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+6)%svpoint%d3(stateinfo(iadd+6)%dims(1),stateinfo(iadd+6)%dims(2),stateinfo(iadd+6)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%partP))THEN
      stateinfo(iadd+7)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+7)%svpoint%d3(stateinfo(iadd+7)%dims(1),stateinfo(iadd+7)%dims(2),stateinfo(iadd+7)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%fastP))THEN
      stateinfo(iadd+8)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+8)%svpoint%d3(stateinfo(iadd+8)%dims(1),stateinfo(iadd+8)%dims(2),stateinfo(iadd+8)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%humusP))THEN
      stateinfo(iadd+9)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+9)%svpoint%d3(stateinfo(iadd+9)%dims(1),stateinfo(iadd+9)%dims(2),stateinfo(iadd+9)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%fastC))THEN
      stateinfo(iadd+10)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+10)%svpoint%d3(stateinfo(iadd+10)%dims(1),stateinfo(iadd+10)%dims(2),stateinfo(iadd+10)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%humusC))THEN
      stateinfo(iadd+11)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+11)%svpoint%d3(stateinfo(iadd+11)%dims(1),stateinfo(iadd+11)%dims(2),stateinfo(iadd+11)%dims(3)))
    ENDIF
    IF(ALLOCATED(soilstate%PPrelpool))THEN
      stateinfo(iadd+12)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+12)%svpoint%d2(stateinfo(iadd+12)%dims(1),stateinfo(iadd+12)%dims(2)))
    ENDIF
    IF(ALLOCATED(soilstate%Srelpool))THEN
      stateinfo(iadd+13)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+13)%svpoint%d2(stateinfo(iadd+13)%dims(1),stateinfo(iadd+13)%dims(2)))
    ENDIF
    IF(ALLOCATED(soilstate%oldgrw))THEN
      stateinfo(iadd+14)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+14)%svpoint%d2(stateinfo(iadd+14)%dims(1),stateinfo(iadd+14)%dims(2)))
    ENDIF
    IF(ALLOCATED(soilstate%icelens))THEN
      stateinfo(iadd+15)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+15)%svpoint%d2i(stateinfo(iadd+15)%dims(1),stateinfo(iadd+15)%dims(2)))
    ENDIF
    IF(ALLOCATED(soilstate%partT1))THEN
      stateinfo(iadd+16)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+16)%svpoint%d3(stateinfo(iadd+16)%dims(1),stateinfo(iadd+16)%dims(2),stateinfo(iadd+16)%dims(3)))
    ENDIF
    stateinfo(iadd+1)%svpoint%d3  => soilstate%water
    stateinfo(iadd+2)%svpoint%d3  => soilstate%temp
    stateinfo(iadd+3)%svpoint%d2  => soilstate%deeptemp
    stateinfo(iadd+4)%svpoint%d4  => soilstate%conc
    stateinfo(iadd+5)%svpoint%d3  => soilstate%humusN
    stateinfo(iadd+6)%svpoint%d3  => soilstate%fastN
    stateinfo(iadd+7)%svpoint%d3  => soilstate%partP
    stateinfo(iadd+8)%svpoint%d3  => soilstate%fastP
    stateinfo(iadd+9)%svpoint%d3  => soilstate%humusP
    stateinfo(iadd+10)%svpoint%d3 => soilstate%fastC
    stateinfo(iadd+11)%svpoint%d3 => soilstate%humusC
    stateinfo(iadd+12)%svpoint%d2 => soilstate%PPrelpool
    stateinfo(iadd+13)%svpoint%d2 => soilstate%Srelpool
    stateinfo(iadd+14)%svpoint%d2 => soilstate%oldgrw
    stateinfo(iadd+15)%svpoint%d2i => soilstate%icelens
    stateinfo(iadd+16)%svpoint%d3  => soilstate%partT1

    !Set state information for aquiferstatetype states
    iadd=nfrozenstates+nsoilstates
    stateinfo(iadd+1:iadd+naquiferstates)%stname = 'aquiferstate'
    stateinfo(iadd+1)%svname = 'water'
    stateinfo(iadd+2)%svname = 'conc'
    stateinfo(iadd+3)%svname = 'lastrecharge'
    stateinfo(iadd+4)%svname = 'clastrecharge'
    stateinfo(iadd+5)%svname = 'nextoutflow'
    stateinfo(iadd+6)%svname = 'cnextoutflow'
    stateinfo(iadd+1)%ndim = 1
    stateinfo(iadd+2)%ndim = 2
    stateinfo(iadd+3)%ndim = 1
    stateinfo(iadd+4)%ndim = 2
    stateinfo(iadd+5)%ndim = 1
    stateinfo(iadd+6)%ndim = 2
    stateinfo(iadd+1)%dims(1) = na
    stateinfo(iadd+2)%dims(1:2) = (/ns,na/)
    stateinfo(iadd+3)%dims(1) = na
    stateinfo(iadd+4)%dims(1:2) = (/ns,na/)
    stateinfo(iadd+5)%dims(1) = na
    stateinfo(iadd+6)%dims(1:2) = (/ns,na/)
    IF(ALLOCATED(aquiferstate%water))THEN
      stateinfo(iadd+1)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+1)%svpoint%d1(stateinfo(iadd+1)%dims(1)))
    ENDIF
    IF(ALLOCATED(aquiferstate%conc))THEN
      stateinfo(iadd+2)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+2)%svpoint%d2(stateinfo(iadd+2)%dims(1),stateinfo(iadd+2)%dims(2)))
    ENDIF
    IF(ALLOCATED(aquiferstate%lastrecharge))THEN
      stateinfo(iadd+3)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+3)%svpoint%d1(stateinfo(iadd+3)%dims(1)))
    ENDIF
    IF(ALLOCATED(aquiferstate%clastrecharge))THEN
      stateinfo(iadd+4)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+4)%svpoint%d2(stateinfo(iadd+4)%dims(1),stateinfo(iadd+4)%dims(2)))
    ENDIF
    IF(ALLOCATED(aquiferstate%nextoutflow))THEN
      stateinfo(iadd+5)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+5)%svpoint%d1(stateinfo(iadd+5)%dims(1)))
    ENDIF
    IF(ALLOCATED(aquiferstate%cnextoutflow))THEN
      stateinfo(iadd+6)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+6)%svpoint%d2(stateinfo(iadd+6)%dims(1),stateinfo(iadd+6)%dims(2)))
    ENDIF
    stateinfo(iadd+1)%svpoint%d1  => aquiferstate%water
    stateinfo(iadd+2)%svpoint%d2  => aquiferstate%conc
    stateinfo(iadd+3)%svpoint%d1  => aquiferstate%lastrecharge
    stateinfo(iadd+4)%svpoint%d2  => aquiferstate%clastrecharge
    stateinfo(iadd+5)%svpoint%d1  => aquiferstate%nextoutflow
    stateinfo(iadd+6)%svpoint%d2  => aquiferstate%cnextoutflow

    !Set state information for riverstatetype states
    iadd=nfrozenstates+nsoilstates+naquiferstates
    stateinfo(iadd+1:iadd+nriverstates)%stname = 'riverstate'
    stateinfo(iadd+1)%svname = 'water'
    stateinfo(iadd+2)%svname = 'temp'
    stateinfo(iadd+3)%svname = 'conc'
    stateinfo(iadd+4)%svname = 'TPmean'
    stateinfo(iadd+5)%svname = 'temp10'
    stateinfo(iadd+6)%svname = 'temp20'
    stateinfo(iadd+7)%svname = 'Psed'
    stateinfo(iadd+8)%svname = 'qqueue'
    stateinfo(iadd+9)%svname = 'cqueue'
    stateinfo(iadd+10)%svname = 'cwetland'
    stateinfo(iadd+11)%svname = 'Qdayacc'
    stateinfo(iadd+12)%svname = 'Q365'
    stateinfo(iadd+13)%svname = 'Qmean'
    stateinfo(iadd+14)%svname = 'T1sed'
    stateinfo(iadd+15)%svname = 'sed'
    stateinfo(iadd+1:iadd+2)%ndim = 2
    stateinfo(iadd+3)%ndim = 3
    stateinfo(iadd+4:iadd+7)%ndim = 2
    stateinfo(iadd+8)%ndim = 3
    stateinfo(iadd+9)%ndim = 4
    stateinfo(iadd+10:iadd+12)%ndim = 3
    stateinfo(iadd+13:iadd+15)%ndim = 2
    stateinfo(iadd+1)%dims(1:2) = (/nr,n/)
    stateinfo(iadd+2)%dims(1:2) = (/nr,n/)
    stateinfo(iadd+3)%dims(1:3) = (/ns,nr,n/)
    DO i = iadd+4,iadd+7
      stateinfo(i)%dims(1:2) = (/nr,n/)
    ENDDO
    stateinfo(iadd+8)%dims(1:3) = (/ml+1,nr,n/)   !(0:ml,nr,n)
    stateinfo(iadd+9)%dims(1:4) = (/ns,ml+1,nr,n/)   !(ns,0:ml,nr,n)
    stateinfo(iadd+10)%dims(1:3) = (/ns,nr,n/)
    stateinfo(iadd+11)%dims(1:3) = (/mt,nr,n/)
    stateinfo(iadd+12)%dims(1:3) = (/366,nr,n/)
    DO i = iadd+13,iadd+15
      stateinfo(i)%dims(1:2) = (/nr,n/)
    ENDDO
    IF(ALLOCATED(riverstate%water))THEN
      stateinfo(iadd+1)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+1)%svpoint%d2(stateinfo(iadd+1)%dims(1),stateinfo(iadd+1)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%temp))THEN
      stateinfo(iadd+2)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+2)%svpoint%d2(stateinfo(iadd+2)%dims(1),stateinfo(iadd+2)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%conc))THEN
      stateinfo(iadd+3)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+3)%svpoint%d3(stateinfo(iadd+3)%dims(1),stateinfo(iadd+3)%dims(2),stateinfo(iadd+3)%dims(3)))
    ENDIF
    IF(ALLOCATED(riverstate%TPmean))THEN
      stateinfo(iadd+4)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+4)%svpoint%d2(stateinfo(iadd+4)%dims(1),stateinfo(iadd+4)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%temp10))THEN
      stateinfo(iadd+5)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+5)%svpoint%d2(stateinfo(iadd+5)%dims(1),stateinfo(iadd+5)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%temp20))THEN
      stateinfo(iadd+6)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+6)%svpoint%d2(stateinfo(iadd+6)%dims(1),stateinfo(iadd+6)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%Psed))THEN
      stateinfo(iadd+7)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+7)%svpoint%d2(stateinfo(iadd+7)%dims(1),stateinfo(iadd+7)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%qqueue))THEN    !0:ml+1-1
      stateinfo(iadd+8)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+8)%svpoint%d3(0:stateinfo(iadd+8)%dims(1)-1,stateinfo(iadd+8)%dims(2),stateinfo(iadd+8)%dims(3)))
    ENDIF
    IF(ALLOCATED(riverstate%cqueue))THEN    !0:ml+1-1
      stateinfo(iadd+9)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+9)%svpoint%d4(stateinfo(iadd+9)%dims(1),0:stateinfo(iadd+9)%dims(2)-1,stateinfo(iadd+9)%dims(3),stateinfo(iadd+9)%dims(4)))
    ENDIF
    IF(ALLOCATED(riverstate%cwetland))THEN
      stateinfo(iadd+10)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+10)%svpoint%d3(stateinfo(iadd+10)%dims(1),stateinfo(iadd+10)%dims(2),stateinfo(iadd+10)%dims(3)))
    ENDIF
    IF(ALLOCATED(riverstate%Qdayacc))THEN
      stateinfo(iadd+11)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+11)%svpoint%d3(stateinfo(iadd+11)%dims(1),stateinfo(iadd+11)%dims(2),stateinfo(iadd+11)%dims(3)))
    ENDIF
    IF(ALLOCATED(riverstate%Q365))THEN
      stateinfo(iadd+12)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+12)%svpoint%d3(stateinfo(iadd+12)%dims(1),stateinfo(iadd+12)%dims(2),stateinfo(iadd+12)%dims(3)))
    ENDIF
    IF(ALLOCATED(riverstate%Qmean))THEN
      stateinfo(iadd+13)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+13)%svpoint%d2(stateinfo(iadd+13)%dims(1),stateinfo(iadd+13)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%T1sed))THEN
      stateinfo(iadd+14)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+14)%svpoint%d2(stateinfo(iadd+14)%dims(1),stateinfo(iadd+14)%dims(2)))
    ENDIF
    IF(ALLOCATED(riverstate%Ssed))THEN
      stateinfo(iadd+15)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+15)%svpoint%d2(stateinfo(iadd+15)%dims(1),stateinfo(iadd+15)%dims(2)))
    ENDIF
    stateinfo(iadd+1)%svpoint%d2  => riverstate%water
    stateinfo(iadd+2)%svpoint%d2  => riverstate%temp
    stateinfo(iadd+3)%svpoint%d3  => riverstate%conc
    stateinfo(iadd+4)%svpoint%d2  => riverstate%TPmean
    stateinfo(iadd+5)%svpoint%d2  => riverstate%temp10
    stateinfo(iadd+6)%svpoint%d2  => riverstate%temp20
    stateinfo(iadd+7)%svpoint%d2  => riverstate%Psed
    stateinfo(iadd+8)%svpoint%d3  => riverstate%qqueue
    stateinfo(iadd+9)%svpoint%d4  => riverstate%cqueue
    stateinfo(iadd+10)%svpoint%d3 => riverstate%cwetland
    stateinfo(iadd+11)%svpoint%d3 => riverstate%Qdayacc
    stateinfo(iadd+12)%svpoint%d3 => riverstate%Q365
    stateinfo(iadd+13)%svpoint%d2 => riverstate%Qmean
    stateinfo(iadd+14)%svpoint%d2 => riverstate%T1sed
    stateinfo(iadd+15)%svpoint%d2 => riverstate%Ssed

    !Set state information for lakestatetype states
    iadd=nfrozenstates+nsoilstates+naquiferstates+nriverstates
    stateinfo(iadd+1:iadd+nlakestates)%stname = 'lakestate'
    stateinfo(iadd+1)%svname = 'water'
    stateinfo(iadd+2)%svname = 'temp'
    stateinfo(iadd+3)%svname = 'conc'
    stateinfo(iadd+4)%svname = 'TPmean'
    stateinfo(iadd+5)%svname = 'slowwater'
    stateinfo(iadd+6)%svname = 'concslow'
    stateinfo(iadd+7)%svname = 'temp10'
    stateinfo(iadd+8)%svname = 'temp20'
    stateinfo(iadd+9)%svname = 'uppertemp'
    stateinfo(iadd+10)%svname = 'lowertemp'
    stateinfo(iadd+1:iadd+2)%ndim = 2
    stateinfo(iadd+3)%ndim = 3
    stateinfo(iadd+4:iadd+5)%ndim = 2
    stateinfo(iadd+6)%ndim = 3
    stateinfo(iadd+7:iadd+10)%ndim = 2
    stateinfo(iadd+1)%dims(1:2) = (/nl,n/)
    stateinfo(iadd+2)%dims(1:2) = (/nl,n/)
    stateinfo(iadd+3)%dims(1:3) = (/ns,nl,n/)
    stateinfo(iadd+4)%dims(1:2) = (/nl,n/)
    stateinfo(iadd+5)%dims(1:2) = (/nl,n/)
    stateinfo(iadd+6)%dims(1:3) = (/ns,nl,n/)
    DO i = iadd+7,iadd+10
      stateinfo(i)%dims(1:2) = (/nl,n/)
    ENDDO
    IF(ALLOCATED(lakestate%water))THEN
      stateinfo(iadd+1)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+1)%svpoint%d2(stateinfo(iadd+1)%dims(1),stateinfo(iadd+1)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%temp))THEN
      stateinfo(iadd+2)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+2)%svpoint%d2(stateinfo(iadd+2)%dims(1),stateinfo(iadd+2)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%conc))THEN
      stateinfo(iadd+3)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+3)%svpoint%d3(stateinfo(iadd+3)%dims(1),stateinfo(iadd+3)%dims(2),stateinfo(iadd+3)%dims(3)))
    ENDIF
    IF(ALLOCATED(lakestate%TPmean))THEN
      stateinfo(iadd+4)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+4)%svpoint%d2(stateinfo(iadd+4)%dims(1),stateinfo(iadd+4)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%slowwater))THEN
      stateinfo(iadd+5)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+5)%svpoint%d2(stateinfo(iadd+5)%dims(1),stateinfo(iadd+5)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%concslow))THEN
      stateinfo(iadd+6)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+6)%svpoint%d3(stateinfo(iadd+6)%dims(1),stateinfo(iadd+6)%dims(2),stateinfo(iadd+6)%dims(3)))
    ENDIF
    IF(ALLOCATED(lakestate%temp10))THEN
      stateinfo(iadd+7)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+7)%svpoint%d2(stateinfo(iadd+7)%dims(1),stateinfo(iadd+7)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%temp20))THEN
      stateinfo(iadd+8)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+8)%svpoint%d2(stateinfo(iadd+8)%dims(1),stateinfo(iadd+8)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%uppertemp))THEN
      stateinfo(iadd+9)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+9)%svpoint%d2(stateinfo(iadd+9)%dims(1),stateinfo(iadd+9)%dims(2)))
    ENDIF
    IF(ALLOCATED(lakestate%lowertemp))THEN
      stateinfo(iadd+10)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+10)%svpoint%d2(stateinfo(iadd+10)%dims(1),stateinfo(iadd+10)%dims(2)))
    ENDIF
    stateinfo(iadd+1)%svpoint%d2  => lakestate%water
    stateinfo(iadd+2)%svpoint%d2  => lakestate%temp
    stateinfo(iadd+3)%svpoint%d3  => lakestate%conc
    stateinfo(iadd+4)%svpoint%d2  => lakestate%TPmean
    stateinfo(iadd+5)%svpoint%d2  => lakestate%slowwater
    stateinfo(iadd+6)%svpoint%d3  => lakestate%concslow
    stateinfo(iadd+7)%svpoint%d2  => lakestate%temp10
    stateinfo(iadd+8)%svpoint%d2  => lakestate%temp20
    stateinfo(iadd+9)%svpoint%d2  => lakestate%uppertemp
    stateinfo(iadd+10)%svpoint%d2 => lakestate%lowertemp

    !Set state information for miscstatetype states
    iadd=nfrozenstates+nsoilstates+naquiferstates+nriverstates+nlakestates
    stateinfo(iadd+1:iadd+nmiscstates)%stname = 'miscstate'
    stateinfo(iadd+1)%svname = 'temp5'
    stateinfo(iadd+2)%svname = 'temp10'
    stateinfo(iadd+3)%svname = 'temp20'
    stateinfo(iadd+4)%svname = 'temp30'
    stateinfo(iadd+5)%svname = 'gdd'
    stateinfo(iadd+6)%svname = 'gsbegin'    !special integer i=65
    stateinfo(iadd+7)%svname = 'nextirrigation'
    stateinfo(iadd+8)%svname = 'cnextirrigation'
    stateinfo(iadd+9)%svname = 'updatestationsarcorr'
    stateinfo(iadd+10)%svname = 'floodwater'
    stateinfo(iadd+11)%svname = 'cfloodwater'
    stateinfo(iadd+12)%svname = 'partT1sf'
    stateinfo(iadd+1:iadd+4)%ndim = 1
    stateinfo(iadd+5:iadd+6)%ndim = 3
    stateinfo(iadd+7)%ndim = 2
    stateinfo(iadd+8)%ndim = 3
    stateinfo(iadd+9)%ndim = 1
    stateinfo(iadd+10)%ndim = 2
    stateinfo(iadd+11)%ndim = 3
    stateinfo(iadd+12)%ndim = 2
    DO i = iadd+1,iadd+4
      stateinfo(i)%dims(1) = n
    ENDDO
    stateinfo(iadd+5)%dims(1:3) = (/2,nc,n/)
    stateinfo(iadd+6)%dims(1:3) = (/2,nc,n/)
    stateinfo(iadd+7)%dims(1:2) = (/nc,n/)
    stateinfo(iadd+8)%dims(1:3) = (/ns,nc,n/)
    stateinfo(iadd+9)%dims(1) = n
    stateinfo(iadd+10)%dims(1:2) = (/2,n/)
    stateinfo(iadd+11)%dims(1:3) = (/ns,2,n/)
    stateinfo(iadd+12)%dims(1:2) = (/nc,n/)
    IF(ALLOCATED(miscstate%temp5))THEN
      stateinfo(iadd+1)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+1)%svpoint%d1(stateinfo(iadd+1)%dims(1)))
    ENDIF
    IF(ALLOCATED(miscstate%temp10))THEN
      stateinfo(iadd+2)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+2)%svpoint%d1(stateinfo(iadd+2)%dims(1)))
    ENDIF
    IF(ALLOCATED(miscstate%temp20))THEN
      stateinfo(iadd+3)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+3)%svpoint%d1(stateinfo(iadd+3)%dims(1)))
    ENDIF
    IF(ALLOCATED(miscstate%temp30))THEN
      stateinfo(iadd+4)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+4)%svpoint%d1(stateinfo(iadd+4)%dims(1)))
    ENDIF
    IF(ALLOCATED(miscstate%gdd))THEN
      stateinfo(iadd+5)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+5)%svpoint%d3(stateinfo(iadd+5)%dims(1),stateinfo(iadd+5)%dims(2),stateinfo(iadd+5)%dims(3)))
    ENDIF
    IF(ALLOCATED(miscstate%gsbegin))THEN
      stateinfo(iadd+6)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+6)%svpoint%d3i(stateinfo(iadd+6)%dims(1),stateinfo(iadd+6)%dims(2),stateinfo(iadd+6)%dims(3)))
    ENDIF
    IF(ALLOCATED(miscstate%nextirrigation))THEN
      stateinfo(iadd+7)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+7)%svpoint%d2(stateinfo(iadd+7)%dims(1),stateinfo(iadd+7)%dims(2)))
    ENDIF
    IF(ALLOCATED(miscstate%cnextirrigation))THEN
      stateinfo(iadd+8)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+8)%svpoint%d3(stateinfo(iadd+8)%dims(1),stateinfo(iadd+8)%dims(2),stateinfo(iadd+8)%dims(3)))
    ENDIF
    IF(ALLOCATED(miscstate%updatestationsarcorr))THEN
      stateinfo(iadd+9)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+9)%svpoint%d1(stateinfo(iadd+9)%dims(1)))
    ENDIF
    IF(ALLOCATED(miscstate%floodwater))THEN
      stateinfo(iadd+10)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+10)%svpoint%d2(stateinfo(iadd+10)%dims(1),stateinfo(iadd+10)%dims(2)))
    ENDIF
    IF(ALLOCATED(miscstate%cfloodwater))THEN
      stateinfo(iadd+11)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+11)%svpoint%d3(stateinfo(iadd+11)%dims(1),stateinfo(iadd+11)%dims(2),stateinfo(iadd+11)%dims(3)))
    ENDIF
    IF(ALLOCATED(miscstate%partT1sf))THEN
      stateinfo(iadd+12)%allok = .TRUE.
      ALLOCATE(stateinfo(iadd+12)%svpoint%d2(stateinfo(iadd+12)%dims(1),stateinfo(iadd+12)%dims(2)))
    ENDIF
    stateinfo(iadd+1)%svpoint%d1  => miscstate%temp5
    stateinfo(iadd+2)%svpoint%d1  => miscstate%temp10
    stateinfo(iadd+3)%svpoint%d1  => miscstate%temp20
    stateinfo(iadd+4)%svpoint%d1  => miscstate%temp30
    stateinfo(iadd+5)%svpoint%d3  => miscstate%gdd
    stateinfo(iadd+6)%svpoint%d3i => miscstate%gsbegin
    stateinfo(iadd+7)%svpoint%d2  => miscstate%nextirrigation
    stateinfo(iadd+8)%svpoint%d3  => miscstate%cnextirrigation
    stateinfo(iadd+9)%svpoint%d1  => miscstate%updatestationsarcorr
    stateinfo(iadd+10)%svpoint%d2 => miscstate%floodwater
    stateinfo(iadd+11)%svpoint%d3 => miscstate%cfloodwater
    stateinfo(iadd+12)%svpoint%d2  => miscstate%partT1sf

  ENDSUBROUTINE set_stateinfo

  !>Allocate state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_model_states(n,ns,nc,na,nsl,nr,nl,ml,mt,isN,isP,isC,isS,  &
                                   isflood,ist1,ist2,iswetl,isglac,islrice, &
                                   isirr,isqar,isgsm,isinf,  &
                                   frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    INTEGER,INTENT(IN) :: ml      !<Maximum river lag steps
    INTEGER,INTENT(IN) :: mt      !<timestep per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: isflood !<Status of flooded area simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN) :: isinf   !<Status of infiltration model == 1
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Soil states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    CALL allocate_frozenstate_variables(n,ns,nc,nr,nl,isglac,islrice,frozenstate)
    CALL allocate_soilstate_variables(n,ns,nc,nsl,isN,isP,isC,isS,ist1,isinf,soilstate)
    CALL allocate_aquiferstate_variables(na,ns,aquiferstate)
    CALL allocate_riverstate_variables(n,ns,nr,ml,mt,isN,isP,isC,isS,ist1,iswetl,riverstate)
    CALL allocate_lakestate_variables(n,ns,nl,isN,isP,isC,isS,ist2,lakestate)
    CALL allocate_miscstate_variables(n,ns,nc,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,miscstate)

  END SUBROUTINE allocate_model_states

  !---------------------------------------------------------------
  !>Deallocate model states 
  !---------------------------------------------------------------
  SUBROUTINE deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    CALL deallocate_frozenstate_variables(frozenstate)
    CALL deallocate_soilstate_variables(soilstate)
    CALL deallocate_aquiferstate_variables(aquiferstate)
    CALL deallocate_riverstate_variables(riverstate)
    CALL deallocate_lakestate_variables(lakestate)
    CALL deallocate_miscstate_variables(miscstate)

  END SUBROUTINE deallocate_model_states

  !---------------------------------------------------------------
  !>Initiate states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_state_zero(ns,na,ist1,ist2,iswetl,isglac,islrice,isirr,isqar,isflood,isgsm,isN,isP,isC,isS,   &
                                     frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isglac    !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar     !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isflood   !<Status of flooded area simulation
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states

    CALL initiate_frozenstate_zero(ns,isglac,islrice,frozenstate)
    CALL initiate_soilstate_zero(ns,soilstate)
    CALL initiate_aquiferstate_zero(ns,na,aquiferstate)
    CALL initiate_riverstate_zero(ns,isP,isN.OR.isP.OR.isC.OR.isS,isS,ist1,iswetl,riverstate)
    CALL initiate_lakestate_zero(ns,isN.OR.isP.OR.isC,ist2,lakestate)
    CALL initiate_miscstate_zero(ns,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,miscstate)
  
  END SUBROUTINE initiate_state_zero

  !>Allocate snow and ice state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_frozenstate_variables(n,ns,nc,nr,nl,isglac,islrice,frozenstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice states

    IF(.NOT.ALLOCATED(frozenstate%snow))      ALLOCATE(frozenstate%snow(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowage))   ALLOCATE(frozenstate%snowage(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowdepth))   ALLOCATE(frozenstate%snowdepth(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowcov))   ALLOCATE(frozenstate%snowcov(nc,n))
    IF(.NOT.ALLOCATED(frozenstate%snowmax))   ALLOCATE(frozenstate%snowmax(nc,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(frozenstate%csnow))   ALLOCATE(frozenstate%csnow(ns,nc,n))
    ENDIF  
    IF(isglac)THEN
      IF(.NOT.ALLOCATED(frozenstate%glacvol))   ALLOCATE(frozenstate%glacvol(n))
    ENDIF
    IF(islrice)THEN
      IF(.NOT.ALLOCATED(frozenstate%lakesnow))     ALLOCATE(frozenstate%lakesnow(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakesnowage))  ALLOCATE(frozenstate%lakesnowage(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakesnowdepth))  ALLOCATE(frozenstate%lakesnowdepth(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakeice))      ALLOCATE(frozenstate%lakeice(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakebice))     ALLOCATE(frozenstate%lakebice(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%lakeicecov))   ALLOCATE(frozenstate%lakeicecov(nl,n))
      IF(.NOT.ALLOCATED(frozenstate%riversnow))    ALLOCATE(frozenstate%riversnow(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riversnowage)) ALLOCATE(frozenstate%riversnowage(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riversnowdepth)) ALLOCATE(frozenstate%riversnowdepth(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riverice))     ALLOCATE(frozenstate%riverice(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%riverbice))    ALLOCATE(frozenstate%riverbice(nr,n))
      IF(.NOT.ALLOCATED(frozenstate%rivericecov))  ALLOCATE(frozenstate%rivericecov(nr,n))
    ENDIF  
    
  END SUBROUTINE allocate_frozenstate_variables

  !>Deallocate snow and ice states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_frozenstate_variables(frozenstate)

    !Argument declarations
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<River states

    IF(ALLOCATED(frozenstate%snow))    DEALLOCATE(frozenstate%snow)
    IF(ALLOCATED(frozenstate%csnow))   DEALLOCATE(frozenstate%csnow)
    IF(ALLOCATED(frozenstate%snowage)) DEALLOCATE(frozenstate%snowage)
    IF(ALLOCATED(frozenstate%snowdepth)) DEALLOCATE(frozenstate%snowdepth)
    IF(ALLOCATED(frozenstate%snowcov)) DEALLOCATE(frozenstate%snowcov)
    IF(ALLOCATED(frozenstate%snowmax)) DEALLOCATE(frozenstate%snowmax)
    IF(ALLOCATED(frozenstate%glacvol)) DEALLOCATE(frozenstate%glacvol)
    IF(ALLOCATED(frozenstate%lakesnow))  DEALLOCATE(frozenstate%lakesnow)
    IF(ALLOCATED(frozenstate%lakesnowage)) DEALLOCATE(frozenstate%lakesnowage)
    IF(ALLOCATED(frozenstate%lakesnowdepth)) DEALLOCATE(frozenstate%lakesnowdepth)
    IF(ALLOCATED(frozenstate%lakeice))   DEALLOCATE(frozenstate%lakeice)
    IF(ALLOCATED(frozenstate%lakebice))  DEALLOCATE(frozenstate%lakebice)
    IF(ALLOCATED(frozenstate%lakeicecov))  DEALLOCATE(frozenstate%lakeicecov)
    IF(ALLOCATED(frozenstate%riversnow)) DEALLOCATE(frozenstate%riversnow)
    IF(ALLOCATED(frozenstate%riversnowage)) DEALLOCATE(frozenstate%riversnowage)
    IF(ALLOCATED(frozenstate%riversnowdepth)) DEALLOCATE(frozenstate%riversnowdepth)
    IF(ALLOCATED(frozenstate%riverice))  DEALLOCATE(frozenstate%riverice)
    IF(ALLOCATED(frozenstate%riverbice)) DEALLOCATE(frozenstate%riverbice)
    IF(ALLOCATED(frozenstate%rivericecov))  DEALLOCATE(frozenstate%rivericecov)
    
  END SUBROUTINE deallocate_frozenstate_variables

  !---------------------------------------------------------------
  !>Initiate submodel snow and ice states by temporary storage of initial state
  !---------------------------------------------------------------
  SUBROUTINE initiate_frozenstate_submodel(ns,n,indexarray,isglac,islrice,frozenstate,frozenstate2)
  
    INTEGER, INTENT(IN) :: ns             !<number of substances
    INTEGER, INTENT(IN) :: n              !<number of subbasins in states to be set
    INTEGER, INTENT(IN) :: indexarray(n)  !<index for basemodel
    LOGICAL,INTENT(IN)  :: isglac         !<Status of glacier simulation
    LOGICAL,INTENT(IN)  :: islrice        !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice state, submodel
    TYPE(snowicestatetype),INTENT(IN) :: frozenstate2    !<Snow and ice state, whole model setup

    frozenstate%snow(:,:) = frozenstate2%snow(:,indexarray(:))
    frozenstate%snowage(:,:) = frozenstate2%snowage(:,indexarray(:))
    frozenstate%snowdepth(:,:) = frozenstate2%snowdepth(:,indexarray(:))
    frozenstate%snowcov(:,:) = frozenstate2%snowcov(:,indexarray(:))
    frozenstate%snowmax(:,:) = frozenstate2%snowmax(:,indexarray(:))
    IF(ns>0)THEN
      frozenstate%csnow(:,:,:) = frozenstate2%csnow(:,:,indexarray(:))
    ENDIF
    IF(isglac) frozenstate%glacvol(:) = frozenstate2%glacvol(indexarray(:))
    IF(islrice)THEN
      frozenstate%lakesnow(:,:)     = frozenstate2%lakesnow(:,indexarray(:))
      frozenstate%lakesnowage(:,:)  = frozenstate2%lakesnowage(:,indexarray(:))
      frozenstate%lakesnowdepth(:,:)  = frozenstate2%lakesnowdepth(:,indexarray(:))
      frozenstate%lakeice(:,:)      = frozenstate2%lakeice(:,indexarray(:))
      frozenstate%lakebice(:,:)     = frozenstate2%lakebice(:,indexarray(:))
      frozenstate%lakeicecov(:,:)   = frozenstate2%lakeicecov(:,indexarray(:))
      frozenstate%riversnow(:,:)    = frozenstate2%riversnow(:,indexarray(:))
      frozenstate%riversnowage(:,:) = frozenstate2%riversnowage(:,indexarray(:))
      frozenstate%riversnowdepth(:,:) = frozenstate2%riversnowdepth(:,indexarray(:))
      frozenstate%riverice(:,:)     = frozenstate2%riverice(:,indexarray(:))
      frozenstate%riverbice(:,:)    = frozenstate2%riverbice(:,indexarray(:))
      frozenstate%rivericecov(:,:)  = frozenstate2%rivericecov(:,indexarray(:))
    ENDIF
  
  END SUBROUTINE initiate_frozenstate_submodel

  !>Initiate snow and ice states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_frozenstate_zero(ns,isglac,islrice,frozenstate)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    LOGICAL,INTENT(IN)  :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN)  :: islrice   !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice states

    frozenstate%snow = 0.
    frozenstate%snowage = 0.
    frozenstate%snowdepth = 0.
    frozenstate%snowcov = 0.
    frozenstate%snowmax = 0.
    IF(ns>0)THEN
      frozenstate%csnow = 0.
    ENDIF 
    IF(isglac) frozenstate%glacvol = 0.
    IF(islrice)THEN
      frozenstate%lakesnow(:,:) = 0.
      frozenstate%lakesnowage(:,:) = 0.
      frozenstate%lakesnowdepth(:,:) = 0.
      frozenstate%lakeice(:,:) = 0.
      frozenstate%lakebice(:,:) = 0.
      frozenstate%lakeicecov(:,:) = 0.
      frozenstate%riversnow(:,:) = 0.
      frozenstate%riversnowage(:,:) = 0.
      frozenstate%riversnowdepth(:,:) = 0.
      frozenstate%riverice(:,:) = 0.
      frozenstate%riverbice(:,:) = 0.
      frozenstate%rivericecov(:,:) = 0.
    ENDIF
  
  END SUBROUTINE initiate_frozenstate_zero

  !>Allocate soil state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_soilstate_variables(n,ns,nc,nsl,isN,isP,isC,isS,ist1,isinf,soilstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1simulation
    LOGICAL,INTENT(IN) :: isinf   !<Status of infiltration model == 1
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states

    IF(.NOT.ALLOCATED(soilstate%water))    ALLOCATE(soilstate%water(nsl,nc,n))
    IF(.NOT.ALLOCATED(soilstate%temp))     ALLOCATE(soilstate%temp(nsl,nc,n))
    IF(.NOT.ALLOCATED(soilstate%deeptemp)) ALLOCATE(soilstate%deeptemp(nc,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(soilstate%conc))      ALLOCATE(soilstate%conc(ns,nsl,nc,n))
      IF(isN)THEN
        IF(.NOT.ALLOCATED(soilstate%humusN)) ALLOCATE(soilstate%humusN(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%fastN))  ALLOCATE(soilstate%fastN(nsl,nc,n))
      ENDIF
      IF(isP)THEN
        IF(.NOT.ALLOCATED(soilstate%humusP)) ALLOCATE(soilstate%humusP(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%fastP))  ALLOCATE(soilstate%fastP(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%partP))  ALLOCATE(soilstate%partP(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%PPrelpool)) ALLOCATE(soilstate%PPrelpool(nc,n))
      ENDIF
      IF(isC)THEN
        IF(.NOT.ALLOCATED(soilstate%humusC)) ALLOCATE(soilstate%humusC(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%fastC))  ALLOCATE(soilstate%fastC(nsl,nc,n))
        IF(.NOT.ALLOCATED(soilstate%oldgrw)) ALLOCATE(soilstate%oldgrw(nc,n))
      ENDIF
      IF(isS)THEN
        IF(.NOT.ALLOCATED(soilstate%Srelpool)) ALLOCATE(soilstate%Srelpool(nc,n))
      ENDIF
      IF(ist1)THEN
        IF(.NOT.ALLOCATED(soilstate%partT1)) ALLOCATE(soilstate%partT1(nsl,nc,n))
      ENDIF   
    ENDIF
    IF(isinf)THEN
      IF(.NOT.ALLOCATED(soilstate%icelens)) ALLOCATE(soilstate%icelens(nc,n))
    ENDIF

  END SUBROUTINE allocate_soilstate_variables

  !>Deallocate soil states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_soilstate_variables(soilstate)

    !Argument declarations
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states

    IF(ALLOCATED(soilstate%water))    DEALLOCATE(soilstate%water)
    IF(ALLOCATED(soilstate%temp))     DEALLOCATE(soilstate%temp)
    IF(ALLOCATED(soilstate%deeptemp)) DEALLOCATE(soilstate%deeptemp)
    IF(ALLOCATED(soilstate%conc))     DEALLOCATE(soilstate%conc)
    IF(ALLOCATED(soilstate%humusN))   DEALLOCATE(soilstate%humusN)
    IF(ALLOCATED(soilstate%fastN))    DEALLOCATE(soilstate%fastN)
    IF(ALLOCATED(soilstate%humusP))   DEALLOCATE(soilstate%humusP)
    IF(ALLOCATED(soilstate%fastP))    DEALLOCATE(soilstate%fastP)
    IF(ALLOCATED(soilstate%partP))    DEALLOCATE(soilstate%partP)
    IF(ALLOCATED(soilstate%PPrelpool))DEALLOCATE(soilstate%PPrelpool)
    IF(ALLOCATED(soilstate%humusC))   DEALLOCATE(soilstate%humusC)
    IF(ALLOCATED(soilstate%fastC))    DEALLOCATE(soilstate%fastC)
    IF(ALLOCATED(soilstate%oldgrw))   DEALLOCATE(soilstate%oldgrw)
    IF(ALLOCATED(soilstate%Srelpool)) DEALLOCATE(soilstate%Srelpool)
    IF(ALLOCATED(soilstate%partT1))   DEALLOCATE(soilstate%partT1)
    IF(ALLOCATED(soilstate%icelens))  DEALLOCATE(soilstate%icelens)

  END SUBROUTINE deallocate_soilstate_variables

  !>Initiate submodel soil states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_soilstate_submodel(ns,n,indexarray,soilstate,soilstate2)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    INTEGER, INTENT(IN) :: n         !<number of subbasins in states to be set
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil state, submodel
    TYPE(soilstatetype),INTENT(IN) :: soilstate2    !<Soil state, whole model setup

    soilstate%water(:,:,:)  = soilstate2%water(:,:,indexarray(:))
    soilstate%temp(:,:,:)   = soilstate2%temp(:,:,indexarray(:))
    soilstate%deeptemp(:,:) = soilstate2%deeptemp(:,indexarray(:))
    IF(ns>0)THEN
      soilstate%conc(:,:,:,:) = soilstate2%conc(:,:,:,indexarray(:))
      IF(ALLOCATED(soilstate%humusN))THEN
        soilstate%humusN(:,:,:) = soilstate2%humusN(:,:,indexarray(:))
        soilstate%fastN(:,:,:)  = soilstate2%fastN(:,:,indexarray(:))
      ENDIF  
      IF(ALLOCATED(soilstate%humusP))THEN
        soilstate%humusP(:,:,:) = soilstate2%humusP(:,:,indexarray(:))
        soilstate%fastP(:,:,:)  = soilstate2%fastP(:,:,indexarray(:))
        soilstate%partP(:,:,:)  = soilstate2%partP(:,:,indexarray(:))
        soilstate%PPrelpool(:,:) = soilstate2%PPrelpool(:,indexarray(:))
      ENDIF  
      IF(ALLOCATED(soilstate%humusC))THEN
        soilstate%humusC(:,:,:) = soilstate2%humusC(:,:,indexarray(:))
        soilstate%fastC(:,:,:)  = soilstate2%fastC(:,:,indexarray(:))
        soilstate%oldgrw(:,:)   = soilstate2%oldgrw(:,indexarray(:))
      ENDIF  
      IF(ALLOCATED(soilstate%Srelpool))THEN
        soilstate%Srelpool(:,:) = soilstate2%Srelpool(:,indexarray(:))
      ENDIF  
      IF(ALLOCATED(soilstate%partT1))THEN
        soilstate%partT1(:,:,:) = soilstate2%partT1(:,:,indexarray(:))
      ENDIF  
    ENDIF
    IF(ALLOCATED(soilstate%icelens)) soilstate%icelens(:,:) = soilstate2%icelens(:,indexarray(:))
  
  END SUBROUTINE initiate_soilstate_submodel

  !>Initiate soil states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_soilstate_zero(ns,soilstate)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states

    soilstate%water = 0.
    soilstate%temp = 0.
    soilstate%deeptemp = 0.
    IF(ns>0)THEN
      soilstate%conc = 0.
      IF(ALLOCATED(soilstate%humusN))THEN
        soilstate%humusN(:,:,:) = 0.
        soilstate%fastN(:,:,:)  = 0.
      ENDIF  
      IF(ALLOCATED(soilstate%humusP))THEN
        soilstate%humusP(:,:,:) = 0.
        soilstate%fastP(:,:,:)  = 0.
        soilstate%partP(:,:,:)  = 0.
        soilstate%PPrelpool = 0.
      ENDIF  
      IF(ALLOCATED(soilstate%humusC))THEN
        soilstate%humusC(:,:,:) = 0.
        soilstate%fastC(:,:,:)  = 0.
        soilstate%oldgrw(:,:)   = 0.
      ENDIF  
      IF(ALLOCATED(soilstate%Srelpool))THEN
        soilstate%Srelpool = 0.
      ENDIF  
      IF(ALLOCATED(soilstate%partT1))THEN
        soilstate%partT1(:,:,:) = 0.
      ENDIF        
    ENDIF
    IF(ALLOCATED(soilstate%icelens)) soilstate%icelens = 0
  
  END SUBROUTINE initiate_soilstate_zero

  !>Allocate aquifer state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_aquiferstate_variables(na,ns,aquiferstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states

    IF(na>0)THEN
      IF(.NOT.ALLOCATED(aquiferstate%water))   ALLOCATE(aquiferstate%water(na))
      IF(.NOT.ALLOCATED(aquiferstate%lastrecharge))   ALLOCATE(aquiferstate%lastrecharge(na))
      IF(.NOT.ALLOCATED(aquiferstate%nextoutflow))   ALLOCATE(aquiferstate%nextoutflow(na))
      IF(ns>0)THEN
        IF(.NOT.ALLOCATED(aquiferstate%conc))  ALLOCATE(aquiferstate%conc(ns,na))
        IF(.NOT.ALLOCATED(aquiferstate%clastrecharge))  ALLOCATE(aquiferstate%clastrecharge(ns,na))
        IF(.NOT.ALLOCATED(aquiferstate%cnextoutflow))  ALLOCATE(aquiferstate%cnextoutflow(ns,na))
      ENDIF
    ENDIF

  END SUBROUTINE allocate_aquiferstate_variables

  !>Deallocate aquifer states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_aquiferstate_variables(aquiferstate)

    !Argument declarations
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states

    IF(ALLOCATED(aquiferstate%water))         DEALLOCATE(aquiferstate%water)
    IF(ALLOCATED(aquiferstate%lastrecharge))  DEALLOCATE(aquiferstate%lastrecharge)
    IF(ALLOCATED(aquiferstate%nextoutflow))   DEALLOCATE(aquiferstate%nextoutflow)
    IF(ALLOCATED(aquiferstate%conc))          DEALLOCATE(aquiferstate%conc)
    IF(ALLOCATED(aquiferstate%clastrecharge)) DEALLOCATE(aquiferstate%clastrecharge)
    IF(ALLOCATED(aquiferstate%cnextoutflow))  DEALLOCATE(aquiferstate%cnextoutflow)

  END SUBROUTINE deallocate_aquiferstate_variables

  !>Initiate submodel aquifer states by temporary storage of initial
  !>state (copy)
  !---------------------------------------------------------------
  SUBROUTINE initiate_aquiferstate_submodel(ns,na,aquiferstate,aquiferstate2)
  
    INTEGER, INTENT(IN) :: ns        !<number of substances
    INTEGER, INTENT(IN) :: na        !<Number of aquifers
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer state, submodel
    TYPE(aquiferstatetype),INTENT(IN) :: aquiferstate2    !<aquifer state, whole model setup

    IF(na>0.)THEN
      aquiferstate%water(:)        = aquiferstate2%water(:)
      aquiferstate%lastrecharge(:) = aquiferstate2%lastrecharge(:)
      aquiferstate%nextoutflow(:)  = aquiferstate2%nextoutflow(:)
      IF(ns>0)THEN
        aquiferstate%conc(:,:)          = aquiferstate2%conc(:,:)
        aquiferstate%clastrecharge(:,:) = aquiferstate2%clastrecharge(:,:)
        aquiferstate%cnextoutflow(:,:)  = aquiferstate2%cnextoutflow(:,:)
      ENDIF
    ENDIF
  
  END SUBROUTINE initiate_aquiferstate_submodel

  !>Initiate aquifer states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_aquiferstate_zero(ns,na,aquiferstate)
  
    INTEGER, INTENT(IN) :: ns      !<number of substances
    INTEGER, INTENT(IN) :: na      !<number of aquifers
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states

    IF(na>0)THEN
      aquiferstate%water = 0.
      aquiferstate%lastrecharge = 0.
      aquiferstate%nextoutflow = 0.
      IF(ns>0)THEN
        aquiferstate%conc = 0.
        aquiferstate%clastrecharge = 0.
        aquiferstate%cnextoutflow = 0.
      ENDIF
    ENDIF
  
  END SUBROUTINE initiate_aquiferstate_zero

  !>Allocate river state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_riverstate_variables(n,ns,nr,ml,mt,isN,isP,isC,isS,ist1,iswetl,riverstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<number of subbasins
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: nr      !<number of river types
    INTEGER,INTENT(IN) :: ml      !<maximum river lag steps
    INTEGER,INTENT(IN) :: mt      !<timestep per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    IF(.NOT.ALLOCATED(riverstate%water)) ALLOCATE(riverstate%water(nr,n))
    IF(.NOT.ALLOCATED(riverstate%temp)) ALLOCATE(riverstate%temp(nr,n))
    IF(.NOT.ALLOCATED(riverstate%qqueue)) ALLOCATE(riverstate%qqueue(0:ml,nr,n))
    IF(.NOT.ALLOCATED(riverstate%Qmean)) ALLOCATE(riverstate%Qmean(nr,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(riverstate%conc)) ALLOCATE(riverstate%conc(ns,nr,n))
      IF(.NOT.ALLOCATED(riverstate%cqueue)) ALLOCATE(riverstate%cqueue(ns,0:ml,nr,n))
    ENDIF
    IF(isP)THEN
      IF(.NOT.ALLOCATED(riverstate%Psed)) ALLOCATE(riverstate%Psed(nr,n))
    ENDIF
    IF(isN.OR.isP.OR.isC.OR.isS)THEN
      IF(.NOT.ALLOCATED(riverstate%TPmean)) ALLOCATE(riverstate%TPmean(nr,n))
      IF(.NOT.ALLOCATED(riverstate%temp10)) ALLOCATE(riverstate%temp10(nr,n))
      IF(.NOT.ALLOCATED(riverstate%temp20)) ALLOCATE(riverstate%temp20(nr,n))
    ENDIF  
    IF(isN.OR.isP.OR.ist1)THEN
      IF(.NOT.ALLOCATED(riverstate%Qdayacc)) ALLOCATE(riverstate%Qdayacc(mt,nr,n))
      IF(.NOT.ALLOCATED(riverstate%Q365)) ALLOCATE(riverstate%Q365(366,nr,n))
    ENDIF
    IF(iswetl)THEN
      IF(.NOT.ALLOCATED(riverstate%cwetland)) ALLOCATE(riverstate%cwetland(ns,nr,n))
    ENDIF
    IF(ist1)THEN
      IF(.NOT.ALLOCATED(riverstate%T1sed)) ALLOCATE(riverstate%T1sed(nr,n))
    ENDIF
    IF(isS)THEN
      IF(.NOT.ALLOCATED(riverstate%Ssed)) ALLOCATE(riverstate%Ssed(nr,n))
    ENDIF

  END SUBROUTINE allocate_riverstate_variables

  !>Deallocate river states
  !---------------------------------------------------------------
  SUBROUTINE deallocate_riverstate_variables(riverstate)

    !Argument declarations
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    IF(ALLOCATED(riverstate%water))    DEALLOCATE(riverstate%water)
    IF(ALLOCATED(riverstate%temp))     DEALLOCATE(riverstate%temp)
    IF(ALLOCATED(riverstate%qqueue))   DEALLOCATE(riverstate%qqueue)
    IF(ALLOCATED(riverstate%Qmean))    DEALLOCATE(riverstate%Qmean)
    IF(ALLOCATED(riverstate%conc))     DEALLOCATE(riverstate%conc)
    IF(ALLOCATED(riverstate%cqueue))   DEALLOCATE(riverstate%cqueue)
    IF(ALLOCATED(riverstate%Psed))     DEALLOCATE(riverstate%Psed)
    IF(ALLOCATED(riverstate%TPmean))   DEALLOCATE(riverstate%TPmean)
    IF(ALLOCATED(riverstate%temp10))   DEALLOCATE(riverstate%temp10)
    IF(ALLOCATED(riverstate%temp20))   DEALLOCATE(riverstate%temp20)
    IF(ALLOCATED(riverstate%Qdayacc))  DEALLOCATE(riverstate%Qdayacc)
    IF(ALLOCATED(riverstate%Q365))     DEALLOCATE(riverstate%Q365)
    IF(ALLOCATED(riverstate%cwetland)) DEALLOCATE(riverstate%cwetland)
    IF(ALLOCATED(riverstate%T1sed))    DEALLOCATE(riverstate%T1sed)
    IF(ALLOCATED(riverstate%Ssed))     DEALLOCATE(riverstate%Ssed)

  END SUBROUTINE deallocate_riverstate_variables

  !>Initiate submodel river states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_riverstate_submodel(ns,n,isN,isP,isC,isS,ist1,iswetl,indexarray,  &
                                          riverstate,riverstate2)
  
    INTEGER,INTENT(IN) :: ns        !<number of substances
    INTEGER,INTENT(IN) :: n         !<number of subbasins in states to be set
    LOGICAL,INTENT(IN) :: isN       !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP       !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1      !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: iswetl    !<Status of wetland simulation
    INTEGER,INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River state, submodel
    TYPE(riverstatetype),INTENT(IN) :: riverstate2    !<River state, whole model setup

    riverstate%water(:,:) = riverstate2%water(:,indexarray(:))
    riverstate%temp(:,:) = riverstate2%temp(:,indexarray(:))
    riverstate%qqueue(:,:,:) = riverstate2%qqueue(:,:,indexarray(:))
    riverstate%Qmean(:,:) = riverstate2%Qmean(:,indexarray(:))
    IF(ns>0)THEN
      riverstate%conc(:,:,:) = riverstate2%conc(:,:,indexarray(:))
      riverstate%cqueue(:,:,:,:) = riverstate2%cqueue(:,:,:,indexarray(:))
    ENDIF
    IF(isP)THEN
      riverstate%Psed(:,:) = riverstate2%Psed(:,indexarray(:))
    ENDIF
    IF(isN.OR.isP.OR.isC.OR.isS)THEN
      riverstate%TPmean(:,:) = riverstate2%TPmean(:,indexarray(:))
      riverstate%temp10(:,:) = riverstate2%temp10(:,indexarray(:))
      riverstate%temp20(:,:) = riverstate2%temp20(:,indexarray(:))
    ENDIF
    IF(isN.OR.isP)THEN
      riverstate%Qdayacc(:,:,:) = riverstate2%Qdayacc(:,:,indexarray(:))
      riverstate%Q365(:,:,:) = riverstate2%Q365(:,:,indexarray(:))
    ENDIF
    IF(iswetl)THEN
      riverstate%cwetland(:,:,:) = riverstate2%cwetland(:,:,indexarray(:))
    ENDIF
    IF(ist1)THEN
      riverstate%T1sed(:,:) = riverstate2%T1sed(:,indexarray(:))
    ENDIF
    IF(isS)THEN
      riverstate%Ssed(:,:) = riverstate2%Ssed(:,indexarray(:))
    ENDIF
  
  END SUBROUTINE initiate_riverstate_submodel

  !>Initiate river states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_riverstate_zero(ns,isP,isNPCS,isS,ist1,iswetl,riverstate)
  
    INTEGER,INTENT(IN) :: ns      !<number of substances
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isNPCS  !<Status of N, P, C or S simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states

    riverstate%water = 0.
    riverstate%temp = 0.
    riverstate%qqueue = 0.
    riverstate%Qmean = 0.
    IF(ns>0)THEN
      riverstate%conc = 0.
      riverstate%cqueue = 0.
    ENDIF
    IF(isP)THEN
      riverstate%Psed = 0.
    ENDIF
    IF(isNPCS)THEN
      riverstate%TPmean = 0.
      riverstate%temp10 = 0.
      riverstate%temp20 = 0.
    ENDIF
!    IF(isN.OR.isP)THEN
!      riverstate%Qdayacc = 0.   !initialized to 0.0001 later in HYPE
!      riverstate%Q365 = 0.   !initialized to 0.0001 later in HYPE
!    ENDIF
    IF(iswetl) riverstate%cwetland = 0.
    IF(ist1) riverstate%T1sed = 0.
    IF(isS) riverstate%Ssed = 0.
  
  END SUBROUTINE initiate_riverstate_zero

  !>Allocate lake state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_lakestate_variables(n,ns,nl,isN,isP,isC,isS,ist2,lakestate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<number of subbasins
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: nl      !<number of lake types
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of sediment simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake states

    IF(.NOT.ALLOCATED(lakestate%water)) ALLOCATE(lakestate%water(nl,n))
    IF(.NOT.ALLOCATED(lakestate%temp)) ALLOCATE(lakestate%temp(nl,n))
    IF(ns>0)THEN
      IF(.NOT.ALLOCATED(lakestate%conc)) ALLOCATE(lakestate%conc(ns,nl,n))
      IF(.NOT.ALLOCATED(lakestate%slowwater)) ALLOCATE(lakestate%slowwater(nl,n))
      IF(.NOT.ALLOCATED(lakestate%concslow)) ALLOCATE(lakestate%concslow(ns,nl,n))
    ENDIF  
    IF(isN.OR.isP.OR.isC.OR.isS)THEN
      IF(.NOT.ALLOCATED(lakestate%TPmean)) ALLOCATE(lakestate%TPmean(nl,n))
      IF(.NOT.ALLOCATED(lakestate%temp10)) ALLOCATE(lakestate%temp10(nl,n))
      IF(.NOT.ALLOCATED(lakestate%temp20)) ALLOCATE(lakestate%temp20(nl,n))
    ENDIF  
    IF(ist2)THEN
      IF(.NOT.ALLOCATED(lakestate%uppertemp)) ALLOCATE(lakestate%uppertemp(nl,n))
      IF(.NOT.ALLOCATED(lakestate%lowertemp)) ALLOCATE(lakestate%lowertemp(nl,n))
    ENDIF  

  END SUBROUTINE allocate_lakestate_variables

  !>Deallocate lake states used for temporary storage of initial
  !>state for submodel simulation
  !---------------------------------------------------------------
  SUBROUTINE deallocate_lakestate_variables(lakestate)

    !Argument declarations
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake states

    IF(ALLOCATED(lakestate%water))     DEALLOCATE(lakestate%water)
    IF(ALLOCATED(lakestate%temp))      DEALLOCATE(lakestate%temp)
    IF(ALLOCATED(lakestate%conc))      DEALLOCATE(lakestate%conc)
    IF(ALLOCATED(lakestate%slowwater)) DEALLOCATE(lakestate%slowwater)
    IF(ALLOCATED(lakestate%concslow))  DEALLOCATE(lakestate%concslow)
    IF(ALLOCATED(lakestate%TPmean))    DEALLOCATE(lakestate%TPmean)
    IF(ALLOCATED(lakestate%temp10))    DEALLOCATE(lakestate%temp10)
    IF(ALLOCATED(lakestate%temp20))    DEALLOCATE(lakestate%temp20)
    IF(ALLOCATED(lakestate%uppertemp)) DEALLOCATE(lakestate%uppertemp)
    IF(ALLOCATED(lakestate%lowertemp)) DEALLOCATE(lakestate%lowertemp)

  END SUBROUTINE deallocate_lakestate_variables

  !>Initiate submodel lake states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_lakestate_submodel(ns,n,isNPCS,ist2,indexarray,lakestate,lakestate2)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns       !<number of substances
    INTEGER, INTENT(IN) :: n        !<number of subbasins in states to be set
    LOGICAL, INTENT(IN) :: isNPCS    !<Status of N, P, C or S simulation
    LOGICAL, INTENT(IN) :: ist2     !<Status of T2 simulation
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state, submodel
    TYPE(lakestatetype),INTENT(IN) :: lakestate2    !<Lake state, whole model setup

    lakestate%water(:,:) = lakestate2%water(:,indexarray(:))
    lakestate%temp(:,:) = lakestate2%temp(:,indexarray(:))
    IF(ns>0)THEN
      lakestate%conc(:,:,:) = lakestate2%conc(:,:,indexarray(:))
      lakestate%slowwater(:,:) = lakestate2%slowwater(:,indexarray(:))
      lakestate%concslow(:,:,:) = lakestate2%concslow(:,:,indexarray(:))
    ENDIF
    IF(isNPCS)THEN
      lakestate%TPmean(:,:) = lakestate2%TPmean(:,indexarray(:))
      lakestate%temp10(:,:) = lakestate2%temp10(:,indexarray(:))
      lakestate%temp20(:,:) = lakestate2%temp20(:,indexarray(:))
    ENDIF
    IF(ist2)THEN
      lakestate%uppertemp(:,:) = lakestate2%uppertemp(:,indexarray(:))
      lakestate%lowertemp(:,:) = lakestate2%lowertemp(:,indexarray(:))
    ENDIF  
  
  END SUBROUTINE initiate_lakestate_submodel

  !>Initiate lake states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_lakestate_zero(ns,isNPCS,ist2,lakestate)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns        !<number of substances
    LOGICAL, INTENT(IN) :: isNPCS    !<Status of N, P, C or S simulation
    LOGICAL, INTENT(IN) :: ist2      !<Status of T2 simulation
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<Lake state, submodel

    lakestate%water = 0.
    lakestate%temp = 0.
    IF(ns>0)THEN
      lakestate%conc = 0.
      lakestate%slowwater = 0.
      lakestate%concslow = 0.
    ENDIF
    IF(isNPCS)THEN
      lakestate%TPmean = 0.
      lakestate%temp10 = 0.
      lakestate%temp20 = 0.
    ENDIF
    IF(ist2)THEN
      lakestate%uppertemp(:,:) = 0.
      lakestate%lowertemp(:,:) = 0.
    ENDIF  

  END SUBROUTINE initiate_lakestate_zero

  !>Allocate miscellaneous state variables for the model
  !
  !>\b Consequences Memory will be allocated
  !---------------------------------------
  SUBROUTINE allocate_miscstate_variables(n,ns,nc,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,miscstate)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<number of subbasins
    INTEGER,INTENT(IN) :: ns      !<number of substances
    INTEGER,INTENT(IN) :: nc      !<number of classes
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isflood !<Status of flooded area simulation
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state

    IF(iswetl)THEN
      IF(.NOT.ALLOCATED(miscstate%temp5))  ALLOCATE(miscstate%temp5(n))
      IF(.NOT.ALLOCATED(miscstate%temp30)) ALLOCATE(miscstate%temp30(n))
    ENDIF
    IF(isC)THEN
      IF(.NOT.ALLOCATED(miscstate%temp10)) ALLOCATE(miscstate%temp10(n))
      IF(.NOT.ALLOCATED(miscstate%temp20)) ALLOCATE(miscstate%temp20(n))
    ENDIF
    IF(isgsm)THEN
      IF(.NOT.ALLOCATED(miscstate%gdd))     ALLOCATE(miscstate%gdd(2,nc,n))
      IF(.NOT.ALLOCATED(miscstate%gsbegin)) ALLOCATE(miscstate%gsbegin(2,nc,n))
    ENDIF
    IF(isirr)THEN
      IF(.NOT.ALLOCATED(miscstate%nextirrigation)) ALLOCATE(miscstate%nextirrigation(nc,n))
      IF(ns>0)THEN
        IF(.NOT.ALLOCATED(miscstate%cnextirrigation)) ALLOCATE(miscstate%cnextirrigation(ns,nc,n))
      ENDIF
    ENDIF
    IF(isqar)THEN
      IF(.NOT.ALLOCATED(miscstate%updatestationsarcorr)) ALLOCATE(miscstate%updatestationsarcorr(n))
    ENDIF
    IF(isflood)THEN
      IF(.NOT.ALLOCATED(miscstate%floodwater)) ALLOCATE(miscstate%floodwater(2,n))
      IF(ns>0)THEN
        IF(.NOT.ALLOCATED(miscstate%cfloodwater)) ALLOCATE(miscstate%cfloodwater(ns,2,n))
      ENDIF
    ENDIF
    IF(ist1)THEN
      IF(.NOT.ALLOCATED(miscstate%partT1sf))ALLOCATE(miscstate%partT1sf(nc,n))
    ENDIF   

  END SUBROUTINE allocate_miscstate_variables

  !>Deallocate miscellaneous states used for temporary storage of initial
  !>state for submodel simulation
  !---------------------------------------------------------------
  SUBROUTINE deallocate_miscstate_variables(miscstate)

    !Argument declarations
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state

    IF(ALLOCATED(miscstate%temp5))  DEALLOCATE(miscstate%temp5)
    IF(ALLOCATED(miscstate%temp30)) DEALLOCATE(miscstate%temp30)
    IF(ALLOCATED(miscstate%temp10)) DEALLOCATE(miscstate%temp10)
    IF(ALLOCATED(miscstate%temp20)) DEALLOCATE(miscstate%temp20)
    IF(ALLOCATED(miscstate%gdd))    DEALLOCATE(miscstate%gdd)
    IF(ALLOCATED(miscstate%gsbegin))         DEALLOCATE(miscstate%gsbegin)
    IF(ALLOCATED(miscstate%nextirrigation))  DEALLOCATE(miscstate%nextirrigation)
    IF(ALLOCATED(miscstate%cnextirrigation)) DEALLOCATE(miscstate%cnextirrigation)
    IF(ALLOCATED(miscstate%updatestationsarcorr)) DEALLOCATE(miscstate%updatestationsarcorr)
    IF(ALLOCATED(miscstate%floodwater)) DEALLOCATE(miscstate%floodwater)
    IF(ALLOCATED(miscstate%cfloodwater)) DEALLOCATE(miscstate%cfloodwater)
    IF(ALLOCATED(miscstate%partT1sf))    DEALLOCATE(miscstate%partT1sf)

  END SUBROUTINE deallocate_miscstate_variables

  !>Initiate submodel miscellaneous states by temporary storage of initial
  !>state
  !---------------------------------------------------------------
  SUBROUTINE initiate_miscstate_submodel(ns,n,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,indexarray,miscstate,miscstate2)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns        !<number of substances
    INTEGER, INTENT(IN) :: n         !<number of subbasins in states to be set
    LOGICAL,INTENT(IN)  :: ist1      !<Status of T1 simulation
    LOGICAL,INTENT(IN)  :: isirr     !<Status of irrigation simulation
    LOGICAL,INTENT(IN)  :: isqar     !<Status of updating q with AR
    LOGICAL,INTENT(IN)  :: isflood   !<Status of flooded area simulation
    LOGICAL,INTENT(IN)  :: isgsm     !<Status of growth season model
    LOGICAL,INTENT(IN)  :: iswetl    !<Status of wetland simulation
    LOGICAL,INTENT(IN)  :: isC       !<Status of organic carbon simulation
    INTEGER, INTENT(IN) :: indexarray(n)            !<index for basemodel
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state, submodel
    TYPE(miscstatetype),INTENT(IN) :: miscstate2    !<misc state, whole model setup

    IF(iswetl)THEN
      miscstate%temp5(:)  = miscstate2%temp5(indexarray(:))
      miscstate%temp30(:) = miscstate2%temp30(indexarray(:))
    ENDIF
    IF(isC)THEN
      miscstate%temp10(:) = miscstate2%temp10(indexarray(:))
      miscstate%temp20(:) = miscstate2%temp20(indexarray(:))
    ENDIF
    IF(isgsm)THEN
      miscstate%gdd(:,:,:) = miscstate2%gdd(:,:,indexarray(:))
      miscstate%gsbegin(:,:,:) = miscstate2%gsbegin(:,:,indexarray(:))
    ENDIF
    IF(isirr)THEN
      miscstate%nextirrigation(:,:) = miscstate2%nextirrigation(:,indexarray(:))
      IF(ns>0)THEN
        miscstate%cnextirrigation(:,:,:) = miscstate2%cnextirrigation(:,:,indexarray(:))
      ENDIF
    ENDIF
    IF(isqar)THEN
      miscstate%updatestationsarcorr(:) = miscstate2%updatestationsarcorr(indexarray(:))
    ENDIF
    IF(isflood)THEN
      miscstate%floodwater(:,:) = miscstate2%floodwater(:,indexarray(:))
      IF(ns>0)THEN
        miscstate%cfloodwater(:,:,:) = miscstate2%cfloodwater(:,:,indexarray(:))
      ENDIF
    ENDIF
    IF(ist1)THEN
      miscstate%partT1sf(:,:) = miscstate2%partT1sf(:,indexarray(:))
    ENDIF
  
  END SUBROUTINE initiate_miscstate_submodel

  !>Initiate miscellaneous states to zero
  !---------------------------------------------------------------
  SUBROUTINE initiate_miscstate_zero(ns,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,miscstate)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: ns      !<number of substances
    LOGICAL,INTENT(IN)  :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN)  :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN)  :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN)  :: isflood !<Status of flooded area simulation
    LOGICAL,INTENT(IN)  :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN)  :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN)  :: isC     !<Status of organic carbon simulation
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc state, submodel

    IF(iswetl)THEN
      miscstate%temp5  = 0.
      miscstate%temp30 = 0.
    ENDIF
    IF(isC)THEN
      miscstate%temp10 = 0.
      miscstate%temp20 = 0.
    ENDIF
    IF(isgsm)THEN
      miscstate%gdd = 0.
      miscstate%gsbegin = 0
    ENDIF
    IF(isirr)THEN
      miscstate%nextirrigation = 0.
      IF(ns>0)THEN
        miscstate%cnextirrigation = 0.
      ENDIF
    ENDIF
    IF(isqar)THEN
      miscstate%updatestationsarcorr(:) = 0.
    ENDIF
    IF(isflood)THEN
      miscstate%floodwater(:,:) = 0.
      IF(ns>0)THEN
        miscstate%cfloodwater = 0.
    ENDIF
    ENDIF
    IF(ist1)THEN
      miscstate%partT1sf(:,:) = 0.
    ENDIF
    
  END SUBROUTINE initiate_miscstate_zero

!------------------------------------
!Subroutines for state file handling
!------------------------------------

  !>Get number of frozen states
  !-------------------------------------------------
  SUBROUTINE get_frozenstate_variables_arraysize(n,ns,nc,nr,nl,isglac,islrice,dim)

    !Argument declarations
    INTEGER,INTENT(IN)  :: n       !<Number of subbasins
    INTEGER,INTENT(IN)  :: ns      !<Number of substances
    INTEGER,INTENT(IN)  :: nc      !<Number of classes
    INTEGER,INTENT(IN)  :: nr      !<Number of river types
    INTEGER,INTENT(IN)  :: nl      !<Number of lake types
    LOGICAL,INTENT(IN)  :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN)  :: islrice !<Status of lakeriver ice model
    INTEGER,INTENT(OUT) :: dim     !<Number of array elements

    dim = 0
    dim = dim + n*nc   !snow
    dim = dim + n*nc   !snowage
    dim = dim + n*nc   !snowdepth
    dim = dim + n*nc   !snowcov
    dim = dim + n*nc   !snowmax
    IF(ns>0) dim = dim + n*nc*ns    !csnow
    IF(isglac) dim = dim + n    !glacvol
    IF(islrice)THEN
      dim = dim + nl*n   !lakesnow
      dim = dim + nl*n   !lakesnowage
      dim = dim + nl*n   !lakesnowdepth
      dim = dim + nl*n   !lakeice
      dim = dim + nl*n   !lakebice
      dim = dim + nl*n   !lakeicecov
      dim = dim + nr*n   !riversnow
      dim = dim + nr*n   !riversnowage
      dim = dim + nr*n   !riversnowdepth
      dim = dim + nr*n   !riverice
      dim = dim + nr*n   !riverbice
      dim = dim + nr*n   !rivericecov
    ENDIF  

  END SUBROUTINE get_frozenstate_variables_arraysize

  !>Write snow and ice state variables to array
  !-------------------------------------------------
  SUBROUTINE set_frozenstate_variables_to_array(n,ns,nc,nr,nl,isglac,islrice,frozenstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(IN) :: frozenstate  !<Snow and ice states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snow(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowage(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowdepth(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowcov(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%snowmax(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%csnow(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isglac)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%glacvol(i)
      ENDDO
    ENDIF
    IF(islrice)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakesnow(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakesnowage(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakesnowdepth(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakeice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakebice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%lakeicecov(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riversnow(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riversnowage(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riversnowdepth(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riverice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%riverbice(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = frozenstate%rivericecov(j,i)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_frozenstate_variables_to_array

  !>Read snow and ice state variables from array
  !-------------------------------------------------
  SUBROUTINE set_frozenstate_variables_from_array(n,ns,nsfil,nc,nr,nl,isglac,islrice,frozenstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances read from file
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isglac  !<Status of glacier simulation
    LOGICAL,INTENT(IN) :: islrice !<Status of lakeriver ice model
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate  !<Snow and ice states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snow(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowage(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowdepth(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowcov(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%snowmax(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%csnow(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n * nc * nsfil  !skip saved concentrations
    ENDIF  
    IF(isglac)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%glacvol(i) = array(iarr-arrshift)
      ENDDO
    ENDIF
    IF(islrice)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakesnow(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakesnowage(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakesnowdepth(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakeice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakebice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%lakeicecov(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riversnow(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riversnowage(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riversnowdepth(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riverice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%riverbice(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) frozenstate%rivericecov(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_frozenstate_variables_from_array

  !>Get number of soil states
  !-------------------------------------------------
  SUBROUTINE get_soilstate_variables_arraysize(n,ns,nc,nsl,isN,isP,isC,isS,ist1,isinf,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isinf   !<Status of infiltration model == 1
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = nsl*nc*n        !water
    dim = dim + nsl*nc*n  !temp
    dim = dim + nc*n      !deeptemp
    IF(ns>0)THEN
      dim = dim + ns*nsl*nc*n  !conc
      IF(isN)THEN
        dim = dim + nsl*nc*n   !humusN
        dim = dim + nsl*nc*n   !fastN
      ENDIF
      IF(isP)THEN
        dim = dim + nsl*nc*n   !humusP
        dim = dim + nsl*nc*n   !fastP
        dim = dim + nsl*nc*n   !partP
        dim = dim + nc*n       !PPrelpool
      ENDIF
      IF(isC)THEN
        dim = dim + nsl*nc*n   !humusC
        dim = dim + nsl*nc*n   !fastC
        dim = dim + nc*n       !oldgrw
      ENDIF
      IF(isS)THEN
        dim = dim + nc*n       !Srelpool
      ENDIF
      IF(ist1)THEN
        dim = dim + nsl*nc*n   !partT1
      ENDIF
    ENDIF
    IF(isinf) dim = dim + nc*n      !icelens

  END SUBROUTINE get_soilstate_variables_arraysize

  !>Write soil state variables to array
  !-------------------------------------------------
  SUBROUTINE set_soilstate_variables_to_array(n,ns,nc,nsl,isN,isP,isC,isS,ist1,isinf,soilstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isinf   !<Status of infiltration model == 1
    TYPE(soilstatetype),INTENT(IN) :: soilstate  !<Soil states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%water(k,j,i)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%temp(k,j,i)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%deeptemp(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,nsl
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%conc(l,k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      IF(iarr>iarrlast) RETURN    !Escape
      IF(isN)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%humusN(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%fastN(k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(iarr>iarrlast) RETURN    !Escape
      IF(isP)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%humusP(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%fastP(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%partP(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%PPrelpool(j,i)
          ENDDO
        ENDDO
      ENDIF
      IF(iarr>iarrlast) RETURN    !Escape
      IF(isC)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%humusC(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%fastC(k,j,i)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%oldgrw(j,i)
          ENDDO
        ENDDO
      ENDIF
      IF(isS)THEN
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%Srelpool(j,i)
          ENDDO
        ENDDO
      ENDIF
      IF(ist1)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = soilstate%partT1(k,j,i)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    IF(isinf)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = REAL(soilstate%icelens(j,i))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_soilstate_variables_to_array

  !>Read soil state variables from array
  !-------------------------------------------------
  SUBROUTINE set_soilstate_variables_from_array(n,ns,nsfil,nc,nsl,isN,isP,isC,isS,ist1,isinf,soilstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances read from file
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    INTEGER,INTENT(IN) :: nsl     !<Number of soil layers
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isinf   !<Status of infiltration model == 1
    TYPE(soilstatetype),INTENT(INOUT) :: soilstate  !<Soil states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    iarr = 0
    arrshift = iarrfirst - 1 
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%water(k,j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        DO k = 1,nsl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%temp(k,j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nc
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%deeptemp(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,nsl
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%conc(l,k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n*nc*nsl*nsfil
    ENDIF
    IF(ns>0)THEN
      IF(isN)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%humusN(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%fastN(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      IF(isP)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%humusP(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%fastP(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%partP(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%PPrelpool(j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO
      ENDIF
      IF(isC)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%humusC(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%fastC(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%oldgrw(j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO
      ENDIF
      IF(isS)THEN
        DO i = 1,n
          DO j = 1,nc
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%Srelpool(j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO
      ENDIF
      IF(isT1)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,nsl
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%partT1(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ELSEIF(nsfil>0)THEN
      IF(isN) iarr = iarr + 2*n*nc*nsl
      IF(isP) iarr = iarr + 3*n*nc*nsl + n*nc
      IF(isC) iarr = iarr + 2*n*nc*nsl + n*nc
      IF(isT1) iarr = iarr + n*nc*nsl
    ENDIF
    IF(isinf)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) soilstate%icelens(j,i) = INT(array(iarr-arrshift))
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_soilstate_variables_from_array

  !>Get number of aquifer states
  !-------------------------------------------------
  SUBROUTINE get_aquiferstate_variables_arraysize(na,ns,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 3*na        !water,lastrecharge,nextoutflow
    IF(ns>0)THEN
      dim = dim + 3*ns*na  !conc,clastrecharge,cnextoutflow
    ENDIF

  END SUBROUTINE get_aquiferstate_variables_arraysize

  !>Write aquifer state variables to array
  !-------------------------------------------------
  SUBROUTINE set_aquiferstate_variables_to_array(na,ns,aquiferstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    TYPE(aquiferstatetype),INTENT(IN) :: aquiferstate  !<aquifer states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)  !<Array of states

    INTEGER i,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%water(i)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%lastrecharge(i)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%nextoutflow(i)
    ENDDO
    IF(ns>0)THEN
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%conc(l,i)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%clastrecharge(l,i)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = aquiferstate%cnextoutflow(l,i)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_aquiferstate_variables_to_array

  !>Read aquifer state variables from array
  !-------------------------------------------------
  SUBROUTINE set_aquiferstate_variables_from_array(na,ns,nsfil,aquiferstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: na      !<Number of aquifers
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances read from file
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<aquifer states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)  !<Array of states

    INTEGER i,l
    INTEGER iarr,arrshift
    
    iarr = 0
    arrshift = iarrfirst - 1 
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%water(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%lastrecharge(i) = array(iarr-arrshift)
    ENDDO
    DO i = 1,na
      iarr = iarr + 1
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%nextoutflow(i) = array(iarr-arrshift)
    ENDDO
    IF(ns>0)THEN
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%conc(l,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%clastrecharge(l,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,na
        DO l = 1,ns
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) aquiferstate%cnextoutflow(l,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + 3*na*nsfil
    ENDIF

  END SUBROUTINE set_aquiferstate_variables_from_array

  !>Get number of river states
  !-------------------------------------------------
  SUBROUTINE get_riverstate_variables_arraysize(n,ns,nr,ml,mt,isN,isP,isC,isS,ist1,iswetl,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nr      !<number of river types
    INTEGER,INTENT(IN) :: ml      !<maximum lag steps
    INTEGER,INTENT(IN) :: mt      !<number of timesteps per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 0
    dim = dim + n*nr        !water
    dim = dim + n*nr        !temp
    dim = dim + (ml+1)*n*nr !qqueue
    dim = dim + n*nr        !Qmean
    IF(ns>0)THEN
      dim = dim + n*nr            !Psed
      dim = dim + ns*n*nr         !conc
      dim = dim + ns*(ml+1)*nr*n  !cqueue
    ENDIF
    IF(isN.OR.isP.OR.isC.OR.isS)THEN
      dim = dim + n*nr            !TPmean
      dim = dim + n*nr            !temp10
      dim = dim + n*nr            !temp20
    ENDIF
    IF(isN.OR.isP.OR.ist1)THEN
      dim = dim + mt*n*nr   !Qdayacc
      dim = dim + 366*n*nr  !Q365
    ENDIF
    IF(iswetl)THEN
      dim = dim + ns*n*nr   !cwetland
    ENDIF
    IF(isT1) dim = dim + n*nr    !T1sed
    IF(isS) dim = dim + n*nr     !sed

  END SUBROUTINE get_riverstate_variables_arraysize

  !>Write river state variables to array
  !-------------------------------------------------
  SUBROUTINE set_riverstate_variables_to_array(n,ns,nr,ml,mt,isN,isP,isC,isS,ist1,iswetl,riverstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: ml      !<maximum lag steps
    INTEGER,INTENT(IN) :: mt      !<number of timesteps per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    TYPE(riverstatetype),INTENT(IN) :: riverstate  !<River states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%water(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%temp(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        DO k = 0,ml
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%qqueue(k,j,i)
        ENDDO
      ENDDO  
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Qmean(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%conc(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 0,ml
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%cqueue(l,k,j,i)
            ENDDO
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isP)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Psed(j,i)
        ENDDO
      ENDDO
    ENDIF
    IF(isN.OR.isP.OR.isC.OR.isS)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%TPmean(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%temp10(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%temp20(j,i)
        ENDDO
      ENDDO
    ENDIF
    IF(isN.OR.isP.OR.ist1)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,mt
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Qdayacc(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,366
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Q365(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(iswetl)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%cwetland(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isT1)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%T1sed(j,i)
        ENDDO
      ENDDO  
    ENDIF  
    IF(isS)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = riverstate%Ssed(j,i)
        ENDDO
      ENDDO  
    ENDIF  

  END SUBROUTINE set_riverstate_variables_to_array

  !>Read river state variables from array
  !-------------------------------------------------
  SUBROUTINE set_riverstate_variables_from_array(n,ns,nsfil,nr,ml,mt,isN,isP,isC,isS,ist1,iswetl, &
                                                 Nfil,Pfil,Cfil,Sfil,T1fil,riverstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances of file
    INTEGER,INTENT(IN) :: nr      !<Number of river types
    INTEGER,INTENT(IN) :: ml      !<Maximum lag steps
    INTEGER,INTENT(IN) :: mt      !<number of timesteps per day
    LOGICAL,INTENT(IN) :: isN     !<Status of nitrogen simulation
    LOGICAL,INTENT(IN) :: isP     !<Status of phosphorus simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    LOGICAL,INTENT(IN) :: isS     !<Status of suspended sediment simulation
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: Nfil    !<Status of nitrogen simulation of file
    LOGICAL,INTENT(IN) :: Pfil    !<Status of phosphorus simulation of file
    LOGICAL,INTENT(IN) :: Cfil    !<Status of organic carbon simulation of file
    LOGICAL,INTENT(IN) :: Sfil    !<Status of suspended sediment simulation of file
    LOGICAL,INTENT(IN) :: T1fil   !<Status of T1 simulation of file
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim         !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)  !<Array of states

    INTEGER i,j,k,l
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1 
    iarr = 0
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%water(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%temp(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        DO k = 0,ml
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%qqueue(k,j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO  
    ENDDO
    DO i = 1,n
      DO j = 1,nr
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Qmean(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%conc(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 0,ml
            DO l = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%cqueue(l,k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n*nr*nsfil+n*nr*(ml+1)*nsfil
    ENDIF  
    IF(isP)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Psed(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(Pfil)THEN
      iarr = iarr + n*nr
    ENDIF
    IF(isN.OR.isP.OR.isC.OR.isS)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%TPmean(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%temp10(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%temp20(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(Nfil.OR.Pfil.OR.Cfil.OR.Sfil)THEN
      iarr = iarr + n*nr*3
    ENDIF  
    IF(isN.OR.isP.OR.ist1)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,mt
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Qdayacc(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,366
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Q365(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(Nfil.OR.Pfil.OR.T1fil)THEN
      iarr = iarr + n*nr*mt
      iarr = iarr + n*nr*366
    ENDIF  
    IF(iswetl)THEN
      DO i = 1,n
        DO j = 1,nr
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%cwetland(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF
    IF(ist1)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%T1sed(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO  
    ELSEIF(T1fil)THEN
      iarr = iarr + n*nr
    ENDIF  
    IF(isS)THEN
      DO i = 1,n
        DO j = 1,nr
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) riverstate%Ssed(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO  
    ELSEIF(Sfil)THEN
      iarr = iarr + n*nr
    ENDIF  

  END SUBROUTINE set_riverstate_variables_from_array

  !>Get number of lake states
  !-------------------------------------------------
  SUBROUTINE get_lakestate_variables_arraysize(n,ns,nl,isNPCS,ist2,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isNPCS   !<Status of N, P, C or S simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 0
    dim = dim + n*nl      !water
    dim = dim + n*nl      !temp
    IF(ns>0)THEN
      dim = dim + ns*n*nl !conc
      dim = dim + n*nl    !slowwater
      dim = dim + ns*n*nl !concslow
    ENDIF
    IF(isNPCS)THEN
      dim = dim + n*nl    !TPmean
      dim = dim + n*nl    !temp10
      dim = dim + n*nl    !temp20
    ENDIF
    IF(ist2)THEN
      dim = dim + n*nl    !uppertemp
      dim = dim + n*nl    !lowertemp
    ENDIF  

  END SUBROUTINE get_lakestate_variables_arraysize

  !>Write lake state variables to array
  !-------------------------------------------------
  SUBROUTINE set_lakestate_variables_to_array(n,ns,nl,isNPCS,ist2,lakestate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n             !<Number of subbasins
    INTEGER,INTENT(IN) :: ns            !<Number of substances
    INTEGER,INTENT(IN) :: nl            !<Number of lake types
    LOGICAL,INTENT(IN) :: isNPCS   !<Status of N, P, C or S simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation
    TYPE(lakestatetype),INTENT(IN) :: lakestate  !<lake states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim           !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%water(j,i)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%temp(j,i)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%conc(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%slowwater(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%concslow(k,j,i)
          ENDDO
        ENDDO  
      ENDDO
    ENDIF
    IF(isNPCS)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%TPmean(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%temp10(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%temp20(j,i)
        ENDDO
      ENDDO
    ENDIF
    IF(ist2)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%uppertemp(j,i)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = lakestate%lowertemp(j,i)
        ENDDO
      ENDDO
    ENDIF  

  END SUBROUTINE set_lakestate_variables_to_array

  !>Read lake state variables from array
  !-------------------------------------------------
  SUBROUTINE set_lakestate_variables_from_array(n,ns,nsfil,nl,isNPCS,ist2,lakestate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances of file
    INTEGER,INTENT(IN) :: nl      !<Number of lake types
    LOGICAL,INTENT(IN) :: isNPCS  !<Status of N, P or C simulation
    LOGICAL,INTENT(IN) :: ist2    !<Status of T2 simulation (of file)
    TYPE(lakestatetype),INTENT(INOUT) :: lakestate  !<lake states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%water(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    DO i = 1,n
      DO j = 1,nl
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%temp(j,i) = array(iarr-arrshift)
      ENDDO
    ENDDO
    IF(ns>0)THEN
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%conc(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%slowwater(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          DO k = 1,ns
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%concslow(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO  
      ENDDO
    ELSEIF(nsfil>0)THEN
      iarr = iarr + n*nl*nsfil
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%water(j,i) = lakestate%water(j,i) + array(iarr-arrshift)
        ENDDO
      ENDDO
      iarr = iarr + n*nl*nsfil
    ENDIF  
    IF(ns==0.AND.isNPCS)THEN
      iarr = iarr + n*nl*3
    ELSEIF(isNPCS)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%TPmean(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%temp10(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%temp20(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDIF  
    IF(ns==0.AND.ist2)THEN
      iarr = iarr + n*nl*2
    ELSEIF(ist2)THEN
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%uppertemp(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nl
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) lakestate%lowertemp(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ENDIF  

  END SUBROUTINE set_lakestate_variables_from_array

  !>Get number of miscellaneous states
  !-------------------------------------------------
  SUBROUTINE get_miscstate_variables_arraysize(n,ns,nc,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,dim)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isflood !<Status of flooded area simulation
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    INTEGER,INTENT(OUT) :: dim    !<Number of array elements

    dim = 0
    IF(iswetl)THEN
      dim = dim + n   !temp5
      dim = dim + n   !temp30
    ENDIF
    IF(isC)THEN
      dim = dim + n   !temp10
      dim = dim + n   !temp20
    ENDIF
    IF(isgsm)THEN
      dim = dim + 2*n*nc   !gdd
      dim = dim + 2*n*nc   !gsbegin
    ENDIF
    IF(isirr)THEN
      dim = dim + n*nc   !nextirrigation
      IF(ns>0)   dim = dim + n*nc*ns  !cnextirrigation
    ENDIF
    IF(isqar)THEN
      dim = dim + n   !updatestationsarcorr
    ENDIF  
    IF(isflood)THEN
      dim = dim + 2*n   !floodwater
      IF(ns>0)   dim = dim + n*2*ns  !cfloodwater
    ENDIF  
    IF(isT1)THEN
      dim = dim + nc*n       !partT1sf
    ENDIF
    
  END SUBROUTINE get_miscstate_variables_arraysize

  !>Write miscellaneous state variables to array
  !-------------------------------------------------
  SUBROUTINE set_miscstate_variables_to_array(n,ns,nc,ist1,isirr,isqar,isflood,isgsm,iswetl,isC,miscstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isflood !<Status of flooded area simulation
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    TYPE(miscstatetype),INTENT(IN) :: miscstate  !<misc states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(OUT)   :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    IF(iswetl)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp5(i)
      ENDDO
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp30(i)
      ENDDO
    ENDIF
    IF(isC)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp10(i)
      ENDDO
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%temp20(i)
      ENDDO
    ENDIF
    IF(isgsm)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%gdd(k,j,i)
          ENDDO
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = REAL(miscstate%gsbegin(k,j,i))
          ENDDO
        ENDDO  
      ENDDO
    ENDIF  
    IF(isirr)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%nextirrigation(j,i)
        ENDDO
      ENDDO
      IF(ns>0)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%cnextirrigation(k,j,i)
            ENDDO
          ENDDO  
        ENDDO
      ENDIF
    ENDIF  
    IF(isqar)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%updatestationsarcorr(i)
      ENDDO
    ENDIF
    IF(isflood)THEN
      DO i = 1,n
        DO j = 1,2
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%floodwater(j,i)
        ENDDO
      ENDDO
      IF(ns>0)THEN
        DO i = 1,n
          DO j = 1,2
            DO k = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%cfloodwater(k,j,i)
            ENDDO
          ENDDO  
        ENDDO
      ENDIF
    ENDIF
    IF(isT1)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) array(iarr-arrshift) = miscstate%partT1sf(j,i)
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE set_miscstate_variables_to_array

  !>Read misc state variables from array
  !-------------------------------------------------
  SUBROUTINE set_miscstate_variables_from_array(n,ns,nsfil,nc,readT1,ist1,isirr,isqar, &
                      isqar2,isflood,isgsm,iswetl,readC,isC,miscstate,iarrfirst,iarrlast,dim,array)

    !Argument declarations
    INTEGER,INTENT(IN) :: n       !<Number of subbasins
    INTEGER,INTENT(IN) :: ns      !<Number of substances
    INTEGER,INTENT(IN) :: nsfil   !<Number of substances of file
    INTEGER,INTENT(IN) :: nc      !<Number of classes
    LOGICAL,INTENT(IN) :: readT1  !<Status of T1 simulation for file
    LOGICAL,INTENT(IN) :: ist1    !<Status of T1 simulation
    LOGICAL,INTENT(IN) :: isirr   !<Status of irrigation simulation
    LOGICAL,INTENT(IN) :: isqar   !<Status of updating q with AR
    LOGICAL,INTENT(IN) :: isqar2  !<Status of updating q with AR for file
    LOGICAL,INTENT(IN) :: isflood  !<Status of flooded area simulation
    LOGICAL,INTENT(IN) :: isgsm   !<Status of growth season model
    LOGICAL,INTENT(IN) :: iswetl  !<Status of wetland simulation
    LOGICAL,INTENT(IN) :: readC   !<Status of organic carbon simulation for file
    LOGICAL,INTENT(IN) :: isC     !<Status of organic carbon simulation
    TYPE(miscstatetype),INTENT(INOUT) :: miscstate  !<misc states
    INTEGER,INTENT(IN) :: iarrfirst   !<Index of first array element
    INTEGER,INTENT(IN) :: iarrlast    !<Index of last array element
    INTEGER,INTENT(IN) :: dim     !<Number of array elements
    REAL,INTENT(IN)    :: array(dim)    !<Array of states

    INTEGER i,j,k
    INTEGER iarr,arrshift
    
    arrshift = iarrfirst - 1
    iarr = 0
    IF(iswetl)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp5(i) = array(iarr-arrshift)
      ENDDO
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp30(i) = array(iarr-arrshift)
      ENDDO
    ENDIF
    IF(isC)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp10(i) = array(iarr-arrshift)
      ENDDO
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%temp20(i) = array(iarr-arrshift)
      ENDDO
    ELSEIF(ns==0 .AND. readC)THEN
      iarr = iarr + 2*n
    ENDIF
    IF(isgsm)THEN
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%gdd(k,j,i) = array(iarr-arrshift)
          ENDDO
        ENDDO
      ENDDO
      DO i = 1,n
        DO j = 1,nc
          DO k = 1,2
            iarr = iarr + 1
            IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%gsbegin(k,j,i) = INT(array(iarr-arrshift))
          ENDDO
        ENDDO  
      ENDDO
    ENDIF
    IF(isirr)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%nextirrigation(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      IF(ns>0)THEN
        DO i = 1,n
          DO j = 1,nc
            DO k = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%cnextirrigation(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO  
        ENDDO
      ELSEIF(nsfil>0)THEN
        iarr = iarr + n*nc*nsfil  
      ENDIF
    ENDIF
    IF(isqar.AND.isqar2)THEN
      DO i = 1,n
        iarr = iarr + 1
        IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%updatestationsarcorr(i) = array(iarr-arrshift)
      ENDDO
    ELSEIF(isqar2)THEN
      iarr = iarr + n
    ELSEIF(isqar)THEN
      IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%updatestationsarcorr = 0.    !no error from previous simulation
    ENDIF    
    IF(isflood)THEN
      DO i = 1,n
        DO j = 1,2
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%floodwater(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
      IF(ns>0)THEN
        DO i = 1,n
          DO j = 1,2
            DO k = 1,ns
              iarr = iarr + 1
              IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%cfloodwater(k,j,i) = array(iarr-arrshift)
            ENDDO
          ENDDO  
        ENDDO
      ELSEIF(nsfil>0)THEN
        iarr = iarr + n*2*nsfil  
      ENDIF
    ENDIF
    IF(isT1)THEN
      DO i = 1,n
        DO j = 1,nc
          iarr = iarr + 1
          IF(iarr>=iarrfirst.AND.iarr<=iarrlast) miscstate%partT1sf(j,i) = array(iarr-arrshift)
        ENDDO
      ENDDO
    ELSEIF(ns==0 .AND. readT1)THEN
      iarr = n*nc
    ENDIF

  END SUBROUTINE set_miscstate_variables_from_array

END MODULE STATETYPE_MODULE
