!> \file assimilation_interface.f90
!> Contains module assimilation_interface, with model specific interface between model and assimilation routines.
!  Authors: D.Gustafsson, C. Pers, R. Piementa (SMHI)
!
!                          ^. .^
!  HYPE DATA ASSIMILATION  ( Y )
!                            `
  
!> Module with HYPE/HYSS specific subroutines for interface between HYPE/HYSS and the EnKF Data Assimilation routines.
MODULE ASSIMILATION_INTERFACE
!Copyright 2016-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------
!Procedures in this module
!-----------------------------------------------------------------------------------------
! 1) interface to the model variables 
!
! 2) code for how to read/write between Ensemble data and Model data (same comment) 
!
! The interface routines are called from MAIN (initialization, ensemble loop within time loop, enkf_analysis, output)
! 
! PLEASE NOTE: The interface routines MUST ALWAYS be updated to the current model version.
!              ALL state variables (MODVARS; HYPEVARS, etc) HAVE TO BE mirrored in the ensemble data structures.
!              On-going development to base the interface on fortran pointers to the hype state variables.
!
! To-do-list (161026):
!           * update INTERFACE to HYPE 4.12.+ (on-going)
!           * clean code from unused stuff    (on-going)     
!           * write ensemble data to binary file instead of keeping in memory for large applications. (on-going, see assimilation_variables)
!           * more clever handling of interface to the model data and to the observations (on-going with Lotta, using pointers)
!              also make sure more than Tobs and Pobs can be used as forcing!!
!           * only save ensemble of OUTVARS for the selected subbasins and selected outputs (on-going with Lotta)
!           * remove dependancy on Intel MKL library (David)
!           * add more DA filters to assimilation_routines (Rafael)
!           * generalize to enable assimilation of Aquiferstates and parameters as well (ie states on other spatial representation than sub-basins)
!               
!
! Version history:
!
! Version 8 (on-going):        Current development version
! Version 7 (20151117):        Version adapted to the HOPE model
! Version 6 (20141210):        Ongoing update to HYPE 4.9.3 and futher implemenation of spatially correlated ensemble generation:
! Version 5 (20140508):        Updates to HYPE verion 4.5.0, including the 2D random field generator (not yet finalized)
! version 4 (131203):          Renamed to the general name assim_interface.f90, to open up for alterntive DA methods other than EnKF
! version 3 (130925):          This is the third version, implemented for HYPE version 4.5.0, accounting for flipped data matrices in HYPE
! version 2 (120510):          This is the second version, application to HYPE version 3.6.2 extended with assimilation of lake and river ice variables
! versions 0-1:(ending 120417) the first version (never in svn), application to the model HYPE version 3.5.0 for snow data assimilation
!
! Questions: David Gustafsson, david.gustafsson@smhi.se, +46(0)114958647
!-----------------------------------------------------------------------------------------
! The subroutines are divided into three blocks:
!     INITIALIZATION
!     READ/WRITE BETWEEN MODEL AND ENSEMBLE
!     ENKF ANALYSIS AND ENSEMBLE GENERATION
!-----------------------------------------------------------------------------------------

! Data assimilation modules
USE ASSIMILATION_VARIABLES
USE ASSIMILATION_ROUTINES

! HYPE/HYSS modules
USE MODVAR
USE WORLDVAR
USE HYPEVARIABLES
USE STATETYPE_MODULE
USE COMPOUT, ONLY: find_variable_index_type
USE READWRITE_ROUTINES, ONLY: read_next_codestr_on_line
USE DATAMODULE, ONLY : set_outvar_for_variable

IMPLICIT NONE


  CONTAINS

!-----------------------------------------------------------------------------------------------
! FIRST BLOCK: ROUTINES FOR INITIALIZATION
!-----------------------------------------------------------------------------------------------

  !----------------------------------------------------------------------------------------------------
  !>Defines the categories for assimilating groups of (state) variables for HYPE
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE initialize_assim_categories_HYPE(assimInfo,stateinfo)

    !Argument declarations
    TYPE(assim_info_type),INTENT(INOUT) :: assimInfo !<information on assimilation simulation
    TYPE(stateinfotype),INTENT(IN)   :: stateinfo(:) !<information on HYPE state variables
    
    !Allocate and set size of array
    assimInfo%nCat = SIZE(stateinfo)
    IF(ALLOCATED(assimInfo%assim_categories)) DEALLOCATE(assimInfo%assim_categories)
    ALLOCATE(assimInfo%assim_categories(assimInfo%nCat))
    
    !The category definition of HYPE model state variables
    !assimInfo%assim_categories(1:6)   = 'snow'   !CP710614 updated to current HYPE
    !assimInfo%assim_categories(7)     = 'glacier'
    !assimInfo%assim_categories(8:12)  = 'riverice'
    !assimInfo%assim_categories(13:17) = 'lakeice'
    !assimInfo%assim_categories(18:30) = 'soil'
    !assimInfo%assim_categories(31:36) = 'aquifer'
    !assimInfo%assim_categories(37:49) = 'riverwt'
    !assimInfo%assim_categories(50:59) = 'lakewt'
    !assimInfo%assim_categories(60:70) = 'misc'
    assimInfo%assim_categories(1:6)   = 'snow'
    assimInfo%assim_categories(7)     = 'glacier'
    assimInfo%assim_categories(8:13)  = 'riverice'
    assimInfo%assim_categories(14:19) = 'lakeice'
    assimInfo%assim_categories(20:35) = 'soil'
    assimInfo%assim_categories(36:41) = 'aquifer'
    assimInfo%assim_categories(42:56) = 'riverwt'
    assimInfo%assim_categories(57:66) = 'lakewt'
    assimInfo%assim_categories(67:78) = 'misc'
    
  END SUBROUTINE initialize_assim_categories_HYPE

  !----------------------------------------------------------------------------------------------------
  !>Get name of state type for a given state variable in a category
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE get_statetype_from_category(assimInfo,stateinfo,catname,svname,stname)
  
    !Argument declaration
    TYPE(assim_info_type),INTENT(IN) :: assimInfo !<information on assimilation simulation
    TYPE(stateinfotype),INTENT(IN)   :: stateinfo(:) !<information on HYPE state variables
    CHARACTER(LEN=*), INTENT(IN)     :: catname   !<category name
    CHARACTER(LEN=*), INTENT(IN)     :: svname    !<state variable name
    CHARACTER(LEN=*), INTENT(OUT)    :: stname    !<state type name
    
    !Local variables
    INTEGER i
    
    stname = ''
    DO i = 1,assimInfo%nCat
      IF(assimInfo%assim_categories(i)==catname .AND.  &
         stateinfo(i)%svname == svname)THEN
        stname = stateinfo(i)%stname
        EXIT
      ENDIF
    ENDDO
    IF(stname=='')THEN
      WRITE(6,*) 'ERROR: statetype not found for category: ',TRIM(catname),', variable: ',TRIM(svname)
      STOP 1
    ENDIF

  END SUBROUTINE get_statetype_from_category

  !----------------------------------------------------------------------------------------------------
  !>Find out if the given state variable is to be assimilated
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE set_assimilate(assimInfo,stateinfo,ivar,astatus)
  
    !Argument declaration
    TYPE(assim_info_type),INTENT(IN) :: assimInfo !<information on assimilation simulation
    TYPE(stateinfotype),INTENT(IN)   :: stateinfo(:) !<information on HYPE state variables
    INTEGER, INTENT(IN)              :: ivar      !<index of current state variable
    LOGICAL, INTENT(OUT)             :: astatus   !<flag if assimilate this variable
    
    !Local variables
    INTEGER i
    CHARACTER(LEN=20) :: catname   !category name
    
    astatus = .TRUE.   !default
    catname = assimInfo%assim_categories(ivar)
    !1) check category assimilation settings
    CALL get_assimilate_for_category(assimInfo,catname,astatus)
    !2) check state variable assimilation settings (specific variable setting has highest priority)
    DO i = 1,assimInfo%nFlag
      IF(assimInfo%assim_flag(i)%category==stateinfo(ivar)%stname .AND. &  
         assimInfo%assim_flag(i)%varname==stateinfo(ivar)%svname)THEN
        astatus = assimInfo%assim_flag(i)%add
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE set_assimilate

  !----------------------------------------------------------------------------------------------------
  !>Find out if the given category is to be assimilated
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE get_assimilate_for_category(assimInfo,catname,astatus)
  
    !Argument declaration
    TYPE(assim_info_type),INTENT(IN) :: assimInfo !<information on assimilation simulation
    CHARACTER(LEN=*),INTENT(IN)      :: catname   !<category name
    LOGICAL, INTENT(OUT)             :: astatus   !<flag if assimilate this category
    
    !Local variables
    INTEGER i
    
    astatus = .TRUE.   !default
    DO i = 1,assimInfo%nFlag
      !1) check category assimilation settings
      IF(assimInfo%assim_flag(i)%category==catname .AND. & 
         assimInfo%assim_flag(i)%varname=='')THEN
        astatus = assimInfo%assim_flag(i)%add
        EXIT
      ENDIF
    ENDDO

  END SUBROUTINE get_assimilate_for_category

  !----------------------------------------------------------------------------------------------------
  !>A subroutine to count the number of states variables in HYPE + some other HYPEVARIABLES
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE checkModelStructure_States_HYPE2(nx,stateinfo)

    !Argument declarations
    INTEGER, INTENT(OUT) :: nx                     !<number of state variables (each variable is size nsubbasin but aquifer states)
    TYPE(stateinfotype),INTENT(IN) :: stateinfo(:) !<information on HYPE state variables
    
    !Local variables
    INTEGER i,nstate
      
    !1) check how many state ensemble matrices we need (each state ensemble matrix is numsubbasins x ensemble_size, so there will be nclass matrices for the standard HYPE state (one state per class))
    nx=0
    nstate = SIZE(stateinfo)
    DO i = 1,nstate
      IF(stateinfo(i)%allok)THEN
        IF(stateinfo(i)%ndim==1)THEN
          nx = nx + 1
        ELSEIF(stateinfo(i)%ndim==2)THEN
          nx = nx + stateinfo(i)%dims(1)
        ELSEIF(stateinfo(i)%ndim==3)THEN
          nx = nx + stateinfo(i)%dims(1)*stateinfo(i)%dims(2)
        ELSEIF(stateinfo(i)%ndim==4)THEN
          nx = nx + stateinfo(i)%dims(1)*stateinfo(i)%dims(2)*stateinfo(i)%dims(3)
        ENDIF
      ENDIF
    ENDDO
    !Note, if aquifer is simulated not all state ensemble matrices will be of dimension nsubbasin. They will be naquifers.
      
    !2) Add variables for river flows (module variables in HYPEVARIABLES)
    IF(ALLOCATED(Qmax))    nx = nx + 2
    IF(ALLOCATED(Q2max))   nx = nx + 2
    IF(ALLOCATED(iQmax))   nx = nx + 2
    IF(ALLOCATED(iQ2max))  nx = nx + 2
    IF(ALLOCATED(accdiff)) nx = nx + 1

    !How about the load variables for source apportionment (also in hypevar)? 
    !Do we want to be able to have assimilation adn load output at the same time?
    
  END SUBROUTINE checkModelStructure_States_HYPE2

  !----------------------------------------------------------------------------------------------------
	!> Allocate and initialize ensembles for HYPE model states
	!----------------------------------------------------------------------------------------------------
  SUBROUTINE allocate_and_initialize_model_state_ensembles_HYPE2(ne,assimVar,varID,assimInfo,stateinfo)
    
    !Argument declarations
    INTEGER,INTENT(IN)                :: ne           !<ensemble size                               
    TYPE(assim_state_ensemble_type)   :: assimVar(:)  !<vector of state ensemble data  
    INTEGER, INTENT(INOUT)            :: varID        !<Ensemble variable ID (from 1 to number of state ensembles)
    TYPE(assim_info_type),INTENT(IN)  :: assimInfo    !<information, settings, etc for assimilation module
    TYPE(stateinfotype),INTENT(IN)    :: stateinfo(:) !<Information about state variables
   
    !Local variables
    LOGICAL :: assimilate
    INTEGER :: locID, coordID
    INTEGER :: ivar, nstates,reclen,dimmax,xdim,idim
    REAL,ALLOCATABLE :: xini(:)
    INTEGER,ALLOCATABLE :: fileIDarray(:)
    
    nstates = SIZE(stateinfo)
    ALLOCATE(fileIDarray(assimInfo%nBinFiles))
    IF(assimInfo%useBinFilesX==1)THEN
      !If bin-files is used open file for state ensembles
     ! DO ivar=1,nstates
     !   IF(stateinfo(ivar)%allok)THEN
     !     IF(stateinfo(ivar)%ndim==1)THEN
     !       ALLOCATE(xini(stateinfo(ivar)%dims(1)))   !OBS !Assume all have nsub as last dimension here
     !       EXIT
     !     ENDIF
     !   ENDIF
     ! ENDDO
      !loop through states and identify the largest calculation unit (either subbasins or aquifers)
      dimmax=0
      DO ivar=1,nstates
        IF(stateinfo(ivar)%allok)THEN
          !find the last non-zero dimension and compare to the current maximum
          xdim=stateinfo(ivar)%dims(1)
          DO idim=2,4
            IF(stateinfo(ivar)%dims(idim).GT.0)THEN
              xdim=stateinfo(ivar)%dims(idim)
            ENDIF
            IF(xdim.GT.dimmax)dimmax=xdim
          ENDDO
        ENDIF  
      ENDDO
      ALLOCATE(xini(dimmax))
      INQUIRE(IOLENGTH=reclen) xini   !determine suitable record length here, because ifort och gfortran have different file storage unit (i.e. RECL)
      fileIDarray = fid_assim_bin_base
      OPEN(UNIT=fid_assim_bin_base,FILE='ensXstates.bin',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=reclen)
      DEALLOCATE(xini)
    ELSE
      !Set local array with fileIds that will be used in the following subroutine calls.
      fileIDarray = fid_assim_bin
    ENDIF

    !3) allocate and initialize each variable ensemble
    !locID=1
    !coordID=1    !CP161213 changed to include aquifers
    !nb.1, varID is initialized to 1 in the calling procedure,
    !nb.2, varID is increased by 1-xx for each call to the "allocate_..." routines below, depending on dimension of the HYPE states.  
    !nb.3, locID and coordID are fixed to 1 here (until we solve how to use more than one of these...)
    !nb.4, coordID is set to 2 for aquiferstates and 1 for all others, locID the same
   
    DO ivar=1,nstates
      !1) Check if the state should be changed by assimilation
      IF(stateinfo(ivar)%allok)THEN
        CALL set_assimilate(assimInfo,stateinfo,ivar,assimilate)
        !2) Allocate an state ensemble for each state
        IF(stateinfo(ivar)%ndim==1)THEN
!          CALL get_HYPEstate_coordID(TRIM(stateinfo(ivar)%stname),coordID) !CP161213 added to include aquifers
          CALL get_HYPEstate_coordID(stateinfo(ivar)%stname,coordID) !CP161213 added to include aquifers
          locID = coordID
          CALL allocate_0dim_state_ensemble(assimVar(varID),ne, &
                stateinfo(ivar)%dims(1),   &
                TRIM(stateinfo(ivar)%stname)//'_'//TRIM(stateinfo(ivar)%svname), &
                stateinfo(ivar)%svpoint%d1,varID,locID,coordID, &
                fileIDarray(varID),assimInfo%useBinFilesX,0.,assim_defmax,assimilate,assimInfo%initializeFromBinFiles)
        ELSEIF(stateinfo(ivar)%ndim==2)THEN
!          CALL get_HYPEstate_coordID(TRIM(stateinfo(ivar)%stname),coordID) !CP161213 added to include aquifers
          CALL get_HYPEstate_coordID(stateinfo(ivar)%stname,coordID) !CP161213 added to include aquifers
          locID = coordID
          IF(ASSOCIATED(stateinfo(ivar)%svpoint%d2))THEN          !both real and integer state of 3dim exist
            CALL allocate_1dim_state_ensemble(assimVar(:),ne,     &
                stateinfo(ivar)%dims(2),   &
                TRIM(stateinfo(ivar)%stname)//'_'//TRIM(stateinfo(ivar)%svname), &
                stateinfo(ivar)%svpoint%d2,varID,locID,coordID, &
                fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,                  &
                stateinfo(ivar)%dims(1),assimilate,assimInfo%initializeFromBinFiles)
          ELSE
            CALL allocate_1dim_state_ensemble(assimVar(:),ne,     &
                stateinfo(ivar)%dims(2),   &
                TRIM(stateinfo(ivar)%stname)//'_'//TRIM(stateinfo(ivar)%svname), &
                REAL(stateinfo(ivar)%svpoint%d2i),varID,locID,coordID, &
                fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,                  &
                stateinfo(ivar)%dims(1),assimilate,assimInfo%initializeFromBinFiles)
          ENDIF
        ELSEIF(stateinfo(ivar)%ndim==3)THEN
!          CALL get_HYPEstate_coordID(TRIM(stateinfo(ivar)%stname),coordID) !CP161213 added to include aquifers
          CALL get_HYPEstate_coordID(stateinfo(ivar)%stname,coordID) !CP161213 added to include aquifers
          locID = coordID
          IF(ASSOCIATED(stateinfo(ivar)%svpoint%d3))THEN          !both real and integer state of 3dim exist
            CALL allocate_2dim_state_ensemble(assimVar(:),ne,   &
                stateinfo(ivar)%dims(3),   &
                TRIM(stateinfo(ivar)%stname)//'_'//TRIM(stateinfo(ivar)%svname), &
                stateinfo(ivar)%svpoint%d3,varID,locID,coordID, &
                fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,                  &
                stateinfo(ivar)%dims(2),stateinfo(ivar)%dims(1),&
                assimilate,assimInfo%initializeFromBinFiles)
          ELSE
            CALL allocate_2dim_state_ensemble(assimVar(:),ne,   &
                stateinfo(ivar)%dims(3),   &
                TRIM(stateinfo(ivar)%stname)//'_'//TRIM(stateinfo(ivar)%svname), &
                REAL(stateinfo(ivar)%svpoint%d3i),varID,locID,  &
                coordID,fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,          &
                stateinfo(ivar)%dims(2),stateinfo(ivar)%dims(1),&
                assimilate,assimInfo%initializeFromBinFiles)
          ENDIF
        ELSEIF(stateinfo(ivar)%ndim==4)THEN
!          CALL get_HYPEstate_coordID(TRIM(stateinfo(ivar)%stname),coordID) !CP161213 added to include aquifers
          CALL get_HYPEstate_coordID(stateinfo(ivar)%stname,coordID) !CP161213 added to include aquifers
          locID = coordID
          CALL allocate_3dim_state_ensemble(assimVar(:),ne,     &
                stateinfo(ivar)%dims(4),   &
                TRIM(stateinfo(ivar)%stname)//'_'//TRIM(stateinfo(ivar)%svname), &
                stateinfo(ivar)%svpoint%d4,varID,locID,coordID, &
                fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,                  &
                stateinfo(ivar)%dims(3),stateinfo(ivar)%dims(2),&
                stateinfo(ivar)%dims(1),assimilate,assimInfo%initializeFromBinFiles)
        ENDIF
      ENDIF
    ENDDO

    !3) Variables for river flows (module variables in HYPEVARIABLES)
    CALL get_assimilate_for_category(assimInfo,'riverwt',assimilate)
    coordID = 1 !CP161213 added to include aquifers
    locID = coordID
    IF(ALLOCATED(Qmax))    CALL allocate_1dim_state_ensemble(assimVar(:),ne,nsub,"Qmax",Qmax,varID,locID,coordID,fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,2,assimilate,assimInfo%initializeFromBinFiles)
    IF(ALLOCATED(Q2max))   CALL allocate_1dim_state_ensemble(assimVar(:),ne,nsub,"Q2max",Q2max,varID,locID,coordID,fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,2,assimilate,assimInfo%initializeFromBinFiles)
    IF(ALLOCATED(iQmax))   CALL allocate_1dim_state_ensemble(assimVar(:),ne,nsub,"iQmax",real(iQ2max),varID,locID,coordID,fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,2,assimilate,assimInfo%initializeFromBinFiles)
    IF(ALLOCATED(iQ2max))  CALL allocate_1dim_state_ensemble(assimVar(:),ne,nsub,"iQ2max",real(iQ2max),varID,locID,coordID,fileIDarray,assimInfo%useBinFilesX,0.,assim_defmax,2,assimilate,assimInfo%initializeFromBinFiles)
    IF(ALLOCATED(accdiff)) CALL allocate_0dim_state_ensemble(assimVar(varID),ne,nsub,"accdiff",accdiff,varID,locID,coordID,fileIDarray(varID),assimInfo%useBinFilesX,0.,assim_defmax,assimilate,assimInfo%initializeFromBinFiles)
 
    DEALLOCATE(fileIDarray)

  END SUBROUTINE allocate_and_initialize_model_state_ensembles_HYPE2

  !----------------------------------------------------------------------------------------------------
  !> Get coordinate ID for states based on number of subbasins or aquifers
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE get_HYPEstate_coordID(stname,coordID)
    !Argument declarations
    CHARACTER(LEN=20),INTENT(IN) :: stname  !<state type name
    INTEGER,INTENT(OUT) :: coordID  !<ID for coordinate system for this
    
    IF(stname(1:12)=='aquiferstate')THEN
      coordID = 2   !aquifer state, coordinates defined for aquifers, possibly introduce parameters here..
    ELSE
      coordID = 1   !coordinates defined for subbasins
    ENDIF
  END SUBROUTINE get_HYPEstate_coordID
    
  !----------------------------------------------------------------------------------------------------
  !> \brief Calculate (middle) coordinates for aquifers.
  !>
  !> The coordinates are calculated as area weighted subbasin coordinates of 
  !> subbasin recieving or recharging to/from the aquifer.
  !----------------------------------------------------------------------------------------------------
  SUBROUTINE calculate_coordinates_for_aquifer(ia,x,y)
  
    !Argument declarations
    INTEGER,INTENT(IN) :: ia  !<index of aquifer
    REAL,INTENT(OUT)   :: x   !<x-coordinate of aquifer
    REAL,INTENT(OUT)   :: y   !<y-coordinate of aquifer
    
    !Local variables
    INTEGER i
    REAL sumx,sumy,suma
    
    sumx=0.;sumy=0.;suma=0.
    DO i=1,nsub
      IF(path(i)%aquid==ia)THEN
        sumx = sumx + basin(i)%xcoord * basin(i)%area
        sumy = sumy + basin(i)%ycoord * basin(i)%area
        suma = suma + basin(i)%area
      ENDIF
    ENDDO
    IF(suma>0.)THEN
      x=sumx/suma
      y=sumy/suma
    ELSE
      WRITE(6,*) 'Warning: Calculating aquifer area.'
      WRITE(6,*) 'Warning: Aquifer with index',ia,'does not include subbasins with area.'
      WRITE(6,*) 'Warning: Check indata.'
    ENDIF
  END SUBROUTINE calculate_coordinates_for_aquifer

  !----------------------------------------------------------------------------------------------------
	!> Allocate and initialize ensembles for HYPE model forcing
  !
	!allocate_and_initialize_model_forcing_ensembles_HYPE   ^_ _^
	!    one more subroutine with self-explaining name      ( Y )
  !                                                          `
  !
  ! todo (20161119): generalize - maybe also forcing variables need a pointer based info-variable, however
  !                  for now we can use the hard-coded stuff, just add more variables if need (TMAX, TMIN etc)
  !                  Btw, the difference TMAX-TMIN is probably a better variable for perturbation, 1) we only need
  !                  to perturb one variable, and 2) we will never end up with TMIN>TMAX, and 3) we mainly use the 
  !                  difference anyway (except for Priestly-Taylor and Penman-Montieth, when we use TMIN for estimating
  !                  actual vapour pressure).
	!---------------------------------------------------------------------------------------------------
  !SUBROUTINE allocate_and_initialize_model_forcing_ensembles_HYPE(ne,nsub,assimVar,varID,nInp,ITYPE,IID,IMIN,IMAX,ISIGMA,ISEMIM,IRESTM,ILSCALE,IGRIDSIZE,ICORRTYPE,XC,YC)    !CP use bin-files
  SUBROUTINE allocate_and_initialize_model_forcing_ensembles_HYPE(ne,nsub,assimVar,varID,assimInfo,nInp,ITYPE,IID,IMIN,IMAX,ISIGMA,ISEMIM,IRESTM,ILSCALE,IGRIDSIZE,ICORRTYPE,XC,YC)
    
    !inputs
    INTEGER :: ne                     !<ensemble size
    INTEGER :: nsub                   !<number of sub-basins
    INTEGER :: varID                  !<ensemble ID
    INTEGER :: nInp                   !<number of forcing variables
    TYPE(assim_input_ensemble_type)  :: assimVar(:) !<ensemble data vector
    TYPE(assim_info_type),INTENT(IN)  :: assimInfo  !<information, settings, etc for assimilation module
    INTEGER :: ITYPE(:)                             !<ensemble generation type
    REAL    :: IMIN(:),IMAX(:),ISIGMA(:)            ! variable minimum, maximum, standard deviation
    INTEGER :: IID(:)                               !<variable ID
    REAL    :: ISEMIM(:),IRESTM(:)                  ! relative sd for semi-restricted, and restricted variables
    REAL    :: ILSCALE(:),IGRIDSIZE(:)              ! correlation length scale, and grid-size for 2D random number generation
    INTEGER :: ICORRTYPE(:)                         !<correlation function for 2D random fields
    REAL XC(:),YC(:)                                ! x,y coordinates
    
    ! locals
    INTEGER i, k                                    ! loop indices 
    real, allocatable :: xini(:)                    ! default initial state vector
    INTEGER,ALLOCATABLE :: fileIDarray(:)           ! local array for setting fileID of all forcing ensembles
    INTEGER reclen
    
    ! start of calculations
    ALLOCATE(xini(nsub))
    xini = 0.
    k=0

    !Open bin-file if separete files are used (useBinFiles==1)
    ALLOCATE(fileIDarray(assimInfo%nBinFiles))
    IF(assimInfo%useBinFilesFA==1)THEN
      !If one bin-files is used open file for forcing ensembles
      INQUIRE(IOLENGTH=reclen) xini   !determine suitable record length here, because ifort och gfortran have different file storage unit (i.e. RECL)
      fileIDarray = fid_assim_bin(varID)    !use first F-variable's fid
      OPEN(UNIT=fid_assim_bin(varID),FILE='ensFstates.bin',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=reclen)
    ELSE
      !Set local array with fileIds that will be used in the following subroutine calls.
      fileIDarray = fid_assim_bin
    ENDIF

    ! loop over forcing variables adn allocate ensemble data
    DO i=1,nInp
      IF(ITYPE(i).GE.1)THEN
        k=k+1
!! TODO: generalize to enable any type of forcing from the model variables, not only Pobs and Tobs (David 20161026)
        IF(IID(i).EQ.1)THEN     !Precipitation
          CALL allocate_assim_forcing_ensemble(assimVar(k),ne,nsub,"Pobs",xini,1,1,1,fileIDarray(varID),IMIN(i),IMAX(i), &
                                               ITYPE(i),ISIGMA(i),ISEMIM(i),IRESTM(i),assim_minsigma,ILSCALE(i),IGRIDSIZE(i),ICORRTYPE(i),XC,YC,assimInfo%useBinFilesFA)
          varID=varID+1
        ELSEIF(IID(i).EQ.2)THEN !Temperature Mean
            CALL allocate_assim_forcing_ensemble(assimVar(k),ne,nsub,"Tobs",xini,2,1,1,fileIDarray(varID),IMIN(i),IMAX(i), &
                                               ITYPE(i),ISIGMA(i),ISEMIM(i),IRESTM(i),assim_minsigma,ILSCALE(i),IGRIDSIZE(i),ICORRTYPE(i),XC,YC,assimInfo%useBinFilesFA)
            varID=varID+1
        ELSEIF(IID(i).EQ.3)THEN !Temperature Min
            CALL allocate_assim_forcing_ensemble(assimVar(k),ne,nsub,"TMINobs",xini,3,1,1,fileIDarray(varID),IMIN(i),IMAX(i), &
                                               ITYPE(i),ISIGMA(i),ISEMIM(i),IRESTM(i),assim_minsigma,ILSCALE(i),IGRIDSIZE(i),ICORRTYPE(i),XC,YC,assimInfo%useBinFilesFA)
            varID=varID+1
        ELSEIF(IID(i).EQ.4)THEN !Temperature Max
            CALL allocate_assim_forcing_ensemble(assimVar(k),ne,nsub,"TMAXobs",xini,4,1,1,fileIDarray(varID),IMIN(i),IMAX(i), &
                                               ITYPE(i),ISIGMA(i),ISEMIM(i),IRESTM(i),assim_minsigma,ILSCALE(i),IGRIDSIZE(i),ICORRTYPE(i),XC,YC,assimInfo%useBinFilesFA)
            varID=varID+1
        ENDIF    
      ENDIF
    ENDDO
    DEALLOCATE(xini)
    DEALLOCATE(fileIDarray)

  END SUBROUTINE allocate_and_initialize_model_forcing_ensembles_HYPE

  !----------------------------------------------------------------------------------------------------
	!> Allocate and initialize ensembles for HYPE model auxiliary variables, i.e. outvar
  !
  !allocate_and_initialize_model_auxiliary_ensembles_HYPE    ^, ,^
  !    yet another subroutine with self-explaining name      ( Y )   this time for auxiliaries, aka outvars
  !                                                             °
  !
  ! todo (20161119): this should also be revised using the new structures from Lotta to only include
  !                  (automatically, settings in assim_info not needed) 1) outvars requested for output
  !                  in info.txt, 2) outvars needed for criteria (info.txt), 3) outvars used in the assimilation
  !                  (observations (from Xobs or Qobs) and modelled observations)
  !---------------------------------------------------------------------------------------------------
  SUBROUTINE allocate_and_initialize_model_auxiliary_ensembles_HYPE(ne,nsub,assimVar,assimInfo,varID,outid,AID,AMIN,AMAX)
  
    ! inputs
    INTEGER :: ne                     !<ensemble size
    INTEGER :: nsub                   !<number of sub-basins
    INTEGER :: varID                  !<ensemble ID
    INTEGER :: outid                  !<number of auxiliaries
!    INTEGER :: ne,nsub,varID,outid                    ! ensemble size, n:o subbasins, ensemble ID, n:o auxiliaries   !CP170619 added doxygen comments
    TYPE(assim_state_ensemble_type)  :: assimVar(:)   !<ensemble data
    TYPE(assim_info_type),INTENT(IN)  :: assimInfo    !<information, settings, etc for assimilation module
    REAL    :: AMIN(:),AMAX(:)                        ! minimum and maximum allowed values
    INTEGER :: AID(:)                                 !<index of outvar

    ! locals
    INTEGER i, k, reclen
    real, allocatable :: xini(:)
    INTEGER,ALLOCATABLE :: fileIDarray(:)
    
    ! start of calculations
    allocate(xini(nsub))
    xini = 0.
    k=0  !changed to let the subroutine increase k instead of in loop
    !k=1
    
    !Open bin-file if separete files are used (useBinFiles==1)
    ALLOCATE(fileIDarray(assimInfo%nBinFiles))
    IF(assimInfo%useBinFilesFA==1)THEN
      !If one bin-files is used open file for forcing ensembles
      INQUIRE(IOLENGTH=reclen) xini   !determine suitable record length here, because ifort och gfortran have different file storage unit (i.e. RECL)
      fileIDarray = fid_assim_bin(varID)    !use first A-variable's fid
      OPEN(UNIT=fid_assim_bin(varID),FILE='ensAstates.bin',FORM='UNFORMATTED',ACCESS='DIRECT',RECL=reclen)
    ELSE
      !Set local array with fileIds that will be used in the following subroutine calls.
      fileIDarray = fid_assim_bin
    ENDIF

    ! loop over auxiliaries (outvars)
    DO i=1,outid
      !IF(AID(i).GE.1)THEN   !CP161207 this is not necessary only AID>0 in outid=nA
        k=k+1  !k=i
!!TODO: enable more flexible way to set the min and max values, see example in assim_interface_HOPE.f90
        !CALL allocate_0dim_state_ensemble(assimVar(k),ne,nsub,"aux",xini,AID(i),1,1,fid_assim_bin(varID),0.,assim_defmax,.false.) !CP161201 added bin-files, set assimilate here (bug correction), AID(i) will change in subroutine (bug?) added extra variable for record and increase
        !CALL allocate_auxiliary_ensemble(assimVar(k),ne,nsub,"aux",xini,AID(i),k,1,1,fileIDarray(varID),assimInfo%useBinFilesFA,0.,assim_defmax,.true.) !CP161218 use AMIN,AMAX arguments! They are defmin and defmax
        CALL allocate_auxiliary_ensemble(assimVar(k),ne,nsub,"aux",xini,AID(i),k,1,1,fileIDarray(varID),assimInfo%useBinFilesFA,AMIN(i),AMAX(i),.true.)
        varID=varID+1  !this the same varID as used for the other ensembles (used for fileID)
      !ENDIF
    ENDDO
    deallocate(xini)
    DEALLOCATE(fileIDarray)
    
  END SUBROUTINE allocate_and_initialize_model_auxiliary_ensembles_HYPE

  !----------------------------------------------------------------------------------------------------
	!> Allocate and initialize ensembles for HYPE model observations
  !
  !allocate_and_initialize_observation_ensembles_HYPE    ^. .^
  !    one more time, now for observations               ( Y )
  !                                                         `
  !  make sure settings are read correctly in the new assim_info-reader
  !
  !  later on, also add functionality to use more than one set of coordinates for observations
  !---------------------------------------------------------------------------------------------------
  SUBROUTINE allocate_and_initialize_observation_ensembles_HYPE(ne,nsub,assimVar,varID,obsid,OIDD,OIDDOV,OIDHXOV,OTYPE,OMIN,OMAX,OSIGMA,OSEMIM,ORESTM,OLSCALE,OGRIDSIZE,OCORRTYPE,XC,YC)

  ! inputs
    INTEGER :: ne                     !<ensemble size
    INTEGER :: nsub                   !<number of sub-basins
    INTEGER :: varID                  !<ensemble ID
    INTEGER :: obsid                  !<number of observations
    !INTEGER :: varID,obsid,ne,nsub                  ! ensemble ID, n:o observations, ensemble size, n:o subbasins  !CP170619 added doxygen comments
    TYPE(assim_input_ensemble_type)  :: assimVar(:) !<ensemble data
    INTEGER :: OTYPE(:)                             !<ensemble generation type
    REAL    :: OMIN(:),OMAX(:),OSIGMA(:)            ! min, max, sd
    !INTEGER :: OIDD(:),OIDHX(:)                     ! model id, observation id (both referring to outvar id)   !CP170202 changed because of reduced size of new outvar
    INTEGER :: OIDD(:)                     !<model id, referring to index in outvarid
    INTEGER :: OIDDOV(:),OIDHXOV(:)                 ! model and observation id, referring to index in outvar
    REAL    :: OSEMIM(:),ORESTM(:)                  ! relative sd of semi-restricted and restricted observations
    REAL    :: OLSCALE(:),OGRIDSIZE(:)              ! correlation length and gridsize for 2d random number generation
    INTEGER :: OCORRTYPE(:)                         !<correlation function for 2d random number generation
    REAL    :: XC(:), YC(:)                         ! X and Y coordinates
    
    ! locals
    INTEGER i, k
    real, allocatable :: xini(:)
    ! start of calculations
    allocate(xini(nsub))
    xini = 0.
    k=0
    ! loop over observations
    DO i=1,obsid
      IF(OTYPE(i).GE.1 .AND. OIDD(i).LE.max_outvar)THEN
        k=k+1
        CALL allocate_assim_observation_ensemble(assimVar(k),ne,nsub,outvarid(OIDD(i))%shortname,xini,OIDDOV(i),OIDHXOV(i),1,fid_assim_bin(varID),OMIN(i),OMAX(i), &
                                               OTYPE(i),OSIGMA(i),OSEMIM(i),ORESTM(i),assim_minsigma,OLSCALE(i),OGRIDSIZE(i),OCORRTYPE(i),XC,YC)
          varID=varID+1
      ENDIF
    ENDDO
    deallocate(xini)
  END SUBROUTINE allocate_and_initialize_observation_ensembles_HYPE
  
  !-----------------------------------------------------------------------------------------------
  !> Initializes a assimilation simulation:
  !>
  !>\li define the basic info needed        (in the code and read from input file)
  !>\li allocate the assimilation variables (model ensemble, input ensemble, etc)
  !>\li initialize assimilation variables
  !
  ! CALLED FROM MAIN!
  !-----------------------------------------------------------------------------------------------
  SUBROUTINE assim_initialize(dir,assimData,ne,resultdir,stateinfo,noutvar,n_Result)
    !input arguments
    CHARACTER(LEN=maxcharpath),INTENT(IN) :: dir       !<path to AssimInfo.txt
    TYPE(assim_data_type), INTENT(INOUT)  :: assimData !<main assimilation variable containing all data
    INTEGER, INTENT(OUT)                  :: ne        !<Ensemble size (read from AssimInfo.txt by this routine)
    CHARACTER(LEN=200),INTENT(IN)         :: resultdir !<path to output directory   !NOT USED
    TYPE(stateinfotype), INTENT(IN)       :: stateinfo(:) !<Information about state variables

    INTEGER, INTENT(INOUT)                :: noutvar   !<number of output variables to be set (so far)  !CP added this
    INTEGER, INTENT(OUT)                  :: n_Result  !<Error status of subroutine !CP add this!

    !Local parameters
    INTEGER, PARAMETER                :: linelen = 1000      
    CHARACTER (LEN=2), PARAMETER      :: errstr = '##'  !local code for error
    INTEGER, PARAMETER                :: maxenkfpar = 300
    INTEGER, PARAMETER                :: strvlen = 20   !length of strvalue and code

    !Local variables
    CHARACTER(LEN=linelen)            :: newline
    CHARACTER (LEN=maxcharpath+8)     :: filename
    CHARACTER (LEN=strvlen)           :: code        !made longer
    CHARACTER (LEN=strvlen)           :: strvalue    !made longer and easier to change in the future
    CHARACTER (LEN=strvlen)           :: stname
    INTEGER :: ncom, ninfo, parid, obsid, inpid, nblank
    !INTEGER :: na,np,nf,nd,nax,npx,nfx,ndt
    INTEGER :: i, newint, outid
    integer :: varID
    INTEGER :: iass  !count rows with information on assimilation variables

    ! parameter info
!		REAL    :: PVAL(maxenkfpar),PMIN(maxenkfpar),PMAX(maxenkfpar),PSIGMA(maxenkfpar)
!		INTEGER :: POPTION(maxenkfpar), PID(maxenkfpar)   
!		REAL    :: PSEMIM(maxenkfpar),PRESTM(maxenkfpar)

    !observation info
    INTEGER :: OTYPE(maxenkfpar),OIDD(maxenkfpar),OIDHX(maxenkfpar)
    INTEGER :: OIDDOV(maxenkfpar),OIDHXOV(maxenkfpar)   !outvar index of oidd and oidhx
    INTEGER :: OIDDAA(maxenkfpar),OIDHXAA(maxenkfpar)   !area aggretate for oidd and oidhx
    REAL    :: OMIN(maxenkfpar),OMAX(maxenkfpar),OSIGMA(maxenkfpar)
    REAL    :: OSEMIM(maxenkfpar),ORESTM(maxenkfpar)
    REAL    :: OLSCALE(maxenkfpar),OGRIDSIZE(maxenkfpar)
    INTEGER :: OCORRTYPE(maxenkfpar)

    !forcing info
    INTEGER :: ITYPE(maxenkfpar)
    REAL    :: IMIN(maxenkfpar),IMAX(maxenkfpar),ISIGMA(maxenkfpar)
    INTEGER :: IID(maxenkfpar)
    REAL    :: ISEMIM(maxenkfpar),IRESTM(maxenkfpar)
    REAL    :: ILSCALE(maxenkfpar),IGRIDSIZE(maxenkfpar)
    INTEGER :: ICORRTYPE(maxenkfpar)

    !auxiliary info (nb. size from max_outvar)
    INTEGER :: AID(max_outvar)
    REAL    :: AMIN(max_outvar),AMAX(max_outvar)

    !additional help variables
    !INTEGER :: NewSwitch
    INTEGER :: io
    INTEGER :: linepos

    LOGICAL :: error
    LOGICAL :: nostrfound
    INTEGER :: iout, flow, aout

    !variables for localization
    REAL rx,ry !,ix,iy,iz,jx,jy,jz !,xy_dist,z_dist !,xy_scalefac,z_scalefac

    ! default state vector for initialization 
    !REAL, allocatable :: xini(:)

    !-----------------------------------------------------------
    ! 0) Important start: initilize assim%info with DEFAULT values
    !-----------------------------------------------------------
    n_Result = 0  !initial return status ok
    CALL initialize_assim_categories_HYPE(assimData%Info,stateinfo)
    CALL initialize_assim_info(assimData%Info)

    ! file unit number for ENKF files (make sure no overlap with HYPE or HYSS)
    fid_assim_info   = 70001
    !fid_assim_final  = 70002 !CP not used
    !fid_assim_min    = 70003
    !fid_assim_max    = 70004
    !fid_assim_median = 70005
    !fid_assim_mean   = 70006
    !fid_assim_min95  = 70007
    !fid_assim_max95  = 70008
    !fid_assim_std    = 70009

    fid_assim_bin_base = 70010

    !Move assim_open_timeseries until after the model structure identification
    !CALL assim_open_timeseriesoutput(resultdir) ! open the assimilation time series outputs

    ! ------------------------------------
    ! 1) read the input file AssimInfo.TXT
    ! ------------------------------------
    parid = 0 ; obsid = 0 ; inpid = 0 ; ncom = 0 ; ninfo = 0 ; outid = 0 ; nblank = 0
    iass = 0
    error = .FALSE.
    filename=TRIM(dir)//'AssimInfo.txt'
    OPEN(UNIT = fid_assim_info,FILE = filename, STATUS = 'old')   !use funit_temp instead
    WRITE(6,*)'Starting to READ AssimInfo.txt ...'
    ! read line by line, skip comment lines
    DO
      READ(fid_assim_info,'(a1000)',ERR=800,END=700,IOSTAT=io) newline
      IF(LEN_TRIM(newline).GT.0)THEN !check for empty lines
        newline = TRIM(newline)
        linepos = 1
        CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,code,nostrfound,errstr)    !Read code
        IF(code(1:2)==errstr)THEN
          error = .TRUE.
          code=newline(1:strvlen)
          EXIT
        ENDIF
        IF(code(1:2)==comment_str)THEN !Comment is !! now
          ncom = ncom+1
          CYCLE  
        ENDIF
        IF(nostrfound)THEN  !empty row
          nblank = nblank + 1
          CYCLE
        ENDIF
        ninfo = ninfo + 1   !number of rows with information in the file
        ! GENERAL SETTINGS (g_...)
        IF(code(1:7)=='g_xyloc')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) assimData%info%xy_scalefac         ! localization scaling length, XY directions
        ENDIF
        IF(code(1:6)=='g_zloc')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) assimData%info%z_scalefac         ! localization scaling length, Z direction
        ENDIF
        IF(code(1:4)=='g_ne')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) assimData%info%nE                  ! Ensemble Size!!
        ENDIF
        IF(code(1:5)=='g_cnc')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint                   !CHANGE TO USE Y or N?
          IF(newint==1) assimData%info%collapsenoncontrolled = .TRUE.         ! Collapse non-controlled variable ensembles to mean?
        ENDIF
        IF(code(1:4)=='g_fa')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint                   !CHANGE TO USE Y or N?
          IF(newint==1) assimData%info%FA = .TRUE.         ! include "A" (outvar) in filter, default is FALSE
        ENDIF
        IF(code(1:4)=='g_ff')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint
          IF(newint==1)assimData%info%FF = .TRUE.         ! include "F" (forcing) in filter (not needed in HYPE if A is used), default is FALSE
        ENDIF
        IF(code(1:4)=='g_mv')THEN                    ! missing value indicator, -9999, NOT USED?? USE missing_value or use this?
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) assimData%info%missing
        ENDIF
        IF(code(1:9)=='g_meanout')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint
          IF(newint==1)assimData%info%meanout = .TRUE.         ! ensemble MEAN in outputs, default
          IF(newint<=0)assimData%info%meanout = .FALSE.        ! ensemble MEDIAN in outputs
        ENDIF
        IF(code(1:9)=='g_statout')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint
          assimData%info%nstatout = newint         ! print out statistical results 2: min and max, 5:min/max, median,+-0.025percentiles, default is 0 (no output)
                                                   !The output files are 002:minimum, 003:maximum, 004:0.025-percentile, 005:median, 006:0.975-percentile 
        ENDIF
        IF(code(1:9)=='g_usebinx')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint
          assimData%info%useBinFilesX = newint         ! use direct access binary files to store state ensemble data instead of allocated memory, default is 0 (no bin-file), 1 is one bin-fil 2 is several bin-files
          !IF(newint==1) assimData%info%useBinFiles = .TRUE.         ! use direct access binary files to store ensemble data instead of allocated memory, default is FALSE
        ENDIF
        IF(code(1:10)=='g_usebinfa')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint
          assimData%info%useBinFilesFA = newint         ! use direct access binary files to store other (F,A) ensemble data instead of allocated memory, default is 0 (no bin-file), 1 is one bin-fil 2 is several bin-files
        ENDIF
        IF(code(1:9)=='g_inibinx')THEN
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) newint                   !CHANGE TO USE Y or N?
          IF(newint==1) assimData%info%initializeFromBinFiles = .TRUE.         ! Initialize State Ensembles from Existing Bin-files?
        ENDIF
        ! Assimilation switches (a_...) (default is true for all these)
        IF(code(1:2)=='a_')THEN
          iass = iass + 1   !count rows
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          !turn assimilate on or off for this category/variable
          IF(code(3:9)=='include')THEN
            assimData%info%assim_flag(iass)%add = .true.
          ELSEIF(code(3:9)=='exclude')THEN
            assimData%info%assim_flag(iass)%add = .false.
          ELSE
            error = .TRUE.
            EXIT
          ENDIF
          !set the category or variable from the name
          IF(code(13:20)=='category')THEN
            assimData%info%assim_flag(iass)%category = strvalue
          ELSEIF(code(13:16)=='name')THEN
            assimData%info%assim_flag(iass)%varname = strvalue
            !find statetype for a variable and save in assim_flag%category
            DO i = iass-1,1,-1
              IF(assimData%info%assim_flag(i)%category/='' .AND. &
                 assimData%info%assim_flag(i)%varname=='')THEN
                CALL get_statetype_from_category(assimData%Info,stateinfo,  &
                        assimData%info%assim_flag(i)%category,  &
                        assimData%info%assim_flag(iass)%varname,stname)
                assimData%info%assim_flag(iass)%category = stname
                EXIT
              ENDIF
            ENDDO
            IF(i==0)THEN
              WRITE(6,*) 'ERROR: Did not find category for variable',strvalue,'.'
              error = .TRUE.
              EXIT
            ENDIF
          ELSE
            error = .TRUE.
            EXIT
          ENDIF
        ENDIF
        !Observation settings (o_...)
        IF(code(1:1)=='o')THEN  
          obsid = obsid + 1
          !------------------------------------------------------------------------------------------------------------------------------------
          ! Format: o_comment idobs(4char), idmod(4char), enstype(0,1,2), min(real), max(real), sigma(real), SemiMeta(real), RestMeta(real), lscale(real), gridsize(real), corrType
          !------------------------------------------------------------------------------------------------------------------------------------
          !OIDD (observation identifier, "HYPE" output variable name)
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !idobs
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          !n_result = find_variable_index_type(strvalue,iout,flow)  !CP170202 added regional outvar
          n_result = find_variable_index_type(strvalue,iout,flow,aout)
          IF(n_Result/=0)THEN
            WRITE(6,*) 'ERROR: No variable found with that name'
            error = .TRUE.
            EXIT
          ENDIF
          OIDD(obsid) = iout
          OIDDAA(obsid) = aout

          !OIDHX (modelled observation identifier, "HYPE" output variable name)
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            WRITE(6,*) 'ERROR: No variable found with that name'
            error = .TRUE.
            EXIT
          ENDIF
          !n_result = find_variable_index_type(strvalue,iout,flow)  !CP170202 added regional outvar
          n_result = find_variable_index_type(strvalue,iout,flow,aout)
          IF(n_Result/=0)THEN
            error = .TRUE.
            EXIT
          ENDIF
          OIDHX(obsid) = iout
          OIDHXAA(obsid) = aout

          !Remaining numerical observation settings
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OTYPE(obsid)   !Code, ensemble type 0=off, 1=unconstrained with constant sd, 2=semi-constrained(min), 3=semi-constrained(max), 4=constrained(min,max)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OMIN(obsid)   !Real, minimum allowed value (EnsType=2,3,4)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OMAX(obsid)   !real, maximum allowed value (EnsType=2,3,4)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OSIGMA(obsid)   !real, constant standard deviation (constrained, EnsType=1) 
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OSEMIM(obsid)   !real, relative standard deviation (semi-constrained, EnsType=2&3)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) ORESTM(obsid)   !real, relative standard deviation (constrained, EnsType=4)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OLSCALE(obsid)   !real, correlation length for spatially correlated perturbations
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OGRIDSIZE(obsid)   !real, cell size(x&y) in the 2D grid used for generating spatially correlated perturbations
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) OCORRTYPE(obsid)   !Code, correlation function: 1 Gaussian, 2 Compact 5th degree polynomial, 3 Power law

        ENDIF
        !Input: meteorological forcing data (on-going: change to use 'F' as prefix to be consistent with the rest of the code)
        IF(code(1:1)=='f')THEN
          inpid = inpid + 1
          !------------------------------------------------------------------------------------------------------------------------------------
          ! Format:     id, EnsType, min, max, sigma, SemiMeta, RestMeta, lscale, gridsize, corrType
          !------------------------------------------------------------------------------------------------------------------------------------
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) IID(inpid)   !Code 1=Pobs,2=Tobs => todo: switch to use "obs" names, T, TMIN, TMAX, P, SW, etc
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) ITYPE(inpid)   !Code, ensemble type 0=off, 1=unconstrained with constant sd, 2=semi-constrained(min), 3=semi-constrained(max), 4=constrained(min,max)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) IMIN(inpid)   !Real, minimum allowed value (EnsType=2,3,4)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) IMAX(inpid)   !real, maximum allowed value (EnsType=2,3,4)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) ISIGMA(inpid)   !real, constant standard deviation (constrained, EnsType=1) 
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) ISEMIM(inpid)   !real, relative standard deviation (semi-constrained, EnsType=2&3)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) IRESTM(inpid)   !real, relative standard deviation (constrained, EnsType=4)
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) ILSCALE(inpid)   !real, correlation length for spatially correlated perturbations
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) IGRIDSIZE(inpid)   !real, cell size(x&y) in the 2D grid used for generating spatially correlated perturbations
          !---
          CALL read_next_codestr_on_line(linelen,strvlen,linepos,newline,strvalue,nostrfound,errstr) !id
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) ICORRTYPE(inpid)   !Code, correlation function: 1 Gaussian, 2 Compact 5th degree polynomial, 3 Power law
        ENDIF
      ENDIF
    ENDDO
    IF(error)THEN
      WRITE(6,*) 'ERROR: reading AssimInfo.txt (',TRIM(filename),'). Code=',code
      n_Result = 1
      CLOSE(fid_assim_info)
      RETURN
    ENDIF
700 CLOSE(fid_assim_info)
    WRITE(6,*)' Finished reading AssimInfo.txt'
    WRITE(6,*)' ... n:o empty lines: ', nblank
    WRITE(6,*)' ... n:o comment lines: ', ncom
    WRITE(6,*)' ... n:o data lines: ', ninfo

    ! ne = n:o ensemble members, just read from AssimInfo.txt
    ne = assimData%info%nE
    ! nFlag = n:o rows with assimilate information, read from AssimInfo.txt
    assimData%info%nFlag = iass
    ! Check number of statistical outputs.
    IF(assimData%info%nstatout>5+ne) assimData%info%nstatout=5+ne

    ! ---------------------------------------------------------------------------------
    ! 2) Check (HYPE) model structure and allocate and initialize the Assimilation information data
    ! ---------------------------------------------------------------------------------

    !STATES, Check model structure and count nX: the number of state variables (per subbasin (or aquifer?))
    CALL checkModelStructure_States_HYPE2(assimData%info%nX,stateinfo)
    
    !FORCING DATA, check number of input variables, and if they are included in Kalman Filter Analysis
    assimData%info%nF      = inpid   ! n:o inputs listed in assimInfo.txt
    DO i=1,inpid
      !check ensemble type and reduce nF if ensType=0
      IF(ITYPE(i).LE.0) assimData%info%nF  = assimData%info%nF - 1
    ENDDO
    IF(assimData%info%nF.LE.0)THEN
      !make sure input data Kalman filtering is FALSE if no input data used
      assimData%info%FF  = .false.
      assimData%info%nF  = 0
    ENDIF

    !AUXILIARIES & OBSERVATIONS, number of auxiliaries (outvars) requested through info.txt and AssimInfo.txt
    !Update outvar to be calculated with AssimInfo-outvars (OTYPE) (outvar from info.txt already set)
    DO i=1,obsid
      IF(OTYPE(i)>=1)THEN
        CALL set_outvar_for_variable(OIDD(i),OIDDAA(i),OIDDOV(i),noutvar)
        CALL set_outvar_for_variable(OIDHX(i),OIDHXAA(i),OIDHXOV(i),noutvar)
      ENDIF
    ENDDO
    !Number of outvar is noutvar, save all and set max/min values
    DO i=1,noutvar
      AID(i)=i
      AMIN(i)=assim_defmin
      AMAX(i)=assim_defmax
    ENDDO
    assimData%info%na = noutvar
    
    !OBSERVATIONS, check number of observation variable types from assimInfo.txt
    assimData%info%nObs = obsid  ! n:o obs. types
    DO i=1,obsid
      IF(OTYPE(i).LE.0) then
        assimData%info%nObs   = assimData%info%nObs - 1 !dont include if OTYPE was set to 0
      ENDIF
    ENDDO
    IF(assimData%info%nObs.LE.0)  assimData%info%nObs = 0
    
    !NCOORD, NLOC, number of spatial representations (sets of coordinates), and number of localization matrices
    assimData%info%NCOORD = 2 !subbasins, for now, we can only handle one set of coordinates (subbasin dimension) CP161213 changed to include aquifers
    assimData%info%NLOC   = 1 !subbasins, for each set of coordinates, we also need one set of localization matrices
    
    !ALLOCATE ensemble vectors with the general routine (X,A,Obs,F), sets fid_assim_bin
    Call allocate_assim_ensemble_vectors(assimData,fid_assim_bin_base)
    
    !---------------------------------------------------------------------------
    !3) Allocate, initialize and assign options and values to all ensemble data
    !---------------------------------------------------------------------------
    !COORDINATES
    ALLOCATE(assimData%Coordinates(assimData%info%NCOORD))
    !First set of coordinates representing the subbasin centre points, x,y, elevation
    ALLOCATE(assimData%Coordinates(1)%x(nsub))
    ALLOCATE(assimData%Coordinates(1)%y(nsub))
    ALLOCATE(assimData%Coordinates(1)%z(nsub))   
    assimData%Coordinates(1)%n = nsub
    DO i=1,nsub
      assimData%Coordinates(1)%x(i) = basin(i)%xcoord
      assimData%Coordinates(1)%y(i) = basin(i)%ycoord
      assimData%Coordinates(1)%z(i) = basin(i)%elev
    ENDDO
    !Second set of coordinates representing the aquifers weighted from subbasin centre points, x,y !CP161214 added aquifers
    IF(naquifers>0)THEN
      ALLOCATE(assimData%Coordinates(2)%x(naquifers))
      ALLOCATE(assimData%Coordinates(2)%y(naquifers))
      ALLOCATE(assimData%Coordinates(2)%z(naquifers))   
      assimData%Coordinates(2)%n = naquifers
      DO i=1,naquifers
        CALL calculate_coordinates_for_aquifer(i,rx,ry)
        assimData%Coordinates(2)%x(i) = rx
        assimData%Coordinates(2)%y(i) = ry
        assimData%Coordinates(2)%z(i) = 0.  !CP: what to do with this one?
      ENDDO
    ENDIF
    !LOCALIZATION
    IF(naquifers>0) assimData%info%NLOC   = 2 !subbasins and aquifers, for each set of coordinates, we also need one set of localization matrices
    IF(ALLOCATED(assimData%locCXY))DEALLOCATE(assimData%locCXY)
    ALLOCATE(assimData%locCXY(assimData%info%NLOC))

    !STATES (nb. the "interface" routine that knows how HYPE is structured)
    varID=1   !varID is a counter that is used in the following subroutines to identify the state ensemble matrices in the pre-defined order (we will improve this system later, using pointers)
    CALL allocate_and_initialize_model_state_ensembles_HYPE2(assimData%info%ne,assimData%X,varID,assimData%info,stateinfo)
    IF(ALLOCATED(assimData%Info%assim_flag)) DEALLOCATE(assimData%Info%assim_flag)  !information transferred to assimData%X%X%assimilate
    
    !FORCING (nb. the varID counter should not be re-initialized here (used for fileID-counter))
    !CALL allocate_and_initialize_model_forcing_ensembles_HYPE(assimData%info%ne,nsub,assimData%F,varID,inpid,ITYPE,IID,IMIN,IMAX,ISIGMA,ISEMIM,IRESTM,ILSCALE,IGRIDSIZE,ICORRTYPE,assimData%Coordinates(1)%x,assimData%Coordinates(1)%y) !CP161202 added assimData%Info for use bin-file
    CALL allocate_and_initialize_model_forcing_ensembles_HYPE(assimData%info%ne,nsub,assimData%F,varID,assimData%info,inpid,ITYPE,IID,IMIN,IMAX,ISIGMA,ISEMIM,IRESTM,ILSCALE,IGRIDSIZE,ICORRTYPE,assimData%Coordinates(1)%x,assimData%Coordinates(1)%y)

    !AUXILIARIES (outvar)
    CALL allocate_and_initialize_model_auxiliary_ensembles_HYPE(assimData%info%ne,nsub,assimData%A,assimData%info,varID,assimData%info%na,AID,AMIN,AMAX)
    
    !PARAMETERS
    !later...
    
    !OBSERVATIONS
    !CALL allocate_and_initialize_observation_ensembles_HYPE(assimData%info%ne,nsub,assimData%Obs,varID,obsid,OIDD,OIDHX,OTYPE,OMIN,OMAX,OSIGMA,OSEMIM,ORESTM,OLSCALE,OGRIDSIZE,OCORRTYPE,assimData%Coordinates(1)%x,assimData%Coordinates(1)%y)  !CP170202 added regional outvar
    CALL allocate_and_initialize_observation_ensembles_HYPE(assimData%info%ne,nsub,assimData%Obs,varID,obsid,OIDD,OIDDOV,OIDHXOV,OTYPE,OMIN,OMAX,OSIGMA,OSEMIM,ORESTM,OLSCALE,OGRIDSIZE,OCORRTYPE,assimData%Coordinates(1)%x,assimData%Coordinates(1)%y)

    ! Finally, assign initial values using the model_to_ensemble function
    DO i=1,assimData%info%nE
      CALL model_to_ensemble2(i,assimData%X,assimData%A,assimData%info%nA,stateinfo,assimData%info%initializeFromBinFiles)
    ENDDO
    RETURN

800 CONTINUE
    n_Result = 1
    WRITE(6,*) 'ERROR: reading AssimInfo.txt. io=',io
    CLOSE(fid_assim_info)
    RETURN
!801 CONTINUE   !CP170505 not used
!    n_Result = 1
!    WRITE(6,*) 'ERROR: reading AssimInfo.txt. information line',ninfo
!    CLOSE(fid_assim_info)
!    RETURN

  END SUBROUTINE assim_initialize 

!-----------------------------------------------------------------------------------------------------------
! SECOND BLOCK: routines for read/write between model and ensemble
!
! TODO (20161101):  1) update to current HYPE structure (see initialization code above)
!                   2) code read/write to binary files (nb. most efficient matrix orientation, probably rows=ensemble member, cols=subbasins)
!-----------------------------------------------------------------------------------------------------------
!> \brief Copy all current model states to the ensemble holding variables or files

!>Purpose: Write model variables (states, parameters, etc) to the ensemble matrices.
!>The user provide the necessary code to write data from the ensemble member i_ens to the ensemble matrix.
  SUBROUTINE model_to_ensemble2(i_ens,assimX,assimA,nA,stateinfo,skipX)
    
    !Argument declarations  
    INTEGER, INTENT(IN) :: i_ens                       !<current ensemble
    TYPE(assim_state_ensemble_type)  :: assimX(:)      !<ensemble data for model states
    TYPE(assim_state_ensemble_type)  :: assimA(:)      !<ensemble data for auxiliaries
    INTEGER, INTENT(IN) :: nA                          !<number of auxiliaries (outvars) requested
    TYPE(stateinfotype), INTENT(IN)   :: stateinfo(:)  !<Information about state variables
    LOGICAL :: skipX                                  !<logical switch to skip updating X (used when initializing from existing binfiles)

    !Local variables
    INTEGER i, k, nstate
    
    !1) Model states to state ensemble X
    IF(.NOT.skipX)THEN
      k=1   !=varID will be increased in subroutine calls
      nstate = SIZE(stateinfo)
      DO i = 1,nstate
        IF(stateinfo(i)%allok)THEN
          IF(stateinfo(i)%ndim==1)THEN
            CALL write1D_to_X(stateinfo(i)%svpoint%d1,stateinfo(i)%dims(1),k,i_ens,assimX)
          ELSEIF(stateinfo(i)%ndim==2)THEN
            IF(ASSOCIATED(stateinfo(i)%svpoint%d2))THEN          !both real and integer state of 3dim exist
              CALL write2D_to_X(stateinfo(i)%svpoint%d2,stateinfo(i)%dims(1),stateinfo(i)%dims(2),k,i_ens,assimX) 
            ELSE
              CALL write2D_to_X(REAL(stateinfo(i)%svpoint%d2i),stateinfo(i)%dims(1),stateinfo(i)%dims(2),k,i_ens,assimX) 
            ENDIF
          ELSEIF(stateinfo(i)%ndim==3)THEN
            IF(ASSOCIATED(stateinfo(i)%svpoint%d3))THEN          !both real and integer state of 3dim exist
              CALL write3D_to_X(stateinfo(i)%svpoint%d3,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),k,i_ens,assimX) 
            ELSE
             CALL write3D_to_X(REAL(stateinfo(i)%svpoint%d3i),stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),k,i_ens,assimX) 
           ENDIF
          ELSEIF(stateinfo(i)%ndim==4)THEN
            CALL write4D_to_X(stateinfo(i)%svpoint%d4,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),stateinfo(i)%dims(4),k,i_ens,assimX) 
          ENDIF
        ENDIF
      ENDDO

      !Variables for river flows (module variables in HYPEVARIABLES)
      IF(ALLOCATED(Qmax))   CALL write2D_to_X(Qmax,2,nsub,k,i_ens,assimX)
      IF(ALLOCATED(Q2max))  CALL write2D_to_X(Q2max,2,nsub,k,i_ens,assimX)
      IF(ALLOCATED(iQmax))  CALL write2D_to_X(REAL(iQmax),2,nsub,k,i_ens,assimX)
      IF(ALLOCATED(iQ2max)) CALL write2D_to_X(REAL(iQ2max),2,nsub,k,i_ens,assimX)
      IF(ALLOCATED(accdiff))CALL write1D_to_X(accdiff,nsub,k,i_ens,assimX)
    ENDIF
    !2) additional model var. to A ensemble (fluxes and auxiliaries)
    !This is done for requested outvars (nA).
    IF(nA.GT.0)THEN
      k=1
      DO i=1,nA
        CALL write1D_to_X(outvar(:,assimA(i)%info%varID),nsub,k,i_ens,assimA)
      ENDDO
    ENDIF 
    
  END SUBROUTINE model_to_ensemble2

  !>\brief Save data for one ensemble to the state ensembles in allocated matrix or file.
  SUBROUTINE write1D_to_X(x1d,ni,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, intent(IN) :: x1d(:)
    INTEGER, intent(IN) :: ni,i_ens
    INTEGER, intent(INOUT) :: k       !=varID

    INTEGER i
    IF(ALLOCATED(assimVar(k)%x%x))THEN
      DO i=1,ni
        assimVar(k)%x%x(i,i_ens) = x1d(i)
      ENDDO
    ELSE
      WRITE(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x1d
    ENDIF
    k=k+1
  END SUBROUTINE write1D_to_X

  !>\brief Save data for one ensemble to the state ensembles in allocated matrix or file.
  SUBROUTINE write2D_to_X(x2d,ni,nj,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, intent(IN) :: x2d(:,:)
    INTEGER, intent(IN) :: ni,nj, i_ens
    INTEGER, intent(INOUT) :: k

    INTEGER i,j
    DO i=1,ni
      IF(ALLOCATED(assimVar(k)%x%x))THEN
        DO j=1,nj
          assimVar(k)%x%x(j,i_ens) = x2d(i,j)
        ENDDO
      ELSE
        WRITE(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x2d(i,:)
      ENDIF
      k=k+1
    ENDDO
 
  END SUBROUTINE write2D_to_X

  !>\brief Save data for one ensemble to the state ensembles in allocated matrix or file.
  SUBROUTINE write3D_to_X(x3d,ni,nj,nl,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, intent(IN) :: x3d(:,:,:)
    INTEGER, intent(IN) :: ni,nj,nl, i_ens
    INTEGER, intent(INOUT) :: k

    INTEGER i,j,l
    DO i=1,ni
      DO j=1,nj
        IF(ALLOCATED(assimVar(k)%x%x))THEN
          DO l=1,nl
            assimVar(k)%x%x(l,i_ens) = x3d(i,j,l)
          ENDDO
        ELSE
          WRITE(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x3d(i,j,:)
        ENDIF
        k=k+1
      ENDDO
    ENDDO
 
  END SUBROUTINE write3D_to_X

  !>\brief Save data for one ensemble to the state ensembles in allocated matrix or file.
  SUBROUTINE write4D_to_X(x4d,ni,nj,nl,nm,k,i_ens,assimVar)
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, intent(IN) :: x4d(:,:,:,:)
    INTEGER, intent(IN) :: ni,nj,nl,nm,i_ens
    INTEGER, intent(INOUT) :: k

    INTEGER i,j,l,m
    DO i=1,ni
      DO j=1,nj
        DO l=1,nl
          IF(ALLOCATED(assimVar(k)%x%x))THEN
            DO m=1,nm    
              assimVar(k)%x%x(m,i_ens) = x4d(i,j,l,m)
            ENDDO
          ELSE
            WRITE(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x4d(i,j,l,:)
          ENDIF
          k=k+1
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE write4D_to_X

!-----------------------------------------------------------------------------------------------------------
!> \brief Copy all model states from the ensemble variables or files.
!
!> Purpose: Write ensemble data to the model variables (states, outvars, forcing).
!>
!> Reverse of subroutine model_to_ensemble
! -------------------------------------------------
  SUBROUTINE ensemble_to_model2(i_ens,nX,nA,nF,assimX,assimA,assimF,stateinfo)

    !Argument declarations  
    INTEGER, INTENT(IN) :: i_ens    !<current ensemble
    INTEGER, INTENT(IN) :: nX       !<number of X (not used)
    INTEGER, INTENT(IN) :: nA       !<number of auxiliaries (outvars) requested
    INTEGER, INTENT(IN) :: nF       !<number of forcing variables
    TYPE(assim_state_ensemble_type)  :: assimX(:)   !<ensemble data for model states
    TYPE(assim_state_ensemble_type)  :: assimA(:)   !<ensemble data for auxiliaries
    TYPE(assim_input_ensemble_type)  :: assimF(:)   !<ensemble data for forcing
    !TYPE(assim_state_ensemble_type)  :: assimP(:) !to come?
    TYPE(stateinfotype), INTENT(INOUT) :: stateinfo(:)  !<Information about state variables

    !Local variables
    INTEGER i, k, nstate
    REAL,ALLOCATABLE :: help2(:,:),help3(:,:,:)  !needed to handle integer pointers

    !1) Model state variables from the X ensemble
    k=1   !=varID will be increased in subroutine calls
    nstate = SIZE(stateinfo)
    DO i = 1,nstate
      IF(stateinfo(i)%allok)THEN
        IF(stateinfo(i)%ndim==1)THEN
          CALL writeX_to_1D(stateinfo(i)%svpoint%d1,stateinfo(i)%dims(1),k,i_ens,assimX)
        ELSEIF(stateinfo(i)%ndim==2)THEN
          IF(ASSOCIATED(stateinfo(i)%svpoint%d2))THEN          !both real and integer state of 3dim exist
            CALL writeX_to_2D(stateinfo(i)%svpoint%d2,stateinfo(i)%dims(1),stateinfo(i)%dims(2),k,i_ens,assimX) 
          ELSE
            ALLOCATE(help2(stateinfo(i)%dims(1),stateinfo(i)%dims(2)))
            CALL writeX_to_2D(help2,stateinfo(i)%dims(1),stateinfo(i)%dims(2),k,i_ens,assimX) 
            stateinfo(i)%svpoint%d2i = NINT(help2)
            DEALLOCATE(help2)
          ENDIF
        ELSEIF(stateinfo(i)%ndim==3)THEN
          IF(ASSOCIATED(stateinfo(i)%svpoint%d3))THEN          !both real and integer state of 3dim exist
            CALL writeX_to_3D(stateinfo(i)%svpoint%d3,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),k,i_ens,assimX) 
          ELSE
            ALLOCATE(help3(stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3)))
            CALL writeX_to_3D(help3,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),k,i_ens,assimX) 
            stateinfo(i)%svpoint%d3i = NINT(help3)
            DEALLOCATE(help3)
          ENDIF
        ELSEIF(stateinfo(i)%ndim==4)THEN
          CALL writeX_to_4D(stateinfo(i)%svpoint%d4,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),stateinfo(i)%dims(4),k,i_ens,assimX) 
        ENDIF
      ENDIF
    ENDDO

    !Variables for river flows (module variables in HYPEVARIABLES)
    IF(ALLOCATED(Qmax))   CALL writeX_to_2D(Qmax,2,nsub,k,i_ens,assimX)
    IF(ALLOCATED(Q2max))  CALL writeX_to_2D(Q2max,2,nsub,k,i_ens,assimX)
    IF(ALLOCATED(iQmax))  CALL writeX_to_2D_real2int(iQmax,2,nsub,k,i_ens,assimX)
    IF(ALLOCATED(iQ2max)) call writeX_to_2D_real2int(iQ2max,2,nsub,k,i_ens,assimX)
    IF(ALLOCATED(accdiff))call writeX_to_1D(accdiff,nsub,k,i_ens,assimX)

    !2) auxiliary variables from A ensemble
    IF(nA.GT.0)THEN
      k=1
      DO i=1,nA
        CALL writeX_to_1D(outvar(:,assimA(i)%info%varID),nsub,k,i_ens,assimA)
      ENDDO
    ENDIF 
        
    !3) Forcing from F ensemble
    IF(nF.GT.0)THEN
      DO i = 1,nF
        IF(assimF(i)%info%varID.EQ.1)THEN
          !varID = 1 >> precipitation
          k=assimF(i)%info%varID    !local variable needed because k is changed in writeX_to_1D
          CALL writeF_to_1D(preci,nsub,k,i_ens,assimF(i))
        ELSEIF(assimF(i)%info%varID.EQ.2)THEN
          !varID = 2 >> temperature
          k=assimF(i)%info%varID    !local variable needed because k is changed in writeX_to_1D
          CALL writeF_to_1D(tempi,nsub,k,i_ens,assimF(i))
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE ensemble_to_model2

  !>\brief Get one ensemble member from the state ensembles.
  !>
  !>The state ensembles may be in an allocated matrix or read from bin-file.
  SUBROUTINE writeF_to_1D(x1d,ni,k,i_ens,assimVar)
    
    !Argument declarations
    REAL, INTENT(OUT) :: x1d(:)
    INTEGER, INTENT(IN) :: ni !<dimension of matrix x1d
    INTEGER, INTENT(INOUT) :: k
    INTEGER, INTENT(IN) :: i_ens
    TYPE(assim_input_ensemble_type),INTENT(IN) :: assimVar

    INTEGER i
    IF(ALLOCATED(assimVar%x%x))THEN
      DO i=1,ni
        x1d(i) = assimVar%x%x(i,i_ens)
      ENDDO
    ELSE
      READ(assimVar%X%fileID,REC=assimVar%X%rec+i_ens) x1d
    ENDIF
    k=k+1 !=varID
  END SUBROUTINE writeF_to_1D

  !>\brief Get one ensemble member from the state ensembles.
  !>
  !>The state ensembles may be in an allocated matrix or read from bin-file.
  SUBROUTINE writeX_to_1D(x1d,ni,k,i_ens,assimVar)
    
    !Argument declarations
    REAL, INTENT(OUT) :: x1d(:)
    INTEGER, INTENT(IN) :: ni !<dimension of matrix x1d
    INTEGER, INTENT(INOUT) :: k
    INTEGER, INTENT(IN) :: i_ens
    TYPE(assim_state_ensemble_type),INTENT(IN) :: assimVar(:)

    INTEGER i

    IF(ALLOCATED(assimVar(k)%x%x))THEN
      DO i=1,ni
        x1d(i) = assimVar(k)%x%x(i,i_ens)
      ENDDO
    ELSE
      READ(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x1d
    ENDIF
    k=k+1 !=varID
  END SUBROUTINE writeX_to_1D

  !>\brief Get one ensemble member from the state ensembles
  !>
  !>The state ensembles may be in an allocated matrix or read from bin-file.
  SUBROUTINE writeX_to_2D(x2d,ni,nj,k,i_ens,assimVar)
    
    !Argument declarations
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, INTENT(OUT) :: x2d(:,:)
    INTEGER, INTENT(IN) :: ni,nj, i_ens
    INTEGER, INTENT(INOUT) :: k

    INTEGER i,j
    DO i=1,ni
      IF(ALLOCATED(assimVar(k)%x%x))THEN
        DO j=1,nj
          x2d(i,j) = assimVar(k)%x%x(j,i_ens)
        ENDDO
      ELSE
        READ(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x2d(i,:)
      ENDIF
      k=k+1 !=varID
    ENDDO

  END SUBROUTINE writeX_to_2D
	
  !>\brief Get one ensemble member from the state ensembles
  !>
  !>The state ensembles may be in an allocated matrix or read from bin-file.
  SUBROUTINE writeX_to_3D(x3d,ni,nj,nl,k,i_ens,assimVar)
   
    !Argument declarations
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, intent(OUT) :: x3d(:,:,:)
    INTEGER, intent(IN) :: ni,nj,nl, i_ens
    INTEGER, intent(INOUT) :: k

    INTEGER i,j,l
    DO i=1,ni
      DO j=1,nj
        IF(ALLOCATED(assimVar(k)%x%x))THEN
          DO l=1,nl
            x3d(i,j,l) = assimVar(k)%x%x(l,i_ens)
          ENDDO
        ELSE
          READ(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x3d(i,j,:)
        ENDIF
        k=k+1
      ENDDO
    ENDDO

  END SUBROUTINE writeX_to_3D

  !>\brief Get one ensemble member from the state ensembles
  !>
  !>The state ensembles may be in an allocated matrix or read from bin-file.
  SUBROUTINE writeX_to_4D(x4d,ni,nj,nl,nm,k,i_ens,assimVar)
    
    !Argument declaration
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    REAL, intent(OUT) :: x4d(:,:,:,:)
    INTEGER, intent(IN) :: ni,nj,nl,nm,i_ens
    INTEGER, intent(INOUT) :: k

    INTEGER i,j,l,m
    DO i=1,ni
      DO j=1,nj
        DO l=1,nl
          IF(ALLOCATED(assimVar(k)%x%x))THEN
            DO m=1,nm    
              x4d(i,j,l,m) = assimVar(k)%x%x(m,i_ens)
            ENDDO
          ELSE
            READ(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x4d(i,j,l,:)
          ENDIF
          k=k+1
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE writeX_to_4D

  !>\brief Get one ensemble member from the state ensembles
  !>
  !>The recieving matrix is an integer, so transformation may be used.
  !>The state ensembles may be in an allocated matrix or read from bin-file.
  SUBROUTINE writeX_to_2D_real2int(x2d,ni,nj,k,i_ens,assimVar)
    
    !Argument declarations
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER, intent(OUT) :: x2d(:,:)
    INTEGER, intent(IN) :: ni,nj, i_ens
    INTEGER, intent(INOUT) :: k

    INTEGER i,j
    DO i=1,ni
      IF(ALLOCATED(assimVar(k)%x%x))THEN
        DO j=1,nj
          x2d(i,j) = INT(assimVar(k)%x%x(j,i_ens))
        ENDDO
      ELSE
        READ(assimVar(k)%X%fileID,REC=assimVar(k)%X%rec+i_ens) x2d(i,:)
      ENDIF
      k=k+1
    ENDDO

  END SUBROUTINE writeX_to_2D_real2int

  !----------------------------------------------------------
	!>\brief Collect wanted statistics of ensemble to output
	!>
  !> - Writes mean or median of the ensemble to outvar for printout
  !> - Writes the min/max etc of the ensemble to outvar for printout to 
  !>the "sequence" model output files
  !> - Writes individual ensembles to outvar for printout to 
  !>the "sequence" model output files
  SUBROUTINE statistics_to_modeloutput(nA,nE,assimA,istat,meanORmedian)
    
    !Argument declarations
    INTEGER, INTENT(IN) :: nA      !<number of auxiliaries (outvars) requested
    INTEGER, INTENT(IN) :: nE      !<number of ensembles
    TYPE(assim_state_ensemble_type),INTENT(IN) :: assimA(:) !<assimilation data on auxiliary variables (outvar)
    INTEGER, INTENT(IN) :: istat   !<choice of statistics, 2=min,3=max,4=0.025q,5=median,6=0.975q,7=1st ensemble, 8=2nd ensemble...
    LOGICAL,INTENT(IN),OPTIONAL :: meanORmedian  !< flag for setting mean (true) or median (false)

    !Local variables
    INTEGER i, j, k
    REAL, ALLOCATABLE:: x1d(:)
    
    ! auxiliary variables from A ensemble
    IF(nA.GT.0)THEN
      IF(PRESENT(meanORmedian))THEN
        !For standard output mean or median is chosen
        IF(meanORmedian)THEN  !write mean
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outmean(j)
            ENDDO
          ENDDO
        ELSE  !write median
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outquant(2,j)
            ENDDO
          ENDDO
        ENDIF
      ELSE
        !For statistical extra output
        IF(istat==2)THEN
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outmin(j)
            ENDDO
          ENDDO
        ELSEIF(istat==3)THEN
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outmax(j)
            ENDDO
          ENDDO
        ELSEIF(istat==4)THEN
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outquant(1,j)
            ENDDO
          ENDDO
        ELSEIF(istat==5)THEN
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outquant(2,j)
            ENDDO
          ENDDO
        ELSEIF(istat==6)THEN
          DO i=1,nA
            DO j=1,nsub
             outvar(j,assimA(i)%info%varID) = assimA(i)%x%outquant(3,j)
            ENDDO
          ENDDO
        ELSEIF(istat>6)THEN
          !For individual ensemble output, istat=6-5+nE
          k=1
          DO i=1,nA
            ALLOCATE(x1d(nsub))
            CALL writeX_to_1D(x1d,nsub,k,istat-6,assimA)
            DO j=1,nsub
              outvar(j,assimA(i)%info%varID) = x1d(j)  !set outvar
            ENDDO
            DEALLOCATE(x1d)
          ENDDO
        ENDIF
      ENDIF
    ENDIF 

  END SUBROUTINE statistics_to_modeloutput

  !-----------------------------------------------------------------------------------------------------------
  !>\brief Collect mean or median of ensemble to model states and auxiliary variables
  !>
  !> Purpose: write the mean or median of the ensemble to the model domain, for the standard model output files.
  ! method: same as ensemble_to_model
  SUBROUTINE meanORmedian_to_model2(nA,assimX,assimA,meanORmedian,stateinfo)
    
    !Argument declarations
    INTEGER, INTENT(IN) :: nA       !<number of auxiliaries (outvars) requested
    TYPE(assim_state_ensemble_type)  :: assimX(:)   !<ensemble data for model states
    TYPE(assim_state_ensemble_type)  :: assimA(:)   !<ensemble data for auxiliaries
    LOGICAL, INTENT(IN) :: meanORmedian   !< flag for getting mean (true) or median (false)
    TYPE(stateinfotype), INTENT(INOUT) :: stateinfo(:)  !<Information about state variables

    !Local variables
    INTEGER i, j, k, nstate
    REAL,ALLOCATABLE :: help2(:,:),help3(:,:,:)  !needed to handle integer pointers

    !1) Model state variables from the X ensemble
    k=1   !=varID will be increased in subroutine calls
    nstate = SIZE(stateinfo)
    DO i = 1,nstate
      IF(stateinfo(i)%allok)THEN
        IF(stateinfo(i)%ndim==1)THEN
          CALL writeXmedian_to_1D(stateinfo(i)%svpoint%d1,stateinfo(i)%dims(1),k,assimX,meanORmedian)
        ELSEIF(stateinfo(i)%ndim==2)THEN
          IF(ASSOCIATED(stateinfo(i)%svpoint%d2))THEN          !both real and integer state of 3dim exist
            CALL writeXmedian_to_2D(stateinfo(i)%svpoint%d2,stateinfo(i)%dims(1),stateinfo(i)%dims(2),k,assimX,meanORmedian) 
          ELSE
            ALLOCATE(help2(stateinfo(i)%dims(1),stateinfo(i)%dims(2)))
            CALL writeXmedian_to_2D(help2,stateinfo(i)%dims(1),stateinfo(i)%dims(2),k,assimX,meanORmedian) 
            stateinfo(i)%svpoint%d2i = NINT(help2)
            DEALLOCATE(help2)
          ENDIF
        ELSEIF(stateinfo(i)%ndim==3)THEN
          IF(ASSOCIATED(stateinfo(i)%svpoint%d3))THEN          !both real and integer state of 3dim exist
            CALL writeXmedian_to_3D(stateinfo(i)%svpoint%d3,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),k,assimX,meanORmedian) 
          ELSE
            ALLOCATE(help3(stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3)))
            CALL writeXmedian_to_3D(help3,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),k,assimX,meanORmedian) 
            stateinfo(i)%svpoint%d3i = NINT(help3)
            DEALLOCATE(help3)
          ENDIF
        ELSEIF(stateinfo(i)%ndim==4)THEN
          CALL writeXmedian_to_4D(stateinfo(i)%svpoint%d4,stateinfo(i)%dims(1),stateinfo(i)%dims(2),stateinfo(i)%dims(3),stateinfo(i)%dims(4),k,assimX,meanORmedian) 
        ENDIF
      ENDIF
    ENDDO

    !Variables for river flows (module variables in HYPEVARIABLES)
    IF(ALLOCATED(Qmax))    CALL writeXmedian_to_2D(Qmax,2,nsub,k,assimX,meanORmedian)
    IF(ALLOCATED(Q2max))   CALL writeXmedian_to_2D(Q2max,2,nsub,k,assimX,meanORmedian)
    IF(ALLOCATED(iQmax))   CALL writeXmedian_to_2D_real2int(iQmax,2,nsub,k,assimX,meanORmedian)
    IF(ALLOCATED(iQ2max))  CALL writeXmedian_to_2D_real2int(iQ2max,2,nsub,k,assimX,meanORmedian)
    IF(ALLOCATED(accdiff)) CALL writeXmedian_to_1D(accdiff,nsub,k,assimX,meanORmedian)
    
    !2) auxiliary variables from A ensemble
    IF(nA.GT.0)THEN
      IF(meanORmedian)THEN
        DO i=1,nA
          DO j=1,nsub
           outvar(j,assimA(i)%info%varID) = assimA(i)%x%outmean(j)
          ENDDO
        ENDDO
      ELSE
        DO i=1,nA
          DO j=1,nsub
           outvar(j,assimA(i)%info%varID) = assimA(i)%x%outquant(2,j)
          ENDDO
        ENDDO
      ENDIF
    ENDIF 

  END SUBROUTINE meanORmedian_to_model2

  !>\brief Get the median ensemble member or the mean of the ensemble members from the state ensembles
  SUBROUTINE writeXmedian_to_1D(x1d,ni,k,assimVar,meanORmedian)
    LOGICAL :: meanORmedian
    TYPE(assim_state_ensemble_type)  :: assimVar(:)   !<ensemble data
    REAL, intent(OUT) :: x1d(:)
    INTEGER, intent(IN) :: ni
    INTEGER, intent(INOUT) :: k
    INTEGER i
    DO i=1,ni
      IF(meanORmedian)THEN
        x1d(i) = assimVar(k)%x%outmean(i) !enkfX%outmean(k+i)
      ELSE
        x1d(i) = assimVar(k)%x%outquant(2,i) !enkfX%outquant(2,k+i)
      ENDIF
    ENDDO
    k=k+1
    !k=k+ni
  END SUBROUTINE writeXmedian_to_1D

  !>\brief Get the median ensemble member or the mean of the ensemble members from the state ensembles
  SUBROUTINE writeXmedian_to_2D(x2d,ni,nj,k,assimVar,meanORmedian)
    LOGICAL :: meanORmedian
    TYPE(assim_state_ensemble_type)  :: assimVar(:)   !<ensemble data
    REAL, intent(OUT) :: x2d(:,:)
    INTEGER, intent(IN) :: ni,nj
    INTEGER, intent(INOUT) :: k
    INTEGER i,j
    DO i=1,ni
      DO j=1,nj
        IF(meanORmedian)THEN
          x2d(i,j) = assimVar(k)%x%outmean(j) !enkfX%outmean(k+(i-1)*nj + j)
        ELSE
          x2d(i,j) = assimVar(k)%x%outquant(2,j) !enkfX%outquant(2,k+(i-1)*nj + j)
        ENDIF
      ENDDO
      k=k+1
    ENDDO
    !k=k+ni*nj
  END SUBROUTINE writeXmedian_to_2D

  SUBROUTINE writeXmedian_to_2D_real2int(x2d,ni,nj,k,assimVar,meanORmedian)
    LOGICAL :: meanORmedian
    TYPE(assim_state_ensemble_type)  :: assimVar(:)   !<ensemble data
    INTEGER, intent(OUT) :: x2d(:,:)
    INTEGER, intent(IN) :: ni,nj
    INTEGER, intent(INOUT) :: k
    INTEGER i,j
    DO i=1,ni
      DO j=1,nj
        IF(meanORmedian)THEN
          x2d(i,j) = int(assimVar(k)%x%outmean(j)) !int(enkfX%outmean(k+(i-1)*nj + j))
        ELSE
          x2d(i,j) = int(assimVar(k)%x%outquant(2,j)) !int(enkfX%outquant(2,k+(i-1)*nj + j))
        ENDIF
      ENDDO
      k=k+1
    ENDDO
    !k=k+ni*nj
  END SUBROUTINE writeXmedian_to_2D_real2int

  !>\brief Get the median ensemble member or the mean of the ensemble members from the state ensembles
  SUBROUTINE writeXmedian_to_3D(x3d,ni,nj,nl,k,assimVar,meanORmedian)
    LOGICAL :: meanORmedian
    TYPE(assim_state_ensemble_type)  :: assimVar(:)   !<ensemble data
    REAL, intent(OUT) :: x3d(:,:,:)
    INTEGER, intent(IN) :: ni,nj,nl
    INTEGER, intent(INOUT) :: k
    INTEGER i,j,l
    DO i=1,ni
      DO j=1,nj
        DO l=1,nl
          IF(meanORmedian)THEN
            x3d(i,j,l) = assimVar(k)%x%outmean(l) !enkfX%outmean(k+(i-1)*nj*nl + (j-1)*nl + l)
          ELSE
            x3d(i,j,l) = assimVar(k)%x%outquant(2,l) !enkfX%outquant(2,k+(i-1)*nj*nl + (j-1)*nl + l)
          ENDIF
        ENDDO
        k=k+1
      ENDDO
    ENDDO
    !k=k+ni*nj*nl
  END SUBROUTINE writeXmedian_to_3D

  !>\brief Get the median ensemble member or the mean of the ensemble members from the state ensembles
  SUBROUTINE writeXmedian_to_4D(x4d,ni,nj,nl,nm,k,assimVar,meanORmedian)
    LOGICAL :: meanORmedian
    TYPE(assim_state_ensemble_type)  :: assimVar(:)   !<ensemble data
    REAL, intent(OUT) :: x4d(:,:,:,:)
    INTEGER, intent(IN) :: ni,nj,nl,nm
    INTEGER, intent(INOUT) :: k
    INTEGER i,j,l,m
    DO i=1,ni
      DO j=1,nj
        DO l=1,nl
          DO m=1,nm
            IF(meanORmedian)THEN
              x4d(i,j,l,m) = assimVar(k)%x%outmean(m) !enkfX%outmean(k+(i-1)*nj*nl*nm + (j-1)*nl*nm + (l-1)*nm + m)
            ELSE
              x4d(i,j,l,m) =  assimVar(k)%x%outquant(2,m) !enkfX%outquant(2,k+(i-1)*nj*nl*nm + (j-1)*nl*nm + (l-1)*nm + m)
            ENDIF
          ENDDO
          k=k+1
        ENDDO
      ENDDO
    ENDDO
    !k=k+ni*nj*nl*nm
  END SUBROUTINE writeXmedian_to_4D

	!-----------------------------------------------------------------------------------------------------------
	!>\brief Generate observation ensemble and calculate its localization.
	!
	!>Purpose: To check if there are available observations, and if so:
	!> - generate the observation data ensembles
	!> - calculate and populate corresponding localization matrices  
	!
	!>The localization is based on x,y, and z coordinates of the subbasins where the observation is located 
	!> - however, for discharge data we should actually take into account the entire upstream area somehow.  
	!
	!Todo: Modify to enable Asynchronic EnKF: that is, to save observations for analysis at a later timestep.
	!      Possible solution: keep several D-ensembles, one for each time step in the analysis period.
	!-----------------------------------------------------------------------------------------------------------
  SUBROUTINE generate_observation_ensemble(assimData)
    !Input arguments
    TYPE(assim_data_type)              :: assimData !<main assimilation variable containing all data
    !Local variables
    INTEGER i,ndata,ndataold,k,j,i2,j2,l
!    INTEGER modobsCol
    LOGICAL obsOK
    REAL    obsValue
    !variables for localization
    REAL ix,iy,iz,jx,jy,jz,xy_dist,z_dist
    !REAL,allocatable :: xini(:)

    !1) Check the available data
    ndata = 0
    assimData%info%nD = 0  
    !ssimData%D%nvar = 0   

    IF(assimData%info%nE.LE.1)return
    !1.1 loop over observation types
    DO k=1,assimData%info%nObs
      !set number of non-missing to 0
      assimData%Obs(k)%gen%ndata = 0
      !1.2 loop over subbasins
      DO i=1,nsub
        ! check available data for the requesed data types defined by the following rules:
        !   a. EnkfDinfo%xi corresponds to the observation's outvar index
        !   b. EnkfDinfo%mi corresponds to the model variable's outvar index
        !   c. EnkfDinfo%qt defines the quality screening type
        !   d. EnkfDinfo%ql defines the quality screening limiting value

        !!initialize missing value      !CP170202 changed so observations are collected from outvar, works for regional observations
        !obsValue = -9999.
        !!Select data source (Xobs or Qobs)
        !IF(assimData%Obs(k)%info%varID.NE.o_rout)THEN      !rout outvar index is o_rout (was 54)
        !!IF(assimData%Obs(k)%info%varID.NE.54)THEN      !COUT OUTVAR INDEX IS 54  CP161128 replaced with above, made safer
        !	!Search Xobs for data
        !	IF(xobsindex(assimData%Obs(k)%info%varID,i)>0)THEN
        !		obsValue = xobsi(xobsindex(assimData%Obs(k)%info%varID,i))
        !	ENDIF
        !ELSE
        !	!Search Qobs for data
        !	IF(ALLOCATED(qobsi))THEN
        !		obsValue = qobsi(i)
        !	ENDIF
        !ENDIF
        !Get observation data from outvar
        obsValue = outvar(i,assimData%Obs(k)%info%varID)
        !Discard missing values, or values out of bound
        IF(obsValue.GT.-9998.)THEN
           obsOK = .TRUE.
           !!Select quality control option
           !SELECT CASE(EnkfDinfo%qt(k))
           !	CASE(1)      ! >=
           !		IF(obsValue.GE.EnkfDinfo%ql(k))obsOK=.TRUE.
           !	CASE(2)      ! >
           !		IF(obsValue.GT.EnkfDinfo%ql(k))obsOK=.TRUE.
           !	CASE(-1)     ! <=
           !		IF(obsValue.LE.EnkfDinfo%ql(k))obsOK=.TRUE.
           !	CASE(-2)     ! <=
           !		IF(obsValue.LT.EnkfDinfo%ql(k))obsOK=.TRUE.
      !       CASE DEFAULT ! non-missing
      !          obsOK = .TRUE.
           !END SELECT
        ELSE
          obsOK = .FALSE.
        ENDIF
        !use the value if quality control has passed
        IF(obsOK)THEN
          ndata = ndata + 1
          assimData%Obs(k)%gen%ndata = assimData%Obs(k)%gen%ndata + 1         
          assimData%Obs(k)%gen%mean(i)=obsValue
          !	EnkfD%mean(ndata)    = obsValue
          !	EnkfD%mi(ndata)      = EnkfDinfo%mi(k)     ! outvar index of the model observation
          !	EnkfD%ensgen(ndata)  = EnkfDinfo%ensgen(k) ! ensemble generation type
        
          !!subbasin index used by modelobservation "H function"
          !EnkfD%xi(ndata)       = i
          !!minimum, maximum, variance, and other ensemble options
          !EnkfD%maximum(ndata)  = EnkfDinfo%maximum(k)
          !EnkfD%minimum(ndata)  = EnkfDinfo%minimum(k)
          !EnkfD%sigma(ndata)    = EnkfDinfo%sigma(k)
          !EnkfD%ensgen(ndata)   = EnkfDinfo%ensgen(k)
          !EnkfD%restmeta(ndata) = EnkfDinfo%restmeta(k)
          !EnkfD%semimeta(ndata) = EnkfDinfo%semimeta(k)
        ELSE
          assimData%Obs(k)%gen%mean(i)=-9999.
        ENDIF
        ndataold = ndata        
      ENDDO
    ENDDO
    assimData%info%nD = ndata  !number of observations this time step
		!assimData%D%nvar = ndata
		!EnkfInfo%nDA = ndata
		!EnkfD%n = ndata
    
    !Check number of observations, continue to ensemble generation if ndata> 0
    IF(assimData%info%nD.GT.0)THEN
      !Ensemble generation with the general function

      !First, re-allocate the D, HX, and R ensemble
      !IF(allocated(xini))deallocate(xini)
      !ALLOCATE(xini(assimData%info%nD))
      !xini = 0.
      !CALL allocate_assim_ensemble(assimData%D,assimData%info%nE,assimData%info%nD, assimData%D%fileID, 1, 1, xini, assim_defmin, assim_defmax, .true., -9999.)
      !CALL allocate_assim_ensemble(assimData%HX,assimData%info%nE,assimData%info%nD, assimData%HX%fileID, 1, 1, xini, assim_defmin, assim_defmax, .true., -9999.)
      !CALL allocate_assim_ensemble(assimData%R,assimData%info%nE,assimData%info%nD, assimData%HX%fileID, 1, 1, xini, assim_defmin, assim_defmax, .true., -9999.)
      IF(ALLOCATED(assimData%D))DEALLOCATE(assimData%D)
      IF(ALLOCATED(assimData%HX))DEALLOCATE(assimData%HX)
      IF(ALLOCATED(assimData%R))DEALLOCATE(assimData%R)
      ALLOCATE(assimData%D(assimData%info%nD,assimData%info%ne))
      ALLOCATE(assimData%HX(assimData%info%nD,assimData%info%ne))
      ALLOCATE(assimData%R(assimData%info%nD))
      !initialize R to 0, then fill diagonal with sigma**2 further down
      assimData%R=0.
      
      !  2nd argument ("vectorvise") = .true. means 
      !  that ensembles will be spatially un-correlated
      CALL generate_input_ensemble(assimData%info%nObs,assimData%Obs)

      !Localization matrix generation
      !   re-allocate the localization matrices used by the enkf analysis routine
      IF(ALLOCATED(assimData%locCYY))deallocate(assimData%locCYY) 
      ALLOCATE(assimData%locCYY(assimData%info%nD,assimData%info%nD))
      DO l=1,assimData%info%nloc
        IF(ALLOCATED(assimData%locCXY(l)%x)) DEALLOCATE(assimData%locCXY(l)%x) 
        ALLOCATE(assimData%locCXY(l)%x(assimData%Coordinates(l)%n,assimData%info%nD))   !ncoord=nloc => l=l
      ENDDO

      !Fill localization matrices with values from the basin2basin localization matrix
      !for all observations Y...
      !loop over observation types
      
      DO k=1,assimData%info%nObs
        i2=0
        !loop over subbasins
        DO i=1,nsub
          !check available data
          IF(assimData%Obs(k)%gen%mean(i).GT.-9998.)THEN
            i2 = i2 + 1
            
            !... sigma*sigma as observation error variance
            assimData%R(i2)=assimData%Obs(k)%gen%sigma(i)**2
            
            !...coordinates of the observation Y(i) (subbasin index given by i)
            ix = assimData%Coordinates(1)%x(i)
            iy = assimData%Coordinates(1)%y(i)
            iz = assimData%Coordinates(1)%z(i)
            !ix = basin(i)%xcoord   !CP161216 use Coordinates(1)(=subbasin) for observations
            !iy = basin(i)%ycoord
            !iz = basin(i)%elev
            
         !   !...loop over the first 1:nsub model states and fill the CXY[1:nsub,1:nD) localization matrix
            !DO j=1,nsub
              !!coordinates of subbasin j
              !jx = basin(j)%xcoord               !CP161214 here was only one locCXY which was based on nsub and basin-coordinates/elevation
              !jy = basin(j)%ycoord               !CP161214 now we need one for nsub and one for aquifers and I loop over coordinate systems
              !jz = basin(j)%elev                 !CP161214 locCYY is still only one, but also observations may need different coordinates systems (eg. cout)
              !!localization value
         !     xy_dist = sqrt((jx-ix)**2 + (jy-iy)**2)
              !z_dist  = abs(jz-iz)
              !assimData%locCXY(j,i2)= exp(-( 2. * (xy_dist/assimData%info%xy_scalefac)**2 )) * exp(-( 2. * (z_dist/assimData%info%z_scalefac)**2 ))
         !   ENDDO
            !...loop over the number of coordinate systems (i.e. subbasins and aquifer) and set localization matrix for each
            DO l=1,assimData%info%nloc
              !...loop over the number of different coordinates in this coordinate system (nsub or naquifers) and fill the CXY localization matrix
              DO j=1,assimData%Coordinates(l)%n
                !get current coordinates
                jx = assimData%Coordinates(l)%x(j)
                jy = assimData%Coordinates(l)%y(j)
                jz = assimData%Coordinates(l)%z(j)
                !calculate localization value
                xy_dist = sqrt((jx-ix)**2 + (jy-iy)**2)
                z_dist  = abs(jz-iz)
                assimData%locCXY(l)%x(j,i2)= exp(-( 2. * (xy_dist/assimData%info%xy_scalefac)**2 )) * exp(-( 2. * (z_dist/assimData%info%z_scalefac)**2 ))
              ENDDO
            ENDDO

            !...loop over all other observations Y and fill the CYY[nD,nD) localization matrix
            j2=0
            DO j=1,nsub
              IF(assimData%Obs(k)%gen%mean(j).GT.-9998.)THEN
                j2=j2+1
                !coordinates of the observation Y(j)
                jx = assimData%Coordinates(1)%x(j)
                jy = assimData%Coordinates(1)%y(j)
                jz = assimData%Coordinates(1)%z(j)
                !jx = basin(j)%xcoord
                !jy = basin(j)%ycoord
                !jz = basin(j)%elev
                !localization value at [i,j]
                xy_dist = sqrt((jx-ix)**2 + (jy-iy)**2)
                z_dist  = abs(jz-iz)
                assimData%locCYY(j2,i2)= exp(-( 2. * (xy_dist/assimData%info%xy_scalefac)**2 )) * exp(-( 2. * (z_dist/assimData%info%z_scalefac)**2 ))
              ENDIF
            ENDDO
            
            !Finally, also write the ensemble to the D-matrix
            assimData%D(i2,:) = assimData%Obs(k)%x%x(i,:)
          ENDIF
        ENDDO
      ENDDO
    ENDIF 
  END SUBROUTINE generate_observation_ensemble

	!---------------------------------------------------------------------------------------
	!>\brief Model equivalents are transfered to the observation ensemble by the observation operator.
	!
	!>Purpose: This is the observation operator, aka the "H function", which fills the 
	!>         modelled observations ensemble (HX) with model data corresponding to the 
	!>         available observations (the D matrix).
	!>
	!>Actions: Based on the observations in the current D matrix (EnkfD) the corresponding
	!>         model variables are derived from the current outvar. 
	!>         So far, no additional transformations have been necessary, but it is straight 
	!>         forward to include whatever operator is needed in this function.
	!---------------------------------------------------------------------------------------
  SUBROUTINE modelobservations_to_ensemble(i_ens,assimData)
  !Input arguments
    TYPE(assim_data_type)   :: assimData !<main assimilation variable containing all data
    INTEGER                 :: i_ens     !<current ensemble
    !Local variables
    INTEGER i, j, k
    ! step 1, prepare the model observation matrix "HX", but only if any obs.data is available
    j=0
    IF(assimData%info%nD.GT.0)THEN ! check the current number of data in the observation matrix D
      !loop over observation types
      DO k=1,assimData%info%nObs
        !loop over subbasins
        DO i=1,nsub
          !check available data
          IF(assimData%Obs(k)%gen%mean(i).GT.-9998.)THEN
            j=j+1
            assimData%HX(j,i_ens) = outvar(i,assimData%Obs(k)%info%modID)
          ENDIF
        ENDDO
      ENDDO
      
      !DO i=1,EnkfInfo%nDA
      !	! EnkfD%mi(i)      tells us what model variable we look for (outvar-column)
      !	! EnkfD%xi(i)      tells us what subid we look for          (outvar-row)
      !	! EnkfHX%ensgen(i) tells us what H function to use
      !	SELECT CASE(EnkfHX%ensgen(i))
      !		CASE(0) ! no transformation, use OUTVAR values directly
      !			enkfHX%y(i,i_ens) = outvar(EnkfD%xi(i),EnkfD%mi(i))
      !		CASE DEFAULT
      !			enkfHX%y(i,i_ens) = outvar(EnkfD%xi(i),EnkfD%mi(i))
      !	END SELECT
      !ENDDO
    ENDIF
  END SUBROUTINE modelobservations_to_ensemble

  !-------------------------------------------------------------------------
  !>\brief Routine to generate forcing data ensembles
  !--------------------------------------------------------------------------
  SUBROUTINE generate_forcing_ensemble(assimData)
    !Input arguments
    TYPE(assim_data_type)   :: assimData !<main assimilation variable containing all data
    !Local
    INTEGER i, j
    ! a) forcing from input vectors to ensemble mean:
    !Loop over forcing data types
    DO i = 1, assimData%info%nF
      !Loop over subbasins
      DO j=1,nsub
        IF(assimData%F(i)%info%varID.EQ.1)THEN     !precipitation
          assimData%F(i)%gen%mean(j) = preci(j)
        ELSEIF(assimData%F(i)%info%varID.EQ.2)THEN !mean temperature
          assimData%F(i)%gen%mean(j) = tempi(j)
        ELSEIF(assimData%F(i)%info%varID.EQ.3)THEN !min temperature
          assimData%F(i)%gen%mean(j) = tmini(j)
        ELSEIF(assimData%F(i)%info%varID.EQ.4)THEN !max temperature
          assimData%F(i)%gen%mean(j) = tmaxi(j)
        ENDIF
      ENDDO
    ENDDO
    ! b) generate input ensembles with general function:
    CALL generate_input_ensemble(assimData%info%nF,assimData%F)
    ! c) add check to make sure tmin<=tmean<=tmax (or better make perturbation of tmean - tmin and tmax-tmean, respectively)

  END SUBROUTINE generate_forcing_ensemble

  END MODULE assimilation_interface