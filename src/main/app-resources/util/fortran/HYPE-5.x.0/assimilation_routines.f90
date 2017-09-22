!> \file assimilation_routines.f90
!> Contains module assimilation_routines, with model independent subroutines and functions used for data assimilation.
!  Author: D.Gustafsson (SMHI)
!  Versions: 
!  2011-2015:  Developed by D.Gustafsson (SMHI/KTH) and J.Ahlberg (KTH). The latest version from 2015-09-20 was used with the HOPE model for a "HUVA" snow data assimilation project 2013-2015.
!  2016.10.20: Module name changed to ASSIMILATION_ROUTINES and additional code convention adaptations for new implementation in HYPE_4_12_+.
!  2017.04.27: Replace all code dependent on Intel MKL library with intrinsic functions and other new/old code.
  
!> Generic subroutines and functions for (Ensemble Kalman Filter) Data Assimilation.
!> These functions are supposed to be general, and should not be changed
!> Functions special for a specific model application are in assimilation_interface.f90

MODULE ASSIMILATION_ROUTINES
!Copyright 2016-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------------

  USE ASSIMILATION_VARIABLES
  USE RANDOM_ROUTINES
  
  IMPLICIT NONE

CONTAINS
!-----------------------------------------------------------------------------------
!VARIOUS ENSEMBLE DATA MANIPULATION ROUTINES
!-----------------------------------------------------------------------------------
  
!<\brief Re-initilize selected ensembles to mean or median to avoid numerical issues.
  SUBROUTINE meanORmedian_to_ensemble(nx,assimX,meanORmedian)

    !Argument declations
    TYPE(assim_state_ensemble_type) :: assimX(:)   !<ensemble data
    INTEGER, INTENT(IN)       :: nx           !<number of variables
    LOGICAL, INTENT(IN)       :: meanORmedian !<flag for setting mean (T) or median (F)

    INTEGER i,j,k
    !loop over state ensembles
    DO k=1,nx
      IF(.not.assimX(k)%x%assimilate)THEN
        !loop over ensemble members
        DO j=1,assimX(k)%x%nens
          IF(ALLOCATED(assimX(k)%x%x))THEN
            !if matrix allocated...
            !loop over model units (often subbasins)
            DO i=1,assimX(k)%x%nvar
              IF(meanORmedian)THEN
                assimX(k)%x%x(i,j) = assimX(k)%x%outmean(i)
              ELSE
                ! switch to "outmedian"
                assimX(k)%x%x(i,j) = assimX(k)%x%outquant(2,i)
              ENDIF
            ENDDO
          ELSE
            !else..  write model units to bin-file, mean or median
            IF(meanORmedian)THEN
              WRITE(assimX(k)%x%fileID,REC=assimX(k)%x%rec+j) assimX(k)%x%outmean(:)
            ELSE
                ! switch to "outmedian"
              WRITE(assimX(k)%x%fileID,REC=assimX(k)%x%rec+j) assimX(k)%x%outquant(2,:)
            ENDIF
         ENDIF
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE meanORmedian_to_ensemble
    
!>Initializes the information variable for assimilation
!---------------------------------------------------------
  SUBROUTINE initialize_assim_info(assimInfo)
    !INPUT ARGUMENT
    TYPE(assim_info_type), INTENT(INOUT) :: assimInfo !<information on assimilation simulation
  
    !Initialize random number generation seed (NOT USED ANYMORE!!!)
    !CALL time(assimInfo%seed)    !the random seed will be different, maybe make it possible for user to decide via the enkfinfo.txt?
  
    !Initialize various info variables
    assimInfo%nE      = 100       !number of ensemble members      (columns in ensemble matrices)
    assimInfo%nX      = 0         !number of state variables       (length of ensemble vector EnkfX)
    assimInfo%nA      = 0         !number of auxiliary variables   (length of ensemble vector EnkfA)
    assimInfo%nF      = 0         !number of forcing variables     (length of ensemble vector EnkfF)
    assimInfo%nObs    = 0         !number of observation variables (rows of ensemble vector EnkfObs(:,:))
    assimInfo%nD      = 0         !number of observations for next Enkf analysis (rows in EnkfD%x%x(:,:))
  !  assimInfo%nDT     = 0         !number of observation timesteps (colss of ensemble vector EnkfObs(:,:))
  !  assimInfo%nP      = 0         !number of parameters            (length of ensemble vector EnkfP)
    assimInfo%nloc    = 0         !number of localization matrices (length of vector EnkfLocCXY)
    assimInfo%ncoord  = 0         !number of spatial domains       (length of vector EnkfCoordinates)
    
    !LOGICAL flags, may be modified...
    assimInfo%FA      = .false.   !include auxil. in kalman filter     (general switch on/off)
    assimInfo%FP      = .false.   !include parameters in kalman filter (general switch on/off) 
    assimInfo%FF      = .false.   !include inputs in kalman filter     (general switch on/off)
    !EnkfInfo%XS      = .false.   !include X states in statistical output (general switch on/off)
    !assimInfo%EC      = .false.   !ECMWF forecasts (general switch on/off)
    assimInfo%meanout = .true.    !ensemble mean(.true.) or median (.false.) in output
    ALLOCATE(assimInfo%assim_flag(assimInfo%nCat+50)) !assimilate state variables (switch for categories and variables) (the size is unecessary large I think)

    !some Ensemble Kalman Filter parameters
    !assimInfo%moradkhani_delta = 0.95 ! coef.[0-1] to retain variance in parameter ensemble (Moradkhani et al, 2004), value typical around 0.95 !CP170615 not used
    !assimInfo%ensgen_minsigma = 1.e-5
    
    !switch on/off read/write ensemble data to binary files
    assimInfo%useBinFilesX = 0 !.FALSE.  CP170505 these are integers
    assimInfo%useBinFilesFA = 0 !.FALSE.
    assimInfo%nBinFiles = 0
  
    !Output of statistical simulation results
    assimInfo%nstatout = 0

    !localization parameters
    assimInfo%xy_scalefac = 1000000.
    assimInfo%z_scalefac  = 100000.
    
    !Options for non-controlled variables
    assimInfo%collapseNonControlled = .false.
    
    !Options for State Ensemble Initialization
    assimInfo%initializeFromBinFiles = .false.
    
  END SUBROUTINE initialize_assim_info

!>Allocate ensemble vectors needed for the data assimilation application
!----------------------------------------------------------------------------------------------  
  SUBROUTINE allocate_assim_ensemble_vectors(assimData,fid_0)
    TYPE(assim_data_type), INTENT(INOUT) :: assimData !<main assimilation variable containing all data
    INTEGER :: fid_0  !<base fileunit for bin files to be used
    !INTEGER :: nx,na,np,nf,nobs,ncoord,nloc,nobsDT
    INTEGER :: i
    !--------------------------------------------------------------------------------------------
    !X state variable ensemble vector
    IF(ALLOCATED(assimData%X))DEALLOCATE(assimData%X)
    IF(assimData%info%nx.GT.0)ALLOCATE(assimData%X(assimData%info%nx))
    
    !A auxiliary variable ensemble vector
    IF(ALLOCATED(assimData%A))DEALLOCATE(assimData%A)
    IF(assimData%info%na.GT.0)ALLOCATE(assimData%A(assimData%info%na))
    
    !  !P parameter ensemble vector
    !  IF(ALLOCATED(assimData%P))DEALLOCATE(assimData%P)
    !  IF(assimData%info%np.GT.0)ALLOCATE(assimData%P(assimData%info%np))
    
    !Obs observation ensemble vector
    IF(ALLOCATED(assimData%Obs))DEALLOCATE(assimData%Obs)
    !  IF(assimData%info%nobs.GT.0)ALLOCATE(assimData%Obs(nobs,nobsDT))
    IF(assimData%info%nobs.GT.0)ALLOCATE(assimData%Obs(assimData%info%nobs))
    
    !F forcing variable ensemble vector
    IF(ALLOCATED(assimData%F))DEALLOCATE(assimData%F)
    IF(assimData%info%nf.GT.0)ALLOCATE(assimData%F(assimData%info%nf))
    
    !Spatial Coordinate Domain vector
    !IF(ALLOCATED(assimData%Coordinates))DEALLOCATE(assimData%Coordinates)
    !IF(assimData%info%ncoord.GT.0)ALLOCATE(assimData%Coordinates(assimData%info%ncoord))
    
    !Localization (locCXY) vector
    !IF(ALLOCATED(assimData%LocCXY))DEALLOCATE(assimData%LocCXY)
    !IF(assimData%info%nloc.GT.0)ALLOCATE(assimData%LocCXY(assimData%info%nloc))

    !D, HX and R and LocCYY are not vectors, and DO not need to be ALLOCATED, the length is decided for each analysis time window
    !Further more, R is not an ensemble but an nD x nD matrix, and if observation errors are assumed uncorrelated, R is in fact a diagonal matrix and can be saved as a nD x 1 vector
  
    !If binary files are used
    !  IF(assimData%info%useBinFiles)THEN
      !nState + nAux + nPar + nForc + nObs*nObsTimeSteps + 3 (D + HX + R)
      !assimData%info%nBinFiles = assimData%info%nX + assimData%info%nA + assimData%info%nP + assimData%info%nF + assimData%info%nObs*assimData%info%nObsDT + 3
 
    assimData%info%nBinFiles = assimData%info%nX + assimData%info%nA + assimData%info%nF + assimData%info%nObs + 3
    
    IF(ALLOCATED(fid_assim_bin))DEALLOCATE(fid_assim_bin)
    ALLOCATE(fid_assim_bin(assimData%info%nBinFiles))
    
    DO i=1,assimData%info%nBinFiles
      fid_assim_bin(i) = fid_0 + i 
    ENDDO
      !ENDIF
    END SUBROUTINE allocate_assim_ensemble_vectors

  !--------------------------------------------------------
  !>Allocate and initialize an assim_ensemble_type variable
  !--------------------------------------------------------
  SUBROUTINE allocate_assim_ensemble(assimVar, nens, nvar, varID, fileID, useBinFile, locID, coordID, xini, mini, maxi, allocateOutput, missing,assimilate,iniFromBin)
    !INPUT ARGUMENTS
    TYPE(assim_ensemble_type) :: assimVar !<ensemble data
    INTEGER :: nens      !<number of ensemble members
    INTEGER :: nvar      !<number of variables (model units, i.e. n:o subbasins in HYPE) 
    INTEGER, INTENT(IN) :: varID  !<to set record for binary file
    INTEGER :: fileID         !<file unit ID for direct access binary file I/O
    INTEGER :: locID
    INTEGER :: coordID  !<ID for coordinate system for this variable
    REAL    :: xini(nvar)     !<initial values
    REAL    :: mini, maxi     !max/min thresholds for truncation of unrealistic values
    LOGICAL :: allocateOutput !<flag for output variable allocation
    REAL    :: missing        !<initial output values
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
    LOGICAL, INTENT(IN) :: iniFromBin

    !LOCAL
    INTEGER :: j,reclen
    CHARACTER(LEN=10) :: filename
    
    !ASSIGN AND ALLOCATE or open bin-files for saving ensembles (CP added bin-files CP161201)
    assimVar%nvar = nvar
    assimVar%nens = nens
    assimVar%minimum = mini
    assimVar%maximum = maxi
    assimVar%fileID = fileID
    assimVar%locID = locID
    assimVar%coordID = coordID
    IF(useBinFile==1)THEN
      !If one bin-file is used, write initial data to already open file (unless we should use existing file for initialization)
      assimVar%rec = (varID-1)*nens
      IF(.NOT.iniFromBin)THEN
        DO j=1,nens
          WRITE(assimVar%fileID,REC=assimVar%rec+j) xini
        ENDDO
      ENDIF
    ELSEIF(useBinFile==2)THEN
      !If several bin-files is used open file and write initial data
      INQUIRE(IOLENGTH=reclen) xini   !determine suitable record length here, because ifort och gfortran have different file storage unit (i.e. RECL)
      WRITE(filename,'(I5.5,a4)') assimVar%fileID,'.bin'
      OPEN(UNIT=assimVar%fileID,FILE=TRIM(filename),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=reclen)
      IF(.NOT.iniFromBin)THEN
        DO j=1,nens
          WRITE(assimVar%fileID,REC=j) xini
        ENDDO
      ENDIF
      assimVar%rec = 0
    ELSE
      !If no bin-files, then allocate a matrix to hold the ensemble and initialize it.
      assimVar%rec = 0
      IF(ALLOCATED(assimVar%x))DEALLOCATE(assimVar%x)
      ALLOCATE(assimVar%x(nvar,nens))
      DO j=1,nens
        assimVar%x(:,j)=xini
      ENDDO
    ENDIF
  
    IF(allocateOutput)THEN
      IF(ALLOCATED(assimVar%outmean))DEALLOCATE(assimVar%outmean)
      ALLOCATE(assimVar%outmean(nvar))
      assimVar%outmean=missing
      IF(ALLOCATED(assimVar%outquant))DEALLOCATE(assimVar%outquant)
      !ALLOCATE(assimVar%outquant(nvar,3))    !CP161208 This is opposit how the array is used in the code!
      ALLOCATE(assimVar%outquant(3,nvar))
      assimVar%outquant=missing
      IF(ALLOCATED(assimVar%outmin))DEALLOCATE(assimVar%outmin)
      ALLOCATE(assimVar%outmin(nvar))
      assimVar%outmin=missing
      IF(ALLOCATED(assimVar%outmax))DEALLOCATE(assimVar%outmax)
      ALLOCATE(assimVar%outmax(nvar))
      assimVar%outmax=missing
      IF(ALLOCATED(assimVar%outsigma))DEALLOCATE(assimVar%outsigma)
      ALLOCATE(assimVar%outsigma(nvar))
      assimVar%outsigma=missing
    ENDIF
  
    assimVar%assimilate = assimilate
  END SUBROUTINE allocate_assim_ensemble

!--------------------------------------------------------
!>Allocate and initialize a assim_interface_type variable
!--------------------------------------------------------
  SUBROUTINE allocate_assim_interface(assimVar,varName,varID,modID,nSubDim,subDimID)
    !INPUT ARGUMENTS
    TYPE(assim_interface_type)  :: assimVar !<interface data
    CHARACTER(LEN=*)     :: varName        !<character string for model variables (used for debugging, and possibly file names and outputs)
    INTEGER              :: varID          !<variable ID (id number used by interface for linking to model variables)  (in HYPE it can be an outvar index, or the order in which the state variables are considered by interface)
    INTEGER              :: modID          !<model ID,  link to the corresponding variables used for H(X)              (in HYPE: outvar index)
    INTEGER              :: nSubDim        !<number of sub-dimensions (if needed, for instance lateral sub-units, vertical layers, substances, etc, in HYPE for instance SLC, substances, landuse, or slc, etc)
    INTEGER              :: subDimID(:)    !<index in the sub-dimensions
    !ALLOCATE AND ASSIGN DATA
    assimVar%varName(1:(MIN(30,LEN_TRIM(varName)))) = TRIM(varName)  
    assimVar%varID   = varID
    assimVar%modID   = modID
    assimVar%nSubDim = nSubDim
    IF(ALLOCATED(assimVar%subDimID))DEALLOCATE(assimVar%subDimID)
    IF(nSubDim.GT.0)THEN
      ALLOCATE(assimVar%subDimID(nSubDim))
      assimVar%subDimID(1:nSubDim)=subDimID(1:nSubDim)
    ENDIF
  END SUBROUTINE allocate_assim_interface

!---------------------------------------------------------
!>Allocate and initialize an enkf_generation_type variable
!---------------------------------------------------------
  SUBROUTINE allocate_assim_generation(assimVar,nvar,ensgen,fixsigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !,ECMWF)
    !INPUT ARGUMENTS
    TYPE(assim_generation_type) :: assimVar !<generation_data variable
    INTEGER                     :: nvar     !<number of variables, ie n:o "model units" (for instance, number of sub-basins)
    INTEGER                     :: ensgen   !<type of ensemble generation        (0 none, 1 unrestricted, 2 [min,+inf], 3 [-inf,max], 4 restricted [min,max])   
    REAL                        :: fixsigma !<fixed standard deviation           (ensgen=1)
    REAL                        :: semimeta !<relative sigma for semi-restricted (ensgen=2,3,4, following Turner et al 2008)
    REAL                        :: restmeta !<relative sigma for restricted      (ensgen=2,3,4)
    REAL                        :: minsigma !<minimum sigma                      (ensgen=2,3,4)
    REAL                        :: lscale   !<correlation length for spatially correlated perturbation
    REAL                        :: gridsize !<2D gridsize for spatially correlated perturbation generation
    INTEGER                     :: corrtype !<spatial correlation FUNCTION option
    REAL                        :: xcoord(:)!<x coordinate, for spatially correlated perturbations
    REAL                        :: ycoord(:)!<y coordinate, for spatially correlated perturbations
    !LOGICAL                     :: ECMWF
    
    !ASSIGN
    assimVar%nvar = nvar
    assimVar%ensgen = ensgen
    assimVar%fixsigma = fixsigma
    assimVar%semimeta = semimeta
    assimVar%restmeta =restmeta
    assimVar%minsigma = minsigma
    
    !ALLOCATE
    !sigma, standard deviation (input to ensemble generation)
    IF(ALLOCATED(assimVar%sigma))DEALLOCATE(assimVar%sigma)
    ALLOCATE(assimVar%sigma(nvar))
    assimVar%sigma = 0.0
    IF(ensgen.EQ.1)assimVar%sigma=fixsigma
    
    !mean, mean value  (input to ensemble generation)
    IF(ALLOCATED(assimVar%mean))DEALLOCATE(assimVar%mean)
    ALLOCATE(assimVar%mean(nvar))
    assimVar%mean = 0.0
    
    !data structure for spatially correlated random perturbations
    IF(lscale.GT.0. .AND. gridsize.GT.0. .AND. corrtype.GT.0 .AND. corrtype.LT.4)THEN
      !DO randomxy
      assimVar%dorandxy = .TRUE.
      !initialize randxy data structure
      CALL init_randxy_data(assimVar%myrandxy_data,nvar,xcoord,ycoord,lscale,gridsize,corrtype)
      !    assimVar%myrandxy_data%lscale     = lscale
      !    assimVar%myrandxy_data%gridsize   = gridsize
      !    assimVar%myrandxy_data%index_corr = corrtype
    ELSE
      assimVar%dorandxy = .FALSE.
    ENDIF
    
    !assimVar%ECMWF = ECMWF
    
  END SUBROUTINE allocate_assim_generation

!-----------------------------------------------------------------
!Allocate and initialize STATE, auxiliary, forcing and obseervation ensemble variables
!-----------------------------------------------------------------
!>Allocate and initialize auxiliary ensemble variables; 0-dimensional
  SUBROUTINE allocate_auxiliary_ensemble(assimVar,nens,nvar,varName,xini,varID,recID,locID,coordID,fileID,useBinFile,minimum,maximum,assimilate)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar
    INTEGER :: nens,nvar,locID,coordID,fileID,dimID(1),ndim
    INTEGER,INTENT(IN) :: varID   !<outvar-variable index
    INTEGER,INTENT(IN) :: recID   !<record of variable in bin-file
    REAL    :: minimum,maximum
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
  
    ndim=0
    dimID(1)=0
    !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,recID,fileID,useBinFile,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate,.false.)
    !ALLOCATE the interface data
    CALL allocate_assim_interface(assimVar%info,varName,varID,-1,ndim,dimID)
    !update varID
    !varID = varID+1 !not for aux this is o_rout etc.
  END SUBROUTINE allocate_auxiliary_ensemble
  
!>Allocate and initialize state ensemble variables; 0-dimensional
  SUBROUTINE allocate_0dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinFile,minimum,maximum,assimilate,iniFromBin)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar
    INTEGER :: nens,nvar,varID,locID,coordID,fileID,dimID(1),ndim
    REAL    :: minimum,maximum
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
    LOGICAL, INTENT(IN) :: iniFromBin
  
    ndim=0
    dimID(1)=0
    !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
    !CALL allocate_assim_ensemble(assimVar%x,nens,nvar,fileID,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,varID,fileID,useBinFile,locID,coordID,xini,minimum,maximum,.true.,-9999.,assimilate,iniFromBin)
    !ALLOCATE the interface data
    CALL allocate_assim_interface(assimVar%info,varName,varID,-1,ndim,dimID)
    !update varID
    varID = varID+1
  END SUBROUTINE allocate_0dim_state_ensemble
  
!>Allocate and initialize state ensemble variables; 1-dimensional
  SUBROUTINE allocate_1dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinFile,minimum,maximum,n1,assimilate,iniFromBin)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER :: nens,nvar,varID,locID,coordID,fileID(:),n1,i,dimID(1),ndim
    REAL    :: minimum,maximum
    REAL    :: xini(n1,nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
    LOGICAL, INTENT(IN) :: iniFromBin

    ndim=1
    DO i=1,n1
      dimID(1)=i
      !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
      !CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,fileID(varID),locID,coordID,xini(i,:),minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
      CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,varID,fileID(varID),useBinFile,locID,coordID,xini(i,:),minimum,maximum,.true.,-9999.,assimilate,iniFromBin)
      !ALLOCATE the interface data
      CALL allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      !update varID
      varID = varID+1
    ENDDO
  END SUBROUTINE allocate_1dim_state_ensemble
  
!>Allocate and initialize state ensemble variables; 2-dimensional
  SUBROUTINE allocate_2dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinfile,minimum,maximum,n1,n2,assimilate,iniFromBin)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER :: nens,nvar,varID,locID,coordID,fileID(:),n1,n2,i,j,ndim,dimID(2)
    REAL    :: minimum,maximum
    REAL    :: xini(n2,n1,nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
    LOGICAL, INTENT(IN) :: iniFromBin

    ndim=1
    DO j=1,n2
      dimID(1)=j
    DO i=1,n1
      dimID(2)=i
      !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
      !CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,fileID(varID),locID,coordID,xini(j,i,:),minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
      CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,varID,fileID(varID),useBinFile,locID,coordID,xini(j,i,:),minimum,maximum,.true.,-9999.,assimilate,iniFromBin)
      !ALLOCATE the interface data
      CALL allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      !update varID
      varID = varID+1
    ENDDO
    ENDDO
  END SUBROUTINE allocate_2dim_state_ensemble
  
!>Allocate and initialize state ensemble variables; 3-dimensional
  SUBROUTINE allocate_3dim_state_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,useBinFile,minimum,maximum,n1,n2,n3,assimilate,iniFromBin)
    !INPUT
    TYPE(assim_state_ensemble_type)  :: assimVar(:)
    INTEGER :: nens,nvar,varID,locID,coordID,fileID(:),n1,n2,n3,i,j,k,ndim,dimID(3)
    REAL    :: minimum,maximum
    REAL    :: xini(n3,n2,n1,nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-files
    LOGICAL, INTENT(IN) :: assimilate   !<flag for including variable in assimilation
    LOGICAL, INTENT(IN) :: iniFromBin

    ndim=3
    DO k=1,n3
      dimID(1)=k
    DO j=1,n2
      dimID(2)=j
    DO i=1,n1
      dimID(3)=i
      !ALLOCATE and initialize the ensemble matrix data or initialize bin-file
      !CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,fileID(varID),locID,coordID,xini(k,j,i,:),minimum,maximum,.true.,-9999.,assimilate) !CP161201 introduced binfiles
      CALL allocate_assim_ensemble(assimVar(varID)%x,nens,nvar,varID,fileID(varID),useBinFile,locID,coordID,xini(k,j,i,:),minimum,maximum,.true.,-9999.,assimilate,iniFromBin)
      !ALLOCATE the interface data
      CALL allocate_assim_interface(assimVar(varID)%info,varName,varID,-1,ndim,dimID)
      !update varID
      varID = varID+1
    ENDDO
    ENDDO
    ENDDO
  END SUBROUTINE allocate_3dim_state_ensemble
  
!>Allocate and initialize forcing ensemble variables
  SUBROUTINE allocate_assim_forcing_ensemble(assimVar,nens,nvar,varName,xini,varID,locID,coordID,fileID,minimum,maximum, &
    ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord,useBinFile)
    !INPUT
    TYPE(assim_input_ensemble_type) :: assimVar
    INTEGER :: nens,nvar,locID,coordID,fileID
    INTEGER,INTENT(IN) :: varID   !<counter for variables in ensembles (is not all variables varID for F-variables but hardkoded)
    REAL    :: minimum,maximum
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    REAL    :: lscale,gridsize,xcoord(:),ycoord(:)
    INTEGER :: corrtype,ensgen
    REAL    :: sigma, semimeta, restmeta, minsigma
    INTEGER, INTENT(IN) :: useBinFile   !<flag for using bin-file
    !LOCAL VARIABLES
    INTEGER :: nDim,dimID(1),modID
  
    nDim=0;dimID(1)=0
    
    !ALLOCATE and initialize the ensemble matrix data (NOT initialize bin-file for now)
    !CALL allocate_assim_ensemble(assimVar%x,nens,nvar,fileID,locID,coordID,xini,minimum,maximum,.true.,-9999.,.false.) !CP161205 for binfiles and assimilate turned on to apply enkf
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,varID,fileID,useBinFile,locID,coordID,xini,minimum,maximum,.true.,-9999.,.true.,.false.)

    !ALLOCATE the interface data
    nDim = 0 ; dimID(1) = 0 ; modID = -9999
    CALL allocate_assim_interface(assimVar%info,varName,varID,modID,nDim,dimID)

    !ALLOCATE and initialize the ensemble generation data
    CALL allocate_assim_generation(assimVar%gen,nvar,ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !, ECMWF)
  END SUBROUTINE allocate_assim_forcing_ensemble

!>Allocate and initialize observation ensemble variables
  SUBROUTINE allocate_assim_observation_ensemble(assimVar,nens,nvar,varName,xini,obsID,modID,coordID,fileID,minimum,maximum, &
    ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord)
    !INPUT
    TYPE(assim_input_ensemble_type) :: assimVar
    INTEGER :: nens,nvar
    REAL    :: xini(nvar)
    CHARACTER(LEN=*) :: varName
    INTEGER :: obsid, modID
    INTEGER coordID,fileID
    REAL    :: minimum,maximum
    REAL    :: lscale,gridsize
    INTEGER :: corrtype,ensgen
    REAL    :: sigma, semimeta, restmeta, minsigma
    REAL    :: xcoord(:),ycoord(:)
    !local
    INTEGER :: nDim,dimID(1)

    !ALLOCATE and initialize the ensemble matrix data (NOT initialize bin-file for now)
    !CALL allocate_assim_ensemble(assimVar%x,nens,nvar,fileID,0,coordID,xini,minimum,maximum,.true.,-9999.,.false.) !CP161205 for binfiles and assimilate turned on to apply enkf
    CALL allocate_assim_ensemble(assimVar%x,nens,nvar,obsid,fileID,0,0,coordID,xini,minimum,maximum,.true.,-9999.,.true.,.false.)
    !ALLOCATE the interface data
    nDim = 0 ; dimID(1) = 0
    CALL allocate_assim_interface(assimVar%info,varName,obsID,modID,nDim,dimID)
    !ALLOCATE and initialize the ensemble generation data
    CALL allocate_assim_generation(assimVar%gen,nvar,ensgen,sigma,semimeta,restmeta,minsigma,lscale,gridsize,corrtype,xcoord,ycoord) !,.false.)

  END SUBROUTINE allocate_assim_observation_ensemble

!>Truncates an ensemble matrix to minimum and maximum allowed values.
!----------------------------------------------------------------------------------------------
SUBROUTINE assim_checkminmax(nx,ne,ensemble,minval,maxval)
  INTEGER, INTENT(IN) :: nx
  INTEGER, INTENT(IN) :: ne
  REAL, INTENT(INOUT) :: ensemble(nx,ne)
  REAL, INTENT(IN)    :: minval, maxval
  INTEGER :: i,j
  DO i=1,ne
    DO j=1,nx
      !make sure missing values stay missing
      IF(ensemble(j,i).GT.-9998.)THEN
        ensemble(j,i)=AMAX1(minval,AMIN1(maxval,ensemble(j,i)))
      ELSE
        ensemble(j,i)=-9999.
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE assim_checkminmax
  
!>General routine for Ensemble generation (forcing and observation data).
!>
!> The ensemble generation is made by adding random numbers to the input data.
!>
!> The basic assumption is that the random perturbations are gaussian with zero mean and standard deviation sigma.
!>
!> However, based on Turner et al(2008), it is assumed that in order to get unbiased input ensembles, 
!> we need to consider two types of perturbations - systematic and random:
!> Thus, the input x at time k for ensemble member j, x_kj = x_k0+eata_kj+chi_j,
!> where eata_kj is regenerated every time step and chi_j is generated only once.
!> Then, eata and chi can be generated with three different types of restrictions on variables.
!> For the moment, we assume that the static error chi_j = 0.
!>
!> In adition, we now also take into account spatial correlation in the data, by generating 
!> spatially correlated random data using a FFT-based method.
!----------------------------------------------------------------------------------------------
SUBROUTINE generate_input_ensemble(n,assimVar)
  INTEGER n
  TYPE(assim_input_ensemble_type) :: assimVar(n)
  INTEGER i,j
  REAL midvalue
  REAL,ALLOCATABLE :: localx(:,:)
    
  !loop over ensemble vector members
  DO j=1,n

    !get data to a local array for manipulation
    ALLOCATE(localx(assimVar(j)%x%nvar,assimVar(j)%x%nens))

    !before ensemble generation, check that number of ensemble members is greater than 1   
    IF(assimVar(j)%x%nens.gt.1)THEN 

      !a) determine standard deviation (sigma) in each model unit (rows in ensemble matrix)
        
      !loop over number of model units (subbasins in HYPE)
      DO i=1,assimVar(j)%x%nvar
!        ! select General case or special "ECMWF" case
!        IF(.not. ensemble(j)%gen%ECMWF)THEN
        
        !Check missing values
        IF(assimVar(j)%gen%mean(i).GT.-9998.)THEN
          SELECT CASE(assimVar(j)%gen%ensgen)
            CASE(4) ! restricted (Turner et al, 2008)
              midvalue = 0.5 * (assimVar(j)%x%maximum+assimVar(j)%x%minimum)
              IF(assimVar(j)%gen%mean(i).gt.midvalue)THEN
                assimVar(j)%gen%sigma(i) = assimVar(j)%gen%restmeta *(assimVar(j)%x%maximum-assimVar(j)%gen%mean(i)) / (assimVar(j)%x%maximum-midvalue);
              ELSE
                assimVar(j)%gen%sigma(i) = assimVar(j)%gen%restmeta*(assimVar(j)%gen%mean(i)-assimVar(j)%x%minimum) / (midvalue - assimVar(j)%x%minimum);
              ENDIF
            CASE(3) ! semirestricted with max (Turner et al, 2008)
              assimVar(j)%gen%sigma(i) = AMAX1(0.0,assimVar(j)%x%maximum-assimVar(j)%gen%mean(i)) * assimVar(j)%gen%semimeta
            CASE(2) ! semirestricted with min (Turner et al, 2008)
              assimVar(j)%gen%sigma(i) = AMAX1(0.0,assimVar(j)%gen%mean(i)-assimVar(j)%x%minimum) * assimVar(j)%gen%semimeta
            CASE DEFAULT ! 0 or 1, ie. unrestricted (Turner et al, 2008)
          END SELECT
!         ENDIF
        ! just make sure we dont have sigma = 0
          assimVar(j)%gen%sigma(i) = AMAX1(assim_minsigma,assimVar(j)%gen%sigma(i))
        ENDIF
      ENDDO
          
      ! b) generate random values with the assigned sigma
      IF(assimVar(j)%gen%dorandxy)THEN
        !Spatially correlated ensemble matrix for a specific variable type
        CALL get_spatially_correlated_random_data2(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%gen,localx)
      ELSE
        ! Non-correlated observation data
        ! loop over records in the observation matrix (irrespective of variable type)
        DO i=1,assimVar(j)%x%nvar
          ! generate ensemble for each record independently, no correlation
          IF(assimVar(j)%gen%mean(i).GT.-9998.)THEN
            !CALL get_random_vector(assimVar(j)%x%nens,assimVar(j)%gen%mean(i),assimVar(j)%gen%sigma(i),assimVar(j)%x%x(i,:)) !CP161205 have to use local array since %x not always allocated (binfiles)
            CALL get_random_vector_gaussian(assimVar(j)%x%nens,assimVar(j)%gen%mean(i),assimVar(j)%gen%sigma(i),localx(i,:)) ! DG20170427 Changed to replace use of MKL dependent code.
          ELSE
            !assimVar(j)%x%x(i,:)=-9999.  !CP161205 same as above
            localx(i,:)=-9999.
          ENDIF
        ENDDO
      ENDIF
    ELSE
      ! in case someone runs with 1 ensemble member, we just use the mean, and skip all the ensemble generation
      !DO i=1,assimVar(j)%x%nvar
      !  assimVar(j)%x%x(i,1)=assimVar(j)%gen%mean(i)   !CP161205 use local x since x not allocated if binfile is used
      !ENDDO    
      DO i=1,assimVar(j)%x%nvar
        localx(i,1)=assimVar(j)%gen%mean(i)
      ENDDO    
    ENDIF
    
    ! save manipulated data from local array (to matrix or bin-file), also truncate all values outside the min/max range to the limits 
    !CALL assim_set_ensemble_data(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%x,localx,Rec,.true.,.false.) !checkminmax true here
    CALL assim_set_ensemble_data(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%x,localx,.true.,.false.) !checkminmax true here
    DEALLOCATE(localx)
    ! Finally, truncate all values outside the min/max range to the limits.
    !CALL assim_checkminmax(assimVar(j)%x%nvar,assimVar(j)%x%nens,assimVar(j)%x%x,assimVar(j)%x%minimum,assimVar(j)%x%maximum)  !CP161205 checked in assim_set_ensemble_data
  ENDDO
  
END SUBROUTINE generate_input_ensemble

!> Function for spatially correlated random variable generation
SUBROUTINE get_spatially_correlated_random_data2(n,nens,assimG,X)
  
  !Argument declarations
  INTEGER, INTENT(IN) :: n      !<number of element (nsub usually)
  INTEGER, INTENT(IN) :: nens   !<number of ensembles
  TYPE(assim_generation_type),INTENT(INOUT) :: assimG !<information on generation of input data
  REAL,INTENT(OUT) :: X(n,nens) !<generated random matrix
  
  !Local variables
  INTEGER :: i,iens   !loop-variables
     
  !>\b Algorithm
  !> Loop over ensemble members:
   DO iens=1,nens
     
     !>\li generate random numbers (with mean 0, standard deviation 1 and correlation length L) for all subbasins' location
     CALL resample_randxy_data(assimG%myrandxy_data)

     !>\li scale perturbations with the required sigma (standard deviation) and appy to the mean
     DO i=1,n
       X(i,iens) = assimG%myrandxy_data%pert_xy(i) * assimG%sigma(i)  + assimG%mean(i)
     ENDDO
   ENDDO
   
END SUBROUTINE get_spatially_correlated_random_data2

!>\brief Get an array of random gaussian values with specified mean and sigma
! USE THE RGAUSS() FROM THE RANDOM_ROUTINES MODULE INSTEAD OF MKL library, DG 2017-04-27
SUBROUTINE get_random_vector_gaussian(n,a,sigma,r)
  REAL(KIND=4) r(:)
  REAL(KIND=4) a,sigma
  INTEGER n, i
  !loop over vector and generate random numbers
  DO i=1,n
    ! random number N[0,1], scaled to the requested sigma, and applied to the mean
    r(i)=rgauss()*sigma+a
  ENDDO
END SUBROUTINE get_random_vector_gaussian

! -----------------------------------------------------------------------------------
! MATRIX OPERATIONS
! Method MATMUL (intrinsic function) /David 2017-04-27
! -----------------------------------------------------------------------------------
!>\brief Subroutine for matrix multiplication.
SUBROUTINE matrixmatrixmultiply(mat1,mat2,matout)
  REAL(KIND=4),INTENT(IN) ::  mat1(:,:),mat2(:,:)
  REAL(KIND=4),INTENT(OUT) :: matout(:,:)
  ! make sure matout is 0
  matout(:,:)=0.0
  ! CALL MATMUL
  matout = matmul(mat1,mat2)
END SUBROUTINE matrixmatrixmultiply

! -----------------------------------------------------------------------------------
! CHOLESKY SOLUTION
! Methods from Numerical Recipes in Fortran 90, 2nd edition, Press et al, 1996. 
! Adopted by David 2017-04-27
! -----------------------------------------------------------------------------------

!> Cholesky decomposition (from Numerical Recipes, adopted by David)
SUBROUTINE choldc(a,p,n)
  !USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
  IMPLICIT NONE
  REAL, DIMENSION(:,:), INTENT(INOUT) :: a
  REAL, DIMENSION(:), INTENT(OUT) :: p
  !Given an N × N positive-definite symmetric matrix a, this routine constructs its Cholesky
  !decomposition, A = L · LT . On input, only the upper triangle of a need be given; it is
  !not modified. The Cholesky factor L is returned in the lower triangle of a, except for its
  !diagonal elements, which are returned in p, a vector of length N.
  INTEGER :: i,n
  REAL :: summ
  !n=assert_eq(size(a,1),size(a,2),size(p),’choldc’)
  do i=1,n
    summ=a(i,i)-dot_product(a(i,1:i-1),a(i,1:i-1))
    if (summ <= 0.0) write(6,*)'choldc failed'    !CP170619 should we not stop and handle this error somehow?
    p(i)=sqrt(summ)
    a(i+1:n,i)=(a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)))/p(i)
end do
END SUBROUTINE choldc

!> Cholesky solver (from Numerical Recipes, adopted by David)
SUBROUTINE cholsl(a,p,b,x,n)
  !USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL, DIMENSION(:,:), INTENT(IN) :: a
  REAL, DIMENSION(:), INTENT(IN) :: p,b
  REAL, DIMENSION(:), INTENT(INOUT) :: x
  !Solves the set of N linear equations A · x = b, where a is a positive-definite symmetric
  !matrix. a (N × N) and p (of length N) are input as the output of the routine choldc.
  !Only the lower triangle of a is accessed. b is the input right-hand-side vector, of length N.
  !The solution vector, also of length N, is returned in x. a and p are not modified and can be
  !left in place for successive calls with different right-hand sides b. b is not modified unless
  !you identify b and x in the calling sequence, which is allowed.
  INTEGER :: i,n
  !n=assert_eq((/size(a,1),size(a,2),size(p),size(b),size(x)/),’cholsl’)
  do i=1,n !Solve L · y = b, storing y in x.
    x(i)=(b(i)-dot_product(a(i,1:i-1),x(1:i-1)))/p(i)
  enddo
  do i=n,1,-1 !Solve LT · x = y.
    x(i)=(x(i)-dot_product(a(i+1:n,i),x(i+1:n)))/p(i)
  enddo
END SUBROUTINE cholsl

!> Cholesky solution
!!
! new routine called from enkf_analysis_prepare, written by David 2017-04-27
SUBROUTINE choleskysolution(ny,ne,P,Y,M)
  !arguments
  INTEGER, INTENT(IN)  :: ny,ne
  REAL,    INTENT(INOUT)  :: P(ny,ny)
  REAL,    INTENT(IN)     :: Y(ny,ne)
  REAL,    INTENT(OUT) :: M(ny,ne)
  !local variables
  REAL                 :: Ldiag(ny)
  integer              :: i 
  !This routine uses Cholesky factorization to solve the set of linear eqations P*M=Y => M=inv(P)*Y
  !The routines for cholesky factorization and solving above are taken from Numerical Recipes for Fortran 90.
  !The result M, is returned to enkf_analysis_prepare.
  !By David, 2017-04-27
  
  ! 1) Cholesky factorization, P = L · LT
  CALL choldc(P,Ldiag,ny)     !Cholesky factors are now in lower triangle of P, except the diagonal elements which are stored in Ldiag
  
  ! 2) Cholesky solution, solving P · M = Y => M = inv(P) · Y, using P and Ldiag from previous steps
  DO i=1,ne
    ! Since Y and M are multi-column ensemble matrices, we apply the solution for each ensemble member in a loop
    CALL  cholsl(P,Ldiag,Y(:,i),M(:,i),ny)
  ENDDO
  
END SUBROUTINE choleskysolution


! ---------------------------------------------------------------------------------------------
! enkf_analysis_prepare
! enkf_analysis_apply
! enkf_analysis_main
!
! Ensemble Kalman filter equations (Evensen) following (a) Mandel, J. Efficient implementation 
! of the ensemble kalman filter, and (b) DeChant, C. Examining the effectiveness and 
! robustness of sequential data assimilation methods for quantification of uncertainty
! in hydrological forecasting. Localization following Magnusson, Gustafsson et al (2014)
!
! Part I: Innovations (Y=D-HX), and inversion of variances M=1/(var(HX) + R) which is later
!         used in the update equation: X_analysis = X_forecast + cov(X,HX)/(var(HX)+R) * (D-HX)
!
!Author: D.Gustafsson, SMHI/KTH.
!-----------------------------------------------------------------------------------------------
!The implementation is divided into two steps according to this pseudo-code (enkf_analysis_main)
!  CALL enkf_analysis_main()
!      
!   1. CALL enkf_analysis_prepare()
!        a. calculation of innovation, Y = D-HX
!        b. calculation and localization of the covariance of predicted observations Cyy = cov(HX) = locCyy * cov(HX)
!        c. combine cov(HX) with observation error variance R into P = (cov(HX)+R)
!        d. derive intermediate update matrix M, by inversion of variance matrice sum  M = P^(-1)Y
!
!   2. DO i=1,num_model_states
!        CALL enkf_analysis_apply()
!             Here in Part 2, M is multiplied with the covariance Cxy = Cov(X,HX) in the final EnKF update equation:
!      ENDDO
! 
! Terminology, EnKF basics from Mandel and DeChant
! ---------------------------------------------------------------------------------------
! Variable   Mandel        Dechant     Description
! --------   ------------  ----------  --------------------------------------------------
! N          N              Nens        ensemble size, n:o of ensemble members
! NX         n              -           state vector length, n:o state variables
! ND         m              
! X          X              X           model state variables
! A          A              e           ensemble deviation (from ensemble mean, A=X-E(X))
! D          D              y+eps       observations (including observation error eps)
! HX         HX             y'          predicted observations
! -          h(x)           h(x)        observation operator y'=h(x) translating model states to observation space
! HA         HA             -           HX ensemble deviation from ensemble mean (HA=HX-E(HX))
! Y          Y=D-HX         y+eps-y'    innovation, deviation between observations and predicted obserations
! -          C              P           model state error covariance, C = AA^T (Mandel) = ee^T = P (DeChant)
! R          R              R           observation error variance (theoretic or sample, assumed to be uncorrelated)
! -          HCH^T          C_yy        variance of predicted observations
! -          A(HA)^T        C_xy        covariance between model states and predicted observations
! PP         P              C_yy+R      sum of predicted observation variance and observation error variance
! M          M=(P^-1)Y      -           intermediate result of the inversion M = (P^-1)Y, with P and Y from Mandel
! locCyy     -              -           localization matrix for the variance of HX
! locCxy     -              -           localization matrix for the covariance between X and HX 
!-----------------------------------------------------------------------------------
!>\brief Part 1 of ensemble Kalman filter analysis
!!
!> Innovations (Y=D-HX), and inversion of variances M=1/(var(HX) + R) which is later
!> used in part 2; the update equation: 
!>
!> X_analysis = X_forecast + cov(X,HX)/(var(HX)+R) * (D-HX)
!>
SUBROUTINE enkf_analysis_prepare(N,ND,D,HX,R,locCyy,M,Y,HA)
  !INPUT ARGUMENTS
  INTEGER, INTENT(IN)  :: N            !<number of ensemble members
  INTEGER, INTENT(IN)  :: ND           !<number of observations
  REAL,    INTENT(IN)  :: D(ND,N)      !<observation ensemble
  REAL,    INTENT(IN)  :: HX(ND,N)     !<predicted observation ensemble
  REAL,    INTENT(IN)  :: R(ND)        !<observation error variance, uncorrelated
  REAL,    INTENT(IN)  :: locCyy(ND,ND)!<localization matrix for cov(HX)
  REAL,    INTENT(OUT) :: M(ND,N)      !<inverse of (R+CXY) multiplied with innovation Y
  REAL,    INTENT(OUT) :: Y(ND,N)      !<innovation ensemble
  REAL,    INTENT(OUT) :: HA(ND,N)     !<HX deviation from mean(HX)
  !LOCAL VARIABLES
  REAL, ALLOCATABLE    :: e_Nx1(:,:)   !unit vector, size Nx1
  REAL, ALLOCATABLE    :: e_1xN(:,:)   !unit vector, size 1xN
  REAL, ALLOCATABLE    :: z_NDx1(:,:)  !intermediate result vector, size NDx1
  REAL, ALLOCATABLE    :: PP(:,:)      !intermediate result vector, size NDxND
  INTEGER              :: I,J          !loop index
  
  !ALLOCATION and INITIALIZATION
  ALLOCATE(e_Nx1(N,1))
  ALLOCATE(e_1xN(1,N))
  ALLOCATE(z_NDx1(ND,1))
  ALLOCATE(PP(ND,ND))
  
  !1) assign unit vectors
  e_Nx1 = 1.
  e_1xN = 1.
    
  !2) HA = HX-E(HX) = HX-1/N * (HX*e_Nx1)*e_1xN,  in three steps:
  ! a. z_NDx1 = HX*e_Nx1
  CALL matrixmatrixmultiply(HX,e_Nx1,z_NDx1)
  ! b. HA = z_NDx1 * e_1xN
  CALL matrixmatrixmultiply(z_NDx1, e_1xN, HA) !HA = z * e_1xN
  ! c. !HA = HX - 1/N * HA
  HA = HX - HA/N                               
    
  !3) Y = D-HX, in one step
  Y = D - HX

    
  !4) PP = (R+cov(HX)), observation error variance + predicted observation covariance, in three step:              
  ! a. PP = 1/(N-1) * HA * (HA)' = cov(predicted_observations)
  CALL matrixmatrixmultiply(HA, TRANSPOSE(HA), PP) ; PP = PP/(N-1)
    
  ! b. Localization by elementwise multiplication, PP = locCyy[nobs,nobs] .* PP
  DO I=1,ND
    DO J=1,ND
      PP(J,I)=PP(J,I)*locCyy(J,I)
    ENDDO
  ENDDO
    
  ! c. PP = PP + R, adding obs error variance, using the theoretical (assumed uncorrelated) variance R 
  !                 rather than the sample covariance
  DO I = 1,ND
    PP(I,I) = PP(I,I) + R(I)
  ENDDO

  ! 5) Cholesky solution... LL' = P, M=inv(P) * Y    ie. M = "innovations" / (cov(obs) + cov(model))
!  CALL choleskysolution(PP, Y, M) ! call to old MKL based function
  CALL choleskysolution(ND, N, PP, Y, M)
  
  !Deallocate temporary results
  DEALLOCATE(PP)
  DEALLOCATE(e_Nx1)
  DEALLOCATE(e_1xN)
  DEALLOCATE(z_NDx1)

END SUBROUTINE enkf_analysis_prepare

! -----------------------------------------------------------------------------------
!>\brief Part 2 of ensemble Kalman filter analysis
!>
!>Apply ENKF filter on a model variable ensemble using the M
!>matrix derived in the enkf_analysis_prepare function.
!>
!> EnKF analysis, following the "basic implementation" in Jan Mandel report modified
!> to include covariance localization: 
!>
!>                (K = Cxy .* loc_Cxy / (Cyy .* loc_Cyy + R))
!>
!> The localization of Cyy is embedded in the M matrix which is 
!> prepared by the enkf_analysis_prepare routine.
!>
!> The localization of Cxy is introduced in the final step of this
!> "apply"-routine, which involves a change in the calculation of the AZ
!> matrix compared to Mandel:
!>
!>                Mandel:    AZ = A * (HA' * M)
!>
!>                Here:      AZ = ((A*HA') .* loc_Cxy) * M
!>
!> Matrix operations are made in an order defined by the parantheses (most inner operations first).
!> This might have consequences on the number of operations compared to Mandel, but I have not checked that.
! -----------------------------------------------------------------------------------
SUBROUTINE enkf_analysis_apply(N,NX,ND,Y,M,HA,locCxy,X)
  !INPUT VARIABLES
  INTEGER, INTENT(IN)  :: N            !<number of ensemble members
  INTEGER, INTENT(IN)  :: NX           !<number of variables in ensemble (rows)
  INTEGER, INTENT(IN)  :: ND           !<number of observations
  REAL,INTENT(INOUT)   :: X(NX,N)      !<ensemble data
!  REAL,INTENT(IN)       :: M(NX,N)
  REAL,INTENT(IN)      :: M(ND,N)      !<inverse of (R+CXY) multiplied with innovation Y
  REAL,INTENT(IN)      :: Y(ND,N)      !<innovation ensemble
  REAL,INTENT(IN)      :: HA(ND,N)     !<HX deviation from mean(HX)
  REAL,INTENT(IN)      :: locCxy(NX,ND)!<localization matrix
        
  !LOCAL VARIABLES
  REAL, ALLOCATABLE :: e_Nx1(:,:)
  REAL, ALLOCATABLE :: e_1xN(:,:)
  REAL, ALLOCATABLE :: z_NXx1(:,:)
  REAL, ALLOCATABLE :: Z_NXxND(:,:)
  REAL, ALLOCATABLE :: A(:,:)
  REAL, ALLOCATABLE :: AZ(:,:)
  INTEGER           :: I, J
  
  !ALLOCATION and INITIALIZATION
  ALLOCATE(A(NX,N))
  ALLOCATE(AZ(NX,N))
  ALLOCATE(e_Nx1(N,1))
  ALLOCATE(e_1xN(1,N))
  ALLOCATE(z_NXx1(NX,1))
  ALLOCATE(z_NXxND(NX,ND))
    
  !0) assign unit vectors
  e_Nx1 = 1.
  e_1xN = 1.

  !1) A = X-E(X) = X-1/N*(X*e_Nx1)*e_1xN
  CALL matrixmatrixmultiply(X, e_Nx1, z_NXx1)
  CALL matrixmatrixmultiply(z_NXx1, e_1xN, A)
  A = X - A / N

  !2) EnKF filtering equation
  !
  ! Modification of Mandel's two-step algoritm, to include localization of 
  ! covariances between innovations (Y=HX-D) and all model states (X) 
  !
  ! The Kalman filter analysis is written as:
  !     X = X + locCxy .* (A*HA') * M / (N-1) ;
  !  
  ! Which is formulate the two-step solution from Mandel, into a three-step solution:
  !  7.1: Z_NXxND = A*HA'
  CALL matrixmatrixmultiply(A, TRANSPOSE(HA), Z_NXxND)
  !    
  !  7.2: Z_NXxND = Z_NXxND .* locCX (element-wise)
  DO I=1,ND
    DO J=1,NX
      Z_NXxND(J,I)=Z_NXxND(J,I) * locCxy(J,I)
    ENDDO
  ENDDO
  !    
  !  7.3: Z_NXxN  = Z_NXxND * M      (we can use AZ for this)
  CALL matrixmatrixmultiply(Z_NXxND, M, AZ)
  !
  !! 7.1 Z = (HA)' * M
  !CALL matrixmatrixmultiply(TRANSPOSE(HA), M, Z_NxN)
  ! 
  !! 7.2) X=X+1/(N-1) * A * Z
  !CALL matrixmatrixmultiply(A, Z_NxN, AZ)
  !    
  ! 8) final update, just as before:
  X = X + AZ / (N-1)

  !Deallocation of local arrays necessary for gfortran
  IF(ALLOCATED(A)) DEALLOCATE(A)
  IF(ALLOCATED(AZ)) DEALLOCATE(AZ)
  IF(ALLOCATED(e_Nx1)) DEALLOCATE(e_Nx1)
  IF(ALLOCATED(e_1xN)) DEALLOCATE(e_1xN)
  IF(ALLOCATED(z_NXx1)) DEALLOCATE(z_NXx1)
  IF(ALLOCATED(z_NXxND)) DEALLOCATE(z_NXxND)
  
END SUBROUTINE enkf_analysis_apply

!------------------------------------------------------
!>Calculates summary statistics for an ensemble matrix
!------------------------------------------------------
SUBROUTINE assim_ensemble_statistics(xin, DIM, NN, xmean, xmins, xmaxs, xsigma, xmedian, domedian) !xquantiles,xcovar, xcor)
  USE compout, only: calculate_median

  !Argument declarations
  INTEGER DIM, NN
  REAL xin(DIM,NN)
  REAL xmean(DIM), xmins(DIM), xmaxs(DIM), xsigma(DIM), xmedian(DIM) !xquantiles(3,DIM), xcovar(DIM,DIM), xcor(DIM,DIM)
!  REAL o_stat(NN,DIM),xr2(DIM)
!  REAL o_quant(3)
  LOGICAL domedian
  
! Data needed for simplified statistical calculation
  REAL s(NN), pp(NN)
  REAL var, ep
  INTEGER i
  
!LOTTA TO DAVID FOR CHECK VSL SUMMARY STATISTICS BEGIN
!Need to keep?
! Data needed to use the VSL Summary statistics routine
!  TYPE(VSL_SS_TASK) task
!  INTEGER p
!  INTEGER n
!  INTEGER x_storage
!  INTEGER o_storage
!  !      INTEGER cov_storage
!  !      INTEGER cor_storage
!  INTEGER i, errcode
!  INTEGER(KIND=8) estimate
!  INTEGER task_method
!LOTTA TO DAVID FOR CHECK VSL SUMMARY STATISTICS END

  
  ! loop over dimension DIM
  DO i=1,DIM
    ! mean value
    xmean(i)=sum(xin(i,:))/NN
    
    ! Min and Max
    xmins(i)=minval(xin(i,:))
    xmaxs(i)=maxval(xin(i,:))
    
    ! standard deviation
    s(:)=xin(i,:)-xmean(i)
    ep=sum(s(:))
    pp(:)=s(:)*s(:)
    var = (sum(pp(:))-ep**2/NN)/(NN-1)
    xsigma(i)=sqrt(var)
    
    ! median
    IF(domedian) CALL calculate_median(NN,xin(i,:),-9999.,xmedian(i))
      
    ! skip the quantiles
  ENDDO
    
!LOTTA TO DAVID FOR CHECK VSL SUMMARY STATISTICS BEGIN
!Need to keep this?
!! Seems like this VSL code is quite slow, or at least when using it many times, it uses a lot of resources.
!! we try something simpler above
!  
!  ! ***** Initializing parameters for Summary Statistics task *****
!  p               = DIM
!  n               = NN
!  x_storage       = VSL_SS_MATRIX_STORAGE_COLS
!  !   cov_storage     = VSL_SS_MATRIX_STORAGE_FULL
!  !     cor_storage     = VSL_SS_MATRIX_STORAGE_FULL
!  task_method     = VSL_SS_METHOD_FAST
!  o_storage   = VSL_SS_MATRIX_STORAGE_ROWS
!      
!  o_quant(1) = 0.025 ; o_quant(2) = 0.5 ; o_quant(3) = 0.975
!      
!  M=3
!
!  !     ***** Create Summary Statistics task *****
!  errcode = vslsssnewtask( task, p, n, x_storage, xin )
!
!  !     ***** Set initial values of the min/max estimates *****
!  DO i = 1, p
!    xmins(i) = xin(i, 1)
!    xmaxs(i) = xin(i, 1)
!  ENDDO
!
!  !     ***** Edit task parameters for min and max computation *****
!  errcode = vslsssedittask( task, VSL_SS_ED_MIN, xmins )
!  errcode = vslsssedittask( task, VSL_SS_ED_MAX, xmaxs )
!
!  !      ***** Minimum and maximum are included in the list of estimates to compute *****
!  estimate = IOR( VSL_SS_MIN, VSL_SS_MAX )
!
!  !     ***** Compute the estimates using FAST method *****
!  errcode = vslssscompute( task, estimate, task_method )
!
!  !     ***** Delete Summary Statistics task *****
!  errcode = vslssdeletetask( task )
!
!  !     ***** Create Summary Statistics task *****
!  errcode = vslsssnewtask( task, p, n, x_storage, xin )
!
!  !     ***** Edit task parameters for computating of mean estimate and 2nd, 3rd
!  !           and 4th raw and central moments estimates *****
!  errcode = vslssseditmoments( task, xmean, xr2, c2m = xsigma)
!
!  !     ***** Mean and 2nd, 3rd and 4th raw and central moments are included
!  !            in the list of estimates to compute *****
!  estimate = IOR( VSL_SS_MEAN,   IOR(VSL_SS_2R_MOM,VSL_SS_2C_MOM))
!
!  !      estimate = IOR( estimate, IOR( VSL_SS_MEAN,                              &
!  !     &           IOR( VSL_SS_2R_MOM, IOR( VSL_SS_3R_MOM,                       &
!  !     &           IOR( VSL_SS_4R_MOM, IOR( VSL_SS_2C_MOM,                       &
!  !     &           IOR( VSL_SS_3C_MOM, VSL_SS_4C_MOM ) ) ) ) ) ) )
!
!  !     ***** Compute the estimates using FAST method *****
!  errcode = vslssscompute( task, estimate, task_method )
!
!  !     ***** Delete Summary Statistics task *****
!  errcode = vslssdeletetask( task )
!
!!!     ***** Create Summary Statistics task *****
!!      errcode = vslsssnewtask( task, p, n, x_storage, xin )
!!
!!!     ***** Initialization of the task parameters using FULL_STORAGE
!!!           for covariance/correlation matrices computation *****
!!      errcode = vslssseditcovcor( task, xmean, xcovar, cov_storage, xcor, cor_storage )
!!
!!!     ***** Mean, Covariance and correlation matrices are included in the list
!!!           of estimates to compute *****
!!      estimate = IOR( VSL_SS_COV, VSL_SS_COR )
!!
!!!     ***** Compute the estimates using FAST method *****
!!      errcode = vslssscompute( task, estimate, task_method )
!!     
!!!     ***** Delete Summary Statistics task *****
!!      errcode = vslssdeletetask( task )
!
!
!  !     ***** Create Summary Statistics task *****
!  errcode = vslsssnewtask( task, p, n, x_storage, xin )
!  errcode = vslssseditquantiles( task, M, o_quant, xquantiles, o_stat, o_storage )
!  estimate = IOR( VSL_SS_QUANTS, VSL_SS_ORDER_STATS )
!    
!  !     ***** Compute the estimates using FAST method *****
!  errcode = vslssscompute( task, estimate, task_method )
!     
!  !     ***** Delete Summary Statistics task *****
!  errcode = vslssdeletetask( task )
!LOTTA TO DAVID FOR CHECK VSL SUMMARY STATISTICS END
  
END SUBROUTINE assim_ensemble_statistics

!---------------------------------------------------------------------
!> Routine that returns the ensemble data in a matrix. 
!> If needed, the data is read from binary file.
!> For speed, records in binary files contain all states in one
!> ensemble member (X transformated).
!---------------------------------------------------------------------
SUBROUTINE assim_get_ensemble_data(NX,NE,ensData,x)
  !input
  INTEGER, INTENT(IN)                    :: NX, NE    !number of states and ensemble members
  REAL, INTENT(OUT)                      :: X(NX,NE)  !ensemble data matrix
  TYPE(assim_ensemble_type),INTENT(IN)   :: ensData   !ensemble data structure
  
  !Local variables
  INTEGER j
  
  !ensemble data ALLOCATED in memory
  IF(ALLOCATED(ensData%x))THEN   !change because needed more options for different kind of bin-files
    !X=ensData%x(:,Rec:(Rec+NE-1))  !CP161202 change Rec to 0 (from 1 in general)
    X=ensData%x(:,ensData%rec+1:ensData%rec+NE)
  ELSE
  !ensemble data read from direct access binary file
    DO j = 1,NE
      READ(ensdata%fileID,REC=ensData%rec+j) X(:,j)
    ENDDO
  ENDIF
END SUBROUTINE assim_get_ensemble_data

!---------------------------------------------------------------------
!> Routine that writes a matrix into an ensemble data structure. 
!> If needed, the data is written to a binary file.
!---------------------------------------------------------------------
SUBROUTINE assim_set_ensemble_data(NX,NE,ensData,x,checkMinMax,doStat)
  !input
  INTEGER, INTENT(IN)                     :: NX, NE
  REAL, INTENT(INOUT)                     :: X(NX,NE)
  TYPE(assim_ensemble_type),INTENT(INOUT) :: ensData
  LOGICAL, INTENT(IN)                     :: checkMinMax !check min and max limits before update
  LOGICAL, INTENT(IN)                     :: doStat      !calculate output statistics

  !Local variables
  INTEGER j
  
  !Check min and max limits
  IF(checkMinMax)THEN
    CALL assim_checkminmax(nx,ne,x,ensdata%minimum,ensData%maximum)
  ENDIF

  !Write x to ensemble data...
  IF(ALLOCATED(ensData%x))THEN
    !... ALLOCATED in memory I don't under stand how Rec works here, removed.
    !ensData%x(:,Rec:(Rec+NE-1))=X  !CP161202 change Rec to 0 (from 1 in general), Rec is always 0 here
    ensData%x(:,ensData%rec+1:ensData%rec+NE)=X
  ELSE
    !... or in direct access binary file.
    DO j = 1,NE
      WRITE(ensdata%fileID,REC=ensData%rec+j) X(:,j)
    ENDDO
  ENDIF
  
  !Update ensemble output statistics
  !IF(doStat)THEN
  !  CALL assim_ensemble_statistics(x, NX, NE, ensData%outmean, ensData%outquant, ensData%outmin, ensData%outmax, ensData%outsigma)
  !ENDIF
END SUBROUTINE assim_set_ensemble_data

!>Calculates ensemble statistics for state, forcing and auxiliary ensembles
SUBROUTINE updateEnsembleStatistics(assimData,total_time)
  TYPE(assim_data_type) :: assimData !<main assimilation variable containing all data
  REAL, OPTIONAL, INTENT(INOUT) :: total_time(4) !<optional timing of the subroutine
  REAL start_time, stop_time
  INTEGER i
  REAL,ALLOCATABLE :: X(:,:)  !intermediate matrix of ensemble
  LOGICAL domedian

  IF(PRESENT(total_time)) CALL cpu_time(start_time)

  ! set domedian if needed for outputs
  IF(assimData%info%meanout)THEN
    domedian=.FALSE. ! median NOT needed for outputs
  ELSE
    domedian=.TRUE.  ! median IS needed for outputs
  ENDIF
  
  !state ensembles
  DO i=1,assimData%info%nX
    !Allocate and collect X from where it is saved before calculating the statistics CP added to handle bin-files.
    ALLOCATE(X(assimData%X(i)%x%nvar,assimData%X(i)%x%nens))

    IF(PRESENT(total_time))THEN
      call cpu_time(stop_time)
      total_time(1)=total_time(1)+stop_time-start_time
      start_time=stop_time
    ENDIF

    CALL assim_get_ensemble_data(assimData%X(i)%x%nvar,assimData%X(i)%x%nens,assimData%X(i)%x,X)   !CP161201 added call to routine to get X-matrix (useful if bin-file)

    IF(PRESENT(total_time))THEN
      call cpu_time(stop_time)
      total_time(2)=total_time(2)+stop_time-start_time
      start_time=stop_time
    ENDIF

!    CALL assim_ensemble_statistics(X, assimData%X(i)%x%nvar, assimData%X(i)%x%nens, assimData%X(i)%x%outmean, assimData%X(i)%x%outquant, assimData%X(i)%x%outmin, assimData%X(i)%x%outmax, assimData%X(i)%x%outsigma)
    CALL assim_ensemble_statistics(X, assimData%X(i)%x%nvar, assimData%X(i)%x%nens, assimData%X(i)%x%outmean, assimData%X(i)%x%outmin, assimData%X(i)%x%outmax, assimData%X(i)%x%outsigma, assimData%X(i)%x%outquant(2,:),domedian)

    IF(PRESENT(total_time))THEN
      call cpu_time(stop_time)
      total_time(3)=total_time(3)+stop_time-start_time
      start_time=stop_time
    ENDIF

    DEALLOCATE(X)

    IF(PRESENT(total_time))THEN
      call cpu_time(stop_time)
      total_time(4)=total_time(4)+stop_time-start_time
      start_time=stop_time
    ENDIF

    !CALL assim_ensemble_statistics(assimData%X(i)%x%x, assimData%X(i)%x%nvar, assimData%X(i)%x%nens, assimData%X(i)%x%outmean, assimData%X(i)%x%outquant, assimData%X(i)%x%outmin, assimData%X(i)%x%outmax, assimData%X(i)%x%outsigma)
  ENDDO

  !forcing ensembles
  DO i=1,assimData%info%nF
    !Allocate and collect X from where it is saved before calculating the statistics, added to handle bin-files.
    ALLOCATE(X(assimData%F(i)%x%nvar,assimData%F(i)%x%nens))
    CALL assim_get_ensemble_data(assimData%F(i)%x%nvar,assimData%F(i)%x%nens,assimData%F(i)%x,X)
!    CALL assim_ensemble_statistics(X, assimData%F(i)%x%nvar, assimData%F(i)%x%nens, assimData%F(i)%x%outmean, assimData%F(i)%x%outquant, assimData%F(i)%x%outmin, assimData%F(i)%x%outmax, assimData%F(i)%x%outsigma)
    CALL assim_ensemble_statistics(X, assimData%F(i)%x%nvar, assimData%F(i)%x%nens, assimData%F(i)%x%outmean, assimData%F(i)%x%outmin, assimData%F(i)%x%outmax, assimData%F(i)%x%outsigma, assimData%F(i)%x%outquant(2,:),domedian)
    DEALLOCATE(X)
    !CALL assim_ensemble_statistics(assimData%F(i)%x%x, assimData%F(i)%x%nvar, assimData%F(i)%x%nens, assimData%F(i)%x%outmean, assimData%F(i)%x%outquant, assimData%F(i)%x%outmin, assimData%F(i)%x%outmax, assimData%F(i)%x%outsigma) !CP161205 use binfiles!
  ENDDO

  !auxiliary ensembles
  
  DO i=1,assimData%info%nA
    !Allocate and collect X from where it is saved before calculating the statistics, added to handle bin-files.
    ALLOCATE(X(assimData%A(i)%x%nvar,assimData%A(i)%x%nens))
    CALL assim_get_ensemble_data(assimData%A(i)%x%nvar,assimData%A(i)%x%nens,assimData%A(i)%x,X)
!    CALL assim_ensemble_statistics(X, assimData%A(i)%x%nvar, assimData%A(i)%x%nens, assimData%A(i)%x%outmean, assimData%A(i)%x%outquant, assimData%A(i)%x%outmin, assimData%A(i)%x%outmax, assimData%A(i)%x%outsigma)
    CALL assim_ensemble_statistics(X, assimData%A(i)%x%nvar, assimData%A(i)%x%nens, assimData%A(i)%x%outmean, assimData%A(i)%x%outmin, assimData%A(i)%x%outmax, assimData%A(i)%x%outsigma, assimData%A(i)%x%outquant(2,:),domedian)
    DEALLOCATE(X)
    !CALL assim_ensemble_statistics(assimData%A(i)%x%x, assimData%A(i)%x%nvar, assimData%A(i)%x%nens, assimData%A(i)%x%outmean, assimData%A(i)%x%outquant, assimData%A(i)%x%outmin, assimData%A(i)%x%outmax, assimData%A(i)%x%outsigma)
  ENDDO

END SUBROUTINE updateEnsembleStatistics

!----------------------------------------------------------------------
!>\brief Main routine in the ensemble Kalman filter analysis.
!!
!> Ensemble Kalman filter equations (Evensen) following (a) Mandel, J. Efficient implementation 
!> of the ensemble kalman filter, and (b) DeChant, C. Examining the effectiveness and 
!> robustness of sequential data assimilation methods for quantification of uncertainty
!> in hydrological forecasting. Localization following Magnusson, Gustafsson et al (2014)
!
! Organizes the calls to enkf_analysis_prepare and enkf_analysis_apply
! following this pseudo-code:
!   1. CALL enkf_analysis_prepare()
!        a. calculation of innovation, Y = D-HX
!        b. calculation and localization of the covariance of predicted observations cov(HX) = locCyy * cov(HX)
!        c. combine cov(HX) with observation error variance R into P = (cov(HX)+R)
!        d. derive intermediate update matrix M, by inversion of variance matrice sum  M = P^(-1)Y
!
!   2. DO i=1,num_model_states
!        CALL enkf_analysis_apply()
!             Here in Part 2, M is multiplied with the covariance Cov(X,HX) in the final EnKF update equation:
!      ENDDO
!----------------------------------------------------------------------------
SUBROUTINE enkf_analysis_main(assimData)
  !ARGUMENT
  TYPE(assim_data_type), INTENT(INOUT) :: assimData !<main assimilation variable containing all data
  !LOCAL VARIABLES
  INTEGER :: N, ND, NVAR, I, NX, J, NVARold, locID
  REAL, ALLOCATABLE :: X(:,:), Y(:,:), HA(:,:), M(:,:), D(:,:), HX(:,:), R(:) !, locCyy(:,:), locCxy(:,:)
  
  !Continue with ENKF analysis only if there is some data to assimilate
  IF(assimData%info%nD==0) RETURN   !No ENKF analysis if no data
  !IF(assimData%info%nD.GT.0)THEN   !CP161205 changed to return for nD=0 to have one less IF-ENDIF to keep track of
  !Part I:  calculate innovations and covariance matrix inversion
  !--------------------------------------------------------------
    !assign local variables
    ND = assimData%info%nD    !N:o observations to assimlate
    N  = assimData%info%NE    !N:o ensemble members
       
    !ALLOCATE necessary matrices
    ALLOCATE(Y(ND,N))
    ALLOCATE(M(ND,N))
    ALLOCATE(HA(ND,N))
    ALLOCATE(HX(ND,N))
    ALLOCATE(D(ND,N))
    ALLOCATE(R(ND))
    !ALLOCATE(locCyy(ND,ND))
      
    !get D, HX, and locCyy data from Enkf data structure (read from binary files if needed)
    !CALL assim_get_ensemble_data(ND,N,assimData%D,D,1)
    !CALL assim_get_ensemble_data(ND,N,assimData%HX,HX,1)
    !CALL assim_get_ensemble_data(ND,1,assimData%R,R,1)
    !CALL assim_get_ensemble_data(ND,ND,assimData%LocCYY,locCyy,1)
    D=assimData%D
    HX=assimData%HX
    R=assimData%R
  
    !locCyy =   assimData%LocCYY
    
    !CALL preparation routine
    CALL enkf_analysis_prepare(N,ND,D,HX,R,assimData%LocCYY,M,Y,HA)
  
    !Deallocate some temporary variables
    IF(ALLOCATED(HX)) DEALLOCATE(HX)
    IF(ALLOCATED(D)) DEALLOCATE(D)
    IF(ALLOCATED(R)) DEALLOCATE(R)
      
    !Part II: enkf update on model variables one-by-one
    !--------------------------------------------------
    !X state ensembles
    NX=assimData%Info%nX              !number of variable ensembles (types)
    IF(NX.GT.0)THEN
      NVAR=assimData%X(1)%x%nvar        !number of variables in ensemble (rows)
      !locID = assimData%X(i)%x%locID !localization ID      
      NVARold=0
      !Loop over number of variable ensembles, if >0
      DO I=1,NX
        !check if this variable is analysed or re-initialized
        IF(assimData%X(I)%x%assimilate)THEN
          !re-ALLOCATE X and locCxy matrices, if needed
          NVAR=assimData%X(I)%x%nvar        !number of variables in ensemble (rows)
          IF(ALLOCATED(X).AND.NVAR.ne.NVARold)DEALLOCATE(X)
          IF(.not.ALLOCATED(X))ALLOCATE(X(NVAR,N))
          !IF(ALLOCATED(locCxy).AND.NVAR.ne.NVARold)DEALLOCATE(locCxy)
          !IF(.not.ALLOCATED(locCxy))ALLOCATE(locCxy(NVAR,ND))
          NVARold = NVAR
   
          !Get ensemble data into X matrix
          CALL assim_get_ensemble_data(NVAR,N,assimData%X(I)%x,X)   !CP161202 for no bin-file or several bin-files

          !get localization matrix locCXY
          locID = assimData%X(i)%x%locID !localization ID      
          !CALL assim_get_ensemble_data(NVAR,ND,assimData%locCXY(locID),locCxy,1)
          !locCxy = assimData%locCXY

!          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY,X)  !CP161214 added aquifer extra coordinate system
          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY(locID)%x,X)
            
          !save updated X matrix back to the ensemble data
          CALL assim_set_ensemble_data(NVAR,N,assimData%X(I)%x,X,.true.,.true.)
        ENDIF
      ENDDO
    ENDIF
    
    !A auxilary ensembles (outvar in HYPE)
    NX=assimData%Info%nA              !number of variable ensembles (types)
    IF(NX.GT.0)THEN
      NVAR=assimData%A(1)%x%nvar        !number of variables in ensemble (rows)
      !locID = assimData%A(i)%x%locID !localization ID
      NVARold=0
        
      !Loop over number of variable ensembles, if >0
      DO I=1,NX
        !check if this variable is analysed or re-initialized
        IF(assimData%A(I)%x%assimilate)THEN
          !re-ALLOCATE X and locCxy matrices, if needed
          NVAR=assimData%A(I)%x%nvar        !number of variables in ensemble (rows)
          IF(ALLOCATED(X).AND.NVAR.ne.NVARold)DEALLOCATE(X)
          IF(.not.ALLOCATED(X))ALLOCATE(X(NVAR,N))
          !IF(ALLOCATED(locCxy).AND.NVAR.ne.NVARold)DEALLOCATE(locCxy)
          !IF(.not.ALLOCATED(locCxy))ALLOCATE(locCxy(NVAR,ND))
          NVARold = NVAR
   
          !Get ensemble data into X matrix
          CALL assim_get_ensemble_data(NVAR,N,assimData%A(I)%x,X)   !CP161206 for no bin-file or several bin-files

          !get localization matrix locCXY
          locID = assimData%X(i)%x%locID !localization ID      
          !CALL assim_get_ensemble_data(NVAR,ND,assimData%locCXY(locID),locCxy,1)
          !locCxy = assimData%locCXY

          !CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY,X)  !CP161214 added aquifer extra coordinate system
          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY(locID)%x,X)
            
          !save updated X matrix back to the ensemble data
          CALL assim_set_ensemble_data(NVAR,N,assimData%A(I)%x,X,.true.,.true.)

        ENDIF
      ENDDO
    ENDIF

    !F forcing ensembles (P, T, SW etc in HYPE)
    NX=assimData%Info%nF              !number of variable ensembles (types)
    IF(NX.GT.0)THEN
      NVAR=assimData%F(1)%x%nvar        !number of variables in ensemble (rows)
      !locID = assimData%F(i)%x%locID !localization ID
      NVARold=0
        
      !Loop over number of variable ensembles, if >0
      DO I=1,NX
        !check if this variable is analysed or re-initialized
        IF(assimData%F(I)%x%assimilate)THEN
          !re-ALLOCATE X and locCxy matrices, if needed
          NVAR=assimData%F(I)%x%nvar        !number of variables in ensemble (rows)
          IF(ALLOCATED(X).AND.NVAR.ne.NVARold)DEALLOCATE(X)
          IF(.not.ALLOCATED(X))ALLOCATE(X(NVAR,N))
          !IF(ALLOCATED(locCxy).AND.NVAR.ne.NVARold)DEALLOCATE(locCxy)
          !IF(.not.ALLOCATED(locCxy))ALLOCATE(locCxy(NVAR,ND))
          NVARold = NVAR
   
          !Get ensemble data into X matrix
          CALL assim_get_ensemble_data(NVAR,N,assimData%F(I)%x,X)

          !get localization matrix locCXY
          locID = assimData%X(i)%x%locID !localization ID      
          !CALL assim_get_ensemble_data(NVAR,ND,assimData%locCXY(locID),locCxy,1)
          !locCxy = assimData%locCXY

          !CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY,X)  !CP161214 added aquifer extra coordinate system
          CALL enkf_analysis_apply(N,NVAR,ND,Y,M,HA,assimData%locCXY(locID)%x,X)
            
          !save updated X matrix back to the ensemble data
          CALL assim_set_ensemble_data(NVAR,N,assimData%F(I)%x,X,.true.,.true.)
        ENDIF
      ENDDO
    ENDIF
    !To come : parameter ensembles

    !Deallocate
    IF(ALLOCATED(X))     DEALLOCATE(X)
    IF(ALLOCATED(Y))     DEALLOCATE(Y)
    IF(ALLOCATED(M))     DEALLOCATE(M)
    IF(ALLOCATED(HA))    DEALLOCATE(HA)
    !IF(ALLOCATED(locCyy))DEALLOCATE(locCyy)
    !IF(ALLOCATED(locCxy))DEALLOCATE(locCxy)
  !ENDIF
END SUBROUTINE enkf_analysis_main
  


!LOTTA TO DAVID FOR CHECK OUTPUT BEGIN
!These routines might write something that today code cannot write.
!Do we need to keep them?
!!-------------------------------------------------------------------------
!! enkf_ensemble_output
!! PURPOSE: print output from the Enkf ensembles:
!!  To be Modified to the new data structures, Feb-2015 /David
!!-------------------------------------------------------------------------
!subroutine enkf_ensemble_output()
!  ! mean, median, 95% quantiles, min, max, stdev in 7 text files:
!  !
!  ! write to these files... fid_enkf_min, fid_enkf_max, fid_enkf_median, fid_enkf_final, fid_enkf_mean, fid_enkf_min95, fid_enkf_max95,fid_enkf_std
!  !
!  ! LOCAL VARIABLES:
!  INTEGER :: errcode,i,j,k,strlen
!  CHARACTER(LEN=30000) :: outputstringMEDIAN,outputstringMIN95,outputstringMAX95
!  CHARACTER(LEN=30000) :: outputstringSTDEV
!  CHARACTER(LEN=20)    :: fmtstring
!
!  ! ADD DATA TO OUTPUTSTRING:
!  outputstringMEDIAN(:) =' '
!  outputstringMIN95(:)  =' '
!  outputstringMAX95(:)  =' '
!  outputstringSTDEV(:)  =' '
!    
!  ! 1) model state variables from the X ensemble
!  IF(EnkfInfo%XS)THEN
!    DO i=1,EnkFInfo%nX
!      WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfX%outquant(2,i)
!      WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfX%outquant(1,i)
!      WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfX%outquant(3,i)
!      !WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfX%outcovar(i,i))**0.5
!      WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfX%outsigma(i))
!      k=k+10
!    ENDDO
!  ENDIF
!  ! 2) model auxiliary variables from A ensemble
!  DO i=1,EnkFInfo%nA
!    WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfA%outquant(2,i)
!    WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfA%outquant(1,i)
!    WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfA%outquant(3,i)
!    !WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfA%outcovar(i,i))**0.5
!    WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfA%outsigma(i))
!    k=k+10
!  ENDDO
!    
!  ! 3) model parameters from P ensemble
!  DO i=1,EnkFInfo%nP
!    WRITE(outputstringMEDIAN(k:k+13),'(f14.6)') EnkfP%outquant(2,i)
!    WRITE(outputstringMIN95(k:k+13), '(f14.6)') EnkfP%outquant(1,i)
!    WRITE(outputstringMAX95(k:k+13), '(f14.6)') EnkfP%outquant(3,i)
!    !WRITE(outputstringSTDEV(k:k+13), '(f14.6)') (EnkfP%outcovar(i,i))**0.5
!    WRITE(outputstringSTDEV(k:k+13), '(f14.6)') (EnkfP%outsigma(i))
!    k=k+14
!  ENDDO
!   
!  ! 4) forcing from F ensemble
!!    DO i=1,EnkFInfo%nF
!!        WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfF%outquant(2,i)
!!        WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfF%outquant(1,i)
!!        WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfF%outquant(3,i)
!!!        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfF%outcovar(i,i))**0.5
!!        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfF%outsigma(i))
!!        k=k+10
!!    ENDDO
!
!  ! 5) "Innovations" (aka mod-obs) from D ensemble (check the enkf analysis scheme...)
!  !j=0
!  !DO i=1,EnkFInfo%nD
!  !    if ()THEN
!  !        j=j+1
!  !    WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') EnkfD%outquant(2,i)
!  !    WRITE(outputstringMIN95(k:k+9), '(f10.4)') EnkfD%outquant(1,i)
!  !    WRITE(outputstringMAX95(k:k+9), '(f10.4)') EnkfD%outquant(3,i)
!  !    WRITE(outputstringSTDEV(k:k+9), '(f10.4)') (EnkfD%outcovar(i,i))**0.5
!  !    k=k+10
!  !    ELSE
!  !        WRITE(outputstringMEDIAN(k:k+9),'(f10.4)') -9999.
!  !        WRITE(outputstringMIN95(k:k+9), '(f10.4)') -9999.
!  !        WRITE(outputstringMAX95(k:k+9), '(f10.4)') -9999.
!  !        WRITE(outputstringSTDEV(k:k+9), '(f10.4)') -9999.
!  !    ENDIF
!  !ENDDO
!
!
!  ! PREPARE OUTPUT FORMAT
!  outputstringMEDIAN = TRIM(adjustr(outputstringMEDIAN(1:k)))
!  outputstringMIN95 = TRIM(adjustr(outputstringMIN95(1:k)))
!  outputstringMAX95 = TRIM(adjustr(outputstringMAX95(1:k)))
!  outputstringSTDEV = TRIM(adjustr(outputstringSTDEV(1:k)))
!
!  fmtstring = enkf_get_fmtstring(outputstringMEDIAN)
!  WRITE(fid_enkf_median,TRIM(fmtstring))TRIM(outputstringMEDIAN)
!
!  fmtstring = enkf_get_fmtstring(outputstringMIN95)
!  WRITE(fid_enkf_min95,TRIM(fmtstring))TRIM(outputstringMIN95)
!
!  fmtstring = enkf_get_fmtstring(outputstringMAX95)
!  WRITE(fid_enkf_max95,TRIM(fmtstring))TRIM(outputstringMAX95)
!
!  fmtstring = enkf_get_fmtstring(outputstringSTDEV)
!  WRITE(fid_enkf_std,TRIM(fmtstring))TRIM(outputstringSTDEV)
!end subroutine enkf_ensemble_output

!!-------------------------------------------------------------------------
!! enkf_open_timeseriesoutput
!!-------------------------------------------------------------------------
!subroutine enkf_open_timeseriesoutput(resultdir)
!  CHARACTER(LEN=200) resultdir    
!  open(fid_enkf_median,file=TRIM(resultdir)//'enkf_output_median.txt',STATUS = 'unknown')
!  open(fid_enkf_min95,file=TRIM(resultdir)//'enkf_output_min95.txt',STATUS = 'unknown')
!  open(fid_enkf_max95,file=TRIM(resultdir)//'enkf_output_max95.txt',STATUS = 'unknown')
!  open(fid_enkf_std,file=TRIM(resultdir)//'enkf_output_std.txt',STATUS = 'unknown')
!  WRITE(fid_enkf_median,*)'% ENKF MEDIAN OUTPUT:'
!  IF(EnkfInfo%XS)THEN
!    WRITE(fid_enkf_median,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_min95,*)'% ENKF 2.5 quantile OUTPUT:'
!    WRITE(fid_enkf_min95,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_max95,*)'% ENKF 97.5 quantile OUTPUT:'
!    WRITE(fid_enkf_max95,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_std,*)'% ENKF standard deviation OUTPUT:'
!    WRITE(fid_enkf_std,*)'% Timestep X(1:nx) A(1:na) P(1:np) F(1:nf)'
!  ELSE
!    WRITE(fid_enkf_median,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_min95,*)'% ENKF 2.5 quantile OUTPUT:'
!    WRITE(fid_enkf_min95,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_max95,*)'% ENKF 97.5 quantile OUTPUT:'
!    WRITE(fid_enkf_max95,*)'% Timestep A(1:na) P(1:np) F(1:nf)'
!    WRITE(fid_enkf_std,*)'% ENKF standard deviation OUTPUT:'
!    WRITE(fid_enkf_std,*)'% Timestep A(1:na) P(1:np) F(1:nf)'    
!  ENDIF
!end subroutine enkf_open_timeseriesoutput
!
!subroutine enkf_close_timeseriesoutput()
!  close(fid_enkf_median)
!  close(fid_enkf_min95)
!  close(fid_enkf_max95)
!  close(fid_enkf_std)
!end subroutine enkf_close_timeseriesoutput

!FUNCTION assim_get_fmtstring(outputstri) RESULT(fmtstri)
!  CHARACTER(LEN=10000) outputstri
!  CHARACTER(LEN=20)    :: fmtstri
!  INTEGER strlen
!
!  fmtstri(1:3)=' (A'
!  strlen = LEN(TRIM(outputstri))
!
!  IF(strlen.lt.10)THEN
!    WRITE(fmtstri(4:4),'(i1)')strlen
!    fmtstri(5:6)=')'
!  ELSEIF(strlen.lt.100)THEN
!    WRITE(fmtstri(4:5),'(i2)')strlen
!    fmtstri(6:7)=')'
!  ELSEIF(strlen.lt.1000)THEN
!    WRITE(fmtstri(4:6),'(i3)')strlen
!    fmtstri(7:8)=')'
!  ELSEIF(strlen.lt.10000)THEN
!    WRITE(fmtstri(4:7),'(i4)')strlen
!    fmtstri(8:9)=')'
!  ELSE
!    WRITE(fmtstri(4:8),'(i5)')strlen
!    fmtstri(9:10)=')'
!  ENDIF
!END FUNCTION assim_get_fmtstring
!LOTTA TO DAVID FOR CHECK OUTPUT END


END MODULE ASSIMILATION_ROUTINES
