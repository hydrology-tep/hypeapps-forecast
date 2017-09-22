!> \file hype_indata.f90
!> Contains module hype_indata.

!>Handles specific input data (files) for HYPE model
MODULE HYPE_INDATA

  !Copyright 2014-2015 SMHI
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
  USE READWRITE_ROUTINES, ONLY : check_station, &
                                 check_obs_timeperiod, &
                                 prepare_read_matrix, &
                                 read_matrix_line,  &
                                 read_parameterline
  USE WORLDVAR, ONLY : fileunit_free, &
                       fileunit_get, &
                       get_seq_filename, &
                       readdaily
  !Subroutine uses modvar.

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: load_xoms_files,    &
            reload_xoms_files,  &
            close_xoms_files,   &
            get_current_xoms,   &
            load_data_for_regression_parameter_estimate, &
            deallocate_regest_input_variables, &
            set_regest_parameter

!>\name Variables for HYPE input files: XobsXOxx.txt
!>\{
! !>Variables to hold information about XOMS input data. These files hold
! !>observation time series of unknown variable xom0..9 and xos0..9
! !>These files use sequence number.
  INTEGER, PARAMETER :: max_xoms = 10     !<Maximum number of xom- resp xos-files (0-9)
  CHARACTER(LEN=20) :: filename(max_xoms,2) = RESHAPE(SOURCE=(/ &
      'XobsXOM0.txt','XobsXOM1.txt','XobsXOM2.txt','XobsXOM3.txt','XobsXOM4.txt', &
      'XobsXOM5.txt','XobsXOM6.txt','XobsXOM7.txt','XobsXOM8.txt','XobsXOM9.txt', &
      'XobsXOS0.txt','XobsXOS1.txt','XobsXOS2.txt','XobsXOS3.txt','XobsXOS4.txt', &
      'XobsXOS5.txt','XobsXOS6.txt','XobsXOS7.txt','XobsXOS8.txt','XobsXOS9.txt'/), &
      SHAPE=(/max_xoms,2/))    !<Filenames for Xoms-files
  INTEGER            :: num_xomsfiles     !<Number of xom- resp xos-files (max filenumber)
  INTEGER, PUBLIC    :: num_xoms(2)       !<Number of xom- resp xos-files each
  INTEGER, ALLOCATABLE :: xoms_funit(:,:) !<File unit of files with observations (nvar,2)
  INTEGER, ALLOCATABLE :: xoms_nstn(:,:)  !<Number of stations (columns) (nvar,2)
  INTEGER, ALLOCATABLE :: xomindex(:,:)   !<Index of observations of variable averaged over period (nvar,subbasin)
  INTEGER, ALLOCATABLE :: xosindex(:,:)   !<Index of observations of variable summed over period (nvar,subbasin)
  TYPE(DateType), ALLOCATABLE :: xoms_fbdate(:,:) !<File starting date within simulation period (nvar,2)
  TYPE(DateType), ALLOCATABLE :: xoms_fedate(:,:) !<File ending date within simulation period (nvar,2)
  REAL, ALLOCATABLE, PUBLIC :: xom(:,:)         !<Observations of current time step of variable averaged over period (nvar,nobs)
  REAL, ALLOCATABLE, PUBLIC :: xos(:,:)         !<Observations of current time step of variable summed over period (nvar,nobs)
!>\}

!> \name Variables for HYPE regression parameter estimates input data
!> \{
  !Variables for regional parameter estimation
  INTEGER              :: npar_regr       !<number of parameters to be estimated using a regression relationship
  INTEGER              :: npc_catdes      !<maximum number of catchment descriptors
  INTEGER              :: ngroup_cat      !<maximum number of catchment groups
  INTEGER, ALLOCATABLE :: indx_par(:)     !<contains info on whether a given parameter is regionally estimated and which parameter order it has (max_par)
  INTEGER, ALLOCATABLE :: ndes_reg(:,:)   !<number of descriptors used to regionalize each of the parameters in each group (ngroup_cat,npar_regr)
  INTEGER, ALLOCATABLE :: idx_catdes(:,:,:)  !<indices of the pca of catchment descriptors used for regionalization of the parameter (ngroup_cat,npar_regr,npc_catdes)
  REAL, ALLOCATABLE    :: wght_catdes(:,:,:) !<weights of the pca used to estimate the model parameter (ngroup_cat,npar_regr,npc_catdes)
  INTEGER, ALLOCATABLE :: clust_group(:)  !<group id of catchment after catchments are grouped based on their charactersistics (nsub)
  REAL, ALLOCATABLE    :: catdes(:,:)     !<catchment descriptors for catchments to be used in parameter regionalization (nsub,npc_catdes)
!> \}
  
CONTAINS

  !>Reads files with different observations (xom0..xom9,xos0..xos9)
  !----------------------------------------------------------------------
  SUBROUTINE load_xoms_files(dir,ns,bdate,edate,status)

    USE MODVAR, ONLY : basin
    
    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(IN)  :: ns            !<Number of subbasins
    TYPE(DateType), INTENT(IN)  :: bdate  !<Begin simulation date
    TYPE(DateType), INTENT(IN)  :: edate  !<End simulation date
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    !Local variables
    INTEGER i,j                     !loop-variables
    INTEGER subnr(ns)
    INTEGER funit   !move into check-subrutiner
    TYPE(DateType) fbdate(max_xoms,2),fedate(max_xoms,2)       !Begin and end date of file
    LOGICAL fex(max_xoms,2)
    LOGICAL notimefound(max_xoms,2)
    INTEGER numobsstn(max_xoms,2)
    INTEGER localindex(ns,max_xoms,2)

    !Initiation of data
    subnr = basin(:)%subid
    localindex = 0
    
    !Find files with observations in overlapping period, and check stations
    DO i = 1,max_xoms
      DO j = 1,2
        CALL get_seq_filename(filename(i,j))
        INQUIRE(FILE=TRIM(dir)//TRIM(filename(i,j)),EXIST=fex(i,j))
        IF(fex(i,j))THEN
          !Check if stations are correct and present, and get their index tables.
          funit = fileunit_get()
          CALL check_station(funit,TRIM(dir)//TRIM(filename(i,j)),ns,1,subnr,   &   !Skip one row
                           numobsstn(i,j),localindex(:,i,j),status)
          CALL fileunit_free(funit)
          IF(status==3)THEN
            fex(i,j) = .FALSE.
            status = 0
            CYCLE
          ELSEIF(status.NE.0)THEN
            RETURN
          ENDIF
        ENDIF
        IF(fex(i,j))THEN
          !Check time periods of observations
          funit = fileunit_get()
          CALL check_obs_timeperiod(funit,TRIM(dir)//TRIM(filename(i,j)),2,   &
                bdate,edate,fbdate(i,j),fedate(i,j),notimefound(i,j),status)
          CALL fileunit_free(funit)
          IF(status==2)THEN  !Shorter time period
            status = 0
            IF(bdate.GT.fedate(i,j) .OR. edate.LT.fbdate(i,j))THEN    !no overlap
              fex(i,j) = .FALSE.
              CYCLE
            ELSE                                      !some overlap
              fbdate(i,j) = MaxDate(bdate,fbdate(i,j))
              fedate(i,j) = MinDate(edate,fedate(i,j))
            ENDIF
          ELSEIF(status.NE.0)THEN
            RETURN
          ENDIF
        ENDIF
      ENDDO
    ENDDO
   
    !Any Xoms-files with relevant data?
    num_xoms = 0
    DO i = 1,max_xoms
      IF(fex(i,1)) num_xoms(1) = i
      IF(fex(i,2)) num_xoms(2) = i
    ENDDO
    IF(num_xoms(1)==0)THEN
      WRITE(6,*) 'XobsXOM: No special observation files found.'
    ENDIF
    IF(num_xoms(2)==0)THEN
      WRITE(6,*) 'XobsXOS: No special observation files found.'
    ENDIF
    num_xomsfiles = MAX(num_xoms(1),num_xoms(2))
    IF(num_xomsfiles==0) RETURN
    
    !Allocate space for reading XobsXOMS-files
    IF(.NOT.ALLOCATED(xoms_funit))  ALLOCATE(xoms_funit(num_xomsfiles,2))
    IF(.NOT.ALLOCATED(xoms_fbdate)) ALLOCATE(xoms_fbdate(num_xomsfiles,2))
    IF(.NOT.ALLOCATED(xoms_fedate)) ALLOCATE(xoms_fedate(num_xomsfiles,2))
    IF(.NOT.ALLOCATED(xom))         ALLOCATE(xom(num_xoms(1),MAXVAL(numobsstn(:,1))))
    IF(.NOT.ALLOCATED(xos))         ALLOCATE(xos(num_xoms(2),MAXVAL(numobsstn(:,2))))
    IF(.NOT.ALLOCATED(xomindex))    ALLOCATE(xomindex(num_xoms(1),ns))
    IF(.NOT.ALLOCATED(xosindex))    ALLOCATE(xosindex(num_xoms(2),ns))
    IF(.NOT.ALLOCATED(xoms_nstn))   ALLOCATE(xoms_nstn(num_xomsfiles,2))
    
    xoms_funit = 0
    !Load data for readdaily (always)
    DO i = 1,num_xomsfiles !Number of files
      DO j = 1,2
        IF(fex(i,j))THEN
          xoms_funit(i,j) = fileunit_get()
          xoms_fbdate(i,j) = fbdate(i,j)
          xoms_fedate(i,j) = fedate(i,j)
          xoms_nstn(i,j) = numobsstn(i,j)
          IF(j==1) xomindex(i,:) = localindex(:,i,j)
          IF(j==2) xosindex(i,:) = localindex(:,i,j)
          CALL prepare_read_matrix(xoms_funit(i,j),TRIM(dir)//TRIM(filename(i,j)),2,fbdate(i,j),status)
          IF(status.NE.0) RETURN
          WRITE(6,*) 'File ready: ', TRIM(dir)//TRIM(filename(i,j))
        ENDIF
      ENDDO
    ENDDO 
     
  END SUBROUTINE load_xoms_files

  !>Prepare files with different observations (xom0..xom9,xos0..xos9) for rereading
  !> No checks.
  !----------------------------------------------------------------------
  SUBROUTINE reload_xoms_files(dir,status)

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    !Local variables
    INTEGER i,j                     !loop-variables

    status = 0
    !Load data for readdaily (always)
    DO j = 1,2
      DO i = 1,num_xoms(j)
        IF(xoms_funit(i,j)>0)THEN
          CALL prepare_read_matrix(xoms_funit(i,j),TRIM(dir)//TRIM(filename(i,j)),2,xoms_fbdate(i,j),status)
          IF(status.NE.0) RETURN
          WRITE(6,*) 'File ready: ', TRIM(dir)//TRIM(filename(i,j))
        ENDIF
      ENDDO
    ENDDO 
     
  END SUBROUTINE reload_xoms_files

  !>Close files with observations
  !--------------------------------------------------------------------
  SUBROUTINE close_xoms_files()

    !Local variables
    INTEGER i,j

    !Load data for readdaily (always)
    DO j = 1,2
      DO i = 1,num_xoms(j)
        IF(xoms_funit(i,j)>0)THEN
          CLOSE(xoms_funit(i,j))    !file unit not freed (may be used again!)
        ENDIF
      ENDDO
    ENDDO 

  END SUBROUTINE close_xoms_files

  !>Get current observations from file
  !---------------------------------------------------------
  SUBROUTINE get_current_xoms(cd,n)

    USE MODVAR, ONLY : missing_value,nsub_basemodel

    TYPE(DateType), INTENT(IN) :: cd  !Current date
    INTEGER, INTENT(IN) :: n          !Number of subbasins
    
    !Local variables
    INTEGER i,j,k  !loop-variables
    TYPE(DateType) d
    REAL    x(nsub_basemodel)   !räcker det? knappast alltid väl

    !Initalize
    xom = missing_value
    xos = missing_value
    
    !Read XOM0..9
    j = 1
    DO i=1,num_xoms(j)
      IF(xoms_funit(i,j)>0)THEN
        IF(cd.LT.xoms_fbdate(i,j) .OR. cd.GT.xoms_fedate(i,j)) CYCLE
        CALL READ_MATRIX_LINE(xoms_funit(i,j),xoms_nstn(i,j),d,x(1:xoms_nstn(i,j)),missing_value,.TRUE.)  !This do not manage files with time in date.
        DO k=1,n
          IF(xomindex(i,k)>0) xom(i,k) = x(xomindex(i,k))
        ENDDO
      ENDIF
    ENDDO

    !Read XOS0..9
    j = 2
    DO i=1,num_xoms(j)
      IF(xoms_funit(i,j)>0)THEN
        IF(cd.LT.xoms_fbdate(i,j) .OR. cd.GT.xoms_fedate(i,j)) CYCLE
        CALL READ_MATRIX_LINE(xoms_funit(i,j),xoms_nstn(i,j),d,x(1:xoms_nstn(i,j)),missing_value,.TRUE.) 
        DO k=1,n
          IF(xosindex(i,k)>0) xos(i,k) = x(xosindex(i,k))
        ENDDO
      ENDIF
    ENDDO
    
  END SUBROUTINE get_current_xoms

  !>Reads three files with basin characteristics, and regression coefficients for parameters
  !----------------------------------------------------------------------
  SUBROUTINE load_data_for_regression_parameter_estimate(dir,nsbase,ns,indexarray,status)

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory (modeldir)
    INTEGER, INTENT(IN)  :: nsbase        !<Number of subbasins (model set-up)
    INTEGER, INTENT(IN)  :: ns            !<Number of subbasins (submodel)
    INTEGER, INTENT(IN)  :: indexarray(ns)   !<index for basemodel
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine

    !Local variables
    INTEGER i,j
    INTEGER funit
    INTEGER,ALLOCATABLE :: temp_clustg(:)                      !help variable for clust_group
    REAL,ALLOCATABLE    :: temp_pcacdes(:,:)                   !help variable for catdes

    !>\b Algorithm \n
    !>Initialise
    status = 0
    funit = fileunit_get()

    !>Read catchment group and catchment descriptor information (CatchGroup.txt,CatchDes.txt)
    OPEN(UNIT=funit,file=TRIM(dir)//'CatchGroup.txt',STATUS='old',ACTION='read',ERR=100)
    IF(.NOT.ALLOCATED(clust_group)) ALLOCATE(clust_group(nsbase))
    ngroup_cat=0
    DO i=1,nsbase
      READ(funit,*,ERR=101,END=102) clust_group(i)
      IF(clust_group(i) > ngroup_cat) ngroup_cat = clust_group(i)
    ENDDO
    CLOSE(funit)
    WRITE(6,*) 'File read: ', TRIM(dir)//'CatchGroup.txt'
    OPEN(UNIT=funit,file=TRIM(dir)//'CatchDes.txt',STATUS='old',ACTION='read',ERR=103)
    READ(funit,*) npc_catdes
    IF(.NOT.ALLOCATED(catdes)) ALLOCATE(catdes(nsbase,npc_catdes))
      DO i=1,nsbase
        READ(funit,*,ERR=104,END=105)(catdes(i,j),j=1,npc_catdes)
      ENDDO
    CLOSE(funit)
    CALL fileunit_free(funit)
    WRITE(6,*) 'File read: ', TRIM(dir)//'CatchDes.txt'
    
    !Reform catchment caracteristics and groups for regional parameter estimation
    IF(nsbase>ns)THEN
      IF(.NOT.ALLOCATED(temp_clustg)) ALLOCATE(temp_clustg(nsbase))
      temp_clustg = clust_group
      DEALLOCATE(clust_group)
      ALLOCATE(clust_group(ns))
      clust_group(:) = temp_clustg(indexarray(:))
      DEALLOCATE(temp_clustg)
    
      IF(.NOT.ALLOCATED(temp_pcacdes)) ALLOCATE(temp_pcacdes(nsbase,npc_catdes))
      temp_pcacdes = catdes
      DEALLOCATE(catdes)
      ALLOCATE(catdes(ns,npc_catdes))
      catdes(:,:) = temp_pcacdes(indexarray(:),:)
      DEALLOCATE(temp_pcacdes)
    ENDIF

    !>Read reg_par.txt
    CALL load_regression_parameters(dir,status) 
    IF(status/=0) RETURN
    
    RETURN
100 WRITE(6,*) 'ERROR: Open file ',TRIM(dir)//'CatchGroup.txt'
    status = 1
    RETURN
101 WRITE(6,*) 'ERROR: Reading file ',TRIM(dir)//'CatchGroup.txt'
    status = 1
    RETURN
102 WRITE(6,*) 'ERROR: File ending while reading: ',TRIM(dir)//'CatchGroup.txt'
    status = 1
    RETURN
103 WRITE(6,*) 'ERROR: Open file ',TRIM(dir)//'CatchDes.txt'
    status = 1
    RETURN
104 WRITE(6,*) 'ERROR: Reading file ',TRIM(dir)//'CatchDes.txt'
    status = 1
    RETURN
105 WRITE(6,*) 'ERROR: File ending while reading: ',TRIM(dir)//'CatchDes.txt'
    status = 1
    RETURN
     
  END SUBROUTINE load_data_for_regression_parameter_estimate

  !>Read parameter values from file
  !>
  !>/b Consequences The module variables npar_regr, ndes_reg, idx_catdes, 
  !> wght_catdes, indx_par will be allocated and set.
  !------------------------------------------------------------
  SUBROUTINE load_regression_parameters(dir,status) 

    USE WORLDVAR, ONLY : fileunit_temp,     &
                         comment_str
    USE MODVAR, ONLY : modparid,  &
                       max_par    

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER, INTENT(OUT) :: status        !<Status of subroutine
    
    !Local variables
    CHARACTER (LEN=10) varstr       !string with variable name
    CHARACTER(LEN=18000) line
    INTEGER i,j,k,l
    INTEGER nvalues                 !number of values read from line
    INTEGER funit
    REAL,ALLOCATABLE :: values(:)               !values of parameter read from file

    !>\b Algorithm \n
    !Initialisations
    status = 0
    funit = fileunit_get()

    !>Allocate and load regression parameter coefficients
      OPEN(unit = funit,FILE = TRIM(dir)//'reg_par.txt', STATUS = 'old', ACTION='read', ERR=100)
      READ(funit,*, ERR=101) npar_regr
     
      IF(.NOT.ALLOCATED(ndes_reg)) ALLOCATE(ndes_reg(ngroup_cat,npar_regr))
      IF(.NOT.ALLOCATED(indx_par)) ALLOCATE(indx_par(max_par))
      IF(.NOT.ALLOCATED(idx_catdes))  ALLOCATE(idx_catdes(ngroup_cat,npar_regr,npc_catdes))
      IF(.NOT.ALLOCATED(wght_catdes)) ALLOCATE(wght_catdes(ngroup_cat,npar_regr,npc_catdes))
      IF(.NOT.ALLOCATED(values))   ALLOCATE(values(npc_catdes))
      idx_catdes=0;ndes_reg=0;indx_par=0;wght_catdes=0.
    
      DO l = 1,ngroup_cat 
        DO i = 1,npar_regr
          READ(funit,'(a)',ERR=101) line
          IF(line(1:2)==comment_str) CYCLE
          CALL read_parameterline(line,npc_catdes,varstr,values,nvalues)    !read one line with parameter values
          IF(varstr=='          ')EXIT    !end of file

          !Find corresponding index
          DO j = 1,max_par
            IF(varstr==modparid(j)%shortname)THEN
              indx_par(j) = i
              ndes_reg(l,i) = nvalues
              wght_catdes(l,i,1:nvalues) = values(1:nvalues)             !!save weights of catchment descriptors to be used in the regional relationship
              READ(funit,*,ERR=101)line,(idx_catdes(l,i,k),k=1,nvalues)  !!read the index (position) of the catchment descriptors used for regionalization 
              EXIT
            ENDIF
          ENDDO 
        ENDDO
      ENDDO
      DEALLOCATE(values)
      CLOSE(funit)
      CALL fileunit_free(funit)
      WRITE(6,*) 'File read: ', TRIM(dir)//'reg_par.txt'

      RETURN
100 WRITE(6,*) 'ERROR: Open file ',TRIM(dir)//'reg_par.txt'
    status = 1
    RETURN
101 WRITE(6,*) 'ERROR: Reading file ',TRIM(dir)//'reg_par.txt'
    status = 1
    RETURN
  END SUBROUTINE load_regression_parameters

  !>Deallocate input data for regression parameter estimates.
  !------------------------------------------------------------
  SUBROUTINE deallocate_regest_input_variables()
  
    IF(ALLOCATED(ndes_reg)) DEALLOCATE(ndes_reg)
    IF(ALLOCATED(indx_par)) DEALLOCATE(indx_par)
    IF(ALLOCATED(idx_catdes))  DEALLOCATE(idx_catdes)
    IF(ALLOCATED(wght_catdes)) DEALLOCATE(wght_catdes)
    IF(ALLOCATED(clust_group)) DEALLOCATE(clust_group)
    IF(ALLOCATED(catdes)) DEALLOCATE(catdes)

  END SUBROUTINE deallocate_regest_input_variables

  !>Set regression parameter estimate from regression coefficients and
  !>catchment characteristics.
  !------------------------------------------------------------
  SUBROUTINE set_regest_parameter(isub,ipar,regpar,term)
  
    !Argument declarations
    INTEGER, INTENT(IN) :: isub   !<current subbasin
    INTEGER, INTENT(IN) :: ipar   !<current regression parameter
    REAL, INTENT(INOUT) :: regpar !<parameter value
    REAL, INTENT(IN),OPTIONAL :: term !<term to be added to parameter
    
    !Local variables
    INTEGER k
    
    IF(indx_par(ipar)>0) THEN
      regpar=0.0
      DO k=1,ndes_reg(clust_group(isub),indx_par(ipar))
        regpar = regpar + wght_catdes(clust_group(isub),indx_par(ipar),k)*catdes(isub,idx_catdes(clust_group(isub),indx_par(ipar),k))
      ENDDO
      IF(PRESENT(term)) regpar = regpar + term
    ENDIF

  END SUBROUTINE set_regest_parameter

END MODULE HYPE_INDATA
