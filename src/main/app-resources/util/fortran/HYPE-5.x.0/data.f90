!> \file data.f90
!> Contains module datamodule.

!> \brief Read and process model input from files, both static data
!> and forcing.
!>
!> Procedures for observations, for reading and preparing input data, 
!> for saving result files, and for loading and saving initial states
MODULE DATAMODULE
  !Copyright 2011-2017 SMHI
  !
  !This file is part of HYPE.
  !
  !HYPE is free software: you can redistribute it and/or modify it under
  !the terms of the Lesser GNU General Public License as published by
  !the Free Software Foundation, either version 3 of the License, or (at
  !your option) any later version. HYPE is distributed in the hope that
  !it will be useful, but WITHOUT ANY WARRANTY; without even the implied
  !warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See
  !the Lesser GNU General Public License for more details. You should
  !have received a copy of the Lesser GNU General Public License along
  !with HYPE. If not, see <http://www.gnu.org/licenses/>.
  !
  !--------------------------------------------------------------------
  USE COMPOUT, ONLY : find_variable_index_type
  USE CONVERT, ONLY : lower_case, &
                      string_convert_to_DateType    
  USE LIBDATE
  USE READWRITE_ROUTINES, ONLY : count_data_cols, &
                                 count_data_rows, &
                                 read_basindata5, &
                                 read_basindata6, &
                                 read_column_headings, &
                                 read_matrix,  &
                                 read_next_codestr_on_line, &
                                 read_parameterline
  USE TIMEROUTINES, ONLY : get_dayno_from_monthday, &
                           period_length, &
                           set_timestep_variables,  &
                           add_date_to_array
  !Subroutines also uses modvar, worldvar, modelmodule, optimization

  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------
  ! Private procedures 
  !----------------------------------------------
  ! load_one_forcing_variabel
  ! save_xoreginformation
  ! get_one_current_forcing_data
  ! get_current_flow
  ! get_current_otherobs
  ! get_current_oregobs
  ! calculate_special_optpar_parameters
  ! load_lakedata
  ! load_damdata
  ! start_lakedata_table
  ! finish_lakedata_table
  ! read_and_calc_basindata
  ! read_pointsourcedata
  ! load_irrigation_data
  ! read_update_data
  ! get_subid_index
  ! open_subbasinfiles
  ! close_subbasinfiles
  ! check_files_for_doubles
  ! open_timefiles
  ! close_timefiles
  ! get_parametervalues
  ! write_ensemble_simulations_heading
  PUBLIC ::  load_observations,&
             load_ascii_forcing,&
             load_ascii_qx_observations,&
             reset_observations,&
             prepare_for_update,&
             get_current_forcing,&
             close_observations,&
             load_basindata,&
             set_model_base_configuration,&
             set_model_configuration,&
             initiate_model_parameters,&
             load_aquiferdata,&
             load_glacierdata,&
             load_branchdata,&
             load_pointsourcedata,&
             get_current_pointsources, &
             calculate_path,&
             reform_inputdata_for_submodel,&
             load_cropdata,&
             get_hyss_arguments,&
             load_coded_info,&
             load_submodel_info,&
             set_outvar_for_variable,&
             prepare_subbasin_output,&
             save_mapfiles,&
             write_subbasin_assessment,&
             write_simulation_assessment,&
             prepare_outputfiles, &
             close_outputfiles, &
             write_subbasinfiles, &
             write_subbasinfiles_in_parallel, &
             write_regionfiles, &
             write_regionfiles_in_parallel, &
             write_timefiles, &
             write_timefiles_in_parallel, &
             load_parameters,&
             load_optpar,&
             save_respar,&
             initiate_output_routines,&
             initiate_outvar,&
             revise_outvar,&
             save_ensemble_simulations,&
             prepare_save_all_simulations,&
             write_simulation_results,&
             save_loadfiles, &
             checkindata_part1, &
             checkindata_part2, &
             checkindata_stop, &
             load_output_regions

CONTAINS
  
  !>Load forcing data files and other observations into the program
  !---------------------------------------------------------------------------------
  SUBROUTINE load_observations(dir,ns,bdate,edate,ndt,status) 
    
    USE WORLDVAR, ONLY : readformat
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir     !<File directory
    INTEGER, INTENT(IN)  :: ns              !<Number of subbasins, basemodel
    TYPE(datetype), INTENT(IN)  :: bdate    !<Begin time for simulation
    TYPE(datetype), INTENT(IN)  :: edate    !<End time for simulation
    INTEGER, INTENT(IN)  :: ndt             !<Number of timesteps
    INTEGER, INTENT(OUT) :: status          !<Error status of subroutine
    
    !> \b Algorithm \n
    !>Load observations from files
    IF(readformat==0)THEN   !ASCII-files
      CALL load_ascii_forcing(dir,ns,bdate,edate,ndt,status) 
      IF(status/=0) RETURN
      CALL load_ascii_qx_observations(dir,ns,ndt,bdate,edate,status) 
      IF(status/=0) RETURN
    ENDIF

  END SUBROUTINE load_observations
  
  !-----------------------------------------------------------------------------------------
  !>Reads forcing data from file. Also checks data.
  !!
  !>\b Consequences Module worldvar variables forcingdata and dates may change.
  !!Module modvar variables preci, tempi, snowfraci, shortwavei, windi, humidi, tmini, tmaxi may be allocated.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_ascii_forcing(dir,ns,bdate,edate,ndt,status) 
    
    USE WORLDVAR, ONLY : dates, &       !OUT
                         readdaily, &
                         get_seq_filename, &
                         forcingdata, &
                         max_forcingdata, &
                         i_pobs,i_tobs, &
                         i_rhobs,i_sfobs, &
                         i_swobs,i_uobs, &
                         i_tminobs,i_tmaxobs
    USE MODVAR, ONLY : preci, &
                       tempi, &
                       snowfraci, &
                       shortwavei, &
                       windi, &
                       humidi, &
                       tmini, &
                       tmaxi
                       
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    TYPE(DateType), INTENT(IN)  :: bdate !<Begin simulation date=bdate
    TYPE(DateType), INTENT(IN)  :: edate !<End simulation date=sdate
    INTEGER, INTENT(IN)  :: ndt          !<Number of timesteps in simulation
    INTEGER, INTENT(OUT) :: status       !<Error number
    
    !Local variables
    INTEGER iforc
    CHARACTER(LEN=16) :: fname

    status = 0

    !Initial allocation of model variables
    IF(forcingdata(i_pobs)%readfile) ALLOCATE(preci(ns))
    IF(forcingdata(i_tobs)%readfile) ALLOCATE(tempi(ns))
    IF(forcingdata(i_sfobs)%readfile) ALLOCATE(snowfraci(ns))
    IF(forcingdata(i_swobs)%readfile) ALLOCATE(shortwavei(ns))
    IF(forcingdata(i_uobs)%readfile)  ALLOCATE(windi(ns))
    IF(forcingdata(i_rhobs)%readfile) ALLOCATE(humidi(ns))
    IF(forcingdata(i_tminobs)%readfile) ALLOCATE(tmini(ns))
    IF(forcingdata(i_tmaxobs)%readfile) ALLOCATE(tmaxi(ns))
    IF(.NOT. readdaily)THEN
      IF(.NOT.ALLOCATED(dates)) ALLOCATE(dates(ndt))
    ENDIF
    
    !Set path and name of files to be used
    DO iforc = 1,max_forcingdata
      IF(forcingdata(iforc)%readfile)THEN
        fname = forcingdata(iforc)%filename
        CALL get_seq_filename(fname)
        forcingdata(iforc)%filepath = TRIM(dir)//TRIM(fname)
      ENDIF
    ENDDO
    
    !Check data and prepare files to be read, possibly read data to memory
    DO iforc = 1,max_forcingdata
      CALL load_one_forcing_variabel(ns,ndt,bdate,edate,forcingdata(iforc),status)
      IF(status/=0) RETURN
    ENDDO

  END SUBROUTINE load_ascii_forcing

  !>Load one observation forcing file; check file, saves data or prepare to read continously
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_one_forcing_variabel(ns,ndt,bdate,edate,forcobs,status) 
  
    USE WORLDVAR, ONLY : fileunit_temp, &
                         forcingdatatype, &
                         set_timeformat, &
                         dates, &       !OUT (ev.)
                         readdaily, &
                         readmatlab, &
                         fileunit_get, &
                         fileunit_free
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : check_obs_timeperiod, &
                                   check_pt_stn, &
                                   prepare_read_matrix
                             
    !Argument declarations
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: ndt          !<Number of timesteps in simulation
    TYPE(DateType), INTENT(IN)  :: bdate !<Begin simulation date=bdate
    TYPE(DateType), INTENT(IN)  :: edate !<End simulation date=sdate
    TYPE(FORCINGDATATYPE),INTENT(INOUT) :: forcobs !<Forcing data information and values
    INTEGER, INTENT(OUT) :: status       !<Error number
    
    !Local variables
    LOGICAL fileexist
    LOGICAL notimefound
    INTEGER nt                         !Number of simulation time steps
    TYPE(DateType) fbdate,fedate       !Begin and end date of file
  
    status = 0
    IF(forcobs%readfile)THEN
      
      !Check that file exist
      INQUIRE(FILE=TRIM(forcobs%filepath),EXIST=fileexist)
      IF(.NOT.fileexist)THEN
        WRITE(6,*) 'ERROR: Forcing data file is missing'
        WRITE(6,*) 'ERROR: ', TRIM(forcobs%filepath)
        status = 1
        RETURN
      ENDIF
      
      !Get and check stations
      IF(.NOT.ALLOCATED(forcobs%basinindex))ALLOCATE(forcobs%basinindex(ns))
      CALL check_pt_stn(fileunit_temp,forcobs%filepath,ns,forcobs%stationid,forcobs%basinindex,forcobs%ncols,status)
      IF(status.NE.0) RETURN
      IF(ALLOCATED(forcobs%stationid)) DEALLOCATE(forcobs%stationid)

      !Check time period
      CALL check_obs_timeperiod(fileunit_temp,forcobs%filepath,1,bdate,   &
             edate,fbdate,fedate,notimefound,status)
      IF(status.NE.0)THEN
        WRITE(6,*) 'ERROR: Forcing data times are missing'
        WRITE(6,*) 'ERROR: ', TRIM(forcobs%filepath)
        RETURN
      ENDIF  
      CALL set_timeformat(notimefound)
    
      !Load forcing data
      forcobs%fileunit = fileunit_get()
      CALL prepare_read_matrix(forcobs%fileunit,forcobs%filepath,1,bdate,status)
      IF(status.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(forcobs%filepath)
      IF(readdaily)THEN
      ELSE
        ALLOCATE(forcobs%allvalues(ndt,forcobs%ncols))
        CALL READ_MATRIX(forcobs%fileunit,ndt,forcobs%ncols,edate,nt,   &
             dates,forcobs%allvalues,missing_value,readmatlab) 
        CLOSE(forcobs%fileunit)
        CALL fileunit_free(forcobs%fileunit)
      ENDIF
    ENDIF
    
  END SUBROUTINE load_one_forcing_variabel

  !>Reads Qobs and Xobs data values in file
  !!
  !>\b Consequences Module worldvar variables qobs, xobs, xoregobs, xcol, xorcol, readqobs,
  !!qobsindex, numqobsstn, bqdate, eqdate, bxdate, exdate, bxrdate and exrdate may change.
  !!Module modvar variables qobsi, xobsi and xoregobsi may be allocated, and 
  !!xobsindex and xoregindex may change.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE load_ascii_qx_observations(dir,ns,ndt,bdate,edate,status) 

    USE WORLDVAR, ONLY : qobs,    & !OUT
         readqobs,     &  !OUT
         qobsindex,    &  !OUT
         numqobsstn,   &  !OUT
         xobs,    &   !OUT
         xoregobs,    &   !OUT
         readoutregion, &
         xcol,    &   !OUT
         xorcol,    &   !OUT
         readdaily,       &
         readmatlab,      &
         fileunit_qobs,   &
         fileunit_xobs,   &
         fileunit_xoreg,   &
         fileunit_temp,   &
         filename_Xobs,   &
         filename_Xoreg, &
         filename_Qobs,   &
         bqdate,  & !OUT
         eqdate,  & !OUT
         bxdate,  & !OUT
         exdate,  &  !OUT
         bxrdate,  & !OUT
         exrdate, &  !OUT
         noutreg
    USE MODVAR, ONLY : basin,     &
         qobsi,     & !OUT (allocated)
         xobsi,     & !OUT (allocated)
         xoregobsi,     & !OUT (allocated)
         xobsindex, & !OUT
         xoregindex, & !OUT
         missing_value,   &
         initiate_xobsinformation,  &
         save_xobsinformation
    USE READWRITE_ROUTINES, ONLY : check_obs_timeperiod, &
                                   check_q_stn, &
                                   check_xobs, &
                                   prepare_read_matrix

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir       !<File directory
    INTEGER, INTENT(IN)  :: ns                !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: ndt               !<Number of timesteps in simulation
    TYPE(DateType), INTENT(IN)   :: bdate     !<Begin simulation date
    TYPE(DateType), INTENT(IN)   :: edate     !<End simulation date
    INTEGER, INTENT(OUT) :: status            !<Error number
    
    !Local variables
    INTEGER nt,nskip
    TYPE(DateType) fbdate,fedate          !Begin and end date of file
    TYPE(DateType) tempdate
    LOGICAL qex
    LOGICAL xex,xorex
    LOGICAL notimefound
    REAL,ALLOCATABLE :: x2(:,:)           !Data, only used when reading all obs to memory
    REAL,ALLOCATABLE :: x3(:,:)           !Data, only used when reading all obs to memory
    INTEGER,ALLOCATABLE :: subnr(:)       !Help variable
    TYPE(DateType),ALLOCATABLE :: localdates(:) !Help variable
    INTEGER,ALLOCATABLE :: xvar(:,:),xorvar(:,:)      !Help variable (varid,subid or regid)

    status = 0
    !Initial allocation of model variables and local variables
    IF(.NOT. readdaily)THEN
      IF(.NOT.ALLOCATED(localdates)) ALLOCATE(localdates(ndt))
    ENDIF
    IF(.NOT.ALLOCATED(subnr)) ALLOCATE(subnr(ns))
    subnr = basin(:)%subid

    !Check if other observation stations are correct and present (only for Q), and get their index tables.
    INQUIRE(FILE=TRIM(dir)//filename_Qobs,EXIST=qex)
    IF(qex)THEN
      IF(.NOT.ALLOCATED(qobsindex)) ALLOCATE(qobsindex(ns))
      CALL check_q_stn(fileunit_temp,TRIM(dir)//filename_Qobs,ns,subnr,   &
            numqobsstn,qobsindex,status)
      IF(status==3)THEN
        qex = .FALSE.
      ELSEIF(status.NE.0)THEN
        RETURN
      ENDIF
    ENDIF
    INQUIRE(FILE=TRIM(dir)//filename_Xobs,EXIST=xex)
    IF(xex)THEN
      CALL count_data_cols(fileunit_temp,TRIM(dir)//filename_Xobs,1, xcol,status)
      IF(status.NE.0) RETURN
      xcol = xcol-1
      IF(xcol>0)THEN
        IF(.NOT.ALLOCATED(xvar)) ALLOCATE(xvar(xcol,2))
        CALL check_xobs(fileunit_temp,TRIM(dir)//filename_Xobs,xcol,    &
             xvar,status) 
        IF(status.NE.0) RETURN
      ENDIF
    ENDIF
    IF(readoutregion)THEN
      INQUIRE(FILE=TRIM(dir)//filename_Xoreg,EXIST=xorex)
      IF(xorex)THEN
        CALL count_data_cols(fileunit_temp,TRIM(dir)//filename_Xoreg,1,xorcol,status)
        IF(status.NE.0) RETURN
        xorcol = xorcol-1
        IF(xorcol>0)THEN
          IF(.NOT.ALLOCATED(xorvar)) ALLOCATE(xorvar(xorcol,2))
          CALL check_xobs(fileunit_temp,TRIM(dir)//filename_Xoreg,xorcol,    &
                 xorvar,status) 
          IF(status.NE.0) RETURN
        ENDIF
      ENDIF      
    ENDIF

    !Check time periods of observations
    IF(qex)THEN
      CALL check_obs_timeperiod(fileunit_qobs,TRIM(dir)//filename_Qobs,1,   &
            bdate,edate,fbdate,fedate,notimefound,status)
      bqdate = bdate
      eqdate = edate
      IF(status==2)THEN  !Shorter time period
        status = 0
        IF(bdate.GT.fedate .OR. edate.LT.fbdate)THEN    !no overlap
          qex = .FALSE.
        ELSE                                      !some overlap
          tempdate = MaxDate(bdate,fbdate)
          bqdate = tempdate
          tempdate = MinDate(edate,fedate)
          eqdate = tempdate
        ENDIF
      ELSEIF(status.NE.0)THEN
        RETURN
      ENDIF
    ENDIF
    IF(xex .AND. xcol>0)THEN
      CALL check_obs_timeperiod(fileunit_xobs,TRIM(dir)//filename_Xobs,3,   &
            bdate,edate,fbdate,fedate,notimefound,status)
      bxdate = bdate
      exdate = edate
      IF(status==2)THEN  !Shorter time period
        status = 0
        IF(bdate.GT.fedate .OR. edate.LT.fbdate)THEN    !no overlap
          xex = .FALSE.
          xcol = 0
        ELSE                                      !some overlap
          bxdate = MaxDate(bdate,fbdate)
          exdate = MinDate(edate,fedate)
        ENDIF
      ELSEIF(status.NE.0)THEN
        RETURN
      ENDIF
    ENDIF
    IF(readoutregion)THEN
      IF(xorex .AND. xorcol>0)THEN
        CALL check_obs_timeperiod(fileunit_xoreg,TRIM(dir)//filename_Xoreg,3,   &
              bdate,edate,fbdate,fedate,notimefound,status)
        bxrdate = bdate
        exrdate = edate
        IF(status==2)THEN  !Shorter time period
          status = 0
          IF(bdate.GT.fedate .OR. edate.LT.fbdate)THEN    !no overlap
            xorex = .FALSE.
            xorcol = 0
          ELSE                                      !some overlap
            bxrdate = MaxDate(bdate,fbdate)
            exrdate = MinDate(edate,fedate)
          ENDIF
        ELSEIF(status.NE.0)THEN
          RETURN
        ENDIF
      ENDIF
    ENDIF

    !Load discharge data
    IF(qex)THEN
      CALL prepare_read_matrix(fileunit_qobs,TRIM(dir)//filename_Qobs,1,bqdate,status)
      IF(status.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(dir)//filename_Qobs
      IF(.NOT.ALLOCATED(qobsi)) ALLOCATE(qobsi(ns))
      IF(readdaily)THEN
        readqobs = .TRUE.
      ELSE
        IF(.NOT.ALLOCATED(x2)) ALLOCATE(x2(ndt,numqobsstn))
        CALL READ_MATRIX(fileunit_qobs,ndt,numqobsstn,eqdate,nt,   &
             localdates,x2(:,:),missing_value,readmatlab) 
        CLOSE(fileunit_qobs)
        nskip = period_length(bdate,bqdate)
        ALLOCATE (qobs(nskip+1:nskip+nt,numqobsstn))
        qobs(:,:)=x2(1:nt,1:numqobsstn)
        DEALLOCATE(x2)
      ENDIF
    ENDIF

    !Load other observation data
    IF(xex .AND. xcol>0)THEN
      CALL prepare_read_matrix(fileunit_xobs,TRIM(dir)//filename_Xobs,3,bxdate,status)
      IF(status.NE.0) RETURN
      WRITE(6,*) 'File ready: ', TRIM(dir)//filename_Xobs
      IF(.NOT.ALLOCATED(xobsi)) ALLOCATE(xobsi(xcol))
      xobsi = missing_value
      IF(readdaily)THEN
      ELSE
        IF(.NOT.ALLOCATED(x3)) ALLOCATE(x3(ndt,xcol))
        CALL READ_MATRIX(fileunit_xobs,ndt,xcol,exdate,nt,localdates,    &
             x3(:,:),missing_value,readmatlab) 
        CLOSE(fileunit_xobs)
        nskip = period_length(bdate,bxdate)
        IF(.NOT.ALLOCATED(xobs)) ALLOCATE(xobs(nskip+1:nskip+nt,xcol))
        xobs(:,:) = x3(1:nt,:)
        DEALLOCATE(x3)
      ENDIF
    ENDIF
    CALL initiate_xobsinformation(xobsindex,ns)
    IF(xcol>0) CALL save_xobsinformation(xcol,xvar,ns)

    !Load outregion observation data
    IF(readoutregion)THEN
      IF(xorex .AND. xorcol>0)THEN
        CALL prepare_read_matrix(fileunit_xoreg,TRIM(dir)//filename_Xoreg,3,bxrdate,status)
        WRITE(6,*) 'File ready: ', TRIM(dir)//filename_Xoreg
        IF(status.NE.0) RETURN
        IF(.NOT.ALLOCATED(xoregobsi)) ALLOCATE(xoregobsi(xorcol))
        xoregobsi = missing_value
        IF(readdaily)THEN
        ELSE
          IF(.NOT.ALLOCATED(x3)) ALLOCATE(x3(ndt,xorcol))
          CALL READ_MATRIX(fileunit_xoreg,ndt,xorcol,exrdate,nt,localdates,    &
                 x3(:,:),missing_value,readmatlab) 
          CLOSE(fileunit_xoreg)
          nskip = period_length(bdate,bxrdate)
          IF(.NOT.ALLOCATED(xoregobs)) ALLOCATE(xoregobs(nskip+1:nskip+nt,xorcol))
          xoregobs(:,:) = x3(1:nt,:)
          DEALLOCATE(x3)
        ENDIF
      ENDIF
      CALL initiate_xobsinformation(xoregindex,noutreg)
      IF(xorcol>0) CALL save_xoreginformation(xorcol,xorvar,noutreg)      
    ENDIF

    IF(ALLOCATED(xvar))       DEALLOCATE(xvar)
    IF(ALLOCATED(subnr))      DEALLOCATE(subnr)
    IF(ALLOCATED(xorvar))     DEALLOCATE(xorvar)
    IF(ALLOCATED(localdates)) DEALLOCATE(localdates)

  END SUBROUTINE load_ascii_qx_observations

  !>Save information about finding other observation series
  !>
  !> \b Consequences Module index array (xoregindex) is set.
  !-------------------------------------------------------
  SUBROUTINE save_xoreginformation(n,var2,nr)
    
    USE MODVAR, ONLY : xoregindex, &
                       max_outvar
    USE WORLDVAR, ONLY : outregion

    INTEGER, INTENT(IN) :: n          !<Number of columns in Xoregobs
    INTEGER, INTENT(IN) :: var2(n,2)  !<Variables and subid for Xoregobs columns
    INTEGER, INTENT(IN) :: nr         !<Number of outregions
    
    !Local variables
    INTEGER ivar,ireg,icol

    DO ivar = 1, max_outvar
      DO ireg = 1, nr
        DO icol = 1,n
          IF(var2(icol,1)==ivar .AND. var2(icol,2)==outregion(ireg)%outregid) THEN
            xoregindex(ivar,ireg) = icol
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE save_xoreginformation

  !>Prepare forcingdata and other observation for another simulation
  !>in case they are read each time step
  !-----------------------------------------------------------------
  SUBROUTINE reset_observations(dir,status) 

    USE WORLDVAR, ONLY : xcol,    &
         xorcol, &
         readdaily,       &
         readqobs,        &
         readformat,      &
         bdate,           &
         bqdate,          &
         bxdate,bxrdate,  &
         fileunit_qobs,   &
         fileunit_xobs,   &
         fileunit_xoreg,  &
         filename_Xobs,   &
         filename_Xoreg,  &
         filename_Qobs,   &
         readoutregion,   &
         forcingdata, &
         max_forcingdata
    USE READWRITE_ROUTINES, ONLY : prepare_read_matrix
    USE MODELMODULE, ONLY : reload_modeldefined_observations

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir         !<File directory
    INTEGER, INTENT(OUT) :: status            !<Error number
    
    !Local variables
    LOGICAL xex
    INTEGER iforc

    IF(readdaily)THEN
      !Close files
      CALL close_observations(dir)
      IF(readformat==0)THEN
        DO iforc = 1, max_forcingdata   !Mandatory and optional forcing observations
          IF(forcingdata(iforc)%readfile)THEN
            CALL prepare_read_matrix(forcingdata(iforc)%fileunit,forcingdata(iforc)%filepath,1,bdate,status)
            IF(status.NE.0) RETURN
          ENDIF
        ENDDO
        IF(readqobs)THEN
          CALL prepare_read_matrix(fileunit_qobs,TRIM(dir)//filename_Qobs,1,bqdate,status)
          IF(status.NE.0) RETURN
        ENDIF
        INQUIRE(FILE=TRIM(dir)//filename_Xobs,EXIST=xex)    !Other observations
        IF(xex .and. xcol>0)THEN
          CALL prepare_read_matrix(fileunit_xobs,TRIM(dir)//filename_Xobs,3,bxdate,status)
          IF(status.NE.0) RETURN
        ENDIF
        IF(readoutregion .AND. xorcol>0)THEN    !Other observations
          CALL prepare_read_matrix(fileunit_xoreg,TRIM(dir)//filename_Xoreg,3,bxrdate,status)
          IF(status.NE.0) RETURN
        ENDIF
      ENDIF   !readformat
    ENDIF   !readdaily
    
    !Reset model defined and read observations
    CALL reload_modeldefined_observations(dir,status)
    IF(status.NE.0) RETURN

  END SUBROUTINE reset_observations

  !>Read information for basins to update
  !!
  !>\b Consequences Module modvar variables doupdate, updatestations, 
  !!updatestationsqar, updatestationwar,updatestationsarfact, updatetpcorr, updatetncorr, 
  !!updatetploccorr, updatetnloccorr, and wendupdstations may change.
  !-----------------------------------------------------------------------------------------
  SUBROUTINE prepare_for_update(dir,wobsvarname,locupall,locupnone,locupnone_qar,locupnone_war,wupall,wupnone,ns) 

    USE MODVAR, ONLY : doupdate,  & !OUT
         updatename,            &
         i_quseobs,             &
         i_qar,                 &
         i_war,                 &
         i_tpcorr,              &
         i_tncorr,              &
         i_tploccorr,           &
         i_tnloccorr,           &
         i_wendupd,             &
         updatestations,        & !OUT
         updatestationsqar,     & !OUT
         updatestationswar,     & !OUT
         updatestationsarfact,  & !OUT
         updatetpcorr,          & !OUT
         updatetncorr,          & !OUT
         updatetploccorr,       & !OUT
         updatetnloccorr,       & !OUT
         set_update_nutrientcorr, &
         basin,                 &
         xobsindex,             &
         wendupdstations,       & !OUT
         wobsvar,               &
         i_sp,i_pp,             &
         i_in,i_on,             &
         outvarid,              &
         max_outvar,            &
         nsub_basemodel
    USE WORLDVAR, ONLY : fileunit_temp,         &
         filename_upd,          &
         qobsindex

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER(LEN=4), INTENT(IN) :: wobsvarname !<Short-name of output variable that waterstage is updated with
    LOGICAL, INTENT(IN) :: locupall             !<Flag for quseobs update on all stations
    LOGICAL, INTENT(IN) :: locupnone            !<Flag for no quseobs update on any stations
    LOGICAL, INTENT(IN) :: locupnone_qar        !<Flag for no q_AR update on any stations !JD2012_AR	  
    LOGICAL, INTENT(IN) :: locupnone_war        !<Flag for no Q update with wAR on any stations
    LOGICAL, INTENT(IN) :: wupall               !<Flag for wendupd update on all stations
    LOGICAL, INTENT(IN) :: wupnone              !<Flag for no wendupd update on any stations
    INTEGER, INTENT(IN) :: ns                   !<Number of subbasins
    
    !Local variables
    INTEGER i, j            !loop-index
    INTEGER qusenum
    INTEGER qarnum
    INTEGER warnum
    INTEGER wendnum 
    INTEGER filenum 
    INTEGER qusetemp(nsub_basemodel)
    INTEGER qartemp(nsub_basemodel)
    INTEGER wartemp(nsub_basemodel)
    REAL    qarfac(nsub_basemodel)
    REAL    warfac(nsub_basemodel)
    INTEGER wendtemp(nsub_basemodel)
    REAL    tptemp(nsub_basemodel,2)
    REAL    tntemp(nsub_basemodel,2)
    REAL    tploctemp(nsub_basemodel,2)
    REAL    tnloctemp(nsub_basemodel,2)
    LOGICAL fileexist

    !Initiation
    qusenum = 0
    wendnum = 0
    filenum = 0
    qarnum  = 0
    warnum  = 0
    wobsvar = 0

    !Read update-file if available
    INQUIRE(FILE=TRIM(dir)//filename_upd,EXIST=fileexist) 
    IF(fileexist) CALL read_update_data(fileunit_temp,TRIM(dir)//filename_upd,    & 
         nsub_basemodel,qusenum,qusetemp,qarnum,warnum,qartemp,wartemp,qarfac,warfac,wendnum,  &
         wendtemp,filenum,tptemp,tntemp,tploctemp,tnloctemp)    

    !Dicharge data to use for update  (quseobs)
    IF(doupdate(i_quseobs))THEN
      IF(locupnone)THEN
        doupdate(i_quseobs) = .FALSE. 
      ELSEIF(locupall .OR. qusenum>0) THEN
        IF(.NOT.ALLOCATED(updatestations)) ALLOCATE(updatestations(ns))
        updatestations = .FALSE.             
        IF(locupall) THEN
          WHERE(qobsindex(1:ns)>0)
            updatestations(:) = .TRUE. 
          ENDWHERE
        ELSEIF(qusenum>0) THEN
          DO i = 1,ns
            DO j = 1, qusenum
              IF(basin(i)%subid == qusetemp(j)) THEN
                updatestations(i)=.TRUE.
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    !Dicharge data to use for AR-update  (qar)
    IF(doupdate(i_qar))THEN
      IF(locupnone_qar)THEN
        doupdate(i_qar) = .FALSE. 
      ELSEIF(qarnum>0) THEN
        IF(.NOT.ALLOCATED(updatestationsqar)) ALLOCATE(updatestationsqar(ns))
        IF(.NOT.ALLOCATED(updatestationsarfact)) ALLOCATE(updatestationsarfact(ns))
        updatestationsqar = .FALSE.             
        updatestationsarfact = 0
        DO i = 1,ns
          DO j = 1, qarnum
            IF(basin(i)%subid == qartemp(j)) THEN
              updatestationsqar(i)=.TRUE.
              updatestationsarfact(i)=qarfac(j)
              EXIT
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDIF

    !Waterstage data to use for AR-update  (war)
    IF(doupdate(i_war))THEN
      IF(locupnone_war)THEN
        doupdate(i_war) = .FALSE. 
      ELSEIF(warnum>0) THEN
        IF(.NOT.ALLOCATED(updatestationswar)) ALLOCATE(updatestationswar(ns))
        IF(.NOT.ALLOCATED(updatestationsarfact))THEN
          ALLOCATE(updatestationsarfact(ns))
          updatestationsarfact = 0
        ENDIF
        updatestationswar = .FALSE.             
        DO i = 1,ns
          DO j = 1, warnum
            IF(basin(i)%subid == wartemp(j)) THEN
              updatestationswar(i)=.TRUE.
              updatestationsarfact(i)=warfac(j)
              EXIT
            ENDIF
          ENDDO
        ENDDO

        !Check for missing input data
        DO wobsvar=1,max_outvar
          IF(outvarid(wobsvar)%shortname==wobsvarname)EXIT
        ENDDO
        IF(wobsvar>max_outvar)THEN
          WRITE(6,*) 'WARNING: observation variable for wAR-updating not found in info.txt file, no wAR updating'
          IF(ALLOCATED(updatestationswar)) DEALLOCATE(updatestationswar)
          doupdate(i_war) = .FALSE. 
        ENDIF
      ENDIF

    ENDIF
    
    !Correction of phosphorus concentration (tpcorr,tploccorr)
    IF(doupdate(i_tpcorr).OR.doupdate(i_tploccorr))THEN
      IF(i_sp==0.OR.i_pp==0)THEN
        WRITE(6,*) 'ERROR: Updating of phsosphorus not possible unless P is simulated'
        STOP 1
      ENDIF
      IF(filenum>0)THEN
        IF(doupdate(i_tpcorr))    CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tpcorr,tptemp,updatetpcorr)
        IF(doupdate(i_tploccorr)) CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tploccorr,tploctemp,updatetploccorr)
      ELSEIF(.NOT.fileexist)THEN
        WRITE(6,*) 'WARNING: missing update.txt file, necessary for update of phosphourus'
        doupdate(i_tpcorr) = .FALSE.
        doupdate(i_tploccorr) = .FALSE.
      ELSE
        WRITE(6,*) 'WARNING: missing data in update.txt file, no phosphourus updating'
        doupdate(i_tpcorr) = .FALSE.
        doupdate(i_tploccorr) = .FALSE.
      ENDIF
    ENDIF

    !Correction of nitrogen concentration out of subbasin (tncorr,tnloccorr)
    IF(doupdate(i_tncorr).OR.doupdate(i_tnloccorr))THEN
      IF(i_in==0.OR.i_on==0)THEN
        WRITE(6,*) 'ERROR: Updating of nitrogen not possible unless N is simulated'
        STOP 1
      ENDIF
      IF(filenum>0)THEN
        IF(doupdate(i_tncorr))    CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tncorr,tntemp,updatetncorr)
        IF(doupdate(i_tnloccorr)) CALL set_update_nutrientcorr(nsub_basemodel,ns,filenum,i_tnloccorr,tnloctemp,updatetnloccorr)
      ELSEIF(.NOT.fileexist)THEN
        WRITE(6,*) 'WARNING: missing update.txt file, necessary for update nitrogen'
        doupdate(i_tncorr) = .FALSE.
        doupdate(i_tnloccorr) = .FALSE.
      ELSE
        WRITE(6,*) 'WARNING: missing data in update.txt file, no nitrogen updating'
        doupdate(i_tncorr) = .FALSE.
        doupdate(i_tnloccorr) = .FALSE.
      ENDIF
    ENDIF

    !Update lake waterstage at end of time step (wendupd)
    IF(doupdate(i_wendupd))THEN
      IF(wupnone)THEN      
        doupdate(i_wendupd) = .FALSE. 
      ELSEIF(wupall .OR. wendnum>0) THEN
        IF(.NOT.ALLOCATED(wendupdstations)) ALLOCATE(wendupdstations(ns))
        wendupdstations = .FALSE.             
        IF(wupall) THEN
          wendupdstations(:) = .TRUE. 
        ELSEIF(wendnum>0) THEN
          DO i = 1,ns
            DO j = 1, wendnum
              IF(basin(i)%subid == wendtemp(j)) THEN
                wendupdstations(i)=.TRUE.
                EXIT
              ENDIF
            ENDDO
          ENDDO
        ENDIF
        !Check for missmatches
        DO wobsvar=1,max_outvar
          IF(outvarid(wobsvar)%shortname==wobsvarname)EXIT
        ENDDO
        IF(wobsvar>max_outvar)THEN
          WRITE(6,*) 'WARNING: missing variable in info.txt file, no wendupd updating'
          IF(ALLOCATED(wendupdstations)) DEALLOCATE(wendupdstations)
          doupdate(i_wendupd) = .FALSE. 
        ELSE
          DO i = 1,ns
            IF(wendupdstations(i).AND. xobsindex(wobsvar,i)==0) wendupdstations(i)=.FALSE.
          ENDDO
        ENDIF
      ELSE
        WRITE(6,*) 'WARNING: missing stations in update.txt, no wendupd updating'
        doupdate(i_wendupd) = .FALSE. 
      ENDIF
    ENDIF

    !Check for incompatible update methods (war not to be used with qar or wendupd)
    IF(doupdate(i_war).AND.doupdate(i_qar))THEN
      DO i = 1,ns
        IF(updatestationswar(i).AND.updatestationsqar(i))THEN
          updatestationswar(i)=.FALSE.
          WRITE(6,*)'WARNING: Update with q-AR and w-AR at the same subbasin is not allowed. w-AR turned off for subbasin',basin(i)%subid
          warnum = warnum - 1
        ENDIF
      ENDDO
      IF(warnum == 0) doupdate(i_war)=.FALSE.
    ENDIF
    IF(doupdate(i_war).AND.doupdate(i_wendupd))THEN
      DO i = 1,ns
        IF(updatestationswar(i).AND.wendupdstations(i))THEN
          updatestationswar(i)=.FALSE.
          WRITE(6,*)'WARNING: Update with wendupd and w-AR at the same subbasin is not allowed. w-AR turned off for subbasin',basin(i)%subid
          warnum = warnum - 1
        ENDIF
      ENDDO
      IF(warnum == 0) doupdate(i_war)=.FALSE.
    ENDIF

    IF(doupdate(i_quseobs).OR.doupdate(i_qar).OR.doupdate(i_war).OR.doupdate(i_wendupd)  &
       .OR.doupdate(i_tpcorr).OR.doupdate(i_tncorr).OR.doupdate(i_tploccorr).OR.doupdate(i_tnloccorr))THEN
      WRITE(6,*)
      WRITE(6,*)'-----Information about the simulation (cont.)----'
      IF(doupdate(i_quseobs))   WRITE(6,*)'Simulation with updated discharge, ',TRIM(updatename(i_quseobs))
      IF(doupdate(i_qar))       WRITE(6,*)'Simulation with qAR-updated discharge, ',TRIM(updatename(i_qar))
      IF(doupdate(i_war))       WRITE(6,*)'Simulation with wAR-updated discharge, ',TRIM(updatename(i_war))
      IF(doupdate(i_wendupd))   WRITE(6,*)'Simulation with updated lake water stage, ',TRIM(updatename(i_wendupd))
      IF(doupdate(i_tpcorr))    WRITE(6,*)'Simulation with updated TP-concentrations, ',TRIM(updatename(i_tpcorr))
      IF(doupdate(i_tncorr))    WRITE(6,*)'Simulation with updated TN-concentrations, ',TRIM(updatename(i_tncorr))
      IF(doupdate(i_tploccorr)) WRITE(6,*)'Simulation with updated TP-concentrations, ',TRIM(updatename(i_tploccorr))
      IF(doupdate(i_tnloccorr)) WRITE(6,*)'Simulation with updated TN-concentrations, ',TRIM(updatename(i_tnloccorr))
      WRITE(6,*)'-------------------------------------------------'
    ENDIF
  END SUBROUTINE prepare_for_update

  !>Get current forcing from file or memory
  !--------------------------------------------------------------------
  SUBROUTINE get_current_forcing(idt,ns,current_date)

    USE WORLDVAR, ONLY : forcingdata, &
                         i_pobs, &
                         i_tobs, &
                         i_tminobs, &
                         i_tmaxobs, &
                         i_rhobs, &
                         i_sfobs, &
                         i_swobs, &
                         i_uobs, &
                         numqobsstn
    USE MODVAR, ONLY :   preci,tempi, & !OUT
                         snowfraci,shortwavei, & !OUT
                         windi,humidi, & !OUT
                         tmini,tmaxi, & !OUT
                         qobsi,xobsi,xoregobsi  !OUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: idt  !<current timestep index
    INTEGER, INTENT(IN)  :: ns   !<number of subbasins
    TYPE(DateType), INTENT(OUT) :: current_date !<current date
    
      IF(forcingdata(i_pobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_pobs),current_date,preci(1:ns))
      IF(forcingdata(i_tobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_tobs),current_date,tempi(1:ns))
      IF(ALLOCATED(qobsi)) qobsi(1:ns) = get_current_flow(idt,current_date,ns,numqobsstn)
      IF(ALLOCATED(xobsi)) xobsi(:) = get_current_otherobs(idt,current_date)
      IF(ALLOCATED(xoregobsi)) xoregobsi(:) = get_current_oregobs(idt,current_date)
      IF(forcingdata(i_sfobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_sfobs),current_date,snowfraci(1:ns))
      IF(forcingdata(i_swobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_swobs),current_date,shortwavei(1:ns))
      IF(forcingdata(i_uobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_uobs),current_date,windi(1:ns))
      IF(forcingdata(i_rhobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_rhobs),current_date,humidi(1:ns))
      IF(forcingdata(i_tminobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_tminobs),current_date,tmini(1:ns))
      IF(forcingdata(i_tmaxobs)%readfile) CALL get_one_current_forcing_data(idt,ns,forcingdata(i_tmaxobs),current_date,tmaxi(1:ns))

  END SUBROUTINE get_current_forcing

  !>Get current forcing data for one variable for selected subbasins from file or memory
  !---------------------------------------
  SUBROUTINE get_one_current_forcing_data(i,n,forcobs,current_date,current_value)
  
    USE WORLDVAR, ONLY : readdaily, &
                         readformat, &
                         readmatlab, &
                         forcingdatatype, &
                         get_current_date_memory, &
                         get_current_forcing_from_memory
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line
    USE LIBDATE, ONLY : DateType
    
    !Argument declaration
    INTEGER, INTENT(IN) :: i  !<current time step
    INTEGER, INTENT(IN) :: n  !<number of subbasins
    TYPE(FORCINGDATATYPE), INTENT(IN) :: forcobs  !<data/information for one forcing variable
    TYPE(DateType), INTENT(OUT) :: current_date !<current date
    REAL, INTENT(OUT) :: current_value(n)     !<current value of the forcing variable
    
    !Variable declaration
    TYPE(DateType) d    !date
    REAL yt(forcobs%ncols)
    
    IF(readformat==0)THEN   !ASCII
      IF(readdaily)THEN
        CALL read_matrix_line(forcobs%fileunit,forcobs%ncols,d,yt,missing_value,readmatlab) 
        current_date = d
      ELSE
        current_date = get_current_date_memory(i)
        yt = get_current_forcing_from_memory(i,forcobs)
      ENDIF
      current_value = yt(forcobs%basinindex)
    ENDIF
    
  END SUBROUTINE get_one_current_forcing_data

  !>Get current discharge from file or memory
  !---------------------------------------------------------
  FUNCTION get_current_flow(i,cd,n,no) RESULT(current_flow)

    USE WORLDVAR, ONLY : readdaily,        &
                         readformat,       &
                         fileunit_qobs,    &
                         readmatlab,       &
                         qobsindex,        &
                         get_current_qobs, &
                         bqdate,           &
                         eqdate
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current time step
    TYPE(DateType), INTENT(IN) :: cd  !<current date
    INTEGER, INTENT(IN) :: n          !<number of subbasins
    INTEGER, INTENT(IN) :: no         !<number of observation stations
    REAL current_flow(n)              !<current flow
    
    !Local variables
    TYPE(DateType) d
    REAL    y(n),x(no)

    IF(cd.LT.bqdate .OR. cd.GT.eqdate)THEN
      current_flow = missing_value
      RETURN
    ENDIF

    IF(readformat==0.OR.readformat==4)THEN   !ASCII
      IF(readdaily)THEN
        y = missing_value
        CALL READ_MATRIX_LINE(fileunit_qobs,no,d,x,missing_value,readmatlab) 
        WHERE(qobsindex(1:n)>0)
          y=x(qobsindex(1:n))
        ENDWHERE
        current_flow = y
        IF(d.NE.cd) WRITE(*,*) 'ERROR read Qobs, date not equal to Pobs date'   !borde ej kunna inträffa
      ELSE
        current_flow = get_current_qobs(i,n)
      ENDIF
    ENDIF

  END FUNCTION get_current_flow

  !>Get current value of other observations from file or memory
  !------------------------------------------------------------
  FUNCTION get_current_otherobs(i,cd) RESULT(current_obs)

    USE WORLDVAR, ONLY : readdaily,        &
                         readformat,       &
                         fileunit_xobs,    &
                         get_current_xobs, &
                         readmatlab,       &
                         xcol,             &
                         bxdate,           &
                         exdate
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current time step
    TYPE(DateType), INTENT(IN) :: cd  !<current date 
    REAL current_obs(xcol)            !<current observations
    
    ! Local variables
    TYPE(DateType) d   !date
    REAL y(xcol)

    IF(xcol==0) RETURN        !no observations

    IF(readformat==0.OR.readformat==4)THEN   !ASCII
      current_obs = missing_value
      IF(cd.LT.bxdate .OR. cd.GT.exdate) RETURN
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_xobs,xcol,d,y,missing_value,readmatlab) 
        current_obs = y
        IF(d.NE.cd) WRITE(*,*) 'ERROR read Xobs, date not equal to Pobs date'   !borde ej kunna inträffa
      ELSE
        current_obs = get_current_xobs(i)
      ENDIF
    ENDIF

  END FUNCTION get_current_otherobs

  !>Get current value of other observations from file or memory
  !------------------------------------------------------------
  FUNCTION get_current_oregobs(i,cd) RESULT(current_obs)

    USE WORLDVAR, ONLY : readdaily,        &
                         readformat,       &
                         fileunit_xoreg,    &
                         get_current_xoregobs, &
                         readmatlab,       &
                         xorcol,             &
                         bxrdate,           &
                         exrdate
    USE MODVAR, ONLY : missing_value
    USE READWRITE_ROUTINES, ONLY : read_matrix_line

    !Argument declarations
    INTEGER, INTENT(IN) :: i          !<current time step
    TYPE(DateType), INTENT(IN) :: cd  !<current date 
    REAL current_obs(xorcol)          !<current observations
    
    ! Local variables
    TYPE(DateType) d   !date
    REAL y(xorcol)

    IF(xorcol==0) RETURN        !no observations

    IF(readformat==0.OR.readformat==4)THEN   !ASCII
      current_obs = missing_value
      IF(cd.LT.bxrdate .OR. cd.GT.exrdate) RETURN
      IF(readdaily)THEN
        CALL READ_MATRIX_LINE(fileunit_xoreg,xorcol,d,y,missing_value,readmatlab) 
        current_obs = y
        IF(d.NE.cd) WRITE(*,*) 'ERROR read Xoregobs, date not equal to Pobs date'   !borde ej kunna inträffa
      ELSE
        current_obs = get_current_xoregobs(i)
      ENDIF
    ENDIF

  END FUNCTION get_current_oregobs

  !>Close files with observations used for forcing and calibration
  !--------------------------------------------------------------------
  SUBROUTINE close_observations(dir)

    USE WORLDVAR, ONLY : fileunit_qobs,  &
                         fileunit_xobs,  &
                         filename_xobs,  &
                         fileunit_xoreg, &
                         xcol,xorcol,    &
                         readqobs,       &
                         readdaily,      &
                         readoutregion,  &
                         forcingdata, &
                         max_forcingdata

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory

    !Local variables
    INTEGER iforc
    LOGICAL xex

    IF(readdaily)THEN
      DO iforc = 1,max_forcingdata
        IF(forcingdata(iforc)%readfile)CLOSE(forcingdata(iforc)%fileunit)
      ENDDO
      IF(readqobs) CLOSE(fileunit_qobs)
      INQUIRE(FILE=TRIM(dir)//filename_Xobs,EXIST=xex)
      IF(xex .AND. xcol>0) CLOSE(fileunit_xobs)
      IF(readoutregion .AND. xorcol>0) CLOSE(fileunit_xoreg)
    ENDIF

  END SUBROUTINE close_observations

  !>Gets the information about soil-landuse-classes and subbasins from GeoClass and GeoData.
  !!
  !!\b Consequences Module worldvar variables vegtypereset may be set.
  !!Module modvar variables basin, classdata, wetland, nclass, nluse, nsoil, 
  !!modeloption, wetlandexist and numrotations are set. 
  !----------------------------------------------------------------
  SUBROUTINE load_basindata(dir,lfile,infile,n,status) 

    USE WORLDVAR, ONLY : fileunit_temp,   &
                         vegtypereset       !OUT
    USE MODVAR, ONLY : basin,             & !OUT
                       classdata,         & !OUT 
                       wetland,           & !OUT
                       nluse,nsoil,       & !OUT
                       nclass,            & !OUT
                       numsubstances,     &
                       set_soildepth,           &
                       set_coded_classes,       &
                       allocate_basinvariables, &
                       conductbasinpetmodel, &
                       petmodel, &
                       modeloption, &  !OUT
                       p_petmodel, &
                       wetlandexist, & !OUT
                       numrotations    !OUT
    USE READWRITE_ROUTINES, ONLY : read_geoclass

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir      !<File directory
    CHARACTER(LEN=*), INTENT(IN) :: lfile    !<Name of GeoClass file to be read
    CHARACTER(LEN=*), INTENT(IN) :: infile   !<Name of GeoData file to be read
    INTEGER, INTENT(OUT) :: n                !<Number of subbasins
    INTEGER, INTENT(OUT) :: status           !<Error status
    
    !Local parameters
    INTEGER,PARAMETER :: rcols = 50         !Number of data columns
    INTEGER,PARAMETER :: maxgccol = 13      !Maximum number of data columns in GeoClass
    INTEGER :: ncols                        !Total number of columns (rcols+number of s-lu-combinations (slc,dhslc,cr2) 
    
    !Local variables
    INTEGER i
    INTEGER mcols                   !Actual number of land-soil combinations in GeoData.txt
    REAL,ALLOCATABLE :: slcdata(:,:)   !data from Geoclass-file
    INTEGER numslc                  !Number of land-soil use combination in GeoClass.txt
    INTEGER temp                    !help variable
    INTEGER,ALLOCATABLE :: lakedataid(:)    !lakedataid from GeoData

    !> \b Algorithm \n
    !Start of subroutine
    status = 0

    !>Read the class characteristics
    CALL read_geoclass(fileunit_temp,TRIM(dir)//lfile,maxgccol,slcdata,numslc)

    nclass = numslc
    IF(.NOT.ALLOCATED(classdata)) ALLOCATE(classdata(nclass))
    nluse = NINT(MAXVAL(slcdata(1,1:numslc)))
    nsoil = NINT(MAXVAL(slcdata(2,1:numslc)))
    classdata(1:numslc)%luse  = NINT(slcdata(1,1:numslc))
    classdata(1:numslc)%soil  = NINT(slcdata(2,1:numslc))
    classdata(1:numslc)%crop  = NINT(slcdata(3,1:numslc))
    classdata(1:numslc)%crop2 = NINT(slcdata(4,1:numslc))
    classdata(1:numslc)%rotation = NINT(slcdata(5,1:numslc))
    numrotations = MAXVAL(classdata(:)%rotation)
    classdata(1:numslc)%vegtype = NINT(slcdata(6,1:numslc))
    classdata(1:numslc)%tiledepth = slcdata(8,1:numslc)
    classdata(1:numslc)%streamdepth = slcdata(9,1:numslc)
    DO i = 1, numslc 
      temp = set_soildepth(NINT(slcdata(10,i)),i,slcdata(11:10+NINT(slcdata(10,i)),i)) !default values from GeoClass.txt
    ENDDO
    !Check vegtype>0
    vegtypereset = .FALSE.
    DO i = 1,nclass
      IF(classdata(i)%vegtype==0)THEN
        classdata(i)%vegtype=1
        vegtypereset = .TRUE.
        WRITE(6,*) 'Warning: vegetation type not given in GeoClass.txt for slc-class',i,' vegtype=1 (open) is used'
      ENDIF
    ENDDO
    
    !Find coded classes
    CALL set_coded_classes(numslc,slcdata(7,1:numslc))
    IF(ALLOCATED(slcdata)) DEALLOCATE(slcdata)

    !>Count number of columns in GeoData
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,0,ncols,status)
    IF(status/=0)RETURN

    !>Count number of subbasins in GeoData
    CALL count_data_rows(fileunit_temp,TRIM(dir)//infile,1,n,status)
    IF(status/=0)RETURN

    !>Initiate subbasin data
    CALL allocate_basinvariables(n,numsubstances)
    basin(1:n)%ilakecatch = -1. !Flag for not read i GeoData
    wetland(1:n,1:2)%area = 0.
    petmodel = 0
    IF(.NOT.ALLOCATED(lakedataid)) ALLOCATE(lakedataid(n))
    lakedataid = 0

    !>Read the description of the basins (GeoData.txt)
    CALL read_and_calc_basindata(fileunit_temp,TRIM(dir)//infile,n,lakedataid,ncols,mcols)

    !Check if wetland is present
    IF(SUM(wetland(1:n,1:2)%area)==0)THEN
      wetlandexist = .FALSE.
      IF(ALLOCATED(wetland)) DEALLOCATE(wetland)
    ENDIF

    !Check if petmodel is present
    IF(SUM(petmodel(1:n))==0)THEN
      conductbasinpetmodel = .FALSE.
      DEALLOCATE(petmodel)
    ELSE
      conductbasinpetmodel = .TRUE.
      IF(modeloption(p_petmodel)>0)THEN
        WRITE(6,*) 'WARNING: Modeloption for petmodel given in info.txt will be ignored. petmodel from GeoData will be used.'
        modeloption(p_petmodel)=0
      ELSE
        WRITE(6,*)'Model option: Different PET model for each subbasin will be used'
      ENDIF
    ENDIF

    !>Read irrigation information of the basins (MgmtData.txt)
    CALL load_irrigation_data(fileunit_temp,TRIM(dir),status)
    IF(status/=0)RETURN

    !>Read observation coupling and elevation from separate file (ForcKey.txt)
    CALL load_forcing_information_data(fileunit_temp,TRIM(dir),n,status)
    IF(status/=0)RETURN

    !>Read additional description of important lakes (LakeData.txt) (possibly overwrite GeoData information)
    CALL load_lakedata(fileunit_temp,dir,n,lakedataid,status)
    IF(status/=0)RETURN
    IF(ALLOCATED(lakedataid)) DEALLOCATE(lakedataid)

    !>Read additional description of important dams (DamData.txt)
    CALL load_damdata(fileunit_temp,dir,n,status)
    IF(status/=0)RETURN

    !>Read additional description of flooding areas (FloodData.txt)
    CALL load_flooddata(fileunit_temp,dir,n,status)
    IF(status/=0)RETURN

  END SUBROUTINE load_basindata

  !>\brief Set model base configuration
  !>
  !>Module modvar variables dimriverlag is calculated for the submodel, 
  !>but the value for the base model is needed for reading statefile.
  !>
  !>\b Consequences Module modvar variables cropdata and cropirrdata 
  !>may be changed.
  !------------------------------------------------------------------
  SUBROUTINE set_model_base_configuration(nsbase,idir,dimrivlag,status)
  
  USE MODVAR, ONLY : conductregest
  USE MODELMODULE, ONLY : set_parameters_region_division, &
                          load_modeldefined_input, &
                          get_special_model_parameters, &
                          calculate_special_model_parameters
                          
  
  !Argument declaration
  INTEGER, INTENT(IN) :: nsbase  !<number of subbasin
  CHARACTER(LEN=*), INTENT(IN) :: idir  !<file directory
  INTEGER, INTENT(OUT) :: dimrivlag  !<max dimension of river queue
  INTEGER, INTENT(OUT) :: status  !<error status of subroutine
  
  !Local variables
  REAL optrivvel,optdamp
  INTEGER rvindex,dpindex   !index of rivvel and damp in modparid (not used here)
  INTEGER helparray(nsbase)
  TYPE(DateType):: date  !not used
  
  !Initialisations
  status = 0
  helparray = 0
  date = DateType(2000,1,1,0,0)

  !Find maximum dimension of river queue variables
  CALL load_modeldefined_input(idir,nsbase,nsbase,helparray,date,date,.FALSE.,conductregest,status) !read reg_par.txt
  IF(status/=0) RETURN
  CALL get_special_model_parameters(rvindex,dpindex)
  CALL calculate_special_optpar_parameters(idir,rvindex,dpindex,optrivvel,optdamp)  !optpar.txt
  CALL calculate_special_model_parameters(nsbase,optrivvel,optdamp,dimrivlag)  !par.txt,reg_par.txt
  
  !Reform crop input data for southern hemisphere
  CALL reform_input_data_for_pseudo_dayno(nsbase)

  END SUBROUTINE set_model_base_configuration

  !>Set model choice configuration
  !>
  !>\b Consequences Module modvar variables glacierexist and doirrigation is set
  !------------------------------------------------------------------
  SUBROUTINE set_model_configuration(glacexist,irrexist)
  
  USE MODVAR, ONLY : nclass,classmodel, &
                     deallocate_glacier
  USE MODELMODULE, ONLY : set_special_models
  
  !Argument declaration
  LOGICAL, INTENT(OUT) :: glacexist   !<status of glacier model
  LOGICAL, INTENT(OUT) :: irrexist    !<status of irrigation model
 
  CALL set_special_models(nclass,classmodel,glacexist,irrexist)

  !Deallocate glacier-variables
  IF(.NOT.glacexist) CALL deallocate_glacier()
  
  END SUBROUTINE set_model_configuration

  !>\brief Allocate model parameters and set their region division coupling.
  !>
  !>\b Consequences Module modvar arrays related to parameters are 
  !>allocated and set.
  !------------------------------------------------------------------
  SUBROUTINE initiate_model_parameters(ns_base,status)
  
  USE MODVAR, ONLY : allocate_modelparameters
  USE MODELMODULE, ONLY : set_parameters_region_division
                          
  
  !Argument declaration
  INTEGER, INTENT(IN) :: ns_base  !<number of subbasin
  INTEGER, INTENT(OUT) :: status  !<error status of subroutine
  
  !Local variable
  INTEGER nregpar  !number of regional parameters
  
  !Initialisations
  status = 0

  !Allocate model parameters and calculate number of regional parameters
  CALL allocate_modelparameters(ns_base,nregpar)

  !Allocate and set region division to the one used for each parameter.
  CALL set_parameters_region_division(nregpar)
  
  END SUBROUTINE initiate_model_parameters

  !>Get the information about branches from BranchData.txt
  !
  !>\b Consequences Module modvar variable branchsubid will be set
  !------------------------------------------------------------------
  SUBROUTINE load_branchdata(dir,status) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp, &
                         maxcharpath
    USE MODVAR, ONLY : branchsubid         !OUT

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(OUT) :: status        !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i
    INTEGER ncols                   !Total number of columns
    INTEGER mcols                   !Help variable
    INTEGER nrows                   !Number of branches in BranchData.txt
    CHARACTER(LEN=maxcharpath) filename
    INTEGER,ALLOCATABLE :: code(:)        !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)      !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)  !Content string

    status = 0
    filename = TRIM(dir)//'BranchData.txt'

    !> \b Algorithm \n
    !>Check if file exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN   

    !Count number of columns in BranchData.txt
    CALL count_data_cols(fileunit_temp,TRIM(filename),0,ncols,status)
    IF(status/=0)RETURN

    !Count number of rows in BranchData.txt
    CALL count_data_rows(fileunit_temp,TRIM(filename),1,nrows,status)
    IF(status/=0)RETURN

    !>Initiate branch data
    IF(.NOT.ALLOCATED(branchsubid)) ALLOCATE(branchsubid(nrows))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    branchsubid%source = 0
    branchsubid%branch = 0
    branchsubid%mainpart  = 1.
    branchsubid%maxQ = 0.
    branchsubid%minQ = 0.
    branchsubid%maxQbranch = 0.

    !>Read BranchData-file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(filename), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading heading in file: ',TRIM(filename)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type for columns
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='sourceid  ') code(i) = i_intg
      IF(str(i)(1:10)=='branchid  ') code(i) = i_intg
      IF(str(i)(1:10)=='mainpart  ') code(i) = i_real
      IF(str(i)(1:10)=='maxqmain  ') code(i) = i_real
      IF(str(i)(1:10)=='maxqbranch') code(i) = i_real
      IF(str(i)(1:10)=='minqmain  ') code(i) = i_real
    ENDDO

    CALL read_basindata5(fileunit_temp,filename,ncols,nrows,ncols,code,rindex,iindex,xi,xr)  !Read all data
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !>Save loaded branch data temporary
    DO i = 1,mcols
      IF(str(i)(1:10)=='sourceid  ')   branchsubid(1:nrows)%source     = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='branchid  ')   branchsubid(1:nrows)%branch     = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='mainpart  ')   branchsubid(1:nrows)%mainpart   = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='maxqmain  ')   branchsubid(1:nrows)%maxQ       = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='minqmain  ')   branchsubid(1:nrows)%minQ       = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='maxqbranch')   branchsubid(1:nrows)%maxQbranch = xr(1:nrows,rindex(i))
    ENDDO
    
    !>Check for illegal combinations of branch data and other illegal values
    DO i = 1,nrows
      IF(branchsubid(i)%maxQbranch>0)THEN
        IF(branchsubid(i)%maxQ>0)THEN
          branchsubid(i)%maxQbranch = 0.
          WRITE(6,*) 'WARNING: Maximum flow in both main channel and branch is given for source: ',branchsubid(i)%source
          WRITE(6,*) 'WARNING: This is not allowed. Maximum flow of branch will be removed.'
          WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
        ENDIF
      ENDIF
      IF(branchsubid(i)%maxQbranch<0)THEN
        branchsubid(i)%maxQbranch = 0.
        WRITE(6,*) 'WARNING: Maximum flow in branch is negative for source: ',branchsubid(i)%source
        WRITE(6,*) 'WARNING: This is not allowed. Maximum flow of branch will be removed.'
        WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
      ENDIF
      IF(branchsubid(i)%maxQ<0)THEN
        branchsubid(i)%maxQ = 0.
        WRITE(6,*) 'WARNING: Maximum flow in main channel is negative for source: ',branchsubid(i)%source
        WRITE(6,*) 'WARNING: This is not allowed. Maximum flow of main channel will be removed.'
        WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
      ENDIF
      IF(branchsubid(i)%minQ<0)THEN
        branchsubid(i)%minQ = 0.
        WRITE(6,*) 'WARNING: Minimum flow in main channel is negative for source: ',branchsubid(i)%source
        WRITE(6,*) 'WARNING: This is not allowed. Minimum flow of main channel will be removed.'
        WRITE(6,*) 'WARNING: You are suggested to correct your BranchData.txt file.'
      ENDIF
      IF(branchsubid(i)%mainpart<0. .OR. branchsubid(i)%mainpart>1.)THEN
        WRITE(6,*) 'ERROR: Main channel flow fraction is not between 0 and 1 for source: ',branchsubid(i)%source
        WRITE(6,*) 'ERROR: This is not allowed.'
        WRITE(6,*) 'ERROR: Fraction is needed for upstream area calculations.'
        WRITE(6,*) 'ERROR: You need to change your BranchData.txt file.'
        STOP 1
      ENDIF
    ENDDO

    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE load_branchdata

  !>Get the information about aquifers from AquiferData.txt
  !!
  !>\b Consequences Module modvar variable aquifer and pathsubid will be set. 
  !!Module modvar varible nregions may be increased.
  !------------------------------------------------------------------
  SUBROUTINE load_aquiferdata(dir,nsub,numaqu,status) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp, &
                         maxcharpath
    USE MODVAR, ONLY : basin,     &
                       pathsubid, &  !OUT
                       aquifer,   &  !OUT
                       nregions,  &  !OUT
                       modeloption, &
                       p_deepgroundwater

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir    !<File directory
    INTEGER, INTENT(IN)  :: nsub            !<Number of subbasins (basemodel)
    INTEGER, INTENT(OUT) :: numaqu          !<Number of aquifers
    INTEGER, INTENT(OUT) :: status          !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i,j
    INTEGER ncols                   !Total number of columns
    INTEGER mcols                   !Help variable
    INTEGER nrows                   !Number of rows in AquiferData.txt
    INTEGER subid_i,aquid_i,recha_i,outfr_i   !data column index for these columns
    INTEGER area_i,depth_i,basdp_i,topdp_i   !data column index for these columns
    INTEGER rate_i,delay_i,por_i,reg_i  !data column index for these columns
    INTEGER temp_i,coIN_i,coSP_i  !data column index for these columns
    INTEGER aqindex(nsub)
    CHARACTER(LEN=maxcharpath) filename
    INTEGER,ALLOCATABLE :: code(:)        !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)      !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)  !Content string

    status = 0
    numaqu = 0
    filename = TRIM(dir)//'AquiferData.txt'

    !Check if file exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(fex)THEN
      IF(modeloption(p_deepgroundwater)/=2)THEN
        WRITE(6,*) 'Warning: Existing AquiferData.txt will not be used, since modeloption not set.'
        RETURN
      ENDIF
    ELSE
      IF(modeloption(p_deepgroundwater)==2)THEN
        WRITE(6,*) 'Warning: Use modeloption with aquifer without AquiferData.txt not possible.'
        WRITE(6,*) 'Warning: Deep groundwater modeloption turned off.'
        modeloption(p_deepgroundwater)=0
        RETURN
      ENDIF
      RETURN   
    ENDIF

    !Count number of columns in AquiferData.txt
    CALL count_data_cols(fileunit_temp,TRIM(filename),0,ncols,status)
    IF(status/=0)RETURN

    !Count number of rows in AquiferData.txt
    CALL count_data_rows(fileunit_temp,TRIM(filename),1,nrows,status)
    IF(status/=0)RETURN

    !Initiate branch data
    IF(.NOT.ALLOCATED(pathsubid)) ALLOCATE(pathsubid(nsub))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    pathsubid%aquid = 0
    pathsubid%rechargebasin = .FALSE.
    pathsubid%recievefraction = 0.

    !Read AquiferData-file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(filename), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading heading in file: ',TRIM(filename)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type for columns
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='aquid     ') code(i) = i_intg
      IF(str(i)(1:10)=='parreg    ') code(i) = i_intg
      IF(str(i)(1:10)=='recharge  ') code(i) = i_intg
      IF(str(i)(1:10)=='area      ') code(i) = i_real
      IF(str(i)(1:10)=='inidepth  ') code(i) = i_real
      IF(str(i)(1:10)=='topdepth  ') code(i) = i_real
      IF(str(i)(1:10)=='basedepth ') code(i) = i_real
      IF(str(i)(1:10)=='porosity  ') code(i) = i_real
      IF(str(i)(1:10)=='retfrac   ') code(i) = i_real
      IF(str(i)(1:10)=='retrate   ') code(i) = i_real
      IF(str(i)(1:10)=='delay     ') code(i) = i_real
      IF(str(i)(1:10)=='temp      ') code(i) = i_real
      IF(str(i)(1:10)=='conc_in   ') code(i) = i_real
      IF(str(i)(1:10)=='conc_sp   ') code(i) = i_real
    ENDDO

    CALL read_basindata5(fileunit_temp,filename,ncols,nrows,ncols,code,rindex,iindex,xi,xr)  !Read all data
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !Find and check data columns
    subid_i=0;aquid_i=0;recha_i=0;area_i=0;depth_i=0;basdp_i=0;reg_i=0
    topdp_i=0;por_i=0;outfr_i=0;rate_i=0;delay_i=0;temp_i=0;coIN_i=0;coSP_i=0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')   subid_i = iindex(i)
      IF(str(i)(1:10)=='aquid     ')   aquid_i = iindex(i)
      IF(str(i)(1:10)=='parreg    ')   reg_i   = iindex(i)
      IF(str(i)(1:10)=='recharge  ')   recha_i = iindex(i)
      IF(str(i)(1:10)=='area      ')   area_i  = rindex(i)
      IF(str(i)(1:10)=='inidepth  ')   depth_i = rindex(i)
      IF(str(i)(1:10)=='basedepth ')   basdp_i = rindex(i)
      IF(str(i)(1:10)=='topdepth  ')   topdp_i = rindex(i)
      IF(str(i)(1:10)=='porosity  ')   por_i   = rindex(i)
      IF(str(i)(1:10)=='retfrac   ')   outfr_i = rindex(i)
      IF(str(i)(1:10)=='retrate   ')   rate_i  = rindex(i)
      IF(str(i)(1:10)=='delay     ')   delay_i = rindex(i)
      IF(str(i)(1:10)=='temp      ')   temp_i  = rindex(i)
      IF(str(i)(1:10)=='conc_in   ')   coIN_i  = rindex(i)
      IF(str(i)(1:10)=='conc_sp   ')   coSP_i  = rindex(i)
    ENDDO
    IF(subid_i==0 .OR. aquid_i==0 .OR. recha_i==0 .OR. area_i==0 .OR. &
       depth_i==0 .OR. basdp_i==0 .OR. por_i==0  .OR. reg_i==0 .OR. &
       outfr_i==0 .OR. rate_i==0  .OR. delay_i==0)THEN
      WRITE(6,*) 'Some columns are missing in AquiferData.txt'
      status=1
      RETURN
    ENDIF
    
    !Count number of aquifers
    numaqu = 0
    DO i = 1,nrows
      IF(xi(i,subid_i)==0)THEN  !Row with aquifer information
        numaqu = numaqu + 1
        aqindex(xi(i,aquid_i))=numaqu
      ENDIF
    ENDDO
    !Save loaded data
    ALLOCATE(aquifer(numaqu))
    aquifer%temperature = 0.
    aquifer%conc_IN = 0.
    aquifer%conc_SP = 0.
    DO i = 1,nrows
      IF(xi(i,subid_i)>0)THEN   !Row with aquifer subbasin coupling information
        DO j = 1,nsub
          IF(basin(j)%subid==xi(i,subid_i)) EXIT
        ENDDO
        IF(j<=nsub)THEN
          pathsubid(j)%aquid = aqindex(xi(i,aquid_i))
          pathsubid(j)%rechargebasin = (xi(i,recha_i)==1)
          pathsubid(j)%recievefraction = xr(i,outfr_i)
        ELSE
          WRITE(6,*) 'Warning: subid in AquiferData not in GeoData', xi(i,subid_i)
        ENDIF
      ELSEIF(xi(i,subid_i)==0)THEN  !Row with aquifer information
        aquifer(aqindex(xi(i,aquid_i)))%area = xr(i,area_i)
        aquifer(aqindex(xi(i,aquid_i)))%parregion(2) = xi(i,reg_i)
        aquifer(aqindex(xi(i,aquid_i)))%retrate = xr(i,rate_i)
        aquifer(aqindex(xi(i,aquid_i)))%percdelay = xr(i,delay_i)
        aquifer(aqindex(xi(i,aquid_i)))%inivol = (xr(i,depth_i)-xr(i,basdp_i))*xr(i,area_i)*xr(i,por_i)  !m3
        IF(topdp_i>0) aquifer(aqindex(xi(i,aquid_i)))%maxvol = (xr(i,topdp_i)-xr(i,basdp_i))*xr(i,area_i)*xr(i,por_i)  !m3
        aquifer(aqindex(xi(i,aquid_i)))%basedepth = xr(i,basdp_i)
        aquifer(aqindex(xi(i,aquid_i)))%porosity = xr(i,por_i)
        IF(temp_i>0) aquifer(aqindex(xi(i,aquid_i)))%temperature = xr(i,temp_i)
        IF(coIN_i>0) aquifer(aqindex(xi(i,aquid_i)))%conc_IN = xr(i,coIN_i)*1.E-3
        IF(coSP_i>0) aquifer(aqindex(xi(i,aquid_i)))%conc_SP = xr(i,coSP_i)*1.E-3
      ENDIF
    ENDDO
    
    !Check/Increase number of parameter regions
    !nregions(1) = MAX(nregions(1),MAXVAL(aquifer(1:numaqu)%parregion))
    nregions(2) = MAXVAL(aquifer(1:numaqu)%parregion(2))

    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE load_aquiferdata

  !>Get the information about glaciers from GlacierData.txt
  !!
  !>\b Consequences Module modvar variables glacier and glacierindex will be set. 
  !------------------------------------------------------------------
  SUBROUTINE load_glacierdata(dir,nsb,status) 

    USE WORLDVAR, ONLY : i_str, &
                         i_intg, &
                         i_real, &
                         bdate, &
                         fileunit_temp, &
                         maxcharpath, &
                         meandaysofyear
    USE MODVAR, ONLY : basin, &
                       glacier, &       !OUT
                       nglaciers, &     !OUT
                       glacierindex, &  !OUT
                       missing_value

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: nsb           !<Number of subbasins (basemodel)
    INTEGER, INTENT(OUT) :: status        !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i,j
    INTEGER ncols                   !Total number of columns
    INTEGER mcols                   !Help variable
    INTEGER nrows                   !Number of rows in AquiferData.txt
    INTEGER subid_i         !data column index for these columns
    INTEGER gyear,gmon,gday
    CHARACTER(LEN=maxcharpath) filename
    INTEGER,ALLOCATABLE :: glaciersubid(:)  !subid from GlacierData
    INTEGER,ALLOCATABLE :: glacdate(:)  !date of glacier slc information from GlacierData
    INTEGER,ALLOCATABLE :: code(:)        !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)      !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)      !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)  !Content string

    status = 0
    filename = TRIM(dir)//'GlacierData.txt'

    !Check if file exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN   

    !Count number of columns in GlacierData.txt
    CALL count_data_cols(fileunit_temp,TRIM(filename),0,ncols,status)
    IF(status/=0)RETURN

    !Count number of rows in GlacierData.txt
    CALL count_data_rows(fileunit_temp,TRIM(filename),1,nrows,status)
    IF(status/=0)RETURN

    !Initiate
    IF(.NOT.ALLOCATED(glaciersubid)) ALLOCATE(glaciersubid(nrows))
    IF(.NOT.ALLOCATED(glacdate)) ALLOCATE(glacdate(nrows))
    IF(.NOT.ALLOCATED(glacier))      ALLOCATE(glacier(nrows))
    IF(.NOT.ALLOCATED(glacierindex)) ALLOCATE(glacierindex(nsb))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    glacier(1:nrows)%gtype = INT(missing_value)
    glacier(1:nrows)%volcorr = 0.
    glacier(1:nrows)%yeardiff = 0.
    glacdate = INT(missing_value)  !date of SLC area (used to rescale initial volume)
    glacier(1:nrows)%glacinimb = 0.                 !annual mass balance for initial glacier volume correction

    !Read GlacierData-file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(filename), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading heading in file: ',TRIM(filename)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type for columns
    code=i_str    !string, ignore
    subid_i=0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        subid_i = i
      ENDIF
      IF(str(i)(1:10)=='glactype  ') code(i) = i_intg
      IF(str(i)(1:10)=='logvolcor ') code(i) = i_real
      IF(str(i)(1:10)=='slcdate   ') code(i) = i_intg
      IF(str(i)(1:10)=='annualmb  ') code(i) = i_real
    ENDDO

    !Check data columns
    IF(subid_i==0)THEN
      WRITE(6,*) 'SUBID columns are missing in GlacierData.txt'
      status=1
      RETURN
    ENDIF

    !Read data
    CALL read_basindata5(fileunit_temp,filename,ncols,nrows,ncols,code,rindex,iindex,xi,xr)  !Read all data
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !Set glacier characheristics array
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ') glaciersubid(1:nrows) = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='glactype  ') glacier(1:nrows)%gtype = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='logvolcor ') glacier(1:nrows)%volcorr = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='slcdate   ') glacdate(1:nrows) = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='annualmb  ') glacier(1:nrows)%glacinimb = xr(1:nrows,rindex(i))
    ENDDO
       
    !Calculate year difference between date of slc area information and simulation begin date
    DO j=1,nrows
      IF(glacdate(j)>0)THEN !years between SLC date and simulation start date
        gyear = INT(DBLE(glacdate(j)*1.d-4 + 0.1d0))
        gmon = INT((DBLE(glacdate(j)-DBLE(gyear)*1.d4)*1.d-2 + 0.1d0))
        gday = INT(DBLE(glacdate(j)-DBLE(gyear)*1.d4-DBLE(gmon)*1.d2 + 0.1d0))
        glacier(j)%yeardiff = TimeLag(bdate,DateType(gyear,gmon,gday,0,0))/meandaysofyear
      ENDIF
    ENDDO

    !Calculate glacier finding index array
    DO j = 1,nrows
      DO i = 1,nsb
        IF(basin(i)%subid==glaciersubid(j))THEN
          glacierindex(i) = j
          EXIT
        ENDIF
      ENDDO
    ENDDO
    
    nglaciers = nrows
        
    !Deallocate local arrays (for gfortran)
    IF(ALLOCATED(glaciersubid)) DEALLOCATE(glaciersubid)
    IF(ALLOCATED(glacdate)) DEALLOCATE(glacdate)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE load_glacierdata

  !>Collects the information about special lakes from LakeData.txt
  !>
  !>\b Consequences Module modvar variables lake, lakebasin, lakeindex and lakebasinindex
  !>may be allocated and set. Module modvar variables basin, nbasinlakes may be changed.
  !------------------------------------------------------------------
  SUBROUTINE load_lakedata(funit,dir,n,lakedataid,status) 

    USE WORLDVAR, ONLY : maxcharpath,  &
                         checkdata, &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &    !OUT
                         classbasin,   &
                         slc_olake,    &
                         missing_value,  &
                         lake,         &    !OUT
                         lakebasin,    &    !OUT
                         lakeindex,    &    !OUT
                         lakeout2index,  &  !OUT
                         lakebasinindex, &  !OUT
                         max_par,      &    
                         modparid,     &    
                         m_ldpar,      &
                         nbasinlakes        !OUT
    USE CONVERT,  ONLY : integer_convert_to_logical

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit         !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: n             !<Number of subbasins (base model)
    INTEGER, INTENT(IN)  :: lakedataid(n) !<lakedataid from GeoData
    INTEGER, INTENT(OUT) :: status        !<Error status
    
    !Local parameters
    INTEGER,PARAMETER :: i_simple = 1            !Lakedata type; simple olake
    INTEGER,PARAMETER :: i_combil = 2            !Lakedata type; lake composed of several lake basins
    INTEGER,PARAMETER :: i_lbasin = 3            !Lakedata type; lake basin
    INTEGER,PARAMETER :: i_llbasin = 4           !Lakedata type; last lake basin
    INTEGER,PARAMETER :: i_ldout1 = 5            !Lakedata type; outlet 1
    INTEGER,PARAMETER :: i_ldout2 = 6            !Lakedata type; outlet 2
    
    !Local variables
    LOGICAL fex
    INTEGER maxcol
    INTEGER i,j,k,ilake,ibasin
    INTEGER ictype                    !Column with lakedata type ('ldtype')
    INTEGER iclakeid                  !Column with lakeid
    INTEGER iclakedataid              !Column with id (lakedataid in GeoData)
    INTEGER ica                       !Column with area
    INTEGER icdepth                   !Column with depth
    INTEGER nrows                     !Number of lakes in file
    INTEGER ncols                     !Number of columns in file
    INTEGER nlake,nbasin,nlakefb,nlakedim
    INTEGER MonthDayDate              !Non-formated date from lakeData.txt
    INTEGER dayNumber                 !Corresponding date in day number format
    REAL    lakearea                  !lakearea calculated from GeoData
    REAL    parea
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)   !Content string
    INTEGER,ALLOCATABLE :: lakeindex2(:)      !row->isub
    INTEGER,ALLOCATABLE :: lakeindex3(:)      !row->ilake
    INTEGER,ALLOCATABLE :: lastlakeindex(:)   !ilake->row of last lakebasin for this lakebasinlake
    INTEGER,ALLOCATABLE :: qprod1date(:)
    INTEGER,ALLOCATABLE :: qprod2date(:)
    INTEGER,ALLOCATABLE :: lakeidindex(:,:)
    INTEGER,ALLOCATABLE :: code(:)            !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)          !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)          !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    REAL,ALLOCATABLE :: regvolume(:)          !Regulation volume (M m3)
    REAL,ALLOCATABLE :: wamp(:)               !Regulation amplitude (real) (m)

    status = 0
    nbasinlakes = 0
    filename = TRIM(dir)//'LakeData.txt'

    !Check if file exist and if lake-class exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(slc_olake==0)THEN
      WRITE(6,*) 'No class for olake to be used with LakeData.txt'
      status = 1
      RETURN
    ENDIF

    !Count number of columns and rows, and allocate variables for reading
    CALL count_data_cols(funit,filename,0,maxcol,status)
    IF(status/=0)THEN
      WRITE(6,*) 'ERROR: counting columns in LakeData'
      RETURN
    ENDIF
    CALL count_data_rows(funit,filename,1,nrows,status)
    IF(status/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No lakes in LakeData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(maxcol))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(maxcol))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(maxcol))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(maxcol))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))

    !Open LakeData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old', ACTION='read')
    CALL read_column_headings(funit,maxcol,str,ncols,status)
    IF(status.NE.0)THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    ica = 0
    ictype = 0
    iclakedataid = 0
    iclakeid = 0
    icdepth = 0

    DO i = 1,ncols
      IF(str(i)(1:10)=='area      ')THEN
        code(i) = i_real
        ica = i
      ENDIF
      IF(str(i)(1:10)=='ldtype    ')THEN
        code(i) = i_intg   !type of info: lake, basin or whole lake (that is divided in basins), outlet 1 or 2
        ictype = i
      ENDIF
      IF(str(i)(1:10)=='lakedataid')THEN
        code(i) = i_intg   !lakedataid in GeoData
        iclakedataid = i
      ENDIF
      IF(str(i)(1:10)=='lake_depth')THEN
        code(i) = i_real   !lake_depth
        icdepth = i
      ENDIF
      IF(str(i)(1:10)=='w0ref     ') code(i) = i_real   !w0
      IF(str(i)(1:10)=='deltaw0   ') code(i) = i_real   !deltaw0
      IF(str(i)(1:10)=='rate      ') code(i) = i_real   !lake_rate
      IF(str(i)(1:10)=='exp       ') code(i) = i_real   !lake_exp
      IF(str(i)(1:10)=='regvol    ') code(i) = i_real   !regvol
      IF(str(i)(1:10)=='wamp      ') code(i) = i_real   !reg. amplitude
      IF(str(i)(1:10)=='qprod1    ') code(i) = i_real   !two date-based production rates, qprod1 and qprod2, now replacing the unique rate qprod
      IF(str(i)(1:10)=='qprod2    ') code(i) = i_real
      IF(str(i)(1:10)=='maxqprod  ') code(i) = i_real   !max qprod
      IF(str(i)(1:10)=='minflow   ') code(i) = i_intg   !flag for minimum flow
      IF(str(i)(1:10)=='datum1    ') code(i) = i_intg   !two dates for the corresponding two dam production flows (date of start)               
      IF(str(i)(1:10)=='datum2    ') code(i) = i_intg
      IF(str(i)(1:10)=='qamp      ') code(i) = i_real   !lake_qamp
      IF(str(i)(1:10)=='qpha      ') code(i) = i_real   !lake_qpha
      IF(str(i)(1:10)=='limqprod  ') code(i) = i_real   !limprod
      IF(str(i)(1:10)=='lakeid    ')THEN
        code(i) = i_intg
        iclakeid = i
      ENDIF
      !All model parameters that can be given in LakeData is real
      DO j=1,max_par
        IF(str(i)(1:10)==modparid(j)%shortname .AND. modparid(j)%deptype==m_ldpar) code(i) = i_real
      ENDDO
    ENDDO

    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    !Count number of laketypes: lakes, lakes formed by basins and lakebasins.
    IF(ictype==0)THEN
      WRITE(6,*) 'ERROR: ldtype not found in LakeData. ',TRIM(filename)
      status = 1
      RETURN
    ELSE
      nlake = 0
      nbasin = 0
      nlakefb = 0
      DO j = 1,nrows
        IF(xi(j,iindex(ictype))==i_simple) nlake = nlake + 1
        IF(xi(j,iindex(ictype))==i_ldout1 .OR. xi(j,iindex(ictype))==i_ldout2) nlake = nlake + 1
        IF(xi(j,iindex(ictype))==i_combil) nlakefb = nlakefb + 1
        IF(xi(j,iindex(ictype))==i_lbasin.OR.xi(j,iindex(ictype))==i_llbasin) nbasin = nbasin + 1
      ENDDO
    ENDIF
    nlakedim = nlake+nlakefb

    !Allocate variables for lake information
    IF(.NOT.ALLOCATED(lakebasin)) ALLOCATE(lakebasin(nbasin))
    IF(.NOT.ALLOCATED(lake)) ALLOCATE(lake(nlakedim))
    IF(.NOT.ALLOCATED(regvolume)) ALLOCATE(regvolume(nlakedim))
    IF(.NOT.ALLOCATED(wamp)) ALLOCATE(wamp(nlakedim))
    IF(.NOT.ALLOCATED(qprod1date)) ALLOCATE(qprod1date(nlakedim)) 
    IF(.NOT.ALLOCATED(qprod2date)) ALLOCATE(qprod2date(nlakedim))
    regvolume = 0.; wamp = missing_value
    qprod1date = 0; qprod2date = 0  

    !Find subbasin-lake coupling: lakeindex(isb),lakebasinindex(isb), lakeindex2(ldrow) and lakeindex3(ldrow)
    IF(iclakedataid==0)THEN
      WRITE(6,*) 'ERROR: id not found in LakeData. ',TRIM(filename)
      status = 1
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(lakeindex)) ALLOCATE(lakeindex(n))
    IF(.NOT.ALLOCATED(lakeout2index)) ALLOCATE(lakeout2index(n))
    IF(nbasin>0 .AND.(.NOT.ALLOCATED(lakebasinindex))) ALLOCATE(lakebasinindex(n))
    IF(.NOT.ALLOCATED(lakeindex2)) ALLOCATE(lakeindex2(nrows))
    IF(.NOT.ALLOCATED(lakeindex3)) ALLOCATE(lakeindex3(nrows))
    IF(.NOT.ALLOCATED(lastlakeindex)) ALLOCATE(lastlakeindex(nrows))
    lakeindex = 0
    lakeout2index = 0
    IF(nbasin>0) lakebasinindex = 0
    lakeindex2 = 0
    lakeindex3 = 0
    lastlakeindex = 0
    DO i = 1,n
      IF(lakedataid(i)==0)CYCLE
      DO j = 1,nrows
        IF(lakedataid(i)==xi(j,iindex(iclakedataid)))THEN
          lakeindex2(j)=i
          IF(xi(j,iindex(ictype))==i_ldout1)THEN
            !Check for second outlet
            DO k=1,nrows
              IF(xi(j,iindex(iclakeid))==xi(k,iindex(iclakeid))) lakeindex2(k)=i
            ENDDO
          ELSEIF(xi(j,iindex(ictype))==i_ldout2)THEN
            !lakedataid should only be given for main outlet
            WRITE(6,*) 'WARNING: lakedataid is given for second outlet of lake in LakeData.txt'
            WRITE(6,*) 'Lake in subbasin ',basin(i)%subid, '.'
            WRITE(6,*) 'This may lead to lakedata-parameters are not read correctly.'
            WRITE(6,*) 'Please change LakeData.txt to have lakedataid=0 for ldtype=6.'
            DO k=1,nrows
              IF(xi(j,iindex(iclakeid))==xi(k,iindex(iclakeid))) lakeindex2(k)=i  !for safety find first outlet
            ENDDO
          ENDIF
          EXIT
        ENDIF
      ENDDO
      IF(j==nrows+1)THEN
        WRITE(6,*) 'ERROR: A lake in GeoData is missing in LakeData. lakedataid: ',lakedataid(i)
        status = 1
        RETURN
      ENDIF  
    ENDDO

    !First take care of the lakes that are formed by lakebasins
    IF(nlakefb>0)THEN
      IF(ica==0.OR.iclakeid==0)THEN             !Check for necessary data
        WRITE(6,*) 'ERROR: area or lakeid not found in LakeData. ',TRIM(filename)
        WRITE(6,*) 'ERROR: These are necessary for lakes formed by basins (ldtype 3).'
        status = 1
        RETURN
      ENDIF
      IF(.NOT.ALLOCATED(lakeidindex)) ALLOCATE(lakeidindex(nlakefb,2)) !lakeid-jrow coupling for later use in subroutine
      ilake = 0
      DO j = 1,nrows
        IF(xi(j,iindex(ictype))==i_combil)THEN
          ilake = ilake + 1
          lakeidindex(ilake,1) = xi(j,iindex(iclakeid))
          lakeidindex(ilake,2) = j
        ENDIF
      ENDDO
      DO ilake = 1,nlakefb                      !Save information permanent or temporary
        DO i = 1,ncols
          IF(str(i)=='area      ')   lake(ilake)%area      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='w0ref     ')   lake(ilake)%w0ref     = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='deltaw0   ')   lake(ilake)%deltaw0   = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qprod1    ')   lake(ilake)%qprod1    = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qprod2    ')   lake(ilake)%qprod2    = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='maxqprod  ')   lake(ilake)%mqprod    = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='minflow   ')   lake(ilake)%minflow   = integer_convert_to_logical(-xi(lakeidindex(ilake,2),iindex(i)))
          IF(str(i)=='datum1    ')   qprod1date(ilake)     = xi(lakeidindex(ilake,2),iindex(i))
          IF(str(i)=='datum2    ')   qprod2date(ilake)     = xi(lakeidindex(ilake,2),iindex(i))
          IF(str(i)=='rate      ')   lake(ilake)%rate      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='exp       ')   lake(ilake)%exp       = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qamp      ')   lake(ilake)%qamp      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='qpha      ')   lake(ilake)%qpha      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='limqprod  ')   lake(ilake)%limprod   = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='regvol    ')   regvolume(ilake)      = xr(lakeidindex(ilake,2),rindex(i))
          IF(str(i)=='wamp      ')   wamp(ilake)           = xr(lakeidindex(ilake,2),rindex(i))
        ENDDO
        lakeindex3(lakeidindex(ilake,2))=ilake
      ENDDO
    ENDIF

    !Simple lakes and lakebasins
    ilake = nlakefb
    ibasin = 0
    DO j = 1,nrows
      IF(lakeindex2(j)==0) CYCLE    !skip rows with lakes formed by basins and lakes not included in simulation
      IF(xi(j,iindex(ictype))==i_simple)THEN     !single lake
        ilake = ilake + 1
        lake(ilake)%area = xr(j,rindex(ica))    !is this used? remove?
        IF(icdepth>0)THEN
          IF(xr(j,rindex(icdepth)).NE.missing_value) basin(lakeindex2(j))%lakedepth(2) = xr(j,rindex(icdepth))
        ENDIF
        DO i = 1,ncols
          IF(str(i)=='w0ref     ')   lake(ilake)%w0ref    = xr(j,rindex(i))
          IF(str(i)=='deltaw0   ')   lake(ilake)%deltaw0  = xr(j,rindex(i))
          IF(str(i)=='qprod1    ')   lake(ilake)%qprod1   = xr(j,rindex(i))
          IF(str(i)=='qprod2    ')   lake(ilake)%qprod2   = xr(j,rindex(i))
          IF(str(i)=='maxqprod  ')   lake(ilake)%mqprod   = xr(j,rindex(i))
          IF(str(i)=='minflow   ')   lake(ilake)%minflow  = integer_convert_to_logical(-xi(j,iindex(i)))
          IF(str(i)=='datum1    ')   qprod1date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='datum2    ')   qprod2date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='rate      ')   lake(ilake)%rate     = xr(j,rindex(i))
          IF(str(i)=='exp       ')   lake(ilake)%exp      = xr(j,rindex(i))
          IF(str(i)=='qamp      ')   lake(ilake)%qamp     = xr(j,rindex(i))
          IF(str(i)=='qpha      ')   lake(ilake)%qpha     = xr(j,rindex(i))
          IF(str(i)=='limqprod  ')   lake(ilake)%limprod  = xr(j,rindex(i))
          IF(str(i)=='regvol    ')   regvolume(ilake)     = xr(j,rindex(i))
          IF(str(i)=='wamp      ')   wamp(ilake)          = xr(j,rindex(i))
        ENDDO
        lakeindex(lakeindex2(j))=ilake
        lakeindex3(j)=ilake
      ELSEIF(xi(j,iindex(ictype))==i_lbasin.OR.xi(j,iindex(ictype))==i_llbasin)THEN     !lake basin
        ibasin = ibasin + 1
        IF(xr(j,rindex(icdepth)).NE.missing_value) basin(lakeindex2(j))%lakedepth(2) = xr(j,rindex(icdepth))
        DO i = 1,nlakefb
          IF(lakeidindex(i,1)==xi(j,iindex(iclakeid))) EXIT
        ENDDO
        IF(i>nlakefb)THEN 
          WRITE(6,*) 'ERROR: lakebasin misses lake. lakeid',xi(j,iindex(iclakeid))
          WRITE(6,*) 'ERROR: can not find row with info on whole lake.'
          status = 1
          RETURN
        ENDIF
        lakebasin(ibasin)%ilk = i
        IF(xi(j,iindex(ictype))==i_llbasin) lastlakeindex(lakebasin(ibasin)%ilk) = j
        IF(xi(j,iindex(ictype))==i_llbasin) lakebasin(ibasin)%last = .TRUE.
        lakebasinindex(lakeindex2(j))=ibasin
        !These have all w0ref=0, since not set
      ELSEIF(xi(j,iindex(ictype))==i_ldout1)THEN     !lake with 2 outlet, outlet 1
        ilake = ilake + 1
        IF(icdepth>0)THEN
          IF(xr(j,rindex(icdepth)).NE.missing_value) basin(lakeindex2(j))%lakedepth(2) = xr(j,rindex(icdepth))
        ENDIF
        DO i = 1,ncols
          IF(str(i)=='w0ref     ')   lake(ilake)%w0ref    = xr(j,rindex(i))
          IF(str(i)=='deltaw0   ')   lake(ilake)%deltaw0  = xr(j,rindex(i))
          IF(str(i)=='qprod1    ')   lake(ilake)%qprod1   = xr(j,rindex(i))
          IF(str(i)=='qprod2    ')   lake(ilake)%qprod2   = xr(j,rindex(i))
          IF(str(i)=='maxqprod  ')   lake(ilake)%mqprod   = xr(j,rindex(i))
          IF(str(i)=='minflow   ')   lake(ilake)%minflow  = integer_convert_to_logical(-xi(j,iindex(i)))
          IF(str(i)=='datum1    ')   qprod1date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='datum2    ')   qprod2date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='rate      ')   lake(ilake)%rate     = xr(j,rindex(i))
          IF(str(i)=='exp       ')   lake(ilake)%exp      = xr(j,rindex(i))
          IF(str(i)=='qamp      ')   lake(ilake)%qamp     = xr(j,rindex(i))
          IF(str(i)=='qpha      ')   lake(ilake)%qpha     = xr(j,rindex(i))
          IF(str(i)=='limqprod  ')   lake(ilake)%limprod  = xr(j,rindex(i))
          IF(str(i)=='regvol    ')   regvolume(ilake)     = xr(j,rindex(i))
          IF(str(i)=='wamp      ')   wamp(ilake)          = xr(j,rindex(i))
        ENDDO
        lakeindex(lakeindex2(j))=ilake
        lakeindex3(j)=ilake
      ELSEIF(xi(j,iindex(ictype))==i_ldout2)THEN     !lake with 2 outlet, outlet 2
        ilake = ilake + 1
        DO i = 1,ncols
          !no w0ref or lake_depth
          IF(str(i)=='deltaw0   ')   lake(ilake)%deltaw0  = xr(j,rindex(i))
          IF(str(i)=='qprod1    ')   lake(ilake)%qprod1   = xr(j,rindex(i))
          IF(str(i)=='qprod2    ')   lake(ilake)%qprod2   = xr(j,rindex(i))
          IF(str(i)=='maxqprod  ')   lake(ilake)%mqprod   = xr(j,rindex(i))
          IF(str(i)=='minflow   ')   lake(ilake)%minflow  = integer_convert_to_logical(-xi(j,iindex(i)))
          IF(str(i)=='datum1    ')   qprod1date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='datum2    ')   qprod2date(ilake)    = xi(j,iindex(i))
          IF(str(i)=='rate      ')   lake(ilake)%rate     = xr(j,rindex(i))
          IF(str(i)=='exp       ')   lake(ilake)%exp      = xr(j,rindex(i))
          IF(str(i)=='qamp      ')   lake(ilake)%qamp     = xr(j,rindex(i))
          IF(str(i)=='qpha      ')   lake(ilake)%qpha     = xr(j,rindex(i))
          IF(str(i)=='regvol    ')   regvolume(ilake)     = xr(j,rindex(i))
          IF(str(i)=='wamp      ')   wamp(ilake)          = xr(j,rindex(i))
          IF(str(i)=='w0ref     ')   lake(ilake)%w0ref    = xr(j,rindex(i))  !w0 relativt w0ref i out1
          IF(str(i)=='limqprod  ')   lake(ilake)%limprod  = xr(j,rindex(i))
        ENDDO
        lakeout2index(lakeindex2(j))=ilake
        lakeindex3(j)=ilake
      ENDIF
    ENDDO
    IF(ilake<nlakedim) WRITE(6,*) ilake,nlakedim
    IF(ibasin>0) nbasinlakes = lakebasin(ibasin)%ilk

    !Calculate other lake variables; wmin, datum, wampcoeff
    lake(:)%wmin = missing_value
    lake(:)%wampcoeff = missing_value
    DO j = 1,nrows
      IF(lakeindex3(j)>0)THEN
        IF(regvolume(lakeindex3(j))>0.)THEN
          IF(xi(j,iindex(ictype))==i_simple)THEN
            lakearea = basin(lakeindex2(j))%area * classbasin(lakeindex2(j),slc_olake)%part
          ELSEIF(xi(j,iindex(ictype))==i_ldout1)THEN
            lakearea = basin(lakeindex2(j))%area * classbasin(lakeindex2(j),slc_olake)%part
          ELSEIF(xi(j,iindex(ictype))==i_ldout2)THEN
            lakearea = basin(lakeindex2(j))%area * classbasin(lakeindex2(j),slc_olake)%part
          ELSEIF(xi(j,iindex(ictype))==i_combil)THEN
            lakearea = lake(lakeindex3(j))%area
          ENDIF
          IF(lakearea>0)THEN
            lake(lakeindex3(j))%wmin = 0. - regvolume(lakeindex3(j)) * 1000000. / lakearea
            IF(wamp(lakeindex3(j))/=missing_value) &
              lake(lakeindex3(j))%wampcoeff = wamp(lakeindex3(j))/(-1.*lake(lakeindex3(j))%wmin)
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    DO ilake = 1,nlakedim 
      MonthDayDate = qprod1date(ilake)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      lake(ilake)%datum1 = dayNumber
      MonthDayDate = qprod2date(ilake)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      lake(ilake)%datum2 = dayNumber
    ENDDO

    !Check lake area against slc-area
    IF(checkdata(1,1))THEN
      fex = .FALSE.
      DO j=1,nrows
        lakearea = 0
        parea = 0.
        IF(lakeindex2(j)>0 .AND. lakeindex3(j)>0)THEN
          lakearea = basin(lakeindex2(j))%area * classbasin(lakeindex2(j),slc_olake)%part
          IF(lakearea>0)THEN
            parea = ABS(lake(lakeindex3(j))%area-lakearea)/lakearea
          ELSEIF(lake(lakeindex3(j))%area>0)THEN
            parea = 0.
          ENDIF
        ENDIF
        IF(parea > 0.001)THEN
          WRITE(6,'(a29,es14.7,a28)') ' INFO: lakearea in LakeData (',lake(lakeindex3(j))%area,') differ more than 0.1% from'
          WRITE(6,'(a29,es14.7,a1)')  ' INFO: lakearea in GeoData  (',lakearea,')'
          WRITE(6,'(a18,I8)') ' INFO: for subid =',basin(lakeindex2(j))%subid 
          fex = .TRUE.
        ENDIF
      ENDDO
      IF(fex) WRITE(6,*)
    ENDIF

    !Build table with lakedata parameter values
    CALL start_lakedata_table(maxcol,ncols,nrows,str,xi,xr,iindex,rindex,n,lakedataid)

    !Deallocate local variables
    IF(ALLOCATED(lakeindex2)) DEALLOCATE(lakeindex2)
    IF(ALLOCATED(lakeindex3)) DEALLOCATE(lakeindex3)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(regvolume)) DEALLOCATE(regvolume)
    IF(ALLOCATED(qprod1date)) DEALLOCATE(qprod1date) 
    IF(ALLOCATED(qprod2date)) DEALLOCATE(qprod2date)
    IF(ALLOCATED(lakeidindex)) DEALLOCATE(lakeidindex)
    IF(ALLOCATED(lastlakeindex)) DEALLOCATE(lastlakeindex)

    WRITE(6,*) 'File read: ', TRIM(filename)

  END SUBROUTINE load_lakedata

  !>Collects the information about dams from DamData.txt
  !>
  !>\b Consequences Module modvar variables dam, damindex may be allocated and set. 
  !>Module modvar variable basin may be changed (lake_depth).
  !------------------------------------------------------------------
  SUBROUTINE load_damdata(funit,dir,n,status)

    USE WORLDVAR, ONLY : maxcharpath,  &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &    !OUT
                         classbasin,   &
                         slc_olake,    &
                         missing_value,  &
                         dam,         &    !OUT
                         damindex        !INOUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit         !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: n             !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: status        !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i,idam
    INTEGER idtype                    !Column with subid
    INTEGER nrows                     !Number of dams in file
    INTEGER ncols                     !Number of columns in file
    INTEGER maxcol
    INTEGER MonthDayDate              !Non-formated date from DamData.txt
    INTEGER dayNumber                 !Corresponding date in day number format
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)   !Content string
    INTEGER,ALLOCATABLE :: code(:)            !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)          !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)          !Index for column integer variables
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    REAL damarea                            ! area of dam reservoir (same as lakearea)
    INTEGER,ALLOCATABLE :: qprod1date(:)
    INTEGER,ALLOCATABLE :: qprod2date(:)
    INTEGER, ALLOCATABLE :: damindex2(:)   ! row where  
    REAL,ALLOCATABLE :: wamp(:)               !Regulation amplitude (real) (m)
    REAL,ALLOCATABLE :: qinflow(:,:)          !Monthly natural flow (m3/s)

    status = 0
    ALLOCATE(damindex(n))
    damindex = 0
    filename = TRIM(dir)//'DamData.txt'

    !Check if file exist and if outlet lake-class exist
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(slc_olake==0)THEN
      WRITE(6,*) 'No class for olake to be used with DamData.txt'
      status = 1
      RETURN
    ENDIF

    !Count number of columns and rows, and allocate variables for reading
    CALL count_data_cols(funit,filename,0,maxcol,status)
    IF(status/=0)THEN
      WRITE(6,*) 'ERROR: counting columns in DamData'
      RETURN
    ENDIF
    CALL count_data_rows(funit,filename,1,nrows,status)
    IF(status/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No dam in DamData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(maxcol))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(maxcol))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(maxcol))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(maxcol))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))
    IF(.NOT.ALLOCATED(qprod1date)) ALLOCATE(qprod1date(nrows)) 
    IF(.NOT.ALLOCATED(qprod2date)) ALLOCATE(qprod2date(nrows))
    IF(.NOT.ALLOCATED(wamp)) ALLOCATE(wamp(nrows))
    IF(.NOT.ALLOCATED(qinflow)) ALLOCATE(qinflow(12,nrows))
    qprod1date = 0; qprod2date = 0  
    wamp = missing_value
    qinflow = missing_value

    !Open DamData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old', ACTION='read')
    CALL read_column_headings(funit,maxcol,str,ncols,status)
    IF(status.NE.0)THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    idtype = 0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        idtype = i                                        !! idtype is column where SUBID is
      ENDIF
      
      IF(str(i)(1:10)=='lake_depth') code(i) = i_real
      IF(str(i)(1:10)=='purpose   ') code(i) = i_intg   !Purpose of dam (1=irrigation,2=water supply, 3=Flood control, 4 =hydroelectricity)
      IF(str(i)(1:10)=='regvol    ') code(i) = i_real   !dam_regvol
      IF(str(i)(1:10)=='datum1    ') code(i) = i_intg   !date for two date based prod rate
      IF(str(i)(1:10)=='datum2    ') code(i) = i_intg   !date for two date based prod rate
      IF(str(i)(1:10)=='rate      ') code(i) = i_real   !dam_rate
      IF(str(i)(1:10)=='exp       ') code(i) = i_real   !dam_exp
      IF(str(i)(1:10)=='w0ref     ') code(i) = i_real   !dam_w0
      IF(str(i)(1:10)=='wamp      ') code(i) = i_real   !wamp (m)
      IF(str(i)(1:10)=='limqprod  ') code(i) = i_real   !dam_limprod
      IF(str(i)(1:10)=='qprod1    ') code(i) = i_real   !two date-based production rates, qprod1 and qprod2, now replacing the unique rate qprod
      IF(str(i)(1:10)=='qprod2    ') code(i) = i_real
      IF(str(i)(1:10)=='qamp      ') code(i) = i_real   !dam_qamp
      IF(str(i)(1:10)=='qpha      ') code(i) = i_real   !dam_qpha
      IF(str(i)(1:10)=='qinfjan   ') code(i) = i_real   !dam_qinfjan !Mean monthly inflow to dam (simulated using natural run)
      IF(str(i)(1:10)=='qinffeb   ') code(i) = i_real   !dam_qinffeb
      IF(str(i)(1:10)=='qinfmar   ') code(i) = i_real   !dam_qinfmar
      IF(str(i)(1:10)=='qinfapr   ') code(i) = i_real   !dam_qinfapr
      IF(str(i)(1:10)=='qinfmay   ') code(i) = i_real   !dam_qinfmay
      IF(str(i)(1:10)=='qinfjun   ') code(i) = i_real   !dam_qinfjun
      IF(str(i)(1:10)=='qinfjul   ') code(i) = i_real   !dam_qinfjul
      IF(str(i)(1:10)=='qinfaug   ') code(i) = i_real   !dam_qinfaug
      IF(str(i)(1:10)=='qinfsep   ') code(i) = i_real   !dam_qinfsep
      IF(str(i)(1:10)=='qinfoct   ') code(i) = i_real   !dam_qinfoct
      IF(str(i)(1:10)=='qinfnov   ') code(i) = i_real   !dam_qinfnov
      IF(str(i)(1:10)=='qinfdec   ') code(i) = i_real   !dam_qinfdec
      IF(str(i)(1:10)=='snowfrac  ') code(i) = i_real   !dam_snowfrac 
    ENDDO

    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    !Allocate variables for dam information
    IF(.NOT.ALLOCATED(dam)) ALLOCATE(dam(nrows))
    ALLOCATE(damindex2(nrows))

    !Find subbasin-dam coupling: damindex(isb)
    IF(idtype==0)THEN
      WRITE(6,*) 'ERROR: subid not found in DamData. ',TRIM(filename)
      status = 1
      RETURN
    ENDIF
    DO idam = 1,nrows
      DO i = 1,n
        IF(basin(i)%subid==xi(idam,iindex(idtype)))THEN
          damindex(i)=idam
          damindex2(idam)=i    !! find row in geodata for dam row
          EXIT
        ENDIF
      ENDDO
      IF(i>n)THEN
        WRITE(6,*) 'ERROR: A dam in DamData is missing in GeoData. subid: ',xi(idam,iindex(idtype))
        status = 1
        RETURN
      ENDIF  
    ENDDO

    !Save dam information
    DO idam = 1,nrows
      DO i = 1,ncols
        IF(str(i)=='regvol    ')   dam(idam)%regvol    = xr(idam,rindex(i))
        IF(str(i)=='rate      ')   dam(idam)%rate      = xr(idam,rindex(i))
        IF(str(i)=='exp       ')   dam(idam)%exp       = xr(idam,rindex(i))
        IF(str(i)=='w0ref     ')   dam(idam)%w0ref     = xr(idam,rindex(i))
        IF(str(i)=='wamp      ')   wamp(idam)          = xr(idam,rindex(i))
        IF(str(i)=='qprod1    ')   dam(idam)%qprod1    = xr(idam,rindex(i))
        IF(str(i)=='qprod2    ')   dam(idam)%qprod2    = xr(idam,rindex(i))
        IF(str(i)=='qamp      ')   dam(idam)%qamp      = xr(idam,rindex(i))
        IF(str(i)=='qpha      ')   dam(idam)%qpha      = xr(idam,rindex(i))
        IF(str(i)=='limqprod  ')   dam(idam)%limprod   = xr(idam,rindex(i))
        IF(str(i)=='qinfjan   ')   qinflow(1,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinffeb   ')   qinflow(2,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfmar   ')   qinflow(3,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfapr   ')   qinflow(4,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfmay   ')   qinflow(5,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfjun   ')   qinflow(6,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfjul   ')   qinflow(7,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfaug   ')   qinflow(8,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfsep   ')   qinflow(9,idam)     = xr(idam,rindex(i))
        IF(str(i)=='qinfoct   ')   qinflow(10,idam)    = xr(idam,rindex(i))
        IF(str(i)=='qinfnov   ')   qinflow(11,idam)    = xr(idam,rindex(i))
        IF(str(i)=='qinfdec   ')   qinflow(12,idam)    = xr(idam,rindex(i))
        IF(str(i)=='snowfrac  ')   dam(idam)%snowfrac  = xr(idam,rindex(i))
        IF(str(i)=='lake_depth')THEN
          IF(xr(idam,rindex(i)).NE.missing_value)   basin(damindex2(idam))%lakedepth(2) = xr(idam,rindex(i))
        ENDIF
        !Integer input
        IF(str(i)=='datum1    ')   qprod1date(idam)    = xi(idam,iindex(i))
        IF(str(i)=='datum2    ')   qprod2date(idam)    = xi(idam,iindex(i))
        IF(str(i)=='purpose   ')   dam(idam)%purpose   = xi(idam,iindex(i))
      ENDDO
    ENDDO
   
    !Calculate other lake/dam variables; w0ref, wmin, datum, wampcoeff, qinflow
    dam(:)%wmin = missing_value
    dam(:)%wampcoeff = missing_value
    DO i = 1,n
      IF(damindex(i)>0)THEN
        damarea = basin(i)%area * classbasin(i,slc_olake)%part
        IF(damarea>0)THEN
          dam(damindex(i))%wmin = 0. - dam(damindex(i))%regvol * 1000000. / damarea
          IF(wamp(damindex(i))/=missing_value) &
              dam(damindex(i))%wampcoeff = wamp(damindex(i))/(-1.*dam(damindex(i))%wmin)
        ENDIF
      ENDIF
    ENDDO
    DO idam = 1,nrows 
      MonthDayDate = qprod1date(idam)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      dam(idam)%datum1 = dayNumber
      MonthDayDate = qprod2date(idam)
      CALL get_dayno_from_monthday(MonthDayDate, dayNumber)
      dam(idam)%datum2 = dayNumber
      dam(idam)%qinfmed = SUM(qinflow(:,idam))/12.
      dam(idam)%qinfmin = MINVAL(qinflow(:,idam))
      dam(idam)%qinfmax = MAXVAL(qinflow(:,idam))
    ENDDO
 
    !Deallocate local variables
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(qprod1date)) DEALLOCATE(qprod1date) 
    IF(ALLOCATED(qprod2date)) DEALLOCATE(qprod2date)
    IF(ALLOCATED(wamp)) DEALLOCATE(wamp)
    IF(ALLOCATED(damindex2)) DEALLOCATE(damindex2)

    WRITE(6,*) 'File read: ', TRIM(filename)

  END SUBROUTINE load_damdata

  !>Collects the information about flooding areas from FloodData.txt
  !>
  !>\b Consequences Module modvar variables flooding, floodindex may be allocated and set. 
  !------------------------------------------------------------------
  SUBROUTINE load_flooddata(funit,dir,n,status)

    USE WORLDVAR, ONLY : maxcharpath,  &
                         i_str,        &
                         i_intg,       &
                         i_real
    USE MODVAR,   ONLY : basin,        &
                         missing_value, &
                         conductflood,  &   !OUT
                         flooding,      &   !OUT
                         floodindex,    &   !OUT
                         modeloption,   &
                         p_floodplain

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit         !<File unit, temporary file open
    CHARACTER (LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN)  :: n             !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: status        !<Error status
    
    !Local variables
    LOGICAL fex
    INTEGER i,idam
    INTEGER maxcol
    INTEGER idtype                    !Column with subid
    INTEGER nrows                     !Number of dams in file
    INTEGER ncols                     !Number of columns in file
    INTEGER,ALLOCATABLE :: code(:)    !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)  !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)  !Index for column integer variables
    CHARACTER(LEN=maxcharpath) filename
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)   !Content string
    INTEGER,ALLOCATABLE :: xi(:,:)            !Integer data read from file
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file

    status = 0
    filename = TRIM(dir)//'FloodData.txt'
    conductflood = .FALSE.

    !Check if file exist and if floodmodel is set
    INQUIRE(FILE=filename,EXIST=fex)
    IF(.NOT.fex) RETURN
    IF(modeloption(p_floodplain)==0)THEN
      WRITE(6,*) 'WARNING: floodmodel not set in info.txt. Simulation will not simulate floodplains or use FloodData.txt'
      RETURN
    ENDIF

    IF(.NOT.ALLOCATED(floodindex)) ALLOCATE(floodindex(n))
    floodindex = 0

    !Count number of data rows and columns, and allocate variables for reading
    CALL count_data_cols(funit,filename,0,maxcol,status)
    IF(status/=0)THEN
      WRITE(6,*) 'ERROR: counting columns in FloodData.txt'
      RETURN
    ENDIF
    CALL count_data_rows(funit,filename,1,nrows,status)
    IF(status/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'WARNING: No flood areas in FloodData.txt ',TRIM(filename)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(flooding)) ALLOCATE(flooding(nrows))
    flooding%fpfmr = 0.
    flooding%fpfol = 0.
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(maxcol))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(maxcol))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(maxcol))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(maxcol))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,maxcol))

    !Open FloodData-file and read headings
    OPEN(UNIT = funit,FILE = filename, STATUS = 'old')
    CALL read_column_headings(funit,maxcol,str,ncols,status)
    IF(status.NE.0)THEN
      WRITE(6,*) 'ERROR: reading file: ',TRIM(filename)
      RETURN
    ENDIF

    !Find variable columns and set variable type
    code = i_str      !string, ignore
    idtype = 0
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        idtype = i                                        !! idtype is column where SUBID is
      ENDIF
      IF(str(i)(1:10)=='fpfol     ') code(i) = i_real
      IF(str(i)(1:10)=='fpfmr     ') code(i) = i_real
      IF(str(i)(1:10)=='floll     ') code(i) = i_real
      IF(str(i)(1:10)=='flolp     ') code(i) = i_real
      IF(str(i)(1:10)=='flmrr     ') code(i) = i_real
      IF(str(i)(1:10)=='flmrp     ') code(i) = i_real
      IF(str(i)(1:10)=='rclfp     ') code(i) = i_real
      IF(str(i)(1:10)=='rcrfp     ') code(i) = i_real
      IF(str(i)(1:10)=='rcfpl     ') code(i) = i_real
      IF(str(i)(1:10)=='rcfpr     ') code(i) = i_real
      IF(str(i)(1:10)=='fymol     ') code(i) = i_real
      IF(str(i)(1:10)=='fymmr     ') code(i) = i_real
    ENDDO

    !Read all data
    CALL read_basindata5(funit,filename,maxcol,nrows,ncols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)

    !Find subbasin-floodplain coupling: floodindex(isb)
    IF(idtype==0)THEN
      WRITE(6,*) 'ERROR: subid not found in FloodData. ',TRIM(filename)
      status = 1
      RETURN
    ENDIF
    DO idam = 1,nrows
      DO i = 1,n
        IF(basin(i)%subid==xi(idam,iindex(idtype)))THEN
          floodindex(i)=idam
          EXIT
        ENDIF
      ENDDO
      IF(i>n)THEN
        WRITE(6,*) 'Error: A floodplain in FloodData in a subbasin missing in GeoData. subid: ',xi(idam,iindex(idtype))
        status = 1
        RETURN
      ENDIF  
    ENDDO

    !Save floodplain information
    DO idam = 1,nrows
      DO i = 1,ncols
        IF(str(i)(1:10)=='fpfol     ') flooding(idam)%fpfol  = xr(idam,rindex(i))  !floodplain
        IF(str(i)(1:10)=='fpfmr     ') flooding(idam)%fpfmr  = xr(idam,rindex(i))
        IF(str(i)(1:10)=='floll     ') flooding(idam)%floll  = xr(idam,rindex(i))
        IF(str(i)(1:10)=='flolp     ') flooding(idam)%flolp  = xr(idam,rindex(i))
        IF(str(i)(1:10)=='flmrr     ') flooding(idam)%flmrr  = xr(idam,rindex(i))
        IF(str(i)(1:10)=='flmrp     ') flooding(idam)%flmrp  = xr(idam,rindex(i))
        IF(str(i)(1:10)=='rclfp     ') flooding(idam)%rcl2fp = xr(idam,rindex(i))
        IF(str(i)(1:10)=='rcrfp     ') flooding(idam)%rcr2fp = xr(idam,rindex(i))
        IF(str(i)(1:10)=='rcfpl     ') flooding(idam)%rcfp2l = xr(idam,rindex(i))
        IF(str(i)(1:10)=='rcfpr     ') flooding(idam)%rcfp2r = xr(idam,rindex(i))
        IF(str(i)(1:10)=='fymol     ') flooding(idam)%fymol  = xr(idam,rindex(i))
        IF(str(i)(1:10)=='fymmr     ') flooding(idam)%fymmr  = xr(idam,rindex(i))
      ENDDO
    ENDDO
   
    !Check if flooded area is simulated
    IF(SUM(flooding%fpfmr) + SUM(flooding%fpfol)>0.) conductflood = .TRUE.

    !Deallocate local variables
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

    WRITE(6,*) 'File read: ', TRIM(filename)

  END SUBROUTINE load_flooddata

  !------------------------------------------------------------------
  !>\brief Set lakedatapar values to those model variable values specified
  !>in LakeData.txt. 
  !!
  !>\b Consequences Module modvar variables lakedatapar, lakedataparindex may be allocated and set. 
  !------------------------------------------------------------------
  SUBROUTINE start_lakedata_table(maxcol,ncols,nrows,str,xi,xr,iindex,rindex,ns,lakedataid)

    USE MODVAR, ONLY : modparid,          &
                       lakedatapar,       &   !OUT
                       lakedataparindex,  &   !OUT
                       max_par,           &
                       m_ldpar,           &
                       nregions,      &
                       missing_value,     &
                       basin

    !Argument declaration
    INTEGER, INTENT(IN)  :: maxcol                    !<Maximum possible amount of columns on LakeData.txt
    INTEGER, INTENT(IN)  :: ncols                     !<Number of columns effectively used in LakeData.txt
    INTEGER, INTENT(IN)  :: nrows                     !<Number of rows effectively used in LakeData.txt
    CHARACTER(LEN=10), INTENT(IN) :: str(maxcol)      !<Character array with labels of ALL columns given in LakeData.txt (column headers)
    INTEGER, INTENT(IN)  :: xi(nrows,maxcol)          !<Integer values aquired by reading LakeData.txt
    REAL, INTENT(IN)     :: xr(nrows,maxcol)          !<Float values aquired by reading LakeData.txt
    INTEGER, INTENT(IN)  :: iindex(maxcol)            !<Index for column real variables
    INTEGER, INTENT(IN)  :: rindex(maxcol)            !<Index for column integer variables
    INTEGER, INTENT(IN)  :: ns                        !<Number of subbasins (base model)
    INTEGER, INTENT(IN)  :: lakedataid(ns)            !<lakedataid from GeoData.txt
    
    !Local variables
    INTEGER rowCounter, colCounter, modparCounter, basinCounter
    INTEGER nlakedatapar, tempVal1, nlakeregions

    nlakedatapar = 0          !Amount of lakedatapar model variables; set to 0, then find the largest one
    nlakeregions = nregions(6)
    DO modparCounter = 1,max_par
      IF(modparid(modparCounter)%deptype==m_ldpar)THEN
        nlakedatapar = MAX(nlakedatapar,modparid(modparCounter)%parno)
      ENDIF
    ENDDO

    IF(.NOT.ALLOCATED(lakedatapar)) ALLOCATE(lakedatapar(nlakeregions+nrows, nlakedatapar)) !This is more rows than needed
    lakedatapar = missing_value       !Initialize the table of parameter values with the missing value (-9999)
    IF(.NOT.ALLOCATED(lakedataparindex)) ALLOCATE(lakedataparindex(ns,2))
    lakedataparindex = 0              !Initialize the line index table with 0 (useful for 2nd colum, with the row number)

    !Parameter by parameter (= column by column), copy those parameter values given in LakeData.txt in the table "lakedatapar" (it's a modvar)
    !Values are put at the bottom of the table. Eventual parameter values coming from par.txt are appended later in the table top rows, in another subroutine.
    !The order in which lakes are given in LakeData.txt is respected linewise
    DO colCounter = 1,ncols

      !Loop through all model variables to look for label match with the currently considered parameter, labelled str(colCounter) in lakeData
      DO modparCounter = 1,max_par
        !When a matching variable shortname is found, determine if the model variable has the parameter value read from lakedata.txt (deptype = m_ldpar), or that from par.txt (deptype = others)
        IF(str(colCounter)==modparid(modparCounter)%shortname .AND. modparid(modparCounter)%deptype==m_ldpar)THEN
          !Parameter values coming from lakedata are stored in xr for float values; take them from there
          IF(rindex(colCounter).NE.0)THEN
            DO rowCounter = 1,nrows
              lakedatapar(nlakeregions+rowCounter, modparid(modparCounter)%parno) = xr(rowCounter, rindex(colCounter))
            ENDDO
          ENDIF
          EXIT
        ENDIF   !End of the check on label match and on the origin of the parameter value (LakeData.txt or par.txt)
      ENDDO   !End of model variable name loop

      !Find corresponding subbasin for every lake, and set pointer to that lake
      IF(str(colCounter)=='lakedataid')THEN
        DO rowCounter = 1,nrows
          tempVal1 = xi(rowCounter, iindex(colCounter))
          IF(tempVal1>0)THEN
            DO basinCounter = 1,ns
              IF(lakedataid(basinCounter) == tempVal1)THEN    
                lakedataparindex(basinCounter, 2) = nlakeregions+rowCounter
                EXIT
              ENDIF
            ENDDO
          ENDIF
        ENDDO
      ENDIF

    ENDDO   !End of lakeData column labels loop

    !Point ilake and non-LakeData lakes to lakeregion-row in lakedatapar
    DO basinCounter = 1,ns
      lakedataparindex(basinCounter,1) = basin(basinCounter)%parregion(6)
      IF(lakedataparindex(basinCounter,2) == 0)THEN
        lakedataparindex(basinCounter,2) = basin(basinCounter)%parregion(6)
      ENDIF
    ENDDO

  END SUBROUTINE start_lakedata_table

  !>Set lakedatapar values to those specified in par.txt.
  !------------------------------------------------------------------
  SUBROUTINE finish_lakedata_table(lakedatafile,ns)

    USE MODVAR, ONLY : modparid, &
                       lakedatapar, &
                       lakedataparindex, &
                       basin, &
                       genpar, &
                       regpar, &    
                       max_par, &
                       m_gpar, &
                       m_rpar, &
                       m_ldpar, &
                       nregions, &  
                       missing_value

    !Argument declaration
    LOGICAL, INTENT(IN) :: lakedatafile   !<status of LakeData.txt
    INTEGER, INTENT(IN) :: ns             !<number of subbasins (nsub_basemodel)
    
    !Local variables
    LOGICAL found         !corresponding parameter found
    INTEGER parCounter    !index for loop to find corresponding par
    INTEGER ldparCounter  !index for loop to find ldpar
    INTEGER lakeCounter   !lake/row in lakedata
    INTEGER isb           !subbasin index/number
    INTEGER dimlakepar    !size of lakedatapar
    INTEGER dimlakepar2   !second dimension of lakedatapar
    INTEGER nlakeregions
    INTEGER,ALLOCATABLE :: subbasinindex(:) !subbasin index for lakedata (rows)

    !Start of subroutine
    IF(lakedatafile)THEN

      !Set index table for finding subbasin index for lakedata (row)
      dimlakepar = SIZE(lakedatapar,1)
      nlakeregions = nregions(6)
      IF(.NOT.ALLOCATED(subbasinindex)) ALLOCATE(subbasinindex(dimlakepar))
      subbasinindex = 0
      DO lakeCounter = nlakeregions+1,dimlakepar
        DO isb = 1,ns
          IF(lakedataparindex(isb,2)==lakeCounter)THEN
            EXIT
          ENDIF
        ENDDO
        IF(isb<=ns) subbasinindex(lakeCounter) = isb
      ENDDO

      DO ldparCounter = 1,max_par
        !For every lakedata-parameter find the corresponding parameter in par.txt
        IF(modparid(ldparCounter)%deptype==m_ldpar)THEN

          !Loop through all model parameters to look for label match with the 
          !currently considered parameter, labelled str(colCounter) in LakeData.txt
          found = .FALSE.
          DO parCounter = 1,max_par
            IF(modparid(parCounter)%deptype==m_gpar .OR. modparid(parCounter)%deptype==m_rpar)THEN
              IF(modparid(parCounter)%shortname==modparid(ldparCounter)%shortname)THEN  !Check label match with model parameter of ldpar-type
                !Corresponding parameter found!
                found = .TRUE.
                !Set general/lakeregion parameter values in lakedatapar
                IF(modparid(parCounter)%deptype==m_gpar)THEN 
                  lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = genpar(modparid(parCounter)%parno)
                ELSEIF(modparid(parCounter)%deptype==m_rpar)THEN
                  lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = regpar(modparid(parCounter)%parno,1:nlakeregions)
                ENDIF

                !Set missing parameters from LakeData file to general/lakeregion parameter value
                DO lakeCounter = nlakeregions+1,dimlakepar
                  IF(lakedatapar(lakeCounter, modparid(ldparCounter)%parno)==missing_value)THEN
                    IF(modparid(parCounter)%deptype==m_gpar)THEN
                      lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = genpar(modparid(parCounter)%parno)
                    ENDIF
                    IF(modparid(parCounter)%deptype==m_rpar)THEN
                      IF(subbasinindex(lakeCounter)>0)THEN
!                        lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = lregpar(modparid(parCounter)%parno,basin(subbasinindex(lakeCounter))%lakeregion(2))
                        lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = regpar(modparid(parCounter)%parno,basin(subbasinindex(lakeCounter))%parregion(6))
                      ENDIF
                    ENDIF
                  ENDIF
                ENDDO

                !Next ldpar
                EXIT
              ENDIF !shortname is equal
            ENDIF !gpar or lrpar
          ENDDO !parCounter

          IF(.NOT.found)THEN  !This variable does not have a corresponing par.txt parameter
            !Set non-LakeData set lake to zero
            lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = 0.
            DO lakeCounter = nlakeregions+1,dimlakepar
              IF(lakedatapar(lakeCounter, modparid(ldparCounter)%parno)==missing_value)THEN
                lakedatapar(lakeCounter, modparid(ldparCounter)%parno) = 0.
              ENDIF
            ENDDO
          ENDIF

        ENDIF !ldpar
      ENDDO !ldparCounter

      IF(ALLOCATED(subbasinindex)) DEALLOCATE(subbasinindex)

    ELSE
      !No LakeData.txt: set lakedatapar with only the general/lakeregion parameters 
      nlakeregions = nregions(6)
      dimlakepar = nlakeregions
      dimlakepar2 = 0          !Amount of lakedatapar model variables; set to 0, then find the largest one
      DO parCounter = 1,max_par
        IF(modparid(parCounter)%deptype==m_ldpar)THEN
          dimlakepar2 = MAX(dimlakepar2,modparid(parCounter)%parno)
        ENDIF
      ENDDO
      IF(.NOT.ALLOCATED(lakedatapar)) ALLOCATE(lakedatapar(nlakeregions,dimlakepar2))
      IF(.NOT.ALLOCATED(lakedataparindex)) ALLOCATE(lakedataparindex(ns,2))

      !Point ilake and non-LakeData olakes to lakeregion-row in lakedatapar
!      lakedataparindex(:,1) = basin(:)%lakeregion(2)
      lakedataparindex(:,1) = basin(:)%parregion(6)
      lakedataparindex(:,2) = lakedataparindex(:,1)

      DO ldparCounter = 1,max_par
          !For every lakedata-parameter find the corresponding parameter in par.txt
        IF(modparid(ldparCounter)%deptype==m_ldpar)THEN

            !Loop through all model parameters to look for label match with the currently considered parameter, labelled str(colCounter) in LakeData.txt
          DO parCounter = 1,max_par
            IF(modparid(parCounter)%deptype==m_gpar .OR. modparid(parCounter)%deptype==m_rpar)THEN
              IF(modparid(parCounter)%shortname==modparid(ldparCounter)%shortname)THEN  !Check label match with model parameter of ldpar-type
                !Corresponding parameter found!

                !Set general/lakeregion parameter values in lakedatapar
                IF(modparid(parCounter)%deptype==m_gpar)THEN 
                  lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = genpar(modparid(parCounter)%parno)
                ELSEIF(modparid(parCounter)%deptype==m_rpar)THEN
                  lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = regpar(modparid(parCounter)%parno,1:nlakeregions)
!                  lakedatapar(1:nlakeregions, modparid(ldparCounter)%parno) = lregpar(modparid(parCounter)%parno,1:nlakeregions)
                ENDIF

                !Next ldpar
                EXIT
              ENDIF !shortname is equal
            ENDIF !gpar or lrpar
          ENDDO !parCounter
        ENDIF !ldpar
      ENDDO !ldparCounter
    ENDIF

  END SUBROUTINE finish_lakedata_table

  !>\brief Gets the information about subbasins from GeoData.
  !!Reads the matrix of basin data values from the file vith mcols columns
  !!
  !>\b Consequences Module variables will be allocated and set
  !------------------------------------------------------------------------------
  SUBROUTINE read_and_calc_basindata(funit,infile,n,lakedataid,maxcol,mcols) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real
    USE MODVAR, ONLY : timesteps_per_day, &
                       i_in,i_on,i_sp,i_pp,  &
                       max_pstype, &
                       basin,        &    !OUT
                       classbasin,   &    !OUT
                       pathsubid,    &    !OUT
                       load,         &    !OUT
                       wetland,      &    !OUT
                       petmodel,     &    !OUT
                       nregiondivisions, &
                       nregions         !OUT  !Number of parameter regions

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile        !<Name of characteristics file to be read
    INTEGER, INTENT(IN)  :: n                      !<Number of subbasins
    INTEGER, INTENT(OUT) :: lakedataid(n)          !<lakedataid from GeoData
    INTEGER, INTENT(IN)  :: maxcol                 !<Maximum number of data columns
    INTEGER, INTENT(OUT) :: mcols                  !<Actual number of columns
    
    !Local variables
    INTEGER i,k
    INTEGER idslc
    INTEGER xi(n,maxcol)               !Integer data read from file
    INTEGER code(maxcol)               !Code for column variable
    INTEGER rindex(maxcol)             !Index for column real variables
    INTEGER iindex(maxcol)             !Index for column integer variables
    INTEGER status                     !Error status
    REAL    help                       !Total load
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL    xr(n,maxcol)               !Real data read from file
    REAL    pointsource(n,5,max_pstype) !Read pointsource info (subbasins,column,pstype)
    REAL    localsource(n,4)           !Read rural source info (subbasins,column)
    CHARACTER(LEN=10) str(maxcol)      !Content string

    !Start of subroutine
    status = 0
    onetspday = 1./REAL(timesteps_per_day)
    pointsource = 0.
    
    OPEN(UNIT = funit,FILE = infile, STATUS = 'old', ACTION='read')     !Open GeoData-file

    !Reads the column headings from file
    CALL read_column_headings(funit,maxcol,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
      IF(str(i)(1:10)=='area      ') code(i) = i_real
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='xcoord    ') code(i) = i_real
      IF(str(i)(1:10)=='ycoord    ') code(i) = i_real
      IF(str(i)(1:10)=='longitude ') code(i) = i_real
      IF(str(i)(1:10)=='latitude  ') code(i) = i_real
      IF(str(i)(1:10)=='elev_mean ') code(i) = i_real
      IF(str(i)(1:10)=='elev_std  ') code(i) = i_real
      IF(str(i)(1:10)=='slope_mean') code(i) = i_real
      IF(str(i)(1:10)=='slope_std ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_in    ') code(i) = i_real
      IF(str(i)(1:10)=='ps1_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_in    ') code(i) = i_real
      IF(str(i)(1:10)=='ps2_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_in    ') code(i) = i_real
      IF(str(i)(1:10)=='ps3_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='loc_tp    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_sp    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_tn    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_in    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_vol   ') code(i) = i_real
      IF(str(i)(1:10)=='region    ') code(i) = i_intg
      IF(str(i)(1:10)=='petmodel  ') code(i) = i_intg
      IF(str(i)(1:10)=='lakeregion') code(i) = i_intg
      IF(str(i)(1:10)=='parreg_4  ') code(i) = i_intg
      IF(str(i)(1:10)=='parreg_5  ') code(i) = i_intg
      IF(str(i)(1:10)=='ilregion  ') code(i) = i_intg
      IF(str(i)(1:10)=='olregion  ') code(i) = i_intg
      IF(str(i)(1:10)=='lake_depth') code(i) = i_real
      IF(str(i)(1:10)=='lake_whigh'.OR.str(i)(1:10)=='lake_qavg '.OR.   &
         str(i)(1:10)=='lake_qamp '.OR.str(i)(1:10)=='lake_rate '.OR.   &
         str(i)(1:10)=='lake_qhigh'.OR.str(i)(1:10)=='lake_exp  '.OR.   &
         str(i)(1:10)=='mq        '.OR.str(i)(1:10)=='regvol    '.OR.   &
         str(i)(1:10)=='lake_qpha '.OR.str(i)(1:10)=='lake_wref ')THEN
        WRITE(6,*) 'WARNING: specific lake parameter no longer read from GeoData',str
        WRITE(6,*)
      ENDIF
      IF(str(i)(1:10)=='icatch    ') code(i) = i_real
      IF(str(i)(1:10)=='mainfl    ') THEN
        WRITE(6,*) 'WARNING: mainfl in GeoData is an column name no longer used' 
        WRITE(6,*)
      ENDIF
      IF(str(i)(1:10)=='branch    ') WRITE(6,*) 'WARNING: branch in GeoData is an column name no longer used'
      IF(str(i)(1:10)=='grwflow1  ') WRITE(6,*) 'WARNING: grwflow1 in GeoData is an column name no longer used'
      IF(str(i)(1:10)=='branchdown') WRITE(6,*) 'WARNING: branchdown in GeoData is an column name no longer used'
      IF(str(i)(1:10)=='maindown  ') code(i) = i_intg
      IF(str(i)(1:10)=='grwdown   ') code(i) = i_intg
      IF(str(i)(1:10)=='grwolake  ') code(i) = i_real
      IF(str(i)(1:3)=='slc')         code(i) = i_real
      IF(str(i)(1:5)=='dhslc')       code(i) = i_real
      IF(str(i)(1:3)=='scr')         code(i) = i_real
      IF(str(i)(1:10)=='parreg    ') code(i) = i_intg
      IF(str(i)(1:10)=='wqparreg  ') code(i) = i_intg
      IF(str(i)(1:10)=='rivlen    ') code(i) = i_real
      IF(str(i)(1:10)=='loc_rivlen') code(i) = i_real
      IF(str(i)(1:10)=='wetdep_n  ') code(i) = i_real
      IF(str(i)(1:10)=='drydep_n1 ') code(i) = i_real
      IF(str(i)(1:10)=='drydep_n2 ') code(i) = i_real
      IF(str(i)(1:10)=='drydep_n3 ') code(i) = i_real
      IF(str(i)(1:10)=='lrwet_area') code(i) = i_real
      IF(str(i)(1:10)=='mrwet_area') code(i) = i_real
      IF(str(i)(1:10)=='lrwet_dep ') code(i) = i_real
      IF(str(i)(1:10)=='mrwet_dep ') code(i) = i_real
      IF(str(i)(1:10)=='lrwet_part') code(i) = i_real
      IF(str(i)(1:10)=='mrwet_part') code(i) = i_real
      IF(str(i)(1:10)=='buffer    ') code(i) = i_real
      IF(str(i)(1:10)=='eroindex  ') code(i) = i_real
      IF(str(i)(1:10)=='close_w   ') code(i) = i_real
      IF(str(i)(1:10)=='lakedataid') code(i) = i_intg   !id i LakeData
      IF(str(i)(1:10)=='pobsid    ') WRITE(6,*) 'WARNING: pobsid etc in GeoData in no longer used, remove confusing columns'
      IF(str(i)(1:10)=='tobsid    ') WRITE(6,*) 'WARNING: tobsid etc in GeoData in no longer used, remove confusing columns'
    ENDDO

    !Read all data
    CALL read_basindata5(funit,infile,maxcol,n,mcols,code,rindex,iindex,xi,xr)
    CLOSE(UNIT=funit)
    WRITE(6,*) 'File read: ', TRIM(infile)

    DO i = 1,mcols
      IF(str(i)(1:10)=='area      ')   basin(1:n)%area      = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='subid     ')   basin(1:n)%subid     = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='xcoord    ')   basin(1:n)%xcoord    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ycoord    ')   basin(1:n)%ycoord    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='longitude ')   basin(1:n)%longitude = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='latitude  ')   basin(1:n)%latitude  = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='elev_mean ')   basin(1:n)%elev      = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='elev_std  ')   basin(1:n)%selev     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='slope_mean')   basin(1:n)%slope     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='slope_std ')   basin(1:n)%sslope    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_tp    ')   pointsource(1:n,1,1) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_sp    ')   pointsource(1:n,2,1) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_tn    ')   pointsource(1:n,3,1) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_in    ')   pointsource(1:n,4,1) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps1_vol   ')   pointsource(1:n,5,1) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_tp    ')   pointsource(1:n,1,2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_sp    ')   pointsource(1:n,2,2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_tn    ')   pointsource(1:n,3,2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_in    ')   pointsource(1:n,4,2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps2_vol   ')   pointsource(1:n,5,2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_tp    ')   pointsource(1:n,1,3) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_sp    ')   pointsource(1:n,2,3) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_tn    ')   pointsource(1:n,3,3) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_in    ')   pointsource(1:n,4,3) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='ps3_vol   ')   pointsource(1:n,5,3) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_tp    ')   localsource(1:n,1) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_sp    ')   localsource(1:n,2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_tn    ')   localsource(1:n,3) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_in    ')   localsource(1:n,4) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='loc_vol   ')  load(1:n)%volloc    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='region    ')  basin(1:n)%region   = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='petmodel  ')  petmodel(1:n) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='lakeregion')  basin(1:n)%parregion(6) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='ilregion  ')  basin(1:n)%parregion(4) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='olregion  ')  basin(1:n)%parregion(5) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='parreg_4  ')  basin(1:n)%parregion(4) = xi(1:n,iindex(i)) !Alternative
      IF(str(i)(1:10)=='parreg_5  ')  basin(1:n)%parregion(5) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='lake_depth')  basin(1:n)%lakedepth(2) = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='icatch    ')  basin(1:n)%ilakecatch = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='maindown  ')  pathsubid(1:n)%main   = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='grwdown   ')  pathsubid(1:n)%grw1   = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='grwolake  ')  pathsubid(1:n)%grwtolake = xr(1:n,rindex(i))
      IF(str(i)(1:3)=='slc') THEN
        idslc = 0
        IF(ICHAR(str(i)(5:5))>=49 .AND. ICHAR(str(i)(5:5))<=57)THEN
          idslc = ICHAR(str(i)(5:5)) - 48
        ENDIF
        IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(6:6)) - 48
        ENDIF
        classbasin(1:n,idslc)%part=xr(1:n,rindex(i))
      ENDIF
      IF(str(i)(1:5)=='dhslc') THEN
        idslc = 0
        IF(ICHAR(str(i)(7:7))>=49 .AND. ICHAR(str(i)(7:7))<=57)THEN
          idslc = ICHAR(str(i)(7:7)) - 48
        ENDIF
        IF(ICHAR(str(i)(8:8))>=48 .AND. ICHAR(str(i)(8:8))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(8:8)) - 48
        ENDIF
        classbasin(1:n,idslc)%deltah = xr(1:n,rindex(i))
      ENDIF
      IF(str(i)(1:3)=='scr') THEN    !area part with secondary crop
        idslc = 0
        IF(ICHAR(str(i)(5:5))>=49 .AND. ICHAR(str(i)(5:5))<=57)THEN
          idslc = ICHAR(str(i)(5:5)) - 48
        ENDIF
        IF(ICHAR(str(i)(6:6))>=48 .AND. ICHAR(str(i)(6:6))<=57)THEN
          IF(idslc>0) idslc = idslc * 10
          idslc = idslc + ICHAR(str(i)(6:6)) - 48
        ENDIF
        classbasin(1:n,idslc)%part2cr = xr(1:n,rindex(i))
      ENDIF
      IF(str(i)(1:10)=='parreg    ')  basin(1:n)%parregion(1) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='wqparreg  ')  basin(1:n)%parregion(3) = xi(1:n,iindex(i))
      IF(str(i)(1:10)=='loc_rivlen')  basin(1:n)%rivlen(1)    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='rivlen    ')  basin(1:n)%rivlen(2)    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='wetdep_n  ')  load(1:n)%inwetdep      = xr(1:n,rindex(i))*1.E-3    !ug/L -> mg/L
      IF(str(i)(1:10)=='drydep_n1 ')  load(1:n)%indrydep(1)   = xr(1:n,rindex(i))*onetspday  !kg/km2/d->kg/km2/ts
      IF(str(i)(1:10)=='drydep_n2 ')  load(1:n)%indrydep(2)   = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:10)=='drydep_n3 ')  load(1:n)%indrydep(3)   = xr(1:n,rindex(i))*onetspday
      IF(str(i)(1:10)=='lrwet_area')  wetland(1:n,1)%area     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='mrwet_area')  wetland(1:n,2)%area     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='lrwet_dep ')  wetland(1:n,1)%depth    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='mrwet_dep ')  wetland(1:n,2)%depth    = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='lrwet_part')  wetland(1:n,1)%part     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='mrwet_part')  wetland(1:n,2)%part     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='buffer    ')  basin(1:n)%buffer       = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='eroindex  ')  basin(1:n)%eroindex     = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='close_w   ')  basin(1:n)%closewater   = xr(1:n,rindex(i))
      IF(str(i)(1:10)=='lakedataid')  lakedataid(1:n)         = xi(1:n,iindex(i))
    ENDDO

    !Check groundwater flow direction      
    IF(SUM(pathsubid(1:n)%grw1)==0) pathsubid(1:n)%grw1 = pathsubid(1:n)%main
    !Calculate number of parameter regions
    DO i=1,nregiondivisions
      nregions(i) = MAXVAL(basin(1:n)%parregion(i))
    ENDDO
    !Calculate point source loads (PS in GeoData) or zeroing the loads.
    DO i = 1,n
      DO k = 1,max_pstype
        load(i)%psvol(k) = pointsource(i,5,k)/86400.  !m3/d->m3/s
      ENDDO
      IF(i_in>0)THEN
        DO k = 1,max_pstype
          help = pointsource(i,3,k)*pointsource(i,5,k)/timesteps_per_day*1.E-3    !kg TN/timestep
          load(i)%psload(k,i_in) = help*pointsource(i,4,k)
          load(i)%psload(k,i_on) = help*(1.-pointsource(i,4,k))
        ENDDO
      ENDIF
      IF(i_sp>0)THEN
        DO k = 1,max_pstype
          help = pointsource(i,1,k)*pointsource(i,5,k)/timesteps_per_day*1.E-3    !kg TP/timestep
          load(i)%psload(k,i_sp) = help*pointsource(i,2,k)
          load(i)%psload(k,i_pp) = help*(1.-pointsource(i,2,k))
        ENDDO
      ENDIF
    ENDDO
    !Calculate rural household loads (locload in GeoData) or zeroing the loads.
    DO i = 1,n
      load(i)%locconc = 0.  !default/initialise
      IF(i_in>0)THEN
        load(i)%locconc(i_in) = localsource(i,3)*localsource(i,4)
        load(i)%locconc(i_on) = localsource(i,3)*(1. - localsource(i,4))
      ENDIF
      IF(i_sp>0)THEN
        load(i)%locconc(i_sp) = localsource(i,1)*localsource(i,2)
        load(i)%locconc(i_pp) = localsource(i,1)*(1. - localsource(i,2))
      ENDIF
    ENDDO


  END SUBROUTINE read_and_calc_basindata

  !>\brief Gets the information about nutrient point sources from PointSourceData.txt
  !>Recalculate into form suitable for model. Set permanent loads.
  !!
  !>\b Consequences Module modvar variable load is set. Module worldvar variable 
  !>psdates may be allocated and set.
  !>
  !> \b Reference ModelDescription Chapter 5 (Point sources)
  !------------------------------------------------------------------------------
  SUBROUTINE load_pointsourcedata(dir,infile,n,status) 

    USE WORLDVAR, ONLY : steplen, &
                         psdates    !OUT
    USE MODVAR, ONLY : timesteps_per_day, &
                       i_in,i_on,i_sp,i_pp, &
                       i_t1,i_t2, &
                       basin,       &
                       load             !OUT
    USE READWRITE_ROUTINES, ONLY : find_reorder_index
 
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of pointsource file to be read (PointSourceData.txt)
    INTEGER, INTENT(IN)  :: n                    !<Number of subbasins
    INTEGER, INTENT(OUT) :: status               !<Error status
    
    !Local constants
    INTEGER, PARAMETER :: dim = 50
    
    !Local variables
    LOGICAL fexist                     !Status of file exist
    INTEGER i,j,k
    INTEGER nrows
    INTEGER maxdates
    INTEGER,ALLOCATABLE :: aindex(:)                  !Index connectingPS-row to basinindex
    INTEGER,ALLOCATABLE :: readsubid(:)
    INTEGER,ALLOCATABLE :: readpstype(:)
    INTEGER,ALLOCATABLE :: readsource(:)
    REAL    help                       !Total load
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL,ALLOCATABLE :: pointsource(:,:) !Read pointsource info (subbasins,column)
    CHARACTER(LEN=50),ALLOCATABLE :: readdate1(:)  !String date read from file
    CHARACTER(LEN=50),ALLOCATABLE :: readdate2(:)  !String date read from file
    TYPE(DateType) :: onedate                 !One date
    TYPE(DateType),ALLOCATABLE :: psdateslocal(:)  !Array of dates from file

    !Start of subroutine
    status = 0
    !Note. If(when) old format is removed from GeoData-reading, zeroing of load%psvol and %psload must be done here.

    !Read PointsourceData.txt, subroutine will allocate readsubid etc.
    CALL read_pointsourcedata(dir,infile,nrows,status,fexist,readsubid,readpstype,readsource,pointsource,readdate1,readdate2) 
    IF(.NOT.fexist) RETURN
    IF(status/=0)RETURN
    
    !Allocate and initiate variables 
    IF(.NOT.ALLOCATED(aindex)) ALLOCATE(aindex(nrows))
    DO i=1,n    !Zeroing the point source loads for accumulation
      load(i)%psvol(:) = 0.
      load(i)%psload(:,:) = 0.
    ENDDO
 
    !Collect dates where point source load changes; fromdata and todate+1
    maxdates = 0
    IF(.NOT.ALLOCATED(psdateslocal)) ALLOCATE(psdateslocal(2*nrows))
    DO i = 1,nrows
      IF(readdate1(i)(1:1)=='0')THEN
      ELSE
        CALL string_convert_to_datetype(readdate1(i),onedate)
        CALL add_date_to_array(onedate,2*nrows,maxdates,psdateslocal)
      ENDIF
      IF(readdate2(i)(1:1)=='0')THEN
      ELSE
        CALL string_convert_to_datetype(readdate2(i),onedate)
        onedate = AddDates(onedate,steplen)
        CALL add_date_to_array(onedate,2*nrows,maxdates,psdateslocal)
      ENDIF
    ENDDO
    IF(maxdates>0)THEN
      IF(.NOT.ALLOCATED(psdates)) ALLOCATE(psdates(maxdates))
      psdates=psdateslocal(1:maxdates)
    ENDIF
    DEALLOCATE(psdateslocal)
    
    IF(.NOT.ALLOCATED(psdates))THEN
      !Calculate permanent point source loads
      onetspday = 1./REAL(timesteps_per_day)
      CALL find_reorder_index(infile,nrows,readsubid,n,basin%subid,.FALSE.,aindex,status) !status=1 is not an error
      DO j=1,nrows
        IF(aindex(j)>0)THEN !Point sources not in model are ignored
          IF(readdate1(j)=='0'.AND.readdate2(j)=='0')THEN
            i=aindex(j)
            IF(pointsource(j,5)<0.)THEN !Negative point source, abstraction
              load(i)%abstrvol = - pointsource(j,5)/86400.  !m3/d->m3/s
              load(i)%abstrind = readsource(j)
            ELSEIF(pointsource(j,5)>0.)THEN
              k=readpstype(j)
              IF(k>0)THEN
                load(i)%psvol(k) = load(i)%psvol(k) + pointsource(j,5)/86400.  !m3/d->m3/s
                IF(i_in>0)THEN
                  help = pointsource(j,3)*pointsource(j,5)*onetspday*1.E-3    !kg TN/timestep
                  load(i)%psload(k,i_in) = load(i)%psload(k,i_in) + help*pointsource(j,4)
                  load(i)%psload(k,i_on) = load(i)%psload(k,i_on) + help*(1.-pointsource(j,4))
                ENDIF
                IF(i_sp>0)THEN
                  help = pointsource(j,1)*pointsource(j,5)*onetspday*1.E-3    !kg TP/timestep
                  load(i)%psload(k,i_sp) = load(i)%psload(k,i_sp) + help*pointsource(j,2)
                  load(i)%psload(k,i_pp) = load(i)%psload(k,i_pp) + help*(1.-pointsource(j,2))
                ENDIF
                IF(i_t1>0)THEN
                  help = pointsource(j,6)*pointsource(j,5)*onetspday*1.E-3    !amount T1/timestep
                  load(i)%psload(k,i_t1) = load(i)%psload(k,i_t1) + help
                ENDIF
                IF(i_t2>0)THEN
                  help = pointsource(j,7)*pointsource(j,5)*onetspday*1.E-3    !amount T2/timestep
                  load(i)%psload(k,i_t2) = load(i)%psload(k,i_t2) + help
                ENDIF
              ELSE
                WRITE(6,*) 'Warning: An permanent point source without ps_type found'
                WRITE(6,*) 'Warning: This is not a recommended use of PointSourceData.txt'
                WRITE(6,*) 'Warning: Point source may be ignored!'
              ENDIF
            ENDIF
          ENDIF
        ELSE
          WRITE(6,*) 'Point source in subid',readsubid(j),'ignored, not in model domain.'
        ENDIF
      ENDDO
    ENDIF
    
    !Deallocate local arrays
    IF(ALLOCATED(aindex)) DEALLOCATE(aindex)
    IF(ALLOCATED(readsubid)) DEALLOCATE(readsubid)
    IF(ALLOCATED(readpstype)) DEALLOCATE(readpstype)
    IF(ALLOCATED(readsource)) DEALLOCATE(readsource)
    IF(ALLOCATED(readdate1)) DEALLOCATE(readdate1)
    IF(ALLOCATED(readdate2)) DEALLOCATE(readdate2)
    IF(ALLOCATED(pointsource)) DEALLOCATE(pointsource)

  END SUBROUTINE load_pointsourcedata

  !>\brief Reads the matrix of point source data from the file. 
  !!
  !> \b Reference ModelDescription Chapter 5 (Point sources)
  !------------------------------------------------------------------------------
  SUBROUTINE read_pointsourcedata(dir,infile,nrows,status,fexist,readsubid,readpstype,readsource,pointsource,readdate1,readdate2) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         fileunit_temp

    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of pointsource file to be read (PointSourceData.txt)
    INTEGER, INTENT(OUT) :: nrows                !<Number of rows in file
    INTEGER, INTENT(OUT) :: status               !<Error status
    LOGICAL, INTENT(OUT) :: fexist               !<Status of file exist
    INTEGER,ALLOCATABLE, INTENT(INOUT) :: readsubid(:)    !<subid if point sources
    INTEGER,ALLOCATABLE, INTENT(INOUT) :: readpstype(:)   !<type of point source
    INTEGER,ALLOCATABLE, INTENT(INOUT) :: readsource(:)   !<abstraction source 
    REAL,ALLOCATABLE, INTENT(INOUT) :: pointsource(:,:)   !<Read pointsource info (subbasins,column)
    CHARACTER(LEN=50),ALLOCATABLE, INTENT(INOUT) :: readdate1(:)  !<From date of point source
    CHARACTER(LEN=50),ALLOCATABLE, INTENT(INOUT) :: readdate2(:)  !<To date of point source
    
    !Local constants
    INTEGER, PARAMETER :: dim = 50
    
    !Local variables
    INTEGER i
    INTEGER ncols,mcols
    INTEGER,ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER,ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER,ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER,ALLOCATABLE :: iindex(:)             !Index for column integer variables
    INTEGER,ALLOCATABLE :: sindex(:)             !Index for column string variables
    REAL,ALLOCATABLE :: xr(:,:)               !Real data read from file
    CHARACTER(LEN=50),ALLOCATABLE :: xs(:,:)  !String data read from file
    CHARACTER(LEN=10),ALLOCATABLE :: str(:)   !Content string

    !Start of subroutine
    status = 0
    
    !Check existence of file
    INQUIRE(FILE=TRIM(dir)//infile,EXIST=fexist)
    IF(.NOT.fexist) RETURN
    
    !Count number of columns and rows in PointSourceData
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,0,ncols,status)
    IF(status/=0)RETURN
    CALL count_data_rows(fileunit_temp,TRIM(dir)//infile,1,nrows,status)
    IF(status/=0)RETURN

    !Allocate variables for holding file data
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(sindex)) ALLOCATE(sindex(ncols))
    IF(.NOT.ALLOCATED(readsubid)) ALLOCATE(readsubid(nrows))
    IF(.NOT.ALLOCATED(readpstype)) ALLOCATE(readpstype(nrows))
    IF(.NOT.ALLOCATED(readsource)) ALLOCATE(readsource(nrows))
    IF(.NOT.ALLOCATED(readdate1)) ALLOCATE(readdate1(nrows))
    IF(.NOT.ALLOCATED(readdate2)) ALLOCATE(readdate2(nrows))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    IF(.NOT.ALLOCATED(xs)) ALLOCATE(xs(nrows,ncols))
    IF(.NOT.ALLOCATED(pointsource)) ALLOCATE(pointsource(nrows,7))
    readpstype = 0  !default point source is first
    readsource = 1  !default abstraction source is main river
    pointsource = 0.
    readdate1 = '0'
    readdate2 = '0'
    
    !Read PointSourceData.txt file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(dir)//infile, STATUS = 'old')     !Open PointSourceData-file

    !Reads the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(TRIM(dir)//infile)
      RETURN
    ENDIF

    !Code column for variable type
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='ps_source ') code(i) = i_intg
      IF(str(i)(1:10)=='ps_type   ') code(i) = i_intg
      IF(str(i)(1:10)=='ps_tpconc ') code(i) = i_real
      IF(str(i)(1:10)=='ps_spfrac ') code(i) = i_real
      IF(str(i)(1:10)=='ps_tnconc ') code(i) = i_real
      IF(str(i)(1:10)=='ps_infrac ') code(i) = i_real
      IF(str(i)(1:10)=='ps_vol    ') code(i) = i_real
      IF(str(i)(1:10)=='ps_t1     ') code(i) = i_real
      IF(str(i)(1:10)=='ps_t2     ') code(i) = i_real
    ENDDO

    !Read all data values
    CALL read_basindata6(fileunit_temp,TRIM(dir)//infile,ncols,nrows,ncols,50,code,rindex,iindex,sindex,xi,xr,xs)
    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(TRIM(dir)//infile)

    !Move data from read matrix to pointsource variables
    DO i = 1,ncols
      IF(str(i)(1:10)=='fromdate  ') readdate1(1:nrows)     = xs(1:nrows,sindex(i)) !String date
      IF(str(i)(1:10)=='todate    ') readdate2(1:nrows)     = xs(1:nrows,sindex(i))
      IF(str(i)(1:10)=='subid     ') readsubid(1:nrows)     = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='ps_type   ') readpstype(1:nrows)    = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='ps_source ') readsource(1:nrows)    = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='ps_tpconc ') pointsource(1:nrows,1) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_spfrac ') pointsource(1:nrows,2) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_tnconc ') pointsource(1:nrows,3) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_infrac ') pointsource(1:nrows,4) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_vol    ') pointsource(1:nrows,5) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_t1     ') pointsource(1:nrows,6) = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ps_t2     ') pointsource(1:nrows,7) = xr(1:nrows,rindex(i))
    ENDDO
    
    !Deallocate local arrays
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(sindex)) DEALLOCATE(sindex)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(xs)) DEALLOCATE(xs)

  END SUBROUTINE read_pointsourcedata

  !>\brief Gets the information about temporary point sources from PointSourceData.txt.
  !>Recalculate into form suitable for model.
  !!
  !>\b Consequences Module modvar variable load and t1load are changed
  !>
  !> \b Reference ModelDescription Chapter 5 (Point sources)
  !------------------------------------------------------------------------------
  SUBROUTINE get_current_pointsources(dir,infile,ns,time,status) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         bdate,  &
                         psdates,     &
                         fileunit_temp
    USE MODVAR, ONLY : timesteps_per_day, &
                       i_in,i_on,i_sp,i_pp, &
                       i_t1,i_t2,   &
                       max_pstype,  &
                       basin,       &
                       load,        & !OUT
                       tload,       & !OUT
                       tloadexist     !OUT
    USE READWRITE_ROUTINES, ONLY : find_reorder_index
 
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir         !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of pointsource file to be read (PointSourceData.txt)
    INTEGER, INTENT(IN)  :: ns                   !<Number of subbasins (submodel)
     TYPE(DateType), INTENT(IN) :: time           !<Current time (e.g. 2005-01-26 18:00)
    INTEGER, INTENT(OUT) :: status               !<Error status
    
    !Local constants
    INTEGER, PARAMETER :: dim = 50
    
    !Local variables
    LOGICAL fexist                     !Status of file exist
    INTEGER i,j,k
    INTEGER nrows
    INTEGER maxdates
    INTEGER numt1
    INTEGER,ALLOCATABLE :: aindex(:)                  !Index connectingPS-row to basinindex
    INTEGER,ALLOCATABLE :: readsubid(:)
    INTEGER,ALLOCATABLE :: readpstype(:)
    INTEGER,ALLOCATABLE :: readsource(:)
    INTEGER,ALLOCATABLE :: tloadindexlocal(:,:)
    REAL    help                       !Total load
    REAL    onetspday                  !Reciprocal of timesteps per day
    REAL,ALLOCATABLE :: pointsource(:,:) !Read pointsource info (subbasins,column)
    REAL,ALLOCATABLE :: tloadlocal(:,:)  !Read pointsource info on T1 sources
    CHARACTER(LEN=50),ALLOCATABLE :: readdate1(:)  !String date read from file
    CHARACTER(LEN=50),ALLOCATABLE :: readdate2(:)  !String date read from file
    TYPE(DateType),ALLOCATABLE :: psdateslocal(:,:)  !Array of dates from file

    !Start of subroutine
    status = 0
    
    !Check if beginning of simulation or current time in change array
    IF(time>bdate)THEN
      maxdates = SIZE(psdates)
      DO i = 1,maxdates
        IF(time==psdates(i)) EXIT
      ENDDO
      IF(i>maxdates) RETURN
    ENDIF
    
    !Read PointsourceData.txt, subroutine will allocate readsubid etc.
    CALL read_pointsourcedata(dir,infile,nrows,status,fexist,readsubid,readpstype,readsource,pointsource,readdate1,readdate2) 
    IF(status/=0)RETURN

    !Allocate and initiate variables 
    IF(.NOT.ALLOCATED(aindex)) ALLOCATE(aindex(nrows))
    IF(.NOT.ALLOCATED(psdateslocal)) ALLOCATE(psdateslocal(nrows,2))
    IF(.NOT.ALLOCATED(tloadindexlocal)) ALLOCATE(tloadindexlocal(2,nrows))
    IF(.NOT.ALLOCATED(tloadlocal)) ALLOCATE(tloadlocal(3,nrows))
    DO i=1,ns    !Zeroing the point source loads for accumulation
      load(i)%psvol(:) = 0.
      load(i)%psload(:,:) = 0.
    ENDDO

    !Transform dates to DateType
    DO j=1,nrows
      IF(readdate1(j)=='0')THEN
        psdateslocal(j,1)=datetype(3000,0,0,0,0)
      ELSE
        CALL string_convert_to_datetype(readdate1(j),psdateslocal(j,1))
      ENDIF
      IF(readdate2(j)=='0')THEN
        psdateslocal(j,2)=datetype(0,0,0,0,0)
      ELSE
        CALL string_convert_to_datetype(readdate2(j),psdateslocal(j,2))
      ENDIF
    ENDDO

    !Calculate and add permanent point source loads and 
    !temporary point source load discharging today
    numt1 = 0
    onetspday = 1./REAL(timesteps_per_day)
    CALL find_reorder_index(infile,nrows,readsubid,ns,basin%subid,.FALSE.,aindex,status) !status=1 is not an error
    DO j=1,nrows
      IF(aindex(j)>0)THEN !Point sources not in (sub)model are ignored
        i=aindex(j)
        IF(readdate1(j)=='0'.AND.readdate2(j)=='0'.OR. &              !permanent point source
          (readdate1(j)=='0'.AND.psdateslocal(j,2)>=time) .OR. &      !ending point source
          (psdateslocal(j,1)<=time.AND.readdate2(j)=='0') .OR. &      !starting point source
          (psdateslocal(j,1)<=time.AND.psdateslocal(j,2)>=time))THEN  !temporary point source
          IF(pointsource(j,5)<0.)THEN   !Negative point source, abstraction
            load(i)%abstrvol = - pointsource(j,5)/86400.  !m3/d->m3/s
            load(i)%abstrind = readsource(j)
          ELSEIF(pointsource(j,5)>0.)THEN
            k=readpstype(j)
            IF(k>0)THEN
              load(i)%psvol(k) = load(i)%psvol(k) + pointsource(j,5)/86400.  !m3/d->m3/s
              IF(i_in>0)THEN
                help = pointsource(j,3)*pointsource(j,5)*onetspday*1.E-3    !kg TN/timestep
                load(i)%psload(k,i_in) = load(i)%psload(k,i_in) + help*pointsource(j,4)
                load(i)%psload(k,i_on) = load(i)%psload(k,i_on) + help*(1.-pointsource(j,4))
              ENDIF
              IF(i_sp>0)THEN
                help = pointsource(j,1)*pointsource(j,5)*onetspday*1.E-3    !kg TP/timestep
                load(i)%psload(k,i_sp) = load(i)%psload(k,i_sp) + help*pointsource(j,2)
                load(i)%psload(k,i_pp) = load(i)%psload(k,i_pp) + help*(1.-pointsource(j,2))
              ENDIF
              IF(i_t1>0)THEN
                help = pointsource(j,6)*pointsource(j,5)*onetspday*1.E-3    !amount T1/timestep
                load(i)%psload(k,i_t1) = load(i)%psload(k,i_t1) + help
              ENDIF
              IF(i_t2>0)THEN
                help = pointsource(j,7)*pointsource(j,5)*onetspday*1.E-3    !amount T1/timestep
                load(i)%psload(k,i_t2) = load(i)%psload(k,i_t2) + help
              ENDIF
            ELSEIF(readsource(j)>0)THEN
              IF(i_t1>0 .OR. i_t2>0)THEN
                numt1 = numt1 + 1
                tloadlocal(1,numt1) = pointsource(j,5)/86400.  !m3/d->m3/s
                tloadlocal(2,numt1) = pointsource(j,6)
                tloadlocal(3,numt1) = pointsource(j,7)
                tloadindexlocal(1,numt1) = readsource(j)
                tloadindexlocal(2,numt1) = i
              ENDIF
            ENDIF
          ENDIF
        ENDIF
      ELSE
        WRITE(6,*) 'Point source in subid',readsubid(j),'ignored, not in model domain.'
      ENDIF
    ENDDO

    !Set T1 loads if present
    IF(numt1>0)THEN
      IF(.NOT.ALLOCATED(tload))THEN
        ALLOCATE(tload(numt1))
        ALLOCATE(tloadexist(ns))
        tloadexist = .FALSE.
      ENDIF
      tload%psvol = tloadlocal(1,1:numt1)
      tload%psconc = tloadlocal(2,1:numt1)
      tload%pstemp = tloadlocal(3,1:numt1)
      tload%sw_code = tloadindexlocal(1,1:numt1)
      tload%subindex = tloadindexlocal(2,1:numt1)
      tloadexist(tloadindexlocal(2,1:numt1)) = .TRUE.
    ELSE
      IF(ALLOCATED(tload))THEN
        DEALLOCATE(tload)
        DEALLOCATE(tloadexist)
      ENDIF
    ENDIF
    
    !Deallocate local variables
    IF(ALLOCATED(readsubid)) DEALLOCATE(readsubid)
    IF(ALLOCATED(readpstype)) DEALLOCATE(readpstype)
    IF(ALLOCATED(readsource)) DEALLOCATE(readsource)
    IF(ALLOCATED(readdate1)) DEALLOCATE(readdate1)
    IF(ALLOCATED(readdate2)) DEALLOCATE(readdate2)
    IF(ALLOCATED(pointsource)) DEALLOCATE(pointsource)
    IF(ALLOCATED(aindex)) DEALLOCATE(aindex)
    IF(ALLOCATED(psdateslocal)) DEALLOCATE(psdateslocal)
    IF(ALLOCATED(tloadindexlocal)) DEALLOCATE(tloadindexlocal)
    IF(ALLOCATED(tloadlocal)) DEALLOCATE(tloadlocal)
    
  END SUBROUTINE get_current_pointsources

  !>\brief Gets the information about irrigation from MgmtData.
  !!Reads the matrix of basin data values from the file vith mcols
  !!columns. Output: data variables, heading variable and actual
  !!number of data columns
  !--------------------------------------------------------------------
  SUBROUTINE load_irrigation_data(funit,dir,status) 

    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         maxcharpath
    USE MODVAR, ONLY : irrigationsystem    !OUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: dir           !<File directory
    INTEGER, INTENT(OUT) :: status                 !<Error status
    
    !Local variables
    INTEGER ncols,nrows                !size of data in file
    INTEGER mcol
    INTEGER i,j
    LOGICAL fileex                     !Status of file
    LOGICAL swcol,gwcol                !Status of sw/gw_part columns
    CHARACTER(LEN=maxcharpath) infile             !Name of irrigation characteristics file 
    INTEGER, ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER, ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER, ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER, ALLOCATABLE :: iindex(:)             !Index for column integer variables
    REAL, ALLOCATABLE    :: xr(:,:)               !Real data read from file
    CHARACTER(LEN=10), ALLOCATABLE :: str(:)      !File column content string

    status = 0
    infile = TRIM(dir)//'MgmtData.txt'

    !Check if file exist
    INQUIRE(FILE=TRIM(infile),EXIST=fileex)
    IF(.NOT.fileex)THEN
       WRITE(6,*) 'No MgmtData.txt file found.'
       RETURN
    ENDIF

    !Count number of columns in MgmtData
    CALL count_data_cols(funit,TRIM(infile),0,ncols,status)
    IF(status/=0)RETURN

    !Count number of subbasins in MgmtData
    CALL count_data_rows(funit,TRIM(infile),1,nrows,status)
    IF(status/=0)RETURN
    IF(nrows==0)THEN
       WRITE(6,*) 'No data in MgmtData.txt file found.'
       RETURN
    ENDIF

    !Allocate and initiate local variables
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))  
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))  
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols)) 
    swcol = .FALSE.; gwcol = .FALSE.

    !Open file
    OPEN(UNIT = funit,FILE = TRIM(infile), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(funit,ncols,str,mcol,status)
    IF(status.NE.0) THEN
       WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
       RETURN
    ENDIF

    !Allocate irrigation variable 
    IF(.NOT.ALLOCATED(irrigationsystem)) ALLOCATE(irrigationsystem(nrows))

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,ncols
       IF(str(i)(1:10)=='subid     ') code(i) = i_intg
       IF(str(i)(1:10)=='regsrcid  ') code(i) = i_intg
       IF(str(i)(1:10)=='local_eff ') code(i) = i_real
       IF(str(i)(1:10)=='region_eff') code(i) = i_real
       IF(str(i)(1:10)=='gw_part   ') code(i) = i_real
       IF(str(i)(1:10)=='demandtype') code(i) = i_intg
       IF(str(i)(1:10)=='irrdam    ') code(i) = i_intg
    ENDDO

    !Read all data
    CALL read_basindata5(funit,infile,ncols,nrows,ncols,code,rindex,iindex,xi,xr)

    CLOSE(UNIT=funit)
    WRITE(6,*) 'File read: ', TRIM(infile)

    DO i = 1,ncols
       IF(str(i)(1:10)=='subid     ') irrigationsystem(1:nrows)%subid       = xi(1:nrows,iindex(i))
       IF(str(i)(1:10)=='regsrcid  ') irrigationsystem(1:nrows)%regsourceid = xi(1:nrows,iindex(i))
       IF(str(i)(1:10)=='local_eff ') irrigationsystem(1:nrows)%local_eff   = xr(1:nrows,rindex(i))
       IF(str(i)(1:10)=='region_eff') irrigationsystem(1:nrows)%reg_eff     = xr(1:nrows,rindex(i))
       IF(str(i)(1:10)=='demandtype') irrigationsystem(1:nrows)%demandtype  = xi(1:nrows,iindex(i))
       IF(str(i)(1:10)=='gw_part   ')THEN
          irrigationsystem(1:nrows)%gw_part     = xr(1:nrows,rindex(i))
          gwcol = .TRUE.
       ENDIF
       IF(str(i)(1:10)=='irrdam    ')THEN
          DO j = 1,nrows
             irrigationsystem(j)%dam     = xi(j,iindex(i))==1
          ENDDO
       ENDIF
    ENDDO

    !Check irrigation used (if file used for other things in the future)
    IF(gwcol)THEN
       irrigationsystem(1:nrows)%sw_part = 1. - irrigationsystem(1:nrows)%gw_part
    ELSE
       WRITE(6,*) 'No irrigation source division found. Irrigation not used.'
       IF(ALLOCATED(irrigationsystem)) DEALLOCATE(irrigationsystem)
    ENDIF

    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)  
    IF(ALLOCATED(code)) DEALLOCATE(code)  
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(str)) DEALLOCATE(str) 

  END SUBROUTINE load_irrigation_data

!--------------------------------------------------------------------
  !>\brief Gets information about forcing data coupling to subbasin 
  !!and other information.
  !!Reads the matrix of basin data values from the file with mcols
  !!columns. Output: data variables, heading variable and actual
  !!number of data columns
  !!
  !>\b Consequences Module worldvar variables tobsid, pobsid
  !> are set. 
  !>Module modvar variable tobselevation is allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE load_forcing_information_data(funit,dir,ns,status) 

    USE WORLDVAR, ONLY : get_seq_filename, &
                         forcingdata, &
                         max_forcingdata, &
                         i_pobs, &
                         i_tobs, &
                         i_tminobs, &
                         i_tmaxobs, &
                         i_rhobs, &
                         i_sfobs, &
                         i_swobs, &
                         i_uobs, &
                         readobsid,   &
                         i_str,       &
                         i_intg,      &
                         i_real,      &
                         maxcharpath
    USE MODVAR, ONLY : basin,   &
                       tobselevation  !OUT

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                  !<Unit number for file
    CHARACTER (LEN=*), INTENT(IN) :: dir           !<File directory
    INTEGER, INTENT(IN) :: ns                      !<Number of subbasins (base model)
    INTEGER, INTENT(OUT) :: status                 !<Error status
    
    !Local variables
    INTEGER ncols                      !size of data in file
    INTEGER mcols
    INTEGER i,j,iforc
    INTEGER isubid,itobselev  !Columns in file
    INTEGER ifobsid(max_forcingdata)   !Columns in file
    LOGICAL fileex                     !Existence of file
    CHARACTER(LEN=maxcharpath) infile  !Name and path of file 
    CHARACTER(LEN=16) filename         !Name of file 
    INTEGER, ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER, ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER, ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER, ALLOCATABLE :: iindex(:)             !Index for column integer variables
    REAL, ALLOCATABLE    :: xr(:,:)               !Real data read from file
    CHARACTER(LEN=10), ALLOCATABLE :: colstr(:)   !File column content string

    status = 0
    !Set default values for obsid
    DO iforc=1,max_forcingdata
      IF(forcingdata(iforc)%readfile)THEN
        IF(.NOT.ALLOCATED(forcingdata(iforc)%stationid))ALLOCATE(forcingdata(iforc)%stationid(ns))
        forcingdata(iforc)%stationid(1:ns) = basin(1:ns)%subid
      ENDIF
    ENDDO

    !Check if file exist and need to be read
    filename = 'ForcKey.txt'
    CALL get_seq_filename(filename)
    infile = TRIM(dir)//TRIM(filename)
    INQUIRE(FILE=TRIM(infile),EXIST=fileex)
    IF(.NOT.fileex)THEN
      IF(readobsid) WRITE(6,*) 'WARNING: No ForcKey-file found, even though readobsid turned on in info.txt.'
      RETURN
    ELSEIF(fileex.AND..NOT.readobsid)THEN
      WRITE(6,*) 'WARNING: Existing ForcKey-file will not be read, because readobsid turned off in info.txt'
      RETURN
    ENDIF

    !Count number of columns in file
    CALL count_data_cols(funit,TRIM(infile),0,ncols,status)
    IF(status/=0)RETURN

    !Allocation of local variables
    ALLOCATE(xi(ns,ncols))  
    ALLOCATE(code(ncols))  
    ALLOCATE(rindex(ncols))
    ALLOCATE(iindex(ncols))
    ALLOCATE(xr(ns,ncols))
    ALLOCATE(colstr(ncols)) 

    !Read file
    OPEN(UNIT=funit,FILE=TRIM(infile),STATUS='old',ACTION='read',ERR=900)
    CALL read_column_headings(funit,ncols,colstr,mcols,status)
    IF(status/=0) RETURN
    code = i_str
    isubid = 0
    itobselev = 0
    ifobsid = 0
    DO i = 1,mcols
      IF(colstr(i)(1:10)=='subid     ')THEN
        code(i) = i_intg
        isubid = i
      ENDIF
      DO iforc=1,max_forcingdata
        IF(colstr(i)(1:10)==forcingdata(iforc)%idcode)THEN
          code(i) = i_intg
          IF(forcingdata(iforc)%readfile) ifobsid(iforc) = i
        ENDIF  
      ENDDO
      IF(colstr(i)(1:10)=='tobselev  ')THEN
        code(i) = i_real
        itobselev = i
      ENDIF  
    ENDDO
    IF(isubid==0)THEN
      WRITE(6,*) 'ERROR: subid missing in file:', TRIM(infile)
      status = 1
      RETURN
    ENDIF
    CALL read_basindata5(funit,infile,ncols,ns,mcols,code,rindex,iindex,xi,xr) 
    CLOSE(funit) 
    WRITE(6,*) 'File read: ', TRIM(infile)
    IF(itobselev==0) WRITE(6,*) 'WARNING: no elevation for tobs found in file:', TRIM(infile)

    !Set forcing data information to subbasin array
    IF(itobselev>0)THEN
      IF(.NOT.ALLOCATED(tobselevation)) ALLOCATE(tobselevation(ns))
      DO i = 1,ns
        IF(basin(i)%subid==xi(i,iindex(isubid)))THEN
          tobselevation(i) = xr(i,rindex(itobselev))
        ELSE
          DO j = 1,ns
            IF(basin(i)%subid==xi(j,iindex(isubid)))THEN
              tobselevation(i) = xr(j,rindex(itobselev))
              EXIT
            ENDIF  
          ENDDO
        ENDIF
      ENDDO
    ENDIF

    !Set subbasin - observation coupling variables
    DO i = 1,ns
      IF(basin(i)%subid==xi(i,iindex(isubid)))THEN
        DO iforc=1,max_forcingdata
          IF(ifobsid(iforc)>0) forcingdata(iforc)%stationid(i) = xi(i,iindex(ifobsid(iforc)))
        ENDDO
      ELSE
        DO j = 1,ns
          IF(basin(i)%subid==xi(j,iindex(isubid)))THEN
            DO iforc=1,max_forcingdata
              IF(ifobsid(iforc)>0) forcingdata(iforc)%stationid(i) = xi(j,iindex(ifobsid(iforc)))
            ENDDO
            EXIT
          ENDIF  
        ENDDO
      ENDIF
    ENDDO

    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)  
    IF(ALLOCATED(code)) DEALLOCATE(code)  
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(colstr)) DEALLOCATE(colstr) 

    RETURN
    
900 WRITE(6,*) 'Error open file: ', TRIM(infile)
    status = 1
    RETURN

  END SUBROUTINE load_forcing_information_data

  !>\brief Gets the information about update from update.txt.
  !-----------------------------------------------------------------
  SUBROUTINE read_update_data(funit,infile,ns,nqsub,qupdsub,nqarsub,nwarsub, &
       qarupdsub,warupdsub,qarupdar,warupdar,nwsub,wendsub,nrow,tpcorrarr,tncorrarr,tploccorrarr,tnloccorrarr)

    USE WORLDVAR, ONLY : i_str,i_intg,i_real

    !Argument declarations
    INTEGER, INTENT(IN)  :: funit                !<Unit for file
    CHARACTER (LEN=*), INTENT(IN) :: infile      !<Name of characteristics file to be read
    INTEGER, INTENT(IN)  :: ns                   !<Number of subbasins (basemodel)
    INTEGER, INTENT(OUT) :: nqsub                !<Number of subbasins with quseobs update
    INTEGER, INTENT(OUT) :: qupdsub(ns)          !<Subid of subbasins with quseobs update
    INTEGER, INTENT(OUT) :: nqarsub              !<Number of subbasins with qAR date
    INTEGER, INTENT(OUT) :: nwarsub              !<Number of subbasins with wAR date
    INTEGER, INTENT(OUT) :: qarupdsub(ns)        !<Subid of subbasins with qAR update
    INTEGER, INTENT(OUT) :: warupdsub(ns)        !<Subid of subbasins with wAR update
    REAL, INTENT(OUT)    :: qarupdar(ns)         !<AR-factor for qAR-update
    REAL, INTENT(OUT)    :: warupdar(ns)         !<AR-factor for wAR-update
    INTEGER, INTENT(OUT) :: nwsub                !<Number of subbasins with wendupd update
    INTEGER, INTENT(OUT) :: wendsub(ns)          !<Subid of subbasins with wendupd update
    INTEGER, INTENT(OUT) :: nrow                 !<Number of data rows in file
    REAL,    INTENT(OUT) :: tpcorrarr(ns,2)      !<Subid,value of subbasins with tpcorr update
    REAL,    INTENT(OUT) :: tncorrarr(ns,2)      !<Subid,value of subbasins with tncorr update
    REAL,    INTENT(OUT) :: tploccorrarr(ns,2)   !<Subid,value of subbasins with tploccorr update
    REAL,    INTENT(OUT) :: tnloccorrarr(ns,2)   !<Subid,value of subbasins with tnloccorr update

    !Local parameters
    INTEGER, PARAMETER :: nskip  = 1      !Number of initial heading rows 

    !Local variables declarations
    INTEGER i,irow                     !Loop counters
    INTEGER maxcol
    INTEGER isub                       !Data column with subid
    INTEGER iquseobs                   !Data column with status for quseobs
    INTEGER iqar,iwar                  !Data column with status for q-AR and w-AR
    INTEGER iqarar                     !Data column with AR-factor for q-AR
    INTEGER iwendupd                   !Data column with status for wendupd
    INTEGER itpcorr                    !Data column with value for tp correction
    INTEGER itncorr                    !Data column with value for tn correction
    INTEGER itploccorr                 !Data column with value for local tp correction
    INTEGER itnloccorr                 !Data column with value for local tn correction
    INTEGER mcols                      !Actual number of columns
    INTEGER status                     !Error status
    INTEGER, ALLOCATABLE ::  code(:)   !Code for column variable
    INTEGER, ALLOCATABLE ::  rindex(:) !Index for column real variables
    INTEGER, ALLOCATABLE ::  iindex(:) !Index for column integer variables
    INTEGER, ALLOCATABLE :: xi(:,:)    !Integer data read from file
    INTEGER, ALLOCATABLE :: temparray(:) !Temporary data array
    REAL, ALLOCATABLE :: temparray2(:) !Temporary data array2 !JD2012_AR	  
    REAL, ALLOCATABLE    :: xr(:,:)    !Real data read from file
    CHARACTER(LEN=10), ALLOCATABLE :: str(:)      !Content string

    !Initiations
    status = 0
    nqsub  = 0
    qupdsub = 0
    nqarsub  = 0 
    nwarsub  = 0 
    qarupdsub = 0
    warupdsub = 0
    qarupdar = 0 
    warupdar = 0 
    nrow = 0
    tpcorrarr = 0.
    tncorrarr = 0.
    tploccorrarr = 0.
    tnloccorrarr = 0.

    !Count the number of columns and rows and allocate variables accordingly
    CALL count_data_cols(funit,infile,0,maxcol,status)
    IF(status/=0)THEN
      WRITE(6,*) 'ERROR: counting columns in update.txt'
      RETURN
    ENDIF
    CALL count_data_rows(funit,infile,nskip,nrow,status)
    IF(status/=0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(maxcol))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(maxcol))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(maxcol))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(maxcol))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrow,maxcol))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrow,maxcol))

    !Open file for reading data
    OPEN(UNIT = funit,FILE = infile, STATUS = 'old', ACTION='read')     

    !Reads the column headings
    CALL read_column_headings(funit,maxcol,str,mcols,status)
    IF(status/=0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
      IF(str(i)(1:10)=='subid     ') code(i) = i_intg
      IF(str(i)(1:10)=='quseobs   ') code(i) = i_intg
      IF(str(i)(1:10)=='wendupd   ') code(i) = i_intg
      IF(str(i)(1:10)=='tpcorr    ') code(i) = i_real
      IF(str(i)(1:10)=='tncorr    ') code(i) = i_real
      IF(str(i)(1:10)=='tploccorr ') code(i) = i_real
      IF(str(i)(1:10)=='tnloccorr ') code(i) = i_real
      IF(str(i)(1:10)=='qarupd    ') code(i) = i_intg !Q-AR on(1), off(0)
      IF(str(i)(1:10)=='warupd    ') code(i) = i_intg !W-AR on(1), off(0)
      IF(str(i)(1:10)=='arfact    ') code(i) = i_real !factor for AR
    ENDDO

    !Read all data
    CALL read_basindata5(funit,infile,maxcol,nrow,mcols,code,rindex,iindex,xi,xr)

    CLOSE(UNIT=funit)
    WRITE(6,*) 'File read: ', TRIM(infile)

    !Find variables for update function
    isub = 0; iquseobs = 0; iwendupd=0; itpcorr=0; itncorr=0; itploccorr=0; itnloccorr=0; iqar=0; iqarar=0; iwar=0
    DO i = 1,mcols
      IF(str(i)(1:10)=='subid     ')   isub = i
      IF(str(i)(1:10)=='quseobs   ')   iquseobs = i
      IF(str(i)(1:10)=='wendupd   ')   iwendupd = i
      IF(str(i)(1:10)=='tpcorr    ')   itpcorr = i
      IF(str(i)(1:10)=='tncorr    ')   itncorr = i
      IF(str(i)(1:10)=='tploccorr ')   itploccorr = i
      IF(str(i)(1:10)=='tnloccorr ')   itnloccorr = i
      IF(str(i)(1:10)=='qarupd    ')   iqar = i  
      IF(str(i)(1:10)=='warupd    ')   iwar = i  
      IF(str(i)(1:10)=='arfact    ')   iqarar = i
    ENDDO
    IF(isub==0)THEN
      WRITE(6,*) 'No subid information in update.txt'
      RETURN
    ENDIF

    !Set variables for update function quseobs      
    IF(iquseobs>0)THEN      !Save subid for stations to be updated
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      temparray = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iquseobs))==1)THEN
          nqsub = nqsub + 1
          temparray(nqsub) = xi(irow,iindex(isub))
        ENDIF
      ENDDO
      IF(nqsub>0)THEN
        qupdsub(1:nqsub)=temparray(1:nqsub)
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
    ENDIF

    !Set variables for update function qAR  
    IF(iqar>0)THEN
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      IF(.NOT.ALLOCATED(temparray2)) ALLOCATE(temparray2(nrow))
      temparray = 0
      temparray2 = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iqar))==1)THEN
          nqarsub = nqarsub + 1
          temparray(nqarsub) = xi(irow,iindex(isub))
          temparray2(nqarsub) = xr(irow,rindex(iqarar))
        ENDIF
      ENDDO
      IF(nqarsub>0)THEN
        qarupdsub(1:nqarsub) = temparray(1:nqarsub)
        qarupdar(1:nqarsub) = temparray2(1:nqarsub)
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
      IF(ALLOCATED(temparray2)) DEALLOCATE(temparray2)
    ENDIF

    !Set variables for update function wAR  
    IF(iwar>0)THEN
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      IF(.NOT.ALLOCATED(temparray2)) ALLOCATE(temparray2(nrow))
      temparray = 0
      temparray2 = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iwar))==1)THEN
          nwarsub = nwarsub + 1
          temparray(nwarsub) = xi(irow,iindex(isub))
          temparray2(nwarsub) = xr(irow,rindex(iqarar))
        ENDIF
      ENDDO
      IF(nwarsub>0)THEN
        warupdsub(1:nwarsub) = temparray(1:nwarsub)
        warupdar(1:nwarsub) = temparray2(1:nwarsub) !Same as for qAR
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
      IF(ALLOCATED(temparray2)) DEALLOCATE(temparray2)
    ENDIF

    !Set variables for update function tpcorr
    IF(itpcorr>0)THEN      !Save subid/value for all stations to be updated
      tpcorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tpcorrarr(1:nrow,2) = xr(1:nrow,rindex(itpcorr))
    ENDIF
    !Set variables for update function tncorr
    IF(itncorr>0)THEN      !Save subid/value for all stations to be updated
      tncorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tncorrarr(1:nrow,2) = xr(1:nrow,rindex(itncorr))
    ENDIF
    !Set variables for update function tploccorr
    IF(itploccorr>0)THEN      !Save subid/value for all stations to be updated
      tploccorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tploccorrarr(1:nrow,2) = xr(1:nrow,rindex(itploccorr))
    ENDIF
    !Set variables for update function tnloccorr
    IF(itnloccorr>0)THEN      !Save subid/value for all stations to be updated
      tnloccorrarr(1:nrow,1) = REAL(xi(1:nrow,iindex(isub)))
      tnloccorrarr(1:nrow,2) = xr(1:nrow,rindex(itnloccorr))
    ENDIF

    !Set variables for update function wendupd      
    IF(iwendupd>0)THEN      !Save subid for stations to be updated
      IF(.NOT.ALLOCATED(temparray)) ALLOCATE(temparray(nrow))
      temparray = 0
      DO irow = 1,nrow
        IF(xi(irow,iindex(iwendupd))==1)THEN
          nwsub = nwsub + 1
          temparray(nwsub) = xi(irow,iindex(isub))
        ENDIF
      ENDDO
      IF(nwsub>0)THEN
        wendsub(1:nwsub)=temparray(1:nwsub)
      ENDIF
      IF(ALLOCATED(temparray)) DEALLOCATE(temparray)
    ENDIF

    IF(ALLOCATED(code)) DEALLOCATE(code)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(str)) DEALLOCATE(str)
    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

  END SUBROUTINE read_update_data

  !>Calculate the coupling between subbasins by index from the existing
  !>coupling by subid
  !>
  !>\b Consequences Module modvar variables path, branchdata and branchindex 
  !>may be allocated and set. Module modvar variables pathsubid and branchsubid 
  !>is deallocated. 
  !--------------------------------------------------------------------
  SUBROUTINE calculate_path(ns) 

    USE MODVAR, ONLY : path,        &   !OUT
                       pathsubid,   &
                       branchdata,  &   !OUT
                       branchsubid, &
                       branchindex      !OUT

    !Argument declarations
    INTEGER, INTENT(IN) :: ns                   !<Number of subbasins (of submodel)
    
    !Local variables
    INTEGER i
    INTEGER dim   !size of branchdata table
    INTEGER mainflow(ns)
    INTEGER grwflow(ns)
    INTEGER,ALLOCATABLE :: source(:),branch(:)

    !>\b Algorithm \n
    !>Set main path and groundwater flow path
    IF(.NOT.ALLOCATED(path)) ALLOCATE(path(ns))
    mainflow = pathsubid%main
    grwflow = pathsubid%grw1
    path%grwtolake = pathsubid%grwtolake
    path%aquid     = pathsubid%aquid
    path%rechargebasin = pathsubid%rechargebasin
    path%recievefraction = pathsubid%recievefraction
    DO i = 1,ns
      path(i)%main   = get_subid_index(mainflow(i),ns)      !uses basin%subid!
      path(i)%grw1   = get_subid_index(grwflow(i),ns)
    ENDDO
    IF(ALLOCATED(pathsubid)) DEALLOCATE(pathsubid)

    !>Set branch path and data, if present
    dim = 0
    IF(ALLOCATED(branchsubid))THEN
      dim = SIZE(branchsubid)
      IF(.NOT.ALLOCATED(branchdata)) ALLOCATE(branchdata(dim))
      IF(.NOT.ALLOCATED(branchindex)) ALLOCATE(branchindex(ns))
      IF(.NOT.ALLOCATED(source)) ALLOCATE(source(dim))
      IF(.NOT.ALLOCATED(branch)) ALLOCATE(branch(dim))
      branchindex = 0
      source = branchsubid%source
      branch = branchsubid%branch
      
      !Copy branch flow data
      branchdata%mainpart   = branchsubid%mainpart
      branchdata%maxQ       = branchsubid%maxQ
      branchdata%minQ       = branchsubid%minQ
      branchdata%maxQbranch = branchsubid%maxQbranch
      
      !Calculate corresponding index
      DO i = 1,dim
        branchdata(i)%source = get_subid_index(source(i),ns)
        IF(branchdata(i)%source>0) branchindex(branchdata(i)%source) = i
        branchdata(i)%branch = get_subid_index(branch(i),ns)
      ENDDO
      IF(ALLOCATED(branchsubid)) DEALLOCATE(branchsubid)
      IF(ALLOCATED(source)) DEALLOCATE(source)
      IF(ALLOCATED(branch)) DEALLOCATE(branch)
      IF(SUM(branchindex)==0)THEN   !No branches in model
        IF(ALLOCATED(branchindex)) DEALLOCATE(branchindex)
        IF(ALLOCATED(branchdata)) DEALLOCATE(branchdata)
        WRITE(6,*) 'No branches in (sub-)model found'
      ENDIF
    ENDIF

  END SUBROUTINE calculate_path

  !>Reform the information about subbasins to arrays for submodel
  !-------------------------------------------------------------------
  SUBROUTINE reform_inputdata_for_submodel(nsfile,ns,indexarray)

    USE MODVAR, ONLY : basin,        &
                       classbasin,   &
                       pathsubid,    &
                       load,         &
                       glacierexist, &
                       wetland,      &
                       wetlandexist,     &
                       basintype,        &
                       classbasintype,   &
                       pathtype,         &
                       npcloadtype,      &
                       wetlandtype,      &
                       tobselevation, &
                       basinpar,         &
                       damindex,         & 
                       floodindex,       & 
                       lakeindex,        & 
                       lakeout2index,    & 
                       lakebasinindex,   & 
                       lakedataparindex, & 
                       glacierindex,     &
                       nclass,           &
                       numsubstances,    &
                       modparid,         &
                       m_bpar,           &
                       max_par,          &
                       allocate_basinvariables, &
                       xobsindex,    &
                       max_outvar,   &
                       max_pstype
    USE WORLDVAR, ONLY : qobsindex,    &
                         forcingdata, &
                         max_forcingdata, &
                         noutreg, &
                         outregion

    !Argument declarations
    INTEGER, INTENT(IN) :: nsfile           !<number of subbasins in file
    INTEGER, INTENT(IN) :: ns               !<number of subbasins to be simulated
    INTEGER, INTENT(IN) :: indexarray(ns)   !<index for basemodel
    
    !Local variables 
    INTEGER i,isub,iforc,k,s,nbasinpar
    TYPE(BASINTYPE),ALLOCATABLE :: temp_basin(:)  !help variables for reindexing
    TYPE(PATHTYPE),ALLOCATABLE  :: temp_pathsubid(:)
    TYPE(NPCLOADTYPE),ALLOCATABLE  :: temp_load(:)
    TYPE(CLASSBASINTYPE),ALLOCATABLE :: temp_classbasin(:,:)
    TYPE(WETLANDTYPE), ALLOCATABLE :: temp_wetland(:,:)
    REAL,ALLOCATABLE    :: temp_tobselev(:)
    REAL,ALLOCATABLE    :: temp_basinpar(:,:)
    INTEGER,ALLOCATABLE :: temp_lakeindex(:)
    INTEGER,ALLOCATABLE :: temp_glacierindex(:)
    INTEGER,ALLOCATABLE :: temp_lakebasinindex(:)
    INTEGER,ALLOCATABLE :: temp_lakedataparindex(:,:)
    INTEGER,ALLOCATABLE :: tmp(:),xtmp(:,:)

    !Copy subbasin information from GeoData to temporary arrays
    IF(.NOT.ALLOCATED(temp_basin)) ALLOCATE(temp_basin(nsfile))
    temp_basin      = basin
    DEALLOCATE(basin)
    IF(.NOT.ALLOCATED(temp_pathsubid)) ALLOCATE(temp_pathsubid(nsfile))
    temp_pathsubid  = pathsubid
    DEALLOCATE(pathsubid)
    IF(.NOT.ALLOCATED(temp_classbasin)) ALLOCATE(temp_classbasin(nsfile,nclass))
    temp_classbasin = classbasin
    DEALLOCATE(classbasin)
    ALLOCATE(temp_load(nsfile))
    DO i = 1,nsfile
      ALLOCATE(temp_load(i)%psvol(max_pstype))
      ALLOCATE(temp_load(i)%psload(max_pstype,numsubstances))
      ALLOCATE(temp_load(i)%locconc(numsubstances))
      DO k=1,max_pstype
        temp_load(i)%psvol(k) = load(i)%psvol(k)
        DO s=1,numsubstances
          temp_load(i)%psload(k,s) = load(i)%psload(k,s)
          IF(k==1) temp_load(i)%locconc(s) = load(i)%locconc(s)
        ENDDO
      ENDDO
    ENDDO
    temp_load%inwetdep    = load%inwetdep
    temp_load%indrydep(1) = load%indrydep(1)
    temp_load%indrydep(2) = load%indrydep(2)
    temp_load%indrydep(3) = load%indrydep(3)
    temp_load%volloc      = load%volloc
    temp_load%abstrvol = load%abstrvol
    temp_load%abstrind = load%abstrind
    DO i = 1,nsfile
      DEALLOCATE(load(i)%psvol)
      DEALLOCATE(load(i)%psload)
      DEALLOCATE(load(i)%locconc)
    ENDDO
    DEALLOCATE(load)
    IF(.NOT.ALLOCATED(temp_wetland)) ALLOCATE(temp_wetland(nsfile,2))
    IF(wetlandexist)THEN  
      temp_wetland = wetland
      DEALLOCATE(wetland)
    ENDIF

    !Allocate for and copy reform subbasin information
    CALL allocate_basinvariables(ns,numsubstances)
    basin(:)        = temp_basin(indexarray(:))
    DEALLOCATE(temp_basin)
    pathsubid(:)    = temp_pathsubid(indexarray(:))
    DEALLOCATE(temp_pathsubid)
    classbasin(:,:) = temp_classbasin(indexarray(:),:)
    DEALLOCATE(temp_classbasin)
    load(:)%inwetdep    = temp_load(indexarray(:))%inwetdep
    load(:)%indrydep(1) = temp_load(indexarray(:))%indrydep(1)
    load(:)%indrydep(2) = temp_load(indexarray(:))%indrydep(2)
    load(:)%indrydep(3) = temp_load(indexarray(:))%indrydep(3)
    load(:)%volloc      = temp_load(indexarray(:))%volloc
    load(:)%abstrvol = temp_load(indexarray(:))%abstrvol
    load(:)%abstrind = temp_load(indexarray(:))%abstrind
    DO i = 1,ns
      DO k=1,max_pstype
        load(i)%psvol(k) = temp_load(indexarray(i))%psvol(k)
        DO s=1,numsubstances
          load(i)%psload(k,s) = temp_load(indexarray(i))%psload(k,s)
          IF(k==1) load(i)%locconc(s) = temp_load(indexarray(i))%locconc(s)
        ENDDO
      ENDDO  
    ENDDO
    DO i = 1,nsfile
      DEALLOCATE(temp_load(i)%psvol)
      DEALLOCATE(temp_load(i)%psload)
      DEALLOCATE(temp_load(i)%locconc)
    ENDDO
    DEALLOCATE(temp_load)
    IF(wetlandexist)THEN  
      wetland(:,:)    = temp_wetland(indexarray(:),:)
      IF(SUM(wetland%area)==0)THEN
        wetlandexist = .FALSE.
        DEALLOCATE(wetland)
      ENDIF
    ELSE
      DEALLOCATE(wetland)
    ENDIF
    DEALLOCATE(temp_wetland)
    ! nregions shall not be changed!

    DO i=1,noutreg
      DO isub=1,outregion(i)%nsubbasin
        DO k=1,ns
          IF(indexarray(k)==outregion(i)%subindex(isub))THEN
            outregion(i)%subindex(isub) = k
            EXIT
          ENDIF
        ENDDO
        IF(k>ns)THEN
          WRITE(6,*) 'Warning: outregion',i,'is not (fully) included in submodel.'
          outregion(i)%nsubbasin = 0  !To check against calculation
          outregion(i)%subindex = 0 !To find error
          EXIT
        ENDIF
      ENDDO
    ENDDO

    !Reform observation indexarray for submodel
    DO iforc = 1, max_forcingdata
      IF(ALLOCATED(forcingdata(iforc)%basinindex))THEN
        IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
        tmp = forcingdata(iforc)%basinindex(indexarray(:))
        DEALLOCATE(forcingdata(iforc)%basinindex)
        CALL MOVE_ALLOC(tmp,forcingdata(iforc)%basinindex)
      ENDIF
    ENDDO
    IF(ALLOCATED(qobsindex))THEN
      IF(.NOT.ALLOCATED(tmp)) ALLOCATE(tmp(ns))
      tmp = qobsindex(indexarray(:))
      DEALLOCATE(qobsindex)
      CALL MOVE_ALLOC(tmp,qobsindex)
    ENDIF
    IF(ALLOCATED(xobsindex))THEN 
      IF(.NOT.ALLOCATED(xtmp)) ALLOCATE(xtmp(max_outvar,ns))
      xtmp = xobsindex(:,indexarray(:))
      DEALLOCATE(xobsindex)
      CALL MOVE_ALLOC(xtmp,xobsindex)
    ENDIF

    !Reform subbasin information about forcing
    IF(ALLOCATED(tobselevation))THEN
      ALLOCATE(temp_tobselev(nsfile))
      temp_tobselev = tobselevation
      DEALLOCATE(tobselevation)
      ALLOCATE(tobselevation(ns))
      tobselevation(:)   = temp_tobselev(indexarray(:))
      DEALLOCATE(temp_tobselev)
    ENDIF
    
    !Reform subbasin information from LakeData and DamData
    !Only index for lakes handled the lakedata is kept as is
    IF(ALLOCATED(lakeindex))THEN
      IF(.NOT.ALLOCATED(temp_lakeindex)) ALLOCATE(temp_lakeindex(nsfile))
      temp_lakeindex    = lakeindex
      DEALLOCATE(lakeindex)
      ALLOCATE(lakeindex(ns))
      lakeindex(:)      = temp_lakeindex(indexarray(:))
      DEALLOCATE(temp_lakeindex)
    ENDIF
    IF(ALLOCATED(lakeout2index))THEN
      IF(.NOT.ALLOCATED(temp_lakeindex)) ALLOCATE(temp_lakeindex(nsfile))
      temp_lakeindex    = lakeout2index
      DEALLOCATE(lakeout2index)
      ALLOCATE(lakeout2index(ns))
      lakeout2index(:)      = temp_lakeindex(indexarray(:))
      DEALLOCATE(temp_lakeindex)
    ENDIF

    IF(ALLOCATED(lakebasinindex))THEN
      IF(.NOT.ALLOCATED(temp_lakebasinindex)) ALLOCATE(temp_lakebasinindex(nsfile))
      temp_lakebasinindex = lakebasinindex
      DEALLOCATE(lakebasinindex)
      ALLOCATE(lakebasinindex(ns))
      lakebasinindex(:) = temp_lakebasinindex(indexarray(:))
      DEALLOCATE(temp_lakebasinindex)
    ENDIF

    IF(ALLOCATED(lakedataparindex))THEN
      IF(.NOT.ALLOCATED(temp_lakedataparindex)) ALLOCATE(temp_lakedataparindex(nsfile,2))
      temp_lakedataparindex = lakedataparindex
      DEALLOCATE(lakedataparindex)
      ALLOCATE(lakedataparindex(ns,2))
      lakedataparindex(:,:) = temp_lakedataparindex(indexarray(:),:)
      DEALLOCATE(temp_lakedataparindex)
    ENDIF

    IF(ALLOCATED(damindex))THEN
      IF(.NOT.ALLOCATED(temp_lakeindex)) ALLOCATE(temp_lakeindex(nsfile))
      temp_lakeindex    = damindex
      DEALLOCATE(damindex)
      ALLOCATE(damindex(ns))
      damindex(:)      = temp_lakeindex(indexarray(:))
      DEALLOCATE(temp_lakeindex)
    ENDIF

    IF(glacierexist)THEN
      IF(ALLOCATED(glacierindex))THEN
        IF(.NOT.ALLOCATED(temp_glacierindex)) ALLOCATE(temp_glacierindex(nsfile))
        temp_glacierindex    = glacierindex
        DEALLOCATE(glacierindex)
        ALLOCATE(glacierindex(ns))
        glacierindex(:)      = temp_glacierindex(indexarray(:))
        DEALLOCATE(temp_glacierindex)
      ENDIF
    ENDIF

    IF(ALLOCATED(floodindex))THEN
      IF(.NOT.ALLOCATED(temp_lakeindex)) ALLOCATE(temp_lakeindex(nsfile))
      temp_lakeindex    = floodindex
      DEALLOCATE(floodindex)
      ALLOCATE(floodindex(ns))
      floodindex(:)      = temp_lakeindex(indexarray(:))
      DEALLOCATE(temp_lakeindex)
    ENDIF

    !Model parameters reformation, only subbasin dependent parameters need to be changed
    nbasinpar = 0
    DO i = 1,max_par
      IF(modparid(i)%deptype==m_bpar)  nbasinpar = MAX(nbasinpar,modparid(i)%parno)
    ENDDO
    IF(.NOT.ALLOCATED(temp_basinpar)) ALLOCATE(temp_basinpar(nbasinpar,nsfile))
    temp_basinpar = basinpar
    IF(ALLOCATED(basinpar)) DEALLOCATE(basinpar)
    ALLOCATE(basinpar(nbasinpar,ns))
    basinpar = temp_basinpar(:,indexarray)
    IF(ALLOCATED(temp_basinpar)) DEALLOCATE(temp_basinpar)

  END SUBROUTINE reform_inputdata_for_submodel
  
  !>Finds the corresponding index for a subid
  !--------------------------------------------
  INTEGER FUNCTION get_subid_index(id,ns)

    USE MODVAR, ONLY : basin

    !Argument declarations
    INTEGER, INTENT(IN)  :: id    !<subid
    INTEGER, INTENT(IN)  :: ns    !<number of subbasins
    
    !Local variables
    INTEGER j

    !>\b Algorithm \n
    !>Loop through subbasin id: to find the index of the current one.
    get_subid_index = 0
    DO j=1,ns
      IF(basin(j)%subid==id)THEN
        get_subid_index = j
        EXIT
      ENDIF
    ENDDO

  END FUNCTION get_subid_index

  !>Reads information about crops from file
  !>
  !>\b Consequences Module modvar variables cropdata, cropirrdata 
  !>and cropdataindex is allocated and set.
  !----------------------------------------------------------
  SUBROUTINE load_cropdata(dir,infile,n,status) 

    USE WORLDVAR, ONLY : fileunit_temp, &
                         i_str,i_intg,i_real
    USE MODVAR, ONLY : cropindex,      &  !OUT
                       cropdata,       &  !OUT
                       cropirrdata,    &  !OUT
                       numsubstances,  &
                       modeloption, &
                       p_erosion, &
                       set_cropdataindex, &
                       initiate_cropdata, &
                       i_in,i_on,i_sp,i_pp,i_t1,i_ss
    
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir     !<File directory
    CHARACTER (LEN=*), INTENT(IN) :: infile  !<Name of file to be read (CropData.txt)
    INTEGER, INTENT(OUT) :: n                !<Number of rows in file
    INTEGER, INTENT(OUT) :: status           !<Error status
    
    !Local parameters
    INTEGER, PARAMETER :: nskip = 1   !heading on row 1
    
    !Local variables
    INTEGER i
    INTEGER mcols                     !Actual number of coulmns in file
    INTEGER sumirr
    INTEGER ncols                     !Number of columns in file
    INTEGER, ALLOCATABLE :: code(:)   !Data type code
    INTEGER, ALLOCATABLE :: iindex(:) !Index integer data from file
    INTEGER, ALLOCATABLE :: rindex(:) !Index real data from file
    INTEGER, ALLOCATABLE :: xi(:,:)   !integer data from file
    INTEGER, ALLOCATABLE :: cropregion(:) !CropData region read from file
    INTEGER, ALLOCATABLE :: cropid(:) !CropData cropid read from file
    REAL, ALLOCATABLE    :: xr(:,:)   !real data from file
    CHARACTER (LEN=10),ALLOCATABLE :: colstr(:)   !Content string of data from file
    LOGICAL calcNP
    LOGICAL fileexist

    !>\b Algorithm \n
    status = 0

    !NP simulation?
    calcNP = .FALSE.
    DO i = 1, numsubstances
      IF(i==i_in .OR. i==i_on .OR. i==i_sp .OR. i==i_pp) calcNP = .TRUE.
    ENDDO

    !>Check if file exist and if it is required
    INQUIRE(FILE=TRIM(dir)//infile,EXIST=fileexist)
    IF(.NOT.fileexist)THEN
      IF(calcNP)THEN
        !N / P modelling require the file
        WRITE(6,*) ' ERROR: Missing CropData.txt for NP simulation'  
        status = 1
      ENDIF
      IF(i_ss>0.AND. modeloption(p_erosion)==0)THEN
        !S modelling with default erosion require the file
        WRITE(6,*) ' ERROR: Missing CropData.txt for S simulation'  
        status = 1
      ENDIF
      RETURN
    ENDIF

    !>Count number of crops and data columns
    CALL count_data_rows(fileunit_temp,TRIM(dir)//infile,nskip,n,status)
    IF(status/=0)RETURN
    CALL count_data_cols(fileunit_temp,TRIM(dir)//infile,nskip,ncols,status)
    IF(status/=0)RETURN
    IF(.NOT.ALLOCATED(colstr)) ALLOCATE(colstr(ncols))
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(n,ncols))
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(n,ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))

    !>Allocate and initiate cropdata variables
    IF(.NOT.ALLOCATED(cropregion)) ALLOCATE(cropregion(n))   
    IF(.NOT.ALLOCATED(cropid)) ALLOCATE(cropid(n))   
    IF(.NOT.ALLOCATED(cropdata)) ALLOCATE(cropdata(n))   
    IF(.NOT.ALLOCATED(cropirrdata)) ALLOCATE(cropirrdata(n))
    CALL initiate_cropdata()

    !>Read column headings from file
    OPEN(UNIT=fileunit_temp,FILE=TRIM(dir)//infile,STATUS='old',ACTION='read')
    CALL read_column_headings(fileunit_temp,ncols,colstr,mcols,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(dir)//infile
      RETURN
    ENDIF

    !>Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,mcols
       IF(colstr(i)(1:10)=='cropid    ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='reg       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='fn1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fn2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fp1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fp2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mn1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mn2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mp1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mp2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fday1     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='fday2     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='mday1     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='mday2     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='fdown1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='fdown2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mdown1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='mdown2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resn      ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resp      ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resc      ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resday    ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='resdown   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='resfast   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='up1       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='up2       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='up3       ')   code(i) = i_real
       IF(colstr(i)(1:10)=='bd1       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd2       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd3       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd4       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='bd5       ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='ccmax1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='ccmax2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='gcmax1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='gcmax2    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='kcbini    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='kcbmid    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='kcbend    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='dlref     ')   code(i) = i_real
       IF(colstr(i)(1:10)=='plantday  ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthini ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthdev ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthmid ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='lengthlate')   code(i) = i_intg
       IF(colstr(i)(1:10)=='imm_start ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='imm_end   ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='upupper   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='pnupr     ')   code(i) = i_real
       IF(colstr(i)(1:10)=='daylength ')   code(i) = i_real
       IF(colstr(i)(1:10)=='gddsow    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='basetemp  ')   code(i) = i_real
       IF(colstr(i)(1:10)=='firstday  ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='tamount   ')   code(i) = i_real
       IF(colstr(i)(1:10)=='tyear     ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='tday      ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='tnumdays  ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='tdaydown  ')   code(i) = i_intg
       IF(colstr(i)(1:10)=='tdown1    ')   code(i) = i_real
       IF(colstr(i)(1:10)=='tdown2    ')   code(i) = i_real
    ENDDO

    !>Read the data of the crops to matrix, first column is string which is skipped
    CALL read_basindata5(fileunit_temp,infile,ncols,n,ncols,code,rindex,iindex,xi,xr) 
    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(dir)//infile

    !>Set cropdata variables according to matrix with read data
    DO i = 1,mcols
      IF(colstr(i)(1:10)=='cropid    ')   cropid(1:n)                = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='reg       ')   cropregion(1:n)            = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='fn1       ')   cropdata(1:n)%fertnamount1 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fn2       ')   cropdata(1:n)%fertnamount2 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fp1       ')   cropdata(1:n)%fertpamount1 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fp2       ')   cropdata(1:n)%fertpamount2 = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mn1       ')   cropdata(1:n)%mannamount1  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mn2       ')   cropdata(1:n)%mannamount2  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mp1       ')   cropdata(1:n)%manpamount1  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='mp2       ')   cropdata(1:n)%manpamount2  = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='fday1     ')   cropdata(1:n)%fertday1     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='fday2     ')   cropdata(1:n)%fertday2     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='mday1     ')   cropdata(1:n)%manday1      = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='mday2     ')   cropdata(1:n)%manday2      = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='fdown1    ')   cropdata(1:n)%fertdown1    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='fdown2    ')   cropdata(1:n)%fertdown2    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='mdown1    ')   cropdata(1:n)%mandown1     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='mdown2    ')   cropdata(1:n)%mandown2     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='resn      ')   cropdata(1:n)%resnamount   = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='resp      ')   cropdata(1:n)%respamount   = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='resc      ')   cropdata(1:n)%rescamount   = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='resday    ')   cropdata(1:n)%resdayno     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='resdown   ')   cropdata(1:n)%resdown      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='resfast   ')   cropdata(1:n)%resfast      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='up1       ')   cropdata(1:n)%uptake1      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='up2       ')   cropdata(1:n)%uptake2      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='up3       ')   cropdata(1:n)%uptake3      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='bd1       ')   cropdata(1:n)%baredayno1   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd2       ')   cropdata(1:n)%baredayno2   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd3       ')   cropdata(1:n)%baredayno3   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd4       ')   cropdata(1:n)%baredayno4   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='bd5       ')   cropdata(1:n)%baredayno5   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='ccmax1    ')   cropdata(1:n)%ccmax1       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='ccmax2    ')   cropdata(1:n)%ccmax2       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='gcmax1    ')   cropdata(1:n)%gcmax1       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='gcmax2    ')   cropdata(1:n)%gcmax2       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='kcbini    ')   cropirrdata(1:n)%kcbini    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='kcbmid    ')   cropirrdata(1:n)%kcbmid    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='kcbend    ')   cropirrdata(1:n)%kcbend    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='dlref     ')   cropirrdata(1:n)%dlref     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='plantday  ')   cropirrdata(1:n)%plantingdayno = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthini ')   cropirrdata(1:n)%lengthini     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthdev ')   cropirrdata(1:n)%lengthdev     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthmid ')   cropirrdata(1:n)%lengthmid     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='lengthlate')   cropirrdata(1:n)%lengthlate    = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='imm_start ')   cropirrdata(1:n)%imm_start     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='imm_end   ')   cropirrdata(1:n)%imm_end       = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='upupper   ')   cropdata(1:n)%uptakeupper      = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='pnupr     ')   cropdata(1:n)%PNuptakeRatio    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='daylength ')   cropdata(1:n)%daylength    = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='gddsow    ')   cropdata(1:n)%gddsow       = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='basetemp  ')   cropdata(1:n)%basetemp     = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='firstday  ')   cropdata(1:n)%firstday     = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='tamount   ')   cropdata(1:n)%T1amount         = xr(1:n,rindex(i))*100.   !kg/ha->kg/km2
      IF(colstr(i)(1:10)=='tyear     ')   cropdata(1:n)%T1year           = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='tday      ')   cropdata(1:n)%T1day            = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='tnumdays  ')   cropdata(1:n)%T1numberofdays   = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='tdaydown  ')   cropdata(1:n)%T1daydown        = xi(1:n,iindex(i))
      IF(colstr(i)(1:10)=='tdown1    ')   cropdata(1:n)%T1down1          = xr(1:n,rindex(i))
      IF(colstr(i)(1:10)=='tdown2    ')   cropdata(1:n)%T1down2          = xr(1:n,rindex(i))    
    ENDDO

    !>Check if cropdata is needed?
    IF(.NOT.(calcNP .OR. i_ss>0 .OR. i_t1>0)) DEALLOCATE(cropdata)
    sumirr = SUM(cropirrdata(1:n)%plantingdayno)
    IF(sumirr==0) DEALLOCATE(cropirrdata)

    IF(calcNP .OR. i_ss>0 .OR. (i_t1>0) .OR. (sumirr>0))THEN
      !>Calculate index variable to find crops
      CALL set_cropdataindex(n,cropid,cropregion,cropindex)
    ENDIF

    !Deallocate local arrays
    DEALLOCATE(xr,xi,iindex,rindex,code,colstr,cropid,cropregion)
    
  END SUBROUTINE load_cropdata
  
  !>Get command line argument: infodir and simsequence, 
  !>or read them from filedir.txt (old variant).
  !-----------------------------------------------------------------------
  SUBROUTINE get_hyss_arguments(dir,iseq)
    
    USE WORLDVAR, ONLY : fileunit_temp,   &
         maxcharpath

    !Argument declarations
    CHARACTER(LEN=maxcharpath), INTENT(OUT) :: dir  !<Directory for information about simulation (infodir)
    INTEGER, INTENT(OUT) :: iseq                    !<Sequence number (simsequence)
    
    !Local variables
    INTEGER i       !Argument number
    INTEGER narg    !Number of command line argument
    INTEGER st      !Status
    INTEGER pos,oldpos  !Position on line
    LOGICAL nostr
    CHARACTER(LEN=maxcharpath) argument,argi  !Command line argument i
    CHARACTER(LEN=500) line

    !>\b Algorithm \n
    iseq = 0
    pos = 1
    dir=''
    narg = COMMAND_ARGUMENT_COUNT()
    !>If no command line argument, read information from filedir.txt
    IF(narg==0)THEN
      OPEN(fileunit_temp,FILE = 'filedir.txt',STATUS = 'old',ACTION='read',ERR=200)
      READ(fileunit_temp,*,END=300,ERR=300) line
      DO
        nostr = .FALSE.
        oldpos = pos
        CALL read_next_codestr_on_line(500,maxcharpath,pos,line,argument,nostr,'err',.TRUE.)
        IF(nostr) EXIT
        argi=TRIM(argument)
        SELECT CASE(argi)
        CASE('-sequence','-s')
          CALL read_next_codestr_on_line(500,maxcharpath,pos,line,argument,nostr,'err')
          READ(argument,*,END=101,ERR=101) iseq
        CASE('-infodir','-i')
          CALL read_next_codestr_on_line(500,maxcharpath,pos,line,argument,nostr,'err',.TRUE.)
          CALL GET_COMMAND_ARGUMENT(i,argument,STATUS=st)
          IF(.NOT.nostr)THEN
            dir=TRIM(argument)
          ELSE
            WRITE(*,*) 'ERROR reading infodir from filedir.txt'
          ENDIF
        CASE DEFAULT
          IF(oldpos==1)THEN
            dir=argi
          ELSE
            WRITE(*,*) 'ERROR: No file path found in filedir.txt.'
            STOP 1
          ENDIF
        END SELECT
      END DO
       CLOSE(fileunit_temp)
       RETURN
    ENDIF
    
    !>Get command line argument
    i = 1
    DO WHILE(i<=narg)
       CALL GET_COMMAND_ARGUMENT(i,argument)
       argi=TRIM(argument)
       SELECT CASE(argi)
       CASE('-sequence','-s')
          i = i + 1
          CALL GET_COMMAND_ARGUMENT(i,argument)
          READ(argument,*,END=100,ERR=100) iseq
       CASE('-infodir','-i')
          i = i + 1
          CALL GET_COMMAND_ARGUMENT(i,argument,STATUS=st)
          IF(st==0)THEN
             dir=TRIM(argument)
          ELSE
             WRITE(*,*) 'ERROR reading infodir command line argument'
          ENDIF
       CASE DEFAULT
          IF(i==1)THEN
             dir=argi
          ELSE
             WRITE(*,*) 'ERROR in command line argument'
             STOP 1
          ENDIF
       END SELECT
       i = i + 1
    END DO
    RETURN

    !Error handling
100 WRITE(*,*) 'ERROR reading sequence command line argument'      
    STOP 1
101 WRITE(*,*) 'ERROR reading sequence from filedir.txt'      
    STOP 1
200 WRITE(*,*) 'ERROR no command line argument and no filedir.txt found' 
    STOP 1
300 WRITE(*,*) 'ERROR reading filedir.txt' 
    STOP 1

  END SUBROUTINE get_hyss_arguments

  !\brief Reads file with information about the current simulation 
  !!
  !!Format of file: first on line is a coded string followed by the value 
  !!The string "!!" is used for comment rows and are not read
  !!
  !>\b Consequences The following worldvar module variables is set: readdaily,
  !> readobsid, outvarallbasin, writeload, readmatlab, writematlab,
  !> ncrit, readformat, simsubmodel, resultseq, steplen, 
  !> forcingdata and checkdata. The following modvar 
  !> module variables is set: i_in,i_on,i_sp,i_pp,i_t1,i_t2,i_oc, 
  !> conductN,conductP,conductC,conductT,conductload,conductregest and modeloption.
  !--------------------------------------------------------------------
  SUBROUTINE load_coded_info(dir,status,date1,date2,date3,skip,step,nsubst,     &
                             stateinput,somax,sonum,bstateoutdate,pstateoutdate,  &
                             mdir,rdir,                           &
                             locupall,locupnone,locupnone_qar,locupnone_war,wupall,wupnone,wupdname,noutvar) 

    USE WORLDVAR, ONLY : fileunit_temp, &
                         maxcharpath, &
                         maxoutbasins, &
                         maxcrit, &
                         max_typeofperiods, &
                         comment_str, &
                         set_calibration, &
                         set_output, &
                         noutput, &   !OUT
                         set_maxmap, &
                         readdaily, &   !OUT
                         readobsid, &   !OUT
                         outvarallbasin, &   !OUT
                         writeload, &   !OUT
                         readmatlab, &   !OUT
                         writematlab, &   !OUT
                         ncrit, &   !OUT
                         readformat, &   !OUT
                         readoutregion, &   !OUT
                         simsubmodel, &   !OUT
                         simsequence, &
                         resultseq, &   !OUT
                         i_d, &
                         steplen, &   !(OUT)
                         forcingdata, & !OUT
                         max_forcingdata, &
                         allocate_and_initialize_forcingdata_structure, &
                         i_sfobs,i_swobs, & 
                         i_tminobs,i_tmaxobs, & 
                         checkdata, &   !OUT
                         doassimilation       !OUT
    USE MODVAR, ONLY : i_in,i_on,i_sp,i_pp, &   !OUT
                       i_t1,i_t2,i_oc,i_ss,i_ae, &   !OUT
                       doupdate, &
                       conductN,conductP,conductC,conductT,conductS, & !OUT
                       conductxoms, & !OUT
                       conductwb, & !OUT
                       conductload, &  !OUT
                       conductregest, &  !OUT
                       soiliniwet, &  !OUT
                       soildepthstretch, &  !OUT
                       max_outvar, &
                       i_quseobs, &
                       i_qar,i_war, &
                       i_wendupd, &
                       i_tpcorr, &
                       i_tncorr, &
                       i_tploccorr, &
                       i_tnloccorr, &
                       irrunlimited, &
                       allocate_initialize_modeloptions, &
                       modeloptionname, &
                       modeloption, &    !OUT
                       num_modelprocess, &   
                       p_snowfall, &
                       p_snowmelt, &
                       p_snowevap, &
                       p_snowdensity, &
                       p_lakeriverice, & 
                       p_petmodel, &
                       p_deepgroundwater, &
                       p_swtemperature, &
                       p_growthstart, &
                       p_floodplain, &
                       p_glacierini, &
                       p_infiltration, &
                       p_erosion
  USE READWRITE_ROUTINES, ONLY : convert_string_to_integer, &
                                 read_next_date_on_line, &
                                 print_output_information_to_logfile

    !Argument declarations
    CHARACTER(LEN=maxcharpath), INTENT(IN) :: dir   !<File directory
    INTEGER, INTENT(OUT) :: status                  !<Status of subroutine
    TYPE(DateType), INTENT(OUT) :: date1            !<Begin date of simulation
    TYPE(DateType), INTENT(OUT) :: date2            !<End date of simulation
    TYPE(DateType), INTENT(OUT) :: date3            !<Begin date for criteria/output date
    INTEGER, INTENT(OUT) :: skip                    !<Length of warm-up period
    INTEGER, INTENT(OUT) :: step                    !<Length of simulation period
    INTEGER, INTENT(OUT) :: nsubst                  !<Number of substances simulated
    LOGICAL, INTENT(OUT) :: stateinput              !<Code for reading input state
    INTEGER, INTENT(IN)  :: somax                   !<Dimension dates for output state
    INTEGER, INTENT(OUT) :: sonum                   !<Number of dates for output state
    TYPE(DateType), INTENT(OUT) :: bstateoutdate(somax) !<Dates for writing state variables
    TYPE(DateType), INTENT(OUT) :: pstateoutdate(somax) !<Dates previous to date for writing 
    CHARACTER(LEN=maxcharpath), INTENT(OUT) :: mdir           !<File directory model
    CHARACTER(LEN=maxcharpath), INTENT(OUT) :: rdir           !<File directory result
    LOGICAL, INTENT(OUT) :: locupall                !<Flag for Q-update for all stations
    LOGICAL, INTENT(OUT) :: locupnone               !<Flag for Q-update on no stations
    LOGICAL, INTENT(OUT) :: locupnone_qar           !<Flag for Q-update on no stations
    LOGICAL, INTENT(OUT) :: locupnone_war           !<Flag for WAR-update on no stations
    LOGICAL, INTENT(OUT) :: wupall                  !<Flag for W-update all stations
    LOGICAL, INTENT(OUT) :: wupnone                 !<Flag for W-update nostations
    CHARACTER(LEN=4), INTENT(OUT) :: wupdname       !<Variable name for wendupd-function
    INTEGER, INTENT(OUT) :: noutvar      !<number of used output variables to be set (so far)
    
    !Local parameters
    INTEGER, PARAMETER :: linelen = 18000      
    CHARACTER (LEN=2), PARAMETER :: errstr = '##'  
    
    !Local variables
    LOGICAL error
    LOGICAL onoff
    LOGICAL calibration
    LOGICAL nostrfound,notimefound
    LOGICAL exitbigloop
    LOGICAL cyclebigloop
    INTEGER io,i,iforc
    INTEGER temp
    INTEGER maxdim
    INTEGER intvalue
    CHARACTER(LEN=3), ALLOCATABLE :: c(:)
    INTEGER, ALLOCATABLE :: critc(:)
    INTEGER, ALLOCATABLE :: critr(:)
    INTEGER, ALLOCATABLE :: flowc(:)
    INTEGER, ALLOCATABLE :: flowr(:)
    INTEGER, ALLOCATABLE :: critareac(:)
    INTEGER, ALLOCATABLE :: critarear(:)
    REAL, ALLOCATABLE :: critweight(:)
    REAL, ALLOCATABLE :: critpar(:)
    INTEGER istate
    INTEGER ioutbasin(max_typeofperiods)
    INTEGER readbasins(max_typeofperiods,maxoutbasins)
    INTEGER ioutregion(max_typeofperiods)
    INTEGER readregions(max_typeofperiods,maxoutbasins)
    INTEGER ibasinoutputvar(max_typeofperiods)
    INTEGER iregionoutputvar(max_typeofperiods)
    INTEGER imapoutputvar(max_typeofperiods)
    INTEGER itimeoutputvar(max_typeofperiods)
    INTEGER ivarindex
    INTEGER varindex,flowtype
    INTEGER basinoutputperiod(max_typeofperiods)
    INTEGER regionoutputperiod(max_typeofperiods)
    INTEGER mapoutputperiod(max_typeofperiods)
    INTEGER timeoutputperiod(max_typeofperiods)
    INTEGER basinoutputdecimal(max_typeofperiods)
    INTEGER regionoutputdecimal(max_typeofperiods)
    INTEGER mapoutputdecimal(max_typeofperiods)
    INTEGER timeoutputdecimal(max_typeofperiods)
    INTEGER basinoutputsignif(max_typeofperiods)
    INTEGER regionoutputsignif(max_typeofperiods)
    INTEGER mapoutputsignif(max_typeofperiods)
    INTEGER timeoutputsignif(max_typeofperiods)
    INTEGER basinoutputvar(max_typeofperiods,max_outvar,2)
    INTEGER regionoutputvar(max_typeofperiods,max_outvar)
    INTEGER mapoutputvar(max_typeofperiods,max_outvar,2)
    INTEGER timeoutputvar(max_typeofperiods,max_outvar,2)
    INTEGER icrit,critperiod,critlimit
    INTEGER tstep,indx !Time step length [stepunit] and string index (used with stepstr)
    INTEGER linepos
    INTEGER icheck,areagg
    CHARACTER (LEN=linelen) :: line
    CHARACTER (LEN=16) :: code,code2, code3
    CHARACTER (LEN=10) :: strvalue
    CHARACTER (LEN=16) :: strdate,strdate1,strdate2,strdate3
    CHARACTER (LEN=16), ALLOCATABLE :: strstates(:)
    CHARACTER (LEN=6)  :: varstr
    CHARACTER (LEN=3)  :: yesnostr
    CHARACTER (LEN=3)  :: critstr
    CHARACTER (LEN=maxcharpath) :: filedir
    CHARACTER (LEN=maxcharpath+8) :: filename
    CHARACTER (LEN=4)  :: stepunit      !Time step unit; mo,d,h or min
    CHARACTER (LEN=10) :: stepstr       !Temporary storage for time step string
    LOGICAL, ALLOCATABLE :: critcond(:) !flag that this criteria should only be used conditionally /DG
    REAL, ALLOCATABLE :: critthres(:)   !acceptance criteria for conditional criteria /DG
    
    status = 0
    error = .FALSE.
    ncrit = maxcrit
    ALLOCATE(c(ncrit),critweight(ncrit),critc(ncrit),flowc(ncrit),    &
             critr(ncrit),flowr(ncrit),critpar(ncrit), &
             critcond(ncrit),critthres(ncrit), &
             critareac(ncrit),critarear(ncrit))
    ncrit = 0
    !initialize the conditional criteria and threshold
    critcond = .FALSE.
    critthres = 0.

    !Open file      
    filename=TRIM(dir)//'info.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')

    !Default values
    calibration    = .FALSE.
    doassimilation = .FALSE.
    stateinput     = .FALSE.
    writematlab    = .FALSE.
    readmatlab     = .FALSE.
    outvarallbasin = .FALSE.
    exitbigloop    = .FALSE.
    simsubmodel    = .FALSE.
    doupdate       = .FALSE.
    locupall       = .FALSE.
    locupnone      = .FALSE.
    locupnone_qar  = .FALSE.
    locupnone_war  = .FALSE.
    wupall         = .FALSE.
    wupnone        = .FALSE.
    readdaily      = .FALSE.
    writeload      = .FALSE.
    readobsid      = .TRUE. 
    irrunlimited   = .FALSE.                !Unlimited irrigation, defaults to False
    soiliniwet     = .FALSE.
    soildepthstretch = .FALSE.
    conductregest  = .FALSE.
    resultseq      = .FALSE.  
    CALL allocate_and_initialize_forcingdata_structure()
    readoutregion  = .FALSE.
    conductxoms    = .FALSE.
    conductwb      = .FALSE.
    CALL allocate_initialize_modeloptions()
    checkdata      = .FALSE.
    conductN       = .FALSE.
    conductP       = .FALSE.
    conductC       = .FALSE.
    conductT       = .FALSE.
    conductS       = .FALSE.
    conductload    = .FALSE.
    IF(simsequence>0) resultseq = .TRUE.  
    date3          = DateType(0,0,0,0,0)
    bstateoutdate  = DateType(0,0,0,0,0)
    ncrit          = 0
    readformat     = 0         !ascii
    nsubst         = 0
    i_in=0;i_on=0;i_sp=0;i_pp=0;i_t1=0;i_t2=0;i_oc=0;i_ss=0;i_ae=0
    mdir           = dir
    rdir           = dir
    ioutbasin          = 0
    ioutregion         = 0
    ibasinoutputvar    = 0
    iregionoutputvar   = 0
    imapoutputvar      = 0
    itimeoutputvar     = 0
    basinoutputvar     = 0
    regionoutputvar    = 0
    mapoutputvar       = 0
    timeoutputvar      = 0
    basinoutputperiod  = 1   !daily
    regionoutputperiod = 1   !daily
    mapoutputperiod    = 5   !mean over simulation period
    timeoutputperiod   = 1   !daily
    basinoutputdecimal = 1
    regionoutputdecimal = 1
    mapoutputdecimal   = 1
    timeoutputdecimal  = 3
    basinoutputsignif  = 0
    regionoutputsignif = 0
    mapoutputsignif    = 0
    timeoutputsignif   = 0
    critperiod         = i_d  !for calvarper
    critlimit          = 3    !for calvarlim
    wupdname=''
    sonum              = 0
    tstep              = 1    ! Default time step length is 1 ..
    stepunit           = 'd'  ! .. day

    !Read content of file      
    DO
      READ(fileunit_temp,600,END=200,ERR=800,IOSTAT=io) line
      linepos = 1
      CALL read_next_codestr_on_line(linelen,16,linepos,line,code,nostrfound,errstr)    !Read code
      IF(code(1:2)==errstr)THEN
        error = .TRUE.
        code=line(1:16)
        EXIT
      ENDIF
      IF(code(1:2)==comment_str) CYCLE
      IF(nostrfound) CYCLE  !empty row
      cyclebigloop = .FALSE.

      !Read date   
      IF(code(1:5)=='bdate' .OR. &
         code(1:5)=='cdate' .OR. &
         code(1:5)=='edate' .OR. &
         code(1:12)=='outstatedate')THEN
        CALL read_next_date_on_line(linelen,16,linepos,line,strdate,nostrfound,notimefound,errstr)
        IF(nostrfound.OR.strdate(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        SELECT CASE(code)
        CASE ('bdate')
          CALL string_convert_to_DateType(strdate,date1)
        CASE ('cdate')
          CALL string_convert_to_DateType(strdate,date3)
        CASE ('edate')
          CALL string_convert_to_DateType(strdate,date2)
        CASE ('outstatedate')
          DO
            sonum = sonum + 1
            IF(sonum>somax) THEN
              WRITE(6,*) 'ERROR: Too many outstatedate i info.txt'
              status = 1
              RETURN
            ENDIF
            CALL string_convert_to_DateType(strdate,bstateoutdate(sonum))
            CALL read_next_date_on_line(linelen,16,linepos,line,strdate,nostrfound,notimefound,errstr)
            IF(nostrfound.OR.strdate(1:2)==errstr)THEN
              EXIT
            ENDIF
          ENDDO             
        END SELECT
      ENDIF

      !Read time step, if present
      IF(code(1:10)=='steplength')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,stepstr,nostrfound,errstr,.TRUE.)
        IF(nostrfound)THEN   ! No error here if not defined, only using default values
          WRITE(6,'(A60,I2,A3)') &
                  'WARNING: No value for steplength found, using default values',tstep,stepunit
        ELSE
          indx=SCAN(stepstr,'dhm') ! Looking for unit: d,h or min
          IF(indx.GT.0)THEN
            READ(stepstr(1:indx-1),'(I4)') tstep
            stepunit=stepstr(indx:LEN(stepstr))
          ELSE
            WRITE(6,'(A65,I2,A3)') &
                   'WARNING: No valid unit for steplength found, using default values',tstep,stepunit
          ENDIF
        ENDIF
      ENDIF

      !Read yes or no coded variables
      DO iforc = 1,max_forcingdata
        IF(code(1:12)==forcingdata(iforc)%infocode) status = set_yesno_variable_from_line(linelen,linepos,line,errstr,forcingdata(iforc)%readfile)
        IF(status/=0)THEN
          error = .TRUE.
          EXIT  
        ENDIF
      ENDDO
      IF(code(1:12)=='assimilation'.OR. &
         code(1:11)=='calibration' .OR. &
         code(1:9) =='readdaily'   .OR. &
         code(1:9) =='printload'   .OR. &
         code(1:8) =='submodel'    .OR. &
         code(1:7) =='instate'     .OR. &
         code(1:8) =='resseqnr'    .OR. &
         code(1:12)=='irrunlimited'.OR. &
         code(1:9) =='readobsid'   .OR. &
         code(1:13)=='printwaterbal' .OR. &
         code(1:10)=='soiliniwet' .OR. &
         code(1:11)=='soilstretch' .OR. &
         code(1:11)=='regestimate' .OR. &
         code(1:13)=='readoutregion' .OR. &
         code(1:13)=='readxomsfiles')THEN
        CALL read_next_codestr_on_line(linelen,3,linepos,line,yesnostr,nostrfound,errstr)
        IF(nostrfound.OR.yesnostr(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(yesnostr(1:1)=='Y' .OR. yesnostr(1:1)=='y' .OR.    &
           yesnostr(1:1)=='J' .OR. yesnostr(1:1)=='j')THEN
          onoff = .TRUE.
        ELSEIF(yesnostr(1:1)=='N' .OR. yesnostr(1:1)=='n')THEN
          onoff = .FALSE.
        ELSE
          error = .TRUE.
          EXIT  
        ENDIF
        SELECT CASE(code)
        CASE ('assimilation')
          doassimilation = onoff
        CASE ('calibration')
          calibration = onoff
        CASE ('readdaily')
          readdaily = onoff
        CASE ('printload')
          writeload = onoff
        CASE ('submodel')
          simsubmodel = onoff
        CASE ('instate')
          stateinput = onoff
        CASE ('resseqnr')
          resultseq = onoff
        CASE ('readobsid')
          readobsid = onoff
        CASE ('readoutregion')
          readoutregion = onoff
        CASE ('readxomsfiles')
          conductxoms = onoff
        CASE ('printwaterbal')
          conductwb = onoff
        CASE ('irrunlimited')
          irrunlimited = onoff
        CASE ('soiliniwet')
          soiliniwet = onoff
        CASE ('soilstretch')
          soildepthstretch = onoff
        CASE ('regestimate')
          conductregest = onoff
        END SELECT
      ENDIF

      !Read integer coded variables
      IF(code(1:10)=='readformat' .OR. &
         code(1:11)=='writeformat')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        READ(strvalue,*) intvalue
        IF(code(1:10)=='readformat')THEN
          IF(intvalue==1) readmatlab = .TRUE.     
        ELSEIF(code(1:11)=='writeformat')THEN
          IF(intvalue==1) writematlab = .TRUE.
        ENDIF
      ENDIF

      !Read optional model structures coded with integers
      IF(code(1:11)=='modeloption')THEN
        CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        READ(strvalue,*) intvalue
        DO ivarindex = 1,num_modelprocess
          IF(code2==modeloptionname(ivarindex))THEN
            modeloption(ivarindex) = intvalue
            IF(ivarindex==p_lakeriverice)THEN
              !Make sure T2 temperature model is switched on if lake and river ice model is used
              IF(i_t2.LE.0 .AND. modeloption(p_lakeriverice).GE.1)THEN
                nsubst=nsubst+1
                i_t2 = nsubst
                conductT = .TRUE.
              ENDIF
            ENDIF
            EXIT
          ENDIF
        ENDDO
      ENDIF

      !Read substances to be simulated
      IF(code(1:9) =='substance')THEN
        DO
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(strvalue(1:2)==errstr)THEN
            error = .TRUE.
            exitbigloop = .TRUE.
            EXIT
          ENDIF
          IF(nostrfound)THEN
            cyclebigloop = .TRUE.
            EXIT
          ENDIF
          SELECT CASE(strvalue(1:1))
          CASE('N','n')
            i_in = nsubst+1
            i_on = nsubst+2
            nsubst=nsubst+2
            conductN = .TRUE.
          CASE('P','p')
            i_sp = nsubst+1
            i_pp = nsubst+2
            nsubst=nsubst+2
            conductP = .TRUE.
          CASE('C','c')
            nsubst=nsubst+1
            i_oc = nsubst
            conductC = .TRUE.
          CASE('T','t')
            SELECT CASE(strvalue(1:2))
            CASE('T1','t1')
              nsubst=nsubst+1
              i_t1 = nsubst
            CASE('T2','t2')
              !Make sure T2 is only initiated once (see lake and river ice setting above)
              IF(i_t2.LE.0)THEN
                nsubst=nsubst+1
                i_t2 = nsubst
                !Is T2 model relevant without lake and river ice model?
              ENDIF
            END SELECT
            conductT = .TRUE.
          CASE('S','s')
            i_ss = nsubst+1
            i_ae = nsubst+2
            nsubst=nsubst+2
            conductS = .TRUE.
          END SELECT
        ENDDO
        IF(exitbigloop)EXIT
        IF(cyclebigloop)CYCLE
      ENDIF

      !Read directory paths
      IF(code(1:8)=='modeldir' .OR. &
         code(1:9)=='resultdir')THEN
        CALL read_next_codestr_on_line(linelen,maxcharpath,linepos,line,filedir,nostrfound,errstr,.TRUE.)
        IF(nostrfound.OR.filedir(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(code(1:8)=='modeldir')  mdir = filedir
        IF(code(1:9)=='resultdir') rdir = filedir
      ENDIF

      !Read output information
      IF(code(1:11)=='basinoutput' .OR. &
         code(1:9) =='mapoutput'   .OR. &
         code(1:10)=='timeoutput'  .OR. &
         code(1:12)=='regionoutput') THEN
        CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
        IF(code2(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(nostrfound) CYCLE
        icheck = 1
        READ(code2,*,ERR=25) icheck
        CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
        IF(code2(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(nostrfound) CYCLE
25      CONTINUE
        IF(code2(1:8)=='allbasin')THEN
          outvarallbasin = .TRUE.
        ENDIF
        IF(code2(1:10)=='meanperiod')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) intvalue
          IF(code(1:11)=='basinoutput')THEN
            basinoutputperiod(icheck) = intvalue
          ELSEIF(code(1:9)=='mapoutput')THEN
            mapoutputperiod(icheck) = intvalue
          ELSEIF(code(1:10)=='timeoutput')THEN
            timeoutputperiod(icheck) = intvalue
          ELSEIF(code(1:12)=='regionoutput')THEN
            regionoutputperiod(icheck) = intvalue
          ENDIF
        ENDIF
        IF(code2(1:8)=='decimals')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) intvalue
          IF(code(1:11)=='basinoutput')THEN
            basinoutputdecimal(icheck) = intvalue
          ELSEIF(code(1:9)=='mapoutput')THEN
            mapoutputdecimal(icheck) = intvalue
          ELSEIF(code(1:10)=='timeoutput')THEN
            timeoutputdecimal(icheck) = intvalue
          ELSEIF(code(1:12)=='regionoutput')THEN
            regionoutputdecimal(icheck) = intvalue
          ENDIF
        ENDIF
        IF(code2(1:11)=='signfigures')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) intvalue
          IF(code(1:11)=='basinoutput')THEN
            basinoutputsignif(icheck) = intvalue
          ELSEIF(code(1:9)=='mapoutput')THEN
            mapoutputsignif(icheck) = intvalue
          ELSEIF(code(1:10)=='timeoutput')THEN
            timeoutputsignif(icheck) = intvalue
          ELSEIF(code(1:12)=='regionoutput')THEN
            regionoutputsignif(icheck) = intvalue
          ENDIF
        ENDIF
        IF(code2(1:8)=='variable')THEN
          DO
            CALL read_next_codestr_on_line(linelen,6,linepos,line,varstr,nostrfound,errstr)
            IF(varstr(1:2)==errstr)THEN
              error = .TRUE.
              WRITE(6,*)  'ERROR: In variable for ',TRIM(code)
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(nostrfound)THEN
              cyclebigloop = .TRUE.
              EXIT
            ENDIF
            status = find_variable_index_type(varstr,varindex,flowtype,areagg)
            IF(status>0)THEN
              WRITE(6,*)  'ERROR: In variable index calculation: ', varstr,' for ',TRIM(code)
              error = .TRUE.
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(code(1:11)=='basinoutput')THEN
              ibasinoutputvar(icheck) = ibasinoutputvar(icheck) + 1
              IF(ibasinoutputvar(icheck)>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many basinoutput variables, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              basinoutputvar(icheck,ibasinoutputvar(icheck),1) = varindex
              basinoutputvar(icheck,ibasinoutputvar(icheck),2) = areagg
            ELSEIF(code(1:9)=='mapoutput')THEN
              imapoutputvar(icheck) = imapoutputvar(icheck) + 1
              IF(imapoutputvar(icheck)>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many mapoutput variable, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              mapoutputvar(icheck,imapoutputvar(icheck),1) = varindex
              mapoutputvar(icheck,imapoutputvar(icheck),2) = areagg
            ELSEIF(code(1:10)=='timeoutput')THEN
              itimeoutputvar(icheck) = itimeoutputvar(icheck) + 1
              IF(itimeoutputvar(icheck)>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many timeoutput variable, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              timeoutputvar(icheck,itimeoutputvar(icheck),1) = varindex
              timeoutputvar(icheck,itimeoutputvar(icheck),2) = areagg
            ELSEIF(code(1:12)=='regionoutput')THEN
              iregionoutputvar(icheck) = iregionoutputvar(icheck) + 1
              IF(iregionoutputvar(icheck)>max_outvar)THEN
                WRITE(6,*) 'ERROR: To many regionoutput variables, max:',max_outvar
                error = .TRUE.
                exitbigloop = .TRUE.
                EXIT
              ENDIF
              regionoutputvar(icheck,iregionoutputvar(icheck)) = varindex
            ENDIF
          ENDDO
          IF(exitbigloop)EXIT
          IF(cyclebigloop)CYCLE
        ENDIF
        IF(code2(1:8)=='subbasin')THEN
          DO
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(strvalue(1:2)==errstr)THEN
              error = .TRUE.
              WRITE(6,*)  'ERROR: In subbasin for ',TRIM(code)
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(nostrfound)THEN
              cyclebigloop = .TRUE.
              EXIT
            ENDIF
            CALL convert_string_to_integer(10,strvalue,intvalue,status) 
            IF(status/=0)THEN
              error = .TRUE.
              exitbigloop = .TRUE.
            ENDIF  
            ioutbasin(icheck) = ioutbasin(icheck) + 1
            IF(ioutbasin(icheck)>maxoutbasins)THEN
              WRITE(6,*) 'Warning: Number of "basinoutput subbasin" is larger than max',maxoutbasins
              WRITE(6,*) 'Warning: Excess subbasins will be ignored'
            ELSE
              readbasins(icheck,ioutbasin(icheck)) = intvalue
            ENDIF
          ENDDO
          IF(exitbigloop)EXIT
          IF(cyclebigloop)CYCLE
        ENDIF
        IF(code2(1:9)=='outregion')THEN
          DO
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(strvalue(1:2)==errstr)THEN
              error = .TRUE.
              WRITE(6,*)  'ERROR: In outregion for ',TRIM(code)
              exitbigloop = .TRUE.
              EXIT
            ENDIF
            IF(nostrfound)THEN
              cyclebigloop = .TRUE.
              EXIT
            ENDIF
            CALL convert_string_to_integer(10,strvalue,intvalue,status) 
            IF(status/=0)THEN
              error = .TRUE.
              exitbigloop = .TRUE.
            ENDIF  
            ioutregion(icheck) = ioutregion(icheck) + 1
            IF(ioutregion(icheck)>maxoutbasins)THEN
              WRITE(6,*) 'Warning: Number of "regionoutput outregion" is larger than max',maxoutbasins
              WRITE(6,*) 'Warning: Excess output regions will be ignored'
            ELSE
              readregions(icheck,ioutregion(icheck)) = intvalue
            ENDIF
          ENDDO
          IF(exitbigloop)EXIT
          IF(cyclebigloop)CYCLE
        ENDIF
      ENDIF

      !Read update information
      IF(code(1:6)=='update') THEN
        CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
        IF(nostrfound.OR.code2(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF

        !Read information for updating with qobs
        IF(code2(1:7)=='quseobs') THEN
          doupdate(i_quseobs) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(code3(1:10)=='allstation')THEN
            locupall = .TRUE.
          ENDIF
          IF(code3(1:9)=='nostation')THEN
            locupnone = .TRUE.
          ENDIF

        ELSEIF(code2(1:3)=='qar') THEN
          doupdate(i_qar) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)          
          IF(code3(1:9)=='nostation')THEN
            locupnone_qar = .TRUE.
          ENDIF

        !Read information for updating with wstr
        ELSEIF(code2(1:7)=='wendupd') THEN
          doupdate(i_wendupd) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(nostrfound.OR.code3(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          wupdname = code3(1:4)
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(code3(1:10)=='allstation')THEN
            wupall = .TRUE.
          ENDIF
          IF(code3(1:9)=='nostation')THEN
            wupnone = .TRUE.
          ENDIF
        ELSEIF(code2(1:3)=='war') THEN
          doupdate(i_war) = .TRUE. 
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)
          IF(nostrfound.OR.code3(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          wupdname = code3(1:4)
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code3,nostrfound,errstr)          
          IF(code3(1:9)=='nostation')THEN
            locupnone_war = .TRUE.
          ENDIF

        !Read information for updating nitrogen and phosphorus concentrations
        ELSEIF(code2(1:6)=='tpcorr') THEN
          doupdate(i_tpcorr) = .TRUE. 
        ELSEIF(code2(1:6)=='tncorr') THEN
          doupdate(i_tncorr) = .TRUE. 
        ELSEIF(code2(1:9)=='tploccorr') THEN
          doupdate(i_tploccorr) = .TRUE. 
        ELSEIF(code2(1:9)=='tnloccorr') THEN
          doupdate(i_tnloccorr) = .TRUE. 
        ENDIF
      ENDIF

      !Read criteria information
      IF(code(1:4)=='crit')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(strvalue(1:10)=='meanperiod')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) critperiod
        ELSEIF(strvalue(1:9)=='datalimit')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) critlimit
        ELSE
          READ(strvalue,*) icrit
          IF(icrit>maxcrit)THEN
            WRITE(6,*)  'ERROR: Too many criteria, max: ',maxcrit
            error = .TRUE.
            EXIT
          ENDIF
          ncrit = MAX(ncrit,icrit)
          CALL read_next_codestr_on_line(linelen,16,linepos,line,code2,nostrfound,errstr)
          IF(nostrfound.OR.code2(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          IF(code2(1:9)=='criterion')THEN
            CALL read_next_codestr_on_line(linelen,3,linepos,line,critstr,nostrfound,errstr)
            IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            c(icrit) = critstr
          ENDIF
          IF(code2(1:9)=='cvariable')THEN
            CALL read_next_codestr_on_line(linelen,6,linepos,line,varstr,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            status = find_variable_index_type(varstr,varindex,flowtype,areagg)
            IF(status>0)THEN
              WRITE(6,*)  'ERROR in variable index calculation: ', varstr,' for ',TRIM(code2)
              error = .TRUE.
              EXIT
            ENDIF
            critc(icrit) = varindex
            flowc(icrit) = flowtype
            critareac(icrit) = areagg
          ENDIF
          IF(code2(1:9)=='rvariable')THEN
            CALL read_next_codestr_on_line(linelen,6,linepos,line,varstr,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            status = find_variable_index_type(varstr,varindex,flowtype,areagg)
            IF(status>0)THEN
              WRITE(6,*)  'ERROR in variable index calculation: ', varstr,' for ',TRIM(code2)
              error = .TRUE.
              EXIT
            ENDIF
            critr(icrit) = varindex
            flowr(icrit) = flowtype
            critarear(icrit) = areagg
          ENDIF
          IF(code2(1:6)=='weight')THEN
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            READ(strvalue,*) critweight(icrit)
          ENDIF
          IF(code2(1:9)=='parameter')THEN
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            READ(strvalue,*) critpar(icrit)
          ENDIF
          IF(code2(1:11)=='conditional')THEN
            critcond(icrit) = .TRUE.
            CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
            IF(nostrfound.OR.varstr(1:2)==errstr)THEN
              error = .TRUE.
              EXIT
            ENDIF
            READ(strvalue,*) critthres(icrit)
          ENDIF
        ENDIF
      ENDIF

      !Read check flag
      IF(code(1:5)=='check')THEN
        CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
        IF(nostrfound.OR.strvalue(1:2)==errstr)THEN
          error = .TRUE.
          EXIT
        ENDIF
        IF(strvalue(1:6)=='indata')THEN
          CALL read_next_codestr_on_line(linelen,10,linepos,line,strvalue,nostrfound,errstr)
          IF(nostrfound.OR.varstr(1:2)==errstr)THEN
            error = .TRUE.
            EXIT
          ENDIF
          READ(strvalue,*) icheck
          checkdata(icheck,1) = .TRUE.
        ENDIF
      ENDIF

    ENDDO
                
    IF(error)THEN
      WRITE(6,*) 'ERROR: reading info.txt (',TRIM(filename),'). Code=',code
      status = 1
      CLOSE(fileunit_temp)
      RETURN
    ENDIF

200 CONTINUE
    ! Check consistency between optional model structure and optional forcing data input
    IF(.NOT.forcingdata(i_sfobs)%readfile .AND. modeloption(p_snowfall).EQ.1)THEN
      WRITE(6,*) 'ERROR: (info.txt) modeloption SNOWFALLMODEL can not be 1 if READSFOBS is N - please provide input file SFobs.txt and change option READSFOBS to Y!'
      status = 1
      RETURN
    ENDIF
    IF(.NOT.forcingdata(i_swobs)%readfile .AND. .NOT.forcingdata(i_tminobs)%readfile .AND. .NOT.forcingdata(i_tmaxobs)%readfile .AND. modeloption(p_petmodel).GE.3)THEN
      WRITE(6,*) 'WARNING: model option PETMODEL', modeloption(p_petmodel), 'requires shortwave radiation data (SWobs.txt) and/or Tmin and Tmax data. Model will run, but derives diurnal temperature range from clear sky turbidity, which might cause strange results....' 
    ENDIF
    IF(.NOT.forcingdata(i_tminobs)%readfile .AND. .NOT.forcingdata(i_tmaxobs)%readfile .AND. modeloption(p_petmodel).GE.4)THEN
      WRITE(6,*) 'WARNING: model option PETMODEL', modeloption(p_petmodel), 'requires Tmin data to estimate relative humidity (assuming actual vapour pressure = saturation pressure at Tmin). Model will run, but derives diurnal temperature range from clear sky turbidity, which might cause strange results.... (btw RHUM data will be optional in future releases)' 
    ENDIF
    !Check compatible simulation settings
    IF((modeloption(p_lakeriverice)>=1.AND.i_t2==0) .OR. (modeloption(p_lakeriverice)==0.AND.i_t2>0))THEN
      WRITE(6,*) 'ERROR: lakerivericemodel requires simulation of water temperature (T2), and vise versa'
      status = 1
      RETURN
    ENDIF
    IF(calibration .AND. conductregest)THEN
      WRITE(6,*) 'ERROR: Automatic calibration is not possible with reg_par.txt.'
      status = 1
      RETURN
    ENDIF
    IF(modeloption(p_floodplain)==1.AND.nsubst>0)THEN
      WRITE(6,*) 'WARNING: floodmodel 1 is not tested for substance simulations'
    ENDIF
    IF(modeloption(p_swtemperature)==1.AND.i_t2==0)THEN
      WRITE(6,*) 'ERROR: model option swtemperature 1 require T2 simulation'
      status = 1
      RETURN
    ENDIF
    IF(soildepthstretch.AND.(stateinput.OR.sonum>0))THEN
      WRITE(6,*) 'WARNING: function soildepthstretch is not recommended to be used with saved state,'
      WRITE(6,*) 'WARNING: because there is no possibility to check if the state is compatible with model set-up (par).'
    ENDIF
    IF(doassimilation)THEN
      IF(calibration)THEN
        WRITE(6,*) 'WARNING: calibration and assimilation are both ON, which is not allowed -> assimilation is kept ON and calibration is set OFF.'
        calibration = .FALSE.
      ENDIF
      IF(doupdate(i_quseobs) .OR. & 
        doupdate(i_wendupd) .OR. &
        doupdate(i_war) .OR. &
        doupdate(i_tpcorr) .OR. & 
        doupdate(i_tncorr) .OR. & 
        doupdate(i_tploccorr) .OR. & 
        doupdate(i_tnloccorr))THEN
        WRITE(6,*) 'WARNING: assimilation and updating are both ON - make sure you dont unintentially try to assimilate the same data twice with different methods!'
      ENDIF
    ENDIF
       
    !Some calculations
    temp = set_timestep_variables(tstep,stepunit,date1)    
    IF(date3%Year==0.AND.date3%Month==0.AND.date3%Day==0.AND.date3%Hour==0.AND.date3%Minute==0) date3=date1      
    skip = period_length(date1,date3)         !Calculate length of warmup period
    step = period_length(date1,AddDates(date2,steplen)) !Calculate length of simulation period
    temp = set_calibration(calibration,ncrit,c(1:ncrit),critperiod,critlimit,   &
         critweight(1:ncrit),critc(1:ncrit),critr(1:ncrit),flowc(1:ncrit),   &   !NOTE only computed variables flow and areagg is saved, combination of different type of variables may go wrong
         critpar(1:ncrit),critcond(1:ncrit),critthres(1:ncrit),critareac(1:ncrit))
    maxdim = MAX(MAXVAL(ibasinoutputvar),MAXVAL(imapoutputvar),MAXVAL(itimeoutputvar))
    noutput = count_number_of_output(max_typeofperiods,ibasinoutputvar,imapoutputvar,itimeoutputvar,iregionoutputvar)
    io = 0
    i=1
    DO
      IF(i>max_typeofperiods)EXIT
      IF(ibasinoutputvar(i)==0)EXIT
      io=io+1
      CALL check_output_variables(max_outvar,ibasinoutputvar(i),basinoutputvar(i,:,1),basinoutputvar(i,:,2))
      temp = set_output(1,io,noutput,max_outvar,basinoutputvar(i,:,1),    &
         ibasinoutputvar(i),basinoutputperiod(i),    &
         basinoutputdecimal(i),basinoutputsignif(i),basinoutputvar(i,:,2),maxoutbasins,ioutbasin(i),readbasins(i,:))
      i=i+1
    ENDDO
    i=1
    DO
      IF(i>max_typeofperiods)EXIT
      IF(imapoutputvar(i)==0)EXIT
      IF(i>1)THEN
        WRITE(6,*) 'ERROR: Only one output period for mapoutput is allowed.'
        WRITE(6,*) 'ERROR: the following outputs are ignored.'
        EXIT
      ENDIF
      io=io+1
      CALL check_output_variables(max_outvar,imapoutputvar(i),mapoutputvar(i,:,1),mapoutputvar(i,:,2))
      temp = set_output(2,io,noutput,max_outvar,mapoutputvar(i,:,1),    &
         imapoutputvar(i),mapoutputperiod(i),    &
         mapoutputdecimal(i),mapoutputsignif(i),mapoutputvar(i,:,2))
      i=i+1
    ENDDO
    i=1
    DO
      IF(i>max_typeofperiods)EXIT
      IF(itimeoutputvar(i)==0)EXIT
      io=io+1
      CALL check_output_variables(max_outvar,itimeoutputvar(i),timeoutputvar(i,:,1),timeoutputvar(i,:,2))
      temp = set_output(3,io,noutput,max_outvar,timeoutputvar(i,:,1),    &
         itimeoutputvar(i),timeoutputperiod(i),    &
         timeoutputdecimal(i),timeoutputsignif(i),timeoutputvar(i,:,2))
      i=i+1
    ENDDO
    i=1
    DO
      IF(i>max_typeofperiods)EXIT
      IF(iregionoutputvar(i)==0)EXIT
      io=io+1
      temp = set_output(4,io,noutput,max_outvar,regionoutputvar(i,:),   &
         iregionoutputvar(i),regionoutputperiod(i),       &
         regionoutputdecimal(i),regionoutputsignif(i),maxa=maxoutbasins,na=ioutregion(i),area=readregions(i,:))
      i=i+1
    ENDDO
    temp = set_maxmap(step-skip+1,mapoutputperiod(1))   !Only one mapoutput allowed
    temp = set_outvar(noutput,noutvar)
    temp = set_outvar_crit(noutvar)
    !outvar for assim is set later
    DO ivarindex = 1,sonum     !Previous time step for output of state
      pstateoutdate(ivarindex)=SubtractDates(bstateoutdate(ivarindex),steplen)
    ENDDO

    !Dates to formatted strings for hyss.log
    CALL format_date(date1,'yyyy-mm-dd HH:MM',strdate1)
    CALL format_date(date2,'yyyy-mm-dd HH:MM',strdate2)
    CALL format_date(date3,'yyyy-mm-dd HH:MM',strdate3)
    IF(sonum>0) ALLOCATE(strstates(sonum))
    DO istate=1,sonum
      CALL format_date(bstateoutdate(istate),'yyyy-mm-dd HH:MM',strstates(istate))
    ENDDO

    !Check consistency between optional model structure and time step (temporary not working)
    IF(.NOT.(steplen==datetype(0,0,1,0,0)))THEN
      IF(modeloption(p_petmodel).GE.1)THEN
        WRITE(6,*)  'ERROR: model option PETMODEL', modeloption(p_petmodel), 'does not work with other time step than day (to come in future release)' 
        status = 1
        RETURN
      ENDIF  
      IF(i_t2.GE.1 .OR. modeloption(p_lakeriverice).GE.1)THEN
        WRITE(6,*) 'ERROR: modeloption lakeriverice and water temperature model does not work with other time step than day (to come in future release)'
        status = 1
        RETURN
      ENDIF
      IF(.NOT.forcingdata(i_swobs)%readfile.AND.modeloption(p_snowmelt)==2)THEN
        WRITE(6,*) 'ERROR: modeloption snowmeltmodel=2 require SWobs.txt to work with other time step than day (for now)'
        status = 1
        RETURN
      ENDIF
    ENDIF    

    !Write information about simulation to hyss.log
    WRITE(6,*)
    WRITE(6,*)'--------Information about the simulation---------'
    WRITE(6,*)'Current model set-up:  ',TRIM(mdir)
    WRITE(6,*)'Simulation sequence:   ',simsequence
    WRITE(6,*)
    WRITE(6,*)'Simulation begin date  ',strdate1
    WRITE(6,*)'Crit/Output start date ',strdate3
    WRITE(6,*)'Simulation end date    ',strdate2
    WRITE(6,'(A23,I2,A3)') 'Simulation steplength ',tstep,stepunit
    IF(calibration) WRITE(6,*)'Calibration run'
    IF(doassimilation) WRITE(6,*)'Assimilation run'
    WRITE(6,*)'Number of substances modelled',nsubst
    WRITE(6,*)'Model options:'
    WRITE(6,*)'               PET model',modeloption(p_petmodel)
    WRITE(6,*)'          Snowfall model',modeloption(p_snowfall)
    WRITE(6,*)'          Snowmelt model',modeloption(p_snowmelt)
    WRITE(6,*)'          Snowevap model',modeloption(p_snowevap)
    WRITE(6,*)'       Snowdensity model',modeloption(p_snowdensity)
    WRITE(6,*)'     Lakeriver ice model',modeloption(p_lakeriverice)
    WRITE(6,*)'  Deep groundwater model',modeloption(p_deepgroundwater)
    WRITE(6,*)'Surface water temp model',modeloption(p_swtemperature)
    WRITE(6,*)'     Growth season model',modeloption(p_growthstart)
    WRITE(6,*)'        Floodplain model',modeloption(p_floodplain)
    WRITE(6,*)'      Infiltration model',modeloption(p_infiltration)
    WRITE(6,*)'           Erosion model',modeloption(p_erosion)
    WRITE(6,*)'  Glacier initialization',modeloption(p_glacierini)
    WRITE(6,*)
    CALL print_output_information_to_logfile(6)
    IF(writeload.AND.(conductN.OR.conductP))THEN
      WRITE(6,*) 'Yearly loads will be saved'
      conductload = .TRUE.
    ELSEIF(writeload)THEN
      WRITE(6,*) 'Yearly loads will not be saved. Need to simulate N or P.'
      writeload = .FALSE.
    ENDIF
    IF(sonum>0)THEN
      WRITE(6,*)'State variables will be saved for: ',strstates(1:sonum)
    ENDIF
    WRITE(6,*)'-------------------------------------------------'
    WRITE(6,*)

    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)
    IF(ALLOCATED(c)) DEALLOCATE(c)
    IF(ALLOCATED(critweight)) DEALLOCATE(critweight)
    IF(ALLOCATED(critc)) DEALLOCATE(critc)
    IF(ALLOCATED(critr)) DEALLOCATE(critr)
    IF(ALLOCATED(flowc)) DEALLOCATE(flowc)
    IF(ALLOCATED(flowr)) DEALLOCATE(flowr)
    IF(ALLOCATED(critareac)) DEALLOCATE(critareac)
    IF(ALLOCATED(critarear)) DEALLOCATE(critarear)
    IF(ALLOCATED(critpar)) DEALLOCATE(critpar)
    IF(ALLOCATED(strstates)) DEALLOCATE(strstates)
    IF(ALLOCATED(critcond)) DEALLOCATE(critcond)
    IF(ALLOCATED(critthres)) DEALLOCATE(critthres)

    IF(.NOT.readdaily)THEN   
      IF(step>26*366)THEN
        WRITE(6,*)
        WRITE(6,*) 'WARNING: Simulation period >26 year'
        WRITE(6,*) 'WARNING: Use readdaily if problem with saving PT to memory'
        WRITE(6,*)
      ENDIF
    ENDIF

600 FORMAT(A18000)    !linelen
    RETURN

800 CONTINUE
    status = 1
    WRITE(6,*) 'ERROR: reading info.txt. io=',io
    CLOSE(fileunit_temp)
    RETURN

  END SUBROUTINE load_coded_info

  !>
  !>
  !> \b Consequences Module variable noutput is set
  !------------------------------------------------------------------------
  INTEGER FUNCTION set_yesno_variable_from_line(linelen,linepos,line,errstr,onoff)
  
    !Argument declarations
    INTEGER,INTENT(IN) :: linelen    !<length of line
    INTEGER,INTENT(INOUT) :: linepos   !<current position on line
    CHARACTER(LEN=linelen), INTENT(IN) :: line  !<line to be read
    CHARACTER(LEN=*), INTENT(IN) :: errstr    !<value of str if error
    LOGICAL,INTENT(OUT) :: onoff  !<value of yesno variable
  
    !Local variables
    INTEGER status
    LOGICAL nostrfound
    CHARACTER(LEN=3) yesnostr
    
      status = 0
      CALL read_next_codestr_on_line(linelen,3,linepos,line,yesnostr,nostrfound,errstr)
      IF(nostrfound.OR.yesnostr(1:2)==errstr)THEN
        status = 1
      ELSEIF(yesnostr(1:1)=='Y' .OR. yesnostr(1:1)=='y' .OR.    &
             yesnostr(1:1)=='J' .OR. yesnostr(1:1)=='j')THEN
        onoff = .TRUE.
      ELSEIF(yesnostr(1:1)=='N' .OR. yesnostr(1:1)=='n')THEN
        onoff = .FALSE.
      ELSE
        status = 1
      ENDIF
      set_yesno_variable_from_line = status
          
  END FUNCTION set_yesno_variable_from_line

  !>Count number of outputs (periods) for region output
  !>
  !> \b Consequences Module worldvar variable noutput is set
  !------------------------------------------------------------------------
  INTEGER FUNCTION count_number_of_output(dim,arr1,arr2,arr3,arr4)

    INTEGER, INTENT(IN) :: dim        !<dimension of array
    INTEGER, INTENT(IN) :: arr1(dim)   !<array with number of variables for output of type 1
    INTEGER, INTENT(IN) :: arr2(dim)   !<array with number of variables for output of type 2
    INTEGER, INTENT(IN) :: arr3(dim)   !<array with number of variables for output of type 3
    INTEGER, INTENT(IN) :: arr4(dim)   !<array with number of variables for output of type 4

    !Local variables 
    INTEGER i,noutput
    
    !> \b Algorithm \n
    !>Count number of outputs with set variables, i.e. used for this simulation
    noutput = 0
    DO i = 1,dim
      IF(arr1(i)>0) noutput = noutput + 1
      IF(arr2(i)>0) noutput = noutput + 1
      IF(arr3(i)>0) noutput = noutput + 1
      IF(arr4(i)>0) noutput = noutput + 1
    ENDDO
    count_number_of_output = noutput

  END FUNCTION count_number_of_output

  !>Saves information on what and how to write to output files
  !>
  !> \b Consequences Module variable output is allocated and set
  !------------------------------------------------------------------------
  SUBROUTINE check_output_variables(dimout,n,varid,aggarea)
  
  USE MODVAR, ONLY : outvarid
    
    !Argument declarations
    INTEGER, INTENT(IN)    :: dimout       !<dimension of variable arrays
    INTEGER, INTENT(INOUT) :: n            !<number of output variables of this type
    INTEGER, INTENT(INOUT) :: varid(dimout)    !<variables for output
    INTEGER, INTENT(INOUT) :: aggarea(dimout)  !<area aggregation type of variable (subbasin=0,upstream=1,region=2)

    !Local variables 
    INTEGER i,nnew,ivar,varidnew(dimout),aggareanew(dimout)
    LOGICAL skipped(n)
    
    !> \b Algorithm \n
    !>Check outvar information for non-calculatable variables
    nnew = n
    skipped=.FALSE.
    DO i = 1,n
      IF(aggarea(i)==1)THEN   !upstream variables checked
        IF(outvarid(varid(i))%basearea==0)THEN
          !Wrong used variable, remove
          nnew = nnew - 1
          skipped(i)=.TRUE.
          WRITE(6,*) 'WARNING: Erronous output variable up'//outvarid(varid(i))%shortname//' removed from output'
        ENDIF
      ENDIF
    ENDDO
    
    IF(n>nnew)THEN
      !>Change output variable information
      ivar = 0
      DO i = 1,n
        IF(.NOT.skipped(i))THEN
          ivar = ivar + 1
          varidnew(ivar)   = varid(i)
          aggareanew(ivar) = aggarea(i)
        ENDIF
      ENDDO
      !>Set output
      n = nnew
      varid = varidnew
      aggarea = aggareanew
    ENDIF

  END SUBROUTINE check_output_variables

  !>Saves information on variables to write to output files
  !>
  !> \b Consequences Module variable outvarinfotemp and outvarindex is allocated.
  !------------------------------------------------------------------------
  INTEGER FUNCTION set_outvar(nout,noutvar)

    USE WORLDVAR, ONLY : output, & !OUT
                         outvarinfotemp !OUT
    USE MODVAR, ONLY : max_noutvar, &
                       max_outvar, &
                       outvarindex  !OUT
    
    !Argument declaration
    INTEGER, INTENT(IN) :: nout         !<number of output file type to be set
    INTEGER, INTENT(OUT) :: noutvar     !<number of output variables (so far)

    INTEGER i,iout
    
    !>\b Algorithm \n
    set_outvar  = 0
    
    !>Initialize and allocate variables holding outvar information
    noutvar = 0 !Initialize here, will be changed in subroutines
    IF(.NOT.ALLOCATED(outvarinfotemp)) ALLOCATE(outvarinfotemp(max_noutvar))  !temporary allocation, big structure
    IF(.NOT.ALLOCATED(outvarindex)) ALLOCATE(outvarindex(max_outvar))
    outvarindex = 0 !default is turned off
    
    !>For every output file type...
    DO iout = 1,nout
      !>and every variable..
      !>\li Set outvar information
      DO i = 1,output(iout)%nvar
        CALL set_outvar_for_variable(output(iout)%variable(i)%idindex,output(iout)%variable(i)%areaagg,output(iout)%variable(i)%ovindex,noutvar)
      ENDDO
    ENDDO

  END FUNCTION set_outvar

  !>Saves information on output variables to calculate for use in criteria calculations
  !>
  !> \b Consequences Module variable acccalvar, outvarinfotemp, noutvar and 
  !>outvarindex is changed.
  !------------------------------------------------------------------------
  INTEGER FUNCTION set_outvar_crit(noutvar)

    USE WORLDVAR, ONLY : nacrit, &
                         acccalvar
    
    INTEGER, INTENT(INOUT) :: noutvar  !<number of output variables to be set (so far)
    INTEGER ia
    
    !> \b Algorithm \n
    set_outvar_crit  = 0
    
    DO ia = 1,nacrit
      CALL set_outvar_for_variable(acccalvar(ia)%comp,acccalvar(ia)%areaagg,acccalvar(ia)%compoutvar,noutvar)
      CALL set_outvar_for_variable(acccalvar(ia)%rec,acccalvar(ia)%areaagg,acccalvar(ia)%recoutvar,noutvar)
    ENDDO
    
  END FUNCTION set_outvar_crit

  !>\brief Reads file with information about submodel for the current
  !run and gets the coupling to the model set-up
  !!
  !!Format of file: first row is number of subbasins, and the
  !!following holds the subid of these subbasins
  !--------------------------------------------------------------------
  SUBROUTINE load_submodel_info(dir,submodelrun,msub,status)

    USE WORLDVAR, ONLY : fileunit_temp,   &
                         ibasemodel
    USE MODVAR, ONLY : nsub_basemodel, &
                       basin

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory for subinfo-file
    LOGICAL, INTENT(IN)  :: submodelrun  !<Submodel is supposed to be run
    INTEGER, INTENT(OUT) :: msub         !<Number of subbasins of submodel
    INTEGER, INTENT(OUT) :: status       !<Error status
    
    !Local variables
    LOGICAL fileexist
    INTEGER i,j
    CHARACTER (LEN=100) filename
    INTEGER, ALLOCATABLE :: selectedaro(:)

    status = 0

    !File exist for simulation with submodel?
    filename=TRIM(dir)//'pmsf.txt'
    INQUIRE(FILE=filename, EXIST=fileexist)
    IF(submodelrun.AND.(.NOT.fileexist))THEN
      WRITE(6,*) 'ERROR: File is missing:',TRIM(filename)
      status = 1
      RETURN
    ENDIF

    !Read subbasins to be simulated
    IF(submodelrun)THEN
      OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
      READ(fileunit_temp,*) msub
      IF(.NOT.ALLOCATED(selectedaro)) ALLOCATE(selectedaro(msub))
      READ(fileunit_temp,*) selectedaro
      CLOSE(fileunit_temp)
      WRITE(6,*) 'File read: ',TRIM(filename)
    ELSE
      msub = nsub_basemodel
    ENDIF

    !Set reordering of subbasins if necessary
    IF(.NOT.ALLOCATED(ibasemodel)) ALLOCATE(ibasemodel(msub))
    ibasemodel(:)=-9999
    IF(submodelrun)THEN
      DO i = 1,nsub_basemodel
        DO j = 1,msub
          IF(basin(i)%subid==selectedaro(j))THEN
            ibasemodel(j)=i
            EXIT
          ENDIF
        ENDDO
      ENDDO
      !Check if all subbasins were found
      j=0 
      DO i=1,msub 
        IF(ibasemodel(i)==-9999) THEN
          IF(j==0) THEN
            WRITE(6,*) ' '
            WRITE(6,*) '---------------------------------------------------'
            WRITE(6,*) 'ERROR Subbasins in pmsf.txt not found i GeoData.txt'
            WRITE(6,*) '       SUBID:       PMSF#'
            j=1
          ENDIF
          WRITE(6,*) selectedaro(i), i
        ENDIF
      ENDDO
      IF(j==1) THEN
        status = 1
        RETURN
      ENDIF
      IF(ALLOCATED(selectedaro)) DEALLOCATE(selectedaro)   

    ELSE
      !Set index straight if no submodel is simulated
      DO i = 1,nsub_basemodel
        ibasemodel(i)=i
      ENDDO
    ENDIF

    !Write information to hyss.log
    IF(submodelrun)THEN
      WRITE(6,*)
      WRITE(6,*)'-----Information about the simulation (cont.)----'
      WRITE(6,*)'Submodel defined by file: ',TRIM(filename)
      WRITE(6,*)'Number of chosen subbasins to simulate',msub
      WRITE(6,*)'-------------------------------------------------'
      WRITE(6,*)
    ENDIF

  END SUBROUTINE load_submodel_info

  !>Preparation of input data for using pseudo dayno (year starting  
  !>July 1) on southern hemisphere regions/subbasins
  !>
  !>\b Consequences Module modvar variables cropdata and cropirrdata 
  !>may be changed.
  !-----------------------------------------------------------------
  SUBROUTINE reform_input_data_for_pseudo_dayno(n)
  
  USE MODVAR, ONLY : basin,  &
                     cropdata,  &   !INOUT
                     cropirrdata, & !INOUT
                     cropindex
                     
    !Argument declarations
    INTEGER, INTENT(IN) :: n  !<number of subbasins (basemodel)
    
    !Local variables
    INTEGER i,j,d1,d2,reg
    INTEGER,ALLOCATABLE :: rows(:)
    LOGICAL,ALLOCATABLE :: southreg(:)  !help variable to keep track of already changed regions
    
    !Check if data exist that needs to be changed
    IF(.NOT.ALLOCATED(cropdata) .AND. .NOT.ALLOCATED(cropirrdata)) RETURN
    
    d1 = SIZE(cropindex,1) !crop
    d2 = SIZE(cropindex,2) !region
    ALLOCATE(rows(d1))
    ALLOCATE(southreg(d2))
    southreg = .FALSE.
    
    DO i = 1,n
      IF(basin(i)%latitude<0.)THEN
        reg = basin(i)%region
        !Change cropdata-components to pseudo dayno for all crops in this region
        IF(.NOT.southreg(reg))THEN  !check if already changed this region
          rows = cropindex(:,reg)   !rows in cropdata
          DO j=1,d1
            IF(rows(j)>0) CALL change_cropdata_to_pseudo_dayno(rows(j),cropdata,cropirrdata)
          ENDDO
          southreg(reg)=.TRUE.
        ENDIF
      ENDIF
    ENDDO

    DEALLOCATE(rows,southreg)  !for gfortran
  
  END SUBROUTINE reform_input_data_for_pseudo_dayno
  
  !>Change dayno-variables to pseudo dayno for crops in southern hemisphere
  !--------------------------------------------------------------------
  SUBROUTINE change_cropdata_to_pseudo_dayno(row,cropdata,cropirrdata)
  
  USE MODVAR, ONLY : CROPDATATYPE,  &
                     CROPIRRDATATYPE, &
                     get_pseudo_dayno
                     
    !Argument declarations
    INTEGER, INTENT(IN) :: row  !<row in cropdata to change for
    TYPE(CROPDATATYPE),ALLOCATABLE,INTENT(INOUT) :: cropdata(:) !<Crop characteristics
    TYPE(CROPIRRDATATYPE),ALLOCATABLE,INTENT(INOUT) :: cropirrdata(:) !<Crop characteristics regarding irrigation
    
    !Get pseudo dayno for cropdata dayno in southern hemisphere
    IF(ALLOCATED(cropdata))THEN
      cropdata(row)%fertday1 = get_pseudo_dayno(cropdata(row)%fertday1,-1.)
      cropdata(row)%fertday2 = get_pseudo_dayno(cropdata(row)%fertday2,-1.)
      cropdata(row)%manday1 = get_pseudo_dayno(cropdata(row)%manday1,-1.)
      cropdata(row)%manday2 = get_pseudo_dayno(cropdata(row)%manday2,-1.)
      IF(cropdata(row)%resdayno>0) cropdata(row)%resdayno = get_pseudo_dayno(cropdata(row)%resdayno,-1.)
      IF(cropdata(row)%baredayno1>0) cropdata(row)%baredayno1 = get_pseudo_dayno(cropdata(row)%baredayno1,-1.)
      IF(cropdata(row)%baredayno2>0) cropdata(row)%baredayno2 = get_pseudo_dayno(cropdata(row)%baredayno2,-1.)
      IF(cropdata(row)%baredayno3>0) cropdata(row)%baredayno3 = get_pseudo_dayno(cropdata(row)%baredayno3,-1.)
      IF(cropdata(row)%baredayno4>0) cropdata(row)%baredayno4 = get_pseudo_dayno(cropdata(row)%baredayno4,-1.)
      IF(cropdata(row)%baredayno5>0) cropdata(row)%baredayno5 = get_pseudo_dayno(cropdata(row)%baredayno5,-1.)
      cropdata(row)%firstday = get_pseudo_dayno(cropdata(row)%firstday,-1.)
    ENDIF
    IF(ALLOCATED(cropirrdata))THEN
      cropirrdata(row)%plantingdayno = get_pseudo_dayno(cropirrdata(row)%plantingdayno,-1.)
      IF(cropirrdata(row)%imm_start>0)THEN
        cropirrdata(row)%imm_start = get_pseudo_dayno(cropirrdata(row)%imm_start,-1.)
        cropirrdata(row)%imm_end = get_pseudo_dayno(cropirrdata(row)%imm_end,-1.)
      ENDIF
    ENDIF
    
  END SUBROUTINE change_cropdata_to_pseudo_dayno
  
  !>Turns on calculation of choosen variable
  !------------------------------------------------------------------------
  SUBROUTINE set_outvar_for_variable(varid,areaagg,ovindex,noutvar)
    
    USE MODVAR, ONLY : outvarindex, &
                       outvarid, &
                       max_noutvar
    USE WORLDVAR, ONLY : OUTVARINFOTYPE, &
                         outvarinfotemp !OUT

    !Argument declaration
    INTEGER, INTENT(IN) :: varid    !<index in outvarid of current variable
    INTEGER, INTENT(IN) :: areaagg  !<aggregation of variable (0=subbasin average, 1=upstream, 2=region)
    INTEGER, INTENT(OUT) :: ovindex !<index of current variable in outvar
    INTEGER, INTENT(INOUT) :: noutvar !<number of output variables to be set (so far)

    !Local variables 
    INTEGER ioutvar,i2,idindex,icurrent,icurrent2
    LOGICAL found,found2,found3,found4
    
    !> \b Algorithm \n
    IF(noutvar>max_noutvar-5) WRITE(6,*) 'WARNING: A lot of output variables has been chosen. Program may crash.'
    ioutvar = noutvar
    
    found = .FALSE.
    !Check if the variable already is in outvar
    DO i2=1,ioutvar
      IF(outvarinfotemp(i2)%idindex == varid)THEN
        IF(outvarinfotemp(i2)%areaagg==areaagg)THEN
          found = .TRUE.
          ovindex = i2
          EXIT
        ENDIF
      ENDIF
    ENDDO
    IF(.NOT.found)THEN
      !If not, add the variable to outvar...
      ioutvar = ioutvar + 1
      icurrent = ioutvar
      idindex = varid
      IF(areaagg==0) outvarindex(idindex) = icurrent
      ovindex = icurrent
      outvarinfotemp(icurrent) = OUTVARINFOTYPE(idindex,areaagg,ioutvar,0)
      !..together with its water-variable if applicable
      IF(outvarid(idindex)%water>0)THEN
        found2 = .FALSE.
        DO i2=1,ioutvar-1
          IF(outvarinfotemp(i2)%idindex == outvarid(idindex)%water)THEN
            IF(outvarinfotemp(i2)%areaagg==areaagg)THEN
              found2 = .TRUE.
              outvarinfotemp(icurrent)%ovindexwater = i2
              EXIT
            ENDIF
          ENDIF
        ENDDO
        IF(.NOT.found2)THEN
          ioutvar = ioutvar + 1
          icurrent2 = ioutvar
          idindex = outvarid(idindex)%water
          IF(areaagg==0) outvarindex(idindex) = icurrent2
          outvarinfotemp(icurrent2) = OUTVARINFOTYPE(idindex,areaagg,ioutvar,0)
          outvarinfotemp(icurrent)%ovindexwater = ioutvar
          !...If the variable is an average of a large area, add the subbasin average of the water variable
          IF(outvarinfotemp(icurrent2)%areaagg>0)THEN
            found3 = .FALSE.
            DO i2=1,icurrent2-1
              IF(outvarinfotemp(i2)%idindex == idindex)THEN
                IF(outvarinfotemp(i2)%areaagg==0)THEN
                  found3 = .TRUE.
                  outvarinfotemp(icurrent2)%ovindex0 = i2
                  EXIT
                ENDIF
              ENDIF
            ENDDO
            IF(.NOT.found3)THEN
              ioutvar = ioutvar + 1
              outvarindex(idindex) = ioutvar
              outvarinfotemp(icurrent2)%ovindex0 = ioutvar
              outvarinfotemp(ioutvar) = OUTVARINFOTYPE(idindex,0,ioutvar,0)
            ENDIF
          ENDIF
        ENDIF
      ENDIF
      !Check if subbasin average of the variable is already in outvar
      IF(outvarinfotemp(icurrent)%areaagg>0)THEN
        found3 = .FALSE.
        DO i2=1,ioutvar-1
          IF(outvarinfotemp(i2)%idindex == varid)THEN
            IF(outvarinfotemp(i2)%areaagg==0)THEN
              found3 = .TRUE.
              outvarinfotemp(icurrent)%ovindex0 = i2
              EXIT
            ENDIF
          ENDIF
        ENDDO
        IF(.NOT.found3)THEN
          !If not, add the subbasin average variable to outvar...
          ioutvar = ioutvar + 1
          icurrent2 = ioutvar
          idindex = varid
          outvarindex(idindex) = ioutvar
          outvarinfotemp(icurrent)%ovindex0  = ioutvar
          outvarinfotemp(icurrent2) = OUTVARINFOTYPE(idindex,0,ioutvar,0)
          IF(outvarid(idindex)%water>0)THEN
            found4 = .FALSE.
            DO i2=1,ioutvar-1
              IF(outvarinfotemp(i2)%idindex == outvarid(idindex)%water)THEN
                IF(outvarinfotemp(i2)%areaagg==0)THEN
                  outvarinfotemp(icurrent2)%ovindexwater = i2
                  found4 = .TRUE.
                  EXIT
                ENDIF
              ENDIF
            ENDDO
            IF(.NOT.found4)THEN   !maybe I already added this one above
              !..together with its water-variable if applicable
              ioutvar = ioutvar + 1
              idindex = outvarid(idindex)%water
              outvarindex(idindex) = ioutvar
              outvarinfotemp(ioutvar) = OUTVARINFOTYPE(idindex,0,ioutvar,0)
              outvarinfotemp(icurrent2)%ovindexwater = ioutvar
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDIF
    
    !>Increase the counter of number of output variables to calculate
    noutvar = ioutvar

  END SUBROUTINE set_outvar_for_variable

  !>Preparation for subbasin and regional output files
  !--------------------------------------------------------------------
  SUBROUTINE prepare_subbasin_output(status)

    USE MODVAR, ONLY : basin,            &
                       nsub, &
                       nsub_basemodel
    USE WORLDVAR, ONLY : output, &
                         noutput, &
                         noutreg, &
                         outregion
    
    !Argument declarations
    INTEGER, INTENT(OUT) :: status        !<Error status
    
    !Local variables
    INTEGER isub,io,ia,ireg,iaccept
    INTEGER na
    INTEGER tempoutareas(nsub_basemodel)

    status = 0
   
    !Find and set index of areas to be written (areas was given as ids)
    !and reduce them to only those in model
    DO io=1,noutput
      na = output(io)%narea
      IF(na==0)THEN   !all out
        IF(output(io)%fileformat==4)THEN
          na = noutreg
          IF(ALLOCATED(output(io)%areaindex)) DEALLOCATE(output(io)%areaindex)
          ALLOCATE(output(io)%areaindex(na))
          tempoutareas(1:na) = outregion(1:na)%outregid
          iaccept = 0
          DO ia = 1,na
            DO ireg = 1,noutreg
              IF(tempoutareas(ia)==outregion(ireg)%outregid .AND. &
                 outregion(ireg)%nsubbasin>0)THEN
                iaccept = iaccept + 1
                output(io)%areaindex(iaccept) = ireg
                EXIT
              ENDIF
            ENDDO
            IF(ireg>noutreg)THEN
              WRITE(6,*) 'Warning: regionoutput region', tempoutareas(ia), 'not found in model'
            ENDIF
          ENDDO
          IF(iaccept<na) output(io)%narea = iaccept
        ENDIF
      ELSE
        tempoutareas(1:na) = output(io)%areaindex !subid or outregid
        IF(output(io)%fileformat==1)THEN
          iaccept = 0
          DO ia = 1,na
            DO isub = 1,nsub
              IF(tempoutareas(ia)==basin(isub)%subid)THEN
                iaccept = iaccept + 1
                output(io)%areaindex(iaccept) = isub
                EXIT
              ENDIF
            ENDDO
            IF(isub>nsub)THEN
              WRITE(6,*) 'Warning: basinoutput subbasin', tempoutareas(ia), 'not in model'
            ENDIF
          ENDDO
          IF(iaccept<na) output(io)%narea = iaccept
        ELSEIF(output(io)%fileformat==4)THEN
          iaccept = 0
          DO ia = 1,na
            DO ireg = 1,noutreg
              IF(tempoutareas(ia)==outregion(ireg)%outregid .AND. &
                 outregion(ireg)%nsubbasin>0)THEN
                iaccept = iaccept + 1
                output(io)%areaindex(iaccept) = ireg
                EXIT
              ENDIF
            ENDDO
            IF(ireg>noutreg)THEN
              WRITE(6,*) 'Warning: regionoutput region', tempoutareas(ia), 'not found in model'
            ENDIF
          ENDDO
          IF(iaccept<na) output(io)%narea = iaccept
        ENDIF
      ENDIF
    ENDDO
    
  END SUBROUTINE prepare_subbasin_output

  !>Load file with information on output regions; Outregions.txt
  !
  !>\b Consequences Module worldvar variable outregion is allocated 
  !>and set, and variable noutreg is set.
  !--------------------------------------------------------------------
  SUBROUTINE load_output_regions(dir,nreg,status)

    USE MODVAR, ONLY : basin,            &
                       nsub
    USE WORLDVAR, ONLY : i_str,       &
                         i_intg,      &
                         i_real,      &
                         maxcharpath, &
                         readoutregion
    USE WORLDVAR, ONLY : noutreg,  &    !OUT
                         outregion, &   !OUT
                         fileunit_temp
    
    !Argument declarations
    CHARACTER (LEN=*), INTENT(IN) :: dir   !<File directory
    INTEGER, INTENT(OUT) :: nreg           !<Number of output regions
    INTEGER, INTENT(OUT) :: status         !<Error status

    !Local variables
    INTEGER ncols,nrows                !size of data in file
    INTEGER mcol,nmax
    INTEGER i,j,isub,ix
    LOGICAL fileex                     !Status of file
    CHARACTER(LEN=maxcharpath) infile             !Name of irrigation characteristics file 
    INTEGER, ALLOCATABLE :: xi(:,:)               !Integer data read from file
    INTEGER, ALLOCATABLE :: code(:)               !Code for column variable
    INTEGER, ALLOCATABLE :: rindex(:)             !Index for column real variables
    INTEGER, ALLOCATABLE :: iindex(:)             !Index for column integer variables
    REAL, ALLOCATABLE    :: xr(:,:)               !Real data read from file
    INTEGER, ALLOCATABLE :: readsubid(:,:)        !Subid in file
    REAL, ALLOCATABLE    :: readweight(:,:)       !Weights in file
    CHARACTER(LEN=10), ALLOCATABLE :: str(:)      !File column content string

    status = 0
    nreg = 0
    infile = TRIM(dir)//'Outregions.txt'
    IF(.NOT.readoutregion) RETURN

    !Check if file exist
    INQUIRE(FILE=TRIM(infile),EXIST=fileex)
    IF(.NOT.fileex)THEN
      WRITE(6,*) 'Error: No Outregions.txt file found.'
      status = 1
      RETURN
    ENDIF

    !Count number of columns in OutregionData
    CALL count_data_cols(fileunit_temp,TRIM(infile),0,ncols,status)
    IF(status/=0)RETURN

    !Count number of regions in OutregionData
    CALL count_data_rows(fileunit_temp,TRIM(infile),1,nrows,status)
    IF(status/=0)RETURN
    IF(nrows==0)THEN
      WRITE(6,*) 'ERROR: No data in Outregions.txt file found.'
      status = 1
      RETURN
    ENDIF

    !Allocate and initiate local variables
    IF(.NOT.ALLOCATED(xi)) ALLOCATE(xi(nrows,ncols))  
    IF(.NOT.ALLOCATED(code)) ALLOCATE(code(ncols))  
    IF(.NOT.ALLOCATED(rindex)) ALLOCATE(rindex(ncols))
    IF(.NOT.ALLOCATED(iindex)) ALLOCATE(iindex(ncols))
    IF(.NOT.ALLOCATED(xr)) ALLOCATE(xr(nrows,ncols))
    IF(.NOT.ALLOCATED(str)) ALLOCATE(str(ncols)) 
    ALLOCATE(readsubid(nsub,nrows))   !nsub = nsub_base
    ALLOCATE(readweight(nsub,nrows))
    readsubid = 0

    !Open file
    OPEN(UNIT = fileunit_temp,FILE = TRIM(infile), STATUS = 'old', ACTION='read')     

    !Read the column headings from file
    CALL read_column_headings(fileunit_temp,ncols,str,mcol,status)
    IF(status.NE.0) THEN
      WRITE(6,*) 'ERROR reading file: ',TRIM(infile)
      RETURN
    ENDIF

    !>If not allocated: allocate and initialize outregion
    noutreg = nrows !module variable
    IF(noutreg>nsub)THEN
      WRITE(6,*) 'ERROR: HYSS-HYPE can not handle more than nsub output regions'
      STOP 1
    ENDIF
    nreg = nrows    !subroutine argument
    IF(.NOT.ALLOCATED(outregion)) ALLOCATE(outregion(noutreg))

    !Code variables for easy finding of variable type
    code=i_str    !string, ignore
    DO i = 1,ncols
      IF(str(i)(1:10)=='outregid  ') code(i) = i_intg
      IF(str(i)(1:10)=='xcoord    ') code(i) = i_real
      IF(str(i)(1:10)=='ycoord    ') code(i) = i_real
      IF(str(i)(1:10)=='zcoord    ') code(i) = i_real
      IF(str(i)(1:6)=='weight')      code(i) = i_real
      IF(str(i)(1:5)=='subid')       code(i) = i_intg
    ENDDO

    !Read all data
    CALL read_basindata5(fileunit_temp,infile,ncols,nrows,ncols,code,rindex,iindex,xi,xr)

    CLOSE(UNIT=fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(infile)

    isub = 0
    DO i = 1,ncols
      IF(str(i)(1:10)=='outregid  ') outregion(1:nrows)%outregid = xi(1:nrows,iindex(i))
      IF(str(i)(1:10)=='xcoord    ') outregion(1:nrows)%xcoord   = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='ycoord    ') outregion(1:nrows)%ycoord   = xr(1:nrows,rindex(i))
      IF(str(i)(1:10)=='zcoord    ') outregion(1:nrows)%zcoord   = xr(1:nrows,rindex(i))
      IF(str(i)(1:5)=='subid')THEN
        isub = isub + 1
        readsubid(isub,1:nrows) = xi(1:nrows,iindex(i))
      ENDIF
      IF(str(i)(1:6)=='weight')THEN
        readweight(isub,1:nrows) = xr(1:nrows,rindex(i))
      ENDIF
    ENDDO
    nmax = isub
    
    !Allocate and find subbasin index and set outregion components
    DO i = 1,noutreg
      outregion(i)%nsubbasin = nmax
      DO j = nmax,1,-1
        IF(readsubid(j,i)==0) outregion(i)%nsubbasin = outregion(i)%nsubbasin-1
      ENDDO
      IF(.NOT.ALLOCATED(outregion(i)%subindex)) ALLOCATE(outregion(i)%subindex(outregion(i)%nsubbasin))
      IF(.NOT.ALLOCATED(outregion(i)%weight)) ALLOCATE(outregion(i)%weight(outregion(i)%nsubbasin))
      outregion(i)%weight = readweight(1:outregion(i)%nsubbasin,i)
      DO isub=1,outregion(i)%nsubbasin
        DO ix=1,nsub
          IF(readsubid(isub,i)==basin(ix)%subid)THEN
            outregion(i)%subindex(isub) = ix
            EXIT
          ENDIF
        ENDDO
      ENDDO
    ENDDO

    !Deallocate local variables
    IF(ALLOCATED(xi)) DEALLOCATE(xi)  
    IF(ALLOCATED(code)) DEALLOCATE(code)  
    IF(ALLOCATED(rindex)) DEALLOCATE(rindex)
    IF(ALLOCATED(iindex)) DEALLOCATE(iindex)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)
    IF(ALLOCATED(str)) DEALLOCATE(str) 
    IF(ALLOCATED(readsubid)) DEALLOCATE(readsubid) 
    IF(ALLOCATED(readweight)) DEALLOCATE(readweight)

  END SUBROUTINE load_output_regions
  

  !>Writes the files with data in format suitable for mapping.
  !--------------------------------------------------------------------
  SUBROUTINE save_mapfiles(dir,numsub,ntime,iens,runens,allens)

    USE WORLDVAR, ONLY : maptime,     &
                         writematlab, &
                         resultseq,   &
                         simsequence, &
                         fileunit_temp, &
                         output,      &
                         noutput,     &
                         i_t,i_d,i_w,i_m,i_s
    USE MODVAR, ONLY : outvarid,      &
                       basin,         &
                       i_sum
    USE READWRITE_ROUTINES, ONLY : write_dataline, &
                                   write_commentline

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
    INTEGER, INTENT(IN) :: numsub        !<Number of subbasins
    INTEGER, INTENT(IN) :: ntime         !<Number of values in vector = number of timesteps
    INTEGER, INTENT(IN) :: iens          !<Current simulation
    LOGICAL, INTENT(IN) :: runens        !<Flag for ensemble simulation
    LOGICAL, INTENT(IN) :: allens        !<Flag for writing all ensemble results

    !Local variables
    INTEGER io,ivar
    INTEGER itime
    INTEGER isub
    INTEGER nn         !File prefix number
    CHARACTER(LEN=20) filename  !Current filename
    CHARACTER(LEN=20) var
    CHARACTER(LEN=40) var2
    CHARACTER(LEN=30) unit
    REAL x(ntime)
    INTEGER lt,lout
    CHARACTER (LEN=16)   t
    CHARACTER (LEN=200) comment
    CHARACTER (LEN=1800000) outtxt
    LOGICAL writeperiod

    !filename=''
    DO io = 1,noutput
      IF(output(io)%fileformat==2)THEN
        DO ivar = 1,output(io)%nvar
          IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. mapCOUT_005.txt
            IF(resultseq.AND.simsequence>0) nn=simsequence
            IF(runens) nn=iens
            CALL create_filename_for_variable(filename,'map',outvarid(output(io)%variable(ivar)%idindex)%shortname,output(io)%variable(ivar)%areaagg,suffix2=nn,suffix2format=allens)
          ELSE  !Ex. mapCOUT.txt
            CALL create_filename_for_variable(filename,'map',outvarid(output(io)%variable(ivar)%idindex)%shortname,output(io)%variable(ivar)%areaagg)
          ENDIF
          var=outvarid(output(io)%variable(ivar)%idindex)%longname
          IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
            var2='regional average of '//var
          ELSEIF(output(io)%variable(ivar)%areaagg==1)THEN    !upstream output
            var2='upstream average of '//var
          ELSE
            var2=var
          ENDIF
          unit=''
          unit=outvarid(output(io)%variable(ivar)%idindex)%longunit
          IF(outvarid(output(io)%variable(ivar)%idindex)%vartype==i_sum)THEN
            SELECT CASE(output(io)%period)
            CASE(i_t)
              unit=TRIM(unit)//' per timestep'
            CASE(i_d)
              unit=TRIM(unit)//' per day'
            CASE(i_w)
              unit=TRIM(unit)//' per week'
            CASE(i_m)
              unit=TRIM(unit)//' per month'
            CASE DEFAULT
              unit=TRIM(unit)//' per year'
            END SELECT
          ENDIF
          OPEN(file=TRIM(dir)//filename,unit=fileunit_temp,status='unknown',form='formatted')

          !Write heading
          IF(writematlab)THEN
            comment = '%This is a file with '//TRIM(var2)//' in '//TRIM(unit)//' for GIS mapping'
          ELSE
            comment = 'This is a file with '//TRIM(var2)//' in '//TRIM(unit)//' for GIS mapping'
          ENDIF
          CALL write_commentline(fileunit_temp,comment)
          outtxt = ''
          IF(writematlab)THEN
            IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
              outtxt(1:10) = '%OUTREGID'//CHAR(44)
              lout = 10
            ELSE
              outtxt(1:7) = '%SUBID'//CHAR(44)
              lout = 7
            ENDIF
          ELSE
            IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
              outtxt(1:11) = '!!OUTREGID'//CHAR(44)
              lout = 11
            ELSE
              outtxt(1:6) = 'SUBID'//CHAR(44)
              lout = 6
            ENDIF
          ENDIF
          IF(output(io)%period==i_s) THEN
            writeperiod=.TRUE.
          ELSE
            writeperiod=.FALSE.
          ENDIF
          DO itime = 1,ntime
            t = ''
            t = TRIM(maptime(itime))
            t = ADJUSTL(t)
            lt = LEN_TRIM(t)
            outtxt(lout+1:lout+lt) = t(1:lt)
            lout = lout+lt
            IF(itime < ntime) THEN
              outtxt(lout+1:lout+1) = CHAR(44)   !Also the last one necessary
              lout = lout+1                      !for reading in "free format"
            ENDIF
          ENDDO
          WRITE(fileunit_temp,'(a)') outtxt(1:lout) 

          !Write data
          DO isub=1,numsub
            x=output(io)%accdata%value(isub,ivar,1:ntime)
            CALL write_dataline(fileunit_temp,ntime,x,output(io)%decimal,output(io)%signfig,&
                  0,',',0,writematlab,id=basin(isub)%subid)
          ENDDO
          CLOSE(fileunit_temp)
        ENDDO
      ENDIF
    ENDDO

    RETURN
  END SUBROUTINE save_mapfiles

  !>Writes subbasin/output region assessment to file subassX.txt for each pair of compared variables
  !--------------------------------------------------------------------
  SUBROUTINE write_subbasin_assessment(filedir,numsub,na,sas2,iens,runens)

    USE WORLDVAR, ONLY : acccalvar,     &
                         calvarper,     &
                         outregion,     &
                         noutreg,       &
                         fileunit_temp, &
                         writematlab,   &
                         resultseq,     &
                         simsequence,   &
                         outperiodname, &
                         maxsubass
    USE MODVAR, ONLY :   outvarid,      &
                         nsub,          &
                         basin,         &
                         missing_value      
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: filedir       !<File directory
    INTEGER, INTENT(IN) :: numsub                 !<dimension of performance measures (areas)
    INTEGER, INTENT(IN) :: na             !<dimension of performance measures (variables)
    REAL, INTENT(IN)    :: sas2(numsub,maxsubass,na)  !<assessment values
    INTEGER, INTENT(IN) :: iens                   !<Current simulation
    LOGICAL, INTENT(IN) :: runens                 !<Flag for ensemble simulation
    
    !Local variables
    INTEGER ia         !Criteria number
    INTEGER isub       !Subbasin
    INTEGER nn         !File suffix number
    CHARACTER(LEN=20) filename
    REAL x(maxsubass)  !Values
    INTEGER ndec       !Number of decimals to be written (rounding to ndec decimals will occurr)
    CHARACTER(LEN=1) tal

    !>\b Algorithm \n
    tal = '0'
    !>For every pair of variables to be compared 
    DO ia = 1, na
      !> \li Determine the filename and open a file
      tal = CHAR(ICHAR(tal)+1)    !tal=i
      IF(tal==':') tal = 'A'      !för fler kriterier än 9 används bokstäver
      IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. subass1_005.txt
        IF(resultseq.AND.simsequence>0) nn=simsequence
        IF(runens) nn=iens
        filename(1:7) = 'subass'//tal
        filename(8:8) = '_'
        WRITE(filename(9:11),601) nn
        filename(12:15) = '.txt'
      ELSE  !Ex. subass1.txt
        filename='subass'//tal//'.txt'
      ENDIF
601   FORMAT(I3.3)
      OPEN(file=TRIM(filedir)//filename,unit=fileunit_temp,status='unknown',form='formatted')
      
      !> \li Write heading to file
      IF(acccalvar(ia)%areaagg==0)THEN
      IF(writematlab)THEN
          WRITE(fileunit_temp,'(a56,a2,a13,a4,a2,a4,a8,a6)')    &
              '%Subbasin assessment. Criteria is calculated for period ',outperiodname(calvarper),'. Variables: ',outvarid(acccalvar(ia)%rec)%shortname,', ',outvarid(acccalvar(ia)%comp)%shortname,' Unit: ',outvarid(acccalvar(ia)%comp)%shortunit
        WRITE(fileunit_temp,'(a)') '%SUBID'//CHAR(9)//'NSE'//CHAR(9)//'CC'//CHAR(9)//'RE(%)'//CHAR(9)//'RSDE(%)'//CHAR(9)//'Sim'//CHAR(9)//'Rec'//CHAR(9)//'SDSim'//CHAR(9)//'SDRec'//CHAR(9)//'MAE'//CHAR(9)//'RMSE'//CHAR(9)//'Bias'//CHAR(9)//'SDE'//CHAR(9)//'KGE'//CHAR(9)//'KGESD'//CHAR(9)//'KGEM'//CHAR(9)//'NRMSE'//CHAR(9)//'NSEW'
      ELSE
          WRITE(fileunit_temp,'(a57,a2,a13,a4,a2,a4,a8,a6)')    &
                '!!Subbasin assessment. Criteria is calculated for period ',outperiodname(calvarper),'. Variables: ',outvarid(acccalvar(ia)%rec)%shortname,', ',outvarid(acccalvar(ia)%comp)%shortname,' Unit: ',outvarid(acccalvar(ia)%comp)%shortunit
        WRITE(fileunit_temp,'(a)') 'SUBID'//CHAR(9)//'NSE'//CHAR(9)//'CC'//CHAR(9)//'RE(%)'//CHAR(9)//'RSDE(%)'//CHAR(9)//'Sim'//CHAR(9)//'Rec'//CHAR(9)//'SDSim'//CHAR(9)//'SDRec'//CHAR(9)//'MAE'//CHAR(9)//'RMSE'//CHAR(9)//'Bias'//CHAR(9)//'SDE'//CHAR(9)//'KGE'//CHAR(9)//'KGESD'//CHAR(9)//'KGEM'//CHAR(9)//'NRMSE'//CHAR(9)//'NSEW'
      ENDIF
      ELSEIF(acccalvar(ia)%areaagg==2)THEN
        IF(writematlab)THEN
          WRITE(fileunit_temp,'(a56,a2,a15,a4,a4,a4,a8,a6)')    &
                '%Subbasin assessment. Criteria is calculated for period ',outperiodname(calvarper),'. Variables: rg',outvarid(acccalvar(ia)%rec)%shortname,', rg',outvarid(acccalvar(ia)%comp)%shortname,' Unit: ',outvarid(acccalvar(ia)%comp)%shortunit
          WRITE(fileunit_temp,'(a)') '%OUTREGID'//CHAR(9)//'NSE'//CHAR(9)//'CC'//CHAR(9)//'RE(%)'//CHAR(9)//'RSDE(%)'//CHAR(9)//'Sim'//CHAR(9)//'Rec'//CHAR(9)//'SDSim'//CHAR(9)//'SDRec'//CHAR(9)//'MAE'//CHAR(9)//'RMSE'//CHAR(9)//'Bias'//CHAR(9)//'SDE'//CHAR(9)//'KGE'//CHAR(9)//'KGESD'//CHAR(9)//'KGEM'//CHAR(9)//'NRMSE'//CHAR(9)//'NSEW'
        ELSE
          WRITE(fileunit_temp,'(a57,a2,a15,a4,a4,a4,a8,a6)')    &
                '!!Subbasin assessment. Criteria is calculated for period ',outperiodname(calvarper),'. Variables: rg',outvarid(acccalvar(ia)%rec)%shortname,', rg',outvarid(acccalvar(ia)%comp)%shortname,' Unit: ',outvarid(acccalvar(ia)%comp)%shortunit
          WRITE(fileunit_temp,'(a)') 'OUTREGID'//CHAR(9)//'NSE'//CHAR(9)//'CC'//CHAR(9)//'RE(%)'//CHAR(9)//'RSDE(%)'//CHAR(9)//'Sim'//CHAR(9)//'Rec'//CHAR(9)//'SDSim'//CHAR(9)//'SDRec'//CHAR(9)//'MAE'//CHAR(9)//'RMSE'//CHAR(9)//'Bias'//CHAR(9)//'SDE'//CHAR(9)//'KGE'//CHAR(9)//'KGESD'//CHAR(9)//'KGEM'//CHAR(9)//'NRMSE'//CHAR(9)//'NSEW'
        ENDIF
      ENDIF
      
      !> \li Write data to file
      ndec=3
      IF(acccalvar(ia)%areaagg==0)THEN
        DO isub=1,nsub
          x=sas2(isub,1:maxsubass,ia)
          IF(x(1)/=missing_value) CALL write_dataline(fileunit_temp,maxsubass,x,ndec,0,0,CHAR(9),0,writematlab,id=basin(isub)%subid)
        ENDDO
      ELSEIF(acccalvar(ia)%areaagg==1)THEN
        !Add upstream variables?
      ELSEIF(acccalvar(ia)%areaagg==2)THEN
        DO isub=1,noutreg
          x=sas2(isub,1:maxsubass,ia)
          IF(x(1)/=missing_value) CALL write_dataline(fileunit_temp,maxsubass,x,ndec,0,0,CHAR(9),0,writematlab,id=outregion(isub)%outregid)
        ENDDO
      ENDIF
      CLOSE(fileunit_temp)
    ENDDO

  END SUBROUTINE write_subbasin_assessment

  !>Writes the simulation assessment to hyss file and separate file
  !--------------------------------------------------------------------
  SUBROUTINE write_simulation_assessment(dir,iens,n,optcrit,performance,runens,ccrit,cthres)

    USE WORLDVAR, ONLY : acccalvar,   &
                         calvarper,   &
                         resultseq,   &
                         maxperf,     &
                         simsequence, &
                         outperiodname, &
                         i_rnse,i_snse,  &
                         i_mnse,         &
                         i_rmae,i_sre,   &  
                         i_rre,i_mre,    &
                         i_rra,i_sra,    &
                         i_mra,i_tau,    &
                         i_mdnse,i_mdra, &
                         i_mstdre,i_mcc, &
                         i_mdkg,i_mabsre, &
                         i_mnrmse,i_mnw, &
                         i_numrc,i_nummc, &
                         fileunit_temp
    USE MODVAR, ONLY :   outvarid

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir           !<File directory
    INTEGER, INTENT(IN)          :: iens          !<Current simulation
    INTEGER, INTENT(IN)          :: n             !<Dimension of performance measures
    REAL, INTENT(IN)             :: optcrit       !<Opimization criterion value
    REAL, INTENT(IN)             :: performance(maxperf,n)  !Perfomance criteria
    LOGICAL, INTENT(IN)          :: runens        !<Flag for ensemble simulation
    REAL, OPTIONAL, INTENT(IN)   :: ccrit         !<conditional criteria
    REAL, OPTIONAL, INTENT(IN)   :: cthres        !<conditional criteria threshold
 
    !Local variables
    INTEGER i
    INTEGER nn    !file suffix number
    CHARACTER(LEN=20) filename

    !Write simulation assessment to log-file
    WRITE(6,*)
    IF(PRESENT(cthres))THEN
       WRITE(6,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: ',ccrit,', threshold: ', cthres
    ELSE
       WRITE(6,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: -9999, threshold: -9999'
    ENDIF
    DO i = 1, n
       WRITE(6,*)
       !WRITE(6,*) 'Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname
       IF(acccalvar(i)%areaagg==0)THEN
         WRITE(6,*) 'Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname
       ELSEIF(acccalvar(i)%areaagg==1)THEN
         WRITE(6,*) 'Variables: up',outvarid(acccalvar(i)%rec)%shortname,', up',outvarid(acccalvar(i)%comp)%shortname
       ELSEIF(acccalvar(i)%areaagg==2)THEN
         WRITE(6,*) 'Variables: rg',outvarid(acccalvar(i)%rec)%shortname,', rg',outvarid(acccalvar(i)%comp)%shortname
       ENDIF
       WRITE(6,*) 'Period: ',outperiodname(calvarper)
       WRITE(6,*) 'Regional NSE:',performance(i_rnse,i)
       WRITE(6,*) 'Regional  RA:',performance(i_rra,i)
       WRITE(6,*) 'Regional  RE:',performance(i_rre,i)
       WRITE(6,*) 'Regional MAE:',performance(i_rmae,i)
       WRITE(6,*) 'Average  NSE:',performance(i_mnse,i)
       WRITE(6,*) 'Average   RA:',performance(i_mra,i)
       WRITE(6,*) 'Average   RE:',performance(i_mre,i)
       WRITE(6,*) 'Average RSDE:',performance(i_mstdre,i)
       WRITE(6,*) 'Average   CC:',performance(i_mcc,i)
       WRITE(6,*) 'Average  ARE:',performance(i_mabsre,i)
       WRITE(6,*) 'Spatial  NSE:',performance(i_snse,i)
       WRITE(6,*) 'Spatial   RA:',performance(i_sra,i)
       WRITE(6,*) 'Spatial   RE:',performance(i_sre,i)
       WRITE(6,*) 'Kendalls Tau:',performance(i_tau,i)
       WRITE(6,*) 'Median   NSE:',performance(i_mdnse,i)
       WRITE(6,*) 'Median    RA:',performance(i_mdra,i)
       WRITE(6,*) 'Median   KGE:',performance(i_mdkg,i)
       WRITE(6,*) 'Median NRMSE:',performance(i_mnrmse,i)
       WRITE(6,*) 'Mean    NSEW:',performance(i_mnw,i)
       WRITE(6,*) 'Number of data for regional criterion:',performance(i_numrc,i)
       WRITE(6,*) 'Number of areas in mean/median criterion:',performance(i_nummc,i)
    ENDDO
    WRITE(6,*)

    !Open file
    IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. simass_005.txt
       IF(resultseq.AND.simsequence>0) nn=simsequence
       IF(runens) nn=iens
       filename(1:6) = 'simass'
       filename(7:7) = '_'
       WRITE(filename(8:10),601) nn
       filename(11:14) = '.txt'
    ELSE
       filename='simass.txt'
    ENDIF
601 FORMAT(I3.3)
    OPEN(FILE=TRIM(dir)//filename,UNIT=fileunit_temp,STATUS='unknown',FORM='formatted')
    !Write heading
    WRITE(fileunit_temp,*) 'Simulation assessment all variable used in criterium'

    !Write data to simass-file
    WRITE(fileunit_temp,*) 'Simulation number: ',iens
    IF(PRESENT(cthres))THEN
      WRITE(fileunit_temp,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: ',ccrit,', threshold: ', cthres
    ELSE
      WRITE(fileunit_temp,*) 'Total criteria value: ',optcrit, ' Conditional criteria value: -9999, threshold: -9999'
    ENDIF
    DO i = 1, n
       WRITE(fileunit_temp,*)
       IF(acccalvar(i)%areaagg==0)THEN
         WRITE(fileunit_temp,*) 'Variables: ',outvarid(acccalvar(i)%rec)%shortname,', ',outvarid(acccalvar(i)%comp)%shortname
       ELSEIF(acccalvar(i)%areaagg==1)THEN
         WRITE(fileunit_temp,*) 'Variables: up',outvarid(acccalvar(i)%rec)%shortname,', up',outvarid(acccalvar(i)%comp)%shortname
       ELSEIF(acccalvar(i)%areaagg==2)THEN
         WRITE(fileunit_temp,*) 'Variables: rg',outvarid(acccalvar(i)%rec)%shortname,', rg',outvarid(acccalvar(i)%comp)%shortname
       ENDIF
       WRITE(fileunit_temp,*) 'Period: ',outperiodname(calvarper)
       WRITE(fileunit_temp,*) 'Regional NSE:',performance(i_rnse,i)
       WRITE(fileunit_temp,*) 'Regional  RA:',performance(i_rra,i)
       WRITE(fileunit_temp,*) 'Regional  RE:',performance(i_rre,i)
       WRITE(fileunit_temp,*) 'Regional MAE:',performance(i_rmae,i)
       WRITE(fileunit_temp,*) 'Average  NSE:',performance(i_mnse,i)
       WRITE(fileunit_temp,*) 'Average   RA:',performance(i_mra,i)
       WRITE(fileunit_temp,*) 'Average   RE:',performance(i_mre,i)
       WRITE(fileunit_temp,*) 'Average RSDE:',performance(i_mstdre,i)
       WRITE(fileunit_temp,*) 'Average   CC:',performance(i_mcc,i)
       WRITE(fileunit_temp,*) 'Average  ARE:',performance(i_mabsre,i)
       WRITE(fileunit_temp,*) 'Spatial  NSE:',performance(i_snse,i)
       WRITE(fileunit_temp,*) 'Spatial   RA:',performance(i_sra,i)
       WRITE(fileunit_temp,*) 'Spatial   RE:',performance(i_sre,i)
       WRITE(fileunit_temp,*) 'Kendalls Tau:',performance(i_tau,i)
       WRITE(fileunit_temp,*) 'Median   NSE:',performance(i_mdnse,i)
       WRITE(fileunit_temp,*) 'Median    RA:',performance(i_mdra,i)
       WRITE(fileunit_temp,*) 'Median   KGE:',performance(i_mdkg,i)
       WRITE(fileunit_temp,*) 'Median NRMSE:',performance(i_mnrmse,i)
       WRITE(fileunit_temp,*) 'Mean    NSEW:',performance(i_mnw,i)
       WRITE(fileunit_temp,*) 'Number of data for regional criterion:',performance(i_numrc,i)
       WRITE(fileunit_temp,*) 'Number of areas in mean/median criterion:',performance(i_nummc,i)
    ENDDO
    WRITE(fileunit_temp,*)
    CLOSE(fileunit_temp)
    
  END SUBROUTINE write_simulation_assessment

  !>Manage output files; opens files and write heading
  !--------------------------------------------------------------------
  SUBROUTINE prepare_outputfiles(dir,n,na,iens,runens,allsim)

    USE MODELMODULE, ONLY : open_modeldefined_outputfiles
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble simulation
    LOGICAL, INTENT(IN) :: allsim       !<Flag for writing all simulation results

    CALL check_files_for_doubles()
    CALL open_subbasinfiles(dir,n,iens,runens,allsim)
    CALL open_timefiles(dir,n,iens,runens,allsim)
    CALL open_modeldefined_outputfiles(dir,n,na,iens,runens)

  END SUBROUTINE prepare_outputfiles
  
  !>Manage output files; close files
  !--------------------------------------------------------------------
  SUBROUTINE close_outputfiles(n,na,iens)

    USE MODELMODULE, ONLY : close_modeldefined_outputfiles
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: na           !<Number of aquifers
    INTEGER, INTENT(IN) :: iens         !<Current simulation

    CALL close_subbasinfiles(n,iens)
    CALL close_timefiles(iens)
    CALL close_modeldefined_outputfiles(n,na,iens)

  END SUBROUTINE close_outputfiles
  
  !>Opens files for subbasin and output region printout
  !--------------------------------------------------------------------
  SUBROUTINE open_subbasinfiles(dir,n,iens,runens,allens)

    USE WORLDVAR, ONLY : output,    &
                         noutput,   &
                         outregion, &
                         noutreg,   &
                         outvarinfo, &
                         outperiodname, &
                         writematlab,     &
                         resultseq,       &
                         simsequence,     &
                         fileunit_add_sub,  &
                         fileunit_add_dabout
    USE MODVAR, ONLY : outvarid,     &
                       basin, &
                       nsub

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<Result file directory
    INTEGER, INTENT(IN) :: n            !<Number of subbasins
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble simulation
    LOGICAL, INTENT(IN) :: allens       !<Flag for writing all ensemble results

    !Local variables
    INTEGER i,io,ivar,isb,ireg
    INTEGER nn                  !file suffix number
    INTEGER nselected           !number of selected subbasins for output
    INTEGER funit               !fileunit for subbasin output
    INTEGER ispecific           !order number of specific output type
    CHARACTER(LEN=23) str       !filename
    CHARACTER(LEN=4) var
    CHARACTER(LEN=6) var2
    CHARACTER(LEN=6) varunit
    INTEGER lt,ls,lout1,lout2
    CHARACTER (LEN=16)   t,s
    CHARACTER (LEN=2000) outtxt,outtxt2           !max_outvar*7


    !Open subbasin files
    ispecific = 0
    DO io=1,noutput
      IF(output(io)%fileformat==1)THEN
        ispecific = ispecific + 1
        !Preparation of heading
        IF(writematlab)THEN
          lout1 = 6
          lout2 = 7
          outtxt(1:lout1) = '%DATE'//CHAR(9)
          outtxt2(1:lout2) = '%UNITS'//CHAR(9)
        ELSE
          lout1 = 5
          lout2 = 6
          outtxt(1:lout1) = 'DATE'//CHAR(9)
          outtxt2(1:lout2) = 'UNITS'//CHAR(9)
        ENDIF
        DO ivar = 1,output(io)%nvar
          var  = outvarid(output(io)%variable(ivar)%idindex)%shortname
          varunit = outvarid(output(io)%variable(ivar)%idindex)%shortunit
          IF(output(io)%variable(ivar)%areaagg==1)THEN    !upstream output
            var2='up'//var
          ELSE
            var2=var
          ENDIF
          WRITE(t,'(a16)') var2
          WRITE(s,'(a16)') varunit
          t = ADJUSTL(t)
          s = ADJUSTL(s)
          lt = LEN_TRIM(t)
          ls = LEN_TRIM(s)
          outtxt(lout1+1:lout1+lt) = t(1:lt)
          outtxt2(lout2+1:lout2+ls) = s(1:ls)
          lout1 = lout1+lt
          lout2 = lout2+ls
          IF(ivar < output(io)%nvar) THEN
            outtxt(lout1+1:lout1+1) = CHAR(9)
            outtxt2(lout2+1:lout2+1) = CHAR(9) !Also the last one necessary
            lout1 = lout1+1                    !for reading in "free format"
            lout2 = lout2+1
          ENDIF
        ENDDO

        nselected = n
        IF(output(io)%narea>0) nselected = output(io)%narea
        DO i = 1,nselected
          isb = i
          IF(output(io)%narea>0) isb = output(io)%areaindex(i)
          IF(resultseq.AND.simsequence>0 .OR. runens)THEN
            IF(resultseq.AND.simsequence>0) nn=simsequence
            IF(runens) nn=iens
            CALL create_filename_for_basin(str,basin(isb)%subid,ispecific>1,outperiodname(output(io)%period),nn,allens)
          ELSE
            CALL create_filename_for_basin(str,basin(isb)%subid,ispecific>1,outperiodname(output(io)%period))
          ENDIF
          funit=fileunit_add_sub+nsub*(io-1)+i+(iens-1)*fileunit_add_dabout
          !Open and write heading to file
          OPEN(UNIT=funit,FILE=TRIM(dir)//TRIM(str),status='unknown',form='formatted')
          WRITE(funit,'(a)') outtxt(1:lout1)
          WRITE(funit,'(a)') outtxt2(1:lout2)
        ENDDO
      ENDIF
    ENDDO
    

    !Output region files
    ispecific = 0
    DO io=1,noutput
      IF(output(io)%fileformat==4)THEN
        ispecific = ispecific + 1
        !Preparation of heading of region files
        outtxt='';outtxt2=''
        IF(writematlab)THEN
          lout1 = 6
          lout2 = 7
          outtxt(1:lout1) = '%DATE'//CHAR(9)
          outtxt2(1:lout2) = '%UNITS'//CHAR(9)
        ELSE
          lout1 = 5
          lout2 = 6
          outtxt(1:lout1) = 'DATE'//CHAR(9)
          outtxt2(1:lout2) = 'UNITS'//CHAR(9)
        ENDIF
        DO ivar = 1,output(io)%nvar
          var  = outvarid(outvarinfo(output(io)%variable(ivar)%ovindex)%idindex)%shortname
          varunit = outvarid(outvarinfo(output(io)%variable(ivar)%ovindex)%idindex)%shortunit
          var2='rg'//var
          WRITE(t,'(a16)') var2
          WRITE(s,'(a16)') varunit
          t = ADJUSTL(t)
          s = ADJUSTL(s)
          lt = LEN_TRIM(t)
          ls = LEN_TRIM(s)
          outtxt(lout1+1:lout1+lt) = t(1:lt)
          outtxt2(lout2+1:lout2+ls) = s(1:ls)
          lout1 = lout1+lt
          lout2 = lout2+ls
          IF(ivar < output(io)%nvar) THEN
            outtxt(lout1+1:lout1+1) = CHAR(9)
            outtxt2(lout2+1:lout2+1) = CHAR(9) !Also the last one necessary
            lout1 = lout1+1                    !for reading in "free format"
            lout2 = lout2+1
          ENDIF
        ENDDO

        nselected = noutreg
        IF(output(io)%narea>0) nselected = output(io)%narea
        DO i = 1,nselected
          ireg = i
          IF(output(io)%narea>0) ireg = output(io)%areaindex(i)
          IF(resultseq.AND.simsequence>0 .OR. runens)THEN
            IF(resultseq.AND.simsequence>0) nn=simsequence
            IF(runens) nn=iens
            CALL create_filename_for_basin(str,outregion(ireg)%outregid,ispecific>1,outperiodname(output(io)%period),nn,allens)
          ELSE
            CALL create_filename_for_basin(str,outregion(ireg)%outregid,ispecific>1,outperiodname(output(io)%period))
          ENDIF
          funit=fileunit_add_sub+nsub*(io-1)+i+(iens-1)*fileunit_add_dabout
          !Open and write heading to file
          OPEN(UNIT=funit,FILE=TRIM(dir)//TRIM(str),status='unknown',form='formatted')
          WRITE(funit,'(a)') outtxt(1:lout1)
          WRITE(funit,'(a)') outtxt2(1:lout2)
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE open_subbasinfiles

  !>Calculate and write subbasin files
  !--------------------------------------------------------------------
  SUBROUTINE write_subbasinfiles(io,idt,ndt,iens,cd)

    USE LIBDATE, ONLY : DateType
    USE WORLDVAR, ONLY : output, &
                         writematlab,     &
                         outstartdate,    &
                         fileunit_add_sub, &
                         fileunit_add_dabout
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_basinoutput
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    INTEGER, INTENT(IN) :: io           !<Current output
    INTEGER, INTENT(IN) :: idt          !<Current time
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: cd    !<timestep date-time

    !Local variables
    INTEGER ia,n   !loop index over selected subbasins and number of selected basins
    LOGICAL pwrite              !Flag for periodend, time to write to file
    REAL,ALLOCATABLE :: x(:)    !Help variable for continous print out

    IF(output(io)%nvar>0)THEN  !Write period mean for selected subbasins
      ALLOCATE(x(output(io)%nvar))
      IF(output(io)%narea>0)THEN
        n = output(io)%narea
      ELSE
        n = nsub
      ENDIF
      DO ia = 1,n
        CALL compute_basinoutput(cd,io,ia,output(io)%nvar,x,pwrite,idt,ndt,1)
        IF(pwrite) CALL write_dataline(fileunit_add_sub+nsub*(io-1)+ia+(iens-1)*fileunit_add_dabout,output(io)%nvar,x,output(io)%decimal,  &
                        output(io)%signfig,output(io)%period,CHAR(9),0,writematlab,d=cd,odate=outstartdate)
       ENDDO
     ENDIF
     IF(ALLOCATED(x)) DEALLOCATE(x)

  END SUBROUTINE write_subbasinfiles

  !>Calculate and write subbasin files (special for assimilation)
  !>Writes sequences' files in parallel
  !--------------------------------------------------------------------
  SUBROUTINE write_subbasinfiles_in_parallel(io,idt,ndt,iens,cd)

    USE LIBDATE, ONLY : DateType
    USE WORLDVAR, ONLY : output,      &
                         writematlab,     &
                         outstartdate,    &
                         fileunit_add_sub, &
                         fileunit_add_dabout
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_basinoutput
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    INTEGER, INTENT(IN) :: io           !<Current output
    INTEGER, INTENT(IN) :: idt          !<Current time
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: cd     !<timestep date-time

    !Local variables
    INTEGER ia,n
    LOGICAL pwrite              !Flag for periodend, time to write to file
    REAL,ALLOCATABLE :: x(:)    !Help variable for continous print out

    IF(output(io)%nvar>0)THEN  !Write period mean for selected subbasins
      ALLOCATE(x(output(io)%nvar))
      IF(output(io)%narea>0)THEN
        n = output(io)%narea
      ELSE
        n = nsub
      ENDIF
      DO ia = 1,n
        CALL compute_basinoutput(cd,io,ia,output(io)%nvar,x,pwrite,idt,ndt,iens)
        IF(pwrite) CALL write_dataline(fileunit_add_sub+nsub*(io-1)+ia+(iens-1)*fileunit_add_dabout,output(io)%nvar,x,output(io)%decimal,  &
                        output(io)%signfig,output(io)%period,CHAR(9),0,writematlab,d=cd,odate=outstartdate)
       ENDDO
     ENDIF
     IF(ALLOCATED(x)) DEALLOCATE(x)

  END SUBROUTINE write_subbasinfiles_in_parallel

  !>Calculate and write output region files
  !--------------------------------------------------------------------
  SUBROUTINE write_regionfiles(io,idt,ndt,iens,cd)

    USE LIBDATE, ONLY : DateType
    USE WORLDVAR, ONLY : output, &
                         noutreg,       &
                         writematlab,     &
                         outstartdate,    &
                         fileunit_add_sub, &
                         fileunit_add_dabout
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_regionoutput
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    INTEGER, INTENT(IN) :: io           !<Current output
    INTEGER, INTENT(IN) :: idt          !<Current time step
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time steps
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: cd    !<Current date-time

    !Local variables
    INTEGER i,ir,nselected
    LOGICAL pwrite              !Flag for period end, time to write to file
    REAL,ALLOCATABLE :: x(:)    !Data to be printed

    ALLOCATE(x(output(io)%nvar))
    nselected = noutreg
    IF(output(io)%narea>0) nselected = output(io)%narea
    DO i = 1,nselected
      ir = i
      IF(output(io)%narea>0) ir = output(io)%areaindex(i)
      CALL compute_regionoutput(cd,io,ir,output(io)%nvar,idt,ndt,1,pwrite,x)
      IF(pwrite) CALL write_dataline(fileunit_add_sub+nsub*(io-1)+i+(iens-1)*fileunit_add_dabout,output(io)%nvar,x,output(io)%decimal,  &
                      output(io)%signfig,output(io)%period,CHAR(9),0,writematlab,d=cd,odate=outstartdate)
    ENDDO
    IF(ALLOCATED(x)) DEALLOCATE(x)

  END SUBROUTINE write_regionfiles

  !>Calculate and write output region files
  !--------------------------------------------------------------------
  SUBROUTINE write_regionfiles_in_parallel(io,idt,ndt,iens,cd)

    USE LIBDATE, ONLY : DateType
    USE WORLDVAR, ONLY : output, &
                         noutreg,       &
                         writematlab,     &
                         outstartdate,    &
                         fileunit_add_sub, &
                         fileunit_add_dabout
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_regionoutput
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    INTEGER, INTENT(IN) :: io           !<Current output
    INTEGER, INTENT(IN) :: idt          !<Current time step
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time steps
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: cd    !<Current date-time

    !Local variables
    INTEGER i,ir,nselected
    LOGICAL pwrite              !Flag for period end, time to write to file
    REAL,ALLOCATABLE :: x(:)    !Data to be printed

    ALLOCATE(x(output(io)%nvar))
    nselected = noutreg
    IF(output(io)%narea>0) nselected = output(io)%narea
    DO i = 1,nselected
      ir = i
      IF(output(io)%narea>0) ir = output(io)%areaindex(i)
      CALL compute_regionoutput(cd,io,ir,output(io)%nvar,idt,ndt,iens,pwrite,x)
      IF(pwrite) CALL write_dataline(fileunit_add_sub+nsub*(io-1)+i+(iens-1)*fileunit_add_dabout,output(io)%nvar,x,output(io)%decimal,  &
                      output(io)%signfig,output(io)%period,CHAR(9),0,writematlab,d=cd,odate=outstartdate)
    ENDDO
    IF(ALLOCATED(x)) DEALLOCATE(x)

  END SUBROUTINE write_regionfiles_in_parallel

  !>Close files for subbasin and output region printout
  !--------------------------------------------------------------------
  SUBROUTINE close_subbasinfiles(n,iens)
  
    USE WORLDVAR, ONLY : output, &
                         fileunit_add_sub, &
                         fileunit_add_dabout, &
                         noutreg,noutput
    
    !Argument declarations
    INTEGER, INTENT(IN) :: n      !<Number of subbasins
    INTEGER, INTENT(IN) :: iens   !<Current simulation

    !Local variables
    INTEGER i,io,funit,nselected

    DO io=1,noutput
      IF(output(io)%fileformat==1)THEN
        nselected = n
        IF(output(io)%narea>0) nselected = output(io)%narea
        DO i = 1,nselected
          funit=fileunit_add_sub+n*(io-1)+i+(iens-1)*fileunit_add_dabout   !denna funit behöver också io när subbasinfiler blir fler
          CLOSE(funit)
        ENDDO
      ELSEIF(output(io)%fileformat==4)THEN
        nselected = noutreg
        IF(output(io)%narea>0) nselected = output(io)%narea
        DO i = 1,nselected
          funit=fileunit_add_sub+n*(io-1)+i+(iens-1)*fileunit_add_dabout
          CLOSE(funit)
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE close_subbasinfiles

  !>Check time- and map-files for multiple variables
  !--------------------------------------------------------------------
  SUBROUTINE check_files_for_doubles()
  
    USE WORLDVAR, ONLY : output,noutput
    USE MODVAR, ONLY : outvarid

    !Local variables
    INTEGER io,ivar,ivar2

    !Check for double output
    DO io = 1, noutput
      IF(output(io)%fileformat==2 .OR. output(io)%fileformat==3)THEN
        DO ivar = 1,output(io)%nvar
          DO ivar2 = ivar+1,output(io)%nvar
            IF(outvarid(output(io)%variable(ivar)%idindex)%shortname==outvarid(output(io)%variable(ivar2)%idindex)%shortname &
                 .AND. output(io)%variable(ivar)%areaagg==output(io)%variable(ivar2)%areaagg)THEN
              WRITE(6,*) 'ERROR: Same variable multiple times in output. Variable: '//TRIM(outvarid(output(io)%variable(ivar)%idindex)%shortname)
              WRITE(6,*) 'Check timeoutput and/or mapoutput'
              STOP 1
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE check_files_for_doubles

  !>Opens files for timeserie printout
  !--------------------------------------------------------------------
  SUBROUTINE open_timefiles(dir,n,iens,runens,allens)
  
    USE WORLDVAR, ONLY : outregion,noutreg, &
                         outperiodname,  &
                         writematlab,       &
                         resultseq,         &
                         simsequence,       &
                         i_t,i_h,i_d,i_w,i_m,i_y,       &
                         fileunit_add_time, &
                         fileunit_add_timeio, &
                         output,noutput
    USE MODVAR, ONLY : outvarid,           &
                       i_sum,              &
                       basin
    USE READWRITE_ROUTINES, ONLY : write_integer_header, &
                                   write_commentline

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir !<File directory
    INTEGER, INTENT(IN) :: n            !<Number of data columns (subbasins)
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    LOGICAL, INTENT(IN) :: runens       !<Flag for ensemble run
    LOGICAL, INTENT(IN) :: allens       !<Flag for writing all ensemble results
    
    !Local variables
    INTEGER io,ivar
    INTEGER nn                  !file suffix2 number
    INTEGER funit
    INTEGER ispecific
    CHARACTER(LEN=22) filename
    CHARACTER(LEN=20) var
    CHARACTER(LEN=40) var2
    CHARACTER(LEN=30) unit
    INTEGER lout
    CHARACTER (LEN=200) commentheading
    CHARACTER (LEN=200) outtxt

    !Open files
    ispecific = 0
    DO io = 1, noutput
      IF(output(io)%fileformat==3)THEN
        ispecific = ispecific + 1
        DO ivar = 1,output(io)%nvar
          !Find variable and file
          IF(resultseq.AND.simsequence>0 .OR. runens)THEN  !Ex. timeCOUT_005.txt
            IF(resultseq.AND.simsequence>0) nn=simsequence
            IF(runens) nn=iens
            IF(ispecific>1)THEN
              CALL create_filename_for_variable(filename,'time',outvarid(output(io)%variable(ivar)%idindex)%shortname,output(io)%variable(ivar)%areaagg,suffix1=outperiodname(output(io)%period),suffix2=nn,suffix2format=allens)
            ELSE
              CALL create_filename_for_variable(filename,'time',outvarid(output(io)%variable(ivar)%idindex)%shortname,output(io)%variable(ivar)%areaagg,suffix2=nn,suffix2format=allens)
            ENDIF
          ELSE  !Ex. timeCOUT.txt
            IF(ispecific>1)THEN
              CALL create_filename_for_variable(filename,'time',outvarid(output(io)%variable(ivar)%idindex)%shortname,output(io)%variable(ivar)%areaagg,suffix1=outperiodname(output(io)%period))
            ELSE
              CALL create_filename_for_variable(filename,'time',outvarid(output(io)%variable(ivar)%idindex)%shortname,output(io)%variable(ivar)%areaagg)
            ENDIF
          ENDIF

          var=outvarid(output(io)%variable(ivar)%idindex)%longname
          IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
            var2='regional average of '//var
          ELSEIF(output(io)%variable(ivar)%areaagg==1)THEN    !upstream output
            var2='upstream average of '//var
          ELSE
            var2=var
          ENDIF
          unit=outvarid(output(io)%variable(ivar)%idindex)%longunit
          IF(outvarid(output(io)%variable(ivar)%idindex)%vartype==i_sum)THEN
            SELECT CASE(output(io)%period)
            CASE(i_h)
              unit=TRIM(unit)//' per hour'
            CASE(i_t)
              unit=TRIM(unit)//' per timestep'
            CASE(i_d)
              unit=TRIM(unit)//' per day'   
            CASE(i_w)
              unit=TRIM(unit)//' per week'
            CASE(i_m)
              unit=TRIM(unit)//' per month'
            CASE(i_y)
              unit=TRIM(unit)//' per year'
            END SELECT
          ENDIF
          funit = iens*fileunit_add_time+ivar+(io-1)*fileunit_add_timeio

          !Open file
          OPEN(UNIT=funit,FILE=TRIM(dir)//TRIM(filename),STATUS='unknown',FORM='formatted',ERR=900)

          !Write headings
          IF(writematlab)THEN
            commentheading = '%This is a file with timeseries of '//TRIM(var2)//' in '//TRIM(unit)
            CALL write_commentline(funit,commentheading)
            outtxt(1:6) = '%DATE'//CHAR(9)
            lout = 6
            IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
              CALL write_integer_header(funit,noutreg,outregion(:)%outregid,outtxt(1:lout))
            ELSE
              CALL write_integer_header(funit,n,basin(:)%subid,outtxt(1:lout))
            ENDIF
          ELSE
            commentheading = '!!This is a file with timeseries of '//TRIM(var2)//' in '//TRIM(unit)
            CALL write_commentline(funit,commentheading)
            outtxt(1:5) = 'DATE'//CHAR(9)
            lout = 5
            IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
              CALL write_integer_header(funit,noutreg,outregion(:)%outregid,outtxt(1:lout))
            ELSE
              CALL write_integer_header(funit,n,basin(:)%subid,outtxt(1:lout))
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
    RETURN
    
900 WRITE(6,*) 'ERROR: Open file '//TRIM(dir)//TRIM(filename)//' for writing.'
    STOP 1 

  END SUBROUTINE open_timefiles

  !>Create file name string for subbasin- or outregion file
  !--------------------------------------------------------------------
  SUBROUTINE create_filename_for_basin(filename,basinid,usesuffix1,suffix1,suffix2,suffix2format)
  
    !Argument declarations
    CHARACTER(LEN=*),INTENT(OUT) :: filename
    INTEGER, INTENT(IN) :: basinid  !<subid or outregid; to be name of file
    LOGICAL,INTENT(IN) :: usesuffix1      !<status of using suffix1
    CHARACTER(LEN=*),INTENT(IN) :: suffix1   !<string file suffix
    INTEGER,OPTIONAL,INTENT(IN) :: suffix2            !<integer file suffix
    LOGICAL,OPTIONAL,INTENT(IN) :: suffix2format      !<use 6 character suffix instead of three
  
    INTEGER iout
  
    filename = ''
    WRITE(filename(1:7),600) basinid
    iout=7
    IF(usesuffix1)THEN
      filename(iout+1:iout+1) = '_'
      filename(iout+2:iout+3) = TRIM(ADJUSTL(suffix1))
      iout=iout+3
    ENDIF
    IF(PRESENT(suffix2))THEN
      filename(iout+1:iout+1)='_'
      IF(PRESENT(suffix2format))THEN
        IF(suffix2format)THEN
          WRITE(filename(iout+2:iout+8),602) suffix2   !6 character sequence number
          iout = iout+8
        ELSE
          WRITE(filename(iout+2:iout+4),601) suffix2
          iout = iout+4
        ENDIF
      ELSE
        WRITE(filename(iout+2:iout+4),601) suffix2
        iout = iout+4
      ENDIF
    ENDIF
    filename(iout+1:iout+4)='.txt'
    iout = iout+4

600 FORMAT(I7.7)
601 FORMAT(I3.3)
602 FORMAT(I7.7)
       
  END SUBROUTINE create_filename_for_basin
  
  !>Create file name string for time- or map file
  !--------------------------------------------------------------------
  SUBROUTINE create_filename_for_variable(filename,outputtype,variablename,areaagg,suffix1,suffix2,suffix2format)
  
  USE CONVERT, ONLY : scalar_upper_case
  
    !Argument declarations
    CHARACTER(LEN=*),INTENT(OUT) :: filename
    CHARACTER(LEN=*),INTENT(IN) :: outputtype  !<time- or map-file
    CHARACTER(LEN=4),INTENT(IN) :: variablename !<code for variable; to be name of file
    INTEGER, INTENT(IN) :: areaagg  !<subbasin, upstream or output region variable
    CHARACTER(LEN=*),OPTIONAL,INTENT(IN) :: suffix1   !<string file suffix
    INTEGER,OPTIONAL,INTENT(IN) :: suffix2            !<integer file suffix
    LOGICAL,OPTIONAL,INTENT(IN) :: suffix2format      !<use 6 character suffix instead of three
  
    INTEGER iout
    CHARACTER(LEN=4) :: capitalname
  
    filename = ''
    iout = LEN(outputtype)
    filename(1:iout)=outputtype(1:iout)
    IF(areaagg==2)THEN
      filename(iout+1:iout+2) = 'RG'
      iout=iout+2
    ELSEIF(areaagg==1)THEN
      filename(iout+1:iout+2) = 'UP'
      iout=iout+2
    ENDIF
    capitalname = variablename
    CALL scalar_upper_case(capitalname)
    filename(iout+1:iout+4) = capitalname
    iout = iout+4
    IF(PRESENT(suffix1))THEN
      filename(iout+1:iout+1)='_'
      filename(iout+2:iout+3) = TRIM(ADJUSTL(suffix1))
      iout = iout+3
    ENDIF
    IF(PRESENT(suffix2))THEN
      filename(iout+1:iout+1)='_'
      IF(PRESENT(suffix2format))THEN
        IF(suffix2format)THEN
          WRITE(filename(iout+2:iout+8),602) suffix2   !6 character sequence number
          iout = iout+8
        ELSE
          WRITE(filename(iout+2:iout+4),601) suffix2
          iout = iout+4
        ENDIF
      ELSE
        WRITE(filename(iout+2:iout+4),601) suffix2
        iout = iout+4
      ENDIF
    ENDIF
    filename(iout+1:iout+4)='.txt'
    iout = iout+4

601    FORMAT(I3.3)
602    FORMAT(I7.7)

  END SUBROUTINE create_filename_for_variable
  
  !>Calculate time aggregate of time-output and writes to file
  !--------------------------------------------------------------------
  SUBROUTINE write_timefiles(io,idt,ndt,iens,cd)
  
    USE WORLDVAR, ONLY : writematlab,       &
                         output, &
                         noutreg, &
                         outstartdate,      &
                         fileunit_add_timeio, &
                         fileunit_add_time
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_timeoutput
    USE READWRITE_ROUTINES, ONLY : write_dataline
    USE LIBDATE, ONLY : DateType

    !Argument declarations
    INTEGER, INTENT(IN) :: io           !<Current output
    INTEGER, INTENT(IN) :: idt          !<Current time
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: cd    !<current date-time
    
    !Local variables
    INTEGER ivar,n
    LOGICAL pwrite              !Flag for periodend, time to write to file
    REAL,ALLOCATABLE :: y(:)    !Help variable for continous print out

    IF(.NOT.ALLOCATED(y)) ALLOCATE(y(nsub))
    DO ivar = 1, output(io)%nvar
      IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
        n = noutreg
      ELSE
        n = nsub
      ENDIF
      CALL compute_timeoutput(cd,io,ivar,n,y,pwrite,idt,ndt,1)
      IF(pwrite) CALL write_dataline(iens*fileunit_add_time+ivar+(io-1)*fileunit_add_timeio,n,y,output(io)%decimal,  &
                      output(io)%signfig,output(io)%period,CHAR(9),0,writematlab,d=cd,  &
                      odate=outstartdate)
    ENDDO 
    IF(ALLOCATED(y)) DEALLOCATE(y)

  END SUBROUTINE write_timefiles

  !>Calculate time aggregate of time-output and writes to file (special 
  !>for data assimilation)
  !--------------------------------------------------------------------
  SUBROUTINE write_timefiles_in_parallel(io,idt,ndt,iens,cd)
  
    USE WORLDVAR, ONLY : output,        &
                         writematlab,       &
                         noutreg,           &
                         outstartdate,      &
                         fileunit_add_timeio, &
                         fileunit_add_time
    USE MODVAR, ONLY : nsub
    USE COMPOUT, ONLY : compute_timeoutput
    USE READWRITE_ROUTINES, ONLY : write_dataline
    USE LIBDATE, ONLY : DateType

    !Argument declarations
    INTEGER, INTENT(IN) :: io           !<Current output
    INTEGER, INTENT(IN) :: idt          !<Current time
    INTEGER, INTENT(IN) :: ndt          !<Maximum simulation time
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(DateType), INTENT(IN) :: cd     !<current date-time
    
    !Local variables
    INTEGER ivar,n
    LOGICAL pwrite              !Flag for period end, time to write to file
    REAL,ALLOCATABLE :: y(:)    !Help variable for continous print out

    IF(.NOT.ALLOCATED(y)) ALLOCATE(y(nsub))
    DO ivar = 1, output(io)%nvar
      IF(output(io)%variable(ivar)%areaagg==2)THEN    !region output
        n = noutreg
      ELSE
        n = nsub
      ENDIF
      CALL compute_timeoutput(cd,io,ivar,n,y,pwrite,idt,ndt,iens)
      IF(pwrite) CALL write_dataline(iens*fileunit_add_time+ivar+(io-1)*fileunit_add_timeio,n,y,output(io)%decimal,  &
                      output(io)%signfig,output(io)%period,CHAR(9),0,writematlab,d=cd,  &
                      odate=outstartdate)
    ENDDO 
    IF(ALLOCATED(y)) DEALLOCATE(y)

  END SUBROUTINE write_timefiles_in_parallel

  !>Close files for timeserie output
  !--------------------------------------------------------------------
  SUBROUTINE close_timefiles(iens)
  
    USE WORLDVAR, ONLY : output,  &
                         noutput, &
                         fileunit_add_timeio, &
                         fileunit_add_time

    !Argument declarations
    INTEGER, INTENT(IN) :: iens  !<Current simulation
    
    !Local variables
    INTEGER io,ivar
    INTEGER funit

    DO io=1,noutput
      IF(output(io)%fileformat==3)THEN
        DO ivar = 1, output(io)%nvar
          funit = iens*fileunit_add_time+ivar+(io-1)*fileunit_add_timeio
          CLOSE(funit)
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE close_timefiles

  !>Read parameter values from file
  !>
  !>\b Consequences The modvar module variables soilpar, landpar, genpar,
  !> basinpar, regpar, lakedatapar and monthpar will 
  !> be allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE load_parameters(dir,ns,name) 

    USE WORLDVAR, ONLY : fileunit_temp,     &
                         comment_str
    USE MODVAR, ONLY : soilpar,    &   !OUT, parameters depending on soil type
                       landpar,    &   !OUT, parameters depending on land use
                       genpar,     &   !OUT, parameters not dependent
                       basinpar,   &   !OUT, parameters depending on subbbasin
                       regpar,     &   !OUT, parameters depending on parregion
                       lakedatapar, &  !OUT, parameters from LakeData 
                       monthpar, &   !OUT, parameters depending on month
                       modparid, &   !definition of model parameters
                       regiondivision, &
                       nregions, &
                       max_par, &
                       m_gpar, &
                       m_bpar, &
                       m_spar, &
                       m_lpar, &
                       m_rpar, &
                       m_ldpar, &
                       m_mpar, &
                       nluse,   &
                       nsoil

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir            !<File directory
    INTEGER,          INTENT(IN) :: ns             !<Number of subbasins (nsub_basemodel)
    CHARACTER(LEN=*), INTENT(IN) :: name           !<File name
    
    !Local variables
    CHARACTER (LEN=220) filename
    CHARACTER (LEN=10) varstr       !string with variable name
    CHARACTER(LEN=18000) line
    INTEGER j
    INTEGER dim
    INTEGER dimcheck                !current parameter supposed size
    INTEGER nskip                   !number of lines with headings
    INTEGER nline                   !line number in file
    INTEGER nvalues                 !number of values read from line
    INTEGER varindex                !parameter index in array
    LOGICAL lakedatafile            !Status of LakeData file
    REAL,ALLOCATABLE :: values(:)               !values of parameter read from file

    !>\b Algorithm \n
    !Initialisations
    lakedatafile=.FALSE.
    IF(ALLOCATED(lakedatapar)) lakedatafile=.TRUE.

    !>Open par.txt file and read heading
    filename=TRIM(dir)//TRIM(name)
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')

    nskip = 0
    DO j = 1,nskip   !Skip heading
      READ(fileunit_temp,*)
    ENDDO

    !>Allocate and initiate parameter variables
    dim = MAX(1,ns,nsoil,nluse,MAXVAL(nregions),12)
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(dim))

    !>Read parameter values from file
    nline = nskip + 1
    DO 
      READ(fileunit_temp,'(a)',END=100) line
      IF(line(1:2)==comment_str)CYCLE
      CALL read_parameterline(line,dim,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')EXIT    !end of file

      !Find corresponding varindex
      varindex = 0
      DO j = 1,max_par
        IF(varstr==modparid(j)%shortname)THEN
          IF(modparid(j)%deptype/=m_ldpar)THEN 
            varindex = j
            EXIT
          ENDIF
        ENDIF
      ENDDO
      IF(varindex==0)THEN
        WRITE(6,*) 'Warning - unknown variable on line. Parameter ',TRIM(varstr),' ignored.'
      ELSE
        !Find and check variable dimension and set parameter values
        SELECT CASE(modparid(varindex)%deptype)
        CASE(m_gpar)
          dimcheck = 1                   !genpar
        CASE(m_bpar)
          dimcheck = ns                  !basinpar
        CASE(m_spar)
          dimcheck = nsoil
        CASE(m_lpar)
          dimcheck = nluse
        CASE(m_rpar)
          dimcheck = nregions(regiondivision(modparid(varindex)%parno))
        CASE(m_mpar)
          dimcheck = 12
        CASE DEFAULT
          WRITE(6,*) 'ERROR in code. Unknown parameter dependence type'
          STOP 1
        END SELECT
        IF(dimcheck /= nvalues)THEN
          IF(modparid(varindex)%deptype==m_rpar .and. dimcheck<nvalues)then
            nvalues = dimcheck     !Ugly solution to bug2650(10716), number of regions reduced by selectaro, but not in par.txt
          ELSE
            WRITE(6,*) 'Error - wrong number of values in parameter file'
            WRITE(6,*) 'Error - check parameter ',varstr
            STOP 1
          ENDIF
        ENDIF
        SELECT CASE(modparid(varindex)%deptype)
        CASE(m_gpar)
          genpar(modparid(varindex)%parno) = values(1)
        CASE(m_bpar)
          basinpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
        CASE(m_spar)
          soilpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
        CASE(m_lpar)
          landpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
        CASE(m_rpar)
          regpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
        CASE(m_mpar)
          monthpar(modparid(varindex)%parno,1:nvalues) = values(1:nvalues)
        END SELECT
      ENDIF

      nline = nline + 1
    ENDDO
100 CONTINUE
    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)
    IF(ALLOCATED(values)) DEALLOCATE(values)

    !>Set parameter values for missing lakes with general values
    CALL finish_lakedata_table(lakedatafile,ns)

  END SUBROUTINE load_parameters
  
  !>Reads parameters to be calibrated from file and saves them in
  !>worldvar-arrays.
  !!
  !>\b Consequences Module worldvar variables optparmin, optparmax, optparprecision,
  !> optparid, parindex, optim, dimpar, numoptimpar will be allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE load_optpar(dir) 

    USE WORLDVAR, ONLY : optparmin, optparmax, & 
                         optparprecision, optparid, &
                         dimpar,maxoptpar, &
                         parindex, &
                         optim, &
                         count_optim_par, &
                         numoptimpar, &
                         fileunit_temp
    USE MODVAR, ONLY : modparid,max_par, &
                       nsoil,nluse,nsub,nregions !,nlakeregions

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir    !<File directory

    !Local variables
    INTEGER j, caseFlag
    INTEGER nskip                   !number of lines with headings
    INTEGER nline                   !line number in file
    INTEGER nvalues                 !number of values read from line
    INTEGER numpar                  !number of different parameters to be calibrated (parameter names in optpar.txt)
    INTEGER varindex                !index for parameter in varstr
    REAL,ALLOCATABLE :: values(:)   !values of parameter read from file (different soil or landuse)
    REAL,ALLOCATABLE :: par(:,:)    !Tempary storage of parameter value 
    CHARACTER (LEN=210) filename
    CHARACTER (LEN=10) varstr       !string with parameter name
    CHARACTER(LEN=18000) line       !line in file
    CHARACTER(LEN=2) taskchar       !code for optimation task
    CHARACTER(LEN=1) flagChar       !character string for one-character flags (Y/N)

    !Local parameters
    INTEGER, PARAMETER :: ninfolines = 20

    !Initiation of variables for optimization parameters
    dimpar = MAX(1,nsub,nsoil,nluse,MAXVAL(nregions))
    IF(.NOT.ALLOCATED(optparmin)) ALLOCATE(optparmin(maxoptpar,dimpar)) !worldvar
    IF(.NOT.ALLOCATED(optparmax)) ALLOCATE(optparmax(maxoptpar,dimpar)) !worldvar
    IF(.NOT.ALLOCATED(optparprecision)) ALLOCATE(optparprecision(maxoptpar,dimpar))   !worldvar
    optparmin = 0.
    optparmax = 0.
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(dimpar))           !local
    IF(.NOT.ALLOCATED(par)) ALLOCATE(par(maxoptpar,dimpar))       !local

    filename=TRIM(dir)//'optpar.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
    nskip = 1
    DO j = 1,nskip   !Skip heading
      READ(fileunit_temp,*)
    ENDDO
    nline = nskip + 1

    DO j = 1, ninfolines     !Read lines with optimisation information
      READ(fileunit_temp,'(a)',END=100) line

      IF(line(1:4)=='task')THEN
        READ(line(6:16000),*,END=102) taskchar
        SELECT CASE(taskchar)
        CASE('MC')
          optim%task_MC        = .TRUE.
        CASE('BP')
          optim%task_boundps   = .TRUE.
        CASE('SC')                        ! Scanning mode for 2-parameter problems, execute an organised sampling of crit(param)  [Fred, 06.10.11]
          optim%task_Scanning  = .TRUE.
        CASE('SM')                        ! Stage-wise centering Monte Carlo routine  [Fred, 26.08.10]
          optim%task_stageMC   = .TRUE.
        CASE('BN')                        ! Brent optimisation routine, with new line search  [Fred, 27.11.10]
          optim%task_BrentNew  = .TRUE.
        CASE('SD')                        ! Steepest descent method  [Fred, 08.10.11]
          optim%task_stpstDesc = .TRUE.
        CASE('Q1')                        ! Quasi-Newton DFP gradient-based optimisation routine  [Fred, 14.07.10]
          optim%task_DFP       = .TRUE.
        CASE('Q2')                        ! Quasi-Newton BFGS gradient-based optimisation routine  [Fred, 09.12.10]
          optim%task_BFGS      = .TRUE.
        CASE('WA') 
          optim%task_writeall  = .TRUE.
        CASE('WS') 
          optim%task_writesim  = .TRUE.  !Write simulation result for all ensembles (MC-methods)
        CASE('DE')
          optim%task_DEMC      = .TRUE.   !Differential-Evolution Markov Chain method [David, 2013.02.11]
        END SELECT
      ENDIF

102   IF(line(1:6)=='num_mc')THEN
        READ(line(8:16000),*,END=103) optim%nruns_MC
      ENDIF

103   IF(line(1:7)=='num_ens')THEN
        READ(line(9:16000),*,END=104) optim%nruns_best
      ENDIF

104   IF(line(1:8)=='num_bpmc')THEN
        READ(line(10:16000),*,END=105) optim%nruns_MCI
      ENDIF

105   IF(line(1:9)=='num_bpmax')THEN
        READ(line(11:16000),*,END=106) optim%nruns_MCImax
      ENDIF

106   IF(line(1:10)=='num_stages')THEN                        ! Read 'num_stages' (amount of successive zooming MC stages) for zoom MC routine  [Fred, 26.08.10]
        READ(line(11:16000),*,END=107) optim%nruns_stages
      ENDIF

107   IF(line(1:8)=='num_zoom')THEN                           ! Read numerical parameter 'num_zoom' (zooming factor) for zoom MC routine  [Fred, 26.08.10]
        READ(line(11:16000),*,END=108) optim%nruns_zoom       ! Checked later (in stage MC procedure) that it is not larger than 1
      ENDIF

108   IF(line(1:10)=='num_dbgMod')THEN                        ! Flag for debug scenario (0 = HYPE, see optim.f90 for alternatives)
        READ(line(11:16000),*,END=109) optim%cal_debugCase
      ENDIF

109   IF(line(1:7)=='cal_log')THEN                            ! Y/N flag to write or not the file "calibration.log", that contains all details on the calibration routine
        READ(line(11:16000),*,END=110) flagChar
        IF(flagChar == 'N')THEN
          optim%cal_log = .FALSE.
        ELSEIF(flagChar .NE. 'N' .AND. flagChar .NE. 'Y')THEN
          WRITE(6,*) 'ERROR: flag to enable/disable writing calibration.log must be either Y or N'
          STOP 1
        ENDIF
      ENDIF

      !Numerical parameter for interruption of non-MC methods [Fred; checked for consistency with optim-type on 29.09.11]
110   IF(line(1:10)=='num_maxItr')THEN                        ! Max amount of iterations allowed
        READ(line(11:16000),*,END=111) optim%cal_maxIterat
      ENDIF

111   IF(line(1:10)=='num_maxTim')THEN                        ! Max amout of time (hours) allowed to calibration routine
        READ(line(11:16000),*,END=112) optim%cal_maxTime
      ENDIF

112   IF(line(1:10)=='num_criItr')THEN                        ! Amount of last optimisation iterations taken into account for criteria improvement monitoring
        READ(line(11:16000),*,END=113) optim%cal_improvCritIter
      ENDIF

113   IF(line(1:10)=='num_criTol')THEN                        ! Tolerance to consider criteria as optimised (delta/mean of criteria over the last 'num_critIt' iterations)
        READ(line(11:16000),*,END=114) optim%cal_improvCritTol
      ENDIF

114   IF(line(1:10)=='num_parItr')THEN                        ! Amount of last optimisation iterations taken into account for parameter improvement monitoring
        READ(line(11:16000),*,END=115) optim%cal_improvParamIter
      ENDIF

      !Numerical parameters for Quasi-Newton method [Fred]
115   IF(line(1:10)=='QN_nrmTol')THEN                         ! Tolerance for gradient norm to be considered zero
        READ(line(11:16000),*,END=116) optim%QN_flatTol
      ENDIF

116   IF(line(1:10)=='QN_pctDerv')THEN                        ! Offset (percentage of parameter value) for numerical derivative
        READ(line(11:16000),*,END=117) optim%QN_factorDeriv
      ENDIF

117   IF(line(1:10)=='QN_stencil')THEN                        ! Stencil type
        READ(line(11:16000),*,END=118) optim%QN_stencil       ! QN optimisation algorithm checkes that it is either 2, 4, 6 or 8
      ENDIF

118   IF(line(1:10)=='QN_lambMax')THEN                        ! Factor containing lambda, to avoid taking points for gradient outside allowed parameter space
        READ(line(11:16000),*,END=119) optim%QN_lambdaMaxFac
      ENDIF

119   IF(line(1:10)=='QN_lambAcc')THEN                        ! Factor increasing the step length proposed by QN algorithms (case lambda = 1 to be replaced by lambdaAccel), in order to allow for faster iteration progression
        READ(line(11:16000),*,END=120) optim%QN_lambdaAccel
      ENDIF

120   IF(line(1:10)=='BR_diagStp')THEN                        ! Flag to disable diagonal step at the end of each Brent iteration
        READ(line(11:16000),*,END=122) flagChar
        IF(flagChar == 'N')THEN
          optim%Brent_diagonalStep = .FALSE.
        ELSEIF(flagChar .NE. 'N' .AND. flagChar .NE. 'Y')THEN
          WRITE(6,*) 'ERROR: flag to enable/disable diagonal Brent step must be either Y or N'
          STOP 1
        ENDIF
      ENDIF

      !Numerical parameters for line search [Fred]
122   IF(line(1:10)=='lnS_maxItr')THEN                        ! Maximum amount of allowed line search iterations
        READ(line(11:16000),*,END=123) optim%lineSearch_maxIter
      ENDIF

123   IF(line(1:7)=='lnS_tol')THEN                            ! Tolerance for numerical line search
        READ(line(11:16000),*,END=124) optim%lineSearch_tol
      ENDIF

124   IF(line(1:9)=='scan_numx')THEN                          ! Amount of scanning points in the dimension of the 1st parameter
        READ(line(11:16000),*,END=125) optim%scan_xpoints
      ENDIF

125   IF(line(1:9)=='scan_numy')THEN                          ! Amount of scanning points in the dimension of the 2nd parameter
        READ(line(11:16000),*,END=126) optim%scan_ypoints
      ENDIF

      !Parameters for DE-MC [David, 2013-02-11]
126   IF(line(1:9)=='DEMC_ngen')THEN                          ! Number of generations in DE-MC simulation
        READ(line(11:16000),*,END=127) optim%DEMC_ngen
      ENDIF

127   IF(line(1:9)=='DEMC_npop')THEN                          ! Number of populations in DE-MC simulation
        READ(line(11:16000),*,END=128) optim%DEMC_npop
      ENDIF
128   IF(line(1:15)=='DEMC_gammascale')THEN                   ! Scaling of the default gamma factor (2.38/sqrt(2*npar)) in DE-MC
        READ(line(17:16000),*,END=129) optim%DEMC_gammascale
      ENDIF
129   IF(line(1:14)=='DEMC_crossover')THEN                     !Crossover probability in DE-MC simulation
        READ(line(16:16000),*,END=130) optim%DEMC_crossover
      ENDIF
130   IF(line(1:10)=='DEMC_sigma')THEN                         !Standard deviation of sample perturbations in DE-MC simulation
        READ(line(12:16000),*,END=131) optim%DEMC_sigma
      ENDIF
131   IF(line(1:12)=='DEMC_accprob')THEN                         !Standard deviation of sample perturbations in DE-MC simulation
        READ(line(14:16000),*,END=132) optim%DEMC_accprob
      ENDIF
132   CONTINUE
    ENDDO

    nline = nline + ninfolines 

    IF(optim%task_MC .OR. optim%task_boundps .OR. optim%task_stageMC) optim%task_runens = .TRUE.
    IF(optim%task_DEMC)                       optim%task_runens = .TRUE.
    IF(optim%task_DEMC)                       optim%nruns_best = optim%DEMC_npop+1 !Number of populations corresponds to the number of best, add 1 to store median of the populations as "the best" NO1
    IF(.NOT. optim%task_runens) optim%nruns_best = 1  
    IF(optim%task_writesim .AND. .NOT.(optim%task_MC .OR. optim%task_boundps .OR. optim%task_DEMC))THEN
      WRITE(6,*) 'Warning: It is only allowed to write all simulation results for some of the '
      WRITE(6,*) 'Warning: MonteCarlo methods. Write all simulation results is turned off'
      optim%task_writesim = .FALSE.
    ENDIF
    
    numpar = 0
    caseFlag = 1
    DO 
      READ(fileunit_temp,'(a)',END=140) line
      CALL read_parameterline(line,dimpar,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')EXIT    !end of file
      varindex = 0
      DO j = 1,max_par                !find index of parameter
        IF(varstr==modparid(j)%shortname)THEN
          varindex = j
          EXIT
        ENDIF
      ENDDO
      IF(varindex==0)THEN             !check if parameter index found
        WRITE(6,*) ' Error - unknown variable on line: ', nline
        WRITE(6,*) ' Error - check file',filename
        STOP 1
      ENDIF
      IF(caseFlag == 1) THEN      !First occurence of parameter: minimum values
        numpar = numpar + 1
        optparmin(numpar,1:nvalues) = values(1:nvalues)
        optparid(numpar)=varindex
        caseFlag = 2
      ELSEIF(caseFlag == 2) THEN  !Second occurence of parameter: maximum values
        optparmax(numpar,1:nvalues) = values(1:nvalues)
        caseFlag = 3
      ELSEIF(caseFlag == 3) THEN  !Third occurence of parameter: calibration decimal precision for parameter
        optparprecision(numpar,1:nvalues) = values(1:nvalues)
        caseFlag = 1
      ENDIF

      nline = nline + 1
    ENDDO

140 CONTINUE
    CLOSE(fileunit_temp)
    WRITE(6,*) 'File read: ', TRIM(filename)

    !Sort min and max values
    WHERE(optparmin>optparmax)
      par=optparmax
      optparmax=optparmin
      optparmin=par
    ENDWHERE

    !Set and check worldvar numoptimpar = amount of model parameters to be optimized
    CALL count_optim_par(maxoptpar,dimpar)
    IF(numoptimpar == 0) THEN
      WRITE(6,*) 'ERROR: only zero-length model parameters intervals found in optpar.txt!'
      STOP 1 
    ENDIF

    !Allocate parindex
    IF(.NOT.ALLOCATED(parindex)) ALLOCATE(parindex(numoptimpar,2))

    IF(.NOT.(optim%task_MC .OR. optim%task_boundps .OR. optim%task_Scanning .OR. optim%task_stageMC .OR. optim%task_BrentNew .OR. optim%task_stpstDesc .OR. optim%task_DFP .OR. optim%task_BFGS .OR. optim%task_DEMC)) THEN     ! [Fred 06.10.11] Added task_Scanning to checklist
      WRITE(6,*) 'ERROR: no task set in optpar'
      STOP 1
    ENDIF
    IF(ALLOCATED(values)) DEALLOCATE(values)
    IF(ALLOCATED(par)) DEALLOCATE(par)
    RETURN

100 WRITE(6,*) 'End of file during read, optpar.txt' 
    CLOSE(fileunit_temp)
    STOP 1

  END SUBROUTINE load_optpar

  !>Reads parameters to be calibrated to find dimension of river states
  !--------------------------------------------------------------------
  SUBROUTINE calculate_special_optpar_parameters(dir,velindex,dampindex,rivvel,damp) 

    USE WORLDVAR, ONLY : fileunit_temp
    USE MODVAR, ONLY : modparid,nsoil,nluse,nsub,nregions !,nlakeregions

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir    !<file directory
    INTEGER, INTENT(IN)  :: velindex   !<index of rivvel in modparid
    INTEGER, INTENT(IN)  :: dampindex  !<index of damp in modparid
    REAL, INTENT(OUT)    :: rivvel   !<lower river velocity boundary
    REAL, INTENT(OUT)    :: damp     !<lower damp boundary

    !Local variables
    INTEGER j, caseFlag
    INTEGER nskip                   !number of lines with headings
    INTEGER nline                   !line number in file
    INTEGER localdimpar             !maximal dimention of parameters
    INTEGER nvalues                 !number of values read from line
    LOGICAL fexist                  !presence of file
    REAL minvalv,maxvalv,minvald,maxvald
    REAL,ALLOCATABLE :: values(:)   !values of parameter read from file (different soil or landuse)
    CHARACTER (LEN=210) filename
    CHARACTER (LEN=10) varstr       !string with parameter name
    CHARACTER(LEN=18000) line       !line in file

    !Local parameters
    INTEGER, PARAMETER :: ninfolines = 20

    !Initiate default value
    rivvel = 9999.; damp = 9999.   !unreasonable large values
    minvald = 0.; maxvald = 0.
    minvalv = 0.; maxvalv = 0.
    
    !Check file presence
    filename=TRIM(dir)//'optpar.txt'
    INQUIRE(FILE=TRIM(filename),EXIST=fexist)
    IF(.NOT.fexist) RETURN
    
    !Initiation of variables for optimization parameters
    localdimpar = MAX(1,nsub,nsoil,nluse,MAXVAL(nregions))
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(localdimpar))           !local

    filename=TRIM(dir)//'optpar.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
    nskip = 1
    DO j = 1,nskip   !Skip heading
      READ(fileunit_temp,*)
    ENDDO
    nline = nskip + 1
    DO j = 1, ninfolines     !Skip lines with optimisation information
      READ(fileunit_temp,'(a)',END=100) line
    ENDDO
    nline = nline + ninfolines 

    caseFlag = 1
    DO 
      READ(fileunit_temp,'(a)',END=140) line
      CALL read_parameterline(line,localdimpar,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')EXIT    !end of file

      IF(varstr==modparid(velindex)%shortname)THEN
        IF(caseFlag == 1) THEN      !First occurence of parameter: minimum values
          minvalv = values(1)
          caseFlag = 2
        ELSEIF(caseFlag == 2) THEN  !Second occurence of parameter: maximum values
          maxvalv = values(1)
          caseFlag = 3
        ELSEIF(caseFlag == 3) THEN  !Third occurence of parameter: calibration decimal precision for parameter
          caseFlag = 1
        ENDIF
      ENDIF
      IF(varstr==modparid(dampindex)%shortname)THEN
        IF(caseFlag == 1) THEN      !First occurence of parameter: minimum values
          minvald = values(1)
          caseFlag = 2
        ELSEIF(caseFlag == 2) THEN  !Second occurence of parameter: maximum values
          maxvald = values(1)
          caseFlag = 3
        ELSEIF(caseFlag == 3) THEN  !Third occurence of parameter: calibration decimal precision for parameter
          caseFlag = 1
        ENDIF
      ENDIF

      nline = nline + 1
    ENDDO

140 CONTINUE
    CLOSE(fileunit_temp)

    !Set output
    IF(minvalv>0. .OR. maxvalv>0.) rivvel = MIN(minvalv,maxvalv)
    IF(minvald>0. .OR. maxvald>0.) damp = MIN(minvald,maxvald)
    IF(ALLOCATED(values)) DEALLOCATE(values)
    RETURN
    
100 WRITE(6,*) 'End of file during read optpar.txt for river dimension' 
    CLOSE(fileunit_temp)
    STOP 1

  END SUBROUTINE calculate_special_optpar_parameters

  !>Save optimal values of optimized parameters to file respar.txt
  !-----------------------------------------------------------------
  SUBROUTINE save_respar(dir,numpar,n) 

    USE WORLDVAR, ONLY : optparid, &
                         parindex, &
                         fileunit_temp
    USE MODVAR, ONLY   : modparid, & 
                         nluse, &
                         nsoil, &
                         nregions

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir    !<File directory
    INTEGER, INTENT(IN)          :: numpar !<number of parameters that has been calibrated
    INTEGER, INTENT(IN)          :: n      !<number of subbasins
    
    !Local variables 
    INTEGER i,j
    INTEGER dim                     !max number of values for parameter
    INTEGER num                     !number of values for parameter
    INTEGER ndec
    CHARACTER (LEN=210) filename
    CHARACTER (LEN=10)  numparstr   !number of parameter values
    CHARACTER (LEN=50)  form_str    !format for writing parameter
    INTEGER varindex                !array-index of parameter
    LOGICAL written                 !status of parameter
    REAL,ALLOCATABLE :: par(:)      !parameter values

    !>\b Algoritm \n
    !>Allocate and initiate local variables
    dim = MAX(1,n,nsoil,nluse,MAXVAL(nregions))
    IF(.NOT.ALLOCATED(par)) ALLOCATE(par(dim))

    !>Open file for writing parameter and write heading
    filename=TRIM(dir)//'respar.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'unknown')
    WRITE(fileunit_temp,*) 'Optimal value of parameters found during automatic calibration'    !heading

    !>For every calibrated parameter:
    ndec = 7    !7 decimals for parameter values
    DO i = 1,numpar
      varindex = optparid(parindex(i,1))
      written = .FALSE.
      DO j = 1, i-1
        IF(varindex==optparid(parindex(j,1))) written = .TRUE.       !Check if parameter already written to file
      ENDDO
      IF(.NOT. written)THEN
        !>\li Get current parameter value
        CALL get_parametervalues(varindex,n,dim,par,num)
        !>\li Write parameter to file
        form_str=''
        WRITE(numparstr,'(I10)') num   !Integer to character
        form_str = '(A10,'//TRIM(ADJUSTL(numparstr))//'(1x,F16.'//CHAR(ndec+48)//'))'
        WRITE(fileunit_temp,form_str) modparid(varindex)%shortname,par(1:num)       !write parameter values to file
      ENDIF
    ENDDO

    CLOSE(fileunit_temp)
    IF(ALLOCATED(par)) DEALLOCATE(par)

  END SUBROUTINE save_respar

  !>Collect parameter values from model variables
  !--------------------------------------------------------------
  SUBROUTINE get_parametervalues(varindex,n,dim,values,num)

    USE MODVAR, ONLY : soilpar,  &
                       landpar,  &
                       genpar,   &
                       regpar,   &
                       basinpar, &
                       monthpar, &
                       modparid, &
                       m_gpar,   &
                       m_bpar,   &
                       m_spar,   &
                       m_lpar,   &
                       m_rpar,   &
                       m_mpar, &
                       regiondivision, &
                       nluse, &
                       nsoil, &
                       nregions

    !Argument declarations
    INTEGER, INTENT(IN)  :: varindex    !<model parameter index
    INTEGER, INTENT(IN)  :: n           !<number of subbasins
    INTEGER, INTENT(IN)  :: dim         !<max number of parameter values (soil types/land uses/subbasins)
    REAL, INTENT(OUT)    :: values(dim) !<parameter values
    INTEGER, INTENT(OUT) :: num         !<number of parameter values (soil types/land uses/subbasins) that are used

    !>\b Algoritm \n
    !>Depending on variable type: Find variable dimension and parameter values
    IF(modparid(varindex)%deptype==m_gpar)THEN
      num = 1                    !genpar
      values(1:num) = genpar(modparid(varindex)%parno)
    ELSEIF(modparid(varindex)%deptype==m_bpar)THEN
      num = n                    !basinpar
      values(1:num) = basinpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_spar)THEN
      num = nsoil
      values(1:num) = soilpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_lpar)THEN
      num = nluse
      values(1:num) = landpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_rpar)THEN
      num = nregions(regiondivision(modparid(varindex)%parno))
      values(1:num) = regpar(modparid(varindex)%parno,1:num)
    ELSEIF(modparid(varindex)%deptype==m_mpar)THEN
      num = 12
      values(1:num) = monthpar(modparid(varindex)%parno,1:num)
    ENDIF

  END SUBROUTINE get_parametervalues

  !>Initiate accumulation variables for printout
  !>
  !>\b Consequences Module worldvar variables accdata_classload, 
  !> accdata_basinload, maptime, tmap, output are set.
  !------------------------------------------------------------------
  SUBROUTINE initiate_output_routines() 

    USE WORLDVAR, ONLY : accdata_classload, &  !OUT
                         accdata_basinload, &  !OUT
                         maptime, &  !OUT
                         tmap, &  !OUT
                         writeload, &
                         output, & !OUT
                         noutput
    USE MODVAR, ONLY : numsubstances                     

    !Local variables
    INTEGER io

    DO io = 1,noutput
      IF(output(io)%nvar>0)THEN
        output(io)%accdata%value = 0.
        output(io)%accdata%help = 0.
        output(io)%accdata%nok = 1 !=.TRUE.
        IF(output(io)%fileformat==2)THEN
          maptime = ''
          tmap = 1
        ENDIF
      ENDIF
    ENDDO
    
    IF(writeload)THEN
      IF(numsubstances>0)THEN
        accdata_classload = 0.
        accdata_basinload = 0.
      ENDIF
    ENDIF

  END SUBROUTINE initiate_output_routines

  !>Initiate output variables to missing value
  !>
  !>/b Consequences Module modvar variables firstoutstep and outvar is set.
  !-------------------------------------------------------------
  SUBROUTINE initiate_outvar(idt)

    USE MODVAR,   ONLY : outvar, & !OUT
                         missing_value, &
                         firstoutstep !OUT
    USE WORLDVAR, ONLY : dtskip      

    INTEGER, INTENT(IN) :: idt  !<current timestep 

    firstoutstep = .FALSE. 
    IF(idt==dtskip+1) firstoutstep = .TRUE.      

    outvar = missing_value

  END SUBROUTINE initiate_outvar

  !>Set concentrations of waters with zero volume to missing_value and
  !>calculate regional and upstream average for output.
  !>
  !>\b Consequences Module modvar variable outvar may be set.
  !------------------------------------------------------------------
  SUBROUTINE revise_outvar()

    USE WORLDVAR, ONLY : noutreg, &
                         outvarinfo
    USE MODVAR, ONLY : outvarid, &
                       missing_value, &
                       nsub, &
                       noutvar, & 
                       outvar, & !OUT
                       xoregobsi,xoregindex
    USE COMPOUT, ONLY : calculate_region_average, &
                        calculate_upstream_average

    !Local variables
    INTEGER isb,ir,iout
    REAL y

    !Calculate new outvar for regional variables
    DO iout = 1,noutvar
      IF(outvarinfo(iout)%areaagg==2)THEN
        DO ir = 1,noutreg
          IF(xoregindex(outvarinfo(iout)%idindex,ir)>0)THEN
            y=xoregobsi(xoregindex(outvarinfo(iout)%idindex,ir))
          ELSE
            CALL calculate_region_average(ir,outvarinfo(iout)%ovindex0,y)
          ENDIF
          outvar(ir,iout) = y    !noutreg assumed less than nsub
        ENDDO
      ENDIF
    ENDDO
  
    !Calculate new outvar for upstream variables
    DO iout = 1,noutvar
      IF(outvarinfo(iout)%areaagg==1)THEN
        CALL calculate_upstream_average(outvarinfo(iout)%ovindex0,outvarid(outvarinfo(iout)%idindex)%basearea,outvar(:,iout))
      ENDIF
    ENDDO

    !Set concentrations of waters with zero volume to missing_value
    DO iout = 1,noutvar
      IF(outvarinfo(iout)%ovindexwater>0)THEN
        DO isb = 1, nsub
          IF(outvar(isb,outvarinfo(iout)%ovindexwater)==0.) outvar(isb,iout) = missing_value
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE revise_outvar

  !>Write to file result from the nsim best MonteCarlo simulations
   !--------------------------------------------------------------------
  SUBROUTINE save_ensemble_simulations(dir,numpar,nperf,numcrit,&
       nsim,optcrit,performance,parameters)

    USE WORLDVAR, ONLY : writematlab,      &
         filename_best, &
         fileunit_temp

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir       !<File directory
    INTEGER, INTENT(IN) :: numpar             !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf              !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit            !<Number of optimised variables
    INTEGER, INTENT(IN) :: nsim               !<Number of simulations to be saved
    REAL, INTENT(IN) :: optcrit(nsim)         !<Optimation criterion for simulations
    REAL, INTENT(IN) :: performance(nsim,nperf,numcrit) !<Performance measures for simulations
    REAL, INTENT(IN) :: parameters(nsim,numpar) !<Parameter values for simulations
   
    !Local variables 
    INTEGER i

    !Open files
    OPEN(file=TRIM(dir)//filename_best,unit=fileunit_temp,status='unknown',form='formatted')
    !Write heading
    CALL write_ensemble_simulations_heading(fileunit_temp,numpar,nperf,numcrit,.FALSE.)
    !Write data
    DO i = 1,nsim
       CALL write_simulation_results(fileunit_temp,i,numpar,numpar,nperf,numcrit,&
            optcrit(i),performance(i,:,:),parameters(i,:),writematlab)
    ENDDO
    !End routine
    CLOSE(fileunit_temp)

  END SUBROUTINE save_ensemble_simulations

  !>Prepare a file to writes result from all MonteCarlo simulations
  !--------------------------------------------------------------------
  SUBROUTINE prepare_save_all_simulations(dir,filename,funit,numpar,nperf,numcrit)

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir       !<File directory
    CHARACTER(LEN=*), INTENT(IN) :: filename  !<Filename
    INTEGER, INTENT(IN) :: funit           !<Unit to be connected to file
    INTEGER, INTENT(IN) :: numpar          !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf           !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit         !<Number of optimised variables

    !Open file
    OPEN(file=TRIM(dir)//filename,unit=funit,status='unknown',form='formatted')
    !Write heading
    CALL write_ensemble_simulations_heading(funit,numpar,nperf,numcrit,.TRUE.)

  END SUBROUTINE prepare_save_all_simulations

  !>Write heading to file for MonteCarlo simulations
  !--------------------------------------------------------------------
  SUBROUTINE write_ensemble_simulations_heading(funit,numpar,nperf,numcrit,popflag)

    USE WORLDVAR, ONLY : performance_name,       &
                         writematlab,      &
                         optparid,      &
                         parindex, &
                         optim  
    USE MODVAR, ONLY :   modparid

    !Argument declarations
    INTEGER, INTENT(IN) :: funit    !<File unit
    INTEGER, INTENT(IN) :: numpar   !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf    !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit  !<Number of optimised variables
    LOGICAL, INTENT(IN) :: popflag  !<Flag for population columns included
    
    !Local variables
    INTEGER i,j,index
    INTEGER lt,lout
    CHARACTER (LEN=16)  t
    CHARACTER (LEN=800) outtxt     !heading
    CHARACTER (LEN=1)   sep        !separator

    !Write heading
    sep = ','
    IF(writematlab)THEN
       outtxt(1:4) = '%NO'//sep       !i
       outtxt(5:9) = 'CRIT'//sep     !optcrit
       lout = 9
    ELSE
       outtxt(1:3) = 'NO'//sep       !i
       outtxt(4:8) = 'CRIT'//sep     !optcrit
       lout = 8
    ENDIF
    DO i = 1,numcrit                !performance measures
       DO j = 1,nperf
          t  = performance_name(j)
          t  = ADJUSTL(t)
          lt = LEN_TRIM(t)
          outtxt(lout+1:lout+lt+1) = t(1:lt)//sep
          lout = lout+lt+1
       ENDDO
    ENDDO
    DO i = 1,numpar                 !parameters
       index = optparid(parindex(i,1))
       t = modparid(index)%shortname
       t  = ADJUSTL(t)
       lt = LEN_TRIM(t)
       outtxt(lout+1:lout+lt+1) = t(1:lt)//sep
       lout = lout+lt+1
    ENDDO
    !IF DE-MC simulation, also print out the generation and population indeces
    IF(popflag .AND. optim%task_DEMC)THEN
       outtxt(lout+1:lout+4+1) = 'jpop'//sep
       lout = lout+4+1
       outtxt(lout+1:lout+4+1) = 'igen'//sep
       lout = lout+4+1
       outtxt(lout+1:lout+4+1) = 'iacc'//sep
       lout = lout+4+1
    ENDIF
    WRITE(funit,'(a)') outtxt(1:lout-1)

  END SUBROUTINE write_ensemble_simulations_heading

  !>Write the result from the last simulations to the file with all
  !>simulations' results.
  !--------------------------------------------------------------------
  SUBROUTINE write_simulation_results(funit,i,mpar,numpar,nperf,&
       numcrit,optcrit,performance,parameters,mlab,jpop,igen,iacc)

    USE READWRITE_ROUTINES, ONLY : write_dataline

    INTEGER, INTENT(IN) :: funit         !<File unit
    INTEGER, INTENT(IN) :: i             !<Simulation number
    INTEGER, INTENT(IN) :: mpar          !<Dimension of parameters
    INTEGER, INTENT(IN) :: numpar        !<Number of optimised parameters
    INTEGER, INTENT(IN) :: nperf         !<Number of performance measures
    INTEGER, INTENT(IN) :: numcrit       !<Number of optimised variables
    REAL, INTENT(IN) :: optcrit          !<Value of optimation criterion for the simulation
    REAL, INTENT(IN) :: performance(nperf,numcrit)  !<performance criteria for the simulation
    REAL, INTENT(IN) :: parameters(mpar)  !<Parameter values for the simulation
    LOGICAL, INTENT(IN) :: mlab           !<MATLAB format for print out
    INTEGER, OPTIONAL, INTENT(IN) :: jpop !<Population Index in DE-MC simulation
    INTEGER, OPTIONAL, INTENT(IN) :: igen !<Generation Index in DE-MC simulation
    INTEGER, OPTIONAL, INTENT(IN) :: iacc !<Acceptance Index in DE-MC simulation

    !Local variables
    REAL, ALLOCATABLE :: x(:)
    INTEGER n          !Number of columns with data to be written
    INTEGER ndec       !Number of decimals to be written (rounding to ndec decimals)

    !Write data
    ndec=8
    n = 1 + nperf*numcrit + numpar
    IF(PRESENT(iacc)) n = n+3
    IF(.NOT.ALLOCATED(x)) ALLOCATE(x(n))
    IF(PRESENT(iacc))THEN
       x = (/optcrit,performance(:,:),parameters(1:numpar),REAL(jpop),REAL(igen),REAL(iacc)/)
    ELSE
       x = (/optcrit,performance(:,:),parameters(1:numpar)/)
    ENDIF
    CALL write_dataline(funit,n,x,ndec,0,0,',',0,mlab,id=i)
    DEALLOCATE(x)

  END SUBROUTINE write_simulation_results

  !>Save the loads for the last year to files. Subroutine is called once per year.
  !>
  !>\b Consequences Module worldvar variables accdata_classload and 
  !>accdata_basinload are zeroed.
  !--------------------------------------------------------------------
  SUBROUTINE save_loadfiles(dir,year)

    USE WORLDVAR, ONLY : accdata_classload, &     !OUT
                         accdata_basinload, &     !OUT
                         fileunit_temp
    USE MODVAR, ONLY : basin,               &
                       numsubstances,       &
                       nclass,              &
                       i_in,i_on,i_sp,i_pp, &
                       nsub,                &
                       max_classoutvar,     &
                       max_basinoutvar,     &
                       loadheadings
    USE READWRITE_ROUTINES, ONLY : write_dataline

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<File directory
    INTEGER, INTENT(IN) :: year           !<Current year

    !Local variables
    INTEGER i       ! substance number       
    INTEGER j       ! class number
    INTEGER k       ! heading number in 'Headings'
    INTEGER l       ! heading number 'Prepheadings'
    INTEGER s       ! subbasin number
    INTEGER subnr   ! subbasin id
    INTEGER ffunit   ! file unit
    INTEGER nvalues   ! number of values in xload
    REAL xload(max_classoutvar*nclass+max_basinoutvar) !input to write_tab_sep
    CHARACTER (LEN=7) str !file name
    CHARACTER (LEN=10),DIMENSION(max_classoutvar*nclass+max_basinoutvar+1) :: headings ! Array of headings for writing to file

    !Preparation of heading
    headings(:)='0'
    WRITE(headings(1),'(a7,a)') loadheadings(1),CHAR(9)                     !subbasin number
    DO l=2,max_classoutvar+1
      DO j=1,nclass
        k=(l-2)*nclass+j+1 ! position in Headings
        WRITE(headings(k),'(a6,a1,I2.2,a)') loadheadings(l),'_',j,CHAR(9)     !class dependent
      ENDDO
    ENDDO
    DO l=max_classoutvar+2,max_classoutvar+max_basinoutvar+1 
      k=max_classoutvar*nclass+(l-max_classoutvar)
      WRITE(headings(k),'(a6,a)') loadheadings(l),CHAR(9)                   !subbasin dependent
    ENDDO

    !Save data to files
    ffunit=fileunit_temp
    DO i=1,numsubstances
      IF (i==i_in.OR.i==i_on.OR.i==i_sp.OR.i==i_pp)THEN
        IF (i==i_in) WRITE(str,'(I4,3a)') year,'_IN'
        IF (i==i_on) WRITE(str,'(I4,3a)') year,'_ON'
        IF (i==i_sp) WRITE(str,'(I4,3a)') year,'_SP'
        IF (i==i_pp) WRITE(str,'(I4,3a)') year,'_PP'
        OPEN(UNIT=ffunit,FILE=TRIM(dir)//TRIM(ADJUSTL(str))//'.txt',status='unknown',form='formatted')
        WRITE(ffunit,'(999a10)') headings     !Write headings to file
        DO s=1,nsub
          !Write data to a row array
          DO k=1,max_classoutvar
            xload(((k-1)*nclass+1):(k*nclass))=accdata_classload(1:nclass,k,i,s)
          ENDDO
          xload((max_classoutvar*nclass+1):(max_classoutvar*nclass)+max_basinoutvar)=accdata_basinload(i,:,s)
          nvalues=max_classoutvar*nclass+max_basinoutvar
          subnr = basin(s)%subid
          CALL write_dataline(ffunit,nvalues,xload,3,0,0,CHAR(9),0,.FALSE.,id=subnr) 
        ENDDO ! end subarea loop (s)
        CLOSE(ffunit)
      ENDIF
    ENDDO !end substance lop (i)

    !Reset accumulations to zero at end of year 
    accdata_classload(:,:,:,:) = 0.0  
    accdata_basinload(:,:,:)   = 0.0 

  END SUBROUTINE save_loadfiles

  !>Check input data from info, GeoData, BranchData, GeoClass, par and
  !>CropData files
  !--------------------------------------------------------------------
  SUBROUTINE checkindata_part1()

    USE WORLDVAR, ONLY : bdate,ndt,   &
                         sdate,       &
                         checkdata,   &
                         modeldir,    &
                         fileunit_temp, &
                         comment_str, &
                         vegtypereset, &
                         output,noutput
    USE MODVAR, ONLY : basin,         &
                       classbasin,    &
                       classdata,     &
                       cropdata,      &
                       classmodel,    &
                       slc_olake,     &
                       load,          &
                       pathsubid,     &
                       branchsubid,   &
                       petmodel, &
                       soildepth,     &
                       nregions,      &
                       nclass,        &
                       nluse,nsoil,   &
                       conductN,conductP,conductC, &
                       nsub_basemodel, &
                       outvarid, &
                       modparid, &
                       modeloption, &
                       max_par, &
                       regiondivision, &
                       m_gpar,m_bpar, &
                       m_spar,m_lpar, &
                       m_rpar, &
                       m_ldpar,m_mpar, &
                       p_growthstart, &
                       p_erosion, &
                       p_petmodel

    !Local variables
    REAL    basinarea
    REAL    classtest
    REAL    basintest
    INTEGER subidproc,classproc
    INTEGER i,k,bdim,ksource,kbranch
    INTEGER j,io,ivar,ivar2
    INTEGER help1,help2,help3,help4,help5
    INTEGER dim
    INTEGER dimcheck                !current parameter supposed size
    INTEGER ncomment                !number of lines with comments
    INTEGER nline                   !line number in file
    INTEGER nvalues                 !number of values read from line
    INTEGER varindex                !parameter index in array
    LOGICAL found
    CHARACTER (LEN=220) filename
    CHARACTER (LEN=10) varstr       !string with variable name
    CHARACTER(LEN=18000) line
    REAL,ALLOCATABLE :: values(:)               !values of parameter read from file

    WRITE(6,*)
    WRITE(6,*) 'Check indata part one'
    WRITE(6,*) '---------------------'
    
    !Check simulation period
    WRITE(6,*) 'Checking info.txt...'
    IF(ndt<=0) THEN
      WRITE(6,*) 'Number of timesteps less than one. Check bdate and edate in info.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    IF(bdate%day<1.OR.bdate%day>31) WRITE(6,*) 'Check bdate day in info.txt'
    IF(bdate%month<1.OR.bdate%month>12) WRITE(6,*) 'Check bdate month in info.txt'
    IF(sdate%day<1.OR.sdate%day>31) WRITE(6,*) 'Check edate day in info.txt'
    IF(sdate%month<1.OR.sdate%month>12) WRITE(6,*) 'Check edate month in info.txt'
    DO io = 1, noutput
      IF(output(io)%fileformat==2 .OR. output(io)%fileformat==3)THEN
        DO ivar = 1,output(io)%nvar
          DO ivar2 = ivar+1,output(io)%nvar
            IF(outvarid(output(io)%variable(ivar)%idindex)%shortname==outvarid(output(io)%variable(ivar2)%idindex)%shortname &
                 .AND. output(io)%variable(ivar)%areaagg==output(io)%variable(ivar2)%areaagg)THEN
              WRITE(6,*) 'Error: Same variable multiple times in output. Variable: '//TRIM(outvarid(output(io)%variable(ivar)%idindex)%shortname)
              WRITE(6,*) 'Check timeoutput and/or mapoutput.'
              checkdata(1,2) = .TRUE.
            ENDIF
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    
    !Check GeoData
    WRITE(6,*) 'Checking GeoData.txt...'
    basinarea = SUM(basin(:)%area)
    IF(basinarea<=0.)THEN 
      WRITE(6,*) 'Error: No area of catchment. Check GeoData.txt.'
      checkdata(1,2) = .TRUE.
    ENDIF
    basinarea = MINVAL(basin(:)%area)
    IF(basinarea<=0.)THEN
      DO i = 1,nsub_basemodel
        IF(basin(i)%area<=0.)THEN
          WRITE(6,*) 'Subbasin',basin(i)%subid,'has area',basin(i)%area,'.'
          IF(SUM(classbasin(i,:)%part)/=0.)THEN
            WRITE(6,*) 'and its sum of slc-fraction is not zero (Error).'
            checkdata(1,2) = .TRUE.
          ENDIF
        ENDIF  
      ENDDO
    ENDIF
    subidproc = MINVAL(basin(:)%subid)
    IF(subidproc==0)THEN
      WRITE(6,*) 'Error: Some subid are zero. Check GeoData.txt.'
      checkdata(1,2)=.TRUE.
    ENDIF
    subidproc = MAXVAL(basin(:)%subid)
    IF(subidproc>=100000000) WRITE(6,*) 'Subid is too large. Maximum 9999999 may be used. Check GeoData.txt.'
    basintest = MINVAL(basin(:)%slope)
    IF(basintest<0)THEN
      WRITE(6,*) 'Error: Some slope is less than zero. Check GeoData.txt.'
      checkdata(1,2)=.TRUE.
    ENDIF
    IF(conductN.OR.conductP.OR.conductC) THEN
      IF(nregions(6)<=0)THEN
        WRITE(6,*) 'Error: Max value of lake regions is zero. Lake regions are necessary for'
        WRITE(6,*) 'NPC simulation. Check GeoData.txt.'
        checkdata(1,2)=.TRUE.
      ELSE
        subidproc = MINVAL(basin(:)%parregion(6))
        IF(subidproc<=0)THEN
          WRITE(6,*) 'Error: Some lake regions are not set. Lake regions are necessary for NPC simulation. Check GeoData.txt.'
          checkdata(1,2)=.TRUE.
        ENDIF
      ENDIF
    ENDIF  
    IF(MAXVAL(nregions)<=0)THEN
      WRITE(6,*) 'Max value of parameter regions is zero. Check GeoData.txt.'
    ENDIF
    IF(conductN.OR.conductP) THEN
      IF(nregions(3)<=0)THEN
        WRITE(6,*) 'Max value of water quality parameter regions is zero. Check GeoData.txt.'
      ENDIF
      subidproc = MINVAL(basin(:)%parregion(3))
      IF(subidproc<=0) WRITE(6,*) 'Some water quality parameter regions are not set. Check GeoData.txt.'
    ENDIF
    subidproc = MINVAL(basin(:)%parregion(1))
    IF(subidproc<=0) WRITE(6,*) 'Some parameter regions are not set. Check GeoData.txt.'
    subidproc = MINVAL(basin(:)%region)
    IF(subidproc<=0)THEN
      WRITE(6,*) 'Some regions are not set. Regions are necessary for NP simulation'
      WRITE(6,*) '(and irrigation). Check GeoData.txt.'
    ENDIF
    basinarea=SUM(basin(:)%rivlen(2))
    IF(basinarea<=0) WRITE(6,*) 'No main rivers? Check GeoData.txt (rivlen).'
    basinarea=MINVAL(basin(:)%rivlen(2))
    IF(basinarea<0) WRITE(6,*) 'Negative length of some rivers. That is not allowed. Check GeoData.txt'
    !Check slc-fractions
    DO i = 1,nsub_basemodel
      IF(SUM(classbasin(i,:)%part)<0.99 .OR. SUM(classbasin(i,:)%part)>1.01)THEN
        WRITE(6,'(A9,I7,A51)') 'Subbasin ',basin(i)%subid,' has sum of class-fractions deviating more than 1%.'
      ELSEIF(SUM(classbasin(i,:)%part)<0.9999 .OR. SUM(classbasin(i,:)%part)>1.0001)THEN
        WRITE(6,'(A9,I7,A54)') 'Subbasin ',basin(i)%subid,' has sum of class-fractions deviating more than 0.01%.'
      ENDIF
      IF(MINVAL(classbasin(i,:)%part)<0.)THEN
        WRITE(6,'(A9,I7,A36)') 'Subbasin ',basin(i)%subid,' has class-fractions less than zero.'
        WRITE(6,'(A38)') 'This is not allowed. Check GeoData.txt'
        checkdata(1,2) = .TRUE.
      ENDIF
    ENDDO
    !Check olake depth
    IF(slc_olake>0)THEN
      DO i = 1,nsub_basemodel
        IF(classbasin(i,slc_olake)%part>0. .AND. basin(i)%lakedepth(2)<=0.)THEN
          WRITE(6,'(A9,I7,A65)') 'Subbasin ',basin(i)%subid,' has olake with depth zero. This is not allowed, model may crash.'
          WRITE(6,*) 'Check GeoData.txt, LakeData.txt and/or DamData.txt'
          WRITE(6,*) 'If lakedepth set in par.txt (gldepo>0) then ignore this warning.'
        ENDIF
      ENDDO
    ENDIF
    !Check abstractions
    DO i = 1,nsub_basemodel
      IF(load(i)%abstrvol>0)THEN
        IF(load(i)%abstrind==0)THEN 
          WRITE(6,*) 'Warning: Abstraction of water without defined source for subbasin ',basin(i)%subid
        ELSEIF(load(i)%abstrind==1)THEN 
          IF(basin(i)%rivlen(2)==0.)THEN
            WRITE(6,*) 'Warning: Abstraction of water from main river of length 0 for subbasin ',basin(i)%subid
          ENDIF
        ELSEIF(load(i)%abstrind==2)THEN 
          IF(slc_olake==0)THEN
            WRITE(6,*) 'Abstraction of water from outlet lake for subbasin ',basin(i)%subid, 'but no slc-class for olake defined.'
            checkdata(1,2) = .TRUE.
          ELSEIF(classbasin(i,slc_olake)%part==0.)THEN
            WRITE(6,*) 'Abstraction of water from outlet lake for subbasin',basin(i)%subid, 'without outlet lake.'
            checkdata(1,2) = .TRUE.
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    !Check linkage and subbasin order
    WRITE(6,*) 'Checking subbasin linkage...'
    DO i = 1,nsub_basemodel
      DO k = i-1,1,-1
        IF(basin(k)%subid==pathsubid(i)%main)THEN
          WRITE(6,*) 'Linking error: main branch (',pathsubid(i)%main,') of subbasin',basin(i)%subid,'comes before subbasin in GeoData.txt.'
          checkdata(1,2) = .TRUE.
        ENDIF
      ENDDO
      found = .FALSE.
      IF(pathsubid(i)%main<=0) found = .TRUE.
      IF(.NOT.found)THEN
        DO k = i+1,nsub_basemodel
          IF(pathsubid(i)%main==basin(k)%subid)THEN
            found=.TRUE.
            EXIT
          ENDIF  
        ENDDO
      ENDIF
      IF(.NOT.found) WRITE(6,*) 'Main branch (',pathsubid(i)%main,') of subbasin',basin(i)%subid,' not found. Check GeoData.txt'
    ENDDO
    IF(ALLOCATED(branchsubid))THEN
      bdim=SIZE(branchsubid(:)%source)
      DO i = 1,bdim
        ksource = 0
        kbranch = 0
        DO k = 1,nsub_basemodel
          IF(branchsubid(i)%source==basin(k)%subid)  ksource = k
          IF(branchsubid(i)%branch==basin(k)%subid)  kbranch = k
        ENDDO
        IF(ksource==0) WRITE(6,*) 'Branch from subid',branchsubid(i)%source,' not in model. Check BranchData.txt.'
        IF(kbranch==0)THEN
          IF(branchsubid(i)%branch>0)  WRITE(6,*) 'Branch (',branchsubid(i)%branch,') of subbasin',branchsubid(i)%source,' not found. Check BranchData.txt.'
        ELSE
          IF(kbranch<ksource)THEN
            WRITE(6,*) 'Linking error: Branch (',branchsubid(i)%branch,') of subbasin'
            WRITE(6,*) branchsubid(i)%source,' lies upstream subbasin. Check BranchData.txt.'
            checkdata(1,2) = .TRUE.
          ENDIF
        ENDIF  
      ENDDO  
    ENDIF
    IF(modeloption(p_growthstart)==1)THEN
      basinarea=SUM(basin(:)%latitude)
      IF(basinarea<=0) WRITE(6,*) 'Latitude is needed for using growth start model 1. No latitude given. Check GeoData.txt (latitude)?'
    ENDIF
    IF(ALLOCATED(petmodel))THEN
      classproc = MINVAL(petmodel)
      IF(classproc<0)THEN
        WRITE(6,*) 'Error: Illegal value of petmodel in GeoData.txt. Must be 0-5.'
        checkdata(1,2) = .TRUE.
      ENDIF
      classproc = MAXVAL(petmodel)
      IF(classproc>5)THEN
        WRITE(6,*) 'Error: Illegal value of petmodel in GeoData.txt. Must be 0-5.'
        checkdata(1,2) = .TRUE.
      ENDIF
    ELSEIF(modeloption(p_petmodel)<0.OR.modeloption(p_petmodel)>5)THEN
      WRITE(6,*) 'Error: Modeloption petmodel must be 0-5. Check info.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    
    !Check GeoClass
    WRITE(6,*) 'Checking GeoClass.txt...'
    IF(nclass<=0)THEN
      WRITE(6,*) 'Error: Zero classes. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    classproc=MINVAL(classdata(:)%luse)
    IF(classproc<=0)THEN
      WRITE(6,*) 'Error: No land use code given. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    classproc=MINVAL(classdata(:)%soil)
    IF(classproc<=0)THEN
      WRITE(6,*) 'Error: No soil type code given. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    help1=0;help2=0;help3=0;help4=0;help5=0
    DO i=1,nclass
      IF(classmodel(i)==0)THEN
      ELSEIF(classmodel(i)==1)THEN
        help1 = help1 + 1
      ELSEIF(classmodel(i)==2)THEN
        help2 = help2 + 1
      ELSEIF(classmodel(i)==3)THEN
        help3 = help3 + 1
      ELSEIF(classmodel(i)==11)THEN
        help4 = help4 + 1
      ELSEIF(classmodel(i)==12)THEN
        help5 = help5 + 1
      ELSE
        WRITE(6,*) 'Unknown information in special class column for slc-class',i,'. Check GeoClass.txt.'
      ENDIF
    ENDDO
    IF(help1>1)THEN
      WRITE(6,*) 'ERROR: More than one olake class. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    IF(help2>1)THEN
      WRITE(6,*) 'ERROR: More than one ilake class. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    IF(help3>1)THEN
      WRITE(6,*) 'ERROR: More than one glacier class. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    IF(help4>1)THEN
      WRITE(6,*) 'ERROR: More than one main river class. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    IF(help5>1)THEN
      WRITE(6,*) 'ERROR: More than one local river class. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    classtest = MINVAL(soildepth(1,:))
    IF(classtest<=0.)THEN
      WRITE(6,*) 'Zero soil depths found for some class. Check GeoClass.txt'
    ENDIF
    classtest=MAXVAL(soildepth(:,:))
    IF(classtest<=0.)THEN
      WRITE(6,*) 'Error: No soil depths found. Check GeoClass.txt'
      checkdata(1,2) = .TRUE.
    ENDIF
    DO i = 1,nclass
      IF(vegtypereset)THEN
        WRITE(6,*) 'Error: vegetation type not given in GeoClass.txt for slc-class vegtype=1 (open) is used (when no indatacheck is off)'
        checkdata(1,2) = .TRUE.
      ENDIF
      IF(MAXVAL(soildepth(:,i))<classdata(i)%streamdepth)THEN
        WRITE(6,*) 'Warning: Stream depth below deepest soillayer for class:',i
      ENDIF
    ENDDO

    !Check parameters
    WRITE(6,*) 'Checking par.txt...'
    filename=TRIM(modeldir)//'par.txt'
    OPEN(UNIT = fileunit_temp,FILE = filename, STATUS = 'old', ACTION='read')
    !Initiate some variables
    dim = MAX(1,nsub_basemodel,nsoil,nluse,MAXVAL(nregions),12)
    IF(.NOT.ALLOCATED(values)) ALLOCATE(values(dim))
    !Read and check heading
    READ(fileunit_temp,'(a)',END=100) line
    READ(line,*,ERR=10) varstr
    varindex = 0
    DO j = 1,max_par
      IF(varstr==modparid(j)%shortname)THEN
        IF(modparid(j)%deptype/=m_ldpar)THEN 
          varindex = j
          EXIT
        ENDIF
      ENDIF
    ENDDO
    IF(varindex>0) WRITE(6,*) 'First line (heading) has parameter name. This parameter will not be used. Ok?'
10  CONTINUE
    !Read and check rest of parameters and their values    
    nline = 1
    ncomment = 0
    DO 
      nline = nline + 1
      READ(fileunit_temp,'(a)',END=100) line
      IF(line(1:2)==comment_str)THEN
        ncomment = ncomment + 1
        CYCLE
      ENDIF
      CALL read_parameterline(line,dim,varstr,values,nvalues)    !read one line with parameter values
      IF(varstr=='          ')THEN
        WRITE(6,*) 'No parameter name found on line',nline,'. End of par.txt file?'
        EXIT    !end of file
      ENDIF
      varindex = 0
      DO j = 1,max_par      !Find corresponding varindex
        IF(varstr==modparid(j)%shortname)THEN
          IF(modparid(j)%deptype/=m_ldpar)THEN 
            varindex = j
            EXIT
          ENDIF
        ENDIF
      ENDDO
      IF(varindex==0)THEN
        WRITE(6,*) 'Unknown variable (',TRIM(varstr),') on line',nline,'.'
      ELSE
        !Find and check variable dimension
        dimcheck = 0
        SELECT CASE(modparid(varindex)%deptype)
        CASE(m_gpar)
          dimcheck = 1  
        CASE(m_bpar)
          dimcheck = nsub_basemodel 
        CASE(m_spar)
          dimcheck = nsoil
        CASE(m_lpar)
          dimcheck = nluse
        CASE(m_rpar)
          dimcheck = nregions(regiondivision(modparid(varindex)%parno))
        CASE(m_mpar)
          dimcheck = 12
        END SELECT
        IF(dimcheck /= nvalues)THEN
          IF(modparid(varindex)%deptype==m_rpar .AND. dimcheck<nvalues)THEN
            WRITE(6,*) 'To many values in parameter file for parameter ',TRIM(varstr),'.'
            WRITE(6,*) 'This parameter is parameter region dependent, and extra values will be skipped.'
          ELSE
            WRITE(6,*) 'Error: wrong number of values in parameter file for parameter ',TRIM(varstr),'.'
            checkdata(1,2) = .TRUE.
          ENDIF
        ENDIF
      ENDIF
    ENDDO
    
100 CONTINUE
    WRITE(6,*) 'Found',ncomment,' lines with comments in par.txt.'
    CLOSE(fileunit_temp)
    IF(ALLOCATED(values)) DEALLOCATE(values)

    !Check CropData.txt
    WRITE(6,*) 'Checking CropData.txt...'
    IF(modeloption(p_growthstart)==1)THEN
      filename=TRIM(modeldir)//'CropData.txt'
      INQUIRE(FILE = filename, EXIST= found)
      IF(.NOT.found)THEN
        WRITE(6,*) 'CropData.txt necessary for using growth start model 1. File not found.'
        checkdata(1,2) = .TRUE.
      ELSE
        basinarea=SUM(cropdata(:)%gddsow)
        IF(basinarea<=0) WRITE(6,*) 'Four columns needed for using growth start model 1, but gddsow not found. Check CropData.txt.'
        basinarea=MINVAL(cropdata(:)%gddsow)
        IF(basinarea<=0) WRITE(6,*) 'Negative or zero gddsow given. Check CropData.txt'
      ENDIF
    ENDIF
    IF(modeloption(p_erosion)==0)THEN
      filename=TRIM(modeldir)//'CropData.txt'
      INQUIRE(FILE = filename, EXIST= found)
      IF(.NOT.found)THEN
        WRITE(6,*) 'CropData.txt necessary for using erosion model 0 (default). Erosion model is used for phosphorus and sediment simulations. File not found.'
        checkdata(1,2) = .TRUE.
      ENDIF
    ENDIF

    WRITE(6,*) '---------------------'

  END SUBROUTINE checkindata_part1

  !>Check forcing data and other observations
  !------------------------------------------------------------
  SUBROUTINE checkindata_part2() 

    USE WORLDVAR, ONLY : modeldir, &
                         fileunit_temp, &
                         filename_Qobs, &
                         maxcharpath, &
                         get_seq_filename, &
                         forcingdata, &
                         max_forcingdata, &
                         i_pobs,i_tobs, &
                         i_rhobs,i_sfobs, &
                         i_swobs,i_uobs, &
                         bdate,sdate, &
                         checkdata, &
                         i_str,i_real,i_intg
    USE MODVAR, ONLY : basin, &
                       tobselevation, &
                       nsub_basemodel  
  USE READWRITE_ROUTINES, ONLY : check_obs_timeperiod, &
                                 check_pt_stn, &
                                 check_q_stn

    !Argument declarations
    
    !Local parameters

    !Local variables
    CHARACTER(LEN=16) filename
    CHARACTER(LEN=maxcharpath) filepath
    TYPE(DateType) :: fbdate,fedate
    REAL    forcproc
    INTEGER i,iforc
    INTEGER forcid
    INTEGER status
    INTEGER obsindex(nsub_basemodel)         !Index of observation stations
    INTEGER nobscol                !Number of columns with observations in file
    INTEGER numneg(nsub_basemodel)
    LOGICAL fileexist,notimefound
    INTEGER,ALLOCATABLE :: xi(:,:)           !Integer data read from file
    REAL,ALLOCATABLE    :: xr(:,:)           !Real data read from file

    WRITE(6,*)
    WRITE(6,*) 'Check indata part two'
    WRITE(6,*) '---------------------'
    
    !Checking ForcKey and coupling between subbasin and observations
    WRITE(6,*) 'Checking ForcKey.txt...'
    IF(ALLOCATED(tobselevation))THEN
      forcproc = MINVAL(tobselevation)
      IF(forcproc<0.) WRITE(6,*) 'Elevation of forcing data (temperature) is below zero. Check ForcKey.txt.'
    ENDIF
    DO iforc = 1, max_forcingdata
      IF(forcingdata(iforc)%readfile)THEN
        forcid = MINVAL(forcingdata(iforc)%stationid)
        IF(forcid<=0)THEN
          WRITE(6,*) 'Identification number of '//forcingdata(iforc)%filename//' forcing is zero for some subbasin. Check ForcKey.txt'
          checkdata(2,2) = .TRUE.
        ENDIF
      ENDIF
    ENDDO
    
    !Checking Pobs, Tobs, TMINobs, TMAXobs, RHobs, SFobs, SWobs, Uobs
    DO iforc = 1,max_forcingdata
      IF(forcingdata(iforc)%readfile)THEN
        WRITE(6,*) 'Checking '//TRIM(forcingdata(iforc)%filename)//'.txt...'
        filename = forcingdata(iforc)%filename
        CALL get_seq_filename(filename)
        filepath = TRIM(modeldir)//TRIM(filename)
        !Check that file exist
        INQUIRE(FILE=TRIM(filepath),EXIST=fileexist)
        IF(.NOT.fileexist)THEN
          WRITE(6,*) 'Error: Data file is missing, while setting is set to read it.'
          WRITE(6,*) 'Error: ', TRIM(filepath)
          checkdata(2,2)=.TRUE.
        ELSE
          CALL check_pt_stn(fileunit_temp,filepath,nsub_basemodel,forcingdata(iforc)%stationid,obsindex,nobscol,status)
          IF(status/=0)THEN
            WRITE(6,*) 'Error: Missmatch in observations between model and '//TRIM(forcingdata(iforc)%filename)//'.txt. Check ForcKey and '//TRIM(forcingdata(iforc)%filename)//'.txt.'
            checkdata(2,2)=.TRUE.
          ENDIF 
          CALL check_obs_timeperiod(fileunit_temp,filepath,1,bdate,   &
               sdate,fbdate,fedate,notimefound,status)
          IF(status.NE.0)THEN
            WRITE(6,*) 'Error: Forcing data has problems with time period. Check '//TRIM(forcingdata(iforc)%filename)//'.txt.'
            checkdata(2,2)=.TRUE.
          ENDIF
          IF(iforc==i_pobs.OR.iforc==i_rhobs.OR.iforc==i_sfobs.OR.iforc==i_swobs.OR.iforc==i_uobs)THEN
            CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,0,numneg,status)
            IF(status==1)THEN
              WRITE(6,*) 'Negative forcing data exist. Check '//TRIM(forcingdata(iforc)%filename)//'.txt.'
              IF(SUM(numneg)>0)THEN 
                checkdata(2,2)=.TRUE. 
                WRITE(6,*) 'The following subbasins have negative values:'
                DO i = 1,nsub_basemodel
                  IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid,'obsid:',forcingdata(iforc)%stationid(i)
                ENDDO
              ENDIF  
            ENDIF
          ENDIF
        ENDIF
      ENDIF
    ENDDO

    !Checking Qobs
    WRITE(6,*) 'Checking Qobs.txt...'
    filepath = TRIM(modeldir)//filename_Qobs
    INQUIRE(FILE=filepath,EXIST=fileexist)
    IF(fileexist)THEN
      CALL check_q_stn(fileunit_temp,filepath,nsub_basemodel,   &
                       basin(:)%subid,nobscol,obsindex,status)
      IF(status==0)THEN
        CALL check_obs_timeperiod(fileunit_temp,filepath,1,   &
            bdate,sdate,fbdate,fedate,notimefound,status)
        IF(status==1)THEN
          WRITE(6,*) 'Error: Discharge data has problems with time period. Check Qobs.txt.'
          checkdata(2,2)=.TRUE.
        ENDIF
        CALL check_data_positive(filepath,notimefound,nobscol,nsub_basemodel,obsindex,1,numneg,status)
        IF(status==1)THEN
          WRITE(6,*) 'Negative discharge data exist. Check Qobs.txt.'
          IF(SUM(numneg)>0)THEN 
            WRITE(6,*) 'The following subbasins have negative discharge:'
            DO i = 1,nsub_basemodel
              IF(numneg(i)>0) WRITE(6,*) 'subid:',basin(i)%subid
            ENDDO
          ENDIF  
        ENDIF
      ENDIF
    ENDIF

    IF(ALLOCATED(xi)) DEALLOCATE(xi)
    IF(ALLOCATED(xr)) DEALLOCATE(xr)

    WRITE(6,*) '---------------------'
    RETURN

  END SUBROUTINE checkindata_part2

  !>Check for negative data in file maybe also for missing values (Pobs.txt, Qobs.txt)
  !------------------------------------------------------------
  SUBROUTINE check_data_positive(filepath,notimefound,ncols,ns,oindex,miss,numneg,negfound)

    USE WORLDVAR, ONLY : readmatlab,      &
                         fileunit_temp   
    USE MODVAR, ONLY : missing_value

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: filepath  !<File
    LOGICAL, INTENT(IN)  :: notimefound  !<Flag for time format: date or date and time
    INTEGER, INTENT(IN)  :: ncols        !<Number of data columns in file
    INTEGER, INTENT(IN)  :: ns           !<Number of subbasins, basemodel
    INTEGER, INTENT(IN)  :: oindex(ns)   !<Index for columns used in model set-up
    INTEGER, INTENT(IN)  :: miss         !<Flag for counting missing values, 0=count, 1=ignore
    INTEGER, INTENT(OUT) :: numneg(ns)   !<Number of negative values per subbasin
    INTEGER, INTENT(OUT) :: negfound     !<Code for found negative data
    
    !Local variables
    INTEGER timeform
    INTEGER i
    INTEGER :: neg(ncols)
    REAL :: y(ncols)                 !Data (one time step)
    CHARACTER(LEN=16)  d2,d3         !Date yyyy-mm-dd[ hh:mm]

    IF(notimefound)THEN
      timeform = 0
    ELSE
      timeform = 1
    ENDIF

    !Load precipitation forcing data
    OPEN(UNIT = fileunit_temp,FILE = filepath, STATUS = 'old', ACTION='read')
    READ(fileunit_temp,*)  !Skip heading
    neg = 0
    numneg = 0
    DO
      y = missing_value
      IF(readmatlab)THEN
        READ(fileunit_temp,*,END=900) d2,y
      ELSEIF(timeform==0)THEN
        READ(fileunit_temp,*,END=900) d2,y
      ELSEIF(timeform==1)THEN
        READ(fileunit_temp,*,END=900) d2,d3,y    
      ENDIF
      DO i = 1,ncols
        IF(miss==0)THEN
          IF(y(i)<0.) neg(i) = neg(i) + 1
        ELSEIF(miss==1)THEN
          IF(y(i)<0..AND.y(i)/=missing_value) neg(i) = neg(i) + 1
        ENDIF
      ENDDO  
    ENDDO
900 IF(SUM(neg)==0.)THEN
      negfound = 0    !no negative data in file
    ELSE
      negfound = 1
      DO i = 1,ns
        IF(oindex(i)>0)THEN
          numneg(i) = neg(oindex(i))
        ENDIF
      ENDDO
    ENDIF

  END SUBROUTINE check_data_positive

  !>Stops simulation if checkindata routine asks for it
  !------------------------------------------------------------
  SUBROUTINE checkindata_stop() 

    USE WORLDVAR, ONLY : checkdata

    IF(checkdata(1,2).OR.checkdata(2,2))THEN
      WRITE(6,*)
      WRITE(6,*) 'Simulation stopped'
      WRITE(6,*) 'Serious errors in indata found'
      STOP 1
    ENDIF

  END SUBROUTINE checkindata_stop

END MODULE
