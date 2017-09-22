!> \file statedata.f90
!> Contains module state_datamodule.

!> \brief Load and save model states.
!>
!> Procedures for loading and saving initial states from file. Also 
!> processing them for submodel.
MODULE STATE_DATAMODULE
  !Copyright 2014-2017 SMHI
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
  USE STATETYPE_MODULE
  USE LibDate
  !Subroutines also uses modvar, worldvar,modelmodule,readwrite_routines
  
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: initiate_state_for_submodel,&
            load_saved_state,&
            finalize_outstate

CONTAINS
  

  !>Initiate state variables for submodel simulation
  !----------------------------------------------------------------------
  SUBROUTINE initiate_state_for_submodel(dir,indexarray,ml,stateinput,   &
                                         frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 

    USE MODELMODULE, ONLY : initiate_model,   &
                            initiate_model_state
    USE MODVAR, ONLY : nsub,                     &
                       nclass, &
                       maxsoillayers, &
                       nsub_basemodel,           &
                       numsubstances,            &
                       naquifers,     &
                       conductN,conductP,conductC,conductS, &
                       i_t1,i_t2, &
                       timesteps_per_day, &
                       doirrigation, &
                       conductflood,  &
                       modeloption,p_lakeriverice, &
                       p_growthstart,p_infiltration, &
                       doupdate,i_qar,i_war, &
                       wetlandexist,glacierexist, &
                       nrivertypes,nlaketypes

    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir      !<file directory
    INTEGER, INTENT(IN)  :: indexarray(nsub) !<index for basemodel
    INTEGER,INTENT(IN)   :: ml               !<max lag in steps
    LOGICAL, INTENT(IN)  :: stateinput       !<code for reading state
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states

    !Local variables
    TYPE(snowicestatetype) :: frozenstate2   !Temporary snow and ice states
    TYPE(soilstatetype)    :: soilstate2   !Temporary soil states
    TYPE(aquiferstatetype) :: aquiferstate2  !<Aquifer states
    TYPE(riverstatetype)   :: riverstate2  !Temporary river states
    TYPE(lakestatetype)    :: lakestate2   !Temporary lake states
    TYPE(miscstatetype)    :: miscstate2   !Temporary misc states

    !>/b Algoritm /n
    CALL deallocate_model_states(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    !>If statefiles exist: read and store states temporary
    IF(stateinput)THEN
      CALL allocate_model_states(nsub_basemodel,numsubstances,nclass,naquifers,maxsoillayers,nrivertypes,nlaketypes,ml,timesteps_per_day,   &
                               conductN,conductP,conductC,conductS,conductflood,i_t1>0,i_t2>0,wetlandexist,glacierexist, &
                               modeloption(p_lakeriverice)>=1,doirrigation,doupdate(i_qar).OR.doupdate(i_war),modeloption(p_growthstart)==1, &
                               modeloption(p_infiltration)==1,frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2) 
      CALL load_saved_state(dir,nsub_basemodel,ml,frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2)
    ENDIF
    !>Reallocate state variables to submodel size
    CALL allocate_model_states(nsub,numsubstances,nclass,naquifers,maxsoillayers,nrivertypes,nlaketypes,ml,timesteps_per_day,  &
                               conductN,conductP,conductC,conductS,conductflood,i_t1>0,i_t2>0,wetlandexist,glacierexist, &
                               modeloption(p_lakeriverice)>=1,doirrigation,doupdate(i_qar).OR.doupdate(i_war),modeloption(p_growthstart)==1, &
                               modeloption(p_infiltration)==1,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    !>If statefiles exist: Initiate state variables from those and deallocate temporary storage
    IF(stateinput)THEN
      CALL initiate_frozenstate_submodel(numsubstances,nsub,indexarray,glacierexist,modeloption(p_lakeriverice)>=1,frozenstate,frozenstate2)
      CALL initiate_soilstate_submodel(numsubstances,nsub,indexarray,soilstate,soilstate2)
      CALL initiate_aquiferstate_submodel(numsubstances,naquifers,aquiferstate,aquiferstate2)
      CALL initiate_riverstate_submodel(numsubstances,nsub,conductN,conductP,conductC,conductS,i_t1>0,wetlandexist,indexarray,riverstate,riverstate2)
      CALL initiate_lakestate_submodel(numsubstances,nsub,conductN.OR.conductP.OR.conductC.OR.conductS,i_t2>0,indexarray,lakestate,lakestate2)
      CALL initiate_miscstate_submodel(numsubstances,nsub,i_t1>0,doirrigation, &
                doupdate(i_qar).OR.doupdate(i_war),conductflood,modeloption(p_growthstart)==1, &
                wetlandexist,conductC,indexarray,miscstate,miscstate2)
      CALL deallocate_model_states(frozenstate2,soilstate2,aquiferstate2,riverstate2,lakestate2,miscstate2)
    !>Else: Initiate state variables with default values
    ELSE
      CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    ENDIF
    !>Initiate other model variables and parameters
    CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)  

  END SUBROUTINE initiate_state_for_submodel

  !>Load starting state from file and initiate state variables
  !------------------------------------------------------------------
  SUBROUTINE load_saved_state(dir,ns,ml,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE MODVAR, ONLY : numsubstances, &
                       maxsoillayers, &
                       nclass, &   
                       naquifers, &
                       seconds_per_timestep, &
                       timesteps_per_day, &
                       conductN,conductP,conductC,conductS, &
                       i_t1, &
                       wetlandexist,doirrigation, &
                       glacierexist,  &
                       doupdate,i_qar,i_war,  &
                       conductflood,  &
                       modeloption,p_lakeriverice, &
                       p_growthstart,p_infiltration, &
                       nrivertypes,nlaketypes
    USE WORLDVAR, ONLY : fileunit_temp, &   
                         bdate
    USE READWRITE_ROUTINES, ONLY : read_array_from_file
 
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: dir   !<file directory
    INTEGER, INTENT(IN)  :: ns            !<number of subbasins
    INTEGER, INTENT(IN)  :: ml            !<max lag in steps
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Local variables
    INTEGER ffunit               
    INTEGER ios
    INTEGER dim,idim
    INTEGER ipiece,npiece
    INTEGER checkstatus,readsubst
    LOGICAL readN,readP,readC,readS,readar,readT1,readT2
    INTEGER, ALLOCATABLE :: sectionlimits(:,:)
    REAL, ALLOCATABLE :: array(:)
    CHARACTER(LEN=28) filename   
    CHARACTER(LEN=16) bdatestr  
     
    !Local parameters
    INTEGER, PARAMETER :: seconds_per_day  = 86400 

    !Beginning of subroutine
    IF(seconds_per_timestep==seconds_per_day)THEN
      CALL format_date(bdate,'yyyymmdd',bdatestr)
    ELSE
      CALL format_date(bdate,'yyyymmddHHMM',bdatestr)
    ENDIF
    filename = 'state_save'//TRIM(ADJUSTL(bdatestr))//'.txt'                
    ffunit = fileunit_temp
    OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=ffunit,STATUS='old',FORM='formatted',IOSTAT=ios,ACTION='read')
    IF(ios/=0) THEN
      WRITE(6,*) 'ERROR: Statefile ', filename, ' not found'
      STOP 1
    ENDIF

    CALL read_and_perform_state_check(ffunit,checkstatus,readsubst,ml,readN,readP,readC,readS,readar,readT1,readT2)

    CALL get_frozenstate_variables_arraysize(ns,readsubst,nclass,nrivertypes,nlaketypes,glacierexist,modeloption(p_lakeriverice)>=1,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_frozenstate_variables_from_array(ns,numsubstances,readsubst, &
                  nclass,nrivertypes,nlaketypes,glacierexist,modeloption(p_lakeriverice)>=1,frozenstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_soilstate_variables_arraysize(ns,readsubst,nclass,maxsoillayers,readN,readP,readC,readS,readT1,modeloption(p_infiltration)==1,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_soilstate_variables_from_array(ns,numsubstances,readsubst,nclass,  &
                  maxsoillayers,conductN,conductP,conductC,conductS,i_t1>0,modeloption(p_infiltration)==1,soilstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_aquiferstate_variables_arraysize(naquifers,readsubst,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_aquiferstate_variables_from_array(naquifers,numsubstances,readsubst,aquiferstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_riverstate_variables_arraysize(ns,readsubst,nrivertypes,ml,timesteps_per_day,  &
              readN,readP,readC,readS,readT1,wetlandexist,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_riverstate_variables_from_array(ns,numsubstances,readsubst,  &
                  nrivertypes,ml,timesteps_per_day,conductN,conductP,conductC,conductS,i_t1>0,wetlandexist,readN,  &
                  readP,readC,readS,readT1,riverstate,sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_lakestate_variables_arraysize(ns,readsubst,nlaketypes,readN.OR.readP.OR.readC.OR.readS,readT2,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_lakestate_variables_from_array(ns,numsubstances,readsubst,nlaketypes, &
                  readN.OR.readP.OR.readC.OR.readS,readT2,lakestate,sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_miscstate_variables_arraysize(ns,numsubstances,nclass,readT1, &
                     doirrigation,readar,conductflood,modeloption(p_growthstart)==1,wetlandexist,readC,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL read_array_from_file(ffunit,100,idim,array)
        CALL set_miscstate_variables_from_array(ns,numsubstances,readsubst, &
                  nclass,readT1,i_t1>0,doirrigation,doupdate(i_qar).OR.doupdate(i_war),   &
                  readar,conductflood,modeloption(p_growthstart)==1,wetlandexist,readC,conductC,miscstate, &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CLOSE(ffunit)
    WRITE(6,*) 'File read: ', TRIM(dir)//TRIM(filename)


  END SUBROUTINE load_saved_state

  !>Saves state values for later use as starting state
  !---------------------------------------------------
  SUBROUTINE finalize_outstate(dir,ml,stateoutdate,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate) 

    USE MODVAR, ONLY : numsubstances, &
                       maxsoillayers, &
                       nsub, nclass, &   
                       naquifers, &
                       i_t1,i_t2, &
                       seconds_per_timestep,    &
                       timesteps_per_day, &
                       conductN,conductP,conductC,conductS,  &
                       wetlandexist,doirrigation, &
                       glacierexist,  &
                       conductflood,  &
                       doupdate,i_qar,i_war, &
                       modeloption,p_lakeriverice, &
                       p_growthstart,p_infiltration, &
                       nrivertypes,nlaketypes
    USE WORLDVAR, ONLY : fileunit_temp
    USE READWRITE_ROUTINES, ONLY : write_array_to_file

    !Argument declaration
    CHARACTER(LEN=*), INTENT(IN) :: dir            !<file directory
    INTEGER,INTENT(IN)           :: ml             !<max lag in steps
    TYPE(DateType), INTENT(IN) :: stateoutdate     !<date for writing state           
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate   !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    
    !Local variables
    INTEGER ipiece,npiece
    INTEGER ffunit
    INTEGER dim,idim
    INTEGER,ALLOCATABLE :: sectionlimits(:,:)
    REAL,ALLOCATABLE :: array(:)
    CHARACTER(LEN=16) stateoutdatestr
    CHARACTER(LEN=28) filename  
    
    !Local parameters
    INTEGER, PARAMETER :: seconds_per_day  = 86400 

    !Save STATETYPE_MODULE state variables 
    ffunit = fileunit_temp
    IF(seconds_per_timestep==seconds_per_day)THEN
      CALL format_date(stateoutdate,'yyyymmdd',stateoutdatestr)
    ELSE
      CALL format_date(stateoutdate,'yyyymmddHHMM',stateoutdatestr)
    ENDIF
    filename = 'state_save'//TRIM(ADJUSTL(stateoutdatestr))//'.txt'
    OPEN(FILE=TRIM(dir)//TRIM(filename),UNIT=ffunit,STATUS='unknown',FORM='formatted')
    CALL write_state_check(ffunit,ml)

    CALL get_frozenstate_variables_arraysize(nsub,numsubstances,nclass,nrivertypes,nlaketypes, &
              glacierexist,modeloption(p_lakeriverice)>=1,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_frozenstate_variables_to_array(nsub,numsubstances,nclass,nrivertypes,nlaketypes,  &
                  glacierexist,modeloption(p_lakeriverice)>=1,frozenstate,sectionlimits(1,ipiece),  &
                  sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_soilstate_variables_arraysize(nsub,numsubstances,nclass, &
              maxsoillayers,conductN,conductP,conductC,conductS,i_t1>0,modeloption(p_infiltration)==1,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_soilstate_variables_to_array(nsub,numsubstances,nclass,maxsoillayers, &
                  conductN,conductP,conductC,conductS,i_t1>0,modeloption(p_infiltration)==1,soilstate,sectionlimits(1,ipiece), &
                  sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_aquiferstate_variables_arraysize(naquifers,numsubstances,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_aquiferstate_variables_to_array(naquifers,numsubstances,aquiferstate,  &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_riverstate_variables_arraysize(nsub,numsubstances,nrivertypes,ml,  &
              timesteps_per_day,conductN,conductP,conductC,conductS,i_t1>0,wetlandexist,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_riverstate_variables_to_array(nsub,numsubstances,nrivertypes,ml, &
                  timesteps_per_day,conductN,conductP,conductC,conductS,i_t1>0,wetlandexist,riverstate,  &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_lakestate_variables_arraysize(nsub,numsubstances,nlaketypes,conductN.OR.conductP.OR.conductC.OR.conductS,i_t2>0,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_lakestate_variables_to_array(nsub,numsubstances,nlaketypes, &
                  conductN.OR.conductP.OR.conductC.OR.conductS,i_t2>0,lakestate,  &
                  sectionlimits(1,ipiece),sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CALL get_miscstate_variables_arraysize(nsub,numsubstances,nclass, &
                     i_t1>0,doirrigation,doupdate(i_qar).OR.doupdate(i_war), &
                     conductflood,modeloption(p_growthstart)==1,wetlandexist,conductC,dim)
    IF(dim>0)THEN
      CALL divide_large_array(dim,npiece,sectionlimits)
      DO ipiece = 1,npiece
        idim = sectionlimits(2,ipiece)-sectionlimits(1,ipiece)+1
        ALLOCATE(array(idim))
        CALL set_miscstate_variables_to_array(nsub,numsubstances,nclass,i_t1>0,  &
                  doirrigation,doupdate(i_qar).OR.doupdate(i_war),conductflood, &
                  modeloption(p_growthstart)==1,wetlandexist,conductC,miscstate,sectionlimits(1,ipiece), &
                  sectionlimits(2,ipiece),idim,array)
        CALL write_array_to_file(ffunit,100,idim,array)
        DEALLOCATE(array)
      ENDDO
      DEALLOCATE(sectionlimits)
    ENDIF

    CLOSE(ffunit)

  END SUBROUTINE finalize_outstate

  !--------------------------------------------------------------------
  !> Calculates appropriate size sections of large array
  !--------------------------------------------------------------------
  SUBROUTINE divide_large_array(dim,npiece,array)

  !Argument declarations
  INTEGER,INTENT(IN) :: dim
  INTEGER,INTENT(OUT) :: npiece
  INTEGER,ALLOCATABLE,INTENT(OUT) :: array(:,:)
  
  !Local varaibles
  INTEGER i
  INTEGER,PARAMETER :: maxchunksize = 12500000  !corresponds to real array of 50MB
  
  npiece = dim/maxchunksize
  IF(dim-maxchunksize*npiece>0) npiece = npiece + 1
  ALLOCATE(array(2,npiece))
  DO i = 1,npiece
    array(1,i) = 1 + (i-1)*maxchunksize
    array(2,i) = i*maxchunksize
  ENDDO
  array(2,npiece) = dim
  
  END SUBROUTINE divide_large_array
  
  !--------------------------------------------------------------------
  !> Saves values for later use as check if starting state is appropriate
  !--------------------------------------------------------------------
  SUBROUTINE write_state_check(ffunitloc,ml) 

    USE MODVAR, ONLY : i_in,i_sp,i_t1,i_t2,i_oc,i_ss,i_ae, &
                       numsubstances, &
                       nsub, &
                       nclass, &
                       maxsoillayers, &
                       timesteps_per_day, &
                       conductN,conductP,conductC,conductS, &
                       wetlandexist, &
                       doirrigation,  &
                       glacierexist,  &
                       doupdate,i_qar,i_war, &
                       modeloption,p_lakeriverice, &
                       p_growthstart,p_infiltration
    USE CONVERT, ONLY : logical_convert_to_integer
         
    !Argument declarations
    INTEGER, INTENT(IN) :: ffunitloc   !<File unit
    INTEGER, INTENT(IN) :: ml          !<dimension river translation variable
    
    !Local variables
    INTEGER log2intvar1,log2intvar2,log2intvar3,log2intvar4,log2intvar5,log2intvar6,log2intvar8,log2intvar9,log2intvar10,log2intvar11

    !Transform logical variables to integer
    log2intvar1 = logical_convert_to_integer(conductN)
    log2intvar2 = logical_convert_to_integer(conductP)
    log2intvar3 = logical_convert_to_integer(conductC)
    log2intvar11 = logical_convert_to_integer(conductS)
    log2intvar4 = logical_convert_to_integer(wetlandexist)
    log2intvar5 = logical_convert_to_integer(doirrigation)
    log2intvar6 = logical_convert_to_integer(glacierexist)
    log2intvar8 = logical_convert_to_integer(doupdate(i_qar).OR.doupdate(i_war))
    log2intvar9 = logical_convert_to_integer(modeloption(p_growthstart)==1)
    log2intvar10 = logical_convert_to_integer(modeloption(p_infiltration)==1)
    
    !Checkwrite for substances, number of subbasins and slc-classes
    WRITE(ffunitloc,'(8I3,I7,15I5)') numsubstances, i_in,i_sp,i_t1, &
            i_t2,i_oc,i_ss,i_ae,nsub,nclass,maxsoillayers,ml,timesteps_per_day,log2intvar1, &
            log2intvar2,log2intvar3,log2intvar11,log2intvar4,log2intvar5,log2intvar6, &
            modeloption(p_lakeriverice),log2intvar8,log2intvar9,log2intvar10

  END SUBROUTINE write_state_check

  !> Check if starting state is appropriate
  !------------------------------------------------------------
  SUBROUTINE read_and_perform_state_check(ffunitloc,status,nsubst,ml,isN,isP,isC,isS,isAR,ist1,ist2)

    USE MODVAR, ONLY : i_in,i_sp,i_t1,i_t2,i_oc,i_ss,i_ae, &
                       numsubstances, &
                       nsub_basemodel, &
                       nclass,maxsoillayers,  &
                       timesteps_per_day, &
                       conductN,conductP,conductC,conductS, &
                       wetlandexist,   &
                       doirrigation,   &
                       glacierexist,  &
                       doupdate,i_qar,i_war, &
                       modeloption,p_lakeriverice, &
                       p_growthstart,p_infiltration
    USE CONVERT, ONLY : logical_convert_to_integer, &
                        integer_convert_to_logical

    !Argument declarations
    INTEGER, INTENT(IN)  :: ffunitloc !<File unit
    INTEGER, INTENT(OUT) :: status    !<status of check
    INTEGER, INTENT(OUT) :: nsubst    !<number of substances in file
    INTEGER, INTENT(IN)  :: ml        !<dimension river translation variable
    LOGICAL, INTENT(OUT) :: isN       !<N simulated for file
    LOGICAL, INTENT(OUT) :: isP       !<P simulated for file
    LOGICAL, INTENT(OUT) :: isC       !<OC simulated for file
    LOGICAL, INTENT(OUT) :: isS       !<SS simulated for file
    LOGICAL, INTENT(OUT) :: isAR      !<AR-updating simulated for file
    LOGICAL, INTENT(OUT) :: ist1      !<T1 simulated for file
    LOGICAL, INTENT(OUT) :: ist2      !<T2 simulated for file
    
    !Local variables
    INTEGER :: statecheck_file(24)  !Checkrow from saved state file
    INTEGER :: statecheck_sim(24)   !Checkrow from current model simulation

    !>\b Algorithm \n
    !>Get model set-up for statefile and current model
    READ(ffunitloc, *) statecheck_file
    statecheck_sim(1) = numsubstances
    statecheck_sim(2) = i_in
    statecheck_sim(3) = i_sp
    statecheck_sim(4) = i_t1
    statecheck_sim(5) = i_t2
    statecheck_sim(6) = i_oc
    statecheck_sim(7) = i_ss
    statecheck_sim(8) = i_ae
    statecheck_sim(9) = nsub_basemodel
    statecheck_sim(10) = nclass
    statecheck_sim(11) = maxsoillayers
    statecheck_sim(12) = ml
    statecheck_sim(13) = timesteps_per_day
    statecheck_sim(14) = logical_convert_to_integer(conductN)
    statecheck_sim(15) = logical_convert_to_integer(conductP)
    statecheck_sim(16) = logical_convert_to_integer(conductC)
    statecheck_sim(17) = logical_convert_to_integer(conductS)
    statecheck_sim(18) = logical_convert_to_integer(wetlandexist)
    statecheck_sim(19) = logical_convert_to_integer(doirrigation)
    statecheck_sim(20) = logical_convert_to_integer(glacierexist)
    statecheck_sim(21) = modeloption(p_lakeriverice)
    statecheck_sim(22) = logical_convert_to_integer(doupdate(i_qar).OR.doupdate(i_war))
    statecheck_sim(23) = logical_convert_to_integer(modeloption(p_growthstart)==1)
    statecheck_sim(24) = logical_convert_to_integer(modeloption(p_infiltration)==1)
    
    !>Set output variables
    nsubst = statecheck_file(1)
    ist1 = statecheck_file(4)>0
    ist2 = statecheck_file(5)>0
    isN = integer_convert_to_logical(statecheck_file(14))
    isP = integer_convert_to_logical(statecheck_file(15))
    isC = integer_convert_to_logical(statecheck_file(16))
    isS = integer_convert_to_logical(statecheck_file(17))
    isAR = integer_convert_to_logical(statecheck_file(22))
    
    !>Compare model set-up with statefile set-up
    status = 2    !missmatch
    IF(ALL(statecheck_sim==statecheck_file,1)) THEN
      !identical model simulation set-up
      status = 0
    ELSEIF(statecheck_sim(22)/=statecheck_file(22))THEN
      !AR-updating in one of the model set-ups only
      statecheck_file(22)=statecheck_sim(22)
      IF(ALL(statecheck_sim==statecheck_file,1))   &
      !similar model simulation set-up
      status = 1
    ELSEIF(statecheck_sim(1)==0)THEN
      !no substances modelled
      IF(statecheck_file(9)==statecheck_sim(9).AND.  &
         statecheck_file(10)==statecheck_sim(10).AND. &
         statecheck_file(11)==statecheck_sim(11).AND. &
         statecheck_file(12)==statecheck_sim(12).AND. &
         statecheck_file(13)==statecheck_sim(13).AND. &
         statecheck_file(18)==statecheck_sim(18).AND. &
         statecheck_file(19)==statecheck_sim(19).AND. &
         statecheck_file(20)==statecheck_sim(20).AND. &
         statecheck_file(21)==statecheck_sim(21).AND. &
         statecheck_file(24)==statecheck_sim(24))  &
      !similar model simulation set-up, but without substance simulation
      status = 1
    ENDIF
    IF(status==2)THEN
      WRITE(6,*) ' '
      WRITE(6,*) 'ERROR: State file not compatible with model set-up and simulation'
      STOP 1
    ENDIF

  END SUBROUTINE read_and_perform_state_check


END MODULE
