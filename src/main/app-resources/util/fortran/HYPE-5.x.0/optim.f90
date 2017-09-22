!> \file optim.f90
!> Contains module optimization.

!> \brief Global optimization procedures for calibration of HYSS.

!  !!@todo Method documentation

MODULE OPTIMIZATION
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

  IMPLICIT NONE
  PRIVATE
!-----------------------------------------------------------------------------------------
!Private procedures in this file
!-----------------------------------------------------------------------------------------
! find_optpar
! run_model_crit
! run_model_perf
! run_model_simout
! DEMC_runmedianparameters
! DEMC_draw_R1R2
! DEMC_proposalgeneration
! DEMC_crossover
! DEMC_controlparameters
! DEMC_acceptreject_proposal
! write_DEMC_calibration_log
! write_DEMC_calibration_log_tail
! get_randompar
! reduce_parameter_space
! bookkeep_result_from_simulation
! stage_MonteCarlo_core_function
! get_randompar_by_radius
! write_stageMC_calibration_log
! write_stageMC_calibration_log_tail
! linesearch_methods_interruptor
! linesearch_methods_interruptor_printandstop
! linesearch_methods_interruptor_printandstop_core_function
! linesearch_HYSS
! linesearch_check_decim
! function_to_minim
! new_Brent_method
! initialize_Brent_calibrationlog
! write_Brent_calibrationlog
! quasiNewton_algorithm
! load_qNstarting_vector
! grad_criterium_function
! direction_multiplier
! get_epsilon
! QN_inv_hessian_update
! DFP_inv_hessian_update
! BFGS_inv_hessian_update
! print_calib_HYSSlog
! print_calib_HYSSlog_noimp
! initialize_QN_calibration_log
! write_QN_calibration_log
!-----------------------------------------------------------------------------------------
  PUBLIC :: set_optim_modpar, &
            DEMC_simulation, &
            MonteCarlo_simulation, &
            bounded_MonteCarlo_simulation, &
            param_scanning, &
            stage_MonteCarlo, &
            linesearch_methods_calibration
  CONTAINS

!>\brief Find the parameters to be optimized and set the brent-routine parameter variables
!>
!>\b Consequences Module worldvar variable parindex is set
!-----------------------------------------------------------------------------------------
SUBROUTINE find_optpar(par,parmin,parmax,parprecision,npar) 

  USE WORLDVAR, ONLY : parindex, &  !OUT
                       maxoptpar, optparmin, optparmax, optparprecision, dimpar, numoptimpar
  
  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(OUT)    :: par(numoptimpar)  !<Current parameter values for calibrated parameter
  REAL, INTENT(OUT)    :: parmin(numoptimpar)  !<Lower parameter limit for calibration
  REAL, INTENT(OUT)    :: parmax(numoptimpar)  !<Upper parameter limit for calibration
  REAL, INTENT(OUT)    :: parprecision(numoptimpar) !<Decimal precision up to which parameter is to be calibrated
  INTEGER, INTENT(OUT) :: npar              !<Actual number of parameters to be looped, = numoptimpar

  !Local variables
  INTEGER i,j
  REAL     :: parx(maxoptpar,dimpar)

  !Find the starting value
  WHERE(optparmin>0.AND.optparmax>0)    !Find mean value to start with
     parx=SQRT(optparmin*optparmax)
  ELSEWHERE
     parx=(optparmin+optparmax)/2.
  ENDWHERE

  !Find the parameters to be calibrated and assign information to used arrays
  npar = 0
  DO i = 1, maxoptpar
    DO j = 1, dimpar
      IF(optparmin(i,j)==parx(i,j).AND.optparmax(i,j)==parx(i,j))THEN 
        !Not calibrated
      ELSE
        npar = npar + 1
        parindex(npar,1) = i
        parindex(npar,2) = j
        par(npar) = parx(i,j)
        parmin(npar) = optparmin(i,j)    !brent variable set to read variable
        parmax(npar) = optparmax(i,j)
        parprecision(npar) = optparprecision(i,1)
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE find_optpar

!>\brief Set the parameter variables to the values to be used for simulation
!>
!>\b Consequences Module modvar variables soilpar, landpar, genpar, basinpar,
!>regpar and monthpar may change.
!-----------------------------------------------------------------------------------------
SUBROUTINE set_optim_modpar(npar,dpar,par) 

  USE WORLDVAR, ONLY : optparid,      &
                       parindex
  USE MODVAR, ONLY : soilpar,m_spar,  &
                     landpar,m_lpar,  &
                     genpar,m_gpar,   &
                     basinpar,m_bpar, &
                     modparid,        &
                     regiondivision, &
                     regpar,m_rpar,   &
                     monthpar,m_mpar, &
                     nregions

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: npar       !<Actual number of parameters
  INTEGER, INTENT(IN) :: dpar       !<dimension of par
  REAL, INTENT(IN)    :: par(dpar)  !<Array with parameter values to be used for this run
  
  !Local variables
  INTEGER i

  !Set current parameter value to correct parameter
  DO i = 1,npar
    IF(modparid(optparid(parindex(i,1)))%deptype==m_gpar)THEN
      genpar(modparid(optparid(parindex(i,1)))%parno) = par(i)
    ELSEIF(modparid(optparid(parindex(i,1)))%deptype==m_bpar)THEN
      basinpar(modparid(optparid(parindex(i,1)))%parno,parindex(i,2)) = par(i)
    ELSEIF(modparid(optparid(parindex(i,1)))%deptype==m_spar)THEN
      soilpar(modparid(optparid(parindex(i,1)))%parno,parindex(i,2)) = par(i)
    ELSEIF(modparid(optparid(parindex(i,1)))%deptype==m_lpar)THEN
      landpar(modparid(optparid(parindex(i,1)))%parno,parindex(i,2)) = par(i)
    ELSEIF(modparid(optparid(parindex(i,1)))%deptype==m_rpar)THEN
      IF(parindex(i,2)<=nregions(regiondivision(modparid(optparid(parindex(i,1)))%parno))) regpar(modparid(optparid(parindex(i,1)))%parno,parindex(i,2)) = par(i)
    ELSEIF(modparid(optparid(parindex(i,1)))%deptype==m_mpar)THEN
      monthpar(modparid(optparid(parindex(i,1)))%parno,parindex(i,2)) = par(i)
    ELSE
      WRITE(6,*) 'ERROR set_optim_modpar'
    ENDIF
  ENDDO

END SUBROUTINE set_optim_modpar

!>Run a simulation of the model with parameters par. 
!>Calculate criterion.
!------------------------------------------------------
SUBROUTINE run_model_crit(npar,mpar,par,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,criterion,status)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : ndt,              &
                       nacrit,            &
                       psdates,    &
                       modeldir
  USE MODELMODULE, ONLY : model, &
                          set_modelconfig_from_parameters, &
                          initiate_model,initiate_model_state
  USE MODVAR, ONLY : nsub
  USE COMPOUT, ONLY : prepare_to_compute_crit,     &
       calculate_criteria
  USE TIMEROUTINES, ONLY : calculate_time_for_model
  USE LIBDATE, ONLY : DateType
  USE DATAMODULE

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: npar      !<Number of parameters
  INTEGER, INTENT(IN) :: mpar      !<Maximum number of parameters
  REAL, INTENT(IN)    :: par(mpar) !<Parameters
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  REAL, INTENT(OUT)   :: criterion !<Value of optimization criterion
  INTEGER, INTENT(OUT):: status  !<Subroutine error status
  
  !Local variables
  REAL optcrit
  TYPE(DateType) d       
  INTEGER idt

  !Initialisations
  status = 0
  optcrit=0.
  CALL set_optim_modpar(npar,mpar,par)
  CALL set_modelconfig_from_parameters()
  CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)  

  !Run model simulation
  time_loop: DO idt = 1,ndt
     
     !Get current input data
     CALL get_current_forcing(idt,nsub,d)
     CALL calculate_time_for_model(idt,d)
     IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,d,status) 
     IF(status.NE.0) RETURN

     CALL initiate_outvar(idt)
     CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
     CALL revise_outvar()
     CALL prepare_to_compute_crit(d,idt,ndt)

  ENDDO time_loop

  !Calculate optimisation criterion     
  IF(nacrit/=0) THEN
     CALL calculate_criteria(optcrit)
  ENDIF

  !Prepare observation for yet another simulation        
  CALL reset_observations(modeldir,status)
  IF(status.NE.0) RETURN

  !Set output variables
  criterion = optcrit

END SUBROUTINE run_model_crit

!>Run a simulation of the model with parameters par. 
!>Calculate criterion and performance measures
!-------------------------------------------------------
SUBROUTINE run_model_perf(npar,mpar,par,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate, &
                          criterion,performance,status,condcrit,condthres)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : ndt,       &
                       nacrit,    &
                       maxperf,   &
                       maxsubass, &
                       psdates,    &
                       modeldir
  USE MODELMODULE, ONLY : model, &
                          set_modelconfig_from_parameters, &
                          initiate_model,initiate_model_state
  USE MODVAR, ONLY : nsub
  USE COMPOUT, ONLY : prepare_to_compute_crit, &
                      calculate_criteria
  USE TIMEROUTINES, ONLY : calculate_time_for_model
  USE LIBDATE, ONLY : DateType
  USE DATAMODULE

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: npar                     !<Number of parameters
  INTEGER, INTENT(IN) :: mpar                     !<Dimension of parameters-array
  REAL, INTENT(IN)    :: par(mpar)                !<Parameter values to be used for this simulation
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  REAL, INTENT(OUT)   :: criterion                !<Value of optimization criterion
  REAL, INTENT(OUT)   :: performance(maxperf,nacrit)    !<Simulation performance criteria
  INTEGER, INTENT(OUT):: status                 !<Subroutine error status
  REAL,OPTIONAL,INTENT(OUT) :: condcrit  !<conditional criterion value
  REAL,OPTIONAL,INTENT(OUT) :: condthres !<threshold of conditional criterion
  
  !Local variables
  REAL optcrit
  TYPE(DateType) d       
  INTEGER idt

  !Initialisations
  status = 0
  optcrit=0.
  CALL set_optim_modpar(npar,mpar,par)
  CALL set_modelconfig_from_parameters()
  CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)  
  
  !Run model simulation
  time_loop: DO idt = 1,ndt

     !Get current input data
     CALL get_current_forcing(idt,nsub,d)
     CALL calculate_time_for_model(idt,d)
     IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,d,status) 
     IF(status.NE.0) RETURN

     CALL initiate_outvar(idt)
     CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
     CALL revise_outvar()
     CALL prepare_to_compute_crit(d,idt,ndt)
  ENDDO time_loop

  !Calculate optimisation criterion and other performance measures
  IF(nacrit/=0) THEN
     IF(PRESENT(condthres))THEN
        CALL calculate_criteria(optcrit,simcrit=performance(:,:),crit1=condcrit,thres1=condthres)
     ELSE
        CALL calculate_criteria(optcrit,simcrit=performance(:,:))
     ENDIF
  ENDIF

  !Prepare observation for yet another simulation        
  CALL reset_observations(modeldir,status)
  IF(status.NE.0) RETURN

  !Set output variables
  criterion = optcrit

  RETURN
END SUBROUTINE run_model_perf

!>Run a simulation of the model with parameters par. 
!>Calculate criterion and performance measures
!>Print out simulation results for all simulations
!-------------------------------------------------------
SUBROUTINE run_model_simout(npar,mpar,par,iens,runens,allens,frozenstate,soilstate,aquiferstate, &
                            riverstate,lakestate,miscstate,criterion,         &
                            performance,status,condcrit,condthres)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : ndt,         &
                       nacrit,      &
                       maxperf,     &
                       maxsubass,   &
                       writematlab, &
                       psdates,    &
                       modeldir,    &
                       resdir,  &
                       noutput,  &
                       output
  USE MODELMODULE, ONLY : model,  &
                          set_modelconfig_from_parameters, &
                          initiate_model, &
                          initiate_model_state
  USE MODVAR, ONLY : nsub,naquifers
  USE COMPOUT, ONLY : compute_mapoutput, &
                      prepare_to_compute_crit,  &
                      calculate_criteria
  USE TIMEROUTINES, ONLY : calculate_time_for_model
  USE LIBDATE, ONLY : DateType
  USE DATAMODULE

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: npar          !<Number of parameters
  INTEGER, INTENT(IN) :: mpar          !<Dimension of parameters-array
  REAL, INTENT(IN)    :: par(mpar)     !<Parameter values to be used for this simulation
  INTEGER, INTENT(IN) :: iens          !<Current simulation
  LOGICAL, INTENT(IN) :: runens        !<Flag for ensemble simulation
  LOGICAL, INTENT(IN) :: allens        !<Flag for writing all ensemble results
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  REAL, INTENT(OUT)   :: criterion       !<Value of optimization criterion
  REAL, INTENT(OUT)   :: performance(maxperf,nacrit)    !<Simulation performance criteria
  INTEGER, INTENT(OUT):: status        !<Subroutine error status
  REAL,OPTIONAL,INTENT(OUT) :: condcrit  !<conditional criterion value
  REAL,OPTIONAL,INTENT(OUT) :: condthres !<threshold of conditional criterion
  
  !Local variables
  REAL optcrit
  TYPE(DateType) d       
  INTEGER idt,iout
  INTEGER :: nmapperiod     !Number of periods for map print out

  !Initialisations
  status = 0
  optcrit = 0.
  nmapperiod = 0 !initialization needed if no mapoutput
  CALL prepare_outputfiles(resdir,nsub,naquifers,iens,runens,allens)
  CALL set_optim_modpar(npar,mpar,par)
  CALL initiate_output_routines()     !All output accumulation variables zeroed
  CALL set_modelconfig_from_parameters()
  CALL initiate_model_state(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  CALL initiate_model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)  
  
  !Run model simulation
  time_loop: DO idt = 1,ndt

     !Get current input data
     CALL get_current_forcing(idt,nsub,d)
     CALL calculate_time_for_model(idt,d)
     IF(ALLOCATED(psdates)) CALL get_current_pointsources(modeldir,'PointSourceData.txt',nsub,d,status) 
     IF(status.NE.0) RETURN

     CALL initiate_outvar(idt)
     CALL model(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
     CALL revise_outvar()
     CALL prepare_to_compute_crit(d,idt,ndt)
     
     !Calculate and write time dependent output
     DO iout = 1,noutput
       IF(output(iout)%fileformat==4) CALL write_regionfiles(iout,idt,ndt,iens,d)
       IF(output(iout)%fileformat==1) CALL write_subbasinfiles(iout,idt,ndt,iens,d)
       IF(output(iout)%fileformat==3) CALL write_timefiles(iout,idt,ndt,iens,d)
       IF(output(iout)%fileformat==2) CALL compute_mapoutput(d,iout,idt,ndt,writematlab,nmapperiod)
     ENDDO

  ENDDO time_loop

  !Calculate optimisation criterion and other performance measures
  IF(nacrit/=0) THEN
    IF(PRESENT(condthres))THEN
      CALL calculate_criteria(optcrit,simcrit=performance(:,:),crit1=condcrit,thres1=condthres)
    ELSE
      CALL calculate_criteria(optcrit,simcrit=performance(:,:))
    ENDIF
  ENDIF

  !Prepare observation for yet another simulation        
  CALL reset_observations(modeldir,status)
  IF(status.NE.0) RETURN

  !Write results to files
  CALL save_mapfiles(resdir,nsub,nmapperiod,iens,runens,allens)

  !Set output variables
  criterion = optcrit

  RETURN
END SUBROUTINE run_model_simout

!Subroutines for Differential Evolution-Markov Chain simulations
!---------------------------------------------------------------------------------------------
  
!>\brief Perform a DE-MC simulation
!>
!> DE-MC introduces uncertainty estimate in the DE optimization algorithm by applying the Metropolis-
!> Hastings acceptance criteria, and by adding a random number in the DE proposal generation. Another
!> way of describing DE-MC is that it is a simplified version of MCMC, whee the DE proposal generation is 
!> used instead of the otherwise cumbersome MCMC jumps based on covariance and multivariate normal distributions. 
!> The genetic DE algorithm overcomes the difficult jump generation by generating the proposal directly from the 
!> current populations, instead of having to tune the covariance matrix and all that stuff.
!> The advantage of DE_MC versus plain DE is both the possibility to get a probability based uncertainty estimate
!> and a better convergence towards the global optima. \n
!
! The following code is a new implementation following the description of DE-MC in Braak (2006). The DE description
! in Ronkkonen et al. (2009) was also useful for the implementation.
!
! The DE-MC implementation was part of the Arctic-HYPE development 2012-2013.
!
! Variables and functions for random number etc are following the previous MonteCarlo routines further below.
!
!>\b References: \n
!> Cajo J. F. Ter Braak (2006). A Markov Chain Monte Carlo version of the genetic algorithm 
!>          Differential Evolution: easy Bayesian computing for real parameter spaces. 
!>          Stat Comput, 16:239-249, DOI 10.1007/s11222-006-8769-1 \n
!> Ronkkonen et al. (2009). 
!--------------------------------------------------------------------------------------------------------------------
SUBROUTINE DEMC_simulation(dir,writeall,frozenstate,soilstate,aquiferstate,riverstate, &
                           lakestate,miscstate,npar)
      USE STATETYPE_MODULE
      USE WORLDVAR, ONLY : nacrit,maxperf,            &
                           writematlab,               &
                           numoptimpar,               &
                           filename_MC,               &
                           fileunit_MC,               &
                           allocate_MCvariables,      &
                           optim
      USE DATAMODULE, ONLY : prepare_save_all_simulations, &
                             write_simulation_results

      IMPLICIT NONE
      
!Argument declaration
      CHARACTER(LEN=*), INTENT(IN) :: dir  !<File directory
      LOGICAL, INTENT(IN)  :: writeall     !<Write result to file
      TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
      TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
      TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
      TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
      TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
      TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
      INTEGER, INTENT(OUT) :: npar         !<Number of parameters to be changed

!Parameter declaration

!Variable declarations
      INTEGER  :: igen,jpop,isim                       !loop counters
      REAL     :: gamma,sigma,crossover
      INTEGER  :: R1, R2, npop, ngen
      INTEGER  :: iacc
      INTEGER  :: status
      REAL     :: parmin(numoptimpar),parmax(numoptimpar) !Parameter interval limits
      REAL     :: parprop(numoptimpar)                    !Proposal parameter vector
      REAL     :: parprec(numoptimpar)                    !Precision sought for parameters
      REAL     :: optcrit,condcrit,condthres              !Optimization criterion, conditional criterion and threshold
      REAL, ALLOCATABLE :: performance(:,:)
      
      
!Initiate MonteCarlo simulation
      IF(.NOT.ALLOCATED(performance)) ALLOCATE(performance(maxperf,nacrit))
      CALL find_optpar(parprop,parmin,parmax,parprec,npar)  !find parameters and limits, sets npar to numoptimpar
      CALL allocate_MCvariables(optim%DEMC_npop+1,npar,maxperf,nacrit) ! +1 row for the median parameters
      CALL random_seed()                                !initiate random numbers
      IF(writeall)THEN
        CALL prepare_save_all_simulations(dir,filename_MC,fileunit_MC,npar,maxperf,nacrit)
      ENDIF

!DE-MC parameters
     npop      = optim%DEMC_npop                       ! number of populations
     gamma     = optim%DEMC_gammascale * 2.38 / sqrt(2.*npar)   ! mutation scaling factor
     crossover = optim%DEMC_crossover                  ! crossover probability
     sigma     = optim%DEMC_sigma                      ! standard deviation of sampling error
     ngen      = optim%DEMC_ngen                       ! number of generations
     IF(npop.LT.npar)THEN
       WRITE(6,*)'WARNING[DEMC]: The number of populations',npop,'is smaller than number of calibrated parameters. You might consider to increase the DEMC_npop parameter in optpar.txt'
     ENDIF
     IF(gamma.GT.0.9)THEN
       WRITE(6,*)'WARNING[DEMC]: gamma = gamma_scale * 2.38/sqrt(2*npar) = ',gamma,' is larger than 0.9, you might consider decreasing the gamma_scale in optpar.txt'
     ENDIF

!Initialize simulation counter
    isim = 0 ! simulation counter
    igen = 0 ! generation counter
    
!Initial population, generation 0-------------------------------------------------------------
     DO jpop=1,npop
       isim = isim + 1
!Generate initial parameter vector 
       CALL get_randompar(numoptimpar,parmin,parmax,npar,parprop)
!Run model, calculate optimization criteria
       IF(optim%task_writesim)THEN
         CALL run_model_simout(npar,numoptimpar,parprop,isim,optim%task_runens,optim%task_writesim,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status,condcrit,condthres)
       ELSE
         CALL run_model_perf(npar,numoptimpar,parprop,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status,condcrit,condthres)
       ENDIF
       IF(status.NE.0)THEN
         WRITE(6,*) 'Error in DE-MC simulation  - while running model for generation 0'
         STOP 1
       ENDIF
!Save results in Population matrices (similar to bookkeepresults...)
       CALL DEMC_acceptreject_proposal(npar,parprop,optcrit,nacrit,performance,jpop,iacc,condcrit,condthres)
!Print results to allsim.txt
       IF(writeall)THEN
         CALL write_simulation_results(fileunit_MC,isim,numoptimpar,npar,maxperf,nacrit,optcrit,performance,parprop,writematlab,jpop,0,iacc)
       ENDIF
!Start printing progress to calibrationlog
       CALL write_DEMC_calibration_log(igen, jpop)
     ENDDO

!Loop over Generations
     DO igen = 1,ngen-1
    
!Loop over populations
       DO jpop = 1,npop
         isim = isim + 1
!Random mutation population R1 and R2; without replacement excluding j
         CALL DEMC_draw_R1R2(jpop,npop,R1,R2)
!Proposal parameter vector, by mutation of R1, R2 and j
         CALL DEMC_proposalgeneration(jpop,R1,R2,gamma,sigma,npar,parprop,parprec)
!Crossover
         CALL DEMC_crossover(parprop,jpop,crossover,npar)
!Control parameter limits
         CALL DEMC_controlparameters(parprop,parmin,parmax,npar)
!Run model, calculate optimization criteria
         IF(optim%task_writesim)THEN
           CALL run_model_simout(npar,numoptimpar,parprop,isim,optim%task_runens,optim%task_writesim,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status,condcrit,condthres)
         ELSE
           CALL run_model_perf(npar,numoptimpar,parprop,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status,condcrit,condthres)
         ENDIF
         IF(status.NE.0)THEN
           WRITE(6,*) 'Error in DE-MC simulation  - while running model for generation', igen,'population',jpop
           STOP 1
         ENDIF
!Accept/Reject proposal - update population matrices   (similar to bookkeepresults...)
         CALL DEMC_acceptreject_proposal(npar,parprop,optcrit,nacrit,performance,jpop,iacc,condcrit,condthres)
!Some printout?
         IF(writeall)THEN
           CALL write_simulation_results(fileunit_MC,isim,numoptimpar,npar,maxperf,nacrit,optcrit,performance,parprop,writematlab,jpop,igen,iacc)
         ENDIF
!Some printout?
         CALL write_DEMC_calibration_log(igen, jpop)
       ENDDO
     ENDDO
!Calculate and write MEDIAN parameterset and results to row 1 in bestMC:s, and the final populations to row 2:npop+1
     CALL DEMC_runmedianparameters(npop,npar,numoptimpar,isim+1,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,condcrit,condthres)

!Finish up
     CLOSE(fileunit_MC)
     IF(ALLOCATED(performance)) DEALLOCATE(performance)
      
END SUBROUTINE DEMC_simulation

!---------------------------------------------------------------------------------------------
SUBROUTINE DEMC_runmedianparameters(npop,npar,numoptimpar,iens,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,condcrit,condthres)

    USE STATETYPE_MODULE
    USE WORLDVAR, ONLY : bestMCparameters, &
                         bestMCoptcrit,          &
                         bestMCperformance,  &
                         bestMCcondcrit, &
                         maxperf,nacrit,  &
                         optim
    USE COMPOUT, ONLY : calculate_median
    IMPLICIT NONE
    
    !Argument declarations
    INTEGER, INTENT(IN) :: npop   !<Number of populations
    INTEGER, INTENT(IN) :: npar   !<Number of parameters to be changed
    INTEGER, INTENT(IN) :: numoptimpar  !<Effective amount of optimization parameters
    INTEGER, INTENT(IN) :: iens         !<Current simulation
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
    REAL, INTENT(INOUT) :: optcrit      !<Optimization criterion
    REAL, INTENT(INOUT) :: performance(maxperf,nacrit)  !<Simulation performance criteria
    REAL, INTENT(INOUT) :: condcrit     !<Conditional criterion
    REAL, INTENT(INOUT) :: condthres    !<Threshold for conditional criterion
    
    !Local variables
    REAL, ALLOCATABLE :: medianpar(:)
    INTEGER           :: i,status

    !Allocate and calculate median vector
    IF(.NOT.ALLOCATED(medianpar)) ALLOCATE(medianpar(npar))
    DO i=1,npar
      CALL calculate_median(npop,bestMCparameters(1:npop,i),-9999.,medianpar(i))
    ENDDO
!Run model, calculate performance and criteria
!Run model with MEDIAN parameter set
    IF(optim%task_writesim)THEN
      CALL run_model_simout(npar,numoptimpar,medianpar,iens,optim%task_runens,optim%task_writesim,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status,condcrit,condthres)
    ELSE
      CALL run_model_perf(npar,numoptimpar,medianpar,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status,condcrit,condthres)
    ENDIF
    IF(status.NE.0)THEN
      WRITE(6,*) 'Error in DE-MC simulation  - while running model for MEDIAN population'
      STOP 1
    ENDIF

!Reshuffle bestMCs
    DO i=1,npop
      bestMCparameters(npop+1-(i-1),:)=bestMCparameters(npop+1-i,:)
      bestMCperformance(npop+1-(i-1),:,:)=bestMCperformance(npop+1-i,:,:)
      bestMCoptcrit(npop+1-(i-1))=bestMCoptcrit(npop+1-i)
      bestMCcondcrit(npop+1-(i-1))=bestMCcondcrit(npop+1-i)
    ENDDO
!Write median parameters and results to Row 1
    bestMCparameters(1,:) = medianpar
    bestMCcondcrit(1) = condcrit
    bestMCoptcrit(1) = optcrit
    bestMCperformance(1,:,:) = performance

!Finish up
    IF(ALLOCATED(medianpar)) DEALLOCATE(medianpar)    
    
END SUBROUTINE DEMC_runmedianparameters

!>Random mutation population R1 and R2; without replacement excluding j
!-----------------------------------------------------------------------
SUBROUTINE DEMC_draw_R1R2(jpop,npop,R1,R2)

  IMPLICIT NONE
  
    !Argument declarations
    INTEGER, INTENT(IN)  :: jpop
    INTEGER, INTENT(IN)  :: npop
    INTEGER, INTENT(OUT) :: R1, R2

    !Local variables
    REAL :: rand1, rand2
    
    !>\b Algorithm \n
    !>Draw a random number R1 and scale it to an integer (1:npop-1)
    CALL RANDOM_NUMBER(rand1)
    rand1 = rand1*(npop-1)
    R1    = min(npop-1,int(rand1-mod(rand1,1.)+1))

    R2 = R1     !Initialize R2 = R1
    !>Draw second random number R2 until R1 not equal to R2    
    DO WHILE(R2.EQ.R1)
      CALL RANDOM_NUMBER(rand2)
      rand2 = rand2*(npop-1)  !Scale it to an integer (1:npop-1)
      R2    = min(npop-1,int(rand2-mod(rand2,1.)+1))
    ENDDO

    !>Re-scale random numbers to integer scale (1:npop) excluding j
    IF(R1.GE.jpop)R1=R1+1
    IF(R2.GE.jpop)R2=R2+1

END SUBROUTINE DEMC_draw_R1R2

!>Proposal parameter vector, by mutation of R1, R2 and j
!---------------------------------------------------------------------------------------------
SUBROUTINE  DEMC_proposalgeneration(jpop,R1,R2,gamma,sigma,npar,parprop,parprec)

    USE WORLDVAR, ONLY : bestMCparameters
    IMPLICIT NONE
    
    !Argument declarations
    INTEGER, INTENT(IN) :: jpop, R1, R2, npar
    REAL, INTENT(IN)    :: gamma, sigma
    REAL, INTENT(INOUT) :: parprop(npar),parprec(npar)
    
    !Local variables
    INTEGER :: kpar
    REAL rand

    !Mutation parameter vector, v = X(j) + gamma*(X(R1)-X(R2)) + epsilon{N(0,sigma)}
    DO kpar=1,npar

      ! Mutation J = J+gamma(R1-R2)
      parprop(kpar) = bestMCparameters(jpop,kpar) + gamma*(bestMCparameters(R1,kpar)-bestMCparameters(R2,kpar))

      ! Add random perturbation (UNIFORM for now, add standard normal distribution later)
      CALL RANDOM_NUMBER(rand) ! random number between 0 and 1
        
      ! scale the random number with the sigma parameter and the precision
      parprop(kpar) = parprop(kpar) + (rand-0.5)*sigma*parprec(kpar)

    ENDDO
END SUBROUTINE DEMC_proposalgeneration

!---------------------------------------------------------------------------------------------
SUBROUTINE DEMC_crossover(parprop,jpop,crossover,npar)

  USE WORLDVAR, ONLY : bestMCparameters
  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) ::  jpop,npar
  REAL,INTENT(INOUT)  ::  parprop(npar)
  REAL, INTENT(IN)    ::  crossover
    
  !Local variables
  REAL                :: rand
  INTEGER             :: fixed_crossover,kpar    

  !Only try crossover in case crossover probability<1
  IF(crossover.LT.1.)THEN
  !Randomly decide ONE parameter that will crossover independent of probability (to get at least one crossover)
    !Draw first random number
    CALL RANDOM_NUMBER(rand)
    !Scale it to an integer (1:npar)
    rand = rand*(npar)
    fixed_crossover    = min(npar,int(rand-mod(rand,1.)+1))
    
    !Loop over paramteters and reject crossover if random number is above probability for crossover
    DO kpar=1,npar
      !Crossover rejection is only needed for the non-fixed parameter
      IF(kpar.NE.fixed_crossover)THEN
        !Draw first random number
        CALL RANDOM_NUMBER(rand)
        ! Reject crossover, replace proposal with the current population
        IF(rand.GE.crossover)parprop(kpar) = bestMCparameters(jpop,kpar)
      ENDIF
    ENDDO
  ENDIF
END SUBROUTINE DEMC_crossover

!>Control parameter limits, and reflect proposal into the range if outside.  
!>If still outside, make a relative reflection
!--------------------------------------------------------------------------
SUBROUTINE DEMC_controlparameters(parprop,parmin,parmax,npar)
    
    IMPLICIT NONE

    !Argument declaration
    INTEGER, INTENT(IN) :: npar
    REAL, INTENT(INOUT) :: parprop(npar)
    REAL, INTENT(IN)    :: parmin(npar),parmax(npar)

    !Local variables
    INTEGER             :: kpar
    REAL                :: ratio
    
!Loop over parameters
    DO kpar = 1, npar
!Reflect proposal into the Min-Max range if outside, following Braak (2006)
        IF(parprop(kpar).LT.parmin(kpar))parprop(kpar) = parmin(kpar) + (parmin(kpar)-parprop(kpar))
        IF(parprop(kpar).GT.parmax(kpar))parprop(kpar) = parmax(kpar) - (parprop(kpar)-parmax(kpar))
!If proposal is still outside (due to the random term?) - make sure its inside
!Use the ratio of the deviation from the range and tha range itself for reflection:
        IF(parprop(kpar).LT.parmin(kpar))THEN
            ratio = (parmin(kpar)-parprop(kpar)) / (parmax(kpar)-parmin(kpar))
            IF(ratio.GT.1)THEN
                parprop(kpar) = parmin(kpar) +  (parmax(kpar)-parmin(kpar))/ratio
            ELSE
                parprop(kpar) = parmin(kpar) +  (parmax(kpar)-parmin(kpar))*ratio
            ENDIF
        ENDIF
        IF(parprop(kpar).GT.parmax(kpar))THEN
            ratio = (parprop(kpar)-parmax(kpar)) / (parmax(kpar)-parmin(kpar))
            IF(ratio.GT.1)THEN
                parprop(kpar) = parmax(kpar) - (parmax(kpar)-parmin(kpar))/ratio
            ELSE
                parprop(kpar) = parmax(kpar) - (parmax(kpar)-parmin(kpar))*ratio
            ENDIF
        ENDIF
    ENDDO

END SUBROUTINE DEMC_controlparameters

!---------------------------------------------------------------------------------------------
SUBROUTINE DEMC_acceptreject_proposal(npar,par,crit,numcrit,performance,jpop,iacc,condcr,condth)

  USE WORLDVAR, ONLY : bestMCparameters,    &
                       bestMCoptcrit,          &
                       bestMCperformance,  &
                       bestMCcondcrit, &
                       maxperf, &
                       optim      
  
  IMPLICIT NONE
  
    !Argument declaration
    INTEGER, INTENT(IN)  :: npar                         !<Actual number of parameters
    REAL, INTENT(IN)     :: par(npar)                    !<Current parameter values
    REAL, INTENT(IN)     :: crit                         !<Current criterium value
    INTEGER, INTENT(IN)  :: numcrit                      !<Number of variables that criteria are calculated for
    REAL, INTENT(IN)     :: performance(maxperf,numcrit) !<Other performance criteria
    INTEGER, INTENT(IN)  :: jpop                         !<Index of Actual Population 
    INTEGER, INTENT(INOUT)  :: iacc                      !<Acceptance index (0)rejected (1)accepted
    REAL, INTENT(IN)  :: condcr      !<Conditional criteria
    REAL, INTENT(IN)  :: condth      !<Conditional acceptance threshold

    !Local variables
    REAL      :: accprob, rand
    LOGICAL   :: conditional_reject
      
    !Conditional rejection if enebled by input arguments
    IF(condcr.GT.bestMCcondcrit(jpop) .AND. condcr.GT.condth)THEN
      conditional_reject = .TRUE.
      iacc = 0
    ELSE
      conditional_reject = .FALSE.
    ENDIF

    !Go on and apply acceptance criteria if not conditionally rejected
    IF(.NOT.conditional_reject)THEN
      !Always accept if criteria is better, conditional on the conditional criteria being also better or above threshold  
      IF(crit<bestMCoptcrit(jpop))THEN
        bestMCcondcrit(jpop) = condcr
        bestMCoptcrit(jpop) = crit
        bestMCparameters(jpop,:) = par(1:npar)
        bestMCperformance(jpop,:,:) = performance
        iacc = 1
      ELSE
         
        !Enable a choice to apply a probabilistic acceptance rule
        IF(optim%DEMC_accprob.GT.0.)THEN
          !For now, we us the CRIT as defined in HYSS, which is always minimized (versus minus infinity)
          !The acceptance probability is created by taking 1-ratio(deviation_from_crit/crit)
          accprob = optim%DEMC_accprob*(1.-(crit-bestMCoptcrit(jpop))/abs(bestMCoptcrit(jpop)))
          !Random number
          CALL RANDOM_NUMBER(rand)
          IF(rand.lt.accprob)then
            bestMCcondcrit(jpop) = condcr
            bestMCoptcrit(jpop) = crit
            bestMCparameters(jpop,:) = par(1:npar)
            bestMCperformance(jpop,:,:) = performance
            iacc = 1
          ELSE
            iacc = 0
          ENDIF
        ELSE     
          iacc = 0
        ENDIF
      ENDIF
    ENDIF      
      
    RETURN
    
END SUBROUTINE DEMC_acceptreject_proposal

!>Write progress of DEMC simulation to Calibration.log
!-------------------------------------------------------------------------------------
SUBROUTINE write_DEMC_calibration_log(genCounter,popCounter)

      USE WORLDVAR, ONLY : optim,                   &
                           modeldir, resdir,        &
                           optimStartTime,          &
                           filename_callog,         &   
                           fileunit_callog

      IMPLICIT NONE

!Argument declaration
      INTEGER, INTENT(IN)  :: genCounter, popCounter

!Variable declarations
      INTEGER  iterationsDone, iterationsTotal
      INTEGER  secondsAvgTime, minuteAvgTime, secondsRmngTime, minuteRmngTime, hourRmngTime, daysRmngTime
      REAL nowTime, averageTime, remainingTime


      iterationsTotal = optim%DEMC_ngen * optim%DEMC_npop
      iterationsDone = genCounter*optim%DEMC_npop + popCounter   !genCounter starts at 0)

! Time stuff
      CALL cpu_time(nowTime)

      averageTime = (nowTime - optimStartTime)/iterationsDone
      IF(averageTime > 60) THEN
        minuteAvgTime = floor(averageTime/60)
        secondsAvgTime = nint(averageTime - 60*minuteAvgTime)
      ENDIF

      remainingTime = (iterationsTotal - iterationsDone) * averageTime
      IF(remainingTime > 60) THEN
        daysRmngTime = floor(remainingTime/86400)
        hourRmngTime = floor((remainingTime - 86400*daysRmngTime)/3600)
        minuteRmngTime = floor((remainingTime - 86400*daysRmngTime - 3600*hourRmngTime)/60)
        secondsRmngTime = nint(remainingTime - 86400*daysRmngTime - 3600*hourRmngTime - 60*minuteRmngTime)
      ENDIF

! Start printing log file
      OPEN(UNIT = fileunit_callog, FILE = TRIM(resdir)//filename_callog, STATUS = 'replace')   ! Replace calibration log file
      WRITE(fileunit_callog,'(A)') 'Calibration by DEMC (Differential Evolution Markov Chain)'
      WRITE(fileunit_callog,'(A)') '---------------------------------------------------------'
      WRITE(fileunit_callog,*) ''
      WRITE(fileunit_callog,'(A,A)') 'Model area: ', TRIM(modeldir)
      WRITE(fileunit_callog,*) ''
      WRITE(fileunit_callog,'(I6,A,I6,A,I7,A)') optim%DEMC_ngen, ' generations, each with ', optim%DEMC_npop, ' populations  =>  total calibration work: ', iterationsTotal, ' iterations'
      WRITE(fileunit_callog,'(I6,A,I6,A)') optim%DEMC_npop, ' populations (parameter sets) are runned in parallel chains which are mutated and crossed over',optim%DEMC_ngen-1,'times following th DEMC algorithm (Braak, 2006).'
      WRITE(fileunit_callog,*) ''
      WRITE(fileunit_callog,*) ''
      WRITE(fileunit_callog,'(A)') 'Current progress:'
      WRITE(fileunit_callog,'(A)') '-----------------'
      WRITE(fileunit_callog,'(A,I6,A,I6)') 'Generation:  ', genCounter, '/', optim%DEMC_ngen
      WRITE(fileunit_callog,'(A,I6,A,I6)') 'Population:  ', popCounter, '/', optim%DEMC_npop
      WRITE(fileunit_callog,'(A,I6,A,I6)') 'Global: ', iterationsDone, '/', iterationsTotal
      WRITE(fileunit_callog,*) ''

      IF(averageTime > 60) THEN
        WRITE(fileunit_callog,'(A,I2,A,I2)') 'Average iteration time: ', minuteAvgTime, ':', secondsAvgTime
      ELSE
        WRITE(fileunit_callog,'(A,F6.3,A)') 'Average iteration time: ', averageTime, ' seconds'
      ENDIF

      IF(remainingTime > 60) THEN
        IF(daysRmngTime > 0) THEN
          WRITE(fileunit_callog, '(A,I2,A,I2,A,I2,A,I2)') 'Remaining calibration time: ', daysRmngTime, ' days and ', hourRmngTime, ':', minuteRmngTime, ':', secondsRmngTime
        ELSE
          WRITE(fileunit_callog, '(A,I2,A,I2,A,I2)') 'Remaining calibration time: ', hourRmngTime, ':', minuteRmngTime, ':', secondsRmngTime
        ENDIF
      ELSE
        WRITE(fileunit_callog, '(A,F6.3,A)') 'Remaining calibration time: ', remainingTime, ' seconds'
      ENDIF
      
      WRITE(fileunit_callog,*) ''
      WRITE(fileunit_callog,*) ''
      WRITE(fileunit_callog, '(A)') 'Current parameter values in each population:'
      WRITE(fileunit_callog, '(A)') '--------------------------------------------'
      WRITE(fileunit_callog,*) ''

      CALL write_DEMC_calibration_log_tail(fileunit_callog)
      CLOSE(fileunit_callog)

END SUBROUTINE write_DEMC_calibration_log

!--------------------------------------------------------
SUBROUTINE write_DEMC_calibration_log_tail(funit)

      USE WORLDVAR, ONLY : bestMCparameters,        &
                           bestMCoptcrit,           &
                           bestMCcondcrit,          &
                           optim,                   &
                           optparid, parindex,      &
                           numoptimpar,             &
                           filename_callog

      USE MODVAR, ONLY : modparid

      IMPLICIT NONE

      !Argument declaration
      INTEGER, INTENT(IN)  :: funit !<file unit of calibration.log

      !Local variables
      INTEGER  columnCounter, lineCounter
      CHARACTER(LEN=20) junkTxt
      CHARACTER(LEN=20*(numoptimpar+2)) fullTxtLine


! Prepare table top line: parameter labels + criteria
      DO columnCounter = 0,numoptimpar-1
        fullTxtLine(20*columnCounter+ 1:20*columnCounter+ 5) = '     '
        fullTxtLine(20*columnCounter+ 6:20*columnCounter+15) = modparid(optparid(parindex(columnCounter+1,1)))%shortname
        fullTxtLine(20*columnCounter+16:20*columnCounter+20) = '     '
      ENDDO
      fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = '     CRIT           '
      fullTxtLine(20*(columnCounter+1)+1:20*(columnCounter+1)+20) = '     CONDCRIT       '
      WRITE(funit,'(A)') fullTxtLine

! Then print all parameter and criteria values
      DO lineCounter = 1,optim%DEMC_npop

        DO columnCounter = 0,numoptimpar-1
          WRITE(junkTxt, '(F20.9)') bestMCparameters(lineCounter, columnCounter+1)  ! Put into "junkTxt" the string conversion of the parameter values
          fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = junkTxt            ! Concatenate in "fullTxtLine"
        ENDDO

        columnCounter = numoptimpar
        WRITE(junkTxt, '(F20.9)') bestMCoptcrit(lineCounter)        ! Add the criteria value
        fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = junkTxt
        
        columnCounter = columnCounter+1
        WRITE(junkTxt, '(F20.9)') bestMCcondcrit(lineCounter)        ! Add the criteria value
        fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = junkTxt

        WRITE(funit,'(A)') fullTxtLine        ! Output as a string to the calibration log
      ENDDO


END SUBROUTINE write_DEMC_calibration_log_tail
  
!---------------------------------------------------------------------------------------------
!Subroutines for MonteCarlo simulations
!---------------------------------------------------------------------------------------------

!>Simple MonteCarlo simulation with random parameter values
!---------------------------------------------------------------------------------------------
SUBROUTINE MonteCarlo_simulation(dir,writeall,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,npar) 

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : nacrit,          &
                       maxperf,         &
                       writematlab,     &
                       numoptimpar,     &
                       filename_MC,     &
                       fileunit_MC,     &
                       allocate_MCvariables,      &
                       optim
  USE DATAMODULE, ONLY : prepare_save_all_simulations, &
                         write_simulation_results
  
  IMPLICIT NONE

  !Argument declarations
  CHARACTER(LEN=*), INTENT(IN) :: dir                 !<File directory
  LOGICAL, INTENT(IN)  :: writeall                    !<Write result to file
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  INTEGER, INTENT(OUT) :: npar                        !<Number of parameters to be changed

  !Local variables
  INTEGER  :: i
  INTEGER  :: status
  REAL     :: parmin(numoptimpar),parmax(numoptimpar) !Parameter interval limits
  REAL     :: par(numoptimpar)                        !Current parameter values
  REAL     :: parprec(numoptimpar)                    !Precision sought for parameters
  REAL     :: optcrit
  REAL, ALLOCATABLE :: performance(:,:)

  !>\b Algorithm \n
  !>Initiate MonteCarlo simulation; find parameters and intervals for them to use
  IF(.NOT.ALLOCATED(performance)) ALLOCATE(performance(maxperf,nacrit))
  CALL find_optpar(par,parmin,parmax,parprec,npar)  !find parameters and limits, sets npar to numoptimpar
  CALL allocate_MCvariables(optim%nruns_best,npar,maxperf,nacrit)
  CALL random_seed()                                !initiate random numbers
  IF(writeall)THEN
    CALL prepare_save_all_simulations(dir,filename_MC,fileunit_MC,npar,maxperf,nacrit)
  ENDIF

  !>Simulate model num_mc times, generating new random parameter values for each simulation
  DO i = 1,optim%nruns_MC
    IF(MOD(i,10000)==0) WRITE(6,*) 'No:',i
    CALL get_randompar(numoptimpar,parmin,parmax,npar,par)
    IF(optim%task_writesim)THEN
      CALL run_model_simout(npar,numoptimpar,par,i,optim%task_runens,optim%task_writesim,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      !Evaluate the function; simulate model
    ELSE
      CALL run_model_perf(npar,numoptimpar,par,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      !Evaluate the function; simulate model
    ENDIF
    IF(status.NE.0)THEN
      WRITE(6,*) 'Error in MC simulation'
      STOP 1
    ENDIF
    !>Meanwhile keep track of the num_ens best simulations for final output
    CALL bookkeep_result_from_simulation(npar,numoptimpar,optim%nruns_best,par,optcrit,nacrit,performance)
    IF(writeall)THEN
      CALL write_simulation_results(fileunit_MC,i,numoptimpar,npar,maxperf,nacrit,optcrit,performance,par,writematlab)
    ENDIF
  ENDDO

  !Finish up
  CLOSE(fileunit_MC)
  IF(ALLOCATED(performance)) DEALLOCATE(performance)

  RETURN
END SUBROUTINE MonteCarlo_simulation

!>Successive MonteCarlo simulation with reduced parameter space
!---------------------------------------------------------------
SUBROUTINE bounded_MonteCarlo_simulation(taskMC,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,npar) 

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : nacrit,maxperf,         &
                       numoptimpar,            &
                       allocate_MCvariables,   &
                       deallocate_MCvariables, &
                       optim

  IMPLICIT NONE
  
  !Argument declarations
  LOGICAL, INTENT(IN)  :: taskMC  !<Flag for MonteCarlo simulation done
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  INTEGER, INTENT(OUT) :: npar    !<Number of parameters to be changed

  !Local variables 
  INTEGER  :: i, j
  INTEGER  :: status
  REAL     :: parmin(numoptimpar),parmax(numoptimpar) !Parameter interval limits
  REAL     :: par(numoptimpar)                        !Current parameter values
  REAL     :: parprec(numoptimpar)                    !Precision sought for parameters
  REAL     :: optcrit
  REAL, ALLOCATABLE :: performance(:,:)

  !>\b Algorithm \n
  IF(.NOT.ALLOCATED(performance)) ALLOCATE(performance(maxperf,nacrit))
  IF(.NOT.taskMC) THEN
    !>Initiate MonteCarlo simulation; find parameters and intervals for them to use
    CALL find_optpar(par,parmin,parmax,parprec,npar)    !find parameters and limits
    CALL allocate_MCvariables(optim%nruns_best,npar,maxperf,nacrit)
    CALL random_seed()                                  !initiate random numbers
    !>Run initial MonteCarlo simulation; num_bpmc simulations, generating new random parameter values for each simulation
    DO i = 1,optim%nruns_MCI
      CALL get_randompar(numoptimpar,parmin,parmax,npar,par)
      IF(optim%task_writesim)THEN
        CALL run_model_simout(npar,numoptimpar,par,i,optim%task_runens,optim%task_writesim,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      !Evaluate the function; simulate model
      ELSE
        CALL run_model_perf(npar,numoptimpar,par,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      !Evaluate the function; simulate model
      ENDIF
      IF(status.NE.0)THEN
        WRITE(6,*) 'Error MonteCarlo reduce parameter space simulation'
        STOP 1
      ENDIF
      !>Meanwhile keep track of the num_ens best simulations for next stage
      CALL bookkeep_result_from_simulation(npar,numoptimpar,optim%nruns_best,par,optcrit,nacrit,performance)
    ENDDO
  ENDIF
  !>Start MonteCarlo simulation with reducing bounds successively; for each repetition (of num_bpmax)...
  DO j = 1,optim%nruns_MCImax              
    !>\li Reducing the parameter space to contain the best simulations only
    CALL reduce_parameter_space(numoptimpar,parmin,parmax,npar)
    !>\li Run a MonteCarlo simulation; num_bpmc simulations, generating new random parameter values for each simulation
    DO i = 1,optim%nruns_MCI
      CALL get_randompar(numoptimpar,parmin,parmax,npar,par)
      IF(optim%task_writesim)THEN
        CALL run_model_simout(npar,numoptimpar,par,j*optim%nruns_MCI+i,optim%task_runens,optim%task_writesim,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      !Evaluate the function; simulate model
      ELSE
        CALL run_model_perf(npar,numoptimpar,par,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      !Evaluate the function; simulate model
      ENDIF
      IF(status.NE.0)THEN
        WRITE(6,*) 'Error MonteCarlo reduce parameter space simulation'
        STOP 1
      ENDIF
      !>\li Meanwhile keep track of the num_ens best simulations for next stage
      CALL bookkeep_result_from_simulation(npar,numoptimpar,optim%nruns_best,par,optcrit,nacrit,performance)
    ENDDO
    !condition to EXIT
  ENDDO
  IF(ALLOCATED(performance)) DEALLOCATE(performance)

END SUBROUTINE bounded_MonteCarlo_simulation

!>Get new random values of parameters to be optimized
!>from a uniform distribution between parmin and parmax
!---------------------------------------------------------------------------------------------
SUBROUTINE get_randompar(mpar,parmin,parmax,npar,par) 

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN)  :: mpar         !<Dimension of argument parameters
  REAL, INTENT(IN)     :: parmin(mpar) !<Parameter interval lower limit
  REAL, INTENT(IN)     :: parmax(mpar) !<Parameter interval upper limit
  INTEGER, INTENT(IN)  :: npar         !<Actual number of parameters to get random number
  REAL, INTENT(OUT)    :: par(mpar)    !<Current parameter values
  
  !Local variables
  REAL     :: rand(npar)

  par = 0.
  CALL RANDOM_NUMBER(rand)
  par(1:npar)=rand*(parmax(1:npar)-parmin(1:npar))+parmin(1:npar)

END SUBROUTINE get_randompar

!>Determine new bounds for the parameter space based on the best simulations sofar
!---------------------------------------------------------------------------------------------
SUBROUTINE reduce_parameter_space(mpar,parmin,parmax,npar) 
  USE WORLDVAR, ONLY : bestMCparameters

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: mpar         !<Dimension of parameters
  REAL, INTENT(INOUT)  :: parmin(mpar) !<Parameter interval lower limit
  REAL, INTENT(INOUT)  :: parmax(mpar) !<Parameter interval upper limit
  INTEGER, INTENT(IN)  :: npar         !<Actual number of parameters

  parmin(1:npar) = MINVAL(bestMCparameters,1)
  parmax(1:npar) = MAXVAL(bestMCparameters,1)

END SUBROUTINE reduce_parameter_space

!>Save the results from the simulation for later use. Save the sofar nbest best and all values.
!---------------------------------------------------------------------------------------------
SUBROUTINE bookkeep_result_from_simulation(npar,mpar,nbest,par,crit,numcrit,performance)

  USE WORLDVAR, ONLY : bestMCparameters,    &
                       bestMCoptcrit, &
                       bestMCperformance, &
                       maxperf

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: npar                         !<Actual number of parameters
  INTEGER, INTENT(IN) :: mpar                         !<Dimension of parameter-array
  INTEGER, INTENT(IN) :: nbest                        !<Number of ensembles saved as best
  REAL, INTENT(IN)    :: par(mpar)                    !<Current parameter values
  REAL, INTENT(IN)    :: crit                         !<Current optimization criterion value
  INTEGER, INTENT(IN) :: numcrit                      !<Number of variables that performance criteria are calculated for
  REAL, INTENT(IN)    :: performance(maxperf,numcrit) !<Performance criteria
  
  !Local variables
  INTEGER i

  !Save best iterations in order
  DO i = 1,nbest
    IF(crit<bestMCoptcrit(i))THEN
      IF(i<nbest)THEN
        bestMCparameters(i+1:nbest,:) = bestMCparameters(i:nbest-1,:)
        bestMCoptcrit(i+1:nbest) = bestMCoptcrit(i:nbest-1)
        bestMCperformance(i+1:nbest,:,:) = bestMCperformance(i:nbest-1,:,:)
      ENDIF
      bestMCoptcrit(i) = crit
      bestMCparameters(i,:) = par(1:npar)
      bestMCperformance(i,:,:) = performance
      EXIT
    ENDIF
  ENDDO

END SUBROUTINE bookkeep_result_from_simulation

!------------------------------------------------------------------------------------------------------------
!                               Quasi-Newton autocal routines
!------------------------------------------------------------------------------------------------------------


!---------------------
!Scanning routine
!---------------------

!>This routine evaluates the criteria for systematic variation of two
!>parameters, with two constant parameter stepsizes
!------------------------------------------------------------------------
SUBROUTINE param_scanning(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : modeldir,        &
                       resdir,          &
                       optim,           &
                       optparid,        &
                       numoptimpar,     &
                       filename_callog, &
                       fileunit_callog
  USE MODVAR, ONLY   : modparid

  IMPLICIT NONE
  
  !Argument declarations
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  
  !Local variables
  INTEGER counterX, counterY, npar
  INTEGER resultat(numoptimpar)
  REAL crit
  REAL par_center(numoptimpar), par_minvalue(numoptimpar), par_maxvalue(numoptimpar), par_precision(numoptimpar)
  REAL delta(numoptimpar), direction(numoptimpar)
  REAL criteria_matrix(optim%scan_xpoints, optim%scan_ypoints)

  !>\b Algorithm \n
  !>Check number parameters is ok, and make sure results will be output to "calibration.log"
  IF(numoptimpar .NE. 2) THEN
    WRITE(6,*) 'ERROR: Problem with scanning mode; scanning can only use 2 parameters'
    STOP 1
  ENDIF
  optim%cal_log = .TRUE.

  !Scanning procedure
  criteria_matrix = 0.0

  !>Initiate MonteCarlo simulation; find parameters and intervals for them to use
  CALL find_optpar(par_center, par_minvalue, par_maxvalue, par_precision, npar)
  !>Calculate the step in each direction for gridded coverage of parameter space
  delta(1) = (par_maxvalue(1) - par_minvalue(1))/(optim%scan_xpoints - 1)
  delta(2) = (par_maxvalue(2) - par_minvalue(2))/(optim%scan_ypoints - 1)

  !>Run simulation for each parameter value in grid
  DO counterX = 1, optim%scan_xpoints
    DO counterY = 1, optim%scan_ypoints
      direction(1) = (counterX-1) * delta(1)
      direction(2) = (counterY-1) * delta(2)
      CALL function_to_minim(1.0,crit,par_minvalue,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
      criteria_matrix(counterX, counterY) = crit
    ENDDO
  ENDDO

  !Find minimum
  resultat = MINLOC(criteria_matrix)

  !>Write result to calibration.log
  OPEN(UNIT = fileunit_callog, FILE = TRIM(resdir)//filename_callog, STATUS = 'replace') !Start calibration log file
  WRITE(fileunit_callog, '(A)') 'Parameter space scanning'
  WRITE(fileunit_callog, '(A)') '------------------------'
  WRITE(fileunit_callog, '(A)') ' '
  WRITE(fileunit_callog, '(A, A)') 'Model area: ', TRIM(modeldir)
  WRITE(fileunit_callog, '(A)') ' '
  WRITE(fileunit_callog, '(A)') ' '

  WRITE(fileunit_callog, '(A47)') modparid(optparid(2))%shortname
  WRITE(fileunit_callog, '(A)', advance='no') '                          '
  DO counterY = 1, optim%scan_ypoints
    WRITE(fileunit_callog, '(F17.9, A)', advance='no') par_minvalue(2) + (counterY-1)*delta(2), '    '
  ENDDO
  WRITE(fileunit_callog, '(A)') ' '
  WRITE(fileunit_callog, '(A22)') modparid(optparid(1))%shortname

  DO counterX = 1, optim%scan_xpoints
    WRITE(fileunit_callog, '(F17.9, A)', advance='no') par_minvalue(1) + (counterX-1)*delta(1), '         '
    DO counterY = 1, optim%scan_ypoints
      WRITE(fileunit_callog, '(F17.9, A)', advance='no') criteria_matrix(counterX, counterY), '    '
    ENDDO
    WRITE(fileunit_callog, '(A)') ' '
  ENDDO
  WRITE(fileunit_callog, '(A)') ' '
  WRITE(fileunit_callog, '(A)') ' '

  WRITE(fileunit_callog, '(A, F17.9)') 'Minimum value: ', criteria_matrix(resultat(1), resultat(2))
  WRITE(fileunit_callog, '(A15, F17.9)') modparid(optparid(1))%shortname, par_minvalue(1) + (resultat(1)-1)*delta(1)
  WRITE(fileunit_callog, '(A15, F17.9)') modparid(optparid(2))%shortname, par_minvalue(2) + (resultat(2)-1)*delta(2)

  STOP 0

END SUBROUTINE param_scanning

!>\brief MonteCarlo optimisation method, stage-wise zooming on top results
!! This is the "outer" part of the function, core part is defined as
!! separate subroutine (see next)
!-------------------------------------------------------------------------
SUBROUTINE stage_MonteCarlo(frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : nacrit, maxperf, &
                       filename_MC, fileunit_MC, &
                       allocate_MCvariables, &
                       resdir, &
                       optim, numoptimpar, bestMCparameters
  USE DATAMODULE, ONLY : prepare_save_all_simulations

  IMPLICIT NONE
  
  !Argument declarations
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  
  !Local variables
  INTEGER stageCounter, amountStageRuns, bestCounter
  INTEGER npar
  REAL    par_minvalue(numoptimpar),par_maxvalue(numoptimpar)   !Lower/upper interval limits, from optpar.txt; those limits are absolute, serve as permanent/reference bounds!!
  REAL    par_precision(numoptimpar)
  REAL    parRad(numoptimpar)       !Vector of "radii" in all optim param dimensions at start (reference, not to be modified!)
  REAL    parRadStage(numoptimpar)  !Vector of "radii" for each stage (zooming in, changing at each stage)
  REAL    par_center(numoptimpar)    !Mean of the boundaries, current (starting) center at begining; becomes copy of best later on
  LOGICAL writeall                  !Write result to file (logical from optim%task_writeall, corresponds to flag 'WA' in optpar.txt)

  !>\b Algorithm \n
  WRITE(6,'(A)') 'Calibration with progressive Monte Carlo method'
  writeall = optim%task_writeall

  !>First of all, check that the numerical zoom factor is not larger than 1 (otherwise the method will be extending the domain, not reducing it!)
  IF(optim%nruns_zoom > 1) THEN         !Leave a message in hyss.log and exit in case of trouble
    WRITE(6,*) 'Problem in optpar for stage MonteCarlo: numerical zooming factor (num_zoom) is larger than 1!'
    STOP 1
  ENDIF

  !>Initiate MonteCarlo simulation; find parameters and intervals for them to use
  CALL find_optpar(par_center, par_minvalue, par_maxvalue, par_precision, npar)  !Get interval limits and center from optpar.txt
  parRad = (par_maxvalue - par_minvalue)/2          !Define reference radii array
  CALL allocate_MCvariables(optim%nruns_best, npar, maxperf, nacrit)
  IF(writeall) THEN                     !If 'WA' flag in optpar.txt, initiate the "allsim.txt" file, where ALL MC simulations results are listed (why not just the best ones, sorted?)
    CALL prepare_save_all_simulations(resdir, filename_MC, fileunit_MC, npar, maxperf, nacrit)
  ENDIF
  CALL random_seed()                    !Initialize random sequence

  !>Proceed with simulations, stage by stage, for each stage update the parameter space by reducing the interval for each parameter centered around the best points
  DO stageCounter = 0,optim%nruns_stages-1                  !Loop through the stages (counter offset by 1 for convenience)

    IF(stageCounter == 0) THEN
      amountStageRuns = optim%nruns_MC * optim%nruns_best       !First stage (0-stage) has 'num_ens' times more runs than later stages
      CALL stage_MonteCarlo_core_function(writeall,1,stageCounter+1,amountStageRuns,  &
        par_center,parRad,par_minvalue,par_maxvalue,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
    ELSE
      amountStageRuns = optim%nruns_MC
      parRadStage = parRad * optim%nruns_zoom**stageCounter     !Update the radii with the zooming factor, in accordance with stage progression
      DO bestCounter = 1,optim%nruns_best                       !Update center for this stage
        par_center(1:npar) = bestMCparameters(bestCounter,1:npar)
        CALL stage_MonteCarlo_core_function(writeall,bestCounter,stageCounter+1, &
           amountStageRuns,par_center,parRadStage,par_minvalue,par_maxvalue,frozenstate, &
           soilstate,aquiferstate,riverstate,lakestate,miscstate)
      ENDDO
    ENDIF

    !>Print progression in hyss.log
    WRITE(6,'(A)') ''
    WRITE(6,'(A, I2)') 'Stage ', stageCounter
    CALL write_stageMC_calibration_log_tail(6)

  ENDDO

  !Finish up
  CLOSE(fileunit_MC)        ! Close the "allsim.txt" file

END SUBROUTINE stage_MonteCarlo

!>\brief Core of loop in "stage_MonteCarlo" subroutine.
! Runs 'num_ens' MC simulations with parameters randomly generated with a given radius
! around a center. The best ones are stored in "bestMCparameters", which is a 
! global variable so that nothing is returned to the calling function.
! ----------------------------------------------------------------
SUBROUTINE stage_MonteCarlo_core_function(writeall,icenter,istage, &
                           nruns,parCenter,parRadStage,parMin,parMax, &
                           frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

    USE STATETYPE_MODULE
    USE WORLDVAR, ONLY : nacrit, maxperf,           &
                         writematlab,               &
                         filename_MC,               &
                         fileunit_MC,               &
                         allocate_MCvariables,      &
                         optim, numoptimpar
    USE DATAMODULE, ONLY : write_simulation_results
  
    IMPLICIT NONE
  
    !Argument declarations
    LOGICAL, INTENT(IN) :: writeall         !<flag for write all result to file
    INTEGER, INTENT(IN) :: icenter          !<current center point
    INTEGER, INTENT(IN) :: istage           !<current stage, runs from 1 to optim%nruns_stages
    INTEGER, INTENT(IN) :: nruns            !<amount of MC runs for the current stage
    REAL, INTENT(IN) :: parCenter(numoptimpar)    !<Current param space center within param space for random param generation
    REAL, INTENT(IN) :: parRadStage(numoptimpar)  !<Current param space radii within param space for random param generation
    REAL, INTENT(IN) :: parMin(numoptimpar)       !<Absolute param space minimum boundary
    REAL, INTENT(IN) :: parMax(numoptimpar)       !<Absolute param space maximum boundary
    TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)    :: soilstate   !<Soil states
    TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
    TYPE(riverstatetype),INTENT(INOUT)   :: riverstate  !<River states
    TYPE(lakestatetype),INTENT(INOUT)    :: lakestate   !<Lake states
    TYPE(miscstatetype),INTENT(INOUT)    :: miscstate   !<Misc states
  
    !Local variables
    INTEGER  irun, status
    INTEGER  iterations_done
    REAL  par(numoptimpar)    !Current set of parameters
    REAL  optcrit
    REAL, ALLOCATABLE :: performance(:,:)

    !>\b Algorithm \n
    IF(.NOT.ALLOCATED(performance)) ALLOCATE(performance(maxperf,nacrit))

    !> Run a MonteCarlo simulation; nn simulations, generating new random parameter values for each simulation
    DO irun = 1,nruns
      IF(MOD(irun,1000)==0) WRITE(6,*) 'MonteCarlo run no:',irun      !For large amounts of runs, display advancement by printing run counter in hyss.log every 1000 runs
      CALL get_randompar_by_radius(parCenter,parRadStage,parMin,parMax,par)     !Insert random parameters there where they should in "par"
      IF(optim%cal_debugCase == 0) THEN
        CALL run_model_perf(numoptimpar,numoptimpar,par,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,optcrit,performance,status)      ! For the randomized par-vector at hand, run the model and evaluate the criteria function
        IF(status.NE.0) THEN
          WRITE(6,*) 'Error stage MC simulation in run',irun,' of stage',istage
          STOP 1
        ENDIF
      ELSE
        CALL function_to_minim(1.,optcrit,parCenter,par-parCenter,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
      ENDIF

      !> Meanwhile keep track of the num_ens best simulations for next stage
      CALL bookkeep_result_from_simulation(numoptimpar,numoptimpar,optim%nruns_best,par,optcrit,nacrit,performance)
      CALL write_stageMC_calibration_log(irun,icenter,istage)
      IF(writeall) THEN
        iterations_done = (istage-1)*optim%nruns_MC*optim%nruns_best + (icenter-1)*optim%nruns_MC + irun
        CALL write_simulation_results(fileunit_MC,iterations_done,numoptimpar,numoptimpar,maxperf,nacrit,optcrit,performance,par,writematlab)
      ENDIF

    ENDDO

  END SUBROUTINE stage_MonteCarlo_core_function


!>\brief Get random values for parameters to be optimized 
!!
!! Parameter space delimitation achieved by specifying center and
!! radius: parCenter, parRadStage (vectors, all parameter dimensions
!! at once)
!----------------------------------------------------------------------------------
SUBROUTINE get_randompar_by_radius(parCenter, parRadStage, parMin, parMax, par) 

  USE WORLDVAR, ONLY : numoptimpar

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)  :: parCenter(numoptimpar)   !<Parameter space center for current stage
  REAL, INTENT(IN)  :: parRadStage(numoptimpar) !<Array with parameter space radii for current stage
  REAL, INTENT(IN)  :: parMin(numoptimpar)      !<Parameter interval lower limit (reference bounds, do not touch!!)
  REAL, INTENT(IN)  :: parMax(numoptimpar)      !<Parameter interval upper limit (reference bounds, do not touch!!)
  REAL, INTENT(OUT) :: par(numoptimpar)         !<Current parameter values, to be send to model
  
  !Local variables
  INTEGER  :: dimCounter
  INTEGER  :: sumAccept         !Used for check (accept/reject) if randomly generated param set is within absolute boundaries
  REAL     :: rand(numoptimpar) !rand is an array to be filled with random numbers in the interval [0, 1]

  sumAccept = 0                 !Initialize to 0, to enter while-check
  DO WHILE(sumAccept < numoptimpar)
    sumAccept = 0                   !Set to 0 again, in case previous proposed set was rejected and sumAccept is between 0 and npar
    par = 0.                        !Clean the "par" variable (important, might contain unwanted stuff!)

    CALL RANDOM_NUMBER(rand)
    par(:) = parCenter(:) + 2*(rand-0.5)*parRadStage(:)     !Generate parameter set with random numbers

    DO dimCounter = 1,numoptimpar   !Check that generated parameter set is within the abolutely accepted parameter space, i.e. between parMin and parMax
      IF(par(dimCounter)>parMin(dimCounter).AND.par(dimCounter)<parMax(dimCounter)) THEN
        sumAccept = sumAccept + 1   !If param in that dimension is ok, sumAccept is incremented
      ENDIF
    ENDDO
  ENDDO     !End of the while condition for accept/reject param set; test passed only if sumAccept is = npar (max value)

END SUBROUTINE get_randompar_by_radius

!-------------------------------------------------------------------------------------
SUBROUTINE write_stageMC_calibration_log(runCounter, centerCounter, stageCounter)

  USE WORLDVAR, ONLY : optim,                   &
                       modeldir, resdir,        &
                       optimStartTime,          &
                       filename_callog,         &
                       fileunit_callog

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN)  :: runCounter, centerCounter, stageCounter
  
  !Local variables
  INTEGER  iterationsDone, iterationsTotal
  INTEGER  secondsAvgTime, minuteAvgTime, secondsRmngTime, minuteRmngTime, hourRmngTime, daysRmngTime
  REAL nowTime, averageTime, remainingTime

  iterationsTotal = optim%nruns_stages * optim%nruns_MC * optim%nruns_best
  iterationsDone = (stageCounter-1)*optim%nruns_MC*optim%nruns_best + (centerCounter-1)*optim%nruns_MC + runCounter

  ! Time stuff
  CALL cpu_time(nowTime)

  averageTime = (nowTime - optimStartTime)/iterationsDone
  IF(averageTime > 60) THEN
    minuteAvgTime = FLOOR(averageTime/60)
    secondsAvgTime = NINT(averageTime - 60*minuteAvgTime)
  ENDIF

  remainingTime = (iterationsTotal - iterationsDone) * averageTime
  IF(remainingTime > 60) THEN
    daysRmngTime = FLOOR(remainingTime/86400)
    hourRmngTime = FLOOR((remainingTime - 86400*daysRmngTime)/3600)
    minuteRmngTime = FLOOR((remainingTime - 86400*daysRmngTime - 3600*hourRmngTime)/60)
    secondsRmngTime = NINT(remainingTime - 86400*daysRmngTime - 3600*hourRmngTime - 60*minuteRmngTime)
  ENDIF

  ! Start printing log file
  OPEN(UNIT = fileunit_callog, FILE = TRIM(resdir)//filename_callog, STATUS = 'replace')   ! Replace calibration log file
  WRITE(fileunit_callog,'(A)') 'Calibration by stage Monte Carlo simulations'
  WRITE(fileunit_callog,'(A)') '--------------------------------------------'
  WRITE(fileunit_callog,*) ''
  WRITE(fileunit_callog,'(A,A)') 'Model area: ', TRIM(modeldir)
  WRITE(fileunit_callog,*) ''
  WRITE(fileunit_callog,'(I2,A,I6,A,I7,A)') optim%nruns_stages, ' successive zooming stages, each with ', optim%nruns_MC * optim%nruns_best, ' runs  =>  total calibration work: ', iterationsTotal, ' iterations'
  WRITE(fileunit_callog,'(I2,A,I5,A)') optim%nruns_best, ' best simulations used as focus center for next stage  =>  ', optim%nruns_MC, ' runs/center'
  WRITE(fileunit_callog,'(A,F7.3,A,F5.2,A)') 'Zoom factor between two stages is ', optim%nruns_zoom, '  =>  total zoom is ', 100.*optim%nruns_zoom**optim%nruns_stages, '%'   !CP120113

  WRITE(fileunit_callog,*) ''
  WRITE(fileunit_callog,*) ''
  WRITE(fileunit_callog,'(A)') 'Current progress:'
  WRITE(fileunit_callog,'(A)') '-----------------'
  WRITE(fileunit_callog,'(A,I2,A,I2)') 'Stage:  ', stageCounter, '/', optim%nruns_stages
  WRITE(fileunit_callog,'(A,I6,A,I6)') 'Run:    ', (centerCounter-1)*optim%nruns_MC + runCounter, '/', optim%nruns_MC*optim%nruns_best
  WRITE(fileunit_callog,'(A,I7,A,I7)') 'Global: ', iterationsDone, '/', iterationsTotal
  WRITE(fileunit_callog,*) ''

  IF(averageTime > 60) THEN
    WRITE(fileunit_callog, '(A,I2,A,I2)') 'Average iteration time: ', minuteAvgTime, ':', secondsAvgTime
  ELSE
    WRITE(fileunit_callog, '(A,F6.3,A)') 'Average iteration time: ', averageTime, ' seconds'
  ENDIF

  IF(remainingTime > 60) THEN
    IF(daysRmngTime > 0) THEN
      WRITE(fileunit_callog,'(A,I2,A,I2,A,I2,A,I2)') 'Remaining calibration time: ', daysRmngTime, ' days and ', hourRmngTime, ':', minuteRmngTime, ':', secondsRmngTime
    ELSE
      WRITE(fileunit_callog,'(A,I2,A,I2,A,I2)') 'Remaining calibration time: ', hourRmngTime, ':', minuteRmngTime, ':', secondsRmngTime
    ENDIF
  ELSE
    WRITE(fileunit_callog,'(A,F6.3,A)') 'Remaining calibration time: ', remainingTime, ' seconds'
  ENDIF

  WRITE(fileunit_callog,*) ''
  WRITE(fileunit_callog,*) ''
  WRITE(fileunit_callog,'(A)') 'Current best parameter values:'
  WRITE(fileunit_callog,'(A)') '------------------------------'
  WRITE(fileunit_callog,*) ''

  CALL write_stageMC_calibration_log_tail(fileunit_callog)
  CLOSE(fileunit_callog)

END SUBROUTINE write_stageMC_calibration_log

!--------------------------------------------------------
SUBROUTINE write_stageMC_calibration_log_tail(funit)

  USE WORLDVAR, ONLY : bestMCparameters,        &
                       bestMCoptcrit,           &
                       optim,                   &
                       optparid, parindex,      &
                       numoptimpar,             &
                       filename_callog

  USE MODVAR, ONLY : modparid

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: funit  !<file unit of calibration.log
  !Local variables
  INTEGER  columnCounter, lineCounter
  CHARACTER(LEN=20) junkTxt
  CHARACTER(LEN=20*(numoptimpar+1)) fullTxtLine


  !Prepare table top line: parameter labels + criteria
  DO columnCounter = 0,numoptimpar-1
     fullTxtLine(20*columnCounter+ 1:20*columnCounter+ 5) = '     '
     fullTxtLine(20*columnCounter+ 6:20*columnCounter+15) = modparid(optparid(parindex(columnCounter+1,1)))%shortname
     fullTxtLine(20*columnCounter+16:20*columnCounter+20) = '     '
  ENDDO
  fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = '     CRIT           '
  WRITE(funit,'(A)') fullTxtLine

  !Then print all parameter and criteria values
  DO lineCounter = 1,optim%nruns_best
    DO columnCounter = 0,numoptimpar-1
      WRITE(junkTxt,'(F20.9)') bestMCparameters(lineCounter, columnCounter+1)  ! Put into "junkTxt" the string conversion of the parameter values
      fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = junkTxt            ! Concatenate in "fullTxtLine"
    ENDDO
    WRITE(junkTxt,'(F20.9)') bestMCoptcrit(lineCounter)        ! Add the criteria value
    fullTxtLine(20*columnCounter+ 1:20*columnCounter+20) = junkTxt
    WRITE(funit,'(A)') fullTxtLine        ! Output as a string to the calibration log
  ENDDO

  END SUBROUTINE write_stageMC_calibration_log_tail

!----------------------------------------------
!Utilities for methods of line search type
!----------------------------------------------
  
!>This routine initializes the variables common to all non-MonteCarlo methods (variables used in the
!>interruptor routine in particular, see next), then proceeds to the required optimisation method
!-----------------------------------------------
SUBROUTINE linesearch_methods_calibration(npar,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,par)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim,numoptimpar,optimFuncCall,lineSearchCallCount

  IMPLICIT NONE

  !Argument declaration
  INTEGER, INTENT(IN)  :: npar        !<Number of opt parameters/dimension of par-variable
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  REAL,INTENT(OUT) :: par(npar)   !<Optimised values of parameters

  !Variable declarations
  INTEGER lineCounter, rowCounter
  REAL critLastVect(optim%cal_improvCritIter)               !criterion saved for last iterations
  REAL parLastTable(optim%cal_improvParamIter, numoptimpar) !parameters saved for last iterations

  !Initialize function call and line search counters
  optimFuncCall = 0
  lineSearchCallCount = 0

  !Initialize critLastVect and parLastTable
  !Done exotically so as to check that content is handled correctly by interruptor (forget one, shift the others, insert last one; see the iterruptor routine)
  DO lineCounter = 1,optim%cal_improvCritIter
    critLastVect(lineCounter) = lineCounter
  ENDDO

  DO lineCounter = 1,optim%cal_improvParamIter
    DO rowCounter = 1,numoptimpar
      parLastTable(lineCounter, rowCounter) = lineCounter
    ENDDO
  ENDDO

  WRITE(6,*) ''

  IF(optim%task_stpstDesc .OR. optim%task_DFP .OR. optim%task_BFGS) THEN  !Start the quasi-Newton gradient-based DFP or BFGS optimisation procedure, or steepest descent (implemented as particular case of quasi-Newton)  [Fred, 15.07.10]
    IF(optim%task_stpstDesc) THEN
      WRITE(6,'(A)') 'Calibration with steepest descent optimisation method'
    ELSEIF(optim%task_DFP) THEN
      WRITE(6,'(A)') 'Calibration with DFP quasi-Newton optimisation method'
    ELSEIF(optim%task_BFGS) THEN
      WRITE(6,'(A)') 'Calibration with BFGS quasi-Newton optimisation method'
    ENDIF
    CALL quasiNewton_algorithm(critLastVect,parLastTable,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  ENDIF

  IF(optim%task_BrentNew) THEN                  !Start optimisation by Brent routine, with new line search method  [Fred, 28.11.10]
    WRITE(6,'(A)') 'Calibration by Brent method'
    CALL new_Brent_method(critLastVect,parLastTable,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  ENDIF

  !Set optimised parameter values to subroutine argument
  par(:) = parLastTable(1,:)

END SUBROUTINE linesearch_methods_calibration

!>\brief Establish statistics on the latest calibration iterations to
!!check for criteria and parameter improvements
!!It also aborts calibration if method timed-out, or max amount of
!!calibration iterations is reached
!------------------------------------------------------------------
SUBROUTINE linesearch_methods_interruptor(iterCounter, critLast, &
                 critLastVect, parLast, parLastTable, parPrecision, finished)

  USE WORLDVAR, ONLY : optim, numoptimpar, optimStartTime

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN)     :: iterCounter
  REAL, INTENT(IN)        :: critLast
  REAL, INTENT(INOUT)     :: critLastVect(optim%cal_improvCritIter)
  REAL, INTENT(IN)        :: parLast(numoptimpar)
  REAL, INTENT(INOUT)     :: parLastTable(optim%cal_improvParamIter, numoptimpar)
  REAL, INTENT(IN)        :: parPrecision(numoptimpar)
  LOGICAL, INTENT(OUT)    :: finished                 !Flag for aborting calibration
  
  !Local variables
  INTEGER lineCounter, rowCounter, precisCounter
  REAL critLastMin, critLastMax, critLastMean, critLastDelta, ratio
  REAL nowTime, ellapsedTime

  !Initialise
  finished = .FALSE.

  !Shift list of parameter sets by one and insert current parameter set on line 1
  DO lineCounter = 1,optim%cal_improvParamIter-1
    DO rowCounter = 1,numoptimpar
      parLastTable(optim%cal_improvParamIter+1-lineCounter, rowCounter) = parLastTable(optim%cal_improvParamIter-lineCounter, rowCounter)
    ENDDO
  ENDDO
  DO rowCounter = 1,numoptimpar
    parLastTable(1,rowCounter) = parLast(rowCounter)
  ENDDO

  !Check for method time out
  CALL cpu_time(nowTime)
  ellapsedTime = nowTime - optimStartTime
  IF(ellapsedTime > 3600*optim%cal_maxTime)THEN
    CALL linesearch_methods_interruptor_printandstop(1,finished)
    RETURN
  ENDIF

  !Check for max allowed amount of iterations
  IF(iterCounter >= optim%cal_maxIterat)THEN
    CALL linesearch_methods_interruptor_printandstop(2,finished)
    RETURN
  ENDIF

  !Check criteria improvement statistics
  !Shift best criterium list of previous iterations by one, and insert latest result at top
  DO lineCounter = 1,optim%cal_improvCritIter-1
    critLastVect(optim%cal_improvCritIter+1-lineCounter) = critLastVect(optim%cal_improvCritIter-lineCounter)
  ENDDO
  critLastVect(1) = critLast

  !Check if criteria changed significantly over the last iterations
  IF(iterCounter >= optim%cal_improvCritIter) THEN
    critLastMin = MINVAL(critLastVect)
    critLastMax = MAXVAL(critLastVect)
    critLastMean = ABS(critLastMin + critLastMax)/2
    critLastDelta = ABS(critLastMin - critLastMax)
    ratio = critLastDelta/critLastMean
    IF(ratio < optim%cal_improvCritTol)THEN
      CALL linesearch_methods_interruptor_printandstop(3,finished)
      RETURN
    ENDIF
  ENDIF

  !Check parameter improvement statistics
  precisCounter = 0

  !Check of parameter values have changed less than required precision over the last iterations
  IF(iterCounter >= optim%cal_improvParamIter) THEN
    DO rowCounter = 1,numoptimpar
      critLastMin = MINVAL(parLastTable(:,rowCounter))
      critLastMax = MAXVAL(parLastTable(:,rowCounter))
      critLastDelta = ABS(critLastMin - critLastMax)
      IF(critLastDelta < parPrecision(rowCounter)/2) precisCounter = precisCounter + 1
    ENDDO
    IF(precisCounter == numoptimpar)THEN
      CALL linesearch_methods_interruptor_printandstop(4,finished)
      RETURN
    ENDIF
  ENDIF

END SUBROUTINE linesearch_methods_interruptor

!>Print the total calibration time and the amount of function calls
!>before ending calibration
!-----------------------------------------------------------------------------
SUBROUTINE linesearch_methods_interruptor_printandstop(stopFlag,finished)

  USE WORLDVAR, ONLY : optim, fileunit_callog

  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: stopFlag
  LOGICAL, INTENT(OUT) :: finished

  CALL linesearch_methods_interruptor_printandstop_core_function(stopFlag, 6) !Print in hyss.log
  IF(optim%cal_log) CALL linesearch_methods_interruptor_printandstop_core_function(stopFlag, fileunit_callog) !Print in calibration.log if requested
  finished = .TRUE.

END SUBROUTINE linesearch_methods_interruptor_printandstop

!>Print the total calibration time and the amount of function calls before stopping
!-----------------------------------------------------------------------------------------
SUBROUTINE linesearch_methods_interruptor_printandstop_core_function(stopFlag,fileunit)

  USE WORLDVAR, ONLY : optim, optimStartTime, optimFuncCall

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: stopFlag, fileunit
  
  !Local variables
  INTEGER ellapsedDays, ellapsedHours, ellapsedMinutes, ellapsedSeconds
  REAL nowTime, ellapsedTime

  !First, print out the reason for the stop
  WRITE(fileunit,*) ''

  SELECT CASE(stopFlag)
  CASE(1)   !Case 1: calibration method timed out
    WRITE(fileunit,'(A)') 'Calibration stopped: method timed out'
  CASE(2)   !Case 2: max allowed amount of calibration iterations reached
    WRITE(fileunit,'(A)') 'Calibration stopped: maximum amount of iterations allowed reached'
  CASE(3)   !Case 3: criteria improved less than relative tolerance over few last iterations
    WRITE(fileunit,'(A,I2,A)') 'Calibration stopped: criteria changed less than specified tolerance over the last ', optim%cal_improvCritIter, ' iterations'
  CASE(4)   !Case 4: all parameter values improved less than given decimals over few last iterations
    WRITE(fileunit,'(A,I2,A)') 'Calibration stopped: all parameter values changed less than specified tolerances over the last ', optim%cal_improvParamIter, ' iterations'
  CASE(5)   !Case 5, quasi-Newton specific: gradient norm reached tolerance
    WRITE(fileunit,'(A)') 'Calibration stopped: gradient norm smaller than specified tolerance'
  CASE(6)   !Case 6, quasi-Newton specific: deltaParam and deltaGrad are perpendicular => optimum reached
    WRITE(fileunit,'(A)') 'Optimum reached, calibration terminated'
  CASE(7)   !Case 7, Brent specific: direction of diagonal step zero => optimum reached
    WRITE(fileunit,'(A)') 'Optimum reached, calibration terminated'
  END SELECT

  ! Second, print the total calibration running time
  CALL CPU_TIME(nowTime)
  ellapsedTime = nowTime - optimStartTime

  ellapsedDays = FLOOR(ellapsedTime/86400)
  ellapsedHours = FLOOR((ellapsedTime - 86400*ellapsedDays)/3600)
  ellapsedMinutes = FLOOR((ellapsedTime - 86400*ellapsedDays - 3600*ellapsedHours)/60)
  ellapsedSeconds = NINT(ellapsedTime - 86400*ellapsedDays - 3600*ellapsedHours - 60*ellapsedMinutes)

  IF(ellapsedDays > 0) THEN
    WRITE(fileunit,'(A,I2,A,I2,A,I2,A,I2)') 'Calibration total running time: ', ellapsedDays, ' days and ', ellapsedHours, ':', ellapsedMinutes, ':', ellapsedSeconds
  ELSE
    WRITE(fileunit,'(A,I2,A,I2,A,I2)') 'Calibration total running time: ', ellapsedHours, ':', ellapsedMinutes, ':', ellapsedSeconds
  ENDIF

  ! Third, print the total amount of function calls
  WRITE(fileunit,'(A,I6)') 'Total amount of function calls: ', optimFuncCall

END SUBROUTINE linesearch_methods_interruptor_printandstop_core_function

!>\brief Line search for minimum 
!-------------------------------------------------------------------------------------------
SUBROUTINE linesearch_HYSS(xMinIn,xMaxIn,refPar,parPrecis,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,xBest,fBest)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim,numoptimpar,fileunit_callog,lineSearchCallCount

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)    :: xMinIn           !<point search interval minimum
  REAL, INTENT(IN)    :: xMaxIn           !<point search interval maximum
  REAL, INTENT(IN)    :: refPar(numoptimpar)      !<parameter values at origin, where x=0
  REAL, INTENT(IN)    :: parPrecis(numoptimpar)   !<precision
  REAL, INTENT(IN)    :: direction(numoptimpar)   !<direction of line search
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  REAL, INTENT(OUT)   :: xBest    !<best point found
  REAL, INTENT(OUT)   :: fBest    !<function value at best point
  
  !Local variables
  INTEGER iterCounter, goldenFlag, maxIter
  REAL a, b, gRatio, xMean, xCurrent, fCurrent, xBefore, fBefore, xBefore2, fBefore2, xNew, fNew
  REAL Term1, Term2, slope, curvature, shiftDist, shiftDir, maxDist, maxDistOld
  REAL eps, tolerance, pointTol, intervalTol
  LOGICAL stopFlag
  CHARACTER(LEN=9) procStrg

  lineSearchCallCount = lineSearchCallCount + 1

  !>\b Algorithm \n
  !>Initialize the numerical parameters
  a = xMinIn
  b = xMaxIn
  maxIter = optim%lineSearch_maxIter
  tolerance = optim%lineSearch_tol
  stopFlag = .FALSE.      !Flag used to stop line search (stopFlag = .TRUE.) in case the required decimal precision is reached

  !>Compute the starting "golden ratio" point, its related function value, and start the output log
  iterCounter = 1                   ! First iteration (initial)
  procStrg = ' initial '

  gRatio = (3. - SQRT(5.))/2.       ! This is = phi_minus + 1 = 0.382, where phi_minus = (1 - SQRT(5))/2 is the negative solution of the golden ratio equation  phi^2 - phi - 1 = 0
  xCurrent = a + gRatio*(b - a)     ! Set the starting x-value at 38.2% within the [a, b] interval (golden ratio like)
  CALL function_to_minim(xCurrent,fCurrent,refPar,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)   ! Get corresponding function value f(xCurrent)

  IF(optim%cal_log) THEN
    WRITE(fileunit_callog, '(A)') 'Calls count         a               b               x             f(x)         Procedure'
    WRITE(fileunit_callog,'(I12, F16.12, F16.12, F16.12, F16.12, A12)') iterCounter, a, b, xCurrent, fCurrent, procStrg
  ENDIF

  !>Define tolerances for point and parabola accept/reject tests
  ! NOTE: tolerances are evolving! They are defined according to xCurrent, which is iteratively updated in the while-loop below, and so are
  !       the tolerances. Since intervalTol is used to enter the loop, there's nothing else to do then write the tolerance definitions twice :(
  eps = 2.2204E-16                    ! Distance from 1.0 to next largest double-precision number, i.e. eps = 1/2^52
  pointTol = SQRT(eps)*ABS(xCurrent) + tolerance/3
  intervalTol = 2*pointTol

  !Initializations/preparations to enter while-loop
  xMean = (a + b)/2           ! Usual mean x-value within the current bracketing interval [a, b] (= interval center value)
  xBefore = xCurrent          ! 2 previous current points and corresp. function values are retained for parabolic approx. Initialize as xCurrent and fCurrent, it's good enough
  fBefore = fCurrent
  xBefore2 = xCurrent
  fBefore2 = fCurrent
  shiftDist = 0               ! xNext = xCurrent + shiftDist. This is the quantity most of the while-loop is all about
  maxDist = 0                 ! maxDist = max{|xCurrent, a|, |xCurrent, b|} = largest distance between xCurrent and the interval boundary. Initialized to 0 here, to avoid entering parabolic fitting (which is useless at 1st
  ! iteration); is correctly updated later, of course. Alternatively, entering fit could also be avoided by initializing pointTol in a funny way; proceeding via maxDist is more natural, though

  !>Here starts the iterative procedure
  ! The whole point is to determine a value for shiftDist (positive or negative), with either the parabolic or golden step method
  ! Then shift the current x-value by that amount, giving xNext, evaluate the function at that point, and accordingly modify the interval boundaries (shrink)
  DO WHILE((intervalTol - (b - a)/2) < ABS(xCurrent - xMean) .AND. iterCounter < maxIter .AND. (.NOT.stopFlag))
     ! Test is subtle. Observe primo that LHS is < 0 unless half interval size |a, b|/2 < intervalTol, which is typically 1e-07 => takes in principle a few iterations to contract [a, b] down to that level; in particular, entering loop is garanteed at first attempt (unless b < a, but this has been tested for)
     ! Secundo, once [a, b] is contracted enough for LHS to be positive at all, then distance |xCurrent, xMean| has to be even smaller. This is however quickly achieved => test is essentially on contraction of [a, b], rather than on position of xCurrent relative to interval center xMean

     goldenFlag = 1                        ! Default assumption is that parabola fit will fail => set golden step flag so as to take golden step; will be reseted if parabola fit is possible and accepted

     !>\li If distance is large enough, do a parabola fit
     IF(ABS(maxDist) > pointTol) THEN      ! Test if |xCurrent, xBoundary| is large enough for a parabola fit. Will fail at 1st and 2nd attempt, and for distances below machine precision

        ! Parabolic fit is of the form  p(x) = a0 + a1*x + a2*x^2. Calculate parabola parameters {a0,a1,a2} from 3 points: xCurrent, xBefore and xBefore2, with p(xi) = fi, i = Current,Before,Before2
        ! Zero of parabola derivative (for min or max!) is at x0 = - a1/(2*a2)  => just need to determine parabola parameters a1 and a2 to get minimum, don't care about a0
        ! a1 and a2 determined by Gaussian elimination; start by eliminating a0, of course, rest is easy and even enjoyable if done elegantly :)
        ! a1 = [(x1 + x3)T1 - (x1 + x2)T2]/D  ,  a2 = (T1 - T2)/(-D)    where    T1 = (x1 - x3)(f1 - f2)  ,  T2 = (x1 - x2)(f1 - f3)  and  D = (x1 - x2)(x1 - x3)(x3 - x2)
        ! D eliminates from the ratio -a1/(2*a2) and is not necessary to compute:  x0 = - a1/(2*a2) = [(x1 + x3)T1 - (x1 + x2)T2] / 2(T1 - T2)
        ! Further, we are not after x0 itself but the shift d such that  x0 = x1 + d  =>  d = - a1/(2*a2) - x1 = ... = [(x1 - x2)T2 - (x1 - x3)T1]/2(T1 - T2)
        ! "slope" and "curvature" here are NOT the actual a1 and a2 parabola parameters, but the numerator and denominator of above equation. They are related to the parabola slope and curvature, though
        ! [for more details, ask Frédéric]
        Term1 = (xCurrent - xBefore2) * (fCurrent - fBefore)                    ! Take x1 = xCurrent, x2 = xBefore, x3 = xBefore2  and corresp. function evaluations for fi
        Term2 = (xCurrent - xBefore) * (fCurrent - fBefore2)
        slope = (xCurrent - xBefore) * Term2 - (xCurrent - xBefore2) * Term1    ! numerator in formula for d above, related to a1, "slope" parabola parameter
        curvature = 2*(Term1 - Term2)                                           ! denominator in formula for d, related to a2, "curvature" parabola parameter

        IF(curvature < 0) slope = -slope    ! If x0 is maximum (= curvature < 0), go in the opposite direction => take opposite sign for numerator
        curvature = ABS(curvature)          ! Curvature set positive, i.e. as minimum parabola in all cases

        maxDistOld = maxDist                ! This is just a dummy for the test bellow
        maxDist = shiftDist                 ! Update, for next iteration

        IF(ABS(slope) < ABS(curvature*maxDistOld/2) .AND. slope > curvature*(a - xCurrent) .AND. slope < curvature*(b - xCurrent)) THEN     ! Test if parabola predicts a point within an acceptable interval. Observe that tests are actually on the shift, d = slope/curvature
           shiftDist = slope/curvature
           xNew = xCurrent + shiftDist       ! New point is "officially" calculated later; this is just a dummy for the test below
           procStrg = 'parabolic'            ! Register for log that step was done with parabolic approx

           IF((xNew - a) < intervalTol .OR. (b - xNew) < intervalTol) THEN     ! Test to avoid taking point too close to interval boundaries. NOTA: this procedure slows down the algorithm treamendously when the minimum is at a boundary point :(
              shiftDir = SIGN(1., xMean - xCurrent)
              shiftDist = shiftDir*pointTol   ! If xNew is too close to any interval boundary, then take new distShift at safe distance
           ENDIF

           goldenFlag = 0                    ! Reset golden step flag, to avoid having parabola approx of shiftDist overruled by golden step method
        ENDIF

     ENDIF     ! End of the "IF" about doing a parabola fit or not. At this point, "goldenFlag" = 0 only if a parabola fit was possible AND acceptable

     !>If no new point was found by parabola fit, find one with golden step method instead
     IF(goldenFlag == 1) THEN      ! Now, if the golden step flag is set, overule the new point determination by parabolic fit with golden step determination of new point
        IF(xCurrent > xMean) THEN
           maxDist = a - xCurrent            ! xCurrent is on the right of xMean  =>  maxDist < 0  =>  shiftDist < 0  =>  new point will be taken on the left
        ELSE
           maxDist = b - xCurrent            ! xCurrent is on the left of xMean  =>  new point to be taken on the right
        ENDIF

        shiftDist = gRatio * maxDist        ! Shrink interval
        procStrg = '  golden '              ! Register for log that step was golden step approx
     ENDIF

     !>Simulate the model for the new point to calculte the optimization criterion
     ! Now that shiftDist has been calculated by either one or the other method, shift
     ! the current x-point by that amount and get the corresponding function value
     shiftDir = SIGN(1., shiftDist)                                ! shiftDir is direction of shift (depends on sign of shiftDist, of course)
     xNew = xCurrent + shiftDir * MAX(ABS(shiftDist), pointTol)    ! This is the new x...
     CALL function_to_minim(xNew,fNew,refPar,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)         ! ...and the corresponding function value
     iterCounter = iterCounter + 1

     !>Update optimal points (theree are saved)
     IF(fNew <= fCurrent) THEN     ! This is the ideal case: the new point was an improvement :)
        IF(xNew >= xCurrent) THEN
           a = xCurrent                  ! If the improved new point is left from the current one, then set the right boundary as the current point
        ELSE
           b = xCurrent                  ! Otherwise the current point becomes the left boundary
        ENDIF
        xBefore2 = xBefore          ! New point + function value becomes current point. Current becomes previous, and previous points are shifted back in memory
        fBefore2 = fBefore
        xBefore = xCurrent
        fBefore = fCurrent
        xCurrent = xNew
        fCurrent = fNew

     ELSE                          ! Tricky case where fNew > fCurrent: the new point was NOT an improvement
        IF(xNew < xCurrent) THEN
           a = xNew                      ! If xNew is left of xCurrent, then xNew becomes the new left boundary...
        ELSE
           b = xNew                      ! ...and vice versa
        ENDIF

        ! Now, test fNew with respect to the older points function values to find whether to retain xNew, and if so then in which "position"
        IF(fNew <= fBefore .OR. xBefore == xCurrent) THEN
           xBefore2 = xBefore
           fBefore2 = fBefore
           xBefore = xNew
           fBefore = fNew
        ELSEIF(fNew <= fBefore2 .OR. xBefore2 == xCurrent .OR. xBefore2 == xBefore) THEN
           xBefore2 = xNew
           fBefore2 = fNew
        ENDIF

     ENDIF   ! End of the point sorting part

     !>Print iteration results in log
     IF(optim%cal_log) WRITE(fileunit_callog,'(I12, F16.12, F16.12, F16.12, F16.12, A12)') iterCounter, a, b, xNew, fNew, procStrg

     xMean = (a + b)/2     ! Now that a and b are updated, update xMean as well

     !>Update tolerances in accordance to new best point
     pointTol = SQRT(eps)*ABS(xCurrent) + tolerance/3   ! (same as orginal tolerance definitions before the while-loop)
     intervalTol = 2*pointTol

     !>Check if parameter precision is reached
     CALL linesearch_check_decim(ABS(a - b), parPrecis, direction, stopFlag)

  ENDDO   ! End of the WHILE loop

  IF(optim%cal_log) THEN
     WRITE(fileunit_callog,*) ''
     WRITE(fileunit_callog,*) ''
  ENDIF

  xBest = xCurrent
  fBest = fCurrent

END SUBROUTINE linesearch_HYSS

!------------------------------------------------------------------------------
SUBROUTINE linesearch_check_decim(intLength, parPrecis, direction, stopFlag)

  USE WORLDVAR, ONLY : optim, numoptimpar, fileunit_callog

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)     :: intLength   !<Current line search interval length
  REAL, INTENT(IN)     :: parPrecis(numoptimpar), direction(numoptimpar)
  LOGICAL, INTENT(OUT) :: stopFlag
  
  !Local variables
  INTEGER dimCounter, summatorio
  REAL    stuff

  stopFlag = .FALSE.
  summatorio = 0

  DO dimCounter = 1, numoptimpar
    stuff = parPrecis(dimCounter)/(2*direction(dimCounter))
    IF(intLength < stuff) summatorio = summatorio + 1
  ENDDO

  IF(summatorio == numoptimpar) THEN
    stopFlag = .TRUE.
    IF(optim%cal_log) WRITE(fileunit_callog,'(A)') 'Line search exited: length of search interval smaller than required decimal precision'
  ENDIF

END SUBROUTINE linesearch_check_decim

!> Choose objective function, HYSS model or one of the simple debug functions.
!-------------------------------------------------------------------
SUBROUTINE function_to_minim(lambda,flambda,refPar,direction,   &
                             frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim,numoptimpar,optimFuncCall
  USE MODVAR, ONLY   : pi

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)    :: lambda       !<length of current step
  REAL, INTENT(OUT)   :: flambda      !<optimization criteria value
  REAL, INTENT(IN)    :: refPar(numoptimpar)    !<reference parameter values
  REAL, INTENT(IN)    :: direction(numoptimpar) !<current direction of step
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake state
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states

  !Local variables
  INTEGER status
  REAL  optcrit, x, y, z, sigma, epsilon, term1, term2, term3
  REAL  newPar(numoptimpar)

  optimFuncCall = optimFuncCall + 1       ! Increment function call counter
  newPar = refPar + lambda*direction      ! Build new parameter vector

  SELECT CASE(optim%cal_debugCase)
  CASE(0)       !Not a debug: run HYPE, return criteria
    CALL run_model_crit(numoptimpar,numoptimpar,newPar, &
                        frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate, &
                        optcrit,status)   !Crit=999
    IF(status.NE.0)THEN
      WRITE(6,*) 'Error in optimization simulation; run_model_crit'
      STOP 1
    ENDIF
    flambda = optcrit
  CASE(1)       !3-param polynomial debug: f(x,y,z) = (x-7)^2 + (y-\sqrt(2))^4 + (z-\pi)^6, has zero and minimum at x = 7, y = 1.4142136, z = 3.1415927
    x = newpar(1)
    y = newpar(2)
    z = newpar(3)
    flambda = (x - 7)**2 + (y - SQRT(2.0))**4 + (z - pi)**6
  CASE(2)           !A classic: purely trigonometrical 2-parameters function; very nice when plotted
    x = newpar(1)   !Recommanded search space: 0.8 =< x,y =< 2
    y = newpar(2)   !Minimum at approx x = 1.3764, y = 1.6787
    flambda = COS(x**2 - 3*y) + SIN(x**2 + y**2)
  CASE(3)           !Lennard-Jones/Buckingham potential-like function, for 2 parameters: f(x,y) = epsilon * sin(y) * [exp(-x/sigma) + (sigma/x)^8  - gamma*(sigma/x)^6]
    x = newpar(1)   !Recommanded search:  1 =< x =< 2.5  ,  0 =< y =< 2*pi  !! CAUTION: saddle point at x = 1.117, y = pi ! Function diverges strongly for x -> 0, toward + infinity for 0 < y < pi (wall), and toward - infinity for pi < y < 2*pi (well)
    y = newpar(2)   !The "infinity wall/well" are however guarded by a local minimum (the target!), resp. a small hill
    epsilon = 1     !It is therefore recommanded to take starting x-values 1.117 < x0, or to combine x0 < 1.117, y0 < pi
    sigma = pi**2/5 !Minimum (local) is at x = 1.2868407, y = pi/2
    flambda = epsilon * SIN(y) * (EXP(-x/sigma) + (sigma/x)**8 - pi * (sigma/x)**6)
  CASE(4)           !An awsome 2-parameters, featuring several max and min, saddle points, and flatness at infinity :)
    x = newpar(1)   !Minimum at  x = 0.228879, y = -1.626176
    y = newpar(2)   !Recommanded search interval:  -3 =< x =< 3  ,  -2.5 =< y =< 2.5
    term1 = 3*(1 - x)**2 * EXP(-(x**2) - (y + 1)**2)
    term2 = - 10 * (x/5 - x**3 - y**5) * EXP(-x**2 - y**2);
    term3 = - 1/3 * EXP(-(x+1)**2 - y**2);
    flambda = term1 + term2 + term3
  END SELECT

END SUBROUTINE function_to_minim

!------------------------------------
!Block of Brent autocal routines
!------------------------------------

!> Brent algorithm for dimension-wise optimisation
!----------------------------------------------------------
SUBROUTINE new_Brent_method(critLastVect,parLastTable,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim, numoptimpar, optparid, parindex, fileunit_callog
  USE MODVAR, ONLY   : modparid

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(INOUT) :: critLastVect(optim%cal_improvCritIter) !<criterion saved for last iterations
  REAL, INTENT(INOUT) :: parLastTable(optim%cal_improvParamIter, numoptimpar) !<parameters saved for last iterations
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  
  !Local variables
  INTEGER   npar, dimCounter, iteratCounter
  INTEGER   stopFlag
  LOGICAL   finished                                           ! Flag for calibration finished
  REAL      stuff, critBest, lambdaBest, lambdaNeg, lambdaPos
  REAL      parMin(numoptimpar), parMax(numoptimpar)           ! Lower/upper interval limits, from optpar.txt; those limits are absolute, serve as permanent/reference bounds!!
  REAL      parPrecis(numoptimpar)
  REAL      par(numoptimpar), parStep(numoptimpar)             ! Current best parameter set, modified parameter set to send to line search
  REAL      newPar(numoptimpar)                                ! New parameter set after Brent step
  REAL      iterStartPar(numoptimpar)                          ! Parameter set at Brent iteration start
  REAL      direction(numoptimpar)                             ! Vectorial direction between current and new best vector
  REAL      lambdaMin(numoptimpar), lambdaMax(numoptimpar)
  REAL      parBestBefore(numoptimpar),critBestBefore          ! Parameters and Criterion for best point before each linesearch, kept if no better found  !CP120524

  !>\b Algorithm \n
  !Check numerical parameters validity
  optim%QN_lambdaMaxFac = 1

  !>Initiate Brent simulation; find parameters and optimization intervals, and read the starting set
  iteratCounter = 0
  finished = .FALSE.
  CALL find_optpar(par,parMin,parMax,parPrecis,npar)  !Get optimization interval limits from optpar.txt; par is the "center"; npar = numoptimpar (useless)
  CALL load_qNstarting_vector(parMin,parMax,par)      !Read qNstartpar.txt

  !Write the calibration log header in hyss.log    
  WRITE(6,'(A9)', advance='no') 'Iteration'
  DO dimCounter = 1, numoptimpar
    WRITE(6, '(A16)', advance='no') modparid(optparid(parindex(dimCounter,1)))%shortname
  ENDDO
  WRITE(6,'(A)') '   Criteria        Nada'

  !>Run simulation for the starting set of parameters
  lambdaBest = 0
  direction = 0.
  CALL function_to_minim(lambdaBest,critBest,par,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

  IF(optim%cal_log) CALL initialize_Brent_calibrationlog()
  CALL print_calib_HYSSlog(iteratCounter, par, critBest, 0.)

  !>Repeatedly perform Brent steps, as long as max criteria improves more then tolerance between 2 iterations,
  !>and as long as max amount of iterations is not reached
  DO
    iteratCounter = iteratCounter + 1
    iterStartPar = par      ! Copy parameter set at iteration start, for ev. diagonal step at iteration end

    IF(optim%cal_log) THEN
      WRITE(fileunit_callog, *) ''
      WRITE(fileunit_callog, '(A, I4)')    'Iteration ', iteratCounter
      WRITE(fileunit_callog, '(A)') '--------------'
      WRITE(fileunit_callog, *) ''
      WRITE(fileunit_callog, *) ''
    ENDIF

    !>\li Perform line search, parameter-wise
    DO dimCounter = 1,numoptimpar
      IF(optim%cal_log) THEN
        WRITE(fileunit_callog, '(A, A)') 'Parameter ', modparid(optparid(parindex(dimCounter,1)))%shortname
        WRITE(fileunit_callog, *) ''
      ENDIF

      parBestBefore = par               !save for check if linesearch gave improvement
      critBestBefore = critBest          

      parStep = par
      parStep(dimCounter) = parMin(dimCounter)

      direction = 0.
      direction(dimCounter) = parMax(dimCounter) - parMin(dimCounter)

      CALL linesearch_HYSS(0.0,1.0,parStep,parPrecis,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,lambdaBest,critBest)

      newPar = parStep + lambdaBest*direction
      IF(optim%cal_log) CALL write_Brent_calibrationlog(par, parStep, direction, newPar, lambdaBest, critBest)

      !Check if linesearch gave improvement
      IF(critBest<critBestBefore)THEN
        par = newPar
        CALL print_calib_HYSSlog(iteratCounter, par, critBest, 0.)
      ELSE
        par = parBestBefore
        critBest = critBestBefore   !keep old best
        CALL print_calib_HYSSlog_noimp(iteratCounter, par, critBest)
      ENDIF
    ENDDO

    !>\li Perform diagonal step at Brent iteration end if required
    IF(optim%Brent_diagonalStep)THEN
      IF(optim%cal_log) THEN
        WRITE(fileunit_callog, '(A)') 'Diagonal step'
        WRITE(fileunit_callog, *) ''
      ENDIF

      parBestBefore = par               !save for check if linesearch gave improvement
      critBestBefore = critBest          

      direction = par - iterStartPar
      IF(SUM(direction)==0)THEN    !check if optimum reached (no changes) !CP120918
        IF(optim%cal_log) WRITE(fileunit_callog,*) 'diagonal step direction 0, skipping diagonal step'
        CALL print_calib_HYSSlog_noimp(iteratCounter, par, critBest)
        stopFlag = 7
        CALL linesearch_methods_interruptor_printandstop_core_function(stopFlag, 6)                               !Print in hyss.log
        IF(optim%cal_log) CALL linesearch_methods_interruptor_printandstop_core_function(stopFlag, fileunit_callog)    !Print in calibration.log if requested
        EXIT
      ELSE

        DO dimCounter = 1,numoptimpar
          lambdaMin(dimCounter) = (parMin(dimCounter) - iterStartPar(dimCounter)) / direction(dimCounter)
          lambdaMax(dimCounter) = (parMax(dimCounter) - iterStartPar(dimCounter)) / direction(dimCounter)

          IF(lambdaMin(dimCounter) > lambdaMax(dimCounter))THEN
            stuff = lambdaMin(dimCounter)
            lambdaMin(dimCounter) = lambdaMax(dimCounter)
            lambdaMax(dimCounter) = stuff
          ENDIF
          IF(lambdaMin(dimCounter) > 0)THEN
            WRITE(6,*) 'ERROR in Brent optimisation: found parameter outside parameter range'
            WRITE(6,*) 'opt parameter no: ',dimCounter, 'lambdaMin: ',lambdaMin(dimCounter)
            STOP 1
          ENDIF
        ENDDO
        lambdaNeg = MINVAL(ABS(lambdaMin))
        lambdaPos = MINVAL(ABS(lambdaMax))

        CALL linesearch_HYSS(-lambdaNeg,lambdaPos,iterStartPar,parPrecis,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,lambdaBest,critBest)
        newPar = iterStartPar + lambdaBest*direction

        IF(optim%cal_log) CALL write_Brent_calibrationlog(par, iterStartPar, direction, newPar, lambdaBest, critBest)

        !Check if linesearch gave improvement
        IF(critBest<critBestBefore)THEN
          par = newPar
          iterStartPar = newPar
          CALL print_calib_HYSSlog(iteratCounter, par, critBest, 0.)
        ELSE
          par = parBestBefore
          critBest = critBestBefore   !keep old best
          CALL print_calib_HYSSlog_noimp(iteratCounter, par, critBest)
        ENDIF
      ENDIF
    ENDIF   !End of diagonal step at Brent iteration end

    !>\li Check if iterations are enough and to be stopped
    CALL linesearch_methods_interruptor(iteratCounter, critBest, critLastVect, par, parLastTable, parPrecis, finished)
    IF(finished) EXIT   

  ENDDO   !End of Brent iteration loop

  CLOSE(fileunit_callog)

END SUBROUTINE new_Brent_method

!------------------------------------------------
SUBROUTINE initialize_Brent_calibrationlog()

  USE WORLDVAR, ONLY : optim,           &
                       modeldir,        &
                       resdir,          & 
                       filename_callog, & 
                       fileunit_callog

  IMPLICIT NONE

  ! Start printing log file
  OPEN(UNIT = fileunit_callog, FILE = TRIM(resdir)//filename_callog, STATUS = 'replace')   ! Start calibration log file

  WRITE(fileunit_callog, '(A)') 'Calibration by Brent algortihm'
  WRITE(fileunit_callog, '(A)') '------------------------------'

  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, '(A, A)') 'Model area: ', trim(modeldir)
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, '(A)') 'Numerical parameters:'
  WRITE(fileunit_callog, '(A, I4)')    'Maximum amount of allowed Brent iterations                        : ', optim%cal_maxIterat
  WRITE(fileunit_callog, '(A, I4)')    'Amount of iterations taken into account to assess method progress : ', optim%cal_improvCritIter
  IF(optim%Brent_diagonalStep)THEN
     WRITE(fileunit_callog, '(A)')      'Diagonal step at the end of each Brent iteration                  : Yes'
  ELSE
     WRITE(fileunit_callog, '(A)')      'Diagonal step at the end of each Brent iteration                  : No'
  ENDIF
  WRITE(fileunit_callog, '(A, F7.3)')  'Squeeze to contain lambda factor (must be 1.00 for Brent!)        : ', optim%QN_lambdaMaxFac
  WRITE(fileunit_callog, '(A, F13.9)') 'Tolerance for line search                                         : ', optim%lineSearch_tol
  WRITE(fileunit_callog, '(A, I4)')    'Max amount of iterations allowed in line search                   : ', optim%lineSearch_maxIter

  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''

END SUBROUTINE initialize_Brent_calibrationlog

!------------------------------------------------------------------------------------------------
SUBROUTINE write_Brent_calibrationlog(par, parStep, direction, newPar, lambdaBest, critBest)

  USE WORLDVAR, ONLY : numoptimpar, fileunit_callog

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN) :: par(numoptimpar), parStep(numoptimpar), direction(numoptimpar), newPar(numoptimpar)
  REAL, INTENT(IN) :: lambdaBest, critBest
  
  !Local variables
  INTEGER counter

  WRITE(fileunit_callog, '(A)') 'Start parameter set:       Modified parameter set:           Direction vector:           New parameter set:'
  DO counter = 1, numoptimpar
     WRITE(fileunit_callog, '(F17.9, A, F17.9, A, F17.9, A, F17.9)', advance='no') par(counter), '            ', parStep(counter), '            ', direction(counter), '            ', newPar(counter)
     WRITE(fileunit_callog, '(A)') ' '
  ENDDO
  WRITE(fileunit_callog, *) ''

  WRITE(fileunit_callog, '(A,F17.9)') 'Optimal lambda:', lambdaBest
  WRITE(fileunit_callog, '(A,F17.9)') 'Criterium value:', critBest
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''

END SUBROUTINE write_Brent_calibrationlog

!-------------------------------------------
!Block of Quasi-Newton autocal routines
!-------------------------------------------

!> Quasi-Newton DFP or BFGS algorithms for gradient-based minimum search
!----------------------------------------------------------------------
SUBROUTINE quasiNewton_algorithm(critLastVect,parLastTable,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim, numoptimpar, optparid, parindex, fileunit_callog
  USE MODVAR, ONLY   : modparid

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(INOUT) :: critLastVect(optim%cal_improvCritIter) !<criterion saved for last iterations
  REAL, INTENT(INOUT) :: parLastTable(optim%cal_improvParamIter, numoptimpar) !<parameters saved for last iterations
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  
  !Local variables
  INTEGER   npar, dimCounter, iteratCounter
  LOGICAL   finished                                           !Flag for calibration finished through linesearch_methods_interruptor
  LOGICAL   Hessianfinished                                    !Flag for calibration finished through QN_inv_hessian_update with denominator equal to zero
  REAL      gradNorm, lambdaBest, critBest
  REAL      parMin(numoptimpar), parMax(numoptimpar)           !Lower/upper interval limits, from optpar.txt; those limits are absolute, serve as permanent/reference bounds!!
  REAL      parPrecis(numoptimpar)                             !Precision, from optpar.txt
  REAL      par(numoptimpar), gradVect(numoptimpar)            !Current best parameter set and gradient of criterium function there
  REAL      parNew(numoptimpar), gradVectNew(numoptimpar)      !Next best attempt for parameter set and corresponding gradient of criterium
  REAL      deltaVector(numoptimpar), deltaGrad(numoptimpar)   !Difference between current parameter vector and corresponding gradient evaluation
  REAL      direction(numoptimpar)                             !Vectorial direction between current and new best vector
  REAL, dimension(numoptimpar,numoptimpar) :: invHessian, invHessianNew

  !>\b Algorithm \n
  !>Check numerical parameters.
  IF(optim%QN_factorDeriv > 0.5) THEN
    WRITE(6,*) 'Problem in optpar with quasi-Newton: it is not numerically safe to use a factor for numerical derivatives > 0.5'
    WRITE(6,*) 'Adviced is a value around 0.02 (2%) for that numerical parameter'
    STOP 1
  ELSEIF(optim%QN_factorDeriv < 0) THEN
    WRITE(6,*) 'Problem in optpar with quasi-Newton: the factor for numerical derivatives must be positive'
    STOP 1
  ENDIF
  IF(optim%QN_stencil.NE.2 .AND. optim%QN_stencil.NE.4 .AND. optim%QN_stencil.NE.6 .AND. optim%QN_stencil.NE.8) THEN
    WRITE(6,*) 'Problem in optpar with quasi-Newton: stencil parameter must be an integer, and can only be 2, 4, 6 or 8'
    STOP 1
  ENDIF
  IF(optim%QN_lambdaMaxFac > 1) THEN
    WRITE(6,*) 'Problem in optpar with quasi-Newton: the factor to contain lambda is larger than 1'
    WRITE(6,*) 'It must be <= 1, otherwise the line search might exceed the parameter space boundaries'
    STOP 1
  ENDIF

  !>Get the optimisation interval boundaries and the starting parameter set
  CALL find_optpar(par, parMin, parMax, parPrecis, npar)    !Get optimization interval limits from optpar.txt; par is the "center", bogus here; npar = numoptimpar, not very useful either
  CALL load_qNstarting_vector(parMin, parMax, par)      

  !Thirdly, write the calibration header in HYSS.log
  WRITE(6,'(A9)', advance='no') 'Iteration'
  DO dimCounter = 1, numoptimpar
     WRITE(6, '(A16)', advance='no') modparid(optparid(parindex(dimCounter,1)))%shortname
  ENDDO
  WRITE(6,'(A)') '   Criteria       Gradient norm'

  !>Initiate quasi-Newton: evaluate function and gradient at current guess, and initialize inverse Hessian as unit matrix
  iteratCounter = 0
  finished = .FALSE.
  Hessianfinished = .FALSE.
  lambdaBest = 0
  direction = 0.
  invHessian = 0.0      ! For minimization problems, inv. Hessian must be positive definite => take positive unit matrix at start
  DO dimCounter = 1, numoptimpar
    invHessian(dimCounter, dimCounter) = 1.0
  ENDDO
  CALL function_to_minim(lambdaBest,critBest,par,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
  CALL grad_criterium_function(par,parMin,parMax,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,gradVect,gradNorm)
  IF(optim%cal_log) CALL initialize_QN_calibration_log(par, gradVect, gradNorm, invHessian)
  CALL print_calib_HYSSlog(iteratCounter,par,critBest,gradNorm)

  !>Repeatedly perform quasi-Newton steps, as long as the numerical gradient norm is larger than the given threshold
  !>and the maximum allowed number of iterations is not exceeded
  DO WHILE(gradNorm > optim%QN_flatTol)
     iteratCounter = iteratCounter + 1

     !>\li Determine a new parameter set, by linear search of minimum in the direction between current and new best parameter set
     direction = - MATMUL(invHessian, gradVect)
     CALL direction_multiplier(par,parMin,parMax,parPrecis,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,lambdaBest,critBest)

     !>\li Compute the new numerical gradient of the gradient function
     parNew = par + lambdaBest*direction
     CALL grad_criterium_function(parNew,parMin,parMax,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,gradVectNew,gradNorm)

     !>\li Calculate the new inverse Hessian
     deltaVector = parNew - par
     deltaGrad = gradVectNew - gradVect
     CALL QN_inv_hessian_update(deltaVector, deltaGrad, invHessian, invHessianNew, Hessianfinished)
     IF(Hessianfinished) EXIT

     !>\li Update all three quantities: current vector, gradient, and inverse Hessian
     par = parNew
     gradVect = gradVectNew
     invHessian = invHessianNew

     IF(optim%cal_log) CALL write_QN_calibration_log(par, gradVect, gradNorm, invHessian, direction, lambdaBest, critBest, iteratCounter)
     CALL print_calib_HYSSlog(iteratCounter, par, critBest, gradNorm)

     !>\li Check if iterations are enough and to be stopped
     CALL linesearch_methods_interruptor(iteratCounter, critBest, critLastVect, par, parLastTable, parPrecis, finished)
     IF(finished) EXIT

  ENDDO

  IF(.NOT.finished .AND. .NOT.Hessianfinished) CALL linesearch_methods_interruptor_printandstop(5,finished)
  RETURN

END SUBROUTINE quasiNewton_algorithm

!>\brief This subroutine loads the starting set of parameter values
!! for quasi-Newton optimization.
!! Note: parameter order according to optpar.txt, but only the
!! parameters with (larger than zero) intervals. It also checks if
!! given parameter values are within optimization interval boundaries
!------------------------------------------------------------------------
SUBROUTINE load_qNstarting_vector(parMin, parMax, par)

  USE WORLDVAR, ONLY : modeldir,    &     !Directory where to find optmization files 
                       numoptimpar, &     !amount of optimization parameters
                       fileunit_temp    
  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)  :: parMin(numoptimpar)  !<Optimization interval lower boundary
  REAL, INTENT(IN)  :: parMax(numoptimpar)  !<Optimization interval upper boundary
  REAL, INTENT(OUT) :: par(numoptimpar)     !<Vector with the optimization parameters
  
  !Local variables
  INTEGER lineCounter
  CHARACTER (LEN=210) filename
  CHARACTER(LEN=18000) charLine

  filename = TRIM(modeldir)//'qNstartpar.txt'
  OPEN(UNIT=fileunit_temp, FILE=filename, STATUS = 'old', ACTION='read')

  DO lineCounter = 1, numoptimpar
    READ(fileunit_temp, '(a)') charLine
    READ(charLine(10:25),*) par(lineCounter)
    IF(par(lineCounter) < parMin(lineCounter) .OR. par(lineCounter) > parMax(lineCounter)) THEN
      WRITE(6,*) 'Problem in qNstartpar: starting parameter value outside optimization boundaries'
      STOP 1
    ENDIF
  ENDDO

  CLOSE(fileunit_temp)

END SUBROUTINE load_qNstarting_vector

!> Computes the numerical gradient of the gradient function, in all dimensions
!-------------------------------------------------------------------------------
SUBROUTINE grad_criterium_function(par,parMin,parMax,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,gradVect,gradNorm)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim, numoptimpar

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)  :: par(numoptimpar)     !<Current best parameter set = algorithm current position in parameter space; to be tested for gradient
  REAL, INTENT(IN)  :: parMin(numoptimpar)  !<Parameter boundary
  REAL, INTENT(IN)  :: parMax(numoptimpar)  !<Parameter boundary
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  REAL, INTENT(OUT) :: gradVect(numoptimpar)    !<gradient
  REAL, INTENT(OUT) :: gradNorm                 !<gradient norm
  
  !Local variables
  INTEGER   dimCounter, stencilCounter
  INTEGER   derivOrder                      ! Numerical derivatives accuracy order (from the selected stencil)
  REAL      optcrit
  REAL      epsilRel, epsilAbs, epsilVal    ! Offset "epsilon" for numerical derivative, in absolute value (for corresponding parameter value)
  REAL      lambdaVal, derivCoeff, gradValContrib
  REAL      directionVect(numoptimpar)      ! Vector d, zero everywhere except epsilon at component under scrutiny
  REAL      derivCoeffStencil(4,9)

  ! Initialize things
  gradNorm = 0.0        ! Initialize gradient norm, for summing
  derivCoeffStencil(:,1) = (/  0.0,    0.0,    0.0,     1.0/280 /)      ! Coefficients for numerical derivatives, according to stencil
  derivCoeffStencil(:,2) = (/  0.0,    0.0,   -1.0/60, -4.0/105 /)
  derivCoeffStencil(:,3) = (/  0.0,    1.0/12, 3.0/20,  1.0/5 /)
  derivCoeffStencil(:,4) = (/ -1.0/2, -2.0/3, -3.0/4,  -4.0/5 /)
  derivCoeffStencil(:,5) = (/  0.0,    0.0,    0.0,     0.0 /)          ! No points for derivative are EVER taken at current position
  derivCoeffStencil(:,6) = -derivCoeffStencil(:,4)                      ! Coefficient matrix is sign anti-symmetric with respect to middle column (5th)
  derivCoeffStencil(:,7) = -derivCoeffStencil(:,3)
  derivCoeffStencil(:,8) = -derivCoeffStencil(:,2)
  derivCoeffStencil(:,9) = -derivCoeffStencil(:,1)
  derivOrder = optim%QN_stencil/2       ! 2-stencil is 1st order accurate, 4-stencil is 2nd order accurate, etc...

  ! Now, for each dimension of parameter vector, calculate numerical derivative of criteria upon variation of parameter value
  DO dimCounter = 1, numoptimpar

     !Initialize direction and gradient vectors
     CALL get_epsilon(par, parMin, parMax, dimCounter, epsilAbs, epsilRel, epsilVal)     !Determine epsilVal, the numerical offset from current parameter value to numerically approximate derivatives with central differences
     directionVect = 0.0                                 ! Initialize the direction vector as 0
     directionVect(dimCounter) = epsilVal                ! Insert epsilon-value at correct component
     gradVect(dimCounter) = 0.0                          ! Initialize the considered gradient vector component to 0

     !THIS DOES NOT WORK par within 2*epsilVal from border for test case !CP120112
     !Maybe one-sided gradient calculation instead of reducing epsilVal?
     !Check if parameter values will exceed parameter boundaries and reduce the step to keep in bound?

     DO stencilCounter = 1,9
       derivCoeff = derivCoeffStencil(derivOrder,stencilCounter)     ! Coefficient for this specific contribution to numerical derivative
       IF(derivCoeff == 0)THEN   ! If derivative coefficient is zero, no contribution to numerical derivative (and pointless to run model)
         gradValContrib = 0
       ELSE
         lambdaVal = stencilCounter-5
         CALL function_to_minim(lambdaVal,optcrit,par,directionVect,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate)
         gradValContrib = optcrit
       ENDIF
       gradVect(dimCounter) = gradVect(dimCounter) + derivCoeff * gradValContrib / epsilVal      ! Add contribution with coefficient weight to compute derivative
     ENDDO

     gradNorm = gradNorm + gradVect(dimCounter)**2
  ENDDO

  gradNorm = SQRT(gradNorm)

END SUBROUTINE grad_criterium_function

!> Determines "lambda", the fraction of direction step to take, by linear search of minimum
!------------------------------------------------------------------------------------------
SUBROUTINE direction_multiplier(par,parMin,parMax,parPrecis,direction,  &
                frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,lambdaBest,critBest)

  USE STATETYPE_MODULE
  USE WORLDVAR, ONLY : optim, numoptimpar

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)     :: par(numoptimpar)      !<Current algorithm position in parameter space
  REAL, INTENT(IN)     :: parMin(numoptimpar)   !<Lower interval limits, from optpar.txt, these limits are absolute, serve as permanent/reference bounds
  REAL, INTENT(IN)     :: parMax(numoptimpar)   !<Upper interval limits, from optpar.txt, these limits are absolute, serve as permanent/reference bounds
  REAL, INTENT(IN)     :: parPrecis(numoptimpar)  !<Precision, from optpar.txt
  REAL, INTENT(IN)     :: direction(numoptimpar)  !<Vectorial direction between current and new best vector
  TYPE(snowicestatetype),INTENT(INOUT) :: frozenstate   !<Snow and ice states
  TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
  TYPE(aquiferstatetype),INTENT(INOUT) :: aquiferstate  !<Aquifer states
  TYPE(riverstatetype),INTENT(INOUT) :: riverstate  !<River states
  TYPE(lakestatetype),INTENT(INOUT)  :: lakestate   !<Lake states
  TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
  REAL, INTENT(OUT)    :: lambdaBest    !<Best step in direction to take
  REAL, INTENT(OUT)    :: critBest      !<Criterion of best step
  
  !Local variables
  INTEGER  dimCounter
  REAL     maxLambda, maxLambdaBoundary, lambdaAccel, fact
  REAL     limit(numoptimpar), maxLambdaVect(numoptimpar)

  limit = 0.
  maxLambdaVect = 0.
  fact = optim%QN_stencil/2*optim%QN_factorDeriv

  DO dimCounter = 1, numoptimpar
     ! Determine first the "limit displacement" vector, that gives the distance from the current point (par) to the boundary in the direction given by the "direction" vector, while leaving a buffer to the parameter space boundary
     ! This is the max step the algorithm can take in that direction, leaving a buffer with the parameter space boundaries so that it is possible to take points in between for numerical gradient calculations
     ! Defined the way they are defined here, components of "limit" always have the same sign as that of "direction" => ratio lamba always > 0 (see bellow)
     SELECT CASE(optim%QN_epsilType)
     CASE(1)               !Case 1 (default): epsilon is constant ("absolute") because calculated as a fraction of the (constant) parameter space width
        IF(direction(dimCounter) < 0)THEN
           limit(dimCounter) = parMin(dimCounter) + fact*ABS(parMin(dimCounter)) - par(dimCounter)
        ELSE
           limit(dimCounter) = parMax(dimCounter) - fact*ABS(parMax(dimCounter)) - par(dimCounter)
        ENDIF
     END SELECT
     ! Determine the ratio limit/direction as well. That's the vector of the largest possible lambda, for each parameter space dimension
     IF(direction(dimCounter) == 0) THEN     ! Unlikely, but direction = 0 might happen => avoid division by 0 by setting ratio to something very large
        maxLambdaVect(dimCounter) = 1000000000.
     ELSE
        maxLambdaVect(dimCounter) = limit(dimCounter)/direction(dimCounter)
     ENDIF
  ENDDO

  ! Now, determine the maximum fraction of step allowed without crossing a parameter space boundary
  maxLambdaBoundary = MINVAL(maxLambdaVect)     ! Maximum allowed value of "lambda" to remain within boundaries
  lambdaAccel = optim%QN_lambdaAccel    !This parameter allows to accelerate the algorithm stepping, at the expense of max one extra line search iteration only; default value consistant with golden ratio line search algorithm
  !Note: an idea for future development is to update this value from step to step, to reach optimal acceleration, instead of using a constant, golden ration based value
  !Current unresolved question: update with respect to the boundaries, or with respect to the value determined by the line search ?

  ! Compare the step length proposed by the algorithm (lambdaAccel) to see if it exceeds or not the limit imposed by the parameter space boundaries
  ! Note: the old way was simply the case lambdaAccel = 1
  IF(lambdaAccel < optim%QN_lambdaMaxFac * maxLambdaBoundary) THEN      !If algorithm shoots within parameter space ("good acceleration"):
     maxLambda = lambdaAccel                     !  take the algorithm limit as upper limit for the line search
  ELSE                                                      !If algorithm shoots out parameter space ("failed" acceleration):
     maxLambda = optim%QN_lambdaMaxFac * maxLambdaBoundary   !  take the limit imposed by parameter space boundaries as the upper limit for line search, and squeeze it with numerical
  ENDIF                                                     !  parameter "num_lamMax" (to ev. account for space to take numerical derivatives within parameter space boundaries as well!)

  ! Then, perform a line search for the value of lambda in the interval [0, maxLambda] to minimize optcrit
  CALL linesearch_HYSS(0.0,maxLambda,par,parPrecis,direction,frozenstate,soilstate,aquiferstate,riverstate,lakestate,miscstate,lambdaBest,critBest)

  END SUBROUTINE direction_multiplier

!>\brief Determines "epsilon", the offset from parameter value used
!!for numerical gradient calculation (central differences)
!!
!!Note: in this numerical scheme, central differences are implemented
!!so that epsilon MUST have the same sign as the parameter value
!------------------------------------------------------------------------------------------
SUBROUTINE get_epsilon(par, parMin, parMax, dimCounter, epsilAbs, epsilRel, epsilVal)

  USE WORLDVAR, ONLY : optim, numoptimpar

  IMPLICIT NONE

  REAL, INTENT(IN)    :: par(numoptimpar)     ! Current algorithm position in parameter space
  REAL, INTENT(IN)    :: parMin(numoptimpar), parMax(numoptimpar)
  INTEGER, INTENT(IN) :: dimCounter
  REAL, INTENT(OUT)   :: epsilAbs, epsilRel, epsilVal
  !Local variables
  REAL  parLim

  !Calculate the absolute epsilon, defined as a fraction of the parameter space width (optim%QN_factorDeriv is the fraction)
  !and the relative epsilon, defines as a fraction of the parameter value (again, optim%QN_factorDeriv is the fraction)
  IF(par(dimCounter) < 0) THEN
     parLim = parMin(dimCounter)             !Logically, parLim < 0 since parMin(dimCounter) < 0
  ELSE
     parLim = parMax(dimCounter)             !parLim > 0 since parMax(dimCounter) > 0
  ENDIF                                     ! => parLim has the sign of par(dimCounter)
  epsilAbs = optim%QN_factorDeriv * parLim  !Since optim%QN_factorDeriv > 0 (this is checked for), epsilAbs has the sign of par(dimCounter)

  epsilRel = optim%QN_factorDeriv * par(dimCounter) !epsilRel aquires trivially the sign of par(dimCounter)

  !Now, depending on the optimisation instruction, determine the epsilon value to be used
  SELECT CASE(optim%QN_epsilType)
  CASE(1)               !Case 1 (default): epsilon is constant ("absolute") because calculated as a fraction of the (constant) parameter space width
     epsilVal = epsilAbs

  CASE(2)               !Case 2: epsilon is relative, calculated as a fraction of the (current, evolving) parameter value
     epsilVal = epsilRel

  CASE(0)               !Case 0 is mixed: a mobile average between relative and absolute epsilon
     epsilVal = (1 - par(dimCounter)/parLim)*epsilAbs + par(dimCounter)/parLim*epsilRel    !Beyond the fact that parLim has the sign of par(dimCounter), note that 0 < par(dimCounter)/parLim < 1
     ! => epsilVal has the same sign as par(dimCounter) even if par(dimCounter) < 0
  END SELECT                                                                                ! Bottom line: epsilVal has the same sign as par(dimCounter) in all cases :)

  END SUBROUTINE get_epsilon

!------------------------------------------------------------------------------------
!>\brief This subroutine computes the update of the inverse Hessian
!!matrix H for quasi-Newton methods as well as the steepest descent method
!!
!!Very clear math. details on both DFP and BFGS update formula (as well as on the quasi-Newton method in general) can be found in:
!!"Numerical optimization" (second edition), J. Nocedal, S.J. Wright, Springer series in operational research, Springer 2006, chapter 6.1 (p. 135-143)
!!
!!DFP:  H_{k+1} = H_k + (s_k * s_k^T)/(y^T * s_k) - (H_k * y_k * y_k^T * H_k)/(y_k^T * H_k * y_k)                           Eq.(6.15) p. 139 in reference
!!              = H_k + (s_k * s_k^T)/denom       - (H_k * y_k * y_k^T * H_k)/(y_k^T * H_k * y_k)
!!
!!BFGS: H_{k+1} = (s_k * s_k^T)/(y^T * s_k) + [I - (s_k * y_k^T)/(y_k^T * s_k)] * H_k * [I - (y_k * s_k^T)/(y_k^T * s_k)]   Eq.(6.17) p. 140 in reference
!!              =       (s_k * s_k^T)/denom + [I - (s_k * y_k^T)/denom] * H_k * [I - (y_k * s_k^T)/denom]
!!
!!where s_k resp. y_k [column vectors!] are the differences in parameters (deltaVector) resp. gradient (deltaGrad) between step k and k+1,
!!      denom = y^T * s_k, and  H_k = inverse Hessian at step K
!!Notation: vectors are assumed column vectors; all matrix/vector multiplications to be performed left-right
!!
!!Implementing quasi-Newton methods in this way, the steepest descent method turns out to be just a particular case
!!of quasi-Newton method, where the inverse Hessian is always kept as a unit matrix. Easy like on a Sunday morning :)
!------------------------------------------------------------------------------------
SUBROUTINE QN_inv_hessian_update(deltaVector, deltaGrad, invHessian, invHessianNew, finished)

  USE WORLDVAR, ONLY : optim, numoptimpar, fileunit_callog

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)     :: deltaVector(numoptimpar)                 ! s_k
  REAL, INTENT(IN)     :: deltaGrad(numoptimpar)                   ! y_k
  REAL, INTENT(IN)     :: invHessian(numoptimpar, numoptimpar)     ! H_k
  REAL, INTENT(OUT)    :: invHessianNew(numoptimpar, numoptimpar)  ! H_{k+1}
  LOGICAL, INTENT(OUT) :: finished                                 ! flag for zero denominator
  
  !Local variables
  INTEGER dimCounter1, dimCounter2
  REAL denom, normDeltaVector, normDeltaGrad
  REAL, DIMENSION(numoptimpar,numoptimpar) :: commonMatrix, DFPmatrix, BFGSmatrix

  finished = .FALSE.

  ! This preliminary part is specific to quasi-Newton methods and could be skiped by steepest descent. Those preliminary computations and checks
  ! are however executed very fast, so that for the sake of code simplicity we execute them for steepest descent as well, and branch at the end
  !
  ! The matrix  (s_k * s_k^T)/(y^T * s_k)  is commun to both DFP and BFGS expressions; further, the denominator
  ! (y^T * s_k)  appears several times in the BFGS expression  =>  compute everything here!

  ! Calculate denominator first: dot product of deltaGrad and deltaVector [=> scalar, of course]
  denom = 0
  DO dimCounter1 = 1, numoptimpar
     denom = denom + deltaGrad(dimCounter1) * deltaVector(dimCounter1)
  ENDDO

  ! If denominator = 0, it means that the method reached optimum => catch and exit
  IF(denom == 0) THEN
     normDeltaVector = 0
     normDeltaGrad = 0

     DO dimCounter1 = 1, numoptimpar
        normDeltaVector = normDeltaVector + deltaVector(dimCounter1)**2
        normDeltaGrad = normDeltaGrad + deltaGrad(dimCounter1)**2
     ENDDO

     normDeltaVector = SQRT(normDeltaVector)
     normDeltaGrad = SQRT(normDeltaGrad)

     WRITE(6,'(A)')          'Quasi-Newton method: denominator of expression to compute inverse Hessian is zero'
     IF(optim%cal_log) WRITE(fileunit_callog,'(A)') 'Quasi-Newton method: denominator of expression to compute inverse Hessian is zero'

     IF(normDeltaVector == 0) THEN
        WRITE(6,'(A)')          'Norm(deltaVector) = 0: new parameter set was taken identical to previous step!'
        IF(optim%cal_log) WRITE(fileunit_callog,'(A)') 'Norm(deltaVector) = 0: new parameter set was taken identical to previous step!'
     ENDIF

     IF(normDeltaGrad == 0) THEN
        WRITE(6,'(A)')          'Norm(deltaGrad) = 0: gradient of target function at new parameter set is identical to gradient at previous parameter set!'
        IF(optim%cal_log) WRITE(fileunit_callog,'(A)') 'Norm(deltaGrad) = 0: gradient of target function at new parameter set is identical to gradient at previous parameter set!'
     ENDIF

     CALL linesearch_methods_interruptor_printandstop(6,finished)
     RETURN
  ENDIF

  ! Compute the matrix  (s_k * s_k^T)/(y^T * s_k)  common to both DFP and BFGS expressions
  DO dimCounter1 = 1, numoptimpar     ! Numerator is matrix obtained from deltaVector * deltaVector^T; divide by denominator directly
     DO dimCounter2 = 1, numoptimpar
        commonMatrix(dimCounter1, dimCounter2) = deltaVector(dimCounter1) * deltaVector(dimCounter2) / denom
     ENDDO
  ENDDO

  ! Method specific contributions and final expression for updated inverse Hessian
  IF(optim%task_stpstDesc) THEN
     invHessianNew = 0.0
     DO dimCounter1 = 1, numoptimpar
        invHessianNew(dimCounter1, dimCounter1) = 1.0           ! Steepest descent: inverse Hessian is just a unit matrix
     ENDDO

  ELSEIF(optim%task_DFP) THEN
     CALL DFP_inv_hessian_update(deltaGrad, invHessian, DFPmatrix)
     invHessianNew = invHessian + commonMatrix - DFPmatrix     ! DFP inverse Hessian update formula

  ELSEIF(optim%task_BFGS) THEN
     CALL BFGS_inv_hessian_update(deltaVector, deltaGrad, invHessian, denom, BFGSmatrix)
     invHessianNew = commonMatrix + BFGSmatrix                 ! BFGS inverse Hessian update formula

  ENDIF

END SUBROUTINE QN_inv_hessian_update

!>\brief This routine computes the last term (DFPmatrix)
!! DFP: H_{k+1} = H_k + [s_k * s_k^T]/[y^T * s_k] - [H_k * y_k * y_k^T * H_k]/[y_k^T * H_k * y_k]
!! Eq. (6.15) p. 139 in "Numerical optimization" (second edition), J. Nocedal, S.J. Wright
!---------------------------------------------------------------------------------
SUBROUTINE DFP_inv_hessian_update(deltaGrad, invHessian, DFPmatrix)

  USE WORLDVAR, ONLY : numoptimpar

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)  :: deltaGrad(numoptimpar)               !<y_k, difference in gradient at step k
  REAL, INTENT(IN)  :: invHessian(numoptimpar, numoptimpar) !<H_k, inverse Hessian at step k
  REAL, INTENT(OUT) :: DFPmatrix(numoptimpar, numoptimpar)  !<Last term of the DFP upgrade formula
  
  !Local variables
  INTEGER dimCounter1, dimCounter2
  REAL denom
  REAL dummyVector(numoptimpar)
  REAL, DIMENSION(numoptimpar,numoptimpar) :: dummyMatrix1, dummyMatrix2

  ! Denominator is given by  y_k^T * H_k * y_k (to be performed left-right, otherwise it's not a scalar!)
  !First compute  deltaGrad^T * H;  it gives an intermediate ("dummy") line vector
  dummyVector = 0.0       ! Initialize to 0
  DO dimCounter1 = 1, numoptimpar
     DO dimCounter2 = 1, numoptimpar
        dummyVector(dimCounter1) = dummyVector(dimCounter1) + deltaGrad(dimCounter2) * invHessian(dimCounter2, dimCounter1)
     ENDDO
  ENDDO

  !Now, the denominator is given by the dot product of  (deltaGrad^T * H)  and  deltaGrad
  denom = 0       ! Initialize to 0
  DO dimCounter1 = 1, numoptimpar
     denom = denom + dummyVector(dimCounter1) * deltaGrad(dimCounter1)     !y_k^T * H_k * y_k
  ENDDO

  ! Numerator is matrix given by  H_k * y_k * y_k^T * H_k  (take left-right order for consistency)
  !First compute  H * deltaGrad  => it gives an intermediate, "dummy" column vector
  dummyVector = 0.0       ! Initialize to 0
  DO dimCounter1 = 1, numoptimpar
     DO dimCounter2 = 1, numoptimpar
        dummyVector(dimCounter1) = dummyVector(dimCounter1) + invHessian(dimCounter1, dimCounter2) * deltaGrad(dimCounter2)
     ENDDO
  ENDDO

  !Secondly, compute the product of the dummy vector with deltaGrad^T; it gives an intermediate, "dummy" matrix
  dummyMatrix1 = 0.0      ! Initialize to 0 (not necessary, right?)
  DO dimCounter1 = 1, numoptimpar
     DO dimCounter2 = 1, numoptimpar
        dummyMatrix1(dimCounter1, dimCounter2) = dummyVector(dimCounter1) * deltaGrad(dimCounter2)
     ENDDO
  ENDDO

  !Thirdly, compute the matrix product between the dummy matrix and the inverse Hessian; that's the numerator
  dummyMatrix2 = 0.0      ! Initialize to 0 (not necessary either)
  dummyMatrix2 = MATMUL(dummyMatrix1, invHessian)     !H_k * y_k * y_k^T * H_k

  DFPmatrix = dummyMatrix2 / denom    ! DFP specific contribution term: [H_k * y_k * y_k^T * H_k]/[y_k^T * H_k * y_k]

  END SUBROUTINE DFP_inv_hessian_update

!>\brief This routine computes the last term (BFGSmatrix)
!>
!!BFGS: H_{k+1} = (s_k * s_k^T)/(y^T * s_k) + [I - (s_k * y_k^T)/(y_k^T * s_k)] * H_k * [I - (y_k * s_k^T)/(y_k^T * s_k)]
!!Eq. (6.17) p. 140 in "Numerical optimization" (second edition), J. Nocedal, S.J. Wright
!------------------------------------------------------------------------------------------
SUBROUTINE BFGS_inv_hessian_update(deltaVector, deltaGrad, invHessian, denom, BFGSmatrix)

  USE WORLDVAR, ONLY : numoptimpar

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)    :: denom        ! deltaGrad^T * deltaVector, calculated in the main quasi-Newton subroutine
  REAL, INTENT(IN)    :: deltaVector(numoptimpar), deltaGrad(numoptimpar)
  REAL, INTENT(IN)    :: invHessian(numoptimpar, numoptimpar)
  REAL, INTENT(OUT)   :: BFGSmatrix(numoptimpar, numoptimpar)
  
  !Local variables
  INTEGER dimCounter1, dimCounter2
  REAL, DIMENSION(numoptimpar,numoptimpar) :: unitMatrix, dummyMatrix, matrix1, matrix3

  ! First put together a unit matrix of "numoptimpar" size
  unitMatrix = 0.0
  DO dimCounter1 = 1, numoptimpar
     unitMatrix(dimCounter1, dimCounter1) = 1.0
  ENDDO

  ! Second, get the matrix given by deltaGrad * deltaVector^T, and divide by denominator deltaGrad^T * deltaVector (denominator passed as argument)
  dummyMatrix = 0.0
  DO dimCounter1 = 1, numoptimpar
     DO dimCounter2 = 1, numoptimpar
        dummyMatrix(dimCounter1, dimCounter2) = deltaGrad(dimCounter1) * deltaVector(dimCounter2) / denom   !(y_k * s_k^T)/(y_k^T * s_k)
     ENDDO
  ENDDO

  ! Third, withdraw from the unit matrix and build the terms 1 and 3 of the BFGS expression
  matrix1 = unitMatrix - TRANSPOSE(dummyMatrix)       ![I - (s_k * y_k^T)/(y_k^T * s_k)];  s_k * y_k^T  is the transpose of  y_k * s_k^T, the "dummy matrix"
  matrix3 = unitMatrix - dummyMatrix                  ![I - (y_k * s_k^T)/(y_k^T * s_k)]

  ! Finally, the BFGS contribution is given by  matrix1 * H * matrix3 (from left to right); it's ok to re-use the "dummy" matrix, since no longer useful
  dummyMatrix = MATMUL(matrix1, invHessian)
  BFGSmatrix = MATMUL(dummyMatrix, matrix3)

END SUBROUTINE BFGS_inv_hessian_update

!-----------------------------------------------------------------------
SUBROUTINE print_calib_HYSSlog(iterCounter, param, critVal, gradNorm)

  USE WORLDVAR, ONLY : numoptimpar

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: iterCounter
  REAL, INTENT(IN)    :: critVal, gradNorm
  REAL, INTENT(IN)    :: param(numoptimpar)
  
  !Local variables
  INTEGER i

  WRITE(6,'(I7)', advance='no') iterCounter
  DO i = 1, numoptimpar
    WRITE(6,'(F16.9)', advance='no') param(i)
  ENDDO
  WRITE(6,'(F16.9, F16.9)') critVal, gradNorm

END SUBROUTINE print_calib_HYSSlog

!-----------------------------------------------------------------------
SUBROUTINE print_calib_HYSSlog_noimp(iterCounter, param, critVal)

  USE WORLDVAR, ONLY : numoptimpar

  IMPLICIT NONE

  !Argument declarations
  INTEGER, INTENT(IN) :: iterCounter
  REAL, INTENT(IN)    :: critVal
  REAL, INTENT(IN)    :: param(numoptimpar)
  
  !Local variables
  INTEGER i

  WRITE(6,'(I7)', advance='no') iterCounter
  DO i = 1, numoptimpar
    WRITE(6,'(F16.9)', advance='no') param(i)
  ENDDO
  WRITE(6,'(F16.9, A19)') critVal, '     no improvement'

END SUBROUTINE print_calib_HYSSlog_noimp

!------------------------------------------------------------------------------
SUBROUTINE initialize_QN_calibration_log(par, gradVect, gradNorm, invHessian)

  USE WORLDVAR, ONLY : optim,           &
                       modeldir,        &
                       resdir,          & 
                       numoptimpar,     &
                       filename_callog, & 
                       fileunit_callog

  IMPLICIT NONE

  !Argument declarations
  REAL, INTENT(IN)    :: gradNorm
  REAL, INTENT(IN)    :: par(numoptimpar), gradVect(numoptimpar)
  REAL, INTENT(IN)    :: invHessian(numoptimpar,numoptimpar)
  !Local variables
  INTEGER counter1, counter2

  ! Start printing log file
  OPEN(UNIT = fileunit_callog, FILE = TRIM(resdir)//filename_callog, STATUS = 'replace')   ! Start calibration log file

  IF(optim%task_stpstDesc) THEN
    WRITE(fileunit_callog, '(A)') 'Calibration by steepest descent method'
    WRITE(fileunit_callog, '(A)') '--------------------------------------'
  ELSEIF(optim%task_DFP) THEN
    WRITE(fileunit_callog, '(A)') 'Calibration by DFP quasi-Newton algortihm'
    WRITE(fileunit_callog, '(A)') '-----------------------------------------'
  ELSEIF(optim%task_BFGS) THEN
    WRITE(fileunit_callog, '(A)') 'Calibration by BFGS quasi-Newton algortihm'
    WRITE(fileunit_callog, '(A)') '------------------------------------------'
  ENDIF

  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, '(A, A)') 'Model area: ', TRIM(modeldir)
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, '(A)') 'Numerical parameters:'
  WRITE(fileunit_callog, '(A, I4)')    'Maximum amount of allowed quasi-Newton iterations             : ', optim%cal_maxIterat
  WRITE(fileunit_callog, '(A, F13.9)') 'Threshold for gradient norm to be considered zero             : ', optim%QN_flatTol
  WRITE(fileunit_callog, '(A, F10.6)')  'Factor of parameter taken as offset for numerical derivatives : ', optim%QN_factorDeriv
  WRITE(fileunit_callog, '(A, I1)')    'Stencil type                                                  : ', optim%QN_stencil
  WRITE(fileunit_callog, '(A, F7.3)')  'Squeeze to contain lambda factor                              : ', optim%QN_lambdaMaxFac
  WRITE(fileunit_callog, '(A, F13.9)') 'Tolerance for line search                                     : ', optim%lineSearch_tol
  WRITE(fileunit_callog, '(A, I4)')    'Max amount of iterations allowed in line search               : ', optim%lineSearch_maxIter

  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''

  WRITE(fileunit_callog, '(A)') 'Iteration 0 (start state)'
  WRITE(fileunit_callog, '(A)') '-------------------------'
  WRITE(fileunit_callog, *) ''

  WRITE(fileunit_callog, '(A, F17.9)') 'Parameter set:              Crit. grad. at param set:              Gradient norm:', gradNorm
  DO counter1 = 1, numoptimpar
    WRITE(fileunit_callog, '(F17.9, A, F17.9)', advance='no') par(counter1), '             ', gradVect(counter1)
    WRITE(fileunit_callog, '(A)') ' '
  ENDDO
  WRITE(fileunit_callog, *) ''

  WRITE(fileunit_callog, '(A)') 'Inverse Hessian:'
  DO counter1 = 1, numoptimpar
    DO counter2 = 1, numoptimpar
      WRITE(fileunit_callog, '(F17.9, A)', advance='no') invHessian(counter1, counter2), '   '
    ENDDO
    WRITE(fileunit_callog, '(A)') ' '
  ENDDO
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''
  WRITE(fileunit_callog, *) ''

END SUBROUTINE initialize_QN_calibration_log

!-----------------------------------------------------------------
  SUBROUTINE write_QN_calibration_log(par, gradVect, gradNorm, invHessian, &
     direction, lambdaBest, critBest, iteratCounter)

    USE WORLDVAR, ONLY : numoptimpar, fileunit_callog

    IMPLICIT NONE

    !Argument declarations
    INTEGER, INTENT(IN) :: iteratCounter
    REAL, INTENT(IN)    :: gradNorm, lambdaBest, critBest
    REAL, INTENT(IN)    :: par(numoptimpar), gradVect(numoptimpar), direction(numoptimpar)
    REAL, INTENT(IN)    :: invHessian(numoptimpar,numoptimpar)

    !Local variables
    INTEGER counter1, counter2

    WRITE(fileunit_callog, '(A, I4)') 'Iteration ', iteratCounter
    WRITE(fileunit_callog, '(A, I4)') '--------------', iteratCounter
    WRITE(fileunit_callog, *) ''

    WRITE(fileunit_callog, '(A, F17.9)') 'Step direction:              Optimal lambda:', lambdaBest
    DO counter1 = 1, numoptimpar
      WRITE(fileunit_callog, '(F17.9)') direction(counter1)
    ENDDO
    WRITE(fileunit_callog, *) ''

    WRITE(fileunit_callog, '(A, F17.9, A, F17.9)') 'Parameter set:              Crit. grad. at param set:          Criterium value:', critBest, '      Gradient norm:', gradNorm
    DO counter1 = 1, numoptimpar
      WRITE(fileunit_callog, '(F17.9, A, F17.9)', advance='no') par(counter1), '             ', gradVect(counter1)
      WRITE(fileunit_callog, '(A)') ' '
    ENDDO
    WRITE(fileunit_callog, *) ''

    WRITE(fileunit_callog, '(A)') 'Inverse Hessian:'
    DO counter1 = 1, numoptimpar
      DO counter2 = 1, numoptimpar
        WRITE(fileunit_callog, '(F17.9, A)', advance='no') invHessian(counter1, counter2), '   '
      ENDDO
      WRITE(fileunit_callog, '(A)') ' '
    ENDDO
    WRITE(fileunit_callog, *) ''
    WRITE(fileunit_callog, *) ''
    WRITE(fileunit_callog, *) ''
    WRITE(fileunit_callog, *) ''

  END SUBROUTINE write_QN_calibration_log
     
END MODULE OPTIMIZATION
