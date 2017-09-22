!> \file worvar.f90
!> Contains module worldvar.

!> \brief Module for declaration of HYSS variables. These are NOT to be used
!>in the model. 
!>
!>NOTE: everything is public.
!>
!>The module worldvar contains constants, variables, and procedures that are used 
!>globally within HYSS (but not within the model). It also holds procedures 
!>for setting simulation information, handling memory, and getting single observation values.
!>
!>The constants define file names and settings for I/O, dimensions of arrays, 
!>code index for accumulation period and input variable type etc.
!>
!>The variables are for model simulation settings, for holding observations both  
!>for forcing and for criteria calculations, for temporary accumulation and calculation 
!>of output, for optimisation setting and results and for calculation of criteria.
MODULE WORLDVAR

  !Copyright 2011-2017 SMHI
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

  !----------------------------------------------------------------------
  !Procedures in this module
  !----------------------------------------------------------------------
  ! da_allocate_accumulation
  ! allocate_accumulation
  ! reallocate_outvar_information
  ! get_current_date_memory
  ! get_current_qobs
  ! get_current_xobs
  ! get_current_xoregobs
  ! get_current_forcing_from_memory
  ! set_outvar
  ! set_output
  ! set_calibration
  ! set_maxmap
  ! deallocate_worldvar
  ! allocate_MCvariables
  ! deallocate_MCvariables
  ! set_timeformat
  ! fileunit_get
  ! fileunit_free
  ! get_seq_filename
  ! count_optim_par
  ! allocate_and_initialize_forcingdata_structure
  !----------------------------------------------------------------------

  USE MODVAR,  ONLY: missing_value,   &
                     seconds_per_timestep
  USE LIBDATE, ONLY: DateType

  IMPLICIT NONE
  SAVE

!> \name Miscellaneous constant parameters and variables
!> \{
  !Dimension parameters
  INTEGER, PARAMETER :: maxcharpath = 200    !<Number of characters in path variables
  INTEGER, PARAMETER :: maxcrit = 20         !<Maximum numbers of criteria
  INTEGER, PARAMETER :: maxoptpar = 100      !<Maximum different parameters in optimization
  !Time related parameters
  REAL,   PARAMETER :: meandaysofyear = 365.25  !<Average number of days in a year
  INTEGER,PARAMETER :: seconds_per_day = 86400  !<Seconds per day

  CHARACTER(LEN=2), PARAMETER :: comment_str ='!!'     !<Comment string code
  !Variable for submodel handling
  INTEGER,ALLOCATABLE :: ibasemodel(:)          !<Index in basemodel as a function of index in submodel
!> \}
!> \name Constant file parameters
!> Files that are frequently used have their own file unit number and file name parameter. 
!> Files opened and closed in sequence use the temporary fileunit 101.
!> \{
  INTEGER, PARAMETER :: fileunit_base = 200    !<Fileunit start range for WB-output and XOMS-files
  INTEGER, PARAMETER :: max_files = 140 !10000000    !<Maximum number of files open at the same time (enough?) TODO: use fileget for all output files
  LOGICAL            :: fileunits(max_files)   !<File units available or occupied (fileunit_base+1--fileunit_base+max_files)
  INTEGER, PARAMETER :: fileunit_temp  = 101   !<temporary fileunit for files that are opened and closed in the same routine, e.g. info,par,optpar,respar,Qobs,Pobs,Tobs,Xobs (check station),GeoClass,GeoData
  INTEGER, PARAMETER :: fileunit_qobs  = 104   !<Fileunit for disckarge observation data
  INTEGER, PARAMETER :: fileunit_xobs  = 105   !<Fileunit for other observations data
  INTEGER, PARAMETER :: fileunit_MC    = 106   !<Fileunit for allsim.txt
  INTEGER, PARAMETER :: fileunit_callog = 107  !<Fileunit for calibration.log
  INTEGER, PARAMETER :: fileunit_xoreg  = 108   !<Fileunit for observations data on output regions
  INTEGER, PARAMETER :: fileunit_add_sub  = 10000000   !<Number series for basin output files (10000001-10099999)
  INTEGER, PARAMETER :: fileunit_add_time = 1000       !<Number series for time output files (1001-999999)
  INTEGER, PARAMETER :: fileunit_add_timeio = 100      !<Number series for time output files (max 100 per output period)
  INTEGER, PARAMETER :: fileunit_add_dabout  = 100000  !<Number series for DA extra basin output files (10100001-10999999)
  CHARACTER(LEN=8), PARAMETER  :: filename_Qobs   = 'Qobs.txt'
  CHARACTER(LEN=8), PARAMETER  :: filename_Xobs   = 'Xobs.txt'
  CHARACTER(LEN=12), PARAMETER :: filename_Xoreg  = 'Xoregobs.txt'
  CHARACTER(LEN=10), PARAMETER :: filename_MC     = 'allsim.txt'
  CHARACTER(LEN=12), PARAMETER :: filename_best   = 'bestsims.txt'
  CHARACTER(LEN=10), PARAMETER :: filename_upd    = 'update.txt'
  CHARACTER(LEN=15), PARAMETER :: filename_callog = 'calibration.log'
!> \}
!> \name Data type parameters
!> These codes are used when reading model run settings and model parameter values from files.
!> \{
  INTEGER,PARAMETER :: i_str=0  !<Data type parameter, code for reading data: string (not read)
  INTEGER,PARAMETER :: i_intg=1 !<Data type parameter, code for reading data: integer
  INTEGER,PARAMETER :: i_real=2 !<Data type parameter, code for reading data: real
!> \}

!> \name Forcing data constants
!> Codes for forcing data and other observations in forcing data structure
!> \{
  INTEGER, PARAMETER :: i_pobs = 1  !<Index of Pobs data in forcingdata structure
  INTEGER, PARAMETER :: i_tobs = 2  !<Index of Tobs data in forcingdata structure
  INTEGER, PARAMETER :: i_tminobs = 3  !<Index of TMINobs data in forcingdata structure
  INTEGER, PARAMETER :: i_tmaxobs = 4  !<Index of TMAXobs data in forcingdata structure
  INTEGER, PARAMETER :: i_rhobs = 5    !<Index of RHobs data in forcingdata structure
  INTEGER, PARAMETER :: i_sfobs = 6    !<Index of SFobs data in forcingdata structure
  INTEGER, PARAMETER :: i_swobs = 7    !<Index of SWobs data in forcingdata structure
  INTEGER, PARAMETER :: i_uobs = 8     !<Index of Uobs data in forcingdata structure
  INTEGER, PARAMETER :: max_forcingdata = 8 !<Maximum number of forcing data variables
!> \}
!> \name Criterion related constant parameters
!> \brief An index code is used for simulation assessment criteria-arrays. Short names for criteria are defined. 
!> \{
  INTEGER,PARAMETER :: i_rnse   = 1   !<Simulation assessment criteria: Regional R2 (NSE)
  INTEGER,PARAMETER :: i_snse   = 2   !<Simulation assessment criteria: Spatial R2 (NSE)
  INTEGER,PARAMETER :: i_mnse   = 3   !<Simulation assessment criteria: Aritmetic mean R2 (NSE)
  INTEGER,PARAMETER :: i_rmae   = 4   !<Simulation assessment criteria: Regional Mean Absolute Error
  INTEGER,PARAMETER :: i_sre    = 5   !<Simulation assessment criteria: spatial RE
  INTEGER,PARAMETER :: i_rre    = 6   !<Simulation assessment criteria: regional RE
  INTEGER,PARAMETER :: i_mre    = 7   !<Simulation assessment criteria: arithmetic mean RE
  INTEGER,PARAMETER :: i_rra    = 8   !<Simulation assessment criteria: regional RA
  INTEGER,PARAMETER :: i_sra    = 9   !<Simulation assessment criteria: spatial RA
  INTEGER,PARAMETER :: i_mra    = 10  !<Simulation assessment criteria: arithmetic mean RA
  INTEGER,PARAMETER :: i_tau    = 11  !<Simulation assessment criteria: Kendalls Tau
  INTEGER,PARAMETER :: i_mdnse  = 12  !<Simulation assessment criteria: Median NSE
  INTEGER,PARAMETER :: i_mdra   = 13  !<Simulation assessment criteria: Median RA
  INTEGER,PARAMETER :: i_mstdre = 14  !<Simulation assessment criteria: Arithmetic mean std relative error
  INTEGER,PARAMETER :: i_mcc    = 15  !<Simulation assessment criteria: Arithmetic mean correlation coefficient
  INTEGER,PARAMETER :: i_mdkg   = 16  !<Simulation assessment criteria: Median Kling-Gupta Efficiency
  INTEGER,PARAMETER :: i_mabsre = 17  !<Simulation assessment criteria: Mean Absolute Relative Error
  INTEGER,PARAMETER :: i_mnrmse = 18  !<Simulation assessment criteria: Median Normalized Root Mean Square Error
  INTEGER,PARAMETER :: i_mnw    = 19  !<Simulation assessment criteria: Mean NSE adjusted for bias
  INTEGER,PARAMETER :: i_numrc  = 20  !<Simulation assessment: Number of data points for regional criteria
  INTEGER,PARAMETER :: i_nummc  = 21  !<Simulation assessment: Number of data points for mean or median criteria
  INTEGER,PARAMETER :: maxperf  = 21  !<Number of simulation assessment criteria
  INTEGER,PARAMETER :: maxsubass = 17 !<Number of subbasin assessment performance criteria
  CHARACTER(LEN=5),PARAMETER :: performance_name(maxperf) = (/'rr2  ','sr2  ',    &
                                'mr2  ','rmae ','sre  ','rre  ','mre  ','rra  ',  &
                                'sra  ','mra  ','tau  ','md2  ','mda  ','mrs  ',  &
                                'mcc  ','mdkg ','mare ','mnr  ','mnw  ','numrc','nummc'/) !<Heading for performance criteria file
!> \}

  !Type declarations
!> \brief Type for holding information about forcing data
  TYPE FORCINGDATATYPE
    LOGICAL :: readfile   !<status of file
    INTEGER :: fileunit   !<file unit number
    CHARACTER(LEN=maxcharpath) :: filepath  !<file path
    CHARACTER(LEN=10) :: filename       !<significant file name (e.g. RHobs)
    CHARACTER(LEN=12) :: infocode       !<code in info.txt to read file
    CHARACTER(LEN=10) :: idcode       !<code in headings for observation id
    INTEGER,ALLOCATABLE :: basinindex(:)   !<index to find data for each subbasin (nsub)
    INTEGER,ALLOCATABLE :: stationid(:)    !<forcing data id which is coupled to subbasin
    INTEGER          :: ncols       !<number of columns of data
    REAL,ALLOCATABLE :: allvalues(:,:)    !<matrix for keeping forcing data in memory (time,ncols)
  END TYPE FORCINGDATATYPE

!> \brief Type for holding information about criteria
  TYPE CALIBRATIONINFO
    CHARACTER(LEN=3) :: crit           !<criterion short name
    INTEGER          :: comp           !<outvar number
    INTEGER          :: rec            !<outvar number
    INTEGER          :: vartype        !<state,flow or conc
    INTEGER          :: areaagg        !<basin,upstream or regional
    REAL             :: weight = 0.    !<weight
    REAL             :: coeff  = 0.    !<parameter for criterion
    LOGICAL          :: cond   = .FALSE. !<To enable conditional acceptence with this criterion
    REAL             :: thres  = 0.      !<Threshold for conditional criterion
  END TYPE CALIBRATIONINFO

!> \brief Type for holding information about optimization
  TYPE OPTIMIZATIONTYPE
    LOGICAL  :: task_MC          = .FALSE.  !<MonteCarlo simulation
    LOGICAL  :: task_boundps     = .FALSE.  !<Repeated MonteCarlo simulation with reduced parameter space
    LOGICAL  :: task_Scanning    = .FALSE.  !<Scan mode, for 2-param situations only
    LOGICAL  :: task_BrentNew    = .FALSE.  !<Calibration using the Brent method (2010 version with new line search)
    LOGICAL  :: task_stpstDesc   = .FALSE.  !<Steepest descent, implemented as a particular case of quasi-Newton gradient-based calibration 
    LOGICAL  :: task_DFP         = .FALSE.  !<Quasi-Newton gradient-based calibration with DFP algorithm for inverse Hessian update 
    LOGICAL  :: task_BFGS        = .FALSE.  !<Quasi-Newton gradient-based calibration with BFGS algorithm for inverse Hessian update
    LOGICAL  :: task_stageMC     = .FALSE.  !<MonteCarlo simulations using stage-wise centering around set of bests MC runs
    LOGICAL  :: task_runens      = .FALSE.  !<Running ensembles and write best result from ensemble runs
    LOGICAL  :: task_writeall    = .FALSE.  !<Write performance for all simulations from ensemble runs
    LOGICAL  :: task_writesim    = .FALSE.  !<Write simulations results for all simulations from ensemble runs (MC-methods only)
    LOGICAL  :: task_DEMC        = .FALSE.  !<DEMC Differential-Evolution Markov Chain (Monte Carlo) simulation

    INTEGER  :: scan_xpoints     = 1        !<Number of points taken for 1st parameter
    INTEGER  :: scan_ypoints     = 1        !<Number of points taken for 2nd parameter

    INTEGER  :: nruns_MC         = 1000     !<Number of runs for Monte Carlo simulation
    INTEGER  :: nruns_MCI        = 200      !<Number of runs for each iteration reducing the parameter space
    INTEGER  :: nruns_MCImax     = 100      !<Maximum number of iteration for reducing the parameter space
    INTEGER  :: nruns_best       = 1        !<Number of ensembles = number of best runs saved from random simulations
    INTEGER  :: nruns_stages     = 1        !<Number of stages to repeat num_mc runs, centering around the best runs at each stage [Fred, 26.08.10]
    REAL     :: nruns_zoom       = 0.9      !<Zooming factor for the centering [Fred, 26.08.10]

    LOGICAL  :: cal_log             = .TRUE. !<Flag to write all calibration details to file "calibration.log" (Y/N flag; Y = .TRUE., N = .FALSE.)
    INTEGER  :: cal_debugCase       = 0     !<Flag to indicate debug case
    INTEGER  :: cal_maxIterat       = 100   !<Max amount of allowed iterations
    REAL     :: cal_maxTime         = 72    !<Max amout of time (hours) allowed to calibration routine
    INTEGER  :: cal_improvParamIter = 10    !<Amount of last optimisation iterations taken into account for parameter improvement monitoring
    INTEGER  :: cal_improvCritIter  = 10    !<Amount of last optimisation iterations taken into account for criteria improvement monitoring
    REAL     :: cal_improvCritTol   = 0.001 !<Tolerance to consider criteria as optimised (delta/mean over the last iterations)

    REAL     :: QN_flatTol          = 0.001 !<Tolerance for gradient norm to be considered zero (quasi-Newton)
    INTEGER  :: QN_epsilType        = 1     !<Specifies whether to use absolute (case 1, default), relative (case 2) or mixed epsilon (case 0) values
    REAL     :: QN_factorDeriv      = 0.02  !<Factor to offset current parameter value for numerical derivative
    INTEGER  :: QN_stencil          = 2     !<Numerical derivative stencil type
    REAL     :: QN_lambdaMaxFac     = 0.90  !<Factor to contain lambda prior to line search (allows to leave space between current point and boundary, for numerical derivatives)
    REAL     :: QN_lambdaAccel      = 1.618 !<Factor increasing the step length proposed by QN algorithms (case lambda = 1 to be replaced by lambdaAccel), in order to allow for faster iteration progression; default value consistant with golden ratio line search algorithm

    LOGICAL  :: Brent_diagonalStep = .TRUE. !<Flag to take a diagonal step at the end of each Brent iteration

    INTEGER  :: lineSearch_maxIter = 500    !<Maximum amount of iterations allowed for line search algorithm
    REAL     :: lineSearch_tol     = 0.001  !<Tolerance for line search of minimum [Fred]
 
    INTEGER  :: DEMC_ngen          = 100    !<Number of generations
    INTEGER  :: DEMC_npop          = 25     !<Number of populations
    REAL     :: DEMC_gammascale    = 1.     !<Multiplicative scaling of the mutation strength gamma = 2.38/(2*num_par^2)^0.5
    REAL     :: DEMC_crossover     = 1.     !<Crossover probability (default=1)
    REAL     :: DEMC_sigma         = 0.1    !<Sample Error Standard deviation
    INTEGER  :: DEMC_accprob       = 0      !<Switch on/off probabilistic acceptance
  END TYPE OPTIMIZATIONTYPE

!> \brief Type for holding information about accumulation of output  
  TYPE ACCUMULATIONINFO
    INTEGER :: comp              !<outvarid number
    INTEGER :: rec               !<outvarid number
    INTEGER :: compoutvar        !<outvar index
    INTEGER :: recoutvar         !<outvar index
    INTEGER :: vartype           !<state,flow or conc
    INTEGER :: areaagg           !<basin,upstream or regional(?)
    LOGICAL :: saveend = .FALSE. !<criteria that need values to be saved until last and then calculated
  END TYPE ACCUMULATIONINFO

!> \brief Types for holding information about wanted output  
!> \{
!> \brief Type for information about content in array with selected output from model (outvar)
  TYPE OUTVARINFOTYPE
    INTEGER        :: idindex        !<variable's index in outvarid
    INTEGER        :: areaagg        !<type of area aggregation; 0=none,1=upstream,2=region output
    INTEGER        :: ovindex0       !<non-aggregated variable's index in outvar
    INTEGER        :: ovindexwater   !<conc variable's index to water in outvar
  END TYPE OUTVARINFOTYPE
!> \brief Type for information about variables for output  
  TYPE OUTPUTVARIABLETYPE
    INTEGER        :: idindex        !<variable's index in outvarid
    INTEGER        :: ovindex        !<variable's index in outvar
    INTEGER        :: areaagg        !<type of area aggregation; 2=region output
  END TYPE OUTPUTVARIABLETYPE
!> \brief Type for accumulated output data
  TYPE OUTPUTDATATYPE
    REAL,ALLOCATABLE    :: value(:,:,:) !<Accumulated data for print out at wperiod (nsub/nreg,nvar,nens)
    REAL,ALLOCATABLE    :: help(:,:,:)  !<Accumulated data of volume for concentrations     
    INTEGER,ALLOCATABLE :: nok(:,:,:)   !<Flag for missing value in accumulation period (0/1) or number of accumulated values in period (maybe)
  END TYPE OUTPUTDATATYPE
!> \brief Type for information about wanted output  
  TYPE OUTPUTTYPE
    INTEGER        :: fileformat        !<code for file format; 4=region output
    INTEGER        :: period            !<time period to calculate average or sum
    INTEGER        :: decimal           !<number of decimals of output
    INTEGER        :: signfig           !<number of significant figures of output
    INTEGER        :: nvar              !<number of variables
    INTEGER        :: narea             !<number of output areas (basin and region files (so far)) (narea=0 means all?)
    TYPE(OUTPUTVARIABLETYPE),ALLOCATABLE :: variable(:) !<variables of output
    INTEGER,ALLOCATABLE  :: areaindex(:) !<index of subbasin or outregion to output
    TYPE(OUTPUTDATATYPE) :: accdata     !<Accumulated output data
  END TYPE OUTPUTTYPE
!> \brief Type for information about output regions
  TYPE OUTREGIONINFOTYPE
    INTEGER :: outregid    !<identification number for output region (must not overlap with subid)
    REAL    :: xcoord      !<x-coordinate
    REAL    :: ycoord      !<y-coordinate
    REAL    :: zcoord      !<z-coordinate, eg. elevation above sea level
    INTEGER :: nsubbasin        !<number of subbasin in region
    INTEGER,ALLOCATABLE :: subindex(:)   !<index of subbasin in region (nsubbasin)    !handle submodel!!?
    REAL,ALLOCATABLE :: weight(:)     !<weight of subbaisn in region average (nsubbasin)
  END TYPE OUTREGIONINFOTYPE
!> \}

!> \name Variables for model simulation setting
!> Information about the model simulation comes from the file info.txt. 
!> \{
  CHARACTER(LEN=maxcharpath) infodir  !<Directory for information about simulation to be run (info.txt, pmsf.txt, hyss.txt)
  CHARACTER(LEN=maxcharpath) modeldir !<Directory for model setup (e.g. GeoData.txt, Tobs.txt,...)
  CHARACTER(LEN=maxcharpath) resdir   !<Directory for result-files (e.g. simass.txt, 10001.txt, ...)
  TYPE(DateType) :: steplen      !<Length of simulation time step, from timestep in info.txt
  TYPE(DateType) :: bdate        !<Begin simulation date
  TYPE(DateType) :: sdate        !<End simulation date
  TYPE(DateType) :: outstartdate !<Date for first output and criteria calculation
  TYPE(DateType),ALLOCATABLE :: psdates(:)    !<Matrix with datetimes of changing point sources
  LOGICAL :: readmatlab          !<Flag for read observations in format suitable for MATLAB to write
  LOGICAL :: readdaily           !<Flag for read forcing data each time step instead of all in the beginning
  LOGICAL :: writematlab         !<Flag for print out in format suitable for MATLAB to read
  LOGICAL :: writeload           !<Flag for print out of yearly loads
  LOGICAL :: simsubmodel         !<Flag for simulation of submodel smaller than basemodel
  LOGICAL :: readobsid           !<Flag for read/use pobsid/tobsid/etc
  LOGICAL :: resultseq           !<Flag for use sequence number on result files, e.g. timeCOUT_007.txt
  LOGICAL :: checkdata(2,2)      !<Flags for checking input data: 1=model set-up, 2=observations, and stopping simulation
  LOGICAL :: vegtypereset        !<Flag for vegetation type is missing and set to 1
  LOGICAL :: readoutregion       !<Flag for reading Outregions.txt
  INTEGER :: simsequence         !<Number of the sequence to be simulated (used for PTobs and possibly result-files).
  INTEGER :: ndt                 !<Number of time steps in simulation
!> \}  

!> \name Forcing data and other observations variables and structures
!> \{
  TYPE(DateType),ALLOCATABLE :: dates(:)    !<Matrix with all datetimes in simulation
  INTEGER,ALLOCATABLE :: qobsindex(:)       !<Index to find correct in qobs-array
  REAL,ALLOCATABLE :: qobs(:,:)             !<Matrix with all observed runoff data (time,stn)
  REAL,ALLOCATABLE :: xobs(:,:)             !<Matrix with all other observation data
  REAL,ALLOCATABLE :: xoregobs(:,:)         !<Matrix with all other outregion observation data
  INTEGER :: numqobsstn                     !<Number of qobs station
  INTEGER :: xcol = 0                       !<Number of columns in xobs, default is zero
  INTEGER :: xorcol = 0                     !<Number of columns in xoregobs, default is zero
  INTEGER :: readformat                     !<File format of obs-files; 0=ASCII, (4=netcdf)
  INTEGER :: timeformat                     !<File format of time; 0=only date, 1=date and time (HHMM)
  TYPE(DateType) :: bqdate                  !<Begin date for Q observations within simulation period
  TYPE(DateType) :: eqdate                  !<End date for Q observations within simulation period
  TYPE(DateType) :: bxdate                  !<Begin date for X observations within simulation period
  TYPE(DateType) :: exdate                  !<End date for X observations within simulation period
  TYPE(DateType) :: bxrdate                 !<Begin date for X outregion observations within simulation period
  TYPE(DateType) :: exrdate                 !<End date for X outregion observations within simulation period
  LOGICAL :: readqobs  = .FALSE.            !<Flag for reading qobs daily (Qobs.txt exist)
  TYPE(FORCINGDATATYPE),ALLOCATABLE :: forcingdata(:)
!> \}
!> \name Constant parameters for output
!> \{
  INTEGER, PARAMETER :: max_typeofoutput = 4 !<Maximum kinds of output; 1=basinoutput,2=mapoutput,3=timeoutput,4=regionoutput
  INTEGER, PARAMETER :: maxoutbasins = 10000 !<Maximum numbers of subbasins for subbasin output
  !Code for accumulation period of output (used for wperiod)
  INTEGER,PARAMETER :: i_t=0  !<Code for accumulation period of output: timesteply
  INTEGER,PARAMETER :: i_d=1  !<Code for accumulation period of output: daily
  INTEGER,PARAMETER :: i_w=2  !<Code for accumulation period of output: weekly
  INTEGER,PARAMETER :: i_m=3  !<Code for accumulation period of output: monthly
  INTEGER,PARAMETER :: i_y=4  !<Code for accumulation period of output: yearly
  INTEGER,PARAMETER :: i_s=5  !<Code for accumulation period of output: simulation period
  INTEGER,PARAMETER :: i_h=6  !<Code for accumulation period of output: hourly (not useable yet)
  INTEGER,PARAMETER :: max_typeofperiods = 6 !<Maximum kinds of output periods
  CHARACTER(LEN=2),PARAMETER :: outperiodname(0:max_typeofperiods) = (/'TS','DD','WK','MO','YR','SP','HR'/)
!> \}

!> \name Variables for output
!> These variables hold information about wanted output, accumulation of output variables 
!> for period output. HYSS handles three kind of output. Basin output gives one file for 
!> each subbasin with time series of selected variables in columns. Time output gives one 
!> file for each selected variable with time series for each subbasin in columns. Map output 
!> gives one file for each selected variable with subbasins in rows and possibly several time periods in columns.
!> \{
  LOGICAL :: outvarallbasin                         !<All subbasins will be printed separately
  TYPE(OUTPUTTYPE),ALLOCATABLE :: output(:)     !<Information about wanted output
  TYPE(OUTREGIONINFOTYPE),ALLOCATABLE :: outregion(:) !<Information about output regions (noutreg)
  TYPE(OUTVARINFOTYPE),ALLOCATABLE :: outvarinfo(:) !<variable to hold information about output for print out in outvar
  TYPE(OUTVARINFOTYPE),ALLOCATABLE :: outvarinfotemp(:) !<temporary variable to hold output for print out in outvar
  INTEGER :: noutput                                !<Number of output 
  INTEGER :: noutreg                                !<Number of output regions
  REAL,ALLOCATABLE    :: accdata_classload(:,:,:,:)!<Accumulated data for yearly loads (class dependent)
  REAL,ALLOCATABLE    :: accdata_basinload(:,:,:)  !<Accumulated data for yearly loads (not class dependent)
  CHARACTER(LEN=20),ALLOCATABLE :: maptime(:)      !<Serie of time for map data
  INTEGER :: tmap                                  !<Current period index for mapdata
  INTEGER :: dtskip                   !Number of timesteps skipped in criteria/output calculation     
  INTEGER :: maxmap                   !Numbers of timesteps for mapping (for array-dimension)
  INTEGER :: idtlag                   !Number of timesteps between 00:00 and first timestep
!> \}

!> \name Variables for optimisation
!> A number of variables are used during calibration. These include information about what calibration to do, 
!> including calibration parameter settings, model parameter space to optimize, and index-variables 
!> to locate them. Also result of on-going optimization, control variables for optimization progress. 
!> \{
  INTEGER :: dimpar                  !<Maximum number of parameter values per parameter
  INTEGER :: numoptimpar             !<Effective amount of optimization parameters
  INTEGER :: optimFuncCall           !<Total amount of model runs during calibration
  INTEGER :: lineSearchCallCount     !<Total amount of line search calls
  INTEGER :: calvarper               !<Average period for criteria
  INTEGER :: calvarlim               !<Minimum number of data for criteria to be calculated
  INTEGER :: ncrit                   !<Number of criteria to be included in calibration
  INTEGER :: nacrit                  !<Number of criteria with unique variables
  INTEGER :: optparid(maxoptpar)     !<Index of optimisation parameter in modparid
  REAL    :: optimStartTime          !<CPU time at calibration start
  LOGICAL :: doopt = .FALSE.         !<Calibration flag
  TYPE(CALIBRATIONINFO), ALLOCATABLE :: calvar(:)    !<Calibration information
  TYPE(ACCUMULATIONINFO), ALLOCATABLE :: acccalvar(:)!<Calibration information for accumulation
  REAL,ALLOCATABLE :: optparmin(:,:)       !<Lower limit of parameter space
  REAL,ALLOCATABLE :: optparmax(:,:)       !<Upper limit of parameter space
  REAL,ALLOCATABLE :: optparprecision(:,:) !<Decimal precision up to which parameter has to be calibrated
  INTEGER,ALLOCATABLE :: parindex(:,:)     !<Parameter index and subbasin/landuse/soiltype
  TYPE(OPTIMIZATIONTYPE) :: optim          !<Variable defining optimisation tasks and number of simulations

  !Variables for MonteCarlo simulation
  REAL, ALLOCATABLE :: bestMCparameters(:,:)   !<Parameter values for the best MonteCarlo simulations (so far)
  REAL, ALLOCATABLE :: bestMCoptcrit(:)        !<Optcrit values for the best MonteCarlo simulations (so far)
  REAL, ALLOCATABLE :: bestMCperformance(:,:,:)!<Performance values for the best MonteCarlo simulations (so far)
  REAL, ALLOCATABLE :: bestMCcondcrit(:)       !<Conditional crit values for the best MonteCarlo simulations (so far)
!> \}

!> \name Criteria calculation variables
!> Several criteria may be calculated during a simulation, but only one 
!> (often a combination of different criteria) is used for optimization. 
!> Criteria are most often calculated from variables accumulated during 
!> the simulation, but some need all values present when the criterion is calculated.
!> \{
  INTEGER nsubCrit                     !<Number of subbasins+outregions used for subbasin criteria calculations
  DOUBLE PRECISION, ALLOCATABLE :: critvec(:,:,:) !<Accumulated errors for later criteria calculation; variable,type,basin
  REAL, ALLOCATABLE    :: rs(:,:)      !<Variable for accumulate recorded values for periodmean
  REAL, ALLOCATABLE    :: cs(:,:)      !<Variable for accumulate observed values for periodmean
  INTEGER, ALLOCATABLE :: ts(:,:)      !<Variable for accumulate number of values for periodmean
  REAL, ALLOCATABLE    :: ktcomp(:,:)  !<Variable for computed values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktcomp2(:,:) !<Variable for computed values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktcomp3(:,:) !<Variable for computed values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktcomp4(:,:) !<Variable for computed values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktrec(:,:)   !<Variable for recorded values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktrec2(:,:)  !<Variable for recorded values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktrec3(:,:)  !<Variable for recorded values for Kendalls Tau or RA calculation
  REAL, ALLOCATABLE    :: ktrec4(:,:)  !<Variable for recorded values for Kendalls Tau or RA calculation
  INTEGER, ALLOCATABLE :: ktnum(:)     !<Variable for number of pairs for Kendalls Tau or RA calculation
  INTEGER, ALLOCATABLE :: ktnum2(:)    !<Variable for number of pairs for Kendalls Tau or RA calculation
  INTEGER, ALLOCATABLE :: ktnum3(:)    !<Variable for number of pairs for Kendalls Tau or RA calculation
  INTEGER, ALLOCATABLE :: ktnum4(:)    !<Variable for number of pairs for Kendalls Tau or RA calculation
!> \}

!> \name Variables for Data Assimilation
!> Most variables needed for the assimilation-routines are defined in the assimilation-modules.
!> However, one basic flag called doassimilation is declared here, and read from info.txt by data.f90 (assimilation y/n).
!> \{  
  LOGICAL :: doassimilation = .FALSE.        !<Data Assimilation on/off flag
!> \}

!Subroutines and functions
CONTAINS

  !---------------------------------------------------------------------
  !>\brief Allocate variables for accumulation of data for printout 
  !> for data assimilation simulation
  !>
  !> \b Consequences Module variables output may be reallocated
  !---------------------------------------------------------------------
  SUBROUTINE da_allocate_accumulation(n,nsta)

    !Argument declarations
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    INTEGER, INTENT(IN) :: nsta   !<number of statistics output

    !Local variables
    INTEGER io,na
    
    DO io=1,noutput
      IF(output(io)%narea>0)THEN
        na=output(io)%narea
      ELSE
        na=n
      ENDIF
      IF(output(io)%fileformat==1.OR.output(io)%fileformat==3.OR.output(io)%fileformat==4)THEN
        IF(.NOT.ALLOCATED(output(io)%accdata%value)) ALLOCATE(output(io)%accdata%value(na,output(io)%nvar,nsta+1))
        IF(.NOT.ALLOCATED(output(io)%accdata%help)) ALLOCATE(output(io)%accdata%help(na,output(io)%nvar,nsta+1))
        IF(.NOT.ALLOCATED(output(io)%accdata%nok)) ALLOCATE(output(io)%accdata%nok(na,output(io)%nvar,nsta+1))
      ELSEIF(output(io)%fileformat==2)THEN
        IF(.NOT.ALLOCATED(output(io)%accdata%value)) ALLOCATE(output(io)%accdata%value(na,output(io)%nvar,maxmap))
        IF(.NOT.ALLOCATED(output(io)%accdata%help)) ALLOCATE(output(io)%accdata%help(na,output(io)%nvar,maxmap))
        IF(.NOT.ALLOCATED(output(io)%accdata%nok)) ALLOCATE(output(io)%accdata%nok(na,output(io)%nvar,maxmap))
        IF(.NOT.ALLOCATED(maptime)) ALLOCATE(maptime(maxmap))
        WRITE(6,*) 'WARNING: Mapfiles not saved for data assimilation run statistics'
      ENDIF
    ENDDO
    
    IF(writeload)THEN
      WRITE(6,*) 'WARNING: Loadfiles not saved for data assimilation run statistics'
      STOP 1
    ENDIF

  END SUBROUTINE da_allocate_accumulation

  !---------------------------------------------------------------------
  !>\brief Allocate variables for accumulation of data for printout 
  !>
  !> \b Consequences Module variables maptime,output and/or 
  !>accdata_classload,accdata_basinload may be allocated
  !---------------------------------------------------------------------
  SUBROUTINE allocate_accumulation(n,nj,ns,nc,nb)

    !Argument declarations
    INTEGER, INTENT(IN) :: n      !<number of subbasins
    INTEGER, INTENT(IN) :: nj     !<number of classes
    INTEGER, INTENT(IN) :: ns     !<number of substances
    INTEGER, INTENT(IN) :: nc     !<maximum class dependent loads
    INTEGER, INTENT(IN) :: nb     !<maximum non-class dependent loads

    !Local variables
    INTEGER io,na

    !Allocate space for accumulation of output to output periods
    DO io = 1,noutput
      IF(output(io)%narea>0)THEN
        na=output(io)%narea
      ELSE
        na=n
      ENDIF
      IF(output(io)%fileformat==1.OR.output(io)%fileformat==3.OR.output(io)%fileformat==4)THEN
        IF(.NOT.ALLOCATED(output(io)%accdata%value)) ALLOCATE(output(io)%accdata%value(na,output(io)%nvar,1))
        IF(.NOT.ALLOCATED(output(io)%accdata%help)) ALLOCATE(output(io)%accdata%help(na,output(io)%nvar,1))
        IF(.NOT.ALLOCATED(output(io)%accdata%nok)) ALLOCATE(output(io)%accdata%nok(na,output(io)%nvar,1))
      ELSEIF(output(io)%fileformat==2)THEN
        IF(.NOT.ALLOCATED(output(io)%accdata%value)) ALLOCATE(output(io)%accdata%value(na,output(io)%nvar,maxmap))
        IF(.NOT.ALLOCATED(output(io)%accdata%help)) ALLOCATE(output(io)%accdata%help(na,output(io)%nvar,maxmap))
        IF(.NOT.ALLOCATED(output(io)%accdata%nok)) ALLOCATE(output(io)%accdata%nok(na,output(io)%nvar,maxmap))
        IF(.NOT.ALLOCATED(maptime)) ALLOCATE(maptime(maxmap))
      ENDIF
    ENDDO
    
    !Allocate space for accumulation of nutrient load output
    IF(writeload)THEN
      IF(ns>0)THEN
        IF(.NOT.ALLOCATED(accdata_classload)) ALLOCATE(accdata_classload(nj,nc,ns,n))
        IF(.NOT.ALLOCATED(accdata_basinload)) ALLOCATE(accdata_basinload(ns,nb,n))
      ENDIF
    ENDIF

  END SUBROUTINE allocate_accumulation

  !>Reallocate output information to suitable size after reading actual output
  !>
  !> \b Consequences Module variables outvarinfo may be reallocated.
  !----------------------------------------------------------------
  SUBROUTINE reallocate_outvar_information(nnew)

    INTEGER, INTENT(IN) :: nnew   !<number of output variables to be calculated

    INTEGER i
    
    IF(.NOT.ALLOCATED(outvarinfo)) ALLOCATE(outvarinfo(nnew))
    DO i = 1,nnew
      outvarinfo(i)=outvarinfotemp(I)
    ENDDO
    DEALLOCATE(outvarinfotemp)

  END SUBROUTINE reallocate_outvar_information

  !>Collects current date from array with dates      
  !----------------------------------------------------------------
  TYPE(DateType) FUNCTION get_current_date_memory(i)

    INTEGER, INTENT(IN) :: i          !<current time step
    !< \retval get_current_date_memory current date
 
    get_current_date_memory = dates(i)

  END FUNCTION get_current_date_memory

  !>Collects current discharge observation
  !--------------------------------------------------------
  FUNCTION get_current_qobs(i,n)

    INTEGER, INTENT(IN) :: i    !<current time step
    INTEGER, INTENT(IN) :: n    !<number of subbasins
    REAL :: get_current_qobs(n) !< \retval get_current_qobs current discharge
    
    !Local variables
    REAL x(n)

    x = missing_value
    IF(ALLOCATED(qobs))THEN
       WHERE(qobsindex(1:n)>0)
          x=qobs(i,qobsindex(1:n))
       ENDWHERE
    ENDIF
    get_current_qobs = x

  END FUNCTION get_current_qobs

  !>Collects current other observation
  !-----------------------------------------
  FUNCTION get_current_xobs(i)

    INTEGER, INTENT(IN) :: i    !< current time step
    REAL get_current_xobs(xcol) !< \retval get_current_xobs current value of other observations
    
    !Local variables
    REAL x(xcol)

    IF(xcol==0) RETURN        !no observations
    x(:)=xobs(i,:)
    get_current_xobs = x

  END FUNCTION get_current_xobs

  !>Collects current other observation
  !-----------------------------------------
  FUNCTION get_current_xoregobs(i)

    INTEGER, INTENT(IN) :: i    !< current time step
    REAL get_current_xoregobs(xorcol) !< \retval get_current_xoregobs current value of other observations
    
    !Local variables
    REAL x(xorcol)

    IF(xorcol==0) RETURN        !no observations
    x(:)=xoregobs(i,:)
    get_current_xoregobs = x

  END FUNCTION get_current_xoregobs

  !>Collects current forcing data from memory
  !--------------------------------------------------------
  FUNCTION get_current_forcing_from_memory(i,obsdata) RESULT(current_value)
  
    INTEGER, INTENT(IN) :: i      !<current time step
    TYPE(FORCINGDATATYPE), INTENT(IN) :: obsdata
    REAL current_value(obsdata%ncols)         !< \retval current forcing data
  
    current_value = obsdata%allvalues(i,1:obsdata%ncols)   
  
  END FUNCTION get_current_forcing_from_memory

  !>Saves information on what and how to write to output files
  !>(except index to locate the values in outvar)
  !>
  !> \b Consequences Module variable output is allocated and set
  !------------------------------------------------------------------------
  INTEGER FUNCTION set_output(itype,iout,maxout,dimout,outv,n,per,ndec,nsig,aggarea,maxa,na,area)

    INTEGER, INTENT(IN) :: itype        !<type of output file
    INTEGER, INTENT(IN) :: iout         !<current output file type to be set
    INTEGER, INTENT(IN) :: maxout       !<number of output file types for this simulation 
    INTEGER, INTENT(IN) :: dimout       !<dimension of variable arrays
    INTEGER, INTENT(IN) :: outv(dimout) !<variables for output
    INTEGER, INTENT(IN) :: n            !<number of output variables of this type
    INTEGER, INTENT(IN) :: per          !<code for print out period
    INTEGER, INTENT(IN) :: ndec         !<number of decimals
    INTEGER, INTENT(IN) :: nsig         !<number of significant figures
    INTEGER, INTENT(IN),OPTIONAL :: aggarea(dimout)  !<area aggregation type of variable (uubbasin=0,upstream=1,region=2)
    INTEGER, INTENT(IN),OPTIONAL :: maxa         !<dimension of area array
    INTEGER, INTENT(IN),OPTIONAL :: na           !<number of areas
    INTEGER, INTENT(IN),OPTIONAL :: area(:)   !<id of areas
    !< \retval set_output error status of function

    !Local variables 
    INTEGER i
    
    !> \b Algorithm \n
    set_output  = 0

    !>If not allocated: allocate and initialize output
    IF(.NOT.ALLOCATED(output)) ALLOCATE(output(maxout))
    ALLOCATE(output(iout)%variable(n))
    
    !>Set outvar information about variables to be written
    DO i = 1,n
      output(iout)%variable(i)%idindex     = outv(i)
      IF(itype==4)THEN
        output(iout)%variable(i)%areaagg   = 2 !Only region output in region output files
      ELSEIF(PRESENT(aggarea))THEN
        output(iout)%variable(i)%areaagg   = aggarea(i) ! (up->1, rg->2, other 0)
      ELSE
        output(iout)%variable(i)%areaagg   = 0
      ENDIF
    ENDDO
    !>Set output information about the file format
    output(iout)%fileformat = itype
    output(iout)%nvar       = n
    output(iout)%period     = per
    output(iout)%decimal    = ndec
    IF(nsig>0)THEN
      output(iout)%signfig  = nsig-1
    ELSE
      output(iout)%signfig  = 0
    ENDIF
    !>Set output information about subbasins or regions to be written (if applicable)
    IF(PRESENT(maxa))THEN
      output(iout)%narea    = na
      ALLOCATE(output(iout)%areaindex(na))
      DO i = 1,na
        output(iout)%areaindex(i) = area(i)   !this is now subid or outregionid, will be changed to index after submodelis read
      ENDDO
    ELSE
      output(iout)%narea    = 0
    ENDIF

  END FUNCTION set_output

  !>Saves information of calibration and criteria to use
  !>
  !> \b Consequences Module variables doopt,calvarper,calvarlim,calvar and nacrit are set.
  !> Module variable acccalvar may be allocated and set.
  !------------------------------------------------------------------------
  INTEGER FUNCTION set_calibration(a,n,c,per,lim,w,cvar,rvar,typ,par,cond,thres,aagg)

    USE MODVAR, ONLY : changecritvar
    
    LOGICAL, INTENT(IN) :: a               !<calibration
    INTEGER, INTENT(IN) :: n               !<number of criteria
    CHARACTER(LEN=*), INTENT(IN) :: c(n)   !<name of criteria
    INTEGER, INTENT(IN) :: per             !<code for accumulation period
    INTEGER, INTENT(IN) :: lim             !<limit
    REAL, INTENT(IN)    :: w(n)            !<weight of criteria
    INTEGER, INTENT(IN) :: cvar(n)         !<computed variable(s)
    INTEGER, INTENT(IN) :: rvar(n)         !<recorded variable(s)
    INTEGER, INTENT(IN) :: typ(n)          !<type of variable(s) (state,flow,conc)
    REAL, INTENT(IN)    :: par(n)          !<coefficient of criteria
    LOGICAL, INTENT(IN) :: cond(n)         !<flag for conditional criteria
    REAL, INTENT(IN)    :: thres(n)        !<threshold for conditional criteria
    INTEGER, INTENT(IN) :: aagg(n)         !<areaaggregation of variable (basin,upstream,region)
    !< \retval set_calibration error status of function
    
    !Local variables
    LOGICAL up,skip,tauorra(n)
    INTEGER i,j,k,d
    CHARACTER(LEN=3) c2

    !> \b Algorithm \n
    !>Set general calibration information variables
    set_calibration = 0
    tauorra = .FALSE.
    IF(.NOT.ALLOCATED(calvar)) ALLOCATE(calvar(n))
    IF(a)THEN
      doopt = .TRUE.
    ENDIF
    calvarper = per
    calvarlim = lim

    !>For every criteria in optimization function: Set calibration information variables
    DO j = 1,n
      c2=ADJUSTL(c(j))
      DO i = 1, LEN_TRIM(c2)
        up = (c2(i:i) .GE. 'a') .AND. (c2(i:i) .LE. 'z')
        IF(up)THEN
          c2(i:i) = CHAR(ICHAR(c2(i:i)) - 32)             !Change to capitals
        ENDIF
      ENDDO
      calvar(j)%crit    = c2
      calvar(j)%rec     = rvar(j)
      calvar(j)%comp    = cvar(j)
       
      !Special case for lake water stage in HYPE
      !Change critera variables to water stage outside w-ref system for better numerical performace
      IF(ALLOCATED(changecritvar))THEN
        IF(calvar(j)%comp==changecritvar(1,1) .AND. calvar(j)%rec==changecritvar(2,1))THEN
          calvar(j)%comp  = changecritvar(1,2)
          calvar(j)%rec   = changecritvar(2,2)
        ENDIF
      ENDIF
       
      calvar(j)%weight  = w(j)
      calvar(j)%vartype = typ(j)
      calvar(j)%areaagg = aagg(j)
      calvar(j)%coeff   = par(j)
      IF(calvar(j)%crit=='TAU'.OR.calvar(j)%crit=='RRA'.OR.  &
         calvar(j)%crit=='MRA'.OR.calvar(j)%crit=='SRA'.OR.  &
         calvar(j)%crit=='MDA') tauorra(j) = .TRUE.
      !Elements for conditional acceptance
      IF(cond(j))THEN
        calvar(j)%cond = .TRUE.
        calvar(j)%thres = thres(j)
      ELSE
        calvar(j)%cond = .FALSE.
        calvar(j)%thres = 0.
      ENDIF
    ENDDO

    !>Count number of unique variables in optimization function
    !Temporary for allocating acccalvar
    k = 1
    DO j = 1,n
      skip = .FALSE.
      DO i = 1,j-1
        IF(calvar(j)%rec.EQ.calvar(i)%rec.AND.calvar(j)%comp.EQ.calvar(i)%comp.AND.calvar(j)%areaagg.EQ.calvar(i)%areaagg)THEN
          skip = .TRUE.
        ENDIF
      ENDDO
      IF(.NOT.skip)THEN
        k = k + 1
      ENDIF
    ENDDO
    k = k - 1
    nacrit = k

    !>For every unique variable: 
    !Set information about accumulation for calibration criteria calculation
    IF(.NOT.ALLOCATED(acccalvar)) ALLOCATE(acccalvar(nacrit))
    k = 1
    DO j = 1,n
      d = 0
      DO i = 1,k-1
        IF(calvar(j)%rec.EQ.acccalvar(i)%rec.AND.calvar(j)%comp.EQ.acccalvar(i)%comp.AND.calvar(j)%areaagg.EQ.acccalvar(i)%areaagg)THEN
          d = i
        ENDIF
      ENDDO
      IF(d==0)THEN
        acccalvar(k)%rec     = calvar(j)%rec
        acccalvar(k)%comp    = calvar(j)%comp
        acccalvar(k)%recoutvar  = 0
        acccalvar(k)%compoutvar = 0
        acccalvar(k)%vartype = calvar(j)%vartype
        acccalvar(k)%areaagg = calvar(j)%areaagg
        acccalvar(k)%saveend = tauorra(j)
        k = k + 1
      ELSE
        IF(tauorra(j)) acccalvar(d)%saveend = .TRUE.
      ENDIF
    ENDDO
    
    !>Check for tau or RA criteria order placement
    DO i = 4,nacrit
      IF(acccalvar(i)%saveend)THEN
        WRITE(6,*) 'ERROR: Only the first four optimisation criteria may be TAU or RA.'
        WRITE(6,*) 'ERROR: Change your criteria (order) in info.txt'
        STOP 1
      ENDIF
    ENDDO

  END FUNCTION set_calibration

  !>Saves information on number of times to write for map file
  !>
  !> \b Consequences Module variable maxmap is set.
  !------------------------------------------------------------------
  INTEGER FUNCTION set_maxmap(tstep,per)

    INTEGER, INTENT(IN) :: tstep !<number of timestep in output period
    INTEGER, INTENT(IN) :: per   !<code for accumulation period
    !< \retval set_maxmap error status of function
    
    !Local variables
    INTEGER num_ts_per_day
    set_maxmap  = 0
    num_ts_per_day = seconds_per_day/seconds_per_timestep

    IF(per==i_t)THEN
      maxmap = tstep
    ELSEIF(per==i_d)THEN
      maxmap = tstep/num_ts_per_day + 2    !safe for start AND end time not beginning/end of day
    ELSEIF(per==i_w)THEN
      maxmap = tstep/(7*num_ts_per_day) + 2
    ELSEIF(per==i_m)THEN
      maxmap = tstep/(28*num_ts_per_day) + 2
    ELSEIF(per==i_y)THEN
      maxmap = tstep/(365*num_ts_per_day) + 2
    ELSEIF(per==i_s)THEN
      maxmap = 1
    ENDIF

  END FUNCTION set_maxmap

  !>Deallocate worldvar-arrays
  !>
  !> \b Consequences A lot of module variables are deallocated
  !-----------------------------------------------------------
  SUBROUTINE deallocate_worldvar()

    INTEGER iforc 
    
    IF(ALLOCATED(dates))        DEALLOCATE(dates)
    IF(ALLOCATED(qobs))         DEALLOCATE(qobs)
    IF(ALLOCATED(xobs))         DEALLOCATE(xobs)
    IF(ALLOCATED(qobsindex))    DEALLOCATE(qobsindex)
    IF(ALLOCATED(calvar))       DEALLOCATE(calvar)
    IF(ALLOCATED(critvec))      DEALLOCATE(critvec)
    IF(ALLOCATED(rs))           DEALLOCATE(rs)
    IF(ALLOCATED(cs))           DEALLOCATE(cs)
    IF(ALLOCATED(ts))           DEALLOCATE(ts)
    IF(ALLOCATED(maptime))      DEALLOCATE(maptime)
    IF(ALLOCATED(optparmin))    DEALLOCATE(optparmin)
    IF(ALLOCATED(optparmax))    DEALLOCATE(optparmax)
    IF(ALLOCATED(parindex))     DEALLOCATE(parindex)
    DO iforc = 1,max_forcingdata
      IF(ALLOCATED(forcingdata(iforc)%allvalues)) DEALLOCATE(forcingdata(iforc)%allvalues)
      IF(ALLOCATED(forcingdata(iforc)%basinindex)) DEALLOCATE(forcingdata(iforc)%basinindex)
    ENDDO

  END SUBROUTINE deallocate_worldvar

  !>Allocate and initiate variables for MonteCarlo simulation
  !>
  !> \b Consequences Module variables bestMCparameters,bestMCoptcrit, 
  !> bestMCperformance and bestMCcondcrit are allocated and set.
  !--------------------------------------------------------------------
  SUBROUTINE allocate_MCvariables(nbest,npar,nperf,numcrit)

    INTEGER, INTENT(IN) :: nbest   !<optimruns_best, number of ensambles saved
    INTEGER, INTENT(IN) :: npar    !<number of parameters
    INTEGER, INTENT(IN) :: nperf   !<number of performance criteria (=maxperf)
    INTEGER, INTENT(IN) :: numcrit !<ncrit, number of variables criteria are calculated for

    IF(.NOT.ALLOCATED(bestMCparameters))  ALLOCATE(bestMCparameters(nbest,npar))
    IF(.NOT.ALLOCATED(bestMCoptcrit))     ALLOCATE(bestMCoptcrit(nbest))
    IF(.NOT.ALLOCATED(bestMCperformance)) ALLOCATE(bestMCperformance(nbest,nperf,numcrit))
    IF(.NOT.ALLOCATED(bestMCcondcrit))       ALLOCATE(bestMCcondcrit(nbest))
    bestMCparameters = 0.
    bestMCoptcrit = 99999999.
    bestMCperformance = 0.
    bestMCcondcrit = 99999999.
    
  END SUBROUTINE allocate_MCvariables

  !>Deallocate variables for MonteCarlo simulation
  !>
  !> \b Consequences Module variables are deallocated
  !-------------------------------------------------------------
  SUBROUTINE deallocate_MCvariables()

    IF(ALLOCATED(bestMCparameters))  DEALLOCATE(bestMCparameters)
    IF(ALLOCATED(bestMCoptcrit))     DEALLOCATE(bestMCoptcrit)
    IF(ALLOCATED(bestMCperformance)) DEALLOCATE(bestMCperformance)
    IF(ALLOCATED(bestMCcondcrit))    DEALLOCATE(bestMCcondcrit)

  END SUBROUTINE deallocate_MCvariables

  !>Set time format used in model set-up (date or date-time)
  !>
  !> \b Consequences Module variable timeformat are set
  !-------------------------------------------------------------
  SUBROUTINE set_timeformat(onlydate)

    LOGICAL, INTENT(IN) :: onlydate   !<Flag for time format with no time, only date

    timeformat = 0 !default, only date
    IF(.NOT.onlydate) timeformat=1 !date and time

  END SUBROUTINE set_timeformat

  !> Find a free file unit for file connecting
  !------------------------------------------------------------
  INTEGER FUNCTION fileunit_get()

    !Local variables
    INTEGER i
    
    DO i = 1,max_files
      IF(.NOT.fileunits(i))THEN
        fileunit_get = i+fileunit_base
        fileunits(i) = .TRUE.
        RETURN
      ENDIF
    ENDDO
    WRITE(6,*) 'Error: No free file unit found. To many files open'
    STOP 1

  END FUNCTION fileunit_get

  !> Release file unit no longer used for file connection
  !------------------------------------------------------------
  SUBROUTINE fileunit_free(funit)

    !Argument declaration
    INTEGER, INTENT(IN) :: funit    !<File unit to be released
    
    fileunits(funit-fileunit_base) = .FALSE.

  END SUBROUTINE fileunit_free

  !> Set sequence filename
  !------------------------------------------------------------
  SUBROUTINE get_seq_filename(fname)

    !Argument declaration
    CHARACTER(LEN=*), INTENT(INOUT) :: fname    !<File name
    
    !Local variable
    INTEGER l
    CHARACTER(LEN=3) seqnr
    CHARACTER(LEN=LEN(fname)) filename1
    
    l=LEN_TRIM(fname)
    filename1 = fname   !Determine first part of filename
    IF(l>4)THEN
      IF(fname(l-3:l)=='.txt')THEN
        filename1 = fname(1:l-4)
      ENDIF
    ENDIF
    fname = ''    !Set filename
    IF(simsequence>0)THEN
      WRITE(seqnr,'(I3.3)') simsequence
      fname = TRIM(filename1)//'_'//seqnr//'.txt'
    ELSE
      fname = TRIM(filename1)//'.txt'
    ENDIF

  END SUBROUTINE get_seq_filename
  
!>\brief This routine determines the value of numoptimpar, i.e.
!! the number of model parameter values to be optimized
!-----------------------------------------------------------------------
SUBROUTINE count_optim_par(dim1,dim2)

  IMPLICIT NONE
  
  !Argument declarations
  INTEGER, INTENT(IN) :: dim1,dim2      !Dimension of optparmin/max-variables
  
  !Local variables
  INTEGER i,j

  numoptimpar = 0
  DO i = 1, dim1
    DO j = 1, dim2
      IF(optparmin(i,j) .NE. optparmax(i,j)) THEN
        numoptimpar = numoptimpar + 1
      ENDIF
    ENDDO
  ENDDO

  END SUBROUTINE count_optim_par

  !>Allocate and set forcingdata structure variable    
  !>
  !> \b Consequences Module variables for model options are allocated.
  !--------------------------------------------------------------
  SUBROUTINE allocate_and_initialize_forcingdata_structure()

    ALLOCATE(forcingdata(max_forcingdata))
    forcingdata(i_pobs)%readfile = .TRUE.
    forcingdata(i_pobs)%filename = 'Pobs'
    forcingdata(i_pobs)%infocode = ''
    forcingdata(i_pobs)%idcode = 'pobsid    '
    forcingdata(i_tobs)%readfile = .TRUE.
    forcingdata(i_tobs)%filename = 'Tobs'
    forcingdata(i_tobs)%infocode = ''
    forcingdata(i_tobs)%idcode = 'tobsid    '
    forcingdata(i_tminobs)%readfile = .FALSE.
    forcingdata(i_tminobs)%filename = 'TMINobs'
    forcingdata(i_tminobs)%infocode = 'readtminobs'
    forcingdata(i_tminobs)%idcode = 'tminobsid '
    forcingdata(i_tmaxobs)%readfile = .FALSE.
    forcingdata(i_tmaxobs)%filename = 'TMAXobs'
    forcingdata(i_tmaxobs)%infocode = 'readtmaxobs'
    forcingdata(i_tmaxobs)%idcode = 'tmaxobsid '
    forcingdata(i_rhobs)%readfile = .FALSE.
    forcingdata(i_rhobs)%filename = 'RHobs'
    forcingdata(i_rhobs)%infocode = 'readrhobs'
    forcingdata(i_rhobs)%idcode = 'rhobsid   '
    forcingdata(i_sfobs)%readfile = .FALSE.
    forcingdata(i_sfobs)%filename = 'SFobs'
    forcingdata(i_sfobs)%infocode = 'readsfobs'
    forcingdata(i_sfobs)%idcode = 'sfobsid   '
    forcingdata(i_swobs)%readfile = .FALSE.
    forcingdata(i_swobs)%filename = 'SWobs'
    forcingdata(i_swobs)%infocode = 'readswobs'
    forcingdata(i_swobs)%idcode = 'swobsid   '
    forcingdata(i_uobs)%readfile = .FALSE.
    forcingdata(i_uobs)%filename = 'Uobs'
    forcingdata(i_uobs)%infocode = 'readuobs'
    forcingdata(i_uobs)%idcode = 'uobsid    '

  END SUBROUTINE allocate_and_initialize_forcingdata_structure


END MODULE WORLDVAR

