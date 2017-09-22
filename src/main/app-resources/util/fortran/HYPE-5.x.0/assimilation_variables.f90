!> \file assimilation_variables.f90
!> Contains module assimilation_variables, with variables needed to run the (EnKF) data assimilation routines

!  Author: D.Gustafsson (SMHI)
!  Versions: 
!  2012.05.10 Original adapted to HYPE ver 3.6
!  2013.09.19 Adapted to HYPE ver 4.5.x
!  2013.12.03 Renamed to assim_data.f90, to open up for more DA methods in addition to EnKF
!  2014.11.09 Reformated these comments to be read  by doxygen
!  2014.12.10 New data for saving ensemble data to direct access binary files instead of keeping in memory.
!  2015.06.22 Clean-up version (not finalized)
!  2015.09.18 Quick and dirty clean-up for the HYPE course 2015-09-25
!  2016.10.20 File and module name changed to assimilation_variables(.f90), plus adaption to HYPE coding style (FORTRAN names in UPPER CASE and variable/routine names in lower case.
  
!> Variables needed to run the (EnKF) data assimilation routines within HYSS

MODULE ASSIMILATION_VARIABLES
!Copyright 2016-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------------------

USE random_routines, ONLY: randxy_data

IMPLICIT NONE
!---------------------------------------------------------------------------------------------------------------------------
! ASSIMILATION APPLICATION INFO VARIABLE (application specific settings)
!---------------------------------------------------------------------------------------------------------------------------

!Type declarations
!>Type for assimilation status of different categories/variables.
TYPE assim_flag_type
  LOGICAL           :: add           !<assimilation status
  CHARACTER(LEN=50) :: category = '' !<name of category or statetype
  CHARACTER(LEN=50) :: varname  = '' !<name of variable, empty string means apply to whole category
END TYPE assim_flag_type

!>Assimilation application information type (application specific settings).
TYPE assim_info_type
  
  !Model structure variables
  INTEGER nX           !<number of model state variables      (dimension of assimData%X vector)
  INTEGER nF           !<number of model forcing variables    (dimension of assimData%F vector)
  INTEGER nA           !<number of model outvar variables     (dimension of assimData%A vector)
  INTEGER nA2          !<number of model outvarbasin variables (dimension of assimData%A2 vector) !not used in HYPE maybe in HOPE?
  INTEGER nObs         !<number of observation variables      (dimension of assimData%Obs vector)
  INTEGER nD           !<number of observations (rows) in D and HX ensemble matrices for the current Analysis
  !INTEGER nP           !number of parameters in P ensemble  (rows in the P ensemble)
  !INTEGER nDT          !number of timesteps in assimilation time window
  
  INTEGER nLoc         !<number of localization matrices     (depend on the number of spatial domains represented by model state variables, 
                       !<                                     in HYPE there are several potential domains: subbasins, aquifers, parameter regions,
                       !<                                     as well as non-spatial states such as general and land use parameters)
  INTEGER ncoord       !<number of spatial domains           (length of EnkfCoordinates(:) vector)
  
  !ENKF settings
  INTEGER nE           !<number of ensemble members
  !INTEGER locFunc      !localization function (0 none, 1 exponential (xy, z) !CP170615 not used
  
  !LOGICAL flags, may be modified...
  LOGICAL FA           !<include auxiliaries in kalman filter (general switch on/off)
  LOGICAL FP           !<include parameters  in kalman filter (general switch on/off) 
  LOGICAL FF           !<include forcing     in kalman filter (general switch on/off)
  LOGICAL meanout      !<print ensemble mean (.true.) or median (.false.) in output files
  
  !assimilation flags for different variable categories
  TYPE(assim_flag_type),ALLOCATABLE :: assim_flag(:) !<assimilation status of different categories/variables.
  INTEGER :: nFlag  !<number of data in assim_flag structure
  INTEGER :: nCat   !<size of assim_categories (number of possible state variables)
  CHARACTER(LEN=20), ALLOCATABLE :: assim_categories(:)
    
  !some parameters
  !REAL :: moradkhani_delta ! coef.[0-1] to retain variance in parameter ensemble (Moradkhani et al, 2004), value typical around 0.95 !CP170615 not used
  
  !random number generation seed
  !INTEGER :: seed  !CP not used
  
  !missing value
  REAL :: missing !<value of missing data 
  
  !minimum allowed sigma in random number generation (general)
  !REAL :: ensgen_minsigma
  
  !option to read/write ensemble data to direct access binary files instead of saving in allocated memory
  !LOGICAL useBinFiles
  INTEGER useBinFilesX   !<flag for X; 0=save in memory, 1=save to one bin-fil, 2=save to several bin-files
  INTEGER useBinFilesFA  !<flag for F and A; 0=save in memory, 1=save to one bin-fil, 2=save to several bin-files
  INTEGER nBinFiles      !<number of direct access binary files to be used
  
  !options regarding outputs
  INTEGER :: nstatout   !<number of statistical extra outputs of simulation (e.g. 2 means min and max)
  
  !Localization parameters
  REAL xy_scalefac,z_scalefac
  
  !option to collapse ensembles of non-controlvariables or not
  LOGICAL :: collapseNonControlled
  
  !option to initialize ensembles from existing binfiles and save final state to model resultdir
  LOGICAL :: initializeFromBinFiles

END TYPE assim_info_type

!---------------------------------------------------------------------------------------------------------------------------
!>Type for holding ensemble data in matrix and assisiated variables
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_ensemble_type
  !number of rows, columns, and the ensemble matrix (y)
  INTEGER           :: nvar           !<number of variables, ie n:o "model units" (for instance, number of sub-basins)
  INTEGER           :: nens           !<number of ensemble members
  !ensemble data matrix (x)
  REAL, ALLOCATABLE :: x(:,:)         !<ensemble data [nvar x nens]
  !variables used for min/max checks, etc
  REAL              :: minimum        !<min allowed 
  REAL              :: maximum        !<max allowed 
  !variables for output statistics
  REAL, ALLOCATABLE :: outmean(:), outquant(:,:), outmin(:), outmax(:), outsigma(:)
  !indices used for I/O
  INTEGER           :: fileID         !<file id used for read/write to direct access binary file (if needed)
  INTEGER           :: rec            !<record of binary file, CP161201
  !indices to link to spatial information
  INTEGER           :: coordID     !<spatial coordinates ID
  INTEGER           :: locID       !<localization ID
  !Flag for assimilation (true) or re-initialization(false)
  LOGICAL           :: assimilate  !<flag for assimilation (true) or re-initialization (false) of data
END TYPE assim_ensemble_type

!---------------------------------------------------------------------------------------------------------------------------
!ENSEMBLE ENKF2MODEL INTERFACE DATA TYPE 
!> Type for holding indices used as interface between ENKF data and Model data
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_interface_type
  !Indices used to interface between ENKF data and Model data
  CHARACTER(LEN=30)    :: varName     !<character string for model variables (used for debugging, and possibly file names and outputs)
  INTEGER              :: varID       !<variable ID (id number used by interface for linking to model variables)  (in HYPE it can be an outvar index, or the order in which the state variables are considered by interface)
  INTEGER              :: modID       !<model ID,  link to the corresponding variables used for H(X)              (in HYPE: outvar index)
  INTEGER              :: selID       !<special index, to select 1st or 2nd set of varID:s and modID:s (ie. outvars or outvarbasinoutvars)  !CP not used in HYPE
  INTEGER              :: nSubDim     !<number of sub-dimensions (if needed, for instance lateral sub-units, vertical layers, substances, etc, in HYPE for instance SLC, substances, landuse, or slc, etc)
  INTEGER, ALLOCATABLE :: subDimID(:) !<index in the sub-dimensions
END TYPE assim_interface_type
  
!---------------------------------------------------------------------------------------------------------------------------
!ENSEMBLE GENERATION INFORMATION DATA TYPE 
!>Type for holding data needed for ensemble generation
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_generation_type
  !number of variables
  INTEGER           :: nvar      !<number of variables, ie n:o "model units" (for instance, number of sub-basins)
  INTEGER           :: ndata     !<number of non-missing data
  !standard deviation (sigma) and mean (mean) used in ensemble generation
  REAL, ALLOCATABLE :: sigma(:)  !<standard deviation [nvar x 1]
  REAL, ALLOCATABLE :: mean(:)   !<mean               [nvar x 1]
  !parameter used for ensemble generation
  INTEGER           :: ensgen   !<type of ensemble generation        (0 none, 1 unrestricted, 2 [min,+inf], 3 [-inf,max], 4 restricted [min,max])   
  REAL              :: fixsigma !<fixed standard deviation           (ensgen=1)
  REAL              :: semimeta !<relative sigma for semi-restricted (ensgen=2,3,4, following Turner et al 2008)
  REAL              :: restmeta !<relative sigma for restricted      (ensgen=2,3,4)
  REAL              :: minsigma !<minimum sigma                      (ensgen=2,3,4)
  !LOGICAL           :: ECMWF   !CP170615 not used
  !variable for generation of spatially correlated random data, if needed
  LOGICAL           :: dorandxy  !<flag for generation of spatially correlated random data
  TYPE(randxy_data) :: myrandxy_data  !<data needed for generation of spatially correlated random data
END TYPE assim_generation_type

!---------------------------------------------------------------------------------------------------------------------------
!STATE VARIABLE ENSEMBLE TYPE (one variable, all model units)
!> Type for holding state ensembles (data and interface information)
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_state_ensemble_type
  TYPE(assim_ensemble_type)  :: x  !<ensemble data
  TYPE(assim_interface_type) :: info  !<interface data
END TYPE assim_state_ensemble_type

!---------------------------------------------------------------------------------------------------------------------------
!ENSEMBLE TYPE USED FOR FORCING, OBSERVATION, and PARAMETER ENSEMBLES 
!>Type for holding other ensembles than states (data, interface and generation information)
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_input_ensemble_type
  TYPE(assim_ensemble_type) :: x  !<ensemble data
  TYPE(assim_interface_type) :: info  !<interface data
  TYPE(assim_generation_type) :: gen  !<generation data
END TYPE assim_input_ensemble_type

!---------------------------------------------------------------------------------------------------------------------------
!SPATIAL COORDINATES MATRIX VARIABLES TYPE 
!> Type for spatial coordinates
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_coordinate_type
  INTEGER n                 !<number of model positions
  REAL, ALLOCATABLE :: x(:) !<x-coordinate
  REAL, ALLOCATABLE :: y(:) !<y-coordinate
  REAL, ALLOCATABLE :: z(:) !<z-coordinate
END TYPE assim_coordinate_type

!---------------------------------------------------------------------------------------------------------------------------
!LOCALIZATION MATRIX VARIABLES TYPE !CP161214 needed several locCXY-matrixes for different coordinate systems
!> Type for 2-dimensional matrix, used for localizations
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_twoDmatrix_type
  INTEGER n1,n2               !<dimensions
  REAL, ALLOCATABLE :: x(:,:) !<data matrix
END TYPE assim_twoDmatrix_type

!---------------------------------------------------------------------------------------------------------------------------
!ASSIMILATION DATA TOP LEVEL DATA TYPE
!>Type for holding all assimilation data
!---------------------------------------------------------------------------------------------------------------------------
TYPE assim_data_type
  !general information about the assimilation application
  TYPE(assim_info_type)                           :: info  !<general information about the assimilation application
  !states, forcing, auxiliaries, parameters and observations ensembles
  TYPE(assim_state_ensemble_type), ALLOCATABLE :: X(:)     !<model state ensemble (vector over variables)
  TYPE(assim_input_ensemble_type), ALLOCATABLE :: F(:)     !<model forcing ensemble (vector over variables)
  TYPE(assim_state_ensemble_type), ALLOCATABLE :: A(:)     !<model auxiliary ensemble (vector over variables)
  TYPE(assim_state_ensemble_type), ALLOCATABLE :: A2(:)    !<model auxiliary ensemble (vector over variables) !not used in HYPE
  !TYPE(assim_input_ensemble_type), ALLOCATABLE :: P(:)     !model parameter ensemble (vector over variables)
  TYPE(assim_input_ensemble_type), ALLOCATABLE :: Obs(:)   !<observation ensemble (vector over variables)
  
  !observations, predicted observations, observation error covariance
     !TYPE(assim_ensemble_type)                    :: D        !observations    -"-      (incl all variables, used for next analysis, may span several timesteps)
     !TYPE(assim_ensemble_type)                    :: HX       !predictions(HX) -"-      (incl all variables, used for next analysis, may span several timesteps)
     !TYPE(assim_ensemble_type)                    :: R        !obs.error covar matrix   (incl all variables, used for next analysis, may span several timesteps)
  REAL, ALLOCATABLE :: D(:,:)   !<observation matrix
  REAL, ALLOCATABLE :: HX(:,:)  !<predictions matrix
  REAL, ALLOCATABLE :: R(:)     !<obs.error covariance matrix
  
 !spatial coordinates and localization data
 TYPE(assim_coordinate_type), ALLOCATABLE :: Coordinates(:) !<spatial coordinates
     !TYPE(assim_ensemble_type)                    :: LocCYY
     !TYPE(assim_ensemble_type),       ALLOCATABLE :: LocCXY(:)
  REAL, ALLOCATABLE :: LocCYY(:,:)  !<localization data
!  REAL, ALLOCATABLE :: LocCXY(:,:) !CP161214 need one LocCXY for every Coordinate system
  TYPE(assim_twoDmatrix_type), ALLOCATABLE :: LocCXY(:)  !<localization data

END TYPE assim_data_type

!------------------------------------------------------------------------------
!Declaration of the variables needed for the assimilation application
!------------------------------------------------------------------------------
!myAssimData with all necessary ensemble data etc
TYPE(assim_data_type), SAVE :: myAssimData  !<variable holding all nescessary ensemble data for assimilation
!---------------------------------------------------
!file id for input/output files
INTEGER fid_assim_info  !<file unit for AssimInfo.txt file
!INTEGER fid_assim_info  !, fid_assim_min, fid_assim_max, fid_assim_median, fid_assim_final  !CP not used
!INTEGER fid_assim_mean, fid_assim_min95, fid_assim_max95, fid_assim_std  !CP not used
!file id for binary direct access files
INTEGER, ALLOCATABLE :: fid_assim_bin(:)    !<file unit for binary direct access files (fid_assim_bin_base-)
INTEGER              :: fid_assim_bin_base  !<first file unit for binary direct access files
!defalt min and max values
REAL, PARAMETER :: assim_defmin = -1.1E38   !<default minimum value in assimilation
REAL, PARAMETER :: assim_defmax = 1.1E38    !<default maximum value in assimilation
REAL, PARAMETER :: assim_minsigma = 1.E-6   !<default minimum value for sigma used in ensemble generation


END MODULE ASSIMILATION_VARIABLES