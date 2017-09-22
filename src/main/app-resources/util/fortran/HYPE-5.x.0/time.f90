!> \file time.f90
!> Contains module timeroutines.

!>Module contains procedures relating to time calculations.
MODULE TIMEROUTINES

  !Copyright 2011-2013,2015-2017 SMHI
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

  !-----------------------------------------------------------------------------------------
  !Private procedures in this module
  !-----------------------------------------------------------------------------------------
  USE LIBDATE
  !The module uses modvar and worldvar
  IMPLICIT NONE
  PRIVATE
  PUBLIC ::  period_length, &
             day_of_month, &
             numdays_of_year, &
             calculate_time_for_model, &
             set_timestep_variables, &
             get_dayno_from_monthday, &
             add_date_to_array

  !Private parameters, global in this module
  INTEGER,PARAMETER :: float = KIND(1.0)

CONTAINS

  !>Calculates the number of time steps between sdate and edate
  !--------------------------------------------------------------------
  INTEGER FUNCTION period_length(sdate,edate)

    USE WORLDVAR, ONLY : steplen

    !Argument declarations
    TYPE(datetype), INTENT(IN) :: sdate !<start date
    TYPE(datetype), INTENT(IN) :: edate !<end date
    
    !Local variables
    TYPE(datetype)       :: period
    TYPE(juliandatetype) :: julstart,julend,tmp
    REAL                 :: nrsteps

    IF(steplen%year.NE.0.OR.steplen%month.NE.0) THEN
      WRITE(6,*) 'ERROR: monthly or yearly time step not implemented.'
      STOP 1
    ENDIF

    julstart = date2julian(sdate)
    julend = date2julian(edate)

    tmp%head = julend%head - julstart%head
    tmp%tail = julend%tail - julstart%tail
    IF(tmp%tail.LT.0) THEN
      tmp%tail = tmp%tail + 1._float
      tmp%head = tmp%head - 1._float
    ENDIF

    period = julian2dhm(tmp)
    nrsteps = (period%day*24.*60. + period%hour*60. + period%minute)/ &
         (steplen%day*24.*60. + steplen%hour*60. + steplen%minute)
    IF(MOD(nrsteps,1.).NE.0) THEN
      WRITE(6,*) 'WARNING! Defined period does not match steplength.'
    ENDIF

    period_length = INT(nrsteps)

  END FUNCTION period_length

  !>Calculates the number of days in month during year
  !--------------------------------------------------------------------
  INTEGER FUNCTION day_of_month(year,month)

    !Argument declaration
    INTEGER, INTENT(IN) :: year   !<year in format YYYY
    INTEGER, INTENT(IN) :: month  !<month
    
    !Local variables
    INTEGER x

    SELECT CASE(month)
    CASE(1,3,5,7,8,10,12)
      x=31
    CASE(4,6,9,11)
      x=30
    CASE(2)
      IF(LeapYear(year))THEN
        x=29
      ELSE
        x=28
      ENDIF
    CASE DEFAULT
      WRITE(6,*) 'ERROR in code, while trying to calculate days of month', month
      STOP 1
    END SELECT

    day_of_month = x

  END FUNCTION day_of_month

  !>Calculates the number of days of year
  !--------------------------------------------------------------------
  INTEGER FUNCTION numdays_of_year(date)

    !Argument declarations
    TYPE(datetype), INTENT(IN) :: date  !<date
    
    IF(leapyear(date%year))THEN
      numdays_of_year=366
    ELSE
      numdays_of_year=365
    ENDIF

  END FUNCTION numdays_of_year

  !>Calculates the time variables available for the model
  !>
  !>\b Consequences Module modvar variables dayno, month, prevdoy, currentdate 
  !>tsofday and endofday are set. 
  !---------------------------------------------------------------------
  SUBROUTINE calculate_time_for_model(idt,date)

    USE MODVAR, ONLY : dayno,month,prevdoy,currentdate,tsofday,endofday, &    !OUT
                       timesteps_per_day
    USE WORLDVAR, ONLY : idtlag   !number of timesteps between 00:00 and first timestep

    !Argument declarations
    INTEGER, INTENT(IN)        :: idt    !<index of current timestep
    TYPE(datetype), INTENT(IN) :: date   !<current date
    
    !Local variables
    TYPE(datetype) yearago

    month = date%month
    dayno = doy(date)

    yearago = datetype(date%year-1,date%month,date%day,date%hour,date%minute)
    prevdoy = numdays_of_year(yearago)
    currentdate = date

    endofday = .FALSE.
    tsofday = MOD(idt+idtlag,timesteps_per_day)
    IF(tsofday==0)THEN
      tsofday = timesteps_per_day 
      endofday = .TRUE.
    ENDIF

  END SUBROUTINE calculate_time_for_model

  !>Set time step variables
  !>
  !>\b Consequences Module modvar variables seconds_per_timestep and 
  !>timesteps_per_day are set. Module worldvar variable steplen and idtlag are set.
  !-------------------------------------------------------------------
  INTEGER FUNCTION set_timestep_variables(tstep,stepunit,fdate)

    USE WORLDVAR, ONLY : steplen, &  !OUT
                         idtlag      !OUT
    USE MODVAR, ONLY : seconds_per_timestep, &  !OUT
                       timesteps_per_day !OUT

    INTEGER, INTENT(IN)  :: tstep               !<time step read from file
    CHARACTER(LEN=*), INTENT(IN) :: stepunit    !<unit of time step read from file
    TYPE(datetype), INTENT(IN) :: fdate         !<first date of simulation

    !Calculate step length of datetype
    SELECT CASE(stepunit)
    CASE ('d')
      steplen=datetype(0,0,tstep,0,0)
    CASE ('h')
      steplen=datetype(0,0,0,tstep,0)
    CASE ('min')
      steplen=datetype(0,0,0,0,tstep)
    CASE DEFAULT
      WRITE(6,*) 'ERROR: Unknown stepunit in info.txt reading steplength'
      set_timestep_variables = 1
    END SELECT

    !Calculate step length in seconds
    seconds_per_timestep = NINT(DelaySeconds(steplen))
    timesteps_per_day = NINT(86400./seconds_per_timestep)
    
    !Calculate number of time steps simulation deviate from 00:00
    idtlag = 0  !for daily time step
    IF(fdate%Hour/=0) idtlag = fdate%Hour/(24/timesteps_per_day)

    set_timestep_variables = 0

  END FUNCTION set_timestep_variables

  !>Calculate daynumber of year from month-day (mmdd)
  !------------------------------------------------------------
  SUBROUTINE get_dayno_from_monthday(datumIn,daysAmount)

    !Argument declarations
    INTEGER, INTENT(IN)  :: datumIn       !<month day (mmdd): 31st of January = 0131
    INTEGER, INTENT(OUT) :: daysAmount    !<Corresponding amount of days since beginning of year
    
    !Local variables
    INTEGER month          !Current month
    INTEGER day,feb        !Current day and number of days in February

    month = datumIn/100
    day = datumIn - 100*month
    feb = 28  !since no year, no leap year is assumed
    SELECT CASE (month)
    CASE (1)
      daysAmount = day
    CASE (2)
      daysAmount = 31 + day
    CASE (3)
      daysAmount = 31 + feb + day
    CASE (4)
      daysAmount = 62 + feb + day
    CASE (5)
      daysAmount = 92 + feb + day
    CASE (6)
      daysAmount = 123 + feb + day
    CASE (7)
      daysAmount = 153 + feb + day
    CASE (8)
      daysAmount = 184 + feb + day
    CASE (9)
      daysAmount = 215 + feb + day
    CASE (10)
      daysAmount = 245 + feb + day
    CASE (11)
      daysAmount = 276 + feb + day
    CASE (12)
      daysAmount = 306 + feb + day
    CASE DEFAULT
      daysAmount = 0
    END SELECT

  END SUBROUTINE get_dayno_from_monthday

!>Add new date to date-array in order. 
!>Already existing dates are not added.
!----------------------------------------------------------
  SUBROUTINE add_date_to_array(onedate,dim,n,datearray)
  
    !Argument declarations
    TYPE(datetype), INTENT(IN) :: onedate  !<date to be added
    INTEGER, INTENT(IN)        :: dim      !<size of array
    INTEGER, INTENT(INOUT)     :: n        !<number of dates in array
    TYPE(datetype), INTENT(INOUT) :: datearray(dim) !<array of dates
    
    !Local varaibles
    INTEGER j
    !Start of subroutine
    IF(n==0)THEN
      datearray(n+1) = onedate
      n = n + 1
    ELSE
      IF(onedate>datearray(n))THEN
        datearray(n+1) = onedate
        n = n + 1
      ELSE
        DO j = 1,n
          IF(onedate==datearray(j))THEN
            EXIT
          ELSEIF(onedate<datearray(j))THEN
            datearray(j+1:n+1) = datearray(j:n)
            datearray(j) = onedate
            n = n + 1
            EXIT
          ENDIF
        ENDDO
      ENDIF
    ENDIF 

  END SUBROUTINE add_date_to_array

END MODULE TIMEROUTINES
