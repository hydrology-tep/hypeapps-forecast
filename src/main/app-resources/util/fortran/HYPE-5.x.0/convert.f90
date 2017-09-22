!> \file convert.f90
!> Contains module convert.

!> Procedures for transforming between different variable types and between upper and lower case
MODULE CONVERT

!Copyright 2011-2013,2015,2017 SMHI
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
  USE LIBDATE, ONLY : DateType

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: lower_case,&
            scalar_lower_case,&
            get_lower_case,&
            scalar_upper_case,&
            string_convert_to_datetype,&
            logical_convert_to_integer,&
            integer_convert_to_logical
  
CONTAINS
  
  !>Transforms a string array of dimension n to all lowercase letters
  !---------------------------------------------------------------------------------------------
  SUBROUTINE lower_case(str,n) 
  
    !Argument declarations
    CHARACTER(LEN=*), INTENT(INOUT) :: str(n) !<string-array to convert
    INTEGER, INTENT(IN) :: n                  !<dimension of str
    
    !Local variables
    INTEGER i,j,l

    l = LEN(str)
    DO i=1,n
      DO j=1,l
        IF(ICHAR(str(i)(j:j)).GE.65 .and. ICHAR(str(i)(j:j)).LE.90)THEN
          str(i)(j:j)=CHAR(ICHAR(str(i)(j:j))+32)
        ENDIF
      ENDDO
    ENDDO
  END SUBROUTINE lower_case

  !>Transforms a string to all lowercase letters
  !---------------------------------------------------------------------------------------------
  SUBROUTINE scalar_lower_case(str) 
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(INOUT) :: str  !<string to convert
    
    !Local variables
    INTEGER j,l

    l = LEN(str)
    DO j=1,l
      IF(ICHAR(str(j:j)).GE.65 .and. ICHAR(str(j:j)).LE.90)THEN
        str(j:j)=CHAR(ICHAR(str(j:j))+32)
      ENDIF
    ENDDO
  END SUBROUTINE scalar_lower_case

  !>Transforms a string to all lowercase letters
  !---------------------------------------------------------------------------------------------
  FUNCTION get_lower_case(length,inputstr) RESULT (lowcasestring)
    
    !Argument declarations
    INTEGER, INTENT(IN) :: length  !<character string length
    CHARACTER(LEN=length), INTENT(IN) :: inputstr  !<character string to convert
    CHARACTER(LEN=length) :: lowcasestring  !<lower case string
    
    !Local variables
    INTEGER i

    lowcasestring = inputstr
    DO i=1,length
      IF(ICHAR(lowcasestring(i:i)).GE.65 .and. ICHAR(lowcasestring(i:i)).LE.90)THEN
        lowcasestring(i:i)=CHAR(ICHAR(lowcasestring(i:i))+32)
      ENDIF
    ENDDO
    
  END FUNCTION get_lower_case

  !>Transforms a string to capital letters
  !---------------------------------------------------------------------------------------------
  SUBROUTINE scalar_upper_case(str) 
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(INOUT) :: str  !<string to convert
    
    !Local variables
    INTEGER j,l

    l = LEN(str)
    DO j=1,l
      IF(ICHAR(str(j:j)).GE.97 .and. ICHAR(str(j:j)).LE.122)THEN
        str(j:j)=CHAR(ICHAR(str(j:j))-32)
      ENDIF
    ENDDO
  END SUBROUTINE scalar_upper_case

  !>Creates DateType object from date stored as string
  !>Handle date formats: yyyy-mm-dd, yyyymmdd, yyyy-mm-dd hh:mm
  !!TODO: Add format yyyymmddhhmm
  !------------------------------------------------------------------------
  SUBROUTINE string_convert_to_datetype(datestr,date)
    
    !Argument declarations
    CHARACTER(LEN=*), INTENT(IN) :: datestr !<string containing a date
    TYPE(datetype), INTENT(OUT)  :: date    !<datetype object corresponding to datestr
    
    ! Local variables    
    INTEGER i,l
    
    l = LEN(TRIM(datestr))
    READ(datestr,'(I4)') date%Year
    i = 5
    IF((datestr(i:i)<=CHAR(47)).OR.(datestr(i:i)>=CHAR(58)))THEN !Not 0-9
      i = i+1
    ENDIF
    READ(datestr(i:i+1),'(I2)') date%Month
    i = i+2
    IF((datestr(i:i)<=CHAR(47)).OR.(datestr(i:i)>=CHAR(58)))THEN !Not 0-9
      i = i+1
    ENDIF
    READ(datestr(i:i+1),'(I2)') date%Day
    i = i+2
    IF(l<i)THEN !String does not contain hours and minutes
      date%Hour = 0    !beginning of last(only) timestep of day
      date%Minute = 0  !not used
      RETURN
    ENDIF
    IF((datestr(i:i)<=CHAR(47)).OR.(datestr(i:i)>=CHAR(58)))THEN !Not 0-9
      i = i+1
    ENDIF
    READ(datestr(i:i+1),'(I2)') date%Hour
    i = i+2
    IF((datestr(i:i)<=CHAR(47)).OR.(datestr(i:i)>=CHAR(58)))THEN !Not 0-9
      i = i+1
    ENDIF
    READ(datestr(i:i+1),'(I2)') date%Minute
  END SUBROUTINE string_convert_to_datetype
  
  !>Turns logical to integer -1 for true and -2 for false
  !-------------------------------------------------------
  INTEGER FUNCTION logical_convert_to_integer(logvar)
  
  !Argument declaration
  LOGICAL, INTENT(IN):: logvar  !<logical value to be converted
  
  !Local variables
  INTEGER intvar
  
  IF(logvar)THEN
    intvar = -1   !=.TRUE.
  ELSE
    intvar = -2   !=.FALSE.
  ENDIF
  logical_convert_to_integer = intvar
  
  END FUNCTION logical_convert_to_integer

  !>Turns integer to logical; 1 is true, other is false
  !-----------------------------------------------------
  LOGICAL FUNCTION integer_convert_to_logical(intvar)
  
  !Argument declaration
  INTEGER, INTENT(IN):: intvar  !<integer value to be converted
  
  !Local variables
  LOGICAL logvar
  
  IF(intvar==-1)THEN
    logvar = .TRUE.
  ELSE
    logvar = .FALSE.
  !ELSEIF(intvar==0)THEN
  !  logvar = .FALSE.
  !ELSE
  !  !ERROR
  !  WRITE(6,*) 'ERROR in code. Trying to transform strange integer to logical variable'
  !  STOP 1
  ENDIF
  integer_convert_to_logical = logvar
  
  END FUNCTION integer_convert_to_logical
  
END MODULE CONVERT
