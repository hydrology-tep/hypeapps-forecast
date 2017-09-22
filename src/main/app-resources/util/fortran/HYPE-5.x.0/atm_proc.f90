!> \file atm_proc.f90
!> Contains module atmospheric_processes.

!>Subroutines for calculating current atmospheric forcing
MODULE ATMOSPHERIC_PROCESSES

  !Copyright 2014-2017 SMHI
  !
  !This file is part of HYPE.
  !HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
  !HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
  !You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.
  !-----------------------------------------------------------------------------------------

  !Uses modvar, hypevariables, hype_indata

  IMPLICIT NONE
  PRIVATE 
  !----------------------------------------------
  ! Private procedures 
  !----------------------------------------------
  !saturationpressure_function
  !netlongwaveradiation_function
  !calculate_vapour_pressures
  !calculate_shortwave_radiation
  !calculate_net_radiation
  !----------------------------------------------
  PUBLIC :: set_atmospheric_parameters_corrections, &
            calculate_class_atmospheric_forcing,  &
            calculate_subbasin_temperature, &
            calculate_subbasin_precipitation, &
            calculate_rain_snow_from_precipitation, & 
            calculate_extraterrestrial_radiation, &
            deltasaturationpressure_function, &
            set_precipitation_concentration,  &
            calculate_class_wind_transformation_factor, &
            calculate_daylength

  CONTAINS

  !>Calculate basin temperature and precipitation corrections
  !>
  !>\b Consequences Module hypevariables variable basintcalt, basintempadd, 
  !>basinpreccorr, basinpcurain and basinpcusnow are allocated and set.
  !>
  !>\b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  SUBROUTINE set_atmospheric_parameters_corrections()
    
    USE MODVAR, ONLY : conductregest, &
                       genpar, &
                       regpar, &
                       regiondivision, &
                       tobselevation, &
                       basin, &
                       classbasin, &
                       nsub, &
                       nclass
    USE HYPE_INDATA, ONLY : set_regest_parameter
    USE HYPEVARIABLES, ONLY : m_tcalt,    &
                              m_pcelevth, &
                              m_pcelevadd,  &
                              m_pcelevmax,  &
                              m_pcelevstd,  &
                              m_tcobselev,   &
                              m_tcadd, m_tcelevadd, &
                              m_pcurain, &
                              m_pcusnow, &
                              n_tcalt,n_tcea,n_tcorr,  &
                              n_pcet,n_pcea,n_pcem,  &
                              n_pcur,n_pcus,  &
                              basintcalt,   &   !OUT
                              basintempadd, &   !OUT
                              basinpreccorr, &  !OUT
                              basinpcurain,  &  !OUT
                              basinpcusnow      !OUT

    !Local variables
    INTEGER i,j,isb
    REAL temppar1,temppar2
    REAL classheight  !masl
    REAL pcelevth,pcelevadd,pcelevmax   !loop-values
    REAL preccorr_hight
     
    !>\b Algoritm \n

    !>Set temperature altitude correction
    IF(.NOT.ALLOCATED(basintcalt)) ALLOCATE(basintcalt(nsub))
    basintcalt = genpar(m_tcalt)
    IF(conductregest)THEN
      DO i = 1,nsub
        CALL set_regest_parameter(i,n_tcalt,basintcalt(i))
      ENDDO
    ENDIF

    !>Set temperature correction
    IF(.NOT.ALLOCATED(basintempadd)) ALLOCATE(basintempadd(nsub))
    basintempadd = 0.
    IF(genpar(m_tcobselev)/=0.AND. ALLOCATED(tobselevation)) basintempadd = basintempadd - genpar(m_tcobselev)*(basin%elev-tobselevation)*0.01    !adjust for Tobs elevation
    temppar1 = genpar(m_tcelevadd)
    DO i = 1,nsub
      IF(basin(i)%parregion(regiondivision(m_tcadd))>0)THEN
        temppar2 = regpar(m_tcadd,basin(i)%parregion(regiondivision(m_tcadd)))
      ELSE
        temppar2  = 0.
      ENDIF
      
      !Replace parameter values with regional parameter estimates
      IF(conductregest)THEN
        CALL set_regest_parameter(i,n_tcea,temppar1)
        CALL set_regest_parameter(i,n_tcorr,temppar2)
      ENDIF
      
      !Adjust subbasin air temperature for observation elevation, subbasin elevation and regional correction.
      basintempadd(i) = basintempadd(i) + temppar2 - temppar1*basin(i)%elev*0.01       !Corrected subbasin temperature
    ENDDO

    !>Set precipitation height correction
    IF(.NOT.ALLOCATED(basinpreccorr)) ALLOCATE(basinpreccorr(nsub,nclass))
    pcelevth  = genpar(m_pcelevth)
    pcelevadd = genpar(m_pcelevadd)
    pcelevmax = genpar(m_pcelevmax)
    DO i = 1,nsub
      DO j= 1,nclass
        !Replace parameter values with regional parameter estimates
        IF(conductregest)THEN
          CALL set_regest_parameter(i,n_pcet,pcelevth)
          CALL set_regest_parameter(i,n_pcea,pcelevadd)
          CALL set_regest_parameter(i,n_pcem,pcelevmax)
        ENDIF
        classheight = basin(i)%elev+classbasin(i,j)%deltah  !masl

        !Calculate height correction
        IF(classheight>pcelevth)THEN
          preccorr_hight = (classheight-pcelevth) / 100. *  pcelevadd
          preccorr_hight = preccorr_hight + basin(i)%selev / 100. * genpar(m_pcelevstd)
          IF(preccorr_hight > pcelevmax) preccorr_hight = pcelevmax !maxmimum for prec hight correction
          basinpreccorr(i,j) = 1. + preccorr_hight
        ELSE
          basinpreccorr(i,j) = 1.
        ENDIF
      ENDDO
    ENDDO

    !>Set subbasin rain and snow parameters
    IF(.NOT.ALLOCATED(basinpcurain))THEN
      ALLOCATE(basinpcurain(nsub))
      ALLOCATE(basinpcusnow(nsub))
    ENDIF
    basinpcurain = genpar(m_pcurain)
    basinpcusnow = genpar(m_pcusnow)
    !>Replace parameter values with regional parameter estimates
    IF(conductregest)THEN
      DO isb=1,nsub
        CALL set_regest_parameter(isb,n_pcur,basinpcurain(isb))
        CALL set_regest_parameter(isb,n_pcus,basinpcusnow(isb))
      ENDDO
    ENDIF

  END SUBROUTINE set_atmospheric_parameters_corrections

  !>Calculate class temperature and precipitation
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  SUBROUTINE calculate_class_atmospheric_forcing(i,j,radext,  &
             temp,prec,tmin,tmax,swrad,rhmin,actvap,satvap,icpevap,netrad,wind,sffrac)
    
    USE MODVAR, ONLY : tempi,     &
                       preci,     &
                       tmini,     &
                       tmaxi,     &
                       windi,     &
                       humidi,    &
                       shortwavei,  &
                       snowfraci, &
                       genpar,    &
                       landpar,   &
                       basin,     &
                       classbasin,  &
                       classdata, &
                       missing_value
    USE HYPEVARIABLES, ONLY : m_pcluse,   &
                              m_alb,      &
                              m_mwind,    &
                              basintcalt, &
                              basinpreccorr, &
                              windtrans,  &
                              calcSWRAD,  &
                              calcVAPOUR, &
                              calcWIND
    
    !Argument declarations
    INTEGER,INTENT(IN) :: i       !<current index of subbasin
    INTEGER,INTENT(IN) :: j       !<current index of class
    REAL, INTENT(IN)  :: radext   !<subbasin extraterrestrial solar radiation (MJ/m2/day)
    REAL, INTENT(OUT) :: temp     !<current class temperature (C)
    REAL, INTENT(OUT) :: prec     !<current class precipitation (mm/timestep)
    REAL, INTENT(OUT) :: tmin     !<current daily min temperature (C)
    REAL, INTENT(OUT) :: tmax     !<current daily max temperature (C)
    REAL, INTENT(OUT) :: swrad    !<daily mean shortwave radiation (MJ/m2/day)
    REAL, INTENT(OUT) :: rhmin    !<daily min relative humidity [fraction 0-1]
    REAL, INTENT(OUT) :: actvap   !<actual vapour pressure [kPa]
    REAL, INTENT(OUT) :: satvap   !<saturated vapour pressure [kPa]
    REAL, INTENT(OUT) :: icpevap  !<interception losses (evaporation) [mm]
    REAL, INTENT(OUT) :: netrad   !<net downward radiation [MJ/m2/day]
    REAL, INTENT(OUT) :: wind     !<current wind speed (m/s)
    REAL, INTENT(OUT) :: sffrac   !<fraction of precipitation as snow (-)
   
    !Local variables
    REAL rhmax,rhmean   !for future use: daily max/mean relative humidity
    REAL albedo         !current albedo (-)  
    REAL relswrad       !daily relative shortwave radiation (MJ/m2/day)
     
    !>\b Algoritm \n
    !Default output values, missing
    tmin = missing_value
    tmax = missing_value
    swrad = missing_value
    rhmin = missing_value
    actvap = missing_value
    satvap = missing_value
    icpevap = 0.
    netrad = missing_value
    wind = missing_value
    sffrac = missing_value
      
    !Initiation of other variables
    rhmax = missing_value
    rhmean = missing_value
    relswrad = missing_value

    !>Calculate class temperature and precipitation adjusted for elevation and landuse bias
    temp = tempi(i)-basintcalt(i)*classbasin(i,j)%deltah*0.01 !Adds SLC-hight dependent correction for temperature
    prec = preci(i) * basinpreccorr(i,j)
    IF(landpar(m_pcluse,classdata(j)%luse)>0.) icpevap = prec * landpar(m_pcluse,classdata(j)%luse)
    prec = prec * (1. - landpar(m_pcluse,classdata(j)%luse))
               
    !>Calculate class minimum and maximum air temperature, if available
    IF(ALLOCATED(tmini)) tmin = tmini(i)-basintcalt(i)*classbasin(i,j)%deltah*0.01 ! Adds SLC-hight dependent correction for minimum temperature
    IF(ALLOCATED(tmaxi)) tmax = tmaxi(i)-basintcalt(i)*classbasin(i,j)%deltah*0.01 ! Adds SLC-hight dependent correction for maximum temperature
      
    !>Calculate class specific shortwave radiation, if needed
    IF(calcSWRAD)THEN
      IF(ALLOCATED(shortwavei)) swrad = shortwavei(i)
      CALL calculate_shortwave_radiation(tmin,tmax,(basin(i)%elev + classbasin(i,j)%deltah),radext,swrad,relswrad)
    ENDIF

    !>Calculate class specific vapor pressures and net radiation, if needed (and possibly tmin and tmax)
    IF(calcVAPOUR)THEN
      IF(ALLOCATED(humidi)) rhmean = humidi(i)
      CALL calculate_vapour_pressures(temp,tmin,tmax,rhmean,rhmin,rhmax,swrad,radext,actvap,satvap)
      albedo = landpar(m_alb,classdata(j)%luse)
      IF(albedo.LE.0.) albedo = 0.23 ! FAO standard albedo
      CALL calculate_net_radiation(temp,tmin,tmax,albedo,actvap,swrad,relswrad,netrad)
    ENDIF
      
    !>Set class specific windspeed
    IF(calcWIND)THEN
      wind = genpar(m_mwind)
      IF(ALLOCATED(windi)) wind = windtrans(j) * windi(i)
    ENDIF

    !>Set snow fall fraction from file
    IF(ALLOCATED(snowfraci)) sffrac = snowfraci(i)

  END SUBROUTINE calculate_class_atmospheric_forcing

  !>Apply corrections to get subbasin average temperature
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !----------------------------------------------------------
  SUBROUTINE calculate_subbasin_temperature(n,month,temparr)
    
    USE HYPEVARIABLES, ONLY : basintempadd, &
                              m_mlapse
    USE MODVAR, ONLY : monthpar,  &
                       basin
       
    !Argument declaration
    INTEGER, INTENT(IN) :: n           !<number of subbasins
    INTEGER, INTENT(IN) :: month       !<current month
    REAL, INTENT(INOUT) :: temparr(n)  !<subbasin average temperature
    
    !>\b Algoritm \n
    
    !>Add temperature correction to current air temperature
    temparr = temparr + basintempadd - monthpar(m_mlapse,month)*basin%elev*0.01
    
  END SUBROUTINE calculate_subbasin_temperature

  !>Apply corrections to get subbasin average precipitation
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !-----------------------------------------------------------
  SUBROUTINE calculate_subbasin_precipitation(n,temppobs,precarr,pcorricep)
    
    USE MODVAR, ONLY : genpar,regpar, &
                       regiondivision, &
                       basin, &
                       snowfraci, &
                       classbasin, &
                       classdata, &
                       nclass
    USE HYPEVARIABLES, ONLY : m_preccorr,   &
                              m_pcaddg, & 
                              basinpcurain,basinpcusnow

    !Argument declaration
    INTEGER, INTENT(IN) :: n                 !<number of subbasins
    REAL, INTENT(IN)    :: temppobs(n)       !<temperature at precipitation input level
    REAL, INTENT(INOUT) :: precarr(n)        !<subbasin average precipitation [mm/timestep]
    REAL, INTENT(OUT)   :: pcorricep(n)      !<interception evaporation due to negative preccorr-parameter [mm/timestep]
    
    !Local variables
    INTEGER i,j
    REAL preccorr
    REAL sffrac, snowfall, rainfall, snowfrac, preccorr_undercatch
 
    pcorricep = 0.    !Initiate default value, no interception evaporation
    !>\b Algoritm \n
    DO i = 1,n
      !>Calculate average snowfall fraction
    !DG 20141113 undercatch correction of precipitation, dependent on snowfall fraction, applied before general and regional corrections:
      IF(ALLOCATED(snowfraci)) sffrac = snowfraci(i) !set snowfall fraction from input data, if existing
      ! loop over classes to derive average snowfall fraction dependent on class landuse and areal fractions in the subbasin.
      !   In the subbasin, the class snowfall fraction will also depend on the classtemperature, however, for the undercatch
      !   correction, we derive the snowfall fraction for each class, as if they were located at the elevation of 
      !   the precipitation input data. Ideally, we should also know what type of landuse the observation represent. But for the moment,
      !   we assume that the subbasin landuse distribution is relevant also for the input data. It could potentially give much more error
      !   to derive the snowfall fraction using the elevation corrected class temperature - so the current solution should be more conservative.
      snowfrac = 0.
      DO j=1,nclass
        IF(classbasin(i,j)%part.GT.0.)THEN
          ! calculate snowfall fraction with existing subroutine, using unit precipitation
          CALL calculate_rain_snow_from_precipitation(classdata(j)%luse,1.,temppobs(i),snowfall,rainfall,sffrac) !form of precipitation
          ! summing areally weighted snowfall fraction
          snowfrac = snowfrac + snowfall * classbasin(i,j)%part
        ENDIF
      ENDDO
      
      !>Undercatch correction factor give initial correction of precarr(i)
      preccorr_undercatch = basinpcurain(i) * (1. - snowfrac) + basinpcusnow(i) * snowfrac
      precarr(i) = precarr(i)*(1.+preccorr_undercatch)
      
      IF(basin(i)%parregion(regiondivision(m_preccorr))>0)THEN
        preccorr = 1. + regpar(m_preccorr,basin(i)%parregion(regiondivision(m_preccorr)))   !Regional correction of precipitation
      ELSE
        preccorr  = 1.
      ENDIF
      !Calculate interception as negative preccorr
      IF(preccorr<1.) pcorricep(i) = precarr(i)*(1.+genpar(m_pcaddg)) * (1.-preccorr)
      !>Adjust subbasin precipitation for undercatch of snow and rain, general bias and regional correction.
      precarr(i) = precarr(i)*(1.+genpar(m_pcaddg)) * preccorr      !Correction for precipitation          
    ENDDO 

  END SUBROUTINE calculate_subbasin_precipitation

  !>Subroutine for calculation amount of rainfall and snowfall from precipitaion
  !>
  !> \b Reference ModelDescription Chapter Processes above ground (Temperature and precipitation)
  !--------------------------------------------------------------------------
  SUBROUTINE calculate_rain_snow_from_precipitation(iluse,prec,temp,snowfall,rainfall,sffrac)

    USE MODVAR, ONLY : genpar,      &
                       landpar,     &
                       modeloption, &
                       p_snowfall
    USE HYPEVARIABLES, ONLY : m_ttpd,   &
                              m_ttpi,   &
                              m_ttmp

    !Argument declarations
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    REAL, INTENT(IN)    :: prec     !<precipitation (mm/timestep)
    REAL, INTENT(IN)    :: temp     !<air temperature (C)
    REAL, INTENT(OUT)   :: snowfall !<Precipitation as snow (mm/timestep)
    REAL, INTENT(OUT)   :: rainfall !<Precipitation as rain (mm/timestep)
    REAL, INTENT(IN)    :: sffrac   !<fraction of precipitation as snow (-), optional input data
   
    !Local variables
    REAL tt       !threshold temperature for rain/snow precip (C)
    REAL dtt      !temperature interval for rain/snow precip (C)
    REAL arain    !fraction of precipitation as rain (-)

    !> \b Algorithm \n
    !Set local parameter values
    dtt = genpar(m_ttpi)              !temperature interval rain/snow precip
    tt  = landpar(m_ttmp,iluse) + genpar(m_ttpd) !threshold temperature for rain/snow precip (ttmp for melt/evap + diff)
   
    !>Select model for rain/snow separation and calculate rain fraction
    SELECT CASE(modeloption(p_snowfall))
      CASE(1)       !Fraction of precipitation as snowfall is given in input data
        arain = 1. - sffrac
      CASE DEFAULT  !Original model based on threshold temperatures
        arain = 0.  !Pure snow
        IF(prec>0)THEN
          IF(temp >= tt + dtt) THEN 
            arain = 1. !Pure rain
          ELSEIF(temp > tt - dtt .AND. temp < tt + dtt) THEN
            arain = (temp - (tt - dtt)) / (2 * dtt) !snow and rain mix
          ENDIF
        ENDIF
    END SELECT
    
    !>Calculate rainfall and snowfall
    rainfall = prec * arain    
    snowfall = prec * (1-arain)

  END SUBROUTINE calculate_rain_snow_from_precipitation

  !>Calculates extraterrestrial solar radiation as a function of
  !>latitude and day of year for all subbasins in the model domain
  !---------------------------------------------------------------
  SUBROUTINE calculate_extraterrestrial_radiation(n,jday,radext)

    USE MODVAR, ONLY : basin, &
                       pi,    & !=3.1415927
                       solar    !=0.0820 MJ/m2/min, solar constant
     
    !Argument declarations
    INTEGER,INTENT(IN) :: n         !<number of subbasins   
    INTEGER,INTENT(IN) :: jday      !<day of year (1-366)
    REAL,INTENT(OUT)   :: radext(n) !<extraterrestrial solar radiation [MJ/m2/day]
     
    !Local variables
    INTEGER i   !loop index
    REAL latrad !latitude in radians
    REAL dr     !inverse relative distance Earth-Sun
    REAL d      !Solar declination
    REAL omega  !Sunset hour angle
     
    !> \b Algorithm \n
    !>Calculate variables needed for extraterrestrial radiation calculation; distance to sun and declination.  
    dr = 1 + 0.033 * cos(2. * pi * jday / 365.)    !Inverse relative distance Earth-Sun (FAO, ekv 23)
    d = 0.409 * sin(2. * pi * jday / 365. - 1.39)  !Solar declination (FAO, ekv 24)

    !>For every subbasin:
    DO i = 1,n     
      latrad = basin(i)%latitude * pi / 180.0     !Transform basin latitude from degrees to radians 

      !>\li Calculate sunset hour angle, with special care for high latitudes
      omega = -tan(latrad) * tan(d) ! first check value of omega=cos(omega)
      IF(omega.GT.1)THEN
        omega = 0.              !Polar night, cos(omega)>1, set omega = 0
      ELSEIF(omega.LT.-1)THEN
        omega = pi              !Midnight sun, cos(omega)<1, set omega = pi
      ELSE
        omega = acos(-tan(latrad) * tan(d))    !Sunset hour angle,(FAO, ekv 25)
      ENDIF
       
      !>\li Calculate extraterrestrial radiation
      IF(omega.GT.0.)THEN
        radext(i) = 1440. / pi * solar * dr * (omega * sin(latrad)*sin(d)+cos(latrad)*cos(d)*sin(omega))
      ELSE
        radext(i) = 0.
      ENDIF
    ENDDO
    
  END SUBROUTINE calculate_extraterrestrial_radiation
   
  !>Calculates downward?? solar radiation at land surface 
  !---------------------------------------------------------------
  SUBROUTINE calculate_shortwave_radiation(tmin,tmax,elev,radext,swrad,relswrad)

    USE MODVAR, ONLY: genpar
    USE HYPEVARIABLES, ONLY: m_krs

    !Argument declarations
    REAL, INTENT(IN)  :: tmin           !<daily minimum air temperature [C]
    REAL, INTENT(IN)  :: tmax           !<daily maximum air temperature [C]
    REAL, INTENT(IN)  :: elev           !<elevation [m.a.s.l]
    REAL, INTENT(IN)  :: radext         !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(INOUT) :: swrad        !<daily mean shortwave radiation [MJ/m2/day]
    REAL, INTENT(OUT) :: relswrad       !<relative shortwave radiation (actual/clearsky) [-]
  
    !Local parameters
    REAL, PARAMETER :: turbmin = 0.25   ! Ångström formula, minimum turbidity
    REAL, PARAMETER :: turbmax = 0.75   ! Ångström formula, maximum turbidity (at sea level)
    
    !Local variables
    REAL turbidity                      ! Turbidity (surface_radiation/top_of_atmosphere_radiation)
    REAL turbidity_clearsky             ! Ångströms clear sky turbidity
    REAL turbidity_hargreaves           ! Hargreaves turbidity krs * (Tmax-Tmin)^0.5
 
    !Initialize turbidity as missing
    turbidity = -9999.
     
    !Clear sky turbidity, Ångström formula, with elevation correction from FAO ekv 37:
    turbidity_clearsky = min(1.,turbmax + elev * 2.E-5)
   
    !Estimate turbidity from the input radiation data - if available - limited by Ångström turbidity
    IF(swrad.GT.-9999. .AND. radext.GT.0.)THEN
      turbidity = max(turbmin,min(turbidity_clearsky, swrad / radext))
    ENDIF     
   
    !Estimate Hargreaves turbidity from tmin and tmax, if available
    IF(tmin.GT.-9999. .AND. tmax.GT.tmin)THEN
      !Hargreaves turbidity function, limited by the Ångströms function
      turbidity_hargreaves = genpar(m_krs) * (tmax - tmin)**0.5
      turbidity_hargreaves = max(turbmin,min(turbidity_clearsky,turbidity_hargreaves))
      !Use Hargreaves turbidity if swrad was missing
      IF(turbidity.LE.-9999.)turbidity = turbidity_hargreaves
    ENDIF
    
    !Assume clear sky turbidity if both swrad and tmin/tmax were missing
    IF(turbidity.LE.-9999.) turbidity = turbidity_clearsky

    !Calculate downward?? shortwave radiation, if missing
    IF(swrad.LE.-9999.) swrad = radext * turbidity
    
    !Calculate relative shortwave radiation (shortwave/clearsky_shortwave)
    IF(turbidity_clearsky.GT.0 .AND. turbidity.GT.0)THEN
      relswrad = turbidity / turbidity_clearsky
    ELSE
      relswrad = 0.
    ENDIF
  END SUBROUTINE calculate_shortwave_radiation
  
  !>Calculates daily mean actual and saturated vapour pressure
  !!depending on data availability, following recommended FAO
  !!procedures.
  !---------------------------------------------------------------
  SUBROUTINE calculate_vapour_pressures(tmean,tmin,tmax,rhmean,rhmin,rhmax,swrad,radext,actvap,satvap)

    USE MODVAR, ONLY: genpar
    USE HYPEVARIABLES, ONLY: m_krs

    !Argument declarations
    REAL, INTENT(IN)    :: tmean     !<daily mean temperature [C]? temperature of timestep 
    REAL, INTENT(INOUT) :: tmin      !<daily min temperature [C]
    REAL, INTENT(INOUT) :: tmax      !<daily max temperature [C]
    REAL, INTENT(IN)    :: rhmean    !<daily mean relative humidity [fraction 0-1]
    REAL, INTENT(IN)    :: rhmin     !<daily min relative humidity [fraction 0-1]
    REAL, INTENT(IN)    :: rhmax     !<daily max relative humidity [fraction 0-1]
    REAL, INTENT(IN)    :: swrad     !<daily mean shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: radext    !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(OUT)   :: actvap    !<actual vapour pressure [kPa]
    REAL, INTENT(OUT)   :: satvap    !<saturated vapour pressure [kPa]
     
    !Local variables
    REAL turbidity, trange
     
    !Initialize output variables as missing
    actvap = -9999. ; satvap = -9999.
     
    !Saturated vapor pressure from Tmin and Tmax (if avalable) or from Tmean
    IF(tmin.GT.-9999. .AND. tmax.GT.-9999.)THEN
      satvap = 0.5 * (saturationpressure_function(tmax)+saturationpressure_function(tmin))
    ELSE
      IF(tmean.GT.-9999) satvap = saturationpressure_function(tmean)
    ENDIF
     
    !Actual vapour pressure, following FAO recommended procedure and function/data priority
    IF(tmin.GT.-9999 .AND. tmax.GT.-9999. .AND. rhmin.GT.-9999 .AND. rhmax.GT.-9999.)THEN
      !FAO, ekv 17, using rhmin, rhmax, tmin, and tmax
      actvap = 0.5 * (saturationpressure_function(tmax) * rhmin + & 
                      saturationpressure_function(tmin) * rhmax)
    ELSEIF(rhmax.GT.-9999. .AND. tmin.GT.-9999.)THEN
      !FAO, ekv 18, using rhmax and tmin
      actvap = saturationpressure_function(tmin) * rhmax
    ELSEIF(rhmean.GT.-9999. .AND. satvap.GT.-9999.)THEN
      !FAO ekv 19, using rhmean and saturation pressure from tmin/max or tmean
      actvap = rhmean * satvap 
    ELSEIF(tmin.GT.-9999)THEN
      !Final solution, taking actual vapour pressure = saturation pressure at tmin
      actvap = saturationpressure_function(tmin)
    ELSEIF(swrad.GT.-9999. .AND. radext.GT.-9999.)THEN
      !As a finalfinal solution, Tmin and Tmax is inferred from 
      ! the Hargreaves turbidity function (if swrad and radext is available), 
      ! and then actvap is estimated from the estimated Tmin
     
      !Turbidity
      turbidity = swrad / radext
      !Diurnal temperature range inferred from Hargreaves turbidity function
      trange = (turbidity/genpar(m_krs))**2
      !Assume min and max is evenly distributed around the mean
      tmin = tmean - 0.5 * trange
      tmax = tmean + 0.5 * trange
      !Actual vapour pressure from tmin
      actvap = saturationpressure_function(tmin)
    ENDIF
     
    !Finally, make sure actvap <= satvap
    IF(actvap.GT.-9999 .AND. satvap.GT.-9999) actvap = min(actvap,satvap)
       
  END SUBROUTINE calculate_vapour_pressures
  
  !>Calculates net radiation at land surface 
  !---------------------------------------------------------------
  SUBROUTINE calculate_net_radiation(tmean,tmin,tmax,albedo,actvap,swrad,relswrad,netrad)
    
    !Argument declarations
    REAL, INTENT(IN)    :: tmean      !<daily mean temperature [C]? temperature of timestep 
    REAL, INTENT(IN)    :: tmin       !<daily min temperature [C]
    REAL, INTENT(IN)    :: tmax       !<daily max temperature [C]
    REAL, INTENT(IN)    :: albedo     !<albedo [fraction, 0-1]
    REAL, INTENT(IN)    :: actvap     !<actual vapour pressure [kPa]
    REAL, INTENT(IN)    :: swrad      !<daily mean shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN)    :: relswrad   !<relative shortwave radiation (actual/clearsky) [-]
    REAL, INTENT(INOUT) :: netrad     !<net downward radiation [MJ/m2/day]
      
    !Local variables
    REAL radnetshort, radnetlong
   
    !Estimate net radiation if missing, following FAO recommended procedure
    IF(netrad.LE.-9999)THEN
      !Net shortwave radiation
      radnetshort = swrad * (1.-albedo)

      !Net downward longwave radiation using tmin/max, actual vapour pressure, and relative shortwave, if available
      IF(tmin.GT.-9999. .AND. tmax.GT.-9999. .AND. relswrad.GT.-9999. .AND. actvap.GT.-9999)THEN
        radnetlong = netlongwaveradiation_function(tmax,tmin,actvap,relswrad)
      ELSEIF(tmean.GT.-9999. .AND. relswrad.GT.-9999. .AND. actvap.GT.-9999)THEN
        !Use Tmean if Tmin and Tmax is missing
        radnetlong = netlongwaveradiation_function(tmean,tmean,actvap,relswrad)
      ELSE
        radnetlong = 0.
      ENDIF

      !net radiation
      netrad = radnetshort - radnetlong
    ENDIF
    
  END SUBROUTINE calculate_net_radiation
 
  !>Saturation pressure (kPa) as a function of temperature in deg C, following FAO
  !-------------------------------------------------------------------------------
  REAL FUNCTION saturationpressure_function(temp)
    
    !Argument declarations
    REAL,INTENT(IN) :: temp
    
    saturationpressure_function = 0.6108 * exp((17.27 * temp)/(temp+237.3))
    
  END FUNCTION saturationpressure_function

  !>Slope of the saturation pressure temperature function (kPa/C)
  !-------------------------------------------------------------------------------
  REAL FUNCTION deltasaturationpressure_function(temp)
    
    !Argument declarations
    REAL,INTENT(IN) :: temp 
    
    deltasaturationpressure_function = 4098 * 0.6108 * exp((17.27 * temp)/(temp+237.3)) / (temp + 237.3)**2
    
  END FUNCTION deltasaturationpressure_function
  
  !>Net (upward) longwave radiation, FAO, ekv 39
  !-------------------------------------------------------------------------------
  REAL FUNCTION netlongwaveradiation_function(tmax,tmin,actvap,relshortwave)
    
    !Argument declarations
    REAL,INTENT(IN) :: tmax,tmin,actvap,relshortwave
     
    !Local parameters
    REAL, PARAMETER :: sigma = 4.903E-9 !Stefan-Boltzmann, MJ K^-4 m^-2 day^-1
    REAL, PARAMETER :: vpar1 = 0.34     !vapour pressure constants
    REAL, PARAMETER :: vpar2 = 0.14
    REAL, PARAMETER :: rpar1 = 1.35     !cloudiness constants
    REAL, PARAMETER :: rpar2 = 0.35
    
    !Net upward longwave radiation
    netlongwaveradiation_function = sigma * 0.5 * ((tmax+273.15)**4 + (tmin+273.15)**4) * &
       (vpar1-vpar2*(actvap)**0.5) * (rpar1 * min(1.,relshortwave)- rpar2)
    
  END FUNCTION netlongwaveradiation_function

  !>Set concentrations of precipitation from time series or parameters
  !>
  !>\b Reference ModelDescription Chapter Processes above ground (Atmospheric deposition of nitrogen and phosphorus)
  !---------------------------------------------------------------------
  SUBROUTINE set_precipitation_concentration(i,ns,cprec)
    
    USE MODVAR, ONLY : i_t1,i_in,i_sp,    &
         xobsi,             &
         xobsindex,         &
         load,              &
         genpar
    USE HYPEVARIABLES, ONLY : o_cprecT1,    &
         o_cprecIN,    &
         o_cprecSP,    &
         m_wetsp
    
    !Argument declarations
    INTEGER, INTENT(IN)    :: i           !<index of current subbasin
    INTEGER, INTENT(IN)    :: ns          !<number of substances, array dimension
    REAL,    INTENT(OUT)   :: cprec(ns)   !<concentration of precipitation 

    !>\b Algorithm \n
    !Initiations
    cprec = 0.

    !>Set T1 concentration
    IF(i_t1>0 .AND. xobsindex(o_cprecT1,i)>0)THEN
      cprec(i_t1) = xobsi(xobsindex(o_cprecT1,i))
    ENDIF
    !>Set wet deposition of inorganic nitrogen and soluble phosphorus
    IF(i_in>0)THEN
      IF(xobsindex(o_cprecIN,i)>0)THEN
        cprec(i_in) = xobsi(xobsindex(o_cprecIN,i))
      ELSE
        cprec(i_in) = load(i)%inwetdep
      ENDIF
    ENDIF
    IF(i_sp>0)THEN
      IF(xobsindex(o_cprecSP,i)>0)THEN
        cprec(i_sp) = xobsi(xobsindex(o_cprecSP,i))
      ELSE
        cprec(i_sp) = genpar(m_wetsp)*1.E-3
      ENDIF
    ENDIF

  END SUBROUTINE set_precipitation_concentration

  !>Calculate transformation factor for wind speed to different height than observations
  !---------------------------------------------------------------------
  SUBROUTINE calculate_class_wind_transformation_factor(windtrans)
  
  USE HYPEVARIABLES, ONLY : m_zwind,    &
                            m_zwish,    &
                            m_zpdh,     &
                            m_roughness
  USE MODVAR, ONLY : nclass,  &
                     genpar

  !Argument declarations
  REAL,ALLOCATABLE, INTENT(INOUT) :: windtrans(:)
  
  !Local parameters
  INTEGER j
  REAL lnz0,d0,zwind,zwish

  !Allocate wind transformation factor  
  IF(.NOT.(ALLOCATED(windtrans))) ALLOCATE(windtrans(nclass))

  !Default: no transformation of wind 
  IF(genpar(m_zwind)==genpar(m_zwish) .OR. &
     (genpar(m_zwind)==0. .OR. genpar(m_zwish)==0.) .AND. &
     genpar(m_roughness)==0. .AND. genpar(m_zpdh)==0.)THEN
    windtrans = 1.
  ELSE
    !Wind is transformed to wanted height
    zwind = genpar(m_zwind)
    zwish = genpar(m_zwish)
    DO j = 1,nclass
      lnz0 = LOG(genpar(m_roughness))   !Could be land use parameters
      d0 = genpar(m_zpdh)
      windtrans(j) = (LOG(zwind-d0)-lnz0)/(LOG(zwish-d0)-lnz0)
    ENDDO
  ENDIF
  
  END SUBROUTINE calculate_class_wind_transformation_factor
  
  !>Calculate daylength based on latitude and julian day
  !---------------------------------------------------------------------
  SUBROUTINE calculate_daylength(dayno,lat,length)
  
  USE MODVAR, ONLY : pi

  !Argument declarations
  INTEGER,INTENT(IN) :: dayno    !<Current julian day number
  REAL, INTENT(IN)   :: lat      !<latitude
  REAL, INTENT(OUT)  :: length   !<day length (hours)
  
  !Local parameters
  REAL dec,a1

    dec = -23.45 * COS(pi*(REAL(dayno) + 10.173)/182.61)
    a1 = MIN(1.,MAX(-1.,(SIN(lat * pi / 180.) * SIN(dec * pi/180)) / (COS(lat * pi / 180.) * COS(dec * pi/180.))))
    length = (1440. - 120./(15.*pi/180.) * ACOS(a1)) /60.
  
  END SUBROUTINE calculate_daylength
  
END MODULE ATMOSPHERIC_PROCESSES
