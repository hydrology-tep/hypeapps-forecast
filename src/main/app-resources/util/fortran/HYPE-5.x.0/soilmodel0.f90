!> \file soilmodel0.f90
!> Contains module soilmodel_default.

!>HYPE default soil model
MODULE SOILMODEL_DEFAULT

  !Copyright 2012-2017 SMHI
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
  !Procedures in this module
  !-----------------------------------------------------------------------------------------
  ! soilmodel_0
  !-----------------------------------------------------------------------------------------

  USE STATETYPE_MODULE
  USE MODVAR,  ONLY: missing_value,pi,   &
                     currentdate, &
                     genpar,landpar,soilpar,  &
                     load,basin,path, &
                     numsubstances,nsub,maxsoillayers, &
                     classdata,soilthick,soildepth, &
                     i_t1,i_t2,i_in,i_on,i_sp,i_pp,i_oc,&
                     modeloption,p_deepgroundwater, &
                     conductN,conductP,conductC,conductS,  &
                     doirrigation
  USE HYPEVARIABLES, ONLY : wpmm,fcmm,epmm,pwmm, &
                            m_ttmp,m_cmlt,m_srrcs, &
                            m_wetsp,m_drypp, m_ponatm,m_cfrost,m_sfrost,m_perc1,m_perc2,  &
                            m_sreroexp,m_filtpbuf,m_filtpinner,m_filtpother,m_macfilt, &
                            m_pprelexp,m_fertdays,m_minerfn,m_minerfp,m_degradhn,  &
                            m_degradhp,m_denitrlu, &
                            m_dissolfN,m_dissolhN,m_dissolfP,m_dissolhP,m_littdays,  &
                            m_crate1,m_crate2,m_crate3,m_crate9,m_crate10,m_minc,m_freuc, &
                            m_freuexp,m_freurate,m_sswcorr,m_immdep,m_iwdfrac,m_wdpar, &
                            m_ripz,m_rips,m_ripe,m_ocsoim,m_ocsmslp,   &
                            soilmem,m_deepmem, &
                            epotdist
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water, &
                                          error_remove_water
  USE ATMOSPHERIC_PROCESSES, ONLY :  calculate_rain_snow_from_precipitation
  USE SOIL_PROCESSES, ONLY : calculate_snow, &
                             calculate_potential_evaporation,          &
                             calculate_actual_soil_evapotranspiration, &
                             add_macropore_flow,    &
                             calculate_tile_drainage,  &
                             calculate_soil_runoff,    &
                             infiltration, &
                             add_infiltration,   &
                             percolation, &
                             calculate_groundwater_table,  &
                             calculate_soiltemp,   &
                             calculate_frostdepth,  &
                             calculate_soil_moisture_deficit
  USE NPC_SOIL_PROCESSES, ONLY : add_dry_deposition_to_landclass,  &
                                 calculate_plant,        &
                                 soil_np_processes,      &
                                 soil_carbon_processes,  &
                                 balance_spsoil,         &
                                 particle_processes_for_runoff, &
                                 local_diffuse_source,   &
                                 class_riparian_zone_processes
  USE TRACER_PROCESSES, ONLY : soil_tracer_processes
  USE IRRIGATION_MODULE, ONLY : apply_irrigation, &
                                calculate_irrigation_water_demand
  USE REGIONAL_GROUNDWATER_MODULE, ONLY : add_regional_groundwater_flow_to_soil

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: soilmodel_0
  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(1)  !error message for location of remove_water call
  PARAMETER (errstring = (/'surface runoff, soillayer 1'/))

CONTAINS

  !----------------------------------------------------------------
  !>\brief Default soilmodel for land classes
  !!Calculate snow and soil processes for a land class. 
  !----------------------------------------------------------------
  SUBROUTINE soilmodel_0(i,j,isoil,iluse,subid,dayno,classarea,classheight,prec,cprec,temp, & 
       daylength,mintemp,maxtemp,sffrac,swrad,  &
       radext,netrad,actvap,satvap,wind,rrcscorr,phoscorr,cevpcorr,incorr,oncorr,  &
       frozenstate,soilstate,miscstate,surfaceflow,csrunoff,crunoffd,   &
       cropuptakein,nitrif,denitrif,epot,gwat,frostdepth,    &
       smdef,evap,cevap,crunoff1,crunoff2,crunoff3,  &
       pwneed,irrappl,irrsources,snowfall,rainfall,cropsources,ruralaload,rgrwload,atmdepload,infiltrationflows,evapflows,  &
       runofflows,crunofflows,verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,cruralflow)  

    INTEGER, INTENT(IN) :: i        !<index for current subbasin
    INTEGER, INTENT(IN) :: j        !<index for current class 
    INTEGER, INTENT(IN) :: isoil    !<index of soil type
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    INTEGER, INTENT(IN) :: dayno    !<pseudo dayno for use in soil model subroutines
    REAL, INTENT(IN) :: classarea   !<class area [km2]
    REAL, INTENT(IN) :: classheight !<elevation [möh]
    REAL, INTENT(IN) :: prec        !<precipitation (mm/timestep)
    REAL, INTENT(IN) :: cprec(numsubstances)        !<concentration of precipitation
    REAL, INTENT(IN) :: temp        !<temperature
    REAL, INTENT(IN) :: daylength   !<day length (hours)
    REAL, INTENT(IN) :: mintemp     !<current daily min temperature (C)
    REAL, INTENT(IN) :: maxtemp     !<current daily max temperature (C)
    REAL, INTENT(IN) :: sffrac      !<snowfall fraction of precipitation [-]
    REAL, INTENT(IN) :: swrad       !<downward shortwave radiation [MJ/m2/day]
    REAL, INTENT(IN) :: radext      !<extraterrestrial solar radiation [MJ/m2/day]
    REAL, INTENT(IN) :: netrad      !<net downward radiation [MJ/m2/day]
    REAL, INTENT(IN) :: actvap      !<actual vapor pressure [kPa]
    REAL, INTENT(IN) :: satvap      !<saturated vapour pressure [kPa]
    REAL, INTENT(IN) :: wind        !<wind speed [m/s]
    REAL, INTENT(IN) :: rrcscorr    !<correction of recession coefficients
    REAL, INTENT(IN) :: phoscorr    !<correction of phosphorus level
    REAL, INTENT(IN) :: cevpcorr    !<correction of potential evaporation
    REAL, INTENT(IN) :: incorr      !<correction of IN
    REAL, INTENT(IN) :: oncorr      !<correction of ON
    TYPE(snowicestatetype),INTENT(INOUT)  :: frozenstate   !<Snow and ice states
    TYPE(soilstatetype),INTENT(INOUT)  :: soilstate   !<Soil states
    TYPE(miscstatetype),INTENT(INOUT)  :: miscstate   !<Misc states
    REAL, INTENT(OUT) :: surfaceflow(2)  !<saturated overflow and surface excess infilt
    REAL, INTENT(OUT) :: csrunoff(numsubstances)   !<concentration surface flow
    REAL, INTENT(OUT) :: crunoffd (numsubstances)  !<concentration tile runoff
    REAL, INTENT(OUT) :: cropuptakein  !<crop uptake of IN      
    REAL, INTENT(OUT) :: nitrif     !<nitrification
    REAL, INTENT(OUT) :: denitrif(maxsoillayers)   !<denitrification
    REAL, INTENT(OUT) :: epot       !<potential evaporation (mm/timestep)
    REAL, INTENT(OUT) :: gwat       !<groundwater table (m)
    REAL, INTENT(OUT) :: frostdepth   !<soil frost depth 
    REAL, INTENT(OUT) :: smdef        !<soil moisture deficit (mm)
    REAL, INTENT(OUT) :: evap       !<evapotranspiration (mm) weighted sum of evap(snowfree)+evapsnow
    REAL, INTENT(OUT) :: cevap(numsubstances)   !<concentration of evapotranspiration
    REAL, INTENT(OUT) :: crunoff1(numsubstances)   !<concentration of runoff from soil layer 1 (mg/L)
    REAL, INTENT(OUT) :: crunoff2(numsubstances)   !<concentration of runoff from soil layer 2 (mg/L)
    REAL, INTENT(OUT) :: crunoff3(numsubstances)   !<concentration of runoff from soil layer 3 (mg/L)
    REAL, INTENT(OUT) :: pwneed                    !<irrigation water demand for this classe (m3)
    REAL, INTENT(INOUT) :: irrappl(2)   !<applied irrigation (mm), for summation basin output
    REAL, INTENT(INOUT) :: irrsources(numsubstances)     !<Load from irrigation to soil (kg/timestep)
    REAL, INTENT(OUT) :: snowfall     !<Precipitation as rain (mm)
    REAL, INTENT(OUT) :: rainfall     !<Precipitation as snow (mm)
    REAL, INTENT(INOUT) :: cropsources(2,numsubstances)  !<Load from fertiliser and resudues (kg/timestep)
    REAL, INTENT(OUT) :: ruralaload(numsubstances)   !<Load from rural households (kg/timestep)
    REAL, INTENT(INOUT) :: rgrwload(numsubstances)     !<Load from regional groundwater flow to soil (kg/timestep)
    REAL, INTENT(INOUT) :: atmdepload(numsubstances)   !<Load of atmospheric dry deposition (kg/timestep)
    REAL, INTENT(OUT) :: infiltrationflows(7)  !<several infiltration flows [mm]
    REAL, INTENT(OUT) :: evapflows(4)  !<evaporation from soillayers, snow and (glacier) [mm]
    REAL, INTENT(OUT) :: runofflows(7) !<different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3,7=saturated surface runoff
    REAL, INTENT(OUT) :: crunofflows(numsubstances,6) !<concentration of different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3
    REAL, INTENT(OUT) :: verticalflows(6) !<vertical flows:1-2=percolation,3-4=upwelling due to rural,5-6=upwelling due to reg. grw flows
    REAL, INTENT(OUT) :: cverticalflows(2,numsubstances) !<concentration of vertical flows:1-2=percolation
    REAL, INTENT(OUT) :: horizontalflows(3)  !<horizontal flows:1-3=recieved rural load flow
    REAL, INTENT(INOUT) :: horizontalflows2(3,nsub) !<horizontal flows:1-3=division of regional groundwater flows to grwdown
    REAL, INTENT(OUT) :: evapsnow !<actual evaporation from snow
    REAL, INTENT(OUT) :: cruralflow(numsubstances) !<concentration of rural flow
    
    !Local variables
    INTEGER nc      !numsubstances
    INTEGER status  !error status of subroutine call
    REAL tt       !threshold temperature for melting (C)
    REAL cm       !coefficient for snow melt (mm/Cday)
    REAL sc       !coefficient for runoff recession surface runoff(no unit)
    REAL oldsnow
    REAL helpmm
    REAL plantuptake(2,2)          !uptake of plant (what they want), kg NP/km2/timestep
    REAL sink(numsubstances), source(numsubstances)

    !Variables for class values
    REAL ginfilt,infilt   !gross infiltration (rain+melt), actual infiltration (after removed surfaceflow and macroporeflow)
    REAL cginfilt(numsubstances),cinfilt(numsubstances)   !concentration of infiltration
    REAL totalsurfaceflow
    REAL satoverflow    !surface flow from saturated soil
    REAL excessinfilt   !infiltration excess surface runoff 
    REAL macroflow
    REAL cmacroflow(numsubstances), cexcessinfilt(numsubstances)          !concentration in infiltration excess runoff and macropore flow
    REAL melt, cmelt(numsubstances)      !snow melt,concentration in snow melt water
    REAL cevapsnow(numsubstances) !concentration of snow evaporation
    REAL soilrunoff(maxsoillayers)    !soil runoff
    REAL csoilrunoff(numsubstances,maxsoillayers)    !concentration of soil runoff
    REAL cweights(maxsoillayers) ! weigths to calculate T2 conc. in tile drainage 
    REAL trunofftemp(maxsoillayers)
    REAL effsnowcov       !effective snow cover with respect to evaporation evap = evap(snowfree)*(1-effsnowcov)+evapsnow
    REAL epotsnow         !potential evaporation from snow
    REAL T1release !release of T1 from surface sources
    REAL runoffd      !tile drainage runoff
   
    !Output, default values
    infiltrationflows = 0.
    runofflows = 0.
    crunofflows = 0.
    evapflows = 0.
    verticalflows=0.
    cverticalflows=0.
    
    !Short notation of parameter values to be used in this subroutine
    tt=landpar(m_ttmp,iluse)       !Threshold temperature for snow melt and evaporation
    cm=landpar(m_cmlt,iluse)       !Coefficient for snow melt
    sc=landpar(m_srrcs,iluse)*rrcscorr      !Runoff coefficient for surface runoff
    IF(sc>1.) sc = 1.

    !Locally defined variables for indata to be used in this subroutine
    T1release = 0.
    nc = numsubstances

    !Atmospheric deposition, add to soil or snow
    CALL add_dry_deposition_to_landclass(i,j,iluse,conductN,conductP,classarea,atmdepload,&
         classdata(j)%vegtype,load(i)%indrydep,landpar(m_drypp,iluse),frozenstate,soilstate)
    
    !Calculate plant growth (for uptake of nutrients)
    CALL calculate_plant(i,j,dayno,temp,daylength,plantuptake,miscstate)

    !Irrigate the soil     
    IF(doirrigation) CALL apply_irrigation(i,j,epotdist(:,j),pwmm(2,j),soilstate,irrappl,irrsources)

    !Potential evapotranspiration (before snow calculations, to calculate snow evaporation)
    CALL calculate_potential_evaporation(i,j,temp,epot,radext,swrad,netrad,actvap,satvap,wind,epotsnow)
    epot = epot * cevpcorr          !PET for snow free soil, regionally adjusted
    epotsnow = epotsnow * cevpcorr  !PET for snow covered soil, regionally adjusted
    
    !Update snow pack; snowfall, melting, evaporation (sublimation)
    CALL calculate_rain_snow_from_precipitation(iluse,prec,temp,snowfall,rainfall,sffrac) !form of precipitation
    oldsnow = frozenstate%snow(j,i)   !for snowdepth calculation
    CALL calculate_snow(i,j,basin(i)%subid,iluse,snowfall,cprec,frozenstate%snow(j,i),   &
              frozenstate%csnow(:,j,i),temp,melt,cmelt,swrad, &
              frozenstate%snowage(j,i),frozenstate%snowmax(j,i),frozenstate%snowdepth(j,i),  &
              frozenstate%snowcov(j,i),epotsnow,evapsnow,cevapsnow,effsnowcov)
    evapflows(3) = evapsnow
    
    !Gross infiltration
    ginfilt  = rainfall + melt 
    IF(numsubstances>0)THEN
      IF(ginfilt>0.)THEN
        cginfilt = (cmelt*melt + cprec*rainfall) / ginfilt
      ELSE 
        cginfilt = 0.
      ENDIF
    ENDIF

    !Calculate and add infiltration to soil, including calculation of surface flow and macropore flow due to limited infiltration capacity 
    CALL infiltration(i,j,isoil,wpmm(:,j),fcmm(:,j),epmm(:,j),ginfilt,cginfilt,temp,mintemp,maxtemp,  &
       infilt,cinfilt,excessinfilt,cexcessinfilt,macroflow,cmacroflow,frozenstate,soilstate)
    CALL add_infiltration(i,j,iluse,infilt,cinfilt,soilstate)
    IF(ginfilt>0.)THEN
      infiltrationflows(1) = melt/ginfilt
    ENDIF
    infiltrationflows(2) = infilt
    infiltrationflows(3) = excessinfilt
    
    !Percolation down through the soil layers, including DOC-reduction
    CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),epmm(:,j),soilthick(:,j),verticalflows(1:2),cverticalflows,soilstate)

    !Calculate and remove evapotranspiration from the soil upper two layers
    CALL calculate_actual_soil_evapotranspiration(i,j,temp,epot,wpmm(:,j),fcmm(:,j),epotdist(:,j),soilstate,evap,evapflows(1:2),cevap,MIN(1.,MAX(1.-effsnowcov,0.)))

    !Total evaporation, and Weighted average concentrations and potential evapotranspiration
    evap = evap + evapsnow
    IF(numsubstances.GT.0 .AND. evap.GT.0.) cevap(:) = (cevap(:)*(evap-evapsnow) + cevapsnow(:)*evapsnow)/evap
    epot = epot * (1.-effsnowcov) + epotsnow * effsnowcov
    
    !Calculate and remove soil runoff
    CALL calculate_soil_runoff(i,j,subid,wpmm(:,j),fcmm(:,j),epmm(:,j), &
              soildepth(:,j),soilthick(:,j),classdata(j)%streamdepth,soilstate,soilrunoff,csoilrunoff)
    runofflows(1:3) = soilrunoff(1:3)
    crunofflows(:,1:3) = csoilrunoff(:,1:3)
    crunoff1 = csoilrunoff(:,1)
    crunoff2 = csoilrunoff(:,2)
    crunoff3 = csoilrunoff(:,3)

    !Calculate and remove runoff by tile or drainage pipe
    CALL calculate_tile_drainage(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),epmm(:,j),   &
              soildepth(:,j),soilthick(:,j),classdata(j)%tiledepth,rrcscorr,soilstate,runoffd, &
              crunoffd,cweights)
    runofflows(4:6) = runoffd*cweights(1:3)
    crunofflows(:,4) = crunoffd*cweights(1)
    crunofflows(:,5) = crunoffd*cweights(2)
    crunofflows(:,6) = crunoffd*cweights(3)

    !Regional lateral groundwater flow from other subbasins
    IF(modeloption(p_deepgroundwater)==1) CALL add_regional_groundwater_flow_to_soil(i,j,classarea,pwmm(:,j),soilstate,rgrwload,verticalflows(5:6),horizontalflows2)
    
    !Load from local diffuse NP-source to the lowest soil layer
    CALL local_diffuse_source(i,j,pwmm(:,j),classarea,soilstate,ruralaload,verticalflows(3:4),horizontalflows(1:3),cruralflow)
    
    !Calculate soil temperature and frost depth
    CALL calculate_soiltemp(maxsoillayers,temp,frozenstate%snowdepth(j,i),genpar(m_deepmem),soilmem(:,j),soilstate%deeptemp(j,i),soilstate%temp(:,j,i))
    helpmm = fcmm(1,j)+fcmm(2,j)+fcmm(3,j)+wpmm(1,j)+wpmm(2,j)+wpmm(3,j)
    IF(soilthick(3,j)>0)THEN
      CALL calculate_frostdepth(helpmm,landpar(m_cfrost,iluse),soilpar(m_sfrost,isoil),   &
            soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i),frostdepth,soilstate%temp(1:2,j,i),soilthick(:,j))
    ELSEIF(soilthick(2,j)>0)THEN
      CALL calculate_frostdepth(helpmm,landpar(m_cfrost,iluse),soilpar(m_sfrost,isoil),   &
            soilstate%water(1,j,i)+soilstate%water(2,j,i),frostdepth,soilstate%temp(1:2,j,i),soilthick(:,j))
    ELSE
      CALL calculate_frostdepth(helpmm,landpar(m_cfrost,iluse),soilpar(m_sfrost,isoil),   &
            soilstate%water(1,j,i),frostdepth,soilstate%temp(1:1,j,i),soilthick(:,j))
    ENDIF

    !Surface runoff from saturated overland flow of uppermost soil layer
    satoverflow = MAX(sc * (soilstate%water(1,j,i)-pwmm(1,j)),0.)
    csrunoff = 0.
    IF(satoverflow > 0.) THEN
      CALL remove_water(soilstate%water(1,j,i),nc,soilstate%conc(:,1,j,i),satoverflow,soilstate%conc(:,1,j,i),status)
      IF(status.NE.0) CALL error_remove_water(errstring(1),subid,i,j)
    ENDIF
    runofflows(7) = satoverflow

    !Total surfaceflow (saturated overland flow and excess infiltration)
    totalsurfaceflow = satoverflow + excessinfilt
    surfaceflow(1) = satoverflow
    surfaceflow(2) = excessinfilt
    IF(totalsurfaceflow > 0. .AND. numsubstances>0) THEN
       csrunoff(:) = (soilstate%conc(:,1,j,i) * satoverflow + cexcessinfilt(:) * excessinfilt) / totalsurfaceflow     !used for satoverflow and excessinfilt
    ENDIF

    !Erosion of particles with fastflow (surface flow and macropore flow) including delay in temporary storage.
    IF(conductP.OR.conductS) CALL particle_processes_for_runoff(i,j,isoil,iluse,dayno,rainfall,totalsurfaceflow,   &
              macroflow,runoffd,runofflows(1)+runofflows(2)+runofflows(3)+runoffd+totalsurfaceflow,phoscorr,  &
              csrunoff,cmacroflow,crunoffd,crunoff1,    &
              crunoff2,crunoff3,frozenstate%snow(j,i),soilstate)  

    !Add macropore water to soil layer with groundwater level (except the PP)
    CALL add_macropore_flow(i,j,macroflow,cmacroflow,wpmm(:,j),fcmm(:,j), &
            epmm(:,j),pwmm(:,j),soildepth(:,j),soilthick(:,j),infiltrationflows(4:6),soilstate)

    !Second percolation down through the soil layers, including DOC-reduction (limited to same maxperc)
    CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),epmm(:,j),soilthick(:,j),verticalflows(1:2),cverticalflows,soilstate)

    !Groundwater level and soil moisture deficit for output variable
    CALL calculate_groundwater_table(soilstate%water(:,j,i),wpmm(:,j),  &
            fcmm(:,j),epmm(:,j),pwmm(:,j),soildepth(:,j),soilthick(:,j),gwat) 
    CALL calculate_soil_moisture_deficit(soilstate%water(:,j,i),wpmm(:,j),  &
            fcmm(:,j),soilthick(:,j),smdef) 

    !Soil transformation processes for substances             
    CALL soil_np_processes(i,j,iluse,conductN,conductP,conductC,dayno,classarea,wpmm(:,j),  &
         fcmm(:,j),epmm(:,j),plantuptake,soilthick(:,j),genpar(m_fertdays),genpar(m_littdays),  &
         source,sink,nitrif,denitrif,cropuptakein,cropsources,   &
         landpar(m_dissolfN,iluse),landpar(m_dissolfP,iluse),   &
         oncorr * landpar(m_dissolhN,iluse),phoscorr * landpar(m_dissolhP,iluse),      &
         landpar(m_minerfn,iluse),landpar(m_minerfp,iluse),incorr * landpar(m_degradhn,iluse),      &
         (2.-incorr) * landpar(m_denitrlu,iluse),landpar(m_degradhp,iluse),soilstate)
    IF(conductC) CALL soil_carbon_processes(i,j,wpmm(:,j),fcmm(:,j),epmm(:,j),pwmm(:,j),soilthick(:,j),    &
         genpar(m_crate1),genpar(m_crate2),genpar(m_crate3),genpar(m_crate9),genpar(m_crate10),genpar(m_minc),  &
         landpar(m_ocsoim,iluse),landpar(m_ocsmslp,iluse),soilstate)
    CALL balance_spsoil(i,j,conductP,soilthick(:,j),soilpar(m_freuc,isoil), &
         soilpar(m_freuexp,isoil),soilpar(m_freurate,isoil),soilstate)
    CALL soil_tracer_processes(i,j,currentdate%year,soilthick(:,j),soilstate,miscstate,ginfilt,T1release)
    
    !Add released T1 to surface runoff or top soil
    IF(i_t1>0)THEN
      IF(T1release>0.)THEN
        IF(totalsurfaceflow>0.)THEN
          csrunoff(i_t1) = csrunoff(i_t1) + (T1release * (totalsurfaceflow/(infilt+totalsurfaceflow))) / totalsurfaceflow
        ENDIF
        soilstate%partT1(1,j,i) = soilstate%partT1(1,j,i) + T1release * (infilt/(infilt+totalsurfaceflow))
      ENDIF
    ENDIF      

    !Calculate irrigation water demand (for next timestep)
    IF(doirrigation) CALL calculate_irrigation_water_demand(i,j,dayno,  &
         classarea,genpar(m_sswcorr),genpar(m_immdep),genpar(m_iwdfrac),  &
         genpar(m_wdpar),soilstate%water(:,j,i),wpmm(:,j),fcmm(:,j),  &
         epmm(:,j),epot,epotdist(:,j),pwneed)

    !Riparian zone for OC
    IF(conductC) CALL class_riparian_zone_processes(numsubstances,classheight,runofflows(1)+runofflows(2)+runofflows(3), &
         crunoff1(:),crunoff2(:),crunoff3(:),soilstate%oldgrw(j,i),landpar(m_ripz,iluse),        &
         soilstate%temp(:,j,i),genpar(m_ripe),genpar(m_rips),gwat,miscstate%temp10(i),miscstate%temp20(i),     &
         soilstate%water(1,j,i)+soilstate%water(2,j,i)+soilstate%water(3,j,i),SUM(wpmm(:,j)),SUM(pwmm(:,j)),soildepth(maxsoillayers,j))
         
    !Runoff Temperature concentrations dependent on soiltemp.calculation [DG/JS Temp.model, May-2013]
    IF(i_t2>0) THEN
        trunofftemp(1) = amax1(0.,soilstate%temp(1,j,i))
        trunofftemp(2) = amax1(0.,soilstate%temp(2,j,i))
        trunofftemp(3) = amax1(0.,soilstate%temp(3,j,i))
        
        !T2 conc. in soil layers
        soilstate%conc(i_t2,1,j,i) = trunofftemp(1)
        soilstate%conc(i_t2,2,j,i) = trunofftemp(2)
        soilstate%conc(i_t2,3,j,i) = trunofftemp(3)
    
        !T2 conc. in runoff from soil layers
        crunoff1(i_t2) = trunofftemp(1)
        crunoff2(i_t2) = trunofftemp(2)
        crunoff3(i_t2) = trunofftemp(3)
        
        !T2 conc. in surface runoff, mixture of excess infiltration and saturated topsoil
        IF(totalsurfaceflow>0)THEN
          csrunoff(i_t2) = (trunofftemp(1) * satoverflow + excessinfilt * cexcessinfilt(i_t2)) / totalsurfaceflow 
        ELSE
          csrunoff(i_t2) = 0.0
        ENDIF    
        !T2 conc. in tile drainage (weigthed average over tile drainage depth)
        crunoffd(i_t2) = cweights(1) * trunofftemp(1) + cweights(2) * trunofftemp(2) + &
                         cweights(3) * trunofftemp(3)
    ENDIF

  END SUBROUTINE soilmodel_0


END MODULE SOILMODEL_DEFAULT
