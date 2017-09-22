!> \file soilmodel4.f90
!> Contains module floodplain_soilmodel.

!>HYPE soil model for floodplain areas
MODULE FLOODPLAIN_SOILMODEL

  !Copyright 2015-2017 SMHI
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
  ! soilmodel_4
  !-----------------------------------------------------------------------------------------

  USE STATETYPE_MODULE
  USE MODVAR, ONLY: missing_value,pi,   &
       currentdate, &
       genpar,landpar,soilpar,  &
       load,basin,path, &
       classdata,classbasin,  &
       numsubstances,nsub,maxsoillayers, &
       soilthick,soildepth, &
       i_t1,i_t2,i_in,i_on,i_sp,i_pp,i_oc,&
       modeloption,p_floodplain, &
       conductN,conductP,conductC,conductS,  &
       flooding,floodindex, &
       slc_mriver,slc_olake, &
       doirrigation
  USE HYPEVARIABLES, ONLY : wpmm,fcmm,epmm,pwmm, &
       m_ttmp,m_cmlt,m_srrcs,m_mactrinf, &
       m_wetsp,m_drypp, m_ponatm,m_cfrost,m_sfrost,m_perc1,m_perc2,  &
       m_soilcoh,m_soilerod,m_sreroexp,m_filtpbuf,m_filtpinner, &
       m_filtpother,m_macfilt,m_pprelmax,m_pprelexp,m_fertdays, &
       m_minerfn,m_minerfp,m_degradhn,m_degradhp,m_denitrlu, &
       m_dissolfN,m_dissolhN,m_dissolfP,m_dissolhP,m_littdays,  &
       m_crate1,m_crate2,m_crate3,m_crate9,m_crate10,m_minc,  &
       m_freuc,m_freuexp,m_freurate,m_sswcorr,m_immdep,m_iwdfrac, &
       m_wdpar,m_ripz,m_rips,m_ripe,m_ocsoim,m_ocsmslp,   &
       m_opt7,m_optonoff,  &
       soilmem,m_deepmem, &
       epotdist
  USE GENERAL_WATER_CONCENTRATION, ONLY : remove_water, &
                                          error_remove_water, &
                                          inflow_lowest_soillayer, &
                                          add_water
  USE ATMOSPHERIC_PROCESSES, ONLY :  calculate_rain_snow_from_precipitation
  USE SOIL_PROCESSES, ONLY : calculate_snow, &
                             calculate_potential_evaporation,          &
                             calculate_actual_soil_evapotranspiration, &
                             add_macropore_flow,    &
                             calculate_tile_drainage,  &
                             calculate_soil_runoff,    &
                             infiltration,   &
                             add_infiltration,   &
                             flood_infiltration,   &
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
  USE REGIONAL_GROUNDWATER_MODULE, ONLY : add_regional_groundwater_flow_to_soil
  USE IRRIGATION_MODULE, ONLY : apply_irrigation,              &
                                calculate_irrigation_water_demand
  USE SURFACEWATER_PROCESSES, ONLY : add_precipitation_to_floodplain, &
                                     calculate_floodplain_evaporation, &
                                     calculate_floodplain_waterlevel

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: soilmodel_4
  
  !Private parameters, global in this module
  CHARACTER(LEN=80) :: errstring(2)  !error message for location of remove_water call
  PARAMETER (errstring = (/'surface runoff, soillayer 1   ', &
                           'infiltration from flooded area'/))

CONTAINS

  !>\brief Soilmodel for floodplain
  !----------------------------------------------------------------
  SUBROUTINE soilmodel_4(i,j,isoil,iluse,subid,dayno,classarea,prec,cprec,temp, & 
       daylength,mintemp,maxtemp,sffrac,swrad,  &
       radext,netrad,actvap,satvap,wind,rrcscorr,phoscorr,cevpcorr,incorr,oncorr,  &
       frozenstate,soilstate,miscstate,surfaceflow,csrunoff,   &
       crunoffd,cropuptakein,nitrif,denitrif,epot,gwat,frostdepth,    &
       smdef,evap,cevap,crunoff1,crunoff2,crunoff3,snowcov,nonflood_part,pwneed,irrappl,irrsources,  &
       snowfall,rainfall,cropsources,ruralaload,atmdepload1,atmdepload,infiltrationflows,floodplainflows,evapflows,  &
       runofflows,verticalflows,cverticalflows,horizontalflows,horizontalflows2,evapsnow,totalsoilrunoff,cruralflow)  

    INTEGER, INTENT(IN) :: i        !<index for current subbasin
    INTEGER, INTENT(IN) :: j        !<index for current class 
    INTEGER, INTENT(IN) :: isoil    !<index of soil type
    INTEGER, INTENT(IN) :: iluse    !<index of landuse
    INTEGER, INTENT(IN) :: subid    !<subbasin id
    INTEGER, INTENT(IN) :: dayno    !<pseudo dayno for use in soil model subroutines
    REAL, INTENT(IN) :: classarea   !<class area [km2]
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
    REAL, INTENT(OUT) :: snowcov      !<snow cover (based on classfparea)
    REAL, INTENT(OUT) :: nonflood_part    !<(old) fraction of floodplain is non-flooded (-)
    REAL, INTENT(OUT) :: pwneed                    !<irrigation water demand for this classe (m3)
    REAL, INTENT(INOUT) :: irrappl(2)              !<applied irrigation (mm), for summation basin output
    REAL, INTENT(INOUT) :: irrsources(numsubstances) !<Load from irrigation to soil (kg/timestep)
    REAL, INTENT(OUT) :: snowfall     !<Precipitation as rain (mm)
    REAL, INTENT(OUT) :: rainfall     !<Precipitation as snow (mm)
    REAL, INTENT(INOUT):: cropsources(2,numsubstances)  !<Load from fertiliser and resudues (kg/timestep)
    REAL, INTENT(OUT):: ruralaload(numsubstances)   !<Load from rural households (kg/timestep)
    REAL, INTENT(INOUT):: atmdepload1(numsubstances)  !<Load of atmospheric wet deposition on flooded flodplain only (kg/timestep)
    REAL, INTENT(INOUT):: atmdepload(numsubstances)   !<Load of atmospheric dry deposition on flooded flodplain only (kg/timestep)
    REAL, INTENT(OUT) :: infiltrationflows(7)  !<several infiltration flows [mm]
    REAL, INTENT(OUT) :: floodplainflows(3)   !<snow incorporation to floodplain, precipitation on floodplain,infiltration [m3]
    REAL, INTENT(OUT) :: evapflows(4)  !<evaporation from soillayers (1-2), snow(3) and flooded floodplain(4) [mm/m3]
    REAL, INTENT(OUT) :: runofflows(7) !<different runoff flows:1-3=soil runoff sl 1-3,4-6=tile runoff sl 1-3,7=saturated surface runoff
    REAL, INTENT(OUT) :: verticalflows(6) !<vertical flows:1-2=percolation,3-4=upwelling due to rural,5-6=upwelling due to reg. grw flows
    REAL, INTENT(OUT) :: cverticalflows(2,numsubstances) !<concentration of vertical flows:1-2=percolation
    REAL, INTENT(OUT) :: horizontalflows(3)  !<horizontal flows:1-3=recieved rural load flow
    REAL, INTENT(INOUT) :: horizontalflows2(3,nsub) !<horizontal flows:1-3=division of regional groundwater flows to grwdown
    REAL, INTENT(OUT) :: evapsnow   !<actual evaporation from snow (mm) for classfparea
    REAL, INTENT(OUT) :: totalsoilrunoff(2)  !<total runoff from flood plain soil to waterbody:1=main river, 2=olake
    REAL, INTENT(OUT) :: cruralflow(numsubstances) !<conc of rural flow
    
    !Local variables
    INTEGER nc      !numsubstances
    INTEGER status  !error status of subroutine call
    INTEGER watertype   !water body in misctype (1=main river, 2= olake)
    REAL tt       !threshold temperature for melting (C)
    REAL cm       !coefficient for snow melt (mm/Cday)
    REAL sc       !coefficient for runoff recession surface runoff(no unit)
    REAL snowonfp, snowmaxonfp,oldsnowonfp, snowdepthonfp  !snow rescaled to non-flooded floodplain area (mm)
    REAL helpmm,helparea
    REAL plantuptake(2,2)          !uptake of plant (what they want), kg NP/km2/timestep
    REAL sink(numsubstances), source(numsubstances)
    REAL fpfrac      !fraction of classarea that is floodplain
    REAL trans

    !Variables for class values
    REAL ginfilt,infilt   !gross infiltration (rain+melt), actual infiltration (after removed surfaceflow and macroporeflow)
    REAL cginfilt(numsubstances),cinfilt(numsubstances)   !concentration of infiltration
    REAL infiltfp   !infiltration from flooded floodplain
    REAL totalsurfaceflow
    REAL satoverflow    !surface flow from saturated soil
    REAL excessinfilt   !infiltration excess surface runoff 
    REAL macroflow
    REAL cmacroflow(numsubstances), cexcessinfilt(numsubstances)          !concentration in infiltration excess runoff and macropore flow
    REAL melt,cmelt(numsubstances)    !snow melt
    REAL soilrunoff(maxsoillayers)    !soil runoff
    REAL csoilrunoff(numsubstances,maxsoillayers)    !concentration of soil runoff
    REAL runoffd    !tile drainage runoff
    REAL cweights(maxsoillayers) ! weigths to calculate T2 conc. in tile drainage 
    REAL trunofftemp(maxsoillayers)
    REAL netaddtosoil,addedtosoil
    REAL effsnowcov       !effective snow cover with respect to evaporation evap = evap(snowfree)*(1-effsnowcov)+evapsnow
    REAL epotsnow         !potential evaporation from snow
    REAL cevapsnow(numsubstances)  !conc in evaporation from snow covered land
    REAL epotsoilwater    !potential evaporation for soil and water surface
    REAL evapsoil    !evaporation from soil [mm]
    REAL classfparea  !Area of floodplain part of class [m2]
    REAL evaprfp,cevaprfp(numsubstances)
    REAL ffpwl    !current water level on floodplain [m]
    REAL ffparea, ffpafrac  !current flooded area of floodplain [m2] and fraction
    REAL fpmaxlevel
    REAL ctotalsoilrunoff(numsubstances)  ! concentration of total runoff from flood plain soil to current waterbody
    REAL T1release !release of T1 from surface sources
    
    !Output, default values
    infiltrationflows = 0.
    floodplainflows = 0.
    runofflows = 0.
    epot = 0.
    evap = 0.
    cevap = 0.
    evapsnow = 0.; cevapsnow = 0.
    runoffd = 0.; crunoffd = 0.
    cropuptakein = 0.
    nitrif = 0.; denitrif=0.
    gwat=0.; smdef=0.
    frostdepth = 0.
    horizontalflows=0.
    verticalflows=0.
    cverticalflows=0.
    evapflows = 0.
    snowcov = 0.
    atmdepload1 = 0.
    
    !Floodplain initiation
    totalsoilrunoff = 0.
    netaddtosoil = 0.
    totalsurfaceflow = 0.
    surfaceflow = 0.
    
    IF(j==slc_mriver)THEN
      watertype = 1
      fpfrac = flooding(floodindex(i))%fpfmr !floodplain fraction of the main river area
    ELSEIF(j==slc_olake)THEN
      watertype = 2
      fpfrac = flooding(floodindex(i))%fpfol !floodplain fraction of the olake area
    ENDIF
    classfparea = classbasin(i,j)%part * basin(i)%area * fpfrac      !Area of floodplain part of class [m2]
    
    !Short notation of parameter values to be used in this subroutine
    tt=landpar(m_ttmp,iluse)       !Threshold temperature for snow melt and evaporation
    cm=landpar(m_cmlt,iluse)       !Coefficient for snow melt
    sc=landpar(m_srrcs,iluse)*rrcscorr      !Runoff coefficient for surface runoff
    IF(sc>1.) sc = 1.
    !Get current interflow parameter value
    IF(genpar(m_optonoff).LE.0)THEN
      IF(watertype==1)THEN
        fpmaxlevel = flooding(floodindex(i))%fymmr
      ELSEIF(watertype==2)THEN
        fpmaxlevel = flooding(floodindex(i))%fymol
      ENDIF
    ELSE
      fpmaxlevel = genpar(m_opt7)
    ENDIF

    !Locally defined variables for indata to be used in this subroutine
    T1release = 0.
    nc = numsubstances

    !Calculation of common processes, flooded or not flooded floodplain soil
    !-----------------------------------------------------------------------
    !Atmospheric deposition, add to soil or snow
    CALL add_dry_deposition_to_landclass(i,j,iluse,conductN,conductP,classfparea*1.E-6,atmdepload,&
         classdata(j)%vegtype,load(i)%indrydep,landpar(m_drypp,iluse),frozenstate,soilstate)
    
    !Calculate plant growth (for uptake of nutrients)
    CALL calculate_plant(i,j,dayno,temp,daylength,plantuptake,miscstate)

    !Irrigate the soil     
    IF(doirrigation) CALL apply_irrigation(i,j,epotdist(:,j),pwmm(2,j),soilstate,irrappl,irrsources)

    !Potential evapotranspiration
    CALL calculate_potential_evaporation(i,j,temp,epotsoilwater,radext,swrad,netrad,actvap,satvap,wind,epotsnow)
    epotsoilwater = epotsoilwater * cevpcorr          !PET for snow free soil and water surface, regionally adjusted
    epotsnow = epotsnow * cevpcorr  !PET for snow covered soil, regionally adjusted
    
    !Calculate current area flooded
    IF(miscstate%floodwater(watertype,i)>0)THEN
      CALL calculate_floodplain_waterlevel(miscstate%floodwater(watertype,i),classfparea,fpmaxlevel,ffpwl,ffparea)
      ffpafrac = ffparea/classfparea
    ELSE
      ffparea = 0.
      ffpafrac = 0.
    ENDIF
    nonflood_part = 1. - ffpafrac
    
    !Calculations for non-flooded floodplain soil
    !--------------------------------------------
    IF(ffparea<classfparea)THEN

      !Scale snow to non-flooded flodplain area for class calculations
      snowonfp = frozenstate%snow(j,i)/nonflood_part
      snowdepthonfp = frozenstate%snowdepth(j,i)/nonflood_part
      snowmaxonfp = frozenstate%snowmax(j,i)/nonflood_part
      oldsnowonfp = snowonfp
      
      !Update snow pack; snowfall, melting, evaporation (sublimation)
      IF(modeloption(p_floodplain)>=2)THEN
        CALL calculate_rain_snow_from_precipitation(iluse,prec,temp,snowfall,rainfall,sffrac) !form of precipitation
        CALL calculate_snow(i,j,basin(i)%subid,iluse,snowfall,cprec,snowonfp,   &
                frozenstate%csnow(:,j,i),temp,melt,cmelt,swrad,frozenstate%snowage(j,i),  &
                snowmaxonfp,snowdepthonfp,frozenstate%snowcov(j,i),epotsnow,evapsnow,cevapsnow,effsnowcov)
        evapflows(3) = evapsnow
        evapsnow = evapsnow*nonflood_part    !rescaled to classfparea
      ELSEIF(modeloption(p_floodplain)==1)THEN
        snowfall =0; rainfall=0; melt=0.; evapsnow=0.; effsnowcov=0.
        !Add more here for snowcov,snowdepth?
      ENDIF
    
      !Gross infiltration
      ginfilt  = rainfall + melt 
      IF(numsubstances>0)THEN
        IF(ginfilt>0.)THEN
          cginfilt(:) = (cmelt(:)*melt + cprec(:)*rainfall) / ginfilt
        ELSE 
          cginfilt(:) = 0.
        ENDIF
      ENDIF

      !Calculate and add infiltration to soil, including calculation of surface flow and macropore flow due to limited infiltration capacity 
      CALL infiltration(i,j,isoil,wpmm(:,j),fcmm(:,j),epmm(:,j),ginfilt,cginfilt,temp,mintemp,maxtemp,  &
         infilt,cinfilt,excessinfilt,cexcessinfilt,macroflow,cmacroflow,frozenstate,soilstate)
      CALL add_infiltration(i,j,iluse,infilt*nonflood_part,cinfilt,soilstate)
      IF(ginfilt>0.)THEN
        infiltrationflows(1) = melt/ginfilt
      ENDIF
      infiltrationflows(2) = infilt*nonflood_part
      infiltrationflows(3) = excessinfilt*nonflood_part

      !Percolation down through the soil layers, including DOC-reduction
      CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),epmm(:,j),soilthick(:,j),verticalflows(1:2),cverticalflows,soilstate)

      !Calculate and remove evapotranspiration from the soil upper two layers
      IF(modeloption(p_floodplain)>=2)THEN
        CALL calculate_actual_soil_evapotranspiration(i,j,temp,epotsoilwater,wpmm(:,j),  &
              fcmm(:,j),epotdist(:,j),soilstate,evapsoil,evapflows(1:2),cevap,  &
              MIN(nonflood_part,MAX(nonflood_part*(1-effsnowcov),0.)))
      ELSEIF(modeloption(p_floodplain)==1)THEN
        evapsoil =0.; cevap = 0.
      ENDIF
      
      !Accumulate evaporation and weighted average concentrations and potential evapotranspiration for non-flooded part
      IF(numsubstances.GT.0) cevap(:) = (cevap(:)*evapsoil + cevapsnow(:)*evapsnow)
      evap = evapsoil + evapsnow
      epot = (epotsoilwater * (1.-effsnowcov) + epotsnow * effsnowcov)*nonflood_part
    
      !Calculate and remove soil runoff
      CALL calculate_soil_runoff(i,j,subid,wpmm(:,j),fcmm(:,j),epmm(:,j), &
              soildepth(:,j),soilthick(:,j),nonflood_part*classdata(j)%streamdepth, &
              soilstate,soilrunoff,csoilrunoff)
      runofflows(1:3) = soilrunoff(1:3)
      crunoff1 = csoilrunoff(:,1)
      crunoff2 = csoilrunoff(:,2)
      crunoff3 = csoilrunoff(:,3)

      !Calculate and remove runoff by tile or drainage pipe
      CALL calculate_tile_drainage(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),epmm(:,j),   &
              soildepth(:,j),soilthick(:,j),nonflood_part*classdata(j)%tiledepth,rrcscorr,soilstate,runoffd, &
              crunoffd,cweights)
      runofflows(4:6) = runoffd*cweights(1:3)
      
      !Add load from local diffuse sources to the lowest soil layer
      CALL local_diffuse_source(i,j,pwmm(:,j),classfparea*1.E-6,soilstate,ruralaload,verticalflows(3:4),horizontalflows(1:3),cruralflow)
      
      !Calculate soil temperature, snow age, snow depth and frost depth
      CALL calculate_soiltemp(maxsoillayers,temp,snowdepthonfp,genpar(m_deepmem),soilmem(:,j),soilstate%deeptemp(j,i),soilstate%temp(:,j,i))
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

      !!Percolation down through the soil layers, including DOC-reduction
      !CALL percolation(i,j,isoil,subid,wpmm(:,j),fcmm(:,j),epmm(:,j),soilthick(:,j),verticalflows(1:2),cverticalflows,soilstate)
      !
      !Surface runoff from saturated overland flow of uppermost soil layer, (sc is set to 1 in par.txt)
      satoverflow = MAX(sc * (soilstate%water(1,j,i)-pwmm(1,j)),0.)
      csrunoff = 0.
      IF(satoverflow > 0.) THEN
        CALL remove_water(soilstate%water(1,j,i),nc,soilstate%conc(:,1,j,i),satoverflow,soilstate%conc(:,1,j,i),status)
        IF(status.NE.0) CALL error_remove_water(errstring(1),subid,i,j)
      ENDIF
      runofflows(7) = satoverflow

      !Total surfaceflow (saturated overland flow and excess infiltration)
      totalsurfaceflow = satoverflow + excessinfilt*nonflood_part
      surfaceflow(1) = satoverflow
      surfaceflow(2) = excessinfilt*nonflood_part
      IF(totalsurfaceflow > 0. .AND. numsubstances>0) THEN
        csrunoff(:) = (soilstate%conc(:,1,j,i) * satoverflow + excessinfilt * cexcessinfilt(:)*nonflood_part) / totalsurfaceflow     !used for satoverflow and excessinfilt
      ENDIF

      !Erosion of particles with fastflow (surface flow and macropore flow) including delay in temporary storage.
      IF(conductP.OR.conductS) CALL particle_processes_for_runoff(i,j,isoil,iluse,dayno,rainfall,totalsurfaceflow,   &
              macroflow,runoffd,runofflows(1)+runofflows(2)+runofflows(3)+runoffd+totalsurfaceflow,phoscorr,  &
              csrunoff,cmacroflow,crunoffd,crunoff1,    &
              crunoff2,crunoff3,snowonfp,soilstate)  

      !Add macropore water to soil layer with groundwater level (except the PP)
      CALL add_macropore_flow(i,j,macroflow*nonflood_part,cmacroflow,wpmm(:,j),fcmm(:,j), &
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
    
      !Add released T1 to surface runoff or top soil, all surfacepool released T1 handled here on nonflooded area
      IF(i_t1>0)THEN
        IF(T1release>0.)THEN
          IF(totalsurfaceflow>0.)THEN
            csrunoff(i_t1) = csrunoff(i_t1) + (T1release * (totalsurfaceflow/(infilt*nonflood_part+totalsurfaceflow))) / totalsurfaceflow
          ENDIF
          soilstate%partT1(1,j,i) = soilstate%partT1(1,j,i) + T1release * (infilt*nonflood_part/(infilt*nonflood_part+totalsurfaceflow))
        ENDIF
      ENDIF      

      !Calculate irrigation water demand (for next timestep) only if no flooded area
      IF(ffpafrac==0.)THEN
        IF(doirrigation) CALL calculate_irrigation_water_demand(i,j,dayno,  &
           classfparea,genpar(m_sswcorr),genpar(m_immdep),genpar(m_iwdfrac),  &
           genpar(m_wdpar),soilstate%water(:,j,i),wpmm(:,j),fcmm(:,j),  &
           epmm(:,j),epot,epotdist(:,j),pwneed)
      ENDIF
      !No riparian zone for OC on floodplain
      
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

      !Flow from floodplain soil to (river/lake or) floodplain water body [m3] ?
      totalsoilrunoff(watertype) = SUM(soilrunoff)*classfparea*1.E-3
      IF(nc>0.) ctotalsoilrunoff = (soilrunoff(1)*crunoff1 + soilrunoff(2)*crunoff2 + soilrunoff(3)*crunoff3)/(SUM(soilrunoff))

      !Scale snow back to class floodplain area (for subbasin output)
      frozenstate%snow(j,i) = snowonfp*nonflood_part
      frozenstate%snowmax(j,i) = snowmaxonfp*nonflood_part
      frozenstate%snowdepth(j,i) = snowdepthonfp*nonflood_part
      snowcov = frozenstate%snowcov(j,i)*nonflood_part
      
    ELSE
    
      !Add any remaining snow to floodwater
      IF(frozenstate%snow(j,i)>0.)THEN
        floodplainflows(1) = frozenstate%snow(j,i)*classfparea*1.E-3
        miscstate%floodwater(watertype,i) = miscstate%floodwater(watertype,i) + floodplainflows(1) !m3
        frozenstate%snow(j,i) = 0.
        frozenstate%snowdepth(j,i) = 0.
        frozenstate%snowcov(j,i) = 0.
        frozenstate%snowmax(j,i) = 0.
        frozenstate%snowage(j,i) = 0.
      ENDIF
      
    ENDIF

    !Calculate flooded floodplain
    !----------------------------
    IF(miscstate%floodwater(watertype,i)>0)THEN
      
      !Add precipitation to flooded floodplain
      IF(prec>0)THEN
        IF(modeloption(p_floodplain)>=2)THEN
          helparea = ffparea
        ELSEIF(modeloption(p_floodplain)==1)THEN
          helparea = classfparea
        ENDIF
        CALL add_precipitation_to_floodplain(i,watertype,helparea,prec,cprec,miscstate,atmdepload1)
        floodplainflows(2) = prec*helparea*1.E-3
      ENDIF
      atmdepload = 0.   !dry dep is added to soil directly
    
      !Evaporation flooded floodplain
      evaprfp = 0.  !Floodplain evaporation
      IF(epotsoilwater>0.)THEN
        IF(modeloption(p_floodplain)>=2)THEN
          helparea = ffparea
        ELSEIF(modeloption(p_floodplain)==1)THEN
          helparea = classfparea
        ENDIF
          CALL calculate_floodplain_evaporation(i,j,watertype,numsubstances,    &
                helparea,temp,epotsoilwater,evaprfp,cevaprfp,miscstate)
        evapflows(4) = evaprfp*helparea*1.E-3   !m3
      !Accumulated for floodplain
      IF(numsubstances.GT.0) cevap = cevap*evap + cevaprfp*evaprfp*helparea/classfparea
      evap = evap + evaprfp * helparea/classfparea
      epot = epot + epotsoilwater * helparea/classfparea
      ENDIF

      !TODO: add T2 (and ice?) processes for floodplains
      !TODO: add NP and C processes in flooded water

      !Infiltration from flooded floodplain into soil (hur ska jag beräkna detta flöde?)
      IF(modeloption(p_floodplain)>=2)THEN
        infiltfp = soilpar(m_mactrinf,isoil)*ffparea*1.E-3   !m3, maximum
        trans=1/ffparea*1.E3  !m3_mm
        IF(miscstate%floodwater(watertype,i)>infiltfp)THEN
          netaddtosoil = infiltfp*trans*ffpafrac !mm floodplain area
        ELSE
          infiltfp = miscstate%floodwater(watertype,i)
          netaddtosoil = miscstate%floodwater(watertype,i)*trans*ffpafrac  !mm
        ENDIF
          CALL flood_infiltration(i,j,pwmm(1,j),netaddtosoil,addedtosoil,miscstate%cfloodwater(:,watertype,i),soilstate) 
          IF(netaddtosoil>addedtosoil) infiltfp = addedtosoil/trans/ffpafrac
        CALL remove_water(miscstate%floodwater(watertype,i),nc,miscstate%cfloodwater(:,watertype,i),infiltfp,miscstate%cfloodwater(:,watertype,i),status)
        IF(status.NE.0) CALL error_remove_water(errstring(2),subid,i,watertype)
          floodplainflows(3) = infiltfp
      ENDIF
    ENDIF
    
    !Add runoff from floodplain soil to flooded floodplain waterbody
    IF(totalsurfaceflow>0.)THEN
      CALL add_water(nc,miscstate%floodwater(watertype,i),miscstate%cfloodwater(:,watertype,i),totalsurfaceflow*1.E-3*classfparea,csrunoff)
    ENDIF
    IF(totalsoilrunoff(watertype)>0.)THEN  !maybe move this one to river/lakevolume instead, it is a loop!
      CALL add_water(nc,miscstate%floodwater(watertype,i),miscstate%cfloodwater(:,watertype,i),totalsoilrunoff(watertype),ctotalsoilrunoff)
    ENDIF
    
    !Finalise output concentration average variables
    IF(numsubstances.GT.0 .AND. evap>0.) cevap = cevap/evap

    
  END SUBROUTINE soilmodel_4


END MODULE FLOODPLAIN_SOILMODEL
