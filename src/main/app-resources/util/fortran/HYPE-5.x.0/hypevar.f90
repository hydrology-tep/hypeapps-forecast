!> \file hypevar.f90
!> Contains module hypevariables.

!>Variables for the HYPE model (HYdrological Predictions for the Environment)
!!
MODULE HYPEVARIABLES

!Copyright 2011-2017 SMHI
!
!This file is part of HYPE.
!HYPE is free software: you can redistribute it and/or modify it under the terms of the Lesser GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!HYPE is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Lesser GNU General Public License for more details.
!You should have received a copy of the Lesser GNU General Public License along with HYPE. If not, see <http://www.gnu.org/licenses/>.

!-----------------------------------------------------------------------------------------
!Used by modules
!-----------------------------------------------------------------------------------------
!modelmodule 
!np_processes_module
!oc_processes_module
!irrigation_module
!soil_npc_processes
!soil_processes
!soilmodel_default
!glacier_soilmodel
!-----------------------------------------------------------------------------------------

  SAVE

!General Fortran parameter for HYPE

!> \name Output variable indices
!> \{
  INTEGER,PARAMETER :: o_crun    = 1    !class
  INTEGER,PARAMETER :: o_rrun    = 2    !<Output variable index
  INTEGER,PARAMETER :: o_prec    = 3    !
  INTEGER,PARAMETER :: o_tobs    = 4
  INTEGER,PARAMETER :: o_crunT1  = 5    !class
  INTEGER,PARAMETER :: o_crunT2  = 6    !class
  INTEGER,PARAMETER :: o_roum    = 7    !<Output variable index
  INTEGER,PARAMETER :: o_cprecT1 = 8    !<Output variable index
  INTEGER,PARAMETER :: o_cevapT1 = 9    !class
  INTEGER,PARAMETER :: o_soim    = 10   !<Output variable index
  INTEGER,PARAMETER :: o_csoilT1 = 11   !class
  INTEGER,PARAMETER :: o_csoilT2 = 12   !class
  INTEGER,PARAMETER :: o_roub    = 13   !<Output variable index
  INTEGER,PARAMETER :: o_snow    = 14   !class
  INTEGER,PARAMETER :: o_evap    = 15   !<Output variable index
  INTEGER,PARAMETER :: o_snowdens = 18  !class
  INTEGER,PARAMETER :: o_grwlevel = 19  !class
  INTEGER,PARAMETER :: o_crunIN  = 20   !class
  INTEGER,PARAMETER :: o_crunON  = 21   !class
  INTEGER,PARAMETER :: o_crunSP  = 22   !class
  INTEGER,PARAMETER :: o_crunPP  = 23   !class
  INTEGER,PARAMETER :: o_reTN    = 28   !subbasin
  INTEGER,PARAMETER :: o_reTP    = 29   !subbasin
  INTEGER,PARAMETER :: o_csoilIN = 36   !class
  INTEGER,PARAMETER :: o_soilfrost = 37   !class
  INTEGER,PARAMETER :: o_soiltmp = 38   !class
  INTEGER,PARAMETER :: o_snowdepth = 39   !class
  INTEGER,PARAMETER :: o_epot    = 40   !class
  INTEGER,PARAMETER :: o_reepot  = 41   !<Output variable index
  INTEGER,PARAMETER :: o_cprc    = 43   !class
  INTEGER,PARAMETER :: o_crunOC  = 44   !class
  INTEGER,PARAMETER :: o_csoilOC = 45   !class
  INTEGER,PARAMETER :: o_wcom    = 51   !outlet lake
  INTEGER,PARAMETER :: o_rewstr  = 52   !<Output variable index
  INTEGER,PARAMETER :: o_cout    = 53   !subbasin
  INTEGER,PARAMETER :: o_rout    = 54   !subbasin
  INTEGER,PARAMETER :: o_crunTN  = 66   !class
  INTEGER,PARAMETER :: o_crunTP  = 67   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_pfN    = (/68,69,70/)   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_phN    = (/71,72,73/)   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_psoilIN    = (/74,75,76/)   !class
  INTEGER,PARAMETER :: o_ro1     = 79   !class
  INTEGER,PARAMETER :: o_ro2     = 80   !class
  INTEGER,PARAMETER :: o_rod     = 82   !class
  INTEGER,PARAMETER :: o_ros     = 83   !class
  INTEGER,PARAMETER :: o_soildenitr = 84   !class
  INTEGER,PARAMETER :: o_cropNupt = 85   !class
  INTEGER,PARAMETER :: o_degrfN   = 86   !class
  INTEGER,PARAMETER :: o_soilNatm = 87   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_psoilSP    = (/91,92,93/)   !class
  INTEGER,PARAMETER :: o_soilPatm = 116   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_pfP = (/94,95,96/)    !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_phP = (/132,133,134/) !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_ppP = (/88,89,90/)    !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_pfC = (/48,50,108/)   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_phC = (/47,49,109/)   !class
  INTEGER,PARAMETER :: o_cprecIN = 97   !<Output variable index
  INTEGER,PARAMETER :: o_cprecSP = 98   !<Output variable index
  INTEGER,PARAMETER :: o_reswe   = 99   !<Output variable index
  INTEGER,PARAMETER :: o_cloc    = 101  !<Output variable index
  INTEGER,PARAMETER :: o_ctmp    = 121  !class
  INTEGER,PARAMETER :: o_coum    = 123  !<Output variable index
  INTEGER,PARAMETER :: o_applirr = 126  !class
  INTEGER,PARAMETER :: o_coub    = 130  !<Output variable index
  INTEGER,PARAMETER :: o_soildef = 131  !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_psoilON = (/136,137,138/)   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_sltmp = (/139,140,141/)   !class
  INTEGER,PARAMETER :: o_rainfall= 142  !class
  INTEGER,PARAMETER :: o_snowfall= 143  !class
  INTEGER,PARAMETER :: o_ro3     = 147  !<Output variable index
  INTEGER,PARAMETER :: o_snowcover = 189   !class
  INTEGER,PARAMETER :: o_snowmax = 191   !class
  INTEGER,PARAMETER :: o_landevap = 195  !class
  INTEGER,PARAMETER :: o_sml1    = 197  !class
  INTEGER,PARAMETER :: o_sml2    = 198  !class
  INTEGER,PARAMETER :: o_sml3    = 199  !class
  INTEGER,PARAMETER :: o_icloss  = 206  !class och subbasin
  INTEGER,PARAMETER :: o_smffc   = 214  !class
  INTEGER,PARAMETER :: o_smfdep  = 215  !class
  INTEGER,PARAMETER :: o_smrzfdep = 216  !class
  INTEGER,PARAMETER :: o_smfpw   = 217  !class
  INTEGER,PARAMETER :: o_smrzfpw = 218  !class
  INTEGER,PARAMETER :: o_xobsm   = 220  !<Output variable index (220-229)
  INTEGER,PARAMETER :: o_xobss   = 230  !<Output variable index (230-239)
  INTEGER,PARAMETER :: o_specificq    = 247  !<Output variable index
  INTEGER,PARAMETER :: o_ros1 = 275 !class
  INTEGER,PARAMETER :: o_ros2 = 276 !class
  INTEGER,PARAMETER :: o_evapsnow = 278 !class
  INTEGER,PARAMETER :: o_cleanwcom = 280  !outlet lake
  INTEGER,PARAMETER :: o_cleanwstr = 281  !outlet lake
  INTEGER,PARAMETER :: o_psim    = 284  !class
  INTEGER,PARAMETER :: o_soilden3  = 303  !class
  INTEGER,PARAMETER :: o_soildenrz = 304  !class
  INTEGER,PARAMETER :: o_sml9    = 330  !class
  INTEGER,PARAMETER :: o_snowmelt= 331  !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_ppT1    = (/333,334,335/)   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_psoilT1 = (/336,337,338/)   !class
  INTEGER,PARAMETER,DIMENSION(3) :: o_csoillayerIN = (/327,328,329/)   !class
  INTEGER,PARAMETER :: o_T1sf    = 341  !class
  INTEGER,PARAMETER :: o_cro1T1  = 343  !class
  INTEGER,PARAMETER :: o_cro2T1  = 344  !class
  INTEGER,PARAMETER :: o_cro3T1  = 345  !class
  INTEGER,PARAMETER :: o_crodT1  = 346  !class
  INTEGER,PARAMETER :: o_crosT1  = 347  !class
  INTEGER,PARAMETER :: o_crunSS  = 348  !class
  INTEGER,PARAMETER :: o_ccSS    = 349  !subbasin
  INTEGER,PARAMETER :: o_reSS    = 350  !subbasin
  INTEGER,PARAMETER :: o_ccAE    = 351  !subbasin
  INTEGER,PARAMETER :: o_ccTS    = 352  !subbasin
  INTEGER,PARAMETER :: o_aowl    = 353  !altimetry  olake water level, m
  INTEGER,PARAMETER :: o_rswa    = 354  !Total surface water area from EO data
  INTEGER,PARAMETER :: o_cswa    = 355  !Total surface water area from model
  INTEGER,PARAMETER :: o_olfv    = 356  !Water volume in olake floodplain water
  INTEGER,PARAMETER :: o_mrfv    = 357  !Water volume in main river floodplain water
  
!>\}

!> \name Model parameter identification variable indices
!> \{
  INTEGER,PARAMETER :: n_lp    = 1    !<Parameter variable index
  INTEGER,PARAMETER :: n_cevpa = 2    !<Parameter variable index
  INTEGER,PARAMETER :: n_cevpp = 3    !<Parameter variable index
  INTEGER,PARAMETER :: n_rivv  = 10   !<Parameter variable index
  INTEGER,PARAMETER :: n_damp  = 11   !<Parameter variable index
  INTEGER,PARAMETER :: n_tcalt = 12   !<Parameter variable index
  INTEGER,PARAMETER :: n_tcea  = 21   !<Parameter variable index
  INTEGER,PARAMETER :: n_pcem  = 22   !<Parameter variable index
  INTEGER,PARAMETER :: n_tcorr = 23   !<Parameter variable index
  INTEGER,PARAMETER :: n_pcea  = 24   !<Parameter variable index
  INTEGER,PARAMETER :: n_pcet  = 29   !<Parameter variable index
  INTEGER,PARAMETER :: n_rrcs3 = 87   !<Parameter variable index
  INTEGER,PARAMETER :: n_cevpc = 117  !<Parameter variable index
  INTEGER,PARAMETER :: n_rrcsc = 140  !<Parameter variable index
  INTEGER,PARAMETER :: n_pcur  = 252  !<Parameter variable index
  INTEGER,PARAMETER :: n_pcus  = 253  !<Parameter variable index
!>\}

!Indices for model parameters
!> \name General model parameter indices
!> \{
  INTEGER,PARAMETER :: m_lp        = 1  !<General model parameter index
  INTEGER,PARAMETER :: m_cevpam    = 2  !<General model parameter index 
  INTEGER,PARAMETER :: m_cevpph    = 3  !<General model parameter index 
  INTEGER,PARAMETER :: m_deadl     = 4  !<General model parameter index 
  INTEGER,PARAMETER :: m_fastN0    = 5  !<General model parameter index  
  INTEGER,PARAMETER :: m_fastP0    = 6  !<General model parameter index  
  INTEGER,PARAMETER :: m_iniT1     = 7  !<General model parameter index
  INTEGER,PARAMETER :: m_iniT2     = 8  !<General model parameter index
  INTEGER,PARAMETER :: m_T1evap    = 9  !<General model parameter index 
  INTEGER,PARAMETER :: m_rivvel    = 10 !<General model parameter index 
  INTEGER,PARAMETER :: m_damp      = 11 !<General model parameter index 
  INTEGER,PARAMETER :: m_tcalt     = 12 !<General model parameter index 
  INTEGER,PARAMETER :: m_grat1     = 13 !<General model parameter index 
  INTEGER,PARAMETER :: m_grat2     = 14 !<General model parameter index 
  INTEGER,PARAMETER :: m_deadm     = 15 !<General model parameter index 
  INTEGER,PARAMETER :: m_denitwr   = 16 !<General model parameter index 
  INTEGER,PARAMETER :: m_limprod   = 17 !<General model parameter index 
  INTEGER,PARAMETER :: m_denitwl   = 18 !<General model parameter index 
  INTEGER,PARAMETER :: m_laketemp  = 19 !<General model parameter index 
  INTEGER,PARAMETER :: m_littdays  = 20 !<General model parameter index 
  INTEGER,PARAMETER :: m_denitwrl  = 21 !<General model parameter index 
  INTEGER,PARAMETER :: m_sndens0   = 22 !<General model parameter index 
  INTEGER,PARAMETER :: m_dsndens   = 23 !<General model parameter index 
  INTEGER,PARAMETER :: m_pprelmax  = 24 !<General model parameter index 
  INTEGER,PARAMETER :: m_pprelexp  = 25 !<General model parameter index 
  INTEGER,PARAMETER :: m_epotdist  = 26 !<General model parameter index 
  INTEGER,PARAMETER :: m_fastlake  = 27 !<General model parameter index 
  INTEGER,PARAMETER :: m_wprodn    = 28 !<General model parameter index 
  INTEGER,PARAMETER :: m_wprodc    = 29 !<General model parameter index 
  INTEGER,PARAMETER :: m_wprodp    = 30 !<General model parameter index 
  INTEGER,PARAMETER :: m_sedon     = 31 !<General model parameter index 
  INTEGER,PARAMETER :: m_sedoc     = 32 !<General model parameter index 
  INTEGER,PARAMETER :: m_sedpp     = 33 !<General model parameter index 
  INTEGER,PARAMETER :: m_sedexp    = 34 !<General model parameter index 
  INTEGER,PARAMETER :: m_rcgrw     = 35 !<General model parameter index 
  INTEGER,PARAMETER :: m_gldepi    = 36 !<General model parameter index 
  INTEGER,PARAMETER :: m_locsoil   = 37 !<General model parameter index 
  INTEGER,PARAMETER :: m_qmean     = 38 !<General model parameter index 
  INTEGER,PARAMETER :: m_deepmem   = 39 !<General model parameter index 
  INTEGER,PARAMETER :: m_wetsp     = 40 !<General model parameter index 
  INTEGER,PARAMETER :: m_sreroexp  = 41 !<General model parameter index 
  INTEGER,PARAMETER :: m_fertdays  = 42 !<General model parameter index 
  INTEGER,PARAMETER :: m_deeplake  = 43 !<General model parameter index 
  INTEGER,PARAMETER :: m_rrcs3     = 44 !<General model parameter index 
  INTEGER,PARAMETER :: m_pcaddg    = 45 !<General model parameter index 
  INTEGER,PARAMETER :: m_pcelevth  = 46 !<General model parameter index 
  INTEGER,PARAMETER :: m_minc      = 47 !<General model parameter index 
  INTEGER,PARAMETER :: m_crate1    = 48 !<General model parameter index 
  INTEGER,PARAMETER :: m_crate2    = 49 !<General model parameter index  
  INTEGER,PARAMETER :: m_crate3    = 50 !<General model parameter index 
  INTEGER,PARAMETER :: m_ripe      = 52 !<General model parameter index 
  INTEGER,PARAMETER :: m_rips      = 53 !<General model parameter index 
  INTEGER,PARAMETER :: m_crate5    = 54 !<General model parameter index   
  INTEGER,PARAMETER :: m_crate6    = 55 !<General model parameter index  
  INTEGER,PARAMETER :: m_crate9    = 56 !<General model parameter index 
  INTEGER,PARAMETER :: m_crate10   = 57 !<General model parameter index 
  INTEGER,PARAMETER :: m_maxwidth  = 58 !<General model parameter index 
  INTEGER,PARAMETER :: m_pcelevadd = 60 !<General model parameter index 
  INTEGER,PARAMETER :: m_pcelevmax = 61 !<General model parameter index 
  INTEGER,PARAMETER :: m_tcelevadd = 62 !<General model parameter index 
  INTEGER,PARAMETER :: m_regirr    = 63 !<General model parameter index 
  INTEGER,PARAMETER :: m_sswcorr   = 64 !<General model parameter index 
  INTEGER,PARAMETER :: m_immdep    = 65 !<General model parameter index 
  INTEGER,PARAMETER :: m_iwdfrac   = 66 !<General model parameter index 
  INTEGER,PARAMETER :: m_wdpar     = 67 !<General model parameter index 
  INTEGER,PARAMETER :: m_ttpd      = 68 !<General model parameter index 
  INTEGER,PARAMETER :: m_ttpi      = 69 !<General model parameter index 
  INTEGER,PARAMETER :: m_irrcomp   = 70 !<General model parameter index 
  INTEGER,PARAMETER :: m_tcobselev = 71 !<General model parameter index 
  INTEGER,PARAMETER :: m_atmload   = 72 !<General model parameter index
  INTEGER,PARAMETER :: m_denitaq   = 73 !<General model parameter index
  INTEGER,PARAMETER :: m_krelflood = 74 !<General model parameter index
!--Lake and river water parameters  
  INTEGER,PARAMETER :: m_t2trriver = 75 !<General model parameter index
  INTEGER,PARAMETER :: m_upper2deep= 76 !<General model parameter index
  INTEGER,PARAMETER :: m_t2trlake  = 77 !<General model parameter index
!---Lake ice parameters
  INTEGER,PARAMETER :: m_licewme   = 78 !<General model parameter index
  INTEGER,PARAMETER :: m_liceTF    = 79 !<General model parameter index
  INTEGER,PARAMETER :: m_licesndens= 80 !<General model parameter index
  INTEGER,PARAMETER :: m_licekika  = 81 !<General model parameter index
  INTEGER,PARAMETER :: m_licekexp  = 82 !<General model parameter index
  INTEGER,PARAMETER :: m_licetmelt = 83 !<General model parameter index
!---River ice parameter
  INTEGER,PARAMETER :: m_ricewme   = 84 !<General model parameter index
  INTEGER,PARAMETER :: m_riceTF    = 85 !<General model parameter index
  INTEGER,PARAMETER :: m_ricesndens= 86 !<General model parameter index
  INTEGER,PARAMETER :: m_ricekika  = 87 !<General model parameter index
  INTEGER,PARAMETER :: m_ricekexp  = 88 !<General model parameter index
  INTEGER,PARAMETER :: m_ricetmelt = 89 !<General model parameter index
  INTEGER,PARAMETER :: m_licewcorr = 90 !<General model parameter index
!--Snow cover fraction, general
  INTEGER,PARAMETER :: m_fscmax    = 91 !<General model parameter index
  INTEGER,PARAMETER :: m_fscmin    = 92 !<General model parameter index
  INTEGER,PARAMETER :: m_fsclim    = 93 !<General model parameter index
  INTEGER,PARAMETER :: m_fsck1     = 94 !<General model parameter index
  INTEGER,PARAMETER :: m_fsckexp   = 95 !<General model parameter index
!--Radiation estimation, and Potential Evapotranspiration parameters
  INTEGER,PARAMETER :: m_krs       = 96
  INTEGER,PARAMETER :: m_jhtadd    = 97
  INTEGER,PARAMETER :: m_jhtscale  = 98
  INTEGER,PARAMETER :: m_alfapt    = 99
  INTEGER,PARAMETER :: m_mwind     = 100
  INTEGER,PARAMETER :: m_grat3     = 101 !<General model parameter index 
!--Glacier parameters
  INTEGER,PARAMETER :: m_glacvcoef = 102
  INTEGER,PARAMETER :: m_glacvexp  = 103
  INTEGER,PARAMETER :: m_glacdens  = 104
  INTEGER,PARAMETER :: m_glacvcoef1 = 105
  INTEGER,PARAMETER :: m_glacvexp1  = 106
  INTEGER,PARAMETER :: m_glac2arlim = 107
!--Lake and river water temperature parameters, eventually replacing the t2trriver and t2trlake parameters
  INTEGER,PARAMETER :: m_tcfriver  = 108
  INTEGER,PARAMETER :: m_scfriver  = 109
  INTEGER,PARAMETER :: m_ccfriver  = 110
  INTEGER,PARAMETER :: m_lcfriver  = 111
  INTEGER,PARAMETER :: m_tcflake   = 112
  INTEGER,PARAMETER :: m_scflake   = 113
  INTEGER,PARAMETER :: m_ccflake   = 114
  INTEGER,PARAMETER :: m_lcflake   = 115
  INTEGER,PARAMETER :: m_stbcorr1  = 116
  INTEGER,PARAMETER :: m_stbcorr2  = 117
  INTEGER,PARAMETER :: m_stbcorr3  = 118
  INTEGER,PARAMETER :: m_zwind     = 119
  INTEGER,PARAMETER :: m_zwish     = 120
  INTEGER,PARAMETER :: m_zpdh      = 121
  INTEGER,PARAMETER :: m_roughness = 122
  INTEGER,PARAMETER :: m_pcelevstd = 123
  INTEGER,PARAMETER :: m_pcurain   = 124
  INTEGER,PARAMETER :: m_pcusnow   = 125
  INTEGER,PARAMETER :: m_cmrefr    = 126 
  INTEGER,PARAMETER :: m_fsceff    = 127
  INTEGER,PARAMETER :: m_glacalb   = 128           !glacier ice albedo (0.35)
  INTEGER,PARAMETER :: m_glacttmp  = 129
  INTEGER,PARAMETER :: m_glaccmlt  = 130
  INTEGER,PARAMETER :: m_glaccmrad = 131
  INTEGER,PARAMETER :: m_glaccmrefr= 132
  INTEGER,PARAMETER :: m_fepotglac = 133

  INTEGER,PARAMETER :: m_kthrflood = 134 !<General model parameter index
  INTEGER,PARAMETER :: m_klowflood = 135 !<General model parameter index
  INTEGER,PARAMETER :: m_limsedon  = 136 !<General model parameter index
  INTEGER,PARAMETER :: m_limsedpp  = 137 !<General model parameter index
!--Optimization parameters
  INTEGER,PARAMETER :: m_opt1     = 138
  INTEGER,PARAMETER :: m_opt2     = 139
  INTEGER,PARAMETER :: m_opt3     = 140
  INTEGER,PARAMETER :: m_opt4     = 141
  INTEGER,PARAMETER :: m_opt5     = 142
  INTEGER,PARAMETER :: m_opt6     = 143
  INTEGER,PARAMETER :: m_opt7     = 144
  INTEGER,PARAMETER :: m_opt8     = 145
  INTEGER,PARAMETER :: m_optonoff = 146

  INTEGER,PARAMETER :: m_sdnsmax  = 147
  INTEGER,PARAMETER :: m_sdnsrate = 148
  INTEGER,PARAMETER :: m_sdnsradd = 149
  INTEGER,PARAMETER :: m_gt2mix   = 150
  INTEGER,PARAMETER :: m_limt2exch = 151
  INTEGER,PARAMETER :: m_hsatINsoil = 152
  INTEGER,PARAMETER :: m_hsatINwater = 153
  INTEGER,PARAMETER :: m_hsatTP   = 154
!--Tracer parameters
  INTEGER,PARAMETER :: m_init1sw   = 155
  INTEGER,PARAMETER :: m_expdec    = 156
  INTEGER,PARAMETER :: m_t1freuc   = 157
  INTEGER,PARAMETER :: m_t1rel     = 158
  INTEGER,PARAMETER :: m_t1sed     = 159
  INTEGER,PARAMETER :: m_t1sedexp  = 160
  INTEGER,PARAMETER :: m_gldepo    = 161
  INTEGER,PARAMETER :: m_gicatch   = 162
  INTEGER,PARAMETER :: m_glacannmb = 163
  INTEGER,PARAMETER :: m_erodslope = 164
  INTEGER,PARAMETER :: m_erodexp   = 165
  INTEGER,PARAMETER :: m_erodindex = 166
  INTEGER,PARAMETER :: m_sedss = 167
  INTEGER,PARAMETER :: m_sedae = 168
  INTEGER,PARAMETER :: m_limsedss = 169
!>\}

  !----basin parameters, none

!> \name Soil type dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_wcfc     = 1  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcwp     = 2  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcep     = 3  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcfc1    = 4  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcwp1    = 5  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcep1    = 6  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcfc2    = 7  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcwp2    = 8  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcep2    = 9  !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcfc3    = 10 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcwp3    = 11 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_wcep3    = 12 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_srrate   = 13 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_trrcs    = 14 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_freuc    = 16 !<Soil type dependent model parameter index   
  INTEGER,PARAMETER :: m_freuexp  = 17 !<Soil type dependent model parameter index  
  INTEGER,PARAMETER :: m_freurate = 18 !<Soil type dependent model parameter index   
  INTEGER,PARAMETER :: m_rrcs1    = 19 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_rrcs2    = 20 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_macrate  = 21 !<Soil type dependent model parameter index  
  INTEGER,PARAMETER :: m_mactrinf = 22 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_mactrsm  = 23 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_soilcoh  = 24 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_soilerod = 25 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_sfrost   = 26 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_macfilt  = 27 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_perc1    = 28 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_perc2    = 29 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_rcgrwst  = 30 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_t1leakst = 31 !<Soil type dependent model parameter index 
  INTEGER,PARAMETER :: m_bfroznsoil = 32
  INTEGER,PARAMETER :: m_erodsoil = 33
!>\}

!> \name Land use dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_ttmp       = 1  !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_cmlt       = 2  !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_cevp       = 3  !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_cfrost     = 4  !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_srrcs      = 5  !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_humusN0    = 6  !<Land use dependent model parameter index   
  INTEGER,PARAMETER :: m_partP0     = 7  !<Land use dependent model parameter index   
  INTEGER,PARAMETER :: m_depthrel   = 8  !<Land use dependent model parameter index   
  INTEGER,PARAMETER :: m_surfmem    = 9  !<Land use dependent model parameter index   
  INTEGER,PARAMETER :: m_hNhalf     = 10 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_pPhalf     = 11 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_drypp      = 12 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_humusP0    = 13 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_hPhalf     = 14 !<Land use dependent model parameter index   
  INTEGER,PARAMETER :: m_denitrlu   = 15 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_ponatm     = 16 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_partP1     = 17 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_partP2     = 18 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_partP3     = 19 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_pcluse     = 20 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_minerfN    = 21 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_minerfP    = 22 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_degradhN   = 23 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_ripz       = 24 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_degradhP   = 25 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_filtPbuf   = 26 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_filtPinner = 27 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_filtPother = 28 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_humusC1    = 29 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_humusC2    = 30 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_humusC3    = 31 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_fastC1     = 32 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_fastC2     = 33 !<Land use dependent model parameter index  
  INTEGER,PARAMETER :: m_fastC3     = 34 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_dissolfN   = 35 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_dissolfP   = 36 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_dissolhN   = 37 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_dissolhP   = 38 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_ocsoim     = 39 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_ocsmslp    = 40 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_onconc0    = 41 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_ppconc0    = 42 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_occonc0    = 43 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_onpercred  = 44 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_pppercred  = 45 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_snalbmin   = 46 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_snalbmax   = 47 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_snalbkexp  = 48 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_cmrad      = 49 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_fscdistmax = 50 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_fscdist0   = 51 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_fscdist1   = 52 !<Land use dependent model parameter index
  INTEGER,PARAMETER,DIMENSION(5) :: m_kc = (/53,62,63,64,65/) !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_alb        = 54 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_fepotsnow  = 55 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_denitr3    = 56 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_t1leaklu   = 57 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_soilstretch = 58 !<Land use dependent model parameter index 
  INTEGER,PARAMETER :: m_ttrig      = 59 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_tredA      = 60 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_tredB      = 61 !<Land use dependent model parameter index
  INTEGER,PARAMETER :: m_erodluse   = 66 !<Land use dependent model parameter index
!>\}

!> \name Ilake region dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_ilrrat1     = 16
  INTEGER,PARAMETER :: m_ilrrat2     = 17
  INTEGER,PARAMETER :: m_ilrldep     = 18
  INTEGER,PARAMETER :: m_ilricatch   = 19
!>\}

!> \name Olake region dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_olrrat1  = 20
  INTEGER,PARAMETER :: m_olrrat2  = 21
  INTEGER,PARAMETER :: m_olrldep  = 22
!>\}

!> \name Surface water (lake) region dependent model parameter indices (some only for olakes)
!> \{
  INTEGER,PARAMETER :: m_velpar1  = 23
  INTEGER,PARAMETER :: m_velpar2  = 24
  INTEGER,PARAMETER :: m_velpar3  = 25
  INTEGER,PARAMETER :: m_widpar1  = 26
  INTEGER,PARAMETER :: m_widpar2  = 27
  INTEGER,PARAMETER :: m_widpar3  = 28
  INTEGER,PARAMETER :: m_tpmean   = 29
  INTEGER,PARAMETER :: m_tnmean   = 30
  INTEGER,PARAMETER :: m_tocmean  = 31
!>\}

!> \name Parameter region dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_cevpcorr   = 1  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_rrcscorr   = 2  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_ratcorr    = 3  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_pirrs      = 4  !<Parameter region parameter index 
  INTEGER,PARAMETER :: m_preccorr   = 5  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_cirrsink   = 6  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_tcadd      = 7  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_pirrg      = 8  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_aqretcorr  = 9  !<Parameter region parameter index
  INTEGER,PARAMETER :: m_aqdelcorr  = 10 !<Parameter region parameter index
  INTEGER,PARAMETER :: m_aqpercorr  = 11 !<Parameter region parameter index
!>\}

!> \name Water quality parameter region dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_incorr     = 12  !<Water quality parameter region parameter index  
  INTEGER,PARAMETER :: m_oncorr     = 13  !<Water quality parameter region parameter index  
  INTEGER,PARAMETER :: m_phoscorr   = 14  !<Water quality parameter region parameter index
  INTEGER,PARAMETER :: m_denit3reg  = 15  !<Water quality parameter region parameter index
  !>\}

!> \name LakeData dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_lddenitwl  = 1  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldqmean    = 2  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldtpmean   = 3  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldtnmean   = 4  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldtocmean  = 5  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldwprodn   = 6  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldwprodp   = 7  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldwprodc   = 8  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldsedon    = 9  !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldsedoc    = 10 !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldsedpp    = 11 !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldprodpp   = 12 !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldprodsp   = 13 !<LakeData parameter index
  INTEGER,PARAMETER :: m_lddeeplake = 14 !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldfastlake = 15 !<LakeData parameter index
  INTEGER,PARAMETER :: m_ldt2mix    = 16 !<LakeData parameter index
!>\}

!> \name Monthly dependent model parameter indices
!> \{
  INTEGER,PARAMETER :: m_mlapse  = 1
!>\}

!> \name Glacier related constants
!> \{
  !Soil models, starting on 3 (0=default,1=olake,2=ilake)
  INTEGER, PARAMETER :: glacier_model    = 3   !<Glacier soil model identifier
!>\}

!> \name Nitrogen, phosphorus and organic carbon process constants
!> \{
    REAL,PARAMETER :: bulkdensity     = 1300.   !<soil density (kg/m3)
    REAL,PARAMETER :: NPratio         = 1.0/7.2 !<NP ratio for production/mineralisation in water 
    REAL,PARAMETER :: dryNratio       = 20.     !<dry weight to N ratio of algae
    REAL,PARAMETER :: maxdenitriwater = 0.5     !<maximum part of IN pool that can be denitrified per timestep
    REAL,PARAMETER :: maxprodwater    = 0.5     !<maximum part of IN/SRP pool that can be used for production per timestep
    REAL,PARAMETER :: maxdegradwater  = 0.5     !<maximum part of ON/PP/OC pool that can be degraded per timestep
    REAL,PARAMETER :: Satact          = 0.6     !<parameter for soil moisture factor for degradation/transformation in soil
    REAL,PARAMETER :: Thetalow        = 8.      !<parameter for soil moisture factor for degradation/transformation in soil
    REAL,PARAMETER :: Thetaupp        = 12.     !<parameter for soil moisture factor for degradation/transformation in soil
    REAL,PARAMETER :: Thetapow        = 1.      !<parameter for soil moisture factor for degradation/transformation in soil
    REAL,PARAMETER :: smfdenitlim     = 0.7     !<parameter for soil moisture factor for denitrification in soil
    REAL,PARAMETER :: smfdenitpow     = 2.5     !<parameter for soil moisture factor for denitrification in soil
    REAL,PARAMETER :: NCratio         = 5.7     !<NC ratio for production/mineralisation in water (värde från BIOLA)
!>\}

!>\name Static variables
!>
!>The static variables are calculated before model time simulation starts and holds their value during simulation. 
!>This includes some status variables, recalculation of model parameters into a form suitable for the model 
!>processes and other useful variables.
!>\{
  LOGICAL :: calcSWRAD                     !<flagging need of shortwave radiation (depending on snowmelt or PET function)
  LOGICAL :: calcVAPOUR                    !<flagging need of vapour pressures and net radiation (depending on PET function)
  LOGICAL :: calcWIND                      !<flagging need of wind speed data (for Penman-Monteith PET function)
  LOGICAL :: T1leakage                     !<T1 is having a parameter leakage source
  REAL, ALLOCATABLE :: epotdist(:,:)       !<relative distribution of potential evaporation between upper two soil layers (-) (soillayer,class)
  REAL, ALLOCATABLE :: transtime(:,:)      !<translation time in river (total)
  REAL, ALLOCATABLE :: ttpart(:,:)         !<translation time in river (part of time step)
  INTEGER, ALLOCATABLE :: ttstep(:,:)      !<translation time in river (whole time step)
  REAL, ALLOCATABLE :: riverrc(:,:)        !<recession parameter for river damping
  REAL, ALLOCATABLE :: riverlength(:,:)    !<river length [m] (rivertype,subbasin)
  REAL, ALLOCATABLE :: deadriver(:,:)      !<volume of dead water in rivers (m3) (rivertype,subbasin)
  REAL, ALLOCATABLE :: deadwidth(:,:)      !<river width based om deadvolume (m) (rivertype,subbasin)
  REAL, ALLOCATABLE :: ratingk(:,:)        !<k-parameter for rating curve (laketype,subbasin)
  REAL, ALLOCATABLE :: wpmm(:,:)           !<water holding parameter (maxsoillayers,nclass)
  REAL, ALLOCATABLE :: fcmm(:,:)           !<water holding parameter (maxsoillayers,nclass)
  REAL, ALLOCATABLE :: epmm(:,:)           !<water holding parameter (maxsoillayers,nclass)
  REAL, ALLOCATABLE :: pwmm(:,:)           !<water holding capacity, pore volume (maxsoillayers,nclass)
  REAL, ALLOCATABLE :: soilrc(:,:,:)       !<soil runoff recession coefficients (maxsoillayers,nclass,nsub)
  REAL, ALLOCATABLE :: soilmem(:,:)        !<soil temperature memory (maxsoillayers,nclass)
  REAL, ALLOCATABLE :: basinrrcscorr(:)    !<subbasin rrcs correction (subbasin)
  REAL, ALLOCATABLE :: basincevpcorr(:)    !<subbasin evaporation correction (subbasin)
  REAL, ALLOCATABLE :: basincevpam(:)      !<subbasin evaporation correction (subbasin)
  REAL, ALLOCATABLE :: basincevpph(:)      !<subbasin evaporation correction (subbasin)
  REAL, ALLOCATABLE :: basinlp(:)          !<subbasin evaporation correction (subbasin)
  REAL, ALLOCATABLE :: basintcalt(:)       !<subbasin temperature correction (subbasin)
  REAL, ALLOCATABLE :: basintempadd(:)     !<subbasin temperature correction (subbasin)
  REAL, ALLOCATABLE :: basinpreccorr(:,:)  !<subbasin precipitation correction (subbasin,nclass)
  REAL, ALLOCATABLE :: basinpcurain(:)     !<subbasin rain correction (subbasin)
  REAL, ALLOCATABLE :: basinpcusnow(:)     !<subbasin snow correction (subbasin)
  REAL, ALLOCATABLE :: windtrans(:)        !<Wind transformation factor (nclass)
  REAL, ALLOCATABLE :: slowlakeini(:,:)    !<target volume of slowlake (mm) (laketype,subbasin)
  REAL :: avertemp(4)                      !<number of timesteps over which meantemp is calculated
!> \brief type for outflow for lakes with two outflows
  TYPE OUTFLOWTYPE
    INTEGER :: otype(2) = 0 !<type of outflow for two outlet lakes (depending on LakeData-parameters)
    INTEGER :: change = 0   !<flag for method to change outlet flows after updating total flow (depent on otype)
  END TYPE OUTFLOWTYPE
  TYPE(OUTFLOWTYPE), ALLOCATABLE :: lakeoutlet(:) !<information of lake outlet type (nbranch)
!>\}

!>\name Variables for river flows
!>
!>Some variables to help with river flow calculations. The variables hold values for each subbasin 
!>and often for both local and main flow. The variables for maximum daily flow during the last 
!>year are relating to the flow values saved in the state variable riverQ365.
!>\{
  !Model variables for river
  REAL, ALLOCATABLE :: Qmax(:,:),Q2max(:,:)      !<max Q and 2nd highest Q (=bankfull flow), from riverQ365 (rivertype,subbasin)
  INTEGER, ALLOCATABLE :: iQmax(:,:),iQ2max(:,:) !<index of max Q and 2nd highest Q in riverQ365 (rivertype,subbasin)
  !Variable for output
  REAL, ALLOCATABLE :: accdiff(:)           !<accumulated flow difference between simulated and observed flow, since last observed flow
!>\}

!>\name Load variables for source apportionment
!>
!>HYPE has variables to hold calculated nutrient loads. These variables contain the current time 
!>step load or transport of nitrogen and phosphorus for each source and subbasin. Two of the 
!>variables hold transport in different points in the river system (Lstream, Lpathway).
!>\{
  REAL,DIMENSION(:,:,:), ALLOCATABLE    :: Latmdep  !< wet & dry atmospheric loads, (kg/timestep)
  REAL,DIMENSION(:,:,:), ALLOCATABLE    :: Lcultiv  !< Cultivation loads, fertilizer and plant decay (kg/timestep)
  REAL,DIMENSION(:,:),   ALLOCATABLE    :: Lirrsoil !< Class source load; irrigation on soil (kg/timestep)
  REAL,DIMENSION(:,:),   ALLOCATABLE    :: Lrurala  !< Class source load; rural a (kg/timestep)
  REAL,DIMENSION(:),     ALLOCATABLE    :: Lruralb  !< River load; rural b (kg/timestep)
  REAL,DIMENSION(:,:),   ALLOCATABLE    :: Lpoints  !< Point Source Loads, urban/indust (kg/timestep)
  REAL,DIMENSION(:,:),   ALLOCATABLE    :: Lstream  !< Loads by runoff to stream (soil, surface, tile), class dependent (kg/timestep)
  REAL,DIMENSION(:,:),   ALLOCATABLE    :: Lpathway !< Loads from local stream to main lake (kg/timestep)
  REAL,DIMENSION(:),     ALLOCATABLE    :: Lbranch  !< Load in branch (kg/timestep)
  REAL,DIMENSION(:,:),   ALLOCATABLE    :: Lgrwsoil !< Soil loads from regional groundwater flow (kg/timestep)
  REAL,DIMENSION(:),     ALLOCATABLE    :: Lgrwmr   !< Main river loads from regional groundwater (kg/timestep)
  REAL,DIMENSION(:),     ALLOCATABLE    :: Lgrwol   !< Outlet lake loads from regional groundwater (kg/timestep)
  REAL,DIMENSION(:,:,:), ALLOCATABLE    :: Lgrwclass !< Loads from soil by regional groundwater flow (kg/timestep) (not printed)
!>\}

CONTAINS


  !>Initiate variables used for calculation of bankful flow
  !-----------------------------------------------------------
  SUBROUTINE set_Qvariables_for_bankfulflow(nsub,riverQ365)

    IMPLICIT NONE

    !Argument declarations
    INTEGER,INTENT(IN) :: nsub  !<number of subbasins
    REAL,INTENT(IN) :: riverQ365(366,2,nsub)  !<daily river flow last year

    !Local variables
    INTEGER i,j,iQ
    REAL riverQ365Temp(366)

      !Allocate 
      IF(.NOT.ALLOCATED(Qmax)) ALLOCATE(Qmax(2,nsub))   !!Qmax och Q2max behövs väl inte sparas?
      IF(.NOT.ALLOCATED(Q2max)) ALLOCATE(Q2max(2,nsub))
      IF(.NOT.ALLOCATED(iQmax)) ALLOCATE(iQmax(2,nsub))
      IF(.NOT.ALLOCATED(iQ2max)) ALLOCATE(iQ2max(2,nsub))
     
      !Set parameter Qmax, iQmax, Q2max and iQ2max from riverstate%Q365
      DO i = 1,nsub
        DO j = 1,2
          riverQ365Temp = riverQ365(:,j,i)
          iQmax(j,i) = MAXLOC(riverQ365Temp,1)
          Qmax(j,i)  = riverQ365Temp(iQmax(j,i))
          riverQ365Temp(iQmax(j,i)) = 0
          iQ = MAXLOC(riverQ365Temp,1)
          riverQ365Temp(iQ) = 0              !!!????varför ta bort två, få 3:e högsta???
          iQ2max(j,i) = MAXLOC(riverQ365Temp,1)
          Q2max(j,i)  = riverQ365Temp(iQ2max(j,i))
        ENDDO
      ENDDO

  END SUBROUTINE set_Qvariables_for_bankfulflow

END MODULE HYPEVARIABLES
