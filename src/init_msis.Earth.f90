subroutine get_msis_temperature(lon, lat, alt, t, h)

 use ModIndicesInterfaces
  use ModTime
  use ModInputs
  use ModPlanet
  use ModGITM

  use EUA_ModMsis90, only: meter6, gtd6

  implicit none

  real, intent(in) :: lon, lat, alt
  real, intent(out) :: t, h

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens

  real :: LonDeg, LatDeg, AltKm, LST
  real,dimension(7)    :: AP  
  real :: nO, nO2, nN2, m, r, g

   integer :: iError
  !-------------------------------------------------------
  
  ap=10.0

  call meter6(.true.)

  LonDeg = lon*180.0/pi
  LatDeg = lat*180.0/pi
  AltKm  = alt/1000.0
  LST = mod(utime/3600.0+LonDeg/15.0,24.0)
  iError = 0
  
  call get_f107(CurrentTime, f107, iError)

  if (iError /= 0) then
     write(*,*) "Error in getting F107 value.  Is this set?"
     write(*,*) "Code : ",iError
     call stop_gitm("Stopping in euv_ionization_heat")
  endif

  call get_f107a(CurrentTime, f107a, iError)
  if (iError /= 0) then
     write(*,*) "Error in getting F107a value.  Is this set?"
     write(*,*) "Code : ",iError
     call stop_gitm("Stopping in euv_ionization_heat")
  endif

  CALL GTD6(iJulianDay,utime,AltKm,LatDeg,LonDeg,LST, &
       F107A,F107,AP,48,msis_dens,msis_temp)

  t = msis_temp(2)
  nO  = msis_dens(2)
  nN2 = msis_dens(3)
  nO2 = msis_dens(4)

  m = (nO * mass(iO_3P_) + nO2 * mass(iO2_) + nN2 * mass(iN2_)) / (nO + nO2 + nN2)

  r = RBody + alt
!  g = Gravitational_Constant * (RBody/r) ** 2
  g = Gravitational_Constant

  h = Boltzmanns_Constant * t / (m*g)

end subroutine get_msis_temperature

!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine init_msis

  use ModGITM
  use ModInputs
  use ModConstants
  use ModPlanet
  use ModTime

  use EUA_ModMsis90, ONLY: meter6, gtd6, tselec

  implicit none

  ! msis variables

  real, dimension(1:2) :: msis_temp
  real, dimension(1:8) :: msis_dens

  real, dimension(25) :: sw

  integer :: iBlock, iAlt, iLat, iLon, iSpecies, iyd
  real :: geo_lat, geo_lon, geo_alt, geo_lst,m,k, ut
  real, dimension(7)  :: ap = 10.0

  real(4) :: hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst
  real(4) :: hwm_f107a, hwm_f107, hwm_ap(2), qw(2)

  call report("init_msis",1)

  !--------------------------------------------------------------------------
  !
  !  From the msis90 library:
  !
  !     OUTPUT:
  !        D(1) - HE NUMBER DENSITY(CM-3)
  !        D(2) - O NUMBER DENSITY(CM-3)
  !        D(3) - N2 NUMBER DENSITY(CM-3)
  !        D(4) - O2 NUMBER DENSITY(CM-3)
  !        D(5) - AR NUMBER DENSITY(CM-3)                       
  !        D(6) - TOTAL MASS DENSITY(GM/CM3)
  !        D(7) - H NUMBER DENSITY(CM-3)
  !        D(8) - N NUMBER DENSITY(CM-3)
  !        T(1) - EXOSPHERIC TEMPERATURE
  !        T(2) - TEMPERATURE AT ALT
  !
  !      TO GET OUTPUT IN M-3 and KG/M3:   CALL METER6(.TRUE.) 
  !
  !      O, H, and N set to zero below 72.5 km
  !      Exospheric temperature set to average for altitudes below 120 km.
  !
  !--------------------------------------------------------------------------

  ! We want units of /m3 and not /cm3

  call meter6(.true.)

  if (UseMsisTides) then
     sw = 1
  ELSE IF (UseMSISOnly) THEN
     ! Diurnal, semidiurnal, and terdiurnal variations are excluded, EYigit:16June09
     CALL report("...Using MSIS without tidal variations...",0)
     sw = 1
     sw(7) = 0
     sw(8) = 0
     sw(14) = 0
  ELSE
     sw = 0
     sw(1) = 1
     sw(9) = 1
  endif

  sw(9) = 0
  sw(2) = 0

  call tselec(sw)

  if (DoRestart) return

  !           The following is for test and special purposes:
  !            TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW)
  !               WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1. 
  !               FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
  !               FOR THE FOLLOWING VARIATIONS
  !               1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
  !               3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
  !               5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
  !               7 - DIURNAL               8 - SEMIDIURNAL
  !               9 - DAILY AP             10 - ALL UT/LONG EFFECTS
  !              11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
  !              13 - MIXED AP/UT/LONG     14 - TERDIURNAL
  !              15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
  !              16 - ALL TINF VAR         17 - ALL TLB VAR
  !              18 - ALL TN1 VAR           19 - ALL S VAR
  !              20 - ALL TN2 VAR           21 - ALL NLB VAR
  !              22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR

  ! Initialize data

  iyd = iTimeArray(1)*1000 + iJulianDay

  do iBlock = 1, nBlocks
     do iAlt = -1, nAlts+2
        do iLon=-1,nLons+2
           do iLat=-1,nLats+2

              geo_lat = Latitude(iLat,iBlock)*180.0/pi
              geo_lon = Longitude(iLon,iBlock)*180.0/pi

              geo_alt = Altitude_GB(iLon, iLat, iAlt, iBlock)/1000.0
              geo_lst = mod(utime/3600.0+geo_lon/15.0,24.0)

              !
              ! Call MSIS (results will be im mks units)
              !

              CALL GTD6(iJulianDay,utime,geo_alt,geo_lat,geo_lon,geo_lst, &
                   F107A,F107,AP,48,msis_dens,msis_temp)

              ! Initialize densities to zero in case msis does not set it
              NDensityS(iLon,iLat,iAlt,:,iBlock) = 0.0

              NDensityS(iLon,iLat,iAlt,iH_,iBlock)          = msis_dens(1)
              NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)          = msis_dens(2)
              NDensityS(iLon,iLat,iAlt,iN2_,iBlock)         = msis_dens(3)
              NDensityS(iLon,iLat,iAlt,iO2_,iBlock)         = msis_dens(4)
              NDensityS(iLon,iLat,iAlt,iAr_,iBlock)         = msis_dens(5)
              NDensityS(iLon,iLat,iAlt,iHe_,iBlock)         = msis_dens(7)
              NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)       = msis_dens(8)
              NDensityS(iLon,iLat,iAlt,iN_2P_,iBlock)       = &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)/10000.0
              NDensityS(iLon,iLat,iAlt,iN_2D_,iBlock)       = &
                   NDensityS(iLon,iLat,iAlt,iN_4S_,iBlock)/100.0
              NDensityS(iLon,iLat,iAlt,iO_1D_,iBlock)       = &
                   NDensityS(iLon,iLat,iAlt,iO_3P_,iBlock)/1000000.0 *0.0 + 1

              MeanMajorMass(iLon,iLat,iAlt)=0

              do iSpecies = 1, nSpecies
                 MeanMajorMass(iLon,iLat,iAlt) = MeanMajorMass(iLon,iLat,iAlt) +   &
                      Mass(iSpecies) * NDensityS(iLon,iLat,iAlt,iSpecies,iBlock)/   &
                      sum(NDensityS(iLon,iLat,iAlt,1:3,iBlock))
              enddo
  
              TempUnit(iLon,iLat,iAlt) = &
                   MeanMajorMass(iLon,iLat,iAlt)/ Boltzmanns_Constant

              Temperature(iLon,iLat,iAlt,iBlock) = &
                   msis_temp(2)/TempUnit(iLon,iLat,iAlt)

              Rho(iLon,iLat,iAlt,iBlock) = msis_dens(6)

              ! The initial profile of [NO] is refered to:
              !  [Charles A. Barth, AGU, 1995]

              if (geo_alt < 120.) then
                 NDensityS(iLon,iLat,iAlt,iNO_,iBlock)=  &
                      1e14-1e10*abs((geo_alt-110.0))**3.5
                      !10**(-0.003*(geo_alt-105.)**2 +14+LOG10(3.))
              else 
                 m = (1e10-3.9e13)/(200)
                 k = 1e10+(-m*300.) 
                 NDensityS(iLon,iLat,iAlt,iNO_,iBlock)=  &
                      MAX(k+(m*geo_alt)-(geo_alt - 120.0)**2,1.0)
                   !   MAX(10**(13.-LOG10(3.)*(geo_alt-165.)/35.),1.0)
              endif

              LogNS(iLon,iLat,iAlt,:,iBlock) = &
                   log(NDensityS(iLon,iLat,iAlt,iNO_,iBlock))

              NDensity(iLon,iLat,iAlt,iBlock) = &
                   sum(NDensityS(iLon,iLat,iAlt,1:nSpecies,iBlock))

              hwm_utime = utime
              hwm_alt = geo_alt
              hwm_lat = geo_lat
              hwm_lon = geo_lon
              hwm_lst = geo_lst
              hwm_f107a = f107a
              hwm_f107 = f107
              hwm_ap(1) = -1.0
              hwm_ap(2) = -1.0
              
              call HWM07(iyd,hwm_utime,hwm_alt,hwm_lat,hwm_lon,hwm_lst,&
                   hwm_f107a,hwm_f107,hwm_ap,qw)

              ! qw is north&east
              Velocity(iLon,iLat,iAlt,iEast_,iBlock) = qw(2)
              Velocity(iLon,iLat,iAlt,iNorth_,iBlock) = qw(1)

           enddo
        enddo
     enddo

     Rho(:,:,:,iBlock) = 0.0
     NDensity(:,:,:,iBlock) = 0.0

     do iSpecies=1,nSpecies

        NDensity(:,:,:,iBlock) = NDensity(:,:,:,iBlock) + &
             NDensityS(:,:,:,iSpecies,iBlock)

        Rho(:,:,:,iBlock) = Rho(:,:,:,iBlock) + &
             Mass(iSpecies)*NDensityS(:,:,:,iSpecies,iBlock)
        
     enddo

  enddo
 
end subroutine init_msis

!--------------------------------------------------------------
!
!--------------------------------------------------------------

subroutine msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
     F107A,F107,AP,LogNS, Temp, LogRho, v)

  use ModTime, only : iTimeArray
  use ModPlanet

  use EUA_ModMsis90, ONLY: gtd6

  implicit none

  integer, intent(in) :: iJulianDay
  real, intent(in) :: uTime, Alt, Lat, Lon, LST, f107a, f107
  real, intent(in):: ap
  real, intent(out) :: LogNS(nSpecies), Temp, LogRho, v(2)

  real :: msis_temp(2)
  real :: msis_dens(8)
  real :: AP_I(7), ffactor, no
  integer :: iyd

  real(4) :: hwm_utime, hwm_alt, hwm_lat, hwm_lon, hwm_lst
  real(4) :: hwm_f107a, hwm_f107, hwm_ap(2), qw(2)

  !----------------------------------------------------------------------------
  AP_I = AP
  CALL GTD6(iJulianDay,uTime,Alt,Lat,Lon,LST, &
       F107A,F107,AP_I,48,msis_dens,msis_temp)

  LogNS(iO_3P_)  = alog(msis_dens(2))
  LogNS(iO2_) = alog(msis_dens(4))
  LogNS(iN2_) = alog(msis_dens(3))
  if (nSpecies >= iN_4S_) LogNS(min(nSpecies,iN_4S_)) = alog(msis_dens(8))

  if (nSpecies >= iNO_) then
     ffactor = 6.36*log(f107)-13.8
     no = ffactor * 1.0e13 + 8.0e13
     LogNS(min(nSpecies,iNO_))   = alog(no)
  endif

  Temp        = msis_temp(2)
  LogRho      = alog(msis_dens(6))

  iyd = iTimeArray(1)*1000 + iJulianDay
  hwm_utime = utime
  hwm_alt = alt
  hwm_lat = lat
  hwm_lon = lon
  hwm_lst = lst
  hwm_f107a = f107a
  hwm_f107 = f107
  hwm_ap(1) = -1.0
  hwm_ap(2) = -1.0

  call HWM07(iyd,hwm_utime,hwm_alt,hwm_lat,hwm_lon,hwm_lst,&
       hwm_f107a,hwm_f107,hwm_ap,qw)

  ! qw is north&east
  V(1) = qw(2)
  V(2) = qw(1)

end subroutine msis_bcs

