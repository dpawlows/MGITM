!\
! ------------------------------------------------------------
! set_boundary
! ------------------------------------------------------------
!/

subroutine set_vertical_bcs(LogRho,LogNS,Vel_GD,Temp, LogINS, iVel, VertVel)

  ! Fill in ghost cells at the top and bottom

  use ModSizeGitm, only: nAlts
  use ModPlanet, only: nSpecies, nIonsAdvect, Mass, nIons, IsEarth,iN2_, iH2_, iCH4_
  use ModGITM, only: TempUnit, iEast_, iNorth_, iUp_
  use ModInputs
  use ModConstants
  use ModTime, only: UTime, iJulianDay,currenttime
  use ModVertical, only: Lat, Lon, Gravity_G, Altitude_G, dAlt_F
  use ModIndicesInterfaces, only: get_HPI

  use EUA_ModMsis90, ONLY: meter6

  implicit none

  real, intent(inout) :: &
       LogRho(-1:nAlts+2), &
       LogNS(-1:nAlts+2,nSpecies), &
       LogINS(-1:nAlts+2,nIonsAdvect), &
       Vel_GD(-1:nAlts+2,3), &
       IVel(-1:nAlts+2,3), &
       Temp(-1:nAlts+2), &
       VertVel(-1:nAlts+2,nSpecies)

  integer :: iSpecies, iAlt
  real    :: InvScaleHeightS, InvScaleHgt, Alt, Lst, Ap = 4.0, dn, dt
  logical :: IsFirstTime = .true., UseMsisBCs = .false.
  real    :: HP
  integer :: ierror

  integer, dimension(25) :: sw

  !-----------------------------------------------------------
  ! Bottom
  !-----------------------------------------------------------

  if (IsEarth) UseMsisBCs = UseMsis

  if (IsFirstTime .and. UseMsisBCs) then
     call meter6(.true.)
     sw = 1
     IsFirstTime = .true.
  endif

  if (UseMsisBCs) then
     call get_HPI(CurrentTime, HP, iError)  
     if (iError > 0) hp = 40.0
     Ap = min(200.,max(-40.72 + 1.3 * HP, 10.))
     do iAlt = -1, 0
        Alt = Altitude_G(iAlt)/1000.0
        Lst = mod(UTime/3600.0+Lon/15.0,24.0)

        call msis_bcs(iJulianDay,UTime,Alt,Lat,Lon,Lst, &
             F107A,F107,AP,LogNS(iAlt,:), Temp(iAlt), &
             LogRho(iAlt))

     enddo
  else

     ! do nothing - which means don't change the initial condition.

  endif


  ! Let the winds blow !!!!
!  Vel_GD(-1:0,iEast_)  = 0.0
!  Vel_GD(-1:0,iNorth_) = 0.0

!! Titan's Winds
   Vel_GD(0,iNorth_) = 2.0*sqrt(cos(Lat*PI/180.0))
   Vel_GD(-1,iNorth_) = 2.0*sqrt(cos(Lat*PI/180.0))

  Vel_GD(-1:0,iUp_)    = 0.0
  VertVel(-1:0,:)      = 0.0
  IVel(-1:0,iUp_)      = 0.0

!  VertVel(0 ,1:nSpecies)     =  -1.0*VertVel(1,1:nSpecies)
!  VertVel(-1,1:nSpecies)     =  -1.0*VertVel(2,1:nSpecies)

  do iSpecies=1,nIonsAdvect
     dn = (LogINS(2,iSpecies) - LogINS(1,iSpecies))
     
     LogINS(0,iSpecies) = LogINS(1,iSpecies) - dn
     LogINS(-1,iSpecies) = LogINS(0,iSpecies) - dn
  enddo

!  ! Lower boundary for NO on Earth
!  if (nSpecies == iNO_) then
!     dn = (LogNS(2,nSpecies) - LogNS(1,nSpecies))
!     if (dn >= 0) then
!        LogNS(0,nSpecies) = LogNS(1,nSpecies) - dn
!        LogNS(-1,nSpecies) = LogNS(0,nSpecies) - dn
!     else
!        LogNS(0,nSpecies) = LogNS(1,nSpecies) + dn
!        LogNS(-1,nSpecies) = LogNS(0,nSpecies) + dn
!     endif
!  endif

!  dn = (LogNS(2,iN_4S_) - LogNS(1,iN_4S_))
!  LogNS(0,iN_4S_) = LogNS(1,iN_4S_) - dn
!  LogNS(-1,iN_4S_) = LogNS(0,iN_4S_) - dn

  ! If you want to slip:
  !  Vel_GD(-1:0,iEast_)  = Vel_GD(1,iEast_)
  !  Vel_GD(-1:0,iNorth_) = Vel_GD(1,iNorth_)

  !-----------------------------------------------------------
  ! Top
  !-----------------------------------------------------------

  ! Slip flow at the top

  Vel_GD(nAlts+1:nAlts+2,iEast_)  = Vel_GD(nAlts,iEast_)
  Vel_GD(nAlts+1:nAlts+2,iNorth_) = Vel_GD(nAlts,iNorth_)

  IVel(nAlts+1:nAlts+2,iEast_)  = IVel(nAlts,iEast_)
  IVel(nAlts+1:nAlts+2,iNorth_) = IVel(nAlts,iNorth_)

  ! Things can go up or down in the ions

  IVel(nAlts+1,iUp_)   =  IVel(nAlts,iUp_)
  IVel(nAlts+2,iUp_)   =  IVel(nAlts-1,iUp_)

  ! We only let stuff flow out in the neutrals

  if(Vel_GD(nAlts,iUp_)>0.)then
     Vel_GD(nAlts+1:nAlts+2,iUp_) = Vel_GD(nAlts,iUp_)
     VertVel(nAlts+1,:) = VertVel(nAlts,:)
     VertVel(nAlts+2,:) = VertVel(nAlts,:)
  else
    ! Vel_GD(nAlts+1:nAlts+2,iUp_) = 0.0 ! -Vel(nAlts)
     Vel_GD(nAlts+1,iUp_) = -Vel_GD(nAlts,iUp_)
     Vel_GD(nAlts+2,iUp_) = -Vel_GD(nAlts-1,iUp_)
     VertVel(nAlts+1,:) = -VertVel(nAlts,:)
     VertVel(nAlts+2,:) = -VertVel(nAlts-1,:)
  endif

     VertVel(nAlts+1,iH2_) = 100.0
     VertVel(nAlts+2,iH2_) = 100.0

     VertVel(nAlts+1,iCH4_) = 9.0
     VertVel(nAlts+2,iCH4_) = 9.0

  ! Constant temperature (zero gradient)

  Temp(nAlts+1) = Temp(nAlts)
  Temp(nAlts+2) = Temp(nAlts)

  dn = (LogRho(nAlts) - LogRho(nAlts-1))
  LogRho(nAlts+1) = LogRho(nAlts) + dn
  LogRho(nAlts+2) = LogRho(nAlts+1) + dn

  ! Limit the slope of the ion density

  do iSpecies=1,nIonsAdvect
     dn = (LogINS(nAlts,iSpecies) - LogINS(nAlts-1,iSpecies))
!     if (dn < 0.75*LogINS(nAlts,iSpecies) .and. dn > 0) &
!          dn = 0.75*LogINS(nAlts,iSpecies)
     if (dn > 0) dn = -0.1*LogINS(nAlts,iSpecies)
     LogINS(nAlts+1,iSpecies) = LogINS(nAlts,iSpecies) + dn
     LogINS(nAlts+2,iSpecies) = LogINS(nAlts+1,iSpecies) + dn
  enddo

  ! Hydrostatic pressure for the neutrals

  do iSpecies=1,nSpecies
     do iAlt = nAlts+1, nAlts+2
        InvScaleHeightS = -Gravity_G(iAlt) * &
             Mass(iSpecies) / (Temp(iAlt)*Boltzmanns_Constant)
        LogNS(iAlt,iSpecies) = &
             LogNS(iAlt-1,iSpecies) - dAlt_F(iAlt)*InvScaleHeightS
        if (LogNS(nAlts+1,iSpecies) > 75.0 .or. &
             LogNS(nAlts+2,iSpecies) > 75.0) then
           write(*,*) "======> bcs : ", iSpecies, 1.0e-3/InvScaleHeightS, &
                Gravity_G(nAlts), Mass(iSpecies), Temp(nAlts), &
                LogNS(nAlts,iSpecies), LogNS(nAlts+1,iSpecies), &
                dAlt_F(nAlts), LogNS(nAlts+2,iSpecies)
        endif
     enddo
  enddo

end subroutine set_vertical_bcs

