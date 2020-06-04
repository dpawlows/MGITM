
subroutine fill_photo(photoion, photoabs, photodis)

  use ModPlanet
  use ModEUV

  implicit none

  real, intent(out) :: photoion(Num_WaveLengths_High, nIons-1)
  real, intent(out) :: photoabs(Num_WaveLengths_High, nSpecies)
  real, intent(out) :: photodis(Num_WaveLengths_High, nSpecies)

  integer :: iSpecies, iWave

  PhotoAbs = 0.0
  PhotoIon = 0.0
  PhotoDis = 0.0

  photoabs(:,iO_3P_)     = PhotoAbs_O
  photoabs(:,iO2_)    = PhotoAbs_O2

  if (nSpecies > 2) then
     iSpecies = iN2_
     photoabs(:,iSpecies)    = PhotoAbs_N2
  endif
  if (nSpecies > 3) then
     iSpecies = iN_4S_
     photoabs(:,min(iSpecies,nSpecies))    = PhotoIon_N
  endif

  ! This may need to be as defined below....
  photoion(:,iN2P_)   = PhotoIon_N2
  photoion(:,iO2P_)   = PhotoIon_O2
  photoion(:,iNP_)    = PhotoIon_N
  photoion(:,iO_4SP_) = PhotoIon_OPlus4S
  photoion(:,iO_2DP_) = PhotoIon_OPlus2D
  photoion(:,iO_2PP_) = PhotoIon_OPlus2P

  do iWave = 1, Num_WaveLengths_High
     if (waves(iWave) >= 1250.0 .and. wavel(iWave) <= 1750.0) then
        PhotoDis(iWave, iO2_) = &
             photoabs(iWave,iO2_) - PhotoIon(iWave, iO2P_)
     endif

     if (waves(iWave) >= 800.0 .and. wavel(iWave) <= 1250.0) then
        PhotoDis(iWave, iN2_) = &
             photoabs(iWave,iN2_) - PhotoIon(iWave, iN2P_)
     endif

  enddo

  !     Ion_Rate_Eff_N2(:,:,:) = Ion_Rate_Eff_N2(:,:,:) +              &
  !          Intensity(:,:,:,N) * (PhotoAbs_N2(N) - PhotoIon_N2(N)) *  &
  !          NDensity(:,:,:,N_N2,index)

end subroutine fill_photo

subroutine calc_planet_sources(iBlock)

  use ModInputs
  use ModSources
  use ModGITM
  use ModTime
  
  implicit none

  integer, intent(in) :: iBlock

  integer :: iAlt, iError, iDir, iLat, iLon

  real :: tmp2(nLons, nLats, nAlts)
  real :: tmp3(nLons, nLats, nAlts)
  real :: Omega(nLons, nLats, nAlts)

  LowAtmosRadRate = 0.0

  !\
  ! Cooling ----------------------------------------------------------
  !/

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> NO cooling", iproc, UseNOCooling

  if (UseNOCooling) then

     !  [NO] cooling 
     ! [Reference: Kockarts,G., G.R.L.,VOL.7, PP.137-140,Feberary 1980 ]
 
     Omega = 6.5e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) /      &
          (6.5e-17 * NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock) + 13.3)

     ! We need to check this out. I don't like the first / sign....

     NOCooling = Planck_Constant * Speed_Light / &
          5.3e-6 * &
          Omega * 13.3 *  &
          exp(- Planck_Constant * Speed_Light / &
          (5.3e-6 * Boltzmanns_Constant * &
          Temperature(1:nLons,1:nLats,1:nAlts,iBlock)* &
          TempUnit(1:nLons,1:nLats,1:nAlts))) * &
          NDensityS(1:nLons,1:nLats,1:nAlts,iNO_,iBlock)

     NOCooling = NOCooling / TempUnit(1:nLons,1:nLats,1:nAlts) / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     NOCooling = 0.0

  endif

  if (UseBarriers) call MPI_BARRIER(iCommGITM,iError)
  if (iDebugLevel > 4) write(*,*) "=====> UseOCooling", iproc, UseOCooling

  if (UseOCooling) then 

     ! [O] cooling 
     ! Reference: Kockarts, G., Plant. Space Sci., Vol. 18, pp. 271-285, 1970
     ! We reduce the LTE 63-um cooling rate by a factor of 2 for 
     ! the non-LTE effects.[Roble,1987]         

     tmp2 = exp(-228./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)))
     tmp3 = exp(-326./(Temperature(1:nLons,1:nLats,1:nAlts,iBlock)*&
          TempUnit(1:nLons,1:nLats,1:nAlts)))

     ! In erg/cm3/s
     OCooling = (1.69e-18*tmp2 + 4.59e-20*tmp3) * &
          (NDensityS(1:nLons,1:nLats,1:nAlts,iO_3P_,iBlock)/1.0e6) / &
          (1.0 + 0.6*tmp2 + 0.2*tmp3)
     ! In w/m3/3
     OCooling = OCooling/10.0
     ! In our special units:
     OCooling = OCooling/ TempUnit(1:nLons,1:nLats,1:nAlts) / &
          (Rho(1:nLons,1:nLats,1:nAlts,iBlock)*cp(:,:,1:nAlts,iBlock))

  else

     OCooling = 0.0

  endif

  RadCooling(1:nLons,1:nLats,1:nAlts,iBlock) = OCooling + NOCooling


!--------------------------------------------------------------------
! GLOW
!--------------------------------------------------------------------

if (UseGlow) then
     if (dt < 10000.) then
        if  (floor((tSimulation-dt)/DtGlow) /= &
             floor(tsimulation/DtGlow)) then   

           call start_timing("glow")
           isInitialGlow = .True.

           if (iDebugLevel > 4) write(*,*) "=====> going into get_glow", iproc

           do iLat = 1, nLats
              do iLon = 1, nLons

                 call get_glow(iLon,iLat,iBlock)
                 
              enddo
           enddo

           call end_timing("glow")

        endif
     endif
     PhotoElectronDensity(:,:,:,:,iBlock) = PhotoElectronRate(:,:,:,:,iBlock) * dt
  endif


end subroutine calc_planet_sources

!---------------------------------------------------------------------
! Initialize Heating Rates
!---------------------------------------------------------------------

subroutine init_heating_efficiency

  use ModGITM, only: nLons, nLats, nAlts, nBlocks, Altitude_GB
  use ModEUV, only: HeatingEfficiency_CB, eHeatingEfficiency_CB

  implicit none

  integer :: iLon, iLat, iAlt
  !------------------------------------------------------------------
  HeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.05
!  max(0.1, &
!       0.40 - &
!       5.56e-5*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 165)**2)

  where(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000. > 150.)
     eHeatingEfficiency_CB(:,:,:,1:nBlocks) = 0.04
!!! min(0.4, &
!!!     0.04 + &
!!!     0.05*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 150)/100)
  elsewhere        
     eHeatingEfficiency_CB(:,:,:,1:nBlocks) = max(0.000001, &
          0.05 + &
          0.07*(Altitude_GB(1:nLons,1:nLats,1:nAlts,1:nBlocks)/1000 - 200)/100)
  end where

end subroutine init_heating_efficiency

!---------------------------------------------------------------------
! Calculate Eddy Diffusion Coefficient
!---------------------------------------------------------------------

subroutine calc_eddy_diffusion_coefficient(iBlock)

  use ModSizeGITM
  use ModGITM, only: pressure
  use ModInputs, only: EddyDiffusionPressure0,EddyDiffusionPressure1, &
       EddyDiffusionCoef
  use ModSources, only: KappaEddyDiffusion

  implicit none

  integer, intent(in) :: iBlock
  integer :: iAlt, iLat, iLon

  KappaEddyDiffusion=0.
  do iAlt = -1, nAlts+2

     do iLat = 1, nLats
        do iLon = 1, nLons

           if (pressure(iLon,iLat,iAlt,iBlock) >EddyDiffusionPressure0) then
              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef
              
           else if (pressure(iLon,iLat,iAlt,iBlock) > &
                EddyDiffusionPressure1) then

              KappaEddyDiffusion(iLon,iLat,iAlt,iBlock) = EddyDiffusionCoef * &
                   (pressure(iLon,iLat,iAlt,iBlock) - &
                   EddyDiffusionPressure1)/&
                   (EddyDiffusionPressure0 - EddyDiffusionPressure1)

           endif
        enddo
     enddo
  enddo

end subroutine calc_eddy_diffusion_coefficient

subroutine set_planet_defaults

  use ModInputs

  return

end subroutine set_planet_defaults

subroutine planet_limited_fluxes(iBlock)
!! Do Nothing
end subroutine planet_limited_fluxes
