
! ----------------------------------------------------------------
! If you want to output some specific variables, then do that here.
! In ModUserGITM, there are three variables defined, UserData1D,UserData2D,UserData3D.
! To output variables:
! 1. Figure out which variable you want to output.
! 2. Go into the code where the variable is set and copy it into
!    UserData3D or UserData2D or UserData1D
! 3. Do this for each variable you want to output.
! 4. Edit output_header_user below, making a list of all the variables
!    that you have added.  Make sure you leave longitude, latitude, and
!    altitude in the list of variables.
! 5. Count the number of variables that you output (including
!    Lon, Lat, and Alt). Change nVarsUser3d or nVarsUser2d or nVarsUser1d in the
!    subroutines towards the top of this file.
! 6. If you add more than 40 variables, you probably should check
!    nUserOutputs in ModUserGITM.f90 and make sure that this number is
!    larger than the number of variables that you added.
! 7. Recompile and run. Debug. Repeat 7.
! ----------------------------------------------------------------

! ----------------------------------------------------------------
! S. W. BOUGHER revision for added UserData3D : 11-10-04, 11-10-26, 11-01-11, 11-07-11
! S. W. BOUGHER revision for added UserData1D : 11-10-25, 11-10-26  11-01-11, 11-07-11
! --3 location Variables
! --11 radiative balance terms for diagnostics
! --7-printed outputs for logfile
! S. W. BOUGHER revision for added UserData1D : 12-11-28
! --3 location Variables
! --11 radiative balance terms for diagnostics
! --8  printed outputs for logfile
! ----------------------------------------------------------------

subroutine set_nVarsUser3d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser3d = 4

  if (nVarsUser3d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser3d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser2d

  use ModUserGITM
  use ModSources, only:ED_N_Energies

  ! Make sure to include Lat, Lon, and Alt

  nVarsUser2d = 10+ED_N_Energies

  if (nVarsUser2d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser2d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine set_nVarsUser1d

  use ModUserGITM

  ! Make sure to include Lat, Lon, and Alt

	!MM adding variables for debugging extra ei processes
  nVarsUser1d = 23

  if (nVarsUser1d-3 > nUserOutputs) &
       call stop_gitm("Too many user outputs!! Increase nUserOutputs!!")

end subroutine set_nVarsUser1d

! ----------------------------------------------------------------
!
! ----------------------------------------------------------------

subroutine output_header_user(cType, iOutputUnit_)

  use ModUserGITM
  use ModSources, only:ED_Energies, ED_N_Energies

  implicit none

  character (len=5), intent(in) :: cType
  integer, intent(in)           :: iOutputUnit_
  integer :: n

  ! ------------------------------------------
  ! 3D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '3D') then

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser3d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats+4, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons+4, " nLongitudes"

     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Eddy Diffusion Coefficient"
     !     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "CO2ImpactIonizationFrequency"
     ! write(iOutputUnit_,"(I7,A1,a)")  6, " ", "EUVHeatingRate"
     ! write(iOutputUnit_,"(I7,A1,a)")  7, " ", "QNIRTotalRate"
     ! write(iOutputUnit_,"(I7,A1,a)")  8, " ", "QNIRLTERate"
     ! write(iOutputUnit_,"(I7,A1,a)")  9, " ", "CIRLTERate"
     ! write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Conduction"
     ! write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Pressure(ubar)"
     ! write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Temperature"
     ! write(iOutputUnit_,"(I7,A1,a)") 13, " ", "VMRO1"
     ! write(iOutputUnit_,"(I7,A1,a)") 14, " ", "VMRCO2"

  endif

  ! ------------------------------------------
  ! 2D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '2D') then

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser2d, " nvars"
     write(iOutputUnit_,"(I7,7A)")     1, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)") nLats, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)") nLons, " nLongitudes"

     write(iOutputUnit_,*) ""
     write(iOutputUnit_,*) "NO GHOSTCELLS"
     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Potential (kV)"
     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "Total Energy (ergs)"
     write(iOutputUnit_,"(I7,A1,a)")  7, " ", "Discrete Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)")  8, " ", "Discrete Total Energy (ergs)"
     write(iOutputUnit_,"(I7,A1,a)")  9, " ", "Wave Average Energy (keV)"
     write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Wave Total Energy (ergs)"
     do n=1,ED_N_Energies
        write(iOutputUnit_,"(I7,A6,1P,E9.3,A11)") 10+n, " Flux@",ED_energies(n), "eV (/cm2/s)"
     enddo
  endif

  ! ------------------------------------------
  ! 1D Output Header
  ! ------------------------------------------

  if (cType(1:2) == '1D') then

     write(iOutputUnit_,*) "NUMERICAL VALUES"
     write(iOutputUnit_,"(I7,6A)") nVarsUser1d, " nvars"
     write(iOutputUnit_,"(I7,7A)") nAlts+4, " nAltitudes"
     write(iOutputUnit_,"(I7,7A)")       1, " nLatitudes"
     write(iOutputUnit_,"(I7,7A)")       1, " nLongitudes"

     write(iOutputUnit_,*) ""

     write(iOutputUnit_,*) "VARIABLE LIST"
     write(iOutputUnit_,"(I7,A1,a)")  1, " ", "Longitude"
     write(iOutputUnit_,"(I7,A1,a)")  2, " ", "Latitude"
     write(iOutputUnit_,"(I7,A1,a)")  3, " ", "Altitude"
     !write(iOutputUnit_,"(I7,A1,a)")  4, " ", "Eddy Diffusion Coefficient"
     !     write(iOutputUnit_,"(I7,A1,a)")  5, " ", "CO2+hv->2O+C"
     !     write(iOutputUnit_,"(I7,A1,a)")  6, " ", "CO+hv->C+O"
     !write(iOutputUnit_,"(I7,A1,a)")  7, " ", "QNIRTotalRate"
     !write(iOutputUnit_,"(I7,A1,a)")  8, " ", "QNIRLTERate"
     !write(iOutputUnit_,"(I7,A1,a)")  9, " ", "CIRLTERate"
     !write(iOutputUnit_,"(I7,A1,a)") 10, " ", "Conduction"
     !write(iOutputUnit_,"(I7,A1,a)") 11, " ", "Pressure(ubar)"
     !write(iOutputUnit_,"(I7,A1,a)") 12, " ", "Temperature"
     !write(iOutputUnit_,"(I7,A1,a)") 13, " ", "VMRO1"
     !write(iOutputUnit_,"(I7,A1,a)") 14, " ", "VMRCO2"
     !MM writing reaction rates
     write(iOutputUnit_,"(I7,A1,a)") 4, " ", "|B|"
     write(iOutputUnit_,"(I7,A1,a)") 5, " ", "Magnetic Elevation Angle"
     write(iOutputUnit_,"(I7,A1,a)") 6, " ", "Magnetic Topology"
     write(iOutputUnit_,"(I7,A1,a)") 7, " ", "CO2_X2PI_G_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 8, " ", "CO2_A2PI_U_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 9, " ", "CO2_B2Sig+u_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 10, " ", "CO2_CO+_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 11, " ", "CO2_O+_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 12, " ", "CO2_C+_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 13, " ", "O1_+4S_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 14, " ", "O1_+2D_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 15, " ", "O1_+2P_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 16, " ", "N2_+X2Sg_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 17, " ", "N2_+A2Pu_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 18, " ", "N2_+B2Su_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 19, " ", "CO_CO+X2Sig_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 20, " ", "CO_COMETTAI_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 21, " ", "CO_1STNEG_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 22, " ", "CO_C+_ionization"
     write(iOutputUnit_,"(I7,A1,a)") 23, " ", "CO_O+_ionization"

  endif


  write(iOutputUnit_,*) ""

end subroutine output_header_user

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_3dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  do iAlt=-1,nAlts+2
     do iLat=-1,nLats+2
        do iLon=-1,nLons+2
           write(iOutputUnit_)       &
                Longitude(iLon,iBlock), &
                Latitude(iLat,iBlock), &
                Altitude_GB(iLon, iLat, iAlt, iBlock),&
                UserData3D(iLon,iLat,iAlt,1,iBlock)
           !                UserData3D(iLon,iLat,iAlt,2,iBlock),&
           !                UserData3D(iLon,iLat,iAlt,3,iBlock)

        enddo
     enddo
  enddo

end subroutine output_3dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_2dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon

  iAlt = 1
  do iLat=1,nLats
     do iLon=1,nLons
        write(iOutputUnit_)       &
             Longitude(iLon,iBlock), &
             Latitude(iLat,iBlock), &
             Altitude_GB(iLon, iLat, iAlt, iBlock),&
             UserData2D(iLon,iLat,iAlt,1:nVarsUser2d-3,iBlock)
     enddo
  enddo

end subroutine output_2dUser

!----------------------------------------------------------------
!
!----------------------------------------------------------------

subroutine output_1dUser(iBlock, iOutputUnit_)

  use ModGITM
  use ModUserGITM

  implicit none

  integer, intent(in) :: iBlock, iOutputUnit_
  integer :: iAlt, iLat, iLon,iiAlt

  iLat = 1
  iLon = 1
  !  Plotter outputs
  do iAlt=-1,nAlts+2
     iiAlt = max(min(iAlt,nAlts),1)
     write(iOutputUnit_)       &
          Longitude(iLon,iBlock), &
          Latitude(iLat,iBlock), &
          Altitude_GB(iLon, iLat, iAlt, iBlock),&
          UserData1D(iLon,iLat,iAlt,1:20)
     !                UserData1D(iLon,iLat,iAlt,2), &
     !                UserData1D(iLon,iLat,iAlt,3)
     !                UserData1D(iLon,iLat,iAlt,10), &
     !                UserData1D(iLon,iLat,iAlt,11)
     !               UserData1D(iLon,iLat,iAlt,1:nVarsUser1d-3)
  enddo

end subroutine output_1dUser
