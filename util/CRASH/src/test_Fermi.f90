program test_Alex
  use CRASH_ModFermiGas
  implicit none

  call test_fermi_function
end program test_Alex

subroutine CON_stop(StringError)
  implicit none
  character (len=*), intent(in) :: StringError
end subroutine CON_stop

subroutine CON_set_do_test(String,DoTest,DoTestMe)
  implicit none
  character (len=*), intent(in)  :: String
  logical,           intent(out) :: DoTest, DoTestMe
end subroutine CON_set_do_test

