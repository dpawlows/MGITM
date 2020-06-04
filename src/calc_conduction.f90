
subroutine calc_conduction(iBlock, Quantity, Diff, MulFac, dTdt_cond)

! New Method, A. Ridley fixes (May 2012). From cvs on 120503.

  use ModSizeGitm
  use ModGITM, only: dAlt_GB, Latitude, Longitude, dt, Altitude_GB, &
       RadialDistance_GB
  use ModConstants

  implicit none

  integer, intent(in) :: iBlock
  real, intent(in) :: Quantity(nLons, nLats, -1:nAlts+2)
  real, intent(in) :: Diff(nLons, nLats, 0:nAlts+1)
  real, intent(in) :: MulFac(nLons, nLats,0:nAlts+1)
  real, intent(out) :: dTdt_cond(nLons, nLats, nAlts)

  real, dimension(0:nAlts+1) :: m, du, r, du12, du22, &
       dl, lou, dlou, di, r2

  real :: tempold(0:nAlts+1), temp(0:nAlts+1)
  real, dimension(0:nAlts+1) :: a,b,c,d, cp, dp
  integer :: iLon, iLat, iAlt

  call start_timing("conduction")
  call report("calc_conduction",3)

  do iLon = 1, nLons
     do iLat = 1, nLats

        tempold = Quantity(iLon, iLat, 0:nAlts+1)

        r2 = RadialDistance_GB(iLon,iLat, 0:nAlts+1,iBlock)**2
        di = diff(iLon,iLat,:)*r2

        m = dt/(MulFac(iLon, iLat, 0:nAlts+1)*r2)
        du = Altitude_GB(iLon,iLat, 1:nAlts+2,iBlock) - &
             Altitude_GB(iLon,iLat, 0:nAlts+1,iBlock)
        dl = Altitude_GB(iLon,iLat, 0:nAlts+1,iBlock) - &
             Altitude_GB(iLon,iLat,-1:nAlts+0,iBlock)
        r = du/dl

        du12 = du*du * (1+r*r)
        du22 = 0.5 * (dl*du + du*du)

        lou  = di/du22
        dlou = di/du22

        dl(1:nAlts) =                            di(2:nAlts+1) - &
                         r(1:nAlts)*r(1:nAlts) * di(0:nAlts-1) - &
                      (1-r(1:nAlts)*r(1:nAlts))* di(1:nAlts  )

        dl(0) = dl(1)
        dl(nAlts+1) = dl(nAlts)

        dl = 0.0

        ! Do a google search for a tri-diagnal solver and you will come up
        ! with this:

        a =  di/du22*r - dl/du12 * r*r

        c =  di/du22 + dl/du12

        b = -1/m - di/du22*(1+r) - dl/du12*(1-r*r) 

        d = -tempold/m

        ! Boundary Conditions:

        d(0) = -tempold(0)
        a(0) = 0.0
        b(0) = -1.0
        c(0) = 0.0

        a(nAlts+1) = 1.0
        b(nAlts+1) = -1.0
        c(nAlts+1) = 0.0
        d(nAlts+1) = 0.0

        cp(0) = c(0)/b(0)
        do iAlt = 1, nAlts+1
           cp(iAlt) = c(iAlt)/(b(iAlt)-cp(iAlt-1)*a(iAlt))
        enddo
        dp(0) = d(0)/b(0)
        do iAlt = 1, nAlts+1
           dp(iAlt) = (d(iAlt)-dp(iAlt-1)*a(iAlt))/(b(iAlt)-cp(iAlt-1)*a(iAlt))
        enddo
        temp(nAlts+1) = dp(nAlts+1)
        do iAlt=nAlts,0,-1
           temp(iAlt) = dp(iAlt)-cp(iAlt)*temp(iAlt+1)
        enddo

        dTdt_cond(iLon,iLat,1:nAlts) = temp(1:nAlts) - &
             Quantity(iLon,iLat,1:nAlts)

     enddo
  enddo

  call end_timing("conduction")

end subroutine calc_conduction
