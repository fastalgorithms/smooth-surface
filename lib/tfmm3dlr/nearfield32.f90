!c
!c**********************************************************************
!      subroutine tpotfld3dhess_dp(iffld,ifhess,source,qwt,dipvec,
!     1           target,sigma,pot,fld,hess)
!c**********************************************************************
!c
!c     This subroutine calculates the potential POT field FLD and
!c     Hessian HESS at the target point TARGET, due to a dipole at
!c     SOURCE. The scaling is that required of the delta function
!c     response: i.e.,
!c
!c                  pot = (dipvec(1) x + dipvec(2) y + dipvec(3) z)/r^3
!c        fld = -grad(pot)
!c        hess = (potxx,potyy,potzz,potxy,potxz,potyz)
!c
!c     INPUT:
!c
!c     iffld        : flag for computing gradient
!c                         ffld = 0 -> dont compute
!c                        ffld = 1 -> do compute
!c     ifhess       : flag for computing Hessian
!c                         ifhess = 0 -> dont compute
!c                        ifhess = 1 -> do compute
!c     source(3)  : location of the source
!c     qwt           : quadrature weight
!c     dipvec(3)   : dipole direction
!c     target(3)    : location of the target
!c     sigma        : Gaussian parameter
!c
!c----------------------------------------------------------------------
!c     OUTPUT:
!c
!c     pot          : calculated potential
!c     fld           : calculated -gradient
!c     hess        : calculated hessian
!c

subroutine tpotfld3d_dp(iffld,source,qwt,dipvec,targ,sigma,grad_sigma,pot,fld)
implicit none
!List of calling arguments
integer ( kind = 8 ), intent(in) :: iffld
real ( kind = 8 ) , intent(in) :: source(3),qwt,dipvec(3),targ(3),sigma,grad_sigma(3)
real ( kind = 8 ), intent(out) ::  pot,fld(3)

!List of local variables
real ( kind = 8 ) R, H, HH, HHH, prod

    R=sqrt((source(1)-targ(1))**2+(source(2)-targ(2))**2+(source(3)-targ(3))**2)
    call my_erf_like_v4(R,sigma,H,HH,HHH)
    prod=(source(1)-targ(1))*dipvec(1)+(source(2)-targ(2))*dipvec(2)+(source(3)-targ(3))*dipvec(3)
    pot=H*prod*qwt
    fld(1)=-((HH*(targ(1)-source(1))*prod+HHH*prod*grad_sigma(1)-H*dipvec(1))*qwt)
    fld(2)=-((HH*(targ(2)-source(2))*prod+HHH*prod*grad_sigma(2)-H*dipvec(2))*qwt)
    fld(3)=-((HH*(targ(3)-source(3))*prod+HHH*prod*grad_sigma(3)-H*dipvec(3))*qwt)

return
end




subroutine my_erf_like_v4(r,sgma,H,HH,HHH)
implicit none

!List of calling arguments
real ( kind = 8 ) , intent(in) :: r,sgma
real ( kind = 8 ), intent(out) ::  H, HH, HHH
double precision :: uuu2, uuu4, uuu6, uuu8

!List of local variables
real ( kind = 8 ) pi, my_exp, denom,x(15),w(15),uuu
integer ( kind = 8) count,nlg

    pi=3.141592653589793238462643383d0
    my_exp=exp(-r**2/(2*sgma**2))
    HHH=-1/(4*pi)*sqrt(2/pi)*(1/sgma**4)*my_exp

    uuu = r/(sqrt(2.0d0)*sgma)
    if (r==0.0d0) then
      H=sqrt(2.0d0)/(12*sgma**3*sqrt(pi**3))
      HH=-1.0d0*sqrt(2.0d0)/(20.0d0*sgma**5.0d0*sqrt(pi**3))
    else if (uuu<0.1d0) then
      uuu2 = uuu*uuu
      uuu4 = uuu2*uuu2
      uuu6 = uuu2*uuu4
      uuu8 = uuu2*uuu6
      H=0.0d0
      HH=0.0d0
      H = uuu8/132.0d0 - uuu6/27.0d0 + uuu4/7.0d0 - 2*uuu2/5.0d0 + 2.0d0/3.0d0
      H = H/(4*sqrt(2.0d0)*sqrt(pi**3)*(sgma**3))
    else
      H=(erf((sqrt(2.0d0)*r)/(2*sgma))/(4*r**2*pi) &
          - (sqrt(2.0d0)*my_exp)/(4*sgma*r*sqrt(pi**3)))/r
      denom=(4*r**5*sgma**3*sqrt(pi**3))
      HH=(sqrt(2.0d0)*r**3*my_exp - 3*sgma**3*sqrt(pi) &
          *erf((sqrt(2.0d0)*r)/(2*sgma)) &
          + 3*sqrt(2.0d0)*r*sgma**2*my_exp)/denom
    end if
    return
  end subroutine my_erf_like_v4

