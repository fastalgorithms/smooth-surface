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

subroutine tpotfld3d_dp(iffld,source,qwt,dipvec,targ,sigma, &
    grad_sigma,pot,fld)
  implicit none
  !List of calling arguments
  integer, intent(in) :: iffld
  double precision , intent(in) :: source(3),qwt,dipvec(3)
  double precision :: targ(3),sigma,grad_sigma(3)
  double precision, intent(out) ::  pot,fld(3)

  !List of local variables
  double precision :: R, H, HH, HHH, prod

  R=sqrt((source(1)-targ(1))**2+(source(2)-targ(2))**2 &
      +(source(3)-targ(3))**2)
  call my_erf_like_v4(R,sigma,H,HH,HHH)
  prod=(source(1)-targ(1))*dipvec(1)+(source(2)-targ(2))*dipvec(2) &
      +(source(3)-targ(3))*dipvec(3)
  pot=H*prod*qwt
  fld(1)=-((HH*(targ(1)-source(1))*prod+HHH*prod*grad_sigma(1) &
      -H*dipvec(1))*qwt)
  fld(2)=-((HH*(targ(2)-source(2))*prod+HHH*prod*grad_sigma(2) &
      -H*dipvec(2))*qwt)
  fld(3)=-((HH*(targ(3)-source(3))*prod+HHH*prod*grad_sigma(3) &
      -H*dipvec(3))*qwt)

  return
end subroutine tpotfld3d_dp




subroutine my_erf_like_v4(r,sgma,H,HH,HHH)
implicit none

!List of calling arguments
double precision , intent(in) :: r,sgma
double precision, intent(out) ::  H, HH, HHH

!List of local variables
double precision pi, my_exp, denom,x(15),w(15)
integer  count

    pi=3.141592653589793238462643383d0
    my_exp=exp(-r**2/(2*sgma**2))
    HHH=-1/(4*pi)*sqrt(2/pi)*(1/sgma**4)*my_exp

    if (r==0.0d0) then
        H=sqrt(2.0d0)/(12*sgma**3*sqrt(pi**3))
        HH=-1.0d0*sqrt(2.0d0)/(20.0d0*sgma**5.0d0*sqrt(pi**3))
    else if (r/sgma<0.9d0) then
        x(1)=6.003740989757256d-03
        x(2)=3.136330379964697d-02
        x(3)=7.589670829478640d-02
        x(4)=1.377911343199150d-01
        x(5)=2.145139136957306d-01
        x(6)=3.029243264612183d-01
        x(7)=3.994029530012828d-01
        x(8)=5.000000000000000d-01
        x(9)=6.005970469987173d-01
        x(10)=6.970756735387817d-01
        x(11)=7.854860863042694d-01
        x(12)=8.622088656800850d-01
        x(13)=9.241032917052137d-01
        x(14)=9.686366962003530d-01
        x(15)=9.939962590102427d-01
        w(1)=1.537662099805856d-02
        w(2)=3.518302374405402d-02
        w(3)=5.357961023358594d-02
        w(4)=6.978533896307718d-02
        w(5)=8.313460290849692d-02
        w(6)=9.308050000778104d-02
        w(7)=9.921574266355579d-02
        w(8)=1.012891209627806d-01
        w(9)=9.921574266355579d-02
        w(10)=9.308050000778104d-02
        w(11)=8.313460290849692d-02
        w(12)=6.978533896307718d-02
        w(13)=5.357961023358594d-02
        w(14)=3.518302374405402d-02
        w(15)=1.537662099805856d-02
        x=x*r/sgma/sqrt(2.0d0)
        w=w*r/sgma/sqrt(2.0d0)
        H=0.0d0
        HH=0.0d0
        do count=1,15
!            x(count)=x(count)*r/sgma/sqrt(2.0d0)
!            w(count)=w(count)*r/sgma/sqrt(2.0d0)
          H=H+(x(count)**2)*exp(-x(count)**2) &
              /(sqrt(pi**3))*w(count)/(r**3)
          HH=HH+(x(count)**4)*exp(-x(count)**2)*-2.0d0*w(count) &
              /(r**5*sqrt(pi**3))
        enddo
    else
      H=(erf((sqrt(2.0d0)*r)/(2*sgma))/(4*r**2*pi) &
          - (sqrt(2.0d0)*my_exp)/(4*sgma*r*sqrt(pi**3)))/r
      denom=(4.*r**5*sgma**3*sqrt(pi**3))
      HH=(sqrt(2.0d0)*r**3*my_exp &
          - 3*sgma**3*sqrt(pi)*erf((sqrt(2.0d0)*r)/(2*sgma)) &
          + 3*sqrt(2.0d0)*r*sgma**2*my_exp)/denom
    end if
    return
  end subroutine my_erf_like_v4

