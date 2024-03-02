      implicit real *8 (a-h,o-z)
      real *8 src(3), dipvec(3), targ(3), sigma
      real *8 grad_sigma(3), pot, grad(3)
      real *8 potex, gradex(3)

      done = 1.0d0
      pi = atan(done)*4.0d0

      rfac = 2.0d0 
      src(1) = hkrand(0)*rfac
      src(2) = hkrand(0)*rfac
      src(3) = hkrand(0)*rfac

      targ(1) = hkrand(0)
      targ(2) = hkrand(0)
      targ(3) = hkrand(0)

      sigma = 0.9d0

      grad_sigma(1) = hkrand(0)
      grad_sigma(2) = hkrand(0)
      grad_sigma(3) = hkrand(0)

      dipvec(1) = hkrand(0)
      dipvec(2) = hkrand(0)
      dipvec(3) = hkrand(0)

      qwt = -4.0d0*pi
      iffld = 1
      call nearkernel(src, dipvec, targ, sigma, grad_sigma, pot, grad)

      call tpotfld3d_dp(iffld, src, qwt, dipvec, targ, sigma, 
     1   grad_sigma, potex, gradex)
      print *, pot, potex
      print *, "error in pot=",abs(pot-potex)
      
      print *, grad(1), gradex(1)
      print *, grad(2), gradex(2)
      print *, grad(3), gradex(3)

      err = abs(grad(1)+gradex(1)) + abs(grad(2) + gradex(2))
      err = err + abs(grad(3)+gradex(3))
      
      print *, "error in gradient =", err
      stop

      r = (src(1) - targ(1))**2 + (src(2)-targ(2))**2
      r = r + (src(3) - targ(3))**2
      r = sqrt(r)
      call surf_smooth_ker(r, sigma, h, hh, hhh)
      call my_erf_like_v4(r, sigma, h2, hh2, hhh2)

      print *, h, h2*4*pi
      print *, hh, hh2*4*pi
      print *, hhh, hhh2*4*pi


      return
      end
c
c
c
c
      subroutine nearkernel(src, dipvec, targ, sigma, grad_sigma, pot,
     1    grad)
      implicit real *8 (a-h,o-z)
      real *8 src(3), targ(3), dipvec(3), sigma, grad_sigma(3), pot
      real *8 grad(3)

      real *8 r, h, hh, hhh, prod
      real *8 dxyz(3)

      dxyz(1) = src(1) - targ(1)
      dxyz(2) = src(2) - targ(2)
      dxyz(3) = src(3) - targ(3)
      r = sqrt(dxyz(1)*dxyz(1) + dxyz(2)*dxyz(2) + dxyz(3)*dxyz(3))
      call surf_smooth_ker(r, sigma, h, hh, hhh)
      prod = dxyz(1)*dipvec(1) + dxyz(2)*dipvec(2) + dxyz(3)*dipvec(3)
      pot = -h*prod

      
      grad(1:3) = -prod*(-hh*dxyz(1:3) + hhh*grad_sigma(1:3)) + 
     1    h*dipvec(1:3)

      return
      end



      subroutine surf_smooth_ker(r, sgma, h, hh, hhh)
      implicit real *8 (a-h,o-z)
      real *8 dmexp
      real *8, parameter :: pi=3.141592653589793238462643383d0
      real *8, parameter :: fourpi = 12.5663706143592d0
      real *8, parameter :: foursqrtpi3 = 22.2733119873268d0
      real *8, parameter :: threesqrtpi = 5.31736155271655d0
      real *8, parameter :: sqpi = 1.7724538509055159d0
      real *8, parameter :: sq2 = 1.4142135623730951d0

    
      dmexp=exp(-r**2/(2*sgma**2))
      hhh=-sqrt(2/pi)*(1/sgma**4)*dmexp

      uuu = r/(sq2*sgma)
      if (r==0.0d0) then
        h = sq2/(3.0d0*sgma**3*sqpi)
        hh = -sq2/(5.0d0*sgma**5.0d0*sqpi)
      else if (uuu<0.01d0) then
        uuu2 = uuu*uuu
        uuu4 = uuu2*uuu2
        uuu6 = uuu2*uuu4
        uuu8 = uuu2*uuu6
        h = 0.0d0
        hh = 0.0d0
        h = uuu8/132.0d0 - uuu6/27.0d0 + uuu4/7.0d0 - 
     1      2*uuu2/5.0d0 + 2.0d0/3.0d0
        h = h/(sq2*sqpi*(sgma**3))
        hh = -uuu8/78.0d0 + 2*uuu6/33.0d0 - 2*uuu4/9 + 4*uuu2/7 -
     1     4.0d0/5.0d0
        hh = hh/((sq2*sgma)**5)
      else
        dmerf = erf(uuu) 
        r2 = r*r
        r3 = r2*r
        r5 = r2*r3
      
        h = dmerf/r3 - sq2*dmexp/(sqpi*sgma*r2)
        denom = sqpi*r5*sgma**3
        hh = (sq2*r3*dmexp - threesqrtpi*sgma**3*dmerf 
     1          + 3*sq2*r*sgma**2*dmexp)/denom
      endif


      return
      end

      subroutine tpotfld3d_dp(iffld,source,qwt,dipvec,targ,sigma, 
     1   grad_sigma,pot,fld)
      implicit none
c      List of calling arguments
      integer, intent(in) :: iffld
      double precision , intent(in) :: source(3),qwt,dipvec(3)
      double precision :: targ(3),sigma,grad_sigma(3)
      double precision, intent(out) ::  pot,fld(3)

      !List of local variables
      double precision :: R, H, HH, HHH, prod

      R=sqrt((source(1)-targ(1))**2+(source(2)-targ(2))**2 
     1     +(source(3)-targ(3))**2)
      call my_erf_like_v4(R,sigma,H,HH,HHH)
      prod=(source(1)-targ(1))*dipvec(1)+(source(2)-targ(2))*dipvec(2) 
     1   +(source(3)-targ(3))*dipvec(3)
      pot=H*prod*qwt
      fld(1)=-((HH*(targ(1)-source(1))*prod+HHH*prod*grad_sigma(1) 
     1    -H*dipvec(1))*qwt)
      fld(2)=-((HH*(targ(2)-source(2))*prod+HHH*prod*grad_sigma(2) 
     1   -H*dipvec(2))*qwt)
      fld(3)=-((HH*(targ(3)-source(3))*prod+HHH*prod*grad_sigma(3) 
     1     -H*dipvec(3))*qwt)

      return
      end subroutine tpotfld3d_dp




      subroutine my_erf_like_v4(r,sgma,H,HH,HHH)
      implicit none

c      List of calling arguments
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
          H=H+(x(count)**2)*exp(-x(count)**2) 
     1         /(sqrt(pi**3))*w(count)/(r**3)
          HH=HH+(x(count)**4)*exp(-x(count)**2)*-2.0d0*w(count) 
     1         /(r**5*sqrt(pi**3))
        enddo
      else
        H=(erf((sqrt(2.0d0)*r)/(2*sgma))/(4*r**2*pi) 
     1     - (sqrt(2.0d0)*my_exp)/(4*sgma*r*sqrt(pi**3)))/r
        denom=(4.*r**5*sgma**3*sqrt(pi**3))
        HH=(sqrt(2.0d0)*r**3*my_exp 
     1     - 3*sgma**3*sqrt(pi)*erf((sqrt(2.0d0)*r)/(2*sgma)) 
     1     + 3*sqrt(2.0d0)*r*sgma**2*my_exp)/denom
      end if
      return
      end subroutine my_erf_like_v4

