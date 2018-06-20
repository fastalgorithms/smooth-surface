
!
!
! This is a collection of routines for computing values, derivatives,
! and second derivatives of Koornwinder polynomials. The relevant
! routines are:
!
!
!     ortho2eva6 - returns values, first and second partial
!         derivatives of the Koornwinder polynomials.
!
!     klegeypols6 - returns values, first and second partial
!         derivatives of the functions P_n(x/y) y^n. These functions
!         make up one part of the Koornwinder polynomials
!
!     kjacopols3 - returns values, first and second derivatives
!         of Jacbobi polynomials. These functions are the other part
!         of Koornwinder polynimials
!
!     
!


subroutine ortho2eva6(mmax, z, pols, dersx, dersy, dersxx, &
    dersxy, dersyy, w)
  implicit real *8 (a-h,o-z)
  real *8 :: z(2), pols(*), dersx(*), dersy(*), w(*)
  real *8 :: dersxx(*), dersxy(*), dersyy(*)
  data c0/.7598356856515925473311877506545453353968d0/
  data c1/1.861209718204199197882437493966468175502d0/
  data c2/1.861209718204199197882437493966468175502d0/
  !
  ! c
  ! c     This subroutine evaluates at the user-supplied point z
  ! c     a collection of polynomials (of x,y) orthogonal on the
  ! c     standard triangle, together with their derivatives with 
  ! c     respect to x and y.
  ! c
  ! c     The "standard" triangle is the triangle with the
  ! c     vertices
  ! c     
  ! c     (0,2/sqrt(3)), (-1,-1/sqrt(3), (1,-1/sqrt(3)).       (1)
  ! c     
  ! c     The polynomials evaluated by this subroutine are all 
  ! c     orthogonal polynomials up to order mmax, arranged in the 
  ! c     increasing order. 
  ! c
  ! c     This subroutine is based on the Koornwinder's representation
  ! c     of the orthogonal polynomials on the right triangle
  ! c
  ! c     (-1,-1), (-1,1), (1,-1)                              (2)
  ! c
  ! c     given by
  ! c
  ! c     K_mn(x,y) = P_m ((2*x+1+y)/(1-y)) ((1-y)/2)^m P_n^{2m+1,0} (y)
  ! c
  ! c     where P_m are the Legendre polynomials or order m
  ! c     and P_n^{2m+1,0} are the Jacobi polynomials of order n with
  ! c     the parameters alpha=2*m+1 and beta=0.
  ! c
  ! c                   Input parameters:
  ! c
  ! c  mmax - the maximum order to which the polynomials are to be evaluated;
  ! c  z - the location in R^2 where the polynomials are to be evaluated;
  ! c       normally, expected to be inside (including boundary) the 
  ! c       standard triangle (1) above.
  ! c  w - the work array. must be at least mmax+1+(mmax+1)**2 long
  ! c     
  ! c                   Output parameters:
  ! c
  ! c  pols - the orthogonal polynomials evaluated at the point z 
  ! c       ( (mmax+1)*(mmax+2)/2 of them things)
  ! c  dersx - the derivatives with respect to x of the polynomials 
  ! c       returned in array pols
  ! c  dersy - the derivatives with respect to y of the polynomials 
  ! c       returned in array pols
  !
  !

  !
  ! if the order is small, then do directly and return
  !

  do i = 1,(mmax+1)*(mmax+2)/2
    dersxx(i) = 0
    dersxy(i) = 0
    dersyy(i) = 0
  end do
  
  
  if( mmax .eq. 0 ) then
    pols(1)=c0
    dersx(1)=0
    dersy(1)=0
    return
  endif

  if( mmax .eq. 1 ) then
    pols(1)=c0
    pols(2)=z(1)*c1
    pols(3)=z(2)*c2
    dersx(1)=0
    dersx(2)=c1
    dersx(3)=0
    dersy(1)=0
    dersy(2)=0
    dersy(3)=c2
    return
  endif
  
  if( mmax .eq. 2 ) then
    pols(1)=c0
    pols(2)=z(1)*c1
    pols(3)=z(2)*c2
    dersx(1)=0
    dersx(2)=c1
    dersx(3)=0
    dersy(1)=0
    dersy(2)=0
    dersy(3)=c2
    x=z(1)
    y=z(2)

    pols(4)=-0.653962434751713743412735927577371370559589d0 &
        +1.13269616323141501386119433259387761572248d0*y &
        -0.490471826063785307559551945683028527919690d0*y**2 &
        +4.41424643457406776803596751114725675127721d0*x**2
    pols(5)=0.759835685651592547331187750654545335396773d0* &
        (2d0+8.66025403784438646763723170752936183471405d0*y)*x
    pols(6)=-0.731152229418051367121788278776110586200037d0 &
        -1.01311424753545672977491700087272711386236d0*y &
        +4.38691337650830820273072967265666351720022d0*y*y

    dersx(4)=4.41424643457406776803596751114725675127721d0*x*2
    dersx(5)=0.759835685651592547331187750654545335396773d0* &
        (2d0+8.66025403784438646763723170752936183471405d0*y)
    dersx(6)=0

    dersy(4)= &
        +1.13269616323141501386119433259387761572248d0 &
        -0.490471826063785307559551945683028527919690d0*y*2
    dersy(5)=0.759835685651592547331187750654545335396773d0* &
        (8.66025403784438646763723170752936183471405d0)*x
    dersy(6)= &
        -1.01311424753545672977491700087272711386236d0 &
        +4.38691337650830820273072967265666351720022d0*y*2

    return
  endif

  iw1=1
  iw2=iw1+mmax+1
  iw3=iw2+(mmax+1)**2
  iw4=iw3+mmax+1
  iw5=iw4+mmax+1
  iw6=iw5+(mmax+1)**2
  iw7=iw6+(mmax+1)**2
  call ortho2eva60(mmax, z, pols, dersx, dersy, dersxx, &
      dersxy, dersyy, w(iw1), w(iw2), w(iw3), w(iw4), &
      w(iw5), w(iw6))
  return
end subroutine ortho2eva6


subroutine ortho2eva60(mmax, z, pols, dersx, dersy, dersxx, &
    dersxy, dersyy, pvals, jvals, pdersx, pdersy, jders, f6)
  implicit real *8 (a-h,o-z)
  real *8 :: z(2), pols(*), dersx(*), dersy(*)
  real *8 :: dersxx(*), dersxy(*), dersyy(*)
  real *8 :: pvals(*), jvals(*), pdersx(*), pdersy(*), jders(*), f6(*)
  real *8 :: pdersxx(10000), pdersxy(10000), pdersyy(10000)
  real *8 :: jders2(10000)
  
  real *8, parameter :: zero = 0

  real *8, &
      parameter :: sqrt2 = 1.414213562373095048801688724209698078570d0
  real *8, &
      parameter :: sqrt3 = 1.732050807568877293527446341505872366943d0
  real *8, &
      parameter :: r11 = -.3333333333333333333333333333333333333333d0
  real *8, &
      parameter :: r12 = -.5773502691896257645091487805019574556476d0
  real *8, &
      parameter :: r21 = -.3333333333333333333333333333333333333333d0
  real *8, &
      parameter :: r22 = 1.154700538379251529018297561003914911295d0

  !
  ! evaluate the orthonormal polynomials on the triangle, 
  ! together with their derivatives with respect to x and y.
  !
  ! ... map the standard triangle to the right
  ! triangle with the vertices (-1,-1), (1,-1), (-1,1)
  !
  x=z(1)
  y=z(2)
  a=r11+r12*y+x
  b=r21+r22*y

  !
  ! first evaluate the scaled Legendre polynomials and the Jacobi
  ! polynomials via recurrence, and then combine into the Koornwinder
  ! polynomials
  !
  !
  par1=(2*a+1+b)/2
  par2=(1-b)/2
  !!!!call klegeypols3(par1, par2, mmax, pvals, pdersx, pdersy)
  call klegeypols6(par1, par2, mmax, pvals, pdersx, pdersy, &
      pdersxx, pdersxy, pdersyy)

  do m = 0,mmax
    par1 = 2*m+1
    ind = 1+m*(mmax+1)
    mtot = mmax - m
    !!!!call kjacopols2(b, par1, zero, mtot, jvals(ind), &
    !!!!    jders(ind))
    call kjacopols3(b, par1, zero, mtot, jvals(ind), &
        jders(ind), jders2(ind))
  end do

  !
  ! now combine
  !
  kk=0
  do m=0,mmax
    do n=0,m
      kk=kk+1

      !
      ! ... evaluate the polynomial (m-n, n), and their derivatives
      !     with respect to x,y
      !
      ip = m-n+1
      ij = n+1+(m-n)*(mmax+1)
      
      pols(kk)=pvals(ip)*jvals(ij)

      dersx(kk) = pdersx(ip)*jvals(ij)

      dersy(kk)= pvals(ip)*jders(ij)*r22 + &
          pdersx(ip)*jvals(ij)*(r12+r22/2) + &
          pdersy(ip)*jvals(ij)*(-r22/2) 

      !
      !  ... and normalize it
      !
      scale = sqrt((1.0d0+(m-n)+n)*(1.0d0+(m-n)+(m-n))/sqrt3) 
      !!!! scale=1

      pols(kk)=pols(kk)*scale
      dersx(kk)=dersx(kk)*scale
      dersy(kk)=dersy(kk)*scale

      dersxx(kk) = dersxx(kk)*scale
      dersxy(kk) = dersxy(kk)*scale
      dersyy(kk) = dersyy(kk)*scale
      
    end do
  end do
  
  return
end subroutine ortho2eva60




      
subroutine klegeypols6(x, y, n, pols, dersx, dersy, dersxx, &
    dersxy, dersyy)
  implicit real *8 (a-h,o-z)
  real *8 :: pols(*), dersx(*), dersy(*), dersxx(*)
  real *8 :: dersxy(*), dersyy(*)

  !
  ! Evaluate a sequence of scaled Legendre polynomials P_n(x/y) y^n,
  ! with the parameter y \in [0..1], together with their first and
  ! second derivatives derivatives with respect to the parameters x
  ! and y.
  !
  ! NOTE: Take note that n+1 values are returned, NOT n
  
  ! Input:
  !   x, y - evaluation location
  !   n - the maximum order polynomial to return, n+1 functions
  !
  ! Output:
  !   pols - the value of the polynomials
  !   dersx, dersy - first partial derivatives
  !   dersxx, dersxy, dersyy - second partial derivatives
  !
  !
  
  if (n .ge. 0) then
    pols(1) = 1
    dersx(1) = 0
    dersy(1) = 0
    dersxx(1) = 0
    dersxy(1) = 0
    dersyy(1) = 0
    if (n .eq. 0) return
  end if

  if (n .ge. 1) then
    pols(2) = x
    dersx(2) = 1
    dersy(2) = 0
    dersxx(2) = 0
    dersxy(2) = 0
    dersyy(2) = 0
    if (n .eq. 1) return
  end if

  !
  ! n is greater than 1. conduct recursion
  !
  pkm1=1
  pk=x
  dkm1=0
  dk=1
  ykm1=0
  yk=0

  dxxm1 = 0
  dxx = 0
  dxxp1 = 0
  
  dxym1 = 0
  dxy = 0
  dxyp1 = 0

  dyym1 = 0
  dyy = 0
  dyyp1 = 0
  
  pk=1
  pkp1=x
  dk=0
  dkp1=1
  yk=0
  ykp1=0
  
  do k=1,n-1

    pkm1 = pk
    pk = pkp1
    dkm1 = dk
    dk = dkp1
    ykm1 = yk
    yk = ykp1

    dxxm1 = dxx
    dxx = dxxp1
    
    dxym1 = dxy
    dxy = dxyp1
    
    dyym1 = dyy
    dyy = dyyp1
    
    pkp1=( (2*k+1)*x*pk-k*pkm1*y*y )/(k+1)
    pols(k+2)=pkp1

    dkp1=( (2*k+1)*(x*dk+pk)-k*dkm1*y*y )/(k+1)
    dersx(k+2)=dkp1

    ykp1=( (2*k+1)*(x*yk)-k*(pkm1*2*y+ykm1*y*y) )/(k+1)
    dersy(k+2)=ykp1

    dxxp1 = ( (2*k+1)*(2*dk + x*dxx) - k*dxxm1*y*y )/(k+1)
    dersxx(k+2) = dxxp1

    dxyp1 = ((2*k+1)*(x*dxy + yk) - k*(y*y*dxym1 + 2*y*dkm1))/(k+1)
    dersxy(k+2) = dxyp1
    
    dyyp1 = ((2*k+1)*x*dyy - k*(2*pkm1 + 4*y*ykm1 + dyym1*y*y))/(k+1)
    dersyy(k+2) = dyyp1
    
  end do
  
  return
end subroutine klegeypols6





subroutine kjacopols3(x,a,b,n,pols,ders,ders2)
  implicit real *8 (a-h,o-z)
  dimension pols(*), ders(*), ders2(*)
  !c
  !c       evaluates a bunch of Jacobi polynomials (together
  !c       with their derivatives) at the user-provided point x
  !c
  !c       ... if n=0 or n=1 - exit
  !c
  if (n .eq. 0) then
    pols(1)=1
    ders(1)=0
    return
  end if

  if (n .eq. 1) then
    pols(1)=1
    ders(1)=0
    pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
    ders(2)=(1+a/2+b/2)
    return
  end if

  !
  ! n > 1, run recursion
  !
  pols(1)=1
  ders(1)=0
  ders2(1) = 0
  
  pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
  ders(2)=(1+a/2+b/2)
  ders2(2) = 0
  
  pk = pols(1)
  dk = ders(1)
  pkp1 = pols(2)
  dkp1 = ders(2)

  dxx = ders2(1)
  dxxp1 = ders2(2)
  
  do k=2,n

    pkm1=pk
    pk=pkp1

    dkm1=dk
    dk=dkp1

    dxxm1 = dxx
    dxx = dxxp1

    alpha1=(2*k+a+b-1)*(a**2-b**2)
    alpha2=(2*k+a+b-1)*((2*k+a+b-2)*(2*k+a+b))
    
    beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
    gamma=(2*k*(k+a+b)*(2*k+a+b-2))
    
    pkp1 = ((alpha1+alpha2*x)*pk - beta*pkm1)/gamma
    pols(k+1) = pkp1
    
    dkp1 = ((alpha1+alpha2*x)*dk-beta*dkm1+alpha2*pk)/gamma
    ders(k+1) = dkp1

    dxxp1 = ((alpha1+alpha2*x)*dxx + 2*alpha2*dk - beta*dxxm1)/gamma
    ders2(k+1) = dxxp1
  end do

  return
end subroutine kjacopols3


