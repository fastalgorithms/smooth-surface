!
!     testing !ode for FMM - tests charges and dipoles against
!     O(N^2) direct method
!
program debug_FMM_Manas
implicit real *8 (a-h,o-z)

integer,parameter :: seed = 86456   !This is for random stuff
integer ( kind = 4 ) ier, iprec, ns, nt, ifcharge,ifdipole,ifpottarg,iffldtarg,count
real ( kind = 8 ) source(3,1000000), dipvec(3,1000000),wts(1000000),trads(1000000),stdev(1000000)
real (kind = 8 ) stdev_grad(3,1000000),targets(3,1000000)
complex ( kind = 8 ) sigma(1000000), mu(1000000), pottarg(1000000),fldtarg(3,1000000)
real ( kind = 8 ) tfmm

    !    real *8 source(3,1000000)
!        complex *16 sigma(1000000)
!        complex *16 mu(1000000)
    !    real *8 wts(1000000)
    !    real *8 dipvec(3,1000000)

     !   real *8 trads(1000000)
      !  real *8 stdev(1000000)
!        real *8 stdev_grad(3,1000000)
!!!        complex *16 pot(1000000)
!!!        complex *16 fld(3,1000000)
!
!!!        complex *16 pot2(1000000)
!!!        complex *16 fld2(3,1000000)
!
!        dimension targets(3,2000000)
!        complex *16 pottarg(2000000)
!        complex *16 fldtarg(3,2000000)

!!!        integer flags(2000000)

!       Tree variables
!        integer idivflag,ndiv,isep,nlmax,nlevels,nboxes,mhung,ltree

!        integer, allocatable :: itree(:)
!        real *8, allocatable :: boxsize(:)
!        real *8, allocatable :: treecenters(:,:)

!        integer ipointer(26)

!
!cc        complex *16 ptemp,ftemp(3)
!!!        real *8 ptemp,ftemp(3)
!
!cc        parameter(lw=120 000 000)

!
!!!        complex *16 ima
!!!        complex *16 zk
!!!        data ima/(0.0d0,1.0d0)/
!
!
!!!        done=1
!!!        pi=4*atan(done)
!
!
!       SET ALL PARAMETERS
!
        call prini(6,13)
!
        print *, 'ENTER nsource'
        read *, nsource
!
!
        call prinf('nsource=*',nsource,1)
!

!
        idist=3
!
        if( idist == 1 ) then
!
!       ... construct randomly located sigma distribution on a unit cube
!
        do i=1,nsource
        source(1,i)=hkrand(0)
        source(2,i)=hkrand(0)
        source(3,i)=hkrand(0)
        source(1,i)=source(1,i)-0.5
        source(2,i)=source(2,i)-0.5
        source(3,i)=source(3,i)-0.5
        enddo
!
        endif
!
        if( idist == 2 ) then
!
!       ... construct sigma distribution on a curve in R^3
!
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        source(1,i)=sin(1.1*a)
        source(2,i)=cos(2.1*a)
        source(3,i)=cos(3.1*a)
        enddo
!
        endif
!
        if( idist == 3 ) then
!
!       ... construct randomly located sigma distribution on a unit sphere
!
        done=1
        pi=4*atan(done)
!
        d=hkrand(0)
        do i=1,nsource
!
!
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        source(1,i)=.5d0*cos(phi)*sin(theta)
        source(2,i)=.5d0*sin(phi)*sin(theta)
        source(3,i)=.5d0*cos(theta)
!
        d=hkrand(0)
        d=hkrand(0)
        enddo
!
        endif
!
        if( idist == 4 ) then
!
!       ... construct the grid of sigmas
!
        ngrid=nsource-1
        kk=0
        do i=1,ngrid+1
        do j=1,ngrid+1
        do k=1,ngrid+1
        kk=kk+1
        source(1,kk)=(i-1.0)/ngrid
        source(2,kk)=(j-1.0)/ngrid
        source(3,kk)=(k-1.0)/ngrid
        source(1,kk)=source(1,kk)+hkrand(0)*.0001
        source(2,kk)=source(2,kk)+hkrand(0)*.0001
        source(3,kk)=source(3,kk)+hkrand(0)*.0001
        enddo
        enddo
        enddo
        nsource=kk
!
        call prinf('after grid build, nsource=*',nsource,1)
!
        endif
!
!        do i=1,nsource
!        source(1,i)=source(1,i)*10
!        source(2,i)=source(2,i)*10
!        source(3,i)=source(3,i)*10
!        enddo


!       ... set up the targets
!
        if( idist == 1 .or. idist == 4 ) then
        do i=1,nsource
        targets(1,i)=source(1,i) + 0.1d0
        targets(2,i)=source(2,i)
        targets(3,i)=source(3,i)
        enddo
        ntarget=nsource
        do i=1,nsource
        targets(1,i+nsource)=source(1,i)
        targets(2,i+nsource)=source(2,i) - 0.2d0
        targets(3,i+nsource)=source(3,i)
        enddo
        ntarget=nsource*2
        endif
!
!
        if( idist == 2 ) then
!
!       ... construct target distribution on a curve in R^3
!
        ntarget=nsource*4
        do i=1,ntarget
        a=2*pi*dble(i)/ntarget
        targets(1,i)=sin(1.1*a)/2
        targets(2,i)=cos(2.1*a)/2
        targets(3,i)=cos(3.1*a)/2
        enddo
        endif
!
        if( idist == 3 ) then
!
!       ... construct target distribution on a unit sphere
!       highly oversampled
!
        ntarget=nsource*1
        do i=1,ntarget
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        targets(1,i)=.51d0*cos(phi)*sin(theta)
        targets(2,i)=.51d0*sin(phi)*sin(theta)
        targets(3,i)=.51d0*cos(theta)
        d=hkrand(0)
        d=hkrand(0)
        enddo
        endif


        do i=1,ntarget
        stdev(i) = hkrand(0)*.01
        stdev_grad(1,i) = hkrand(0)
        stdev_grad(2,i) = hkrand(0)
        stdev_grad(3,i) = hkrand(0)
        trads(i) = 6*stdev(i)
!cc        trads(i) = stdev(i)*0.01
        enddo

        do i=1,10
!cc        trads(i) = 1.0d-1
        enddo

!
        call prinf('ntarget=*',ntarget,1)
!
!
        iprec=2
!
!
!
        ifcharge=0
        ifdipole=1
!
        ifpottarg=1
        iffldtarg=1
!
        do i=1,nsource
        sigma(i)=1
        enddo
!
!
        do i=1,nsource
           wts(i)=hkrand(0)
!c           wts(i)=1
           mu(i)=1
           dipvec(1,i)=1
           dipvec(2,i)=2
           dipvec(3,i)=3
        enddo
!







    write (*,*) ifcharge,ifdipole,ifpottarg,iffldtarg,iprec,ier
    ifcharge=0
    ifdipole=1
    ifpottarg=1
    iffldtarg=1
    iprec=2
    ier=0



        call tfmm3dwrap(ier,iprec,source,dipvec,nsource,&
          &wts,ifcharge,sigma,ifdipole,mu,targets,&
          &ntarget,trads,stdev,stdev_grad,&
          &ifpottarg,pottarg,iffldtarg,fldtarg,tfmm)

!

        stop
        end program

!
