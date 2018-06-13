!
!     testing !ode for FMM - tests charges and dipoles against
!     O(N^2) direct method
!
!
        implicit real *8 (a-h,o-z)
        parameter(lw=15000000)
        dimension w(lw)
!
        call test_part(w,lw)
!
        stop
        end
!
!
!
!
!
        subroutine test_part(w,lw)
        implicit real *8 (a-h,o-z)

        real *8 expc(3,1)
        complex *16 texps(10000)
        real *8 scales
        integer ntj

        real *8 rad(1000000)
        real *8 source(3,1000000)
        complex *16 charge(1000000)
        complex *16 dipstr(1000000)
        real *8 wts(1000000)
        real *8 dipvec(3,1000000)

        real *8 trads(1000000)
        real *8 stdev(1000000)
        real *8 stdev_grad(3,1000000)
        complex *16 pot(1000000)
        complex *16 fld(3,1000000)
!
        complex *16 pot2(1000000)
        complex *16 fld2(3,1000000)
!
        dimension targets(3,2000000)
        complex *16 pottarg(2000000)
        complex *16 fldtarg(3,2000000)

        integer flags(2000000)

!       Tree variables
        integer idivflag,ndiv,isep,nlmax,nlevels,nboxes,mhung,ltree

        integer, allocatable :: itree(:)
        real *8, allocatable :: boxsize(:)
        real *8, allocatable :: treecenters(:,:)

        integer ipointer(26)

!
!cc        complex *16 ptemp,ftemp(3)
        real *8 ptemp,ftemp(3)
!
!cc        parameter(lw=120 000 000)
        dimension w(1)
!
        complex *16 ima
        complex *16 zk
        data ima/(0.0d0,1.0d0)/
!
!
        done=1
        pi=4*atan(done)
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
        do i=1,lw
        w(i)=0
        enddo
!
        idist=3
!
        if( idist == 1 ) then
!
!       ... construct randomly located charge distribution on a unit cube
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
!       ... construct charge distribution on a curve in R^3
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
!       ... construct randomly located charge distribution on a unit sphere
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
!       ... construct the grid of charges
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
        charge(i)=1
        enddo
!
!
        do i=1,nsource
           wts(i)=hkrand(0)
!c           wts(i)=1
           dipstr(i)=1
           dipvec(1,i)=1
           dipvec(2,i)=2
           dipvec(3,i)=3
        enddo
!
        call tfmm3dwrap(ier,iprec,source,dipvec,nsource,&
          &wts,ifcharge,charge,ifdipole,dipstr,targets,&
          &ntarget,trads,stdev,stdev_grad,&
          &ifpottarg,pottarg,iffldtarg,fldtarg,tfmm)

!
        t2=second()
!
!


        call prinf('ier=*',ier,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        call prin2('after fmm, time (sec)=*',t2-t1,1)
!cc        call prin2('after fmm, speed (points/sec)=*',nsource/(t2-t1),1)
        call prin2('after fmm, speed (points+targets/sec)=*',&
          &(nsource+ntarget)/(t2-t1),1)

!
!cc        m=nsource
        m=min(ntarget,100)
!
!
        ifprint=0
        if (ifprint == 1) then
        call prin2('source=*',source,3*nsource)
        endif

        ifprint=0
        if (ifprint == 1) then
        call prin2('after fmm, pot=*',pot,2*m)
        if( iffld==1 ) call prin2('after fmm, fld=*',fld,3*2*m)
        endif
!
!
        if( ntarget == 0 ) stop
!
        do i=1,ntarget
        if (ifpottarg==1) pot2(i)=0
        if (iffldtarg==1) then
           fld2(1,i)=0
           fld2(2,i)=0
           fld2(3,i)=0
        endif

        enddo
!
        t1=second()
!
        call prin2('pot2=*',pot2,12)
        call prin2('pottarg=*',pottarg,12)
        if (ifprint == 1) then
        if (ifpottarg==1) call prin2('directly, pottarg=*',pot2,2*m)
        if( iffldtarg==1 ) call prin2('directly, fldtarg=*',fld2,3*2*m)
        endif
!
        call prin2('directly, estimated time (sec)=*',(t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',m/(t2-t1),1)
!
        if (ifpottarg == 1) then
             do i=1,m
                write(27,*) dreal(pottarg(i)),dreal(pot2(i)),&
               &dabs(dreal(pottarg(i)-pot2(i)))/dabs(dreal(pot2(i)))
             enddo
        call h3derror(pottarg,pot2,m,aerr,rerr)
!cc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in target potential=*',rerr,1)
        endif
!
        if (iffldtarg == 1) then
        call h3derror(fldtarg,fld2,3*m,aerr,rerr)
!cc         call prin2('absolute L2 error in field=*',aerr,1)
        call prin2('relative L2 error in target field=*',rerr,1)
        endif

        return
        end
!
!
!
!
!
        subroutine h3dmperr(mpole1,mpole2,nterms,d)
        implicit real *8 (a-h,o-z)
!
        complex *16 mpole1(0:nterms,-nterms:nterms)
        complex *16 mpole2(0:nterms,-nterms:nterms)
!
        d=0
!
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole1(n,m)-mpole2(n,m))**2
        enddo
        enddo
!
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
!
        return
        end
!
!
!
!
!
        subroutine h3dmpnorm(mpole,nterms,d)
        implicit real *8 (a-h,o-z)
!
        complex *16 mpole(0:nterms,-nterms:nterms)
!
        d=0
!
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole(n,m))**2
        enddo
        enddo
!
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
!
        return
        end
!
!
!
!
!
        subroutine h3derror(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
!
!       evaluate absolute and relative errors
!
        complex *16 pot1(n),pot2(n)
!
        d=0
        a=0
!
        do i=1,n
        d=d+dabs(dreal(pot1(i))-dreal(pot2(i)))**2
        a=a+cdabs(pot1(i))**2
        enddo
!
        d=d/n
        d=dsqrt(d)
        a=a/n
        a=dsqrt(a)
!
        ae=d
        re=d/a
!
        return
        end
!
!
!

