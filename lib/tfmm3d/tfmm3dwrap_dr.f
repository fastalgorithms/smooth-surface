c
c     testing code for FMM - tests charges and dipoles against
c     O(N^2) direct method
c
c
        implicit real *8 (a-h,o-z)
        parameter(lw=15 000 000)
        dimension w(lw)
c
        call test_part(w,lw)
c
        stop
        end
c
c
c
c
c0,0, 0.25
        subroutine test_part(w,lw)
        implicit real *8 (a-h,o-z)

        real *8 expc(3,1)
        complex *16 texps(10000)
        real *8 scales
        integer ntj

        real *8 rad(1 000 000)
        real *8 source(3,1 000 000)
        complex *16 charge(1 000 000)
        complex *16 dipstr(1 000 000)
        real *8 wts(1 000 000)
        real *8 dipvec(3,1 000 000)

        real *8 trads(1 000 000)
        real *8 stdev(1 000 000)
        real *8 stdev_grad(3,1 000 000)
        complex *16 pot(1 000 000)
        complex *16 fld(3,1 000 000)
c
        complex *16 pot2(1 000 000)
        complex *16 fld2(3,1 000 000)
c
        dimension targets(3,2 000 000)
        complex *16 pottarg(2 000 000)
        complex *16 fldtarg(3,2 000 000)

        integer flags(2 000 000)

c       Tree variables
        integer idivflag,ndiv,isep,nlmax,nlevels,nboxes,mhung,ltree

        integer, allocatable :: itree(:)
        real *8, allocatable :: boxsize(:)
        real *8, allocatable :: treecenters(:,:)
 
        integer ipointer(26)

c
ccc        complex *16 ptemp,ftemp(3)
        real *8 ptemp,ftemp(3)
c
ccc        parameter(lw=120 000 000)
        dimension w(1)
c
        complex *16 ima
        complex *16 zk
        data ima/(0.0d0,1.0d0)/
c
c
        done=1
        pi=4*atan(done)
c
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        print *, 'ENTER nsource'
        read *, nsource
c
c
        call prinf('nsource=*',nsource,1)
c
        do i=1,lw
        w(i)=0
        enddo
c
        idist=3
c
        if( idist .eq. 1 ) then
c
c       ... construct randomly located charge distribution on a unit cube
c
        do i=1,nsource
        source(1,i)=hkrand(0)
        source(2,i)=hkrand(0)
        source(3,i)=hkrand(0)
        source(1,i)=source(1,i)-0.5
        source(2,i)=source(2,i)-0.5
        source(3,i)=source(3,i)-0.5
        enddo
c
        endif
c
        if( idist .eq. 2 ) then
c
c       ... construct charge distribution on a curve in R^3
c
        do i=1,nsource
        a=2*pi*dble(i)/nsource
        source(1,i)=sin(1.1*a)
        source(2,i)=cos(2.1*a)
        source(3,i)=cos(3.1*a)
        enddo
c
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct randomly located charge distribution on a unit sphere
c
        done=1
        pi=4*atan(done)
c
        d=hkrand(0)
        do i=1,nsource
c
c
        theta=hkrand(0)*pi
        phi=hkrand(0)*2*pi
        source(1,i)=.5d0*cos(phi)*sin(theta)
        source(2,i)=.5d0*sin(phi)*sin(theta)
        source(3,i)=.5d0*cos(theta)
c
        d=hkrand(0)
        d=hkrand(0)
        enddo
c
        endif
c
        if( idist .eq. 4 ) then
c
c       ... construct the grid of charges
c
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
c
        call prinf('after grid build, nsource=*',nsource,1)
c
        endif
c
c        do i=1,nsource
c        source(1,i)=source(1,i)*10
c        source(2,i)=source(2,i)*10
c        source(3,i)=source(3,i)*10
c        enddo


c       ... set up the targets
c
        if( idist .eq. 1 .or. idist .eq. 4 ) then
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
c
c
        if( idist .eq. 2 ) then
c
c       ... construct target distribution on a curve in R^3
c
        ntarget=nsource*4
        do i=1,ntarget
        a=2*pi*dble(i)/ntarget
        targets(1,i)=sin(1.1*a)/2
        targets(2,i)=cos(2.1*a)/2
        targets(3,i)=cos(3.1*a)/2
        enddo
        endif
c
        if( idist .eq. 3 ) then
c
c       ... construct target distribution on a unit sphere
c       highly oversampled
c
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
ccc        trads(i) = stdev(i)*0.01
        enddo

        do i=1,10
ccc        trads(i) = 1.0d-1
        enddo

c
        call prinf('ntarget=*',ntarget,1)
c
c
        iprec=2
c
c
c
        ifcharge=0
        ifdipole=1
c
        ifpottarg=1
        iffldtarg=1
c
        do i=1,nsource
        charge(i)=1
        enddo
c
c
        do i=1,nsource
           wts(i)=hkrand(0)
cc           wts(i)=1
           dipstr(i)=1
           dipvec(1,i)=1
           dipvec(2,i)=2
           dipvec(3,i)=3
        enddo
c
        call tfmm3dwrap(ier,iprec,source,dipvec,nsource,
     $     wts,ifcharge,charge,ifdipole,dipstr,targets,
     $     ntarget,trads,stdev,stdev_grad,
     $     ifpottarg,pottarg,iffldtarg,fldtarg,tfmm)

c
        t2=second()
C$        t2=omp_get_wtime()
c
c


        call prinf('ier=*',ier,1)
        call prinf('nsource=*',nsource,1)
        call prinf('ntarget=*',ntarget,1)
        call prin2('after fmm, time (sec)=*',t2-t1,1)
ccc        call prin2('after fmm, speed (points/sec)=*',nsource/(t2-t1),1)
        call prin2('after fmm, speed (points+targets/sec)=*',
     $     (nsource+ntarget)/(t2-t1),1)

c
ccc        m=nsource
        m=min(ntarget,100)
c
c
        ifprint=0
        if (ifprint .eq. 1) then
        call prin2('source=*',source,3*nsource)
        endif

        ifprint=0
        if (ifprint .eq. 1) then
        call prin2('after fmm, pot=*',pot,2*m)
        if( iffld.eq.1 ) call prin2('after fmm, fld=*',fld,3*2*m)
        endif
c
c
        if( ntarget .eq. 0 ) stop
c
        do i=1,ntarget
        if (ifpottarg .eq. 1) pot2(i)=0
        if (iffldtarg .eq. 1) then
           fld2(1,i)=0
           fld2(2,i)=0
           fld2(3,i)=0
        endif

        enddo
c
        t1=second()
C$        t1=omp_get_wtime()
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i,j,ptemp,ftemp)
cccC$OMP$SCHEDULE(DYNAMIC)
cccC$OMP$NUM_THREADS(4)
        do 8160 j=1,m
        do 8150 i=1,nsource
        if( ifcharge .eq. 1 ) then
        call lpotfld3d(iffldtarg,
     1   source(1,i),charge(i),
     $     targets(1,j),
     1     ptemp,ftemp)
        if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
        if (iffldtarg .eq. 1) then
           fld2(1,j)=fld2(1,j)+ftemp(1)
           fld2(2,j)=fld2(2,j)+ftemp(2)
           fld2(3,j)=fld2(3,j)+ftemp(3)
        endif
        endif
        if (ifdipole .eq. 1) then
           call tpotfld3d_dp(iffldtarg,source(1,i),
     $     wts(i),dipvec(1,i),
     $     targets(1,j),stdev(j),stdev_grad(1,j),ptemp,ftemp)
           if (ifpottarg .eq. 1) pot2(j)=pot2(j)+ptemp
           if (iffldtarg .eq. 1) then
              fld2(1,j)=fld2(1,j)+ftemp(1)
              fld2(2,j)=fld2(2,j)+ftemp(2)
              fld2(3,j)=fld2(3,j)+ftemp(3)
           endif
        endif
c
 8150   continue
 8160   continue
C$OMP END PARALLEL DO
c
        t2=second()
C$        t2=omp_get_wtime()
c
        call prin2('pot2=*',pot2,12)
        call prin2('pottarg=*',pottarg,12)
        if (ifprint .eq. 1) then
        if (ifpottarg .eq. 1)
     $     call prin2('directly, pottarg=*',pot2,2*m)
        if( iffldtarg.eq.1 )
     $     call prin2('directly, fldtarg=*',fld2,3*2*m)
        endif
c
        call prin2('directly, estimated time (sec)=*',
     $     (t2-t1)*dble(ntarget)/dble(m),1)
        call prin2('directly, estimated speed (targets/sec)=*',
     $     m/(t2-t1),1)
c
        if (ifpottarg .eq. 1) then
             do i=1,m
                write(27,*) dreal(pottarg(i)),dreal(pot2(i)),
     1          dabs(dreal(pottarg(i)-pot2(i)))/dabs(dreal(pot2(i)))
             enddo
        call h3derror(pottarg,pot2,m,aerr,rerr)
ccc        call prin2('absolute L2 error in potential=*',aerr,1)
        call prin2('relative L2 error in target potential=*',rerr,1)
        endif
c
        if (iffldtarg .eq. 1) then
        call h3derror(fldtarg,fld2,3*m,aerr,rerr)
ccc         call prin2('absolute L2 error in field=*',aerr,1)
        call prin2('relative L2 error in target field=*',rerr,1)
        endif

        return
        end
c
c
c
c
c
        subroutine h3dmperr(mpole1,mpole2,nterms,d)
        implicit real *8 (a-h,o-z)
c
        complex *16 mpole1(0:nterms,-nterms:nterms)
        complex *16 mpole2(0:nterms,-nterms:nterms)
c
        d=0
c
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole1(n,m)-mpole2(n,m))**2
        enddo
        enddo
c
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
c
        return
        end
c
c
c
c
c
        subroutine h3dmpnorm(mpole,nterms,d)
        implicit real *8 (a-h,o-z)
c
        complex *16 mpole(0:nterms,-nterms:nterms)
c
        d=0
c
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole(n,m))**2
        enddo
        enddo
c
        d=d/(nterms+1)
        d=d/(2*nterms+1)
        d=sqrt(d)
c
        return
        end
c
c
c
c
c
        subroutine h3derror(pot1,pot2,n,ae,re)
        implicit real *8 (a-h,o-z)
c
c       evaluate absolute and relative errors
c
        complex *16 pot1(n),pot2(n)
c
        d=0
        a=0
c
        do i=1,n
        d=d+dabs(dreal(pot1(i))-dreal(pot2(i)))**2
        a=a+cdabs(pot1(i))**2
        enddo
c
        d=d/n
        d=dsqrt(d)
        a=a/n
        a=dsqrt(a)
c
        ae=d
        re=d/a
c
        return
        end
c
c
c
