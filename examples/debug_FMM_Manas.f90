program debug_FMM_Manas
implicit real *8 (a-h,o-z)

integer,parameter :: seed = 86456   !This is for random stuff
integer ( kind = 4 ) ier, iprec, ns, nt, ifcharge,ifdipole,ifpottarg,iffldtarg,count
real ( kind = 8 ) source(3,1000000), dipvec(3,1000000),wts(1000000),trads(1000000),stdev(1000000)
real (kind = 8 ) stdev_grad(3,1000000),targets(3,1000000)
complex ( kind = 8 ) sigma(1000000), mu(1000000), pottarg(1000000),fldtarg(3,1000000)
real ( kind = 8 ) tfmm

    nt=1000
    ns=1000
!    allocate(sigma(ns))
!    allocate(mu(ns))
!    allocate(source(3,ns))
!    allocate(dipvec(3,ns))
!    allocate(wts(ns))
!    allocate(stdev(nt))
!    allocate(stdev_grad(3,nt))
!    allocate(pottarg(nt))
!    allocate(fldtarg(3,nt))
!    allocate(trads(nt))
!    allocate(targets(3,nt))



    do count=1,ns
        sigma(count)=(1.0d0,0.0d0)*rand(0)+(0.0d0,1.0d0)*rand(0)
        mu(count)=(1.0d0,0.0d0)*rand(0)+(0.0d0,1.0d0)*rand(0)
        source(1,count)=rand(0)-0.5d0
        source(2,count)=rand(0)-0.5d0
        source(3,count)=rand(0)-0.5d0
!        source(:,count)=source(:,count)/norm2(source(:,count))
        dipvec(:,count)=source(:,count)
        wts(count)=1.0d0
    enddo
    do count=1,nt
        stdev(count)=0.0001d0
        stdev_grad(1,count)=0.00d0
        stdev_grad(2,count)=0.00d0
        stdev_grad(3,count)=0.00d0
        trads(count)=stdev(count)*0.0d0
!        targets(1,count)=rand(0)-0.5d0
!        targets(2,count)=rand(0)-0.5d0
!        targets(3,count)=rand(0)-0.5d0
!        targets(:,count)=targets(:,count)/norm2(targets(:,count))
        targets(1,count)=source(1,count)
        targets(2,count)=source(2,count)
        targets(3,count)=source(3,count)
    enddo

    do count=1,nt
write (*,*) norm2(targets(:,count)),norm2(source(:,count))
    enddo

    ifcharge=0
    ifdipole=1
    ifpottarg=1
    iffldtarg=1
    iprec=2
    ier=1
    do count=1,ns
!        write (*,*) source(1,count),source(2,count),source(3,count)
!        write (*,*) dipvec(1,count),dipvec(2,count),dipvec(3,count)
!        write (*,*) wts(count)
    enddo
    write (*,*) 'hola'

        call tfmm3dwrap(ier,iprec,source,dipvec,ns,wts,ifcharge,&
    &sigma,ifdipole,mu,targets,nt,trads,stdev,&
    &stdev_grad,ifpottarg,pottarg,iffldtarg,fldtarg,tfmm)


stop
end program


