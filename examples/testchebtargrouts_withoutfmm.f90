program main
    use prefunrouts_withoutfmm
    implicit real *8 (a-h,o-z)

    real *8, allocatable :: sources(:,:)
    real *8, allocatable :: dumtarg(:,:)

    integer, allocatable :: itree(:),iptr(:)
    real *8, allocatable :: treecenters(:,:),boxsize(:)
    real *8, allocatable :: fcoeffs(:)
    real *8, allocatable :: xmat(:,:)
    real *8 targtest(3)
    integer ltree

    external fun1

    call prini(6,13)
    call prin2('Enter n*',n,0)
    read *, n

    ns = 9
    nt = 9

    allocate(sources(3,ns),dumtarg(3,nt))
    
      sources(1,1) = 0
      sources(2,1) = 0
      sources(3,1) = 0

      sources(1,2) = 1
      sources(2,2) = 0
      sources(3,2) = 0


      sources(1,3) = 1
      sources(2,3) = 1
      sources(3,3) = 0

      sources(1,4) = 0
      sources(2,4) = 1
      sources(3,4) = 0

      sources(1,5) = 0
      sources(2,5) = 0
      sources(3,5) = 1

      sources(1,6) = 1
      sources(2,6) = 0
      sources(3,6) = 1


      sources(1,7) = 1
      sources(2,7) = 1
      sources(3,7) = 1

      sources(1,8) = 0
      sources(2,8) = 1
      sources(3,8) = 1

      sources(1,9) = 1 - 0.124d0
      sources(2,9) = 1
      sources(3,9) = 1

      do i=1,9
         do j=1,3
            dumtarg(j,i) = sources(j,i)
         enddo
      enddo

      eps = 1.0d-9
      norder = 8

      call precompphi_fun1(eps,ns,sources,nt,dumtarg,norder, &
         itree,ltree, &
         nlevels,nboxes,nlbox,iptr,treecenters,boxsize,nt2,fcoeffs, &
         fun1)
       call prinf('nboxes=*',nboxes,1)
       call prinf('nt2=*',nt2,1)
       
!
!c       generate random targets in a particular leaf ndoes
!
      ibox = 6
      xtest =treecenters(1,ibox) + (hkrand(0)-0.5d0)*boxsize(1)
      ytest =treecenters(2,ibox) + (hkrand(0)-0.5d0)*boxsize(1)
      ztest =treecenters(3,ibox) + (hkrand(0)-0.5d0)*boxsize(1)

      targtest(1) = xtest
      targtest(2) = ytest
      targtest(3) = ztest

      call prin2('targtest=*',targtest,3)

      call findboxtarg(targtest,iboxtarg,ilevel,itree,iptr, &
            treecenters,nboxes,nlevels)





!
!c      evaluate expansions to get
!

      istart = itree(iptr(6)+iboxtarg-1)
      call prinf('istart=*',istart,1)
      call cheb3deval(targtest,boxsize(ilevel),treecenters(1,iboxtarg), &
          norder,fcoeffs(istart),f)

      call fun1(targtest(1),targtest(2),targtest(3),fex)

      call prin2('f=*',f,1)
      call prin2('fex=*',fex,1)

      err = abs(fex-f)

      call prin2('error=*',err,1)


end program main


subroutine fun1(x,y,z,f)
implicit real *8 (a-h,o-z)

    f = sin(x*y*z)


end subroutine fun1

