program main
    use prefunrouts
    implicit real *8 (a-h,o-z)

    real *8, allocatable :: sources(:,:),rn(:,:)
    real *8, allocatable :: dumtarg(:,:)

    integer, allocatable :: itree(:),iptr(:)
    real *8, allocatable :: treecenters(:,:),boxsize(:)
    real *8, allocatable :: fcoeffs(:),fcoeffsx(:),fcoeffsy(:)
    real *8, allocatable :: fcoeffsz(:)
    real *8, allocatable :: xmat(:,:)
    real *8 targtest(3),fexgrad(3)
    integer ltree

    external sigma_eval1

    call prini(6,13)
    call prin2('Enter n*',n,0)
!    read *, n

!    ns = 100
!    nt = 100
!    allocate(sources(3,ns),dumtarg(3,nt))
!
!    do i=1,ns
!       sources(1,i) = hkrand(0)
!       sources(2,i) = hkrand(0)
!       sources(3,i) = hkrand(0)
!    enddo
!
!    do i=1,nt
!       dumtarg(1,i) = hkrand(0)
!       dumtarg(2,i) = hkrand(0)
!       dumtarg(3,i) = hkrand(0)
!    enddo

    ns = 9
    nt = 9

    allocate(sources(3,ns),dumtarg(3,nt),rn(3,ns))
    
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

      do i=1,ns
         rn(1,i) = hkrand(0)
         rn(2,i) = hkrand(0)
         rn(3,i) = hkrand(0)

         rn(1,i) = 1
         rn(2,i) = 2
         rn(3,i) = 3

         rrr = sqrt(rn(1,i)**2 + rn(2,i)**2 + rn(3,i)**2)
          
         rn(1,i) = rn(1,i)/rrr
         rn(2,i) = rn(2,i)/rrr
         rn(3,i) = rn(3,i)/rrr
      enddo

      do i=1,9
         do j=1,3
            dumtarg(j,i) = sources(j,i)
         enddo
      enddo


      eps = 1.0d-6
      norder = 6

      call precompphi(eps,ns,sources,rn,nt,dumtarg,norder, &
         itree,ltree, &
         nlevels,nboxes,iptr,treecenters,boxsize,nt2,fcoeffs, &
         fcoeffsx, fcoeffsy, fcoeffsz, sigma_eval1)
       call prinf('nboxes=*',nboxes,1)
       call prinf('nt2=*',nt2,1)

       targtest(1) = hkrand(0)
       targtest(2) = hkrand(0)
       targtest(3) = hkrand(0)
       
!
!c       generate random targets in a particular leaf ndoes
!

      call prinf('iptr(6) arr=*',itree(iptr(6)),nboxes)
!      in testing code make sure itree(iptr(6)+ibox-1) .ne. -1     
!        i.e. make sure there is a dummy target in the box you
!        are testing and that it is the leaf box containing the
!        dummy target
!
      ibox = 56
      ilev = itree(iptr(2) + ibox-1)
      call prin2('treecenters=*',treecenters(1,ibox),3)
      call prinf('ilev=*',ilev,1)

      xtest =treecenters(1,ibox) + (hkrand(0)-0.5d0)*boxsize(ilev)
      ytest =treecenters(2,ibox) + (hkrand(0)-0.5d0)*boxsize(ilev)
      ztest =treecenters(3,ibox) + (hkrand(0)-0.5d0)*boxsize(ilev)

      targtest(1) = xtest
      targtest(2) = ytest
      targtest(3) = ztest

      call findboxtarg(targtest,iboxtarg,ilevel,itree,iptr, &
            treecenters,nboxes,nlevels)
!
!c      evaluate expansions to get
!

      istart = itree(iptr(6)+iboxtarg-1)
      call prinf('istart=*',istart,1)
!      call prinf('itree(iptr(6)=*',itree(iptr(6)),nboxes)
!      call prinf('nlevels=*',nlevels,1)
      call cheb3deval(targtest,boxsize(ilevel),treecenters(1,iboxtarg), &
          norder,fcoeffs(istart),f)
!
!!       evaluate function manually
!
      call dirfuneval(ns,sources,rn,targtest,fex,fexgrad) 

      err = abs(f-fex)

      call prin2('f=*',f,1)
      call prin2('fex=*',fex,1)
      call prin2('err=*',err,1)

end program main


!---------------------------------------------

subroutine sigma_eval1(xyz,sigma,sigma_grad)
   implicit real *8 (a-h,o-z)
   real *8 sigma_grad(3),xyz(3)

   rr2 = xyz(1)**2 + xyz(2)**2 + xyz(3)**2

   sigma = 1+(2*rr2)**2

   sigma_grad(1) = 16*rr2*xyz(1)
   sigma_grad(2) = 16*rr2*xyz(2)
   sigma_grad(3) = 16*rr2*xyz(3)

   rscale = 10
   sigma = (0.1 + xyz(1) + xyz(2) + xyz(3))/rscale
   sigma_grad(1) = xyz(1)/rscale
   sigma_grad(2) = xyz(2)/rscale
   sigma_grad(3) = xyz(3)/rscale


end subroutine sigma_eval1

!------------------------------------------------
  
subroutine dirfuneval(ns,sources,rn,targ,fex,fexgrad)
    implicit real *8 (a-h,o-z)
    real *8 sources(3,*),targ(3),fex,fexgrad(3),rn(3,*)
    real *8 std,std_grad(3),gradtmp(3)

    call sigma_eval1(targ,std,std_grad)

    dipstr = 1
    ifgradtarg = 1
    fex = 0
    fexgrad(1) = 0
    fexgrad(2) = 0
    fexgrad(3) = 0

    do i=1,ns
       call tpotfld3d_dp(ifgradtarg,sources(1,i),dipstr,rn(1,i),&
          targ,std,std_grad,pottmp,gradtmp)
       fex = fex + pottmp
       fexgrad(1) = fexgrad(1) + gradtmp(1)
       fexgrad(2) = fexgrad(2) + gradtmp(2)
       fexgrad(3) = fexgrad(3) + gradtmp(3)
    enddo
end subroutine
