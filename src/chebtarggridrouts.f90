!-----------------------------------

module prefunrouts

contains

subroutine precompphi(eps,ns,sources,rn,nt,dumtarg,norder, &
       itree,ltree,nlevels,nboxes,nlbox,iptr,treecenters,boxsize, &
       nt2,fcoeffs,fcoeffsx,fcoeffsy,fcoeffsz,sigma_eval)

!       this subroutine tree sorts a collection of sources
!        and targets, and returns the tree  
!        
!
!
!        input:
!        eps - precision requested
!        ns - number of sources
!        sources(3,*) - xyz locations of the sources
!        rn(3,*) - normal orientations at the sources
!        nt - number of targets
!        dumtarg(3,*) - xyz locations of the dummy targets
!        sigma_eval - function handle for evaluating sigma 
!      
!        output
!        itree(ltree) - tree array
!        ltree - length of tree array
!        iptr(7) - pointer for tree array
!        nlevels - number of levels
!        nboxes - number of boxes
!        nlbox - number of leaf boxes
!        treecenters(3,nboxes) - locations of centers of trees
!        boxsize(0:nlevels) - size of boxes at various levels
!        nt2 - number of dummy targets
!        fcoeffs(nlbox*norder**3) - chebyshev expansion coefficients
!                                   of the function to be evaluated
!                                   for the revelant leaf boxes
!        
!        fcoeffsx(nlbox*norder**3) - chebyshev expansion coefficients
!                                   of the x component of the gradient of
!                                   the function to be evaluated
!                                   for the revelant leaf boxes
!
!        fcoeffsy(nlbox*norder**3) - chebyshev expansion coefficients
!                                   of the y component of the gradient of
!                                   the function to be evaluated
!                                   for the revelant leaf boxes
!
!        fcoeffsz(nlbox*norder**3) - chebyshev expansion coefficients
!                                   of the z component of the gradient of
!                                   the function to be evaluated
!                                   for the revelant leaf boxes
!
      
      implicit real *8 (a-h,o-z)

!      calling sequence variables
      real *8, intent(in) :: sources(3,ns)
      real *8, intent(in) :: rn(3,ns)
      real *8, intent(in) :: dumtarg(3,nt)

      integer, allocatable, intent(inout) :: itree(:)
      integer, allocatable, intent(inout) :: iptr(:)
      integer, intent(inout) :: nboxes,nlevels,ltree,nlbox,nt2

      real *8, allocatable, intent(inout) :: treecenters(:,:)
      real *8, allocatable, intent(inout) :: boxsize(:)
      real *8, allocatable, intent(inout) :: fcoeffs(:)
      real *8, allocatable, intent(inout) :: fcoeffsx(:)
      real *8, allocatable, intent(inout) :: fcoeffsy(:)
      real *8, allocatable, intent(inout) :: fcoeffsz(:)

      real *8, allocatable :: radsrc(:)
      integer, allocatable :: ntargbox(:)

      real *8, allocatable :: targets(:,:)
      real *8, allocatable :: fvals(:),fvalsx(:),fvalsy(:),fvalsz(:)
      real *8, allocatable :: fvalstmp(:)
      real *8, allocatable :: fvalsxtmp(:),fvalsytmp(:),fvalsztmp(:)
      

      real *8, allocatable :: sigma(:),sigma_grad(:,:),trads(:)
      real *8 sigmatmp,sigma_gradtmp(3)

!      fmm variables
      integer, allocatable :: itreefmm(:)
      real *8, allocatable :: treecentersfmm(:,:),boxsizefmm(:)
      integer ipointer(32)
      integer, allocatable :: iaddr(:,:),nterms(:)
      real *8, allocatable :: rmlexp(:)
      complex *16, allocatable :: charges(:),dipole(:)
      complex *16, allocatable :: pottarg(:),fldtarg(:,:)

      real *8 xpts(norder),wts(norder)
      real *8, allocatable :: xmatc(:,:)
      real *8, allocatable :: tails(:)
      real *8, allocatable :: tailsx(:)
      real *8, allocatable :: tailsy(:)
      real *8, allocatable :: tailsz(:)
      integer, allocatable :: irefineflag(:)
      integer, allocatable :: ibcompflag(:)
      integer, allocatable :: itcompflag(:)

      integer, allocatable :: itreeout(:)
      integer, allocatable :: ntargboxout(:)
      integer iptrout(7)
      real *8, allocatable :: targetsout(:,:),fcoeffsout(:), &
            boxsizeout(:),treecentersout(:,:),fcoeffsxout(:), &
            fcoeffsyout(:), fcoeffszout(:)

!
!c      temporary variables
!
      real *8 expc(3),radt,wlege(40000),pottmp,gradtmp(3)
      real *8, allocatable :: sourcesort(:,:),dipvecsort(:,:)
      complex *16, allocatable :: dipstrsort(:),chargesort(:)
      integer, allocatable :: ilevel(:)

      external sigma_eval

      nbmax = 0
      nlmax = 200
      mnlist1 = 0
      mnlist2 = 0
      mnlist3 = 0
      mnlist4 = 0
      mnbors = 0
      ltree = 0

      idivflag = 0
      ndiv = 1

      isep = 1

      allocate(radsrc(ns))
      do i=1,ns
         radsrc(i) = 0
      enddo

      ier = 0
      nexpc = 0
      radexp = 0

      call mklraptreemem(ier,sources,ns,radsrc,dumtarg,nt,expc, &
            nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,nlevelsfmm, &
            nboxesfmm,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung , &
            ltreefmm)


       allocate(itreefmm(ltreefmm))
       allocate(boxsizefmm(0:nlevelsfmm))
       allocate(treecentersfmm(3,nboxesfmm))

       call mklraptree(sources,ns,radsrc,dumtarg,nt,expc,nexpc, &
             radexp,idivflag,ndiv,isep,mhung,mnbors,mnlist1, &
             mnlist2,mnlist3,mnlist4,nlevelsfmm,nboxesfmm, &
             treecentersfmm,boxsizefmm,itreefmm,ltreefmm,ipointer)
!
!c      extract relevant bits of the tree needed from this point on
!
        nboxes = nboxesfmm
        nlevels = nlevelsfmm

        allocate(boxsize(0:levels),treecenters(3,nboxes))


        treecenters = treecentersfmm
        boxsize = boxsizefmm

        ltree = 2*(nlevels+1)+12*nboxes
        allocate(iptr(7),itree(ltree),ntargbox(nboxesfmm))


        call extracttreesubinfo(itreefmm,ipointer,nboxes,nlevels, &
                itree,iptr,ntargbox)
!        call prinf('nboxes=*',nboxes,1)
!        call prinf('nlevel=*',nlevels,1)
!        call prinf('itree=*',itree,ltree)
!        call prinf('ntargbox=*',ntargbox,nboxes)
!        call prinf('iptr=*',iptr,7)



!
!c        generate initial set of targets
!
        itype = 1
        call chebexps(itype,norder,xpts,utmp,vtmp,wts)

!
!c        count the number of leaf boxes

        nlbox = 0
        do i=1,nboxes
           nchild = itree(iptr(4)+i-1)
           if(nchild.eq.0.and.ntargbox(i).gt.0) nlbox = nlbox+1
        enddo
        
        itarg = 0
        ictr = 1
        norder3 = norder*norder*norder
        
        call prinf('norder3=*',norder3,1)
        call prinf('nlbox=*',nlbox,1)
        nt2 = norder3*nlbox
        call prinf('nt2=*',nt2,1)


        allocate(targets(3,nt2),fvals(nt2),fcoeffs(nt2))
        allocate(fvalsx(nt2),fvalsy(nt2),fvalsz(nt2))
        allocate(fcoeffsx(nt2),fcoeffsy(nt2),fcoeffsz(nt2))
        allocate(fvalstmp(norder3),fvalsxtmp(norder3))
        allocate(fvalsytmp(norder3),fvalsztmp(norder3))

        do ibox = 1,nboxes 
           nchild = itree(iptr(4)+ibox-1)
           if(nchild.eq.0.and.ntargbox(ibox).gt.0)  then
              itree(iptr(6)+ibox-1) = ictr
              ictr = ictr + norder3
              ilev = itree(iptr(2)+ibox-1)

              do ii =1,norder
                 do jj= 1,norder
                    do kk =1,norder
                       itarg = itarg + 1
                       targets(1,itarg) = treecentersfmm(1,ibox)+ &
                                 xpts(kk)*boxsizefmm(ilev)/2
                       targets(2,itarg) = treecentersfmm(2,ibox)+ &
                                 xpts(jj)*boxsizefmm(ilev)/2
                       targets(3,itarg) = treecentersfmm(3,ibox)+ &
                                 xpts(ii)*boxsizefmm(ilev)/2
                     enddo
                  enddo
              enddo
          endif
       enddo

!
!c       evaluate sigma, trads, stdev initially at the targets
!
!
        allocate(sigma(nt2),sigma_grad(3,nt2),trads(nt2))

        do i=1,nt2
           call sigma_eval(targets(1,i),sigma(i),sigma_grad(1,i))
              trads(i) = 12*sigma(i)
        enddo


        call prinf('nboxes=*',nboxes,1)
        call prinf('nlevels=*',nlevels,1)
        call prinf('ltreefmm=*',ltreefmm,1)

        call prinf('ipointer=*',ipointer,32)

!
!        recompute tree
!

      call mklraptreemem(ier,sources,ns,radsrc,dumtarg,nt,targets, &
            nt2,trads,idivflag,ndiv,isep,nlmax,nbmax,nlevelsfmm, &
            nboxesfmm,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung , &
            ltreefmm)

      call prinf('ltreefmm=*',ltreefmm,1)
      deallocate(itreefmm)
      allocate(itreefmm(ltreefmm))



       call mklraptree(sources,ns,radsrc,dumtarg,nt,targets,nt2, &
             trads,idivflag,ndiv,isep,mhung,mnbors,mnlist1, &
             mnlist2,mnlist3,mnlist4,nlevelsfmm,nboxesfmm, &
             treecentersfmm,boxsizefmm,itreefmm,ltreefmm,ipointer)

        call prinf('nboxes=*',nboxes,1)
        call prinf('nlevels=*',nlevels,1)
        call prinf('ipointer=*',ipointer,32)
!
!!      call the fmm
!
        iprec = 1
        if(eps.le.0.5d-3.and.eps.gt.0.5d-6) iprec = 2
        if(eps.le.0.5d-6.and.eps.gt.0.5d-9) iprec = 3
        if(eps.le.0.5d-9) iprec = 4
        ier = 0

!
!!        set eps fmm


        allocate(iaddr(2,nboxes))
        allocate(nterms(0:nlevels))

      if( iprec .eq. -2 ) epsfmm=.5d-0 
      if( iprec .eq. -1 ) epsfmm=.5d-1
      if( iprec .eq. 0 ) epsfmm=.5d-2
      if( iprec .eq. 1 ) epsfmm=.5d-3
      if( iprec .eq. 2 ) epsfmm=.5d-6
      if( iprec .eq. 3 ) epsfmm=.5d-9
      if( iprec .eq. 4 ) epsfmm=.5d-12
      if( iprec .eq. 5 ) epsfmm=.5d-15
      if( iprec .eq. 6 ) epsfmm=0

      do i=0,nlevels
         if(isep.eq.1) call l3dterms(epsfmm,nterms(i),ier)
         if(isep.eq.2) call l3dterms_far(epsfmm,nterms(i),ier)
      enddo

      call l3dmpalloc_newtree(itree(1),iaddr,nlevels,lmptot,nterms)
      call prinf(' lmptot is *',lmptot,1)


      allocate(rmlexp(lmptot),stat=ier)
      if(ier.ne.0) then
         call prinf('Cannot allocate mpole expansion workspace, &
                   lmptot is *', lmptot,1)
         ier = 16
         return
      endif

      ifcharge = 0
      ifdipole = 1
      ifpottarg = 1
      iffldtarg = 1
      ier = 0

      allocate(charges(ns),dipole(ns),pottarg(nt2),fldtarg(3,nt2))
      allocate(sourcesort(3,ns))
      allocate(chargesort(ns),dipstrsort(ns),dipvecsort(3,ns))

      do i=1,ns
         dipole(i) = 1
         charges(i) = 0 
      enddo

      call prinf('ns=*',ns,1)
      call prin2('dipole=*',dipole,2*ns)
      call prin2('rn=*',rn,3*ns)
      call prin2('sigma=*',sigma,12)
      call prin2('sigma_grad=*',sigma_grad,36)

      call prinf('iprec=*',iprec,1)
      call prinf('ifpottarg=*',ifpottarg,1)
      call prinf('iffldtarg=*',iffldtarg,1)


      call tfmm3dlr(ier,iprec,ns,sources,ifcharge,charges,ifdipole, &
            dipole,rn,nt2,targets,trads,sigma,sigma_grad,itreefmm, &
            ltreefmm,ipointer,isep,ndiv,nlevels,nboxes,boxsizefmm, &
            mnbors,mnlist1,mnlist2,mnlist3,mnlist4,treecentersfmm, &
            iaddr,nterms,lmptot,ifpottarg,pottarg,iffldtarg,fldtarg, &
            rmlexp)

!
!!      store the sorted arrays
!
      call l3dreorder_newtree(ns,sources,ifcharge,charges, &
        itreefmm(ipointer(5)), ifdipole,dipole,rn,sourcesort, &
        chargesort,dipstrsort,dipvecsort)


!
!!       extract the ilevel array
!

      allocate(ilevel(nboxesfmm))
      do ilev=0,nlevelsfmm
          do ibox = itreefmm(2*ilev+1),itreefmm(2*ilev+2)
              ilevel(ibox) = ilev
          enddo
      enddo

!
!!     extract fvals and the gradient values
!
      do i=1,nt2
         fvals(i) = real(pottarg(i))
         fvalsx(i) = real(fldtarg(1,i))
         fvalsy(i) = real(fldtarg(2,i))
         fvalsz(i) = real(fldtarg(3,i))
      enddo

!      call prin2('fvals=*',fvals,512)

     
!
!c        get the interpolation matrix
!
      allocate(xmatc(norder3,norder3))
      call getxmatc(norder,norder3,xmatc)
      do i=1,nboxes
         nchild = itree(iptr(4)+i-1)
         if(nchild.eq.0.and.ntargbox(i).gt.0) then
            ii = itree(iptr(6)+i-1)
            call matvec(norder3,norder3,xmatc,fvals(ii),fcoeffs(ii))
            call matvec(norder3,norder3,xmatc,fvalsx(ii),fcoeffsx(ii))
            call matvec(norder3,norder3,xmatc,fvalsy(ii),fcoeffsy(ii))
            call matvec(norder3,norder3,xmatc,fvalsz(ii),fcoeffsz(ii))
        endif
      enddo

!      call prin2('fcoeffs=*',fcoeffs,512)

      iflag = 1

      maxit = 100

      nlege = 100
      lw7 = 40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)

      do iter = 1,maxit
        call prinf('iter=*',iter,1)
        rmaxerr = 0
        iflag = 0 
        allocate(tails(nboxes),irefineflag(nboxes))
        allocate(tailsx(nboxes),tailsy(nboxes),tailsz(nboxes))
        do i=1,nboxes
           nchild = itree(iptr(4)+i-1)
           tails(i) = 0
           tailsx(i) = 0
           tailsy(i) = 0
           tailsz(i) = 0
           irefineflag(i) = 0
           if(nchild.eq.0.and.ntargbox(i).gt.0) then
              ii = itree(iptr(6)+i-1)

!              call prinf('ibox=*',i,1)
!              read *, itmp

              call comptail(norder,norder3,fcoeffs(ii),tails(i))
              call comptail(norder,norder3,fcoeffsx(ii),tailsx(i))
              call comptail(norder,norder3,fcoeffsy(ii),tailsy(i))
              call comptail(norder,norder3,fcoeffsz(ii),tailsz(i))
!              call prin2('fcoeffs=*',fcoeffs(ii),norder3)
!              call prin2('fcoeffsx=*',fcoeffsx(ii),norder3)
!              call prin2('fcoeffsy=*',fcoeffsy(ii),norder3)
!              call prin2('fcoeffsz=*',fcoeffsz(ii),norder3)

              if(tails(i).gt.eps.or.tailsx(i).gt.eps.or. &
                   tailsy(i).gt.eps.or.tailsz(i).gt.eps) then
                 irefineflag(i) = 1
                 iflag = 1
              endif
              if(tails(i).gt.rmaxerr) rmaxerr = tails(i)
              if(tailsx(i).gt.rmaxerr) rmaxerr = tailsx(i)
              if(tailsy(i).gt.rmaxerr) rmaxerr = tailsy(i)
              if(tailsz(i).gt.rmaxerr) rmaxerr = tailsz(i)
           endif
         enddo
         call prin2('rmaxerr=*',rmaxerr,1)

         if(iflag.eq.0) goto 2000
      
!         call prin2('tails=*',tails,nboxes)
!         call prinf('iflag=*',iflag,1)


!      do i=1,nboxes
!         irefineflag(i) = 0
!      enddo
!  
!      irefineflag(6) = 1
!      call prinf('irefineflag=*',irefineflag,nboxes)
!         read *, ii

!         call prinf('nboxes=*',nboxes,1)
!         call prinf('itree(iptr(6))=*',itree(iptr(6)),nboxes)
         call reorganizechebtree_withgrad(itree,ltree,iptr, &
            treecenters, &
            boxsize,nboxes, &
            nlevels,ntargbox,norder,nt2, targets, fcoeffs, fcoeffsx, &
            fcoeffsy,fcoeffsz, &
            irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
            iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
            ntargboxout,nt2out,targetsout,fcoeffsout,fcoeffsxout, &
            fcoeffsyout,fcoeffszout)


            deallocate(itree,treecenters,boxsize,ntargbox,targets, &
               fcoeffs,fcoeffsx,fcoeffsy,fcoeffsz)

            ltree = ltreeout
            nboxes = nboxesout
            nlevels = nlevelsout
            nt2 = nt2out

            call prinf('ltree=*',ltree,1)
            call prinf('nboxes=*',nboxes,1)
            call prinf('nlevels=*',nlevels,1)
            call prinf('nt2=*',nt2,1)
            
            allocate(itree(ltree),treecenters(3,nboxes), &
               boxsize(0:nlevels),targets(3,nt2),fcoeffs(nt2), &
               ntargbox(nboxes),fcoeffsx(nt2),fcoeffsy(nt2), &
               fcoeffsz(nt2))

            iptr = iptrout
            itree = itreeout
            treecenters = treecentersout
            boxsize(0:nlevels) = boxsizeout(0:nlevels)
            targets = targetsout
            fcoeffs = fcoeffsout
            fcoeffsx = fcoeffsxout
            fcoeffsy = fcoeffsyout
            fcoeffsz = fcoeffszout
            ntargbox = ntargboxout
!            call prinf('iptr=*',iptr,7)
!            call prinf('itree(iptr(6))=*',itree(iptr(6)),nboxes)



!            call prinf('ibcompflag=*',ibcompflag,nboxes)
!            call prinf('ifcoeffs=*',itree(iptr(6)),nboxes)

            do i=1,nboxes
               if(ibcompflag(i).eq.1) then
                  ii = itree(iptr(6)+i-1)
!                  call prinf('i=*',i,1)
!                  call prinf('ii=*',ii,1)
                  do j=1,norder3
                     iii = ii+j-1
                     call sigma_eval(targets(1,iii),sigmatmp, &
                        sigma_gradtmp)

!                     call prin2('sigmatmp=*',sigmatmp,1)
                     tradtmp = 12*sigmatmp
                 
                     call findboxtarg_fulltree(targets(1,iii),tradtmp,&
                      iboxtarg,ilev,itreefmm,ipointer,boxsizefmm, &
                      treecentersfmm,nboxesfmm,nlevelsfmm)

!                      call prinf('iboxtarg=*',iboxtarg,1)
!                      call prinf('ilev=*',ilev,1)
!                      call prin2('targets=*',targets(1,iii),3)
!                      call prin2('tradtmp=*',tradtmp,1)
!                      call prin2('dipvecsort=*',dipvecsort,3*ns)

!                     call prinf('mnbors=*',mnbors,1)
!                     call prinf('mnlist1=*',mnlist1,1)
!                     call prinf('nlege=*',nlege,1)

                     if(i.eq.18) then
!                        call prinf('iii=*',iii,1)
!                        call prin2('targets=*',targets(1,iii),3)
!                        call prin2('sigmatmp=*',sigmatmp,1)
!                        call prin2('sigma_gradtmp=*',sigma_gradtmp,3)

                     endif

                     call newtargeval(targets(1,iii),sigmatmp, &
                        sigma_gradtmp, iboxtarg,ilev,nboxesfmm, &
                        nlevelsfmm,boxsizefmm,treecentersfmm,iaddr, &
                        itreefmm,ipointer,mnbors,mnlist1,ilevel, &
                        rmlexp,ns,sourcesort,ifcharge,chargesort,&
                        ifdipole,dipstrsort,dipvecsort, &
                        nterms,nlege,wlege,pottmp,gradtmp)


                      fvalstmp(j) = pottmp
                      fvalsxtmp(j) = gradtmp(1)
                      fvalsytmp(j) = gradtmp(2)
                      fvalsztmp(j) = gradtmp(3)
                  enddo
                  if(i.eq.18) then
!                     call prinf('ibox=*',i,1)
!                     call prin2('fvalstmp=*',fvalstmp,norder3)
                  endif
!                  call prin2('fvalstmp=*',fvalstmp,norder3)
!                  call prin2('fvalstmp=*',fvalsxtmp,norder3)
!                  call prin2('fvalstmp=*',fvalsytmp,norder3)
!                  call prin2('fvalstmp=*',fvalsztmp,norder3)
                  call matvec(norder3,norder3,xmatc,fvalstmp,& 
                     fcoeffs(ii))
                  call matvec(norder3,norder3,xmatc,fvalsxtmp,& 
                     fcoeffsx(ii))
                  call matvec(norder3,norder3,xmatc,fvalsytmp,& 
                     fcoeffsy(ii))
                  call matvec(norder3,norder3,xmatc,fvalsztmp,& 
                     fcoeffsz(ii))

!                  call prin2('fcoeffs=*',fcoeffs(ii),norder3)
               endif
            enddo


            deallocate(tails,irefineflag,itreeout,treecentersout, &
             boxsizeout,ntargboxout,targetsout,fcoeffsout)
            deallocate(fcoeffsxout,fcoeffsyout,fcoeffszout)
            deallocate(tailsx,tailsy,tailsz)
            deallocate(itcompflag,ibcompflag)

      enddo
 2000 continue      


end subroutine precompphi 
!-------------------------------------------     

       subroutine getxmatc(k,k3,xmatc)
       implicit real *8 (a-h,o-z)
       real *8 xmatc(k3,*)
       real *8 uk(k,k),ts(k),wts(k),vk(k,k)

       itype = 2
       call chebexps(itype,k,ts,uk,vk,wts)


       do icw = 1,k
          do icv = 1,k
             do icu = 1,k
                ii = (icw-1)*k*k + (icv-1)*k + icu
                do iw = 1,k
                   do iv = 1,k
                      do iu = 1,k
                       jj = (iw-1)*k*k + (iv-1)*k + iu
                       xmatc(ii,jj) = uk(icw,iw)*uk(icv,iv)*uk(icu,iu)
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo

       return
end subroutine getxmatc
!---------------------------------------------------------
      subroutine findboxtarg(targtest,iboxtarg,ilevel,itree,iptr, &
          treecenters,nboxes,nlevels)
!
!!      the only difference between findboxtarg_fulltree and 
!       and findboxtarg is the pointer array and where the children
!      information is stored in the tree
!
      implicit real *8 (a-h,o-z)
      integer iptr(*),itree(*)
      real *8 treecenters(3,*),targtest(3)

      iboxtarg = 1
      ilevel = 0
      do ilev=0,nlevels-1
!
!c       compute relartive corordinates for teh target with
!        respect to the box
!
         nchild = itree(iptr(4)+iboxtarg-1)

         if(nchild.eq.0) goto 3000
         xtmp = targtest(1) - treecenters(1,iboxtarg)
         ytmp = targtest(2) - treecenters(2,iboxtarg)
         ztmp = targtest(3) - treecenters(3,iboxtarg)
         if(xtmp.le.0.and.ytmp.le.0.and.ztmp.le.0) ichild=1
         if(xtmp.gt.0.and.ytmp.le.0.and.ztmp.le.0) ichild=2
         if(xtmp.le.0.and.ytmp.gt.0.and.ztmp.le.0) ichild=3
         if(xtmp.gt.0.and.ytmp.gt.0.and.ztmp.le.0) ichild=4
         if(xtmp.le.0.and.ytmp.le.0.and.ztmp.gt.0) ichild=5
         if(xtmp.gt.0.and.ytmp.le.0.and.ztmp.gt.0) ichild=6
         if(xtmp.le.0.and.ytmp.gt.0.and.ztmp.gt.0) ichild=7
         if(xtmp.gt.0.and.ytmp.gt.0.and.ztmp.gt.0) ichild=8

         iboxtarg = itree(iptr(5)+(iboxtarg-1)*8+ichild-1)
         ilevel = ilev + 1

      enddo
 3000 continue      

      return
end subroutine findboxtarg      
!------------------------------------------------------------      
      subroutine findboxtarg_fulltree(targtest,trad,iboxtarg,ilevel, &
          itree,iptr,boxsize,treecenters,nboxes,nlevels)
!
!!      the only difference between findboxtarg_fulltree and 
!       and findboxtarg is the pointer array and where the children
!      information is stored in the tree
          
      implicit real *8 (a-h,o-z)
      integer iptr(*),itree(*)
      real *8 treecenters(3,*),targtest(3),boxsize(0:nlevels)

!
!c       compute relartive corordinates for teh target with
!        respect to the box
!
      iboxtarg = 1
      ilevel = 0
      do ilev=0,nlevels-1

!          exit loop if target is hung at this level
!
         if(trad.gt.boxsize(ilev)) goto 3000
         nchild = itree(iptr(3)+iboxtarg-1)

         if(nchild.eq.0) goto 3000
         xtmp = targtest(1) - treecenters(1,iboxtarg)
         ytmp = targtest(2) - treecenters(2,iboxtarg)
         ztmp = targtest(3) - treecenters(3,iboxtarg)
         if(xtmp.le.0.and.ytmp.le.0.and.ztmp.le.0) ichild=1
         if(xtmp.gt.0.and.ytmp.le.0.and.ztmp.le.0) ichild=2
         if(xtmp.le.0.and.ytmp.gt.0.and.ztmp.le.0) ichild=3
         if(xtmp.gt.0.and.ytmp.gt.0.and.ztmp.le.0) ichild=4
         if(xtmp.le.0.and.ytmp.le.0.and.ztmp.gt.0) ichild=5
         if(xtmp.gt.0.and.ytmp.le.0.and.ztmp.gt.0) ichild=6
         if(xtmp.le.0.and.ytmp.gt.0.and.ztmp.gt.0) ichild=7
         if(xtmp.gt.0.and.ytmp.gt.0.and.ztmp.gt.0) ichild=8

         iboxtarg = itree(iptr(4)+(iboxtarg-1)*8+ichild-1)
         ilevel = ilev + 1

      enddo
 3000 continue      

      return
end subroutine findboxtarg_fulltree      
!-----------------------------------------------------

      subroutine cheb3deval(targtest,rscale,center,k,fcoeffs,f)
      implicit real *8 (a-h,o-z)
      real *8 xpols(k+10),ypols(k+10),zpols(k+10)
      real *8 fcoeffs(*),targtest(3),center(3)

      x = (targtest(1) - center(1))/rscale*2
      y = (targtest(2) - center(2))/rscale*2
      z = (targtest(3) - center(3))/rscale*2

      f = 0
      call chebpols(x,k-1,xpols)
      call chebpols(y,k-1,ypols)
      call chebpols(z,k-1,zpols)

      do iz = 1,k
         do iy = 1,k
            do ix = 1,k
               i = (iz-1)*k*k  + (iy-1)*k + ix
               f = f + fcoeffs(i)*xpols(ix)*ypols(iy)*zpols(iz)
            enddo
         enddo
      enddo


      return
end subroutine cheb3deval      
!--------------------------------------------      
subroutine extracttreesubinfo(itreetmp,ipointer,nboxes,nlevels, &
     itree,iptr,ntargbox)

!
!c      itree(iptr(1)) - laddr
!       itree(iptr(2)) - ilevel
!       itree(iptr(3)) - iparent
!       itree(iptr(4)) - nchild
!       itree(iptr(5)) - ichild
!       itree(iptr(6)) - pointer to fcoeffs array

     implicit real *8 (a-h,o-z)
     integer itree(*),itreetmp(*),iptr(*),ipointer(*),ntargbox(*)

     
        iptr(1) = 1
        iptr(2) = iptr(1)+2*(nlevels+1)
        iptr(3) = iptr(2)+nboxes
        iptr(4) = iptr(3)+nboxes
        iptr(5) = iptr(4)+nboxes
        iptr(6) = iptr(5)+8*nboxes
        iptr(7) = iptr(6)+nboxes


        do i=1,2*nlevels+2
           itree(i) = itreetmp(i)
        enddo

        do ilev=0,nlevels
           do ibox = itree(2*ilev+1),itree(2*ilev+2)
              itree(iptr(2)+ibox-1) = ilev
           enddo
        enddo

        do ibox=1,nboxes
           itree(iptr(3)+ibox-1) = itreetmp(ipointer(2)+ibox-1)
           itree(iptr(4)+ibox-1) = itreetmp(ipointer(3)+ibox-1)
           do j=1,8
              itree(iptr(5)+ 8*(ibox-1)+j-1) =  &
                  itreetmp(ipointer(4)+8*(ibox-1)+j-1)
           enddo
           itree(iptr(6)+ibox-1) = -1

           istart = itreetmp(ipointer(12)+ibox-1)
           iend = itreetmp(ipointer(13)+ibox-1)

           ntargbox(ibox) = iend-istart+1
        enddo
end subroutine extracttreesubinfo
!-----------------------------------------------
subroutine comptail(norder,norder3,fcoeffs,tail)
    implicit real *8 (a-h,o-z)
    real *8 fcoeffs(*),tail

    tail = 0
    rnorder = norder + 0.0d0

    do ii = 1,norder
      do jj = 1,norder
         do kk = 1,norder

            i = ii+jj+kk-3
            ri = sqrt((ii-1.0d0)**2 + (jj-1.0d0)**2 + (kk-1.0d0)**2)
            if(ri.ge.rnorder) then
!            if(ii.eq.norder.or.jj.eq.norder.or.kk.eq.norder) then
!            if(i.ge.norder) then
                iind = (ii-1)*norder*norder+(jj-1)*norder+kk
                tail = tail + abs(fcoeffs(iind))**2
            endif
         enddo
      enddo
    enddo

    tail = sqrt(tail)
end subroutine

subroutine matvec(m,n,a,x,y)
implicit real *8 (a-h,o-z)
   real *8 a(m,n),x(n),y(m)
   
   do i=1,m
     y(i) = 0
     do j=1,n
        y(i) = y(i) + a(i,j)*x(j)
     enddo
   enddo


end subroutine
!-----------------------------------
      subroutine reorganizechebtree_withgrad(itree,ltree,iptr, &
         treecenters, &
         boxsize,nboxes, &
         nlevels,ntargbox,norder,nt2, targets,fcoeffs, fcoeffsx, &
         fcoeffsy,fcoeffsz, &
         irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
         iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
         ntargboxout,nt2out,targetsout,fcoeffsout,fcoeffsxout, &
         fcoeffsyout,fcoeffszout)

       implicit real *8 (a-h,o-z)
       integer, intent(inout) :: nboxes,nlevels,nt2
       integer, dimension(:), intent(in) ::  itree
       integer iptr(7)
       integer, dimension(:), intent(in) :: ntargbox
       real *8, dimension(:,:), intent(in) :: treecenters
       real *8, intent(in) :: boxsize(0:nlevels)
       real *8, dimension(:,:), intent(in) :: targets
       real *8, dimension(:), intent(in) :: fcoeffs
       real *8, dimension(:), intent(in) :: fcoeffsx
       real *8, dimension(:), intent(in) :: fcoeffsy
       real *8, dimension(:), intent(in) :: fcoeffsz
       integer norder
       integer, intent(out), allocatable :: ibcompflag(:)
       integer, intent(out), allocatable :: itcompflag(:)
       integer, intent(in), dimension(:) :: irefineflag


       integer, intent(out), allocatable :: itreeout(:)
       integer, intent(out) :: ltreeout,iptrout(7),nboxesout, &
              nlevelsout,nt2out
       integer, intent(out), allocatable :: ntargboxout(:)
       real *8, intent(out), allocatable :: targetsout(:,:), &
              fcoeffsout(:), boxsizeout(:),treecentersout(:,:)
       real *8, intent(out), allocatable :: fcoeffsxout(:), &
                  fcoeffsyout(:),fcoeffszout(:)


       integer, allocatable :: tntargbox(:),tibcompflag(:)
       integer, allocatable :: laddrtail(:,:),tilev(:),tladdr(:,:)
       integer, allocatable :: tiparent(:),tnchild(:),tichild(:,:)
       integer tiptr(7)
       real *8, allocatable :: ttargets(:,:),tboxsize(:)
       real *8, allocatable :: tcenters(:,:)
       integer curbox
       integer, allocatable :: iboxtocurbox(:)
       real *8 xpts(norder),wts(norder),utmp,vtmp

!
!!       count number of additional boxes to be formed
!
      
      nextra = 0
      nlevelsout = nlevels
      do i=1,nboxes
         if(irefineflag(i).eq.1) nextra = nextra + 8
         nchild = itree(iptr(4)+i-1)
         ilev = itree(iptr(2)+i-1)
         if(ilev.eq.nlevels.and.irefineflag(i).eq.1) &
             nlevelsout = nlevels+1
      enddo


      nboxesout = nboxes + nextra

      ltreeout = 2*(nlevelsout+1)+12*(nboxesout)



      norder3 = norder*norder*norder

      allocate(tntargbox(nboxesout),tibcompflag(nboxesout))
      allocate(laddrtail(2,0:nlevels+1),tladdr(2,0:nlevels+1))
      allocate(tilev(nboxesout),tnchild(nboxesout),tichild(8,nboxesout))
      allocate(tiparent(nboxesout))
      allocate(tcenters(3,nboxesout))

!
!!      copy everything into temporary arrays
!

      do i=0,nlevels+1
         laddrtail(1,i) = 0
         laddrtail(2,i) = -1
      enddo

      do i=0,nlevels
         tladdr(1,i) = itree(2*i+1)
         tladdr(2,i) = itree(2*i+2)
      enddo

      tladdr(1,nlevels+1) = nboxes+1
      tladdr(2,nlevels+1) = nboxes

!      call prinf('itree=*',itree,ltree)


      do i=1,nboxes
         tilev(i) = itree(iptr(2)+i-1)
         tiparent(i) = itree(iptr(3)+i-1)
         tnchild(i) = itree(iptr(4)+i-1)
         do j=1,8
            tichild(j,i) = itree(iptr(5)+8*(i-1)+j-1)
         enddo
         tntargbox(i) = ntargbox(i)
         tibcompflag(i) = 0
         if(tnchild(i).eq.0.and.ntargbox(i).gt.0) tibcompflag(i) = 2
         do j=1,3
            tcenters(j,i) = treecenters(j,i)
         enddo
      enddo

      nboxesout = nboxes
      do ilev = 0,nlevels
         laddrtail(1,ilev+1) = nboxesout+1
         do ibox = itree(2*ilev+1),itree(2*ilev+2)
            if(irefineflag(ibox).eq.1) then
               tnchild(ibox) = 8
               tibcompflag(ibox) = -1
               do i=1,8
                  ii = 2
                  jj = 2
                  if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
                  if(i.lt.5) jj = 1
                  nboxesout = nboxesout + 1
                  tichild(i,ibox) = nboxesout
                  tibcompflag(nboxesout) = 1
                  tcenters(1,nboxesout)=treecenters(1,ibox)+(-1)**i*&
                                        boxsize(ilev)/8.0d0
                  tcenters(2,nboxesout)=treecenters(2,ibox)+(-1)**ii*&
                                        boxsize(ilev)/8.0d0
                  tcenters(3,nboxesout)=treecenters(3,ibox)+(-1)**jj*&
                                        boxsize(ilev)/8.0d0
                  tnchild(nboxesout) = 0
                  do j=1,8
                     tichild(j,nboxesout) = -1
                  enddo
                  tiparent(nboxesout) = ibox
                  tntargbox(nboxesout) = 1
                  tilev(nboxesout) = ilev+1
                  tibcompflag(nboxesout) = 1
               enddo
            endif
         enddo
         laddrtail(2,ilev+1) = nboxesout
      enddo

!      call prinf('nboxesout=*',nboxesout,1)
!      call prinf('ltreeout=*',ltreeout,1)


      allocate(itreeout(ltreeout),treecentersout(3,nboxesout))
      allocate(ntargboxout(nboxesout),iboxtocurbox(nboxesout))
      allocate(boxsizeout(0:nlevelsout),ibcompflag(nboxesout))

!
!!      reset the iptr array
!
        iptrout(1) = 1
        iptrout(2) = iptrout(1)+2*(nlevelsout+1)
        iptrout(3) = iptrout(2)+nboxesout
        iptrout(4) = iptrout(3)+nboxesout
        iptrout(5) = iptrout(4)+nboxesout
        iptrout(6) = iptrout(5)+8*nboxesout
        iptrout(7) = iptrout(6)+nboxesout

!        call prinf('nlevels=*',nlevels,1)
!
!        call prinf('nboxesout=*',nboxesout,1)
!        call prinf('nlevelsout=*',nlevelsout,1)
!        call prinf('iptrout=*',iptrout,7)



!      call prinf('laddrtail=*',laddrtail,2*(nlevels+2))
!      call prinf('tladdr=*',tladdr,2*(nlevels+1))

      curbox = 1
      do ilev=0,nlevelsout
         itreeout(2*ilev+1) = curbox
         do ibox = tladdr(1,ilev),tladdr(2,ilev)
!            call prinf('ibox=*',ibox,1)
!            call prinf('curbox=*',curbox,1)
            itreeout(iptrout(2)+curbox-1) = tilev(ibox)
            itreeout(iptrout(4)+curbox-1) = tnchild(ibox)
            itreeout(iptrout(6)+curbox-1) = -1
            ibcompflag(curbox) = tibcompflag(ibox)
            treecentersout(1,curbox) = tcenters(1,ibox)
            treecentersout(2,curbox) = tcenters(2,ibox)
            treecentersout(3,curbox) = tcenters(3,ibox)
            ntargboxout(curbox) = tntargbox(ibox)
            iboxtocurbox(ibox) = curbox
            curbox = curbox + 1
         enddo

         do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
!            call prinf('ibox=*',ibox,1)
!            call prinf('curbox=*',curbox,1)
            itreeout(iptrout(2)+curbox-1) = tilev(ibox)
            itreeout(iptrout(4)+curbox-1) = tnchild(ibox)
            itreeout(iptrout(6)+curbox-1) = -1
            ibcompflag(curbox) = tibcompflag(ibox)
            treecentersout(1,curbox) = tcenters(1,ibox)
            treecentersout(2,curbox) = tcenters(2,ibox)
            treecentersout(3,curbox) = tcenters(3,ibox)
            ntargboxout(curbox) = tntargbox(ibox)
            iboxtocurbox(ibox) = curbox
            curbox = curbox + 1
         enddo
         itreeout(2*ilev+2) = curbox-1
      enddo
!      call prinf('itreeout=*',itreeout(iptrout(4)),nboxesout)
!      call prinf('tichild=*',tichild,8*nboxesout)
!      call prinf('iboxtocurbox=*',iboxtocurbox,nboxesout)


!
!!         handle the parent child part of the tree using the
!          mapping iboxtocurbox

       ictr = 1
       do ibox=1,nboxesout
          if(tiparent(ibox).eq.-1) &
            itreeout(iptrout(3)+iboxtocurbox(ibox)-1) = -1
          if(tiparent(ibox).gt.0) &
               itreeout(iptrout(3)+iboxtocurbox(ibox)-1)=&
               iboxtocurbox(tiparent(ibox))
          do i=1,8
             if(tichild(i,ibox).eq.-1) &
                itreeout(iptrout(5)+(iboxtocurbox(ibox)-1)*8+i-1)=-1
             if(tichild(i,ibox).gt.0) &
                itreeout(iptrout(5)+(iboxtocurbox(ibox)-1)*8+i-1)= &
                     iboxtocurbox(tichild(i,ibox))
          enddo
          
          nchild = itreeout(iptrout(4)+ibox-1)
          if(nchild.eq.0.and.ntargboxout(ibox).gt.0) then
             itreeout(iptrout(6)+ibox-1) = ictr
             ictr = ictr + norder3
          endif
       enddo

       nt2out = ictr-1
       allocate(fcoeffsout(nt2out),targetsout(3,nt2out), &
          itcompflag(nt2out),fcoeffsxout(nt2out),fcoeffsyout(nt2out), &
          fcoeffszout(nt2out))

!       call prinf('iptrout=*',iptrout,7)
!
!       call prinf('itreeout=*',itreeout,ltreeout)
!       call prinf('nchild=*',itreeout(iptrout(4)),nboxesout)


!        call prinf('ticompflag=*',tibcompflag,nboxesout)
!
!!        fix targets and fcoeffs
!
       itype = 1
       call chebexps(itype,norder,xpts,utmp,vtmp,wts)
       
       do i=0,nlevels
          boxsizeout(i) = boxsize(i)
       enddo

       if(nlevelsout.eq.nlevels+1) &
          boxsizeout(nlevelsout) = boxsizeout(nlevels)/2.0d0

!       call prin2('boxsizeout=*',boxsizeout,nlevelsout+1)   
       do ibox=1,nboxesout
          if(tibcompflag(ibox).eq.2) then
              curbox = iboxtocurbox(ibox)
!              call prinf('ibox=*',ibox,1)
!              call prinf('curbox=*',curbox,1)
              ii = itree(iptr(6)+ibox-1)
              jj = itreeout(iptrout(6)+curbox-1)
!              call prinf('ii=*',ii,1)
!              call prinf('jj=*',jj,1)
              do i=1,norder3
                 fcoeffsout(jj+i-1) = fcoeffs(ii+i-1)
                 fcoeffsxout(jj+i-1) = fcoeffsx(ii+i-1)
                 fcoeffsyout(jj+i-1) = fcoeffsy(ii+i-1)
                 fcoeffszout(jj+i-1) = fcoeffsz(ii+i-1)
                 targetsout(1,jj+i-1) = targets(1,ii+i-1)
                 targetsout(2,jj+i-1) = targets(2,ii+i-1)
                 targetsout(3,jj+i-1) = targets(3,ii+i-1)
                 itcompflag(jj+i-1) = 0
              enddo
          endif

          if(tibcompflag(ibox).eq.1) then
             curbox = iboxtocurbox(ibox)
             itarg = itreeout(iptrout(6)+curbox-1)
             ilev = itreeout(iptrout(2)+curbox-1)
             do ii=1,norder
                do jj=1,norder
                   do kk=1,norder
                       targetsout(1,itarg)=treecentersout(1,curbox)+&
                                 xpts(kk)*boxsizeout(ilev)/2
                       targetsout(2,itarg)=treecentersout(2,curbox)+&
                                 xpts(jj)*boxsizeout(ilev)/2
                       targetsout(3,itarg)=treecentersout(3,curbox)+&
                                 xpts(ii)*boxsizeout(ilev)/2
                       itcompflag(itarg) = 1
                       itarg = itarg + 1
                   enddo
                enddo
             enddo
          endif
       enddo
end subroutine


subroutine newtargeval(targets,stdev,stdev_grad,ibox,ilev, &
             nboxes,&
             nlevels,boxsize,centers, iaddr,itree, &
             ipointer,mnbors,mnlist1,ilevel,rmlexp,ns,sourcesort, &
             ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
             nterms,nlege,wlege,pot,grad)
     implicit real *8 (a-h,o-z)
     real *8 targets(*),boxsize(0:nlevels),rmlexp(*),wlege(*)
     real *8 sourcesort(3,*),dipvecsort(3,*),centers(3,*)
     complex *16 chargesort(*), dipstrsort(*)
     integer iaddr(2,*),itree(*),ipointer(32),nterms(0:nlevels)
     integer ilevel(*)
     real *8 pot,grad(3),stdev_grad(3)
     complex *16 pottmp,gradtmp(3)

!
!!        evaluate local expansion 
!

     npts = 1
     done = 1
     pi = atan(done)*4

     ifpottarg = 1
     ifgradtarg = 1
     npts = 1

!     call prinf('ilev=*',ilev,1)

     pottmp = 0
     gradtmp(1) = 0
     gradtmp(2) = 0
     gradtmp(3) = 0

     call l3dtaevalall_trunc(boxsize(ilev),centers(1,ibox), &
       rmlexp(iaddr(2,ibox)),nterms(ilev),nterms(ilev),targets, &
       npts,ifpottarg,pottmp,ifgradtarg,gradtmp,wlege,nlege,ier)

     pottmp = -pottmp/(4*pi)
     gradtmp(1) = -gradtmp(1)/(4*pi)
     gradtmp(2) = -gradtmp(2)/(4*pi)
     gradtmp(3) = -gradtmp(3)/(4*pi)




     itstart = 1
     itend = 1

     nnbors = itree(ipointer(18)+ibox-1)
     do i=1,nnbors
        jbox = itree(ipointer(19)+mnbors*(ibox-1)+i-1)
        jstart = itree(ipointer(10)+jbox-1)
        jend = itree(ipointer(11)+jbox-1)

        call tfmm3dparttarg_direct_newtree(jstart,jend,itstart,itend, &
           sourcesort,ifcharge,chargesort,ifdipole,dipstrsort, &
           dipvecsort, targets, stdev, stdev_grad, ifpottarg, &
           pottmp, ifgradtarg, gradtmp)


     enddo


     nlist1 = itree(ipointer(20)+ibox-1)

     do i=1,nlist1
        jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
        if(ilevel(ibox).gt.ilevel(jbox)) then
           jstart = itree(ipointer(10)+jbox-1)
           jend = itree(ipointer(11)+jbox-1)
           
           call tfmm3dparttarg_direct_newtree(jstart, &
             jend, itstart,itend, &
             sourcesort,ifcharge,chargesort,ifdipole,dipstrsort, &
             dipvecsort, targets, stdev, stdev_grad, ifpottarg, &
             pottmp, ifgradtarg, gradtmp)
        endif
     enddo

     pot = real(pottmp)
     grad(1) = real(gradtmp(1))
     grad(2) = real(gradtmp(2))
     grad(3) = real(gradtmp(3))

   
end subroutine newtargeval

end module
