!-----------------------------------
subroutine precompphi(eps,ns,sources,nt,dumtarg,norder, &
       itree,ltree,nlevels,nboxes,nlbox,iptr,treecenters,boxsize, &
       nt2,fcoeffs,fun1)

!       this subroutine tree sorts a collection of sources
!        and targets, and returns the tree  
!        
!
!
!        input:
!        eps - precision requested
!        ns - number of sources
!        sources(3,*) - xyz locations of the sources
!        nt - number of targets
!        dumtarg(3,*) - xyz locations of the dummy targets
!        fun1 - function handle of the function to be evaluated
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
!
      
      implicit real *8 (a-h,o-z)
interface      
      subroutine reorganizechebtree(itree,ltree,iptr,treecenters, &
         boxsize,nboxes, &
         nlevels,ntargbox,norder,nt2, targets, fcoeffs,&
         irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
         iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
         ntargboxout,nt2out,targetsout,fcoeffsout)

       implicit real *8 (a-h,o-z)
       integer, intent(inout) :: nboxes,nlevels,nt2
       integer, dimension(:), intent(in) ::  itree
       integer iptr(7)
       integer, dimension(:), intent(in) :: ntargbox
       real *8, dimension(:,:), intent(in) :: treecenters
       real *8, intent(in) :: boxsize(0:nlevels)
       real *8, dimension(:,:), intent(in) :: targets
       real *8, dimension(:), intent(in) :: fcoeffs
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
      end subroutine
end interface         

!      calling sequence variables
      real *8, intent(in) :: sources(3,ns)
      real *8, intent(in) :: dumtarg(3,nt)

      integer, allocatable, intent(inout) :: itree(:)
      integer, allocatable, intent(inout) :: iptr(:)
      integer, intent(inout) :: nboxes,nlevels,ltree,nlbox,nt2

      real *8, allocatable, intent(inout) :: treecenters(:,:)
      real *8, allocatable, intent(inout) :: boxsize(:)
      real *8, allocatable, intent(inout) :: fcoeffs(:)

      real *8, allocatable :: radsrc(:)
      integer, allocatable :: ntargbox(:)

      real *8, allocatable :: targets(:,:)
      real *8, allocatable :: fvals(:)
      real *8, allocatable :: fvalstmp(:)

      integer, allocatable :: itreetmp(:)
      integer ipointer(32)

      real *8 xpts(norder),wts(norder)
      real *8, allocatable :: xmatc(:,:)
      real *8, allocatable :: tails(:)
      integer, allocatable :: irefineflag(:)
      integer, allocatable :: ibcompflag(:)
      integer, allocatable :: itcompflag(:)

      integer, allocatable :: itreeout(:)
      integer, allocatable :: ntargboxout(:)
      integer iptrout(7)
      real *8, allocatable :: targetsout(:,:),fcoeffsout(:), &
            boxsizeout(:),treecentersout(:,:)
      


!
!c      temporary variables
!
      real *8 expc(3),radt

      external fun1

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
            nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,nlevels, &
            nboxes,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung , &
            ltreetmp)


       allocate(itreetmp(ltreetmp))
       allocate(boxsize(0:nlevels))
       allocate(treecenters(3,nboxes))

!       allocate(itreetmp(ltreetmp),boxsize(0:nlevels), &
!          treecenters(3,nboxes))



       call mklraptree(sources,ns,radsrc,dumtarg,nt,expc,nexpc, &
             radexp,idivflag,ndiv,isep,mhung,mnbors,mnlist1, &
             mnlist2,mnlist3,mnlist4,nlevels,nboxes,treecenters, &
             boxsize,itreetmp,ltreetmp,ipointer)



!
!c      extract relevant bits of the tree needed from this point on
!
        ltree = 2*(nlevels+1)+12*nboxes
        allocate(iptr(7),itree(ltree),ntargbox(nboxes))


        call extracttreesubinfo(itreetmp,ipointer,nboxes,nlevels, &
                itree,iptr,ntargbox)
        deallocate(itreetmp)
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
        allocate(fvalstmp(norder3))

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
                       targets(1,itarg) = treecenters(1,ibox)+ &
                                 xpts(kk)*boxsize(ilev)/2
                       targets(2,itarg) = treecenters(2,ibox)+ &
                                 xpts(jj)*boxsize(ilev)/2
                       targets(3,itarg) = treecenters(3,ibox)+ &
                                 xpts(ii)*boxsize(ilev)/2
                     enddo
                  enddo
              enddo
          endif
       enddo

!
!c       evaluate the function intially at the targets
!

        do i=1,nt2
           call fun1(targets(1,i),targets(2,i),targets(3,i),fvals(i))
        enddo

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
        endif
      enddo

!
!c         compute the coefficients
!
!
!

      iflag = 1

      maxit = 100

      do iter = 1,maxit
        call prinf('iter=*',iter,1)
        rmaxerr = 0
        iflag = 0 
        allocate(tails(nboxes),irefineflag(nboxes))
        do i=1,nboxes
           nchild = itree(iptr(4)+i-1)
           tails(i) = 0
           irefineflag(i) = 0
           if(nchild.eq.0.and.ntargbox(i).gt.0) then
              ii = itree(iptr(6)+i-1)
              call comptail(norder,norder3,fcoeffs(ii),tails(i))
              if(tails(i).gt.eps) then
                 irefineflag(i) = 1
                 iflag = 1
              endif
              if(tails(i).gt.rmaxerr) rmaxerr = tails(i)
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


         call reorganizechebtree(itree,ltree,iptr,treecenters, &
            boxsize,nboxes, &
            nlevels,ntargbox,norder,nt2, targets, fcoeffs,&
            irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
            iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
            ntargboxout,nt2out,targetsout,fcoeffsout)


            deallocate(itree,treecenters,boxsize,ntargbox,targets, &
               fcoeffs)

            ltree = ltreeout
            nboxes = nboxesout
            nlevels = nlevelsout
            nt2 = nt2out

!            call prinf('ltree=*',ltree,1)
!            call prinf('nboxes=*',nboxes,1)
!            call prinf('nlevels=*',nlevels,1)
!            call prinf('nt2=*',nt2,1)
            

            allocate(itree(ltree),treecenters(3,nboxes), &
               boxsize(0:nlevels),targets(3,nt2),fcoeffs(nt2), &
               ntargbox(nboxes))

            iptr = iptrout
            itree = itreeout
            treecenters = treecentersout
            boxsize(0:nlevels) = boxsizeout(0:nlevels)
            targets = targetsout
            fcoeffs = fcoeffsout
            ntargbox = ntargboxout

!            call prinf('ibcompflag=*',ibcompflag,nboxes)
!            call prinf('ifcoeffs=*',itree(iptr(6)),nboxes)

            do i=1,nboxes
               if(ibcompflag(i).eq.1) then
!                  call prinf('i=*',i,1)
                  ii = itree(iptr(6)+i-1)
!                  call prin2('targets=*',targets(1,ii),3*norder3)
                  do j=1,norder3
                     iii = ii+j-1
                     call fun1(targets(1,iii),targets(2,iii),&
                        targets(3,iii),fvalstmp(j))
                  enddo
!                  call prin2('fvalstmp=*',fvalstmp,norder3)
                  call matvec(norder3,norder3,xmatc,fvalstmp,& 
                     fcoeffs(ii))
!                  call prin2('fcoeffs=*',fcoeffs(ii),norder3)
               endif
            enddo

            deallocate(tails,irefineflag,itreeout,treecentersout, &
             boxsizeout,ntargboxout,targetsout,fcoeffsout) 
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
      subroutine reorganizechebtree(itree,ltree,iptr,treecenters, &
         boxsize,nboxes, &
         nlevels,ntargbox,norder,nt2, targets,fcoeffs,&
         irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
         iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
         ntargboxout,nt2out,targetsout,fcoeffsout)

       implicit real *8 (a-h,o-z)
       integer, intent(inout) :: nboxes,nlevels,nt2
       integer, dimension(:), intent(in) ::  itree
       integer iptr(7)
       integer, dimension(:), intent(in) :: ntargbox
       real *8, dimension(:,:), intent(in) :: treecenters
       real *8, intent(in) :: boxsize(0:nlevels)
       real *8, dimension(:,:), intent(in) :: targets
       real *8, dimension(:), intent(in) :: fcoeffs
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
          itcompflag(nt2out))

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
