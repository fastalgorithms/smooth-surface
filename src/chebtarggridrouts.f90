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
      subroutine reorganizechebtree(nboxnew,nlevnew,ltreenew, &
         itree,iptr,treecenters,boxsize,nboxes, &
         nlevels,nlbox,norder,nt2, &
         targets,fcoeffs,irefineflag,icompflag)

       implicit real *8 (a-h,o-z)
       integer, intent(in) ::  nboxnew,nlevnew,ltreenew
       integer, dimension(:), intent(inout) ::  itree
       integer, iptr(7)
       real *8, dimension(:,:), intent(inout) :: treecenters
       real *8, dimension(:), intent(inout) :: boxsize
       integer nboxes,nlevels,norder,nt2
       real *8, dimension(:,:), intent(inout) :: targets
       real *8, dimension(:), intent(inout) :: fcoeffs
       integer, icompflag(nboxnew),irefineflag(nboxes)
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

      integer, allocatable :: itreetmp(:)
      integer ipointer(32)

      real *8 xpts(norder),wts(norder)
      real *8, allocatable :: xmatc(:,:)
      real *8, allocatable :: tails(:)
      integer, allocatable :: irefineflag(:)

      integer, allocatable :: laddrtail(:,:)

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

        

        ictr = 1
        norder3 = norder*norder*norder
        
        call prinf('norder3=*',norder3,1)
        call prinf('nlbox=*',nlbox,1)
        nt2 = norder3*nlbox
        call prinf('nt2=*',nt2,1)


        allocate(targets(3,nt2),fvals(nt2),fcoeffs(nt2))

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

!
!c         compute the coefficients
!
!
!
      allocate(tails(nboxes),irefineflag(nboxes))
      do i=1,nboxes
        nchild = itree(iptr(4)+i-1)
        tails(i) = 0
        irefineflag(i) = 0
        if(nchild.eq.0.and.ntargbox(i).gt.0) then
           ii = itree(iptr(6)+i-1)
           call matvec(norder3,norder3,xmatc,fvals(ii),fcoeffs(ii))
           call comptail(norder,norder3,fcoeffs,tails(i))
           if(tails(i).gt.eps) irefineflag(i) = 1
        endif
      enddo
      
      call prin2('tails=*',tails,nboxes)
      call prinf('irefineflag=*',irefineflag,nboxes)


      do i=1,nboxes

         irefineflag(i) = 0
      enddo
  
      irefineflag(6) = 1

      allocate(laddrtail(2,0:nlevels))

!
!!       count number of additional boxes to be formed
!
      
      nextra = 0
      nlevnew = nlevels
      do i=1,nboxes
         if(irefineflag(i).eq.1) nextra = nextra + 1
         nchild = itree(iptr(4)+i-1)
         if(nchild.eq.0.and.ntargbox(i).gt.0) nlevnew = nlevels+1
      enddo

      nboxnew = nboxes + nextra

      ltreenew = 2*(nlevnew+1)+12*(nboxnew)

      allocate(icompflag(nboxnew))

      call reorganizechebtree(nextra,nlevnew,ltreenew, &
         itree,iptr,treecenters,boxsize,nboxes, &
         nlevels,nlbox,norder,nt2, &
         targets,fcoeffs,irefineflag,icompflag)





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

    do ii = 1,norder
      do jj = 1,norder
         do kk = 1,norder

            i = ii+jj+kk
            if(i.ge.norder) then
                iind = (ii-1)*norder*norder + (jj-1)*norder+kk
                tail = tail + abs(fcoeffs(iind))**2
            endif
         enddo
      enddo
    enddo

    tail = sqrt(tail)
end subroutine

subroutine matvec(m,n,a,x,y)
implicit real *8 (a-h,o-z)
   real *8 a(m,*),x(*),y(*)
   
   do i=1,m
     y(i) = 0
     do j=1,n

        y(i) = y(i) + a(i,j)*x(j)

     enddo

   enddo


end subroutine
!-----------------------------------
      subroutine reorganizechebtree(nboxnew,nlevnew,ltreenew, &
         itree,iptr,treecenters,boxsize,nboxes, &
         nlevels,nlbox,norder,nt2, &
         targets,fcoeffs,irefineflag,icompflag)

       implicit real *8 (a-h,o-z)
       integer, intent(in) ::  nboxnew,nlevnew,ltreenew
       integer, dimension(:), intent(inout) ::  itree
       integer, iptr(7)
       real *8, dimension(:,:), intent(inout) :: treecenters
       real *8, dimension(:), intent(inout) :: boxsize
       integer nboxes,nlevels,norder,nt2
       real *8, dimension(:,:), intent(inout) :: targets
       real *8, dimension(:), intent(inout) :: fcoeffs
       integer, icompflag(nboxnew),irefineflag(nboxes)



end subroutine
