
      subroutine gettargandtreemem(ns,sources,nt,targets,lmem,nlbox,
     1             nlevels,nboxes)
c
c       this subroutine tree sorts a collection of sources
c        and targets, and returns the length of the tree used
c        and the number of leaf boxes, number of levels in the tree,
c        and the number of boxes in the tree
c        
c
c
c        input: 
c        ns - number of sources
c        sources(3,*) - xyz locations of the sources
c        nt - number of targets
c        targets(3,*) - xyz locations of the targets
c      
c        output
c        lmem - memory required for allocating the tree
c        nlbox - number of leaf boxes
c        nlevels - number of levels in the tree
c        nboxes - total number of boxes
c
      
      implicit real *8 (a-h,o-z)
      real *8 sources(3,*),targets(3,*)
      real *8, allocatable :: radsrc(:)
c
cc      temporary variables
c
      real *8 expc(3),radt
      integer ipointer(32)
      integer, allocatable :: itree(:)
      real *8, allocatable :: boxsize(:),treecenters(:,:)

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

      call mklraptreemem(ier,sources,ns,radsrc,targets,nt,expc,
     1       nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,nlevels,
     2       nboxes,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung,
     3       ltree)

       allocate(itree(ltree),boxsize(0:nlevels),treecenters(3,nboxes)) 

       call mklraptree(sources,ns,radsrc,targets,nt,expc,nexpc,
     1        radexp,idivflag,ndiv,isep,mhung,mnbors,mnlist1,
     2        mnlist2,mnlist3,mnlist4,nlevels,nboxes,treecenters,
     3        boxsize,itree,ltree,ipointer)


        lmem = ltree
        nlbox = 0
        do i=1,nboxes
           nchild = itree(ipointer(3)+i-1)
           if(nchild.eq.0) then
              istart = itree(ipointer(12)+i-1)
              iend = itree(ipointer(13)+i-1)
              ntarg = iend - istart + 1
              if(ntarg.gt.0)  nlbox = nlbox + 1
           endif
        enddo

      return
      end
c-----------------------------------
      subroutine gettargandtree(ns,sources,nt,targets,
     1      itree,ipointer,treecenters,boxsize,norder,targets2,iptr)
c
c       this subroutine tree sorts a collection of sources
c        and targets, and returns the tree  
c        
c
c
c        input: 
c        ns - number of sources
c        sources(3,*) - xyz locations of the sources
c        nt - number of targets
c        targets(3,*) - xyz locations of the targets
c      
c        output
c        lmem - memory required for allocating the tree
c        nlbox - number of leaf boxes
c        nlevels - number of levels in the tree
c        nboxes - total number of boxes
c
      
      implicit real *8 (a-h,o-z)
      real *8 sources(3,*),targets(3,*)
      real *8, allocatable :: radsrc(:)

      integer itree(*),ipointer(32),iptr(*)
      real *8 boxsize(0:*),treecenters(3,*)
      real *8 xpts(norder),wts(norder)
      real *8 targets2(3,*)
c
cc      temporary variables
c
      real *8 expc(3),radt

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

      call mklraptreemem(ier,sources,ns,radsrc,targets,nt,expc,
     1       nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,nlevels,
     2       nboxes,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung,
     3       ltree)

       call mklraptree(sources,ns,radsrc,targets,nt,expc,nexpc,
     1        radexp,idivflag,ndiv,isep,mhung,mnbors,mnlist1,
     2        mnlist2,mnlist3,mnlist4,nlevels,nboxes,treecenters,
     3        boxsize,itree,ltree,ipointer)

       call prin2('boxsize=*',boxsize,nlevels+1)

        itype = 1
        call chebexps(itype,norder,xpts,utmp,vtmp,wts)

        call prin2('xpts=*',xpts,norder)

        do i=1,nboxes
           iptr(i) = -1
        enddo

        itarg = 0

        ictr = 1
        do ilev=0,nlevels
           do ibox = itree(2*ilev+1),itree(2*ilev+2)
              nchild = itree(ipointer(3)+ibox-1)
              if(nchild.eq.0) then
                 istart = itree(ipointer(12)+ibox-1)
                 iend = itree(ipointer(13)+ibox-1)
                 ntarg = iend-istart+1
                 if(ntarg.gt.0) then

                     iptr(ibox) = ictr
                     ictr = ictr + norder*norder*norder

                     do ii =1,norder
                        do jj= 1,norder
                           do kk =1,norder
                               itarg = itarg + 1
                               targets2(1,itarg) = treecenters(1,ibox)+
     1                            xpts(kk)*boxsize(ilev)/2
                               targets2(2,itarg) = treecenters(2,ibox)+
     1                            xpts(jj)*boxsize(ilev)/2
                               targets2(3,itarg) = treecenters(3,ibox)+
     1                            xpts(ii)*boxsize(ilev)/2
                           enddo
                        enddo
                     enddo
                 endif
              endif
           enddo
        enddo

      return
      end
c-------------------------------------------     

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
       end
c---------------------------------------------------------

      subroutine findboxtarg(target,iboxtarg,ilevel,itree,ipointer,
     1     boxsize,treecenters,nboxes,nlevels)
      implicit real *8 (a-h,o-z)
      integer ipointer(32),itree(*)
      real *8 boxsize(0:nlevels),treecenters(3,*),target(3)

      iboxtarg = 1
      ilevel = 0
      do ilev=0,nlevels-1
c
cc       compute relartive corordinates for teh target with
c        respect to the box
c
         nchild = itree(ipointer(3)+iboxtarg-1)
         if(nchild.eq.0) goto 3000
         xtmp = target(1) - treecenters(1,iboxtarg)
         ytmp = target(2) - treecenters(2,iboxtarg)
         ztmp = target(3) - treecenters(3,iboxtarg)
         if(xtmp.le.0.and.ytmp.le.0.and.ztmp.le.0) ichild=1
         if(xtmp.gt.0.and.ytmp.le.0.and.ztmp.le.0) ichild=2
         if(xtmp.le.0.and.ytmp.gt.0.and.ztmp.le.0) ichild=3
         if(xtmp.gt.0.and.ytmp.gt.0.and.ztmp.le.0) ichild=4
         if(xtmp.le.0.and.ytmp.le.0.and.ztmp.gt.0) ichild=5
         if(xtmp.gt.0.and.ytmp.le.0.and.ztmp.gt.0) ichild=6
         if(xtmp.le.0.and.ytmp.gt.0.and.ztmp.gt.0) ichild=7
         if(xtmp.gt.0.and.ytmp.gt.0.and.ztmp.gt.0) ichild=8

         iboxtarg = itree(ipointer(4)+(iboxtarg-1)*8+ichild-1)
         ilevel = ilev + 1

      enddo
 3000 continue      

      return
      end
c------------------------------------------------------------      


      subroutine cheb3deval(target,rscale,center,k,fcoeffs,f)
      implicit real *8 (a-h,o-z)
      real *8 xpols(k+10),ypols(k+10),zpols(k+10)
      real *8 fcoeffs(*),target(3),center(3)

      x = (target(1) - center(1))/rscale*2
      y = (target(2) - center(2))/rscale*2
      z = (target(3) - center(3))/rscale*2

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
      end
c--------------------------------------------      
