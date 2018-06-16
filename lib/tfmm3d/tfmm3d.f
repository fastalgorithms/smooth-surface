cc Copyright (C) 2010-2011: Leslie Greengard and
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$

        subroutine tfmm3d(ier,iprec,ns,source,   
     $     ifcharge,charge,ifdipole,dipstr,dipvec,
     $     nt,targ,trads,stdev,stdev_grad,
     $     itree,ltree,ipointer,isep,ndiv,nlevels,nboxes,boxsize,
     $     centers,ifpottarg,pottarg,
     $     ifgradtarg,gradtarg)
c       
c       
c       FMM in R^3 for Surface Smoothing kernel 
c
c       erf(|x-y|/(sqrt(2)*stdev(x)))/(4* pi * |x-y|)
c       where y is a source point and x is a target.
c
c       On input, trads(j) for target x_j ix assumed to be 
c       so large that 
c       erf(|x-y|/(sqrt(2)*stdev(x))) = 1, so that the
c       smoothing kernel is the usual Coulomb source outside
c       a disk of radius trad(j).
c
c
c    The tree is built so that targets are hung at a level in the
c    tree so that the disk of radius trads(j) is contained within the
c    nearest neighbor boxes. 
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   ier:    error return code [NOT THOUGHT THROUGH YET].
c   iprec:  FMM precision flag
c
c                 -2 => tolerance =.5d0
c                 -1 => tolerance =.5d-1
c                  0 => tolerance =.5d-2
c                  1 => tolerance =.5d-3
c                  2 => tolerance =.5d-6
c                  3 => tolerance =.5d-9
c                  4 => tolerance =.5d-12
c
c   ns:     integer:          number of sources
c   source: real *8 (3,ns):   source locations
c   ifcharge:                 charge computation flag
c                             ifcharge = 1 => include charge contribution
c                             otherwise do not
c   charge: complex *16 (ns): charge strengths
c   ifdipole:  		      dipole computation flag
c              		      ifdipole = 1 => include dipole contribution
c                             otherwise do not
c   dipstr: complex *16 (ns): dipole strengths
c   dipvec: real *8 (3,ns):   dipole orientation vectors
c   nt: integer:  	      number of targets == expansion centers
c   targ: real *8 (3,nt):     target locations
c   trads: real *8 (nt):      near field distance -> halts sorting
c                             of target to finer levels.
c   stdev: real *8 (nt):      smoothing kernel parameter 
c   stdev_grad: real *8 (nt): grad of stdev 
c
c   itree    in: integer (ltree)
c            		      This array contains all the information
c             		      about the tree
c   ltree    in: integer
c            length of tree
c
c    ipointer in: integer(26)
c             ipointer is a collection of pointers 
c             which points to where different elements 
c             of the tree are stored in the itree array
c
c             itree(ipointer(1):ipointer(2)-1) == laddr(2,0:nlevels)
c             laddr(1,i) is the first box on level i
c             laddr(2,i) is the last box on level i
c
c             itree(ipointer(2):ipointer(3)-1) == iparent(nboxes)
c             iparent(i) is the parent of box i
c
c             itree(ipointer(3):ipointer(4)-1) == nchild(nboxes)
c             nchild(i) is the number of children of box i
c
c             itree(ipointer(4):ipointer(5)-1) == child(4,nboxes)
c             child(i,j) is the ith child of box j
c
c             itree(ipointer(5):ipointer(6)-1) == nnbor(nboxes)
c             nnbor(i) is the number of nearest neighbors of box i
c
c             itree(ipointer(6):ipointer(7)-1) == nbors(mnbors,nboxes)
c             nnbor(i,j) is the ith neighbor of box j
c
c             itree(ipointer(7):ipointer(8)-1) == nlist2(nboxes)
c             nlist2(i) is the number of boxes in the interaction list
c             of box i
c 
c             itree(ipointer(8):ipointer(9)-1) == list2(mnlist2,nboxes)
c             list2(i,j) is the ith box in the interaction list of
c             box j
c
c             itree(ipointer(9):ipointer(10)-1) == isource(ns)
c             isource is the mapping to tree sort source chunks
c
c             itree(ipointer(10):ipointer(11)-1) == itarget(nt)
c             itarget is the mapping to tree sort targets
c
c             itree(ipointer(11):ipointer(12)-1) == iexpc(nexpc)
c             iexpc is the mapping to tree sort expansion centers
c
c             itree(ipointer(12):ipointer(13)-1) == ihsfirst(nboxes)
c             ihsfirst(i) is the location in isource of the first
c             hung chunk in box i
c 
c             itree(ipointer(13):ipointer(14)-1) == ihslast(nboxes)
c             ihslast(i) is the location in isource of the last hung
c             chunk in box i
c
c             itree(ipointer(14):ipointer(15)-1) == isfirst(nboxes)
c             isfirst(i) is the location in isource of the first
c             source in box i
c 
c             itree(ipointer(15):ipointer(16)-1) == islast(nboxes)
c             islast(i) is the location in isource of the last
c             source in box i
c
c             itree(ipointer(16):ipointer(17)-1) == itfirst(nboxes)
c             itfirst(i) is the location in itarget of the first
c             target in box i
c 
c             itree(ipointer(17):ipointer(18)-1) == itlast(nboxes)
c             itlast(i) is the location in itarget of the last
c             target in box i
c
c             itree(ipointer(18):ipointer(19)-1) == ihefirst(nboxes)
c             ihsfirst(i) is the location in iexpc of the first
c             hung expansion center in box i
c 
c             itree(ipointer(19):ipointer(20)-1) == ihelast(nboxes)
c             ihslast(i) is the location in iexpc of the last hung
c             expansion center in box i
c
c             itree(ipointer(20):ipointer(21)-1) == iefirst(nboxes)
c             iefirst(i) is the location in iexpc of the first
c             expansion center in box i
c 
c             itree(ipointer(21):ipointer(22)-1) == ielast(nboxes)
c             ielast(i) is the location in iexpc of the last
c             expansion center in box i
c
c             itree(ipointer(22):ipointer(23)-1) == nhungsrc(nboxes)
c             nhungsrc(i) is the number of hung chunks in box i
c
c             itree(ipointer(23):ipointer(24)-1) == nhungexp(nboxes)
c             nhungexp(i) is the number of hung expansion centers 
c             in box i
c 
c             NOTE: The next two parameters are irrelevant for
c             this code
c             itree(ipointer(24):ipointer(25)-1) == nhunglistsrc(nboxes)
c             nhunglistsrc(i) is the number of hung sources relevant
c             to box i
c 
c             itree(ipointer(25):ipointer(26)-1) == 
c             ihunglistsrc(mhung,nboxes)
c             ihunglistsrc(i,j) = src id of ith hung source relevant
c             to box j
c
c     ndiv    in: integer
c             Max number of chunks per box
c
c     isep    in: integer
c             separation parameter. isep=1/2, list1 comprises
c             of 1/2 nearest neighbors respectively. If isep=1,
c             mnbors = 27 and mnlist2 = 27*7. If isep=2,
c             mnbors = 125 and mnlist2 = 125*7
c
c     nlevels in: integer
c             number of levels in the tree
c
c     
c     nboxes  in: integer
c             number of boxes in the tree
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c
c     centers in: real *8(3,nboxes)
c                 array containing the centers of all the boxes
c
c     ifpottarg   in: integer
c              flag for evaluating potential at the targets
c              The potential will be evaluated if ifpottarg = 1
c
c     ifgradtarg   in:integer
c                flag for evaluating the gradient at the targets
c                gradient will be evaluated at the targets if
c                ifgradtarg = 1
c
c   OUTPUT
c
c   ier   =  error return code
c
c   pottarg:    out: complex *16(nt) 
c               potential at the target locations
c
c   gradtarg:   out: complex *16(3,nt)
c               gradient at the target locations
c
c------------------------------------------------------------------

      implicit none

      integer ier,iprec,ns,nt
      integer ifcharge,ifdipole
      integer ifpottarg,ifgradtarg,ifgradaatarg

      integer i,j,k,l
      integer ibox,jbox,ilev,nmax


      real *8 source(3,1)
      real *8 targ(3,1)

      complex *16 charge(1)
      complex *16 dipstr(1)
      real *8 dipvec(3,1)

      real *8 trads(*)
      real *8 stdev(*)
      real *8 stdev_grad(3,*)

      complex *16 pottarg(1)
      complex *16 gradtarg(3,1)

c       Tree variables
      integer itree(1),ltree,ipointer(26),ndiv,mhung,nboxes
      integer nlevels, isep
      real *8 centers(3,1), boxsize(0:nlevels)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)

      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: targsort(:,:)

      complex *16, allocatable :: chargesort(:)
      complex *16, allocatable :: dipstrsort(:)
      real *8, allocatable :: dipvecsort(:,:)

      real *8, allocatable :: tradsort(:)
      real *8, allocatable :: stdevsort(:)
      real *8, allocatable :: stdev_gradsort(:,:)

      complex *16, allocatable :: potsort(:)
      complex *16, allocatable :: gradsort(:,:)

      real *8, allocatable :: rmlexp(:)
      complex *16, allocatable :: mptemp(:,:)
      complex *16, allocatable :: mptemp2(:,:)
      integer, allocatable :: iaddr(:,:)

      real *8 time1,time2
c
      real *8 timeinfo(10),epsfmm
      real *8 scales(0:nlevels)
      integer ifprint,lmptot, lmptemp


      ier=0
c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
      ifprint=0
c
c     set fmm tolerance based on iprec flag.
c       
      if( iprec .eq. -2 ) epsfmm=.5d-0 
      if( iprec .eq. -1 ) epsfmm=.5d-1
      if( iprec .eq. 0 ) epsfmm=.5d-2
      if( iprec .eq. 1 ) epsfmm=.5d-3
      if( iprec .eq. 2 ) epsfmm=.5d-6
      if( iprec .eq. 3 ) epsfmm=.5d-9
      if( iprec .eq. 4 ) epsfmm=.5d-12
      if( iprec .eq. 5 ) epsfmm=.5d-15
      if( iprec .eq. 6 ) epsfmm=0
c      
      if(ifprint .eq. 1) call prin2('epsfmm=*',epsfmm,1)

c     Allocate sorted source and target arrays      

      allocate(sourcesort(3,ns))
      allocate(targsort(3,nt))
      allocate(chargesort(ns))
      allocate(dipstrsort(ns))
      allocate(dipvecsort(3,ns))

c     Allocate sorted expansions      
      allocate(stdevsort(nt),tradsort(nt))
      allocate(stdev_gradsort(3,nt))

      allocate(potsort(nt))
      allocate(gradsort(3,nt))


c     Figure out scales of expansions at all levels
      do ilev = 0,nlevels
          scales(ilev) = boxsize(ilev)
      enddo


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,nt
         potsort(i) = 0.0d0
         gradsort(1,i) = 0.0d0
         gradsort(2,i) = 0.0d0
         gradsort(3,i) = 0.0d0

      enddo
C$OMP END PARALLEL DO

c     Compute length of expansions at each level      
      nmax = 0
      do i=0,nlevels
         if(isep.eq.1) call l3dterms(epsfmm,nterms(i),ier)
         if(isep.eq.2) call l3dterms_far(epsfmm,nterms(i),ier)
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))
      allocate(mptemp(0:nmax,-nmax:nmax))
      allocate(mptemp2(0:nmax,-nmax:nmax))
      lmptemp = (nmax+1)*(2*nmax+1)

c     Unpack laddr array
      do i=0,nlevels
         laddr(1,i) = itree(2*i+1)
         laddr(2,i) = itree(2*i+2)
      enddo
c     End of unpacking laddr array      

c     Sort sources and targets
      call l3dreorder_newtree(ns,source,ifcharge,charge,
     $    itree(ipointer(9)),ifdipole,dipstr,dipvec,sourcesort,
     $    chargesort,dipstrsort,dipvecsort)

c     Sort expansion centers/targets
      call l3dreordertarg_newtree(nt,targ,trads,stdev,stdev_grad,
     $  itree(ipointer(11)),targsort,tradsort,stdevsort,stdev_gradsort)


c
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call l3dmpalloc_newtree(laddr,iaddr,nlevels,lmptot,nterms)
      if(1 .eq. 1) call prinf(' lmptot is *',lmptot,1)

      allocate(rmlexp(lmptot),stat=ier)
      if(ier.ne.0) then
         call prinf('Cannot allocate mpole expansion workspace,
     1              lmptot is *', lmptot,1)
         ier = 16
         return
      endif

c     Memory allocation is complete. 
c     Call main fmm routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call tfmm3dmain(ier,iprec,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipstrsort,dipvecsort,
     $   nt,targsort,tradsort,stdevsort,stdev_gradsort,
     $   epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,isep,ndiv,nlevels,
     $   nboxes,boxsize,scales,
     $   centers,laddr,nterms,
     $   ifpottarg,potsort,ifgradtarg,gradsort)

cc      call prin2('texpssort=*',texpssort,180)

      print *, "Have exited main fmm routine"

cccc      time2=second()
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)

c
c     parameter ier from targmain routine is currently 
c     meaningless, reset to 0
      if( ier .ne. 0 ) ier = 0


      if(ifpottarg.eq.1) call l3dpsort_newtree(nt,
     1                   itree(ipointer(11)),potsort,pottarg)
      if(ifgradtarg.eq.1) call l3dfsort_newtree(nt,
     1                   itree(ipointer(11)),gradsort,gradtarg)
      return
      end
c
cc-----------------------------------------------------------------------
cc
      subroutine tfmm3dmain(ier,iprec,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ntarget,targetsort,
     $     tradsort,stdevsort,stdev_gradsort,
     $     epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $     itree,ltree,ipointer,isep,ndiv,nlevels, 
     $     nboxes,boxsize,scales,
     $     centers,laddr,nterms,
     $     ifpottarg,pottarg,ifgradtarg,gradtarg)
      implicit none

      integer ier,iprec
      integer nsource,ntarget,nexpc
      integer ndiv,nlevels

      integer ifcharge,ifdipole
      integer ifpottarg,ifgradtarg
      real *8 epsfmm

      real *8 sourcesort(3,*)

      complex *16 chargesort(*)
      complex *16 dipstrsort(*)
      real *8 dipvecsort(3,*)

      real *8 targetsort(3,*)
      real *8 tradsort(*),stdevsort(*)
      real *8 stdev_gradsort(3,*)


      complex *16 pottarg(*)
      complex *16 gradtarg(3,*)

      integer iaddr(2,*), lmptot, lmptemp
      real *8 rmlexp(*)
      complex *16 mptemp(*)
      complex *16 mptemp2(*)
       
      real *8 timeinfo(10)
      real *8 centers(3,*)

      integer isep, ltree
      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer ipointer(26)
      integer itree(*)
      integer nboxes
      real *8 scales(0:nlevels)
      real *8 boxsize(0:nlevels)

      integer nuall,ndall,nnall,nsall,neall,nwall
      integer nu1234,nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8

      integer uall(200),dall(200),nall(120),sall(120),eall(72),wall(72)
      integer u1234(36),d5678(36),n1256(24),s3478(24)
      integer e1357(16),w2468(16),n12(20),n56(20),s34(20),s78(20)
      integer e13(20),e57(20),w24(20),w68(20)
      integer e1(20),e3(5),e5(5),e7(5),w2(5),w4(5),w6(5),w8(5)

c     temp variables
      integer i,j,k,l,ii,jj,kk,ll
      integer ibox,jbox,ilev,npts
      integer nchild,nlist1,nlist2

      integer istart,iend
      integer isstart,isend,jsstart,jsend
      integer jstart,jend

      integer ifprint

      real *8 d,time1,time2
      complex *16 pottmp,gradtmp(3)

c     PW variables
      integer ntmax, nexpmax, nlams, nmax, nthmax, nphmax
      parameter (ntmax = 1000)
      parameter (nexpmax = 16000)
      real *8, allocatable :: carray(:,:), dc(:,:), rdplus(:,:,:)
      real *8, allocatable :: rdminus(:,:,:), rdsq3(:,:,:)
      real *8, allocatable :: rdmsq3(:,:,:)
  
      real *8 rlams(ntmax), whts(ntmax)

      real *8, allocatable :: rlsc(:,:,:)
      integer nfourier(ntmax), nphysical(ntmax)
      integer nexptot, nexptotp
      complex *16, allocatable :: xshift(:,:)
      complex *16, allocatable :: yshift(:,:)
      real *8, allocatable :: zshift(:,:)

      complex *16 fexpe(100000), fexpo(100000), fexpback(400000)
      complex *16, allocatable :: lexp(:,:,:)
      complex *16, allocatable :: mexpf1(:),mexpf2(:)
      complex *16, allocatable :: mexpp1(:),mexpp2(:),mexppall(:,:)

      complex *16, allocatable :: tmp(:,:)
      
      real *8 sourcetmp(3)
      complex *16 chargetmp

      integer ix,iy,iz
      real *8 rtmp
      complex *16 zmul

      integer nlege, lw7, lused7, itype
      real *8 wlege(40000)
      integer nterms_eval(4,0:nlevels)

      integer mnlist1, mnlist2
      complex *16 eye, ztmp
      real *8 alphaj
      integer ctr
      real *8 pi,errtmp

      pi = 4.0d0*datan(1.0d0)
c     Max number of elements in list1 and list2
      if(isep.eq.1) mnlist1 = 27
      if(isep.eq.2) mnlist1 = 125

      mnlist2 = 7*mnlist1

c     Initialize routines for plane wave mp loc translation

      if(isep.eq.1) then
         if(iprec.le.1) nlams = 14
         if(iprec.eq.2) nlams = 23
         if(iprec.eq.3) nlams = 32
         if(iprec.eq.4) nlams = 37
      endif
      if(isep.eq.2) then
         if(iprec.le.1) nlams = 9
         if(iprec.eq.2) nlams = 15
         if(iprec.eq.3) nlams = 22
         if(iprec.eq.4) nlams = 29
      endif
c      if(isep.eq.2) then
c         if(iprec.le.1) nlams = 10
c         if(iprec.eq.2) nlams = 20
c         if(iprec.eq.3) nlams = 30
c         if(iprec.eq.4) nlams = 40
c      endif

      nmax = 0
      do i=0,nlevels
         if(nmax.lt.nterms(i)) nmax = nterms(i)
      enddo
      allocate(carray(4*nmax+1,4*nmax+1))
      allocate(dc(0:4*nmax,0:4*nmax))
      allocate(rdplus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdminus(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rdmsq3(0:nmax,0:nmax,-nmax:nmax))
      allocate(rlsc(1:nlams,0:nmax,0:nmax))


c     generate rotation matrices and carray
      call rotgen(nterms,carray,rdplus,rdminus,rdsq3,rdmsq3,dc)

c     generate rlams and weights (these are the nodes
c     and weights for the lambda integral)

      if(isep.eq.1) call vwts(rlams,whts,nlams)
      if(isep.eq.2) call lwtsexp3sep2(nlams,rlams,whts,errtmp)

c     generate the number of fourier modes required to represent the
c     moment function in fourier space

      if(isep.eq.1) call numthetahalf(nfourier,nlams)
      if(isep.eq.2) call numthetahalf_isep2(nfourier,nlams)
 
c     generate the number of fourier modes in physical space
c     required for the exponential representation
      if(isep.eq.1) call numthetafour(nphysical,nlams)
      if(isep.eq.2) call numthetasix(nphysical,nlams)

c     Generate powers of lambda for the exponential basis
      call rlscini(rlsc,nlams,rlams,nmax)

c     Compute total number of plane waves
      nexptotp = 0
      nexptot = 0
      nthmax = 0
      nphmax = 0
      do i=1,nlams
         nexptot = nexptot + nfourier(i)
         nexptotp = nexptotp + nphysical(i)
         if(nfourier(i).gt.nthmax) nthmax = nfourier(i)
         if(nphysical(i).gt.nphmax) nphmax = nphysical(i)
      enddo
      allocate(tmp(0:nmax,-nmax:nmax))

      allocate(xshift(5,nexptotp))
      allocate(yshift(5,nexptotp))
      allocate(zshift(5,nexptotp))
 
      allocate(mexpf1(nexptot),mexpf2(nexptot),mexpp1(nexptotp))
      allocate(mexpp2(nexptotp),mexppall(nexptotp,16))

      allocate(lexp(nexptotp,nboxes,6))


c     Precompute table for shifting exponential coefficients in 
c     physical domain
      call mkexps(rlams,nlams,nphysical,nexptotp,xshift,yshift,zshift)

c     Precompute table of exponentials for mapping from
c     fourier to physical domain
      call mkfexp(nlams,nfourier,nphysical,fexpe,fexpo,fexpback)
      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c       
        ifprint=1
c
c       
        do i=1,10
          timeinfo(i)=0
        enddo

c
c       ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            call l3dzero(rmlexp(iaddr(1,ibox)),nterms(ilev))
            call l3dzero(rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO         
       enddo

c    initialize legendre function evaluation routines
      nlege = 100
      lw7 = 40000
      call ylgndrfwini(nlege,wlege,lw7,lused7)

c    initialize nterms_eval
      do ilev=0,nlevels
         do itype = 1,4
            call l3dterms_eval(itype,epsfmm,
     1           nterms_eval(itype,ilev),ier)
         enddo
      enddo
       
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
      call cpu_time(time1)
cccc      time1=second()
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions



C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts)

      do ibox=laddr(1,nlevels),laddr(2,nlevels)
         istart = itree(ipointer(14)+ibox-1)
         iend = itree(ipointer(15)+ibox-1)
         npts = iend-istart+1

         if(npts.gt.0) then
            if(ifcharge.eq.1) then
               call l3dformmp_add_trunc(ier,scales(nlevels),
     1         sourcesort(1,istart),chargesort(istart),npts,
     2         centers(1,ibox),nterms(nlevels), nterms(nlevels),
     3         rmlexp(iaddr(1,ibox)),wlege,nlege)          
            endif
            if(ifdipole.eq.1) then
               call l3dformmp_dp_add_trunc(ier,scales(nlevels),
     1         sourcesort(1,istart),dipstrsort(istart),
     2         dipvecsort(1,istart),          
     3         npts,centers(1,ibox),
     4         nterms(nlevels),nterms(nlevels),
     5         rmlexp(iaddr(1,ibox)),wlege,nlege)
            endif
         endif
      enddo

C$OMP END PARALLEL DO        

      call cpu_time(time2)
cccc      time2=second()
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1
c       
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 2 (merge mp) ====*',i,0)
cccc      time1=second()
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,0,-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i)

         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  istart = itree(ipointer(14)+jbox-1)
                  iend = itree(ipointer(15)+jbox-1)
                  npts = iend-istart+1
                  if(npts.gt.0) then
                     call l3dmpmpquadu_add(scales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),scales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev),nterms(ilev),
     4               ier)
                  endif
               endif
            enddo
         enddo
C$OMP END PARALLEL DO        
      enddo

      call cpu_time(time2)
cccc      time2=second()
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 3 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
cccc      time1 = second()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(i)
      do i=1,nboxes
         do j=1,nexptotp
            do k=1,6
               lexp(j,i,k) = 0.0d0
            enddo
         enddo
      enddo
C$OMP END PARALLEL DO        

      do ilev = 1,nlevels-1
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istart = itree(ipointer(14)+ibox-1)
            iend = itree(ipointer(15)+ibox-1)
            npts = iend-istart+1
            nchild = itree(ipointer(3)+ibox-1)
   
            if(npts.gt.0.and.nchild.gt.0) then
               call getpwlistall(ibox,boxsize(ilev+1),nboxes,
     1         itree(ipointer(5)+ibox-1),itree(ipointer(6)+
     2         mnlist1*(ibox-1)),nchild,itree(ipointer(4)),centers,
     3         isep,nuall,uall,ndall,dall,nnall,nall,nsall,sall,neall,
     4         eall,nwall,wall,nu1234,u1234,nd5678,d5678,nn1256,n1256,
     5         ns3478,s3478,ne1357,e1357,nw2468,w2468,nn12,n12,nn56,n56,
     6         ns34,s34,ns78,s78,ne13,e13,ne57,e57,nw24,w24,nw68,w68,
     7         ne1,e1,ne3,e3,ne5,e5,ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)
 
               do i=1,nexptot
                  mexpf1(i) = 0.0d0
                  mexpf2(i) = 0.0d0
               enddo

               do i=1,nexptotp
                  mexpp1(i) = 0.0d0
                  mexpp2(i) = 0.0d0
                  do j=1,16
                     mexppall(i,j) = 0.0d0
                  enddo
               enddo

               call mkudexp(ibox,ilev,nboxes,itree(ipointer(4)),
     1         scales(ilev+1),nterms(ilev+1),nmax,rmlexp,iaddr,rlams,
     2         nlams,nfourier,nphysical,nthmax,nexptot,nexptotp,mexpf1,
     3         mexpf2,mexpp1,mexpp2,mexppall(1,1),mexppall(1,2),
     4         mexppall(1,3),mexppall(1,4),xshift,yshift,zshift,fexpe,
     5         fexpo,rlsc)

               call processup(ibox,scales(ilev+1),lexp,centers,nboxes,
     1         nuall,uall,nu1234,u1234,mexppall(1,1),mexppall(1,2),
     2         xshift,yshift,zshift,nexptotp)

               call processdn(ibox,scales(ilev+1),lexp,centers,nboxes,
     1         ndall,dall,nd5678,d5678,mexppall(1,3),mexppall(1,4),
     2         xshift,yshift,zshift,nexptotp) 

c             Process north and south interactions
              call mknsexp(ibox,ilev,nboxes,itree(ipointer(4)),
     1        scales(ilev+1),nterms(ilev+1),nmax,rmlexp,iaddr,mptemp,
     2        mptemp2,rlams,nlams,nfourier,nphysical,
     1        nthmax,nexptot,nexptotp,mexpf1,mexpf2,mexpp1,
     2        mexpp2,rdminus,mexppall(1,1),mexppall(1,2),
     3        mexppall(1,3),mexppall(1,4),mexppall(1,5),mexppall(1,6),
     5        mexppall(1,7),mexppall(1,8),xshift,yshift,zshift,
     6        fexpe,fexpo,rlsc)
   
              call processno(ibox,scales(ilev+1),lexp,centers,nboxes,
     1        nnall,nall,nn1256,n1256,nn12,n12,nn56,n56,mexppall(1,1),
     2        mexppall(1,2),mexppall(1,3),mexppall(1,4),
     3        xshift,yshift,zshift,nexptotp)
            
              call processso(ibox,scales(ilev+1),lexp,centers,nboxes,
     1        nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexppall(1,5),
     2        mexppall(1,6),mexppall(1,7),mexppall(1,8),xshift,
     3        yshift,zshift,nexptotp)

              call mkewexp(ibox,ilev,nboxes,itree(ipointer(4)),
     1        scales(ilev+1),nterms(ilev+1),nmax,rmlexp,iaddr,mptemp,
     2        mptemp2,rlams,nlams,nfourier,nphysical,
     1        nthmax,nexptot,nexptotp,mexpf1,mexpf2,mexpp1,
     2        mexpp2,rdplus,mexppall(1,1),mexppall(1,2),
     3        mexppall(1,3),mexppall(1,4),mexppall(1,5),mexppall(1,6),
     5        mexppall(1,7),mexppall(1,8),mexppall(1,9),mexppall(1,10),
     6        mexppall(1,11),mexppall(1,12),mexppall(1,13),
     7        mexppall(1,14),mexppall(1,15),mexppall(1,16),xshift,
     8        yshift,zshift,fexpe,fexpo,rlsc)
   
              call processea(ibox,scales(ilev+1),lexp,centers,nboxes,
     1        neall,eall,ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,ne3,e3,
     2        ne5,e5,ne7,e7,mexppall(1,1),mexppall(1,2),mexppall(1,3),
     2        mexppall(1,4),mexppall(1,5),mexppall(1,6),
     2        mexppall(1,7),mexppall(1,8),xshift,yshift,zshift,nexptotp)
            
              call processwe(ibox,scales(ilev+1),lexp,centers,nboxes,
     1        nwall,wall,nw2468,w2468,nw24,w24,nw68,w68,nw2,w2,nw4,w4,
     2        nw6,w6,nw8,w8,mexppall(1,9),mexppall(1,10),mexppall(1,11),
     2        mexppall(1,12),mexppall(1,13),mexppall(1,14),
     2        mexppall(1,15),mexppall(1,16),xshift,yshift,zshift,
     3        nexptotp)
            endif
         enddo
      enddo

      do ilev=2,nlevels
         do ibox=laddr(1,ilev),laddr(2,ilev)
            istart = itree(ipointer(16)+ibox-1)
            iend = itree(ipointer(17)+ibox-1)
            npts = iend-istart+1

            istart = itree(ipointer(18)+ibox-1)
            iend = itree(ipointer(19)+ibox-1)

            npts = npts + iend-istart+1

            istart = itree(ipointer(20)+ibox-1)
            iend = itree(ipointer(21)+ibox-1)

            npts = npts + iend-istart+1
            if(npts.gt.0) then
c              process up down plane waves            
               call phystof(mexpf1,nlams,rlams,nfourier,nphysical,
     1                     nthmax,lexp(1,ibox,1),fexpback)
 
               call phystof(mexpf2,nlams,rlams,nfourier,nphysical,
     1                     nthmax,lexp(1,ibox,2),fexpback)

               call exptolocal(mptemp,nterms(ilev),rlams,whts,nlams,
     1              nfourier,nthmax,nexptot,mexpf1,mexpf2,scales(ilev))

               do ii = -nterms(ilev),nterms(ilev)
                  do jj = 0,nterms(ilev)
                     mptemp(jj+1+(ii+nterms(ilev))*(nterms(ilev)+1))
     1               =mptemp(jj+1+(ii+nterms(ilev))*(nterms(ilev)+1))
     1                  /scales(ilev)**(2*jj)
                  enddo
               enddo

               call l3dadd(mptemp,rmlexp(iaddr(2,ibox)),nterms(ilev))
c              Process north south plane waves  
               call phystof(mexpf1,nlams,rlams,nfourier,nphysical,
     1                      nthmax,lexp(1,ibox,3),fexpback)
               call phystof(mexpf2,nlams,rlams,nfourier,nphysical,
     1                      nthmax,lexp(1,ibox,4),fexpback)

               call exptolocal(mptemp,nterms(ilev),rlams,whts,nlams,
     1              nfourier,nthmax,nexptot,mexpf1,mexpf2,scales(ilev))

               call rotytoz(nterms(ilev),mptemp,mptemp2,rdplus)
               do ii = -nterms(ilev),nterms(ilev)
                  do jj = 0,nterms(ilev)
                     mptemp2(jj+1+(ii+nterms(ilev))*(nterms(ilev)+1))= 
     1               mptemp2(jj+1+(ii+nterms(ilev))*(nterms(ilev)+1))/
     1               scales(ilev)**(2*jj)
                  enddo
               enddo

               call l3dadd(mptemp2,rmlexp(iaddr(2,ibox)),nterms(ilev))
c              Process the east west plane waves  
               call phystof(mexpf1,nlams,rlams,nfourier,nphysical,
     1                      nthmax,lexp(1,ibox,5),fexpback)
               call phystof(mexpf2,nlams,rlams,nfourier,nphysical,
     1                      nthmax,lexp(1,ibox,6),fexpback)

               call exptolocal(mptemp,nterms(ilev),rlams,whts,nlams,
     1              nfourier,nthmax,nexptot,mexpf1,mexpf2,scales(ilev))

               call rotztox(nterms(ilev),mptemp,mptemp2,rdminus)
               do ii = -nterms(ilev),nterms(ilev)
                  do jj = 0,nterms(ilev)
                     mptemp2(jj+1+(ii+nterms(ilev))*(nterms(ilev)+1))= 
     1               mptemp2(jj+1+(ii+nterms(ilev))*(nterms(ilev)+1))/
     1               scales(ilev)**(2*jj)
                  enddo
               enddo

               call l3dadd(mptemp2,rmlexp(iaddr(2,ibox)),nterms(ilev))
            endif
         enddo
      enddo

      call cpu_time(time2)
cccc      time2 = second()
      timeinfo(3) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (split loc) ===*',i,0)

      call cpu_time(time1)
cccc      time1 = second()
      do ilev = 0,nlevels-1
         do ibox = laddr(1,ilev),laddr(2,ilev)
            do i=1,8
               jbox = itree(ipointer(4)+8*(ibox-1)+i-1)
               if(jbox.gt.0) then
                  call l3dzero(mptemp,nterms(ilev+1))

                  call l3dloclocquadu_add(scales(ilev),centers(1,ibox),
     1            rmlexp(iaddr(2,ibox)),nterms(ilev),
     2            scales(ilev+1),centers(1,jbox),mptemp,nterms(ilev+1),
     3            nterms(ilev+1),ier)

                  call l3dadd(mptemp,rmlexp(iaddr(2,jbox)),
     1            nterms(ilev+1))
               endif
            enddo
         enddo
      enddo
      call cpu_time(time2)
cccc      time2 = second()
      timeinfo(4) = time2-time1

c
      if(ifprint.ge.1)
     $    call prinf('=== step 5 (eval lo) ===*',i,0)

c     ... step 5, evaluate all local expansions
      call cpu_time(time1)
cccc      time1 = second()
      do ilev = 1,nlevels

         do ibox = laddr(1,ilev),laddr(2,ilev)

            istart = itree(ipointer(18)+ibox-1)
            iend = itree(ipointer(19)+ibox-1)

            if(ilev.eq.nlevels) iend = itree(ipointer(21)+ibox-1)

            do i=istart,iend
               pottmp = 0.0d0
               gradtmp(1) = 0.0d0
               gradtmp(2) = 0.0d0
               gradtmp(3) = 0.0d0

               call l3dtaevalall_trunc(scales(ilev),centers(1,ibox),
     1         rmlexp(iaddr(2,ibox)),nterms(ilev),
     2         nterms(nlevels),          
     2         targetsort(1,i),1,ifpottarg,pottmp,ifgradtarg,gradtmp,
     3         wlege,nlege,ier)

               if(ifpottarg.eq.1) pottarg(i) = -pottmp/(4*pi)
               if(ifgradtarg.eq.1) then 
                  gradtarg(1,i) = -gradtmp(1)/(4*pi)
                  gradtarg(2,i) = -gradtmp(2)/(4*pi)
                  gradtarg(3,i) = -gradtmp(3)/(4*pi)
               endif
            enddo
         enddo
      enddo


      call cpu_time(time2)
cccc      time2 = second()
      timeinfo(5) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 6 (direct) =====*',i,0)
      call cpu_time(time1)
cccc      time1=second()
      do ilev = 0,nlevels
      do ibox = laddr(1,ilev),laddr(2,ilev)

         istart = itree(ipointer(14)+ibox-1)
         iend = itree(ipointer(15)+ibox-1)
         nlist1 = itree(ipointer(5)+ibox-1)
         do i =1,nlist1
            jbox = itree(ipointer(6)+mnlist1*(ibox-1)+i-1)

            jstart = itree(ipointer(18)+jbox-1)
            jend = itree(ipointer(19)+jbox-1)
 
            if(ilev.eq.nlevels) jend = itree(ipointer(21)+jbox-1)

            call tfmm3dparttarg_direct_newtree(istart,iend,jstart,jend,
     1      sourcesort,ifcharge,chargesort,ifdipole,dipstrsort,
     2      dipvecsort,targetsort,stdevsort,stdev_gradsort,ifpottarg,
     3      pottarg,ifgradtarg,gradtarg)
         enddo   
      enddo
      enddo
      call cpu_time(time2)
cccc      time2 = second()
      timeinfo(6) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,6)
      d = 0
      do i = 1,6
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)

      return
      end
c
      subroutine tfmm3dparttarg_direct_newtree(istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     targ,stdev,stdev_grad,ifpottarg,pottarg,ifgradtarg,gradtarg)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the computed velocities,
c     gradients 
c
c     INPUT arguments
c-------------------------------------------------------------------
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: real *8(3,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                 dipole strengths at the source locations
c
c     dipvec      in: real *8(3,ns)
c
c     targ        in: real *8(3,nt)
c                 target locations
c
c     stdev:      real *8 (nt): smoothing kernel parameter 
c     stdev_grad: real *8 (nt): grad of stdev 
c
c     ifpottarg    in: Integer
c                  Flag for computing the potential. The potential
c                  will be updated if ifpottarg = 1
c
c     ifgradtarg  in: Integer
c                  Flag for computing the gradient.
c                  The gradient will be updated if
c                  ifgradtarg == 1
c
c------------------------------------------------------------
c     OUTPUT
c
c   pottarg : potential at the targets
c   gradtarg: gradient at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        real *8 source(3,*)
        complex *16 charge(*),dipstr(*)
        real *8 dipvec(3,*)
        real *8 stdev(*)
        real *8 stdev_grad(3,*)
        real *8 qwt

        integer ifpottarg,ifgradtarg
        real *8 targ(3,*)
c
        complex *16 pottarg(*)
        complex *16 gradtarg(3,*)

ccc        complex *16 pottmp,gradtmp(3)
        real *8 pottmp,gradtmp(3)

c
        ns = iend - istart + 1
        do j=jstart,jend
              do i=istart,iend
                 if((targ(1,j)-source(1,i))**2 +
     1              (targ(2,j)-source(2,i))**2.gt.1.0d-28) then
                    if(ifcharge.eq.1) then
                       call lpotfld3d(ifgradtarg,
     1                 source(1,i),
     1                 charge(i),targ(1,j),pottmp,gradtmp)
                       if(ifpottarg.eq.1) pottarg(j) = pottarg(j)+
     1                                    pottmp
                       if(ifgradtarg.eq.1) then
                          gradtarg(1,j) = gradtarg(1,j)+gradtmp(1)
                          gradtarg(2,j) = gradtarg(2,j)+gradtmp(2)
                          gradtarg(3,j) = gradtarg(3,j)+gradtmp(3)
                       endif
                    endif
                    if(ifdipole.eq.1) then
ccc                       call lpotfld3d_dp(ifgradtarg,
ccc     1                  source(1,i),dipstr(i),dipvec(1,i),targ(1,j),
ccc     2                  pottmp,gradtmp)
                       qwt =  dreal(dipstr(i))
                       call tpotfld3d_dp(ifgradtarg,
     1                  source(1,i),qwt,dipvec(1,i),targ(1,j),
     2                  stdev(j),stdev_grad(1,j),pottmp,gradtmp)
                       if(ifpottarg.eq.1) pottarg(j) = pottarg(j)+
     1                                    pottmp
                       if(ifgradtarg.eq.1) then
                          gradtarg(1,j) = gradtarg(1,j)+gradtmp(1)
                          gradtarg(2,j) = gradtarg(2,j)+gradtmp(2)
                          gradtarg(3,j) = gradtarg(3,j)+gradtmp(3)
                       endif
                    endif
                 endif
              enddo
        enddo
c
        return
        end
c------------------------------------------------------------------    
c
        subroutine l2dadd2(mpole,ntbig,mpole2,nterms)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:ntbig)
        complex *16 mpole2(0:nterms)
c       
        do n=0,nterms
        mpole2(n)=mpole2(n)+mpole(n)
        enddo
c
        return
        end
c-----------------------------------------------------        
 
      subroutine l3dreorderint(n,arr,isort,arrsort)
c    This subroutine reorders the array of integers
c    to the array arrsort via the mapping isort.
c    The ordering is reset to arrsort(i) = arr(isort(i))
c
c    INPUT arguments
c    n       in:Integer
c            Number of elements in the array
c
c    arr     in:Integer(n)
c            The array to be sorted
c
c    isort   in:Integer(n)
c            Mapping used to sort the array
c
c    OUTPUT arguments
c    arrsort  in:Integer(n)
c             The sorted array
c -------------------------------------------------------
      implicit none
      integer n,i,arr(1),arrsort(1),isort(1)

      do i=1,n
         arrsort(i) = arr(isort(i))
      enddo

      return
      end
c----------------------------------------------------------      
      subroutine l3dmpalloc_newtree(laddr,iaddr,nlevels,lmptot,
     1                          nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
c  
c     Input arguments
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array provinding access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels)
      integer iaddr(2,1), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr

      iptr = 1
      do i = 0,nlevels
         do ibox = laddr(1,i),laddr(2,i)
c            Allocate memory for the multipole expansion         
             iaddr(1,ibox) = iptr
             iptr = iptr + (2*nterms(i)+1)*2*(nterms(i)+1)

c            Allocate memory for the local expansion
             iaddr(2,ibox) = iptr
             iptr = iptr + (2*nterms(i)+1)*2*(nterms(i)+1)
         enddo
      enddo
      lmptot = iptr

      return
      end
c----------------------------------------------------------------     
c
      subroutine l3dreorder_newtree(ns,source,ifcharge,charge,
     1           isource,ifdipole,dipstr,dipvec,sourcesort,
     2           chargesort,dipstrsort,dipvecsort)

c     This subroutine sorts the sources, charges and dipole
c     vectors into tree ordered arrays given by the mapping
c     isource
c     arrsort(i) = arr(isource(i))
c
c     Input arguments
c     ns       in:Integer
c              Number of sources
c
c     source   in: real *8(3,ns)
c              x and y co-ordinates of the sources
c
c     ifcharge in: Integer
c              flag for sorting charges
c
c     charges  in: complex *16(ns)
c              charge strengths at the source locations
c
c     ifdipole in: Integer
c              flag for sorting dipoles
c
c     dipstr   in: complex *16(ns)
c              dipole strength at the source locations
c
c     dipvec   in: real *8(3,ns)
c              dipole vector orientation at the source locations
c 
c     isource  in: Integer(ns)
c              mapping for sorting the relevant quantities
c
c
c     OUTPUT arguments
c     sourcesort   out: real *8(3,ns)
c                  Tree sorted sources
c
c     chargesort   out: complex *16(ns)
c                  Tree sorted charge strengths
c
c     dipstrsort   out: complex *16(ns)
c                  Tree sorted dipole strengths
c
c     dipvecsort   out: real *8(3,ns)
c                  Tree sorted dipole orientation vectors
c
c---------------------------------------------------------------
      implicit none
      real *8 source(3,*), sourcesort(3,*)
      integer isource(*),ns,i
      integer ifcharge,ifdipole
      complex *16 charge(*),chargesort(*)
      complex *16 dipstr(*),dipstrsort(*)
      real *8 dipvec(3,*), dipvecsort(3,*)

      do i=1,ns
         sourcesort(1,i) = source(1,isource(i))
         sourcesort(2,i) = source(2,isource(i))
         sourcesort(3,i) = source(3,isource(i))
         if(ifcharge.eq.1) chargesort(i) = charge(isource(i))
         if(ifdipole.eq.1) then
            dipstrsort(i) = dipstr(isource(i))
            dipvecsort(1,i) = dipvec(1,isource(i))
            dipvecsort(2,i) = dipvec(2,isource(i))
            dipvecsort(3,i) = dipvec(3,isource(i))
         endif
      enddo

      return
      end
c------------------------------------------------------------------

      subroutine l3dreordertarg_newtree(nt,targ,trads,stdev,stdev_grad,
     1   itarg,targsort,tradsort,stdevsort,stdev_gradsort)
c     This subroutine tree sorts the target locations
c     given by the mapping itarg
c
c     arrsort(i) = arr(itarg(i))
c
c     INPUT arguments
c     nt       in: Integer
c              number of targets
c
c     targ     in: real *8(3,nt)
c              x, y, and z co-ordinates of targets
c
c     trads    in: real *8 (nt)
c              extent of targets
c
c     stdev      in: real *8 (nt)
c                smoothing kernel parameter for targets
c
c     stdevgrad      in: real *8 (nt)
c              grad of smoothing kernel parameter for targets
c
c     itarg    in: Integer(nt)
c              mapping for tree sorting the targets
c
c     OUTPUT
c     targsort   out: real *8(3,nt)
c                tree sorted target locations
c
c     tradsort   out: real *8 (nt)
c                tree sorted extent of targets
c
c     stdevsort    out: real *8 (nt)
c                tree sorted stdev
c     stdev_gradsort    out: real *8 (nt)
c                tree sorted stdev
c------------------------------------------------------------------

      implicit none
      integer nt,i,itarg(*)
      real *8 targ(3,*),targsort(3,*)
      real *8 trads(*),tradsort(*)
      real *8 stdev(*), stdevsort(*)
      real *8 stdev_grad(3,*), stdev_gradsort(3,*)

      do i=1,nt
         targsort(1,i) = targ(1,itarg(i))
         targsort(2,i) = targ(2,itarg(i))
         targsort(3,i) = targ(3,itarg(i))

         tradsort(i) = trads(itarg(i))
         stdevsort(i) = stdev(itarg(i))
         stdev_gradsort(1,i) = stdev_grad(1,itarg(i))
         stdev_gradsort(2,i) = stdev_grad(2,itarg(i))
         stdev_gradsort(3,i) = stdev_grad(3,itarg(i))
      enddo

      return
      end
c-----------------------------------------------------------------      

      subroutine l3dtexpsort_newtree(ns,isource,ntj,jsort,jexps,
     1                   stdevsort,stdev)

c     This subroutine reorders the tree sorted expansions and
c     expansion scaling parameters to the user prescribed ordering.
c     The mapping isource is the mapping from the user prescribed
c     ordering to the tree sorted ordering.
c     arrsort(i) = arr(isource(i))
c
c     INPUT arguments
c     ns        in: Integer
c               Number of sources
c
c     isource   in: Integer (ns)
c               Tree sort map
c
c     ntj       in: integer
c               number of terms in expansion
c
c     jsort     in: complex *16 (0:ntj,ns)
c               Tree sorted expansions
c
c    stdevsort    in: real *8 (ns)
c               Tree sorted smoothing kernel parameters
c
c    OUTPUT arguments
c    jexps      out: complex *16(0:ntj,ns)
c               Expansions in the user prescribed ordering
c
c    stdev        out: real *8(ns)
c               smoothing kernel parameters in the user
c               prescribed ordering
c----------------------------------------------------------

      implicit none
      integer i,j,k,ns,ntj
      integer isource(*)
      complex *16 jsort(0:ntj,-ntj:ntj,*),jexps(0:ntj,-ntj:ntj,*)
      real *8 stdevsort(*),stdev(*)

      do i=1,ns
         stdev(isource(i)) = stdevsort(i)
         do j=0,ntj
            do k=-ntj,ntj
               jexps(j,k,isource(i)) = jsort(j,k,i)
            enddo
         enddo
      enddo

      return
      end
c------------------------------------------------------------
      subroutine l3dpsort_newtree(nt,itarg,potsort,pottarg)
c     This subroutine reorders the tree sorted target
c     potentials to the user prescribed ordering. itarg
c     is a mapping from the user prescribed ordering to
c     the tree sorted ordering for targets
c     arrsort(i) = arr(itarg(i))
c
c     INPUT arguments
c     nt          in:Integer
c                 Number of targets
c
c     itarg       in:Integer(nt)
c                 Tree sorting map
c
c     potsort     in: complex *16 (nt)
c                 Tree sorted potential at target locations
c
c     OUTPUT
c     pottarg     out:complex *16(nt)
c                 Potential at targets in user prescribed ordering
c-----------------------------------------------------------------
      implicit none
      integer nt,i,itarg(*)
      complex *16 potsort(*),pottarg(*)

      do i=1,nt
         pottarg(itarg(i)) = potsort(i)
      enddo

      return
      end
c----------------------------------------------------------------
      subroutine l3dfsort_newtree(nt,itarg,gradsort,gradtarg)
c     This subroutine reorders the tree sorted target
c     gradients to the user prescribed ordering. itarg
c     is a mapping from the user prescribed ordering to
c     the tree sorted ordering for targets
c     arrsort(i) = arr(itarg(i))
c
c     INPUT arguments
c     nt          in:Integer
c                 Number of targets
c
c     itarg       in:Integer(nt)
c                 Tree sorting map
c
c     gradsort     in: complex *16 (3,nt)
c                 Tree sorted gradient at target locations
c
c     OUTPUT
c     gradtarg     out:complex *16(3,nt)
c                 gradients at targets in user prescribed ordering
c-----------------------------------------------------------------
      implicit none
      integer nt,i,itarg(*)
      complex *16 gradsort(3,*),gradtarg(3,*)

      do i=1,nt
         gradtarg(1,itarg(i)) = gradsort(1,i)
         gradtarg(2,itarg(i)) = gradsort(2,i)
         gradtarg(3,itarg(i)) = gradsort(3,i)
      enddo

      return
      end
c----------------------------------------------------------------
