c
c    Uniform PTree - generate quad tree refined to finest level pruning
c    empty boxes. Features include 
c
c     radsrc parameter (preventing sources with large radius 
c                         from being pushed to the finest level)
c     radexp parameter (preventing exp centers with large radius 
c                         from being pushed to the finest level)
c     ndiv parameter (usual refinement parameter extending tree
c                         depth until all leaf nodes hanve fewer than
c                         ndiv particles)
c     idivflag       (allows use of ndiv parameter for sources OR
c                         targets OR sources+targets OR sources+targets
c                         + expansion centers)
c   
c     isep           (separation parameter: isep = 1 or 2
c                     based on whether to include neighbors
c                     and self in list 1 or two neighbors
c                     and self in list 2)
c                     
c
c     mhung is the max number of hung chunks in a given box
c
c     This tree code MUST be used along with its
c     memomry code mkptreemem located at the end of this file.
c     The memory code precomputes the number of levels, the number
c     of boxes, the max number of hung sources and the length
c     of the tree. 
c
c     subroutines in this file
c     mkptree - forms the tree. Note this must be preceeded by 
c     a call to mkptreemem
c
c     mkptreemem - determination of amount of memory needed
c     for forming tree 
c
c     getaq - perform area query for a given tree structure.
c     Note this must be preceeded by a call to getaqmem
c
c     getaqmem - determine amount of memory needed for 
c     area query structures
c
c
c
      subroutine mkptree(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,isep,
     $                   mhung,nlevels,nboxes,
     $                   centers,boxsize,itree,ltree,
     $                   ipointer)
      implicit none
      integer ier,ns,nt,nexpc,idivflag,ndiv,isep,mhung
      integer nlevels,lcenters,ltree
      integer nlmax,nbmax,nboxes, nlevtmp,nbtmp, mhungtmp
      integer itree(ltree)
      integer i,j
      integer, allocatable :: laddr(:,:)
      integer, allocatable :: iparenttemp(:)
      integer, allocatable :: nchild(:)
      integer, allocatable :: ichildtemp(:,:)
      integer, allocatable :: nnbor(:)
      integer, allocatable :: nbors(:,:)
      integer, allocatable :: nlist2(:)
      integer, allocatable :: list2(:,:)
      integer, allocatable :: isourcetemp(:)
      integer, allocatable :: itargettemp(:)
      integer, allocatable :: iexpctemp(:)
      integer, allocatable :: ihsfirsttemp(:)
      integer, allocatable :: ihslasttemp(:)
      integer, allocatable :: isfirsttemp(:)
      integer, allocatable :: islasttemp(:)
      integer, allocatable :: itfirsttemp(:)
      integer, allocatable :: itlasttemp(:)
      integer, allocatable :: ihefirsttemp(:)
      integer, allocatable :: ihelasttemp(:)
      integer, allocatable :: iefirsttemp(:)
      integer, allocatable :: ielasttemp(:)
      integer, allocatable :: nhungsrc(:)
      integer, allocatable :: nhunglistsrc(:)
      integer, allocatable :: nhungexp(:)
      integer, allocatable :: ihunglistsrc(:,:)

      integer ipointer(26)
      real *8 boxsize(0:nlevels)
      real *8 src(3,ns),radsrc(ns)
      real *8 trg(3,nt)
      real *8 centers(3,nboxes)
      real *8 expc(3,nexpc)
      real *8 radexp(nexpc)
c
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez
      integer ictr,ih,irefine,is,ie
      integer nss,nee

      integer mnbors,mnlist2
      integer nupm,ndownm,nnorthm,nsouthm,neastm,nwestm
      integer, allocatable :: uplist(:,:)
      integer, allocatable :: downlist(:,:)
      integer, allocatable :: northlist(:,:)
      integer, allocatable :: southlist(:,:)
      integer, allocatable :: eastlist(:,:)
      integer, allocatable :: westlist(:,:)
      integer, allocatable :: nup(:)
      integer, allocatable :: ndown(:)
      integer, allocatable :: nnorth(:)
      integer, allocatable :: nsouth(:)
      integer, allocatable :: neast(:)
      integer, allocatable :: nwest(:)

c
c     INPUT:
c
c     src           source locations        
c     ns            number of sources 
c     rads          source radii (determines deepest level that
c                   the source can reach) 
c     trg           target locations        
c     nt            number of targets
c
c     expc          expansion center locations
c     nexpc         number of expansion centers
c
c     idivflag      0 => divide on sources
c                   1 => divide on targets
c                   2 => divide on sources+targets
c                   3 => divide on sources+targets+expansion centers
c 
c     ndiv          refinement criterion - extend tree until all
c                   nodes at finest level have less than ndiv 
c                   source/targets/sources+targets depending on
c                   idivflag
c
c     isep          separation parameter for determining 
c                   nearest neighbors
c
c     ltree         length of ltree
c     nlevels       number of levels determined using the memory code
c     nboxes        number of boxes determined using mem code
c                   refer to mkptreemem.
c
c
c     OUTPUT:
c     centers       array of box centers
c     boxsize       box dimensions at all levels
c
c     itree         tree array
c     
c     iladdr = 1
c     itree(iladdr) <-> laddr
c                   indexing array providing access to boxes at
c                   each level. 
c                   the first box on level i is laddr(1,i)
c                   the last box on level i is laddr(2,i)
c
c     iparent = iladdr + 2*(nlevels+1)
c     itree(iparent) <-> parent
c                   parent of box (set to -1 for root node)
c
c     inchild = iparent + nboxes
c     itree(inchild) <-> nchild
c                   number of children for each box 
c
c     ichild = inchild + nboxes
c     itree(ichild) <-> child
c                   only first nchild entries are defined 
c                   (others set to -1). The array should be viewed 
c                   as dimensioned (8,nboxes)
c
c     innbor = ichild + 8*nboxes
c     itree(innbor) <-> nnbor
c                   number of neighbors
c
c     inbors = innbor + nboxes
c     itree(inbors) <-> nbors
c                   only first nnbor entries are set 
c                   (others set to -1). The array should be viewed 
c                   as dimensioned (27,nboxes) if isep = 1
c                   and (125, nboxes) if isep = 2
c                   mnnbors = 27/125 if isep=1/2
c
c     inlist2 = inbors + mnnbors*nboxes
c     itree(inlist2) <-> nlist2
c                   number of elements in the interaction list(list2)
c
c     ilist2 = inlist2 + nboxes
c     itree(ilist2) <-> list2
c                   only first nlist2 entries are set 
c                   (others set to -1). The array should be viewed 
c                   as dimensioned (189,nboxes) if isep = 1
c                   and (875,nboxes) if isep = 2
c                   mnlist2 = 189/875 if isep=1/2
c
c     isource = ilist2 + mnlist2*nboxes
c     itree(isource) <-> isource
c                   tree-ordered array of sources
c
c     itarget = isource + ns
c     itree(itarget) <-> itarget
c                   tree-ordered array of targets
c
c     iexpc = itarget + nt
c     itree(iexpc) <-> iexpc
c                      tree-ordered array of expansion centers
c
c     ihsfirst = itarget + nt
c     itree(ihsfirst) <-> ihsfirst
c                   ihsfirst(j) = location in isource of first hung 
c                                source for box j 
c
c     ihslast = ihsfirst + nboxes
c     itree(ihslast) <-> ihslast
c                   ihslast(j)  = location in isource of last hung
c                                source for box j 
c
c     isfirst = ihlast + nboxes
c     itree(isfirst) <-> isfirst
c                   isfirst(j) = location in isource of first 
c                                source for box j 
c
c     islast = isfirst + nboxes
c     itree(islast) <-> islast
c                   islast(j)  = location in isource of last 
c                                source for box j 
c
c     itfirst = islast + nboxes
c     itree(itfirst) <-> itfirst
c                   itfirst(j) = location in itarget of first 
c                                target for box j 
c
c     itlast = itfirst + nboxes
c     itree(itlast) <-> itlast
c                   itlast(j)  = location in itarget of last 
c                                target for box j
c
c     ihefirst = itlast + nboxes
c     itree(ihefirst) <-> ihefirst
c                   ihefirst(j) = location in iexpc of first hung 
c                                 expansion center for box j 
c
c     ihelast = ihefirst + nboxes
c     itree(ihelast) <-> ihelast
c                   ihelast(j)  = location in iexpc of last hung
c                                expansion center for box j 
c
c     iefirst = ihelast + nboxes
c     itree(iefirst) <-> iefirst
c                    iefirst(j) = location in iexpc of first target
c                                 in box j
c
c     ielast = ielast + nboxes
c     itree(ielast) <-> ielast
c                       ielast(j) = location in iexpc of last target
c                                   in box j
c
c     inhungsrc = ielast + nboxes
c     itree(inhungsrc) <-> inhungsrc
c                inhungsrc(j)  = number of hung chunks in box j
c
c     inhungexp = inhungsrc + nboxes
c     itree(inhungexp) <-> inhungexp
c                inhungexp(j) = number of hung expansion centers
c                               in box j
c
c     nhunglistsrc = inhungexp + nboxes
c     itree(inhunglistsrc) <-> inhunglistsrc
c           inhunglistsrc(j) = Total number of hung sources
c                              relevant for box j
c
c     ihunglistsrc = inhungexp + nboxes
c     itree(ihunglistsrc) <-> ihunglistsrc
c                   ihunglistsrc(m,j) = src id of  mth hung src
c                   relevant to box j. A hung src is relevant to
c                   a box if the box is a descendant of the
c                   neighbor of the box in which the src is hung
c
c     ltree = ihunglistsrc + mhung*nboxes
c
c     ipointers is the collection of pointers
c     ipointer(1) = iladdr
c     ipointer(2) = iparent = ipointer(1) + 2*(nlevels+1)
c     ipointer(3) = inchild = ipointer(2) + nboxes
c     ipointer(4) = ichild = ipointer(3) + nboxes
c     ipointer(5) = nnbor = ipointer(4) + 8*nboxes
c     ipointer(6) = nbors = ipointer(5) + nboxes
c     ipointer(7) = nlist2 = ipointer(6) + mnbors*nboxes
c     ipointer(8) = list2 = ipointer(7) + nboxes
c     ipointer(9) = isource = ipointer(8) + mnlist2*nboxes
c     ipointer(10) = itarget = ipointer(9) + ns
c     ipointer(11) = iexpc = ipointer(10) + nt
c     ipointer(12) = ihsfirst = ipointer(11) + nexpc
c     ipointer(13) = ihslast = ipointer(12) + nboxes
c     ipointer(14) = isfirst = ipointer(13) + nboxes
c     ipointer(15) = islast = ipointer(14) + nboxes
c     ipointer(16) = itfirst = ipointer(15) + nboxes
c     ipointer(17) = itlast = ipointer(16) + nboxes
c     ipointer(18) = ihefirst = ipointer(17) + nboxes
c     ipointer(19) = ihelast = ipointer(18) + nboxes
c     ipointer(20) = iefirst = ipointer(19) + nboxes
c     ipointer(21) = ielast = ipointer(20) + nboxes
c     ipointer(22) = nhungsrc = ipointer(21) + nboxes
c     ipointer(23) = nhungexp = ipointer(22) + nboxes
c     ipointer(24) = nhunglistsrc = ipointer(23) + nboxes
c     ipointer(25) = hunglistsrc = ipointer(24) + nboxes
c     ipointer(26) = lentree = ipointer(25) + mhung*nboxes
c
c     mhung is the max number of hung sources 
c
c---------------------------------------------------------------------
c
c     Other notes:
c
c     refinement criterion w.r.t rads is:
c
c     hang if (rads .geq. boxsize)
c
c
      lcenters = nboxes

      allocate(laddr(2,0:nlevels))
      allocate(iparenttemp(nboxes))
      allocate(nchild(nboxes))
      allocate(ichildtemp(8,nboxes))
      allocate(isourcetemp(ns))
      allocate(itargettemp(nt))
      allocate(iexpctemp(nexpc))
      allocate(ihsfirsttemp(nboxes))
      allocate(ihslasttemp(nboxes))
      allocate(isfirsttemp(nboxes))
      allocate(islasttemp(nboxes))
      allocate(itfirsttemp(nboxes))
      allocate(itlasttemp(nboxes))
      allocate(ihefirsttemp(nboxes))
      allocate(ihelasttemp(nboxes))
      allocate(iefirsttemp(nboxes))
      allocate(ielasttemp(nboxes))
      allocate(nhungsrc(nboxes))
      allocate(nhungexp(nboxes))

      if(isep.eq.1) then
         mnbors = 27
         mnlist2 = 27*7
      endif

      if(isep.eq.2) then
         mnbors = 125
         mnlist2 = 125*7
      endif

      if(isep.ne.1.and.isep.ne.2) then
         ier = 4
         call prinf('Error tree not allocated, isep.ne.1,2*',ier,0)
         return
      endif

c     Step 1: Find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo

      do i=1,nexpc
        if(expc(1,i) .lt. xmin) xmin=expc(1,i)
        if(expc(1,i) .gt. xmax) xmax=expc(1,i)
        if(expc(2,i) .lt. ymin) ymin=expc(2,i)
        if(expc(2,i) .gt. ymax) ymax=expc(2,i)
        if(expc(3,i) .lt. zmin) zmin=expc(3,i)
        if(expc(3,i) .gt. zmax) zmax=expc(3,i)
      enddo
      boxsize(0)=xmax-xmin
      sizey=ymax-ymin
      sizez=zmax-zmin
      if(sizey .gt. boxsize(0)) boxsize(0)=sizey
      if(sizez .gt. boxsize(0)) boxsize(0)=sizez
c
c     initialize arrays at level 0
c
      centers(1,1)=(xmin+xmax)/2
      centers(2,1)=(ymin+ymax)/2
      centers(3,1)=(zmin+zmax)/2
      laddr(1,0) = 1
      laddr(2,0) = 1
      iparenttemp(1) = -1
      isfirsttemp(1) = 1
      nhungsrc(1) = 0
      nhungexp(1) = 0
c
c     count number of hung sources
c     and hang up "big" sources
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) nhungsrc(1)=nhungsrc(1)+1 
      enddo
      isfirsttemp(1) = nhungsrc(1)+1
      islasttemp(1) = ns
      if (nhungsrc(1) .gt. 0) then 
         ihsfirsttemp(1) = 1
         ihslasttemp(1) = nhungsrc(1)
      else
         ihsfirsttemp(1) = 0
         ihslasttemp(1) = -1
      endif

c     Count number of hung expansion centers      
c     and hang up "big" expansion centers
      do i=1,nexpc
         if (radexp(i).gt.boxsize(0)) nhungexp(1)=nhungexp(1)+1
      enddo
      iefirsttemp(1) = nhungexp(1)+1
      ielasttemp(1) = nexpc
      if (nhungexp(1).gt.0) then
          ihefirsttemp(1) = 1
          ihelasttemp(1) = nhungexp(1)
      else
         ihefirsttemp(1) = 0
         ihelasttemp(1) = -1
      endif

c
c     reorder isourcetemp to put hung sources in beginning
      ih = 0
      is = nhungsrc(1)
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) then
            ih = ih+1
            isourcetemp(ih) = i
         else
            is = is+1
            isourcetemp(is) = i
         endif
      enddo
c     reorder iexptemp to put hung expansion centers in beginning
      ih = 0
      ie = nhungexp(1)
      do i= 1,nexpc
         if(radexp(i).gt.boxsize(0)) then
            ih = ih+1
            iexpctemp(ih) = i
         else
            ie = ie+1
            iexpctemp(ie) = i
         endif
      enddo

c     initialize itargettemp 
      do i = 1,nt
         itargettemp(i) = i
      enddo
      itfirsttemp(1) = 1
      itlasttemp(1) = nt

      nchild(1) = 0
      ichildtemp(1,1) = -1
      ichildtemp(2,1) = -1
      ichildtemp(3,1) = -1
      ichildtemp(4,1) = -1
      ichildtemp(5,1) = -1
      ichildtemp(6,1) = -1
      ichildtemp(7,1) = -1
      ichildtemp(8,1) = -1

c
      irefine = 0
      nss = ns - nhungsrc(1)
      nee = nexpc - nhungexp(1)
      if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
      if ((idivflag .eq.1).and.(nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.2).and.(nss+nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv)) irefine=1

c     Reset nlevels, nboxes      
      nlevtmp = 0
      nbtmp = 1
      do i = 1,nlevels
         if (irefine.eq.1) then
            call subdivide(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevtmp,nbtmp,
     $                   centers,lcenters,boxsize,nboxes,nlevels,
     $                   laddr,iparenttemp,nchild,ichildtemp,
     $                   isourcetemp,itargettemp,iexpctemp,
     $                   ihsfirsttemp,ihslasttemp,
     $                   isfirsttemp,islasttemp,
     $                   itfirsttemp,itlasttemp,
     $                   ihefirsttemp,ihelasttemp,
     $                   iefirsttemp,ielasttemp,nhungsrc,
     $                   nhungexp,irefine)

         else
            exit
         endif
      enddo
      call prinf('At the end of main tree code=*',i,0)
c     Set up computation of list1 and list2      
      allocate(nnbor(nboxes))
      allocate(nbors(mnbors,nboxes))
      allocate(nlist2(nboxes))
      allocate(list2(mnlist2,nboxes))

c     Initialize nnbor,nlist2,nbors,list2 arrays
      do i=1,nboxes
          nnbor(i) = 0
          nlist2(i) = 0
          do j=1,mnbors
              nbors(j,i) = -1
              list2(j,i) = -1
          enddo
          do j=mnbors+1,mnlist2
              list2(j,i) = -1
          enddo
      enddo

c     compute list1 and list2
      call computelists(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,isep,nnbor,nbors,mnbors,nlist2,
     3                   list2,mnlist2)

c     Initialize hung variables
      allocate(nhunglistsrc(nboxes))
      allocate(ihunglistsrc(mhung,nboxes))
      do i=1,nboxes
         nhunglistsrc(i) = 0
         do j=1,mhung
            ihunglistsrc(j,i) = -1
         enddo
      enddo

c     compute hunglist      
      if(mhung.gt.0) then
         call computehunglist(mhung,nlevels,nboxes,laddr,
     1                        ns,isourcetemp,iparenttemp,
     1                        ihsfirsttemp, ihslasttemp,
     2                        nnbor,nbors,mnbors,nhungsrc,
     3                        nhunglistsrc,ihunglistsrc)
      endif

c      call findnudnsew(isep,nlevels,laddr,boxsize,
c     1                 nboxes,centers,mnlist2,
c     2                 nlist2,list2,nupm,ndownm,nnorthm,
c     3                 nsouthm,neastm,nwestm)

c      allocate(uplist(nupm,nboxes))
c      allocate(downlist(ndownm,nboxes))
c      allocate(northlist(nnorthm,nboxes))
c      allocate(southlist(nsouthm,nboxes))
c      allocate(eastlist(neastm,nboxes))
c      allocate(westlist(nwestm,nboxes))

c      allocate(nup(nboxes))
c      allocate(ndown(nboxes))
c      allocate(nnorth(nboxes))
c      allocate(nsouth(nboxes))
c      allocate(neast(nboxes))
c      allocate(nwest(nboxes))

c      call mkpwlists(isep,nlevels,laddr,boxsize,
c     1               nboxes,centers,mnlist2,nlist2,list2,
c     2               nupm,ndownm,nnorthm,nsouthm,
c     3               neastm,nwestm,nup,uplist,ndown,downlist,
c     3               nnorth,northlist,nsouth,southlist,
c     4               neast,eastlist,nwest,westlist)

c      call prinf('nupm=*',nupm,1)
c      call prinf('ndownm=*',ndownm,1)
c      call prinf('neastm=*',neastm,1)
c      call prinf('nwestm=*',nwestm,1)
c      call prinf('nnorthm=*',nnorthm,1)
c      call prinf('nsouthm=*',nsouthm,1)
c      do i=1,nboxes
c         call prinf('ibox=*', i,1)
c         call prinf('nuplist ibox=*',neast(i),1)
c      enddo

c     Move the output to itree
c     Store iladdr - first and last box at level ilev
      ictr = 1
      ipointer(1) = ictr
      do i=0,nlevels
         itree(ictr) = laddr(1,i)
         ictr = ictr + 1
         itree(ictr) = laddr(2,i)
         ictr = ictr + 1
      enddo

c     Store iparent - parent of box i
      ipointer(2)  = ictr
      do i = 1,nboxes
          itree(ictr) = iparenttemp(i)
          ictr = ictr + 1
      enddo

c     Store inchild - number of children of given box
      ipointer(3)  = ictr
      do i = 1,nboxes
         itree(ictr) = nchild(i)
         ictr = ictr + 1
      enddo 

c     Store ichild temp - children of given box i
      ipointer(4)  = ictr
      do i = 1,nboxes
         do j=1,8
             itree(ictr) = ichildtemp(j,i)
             ictr = ictr + 1
         enddo
      enddo

c     Store nnbor - number of elements in nearest neighbors
      ipointer(5)  = ictr
      do i =1,nboxes
         itree(ictr) = nnbor(i)
         ictr = ictr + 1
      enddo

c     Store nbors - nearest neighbors
      ipointer(6)  = ictr
      do i =1,nboxes
         do j=1,mnbors
            itree(ictr) = nbors(j,i)
            ictr = ictr + 1
         enddo
      enddo

c     Store nlist2 - number of elements in the interaction list
      ipointer(7)  = ictr
      do i=1,nboxes
         itree(ictr) = nlist2(i)
         ictr = ictr + 1
      enddo

c     list2 - interaction list     
      ipointer(8)  = ictr
      do i=1,nboxes
         do j = 1,mnlist2
            itree(ictr) = list2(j,i)
            ictr = ictr + 1
         enddo
      enddo

c     isource - tree sorted array of sources
      ipointer(9)  = ictr
      do i =1,ns
         itree(ictr) = isourcetemp(i)
         ictr = ictr + 1
      enddo

c     itarget - tree sorted array of targets
      ipointer(10)  = ictr
      do i = 1,nt
         itree(ictr) = itargettemp(i)
         ictr = ictr + 1
      enddo

      ipointer(11) = ictr
      do i= 1,nexpc
         itree(ictr) = iexpctemp(i)
         ictr = ictr + 1
      enddo

c     ihfirst - first hung chunk in box i
      ipointer(12)  = ictr
      do i=1,nboxes
         itree(ictr)  = ihsfirsttemp(i)
         ictr = ictr + 1
      enddo

c     ihlast - last hung chunk in box i
      ipointer(13)  = ictr
      do i=1,nboxes
         itree(ictr) = ihslasttemp(i)
         ictr = ictr + 1
      enddo

c     isfirst - first source in box i
      ipointer(14)  = ictr
      do i=1,nboxes
         itree(ictr)  = isfirsttemp(i)
         ictr = ictr + 1
      enddo

c     islast - last source in box i
      ipointer(15)  = ictr
      do i=1,nboxes
         itree(ictr) = islasttemp(i)
         ictr = ictr + 1
      enddo

c     itfirst - first target in box i
      ipointer(16)  = ictr
      do i=1,nboxes
         itree(ictr)  = itfirsttemp(i)
         ictr = ictr + 1
      enddo

c     itlast - last target in box i
      ipointer(17)  = ictr
      do i=1,nboxes
         itree(ictr) = itlasttemp(i)
         ictr = ictr + 1
      enddo

c     ihefirst - first hung exansion center in
c     box i
      ipointer(18) = ictr
      do i=1,nboxes
          itree(ictr) = ihefirsttemp(i)
          ictr = ictr + 1
      enddo

c     ihelast - last hung expansion center in
c     box i
      ipointer(19) = ictr
      do i=1,nboxes
         itree(ictr) = ihelasttemp(i)
         ictr = ictr + 1
      enddo

c     iefirst - first exansion center in
c     box i
      ipointer(20) = ictr
      do i=1,nboxes
          itree(ictr) = iefirsttemp(i)
          ictr = ictr + 1
      enddo

c     ielast - last hung expansion center in
c     box i
      ipointer(21) = ictr
      do i=1,nboxes
         itree(ictr) = ielasttemp(i)
         ictr = ictr + 1
      enddo

c     nhungsrc  - Number of hung chunks in a given box
      ipointer(22)  = ictr
      do i=1,nboxes
         itree(ictr) = nhungsrc(i)
         ictr = ictr + 1
      enddo
 
c     nhungexp - Number of hung expansion centers
c                in box i
      ipointer(23) = ictr
      do i=1,nboxes
         itree(ictr) = nhungexp(i)
         ictr = ictr + 1
      enddo

c     nhunglistsrc - Number of hung chunks for box i
c        (ordered by given ordering of targets, not tree sorted)
      ipointer(24)  = ictr
      do i=1,nboxes
         itree(ictr) = nhunglistsrc(i)
         ictr = ictr + 1
      enddo

c     hunglistsrc - chunk hunglist for each target
      ipointer(25)  = ictr
      do i=1,nboxes
         do j=1,mhung
            itree(ictr) = ihunglistsrc(j,i)
            ictr=ictr+1
         enddo
      enddo
      ipointer(26) = ictr

      return
      end
c--------------------------------------------------------------
c
      subroutine subdivide(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,iparent,nchild,ichild,
     $                   isource,itarget,iexpc,
     $                   ihsfirst,ihslast,
     $                   isfirst,islast,
     $                   itfirst,itlast,
     $                   ihefirst,ihelast,
     $                   iefirst,ielast,nhungsrc,
     $                   nhungexp,irefine)
      implicit none
      integer ier,ns,nt,nexpc,idivflag,ndiv,mhung
      integer nlevels,nboxes,lcenters,nbmax,nlmax
      integer irefine
      real *8 src(3,ns),radsrc(ns)
      real *8 trg(3,nt)
      real *8 expc(3,nexpc),radexp(nexpc)
      real *8 centers(3,lcenters)
      real *8 boxsize(0:nlmax)
      integer laddr(2,0:nlmax)
      integer iparent(nbmax)
      integer nchild(nbmax)
      integer ichild(8,nbmax)
      integer isource(ns)
      integer itarget(nt)
      integer iexpc(nexpc)
      integer ihsfirst(nbmax)
      integer ihslast(nbmax)
      integer isfirst(nbmax)
      integer islast(nbmax)
      integer itfirst(nbmax)
      integer itlast(nbmax)
      integer ihefirst(nbmax)
      integer ihelast(nbmax)
      integer iefirst(nbmax)
      integer ielast(nbmax)
      integer nhungsrc(nbmax)
      integer nhungexp(nbmax)
c     Temporary variables
      integer isrctmp(ns),itargtmp(nt),iexpctmp(nexpc)
      integer i,j,i12,i34,istart,jstart,kstart,ii,iii,nss,nee
      integer jj
      integer i56, i78, i1234, i5678
      integer ibox,ifirstbox,ilastbox,nbfirst
      integer is,it,ie
      integer nsc(8),ntc(8),nh(8),nexpcc(8),nhc(8)
c
c     for every box at level nlevels,
c     sort into children, updating various arrays 
c     perhaps just build tree here paren/child/particle sorting...
c     lists in second call ???
c     
c     allocate temp array for isourcetemp2 itargtemp2
c     after all done, write back to isourcetemp, itargtemp
c     this is O(N) * nlevels work for rewriting.
c     can be fancier I suppose.
c     
      irefine = 0
      ifirstbox = laddr(1,nlevels)
      ilastbox =  laddr(2,nlevels)
c
      nbfirst = nboxes+1
      boxsize(nlevels+1) = boxsize(0)/2.0d0**(nlevels+1)

      do ibox = ifirstbox,ilastbox
c        Allocate temporary array to figure out which child you
c        belong to
c        which child?  1,2,3,4,5,6,7,8? counter ns1,ns2,ns3,ns4,
c        ns5,ns6,ns7,ns8
c        The box nomenclature is as follows
c        3   4       7  8       <--- Looking down in z direction
c        1   2       5  6
c
c        If the parent box center is at the origin then the centers
c        of box i have co-ordinates
c        1     x<0,y<0,z<0
c        2     x>0,y<0,z<0
c        3     x<0,y>0,z<0
c        4     x>0,y>0,z<0
c        5     x<0,y<0,z>0
c        6     x>0,y<0,z>0
c        7     x<0,y>0,z>0
c        8     x>0,y>0,z>0

         i1234 = isfirst(ibox)-1
         i5678 = 0
         do is = isfirst(ibox),islast(ibox)
            if(src(3,isource(is)) - centers(3,ibox).lt.0) then
                i1234 = i1234+1
                isource(i1234) = isource(is)
            else
                i5678 = i5678 + 1
                isrctmp(i5678) = isource(is)
            endif
         enddo
c        Note at the end of the loop, i1234 is where the particles
c        in part 1234 of the box end

c        Reorder sources to include sources in 5678 in the array
         do i=1,i5678
             isource(i1234+i) = isrctmp(i)
         enddo

c        Sort i1234 into i12 and i34         
         i12 = isfirst(ibox)-1
         i34 = 0
         do is = isfirst(ibox),i1234
            if(src(2,isource(is))-centers(2,ibox).lt.0) then
                i12 = i12 + 1
                isource(i12) = isource(is)
            else
               i34 = i34 + 1
               isrctmp(i34) = isource(is)
            endif
         enddo
c        Note at the end of the loop, i12 is where the particles
c        in part 12 of the box end
c
c        Reorder sources to include 34 in the array
         do i=1,i34
             isource(i12+i) = isrctmp(i)
         enddo

c        sort i5678 into i56 and i78
         i56 = i1234
         i78 = 0
         do is=i1234+1,islast(ibox)
            if(src(2,isource(is))-centers(2,ibox).lt.0) then
               i56 = i56 + 1
               isource(i56) = isource(is)
            else
               i78 = i78 + 1
               isrctmp(i78) = isource(is)
            endif
         enddo

c        Reorder sources to include 78 in the array
         do i=1,i78
             isource(i56+i) = isrctmp(i)
         enddo
c        End of reordering i5678         

         nsc(1) = 0
         nsc(2) = 0
         nsc(3) = 0
         nsc(4) = 0
         nsc(5) = 0
         nsc(6) = 0
         nsc(7) = 0
         nsc(8) = 0
c        Sort into boxes 1 and 2
         do is = isfirst(ibox),i12
            if(src(1,isource(is))-centers(1,ibox).lt.0) then
               isource(isfirst(ibox)+nsc(1)) = isource(is)
               nsc(1) = nsc(1) + 1
            else
               nsc(2) = nsc(2) + 1
               isrctmp(nsc(2)) = isource(is)
            endif
         enddo
c        Reorder sources so that sources in 2 are at the
c        end of this part of the array
         do i=1,nsc(2)
             isource(isfirst(ibox)+nsc(1)+i-1) = isrctmp(i)
         enddo

c        Sort into boxes 3 and 4
         do is = i12+1, i1234
             if(src(1,isource(is))-centers(1,ibox).lt.0) then
                 isource(i12+1+nsc(3)) = isource(is)
                 nsc(3) = nsc(3) + 1
             else
                 nsc(4) = nsc(4)+1
                 isrctmp(nsc(4)) = isource(is)
             endif
         enddo
c        Reorder sources so that sources in 4 are at the
c        end of this part of the array
         do i=1,nsc(4)
             isource(i12+nsc(3)+i) = isrctmp(i)
         enddo

c        Sort into boxes 5 and 6
         do is = i1234+1,i56
            if(src(1,isource(is))-centers(1,ibox).lt.0) then
               isource(i1234+1+nsc(5)) = isource(is)
               nsc(5) = nsc(5) + 1
            else
               nsc(6) = nsc(6) + 1
               isrctmp(nsc(6)) = isource(is)
            endif
         enddo
c        Reorder sources so that sources in 6 are at the
c        end of this part of the array
         do i=1,nsc(6)
            isource(i1234+nsc(5)+i) = isrctmp(i)
         enddo
c        End of sorting sources into boxes 5 and 6

c        Sort into boxes 7 and 8
         do is=i56+1,islast(ibox)
            if(src(1,isource(is))-centers(1,ibox).lt.0) then
               isource(i56+1+nsc(7)) = isource(is)
               nsc(7) = nsc(7) + 1
            else
               nsc(8) = nsc(8) + 1
               isrctmp(nsc(8)) = isource(is)
            endif
         enddo
c        Reorder sources so that sources in 8 are at the
c        end of the array
         do i=1,nsc(8)
            isource(i56+nsc(7)+i) = isrctmp(i)
         enddo

         istart = isfirst(ibox)-1
         do j=1,8
c        check hung -> counter nh1,nh2,nh3,nh4,nh5,nh6,nh7,nh8
             ii = 0
             nh(j) = 0
             do i=1,nsc(j)
                if(radsrc(isource(istart+i)).gt.boxsize(nlevels+1))
     1          then     
                     nh(j) = nh(j) + 1
                     isource(istart+nh(j)) = isource(istart+i)
                 else
                     ii = ii+1
                     isrctmp(ii) = isource(istart+i)
                 endif
             enddo
c        Reorder sources to have hung chunks at the star
c        of the sorted sources in the box ibox
             do i=1,ii
                isource(istart+nh(j)+i) = isrctmp(i)
             enddo
             istart = istart + nsc(j)
         enddo
c        End of sorting sources

c        Sort targets
c        which child?  1,2,3,4,5,6,7,8? counter nt1,nt2,nt3,nt4,
c        nt5,nt6,nt7,nt8
c        The box nomenclature is as follows
c        3   4       7  8       <--- Looking down in z direction
c        1   2       5  6
c
c        If the parent box center is at the origin then the centers
c        of box i have co-ordinates
c        1     x<0,y<0,z<0
c        2     x>0,y<0,z<0
c        3     x<0,y>0,z<0
c        4     x>0,y>0,z<0
c        5     x<0,y<0,z>0
c        6     x>0,y<0,z>0
c        7     x<0,y>0,z>0
c        8     x>0,y>0,z>0
c
         i1234 = itfirst(ibox)-1
         i5678 = 0
         do it = itfirst(ibox),itlast(ibox)
            if(trg(3,itarget(it)) - centers(3,ibox).lt.0) then
                i1234 = i1234+1
                itarget(i1234) = itarget(it)
            else
                i5678 = i5678 + 1
                itargtmp(i5678) = itarget(it)
            endif
         enddo
c        Reorder sources to include targets in 5678 in the array
         do i=1,i5678
             itarget(i1234+i) = itargtmp(i)
         enddo

c        Sort i1234 into i12 and i34         
         i12 = itfirst(ibox)-1
         i34 = 0
         do it = itfirst(ibox),i1234
            if(trg(2,itarget(it))-centers(2,ibox).lt.0) then
                i12 = i12 + 1
                itarget(i12) = itarget(it)
            else
               i34 = i34 + 1
               itargtmp(i34) = itarget(it)
            endif
         enddo
c        Note at the end of the loop, i12 is where the particles
c        in part 12 of the box end
c
c        Reorder targets to include 34 in the array
         do i=1,i34
             itarget(i12+i) = itargtmp(i)
         enddo

c        sort i5678 into i56 and i78
         i56 = i1234
         i78 = 0
         do it=i1234+1,itlast(ibox)
            if(trg(2,itarget(it))-centers(2,ibox).lt.0) then
               i56 = i56 + 1
               itarget(i56) = itarget(it)
            else
               i78 = i78 + 1
               itargtmp(i78) = itarget(it)
            endif
         enddo

c        Reorder sources to include 78 in the array
         do i=1,i78
             itarget(i56+i) = itargtmp(i)
         enddo
c        End of reordering i5678         

         ntc(1) = 0
         ntc(2) = 0
         ntc(3) = 0
         ntc(4) = 0
         ntc(5) = 0
         ntc(6) = 0
         ntc(7) = 0
         ntc(8) = 0

c        Sort into boxes 1 and 2
         do it = itfirst(ibox),i12
            if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
               itarget(itfirst(ibox)+ntc(1)) = itarget(it)
               ntc(1) = ntc(1) + 1
            else
               ntc(2) = ntc(2) + 1
               itargtmp(ntc(2)) = itarget(it)
            endif
         enddo
c        Reorder targets so that sources in 2 are at the
c        end of this part of the array
         do i=1,ntc(2)
             itarget(itfirst(ibox)+ntc(1)+i-1) = itargtmp(i)
         enddo
c        Sort into boxes 3 and 4
         do it = i12+1, i1234
             if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
                 itarget(i12+1+ntc(3)) = itarget(it)
                 ntc(3) = ntc(3) + 1
             else
                 ntc(4) = ntc(4)+1
                 itargtmp(ntc(4)) = itarget(it)
             endif
         enddo
c        Reorder targets so that sources in 4 are at the
c        end of this part of the array
         do i=1,ntc(4)
             itarget(i12+ntc(3)+i) = itargtmp(i)
         enddo

c        Sort into boxes 5 and 6
         do it = i1234+1,i56
            if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
               itarget(i1234+1+ntc(5)) = itarget(it)
               ntc(5) = ntc(5) + 1
            else
               ntc(6) = ntc(6) + 1
               itargtmp(ntc(6)) = itarget(it)
            endif
         enddo
c        Reorder targets so that sources in 6 are at the
c        end of this part of the array
         do i=1,ntc(6)
            itarget(i1234+ntc(5)+i) = itargtmp(i)
         enddo
c        End of sorting sources into boxes 5 and 6

c        Sort into boxes 7 and 8
         do it=i56+1,itlast(ibox)
            if(trg(1,itarget(it))-centers(1,ibox).lt.0) then
               itarget(i56+1+ntc(7)) = itarget(it)
               ntc(7) = ntc(7) + 1
            else
               ntc(8) = ntc(8) + 1
               itargtmp(ntc(8)) = itarget(it)
            endif
         enddo
c        Reorder targets so that sources in 8 are at the
c        end of the array
         do i=1,ntc(8)
            itarget(i56+ntc(7)+i) = itargtmp(i)
         enddo
c        End of sorting targets

c        Sort expansion centers
c        which child?  1,2,3,4,5,6,7,8? counter
c        nexpcc1, nexpcc2, nexpcc3, nexpcc4,
c        nexpcc5, nexpcc6, nexpcc7, nexpcc8,
c        The box nomenclature is as follows
c        3   4       7  8       <--- Looking down in z direction
c        1   2       5  6
c
c        If the parent box center is at the origin then the centers
c        of box i have co-ordinates
c        1     x<0,y<0,z<0
c        2     x>0,y<0,z<0
c        3     x<0,y>0,z<0
c        4     x>0,y>0,z<0
c        5     x<0,y<0,z>0
c        6     x>0,y<0,z>0
c        7     x<0,y>0,z>0
c        8     x>0,y>0,z>0
c
         i1234 = iefirst(ibox)-1
         i5678 = 0
         do ie = iefirst(ibox),ielast(ibox)
            if(expc(3,iexpc(ie)) - centers(3,ibox).lt.0) then
                i1234 = i1234+1
                iexpc(i1234) = iexpc(ie)
            else
                i5678 = i5678 + 1
                iexpctmp(i5678) = iexpc(ie)
            endif
         enddo

c        Reorder sources to include sources in 5678 in the array
         do i=1,i5678
             iexpc(i1234+i) = iexpctmp(i)
         enddo

c        Sort i1234 into i12 and i34         
         i12 = iefirst(ibox)-1
         i34 = 0
         do ie = iefirst(ibox),i1234
            if(expc(2,iexpc(ie))-centers(2,ibox).lt.0) then
                i12 = i12 + 1
                iexpc(i12) = iexpc(ie)
            else
               i34 = i34 + 1
               iexpctmp(i34) = iexpc(ie)
            endif
         enddo
c        Note at the end of the loop, i12 is where the particles
c        in part 12 of the box end
c
c        Reorder sources to include 34 in the array
         do i=1,i34
             iexpc(i12+i) = iexpctmp(i)
         enddo

c        sort i5678 into i56 and i78
         i56 = i1234
         i78 = 0
         do ie=i1234+1,ielast(ibox)
            if(expc(2,iexpc(ie))-centers(2,ibox).lt.0) then
               i56 = i56 + 1
               iexpc(i56) = iexpc(ie)
            else
               i78 = i78 + 1
               iexpctmp(i78) = iexpc(ie)
            endif
         enddo

c        Reorder sources to include 78 in the array
         do i=1,i78
             iexpc(i56+i) = iexpctmp(i)
         enddo
c        End of reordering i5678         

         nexpcc(1) = 0
         nexpcc(2) = 0
         nexpcc(3) = 0
         nexpcc(4) = 0
         nexpcc(5) = 0
         nexpcc(6) = 0
         nexpcc(7) = 0
         nexpcc(8) = 0

c        Sort into boxes 1 and 2
         do ie = iefirst(ibox),i12
            if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
               iexpc(iefirst(ibox)+nexpcc(1)) = iexpc(ie)
               nexpcc(1) = nexpcc(1) + 1
            else
               nexpcc(2) = nexpcc(2) + 1
               iexpctmp(nexpcc(2)) = iexpc(ie)
            endif
         enddo
c        Reorder expc so that sources in 2 are at the
c        end of this part of the array
         do i=1,nexpcc(2)
             iexpc(iefirst(ibox)+nexpcc(1)+i-1) = iexpctmp(i)
         enddo
c        Sort into boxes 3 and 4
         do ie = i12+1, i1234
             if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
                 iexpc(i12+1+nexpcc(3)) = iexpc(ie)
                 nexpcc(3) = nexpcc(3) + 1
             else
                 nexpcc(4) = nexpcc(4)+1
                 iexpctmp(nexpcc(4)) = iexpc(ie)
             endif
         enddo
c        Reorder expc so that sources in 4 are at the
c        end of this part of the array
         do i=1,nexpcc(4)
             iexpc(i12+nexpcc(3)+i) = iexpctmp(i)
         enddo

c        Sort into boxes 5 and 6
         do ie = i1234+1,i56
            if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
               iexpc(i1234+1+nexpcc(5)) = iexpc(ie)
               nexpcc(5) = nexpcc(5) + 1
            else
               nexpcc(6) = nexpcc(6) + 1
               iexpctmp(nexpcc(6)) = iexpc(ie)
            endif
         enddo
c        Reorder expc so that sources in 6 are at the
c        end of this part of the array
         do i=1,nexpcc(6)
            iexpc(i1234+nexpcc(5)+i) = iexpctmp(i)
         enddo
c        End of sorting sources into boxes 5 and 6

c        Sort into boxes 7 and 8
         do ie=i56+1,ielast(ibox)
            if(expc(1,iexpc(ie))-centers(1,ibox).lt.0) then
               iexpc(i56+1+nexpcc(7)) = iexpc(ie)
               nexpcc(7) = nexpcc(7) + 1
            else
               nexpcc(8) = nexpcc(8) + 1
               iexpctmp(nexpcc(8)) = iexpc(ie)
            endif
         enddo
c        Reorder expc so that sources in 8 are at the
c        end of the array
         do i=1,nexpcc(8)
            iexpc(i56+nexpcc(7)+i) = iexpctmp(i)
         enddo
c        End of sorting expanison centers

         istart = iefirst(ibox)-1
         do j=1,8
c        check hung -> counter nhc1,nhc2,nhc3,nhc4,nhc5,nhc6,nhc7,nhc8
             ii = 0
             nhc(j) = 0
             do i=1,nexpcc(j)
                if(radexp(iexpc(istart+i)).gt.boxsize(nlevels+1)/4)
     1          then     
                     nhc(j) = nhc(j) + 1
                     iexpc(istart+nhc(j)) = iexpc(istart+i)
                 else
                     ii = ii+1
                     iexpctmp(ii) = iexpc(istart+i)
                 endif
             enddo
c        Reorder sources to have hung chunks at the star
c        of the sorted sources in the box ibox
             do i=1,ii
                iexpc(istart+nhc(j)+i) = iexpctmp(i)
             enddo
             istart = istart + nexpcc(j)
         enddo

         nchild(ibox) = 0
c        Create the required boxes
         istart = isfirst(ibox)
         jstart = itfirst(ibox)
         kstart = iefirst(ibox)
         do i=1,8
             ii = 2
             jj = 2
             if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
             if(i.lt.5) jj = 1
             if(nsc(i)+ntc(i)+nexpcc(i).gt.0) then
c               Increment total number of boxes               
                nboxes = nboxes + 1
c               Increment number of children for the current box
                nchild(ibox) = nchild(ibox)+1
c               Update the array of children for the current box
                ichild(i,ibox) = nboxes
c               Update the array of parents for the child box
                iparent(nboxes) = ibox
c               Compute center for the child box
                centers(1,nboxes) = centers(1,ibox)+(-1)**i*
     1                              boxsize(nlevels+1)/2.0
                centers(2,nboxes) = centers(2,ibox)+(-1)**ii*
     1                              boxsize(nlevels+1)/2.0
                centers(3,nboxes) = centers(3,ibox)+(-1)**jj*
     1                              boxsize(nlevels+1)/2.0

                nchild(nboxes) = 0
                ichild(1,nboxes) = -1
                ichild(2,nboxes) = -1
                ichild(3,nboxes) = -1
                ichild(4,nboxes) = -1
                ichild(5,nboxes) = -1
                ichild(6,nboxes) = -1
                ichild(7,nboxes) = -1
                ichild(8,nboxes) = -1
c               Update arrays ihsfirst,ihslast,isfirst,islast
                ihsfirst(nboxes) = istart
                ihslast(nboxes) = istart + nh(i) - 1
                nhungsrc(nboxes) = nh(i)

                isfirst(nboxes) = istart + nh(i)
                islast(nboxes) = istart + nsc(i) - 1

c               Update arrays itfirst, itlast
                itfirst(nboxes) = jstart
                itlast(nboxes) = jstart + ntc(i) - 1

c               Update arrays ihefirst,ihelast,iefirst,ielast
                ihefirst(nboxes) = kstart
                ihelast(nboxes) = kstart + nhc(i)-1
                nhungexp(nboxes) = nhc(i)

                iefirst(nboxes) = kstart + nhc(i)
                ielast(nboxes) = kstart + nexpcc(i) - 1

                nss = islast(nboxes) - isfirst(nboxes)+1
                nee = ielast(nboxes) - iefirst(nboxes)+1 
c               Check if further refinement required
                if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
                if ((idivflag .eq.1).and.(ntc(i).gt.ndiv)) irefine=1
                if ((idivflag .eq.2).and.(nss+ntc(i).gt.ndiv)) irefine=1
                if ((idivflag .eq.3).and.(nss+ntc(i)+nee.gt.
     1                                    ndiv)) irefine=1  
             endif
             istart = istart + nsc(i)
             jstart = jstart + ntc(i)
             kstart = kstart + nexpcc(i)
         enddo
      enddo
      nlevels = nlevels+1
      laddr(1,nlevels) = nbfirst
      laddr(2,nlevels) = nboxes

      return
      end
c------------------------------------------------------------------      
      subroutine computelists(nlevels,nboxes,laddr,boxsize,
     1                        centers,iparent,nchild,ichild,isep,
     2                        nnbors,nbors,mnbors,nlist2,list2,mnlist2)

c     This subroutine computes the list1 and list2 for the
c     tree structure
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: real *8(0:;nlevels)
c                 Array of boxsizes
c 
c     centers     in: real *8(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(8,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     isep        in: integer
c                 separation parameter defining list1 of the box.
c                 if isep=1, then list1 includes self and nearest
c                 neighbors. If isep=2, then list1 includes
c                 two neighbors in each direction
c
c     mnbors      in: integer
c                 max number of boxes in list1. If isep=1, 
c                 mnbors should be 27 and if isep=2, mnbors
c                 should be 125
c
c     mnlist2     in: integer
c                 max number of boxes in list2. If isep=1,
c                 mnlist2 should be 189 and if isep=2, 
c                 mnlist2 should be 875
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       out: integer(mnbors,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     nlist2      out: integer(nboxes)
c                 nlist2(i) is the number of boxes in the list 2
c                 of box i
c 
c     list2       out: integer(mnlist2,nboxes)
c                 list2(j,i) is the box id of the jth box in 
c                 list2 of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer mnbors, mnlist2, isep
      integer laddr(2,0:nlevels)
      real *8 boxsize(0:nlevels)
      real *8 centers(3,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(8,nboxes)
      integer nnbors(nboxes), nlist2(nboxes)
      integer nbors(mnbors,nboxes), list2(mnlist2,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox

c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      nlist2(1) = 0
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,8
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev)*isep).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev)*isep).and.
     4                   (abs(centers(3,kbox)-centers(3,ibox)).le.
     5                   1.05*boxsize(ilev)*isep)) then
                     
                      nnbors(ibox) = nnbors(ibox)+1
                      nbors(nnbors(ibox),ibox) = kbox
                      else
                         nlist2(ibox) = nlist2(ibox)+1
                         list2(nlist2(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing list1 and list2 for ibox
         enddo
      enddo

      return
      end
c------------------------------------------------------------------

      subroutine computemhung(nlevels,nboxes,laddr,iparent,nnbors,
     1                        nbors,mnbors,nhungsrc,nhunglistsrc,mhung)
c     This subroutine computes mhung for a given tree 
c     and the number of hung sources per box
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c    
c     iparent     in: integer(nboxes)
c                 iparent(i) is the parent of box i
c
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     mnbors      in: integer
c                 max number of boxes in list1. If isep=1, 
c                 mnbors should be 27 and if isep=2, mnbors
c                 should be 125
c
c     nhungsrc   in: integer(nboxes)
c                 nhung(i) is the number of hung sources in box i
c 
c      OUTPUT
c      nhunglistsrc  out:integer(nboxes)
c                    nhunglistsrc(i) is the number of hung sources
c                    relevant to box i, nhunglistsrc(i) = 
c                    nhunglistsrc(idad) + \sum_{j=1}{nnbors(i)}
c                    nhungsrc(nbors(j,i))
c
c      mhung         out: integer
c                    max(nhunglistsrc)
      implicit none
      integer nlevels,nboxes, mnbors
      integer laddr(2,0:nlevels)
      integer iparent(nboxes),nnbors(nboxes), nbors(mnbors,nboxes)
      integer nhungsrc(nboxes), nhunglistsrc(nboxes)
      integer mhung

c     Temporary variables
      integer ilev,ibox,i,dad,jbox

c     initialize for level 0      
      nhunglistsrc(1) = nhungsrc(1)
      do ilev = 1,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            dad = iparent(ibox)
            nhunglistsrc(ibox) = nhunglistsrc(dad)
            do i=1,nnbors(ibox)
               jbox = nbors(i,ibox)
               nhunglistsrc(ibox) = nhunglistsrc(ibox) +
     1                              nhungsrc(jbox)
            enddo
         enddo
      enddo

      mhung = 0
      do ibox = laddr(1,nlevels),laddr(2,nlevels)
         if(nhunglistsrc(ibox).gt.mhung) mhung = nhunglistsrc(ibox)
      enddo


      return
      end
c------------------------------------------------------------------
      subroutine computehunglist(mhung,nlevels,nboxes,laddr,ns,
     1           isource,iparent,ihsfirst,ihslast,nnbors,nbors,mnbors,
     2           nhungsrc,nhunglistsrc,ihunglistsrc)
c     This subroutine computes the hung list for each box 
c     and stores it. The hung list of sources of a box is
c     the hunglist of the parent + the sources hung in list1
c     of the box
c 
c     INPUT arguments
c     mhung         in: Integer
c                   max number of hung chunks for a box
c
c     nlevels       in: Integer
c                   number of levels in the tree
c
c     nboxes        in: Integer
c                   number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     ns          in: integer
c                 number of sources
c
c     isource     in: integer(ns)
c                 Mapping from tree sorted array of sources
c                 to the user prescribed ordering
c  
c     iparent     in: integer(nboxes)
c                 iparent(i) is the parent of box i
c     
c     ihsfirst    in: integer(nboxes)
c                 ihsfirst(i) points to the first hung source
c                 in box i in the sorted array
c
c      ihslast    in: integer(nboxes)
c                 ihslast(i) points to the last hung source
c                 in box i in the sorted array
c
c     nnbors      in: integer(nboxes)
c                 nnbors(i) is the number of boxes in list 1 of box i
c
c     nbors       in: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth box in
c                 list 1 of box i
c
c     mnbors      in: integer
c                 max number of boxes in list1. If isep=1, 
c                 mnbors should be 27 and if isep=2, mnbors
c                 should be 125
c
c     nhungsrc   in: integer(nboxes)
c                 nhung(i) is the number of hung sources in box i
c 
c      OUTPUT
c      nhunglistsrc  out:integer(nboxes)
c                    nhunglistsrc(i) is the number of hung sources
c                    relevant to box i 
c
c      ihunglistsrc  out: integer(mhung,nboxes)
c                    ihunglistsrc(j,i) is the id of the jth
c                    hung source relevant to box i
c-----------------------------------------------------------------
       implicit none
       integer mhung,nlevels,nboxes,ns, mnbors
       integer laddr(2,0:nlevels)
       integer iparent(nboxes)
       integer isource(ns), ihsfirst(nboxes),ihslast(nboxes)
       integer nnbors(nboxes), nbors(mnbors,nboxes)
       integer nhungsrc(nboxes),nhunglistsrc(nboxes)
       integer ihunglistsrc(mhung,nboxes)
c      Temp variables
       integer i,j,ibox,jbox,ilev,dad


c      Compute hung list for root node       
       nhunglistsrc(1) = nhungsrc(1)
       do i=1,nhungsrc(1)
          ihunglistsrc(i,1) = isource(ihsfirst(1)+i-1)
       enddo

       do ilev=1,nlevels
          do ibox = laddr(1,ilev),laddr(2,ilev)
             dad = iparent(ibox)
             do i=1,nhunglistsrc(dad)
                ihunglistsrc(i,ibox) = ihunglistsrc(i,dad)
             enddo
             nhunglistsrc(ibox) = nhunglistsrc(dad)
             do i=1,nnbors(ibox)
                jbox = nbors(i,ibox)
                do j=1,nhungsrc(jbox)
                   ihunglistsrc(nhunglistsrc(ibox)+j,ibox)=
     1             isource(ihsfirst(jbox)+j-1)
                enddo
                nhunglistsrc(ibox) = nhunglistsrc(ibox)+
     1          nhungsrc(jbox)
             enddo
          enddo
       enddo

       return
       end
c------------------------------------------------------------------

      subroutine mkptreemem(src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,isep,
     $                   nlmax,nlevels,nboxes,mhung,ltree)

      implicit none
      integer ier
      integer ns,nt,nexpc,idivflag,ndiv,mhung,isep
      integer nlevels,nboxes,lcenters,ltree
      integer i,j,nbmax,nlmax
      integer, allocatable :: laddr(:,:)
      integer, allocatable :: iparenttemp(:)
      integer, allocatable :: nchild(:)
      integer, allocatable :: ichildtemp(:,:)
      integer, allocatable :: nnbor(:)
      integer, allocatable :: nbors(:,:)
      integer, allocatable :: nlist2(:)
      integer, allocatable :: list2(:,:)
      integer, allocatable :: isourcetemp(:)
      integer, allocatable :: itargettemp(:)
      integer, allocatable :: iexpctemp(:)
      integer, allocatable :: ihsfirsttemp(:)
      integer, allocatable :: ihslasttemp(:)
      integer, allocatable :: isfirsttemp(:)
      integer, allocatable :: islasttemp(:)
      integer, allocatable :: itfirsttemp(:)
      integer, allocatable :: itlasttemp(:)
      integer, allocatable :: ihefirsttemp(:)
      integer, allocatable :: ihelasttemp(:)
      integer, allocatable :: iefirsttemp(:)
      integer, allocatable :: ielasttemp(:)
      integer, allocatable :: nhungsrc(:)
      integer, allocatable :: nhungexp(:)
      integer, allocatable :: nhunglistsrc(:)

      real *8 boxsize(0:nlmax)
      real *8 src(3,ns),radsrc(ns)
      real *8 trg(3,nt)
      real *8, allocatable :: centers(:,:)
      real *8 expc(3,nexpc)
      real *8 radexp(nexpc)
c
      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey, sizez
      integer ictr,ih,irefine,is,ie
      integer nss,nee

      integer mnbors, mnlist2

c     This code is a memory code for mkptree. This
c     subroutine determines the number of boxes (lcenters) required,
c     the length of the tree array (ltree) and the maximum number
c     of hung chunks (mhung)
c     
c     INPUT:
c
c     src           source locations        
c     ns            number of sources 
c     rads          source radii (determines deepest level that
c                   the source can reach) 
c     trg           target locations        
c     nt            number of targets
c
c     expc          expansion center locations
c     nexpc         number of expansion centers
c
c     radexp        expansion center radius
c
c     idivflag      0 => divide on sources
c                   1 => divide on targets
c                   2 => divide on sources+targets
c                   3 => divide on sources+targets+expansion centers
c 
c     ndiv          refinement criterion - extend tree until all
c                   nodes at finest level have less than ndiv 
c                   source/targets/sources+targets depending on
c                   idivflag
c
c     nlmax         max number of levels
c
c     OUTPUT:
c     nlevels       number of levels
c     nboxes        number of boxes
c     mhung         max number of hung chunks
c     ltree         length of the tree
c
c     mhung is the max number of hung sources on output
c
c     This subroutine essentially creates the quad tree
c     and in process determines the memory required for the tree.
c
c     It can be shown that the maximum number of boxes is less
c     than 2*(nexpc + ns + nt)
c
c     For details of the variables used in the tree, refer
c     to the main routine
c---------------------------------------------------------------------
c
c     initialize temporary arrays
c  
c     assumes nlevels lt 200
c     assumes nboxes lt lcenters
c
c     Other notes:
c
c     refinement criterion w.r.t rads is:
c
c     hang if (rads .geq. boxsize)
c
c
      ier = 0
      nbmax = 16*(ns+nt+nexpc) 
      lcenters = nbmax
c
      allocate(laddr(2,0:nlmax))
      allocate(iparenttemp(nbmax))
      allocate(nchild(nbmax))
      allocate(ichildtemp(8,nbmax))
      allocate(isourcetemp(ns))
      allocate(itargettemp(nt))
      allocate(iexpctemp(nexpc))
      allocate(ihsfirsttemp(nbmax))
      allocate(ihslasttemp(nbmax))
      allocate(isfirsttemp(nbmax))
      allocate(islasttemp(nbmax))
      allocate(itfirsttemp(nbmax))
      allocate(itlasttemp(nbmax))
      allocate(ihefirsttemp(nbmax))
      allocate(ihelasttemp(nbmax))
      allocate(iefirsttemp(nbmax))
      allocate(ielasttemp(nbmax))
      allocate(nhungsrc(nbmax))
      allocate(nhungexp(nbmax))
      allocate(centers(3,lcenters))

      if(isep.eq.1) then
         mnbors = 27
         mnlist2 = 27*7
      endif

      if(isep.eq.2) then
         mnbors = 125
         mnlist2 = 125*7
      endif

      if(isep.ne.1.and.isep.ne.2) then
         ier = 4
         call prinf('Error tree not allocated, isep.ne.1,2*',ier,0)
         return
      endif

c     Step 1: Find enclosing box
c
      xmin = src(1,1)
      xmax = src(1,1)
      ymin = src(2,1)
      ymax = src(2,1)
      zmin = src(3,1)
      zmax = src(3,1)
c
      do i=1,ns
        if(src(1,i) .lt. xmin) xmin=src(1,i)
        if(src(1,i) .gt. xmax) xmax=src(1,i)
        if(src(2,i) .lt. ymin) ymin=src(2,i)
        if(src(2,i) .gt. ymax) ymax=src(2,i)
        if(src(3,i) .lt. zmin) zmin=src(3,i)
        if(src(3,i) .gt. zmax) zmax=src(3,i)
      enddo
      do i=1,nt
        if(trg(1,i) .lt. xmin) xmin=trg(1,i)
        if(trg(1,i) .gt. xmax) xmax=trg(1,i)
        if(trg(2,i) .lt. ymin) ymin=trg(2,i)
        if(trg(2,i) .gt. ymax) ymax=trg(2,i)
        if(trg(3,i) .lt. zmin) zmin=trg(3,i)
        if(trg(3,i) .gt. zmax) zmax=trg(3,i)
      enddo

      do i=1,nexpc
        if(expc(1,i) .lt. xmin) xmin=expc(1,i)
        if(expc(1,i) .gt. xmax) xmax=expc(1,i)
        if(expc(2,i) .lt. ymin) ymin=expc(2,i)
        if(expc(2,i) .gt. ymax) ymax=expc(2,i)
        if(expc(3,i) .lt. zmin) zmin=expc(3,i)
        if(expc(3,i) .gt. zmax) zmax=expc(3,i)
      enddo
      boxsize(0)=xmax-xmin
      sizey=ymax-ymin
      sizez=zmax-zmin
      if(sizey .gt. boxsize(0)) boxsize(0)=sizey
      if(sizez .gt. boxsize(0)) boxsize(0)=sizez
c
c     initialize arrays at level 0
c
      centers(1,1)=(xmin+xmax)/2
      centers(2,1)=(ymin+ymax)/2
      centers(3,1)=(zmin+zmax)/2
      laddr(1,0) = 1
      laddr(2,0) = 1
      iparenttemp(1) = -1
      isfirsttemp(1) = 1
      nhungsrc(1) = 0
      nhungexp(1) = 0
c
c     count number of hung sources
c     and hang up "big" sources
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) nhungsrc(1)=nhungsrc(1)+1 
      enddo
      isfirsttemp(1) = nhungsrc(1)+1
      islasttemp(1) = ns
      if (nhungsrc(1) .gt. 0) then 
         ihsfirsttemp(1) = 1
         ihslasttemp(1) = nhungsrc(1)
      else
         ihsfirsttemp(1) = 0
         ihslasttemp(1) = -1
      endif

c     Count number of hung expansion centers      
c     and hang up "big" expansion centers
      do i=1,nexpc
         if (radexp(i).gt.boxsize(0)) nhungexp(1)=nhungexp(1)+1
      enddo
      iefirsttemp(1) = nhungexp(1)+1
      ielasttemp(1) = nexpc
      if (nhungexp(1).gt.0) then
          ihefirsttemp(1) = 1
          ihelasttemp(1) = nhungexp(1)
      else
         ihefirsttemp(1) = 0
         ihelasttemp(1) = -1
      endif

c
c     reorder isourcetemp to put hung sources in beginning
      ih = 0
      is = nhungsrc(1)
      do i = 1,ns
         if (radsrc(i).gt.boxsize(0)) then
            ih = ih+1
            isourcetemp(ih) = i
         else
            is = is+1
            isourcetemp(is) = i
         endif
      enddo
c     reorder iexptemp to put hung expansion centers in beginning
      ih = 0
      ie = nhungexp(1)
      do i= 1,nexpc
         if(radexp(i).gt.boxsize(0)) then
            ih = ih+1
            iexpctemp(ih) = i
         else
            ie = ie+1
            iexpctemp(ie) = i
         endif
      enddo

c     initialize itargettemp 
      do i = 1,nt
         itargettemp(i) = i
      enddo
      itfirsttemp(1) = 1
      itlasttemp(1) = nt

      nlevels = 0
      nboxes = 1
      nchild(1) = 0
      ichildtemp(1,1) = -1
      ichildtemp(2,1) = -1
      ichildtemp(3,1) = -1
      ichildtemp(4,1) = -1
      ichildtemp(5,1) = -1
      ichildtemp(6,1) = -1
      ichildtemp(7,1) = -1
      ichildtemp(8,1) = -1
c
      irefine = 0
      nss = ns - nhungsrc(1)
      nee = nexpc - nhungexp(1)
      if ((idivflag .eq.0).and.(nss.gt.ndiv)) irefine=1
      if ((idivflag .eq.1).and.(nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.2).and.(nss+nt.gt.ndiv)) irefine=1
      if ((idivflag .eq.3).and.(nss+nt+nee.gt.ndiv)) irefine=1


      do i = 1,nlmax
         if (irefine.eq.1) then
            call subdivide(ier,src,ns,radsrc,trg,nt,expc,nexpc,
     $                   radexp,idivflag,ndiv,
     $                   nlevels,nboxes,
     $                   centers,lcenters,boxsize,nbmax,nlmax,
     $                   laddr,iparenttemp,nchild,ichildtemp,
     $                   isourcetemp,itargettemp,iexpctemp,
     $                   ihsfirsttemp,ihslasttemp,
     $                   isfirsttemp,islasttemp,
     $                   itfirsttemp,itlasttemp,
     $                   ihefirsttemp,ihelasttemp,
     $                   iefirsttemp,ielasttemp,nhungsrc,
     $                   nhungexp,irefine)
         else
            exit
         endif
      enddo
      call prinf('At end of refining in memory allocation code=*',i,0)

c     Set up computation of list1 and list2      
      allocate(nnbor(nboxes))
      allocate(nbors(mnbors,nboxes))
      allocate(nlist2(nboxes))
      allocate(list2(mnlist2,nboxes))

c     Initialize nnbor,nlist2,nbors,list2 arrays
      do i=1,nboxes
          nnbor(i) = 0
          nlist2(i) = 0
          do j=1,mnbors
              nbors(j,i) = -1
              list2(j,i) = -1
          enddo
          do j=mnbors+1,mnlist2
              list2(j,i) = -1
          enddo
      enddo

c     compute list1 and list2
      call computelists(nlevels,nboxes,laddr, boxsize,
     1                   centers,iparenttemp,nchild,
     2                   ichildtemp,isep,nnbor,nbors,mnbors,nlist2,
     3                   list2,mnlist2)

c     compute mhung
      allocate(nhunglistsrc(nboxes))
      do i=1,nboxes
         nhunglistsrc(i) = 0
      enddo
      mhung = 0
      call computemhung(nlevels,nboxes,laddr,iparenttemp,nnbor,nbors,
     1                  mnbors,nhungsrc,nhunglistsrc,mhung)


      ltree = (25+mhung+mnbors+mnlist2)*nboxes + 
     1        2*(nlevels+1) +ns+nt+nexpc

      return
      end

      subroutine mkpwlists(isep,nlevels,laddr,boxsize,
     1                   nboxes,centers,mnlist2,nlist2,list2,
     2                   nupm,ndownm,nnorthm,nsouthm,
     3                   neastm,nwestm,nup,uplist,ndown,downlist,
     3                   nnorth,northlist,nsouth,southlist,
     4                   neast,eastlist,nwest,westlist)

c     This subroutine computes the uplist, downlist, eastlist,
c     west list, north list and south list based on the list2 of
c     the tree structure. The up and down list for a box B
c     are the set of boxes in list 2 
c     that are separated by at least isep boxes
c     in the +z direction and -z direction respectively. The north
c     and south list of box B are those boxes which are not in
c     the up and down list and are well separated from the box B
c     by at least isep boxes in the +y and -y direction
c     respectively. The east and west list of box B are those
c     boxes which are in none of up,down,north and south list of
c     the box B and are well separated from box B by at least
c     isep boxes in the +x and -x direction respectively.
c
c     NOTE: to determine nupm,ndownm,nnorthm,nsouthm,neastm and
c           nwestm use the subroutine findnudnsew
c
c     INPUT arguments:
c     isep        in: Integer
c                 separation parameter
c
c     nlevels     in: integer
c                 number of levels
c
c     laddr       in: integer (2,0:nlevels)
c                 laddr(1,i) is the first box at level i
c                 laddr(2,i) is the last box at level i
c
c     boxsize     in: real *8(0:nlevels)
c                 boxsize(i) is the size of a box at level i
c
c     nboxes      in: Integer
c                 number of boxes
c
c     centers     in: real *8(3,nboxes)
c                 co-ordinates of the box centers in the
c                 tree structure
c
c     mnlist2     in: integer
c                 maximum number of elements in list2 of a box
c  
c     nlist2      in: integer(nboxes)
c                 nlist2(i) is the number of boxes in list2
c                 of box i
c
c     list2       in: integer(mnlist2,nboxes)
c                 list2(j,i) is the box number of the jth box
c                 in the list2 of box i
c
c     nupm        in: integer
c                 max number of boxes in uplist of any box
c
c     ndownm      in: integer
c                 max number of boxes in downlist of any box
c
c     nnorthm     in: integer
c                 max number of boxes in the north list of any box
c
c     nsouthm     in: integer
c                 max number of boxes in the south list of any box
c
c     neastm      in: integer
c                 max number of boxes in the east list of any box
c
c     nwestm      in: integer
c                 max number of boxes in the west list of any box
c
c     OUTPUT
c     nup         out: integer(nboxes)
c                 nup(i) is the number of boxes in the uplist
c                 of box i
c
c     uplist      out: integer(nupm,nboxes)
c                 uplist(j,i) is the jth box in the uplist
c                 of box i
c
c     ndown       out: integer(nboxes)
c                 ndown(i) is the number of boxes in the downlist
c                 of box i
c
c     downlist    out: integer(ndownm,nboxes)
c                 downlist(j,i) is the jth box in the downlist
c                 of box i
c
c     nnorth      out: integer(nboxes)
c                 nnorth(i) is the number of boxes in the northlist
c                 of box i
c
c     northlist   out: integer(nnorthm,nboxes)
c                 northlist(j,i) is the jth box in the northlist
c                 of box i
c
c     nsouth      out: integer(nboxes)
c                 nsouth(i) is the number of boxes in the southlist
c                 of box i
c
c     southlist   out: integer(ndownm,nboxes)
c                 southlist(j,i) is the jth box in the southlist
c                 of box i
c
c     neast       out: integer(nboxes)
c                 neast(i) is the number of boxes in the eastlist
c                 of box i
c
c     eastlist    out: integer(neastm,nboxes)
c                 eastlist(j,i) is the jth box in the eastlist
c                 of box i
c
c     nwest       out: integer(nboxes)
c                 nwest(i) is the number of boxes in the westlist
c                 of box i
c
c     westlist    out: integer(nwestm,nboxes)
c                 westlist(j,i) is the jth box in the westlist
c                 of box i
c---------------------------------------------------------------

      implicit none
      integer isep,nlevels,nboxes,mnlist2,nlist2(nboxes)
      integer laddr(2,0:nlevels)
      real *8 boxsize(0:nlevels)
      integer list2(mnlist2,nboxes)
      integer nup(nboxes),ndown(nboxes)
      integer nnorth(nboxes),nsouth(nboxes)
      integer neast(nboxes),nwest(nboxes)
      integer nupm,ndownm,nnorthm,nsouthm,neastm,nwestm
      integer uplist(nupm,nboxes),downlist(ndownm,nboxes)
      integer northlist(nnorthm,nboxes),southlist(nsouthm,nboxes)
      integer eastlist(neastm,nboxes),westlist(nwestm,nboxes)
      integer ilev,i,ibox,jbox
      real *8 centers(3,nboxes)

      do ilev = 0,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nup(ibox) = 0
            ndown(ibox) = 0
            nnorth(ibox) = 0
            nsouth(ibox) = 0
            neast(ibox) = 0
            nwest(ibox) = 0
            do i=1,nlist2(ibox)
              jbox = list2(i,ibox)
              if(centers(3,jbox) - centers(3,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nup(ibox) = nup(ibox)+1
                 uplist(nup(ibox),ibox) = jbox
                 goto 1111
              endif
              if(centers(3,jbox) - centers(3,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 ndown(ibox) = ndown(ibox)+1
                 downlist(ndown(ibox),ibox) = jbox
                 goto 1111
              endif

              if(centers(2,jbox) - centers(2,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nnorth(ibox) = nnorth(ibox)+1
                 northlist(nnorth(ibox),ibox) = jbox
                 goto 1111
              endif
              if(centers(2,jbox) - centers(2,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nsouth(ibox) = nsouth(ibox)+1
                 southlist(nsouth(ibox),ibox) = jbox
                 goto 1111
              endif

              if(centers(1,jbox) - centers(1,ibox).ge.
     1          1.01d0*isep*boxsize(ilev)) then
                 neast(ibox) = neast(ibox)+1
                 eastlist(neast(ibox),ibox) = jbox
                 goto 1111
              endif
              if(centers(1,jbox) - centers(1,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nwest(ibox) = nwest(ibox)+1
                 westlist(nwest(ibox),ibox) = jbox
                 goto 1111
              endif
1111          continue                 
           enddo
        enddo
      enddo

      return
      end

      subroutine findnudnsew(isep,nlevels,laddr,boxsize,
     1                       nboxes,centers,mnlist2,
     2                       nlist2,list2,nup,ndown,nnorth,
     3                       nsouth,neast,nwest)
  
c     This subroutine computes the maximum number of elements
c     in the uplist,downlist, eastlist, west list, north list
c     and south list based on the list2 of the tree structure.
c     For detailed description on the lists, refer to mkpwlists
c
c     INPUT arguments:
c     isep        in: Integer
c                 separation parameter
c
c     nlevels     in: integer
c                 number of levels
c
c     laddr       in: integer (2,0:nlevels)
c                 laddr(1,i) is the first box at level i
c                 laddr(2,i) is the last box at level i
c
c     boxsize     in: real *8(0:nlevels)
c                 boxsize(i) is the size of a box at level i
c
c     nboxes      in: Integer
c                 number of boxes
c
c     centers     in: real *8(3,nboxes)
c                 co-ordinates of the box centers in the
c                 tree structure
c
c     mnlist2     in: integer
c                 maximum number of elements in list2 of a box
c  
c     nlist2      in: integer(nboxes)
c                 nlist2(i) is the number of boxes in list2
c                 of box i
c
c     list2       in: integer(mnlist2,nboxes)
c                 list2(j,i) is the box number of the jth box
c                 in the list2 of box i
c    
c     OUTPUT
c     nup         out: integer
c                 max number of boxes in the uplist of any box
c
c     ndown       out: integer
c                 max number of boxes in the downlist of any box
c
c     nnorth      out: integer
c                 max number of boxes in the northlist of any box
c
c     nsouth      out: integer
c                 max number of boxes in the southlist of any box
c
c     neast       out: integer
c                 max number of boxes in the eastlist of any box
c
c     nwest       out: integer
c                 max number of boxes in the westlist of any box
c
c----------------------------------------------------------------
      implicit none
      integer isep,nlevels,nboxes,mnlist2,nlist2(nboxes)
      integer laddr(2,0:nlevels)
      real *8 boxsize(0:nlevels)
      integer list2(mnlist2,nboxes)
      integer nup,ndown,nnorth,nsouth,neast,nwest
      integer nupt,ndownt,nnortht,nsoutht,neastt,nwestt
      integer ilev,i,ibox,jbox
      real *8 centers(3,nboxes)

      nup = 0
      ndown = 0
      nnorth = 0
      nsouth = 0
      neast = 0
      nwest = 0
      do ilev = 2,nlevels
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nupt = 0
            ndownt = 0
            nnortht = 0
            nsoutht = 0
            neastt = 0
            nwestt = 0
            do i=1,nlist2(ibox)
              jbox = list2(i,ibox)
              if(centers(3,jbox) - centers(3,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nupt = nupt+1
                 goto 1111
              endif
              if(centers(3,jbox) - centers(3,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 ndownt = ndownt+1
                 goto 1111
              endif

              if(centers(2,jbox) - centers(2,ibox).ge.
     1           1.01d0*isep*boxsize(ilev)) then
                 nnortht = nnortht+1
                 goto 1111
              endif
              if(centers(2,jbox) - centers(2,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nsoutht = nsoutht+1
                 goto 1111
              endif

              if(centers(1,jbox) - centers(1,ibox).ge.
     1          1.01d0*isep*boxsize(ilev)) then
                 neastt = neastt+1
                 goto 1111
              endif
              if(centers(1,jbox) - centers(1,ibox).le.
     1           -1.01d0*isep*boxsize(ilev)) then
                 nwestt = nwestt+1
                 goto 1111
              endif
1111          continue                 
           enddo
           if(nupt.gt.nup) nup = nupt
           if(ndownt.gt.ndown) ndown = ndownt
           if(nnortht.gt.nnorth) nnorth = nnortht
           if(nsoutht.gt.nsouth) nsouth = nsoutht
           if(neastt.gt.neast) neast = neastt
           if(nwestt.gt.nwest) nwest = nwestt
        enddo
      enddo

      return
      end
c--------------------------------------------------------------------
      subroutine getpwlistall(ibox,bs,nboxes,nnbors,nbors,
     1           nchild,ichild,centers,isep,nuall,uall,ndall,dall,nnall,
     2           nall,nsall,sall,neall,eall,nwall,wall,nu1234,u1234,
     3           nd5678,d5678,nn1256,n1256,ns3478,s3478,ne1357,e1357,
     4           nw2468,w2468,nn12,n12,nn56,n56,ns34,s34,ns78,s78,ne13,
     5           e13,ne57,e57,nw24,w24,nw68,w68,ne1,e1,ne3,e3,ne5,e5,
     6           ne7,e7,nw2,w2,nw4,w4,nw6,w6,nw8,w8)
c-------------------------------------------------------------------
      implicit none
      integer ibox
      real *8 boxsize,bs
      integer nboxes,nnbors,nbors(nnbors)
      integer nchild, ichild(8,nboxes)
      real *8 centers(3,nboxes)
      integer isep
      integer nuall,ndall,nnall,nsall,neall,nwall,nu1234
      integer nd5678,nn1256,ns3478,ne1357,nw2468
      integer nn12,nn56,ns34,ns78,ne13,ne57,nw24,nw68
      integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
      integer uall(1),dall(1),nall(1),sall(1),eall(1),wall(1)
      integer u1234(1),d5678(1),n1256(1),s3478(1),e1357(1),w2468(1)
      integer n12(1),n56(1),s34(1),s78(1),e13(1),e57(1),w24(1),w68(1)
      integer e1(1),e3(1),e5(1),e7(1)
      integer w2(1),w4(1),w6(1),w8(1)

      integer jbox,kbox
      integer c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c16
      integer i,j

      nuall = 0
      ndall = 0
      nnall = 0
      nsall = 0
      neall = 0
      nwall = 0
      nu1234 = 0
      nd5678 = 0
      nn1256 = 0
      ns3478 = 0
      ne1357 = 0
      nw2468 = 0
      nn12 = 0
      nn56 = 0
      ns34 = 0
      ns78 = 0
      ne13 = 0
      ne57 = 0
      nw24 = 0
      nw68 = 0
      ne1 = 0
      ne3 = 0
      ne5 = 0
      ne7 = 0
      nw2 = 0
      nw4 = 0
      nw6 = 0
      nw8 = 0
      do i=1,nnbors
         jbox = nbors(i)
         do j=1,8
            kbox = ichild(j,jbox)
            if(kbox.gt.0) then
               c1 = 0
               c2 = 0
               c3 = 0
               c4 = 0
               c5 = 0
               c6 = 0
               c7 = 0
               c8 = 0
               c9 = 0
               c10 = 0
               c11 = 0
               c12 = 0
               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c1 = 1

               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c2 = 1

               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c3 = 1

               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c4 = 1

               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs+bs/2.0d0) c5 = 1

               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs-bs/2.0d0) c6 = 1

               if((centers(3,kbox)-centers(3,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c7 = 1

               if((centers(3,kbox)-centers(3,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c8 = 1

               if((centers(2,kbox)-centers(2,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c9 = 1

               if((centers(2,kbox)-centers(2,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c10 = 1

               if((centers(1,kbox)-centers(1,ibox)).ge.
     1              1.01d0*isep*bs-bs/2.0d0) c11 = 1

               if((centers(1,kbox)-centers(1,ibox)).le.
     1              -1.01d0*isep*bs+bs/2.0d0) c12 = 1

               if(c1.eq.1) then
                  nuall = nuall + 1
                  uall(nuall) = kbox
               endif

               if(c2.eq.1) then
                  ndall = ndall + 1
                  dall(ndall) = kbox
               endif

               if(c3.eq.1.and.c1.eq.0.and.c2.eq.0) then
                  nnall = nnall + 1
                  nall(nnall) = kbox
               endif

               if(c4.eq.1.and.c1.eq.0.and.c2.eq.0) then   
                  nsall = nsall + 1
                  sall(nsall) = kbox
               endif

               if(c5.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1             c4.eq.0) then
                  neall = neall + 1
                  eall(neall) = kbox
               endif

               if(c6.eq.1.and.c1.eq.0.and.c2.eq.0.and.c3.eq.0.and.
     1            c4.eq.0) then
                  nwall = nwall + 1
                  wall(nwall) = kbox
               endif

               c16 = c1 + c2 + c3 + c4 + c5 +c6
               if(c16.eq.0.and.c7.eq.1) then
                  nu1234 = nu1234 + 1
                  u1234(nu1234) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1) then
                  nd5678 = nd5678 + 1
                  d5678(nd5678) = kbox
               endif

               if(c16.eq.0.and.c9.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  nn1256 = nn1256 + 1
                  n1256(nn1256) = kbox
               endif

               if(c16.eq.0.and.c10.eq.1.and.c7.eq.0.and.c8.eq.0) then
                  ns3478 = ns3478 + 1
                  s3478(ns3478) = kbox
               endif

               if(c16.eq.0.and.c11.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  ne1357 = ne1357 + 1
                  e1357(ne1357) = kbox
               endif

               if(c16.eq.0.and.c12.eq.1.and.(c7+c8+c9+c10).eq.0) then
                  nw2468 = nw2468 + 1
                  w2468(nw2468) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c9.eq.1) then
                  nn12 = nn12 + 1
                  n12(nn12) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c9.eq.1) then
                  nn56 = nn56 + 1
                  n56(nn56) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c10.eq.1) then
                  ns34 = ns34 + 1
                  s34(ns34) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c10.eq.1) then
                  ns78 = ns78 + 1
                  s78(ns78) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne13 = ne13 + 1
                  e13(ne13) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c11.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  ne57 = ne57 + 1
                  e57(ne57) = kbox
               endif

               if(c16.eq.0.and.c8.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw24 = nw24 + 1
                  w24(nw24) = kbox
               endif

               if(c16.eq.0.and.c7.eq.1.and.c12.eq.1
     1           .and.c9.eq.0.and.c10.eq.0) then
                  nw68 = nw68 + 1
                  w68(nw68) = kbox
               endif

               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne1 = ne1 + 1
                  e1(ne1) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne3 = ne3 + 1
                  e3(ne3) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c11.eq.1) then
                  ne5 = ne5 + 1
                  e5(ne5) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c11.eq.1) then
                  ne7 = ne7 + 1
                  e7(ne7) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw2 = nw2 + 1
                  w2(nw2) = kbox
               endif
               if(c16.eq.0.and.c7.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw4 = nw4 + 1
                  w4(nw4) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c10.eq.1.and.c12.eq.1) then
                  nw6 = nw6 + 1
                  w6(nw6) = kbox
               endif
               if(c16.eq.0.and.c8.eq.0.and.c9.eq.1.and.c12.eq.1) then
                  nw8 = nw8 + 1
                  w8(nw8) = kbox
               endif
            endif
         enddo
      enddo
      

      return
      end
c--------------------------------------------------------------------      

       subroutine getaqmem(imode,nboxes,nlevels,boxsize,itree,ltree,
     1     ipointer,ilevel,mnbors,n,iarr,rad,flags,laqmem1,laqmem2)
c
c        Memory determination code for performing area queries
c        centered at particles with radius given by array rad.
c        
c        INPUT parameters
c        imode  - if imode = 1, area query is being performed
c                 on sources
c                 if imode = 2, area query is being performed
c                 on targets
c                 if imode = 3, area query is being performed
c                 on expansion centers
c
c        nboxes - number of boxes
c        nlevels - number of levels
c        boxsize - real *8 (0:nlevels)
c                  boxsize at each level
c        itree -  integer (ltree)
c                 Treeinfo
c        ltree - integer, length of tree
c        ipointer - integer(26)
c                   pointer for array itree
c
c        ilevel - integer(nboxes)
c                 ilevel(i) tells the level number of box number i
c
c        mnbors - maximum number of boxes in list1 aka nbors
c        n - integer
c            n=ns if imode =1 (number of sources)
c            n=nt if imode =2 (number of targets)
c            n=nexpc if imode=3 (number of expansion centers)
c
c        iarr - integer(n)
c               resorting array into tree boxes 
c        
c        rad - real *8(n)
c              rad(i) - area query parameter for particle i
c
c        flags - integer(n)
c                perform area query for particle i if 
c                flags(i).eq.1
c
c        Output
c        laqmem1 - number of integers required for number of hung 
c                  particles in each box
c
c        laqmem2 - number of integers required for number of
c                  relevant hung particles in each box

       implicit real *8 (a-h,o-z)
       integer nboxes,nlevels,itree(1),ltree,ipointer(1),ilevel(1)
       real *8 boxsize(0:nlevels),rad(1)
       integer, allocatable :: nlist(:),nlistr(:)
       integer nemax,nermax,dad,iarr(1),flags(*)

       allocate(nlist(nboxes),nlistr(nboxes))


       do 1000 i=1,nboxes

       nlist(i) = 0
       nlistr(i) = 0

 1000  continue       

       do 1050 i=1,n

       rad(i) = max(rad(i),1.0d-16*boxsize(nlevels))

 1050  continue       

       do 2000 ilev = 0,nlevels

       do 1900 ibox = itree(2*ilev+1),itree(2*ilev+2)

       nchild = itree(ipointer(3)+ibox-1)

       if(nchild.eq.0) then

       if(imode.eq.1) then

       istart = itree(ipointer(14)+ibox-1)
       iend = itree(ipointer(15)+ibox-1)

       endif


       if(imode.eq.2) then

       istart = itree(ipointer(16)+ibox-1)
       iend = itree(ipointer(17)+ibox-1)

       endif

       if(imode.eq.3) then

       istart = itree(ipointer(20)+ibox-1)
       iend = itree(ipointer(21)+ibox-1)

       endif

       do 1800 ie=istart,iend

       if(flags(iarr(ie)).eq.1) then

       rr = rad(iarr(ie))/boxsize(ilev)

       ilevup = int(log(rr)/log(2.0))+1
       ilevup = min(ilevup,ilev)
       itmp = ibox

       do 1700 i=1,ilevup

       itmp = itree(ipointer(2)+itmp-1)

 1700  continue     

       nlist(itmp) = nlist(itmp) + 1

       endif

 1800  continue

       endif

 1900  continue
 2000  continue       

       nlistr(1) = nlist(1)
       do 3000 ilev = 1,nlevels
       do 2900 ibox = itree(2*ilev+1),itree(2*ilev+2)

       dad = itree(ipointer(2)+ibox-1)
       nlistr(ibox) = nlistr(dad)

       nnbors = itree(ipointer(5)+ibox-1)

       do 2200 i=1,nnbors

       jbox = itree(ipointer(6)+mnbors*(ibox-1)+i-1)
       nlistr(ibox) = nlistr(ibox) + nlist(jbox)

 2200  continue

 2900  continue
 3000  continue       

       laqmem1 = 0
       laqmem2 = 0

       do 4000 i=1,nboxes

       laqmem1 = laqmem1 + nlist(i)
       laqmem2 = laqmem2 + nlistr(i)

 4000  continue       
       
       return
       end
    
c-----------------------------------------------------------------------    

       subroutine getaq(imode,nboxes,nlevels,boxsize,itree,ltree,
     1     ipointer,ilevel,mnbors,n,iarr,rad,flags,laqmem1,laqmem2,
     2     iaq,nlistr,iaqptr)
c
c        Memory determination code for performing area queries
c        centered at particles with radius given by array rad.
c        
c        INPUT parameters
c        imode  - if imode = 1, area query is being performed
c                 on sources
c                 if imode = 2, area query is being performed
c                 on targets
c                 if imode = 3, area query is being performed
c                 on expansion centers
c
c        nboxes - number of boxes
c        nlevels - number of levels
c        boxsize - real *8 (0:nlevels)
c                  boxsize at each level
c        itree -  integer (ltree)
c                 Treeinfo
c        ltree - integer, length of tree
c        ipointer - integer(26)
c                   pointer for array itree
c
c        ilevel - integer(nboxes)
c                 ilevel(i) tells the level number of box number i
c
c        mnbors - maximum number of boxes in list1 aka nbors
c        n - integer
c            n=ns if imode =1 (number of sources)
c            n=nt if imode =2 (number of targets)
c            n=nexpc if imode=3 (number of expansion centers)
c
c        iarr - integer(n)
c               resorting array into tree boxes 
c        
c        rad - real *8(n)
c              rad(i) - area query parameter for particle i
c
c        flags - integer(n)
c                perform area query for particle i if 
c                flags(i).eq.1
c
c        laqmem1 - length of temp array needed in the
c                  routine
c        laqmem2 - length of iaq array in output
c
c        Output
c        iaq - integer (laqmem)
c                 iaq(iaqptr(i):iaqptr(i)+nlistr(i)) is the list
c                 of particles relevant for box i
c
c        iaqptr - integer (nboxes)
c                 pointer array to iaq array
c
c        nlistr - integer (nboxes)
c                 nlistr(i) is the number of relevant particles
c                 in box i
c
c

       implicit real *8 (a-h,o-z)
       integer nboxes,nlevels,itree(1),ltree,ipointer(1),ilevel(1)
       real *8 boxsize(0:nlevels),rad(1)
       integer, allocatable :: nlist(:),iptr(:),iaqlist(:)
       integer nlistr(*),iaq(*),iaqptr(*)
       integer nemax,nermax,dad,iarr(1),flags(*)

       allocate(nlist(nboxes),iptr(nboxes),iaqlist(laqmem1))


cc       call prinf('laqmem1=*',laqmem1,1)
cc       call prinf('laqmem=*',laqmem2,1)

       do 1000 i=1,nboxes

       nlist(i) = 0
       nlistr(i) = 0

 1000  continue       


cc       call prinf('n=*',n,1)
cc       call prinf('itree=*',itree,12)
cc       call prinf('ltree=*',ltree,1)
cc       call  prin2('boxsize=*',boxsize,nlevels+1)

       do 1050 i=1,n

       rad(i) = max(rad(i),1.0d-16*boxsize(nlevels))

 1050  continue       

       do 2000 ilev = 0,nlevels
       do 1900 ibox = itree(2*ilev+1),itree(2*ilev+2)

       nchild = itree(ipointer(3)+ibox-1)

       if(nchild.eq.0) then

       if(imode.eq.1) then

       istart = itree(ipointer(14)+ibox-1)
       iend = itree(ipointer(15)+ibox-1)

       endif


       if(imode.eq.2) then

       istart = itree(ipointer(16)+ibox-1)
       iend = itree(ipointer(17)+ibox-1)

       endif

       if(imode.eq.3) then

       istart = itree(ipointer(20)+ibox-1)
       iend = itree(ipointer(21)+ibox-1)

       endif

       do 1800 ie=istart,iend

       if(flags(iarr(ie)).eq.1) then

       rr = rad(iarr(ie))/boxsize(ilev)
       ilevup = int(log(rr)/log(2.0))+1
       ilevup = min(ilevup,ilev)
       itmp = ibox

       do 1700 i=1,ilevup

       itmp = itree(ipointer(2)+itmp-1)

 1700  continue     

       nlist(itmp) = nlist(itmp) + 1

       endif

 1800  continue

       endif

 1900  continue
 2000  continue       

       nlistr(1) = nlist(1)
       do 3000 ilev = 1,nlevels
       do 2900 ibox = itree(2*ilev+1),itree(2*ilev+2)

       dad = itree(ipointer(2)+ibox-1)
       nlistr(ibox) = nlistr(dad)

       nnbors = itree(ipointer(5)+ibox-1)

       do 2200 i=1,nnbors

       jbox = itree(ipointer(6)+mnbors*(ibox-1)+i-1)
       nlistr(ibox) = nlistr(ibox) + nlist(jbox)

 2200  continue

 2900  continue
 3000  continue       


       iaqptr(1) = 1
       iptr(1) = 1

       do 4000 i=2,nboxes

       iptr(i) = iptr(i-1) + nlist(i-1)
       iaqptr(i) = iaqptr(i-1) + nlistr(i-1)

 4000  continue      

       do 4100 i=1,nboxes

       nlist(i) = 0
       nlistr(i) = 0

 4100  continue

       do 5000 ilev = 0,nlevels
       do 4900 ibox = itree(2*ilev+1),itree(2*ilev+2)

       nchild = itree(ipointer(3)+ibox-1)

       if(nchild.eq.0) then

       if(imode.eq.1) then

       istart = itree(ipointer(14)+ibox-1)
       iend = itree(ipointer(15)+ibox-1)

       endif


       if(imode.eq.2) then

       istart = itree(ipointer(16)+ibox-1)
       iend = itree(ipointer(17)+ibox-1)

       endif

       if(imode.eq.3) then

       istart = itree(ipointer(20)+ibox-1)
       iend = itree(ipointer(21)+ibox-1)

       endif

       do 4800 ie=istart,iend

       if(flags(iarr(ie)).eq.1) then

       rr = rad(iarr(ie))/boxsize(ilev)
       ilevup = int(log(rr)/log(2.0))+1
       ilevup = min(ilevup,ilev)
       itmp = ibox

       do 4700 i=1,ilevup

       itmp = itree(ipointer(2)+itmp-1)

 4700  continue     

       nlist(itmp) = nlist(itmp) + 1
       ictr = iptr(itmp)+nlist(itmp)-1
       iaqlist(ictr) = iarr(ie)


       endif

 4800  continue

       endif

 4900  continue
 5000  continue      

cc       call prinf('iaqlist=*',iaqlist,laqmem1)
   
cc       call prinf('nlist=*',nlist,nboxes)

       nlistr(1) = nlist(1)
       do 5050 i=1,nlist(1)

       iaq(i) = iaqlist(i)

 5050  continue
  

       do 6000 ilev = 1,nlevels
       do 5900 ibox = itree(2*ilev+1),itree(2*ilev+2)

       dad = itree(ipointer(2)+ibox-1)
       nlistr(ibox) = nlistr(dad)

       do 5075 i=1,nlistr(dad)

       ictr = iaqptr(ibox)+i-1
       ictr2 = iaqptr(dad)+i-1

       iaq(ictr) = iaq(ictr2)

 5075  continue      

       nnbors = itree(ipointer(5)+ibox-1)

       do 5200 i=1,nnbors

       jbox = itree(ipointer(6)+mnbors*(ibox-1)+i-1)
       do 5190 j=1,nlist(jbox)

       ictr = iaqptr(ibox)+nlistr(ibox)+j-1
       ictr2 = iptr(jbox)+j-1

       iaq(ictr) = iaqlist(ictr2)

 5190  continue

       nlistr(ibox) = nlistr(ibox) + nlist(jbox)

 5200  continue


 5900  continue
 6000  continue      

cc       call prinf('nlist=*',nlist,nboxes)
cc       call prinf('nlistr=*',nlistr,nboxes)
 


       return
       end
    
c-----------------------------------------------------------------------    
