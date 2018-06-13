      subroutine tfmm3dwrap(ier,iprec,source,dipvec,ns,wts,ifcharge,
     1           sigma,ifdipole,mu,targets,nt,trads,stdev,
     2           stdev_grad,ifpottarg,pottarg,iffldtarg,fldtarg,tfmm)

      implicit real *8(a-h,o-z)
c-----------------------------------------------------------------
c     This subroutine computes the potential and gradient due to 
c     single and double layer sources for the surface smoothing kernel
c     in 3D. 
c
c     Input:
c     ier     - error return code
c     iprec   - precision flag
c     source  - locations of the points
c     dipvec  - normal vectors
c     ns      - number of discretization points 
c     wts     - quad weights
c     ifcharge - flag =1 if SLP component is present, 0 otherwise
c     sigma   - single layer density at the corresponding points
c     ifdipole - flag =1 if DLP component is present, 0 otherwise
c     mu      - double layer density at the corresponding points
c     targets - target locations
c     nt      - number of targets
c     trads   - target radii (where kernel differs from Coulomb kernel)
c     stdev   - smoothing kernel parameter
c     stdev_grad   - grad of smoothing kernel parameter
c     ifpottarg  - flag=1 means compute pot at targets, 0 otherwise
c     iffldtarg - flag=1 means compute grad at targets, 0 otherwise
c              
c     Output:
c
c     pottarg     - potential at target points
c     fldtarg    - gradient at target points
c     tfmm    - time for FMM call
c--------------------------------------------------------------------
c
      integer ns,nt
      integer ifpottarg,ifgradtarg
      real *8 source(3,ns)
      real *8 dipvec(3,ns)
      real *8 wts(ns)
      real *8 targets(3,nt)
      real *8 trads(nt)
      complex *16 sigma(ns),mu(ns)
      complex *16 pottarg(nt),fldtarg(3,nt)

c     Tree variables

      integer, allocatable :: itree(:)
      real *8, allocatable :: boxsize(:)
      real *8, allocatable :: treecenters(:,:)

      real *8, allocatable :: radsrc(:)

      integer ipointer(26)

c     Temp variables
c
      real *8 expc(3,1)
      integer flags(10)
c
      pi = 4.0d0*datan(1.0d0)
      eye = dcmplx(0.0d0,1.0d0)
c
c     Set up variables of tree code
c
      idivflag = 0
      isep = 1
      if(isep.eq.1) ndiv = 1000
      if(isep.eq.2) ndiv = 200

      nlmax = 200
      nlevels = 0
      nboxes = 0
      mhung = 0
      ltree = 0

      allocate(radsrc(ns))

      do i=1,ns
         radsrc(i) = 0
      enddo

      nexpc = 0

      call mkptreemem(source,ns,radsrc,expc,nexpc,targets,nt,
     1                  trads,idivflag,ndiv,isep,nlmax,nlevels,nboxes,
     2                  mhung,ltree)

      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(treecenters(3,nboxes))

cc      call prinf('ltree=*',ltree,1)
cc      call prinf('nboxes=*',nboxes,1)

c       Call tree code
        call mkptree(ier,source,ns,radsrc,expc,
     1               nexpc,targets,nt,trads,idivflag,ndiv,isep,
     2               mhung,nlevels,nboxes,treecenters,boxsize,itree,
     3               ltree,ipointer)

c       Find nudnsew
        if(isep.eq.1) mnlist2 = 27*7
        if(isep.eq.2) mnlist2 = 125*7
c
      if (ifdipole.eq.1) then
         do i = 1,ns
            mu(i) = mu(i)*wts(i)
         enddo
      endif
      if (ifcharge.eq.1) then
         do i = 1,ns
            sigma(i) = sigma(i)*wts(i)
         enddo
      endif
c
        t1=second()
C$        t1=omp_get_wtime()
c
        call tfmm3d(ier,iprec,
     $     ns,source,
     $     ifcharge,sigma,ifdipole,mu,dipvec,
     $     nt,targets,trads,stdev,stdev_grad,
     $     itree,ltree,ipointer,isep,ndiv,nlevels,nboxes,boxsize,
     $     treecenters,
     $     ifpottarg,pottarg,iffldtarg,fldtarg)

        t2=second()
        tfmm=t2-t1
c
      return
      end
