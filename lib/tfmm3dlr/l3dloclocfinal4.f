cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-01-05 12:57:52 -0500 (Tue, 05 Jan 2010) $
c    $Revision: 782 $
c
c
c     Local to local shift routines, f90 version using allocate
c
C***********************************************************************
      subroutine l3dloclocquadu(sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dloclocquad0 (below).
C
C     Usage:
C
C     Shift center of a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts along
C     the Z-axis, and then rotates back to the original
C     coordinates.
C
C     INPUT:
C
C     sc1       scaling parameter for locold expansion
C     x0y0z0    center of original expansion
C     locold    coefficients of original expansion
C     nterms    order of original expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted expansion
C     nterms2   order of shifted expansion
C
C     OUTPUT:
C
C     local = coefficients of shifted expansion
C     ier   = error return flag
C             CURRENTLY NOT USED
C
C***********************************************************************
      implicit none
      integer nterms,nterms2,ier,l,m,jnew,knew
      integer  ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      integer  lused
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2,d,theta,phi,ctheta
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
      complex *16, allocatable :: cw(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      imarray = 1
      lmarray = (ldc+1)*(2*ldc+1) 
      imarray1 = imarray+lmarray
      lmarray1 = (ldc+1)*(2*ldc+1) 
      iephi = imarray1+lmarray1
      lephi = (2*ldc+3)
      lused = iephi+ lephi
      allocate (cw(lused))
      allocate (w(2*ldc+3))
c
      call l3dloclocquad0(sc1,x0y0z0,locold,nterms,sc2,xnynzn,
     1           local,nterms2,cw(imarray),cw(imarray1),ldc,
     2           cw(iephi),w,ier)
      return
      end
c
c
C***********************************************************************
      subroutine l3dloclocquadu_add(sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,local,ldc,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dloclocquad0 (below).
C
C     Usage:
C
C     Shift center of a local expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts along
C     the Z-axis, and then rotates back to the original
C     coordinates.
C
C     INPUT:
C
C     sc1       scaling parameter for locold expansion
C     x0y0z0    center of original expansion
C     locold    coefficients of original expansion
C     nterms    order of original expansion
C     sc2       scaling parameter for local expansion
C     xnynzn    center of shifted expansion
C     ldc       dimension parameter for local expansion
C     nterms2   order of shifted expansion
C
C     OUTPUT:
C
C     local = coefficients of shifted expansion
C     ier   = error return flag
C                   CURRENTLY NOT USED
C***********************************************************************
      implicit none
      integer nterms,nterms2,ldc,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2,d,theta,phi,ctheta
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 local(0:ldc,-ldc:ldc)
      complex *16 imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dloclocquadu(sc1,x0y0z0,locold,nterms,
     1           sc2,xnynzn,mptemp,nterms2,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            local(l,m) = local(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dloclocquad0(sc1,x0y0z0,locold,nterms,sc2,
     1           xnynzn,local,nterms2,marray,marray1,
     2           ldc,ephi,fr,ier) 
C***********************************************************************
C
C     Usage:
C
C           Shifts center of a local expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for locold expansion
C     x0y0z0  : center of original multiple expansion
C     locold  : coefficients of original multiple expansion
C     nterms  : order of original local expansion
C     sc2     : scaling parameter for local expansion
C     xnynzn  : center of shifted local expansion
C     nterms2 : order of new local expansion
c     marray  : work array
c     marray1 : work array
c     ldc     : dimension parameter for work arrays
c     ephi    : work array 
c     fr      : work array 
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     local   : coefficients of shifted local expansion
c     ier     : error return code 
c               CURRENTLY NOT USED
C
C***********************************************************************
C
      implicit none
      integer nterms,ier,l,m,jnew,knew,nterms2,ldc,mp
      real *8 x0y0z0(3),xnynzn(3),rvec(3)
      real *8 d,theta,ctheta,phi,sc1,sc2
      real *8 fr(0:ldc+1)
      real *8 rshift
      complex *16 locold(0:nterms,-nterms:nterms)
      complex *16 marray1(0:ldc,-ldc:ldc)
      complex *16 local(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 imag,ephi1
      complex *16 ephi(-ldc-1:ldc+1)
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
c
      ephi1 = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c      a rotation of THETA radians about the Yprime-axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      locold and ephi inside the following loop. 
c
      do l=0,nterms
         do mp=-l,l
            marray1(l,mp) = locold(l,mp)*ephi(mp)
         enddo
      enddo
      do l=0,nterms2
         do mp=-nterms2,nterms2
            marray(l,mp) = 0.0d0
         enddo
      enddo
ccc      t1 = second()
      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,
     1      ldc,marray,ldc)
ccc      t2 = second()
c
c----- shift the local expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
ccc      t1 = second()
       call l3dlocloczshift(sc1,marray,ldc,nterms,sc2,local,
     1           nterms2,nterms2,rshift,fr,ier) 
ccc      t2 = second()
c
c      reverse THETA rotation.
c      I.e. rotation of -THETA about Yprime axis.
c
ccc      t1 = second()
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,local,
     1      nterms2,marray,ldc)
ccc      t2 = second()
ccc      call prin2(' time for second rot is *',t2-t1,1)
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            local(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine l3dlocloczshift(scale,locold,lmp,nterms,scale2,
     1  local,lmpn,nterms2,zshift,fr,ier) 
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a local expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the local expansion.
c
c     INPUT:
c
c     scale    : scaling parameter for locold
c     locold   : coefficients of original multipole exp.
c     lmp      : leading dim of locold (may be a work array)
c     nterms   : number of terms in the orig. expansion
c
c     scale2   : scaling parameter for output expansion (local)
c     lmpn     : leading dim of local (may be a work array)
c     nterms2  : number of terms in output local exp.
c     zshift   : shifting distance along z-axis (assumed positive)
c     fr       : work array
c
c     OUTPUT:
c
c     local    : coefficients of shifted local exp.
c     ier      : error return code
c                 CURRENTLY NOT USED
c-----------------------------------------------------------------------
      implicit none
      integer nmax,nterms,nterms2,nquad,ier,l,lw,m,jnew,knew
      integer lmp,lmpn,ll
      real *8  zshift,d
      real *8  scale,scale2
      real *8 fr(0:nterms+1)
      complex *16 locold(0:lmp,-lmp:lmp)
      complex *16 local(0:lmpn,-lmpn:lmpn)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: carray(:,:)
C
C----- shift along z-axis 
C
      nmax = nterms+nterms2
c
      allocate( dc(0:2*nmax,0:2*nmax) )
      allocate( carray(0:2*nmax,0:2*nmax) )
c
      do l = 0,2*nmax
         carray(l,0) = 1.0d0
         dc(l,0) = 1.0d0
      enddo
      do m=1,2*nmax
         carray(m,m) = 1.0d0
         dc(m,m) = 1.0d0
         do l=m+1,2*nmax
	    carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
	    dc(l,m)=sqrt(carray(l,m))
         enddo
      enddo
c
      d = zshift
      fr(0) = 1.0d0
      d = d*scale
      fr(1) = d
      do l=2,nterms+1
         fr(l) = fr(l-1)*d
      enddo
c
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew) = locold(jnew,knew)
	    do l=1,nterms-jnew
	       ll = l+jnew
	       local(jnew,knew)=local(jnew,knew)+locold(ll,knew)*
     1              fr(l)*dc(ll+knew,l)*dc(ll-knew,l)
            enddo
         enddo
      enddo
      do jnew = 0,nterms2
         do knew = -jnew,jnew
	    local(jnew,knew)=local(jnew,knew)*(scale/scale2)**jnew
         enddo
      enddo
      return
      end
C
