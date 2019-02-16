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
c     multipole shift routines   F95 versions, using allocate
c
C***********************************************************************
      subroutine l3dmpmpquadu(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dmpmpquad0 (below).
C
C     Usage:
C
C     Shift center of multipole expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then shifts
C     along the Z-axis, and then rotates back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     :  scaling parameter for mpole expansion
C     x0y0z0  :  center of original multiple expansion
C     mpole   :  coefficients of original multiple expansion
C     nterms  :  order of multipole expansion
C     sc2     :  scaling parameter for shifted expansion
C     xnynzn  :  center of shifted expansion
C     nterms2 :  order of shifted expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen = coefficients of shifted mpole expansion
C     ier   = error return flag
C             CURRENTLY UNUSED
C
C     Notes: Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
C                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           fr      = work array 
C
C***********************************************************************
      implicit none
      integer  nterms,nterms2,ier,l,m,jnew,knew
      integer  ldc,imarray,lmarray,imarray1,lmarray1,iephi,lephi
      integer  ifr,lused
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
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
      lused = iephi + lephi
      allocate (cw(lused))
      allocate (w(2*ldc+3))
c
      call l3dmpmpquad0(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,cw(imarray),cw(imarray1),
     2           ldc,cw(iephi),w,ier)
      return
      end
c
c
C***********************************************************************
      subroutine l3dmpmpquadu_add(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,ldc,nterms2,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine l3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for mpole expansion
C     x0y0z0  : center of original multiple expansion
C     mpole   : coefficients of original multiple expansion
C     nterms  : order of multipole expansion
C     sc2     : scaling parameter for shifted expansion
C     xnynzn  : center of shifted expansion
C     ldc     : dimension parameter of shifted expansion
C     nterms2 : order of shifted expansion
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen  : coefficients of shifted mpole expansion
C     ier     : error return flag
C                     CURRENTLY UNUSED
C
C     Notes: Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           fr      = work array 
C
C***********************************************************************
      implicit none
      integer  nterms,nterms2,ier,l,m,jnew,knew,ldc
      real *8 x0y0z0(3),xnynzn(3)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:ldc,-ldc:ldc)
      complex *16 imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call l3dmpmpquadu(sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mptemp,nterms2,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            mpolen(l,m) = mpolen(l,m)+mptemp(l,m)
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l3dmpmpquad0(sc1,x0y0z0,mpole,nterms,sc2,
     1           xnynzn,mpolen,nterms2,marray,marray1,ldc,ephi,
     2           fr,ier)
C***********************************************************************
C
C     Usage:
C
C     Shift multipole expansion.
C     This is a reasonably fast "point and shoot" version which
C     first rotates the coordinate system, then doing the shifting
C     along the Z-axis, and then rotating back to the original
C     coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     sc1     : scaling parameter for mpole expansion
C     x0y0z0  : center of original multiple expansion
C     mpole   : coefficients of original multiple expansion
C     nterms  : order of multipole expansion
C     sc2     : scaling parameter for shifted expansion
C     xnynzn  : center of shifted expansion
C     nterms2 : order of shifted expansion
C     marray  : work array
C     marray1 : work array
C     ldc     : dimension parameter for work arrays
C     ephi    : work array
C     fr      : work array
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     mpolen  : coefficients of shifted expansion
C     ier     : error flag - UNUSED.
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           ldc      determines dimension of dc
c                    must exceed max(nterms,nterms2).
C           rd     = work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer  nterms, lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 rshift
      real *8 d,theta,ctheta,phi,sc1,sc2,rvec(3)
      real *8 fr(0:nterms+1)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 marray1(0:ldc,-ldc:ldc)
      complex *16 marray(0:ldc,-ldc:ldc)
c
      complex *16 ephi(-ldc-1:ldc+1),imag
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polarl(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c----- a rotation of THETA radians about the Yprime axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      mpole and ephi inside the following loop. 
c
      do l=0,nterms
         do m=-l,l
            marray1(l,m)=mpole(l,m)*ephi(m)
         enddo
      enddo
      do l=0,nterms2
         do m=-nterms2,nterms2
            marray(l,m)= 0.0d0
         enddo
      enddo
      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,
     1        ldc,marray,ldc)
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call l3dmpmpzshift(sc1,marray,ldc,nterms,sc2,mpolen,
     1           nterms2,nterms2,rshift,fr)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine l3dmpmpzshift(scale,mpole,lmp,nterms,scale2,mpolen,
     1      lmpn,nterms2,zshift,fr)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a multipole expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the shifted expansion.
c
c     INPUT:
c
c     scale    : scale parameter for mpole
c     mpole    : coefficients of original multipole exp.
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in the orig. expansion
c     scale2   : scale parameter for new expansion (mpolen)
c     lmpn     : leading dim of shifted (may be work array)
c     nterms2  : number of terms in output expansion
c     zshift   : shifting distance along z-axis
c                              (always assumed positive)
c     fr       : work array
c
c     OUTPUT:
c
c     mpolen  (complex *16)  : coefficients of shifted exp.
c
c-----------------------------------------------------------------------
      implicit none
      integer l0,nmax,nterms,nterms2,nquad,ier,lmp,lmpn,ldc
      integer l,m,jnew,knew
      real *8 d,zshift,scale,scale2,rat
      real *8 fr(0:*)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 mpolen(0:lmpn,-lmpn:lmpn)
      real *8, allocatable :: dc(:,:)
      real *8, allocatable :: carray(:,:)
C
C----- shift along z-axis
C
      nmax = max(nterms,nterms2)
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
      do l=2,nmax+1
         fr(l) = fr(l-1)*d
      enddo
c
      do jnew = 0,nterms2
         l0 = max(0,jnew-nterms)
         do knew = -jnew,jnew
	    mpolen(jnew,knew) = 0.0d0
	    do l=l0,jnew-iabs(knew)
	       mpolen(jnew,knew)=mpolen(jnew,knew)+mpole(jnew-l,knew)*
     1             fr(l)*dc(jnew-knew,l)*dc(jnew+knew,l)*(-1)**l
            enddo
         enddo
      enddo
      do jnew = 1,nterms2
         do knew = -jnew,jnew
	    mpolen(jnew,knew)=mpolen(jnew,knew)*(scale2/scale)**jnew
         enddo
      enddo
      return
      end
C
