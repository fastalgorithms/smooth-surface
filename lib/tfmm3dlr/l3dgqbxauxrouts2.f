c***********************************************************************
      subroutine rlscini2(rlsc,nlambs,rlams,nterms)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      real *8 rlsc(0:nterms,0:nterms,nlambs)
      real *8     rlams(nlambs),rlampow(0:100)
      real *8     facts(0:200)
c
      facts(0) = 1.0d0
      do 100 i = 1,100
	 facts(i) = facts(i-1)*dsqrt(i+0.0d0)
100   continue
c
      do 1000 nl = 1,nlambs
c
c     compute powers of lambda_nl
c
         rlampow(0) = 1.0d0
         rmul = rlams(nl)
         do 200 j = 1,nterms
            rlampow(j) = rlampow(j-1)*rmul
200      continue
         do 600 j = 0,nterms
            do 400 k = 0,j
               rlsc(j,k,nl) = rlampow(j)/(facts(j-k)*facts(j+k))
400         continue
600      continue
1000  continue
      return
      end
c***********************************************************************
      subroutine mkexps2(rlams,nlambs,numphys,nexptotp,xs,ys,zs)
      implicit real *8 (a-h,o-z)
      complex *16 ima
      complex *16 xs(-5:5,nexptotp)
      complex *16 ys(-5:5,nexptotp)
      real *8 zs(5,nexptotp)
      real *8     rlams(nlambs),u
      integer *4  nlambs,numphys(nlambs),nexptotp
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine computes the tables of exponentials needed
c     for translating exponential representations of harmonic
c     functions, discretized via normans quadratures.
c
c     u   = \int_0^\infty e^{-\lambda z}
c           \int_0^{2\pi} e^{i\lambda(x cos(u)+y sin(u))}
c           mexpphys(lambda,u) du dlambda
c
c     mexpphys(*):  discrete values of the moment function 
c                   m(\lambda,u), ordered as follows.
c
c         mexpphys(1),...,mexpphys(numphys(1)) = m(\lambda_1,0),..., 
c              m(\lambda_1, 2*pi*(numphys(1)-1)/numphys(1)).
c         mexpphys(numphys(1)+1),...,mexpphys(numphys(2)) = 
c              m(\lambda_2,0),...,
c                  m(\lambda_2, 2*pi*(numphys(2)-1)/numphys(2)).
c         etc.
c
c     on input:
c
c     rlams(nlambs)   discretization points in lambda integral 
c     nlambs          number of discret. pts. in lambda integral
c     numphys(j)     number of nodes in u integral needed 
c                    for corresponding lambda =  lambda_j. 
c     nexptotp        sum_j numphys(j)
c
c     on output:
c
c     xs(1,nexptotp)   e^{i \lambda_j (cos(u_k)}  in above ordering
c     xs(2,nexptotp)   e^{i \lambda_j (2 cos(u_k)}  in above ordering.
c     xs(3,nexptotp)   e^{i \lambda_j (3 cos(u_k)}  in above ordering.
c     xs(4,nexptotp)   e^{i \lambda_j (4 cos(u_k)}  in above ordering.
c     xs(5,nexptotp)   e^{i \lambda_j (5 cos(u_k)}  in above ordering.
c     ys(1,nexptotp)   e^{i \lambda_j (sin(u_k)}  in above ordering.
c     ys(2,nexptotp)   e^{i \lambda_j (2 sin(u_k)}  in above ordering.
c     ys(3,nexptotp)   e^{i \lambda_j (3 sin(u_k)}  in above ordering.
c     ys(4,nexptotp)   e^{i \lambda_j (4 sin(u_k)}  in above ordering.
c     ys(5,nexptotp)   e^{i \lambda_j (5 sin(u_k)}  in above ordering.
c     zs(1,nexptotp)   e^{-\lambda_j}     in above ordering.
c     zs(2,nexptotp)    e^{-2 \lambda_j}   in above ordering. 
c     zs(3,nexptotp)    e^{-3 \lambda_j}   in above ordering. 
c     zs(4,nexptotp)    e^{-2 \lambda_j}   in above ordering. 
c     zs(5,nexptotp)    e^{-3 \lambda_j}   in above ordering. 
c------------------------------------------------------------
c      
c     loop over each lambda value 
c
      pi = 4*datan(1.0d0)
      ntot = 0
      do 400 nl = 1,nlambs
         hu=2*pi/numphys(nl)
         do 200 mth = 1,numphys(nl)
            u = (mth-1)*hu
            ncurrent = ntot+mth
            zs(1,ncurrent) = dexp( -rlams(nl) )
            zs(2,ncurrent) = dexp( - 2.0d0*rlams(nl) )
            zs(3,ncurrent) = dexp( - 3.0d0*rlams(nl) )
            zs(4,ncurrent) = dexp( - 4.0d0*rlams(nl) )
            zs(5,ncurrent) = dexp( - 5.0d0*rlams(nl) )
            xs(1,ncurrent) = cdexp(ima*rlams(nl)*dcos(u))
            xs(2,ncurrent) = cdexp(ima*rlams(nl)*2.0d0*dcos(u))
            xs(3,ncurrent) = cdexp(ima*rlams(nl)*3.0d0*dcos(u))
            xs(4,ncurrent) = cdexp(ima*rlams(nl)*4.0d0*dcos(u))
            xs(5,ncurrent) = cdexp(ima*rlams(nl)*5.0d0*dcos(u))
            ys(1,ncurrent) = cdexp(ima*rlams(nl)*dsin(u))
            ys(2,ncurrent) = cdexp(ima*rlams(nl)*2.0d0*dsin(u))
            ys(3,ncurrent) = cdexp(ima*rlams(nl)*3.0d0*dsin(u))
            ys(4,ncurrent) = cdexp(ima*rlams(nl)*4.0d0*dsin(u))
            ys(5,ncurrent) = cdexp(ima*rlams(nl)*5.0d0*dsin(u))

            xs(0,ncurrent)  = 1.0d0
            xs(-1,ncurrent) = cdexp(-ima*rlams(nl)*dcos(u))
            xs(-2,ncurrent) = cdexp(-ima*rlams(nl)*2.0d0*dcos(u))
            xs(-3,ncurrent) = cdexp(-ima*rlams(nl)*3.0d0*dcos(u))
            xs(-4,ncurrent) = cdexp(-ima*rlams(nl)*4.0d0*dcos(u))
            xs(-5,ncurrent) = cdexp(-ima*rlams(nl)*5.0d0*dcos(u))

            ys(0,ncurrent) = 1.0d0
            ys(-1,ncurrent) = cdexp(-ima*rlams(nl)*dsin(u))
            ys(-2,ncurrent) = cdexp(-ima*rlams(nl)*2.0d0*dsin(u))
            ys(-3,ncurrent) = cdexp(-ima*rlams(nl)*3.0d0*dsin(u))
            ys(-4,ncurrent) = cdexp(-ima*rlams(nl)*4.0d0*dsin(u))
            ys(-5,ncurrent) = cdexp(-ima*rlams(nl)*5.0d0*dsin(u))
200      continue
         ntot = ntot+numphys(nl)
400   continue
      return
      end
c***********************************************************************
      subroutine mkfexp2(nlambs,numfour,numphys,fexpe,fexpo,fexpback)
      implicit real *8 (a-h,o-z)
      complex *16 ima
      complex *16 fexpe(1)
      complex *16 fexpo(1)
      complex *16 fexpback(1)
      integer *4  nlambs,numphys(nlambs),numfour(nlambs)
      data ima/(0.0d0,1.0d0)/
c
c     this subroutine computes the tables of exponentials needed
c     for mapping from fourier to physical domain. 
c     in order to minimize storage, they are organized in a 
c     one-dimenional array corresponding to the order in which they
c     are accessed by subroutine ftophys.
c    
c     size of fexpe, fexpo =          40000   for nlambs = 39
c     size of fexpe, fexpo =          15000   for nlambs = 30
c     size of fexpe, fexpo =           4000   for nlambs = 20
c     size of fexpe, fexpo =            400   for nlambs = 10
c
c***********************************************************************
      pi = 4*datan(1.0d0)
      nexte = 1
      nexto = 1
      do 600 i=1,nlambs
	 nalpha = numphys(i)
         halpha=2*pi/nalpha
         do 400 j=1,nalpha
            alpha=(j-1)*halpha
	    do 200 mm = 2,numfour(i),2
               fexpe(nexte)  = cdexp(ima*(mm-1)*alpha)
	       nexte = nexte + 1
200         continue
	    do 300 mm = 3,numfour(i),2
               fexpo(nexto)  = cdexp(ima*(mm-1)*alpha)
	       nexto = nexto + 1
300         continue
400      continue
600   continue
      next = 1
      do 1600 i=1,nlambs
	 nalpha = numphys(i)
         halpha=2*pi/nalpha
	 do 1400 mm = 2,numfour(i)
            do 1200 j=1,nalpha
               alpha=(j-1)*halpha
               fexpback(next)  = cdexp(-ima*(mm-1)*alpha)
	       next = next + 1
1200        continue
1400     continue
1600  continue
      return
      end
c***********************************************************************

      subroutine mkudexp2(ibox,ilev,nboxes,ichild,rscale,nterms,nmax,
     1           rmlexp,iaddr,rlams,nlams,nfourier,nphysical,nthmax,
     2           nexptot,nexptotp,mexpup,mexpdown,mexpupphys,
     3           mexpdownphys,mexpuall,mexpu1234,mexpdall,mexpd5678,xs,
     4           ys,zs,fexpe,fexpo,rlsc)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       real *8 rscale,rmlexp(*),rlams(nlams)
       complex *16 tmp(0:nmax,-nmax:nmax)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpuall(nexptotp),mexpdall(nexptotp)
       complex *16 mexpu1234(nexptotp),mexpd5678(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       real *8 zs(5,nexptotp),rlsc(0:nmax,0:nmax,nlams)
       complex *16 fexpe(1),fexpo(1)

c      temp variables
       integer jbox,ctr,ii,jj,i
       complex *16 ztmp

c      add contributions due to child 1
       jbox = ichild(1,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif

       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif
       do i=1,nexptotp
          mexpu1234(i) = mexpupphys(i)
          mexpdall(i) = mexpdownphys(i)
       enddo

c      add contributions due to child 1
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpupphys(i) = mexpupphys(i)*dconjg(xs(1,i))
          mexpu1234(i) = mexpu1234(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*xs(1,i)
          mexpdall(i) = mexpdall(i)+mexpdownphys(i)
       enddo
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpupphys(i) = mexpupphys(i)*dconjg(ys(1,i))
          mexpu1234(i) = mexpu1234(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*ys(1,i)
          mexpdall(i) = mexpdall(i)+mexpdownphys(i)
       enddo

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*ys(1,i)
          mexpupphys(i) = mexpupphys(i)*dconjg(ztmp)
          mexpu1234(i) = mexpu1234(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*ztmp
          mexpdall(i) = mexpdall(i)+mexpdownphys(i)
       enddo

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif

       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpupphys(i) = mexpupphys(i)/zs(1,i)
          mexpuall(i) = mexpu1234(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*zs(1,i)
          mexpd5678(i) = mexpdownphys(i)
       enddo

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          mexpupphys(i) = mexpupphys(i)/ztmp
          mexpuall(i) = mexpuall(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*ztmp
          mexpd5678(i) = mexpd5678(i)+mexpdownphys(i)
       enddo

c      add contributions due to child 7
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = ys(1,i)*zs(1,i)
          mexpupphys(i) = mexpupphys(i)/ztmp
          mexpuall(i) = mexpuall(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*ztmp
          mexpd5678(i) = mexpd5678(i)+mexpdownphys(i)
       enddo

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo

c         convert multipole to plane waves
          call mpoletoexp2(tmp,nterms,nlams,nfourier,nexptot,
     1         mexpup,mexpdown,rlsc)

          call ftophys(mexpup,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpupphys,fexpe,fexpo)

          call ftophys(mexpdown,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpdownphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpupphys(i) = 0.0d0
             mexpdownphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*ys(1,i)*zs(1,i)
          mexpupphys(i) = mexpupphys(i)/ztmp
          mexpuall(i) = mexpuall(i)+mexpupphys(i)
          mexpdownphys(i) = mexpdownphys(i)*ztmp
          mexpd5678(i) = mexpd5678(i)+mexpdownphys(i)
          mexpdall(i) = mexpdall(i)+mexpd5678(i)
       enddo

      return
      end
c--------------------------------------------------------------------
      subroutine processup2(ibox,rscale,lexp,centers,nboxes,nuall,uall,
     1           nu1234,u1234,mexpuall,mexpu1234,xs,ys,zs,nexptotp)
c     shift plane waves for uall and u1234

      integer ibox
      real *8 rscale,zs(5,nexptotp),centers(3,nboxes)
      integer nuall,nu1234,uall(*),u1234(*)
      complex *16 lexp(nexptotp,nboxes,6)
      complex *16 mexpuall(nexptotp),mexpu1234(nexptotp)
      complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      complex *16 zmul

c     temp variables
      integer i,jbox,j,ix,iy,iz
      real *8 ctmp(3)
      
   
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
      do i=1,nuall
         jbox = uall(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(iz,j)
            if(ix.gt.0) zmul = zmul*xs(ix,j)
            if(ix.lt.0) zmul = zmul*dconjg(xs(-ix,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,1) = lexp(j,jbox,1)+mexpuall(j)*zmul
         enddo
      enddo

      do i=1,nu1234
         jbox = u1234(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(iz,j)
            if(ix.gt.0) zmul = zmul*xs(ix,j)
            if(ix.lt.0) zmul = zmul*dconjg(xs(-ix,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,1) = lexp(j,jbox,1)+mexpu1234(j)*zmul
         enddo
      enddo

      return
      end
c------------------------------------------------------------------      
      subroutine processdn2(ibox,rscale,lexp,centers,nboxes,ndall,dall,
     1           nd5678,d5678,mexpdall,mexpd5678,xs,ys,zs,nexptotp)
c     shift plane waves for uall and u1234

      integer ibox
      real *8 rscale,zs(5,nexptotp),centers(3,*)
      integer ndall,nd5678,dall(*),d5678(*)
      complex *16 lexp(nexptotp,nboxes,6)
      complex *16 mexpdall(nexptotp),mexpd5678(nexptotp)
      complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      complex *16 zmul

c     temp variables
      integer i,jbox,j,ix,iy,iz
      real *8 ctmp(3)
      
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
      do i=1,ndall
         jbox = dall(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(-iz,j)
            if(ix.lt.0) zmul = zmul*xs(-ix,j)
            if(ix.gt.0) zmul = zmul*dconjg(xs(ix,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            lexp(j,jbox,2) = lexp(j,jbox,2)+mexpdall(j)*zmul
         enddo
      enddo

      do i=1,nd5678
         jbox = d5678(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(-iz,j)
            if(ix.lt.0) zmul = zmul*xs(-ix,j)
            if(ix.gt.0) zmul = zmul*dconjg(xs(ix,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            lexp(j,jbox,2) = lexp(j,jbox,2)+mexpd5678(j)*zmul
         enddo
      enddo

      return
      end
c------------------------------------------------------------------      
      subroutine mknsexp2(ibox,ilev,nboxes,ichild,rscale,nterms,nmax,
     1           rmlexp,iaddr,mptemp,mptemp2,rlams,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexpnof,mexpsof,
     2           mexpnphys,mexpsphys,rdminus,mexpnall,mexpn1256,mexpn12,
     4           mexpn56,mexpsall,mexps3478,mexps34,mexps78,xs,
     4           ys,zs,fexpe,fexpo,rlsc)
c--------------------------------------------------------------------
c      create north south expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       real *8 rscale,rmlexp(*),rlams(nlams)
       complex *16 tmp(0:nmax,-nmax:nmax),mptemp(0:nmax,-nmax:nmax)
       complex *16 mptemp2(0:nmax,-nmax:nmax)
       real *8 rdminus(0:nmax,0:nmax,-nmax:nmax)
       complex *16 mexpnof(nexptot),mexpsof(nexptot)
       complex *16 mexpnphys(nexptotp),mexpsphys(nexptotp)
       complex *16 mexpnall(nexptotp),mexpsall(nexptotp)
       complex *16 mexpn1256(nexptotp),mexps3478(nexptotp)
       complex *16 mexpn12(nexptotp),mexpn56(nexptotp)
       complex *16 mexps34(nexptotp),mexps78(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       real *8 zs(5,nexptotp),rlsc(0:nmax,0:nmax,nmax)
       complex *16 fexpe(1),fexpo(1)

c      temp variables
       integer jbox,ctr,ii,jj,i
       complex *16 ztmp

c      add contributions due to child 1
       jbox = ichild(1,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpn12(i) = mexpnphys(i)
          mexpsall(i) = mexpsphys(i)
       enddo

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpnphys(i) = mexpnphys(i)*dconjg(ys(1,i))
          mexpn12(i) = mexpn12(i)+mexpnphys(i)
          mexpsphys(i) = mexpsphys(i)*ys(1,i)
          mexpsall(i) = mexpsall(i)+mexpsphys(i)
       enddo
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpnphys(i) = mexpnphys(i)/zs(1,i)
          mexpnall(i) = mexpnphys(i)
          mexpsphys(i) = mexpsphys(i)*zs(1,i)
          mexps34(i) = mexpsphys(i)
       enddo

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          mexpnphys(i) = mexpnphys(i)/ztmp
          mexpnall(i) = mexpnall(i)+mexpnphys(i)
          mexpsphys(i) = mexpsphys(i)*ztmp
          mexps34(i) = mexps34(i)+mexpsphys(i)
       enddo

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpnphys(i) = mexpnphys(i)*dconjg(xs(1,i))
          mexpn56(i) = mexpnphys(i)
          mexpsphys(i) = mexpsphys(i)*xs(1,i)
          mexpsall(i) = mexpsall(i)+mexpsphys(i)
       enddo

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*ys(1,i)
          mexpnphys(i) = mexpnphys(i)*dconjg(ztmp)
          mexpn56(i) = mexpn56(i)+mexpnphys(i)
          mexpn1256(i) = mexpn56(i) + mexpn12(i)
          mexpsphys(i) = mexpsphys(i)*ztmp
          mexpsall(i) = mexpsall(i)+mexpsphys(i)
       enddo

c      add contributions due to child 7
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*zs(1,i)
          mexpnphys(i) = mexpnphys(i)/ztmp
          mexpnall(i) = mexpnall(i)+mexpnphys(i)
          mexpsphys(i) = mexpsphys(i)*ztmp
          mexps78(i) = mexpsphys(i)
       enddo

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztoy(nterms,tmp,mptemp,rdminus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpnof,mexpsof,rlsc)

          call ftophys(mexpnof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpnphys,fexpe,fexpo)

          call ftophys(mexpsof,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpsphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpnphys(i) = 0.0d0
             mexpsphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*ys(1,i)*zs(1,i)
          mexpnphys(i) = mexpnphys(i)/ztmp
          mexpnall(i) = mexpnall(i)+mexpnphys(i)+mexpn1256(i)
          mexpsphys(i) = mexpsphys(i)*ztmp
          mexps78(i) = mexps78(i) + mexpsphys(i)
          mexps3478(i) = mexps78(i) + mexps34(i)
          mexpsall(i) = mexpsall(i) + mexps3478(i)
       enddo

      return
      end
c--------------------------------------------------------------------      
      subroutine processno2(ibox,rscale,lexp,centers,nboxes,nnall,nall,
     1           nn1256,n1256,nn12,n12,nn56,n56,mexpnall,mexpn1256,
     2           mexpn12,mexpn56,xs,ys,zs,nexptotp)
c     shift plane waves for uall and u1234

      integer ibox
      real *8 rscale,zs(5,nexptotp),centers(3,*)
      integer nnall,nn1256,nn12,nn56,nall(*),n1256(*),n12(*),n56(*)
      complex *16 lexp(nexptotp,nboxes,6)
      complex *16 mexpnall(nexptotp),mexpn1256(nexptotp)
      complex *16 mexpn12(nexptotp),mexpn56(nexptotp)
      complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      complex *16 zmul

c     temp variables
      integer i,jbox,j,ix,iy,iz
      real *8 ctmp(3)
      
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
      do i=1,nnall
         jbox = nall(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(iy,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(ix.gt.0) zmul = zmul*ys(ix,j)
            if(ix.lt.0) zmul = zmul*dconjg(ys(-ix,j))
            lexp(j,jbox,3) = lexp(j,jbox,3)+mexpnall(j)*zmul
         enddo
      enddo

      do i=1,nn1256
         jbox = n1256(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(iy,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(ix.gt.0) zmul = zmul*ys(ix,j)
            if(ix.lt.0) zmul = zmul*dconjg(ys(-ix,j))
            lexp(j,jbox,3) = lexp(j,jbox,3)+mexpn1256(j)*zmul
         enddo
      enddo

      do i=1,nn12
         jbox = n12(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(iy,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(ix.gt.0) zmul = zmul*ys(ix,j)
            if(ix.lt.0) zmul = zmul*dconjg(ys(-ix,j))
            lexp(j,jbox,3) = lexp(j,jbox,3)+mexpn12(j)*zmul
         enddo
      enddo

      do i=1,nn56
         jbox = n56(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(iy,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(ix.gt.0) zmul = zmul*ys(ix,j)
            if(ix.lt.0) zmul = zmul*dconjg(ys(-ix,j))
            lexp(j,jbox,3) = lexp(j,jbox,3)+mexpn56(j)*zmul
         enddo
      enddo

      return
      end
c------------------------------------------------------------------      
      subroutine processso2(ibox,rscale,lexp,centers,nboxes,nsall,sall,
     1           ns3478,s3478,ns34,s34,ns78,s78,mexpsall,mexps3478,
     2           mexps34,mexps78,xs,ys,zs,nexptotp)
c     shift plane waves for uall and u1234

      integer ibox
      real *8 rscale,zs(5,nexptotp),centers(3,*)
      integer nsall,ns3478,ns34,ns78,sall(*),s3478(*),s34(*),s78(*)
      complex *16 lexp(nexptotp,nboxes,6)
      complex *16 mexpsall(nexptotp),mexps3478(nexptotp)
      complex *16 mexps34(nexptotp),mexps78(nexptotp)
      complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      complex *16 zmul

c     temp variables
      integer i,jbox,j,ix,iy,iz
      real *8 ctmp(3)
      
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
      do i=1,nsall
         jbox = sall(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(-iy,j)
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(ix.lt.0) zmul = zmul*ys(-ix,j)
            if(ix.gt.0) zmul = zmul*dconjg(ys(ix,j))
            lexp(j,jbox,4) = lexp(j,jbox,4)+mexpsall(j)*zmul
         enddo
      enddo

      do i=1,ns3478
         jbox = s3478(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(-iy,j)
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(ix.lt.0) zmul = zmul*ys(-ix,j)
            if(ix.gt.0) zmul = zmul*dconjg(ys(ix,j))
            lexp(j,jbox,4) = lexp(j,jbox,4)+mexps3478(j)*zmul
         enddo
      enddo

      do i=1,ns34
         jbox = s34(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(-iy,j)
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(ix.lt.0) zmul = zmul*ys(-ix,j)
            if(ix.gt.0) zmul = zmul*dconjg(ys(ix,j))
            lexp(j,jbox,4) = lexp(j,jbox,4)+mexps34(j)*zmul
         enddo
      enddo

      do i=1,ns78
         jbox = s78(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         do j=1,nexptotp
            zmul = zs(-iy,j)
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(ix.lt.0) zmul = zmul*ys(-ix,j)
            if(ix.gt.0) zmul = zmul*dconjg(ys(ix,j))
            lexp(j,jbox,4) = lexp(j,jbox,4)+mexps78(j)*zmul
         enddo
      enddo

      return
      end
c------------------------------------------------------------------      
      subroutine mkewexp2(ibox,ilev,nboxes,ichild,rscale,nterms,nmax,
     1           rmlexp,iaddr,mptemp,mptemp2,rlams,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexpef,mexpwf,
     2           mexpephys,mexpwphys,rdplus,mexpeall,mexpe1357,mexpe13,
     4           mexpe57,mexpe1,mexpe3,mexpe5,mexpe7,mexpwall,mexpw2468,
     5           mexpw24,mexpw68,mexpw2,mexpw4,mexpw6,mexpw8,xs,
     4           ys,zs,fexpe,fexpo,rlsc)
c--------------------------------------------------------------------
c      create north south expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       real *8 rscale,rmlexp(*),rlams(nlams)
       complex *16 tmp(0:nmax,-nmax:nmax),mptemp(0:nmax,-nmax:nmax)
       complex *16 mptemp2(0:nmax,-nmax:nmax)
       real *8 rdplus(0:nmax,0:nmax,-nmax:nmax)
       complex *16 mexpef(nexptot),mexpwf(nexptot)
       complex *16 mexpephys(nexptotp),mexpwphys(nexptotp)
       complex *16 mexpeall(nexptotp),mexpwall(nexptotp)
       complex *16 mexpe1357(nexptotp),mexpw2468(nexptotp)
       complex *16 mexpe13(nexptotp),mexpe57(nexptotp)
       complex *16 mexpw24(nexptotp),mexpw68(nexptotp)
       complex *16 mexpe1(nexptotp),mexpe3(nexptotp)
       complex *16 mexpe5(nexptotp),mexpe7(nexptotp)
       complex *16 mexpw2(nexptotp),mexpw4(nexptotp)
       complex *16 mexpw6(nexptotp),mexpw8(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       real *8 zs(5,nexptotp),rlsc(0:nmax,0:nmax,nlams)
       complex *16 fexpe(1),fexpo(1)

c      temp variables
       integer jbox,ctr,ii,jj,i
       complex *16 ztmp

c      add contributions due to child 1
       jbox = ichild(1,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpe1(i) = mexpephys(i)
          mexpwall(i) = mexpwphys(i)
       enddo

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)

       endif

       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpeall(i) = mexpephys(i)/zs(1,i)
          mexpw2(i) = mexpwphys(i)*zs(1,i)
       enddo
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif
       
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpe3(i) = mexpephys(i)*dconjg(ys(1,i))
          mexpe13(i) = mexpe1(i) + mexpe3(i)
          mexpwall(i) = mexpwall(i) + mexpwphys(i)*ys(1,i)
       enddo

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif       
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = zs(1,i)*ys(1,i)
          mexpephys(i) = mexpephys(i)/ztmp
          mexpeall(i) = mexpeall(i)+mexpephys(i)
          mexpw4(i) = mexpwphys(i)*ztmp
          mexpw24(i) = mexpw2(i) + mexpw4(i)
       enddo

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif

       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          mexpe5(i) = mexpephys(i)*xs(1,i)
          mexpwall(i) = mexpwall(i) + mexpwphys(i)*dconjg(xs(1,i))
       enddo

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)/zs(1,i)
          mexpeall(i) = mexpeall(i) + mexpephys(i)*ztmp
          mexpw6(i) = mexpwphys(i)/ztmp
       enddo

c      add contributions due to child 7
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*dconjg(ys(1,i))
          mexpe7(i) = mexpephys(i)*ztmp
          mexpe57(i) = mexpe5(i)+mexpe7(i)
          mexpe1357(i) = mexpe13(i)+mexpe57(i)
          mexpwall(i) = mexpwall(i) + mexpwphys(i)*dconjg(ztmp)
       enddo

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then
c         rescale multipole expansion to be in insync with
c         scaling for pw expansion
          ctr = 0
          do ii= -nterms,nterms
             do jj = 0,nterms
                tmp(jj,ii) = dcmplx(rmlexp(iaddr(1,jbox)+ctr),
     1                       rmlexp(iaddr(1,jbox)+ctr+1))
                tmp(jj,ii) = tmp(jj,ii)/rscale**(2*jj+1)
                ctr = ctr + 2
             enddo
          enddo
c         rotate multipole expansion
          call rotztox(nterms,tmp,mptemp,rdplus)

c         convert multipole to plane waves
          call mpoletoexp2(mptemp,nterms,nlams,nfourier,nexptot,
     1         mexpef,mexpwf,rlsc)

          call ftophys(mexpef,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpephys,fexpe,fexpo)

          call ftophys(mexpwf,nlams,rlams,nfourier,nphysical,
     1         nthmax,mexpwphys,fexpe,fexpo)
       endif
       if(jbox.le.0) then
          do i=1,nexptotp
             mexpephys(i) = 0.0d0
             mexpwphys(i) = 0.0d0
          enddo
       endif

       do i=1,nexptotp
          ztmp = xs(1,i)*dconjg(ys(1,i))/zs(1,i)
          mexpeall(i) = mexpeall(i)+mexpephys(i)*ztmp+mexpe1357(i)
          mexpw8(i) = mexpwphys(i)/ztmp
          mexpw68(i) = mexpw8(i) + mexpw6(i)
          mexpw2468(i) = mexpw24(i) + mexpw68(i)
          mexpwall(i) = mexpwall(i) + mexpw2468(i)
       enddo

      return
      end
c--------------------------------------------------------------------      
      subroutine processea2(ibox,rscale,lexp,centers,nboxes,neall,eall,
     1           ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,ne3,e3,ne5,e5,
     2           ne7,e7,mexpeall,mexpe1357,mexpe13,mexpe57,mexpe1,
     2           mexpe3,mexpe5,mexpe7,xs,ys,zs,nexptotp)
c     shift plane waves for uall and u1234

      integer ibox
      real *8 rscale,zs(5,nexptotp),centers(3,*)
      integer neall,ne1357,ne13,ne57,eall(*),e1357(*),e13(*),e57(*)
      integer ne1,ne3,ne5,ne7,e1(*),e3(*),e5(*),e7(*)
      complex *16 lexp(nexptotp,nboxes,6)
      complex *16 mexpeall(nexptotp),mexpe1357(nexptotp)
      complex *16 mexpe13(nexptotp),mexpe57(nexptotp),mexpe7(nexptotp)
      complex *16 mexpe1(nexptotp),mexpe3(nexptotp),mexpe5(nexptotp)
      complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      complex *16 zmul

c     temp variables
      integer i,jbox,j,ix,iy,iz
      real *8 ctmp(3)
      
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
      do i=1,neall
         jbox = eall(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         iz = -iz
         
         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpeall(j)*zmul
         enddo
      enddo

      do i=1,ne1357
         jbox = e1357(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         iz = -iz
         
         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe1357(j)*zmul
         enddo
      enddo

      do i=1,ne13
         jbox = e13(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         iz = -iz
         
         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe13(j)*zmul
         enddo
      enddo

      do i=1,ne57
         jbox = e57(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz 

         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe57(j)*zmul
         enddo
      enddo

      do i=1,ne1
         jbox = e1(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe1(j)*zmul
         enddo
      enddo

      do i=1,ne3
         jbox = e3(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe3(j)*zmul
         enddo
      enddo

      do i=1,ne5
         jbox = e5(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe5(j)*zmul
         enddo
      enddo

      do i=1,ne7
         jbox = e7(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(ix,j)
            if(iz.gt.0) zmul = zmul*xs(iz,j)
            if(iz.lt.0) zmul = zmul*dconjg(xs(-iz,j))
            if(iy.gt.0) zmul = zmul*ys(iy,j)
            if(iy.lt.0) zmul = zmul*dconjg(ys(-iy,j))
            lexp(j,jbox,5) = lexp(j,jbox,5)+mexpe7(j)*zmul
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------      
      subroutine processwe2(ibox,rscale,lexp,centers,nboxes,nwall,wall,
     1           nw2468,w2468,nw24,w24,nw68,w68,nw2,w2,nw4,w4,nw6,w6,
     2           nw8,w8,mexpwall,mexpw2468,mexpw24,mexpw68,mexpw2,
     2           mexpw4,mexpw6,mexpw8,xs,ys,zs,nexptotp)
c     shift plane waves for uall and u1234

      integer ibox
      real *8 rscale,zs(5,nexptotp),centers(3,*)
      integer nwall,nw2468,nw24,nw68,wall(*),w2468(*),w24(*),w68(*)
      integer nw2,nw4,nw6,nw8,w2(*),w4(*),w6(*),w8(*)
      complex *16 lexp(nexptotp,nboxes,6)
      complex *16 mexpwall(nexptotp),mexpw2468(nexptotp)
      complex *16 mexpw24(nexptotp),mexpw68(nexptotp),mexpw8(nexptotp)
      complex *16 mexpw2(nexptotp),mexpw4(nexptotp),mexpw6(nexptotp)
      complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
      complex *16 zmul

c     temp variables
      integer i,jbox,j,ix,iy,iz
      real *8 ctmp(3)
      
      ctmp(1) = centers(1,ibox) - rscale/2.0d0
      ctmp(2) = centers(2,ibox) - rscale/2.0d0
      ctmp(3) = centers(3,ibox) - rscale/2.0d0
      do i=1,nwall
         jbox = wall(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
        
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpwall(j)*zmul
         enddo
      enddo

      do i=1,nw2468
         jbox = w2468(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw2468(j)*zmul
         enddo
      enddo

      do i=1,nw24
         jbox = w24(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw24(j)*zmul
         enddo
      enddo

      do i=1,nw68
         jbox = w68(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw68(j)*zmul
         enddo
      enddo

      do i=1,nw2
         jbox = w2(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw2(j)*zmul
         enddo
      enddo

      do i=1,nw4
         jbox = w4(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw4(j)*zmul
         enddo
      enddo

      do i=1,nw6
         jbox = w6(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw6(j)*zmul
         enddo
      enddo

      do i=1,nw8
         jbox = w8(i)
         ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
         iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
         iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
         iz = -iz

         do j=1,nexptotp
            zmul = zs(-ix,j)
            if(iz.gt.0) zmul = zmul*dconjg(xs(iz,j))
            if(iz.lt.0) zmul = zmul*xs(-iz,j)
            if(iy.gt.0) zmul = zmul*dconjg(ys(iy,j))
            if(iy.lt.0) zmul = zmul*ys(-iy,j)
            lexp(j,jbox,6) = lexp(j,jbox,6)+mexpw8(j)*zmul
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------      

      subroutine processudexp(ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nuall,uall,
     3           nu1234,u1234,ndall,dall,nd5678,d5678,mexpup,mexpdown,
     4           mexpupphys,mexpdownphys,mexpuall,mexpu5678,mexpdall,
     5           mexpd1234,xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      process up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nuall,ndall,nu1234,nd5678
       integer uall(*),dall(*),u1234(*),d5678(*)
       real *8 rscale
       real *8 rlams(*),whts(*)
       complex *16 tloc(0:nterms,-nterms:nterms)
       complex *16 mexp(nexptotp,nboxes,6)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpuall(nexptotp),mexpdall(nexptotp)
       complex *16 mexpd1234(nexptotp),mexpu5678(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       real *8 zs(5,nexptotp)
       real *8 rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       complex *16 ztmp,zmul
     
       real *8 ctmp(3)


       do i=1,nexptotp

       mexpuall(i) = 0
       mexpdall(i) = 0
       mexpu5678(i) = 0
       mexpd1234(i) = 0

       enddo
      
   
       ctmp(1) = centers(1,ibox) - rscale/2.0d0
       ctmp(2) = centers(2,ibox) - rscale/2.0d0
       ctmp(3) = centers(3,ibox) - rscale/2.0d0
       
       do i=1,nuall
          jbox = uall(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
             mexpdall(j) = mexpdall(j) + mexp(j,jbox,2)*zmul
          enddo
       enddo

       do i=1,nu1234
          jbox = u1234(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

          do j=1,nexptotp
             zmul = zs(iz,j)*xs(ix,j)*ys(iy,j)
             mexpd1234(j) = mexpd1234(j) + mexp(j,jbox,2)*zmul
          enddo
       enddo


       do i=1,ndall
          jbox = dall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

          do j=1,nexptotp
             zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
             mexpuall(j) = mexpuall(j) + mexp(j,jbox,1)*zmul
          enddo
       enddo




       do i=1,nd5678
          jbox = d5678(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

          do j=1,nexptotp
             zmul = zs(-iz,j)*xs(-ix,j)*ys(-iy,j)
             mexpu5678(j) = mexpu5678(j) + mexp(j,jbox,1)*zmul
          enddo
       enddo



c
cc       move contributions to the children
c


c      add contributions due to child 1

       jbox = ichild(1,ibox)

       

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)
          mexpdownphys(i) = mexpdall(i) + mexpd1234(i)
       enddo


       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)*xs(1,i)
          mexpdownphys(i) = (mexpdall(i) + mexpd1234(i))*xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)*ys(1,i)
          mexpdownphys(i) = (mexpdall(i) + mexpd1234(i))*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpuall(i)*ys(1,i)*xs(1,i)
          mexpdownphys(i) = (mexpdall(i) + mexpd1234(i))*ys(-1,i)*
     1                        xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)
          mexpdownphys(i) = mexpdall(i)/zs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)*xs(1,i)
          mexpdownphys(i) = mexpdall(i)/zs(1,i)*xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo


        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)*ys(1,i)
          mexpdownphys(i) = mexpdall(i)/zs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpuall(i)+mexpu5678(i))*zs(1,i)*ys(1,i)*
     1                       xs(1,i) 
          mexpdownphys(i) = mexpdall(i)/zs(1,i)*ys(-1,i)*xs(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc(jj,ii) = tloc(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processnsexp(ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,nnall,nall,
     3           nn1256,n1256,nn12,n12,nn56,n56,
     4           nsall,sall,ns3478,s3478,ns34,s34,ns78,s78,mexpup,
     5           mexpdown,mexpupphys,mexpdownphys,
     6           mexpnall,mexpn3478,mexpn34,mexpn78,mexpsall,
     7           mexps1256,mexps12,mexps56,rdplus,
     8           xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer nnall,nsall,nn1256,ns3478,nn12,nn56,ns34,ns78
       integer nall(*),sall(*),n1256(*),s3478(*)
       integer n12(*),n56(*),s34(*),s78(*)
       real *8 rscale
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(0:nterms,-nterms:nterms)
       complex *16 tloc2(0:nterms,-nterms:nterms)
       complex *16 mexp(nexptotp,nboxes,6)
       real *8 rdplus(0:nterms,0:nterms,-nterms:nterms)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpnall(nexptotp),mexpsall(nexptotp)
       complex *16 mexps1256(nexptotp),mexpn3478(nexptotp)
       complex *16 mexps12(nexptotp),mexps56(nexptotp)
       complex *16 mexpn34(nexptotp),mexpn78(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       real *8 zs(5,nexptotp)
       real *8 rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j
       complex *16 ztmp,zmul
     
       real *8 ctmp(3)


       do i=1,nexptotp

       mexpnall(i) = 0
       mexpsall(i) = 0
       mexpn3478(i) = 0
       mexpn34(i) = 0
       mexpn78(i) = 0
       mexps1256(i) = 0
       mexps12(i) = 0
       mexps56(i) = 0

       enddo
      
   
       ctmp(1) = centers(1,ibox) - rscale/2.0d0
       ctmp(2) = centers(2,ibox) - rscale/2.0d0
       ctmp(3) = centers(3,ibox) - rscale/2.0d0
       
       do i=1,nnall
          jbox = nall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexpsall(j) = mexpsall(j) + mexp(j,jbox,4)*zmul
          enddo

       enddo

       do i=1,nn1256
          jbox = n1256(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexps1256(j) = mexps1256(j) + mexp(j,jbox,4)*zmul
          enddo
       enddo

       do i=1,nn12
          jbox = n12(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexps12(j) = mexps12(j) + mexp(j,jbox,4)*zmul
          enddo
       enddo


       do i=1,nn56
          jbox = n56(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(iy,j)*xs(iz,j)*ys(ix,j)
             mexps56(j) = mexps56(j) + mexp(j,jbox,4)*zmul
          enddo
       enddo


       do i=1,nsall
          jbox = sall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpnall(j) = mexpnall(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

       do i=1,ns3478
          jbox = s3478(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpn3478(j) = mexpn3478(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

       do i=1,ns34
          jbox = s34(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpn34(j) = mexpn34(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

       do i=1,ns78
          jbox = s78(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-iy,j)*xs(-iz,j)*ys(-ix,j)
             mexpn78(j) = mexpn78(j) + mexp(j,jbox,3)*zmul
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)
          mexpdownphys(i) = mexpsall(i) + mexps1256(i) + mexps12(i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)*ys(1,i)
          mexpdownphys(i) = (mexpsall(i) + mexps1256(i) + mexps12(i))*
     1                        ys(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn34(i)+mexpn3478(i))*
     1         zs(1,i)
          mexpdownphys(i) = mexpsall(i)/zs(1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn34(i)+mexpn3478(i))*
     1         ys(1,i)*zs(1,i)
          mexpdownphys(i) = mexpsall(i)*ys(-1,i)/zs(1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)*xs(1,i)
          mexpdownphys(i) = (mexpsall(i) + mexps1256(i) + mexps56(i))*
     1                        xs(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpnall(i)*ys(1,i)*xs(1,i)
          mexpdownphys(i) = (mexpsall(i) + mexps1256(i) + mexps56(i))*
     1                        ys(-1,i)*xs(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn78(i)+mexpn3478(i))*
     1         xs(1,i)*zs(1,i)
          mexpdownphys(i) = mexpsall(i)*xs(-1,i)/zs(1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpnall(i)+mexpn78(i)+mexpn3478(i))*
     1         ys(1,i)*zs(1,i)*xs(1,i)
          mexpdownphys(i) = mexpsall(i)*ys(-1,i)/zs(1,i)*xs(-1,i)      
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotytoz(nterms,tloc,tloc2,rdplus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      

      subroutine processewexp(ibox,ilev,nboxes,centers,ichild,
     1           rscale,nterms,iaddr,rmlexp,rlams,whts,nlams,nfourier,
     2           nphysical,nthmax,nexptot,nexptotp,mexp,neall,eall,
     3           ne1357,e1357,ne13,e13,ne57,e57,ne1,e1,ne3,e3,ne5,e5,
     4           ne7,e7,nwall,wall,nw2468,w2468,nw24,w24,nw68,w68,
     5           nw2,w2,nw4,w4,nw6,w6,nw8,w8,
     6           mexpup,mexpdown,mexpupphys,mexpdownphys,
     7           mexpeall,mexpe2468,mexpe24,mexpe68,mexpe2,mexpe4,
     8           mexpe6,mexpe8,mexpwall,mexpw1357,mexpw13,mexpw57,
     9           mexpw1,mexpw3,mexpw5,mexpw7,rdminus,
     9           xs,ys,zs,fexpback,rlsc,rscpow)
c--------------------------------------------------------------------
c      create up down expansions for box ibox
c-------------------------------------------------------------------
       implicit none
       integer ibox,ilev,nboxes,nterms,nlams,nthmax
       integer nphysical(nlams),nfourier(nlams)
       integer iaddr(2,nboxes),ichild(8,nboxes)
       integer nexptot,nexptotp,nmax
       integer neall,nwall,ne1357,nw2468,ne13,ne57,nw24,nw68
       integer ne1,ne3,ne5,ne7,nw2,nw4,nw6,nw8
       integer eall(*),wall(*),e1357(*),w2468(*)
       integer e13(*),e57(*),w24(*),w68(*)
       integer e1(*),e3(*),e5(*),e7(*),w2(*),w4(*),w6(*),w8(*)
       real *8 rscale
       complex *16 zk2
       complex *16 rlams(*),whts(*)
       complex *16 tloc(0:nterms,-nterms:nterms)
       complex *16 tloc2(0:nterms,-nterms:nterms)
       complex *16 mexp(nexptotp,nboxes,6)
       real *8 rdminus(0:nterms,0:nterms,-nterms:nterms)
       real *8 rmlexp(*),centers(3,*)
       complex *16 mexpup(nexptot),mexpdown(nexptot)
       complex *16 mexpupphys(nexptotp),mexpdownphys(nexptotp)
       complex *16 mexpeall(nexptotp),mexpwall(nexptotp)
       complex *16 mexpw1357(nexptotp),mexpe2468(nexptotp)
       complex *16 mexpw13(nexptotp),mexpw57(nexptotp)
       complex *16 mexpe24(nexptotp),mexpe68(nexptotp)
       complex *16 mexpw1(nexptotp),mexpw3(nexptotp)
       complex *16 mexpw5(nexptotp),mexpw7(nexptotp)
       complex *16 mexpe2(nexptotp),mexpe4(nexptotp)
       complex *16 mexpe6(nexptotp),mexpe8(nexptotp)
       complex *16 xs(-5:5,nexptotp),ys(-5:5,nexptotp)
       real *8 zs(5,nexptotp)
       real *8 rlsc(0:nterms,0:nterms,nlams),rscpow(0:nterms)
       complex *16 fexpback(*)

c      temp variables
       integer jbox,ctr,ii,jj,i,ix,iy,iz,j,l
       complex *16 ztmp,zmul
     
       real *8 ctmp(3)


       do i=1,nexptotp

       mexpeall(i) = 0
       mexpwall(i) = 0
       mexpe2468(i) = 0
       mexpe24(i) = 0
       mexpe68(i) = 0
       mexpe2(i) = 0
       mexpe4(i) = 0
       mexpe6(i) = 0
       mexpe8(i) = 0
       mexpw1357(i) = 0
       mexpw13(i) = 0
       mexpw57(i) = 0
       mexpw1(i) = 0
       mexpw3(i) = 0
       mexpw5(i) = 0
       mexpw7(i) = 0

       enddo
      
   
       ctmp(1) = centers(1,ibox) - rscale/2.0d0
       ctmp(2) = centers(2,ibox) - rscale/2.0d0
       ctmp(3) = centers(3,ibox) - rscale/2.0d0
       
       do i=1,neall
          jbox = eall(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpwall(j) = mexpwall(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne1357
          jbox = e1357(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw1357(j) = mexpw1357(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne13
          jbox = e13(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw13(j) = mexpw13(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo


       do i=1,ne57
          jbox = e57(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw57(j) = mexpw57(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne1
          jbox = e1(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw1(j) = mexpw1(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo


       do i=1,ne3
          jbox = e3(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw3(j) = mexpw3(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,ne5
          jbox = e5(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw5(j) = mexpw5(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo


       do i=1,ne7
          jbox = e7(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(ix,j)*xs(-iz,j)*ys(iy,j)
             mexpw7(j) = mexpw7(j) + mexp(j,jbox,6)*zmul
          enddo
       enddo

       do i=1,nwall
          jbox = wall(i)

          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale

         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpeall(j) = mexpeall(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw2468
          jbox = w2468(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe2468(j) = mexpe2468(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw24
          jbox = w24(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe24(j) = mexpe24(j) + mexp(j,jbox,5)*zmul
          enddo

       enddo

       

       do i=1,nw68
          jbox = w68(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe68(j) = mexpe68(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw2
          jbox = w2(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe2(j) = mexpe2(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo


       do i=1,nw4
          jbox = w4(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe4(j) = mexpe4(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw6
          jbox = w6(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe6(j) = mexpe6(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

       do i=1,nw8
          jbox = w8(i)
          ix = 1.05d0*(centers(1,jbox)-ctmp(1))/rscale
          iy = 1.05d0*(centers(2,jbox)-ctmp(2))/rscale
          iz = 1.05d0*(centers(3,jbox)-ctmp(3))/rscale
         
          do j=1,nexptotp
             zmul = zs(-ix,j)*xs(iz,j)*ys(-iy,j)
             mexpe8(j) = mexpe8(j) + mexp(j,jbox,5)*zmul
          enddo
       enddo

c
cc       move contributions to the children
c


c      add contributions due to child 1
       jbox = ichild(1,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)
          mexpdownphys(i) = mexpwall(i)+mexpw1357(i)+mexpw13(i)+
     1                        mexpw1(i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 2
       jbox = ichild(2,ibox)

       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe24(i)+
     1                       mexpe2(i))*zs(1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)


        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)


        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif
  
c      add contributions due to child 3
       jbox = ichild(3,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)*ys(1,i)
          mexpdownphys(i) = (mexpwall(i)+mexpw1357(i)+mexpw13(i)+
     1                        mexpw3(i))*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 4
       jbox = ichild(4,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe24(i)+
     1                       mexpe4(i))*zs(1,i)*ys(1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 5
       jbox = ichild(5,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)*xs(-1,i)
          mexpdownphys(i) = (mexpwall(i)+mexpw1357(i)+mexpw57(i)+
     1                        mexpw5(i))*xs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 6
       jbox = ichild(6,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe68(i)+
     1                       mexpe6(i))*zs(1,i)*xs(-1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)*xs(1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)

        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 7
 
       jbox = ichild(7,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = mexpeall(i)*xs(-1,i)*ys(1,i)
          mexpdownphys(i) = (mexpwall(i)+mexpw1357(i)+mexpw57(i)+
     1                        mexpw7(i))*xs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

c      add contributions due to child 8
       jbox = ichild(8,ibox)
       if(jbox.gt.0) then

       do i=1,nexptotp
          mexpupphys(i)  = (mexpeall(i)+mexpe2468(i)+mexpe68(i)+
     1                       mexpe8(i))*zs(1,i)*xs(-1,i)*ys(1,i)      
          mexpdownphys(i) = mexpwall(i)/zs(1,i)*xs(1,i)*ys(-1,i)
       enddo

       call phystof(mexpup,nlams,nfourier,nphysical,
     1               mexpupphys,fexpback)
 
       call phystof(mexpdown,nlams,nfourier,nphysical,
     1              mexpdownphys,fexpback)

        call exptolocal2(tloc,nterms,rlams,whts,
     1         nlams,nfourier,nthmax,nexptot,mexpup,mexpdown,
     2         rscale,rlsc)

        call rotztox(nterms,tloc,tloc2,rdminus)


        do ii=-nterms,nterms
           do jj=0,nterms
              tloc2(jj,ii) = tloc2(jj,ii)/rscpow(jj)
           enddo
        enddo

        call l3dadd(tloc2,rmlexp(iaddr(2,jbox)),nterms)

       endif

      return
      end
c--------------------------------------------------------------------      
