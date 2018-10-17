c     Plane wave routines for 3D FMM

      subroutine mpoletoexp2(mpole,nterms,nlambs,numtets,nexptot,
     1                mexpupf,mexpdownf,rlsc)

c     This subroutine converts a multipole expansion into the
c     corresponding exponential moment function mexp for
c     both the +z direction and the -z direction
c
c     Note: this subroutine is the same as mpoletoexp
c     in mtxbothnew.f but just has a different data structure
c     for mpole
c
c     U(x,y,z) = \sum_{n=0}^{nterms} \sum_{m=-n,n} mpole(n,m)
c                Y_n^m (\cos(\theta)) e^{i m \phi}/r^{n+1}
c  
c              = (1/2\pi) \int_{0}^{\infty} e^{-\lambda z}
c                \int_{0}^{2\pi} e^{i \lambda (x \cos(\alpha) +
c                y \sin(\alpha))} mexpup(\lambda,\alpha)
c                d\alpha d\lambda
c 
c     for the +z direction and
c
c              = (1/2\pi) \int_{0}^{\infty} e^{\lambda z}
c                \int_{0}^{2\pi} e^{-i \lambda (x \cos(\alpha) +
c                y \sin(\alpha))} mexpdown(\lambda,\alpha)
c                d\alpha d\lambda
c 
c     for the -z direction
c
c     NOTE: The expression for -z corresponds to the mapping
c     (x,y,z) -> (-x,-y,-z), ie reflection through
c     the origin.
c
c     NOTE 2: The multipole expansion is assumed to have been
c     rescaled so that the box containing the sources has unit
c     dimension
c
c     NOTE 3: Since we store the exponential moment function in
c     Fourier domain (w.r.t the \alpha variable) we compute
c 
c     M_\lambda(m) = (i)**m \sum_{n=m}^{N} c(n,m) mpole(n,m)
c     lambda^n
c 
c     for m >=0 only, where c(n,m) = 1/sqrt((n+m)!(n-m)!)
c
c     For possible future reference, it should be noted that it
c     is NOT true that M_\lambda(-m) = dconjg(M_\lambda(m))
c
c     Inspection of the integral formula for Y_n^{-m} shows
c     that M_\lambda(-m) = dconjg(M_\lambda) * (-1)**m
c
c     INPUT arguments
c     mpole       in: complex *16 (0:nterms, -nterms:nterms)
c                 The multipole expansion 
c  
c     nterms:     in: integer
c                 Order of the multipole expansion
c
c     nlambs      in: integer
c                 number of discretization points in the \lambda
c                 integral
c
c     numtets     in: integer(nlambs)
c                 number of fourier modes needed in expansion
c                 of \alpha variable for each \lambda variable
c
c     nexptot     in: integer
c                 nexptot = \sum_{j} numtets(j)
c
c     rlsc        in: real *8(0:nterms, 0:nterms,nlambs)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT 
c     mexpupf     out: complex *16 (nexptot)
c                 Fourier coefficients of the function
c                 mexpup(\lambda,\alpha) for successive
c                 discrete lambda values. They are ordered as
c                 follows
c
c                 mexpupf(1,...., numtets(1)) = fourier modes
c                             for \lambda_1
c
c                 mexpupf(numtets(1)+1,...., numters(2) = fourier
c                 modes for \lambda_2
c
c                 ETC
c
c     mexpdownf   out: complex *16 (nexptot)
c                 Fourier coefficients of the function 
c                 mexpdown(\lambda,\alpha) for successive
c                 discrete \lambda values
c---------------------------------------------------------------

      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mexpupf(*)
      complex *16 mexpdownf(*)
      complex *16 zeyep,ztmp1,ztmp2
      real *8 rlsc(0:nterms,0:nterms,nlambs)

c     Temp variables
      real *8 sgn
      integer ntot,ncurrent,nl,mth,nm

      ntot = 0
      do nl=1,nlambs
         sgn = -1.0d0
         zeyep = 1.0d0
         do mth = 0,numtets(nl)-1
            ncurrent = ntot + mth + 1
            ztmp1 = 0.0d0
            ztmp2 = 0.0d0
            sgn = -sgn
            do nm = mth,nterms,2
               ztmp1 = ztmp1 + rlsc(nm,mth,nl)*mpole(nm,mth)
            enddo

            do nm=mth+1,nterms,2
               ztmp2 = ztmp2 + rlsc(nm,mth,nl)*mpole(nm,mth)
            enddo
            mexpupf(ncurrent) = (ztmp1+ztmp2)*zeyep
            mexpdownf(ncurrent) = sgn*(ztmp1-ztmp2)*zeyep
            zeyep = zeyep*dcmplx(0.0d0,1.0d0)
         enddo
         ntot = ntot + numtets(nl)
      enddo

      return
      end

c -----------------------------------------------------------------
      subroutine exptolocalup2(local,nterms,rlambs,whts,nlambs,numtets,
     1                         nthmax,nexptot,lexp1f,scale,rlsc)
c-----------------------------------------------------------------
c     This sburoutine converts the Fourier representation of
c     an exponential moment function into a local
c     multipole expansion (with respect to the same box center)
c
c     u(x,y,z) = \int_{0}^{\infty} e^{-\lambda z}
c                \int_{0}^{2\pi} e^{i\lambda(x\cos(\alpha) +
c                y\sin(\alpha))} lexp1 (\lambda,\alpha) d\lambda
c                d\alpha
c
c              = \sum_{0}^{nterms} \sum_{m=-n,n} local(n,m) Y_n^m
c                 (\cos(\theta)) e^{i m \phi} r^{n}
c
c     INPUT arguments
c     nterms           in: integer
c                      Order of local expansion
c
c     rlambs           in: real *8(nlambs)
c                      discretization points in the \lambda integral
c
c     whts             in: real *8(nlambs)
c                      quadrature weights in \lambda integral
c
c     nlambs           in: integer
c                      number of discretization points in \lambda
c                      integral
c
c     numtets          in: integer(nlambs)
c                      number of fourier modes in expansion of
c                      \alpha variable for \lambda_j
c
c     nthmax           in: integer
c                      max_j numtets(j)
c
c     nexptot          in: integer
c                      sum_j numtets(j)
c                      
c
c     lexp1f(nexptot)  complex *16(nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values.
c                      They are ordered as follows:
c
c                      lexp1f(1,...,numtets(1)) = Fourier modes
c                      for \lambda_1
c                      lexp1f(numtets(1)+1,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c     scale            in: real *8
c                      scaling parameter for local expansion
c
c     rlsc        in: real *8(nlambs, 0:nterms, 0:nterms)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT
c     local(0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot,nthmax
      integer ncurrent,ntot,nl
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 lexp1f(nexptot)
      complex *16 zeye(0:nterms),ima
      real *8 rlambs(nlambs), rlambpow(0:nterms), whts(nlambs)
      real *8 rlsc(0:nterms,0:nterms,nlambs)
      real *8 rmul
      real *8 scale, rscale(0:nterms)
    
c     Temporary variables
      integer i, nm, mth, j, mmax

      data ima/(0.0d0,1.0d0)/

      

      zeye(0) = 1.0d0
      do i=1,nterms
         zeye(i) = zeye(i-1)*ima
      enddo

      rscale(0) = 1
      do nm=0,nterms
         if(nm.gt.0) rscale(nm) = rscale(nm-1)*scale

         do mth = -nterms,nterms
            local(nm,mth) = 0.0d0
         enddo
      enddo

      ntot = 1
      do nl=1,nlambs
c        Add contributions to local expansion
         do nm=0,nterms,2
            mmax = numtets(nl)-1
            if(mmax.gt.nm) mmax = nm
            do mth=0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth)+rlsc(nm,mth,nl)*
     1              lexp1f(ncurrent)*whts(nl)
            enddo
         enddo
         do nm=1,nterms,2
            mmax = numtets(nl) - 1
            if(mmax.gt.nm) mmax = nm
            do mth =0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth) - rlsc(nm,mth,nl)*
     1             lexp1f(ncurrent)*whts(nl)
            enddo
         enddo
         ntot = ntot + numtets(nl)
      enddo

      do nm=0,nterms
         do mth = 0,nm
            local(nm,mth) = local(nm,mth)*zeye(mth)
            if(mth.gt.0) local(nm,-mth) = dconjg(local(nm,mth))
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------      

      subroutine exptolocaldown2(local,nterms,rlambs,whts,nlambs,
     1    numtets,nthmax,nexptot,lexp1f,scale,rlsc)
c-----------------------------------------------------------------
c     This sburoutine converts the Fourier representation of
c     an exponential moment function into a local
c     multipole expansion (with respect to the same box center)
c
c     u(x,y,z) = \int_{0}^{\infty} e^{\lambda z}
c                \int_{0}^{2\pi} e^{-i\lambda(x\cos(\alpha) +
c                y\sin(\alpha))} lexp1 (\lambda,\alpha) d\lambda
c                d\alpha
c
c              = \sum_{0}^{nterms} \sum_{m=-n,n} local(n,m) Y_n^m
c                 (\cos(\theta)) e^{i m \phi} r^{n}
c
c     INPUT arguments
c     nterms           in: integer
c                      Order of local expansion
c
c     rlambs           in: real *8(nlambs)
c                      discretization points in the \lambda integral
c
c     whts             in: real *8(nlambs)
c                      quadrature weights in \lambda integral
c
c     nlambs           in: integer
c                      number of discretization points in \lambda
c                      integral
c
c     numtets          in: integer(nlambs)
c                      number of fourier modes in expansion of
c                      \alpha variable for \lambda_j
c
c     nthmax           in: integer
c                      max_j numtets(j)
c
c     nexptot          in: integer
c                      sum_j numtets(j)
c                      
c
c     lexp1f(nexptot)  complex *16(nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values.
c                      They are ordered as follows:
c
c                      lexp1f(1,...,numtets(1)) = Fourier modes
c                      for \lambda_1
c                      lexp1f(numtets(1)+1,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c     scale            in: real *8
c                      scaling parameter for local expansion
c
c     rlsc        in: real *8(nlambs, 0:nterms, 0:nterms)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT
c     local(0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot,nthmax
      integer ncurrent,ntot,nl
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 lexp1f(nexptot)
      complex *16 zeye(0:nterms),ima
      real *8 rlambs(nlambs), rlambpow(0:nterms) ,whts(nlambs)
      real *8 rmul
      real *8 scale, rscale(0:nterms)
      real *8 rlsc(0:nterms,0:nterms,nlambs)
    
c     Temporary variables
      integer i, nm, mth, j, mmax

      data ima/(0.0d0,1.0d0)/

      

      zeye(0) = 1.0d0
      do i=1,nterms
         zeye(i) = zeye(i-1)*ima
      enddo

      rscale(0) = 1.0d0
      do nm=0,nterms

         if(nm.gt.0) rscale(nm) = rscale(nm-1)*scale
         do mth = -nterms,nterms
            local(nm,mth) = 0.0d0
         enddo
      enddo

      ntot = 1
      do nl=1,nlambs
c        Add contributions to local expansion
         do nm=0,nterms
            mmax = numtets(nl)-1
            if(mmax.gt.nm) mmax = nm
            do mth=0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth)+rlsc(nm,mth,nl)*
     1              lexp1f(ncurrent)*whts(nl)
            enddo
         enddo
         ntot = ntot + numtets(nl)
      enddo

      do nm=0,nterms
         do mth = 0,nm
            local(nm,mth) = local(nm,mth)*zeye(mth)
            if(mth.gt.0) local(nm,-mth) = dconjg(local(nm,mth))
         enddo
      enddo

      return
      end
c------------------------------------------------------------------
      subroutine exptolocal2(local,nterms,rlambs,whts,nlambs,numtets,
     1                         nthmax,nexptot,lexp1f,lexp2f,scale,rlsc)
c-----------------------------------------------------------------
c     INPUT arguments
c     nterms           in: integer
c                      Order of local expansion
c
c     rlambs           in: real *8(nlambs)
c                      discretization points in the \lambda integral
c
c     whts             in: real *8(nlambs)
c                      quadrature weights in \lambda integral
c
c     nlambs           in: integer
c                      number of discretization points in \lambda
c                      integral
c
c     numtets          in: integer(nlambs)
c                      number of fourier modes in expansion of
c                      \alpha variable for \lambda_j
c
c     nthmax           in: integer
c                      max_j numtets(j)
c
c     nexptot          in: integer
c                      sum_j numtets(j)
c                      
c
c     lexp1f(nexptot)  complex *16(nexptot)
c                      Fourier coefficients of the function 
c                      lexp1 for discrete \lambda values.
c                      They are ordered as follows:
c
c                      lexp1f(1,...,numtets(1)) = Fourier modes
c                      for \lambda_1
c                      lexp1f(numtets(1)+1,...,numtets(2) = Fourier
c                      modes for \lambda_2 etc
c
c     scale            in: real *8
c                      scaling parameter for local expansion
c
c     rlsc        in: real *8(nlambs, 0:nterms, 0:nterms)
c                 scaled discretization points in the \lambda
c                 integral
c
c     OUTPUT
c     local(0:nterms,-nterms:nterms): output local expansion of order
c                                     nterms
        
      implicit none
      integer nterms,nlambs,numtets(nlambs),nexptot,nthmax
      integer ncurrent,ntot,nl
      complex *16 local(0:nterms,-nterms:nterms)
      complex *16 lexp1f(nexptot),lexp2f(nexptot)
      complex *16 zeye(0:nterms)
      real *8 rlambs(nlambs), rlambpow(0:nterms) ,whts(nlambs)
      real *8 rmul,rlsc(0:nterms,0:nterms,nlambs)
      real *8 scale, rscale(0:nterms)
      complex *16 ima
    
c     Temporary variables
      integer i, nm, mth, j, mmax

      data ima/(0.0d0,1.0d0)/


      zeye(0) = 1.0d0
      do i=1,nterms
         zeye(i) = zeye(i-1)*ima
      enddo

      rscale(0) = 1
      do nm=0,nterms
         if(nm.gt.0) rscale(nm) = rscale(nm-1)*scale
         do mth = -nterms,nterms
            local(nm,mth) = 0.0d0
         enddo
      enddo

      ntot = 1
      do nl=1,nlambs
c        Add contributions to local expansion
         do nm=0,nterms,2
            mmax = numtets(nl)-1
            if(mmax.gt.nm) mmax = nm
            do mth=0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth)+rlsc(nm,mth,nl)*
     1          (lexp1f(ncurrent)+lexp2f(ncurrent))*whts(nl)
            enddo
         enddo
         do nm=1,nterms,2
            mmax = numtets(nl) - 1
            if(mmax.gt.nm) mmax = nm
            rmul = rlambpow(nm)
            do mth =0,mmax
               ncurrent = ntot+mth
               local(nm,mth) = local(nm,mth)-rlsc(nm,mth,nl)*
     1          (lexp1f(ncurrent)-lexp2f(ncurrent))*whts(nl)
            enddo
         enddo
         ntot = ntot + numtets(nl)
      enddo

      do nm=0,nterms
         do mth = 0,nm
            local(nm,mth) = local(nm,mth)*zeye(mth)
            if(mth.gt.0) local(nm,-mth) = dconjg(local(nm,mth))
         enddo
      enddo

      return
      end
c------------------------------------------------
