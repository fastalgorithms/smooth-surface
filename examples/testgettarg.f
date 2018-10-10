      implicit real *8 (a-h,o-z)
      real *8 sources(3,1000),targ(3,10000)
      real *8, allocatable :: targ2(:,:)

      integer, allocatable :: itree(:),iptr(:)
      integer ipointer(32)
      real *8, allocatable :: boxsize(:),treecenters(:,:)
      real *8, allocatable :: fvals(:),fcoeffs(:),xmatc(:,:)
      real *8 target(3)

      call prini(6,13)
      call prin2('Enter n*',n,0)
      read *, n 

      
      ns = 9
      nt = 9

      sources(1,1) = 0
      sources(2,1) = 0
      sources(3,1) = 0

      sources(1,2) = 1
      sources(2,2) = 0
      sources(3,2) = 0


      sources(1,3) = 1
      sources(2,3) = 1
      sources(3,3) = 0

      sources(1,4) = 0
      sources(2,4) = 1
      sources(3,4) = 0

      sources(1,5) = 0
      sources(2,5) = 0
      sources(3,5) = 1

      sources(1,6) = 1
      sources(2,6) = 0
      sources(3,6) = 1


      sources(1,7) = 1
      sources(2,7) = 1
      sources(3,7) = 1

      sources(1,8) = 0
      sources(2,8) = 1
      sources(3,8) = 1

      sources(1,9) = 1 - 0.124d0
      sources(2,9) = 1
      sources(3,9) = 1

      do i=1,9
         do j=1,3
            targ(j,i) = sources(j,i)
         enddo
      enddo

      call gettargandtreemem(ns,sources,nt,targ,lmem,nlbox,nlevels,
     1      nboxes)
      call prinf('lmem=*',lmem,1)
      call prinf('nlbox=*',nlbox,1)
      call prinf('nboxes=*',nboxes,1)
      call prinf('nlevels=*',nlevels,1)

      norder = 4
      norder3 = norder*norder*norder
      nt2 = nlbox*norder3

      allocate(targ2(3,nt2))
      allocate(itree(lmem),boxsize(0:nlevels),treecenters(3,nboxes))
      allocate(iptr(nboxes))

      call gettargandtree(ns,sources,nt,targ,itree,ipointer,
     1      treecenters,boxsize,norder,targ2,iptr)

       call prinf('iptr=*',iptr,nboxes)

 2100 format(3(2x,e11.5))

      do i=1,nt2
         write(23,2100) targ2(1,i),targ2(2,i),targ2(3,i)
      enddo

c
cc      evaluate predefined function at these nodes
c

      allocate(fvals(nt2),fcoeffs(nt2))

      do i=1,nt2
         call fun1(targ2(1,i),targ2(2,i),targ2(3,i),fvals(i))
      enddo

      allocate(xmatc(norder3,norder3))

      call getxmatc(norder,norder3,xmatc)

c
cc     compute the hebyshev expansion coefficients 
c      at the nodes
      do i=1,nlbox
        istart = (i-1)*norder3+1
        call matvec(norder3,xmatc,fvals(istart),fcoeffs(istart))
      enddo

c
cc       generate random targets in a particular leaf ndoes
c
      ibox = 6
      xtest =treecenters(1,ibox) + (hkrand(0)*-0.5d0)*boxsize(1)
      ytest =treecenters(2,ibox) + (hkrand(0)*-0.5d0)*boxsize(1)
      ztest =treecenters(3,ibox) + (hkrand(0)*-0.5d0)*boxsize(1)

      target(1) = xtest
      target(2) = ytest
      target(3) = ztest


      call findboxtarg(target,iboxtarg,ilevel,itree,ipointer,boxsize,
     1       treecenters,nboxes,nlevels)

      call prinf('iboxtarg=*',iboxtarg,1)
      call prinf('ilevel=*',ilevel,1)

c
cc      evaluate expansions to get
c
      call cheb3deval(target,boxsize(ilevel),treecenters(1,iboxtarg),
     1     norder,fcoeffs(iptr(iboxtarg)),f)

      call fun1(target(1),target(2),target(3),fex)

      call prin2('f=*',f,1)
      call prin2('fex=*',fex,1)

      err = abs(fex-f)

      call prin2('error=*',err,1)


      stop
      end
c-----------------------      
      subroutine fun1(x,y,z,f)
      implicit real *8 (a-h,o-z)
      
      f = sin(x*y*z)
cc      f = 1

      return
      end
c------------------------

       subroutine matvec(n,a,x,y)
       implicit real *8 (a-h,o-z)
       real *8 a(n,n),x(n),y(n)

       do i=1,n
          y(i) = 0 
          do j=1,n
            y(i) =y(i)  +a(i,j)*x(j)
          enddo
       enddo
       
       return
       end
      
