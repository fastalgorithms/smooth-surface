!-----------------------------------

module prefunrouts
!use Mod_Fast_Sigma

contains

subroutine initialize_feval(eps,ns,sources,skeleton_w,rn,nt,dumtarg,norder, &
    itree,ltree,nlevels,nboxes,iptr,treecenters,boxsize, &
    nt2,fcoeffs,fcoeffsx,fcoeffsy,fcoeffsz,FSS_1,adapt_flag)

  use Mod_Fast_Sigma


  !       this subroutine tree sorts a collection of sources
  !        and targets, and returns the tree  
  !        
  !
  !
  !        input:
  !        eps - precision requested
  !        ns - number of sources
  !        sources(3,*) - xyz locations of the sources
  !        rn(3,*) - normal orientations at the sources
  !        nt - number of targets
  !        dumtarg(3,*) - xyz locations of the dummy targets
  !        sigma_eval - function handle for evaluating sigma 
  !      
  !        output
  !        itree(ltree) - tree array
  !        ltree - length of tree array
  !        iptr(7) - pointer for tree array
  !        nlevels - number of levels
  !        nboxes - number of boxes
  !        treecenters(3,nboxes) - locations of centers of trees
  !        boxsize(0:nlevels) - size of boxes at various levels
  !        nt2 - number of dummy targets
  !        fcoeffs(nlbox*norder**3) - chebyshev expansion coefficients
  !                                   of the function to be evaluated
  !                                   for the revelant leaf boxes
  !        
  !        fcoeffsx(nlbox*norder**3) - chebyshev expansion coefficients
  !                                   of the x component of the gradient of
  !                                   the function to be evaluated
  !                                   for the revelant leaf boxes
  !
  !        fcoeffsy(nlbox*norder**3) - chebyshev expansion coefficients
  !                                   of the y component of the gradient of
  !                                   the function to be evaluated
  !                                   for the revelant leaf boxes
  !
  !        fcoeffsz(nlbox*norder**3) - chebyshev expansion coefficients
  !                                   of the z component of the gradient of
  !                                   the function to be evaluated
  !                                   for the revelant leaf boxes
  !
  !      use Mod_TreeLRD
  implicit double precision (a-h,o-z)

  !      calling sequence variables
  integer, intent(in) :: ns,nt
  double precision, intent(in) :: eps
  double precision, intent(in) :: sources(3,ns)
  double precision, intent(in) :: rn(3,ns)
  double precision, intent(in) :: skeleton_w(ns)
  double precision, intent(inout) :: dumtarg(3,nt)

  integer, allocatable, intent(inout) :: itree(:)
  integer, allocatable, intent(inout) :: iptr(:)
  integer, intent(inout) :: nboxes,nlevels,ltree,nt2
  integer, intent(in) :: norder

  double precision, allocatable, intent(inout) :: treecenters(:,:)
  double precision, allocatable, intent(inout) :: boxsize(:)
  double precision, allocatable, intent(inout) :: fcoeffs(:)
  double precision, allocatable, intent(inout) :: fcoeffsx(:)
  double precision, allocatable, intent(inout) :: fcoeffsy(:)
  double precision, allocatable, intent(inout) :: fcoeffsz(:)
  type ( Fast_Sigma_stuff ), pointer :: FSS_1
  integer, intent(in) :: adapt_flag
  !        temporary variables      

  double precision, allocatable :: radsrc(:)
  integer, allocatable :: ntargbox(:)
  integer, allocatable :: itstartbox(:)
  double precision, allocatable :: dumtargsort(:,:)

  double precision, allocatable :: targets(:,:)
  double precision, allocatable :: fvals(:),fvalsx(:),fvalsy(:),fvalsz(:)
  double precision fvalstmp(norder*norder*norder)
  double precision fvalsxtmp(norder*norder*norder)
  double precision fvalsytmp(norder*norder*norder)
  double precision fvalsztmp(norder*norder*norder)
  !      double precision, allocatable :: fvalstmp(:)
  !      double precision, allocatable :: fvalsxtmp(:),fvalsytmp(:),fvalsztmp(:)

  integer nlbox      

  double precision, allocatable :: sigma(:),sigma_grad(:,:),trads(:)
  double precision sigmatmp,sigma_gradtmp(3)

  !      fmm variables
  integer, allocatable :: itreefmm(:)
  double precision, allocatable :: treecentersfmm(:,:),boxsizefmm(:)
  integer ipointer(32)
  integer, allocatable :: iaddr(:,:),nterms(:)
  double precision, allocatable :: rmlexp(:)
  complex *16, allocatable :: charges(:),dipole(:)
  complex *16, allocatable :: pottarg(:),fldtarg(:,:)

  double precision xpts(norder),wts(norder)
  double precision, allocatable :: xmatc(:,:)
  double precision, allocatable :: tails(:)
  double precision, allocatable :: tailsx(:)
  double precision, allocatable :: tailsy(:)
  double precision, allocatable :: tailsz(:)
  integer, allocatable :: irefineflag(:)
  integer, allocatable :: ibcompflag(:)
  integer, allocatable :: itcompflag(:)

  integer, allocatable :: itreeout(:)
  integer, allocatable :: ntargboxout(:)
  integer, allocatable :: itstartboxout(:)
  integer iptrout(7)
  double precision, allocatable :: targetsout(:,:),fcoeffsout(:), &
      boxsizeout(:),treecentersout(:,:),fcoeffsxout(:), &
      fcoeffsyout(:), fcoeffszout(:)
  integer nnnn, nadapt_flag
  !
  !c      more temporary variables
  !
  double precision expc(3),radt,wlege(40000),pottmp,gradtmp(3)
  double precision, allocatable :: sourcesort(:,:),dipvecsort(:,:)
  complex *16, allocatable :: dipstrsort(:),chargesort(:)
  integer, allocatable :: ilevel(:)
  real ( kind = 8 ) sigmatmp_v(1), sigma_gradtmp_v_x(1)
  real ( kind = 8 ) sigma_gradtmp_v_y(1),sigma_gradtmp_v_z(1)
  double precision tradtmp


  !
  !!      sigma tree variables
  !
  !      type (TreeLRD) :: TreeLRD_1

  !      external sigma_eval

  nbmax = 0
  nlmax = 200
  mnlist1 = 0
  mnlist2 = 0
  mnlist3 = 0
  mnlist4 = 0
  mnbors = 0
  ltree = 0

  idivflag = 0
  ndiv = 45

  isep = 1

  allocate(radsrc(ns))
  do i=1,ns
    radsrc(i) = 0
  enddo

  ier = 0
  nexpc = 0
  radexp = 0

  !
  !!         allocate initial tree with sources and dummy targets
  !


  write (*,*) 'Im inside'

  call mklraptreemem(ier,sources,ns,radsrc,dumtarg,nt,expc, &
      nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,nlevelsfmm, &
      nboxesfmm,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung , &
      ltreefmm)
  write (*,*) 'Im inside 2',ltreefmm,nlevelsfmm,nboxesfmm

  allocate(itreefmm(ltreefmm))
  allocate(boxsizefmm(0:nlevelsfmm))
  allocate(treecentersfmm(3,nboxesfmm))

  call mklraptree(sources,ns,radsrc,dumtarg,nt,expc,nexpc, &
      radexp,idivflag,ndiv,isep,mhung,mnbors,mnlist1, &
      mnlist2,mnlist3,mnlist4,nlevelsfmm,nboxesfmm, &
      treecentersfmm,boxsizefmm,itreefmm,ltreefmm,ipointer)


  nboxes = nboxesfmm
  nlevels = nlevelsfmm

  write (*,*) 'Im inside 3',nboxes,nlevels

  allocate(boxsize(0:levels),treecenters(3,nboxes))


  treecenters = treecentersfmm
  boxsize = boxsizefmm

  ltree = 2*(nlevels+1)+12*nboxes

  allocate(iptr(7),itree(ltree),ntargbox(nboxesfmm))
  allocate(itstartbox(nboxesfmm),dumtargsort(3,nt))

  !
  !c      extract relevant bits of the tree needed from this point on
  !

  call extracttreesubinfo(itreefmm,ipointer,nboxes,nlevels, &
      itree,iptr,itstartbox,ntargbox)
  !
  !!        treesort dumtarg
  !

  do i=1,nt
    ii = itreefmm(ipointer(6)+i-1)
    dumtargsort(1,i) = dumtarg(1,ii)
    dumtargsort(2,i) = dumtarg(2,ii)
    dumtargsort(3,i) = dumtarg(3,ii)
  enddo

  !
  !c        generate initial set of targets
  !
  itype = 1
  call chebexps(itype,norder,xpts,utmp,vtmp,wts)

  nlbox = 0
  do i=1,nboxes
    nchild = itree(iptr(4)+i-1)
    if(nchild.eq.0.and.ntargbox(i).gt.0) nlbox = nlbox+1
    !          write (*,*) nchild,nlbox,iptr(4)+i-1
  enddo

  itarg = 0
  ictr = 1
  norder3 = norder*norder*norder

  nt2 = norder3*nlbox


  allocate(targets(3,nt2),fvals(nt2),fcoeffs(nt2))
  allocate(fvalsx(nt2),fvalsy(nt2),fvalsz(nt2))
  allocate(fcoeffsx(nt2),fcoeffsy(nt2),fcoeffsz(nt2))
  !        allocate(fvalstmp(norder3),fvalsxtmp(norder3))
  !        allocate(fvalsytmp(norder3),fvalsztmp(norder3))

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
            targets(1,itarg) = treecentersfmm(1,ibox)+ &
                xpts(kk)*boxsizefmm(ilev)/2
            targets(2,itarg) = treecentersfmm(2,ibox)+ &
                xpts(jj)*boxsizefmm(ilev)/2
            targets(3,itarg) = treecentersfmm(3,ibox)+ &
                xpts(ii)*boxsizefmm(ilev)/2
          enddo
        enddo
      enddo
    endif
  enddo

  !
  !c       evaluate sigma, trads, stdev initially at the targets
  !
  !
  allocate(sigma(nt2),sigma_grad(3,nt2),trads(nt2))

!!! ESTO LO HE CAMBIADO YO...
  !        npts = 1
  !        do i=1,nt2
  !            call sigma_eval1(targets(1,i),sigma(i),sigma_grad(:,i))
  !
  !              trads(i) = 12*sigma(i)
  !        enddo

  !write (*,*) 'nt2', nt2

  call function_eval_sigma(FSS_1,targets,int(nt2),sigma,sigma_grad(1,:),&
      &sigma_grad(2,:),sigma_grad(3,:),int(adapt_flag))

  trads(1:nt2) = 12*sigma(1:nt2)


  !
  !        recompute tree with sources dummy targets and initial 
  !        set of targets
  !



  call mklraptreemem(ier,sources,ns,radsrc,dumtarg,nt,targets, &
      nt2,trads,idivflag,ndiv,isep,nlmax,nbmax,nlevelsfmm, &
      nboxesfmm,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,mhung , &
      ltreefmm)

  deallocate(itreefmm)
  allocate(itreefmm(ltreefmm))



  call mklraptree(sources,ns,radsrc,dumtarg,nt,targets,nt2, &
      trads,idivflag,ndiv,isep,mhung,mnbors,mnlist1, &
      mnlist2,mnlist3,mnlist4,nlevelsfmm,nboxesfmm, &
      treecentersfmm,boxsizefmm,itreefmm,ltreefmm,ipointer)

  !
  !!      set fmm flags
  !
  iprec = 1
  if(eps.le.0.5d-3.and.eps.gt.0.5d-6) iprec = 2
  if(eps.le.0.5d-6.and.eps.gt.0.5d-9) iprec = 3
  if(eps.le.0.5d-9) iprec = 4
  ier = 0


  allocate(iaddr(2,nboxes))
  allocate(nterms(0:nlevels))

  if( iprec .eq. -2 ) epsfmm=.5d-0 
  if( iprec .eq. -1 ) epsfmm=.5d-1
  if( iprec .eq. 0 ) epsfmm=.5d-2
  if( iprec .eq. 1 ) epsfmm=.5d-3
  if( iprec .eq. 2 ) epsfmm=.5d-6
  if( iprec .eq. 3 ) epsfmm=.5d-9
  if( iprec .eq. 4 ) epsfmm=.5d-12
  if( iprec .eq. 5 ) epsfmm=.5d-15
  if( iprec .eq. 6 ) epsfmm=0


  do i=0,nlevels
    if(isep.eq.1) call l3dterms(epsfmm,nterms(i),ier)
    if(isep.eq.2) call l3dterms_far(epsfmm,nterms(i),ier)
  enddo

  call l3dmpalloc_newtree(itree(1),iaddr,nlevels,lmptot,nterms)
  call prinf(' lmptot is *',lmptot,1)


  allocate(rmlexp(lmptot),stat=ier)
  if(ier.ne.0) then
    call prinf('Cannot allocate mpole expansion workspace, &
        lmptot is *', lmptot,1)
    ier = 16
    return
  endif

  ifcharge = 0
  ifdipole = 1
  ifpottarg = 1
  iffldtarg = 1
  ier = 0

  allocate(charges(ns),dipole(ns),pottarg(nt2),fldtarg(3,nt2))
  allocate(sourcesort(3,ns))
  allocate(chargesort(ns),dipstrsort(ns),dipvecsort(3,ns))

  do i=1,ns
    dipole(i) = skeleton_w(i)
    charges(i) = 0 
  enddo

  !
  !!      call the fmm

  call tfmm3dlr(ier,iprec,ns,sources,ifcharge,charges,ifdipole, &
      dipole,rn,nt2,targets,trads,sigma,sigma_grad,itreefmm, &
      ltreefmm,ipointer,isep,ndiv,nlevels,nboxes,boxsizefmm, &
      mnbors,mnlist1,mnlist2,mnlist3,mnlist4,treecentersfmm, &
      iaddr,nterms,lmptot,ifpottarg,pottarg,iffldtarg,fldtarg, &
      rmlexp)

  !
  !!      store the sorted arrays
  !
  call l3dreorder_newtree(ns,sources,ifcharge,charges, &
      itreefmm(ipointer(5)),ifdipole,dipole,rn,sourcesort, &
      chargesort,dipstrsort,dipvecsort)


  !
  !!       extract the ilevel array
  !
  allocate(ilevel(nboxesfmm))
  do ilev=0,nlevelsfmm
    do ibox = itreefmm(2*ilev+1),itreefmm(2*ilev+2)
      ilevel(ibox) = ilev
    enddo
  enddo

  !
  !!     extract fvals and the gradient values
  !
  do i=1,nt2
    fvals(i) = real(pottarg(i))
    fvalsx(i) = real(fldtarg(1,i))
    fvalsy(i) = real(fldtarg(2,i))
    fvalsz(i) = real(fldtarg(3,i))
  enddo


  !
  !c        get the interpolation matrix
  !
  allocate(xmatc(norder3,norder3))
  call getxmatc(norder,norder3,xmatc)

  !
  !!      compute the chebyshev expansions of the function
  !       and the gradient values
  do i=1,nboxes
    nchild = itree(iptr(4)+i-1)
    if(nchild.eq.0.and.ntargbox(i).gt.0) then
      ii = itree(iptr(6)+i-1)
      call matvec(norder3,norder3,xmatc,fvals(ii),fcoeffs(ii))
      call matvec(norder3,norder3,xmatc,fvalsx(ii),fcoeffsx(ii))
      call matvec(norder3,norder3,xmatc,fvalsy(ii),fcoeffsy(ii))
      call matvec(norder3,norder3,xmatc,fvalsz(ii),fcoeffsz(ii))
    endif
  enddo

  iflag = 1

  maxit = 100
  !
  !!      compute coefficients for fast pnm evaluations required
  !       for evaluating the potential at new targets
  !
  nlege = 100
  lw7 = 40000
  call ylgndrfwini(nlege,wlege,lw7,lused7)

  !$      time1 = omp_get_wtime()

  do iter = 1,maxit
    write (*,*) 'Levels refinement: ', iter
    rmaxerr = 0
    iflag = 0 
    allocate(tails(nboxes),irefineflag(nboxes))
    allocate(tailsx(nboxes),tailsy(nboxes),tailsz(nboxes))

    !
    !!       compute the tail of all chebyshev expansions
    !        and decide which boxes need to be refined further
    !

    !$OMP PARALLEL DO DEFAULT(SHARED), &
    !$OMP& PRIVATE(i,nchild,ii)
    do i=1,nboxes
      nchild = itree(iptr(4)+i-1)
      tails(i) = 0
      tailsx(i) = 0
      tailsy(i) = 0
      tailsz(i) = 0
      irefineflag(i) = 0
      if(nchild.eq.0.and.ntargbox(i).gt.0) then
        ii = itree(iptr(6)+i-1)

        call comptail(norder,norder3,fcoeffs(ii),tails(i))
        call comptail(norder,norder3,fcoeffsx(ii),tailsx(i))
        call comptail(norder,norder3,fcoeffsy(ii),tailsy(i))
        call comptail(norder,norder3,fcoeffsz(ii),tailsz(i))

        if(tails(i).gt.eps.or.tailsx(i).gt.eps.or. &
            tailsy(i).gt.eps.or.tailsz(i).gt.eps) then
          irefineflag(i) = 1
          iflag = 1
        endif
        !              if(tails(i).gt.rmaxerr) rmaxerr = tails(i)
        !              if(tailsx(i).gt.rmaxerr) rmaxerr = tailsx(i)
        !              if(tailsy(i).gt.rmaxerr) rmaxerr = tailsy(i)
        !              if(tailsz(i).gt.rmaxerr) rmaxerr = tailsz(i)
      endif
    enddo
    !$OMP END PARALLEL DO         


    if(iflag.eq.0) goto 2000
    !
    !!          refine the flagged boxes and output the new refined tree
    !           
    call reorganizechebtree_withgrad(itree,ltree,iptr, &
        treecenters,boxsize,nboxes,&
        nlevels,nt,dumtargsort,itstartbox,ntargbox,norder,nt2,&
        targets, fcoeffs,fcoeffsx, &
        fcoeffsy,fcoeffsz, &
        irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
        iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
        itstartboxout,ntargboxout,nt2out,targetsout, &
        fcoeffsout,fcoeffsxout,fcoeffsyout,fcoeffszout)

    !
    !!           deallocate the old tree arrays and point them to the
    !            new ones
    !
    deallocate(itree,treecenters,boxsize,ntargbox,targets, &
        fcoeffs,fcoeffsx,fcoeffsy,fcoeffsz,itstartbox)

    ltree = ltreeout
    nboxes = nboxesout
    nlevels = nlevelsout
    nt2 = nt2out

    allocate(itree(ltree),treecenters(3,nboxes), &
        boxsize(0:nlevels),targets(3,nt2),fcoeffs(nt2), &
        ntargbox(nboxes),fcoeffsx(nt2),fcoeffsy(nt2), &
        fcoeffsz(nt2),itstartbox(nboxes))

    iptr = iptrout
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    do i=1,ltree
      itree(i) = itreeout(i)
    enddo
    !$OMP END PARALLEL DO            

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    do i=1,nboxes
      treecenters(1,i) = treecentersout(1,i)
      treecenters(2,i) = treecentersout(2,i)
      treecenters(3,i) = treecentersout(3,i)
      itstartbox(i) = itstartboxout(i)
      ntargbox(i) = ntargboxout(i)
    enddo
    !$OMP END PARALLEL DO            
    boxsize(0:nlevels) = boxsizeout(0:nlevels)


    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
    do i=1,nt2
      targets(1,i) = targetsout(1,i)
      targets(2,i) = targetsout(2,i)
      targets(3,i) = targetsout(3,i)

      fcoeffs(i) = fcoeffsout(i)
      fcoeffsx(i) = fcoeffsxout(i)
      fcoeffsy(i) = fcoeffsyout(i)
      fcoeffsz(i) = fcoeffszout(i)
    enddo
    !$OMP END PARALLEL DO            

    !            call prinf('after reorg ntargbox=*',ntargbox,nboxes)
    !            call prinf('after reorg itstartbox=*',itstartbox,nboxes)
    !            call prinf('after reorg ilev=*',itree(iptr(2)),nboxes)

    !
    !!            compute the function and the gradient at the new
    !             targets
    !
    npts = 1

    nnnn = int(npts)
    nadapt_flag = int(adapt_flag)
    !$OMP PARALLEL DO DEFAULT(SHARED), &
    !$OMP& PRIVATE(i,ii,j,iii,sigmatmp), &
    !$OMP& PRIVATE(sigma_gradtmp,tradtmp,iboxtarg,ilev,pottmp), &
    !$OMP& PRIVATE(gradtmp,fvalstmp,fvalsxtmp,fvalsytmp,fvalsztmp), &
    !$OMP& PRIVATE(sigmatmp_v,sigma_gradtmp_v_x,sigma_gradtmp_v_y), &
    !$OMP& PRIVATE(sigma_gradtmp_v_z)
    do i=1,nboxes
      if(ibcompflag(i).eq.1) then
        ii = itree(iptr(6)+i-1)
        do j=1,norder3
          iii = ii+j-1

          call function_eval_sigma(FSS_1,targets(1,iii),nnnn,&
              &sigmatmp_v,sigma_gradtmp_v_x,&
              &sigma_gradtmp_v_y,sigma_gradtmp_v_z,nadapt_flag)

          sigmatmp=sigmatmp_v(1)
          sigma_gradtmp(1)=sigma_gradtmp_v_x(1)
          sigma_gradtmp(2)=sigma_gradtmp_v_y(1)
          sigma_gradtmp(3)=sigma_gradtmp_v_z(1)


          tradtmp = 12*sigmatmp
          iboxtarg = 1
          ilev = 0

          call findboxtarg_fulltree(targets(1,iii),tradtmp,&
              iboxtarg,ilev,itreefmm,ipointer,boxsizefmm, &
              treecentersfmm,nboxesfmm,nlevelsfmm)

          call newtargeval(targets(1,iii),sigmatmp, &
              sigma_gradtmp,iboxtarg,ilev,nboxesfmm, &
              nlevelsfmm,boxsizefmm,treecentersfmm,iaddr, &
              itreefmm,ipointer,mnbors,mnlist1,ilevel, &
              rmlexp,ns,sourcesort,ifcharge,chargesort,&
              ifdipole,dipstrsort,dipvecsort, &
              nterms,nlege,wlege,pottmp,gradtmp)


          fvalstmp(j) = pottmp
          fvalsxtmp(j) = gradtmp(1)
          fvalsytmp(j) = gradtmp(2)
          fvalsztmp(j) = gradtmp(3)

        enddo
        call matvec(norder3,norder3,xmatc,fvalstmp,& 
            fcoeffs(ii))
        call matvec(norder3,norder3,xmatc,fvalsxtmp,& 
            fcoeffsx(ii))
        call matvec(norder3,norder3,xmatc,fvalsytmp,& 
            fcoeffsy(ii))
        call matvec(norder3,norder3,xmatc,fvalsztmp,& 
            fcoeffsz(ii))
      endif
    enddo
    !$OMP END PARALLEL DO            
    !
    !!            deallocate the new arrays for the next iteration
    !

    deallocate(tails,irefineflag,itreeout,treecentersout, &
        boxsizeout,ntargboxout,targetsout,fcoeffsout)
    deallocate(itstartboxout)
    deallocate(fcoeffsxout,fcoeffsyout,fcoeffszout)
    deallocate(tailsx,tailsy,tailsz)
    deallocate(itcompflag,ibcompflag)


  enddo
2000 continue   

  !$      time2 = omp_get_wtime()

  call prin2('time in iterating=*',time2-time1,1)
  open(unit=23,file='boxsizes.txt')
  write(23,*) "ilevel, boxsize"       
2100 format(2x,i4,2x,e11.5)       
  do i=0,nlevels
    write(23,2100) i,boxsize(i)
  enddo
  close(23)



end subroutine initialize_feval
!-------------------------------------------     

subroutine f_eval(nt,targets,norder,itree,ltree,nlevels,nboxes,iptr, &
    treecenters,boxsize,nt2,fcoeffs,fcoeffsx,fcoeffsy, &
    fcoeffsz, f, fx, fy, fz, flags)

  implicit double precision (a-h,o-z)
  double precision targets(3,nt)
  integer norder,ltree,nlevels,nboxes,iptr(7),nt2
  integer flags(nt)
  integer itree(ltree)
  double precision treecenters(3,nboxes),boxsize(0:nlevels)
  double precision fcoeffs(nt2),fcoeffsx(nt2),fcoeffsy(nt2),fcoeffsz(nt2)

  double precision f(nt),fx(nt),fy(nt),fz(nt)

  !
  !!       find which box the target lives in
  !
  !
  do i=1,nt
    !write (*,*) 'inside loop: ',i
    call findboxtarg(targets(:,i),iboxtarg,ilevel,itree,iptr,&
        treecenters,nboxes,nlevels)

    if(iboxtarg.eq.-1) istart = -1
    if(iboxtarg.ge.0) istart = itree(iptr(6)+iboxtarg-1)

    if(istart.le.0) then
      flags(i) = 1    
    endif

    if(istart.gt.0) then
      flags(i) = 0
      !
      !!      now evaluate the chebyshev expansion and the gradient
      !       at the new target
      !write (*,*) 'about to start cheb: ',i

      call cheb3deval_withgrad(targets(:,i),boxsize(ilevel), &
          treecenters(1,iboxtarg),norder,fcoeffs(istart),&
          fcoeffsx(istart),fcoeffsy(istart),fcoeffsz(istart),&
          f(i),fx(i),fy(i),fz(i))
      !            write (*,*) i, targets(:,i)
      !            write (*,*) f(i),fx(i),fy(i),fz(i)
    endif
  enddo

  ! 9100 format(3(2x,e11.5))
  !          if(istart.le.0) then
  !             write(*,*) "target in region where no box exists"
  !          write(*,*) "Crashing now"
  !         write(*,*) "error report in targcrashreport.txt"
  !         open(unit=23,file="targcrashreport.txt")
  !         write(23,*) "Target location"
  !         write(23,9100) targets(1),targets(2),targets(3)
  !         write(23,*) "iboxtarg = ",iboxtarg
  !         write(23,*) "ilevel =",ilevel
  !         close(23)
  !         stop
  !      endif
  !    write (*,*) 'detectada box', istart
end subroutine f_eval

!--------------------------------
subroutine f_eval_slow(n,targets,iflag,ns,sources,wts,dipvec,stdev, &
    stdev_grad,pot,grad)
  !
  !!     slow function evaluator
  !
  !      n - number of targets
  !      targets(3,nt) - target location
  !      iflag - targets for which potential needs to be computed
  !      ns - number of sources
  !      sources(3,ns) - source locations
  !      wts(ns) - quadrature weights at the source locations
  !      dipvec(3,ns) - dipole orientation vectors
  !      stdev - sigma value at target location
  !      stdev_grad - gradient of sigma at target location
  !
  !      output
  !      pot(n) - value of the function
  !      grad(3,n) - gradient of the function


  implicit double precision (a-h,o-z)
  integer :: n
  double precision targets(3,n), sources(3,ns),wts(ns),dipvec(3,ns),stdev(n)
  double precision stdev_grad(3,n),pot(n),grad(3,n),ptmp,ftmp(3)
  integer iflag(n)


  iffld = 1
  !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j,i,ptmp,ftmp)   
  do j=1,n

    if(iflag(j).eq.1) then
      pot(j) = 0
      grad(1,j) = 0
      grad(2,j) = 0
      grad(3,j) = 0

      do i=1,ns
        call tpotfld3d_dp(iffld,sources(1,i),wts(i),dipvec(1,i),&
            targets(1,j),stdev(j),stdev_grad(1,j),ptmp,ftmp)
        pot(j) = pot(j) + ptmp
        grad(1,j) = grad(1,j) - ftmp(1)
        grad(2,j) = grad(2,j) - ftmp(2)
        grad(3,j) = grad(3,j) - ftmp(3)
      enddo

      pot(j) = pot(j) - 0.5d0
    endif
  enddo
  !$OMP END PARALLEL DO   

end subroutine f_eval_slow
!-------------------------------------------

subroutine getxmatc(k,k3,xmatc)
  implicit double precision (a-h,o-z)
  double precision xmatc(k3,*)
  double precision uk(k,k),ts(k),wts(k),vk(k,k)

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
  !
  !!      the only difference between findboxtarg_fulltree and 
  !       and findboxtarg is the pointer array and where the children
  !      information is stored in the tree
  !
  implicit double precision (a-h,o-z)
  integer iptr(*),itree(*)
  double precision treecenters(3,*),targtest(3)


  iboxtarg = 1
  ilevel = 0

9100 format(3(2x,e11.5))

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

    iboxtargc = itree(iptr(5)+(iboxtarg-1)*8+ichild-1)

    if(iboxtargc.eq.-1) then
      iboxtarg = -1
      return
      !            write(*,*) "target in region where no box exists"
      !            write(*,*) "Crashing now"
      !            write(*,*) "error report in targcrashreport.txt"
      !            open(unit=23,file='targcrashreport.txt')
      !            write(23,*) "Target location"
      !            write(23,9100) targtest(1),targtest(2),targtest(3)
      !            write(23,*) "last known box where targ found = ",iboxtarg
      !            write(23,*) "Tree centers = "
      !            write(23,9100) treecenters(1,iboxtarg), &
      !                   treecenters(2,iboxtarg),treecenters(3,iboxtarg)
      !            write(23,*) "ilevel =",ilevel
      !            close(23)
      !            stop
    endif

    iboxtarg = iboxtargc

    ilevel = ilev + 1

  enddo
3000 continue      

  return
end subroutine findboxtarg
!------------------------------------------------------------      
subroutine findboxtarg_fulltree(targtest,trad,iboxtarg,ilevel, &
    itree,iptr,boxsize,treecenters,nboxes,nlevels)
  !
  !!      the only difference between findboxtarg_fulltree and 
  !       and findboxtarg is the pointer array and where the children
  !      information is stored in the tree

  implicit double precision (a-h,o-z)
  integer iptr(*),itree(*)
  double precision treecenters(3,nboxes),targtest(3),boxsize(0:nlevels)

  !
  !c       compute relartive corordinates for teh target with
  !        respect to the box
  !
  !
  iboxtarg = 1
  ilevel = 0
  do ilev=0,nlevels-1

    !          exit loop if target is hung at this level
    !
    if(trad.gt.boxsize(ilev)) goto 3000
    nchild = itree(iptr(3)+iboxtarg-1)

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

    iboxtarg = itree(iptr(4)+(iboxtarg-1)*8+ichild-1)
    ilevel = ilev + 1

  enddo
3000 continue      

  return
end subroutine findboxtarg_fulltree
!-----------------------------------------------------

subroutine cheb3deval(targtest,rscale,center,k,fcoeffs,f)
  implicit double precision (a-h,o-z)
  double precision xpols(k+10),ypols(k+10),zpols(k+10)
  double precision fcoeffs(*),targtest(3),center(3)

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

subroutine cheb3deval_withgrad(targtest,rscale,center,k, &
    fcoeffs,fcoeffsx,fcoeffsy,fcoeffsz,f,fx,fy,fz)
  implicit double precision (a-h,o-z)
  double precision xpols(k+10),ypols(k+10),zpols(k+10)
  double precision fcoeffs(*),fcoeffsx(*),fcoeffsy(*),fcoeffsz(*)
  double precision targtest(3),center(3)
  double precision f,fx,fy,fz

  x = (targtest(1) - center(1))/rscale*2
  y = (targtest(2) - center(2))/rscale*2
  z = (targtest(3) - center(3))/rscale*2

  f = 0
  fx = 0
  fy = 0
  fz = 0
  call chebpols(x,k-1,xpols)
  call chebpols(y,k-1,ypols)
  call chebpols(z,k-1,zpols)

  do iz = 1,k
    do iy = 1,k
      do ix = 1,k
        i = (iz-1)*k*k  + (iy-1)*k + ix
        f = f + fcoeffs(i)*xpols(ix)*ypols(iy)*zpols(iz)
        fx = fx + fcoeffsx(i)*xpols(ix)*ypols(iy)*zpols(iz)
        fy = fy + fcoeffsy(i)*xpols(ix)*ypols(iy)*zpols(iz)
        fz = fz + fcoeffsz(i)*xpols(ix)*ypols(iy)*zpols(iz)
      enddo
    enddo
  enddo


  return
end subroutine cheb3deval_withgrad
!--------------------------------------------      
subroutine extracttreesubinfo(itreetmp,ipointer,nboxes,nlevels, &
    itree,iptr,itstartbox,ntargbox)

  !
  !c      itree(iptr(1)) - laddr
  !       itree(iptr(2)) - ilevel
  !       itree(iptr(3)) - iparent
  !       itree(iptr(4)) - nchild
  !       itree(iptr(5)) - ichild
  !       itree(iptr(6)) - pointer to fcoeffs array

  implicit double precision (a-h,o-z)
  integer itree(*),itreetmp(*),iptr(*),ipointer(*),ntargbox(*)
  integer itstartbox(*)


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

    itstartbox(ibox) = istart
    ntargbox(ibox) = iend-istart+1
  enddo
end subroutine extracttreesubinfo
!-----------------------------------------------
subroutine comptail(norder,norder3,fcoeffs,tail)
  implicit double precision (a-h,o-z)
  double precision fcoeffs(*),tail

  tail = 0
  rnorder = norder + 0.0d0

  do ii = 1,norder
    do jj = 1,norder
      do kk = 1,norder

        i = ii+jj+kk-3
        ri = sqrt((ii-1.0d0)**2 + (jj-1.0d0)**2 + (kk-1.0d0)**2)
        if(ri.ge.rnorder) then
          !            if(ii.eq.norder.or.jj.eq.norder.or.kk.eq.norder) then
          !            if(i.ge.norder) then
          iind = (ii-1)*norder*norder+(jj-1)*norder+kk
          tail = tail + abs(fcoeffs(iind))**2
        endif
      enddo
    enddo
  enddo

  tail = sqrt(tail)
end subroutine comptail

subroutine matvec(m,n,a,x,y)
  implicit double precision (a-h,o-z)
  double precision a(m,n),x(n),y(m)

  do i=1,m
    y(i) = 0
    do j=1,n
      y(i) = y(i) + a(i,j)*x(j)
    enddo
  enddo


end subroutine matvec
!-----------------------------------




subroutine reorganizechebtree_withgrad(itree,ltree,iptr, &
    treecenters,boxsize,nboxes, &
    nlevels,nt,dumtarg,itstartbox,ntargbox,norder,nt2,&
    targets,fcoeffs, fcoeffsx, &
    fcoeffsy,fcoeffsz, &
    irefineflag,ibcompflag,itcompflag,itreeout,ltreeout, &
    iptrout,treecentersout,boxsizeout,nboxesout,nlevelsout, &
    itstartboxout,ntargboxout,nt2out,targetsout,fcoeffsout, &
    fcoeffsxout,fcoeffsyout,fcoeffszout)

  implicit double precision (a-h,o-z)
  integer :: nboxes,nlevels,nt2
  integer ::  itree(*)
  integer iptr(7)
  integer :: itstartbox(*)
  integer :: ntargbox(*)
  double precision :: treecenters(3,*)
  double precision :: boxsize(0:nlevels)
  double precision :: targets(3,*)
  double precision :: dumtarg(3,nt)
  double precision :: fcoeffs(*)
  double precision :: fcoeffsx(*)
  double precision :: fcoeffsy(*)
  double precision :: fcoeffsz(*)
  integer norder
  integer :: irefineflag(*)

  integer, intent(out), allocatable :: ibcompflag(:)
  integer, intent(out), allocatable :: itcompflag(:)
  integer, intent(out), allocatable :: itreeout(:)

  integer:: ltreeout,iptrout(7),nboxesout, nlevelsout, nt2out

  integer, intent(out), allocatable :: itstartboxout(:)       
  integer, intent(out), allocatable :: ntargboxout(:)
  double precision, intent(out), allocatable :: targetsout(:,:), &
      fcoeffsout(:), boxsizeout(:),treecentersout(:,:)
  double precision, intent(out), allocatable :: fcoeffsxout(:), &
      fcoeffsyout(:),fcoeffszout(:)

  ! local variables

  integer, allocatable :: titstartbox(:)
  integer, allocatable :: tntargbox(:),tibcompflag(:)
  integer, allocatable :: laddrtail(:,:),tilev(:),tladdr(:,:)
  integer, allocatable :: tiparent(:),tnchild(:),tichild(:,:)
  integer tiptr(7)
  double precision, allocatable :: ttargets(:,:),tboxsize(:)
  double precision, allocatable :: tcenters(:,:)
  double precision, allocatable :: dumtargsort(:,:)
  integer curbox
  integer, allocatable :: iboxtocurbox(:)
  double precision xpts(norder),wts(norder),utmp,vtmp
  integer ntc(8)


  !
  !!       count number of additional boxes to be formed
  !

  nextra = 0
  nlevelsout = nlevels

  !      call prinf('irefineflag=*',irefineflag,nboxes)
  do i=1,nboxes
    if(irefineflag(i).eq.1) nextra = nextra + 8
    nchild = itree(iptr(4)+i-1)
    ilev = itree(iptr(2)+i-1)
    if(ilev.eq.nlevels.and.irefineflag(i).eq.1) &
        nlevelsout = nlevels+1
  enddo


  nboxesout = nboxes + nextra
  !      call prinf('nboxesout max=*',nboxesout,1)



  norder3 = norder*norder*norder

  allocate(dumtargsort(3,nt))
  allocate(titstartbox(nboxesout))
  allocate(tntargbox(nboxesout),tibcompflag(nboxesout))
  allocate(laddrtail(2,0:nlevels+1),tladdr(2,0:nlevels+1))
  allocate(tilev(nboxesout),tnchild(nboxesout),tichild(8,nboxesout))
  allocate(tiparent(nboxesout))
  allocate(tcenters(3,nboxesout))

  dumtargsort = dumtarg

  !
  !!      copy everything into temporary arrays
  !
  do i=0,nlevels+1
    laddrtail(1,i) = 0
    laddrtail(2,i) = -1
  enddo

  do i=0,nlevels
    tladdr(1,i) = itree(2*i+1)
    tladdr(2,i) = itree(2*i+2)
  enddo

  tladdr(1,nlevels+1) = nboxes+1
  tladdr(2,nlevels+1) = nboxes

  !      call prinf('itree=*',itree,ltree)


  do i=1,nboxes
    tilev(i) = itree(iptr(2)+i-1)
    tiparent(i) = itree(iptr(3)+i-1)
    tnchild(i) = itree(iptr(4)+i-1)
    titstartbox(i) = itstartbox(i)
    do j=1,8
      tichild(j,i) = itree(iptr(5)+8*(i-1)+j-1)
    enddo
    tntargbox(i) = ntargbox(i)
    tibcompflag(i) = 0
    if(tnchild(i).eq.0.and.ntargbox(i).gt.0) tibcompflag(i) = 2
    do j=1,3
      tcenters(j,i) = treecenters(j,i)
    enddo
  enddo


  nboxesout = nboxes

  !      call prinf('tntargbox=*',tntargbox,nboxesout)
  !      call prinf('titstartbox=*',titstartbox,nboxesout)
  !      call prin2('dumtarg=*',dumtarg,3*nt)


  do ilev = 0,nlevels
    laddrtail(1,ilev+1) = nboxesout+1
    do ibox = itree(2*ilev+1),itree(2*ilev+2)
      if(irefineflag(ibox).eq.1) then

        call targsorttochildren(ibox,nt,dumtargsort, &
            tcenters(1,ibox),titstartbox,tntargbox,ntc)
        tnchild(ibox) = 0
        tibcompflag(ibox) = -1
        kstart = titstartbox(ibox)
        !               call prinf('ibox=*',ibox,1)
        !               call prin2('tcenters=*',tcenters(1,ibox),3)
        !               call prinf('ntc=*',ntc,8)
        !
        do i=1,8 
          ii = 2
          jj = 2
          if(i.eq.1.or.i.eq.2.or.i.eq.5.or.i.eq.6) ii = 1
          if(i.lt.5) jj = 1
          if(ntc(i).gt.0) then
            nboxesout = nboxesout + 1
            tnchild(ibox) = tnchild(ibox)+1
            tichild(i,ibox) = nboxesout
            tibcompflag(nboxesout) = 1
            tcenters(1,nboxesout)=treecenters(1,ibox)+(-1)**i*&
                boxsize(ilev)/4.0d0
            tcenters(2,nboxesout)=treecenters(2,ibox)+(-1)**ii*&
                boxsize(ilev)/4.0d0
            tcenters(3,nboxesout)=treecenters(3,ibox)+(-1)**jj*&
                boxsize(ilev)/4.0d0
            tnchild(nboxesout) = 0
            do j=1,8
              tichild(j,nboxesout) = -1
            enddo
            tiparent(nboxesout) = ibox
            titstartbox(nboxesout) = kstart 
            tntargbox(nboxesout) = ntc(i)
            tilev(nboxesout) = ilev+1
            tibcompflag(nboxesout) = 1
            kstart = kstart + ntc(i)
          endif
        enddo
      endif
    enddo
    laddrtail(2,ilev+1) = nboxesout
  enddo
  ltreeout = 2*(nlevelsout+1)+12*(nboxesout)

  !      call prinf('nboxesout=*',nboxesout,1)
  !      call prinf('tntargbox=*',tntargbox,nboxesout)
  !      call prinf('titstartbox=*',titstartbox,nboxesout)

  !      call prinf('ltreeout=*',ltreeout,1)


  allocate(itreeout(ltreeout),treecentersout(3,nboxesout))
  allocate(ntargboxout(nboxesout),iboxtocurbox(nboxesout))
  allocate(boxsizeout(0:nlevelsout),ibcompflag(nboxesout))
  allocate(itstartboxout(nboxesout))

  !
  !!      reset the iptr array
  !
  iptrout(1) = 1
  iptrout(2) = iptrout(1)+2*(nlevelsout+1)
  iptrout(3) = iptrout(2)+nboxesout
  iptrout(4) = iptrout(3)+nboxesout
  iptrout(5) = iptrout(4)+nboxesout
  iptrout(6) = iptrout(5)+8*nboxesout
  iptrout(7) = iptrout(6)+nboxesout

  !        call prinf('nlevels=*',nlevels,1)
  !
  !        call prinf('nboxesout=*',nboxesout,1)
  !        call prinf('nlevelsout=*',nlevelsout,1)
  !        call prinf('iptrout=*',iptrout,7)



  !      call prinf('laddrtail=*',laddrtail,2*(nlevels+2))
  !      call prinf('tladdr=*',tladdr,2*(nlevels+1))

  curbox = 1
  do ilev=0,nlevelsout
    itreeout(2*ilev+1) = curbox
    do ibox = tladdr(1,ilev),tladdr(2,ilev)
      !            call prinf('ibox=*',ibox,1)
      !            call prinf('curbox=*',curbox,1)
      itreeout(iptrout(2)+curbox-1) = tilev(ibox)
      itreeout(iptrout(4)+curbox-1) = tnchild(ibox)
      itreeout(iptrout(6)+curbox-1) = -1
      ibcompflag(curbox) = tibcompflag(ibox)
      treecentersout(1,curbox) = tcenters(1,ibox)
      treecentersout(2,curbox) = tcenters(2,ibox)
      treecentersout(3,curbox) = tcenters(3,ibox)
      ntargboxout(curbox) = tntargbox(ibox)
      itstartboxout(curbox) = titstartbox(ibox)
      iboxtocurbox(ibox) = curbox
      curbox = curbox + 1
    enddo

    do ibox = laddrtail(1,ilev),laddrtail(2,ilev)
      !            call prinf('ibox=*',ibox,1)
      !            call prinf('curbox=*',curbox,1)
      itreeout(iptrout(2)+curbox-1) = tilev(ibox)
      itreeout(iptrout(4)+curbox-1) = tnchild(ibox)
      itreeout(iptrout(6)+curbox-1) = -1
      ibcompflag(curbox) = tibcompflag(ibox)
      treecentersout(1,curbox) = tcenters(1,ibox)
      treecentersout(2,curbox) = tcenters(2,ibox)
      treecentersout(3,curbox) = tcenters(3,ibox)
      ntargboxout(curbox) = tntargbox(ibox)
      itstartboxout(curbox) = titstartbox(ibox)
      iboxtocurbox(ibox) = curbox
      curbox = curbox + 1
    enddo
    itreeout(2*ilev+2) = curbox-1
  enddo
  !      call prinf('itreeout=*',itreeout(iptrout(4)),nboxesout)
  !      call prinf('tichild=*',tichild,8*nboxesout)
  !      call prinf('iboxtocurbox=*',iboxtocurbox,nboxesout)


  !
  !!         handle the parent child part of the tree using the
  !          mapping iboxtocurbox

  ictr = 1
  do ibox=1,nboxesout
    if(tiparent(ibox).eq.-1) &
        itreeout(iptrout(3)+iboxtocurbox(ibox)-1) = -1
    if(tiparent(ibox).gt.0) &
        itreeout(iptrout(3)+iboxtocurbox(ibox)-1)=&
        iboxtocurbox(tiparent(ibox))
    do i=1,8
      if(tichild(i,ibox).eq.-1) &
          itreeout(iptrout(5)+(iboxtocurbox(ibox)-1)*8+i-1)=-1
      if(tichild(i,ibox).gt.0) &
          itreeout(iptrout(5)+(iboxtocurbox(ibox)-1)*8+i-1)= &
          iboxtocurbox(tichild(i,ibox))
    enddo

    nchild = itreeout(iptrout(4)+ibox-1)
    if(nchild.eq.0.and.ntargboxout(ibox).gt.0) then
      itreeout(iptrout(6)+ibox-1) = ictr
      ictr = ictr + norder3
    endif
  enddo

  nt2out = ictr-1
  allocate(fcoeffsout(nt2out),targetsout(3,nt2out), &
      itcompflag(nt2out),fcoeffsxout(nt2out),fcoeffsyout(nt2out), &
      fcoeffszout(nt2out))

  !       call prinf('iptrout=*',iptrout,7)
  !
  !       call prinf('itreeout=*',itreeout,ltreeout)
  !       call prinf('nchild=*',itreeout(iptrout(4)),nboxesout)


  !        call prinf('ticompflag=*',tibcompflag,nboxesout)
  !
  !!        fix targets and fcoeffs
  !
  itype = 1
  call chebexps(itype,norder,xpts,utmp,vtmp,wts)

  do i=0,nlevels
    boxsizeout(i) = boxsize(i)
  enddo

  if(nlevelsout.eq.nlevels+1) &
      boxsizeout(nlevelsout) = boxsizeout(nlevels)/2.0d0

  !       call prin2('boxsizeout=*',boxsizeout,nlevelsout+1)   
  do ibox=1,nboxesout
    if(tibcompflag(ibox).eq.2) then
      curbox = iboxtocurbox(ibox)
      !              call prinf('ibox=*',ibox,1)
      !              call prinf('curbox=*',curbox,1)
      ii = itree(iptr(6)+ibox-1)
      jj = itreeout(iptrout(6)+curbox-1)
      !              call prinf('ii=*',ii,1)
      !              call prinf('jj=*',jj,1)
      do i=1,norder3
        fcoeffsout(jj+i-1) = fcoeffs(ii+i-1)
        fcoeffsxout(jj+i-1) = fcoeffsx(ii+i-1)
        fcoeffsyout(jj+i-1) = fcoeffsy(ii+i-1)
        fcoeffszout(jj+i-1) = fcoeffsz(ii+i-1)
        targetsout(1,jj+i-1) = targets(1,ii+i-1)
        targetsout(2,jj+i-1) = targets(2,ii+i-1)
        targetsout(3,jj+i-1) = targets(3,ii+i-1)
        itcompflag(jj+i-1) = 0
      enddo
    endif

    if(tibcompflag(ibox).eq.1) then
      curbox = iboxtocurbox(ibox)
      itarg = itreeout(iptrout(6)+curbox-1)
      ilev = itreeout(iptrout(2)+curbox-1)
      do ii=1,norder
        do jj=1,norder
          do kk=1,norder
            targetsout(1,itarg)=treecentersout(1,curbox)+&
                xpts(kk)*boxsizeout(ilev)/2
            targetsout(2,itarg)=treecentersout(2,curbox)+&
                xpts(jj)*boxsizeout(ilev)/2
            targetsout(3,itarg)=treecentersout(3,curbox)+&
                xpts(ii)*boxsizeout(ilev)/2
            itcompflag(itarg) = 1
            itarg = itarg + 1
          enddo
        enddo
      enddo
    endif
  enddo
end subroutine reorganizechebtree_withgrad


subroutine newtargeval(targets,stdev,stdev_grad,ibox,ilev, &
    nboxes,&
    nlevels,boxsize,centers, iaddr,itree, &
    ipointer,mnbors,mnlist1,ilevel,rmlexp,ns,sourcesort, &
    ifcharge,chargesort,ifdipole,dipstrsort,dipvecsort, &
    nterms,nlege,wlege,pot,grad)
  implicit double precision (a-h,o-z)
  double precision targets(*),boxsize(0:nlevels),rmlexp(*),wlege(*)
  double precision sourcesort(3,*),dipvecsort(3,*),centers(3,*)
  complex *16 chargesort(*), dipstrsort(*)
  integer iaddr(2,*),itree(*),ipointer(32),nterms(0:nlevels)
  integer ilevel(*)
  double precision pot,grad(3),stdev_grad(3)
  complex *16 pottmp,gradtmp(3)

  !
  !!        evaluate local expansion 
  !

  npts = 1
  done = 1
  pi = atan(done)*4

  ifpottarg = 1
  ifgradtarg = 1
  npts = 1

  !     call prinf('ilev=*',ilev,1)

  pottmp = 0
  gradtmp(1) = 0
  gradtmp(2) = 0
  gradtmp(3) = 0

  call l3dtaevalall_trunc(boxsize(ilev),centers(1,ibox), &
      rmlexp(iaddr(2,ibox)),nterms(ilev),nterms(ilev),targets, &
      npts,ifpottarg,pottmp,ifgradtarg,gradtmp,wlege,nlege,ier)

  pottmp = -pottmp/(4*pi)
  gradtmp(1) = -gradtmp(1)/(4*pi)
  gradtmp(2) = -gradtmp(2)/(4*pi)
  gradtmp(3) = -gradtmp(3)/(4*pi)




  itstart = 1
  itend = 1

  nnbors = itree(ipointer(18)+ibox-1)
  do i=1,nnbors
    jbox = itree(ipointer(19)+mnbors*(ibox-1)+i-1)
    jstart = itree(ipointer(10)+jbox-1)
    jend = itree(ipointer(11)+jbox-1)

    call tfmm3dparttarg_direct_newtree(jstart,jend,itstart,itend, &
        sourcesort,ifcharge,chargesort,ifdipole,dipstrsort, &
        dipvecsort, targets, stdev, stdev_grad, ifpottarg, &
        pottmp, ifgradtarg, gradtmp)
  enddo


  nlist1 = itree(ipointer(20)+ibox-1)

  do i=1,nlist1
    jbox = itree(ipointer(21)+mnlist1*(ibox-1)+i-1)
    if(ilevel(ibox).gt.ilevel(jbox)) then
      jstart = itree(ipointer(10)+jbox-1)
      jend = itree(ipointer(11)+jbox-1)

      call tfmm3dparttarg_direct_newtree(jstart, &
          jend, itstart,itend, &
          sourcesort,ifcharge,chargesort,ifdipole,dipstrsort, &
          dipvecsort, targets, stdev, stdev_grad, ifpottarg, &
          pottmp, ifgradtarg, gradtmp)
    endif
  enddo

  pot = real(pottmp)
  grad(1) = real(gradtmp(1))
  grad(2) = real(gradtmp(2))
  grad(3) = real(gradtmp(3))


end subroutine newtargeval

subroutine targsorttochildren(ibox,nt,trg,centers, &
    itstart,ntarg,ntc)
  implicit double precision (a-h,o-z)
  double precision trg(3,*), centers(3)
  integer itstart(*),ntarg(*),ntc(8)
  integer, allocatable :: itarget(:),itargtmp(:)
  double precision, allocatable :: trgtmp(:,:)


  allocate(itarget(nt),itargtmp(nt))
  allocate(trgtmp(3,nt))

  do i=1,nt
    itarget(i) = i
  enddo

  itlast = itstart(ibox) + ntarg(ibox)-1
  !
  !           which child?  1,2,3,4,5,6,7,8? counter nt1,nt2,nt3,nt4,
  !           nt5,nt6,nt7,nt8
  !           The box nomenclature is as follows
  !           3   4       7  8       <--- Looking down in z direction
  !           1   2       5  6
  !
  !           If the parent box center is at the origin then the centers
  !           of box i have co-ordinates
  !           1     x<0,y<0,z<0
  !           2     x>0,y<0,z<0
  !           3     x<0,y>0,z<0
  !           4     x>0,y>0,z<0
  !           5     x<0,y<0,z>0
  !           6     x>0,y<0,z>0
  !           7     x<0,y>0,z>0
  !           8     x>0,y>0,z>0
  !
  i1234 = itstart(ibox)-1
  i5678 = 0
  do ii = 1,ntarg(ibox)
    it = ii+itstart(ibox)-1
    if(trg(3,itarget(it)) - centers(3).lt.0) then
      i1234 = i1234+1
      itarget(i1234) = itarget(it)
    else
      i5678 = i5678 + 1
      itargtmp(i5678) = itarget(it)
    endif
  enddo
  !           Reorder targets to include targets in 5678 in the array
  do i=1,i5678
    itarget(i1234+i) = itargtmp(i)
  enddo

  !           Sort i1234 into i12 and i34         
  i12 = itstart(ibox)-1
  i34 = 0
  do it = itstart(ibox),i1234
    if(trg(2,itarget(it))-centers(2).lt.0) then
      i12 = i12 + 1
      itarget(i12) = itarget(it)
    else
      i34 = i34 + 1
      itargtmp(i34) = itarget(it)
    endif
  enddo
  !           Note at the end of the loop, i12 is where the particles
  !           in part 12 of the box end
  !
  !           Reorder targets to include 34 in the array
  do i=1,i34
    itarget(i12+i) = itargtmp(i)
  enddo

  !           sort i5678 into i56 and i78
  i56 = i1234
  i78 = 0
  do it=i1234+1,itlast
    if(trg(2,itarget(it))-centers(2).lt.0) then
      i56 = i56 + 1
      itarget(i56) = itarget(it)
    else
      i78 = i78 + 1
      itargtmp(i78) = itarget(it)
    endif
  enddo

  !           Reorder sources to include 78 in the array
  do i=1,i78
    itarget(i56+i) = itargtmp(i)
  enddo
  !           End of reordering i5678         

  ntc(1) = 0
  ntc(2) = 0
  ntc(3) = 0
  ntc(4) = 0
  ntc(5) = 0
  ntc(6) = 0
  ntc(7) = 0
  ntc(8) = 0

  !           Sort into boxes 1 and 2
  do it = itstart(ibox),i12
    if(trg(1,itarget(it))-centers(1).lt.0) then
      itarget(itstart(ibox)+ntc(1)) = itarget(it)
      ntc(1) = ntc(1) + 1
    else
      ntc(2) = ntc(2) + 1
      itargtmp(ntc(2)) = itarget(it)
    endif
  enddo
  !           Reorder targets so that sources in 2 are at the
  !           end of this part of the array
  do i=1,ntc(2)
    itarget(itstart(ibox)+ntc(1)+i-1) = itargtmp(i)
  enddo
  !           Sort into boxes 3 and 4
  do it = i12+1, i1234
    if(trg(1,itarget(it))-centers(1).lt.0) then
      itarget(i12+1+ntc(3)) = itarget(it)
      ntc(3) = ntc(3) + 1
    else
      ntc(4) = ntc(4)+1
      itargtmp(ntc(4)) = itarget(it)
    endif
  enddo
  !           Reorder targets so that sources in 4 are at the
  !           end of this part of the array
  do i=1,ntc(4)
    itarget(i12+ntc(3)+i) = itargtmp(i)
  enddo

  !           Sort into boxes 5 and 6
  do it = i1234+1,i56
    if(trg(1,itarget(it))-centers(1).lt.0) then
      itarget(i1234+1+ntc(5)) = itarget(it)
      ntc(5) = ntc(5) + 1
    else
      ntc(6) = ntc(6) + 1
      itargtmp(ntc(6)) = itarget(it)
    endif
  enddo
  !           Reorder targets so that sources in 6 are at the
  !           end of this part of the array
  do i=1,ntc(6)
    itarget(i1234+ntc(5)+i) = itargtmp(i)
  enddo
  !           End of sorting sources into boxes 5 and 6

  !           Sort into boxes 7 and 8
  do it=i56+1,itlast
    if(trg(1,itarget(it))-centers(1).lt.0) then
      itarget(i56+1+ntc(7)) = itarget(it)
      ntc(7) = ntc(7) + 1
    else
      ntc(8) = ntc(8) + 1
      itargtmp(ntc(8)) = itarget(it)
    endif
  enddo
  !           Reorder targets so that sources in 8 are at the
  !           end of the array
  do i=1,ntc(8)
    itarget(i56+ntc(7)+i) = itargtmp(i)
  enddo
  !           End of sorting targets

  !
  !!       reorder the relevant set of targets

  do it=itstart(ibox),itlast
    trgtmp(1,it) = trg(1,itarget(it))
    trgtmp(2,it) = trg(2,itarget(it))
    trgtmp(3,it) = trg(3,itarget(it))
  enddo

  do it=itstart(ibox),itlast
    trg(1,it) = trgtmp(1,it)
    trg(2,it) = trgtmp(2,it)
    trg(3,it) = trgtmp(3,it)
  enddo

  deallocate(itarget,itargtmp,trgtmp)

end subroutine targsorttochildren



subroutine sigma_eval1(xyz,sigma,sigma_grad)
  implicit double precision (a-h,o-z)
  double precision sigma_grad(3),xyz(3)

  rr2 = xyz(1)**2 + xyz(2)**2 + xyz(3)**2

  sigma = 1+(2*rr2)**2

  sigma_grad(1) = 16*rr2*xyz(1)
  sigma_grad(2) = 16*rr2*xyz(2)
  sigma_grad(3) = 16*rr2*xyz(3)

  rscale = 10
  sigma = (0.1 + xyz(1) + xyz(2) + xyz(3))/rscale
  sigma_grad(1) = xyz(1)/rscale
  sigma_grad(2) = xyz(2)/rscale
  sigma_grad(3) = xyz(3)/rscale

  !   sigma=0.40820712299565159d0
  !   sigma_grad(1)=0.0d0
  !   sigma_grad(2)=0.0d0
  !   sigma_grad(3)=0.0d0
end subroutine sigma_eval1




end module
