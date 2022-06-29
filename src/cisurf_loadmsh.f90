
!
! (c) 2019 Felipe Vico and Michael O'Neil
! oneil@cims.nyu.edu
!
! This collection of subroutines loads a msh file, puts info into the
! Geometry variable, and records the resulting smooth surface in a
! *.gov file
!

subroutine readgeometry(Geometry1, filename, norder_skel, &
    norder_smooth)
  use ModType_Smooth_Surface
  implicit none
  type (Geometry) :: Geometry1
  character(len=100) :: filename
  integer :: n_order_sf, norder_skel, norder_smooth
  !
  ! This subroutine open a msh file and load the information in a
  ! variable of type Geometry
  !
  ! Input
  !   filename - name of .msh file
  !   n_order_sf - number of points to put on a smooth triangle
  !

  if (norder_skel .gt. 20) then
    call prinf('norder_skel too large = *', norder_skel, 1)
    stop
  end if

  if (norder_smooth .gt. 20) then
    call prinf('norder_smooth too large = *', norder_smooth, 1)
    stop
  end if

  
  if (index(filename,'.msh') > 0) then
    call readmsh(Geometry1, filename, norder_skel, norder_smooth)

  elseif (index(filename,'.gidmsh') > 0) then
    call readgidmsh(Geometry1, filename, norder_skel, norder_smooth)

  elseif (index(filename,'.tri')>0) then
    print *, 'filename = ', trim(filename)
    print *, 'order not converted in .tri files, double check!!'
    call readtri(Geometry1, filename, norder_skel, norder_smooth)

  else
    write (*,*) 'Geometry type not recognized'
    stop
  endif

  return
end subroutine readgeometry





subroutine readmsh(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  !
  ! Input
  !   filename - the file to read
  !   norder_skel - order to discretize the skeleton patches
  !   norder_smooth - order to discretize the smoothed patches
  !
  ! Output
  !   Geometry1 - data structure for geometry
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=100), intent(in) :: filename         !! name of the msh file
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: ierror


  Geometry1%ifflat = 0

  open(UNIT=8, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  read(8,*) aux1,aux2,aux3,m, N

  print *
  print *
  
  write (6,*) 'loading file ', trim(filename)
  write (13,*) 'loading file ', trim(filename)

  call prinf('ntri = *', n, 1)
  !call prinf('npoints = *', m, 1)

  

  nsk = (norder_skel+1)*(norder_skel+2)/2
  !call prinf('nsk = *', nsk, 1)
  !stop

  call prinf('num points on skeleton mesh = *', nsk*n, 1)
  
  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk

  !Geometry1%norder_smooth = norder_smooth
  !Geometry1%nsmooth = n_order_sf

  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  call prinf('num points on smooth mesh = *', nsf*n, 1)
  
  Geometry1%n_order_sf = nsf
  Geometry1%npoints=m
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif

  allocate(Geometry1%Points(3,m))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    read(8,*) Geometry1%Points(1,j),Geometry1%Points(2,j),Geometry1%Points(3,j)
  enddo

  read(8,*) aux1

  do j=1,N
    read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    read(8,*) aux1,aux2,aux3,aux4,aux5,aux6
    Geometry1%Tri(1,j)=aux1
    Geometry1%Tri(2,j)=aux2
    Geometry1%Tri(3,j)=aux3
    Geometry1%Tri(4,j)=aux4
    Geometry1%Tri(5,j)=aux5
    Geometry1%Tri(6,j)=aux6
  enddo
  
  close (8)

  return
end subroutine readmsh





subroutine readgidmsh(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  !
  ! Input
  !   filename - the file to read
  !   norder_skel - order to discretize the skeleton patches
  !   norder_smooth - order to discretize the smoothed patches
  !
  ! Output
  !   Geometry1 - data structure for geometry
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=100), intent(in) :: filename         !! name of the msh file
  character(len=100) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer :: umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: node, nnodes, maxnodes
  integer, allocatable :: elems(:,:)
  integer :: ielem, nelems, maxelems
  double precision :: x, y, z, d, dmin
  double precision, allocatable :: xs(:), ys(:), zs(:)
  integer :: ierror


  Geometry1%ifflat = 0

  open(UNIT=8, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

  write (6,*) 'loading file ', trim(filename)
  write (13,*) 'loading file ', trim(filename)
  print *

  read(8,*) tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  read(8,*) tmp1

  maxnodes = 1000000
  allocate( xs(maxnodes), ys(maxnodes), zs(maxnodes) )

  !
  ! load in the nodes
  !
  
  do i = 1,maxnodes
    read (8,'(a)') tmp1
    if (index(tmp1, 'End Coordinates') > 0) then
      nnodes = i-1
      print *, 'nnodes = ', nnodes
      exit
    end if

    read (tmp1, *) node, xs(i), ys(i), zs(i)
    !print *, 'echoing ', node, x, y, z
    if (i .eq. maxnodes) then
      print *, 'not all nodes loaded, maxnodes = ', maxnodes
      stop
    end if
    
  end do


  ! find the minimun distance between nodes
  if (1 .eq. 0) then

    dmin = 1000

    do i = 1,nnodes
      do j = 1,nnodes
        
        if (i .ne. j) then
          d = (xs(i)-xs(j))**2 + (ys(i)-ys(j))**2 + (zs(i)-zs(j))**2 
          d = sqrt(d)
          if (d .lt. dmin) then
            !print *
            !print *, xs(i), ys(i), zs(i)
            !print *, xs(j), ys(j), zs(j)
            dmin = d
          end if
        end if
        
        
      end do
    end do
    
    !call prin2('minimum dist between nodes = *', dmin, 1)
    !stop
  end if

  
  
  do i = 1,100
    read(8,*) tmp1
    if (index(tmp1, 'Elements') > 0) exit
  end do
  
  maxelems = 1000000
  allocate( elems(6,maxelems) )
  
  do i = 1,maxelems    
    read (8,'(a)') tmp1
    if (index(tmp1, 'End Elements') > 0) then
      nelems = i-1
      print *, 'nelems = ', nelems
      exit
    end if

    read (tmp1, *) ielem, elems(1,i), elems(2,i), elems(3,i), &
        elems(4,i), elems(5,i), elems(6,i)

    !call prinf('elem = *', elems(1,i), 6)
    
  end do


  


  
  n = nelems
  call prinf('ntri = *', n, 1)
  !call prinf('npoints = *', m, 1)

  

  nsk = (norder_skel+1)*(norder_skel+2)/2
  !call prinf('nsk = *', nsk, 1)
  !stop

  call prinf('num points on skeleton mesh = *', nsk*n, 1)
  
  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk

  !Geometry1%norder_smooth = norder_smooth
  !Geometry1%nsmooth = n_order_sf

  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  call prinf('num points on smooth mesh = *', nsf*n, 1)
  
  Geometry1%n_order_sf = nsf

  m = nnodes
  Geometry1%npoints=m
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif

  allocate(Geometry1%Points(3,m))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    !read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    !read(8,*) &
     !Geometry1%Points(1,j),Geometry1%Points(2,j),Geometry1%Points(3,j)
    Geometry1%Points(1,j) = xs(j)
    Geometry1%Points(2,j) = ys(j)
    Geometry1%Points(3,j) = zs(j)
  enddo

  !read(8,*) aux1

  do j=1,N
    !read(8,*) aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    !read(8,*) aux1,aux2,aux3,aux4,aux5,aux6
    do i = 1,6
      Geometry1%Tri(i,j) = elems(i,j)
    end do  
  enddo

  !stop
  
  close(8)

  return
end subroutine readgidmsh








subroutine read_q_gmsh(Geometry1, filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine opens a v2 gmsh file and load the information in a
  !! variable of type Geometry

  !
  ! Input
  !   filename - the file to read
  !   norder_skel - order to discretize the skeleton patches
  !   norder_smooth - order to discretize the smoothed patches
  !
  ! Output
  !   Geometry1 - data structure for geometry
  !

  !List of calling arguments
  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=100), intent(in) :: filename         !! name of the msh file
  character(len=100) :: tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7
  character(len=1000) :: cline
  integer :: n_order_sf, nsk, nsf
  integer  :: norder_skel, norder_smooth
  
  integer :: umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
  integer :: node, nnodes, maxnodes
  integer, allocatable :: elements(:,:), element(:)
  integer :: ielem, nelems, maxelems
  double precision :: x, y, z, d, dmin
  double precision, allocatable :: xs(:), ys(:), zs(:)
  integer :: ierror,iunit,korder,kpols,itype
  integer :: io,numnodes,ind,numelem,nel,ntri,ntag,lll


  Geometry1%ifflat = 0

  iunit = 899

  open(UNIT=iunit, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)

  write (6,*) 'loading file ', trim(filename)
  write (13,*) 'loading file ', trim(filename)
  print *

  korder = 2
  kpols = (korder+1)*(korder+2)/2

  itype = -1
  if (korder .eq. 1) itype = 2
  if (korder .eq. 2) itype = 9
  if (korder .eq. 3) itype = 21
  if (korder .eq. 4) itype = 23
  if (korder .eq. 5) itype = 25
  if (korder .eq. 6) itype = -1



  
  

  do
    read(iunit, *, iostat=io) cline

    if (io .ne. 0) exit

    if (trim(cline) .eq. '$Nodes') then
      print *, 'Reading nodes . . . '
      read(iunit,*) numnodes
      print *, 'Number of nodes = ', numnodes
      print *
      
      allocate(xs(numnodes),ys(numnodes),zs(numnodes))
      do i = 1,numnodes
        read (iunit,*) ind, x, y, z
        xs(i) = x
        ys(i) = y
        zs(i) = z
      end do

    end if

    if (trim(cline) .eq. '$Elements') then
      print *, 'Reading elements . . . '
      read(iunit,*) numelem
      print *, 'Number of elements = ', numelem
      print *

      nel = (korder+1)*(korder+2)/2
      allocate(elements(nel,numelem))

      ntri = 0
      do i = 1,numelem
        
        read(iunit, '(a)', iostat=io) cline
        read (cline,*) ind, ielem, ntag

        if (ielem .eq. itype) then

          ntri = ntri + 1

          lll= 1+1+1+ntag+nel
          allocate(element(lll))
          read(cline,*) element

          do j = 1,kpols
            elements(j,ntri) = element(j+3+ntag)
          end do
          deallocate(element)

        end if
        
      end do

    end if
    
  end do


  
  n = ntri
  call prinf('ntri = *', n, 1)

  nsk = (norder_skel+1)*(norder_skel+2)/2

  call prinf('num points on skeleton mesh = *', nsk*n, 1)
  
  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk


  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  Geometry1%nsmooth = nsf

  call prinf('num points on smooth mesh = *', nsf*n, 1)
  
  Geometry1%n_order_sf = nsf

  m = numnodes
  Geometry1%npoints=m
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif

  allocate(Geometry1%Points(3,m))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    Geometry1%Points(1,j) = xs(j)
    Geometry1%Points(2,j) = ys(j)
    Geometry1%Points(3,j) = zs(j)
  enddo

  do j=1,N
    do i = 1,6
      Geometry1%Tri(i,j) = elements(i,j)
    end do  
  enddo

  
  close(iunit)

  return
end subroutine read_q_gmsh







subroutine readtri(Geometry1,filename, norder_skel, norder_smooth)
  use ModType_Smooth_Surface
  implicit none

  !! This subroutine open a msh file and load the information in a
  !! variable of type Geometry

  type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
  character(len=100), intent(in) :: filename         !! name of the msh file
  integer :: norder_smooth, norder_skel

  integer umio,i,m,N,j,aux1,aux2,aux3,ipointer
  integer :: ierror, nsk,nsf

  ! set the flag for flat vs quadratic
  Geometry1%ifflat = 1

  
  open(UNIT=8, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
  read(8,*) m, N


  print *
  print *
  write (6,*) 'loading file ', trim(filename)
  write (13,*) 'loading file ', trim(filename)
  call prinf('npoints = *', m, 1)
  call prinf('ntri = *', n, 1)


  nsk = (norder_skel+1)*(norder_skel+2)/2
  call prinf('nsk = *', nsk, 1)
  !stop

  Geometry1%norder_skel = norder_skel
  Geometry1%nskel = nsk
  
  Geometry1%norder_smooth = norder_smooth
  nsf = (norder_smooth+1)*(norder_smooth+2)/2
  print *, "nsf=",nsf
  
  Geometry1%n_order_sf=nsf
  Geometry1%nsmooth = nsf
  Geometry1%npoints=m+N*3
  Geometry1%ntri=N
  Geometry1%ntri_sk=N
  Geometry1%n_Sf_points=N*nsf
  Geometry1%n_Sk_points=N*nsk

  if (allocated(Geometry1%Points)) then
    deallocate(Geometry1%Points)
  endif
  if (allocated(Geometry1%Tri)) then
    deallocate(Geometry1%Tri)
  endif
  allocate(Geometry1%Points(3,Geometry1%npoints))
  allocate(Geometry1%Tri(6,N))

  do j=1,m
    read(8,*) Geometry1%Points(1,j),Geometry1%Points(2,j),Geometry1%Points(3,j)
    Geometry1%Points(1:3,j) = Geometry1%Points(1:3,j)*1.0d7
  enddo



  ipointer=m+1

  do j=1,N
    read(8,*) aux1,aux2,aux3
    Geometry1%Tri(1,j)=aux1
    Geometry1%Tri(2,j)=aux2
    Geometry1%Tri(3,j)=aux3
    Geometry1%Tri(4,j)=ipointer
    Geometry1%Tri(5,j)=ipointer+1
    Geometry1%Tri(6,j)=ipointer+2
    Geometry1%Points(:,ipointer)=(Geometry1%Points(:,Geometry1%Tri(1,j))+Geometry1%Points(:,Geometry1%Tri(2,j)))/2.0d0
    Geometry1%Points(:,ipointer+1)=(Geometry1%Points(:,Geometry1%Tri(2,j))+Geometry1%Points(:,Geometry1%Tri(3,j)))/2.0d0
    Geometry1%Points(:,ipointer+2)=(Geometry1%Points(:,Geometry1%Tri(1,j))+Geometry1%Points(:,Geometry1%Tri(3,j)))/2.0d0
    ipointer=ipointer+3
  enddo
  close (8)

  return
end subroutine readtri




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine record_Geometry(Geometry1,filename)
  use ModType_Smooth_Surface
  implicit none

  !List of calling arguments
  type (Geometry), intent(in) :: Geometry1
  character (len=100) filename

  !List of local variables
  integer umio,count1,count2,flag,norder_smooth
  integer :: ierror

  open(8, FILE=filename,STATUS='REPLACE')
  norder_smooth = Geometry1%norder_smooth
  
  write(8,*) norder_smooth
  write(8,*) Geometry1%ntri
!  write(8,*) Geometry1%n_Sf_points
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%S_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%S_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%S_smooth(3,count1)
  enddo

  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%du_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%du_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%du_smooth(3,count1)
  enddo
  
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%dv_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%dv_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
	write(8,*) Geometry1%dv_smooth(3,count1)
  enddo

  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%N_smooth(1,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%N_smooth(2,count1)
  enddo
  do count1=1,Geometry1%n_Sf_points
    write(8,*) Geometry1%N_smooth(3,count1)
  enddo

  close (8)

  return
end subroutine record_Geometry
