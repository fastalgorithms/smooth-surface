!
! (c) 2019 Felipe Vico, Manas Rachh, Mike O'Neil
! oneil@cims.nyu.edu
!
!

program smoother

  use ModType_Smooth_Surface
  use Mod_TreeLRD
  use Mod_Plot_Tools_sigma
  use Mod_Fast_Sigma
  use Mod_Feval
  use Mod_Smooth_Surface

  implicit none

  type ( Geometry ), pointer :: Geometry1
  type ( Feval_stuff ), pointer :: Feval_stuff_1

  integer :: N, count,n_refinement
  integer :: adapt_flag
  integer :: interp_flag,fmm_flag
  integer :: norder_skel, norder_smooth

  character (len=100) :: nombre, filename,plot_name,name_aux
  character (len=8) :: istr1
  double precision :: x0,y0,z0
  double precision, allocatable :: time_report(:), error_report(:)


  call prini(6,13)
  
  ! order with which to discretize the skeleton patches (pick
  ! something high-order)
  norder_skel = 16
  
  ! order with which to discretize the smooth patches, choose
  ! something reasonable: 4, 6, 8, 10, etc.
  norder_smooth = 8

  
  n_refinement=1  ! Specify the numnber of refinements to do starting from 0
  adapt_flag=1    ! this is to enable adaptativity (otherwise sigma is constant)

  ! this is to enable FMM (if =1) otherwise ( =0) iterates with stokes
  ! identity (local surface integral + contour integral)
  fmm_flag=1      

  call prinf('. . . printing flags and options*', norder_skel, 0)
  call prinf('norder_skel = *', norder_skel, 1)
  call prinf('norder_smooth = *', norder_smooth, 1)
  call prinf('n_refinement = *', n_refinement, 1)
  call prinf('adapt_flag = *', adapt_flag, 1)

  

  allocate(Geometry1)
  allocate(error_report(n_refinement+100))

  !
  ! specify the msh file to read in
  !
  nombre='./msh_files/sphere_subtract.msh'
  filename='./plot_files/sphere_subtract'
  ! point inside to check Gauss integral
  x0 = 0
  y0 = 0
  z0 = 0


  ! load in the msh file
  call readgeometry(Geometry1, nombre, norder_skel, &
      norder_smooth)

  ! dump out discretization points on the skeleton mesh
  call funcion_skeleton(Geometry1)
  call funcion_normal_vert(Geometry1)
  
  call start_Feval_tree(Feval_stuff_1, Geometry1)
  call funcion_Base_Points(Geometry1)


  !! Esto es para el modo sin FMM
  if (fmm_flag.eq.0) then
    print *, 'do not run with fmm_flag = 0 !!!'
    stop
    call start_Feval_local(Feval_stuff_1,Geometry1)
  endif

  
  call find_smooth_surface(Geometry1, Feval_stuff_1, adapt_flag)

  print *
  name_aux=trim(filename)// '_r00.gov'
  print *, '. . . saving *.gov file: ', trim(name_aux)
  call record_Geometry(Geometry1,name_aux)

  print *
  print *, '. . . checking gauss identity'
  call check_Gauss(Geometry1,x0,y0,z0,error_report(1))

  
  !
  ! plot the smoothed surface
  !
  print *
  print *
  print *, '. . . plotting vtk smoothed geometry'
  plot_name = 'smoothed.vtk'
  call plotsmoothgeometryvtk(Geometry1, plot_name)
  print *, '. . . finished plotting vtk smoothed geometry'


  !
  ! refinement not working properly, must rewrite
  !
  

  !
  ! do some refinement and experiment
  !
  
  ! do count=1,n_refinement
  !   write (*,*) 'Refinement num: ',count
  !   call refine_geometry_smart(Geometry1)
  !   call funcion_Base_Points(Geometry1)
  !   call find_smooth_surface(Geometry1,Feval_stuff_1,adapt_flag)
  !   !write (*,*) 'SAVING .GOV FILE'
  !   !write(istr1,"(I2.2)") count
  !   !name_aux = trim(filename)// '_r'//trim(istr1)//'.gov'
  !   !call record_Geometry(Geometry1,name_aux)
  !   call check_Gauss(Geometry1,x0,y0,z0,error_report(count+1))
  ! enddo


  ! write (*,*) 'FINAL REPORT'
  ! do count=0,n_refinement
  !   write (*,*) 'Refinement num: ',int(count), '  Error: ', &
  !       &error_report(count+1)
  ! enddo

end program smoother








subroutine plotsmoothgeometryvtk(Geometry1, filename)
  use Mod_Smooth_Surface
  implicit none

  type (Geometry) :: Geometry1
  character (len=*) filename

  integer :: umio,count1,count2,flag,n_order_sf, norder_smooth
  integer :: ierror, id, norder, nover, nsub, k, ntri, i, j, ictr
  integer :: ntot, ltot, npols7, info, iii, n, l, nnn, iw
  real (kind = 8) :: us(100000), vs(100000), ws(100000), dcond
  real (kind = 8) :: uv1(10), uv2(10), uv3(10), uv(10), pols(100000)
  real (kind = 8) :: xcoefs(10000), xrhs(10000)
  real (kind = 8) :: ycoefs(10000), yrhs(10000)
  real (kind = 8) :: zcoefs(10000), zrhs(10000)
  real (kind = 8) :: xval, yval, zval, pinv(1000000)

  real (kind = 8), allocatable :: xyzs(:,:,:), uvs(:,:,:)
  real (kind = 8), allocatable :: pmat(:,:), triout(:,:,:)

  double precision :: umatr(100000), vmatr(100000)
  integer :: itype, npols
  
  !
  ! This routien dumps out smoothed geometry into a vtk file,
  ! oversampling the triangles as necessary to show the smoothness
  !
  ! Input:
  !   Geometry1 - the structure containing all info
  !   filename - VTK ASCII filename, should end in .vtk
  !
  ! Output:
  !   the file 'filename' is created and contains vtk info
  !

  !id = 888
  !open(id, FILE=trim(filename),STATUS='REPLACE')

  norder_smooth = Geometry1%norder_smooth
  n_order_sf = Geometry1%n_order_sf

  !
  ! get the nodes here
  !

  call ortho2siexps(itype, norder_smooth, npols, us, vs, &
      umatr, vmatr, ws)

  norder = norder_smooth
  k = npols

  if (n_order_sf .gt. 4**0) nover = 1
  if (n_order_sf .gt. 4**1) nover = 2
  if (n_order_sf .gt. 4**2) nover = 3
  if (n_order_sf .gt. 4**3) nover = 4
  if (n_order_sf .gt. 4**4) nover = 5

  nover = nover + 1
  nsub = 4**nover

  !
  ! now dump out all the info needed for the triangles, compute xtri
  ! coefficients, and resample and plot
  !
  ntri = Geometry1%ntri
  call prinf('in vtk plotter, original ntri = *', ntri, 1)

  
  allocate(xyzs(3,k,ntri))


  ictr = 0
  do i = 1,ntri
    do j = 1,k
      ictr = ictr + 1
      xyzs(1,j,i) = Geometry1%S_smooth(1,ictr)
      xyzs(2,j,i) = Geometry1%S_smooth(2,ictr)
      xyzs(3,j,i) = Geometry1%S_smooth(3,ictr)
    end do
  end do

  
  allocate(uvs(2,3,nsub))
  uvs(1,1,1) = 0
  uvs(2,1,1) = 0
  uvs(1,2,1) = 1
  uvs(2,2,1) = 0
  uvs(1,3,1) = 0
  uvs(2,3,1) = 1


  !
  ! if necessary, recursively subdivide the triangle - first construct
  ! all the uv points
  !
  if (nover .gt. 0) then

    ntot = 1
    do i = 1,nover

      ltot = ntot

      do j = 1,ltot
        uv1(1) = uvs(1,1,j)
        uv1(2) = uvs(2,1,j)
        uv2(1) = uvs(1,2,j)
        uv2(2) = uvs(2,2,j)
        uv3(1) = uvs(1,3,j)
        uv3(2) = uvs(2,3,j)

        uvs(1,1,j) = uv1(1)
        uvs(2,1,j) = uv1(2)
        uvs(1,2,j) = (uv1(1) + uv2(1))/2
        uvs(2,2,j) = (uv1(2) + uv2(2))/2
        uvs(1,3,j) = (uv1(1) + uv3(1))/2
        uvs(2,3,j) = (uv1(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv2(2))/2
        uvs(1,2,ntot) = uv2(1)
        uvs(2,2,ntot) = uv2(2)
        uvs(1,3,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,3,ntot) = (uv2(2) + uv3(2))/2

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,3,ntot) = uv3(1)
        uvs(2,3,ntot) = uv3(2)

        ntot = ntot + 1
        uvs(1,1,ntot) = (uv2(1) + uv3(1))/2
        uvs(2,1,ntot) = (uv2(2) + uv3(2))/2
        uvs(1,2,ntot) = (uv1(1) + uv3(1))/2
        uvs(2,2,ntot) = (uv1(2) + uv3(2))/2
        uvs(1,3,ntot) = (uv1(1) + uv2(1))/2
        uvs(2,3,ntot) = (uv1(2) + uv2(2))/2

      end do
    end do

  end if

  !call prinf('total triangles, ntot = *', ntot, 1)
  !call prin2('uvs = *', uvs, 6*ntot)

  !
  ! now compute the koornwinder expansion of the triangle
  !
  allocate(pmat(k,k))

  !call prin2('us = *', us, k)
  !call prin2('vs = *', vs, k)
  !print *
  !print *

  do i = 1,k
    uv(1) = us(i)
    uv(2) = vs(i)
    call koorn_pols(uv, norder, npols, pols)
    !call prin2('uv = *', uv, 2)
    !call prinf('npols = *', npols, 1)
    if (npols .ne. n_order_sf) then
      call prinf('npols = *', npols, 1)
      call prinf('n_order_sf = *', n_order_sf, 1)
      stop
    end if
    
    !call prin2('pols = *', pols, k)
    !stop
    do j = 1,npols
      pmat(i,j) = pols(j)
    end do
  end do


  call dinverse(npols, pmat, info, pinv)
  !call prinf('after inverse, info = *', info, 1)
  
  !
  ! loop over each triangle, solve for each of the sets of
  ! coefficients, and then evaluate the subsampled triangles
  !

  allocate(triout(3,3,nsub*ntri))
  nnn = 0

  do i = 1,ntri

    do j = 1,k
      xrhs(j) = xyzs(1,j,i)
      yrhs(j) = xyzs(2,j,i)
      zrhs(j) = xyzs(3,j,i)
    end do

    call dmatvec(k, k, pinv, xrhs, xcoefs)
    call dmatvec(k, k, pinv, yrhs, ycoefs)
    call dmatvec(k, k, pinv, zrhs, zcoefs)

    !
    ! now evaluate the new triangle nodes
    !
    do j = 1,nsub
      
      nnn = nnn + 1
      
      do iii = 1,3
        uv(1) = uvs(1,iii,j)
        uv(2) = uvs(2,iii,j)
        call koorn_pols(uv, norder, npols, pols)
        xval = 0
        yval = 0
        zval = 0
        do l = 1,k
          xval = xval + xcoefs(l)*pols(l)
          yval = yval + ycoefs(l)*pols(l)
          zval = zval + zcoefs(l)*pols(l)
        end do
        
        triout(1,iii,nnn) = xval
        triout(2,iii,nnn) = yval
        triout(3,iii,nnn) = zval
        
      end do

      !call prin2('tri = *', triout(1,1,nnn), 9)
      
    end do
    
    
  end do

  call prinf('num triangles plotted = *', nnn, 1)
  
  
  call xtri_vtk_flat(nnn, triout, 'smoothed geometry', filename)

  !close (id)
  return
end subroutine plotsmoothgeometryvtk



subroutine xtri_vtk_flat(ntri, xtri1s, title, filename)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,3,ntri)
  character(len=*) :: title, filename

  character(len=1024) :: dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !
  ! Output:
  !   files which can be executed in matlab to plot the surface
  !
  !

  !
  ! write the vtk plotting script
  !
  iunit1 = 877
  open(unit = iunit1, file=trim(filename), status='replace')

  write(iunit1,'(a)') "# vtk DataFile Version 3.0"
  write(iunit1,'(a)') "vtk output"
  write(iunit1,'(a)') "ASCII"
  !write(iunit1,'(a)') "DATASET POLYDATA"
  write(iunit1,'(a)') "DATASET UNSTRUCTURED_GRID"
  write(iunit1,'(a,i8,a)') "POINTS ", ntri*3, " float"

  fmt2 = "(E11.5,2X,E11.5,2X,E11.5)"
  do i = 1,ntri
    write(iunit1,fmt2) xtri1s(1,1,i), xtri1s(2,1,i), xtri1s(3,1,i)
    write(iunit1,fmt2) xtri1s(1,2,i), xtri1s(2,2,i), xtri1s(3,2,i)
    write(iunit1,fmt2) xtri1s(1,3,i), xtri1s(2,3,i), xtri1s(3,3,i)
  end do


  write(iunit1,'(a,i8,i8)') "CELLS ", ntri, ntri*4

  do i = 1,ntri
    i1 = 3*(i-1) + 1
    write(iunit1,'(a,i8,i8,i8)') "3 ", i1-1, i1, i1+1
  end do

  write(iunit1,'(a,i8)') "CELL_TYPES ", ntri
  do i = 1,ntri
    write(iunit1,'(a)') "5"
  end do

  write(iunit1,'(a,i8)') "POINT_DATA ", ntri*3
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    do j = 1,3
      write(iunit1,'(E11.5)') xtri1s(3,j,i)
    end do
  end do



  write(iunit1,'(a)') ""
  write(iunit1,'(a)') ""
  write(iunit1,'(a,i8)') "CELL_DATA ", ntri
  write(iunit1,'(a)') "SCALARS scalars float 1"
  write(iunit1,'(a)') "LOOKUP_TABLE default"
  do i = 1,ntri
    write(iunit1,'(E13.5)') (xtri1s(3,1,i) + &
        xtri1s(3,2,i) + xtri1s(3,3,i))/3
  end do

  close(iunit1)

  return
end subroutine xtri_vtk_flat
