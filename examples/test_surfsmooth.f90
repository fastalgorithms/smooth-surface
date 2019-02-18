program Test_6
use Mod_InOut
use ModType_Smooth_Surface
use Mod_Tri_Tools
use Mod_TreeLRD
use Mod_Plot_Tools_sigma
use Mod_Fast_Sigma
use Mod_Feval
use Mod_Smooth_Surface
!use Mod_Plot_Tools_Feval
!use Mod_Plot_Tools_Smooth_Surface

implicit none

!List of Variables
integer,parameter :: seed = 86456   !This is for random stuff
type ( Geometry ), pointer :: Geometry1
type ( Feval_stuff ), pointer :: Feval_stuff_1
!type ( TreeLRD ), pointer :: TreeLRD_1
integer ( kind = 8 ) N,n_order_sk,n_order_sf, count,n_refinement
character ( len=100 ) nombre,filename,plot_name,name_aux
character(len=100) nombre1,nombre2
CHARACTER(LEN=8) :: istr1
real ( kind = 8 ) U(78),V(78),w(78),x0,y0,z0
real ( kind = 8 ), allocatable :: Pts(:,:), sgmas(:)
real ( kind = 8 ) x_min,x_max,y_min,y_max,z_min,z_max
real ( kind = 8 ), allocatable :: F_plot(:,:),targ_vect(:,:),sgma(:),sgma_x(:),sgma_y(:),sgma_z(:)
integer ( kind = 8 ) N_plot,M_plot,count1,count2,icount,adapt_flag,n_targ,n_targets,interp_flag,fmm_flag
INTEGER :: t1, t2,clock_rate, clock_max
real ( kind = 8 ), allocatable :: time_report(:),error_report(:)

    n_order_sk=78   !This is usually 78 always
!Make your choice (number of points per smooth triangle)
    n_order_sf=45
!    n_order_sf=78
    n_refinement=1  ! Specify the numnber of refinements to do starting from 0
    adapt_flag=1    ! this is to enable adaptativity (otherwise sigma is constant)
    interp_flag=1   ! this is to enable the interpolation machinery (otherwise iterates with FMM every time or with stokes identity)
    fmm_flag=1      ! this is to enable FMM (if =1) otherwise ( =0) iterates with stokes identity (local surface integral + contour integral)


    allocate(Geometry1)
    allocate(time_report(n_refinement+1))
    allocate(error_report(n_refinement+1))
    time_report(0)=0.0d0

!Uncomment one of the following geometries
!    nombre='./msh_files/Round_1.msh'
!    filename='./plot_files/Round_1'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./msh_files/Round_2.msh'
!    filename='./plot_files/Round_2'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./msh_files/Cube_substraction.msh'
!    filename='./plot_files/Cube_substraction'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0


!    nombre='./msh_files/Multiscale_1.msh'
!    filename='./plot_files/Multiscale_1'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0

!    nombre='./msh_files/Multiscale_2.msh'
!    filename='./plot_files/Multiscale_2'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0


!    nombre='./msh_files/Round_genus_2.msh'
!    filename='./plot_files/Round_genus_2'
!!!  point inside to check Gauss integral
!    x0=3.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='./msh_files/Multires_cavidad_esfera.msh'   !Not working. Check out the msh file
!    filename='./plot_files/Multires_cavidad_esfera'
!!!  point inside to check Gauss integral
!    x0=-1.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='./msh_files/sphere_union.msh'
!    filename='./plot_files/sphere_union'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./msh_files/sphere_substraction.msh'
!    filename='./plot_files/sphere_substraction'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='./msh_files/substraction2_v2.msh'
!    filename='./plot_files/substraction2_v2'
!!!  point inside to check Gauss integral
!    x0=0.7d0
!    y0=0.3d0
!    z0=1.0d0


!    nombre='cunya1.msh'
!    filename='cunya1'
!!!  point inside to check Gauss integral
!    x0=-0.3d0
!    y0=-0.3d0
!    z0=1.0d0


!    nombre='cunya_local.msh'
!    filename='cunya_local'
!!!  point inside to check Gauss integral
!    x0=0.6d0
!    y0=0.2d0
!    z0=1.0d0


!    nombre='cunya_2_1.msh'
!    filename='cunya_2_1'
!!!  point inside to check Gauss integral
!    x0=-0.2d0
!    y0=0.2d0
!    z0=1.0d0

!    nombre='./msh_files/torus_box.msh'
!    filename='./plot_files/torus_box'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=1.0d0
!    z0=0.0d0



!    nombre='./msh_files/open_cavity_30deg_v2.msh'
!    filename='./plot_files/open_cavity_30deg_v2'
!!!  point inside to check Gauss integral
!    x0=1.5d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='cubo_esfera_multires.msh'
!    filename='cubo_esfera_multires'
!!!  point inside to check Gauss integral
!    x0=-.6d0
!    y0=-.6d0
!    z0=0.6d0


!    nombre='./msh_files/esfera_esfera.msh'
!    filename='./plot_files/esfera_esfera'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=-1.5d0


!    nombre='./msh_files/capsule_multiscale.msh'
!    filename='./plot_files/capsule_multiscale'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.1d0
!    z0=.2d0

!    nombre='sci_fi_3.msh'
!    filename='sci_fi_3'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.01d0
!    z0=.02d0

!    nombre='high_genus_3.msh'
!    filename='high_genus_3'
!!!  point inside to check Gauss integral
!    x0=0.5d0
!    y0=0.5d0
!    z0=0.5d0


!    nombre='parabolic_antenna_v6.msh'
!    filename='parabolic_antenna_v6'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=-0.25d0
!    z0=.00d0


!    nombre='./msh_files/warship3.msh'
!    filename='./plot_files/warship3'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=-.50d0

!    nombre='huge_genus_4.msh'
!    filename='huge_genus_4'
!!!  point inside to check Gauss integral
!    x0=1.5d0
!    y0=1.5d0
!    z0=1.5d0


    nombre='./msh_files/pico_2.msh'
    filename='./plot_files/pico_2'
!    point inside to check Gauss integral
    x0=0d0
    y0=0d0
    z0=2.0d0

!    nombre='two_cavity_filter.msh'
!    filename='two_cavity_filter'
!!!  point inside to check Gauss integral
!    x0=0d0
!    y0=0d0
!    z0=-0.5d0



!    nombre='prueba_cilindro.msh'
!    filename='prueba_cilindro'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=2.0d0

!    nombre='./tri_files/cube0.tri'
!    filename='./plot_files/cube0'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.2d0

!    nombre='./tri_files/Avion_43000_triangles.a.tri'
!    filename='./plot_files/Avion'
!!!  point inside to check Gauss integral
!!    x0=0.0d0
!!    y0=0.0d0
!!    z0=0.2d0

!    nombre='./tri_files/cube_pyramid_conformal_0.a.tri'
!    filename='./plot_files/cube_pyramid_conformal_0'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.2d0


!    nombre='./tri_files/conesphere0.a.tri'
!    filename='./plot_files/conesphere0'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.2d0

!    nombre='./msh_files/simplest_cube_quadratic.msh'
!    filename='./plot_files/simplest_cube_quadratic'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.5d0

!    nombre='./msh_files/antenna_plane_1.msh'
!    filename='./plot_files/antenna_plane_1'
!!!  point inside to check Gauss integral
!    x0=-2.0d0
!    y0=0.0d0
!    z0=0.333d0/2.0d0


!    nombre='./msh_files/antenna_plane_2_v2.msh'
!    filename='./plot_files/antenna_plane_1'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='./msh_files/antenna_plane_1_v2.msh'
!    filename='./plot_files/antenna_plane_1_v2'
!!!  point inside to check Gauss integral
!    x0=-2.0d0
!    y0=0.0d0
!    z0=2.0d0

!    nombre='./msh_files/antenna_plane_2_v2.msh'
!    filename='./plot_files/antenna_plane_2_v2'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./msh_files/antenna_plane_3_v2.msh'
!    filename='./plot_files/antenna_plane_3_v2'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0

!    nombre='./msh_files/Smooth_plane_v0.msh'
!    filename='./plot_files/Smooth_plane_v0'
!!!  point inside to check Gauss integral
!    x0=-5.0d0
!    y0=0.0d0
!    z0=-2.0d0



    call readgeometry(Geometry1,nombre,n_order_sk,n_order_sf)
    call funcion_skeleton(Geometry1,n_order_sk)
    call funcion_normal_vert(Geometry1)
    call start_Feval_tree(Feval_stuff_1,Geometry1)

    call funcion_Base_Points(Geometry1)

!! Esto es para el modo sin FMM
if (fmm_flag.eq.0) then
    call start_Feval_local(Feval_stuff_1,Geometry1)
endif
!! Esto es para el modo sin FMM

    !call plot_sigma(Feval_stuff_1%FSS_1, Geometry1,adapt_flag)
!    call plot_tree_tool(Feval_stuff_1%Tree_local)
!    stop
    call system_clock ( t1, clock_rate, clock_max )
        call find_smooth_surface(Geometry1,Feval_stuff_1,adapt_flag)
    call system_clock ( t2, clock_rate, clock_max )
    time_report(1)=real ( t2 - t1 ) / real ( clock_rate )

    write (*,*) 'SAVING .GOV FILE'
    name_aux=trim(filename)// '_r00.gov'
    call record_Geometry(Geometry1,name_aux)
    call check_Gauss(Geometry1,x0,y0,z0,error_report(1))

!!!THIS COMMAND TRIGGERS THE NEW INTERPOLATION METHOD
    if (interp_flag==1) then
        call system_clock ( t1, clock_rate, clock_max )
        call start_Feval(Feval_stuff_1,Geometry1,adapt_flag)
        call system_clock ( t2, clock_rate, clock_max )
        time_report(0)=real ( t2 - t1 ) / real ( clock_rate )
    endif
!!!IF IT IS COMMENTED, THE OLD FMM NEWTON ITERATION METHOD IS USED

    do count=1,n_refinement
!        read (*,*)
        write (*,*) 'Refinement nº: ',count
        call refine_geometry_smart(Geometry1)
        call funcion_Base_Points(Geometry1)
        call system_clock ( t1, clock_rate, clock_max )
            call find_smooth_surface(Geometry1,Feval_stuff_1,adapt_flag)
        call system_clock ( t2, clock_rate, clock_max )
        time_report(count+1)=real ( t2 - t1 ) / real ( clock_rate )
        write (*,*) 'SAVING .GOV FILE'
        write(istr1,"(I2.2)") count
        name_aux = trim(filename)// '_r'//trim(istr1)//'.gov'
        call record_Geometry(Geometry1,name_aux)
        call check_Gauss(Geometry1,x0,y0,z0,error_report(count+1))
    enddo
    write (*,*) 'FINAL REPORT'
    do count=0,n_refinement
        write (*,*) 'Refinement nº: ',int(count,4), '  Error: ', &
        &real(error_report(count+1),4), '  Time: ',real(time_report(count+1),4),'sec'
    enddo
    write (*,*) 'Interpolation time: ',real(time_report(0),4),'sec'

    nombre1 = trim(filename)//'.rec'

    call record_results(n_refinement,time_report,error_report,nombre1)



  filename='smoothed.vtk'
  call plotSmoothGeometryVTK(Geometry1, filename)


    
stop
end program




subroutine record_results(n_refinement,time_report,error_report,nombre)

integer umio,i,m,n,j
character(len=100) nombre
real *8 time_report(n_refinement+1)
real *8 error_report(n_refinement+1)

    umio=1
    write (*,*) umio
    open(umio, FILE=nombre,STATUS='REPLACE')
    do i=0,n_refinement
        write (umio,*) i,error_report(i+1),time_report(i+1)
    enddo
    close (umio)

return
end







subroutine plotSmoothGeometryVTK(Geometry1, filename)
  use Mod_Smooth_Surface
  implicit none

  type (Geometry) :: Geometry1
  character (len=*) filename

  integer ( kind = 8 ) :: umio,count1,count2,flag,n_order_sf
  integer :: ierror, id, norder, nover, nsub, k, ntri, i, j, ictr
  integer :: ntot, ltot, npols7, npols, info, iii, n, l, nnn, iw
  real (kind = 8) :: us(1000), vs(1000), ws(1000), dcond
  real (kind = 8) :: uv1(10), uv2(10), uv3(10), uv(10), pols(10000)
  real (kind = 8) :: xcoefs(10000), xrhs(10000)
  real (kind = 8) :: ycoefs(10000), yrhs(10000)
  real (kind = 8) :: zcoefs(10000), zrhs(10000)
  real (kind = 8) :: xval, yval, zval, pinv(100000)

  real (kind = 8), allocatable :: xyzs(:,:,:), uvs(:,:,:)
  real (kind = 8), allocatable :: pmat(:,:), triout(:,:,:)
  
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

  n_order_sf = (Geometry1%n_Sf_points)/(Geometry1%ntri)

  if (n_order_sf .eq. 45) then
    norder=8
    nover = 2
    nover = 4
    nsub = 4**nover
    k = 45
    call GaussTri45(us, vs, ws)
  end if
  
  if (n_order_sf .eq. 78) then
    norder=11
    nover = 3
    nover = 5
    nsub = 4**nover
    k = 78
    call GaussTri78(us, vs, ws)
  end if

  !
  ! now dump out all the info needed for the triangles, compute xtri
  ! coefficients, and resample and plot
  !
  ntri = Geometry1%ntri
  call prinf('in vtk plotter, ntri = *', ntri, 1)
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

  call prinf('total triangles, ntot = *', ntot, 1)
  call prin2('uvs = *', uvs, 6*ntot)

  !
  ! now compute the koornwinder expansion of the triangle
  !
  allocate(pmat(k,k))

  call prin2('us = *', us, k)
  call prin2('vs = *', vs, k)
  print *
  print *

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
  call prinf('after inverse, info = *', info, 1)
  
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

    
    !call dgausselim(k, pmat, xrhs, info, xcoefs, dcond)
    !call dgausselim(k, pmat, yrhs, info, ycoefs, dcond)
    !call dgausselim(k, pmat, zrhs, info, zcoefs, dcond)

    call dmatvec(k, k, pinv, xrhs, xcoefs)
    call dmatvec(k, k, pinv, yrhs, ycoefs)
    call dmatvec(k, k, pinv, zrhs, zcoefs)
    
    !call prin2('xcoefs = *', xcoefs, k)
    !call prin2('ycoefs = *', ycoefs, k)
    !call prin2('zcoefs = *', zcoefs, k)
    !stop

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

        !call prin2('xcoefs = *', xcoefs, k)
        !call prin2('ycoefs = *', ycoefs, k)
        !call prin2('zcoefs = *', zcoefs, k)
        
        triout(1,iii,nnn) = xval
        triout(2,iii,nnn) = yval
        triout(3,iii,nnn) = zval
        
      end do

      !call prin2('tri = *', triout(1,1,nnn), 9)
      
    end do
    
    
  end do

  call prinf('total number of triangles = *', nnn, 1)
  call prinf('ntri times nsub = *', ntri*nsub, 1)
  
  
  iw = 33
  call xtri_vtk_flat(iw, nnn, triout, 'smoothed geometry')



  !close (id)
  return
end subroutine plotSmoothGeometryVTK



subroutine xtri_vtk_flat(iw, ntri, xtri1s, title)
  implicit real *8 (a-h,o-z)
  real *8 :: xtri1s(3,3,ntri)
  character(len=*) :: title

  character(len=1024) :: filename, dataname, valsname, imgname
  character(len=1024) :: trisname, vecsname, centname
  character(len=12) :: fmt, fmt3, fmt4
  character(len=25) :: fmt2

  !
  ! This routine plots a sequence of FLAT TRIANGLES with surface
  ! color vals.
  !
  ! Input:
  !   iw - plot number, controls the filenames
  !   ntri - number of flat triangles
  !   xtri1s - full triangle information
  !
  ! Output:
  !   files which can be executed in matlab to plot the surface
  !
  !

  if (iw .lt. 10) then
    fmt = "(A4,I1,A4)"
    fmt3 = "(A8,I1,A4)"
    fmt4 = "(A5,I1,A4)"
  elseif (iw .lt. 100) then
    fmt = "(A4,I2,A4)"
    fmt3 = "(A8,I2,A4)"
    fmt4 = "(A5,I2,A4)"
  elseif (iw .lt. 1000) then
    fmt = "(A4,I3,A4)"
    fmt3 = "(A8,I3,A4)"
    fmt4 = "(A5,I3,A4)"
  elseif (iw .lt. 10000) then
    fmt = "(A4,I4,A4)"
    fmt3 = "(A8,I4,A4)"
    fmt4 = "(A5,I4,A4)"
  end if

  write(filename, fmt) 'plot', iw, '.vtk'

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
