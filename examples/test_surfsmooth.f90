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

  type ( Geometry ), pointer :: Geometry1 => null ()
  type ( Feval_stuff ), pointer :: Feval_stuff_1 => null ()

  integer :: N, count,nrefine, ifplot
  integer :: adapt_flag, ifflatten
  integer :: interp_flag,fmm_flag
  integer :: norder_skel, norder_smooth

  character (len=100) :: nombre, filename,name_aux
  character (len=21) :: plot_name
  character (len=8) :: istr1
  character (len=2) :: arg_comm
  double precision :: x0,y0,z0
  double precision, allocatable :: time_report(:), error_report(:)
  double precision :: err_skel,rlam,t1,t2,omp_get_wtime


  call prini(6,13)

  
  ! order with which to discretize the skeleton patches (pick
  ! something high-order)
  norder_skel = 6
  
  ! order with which to discretize the smooth patches, choose
  ! something reasonable: 4, 6, 8, 10, etc.
  norder_smooth = 4

  nrefine = 3
  ! nrefine=1  

  ! this is to enable adaptativity (otherwise sigma is constant)
  ! adapt_flag = 0  ->  no adaptivity, mean triangle size
  ! adapt_flag = 1  ->  some adaptivity, alpha form
  ! adapt_flag = 2  ->  full recursive definition, slightly slower
  adapt_flag = 1


  !
  !
  ! rlam flag decides the proportionality value for \sigma_{j}
  ! in relation to triangle diameter
  ! \sigma_{j} = D_{j}/rlam
  !
  rlam = 10
  rlam = 2.5d0
  !rlam = .5d0
  !rlam = 1
  !rlam = 2.5d0

  ! this is to enable FMM (if =1) otherwise ( =0) iterates with stokes
  ! identity (local surface integral + contour integral)
  fmm_flag=1      

  call prinf('. . . printing flags and options*', norder_skel, 0)
  call prinf('norder_skel = *', norder_skel, 1)
  call prinf('norder_smooth = *', norder_smooth, 1)
  ! call prinf('nrefine = *', nrefine, 1)
  call prinf('adapt_flag = *', adapt_flag, 1)


  allocate(Geometry1)
  allocate(error_report(nrefine+1))

  !
  ! specify the msh file to read in
  !

  !nombre='./geometries/sphere.msh'
  !nombre='./geometries/sphere128.gidmsh'
  nombre='./geometries/rcube_refined.gidmsh'
  !nombre='./geometries/rcube.gidmsh'
  !nombre='./geometries/prism_3368.gidmsh'
  !filename='./plot_files/high_genus'

  ! point inside to check Gauss integral
  x0 = 4.5d0
  y0 = 4.5d0
  z0 = 5

  !x0 = .1d0
  !y0 = 0d0
  !z0 = 0



!    nombre='./geometries/msh_files/Round_1.msh'
!    filename='./plot_tools/Round_1'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./geometries/msh_files/simplest_cube_quadratic.msh'
!    filename='./plot_tools/simplest_cube_quadratic'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.5d0

!    nombre='./geometries/msh_files/Round_1.msh'
!    filename='./plot_tools/Round_1'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./geometries/msh_files/Round_2.msh'
!    filename='./plot_tools/Round_2'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='./geometries/msh_files/Cube_substraction.msh'
!    filename='./plot_tools/Cube_substraction'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0

!    nombre='./geometries/msh_files/Multiscale_1.msh'
!    filename='./plot_tools/Multiscale_1'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0

!    nombre='./geometries/msh_files/Multiscale_2.msh'
!    filename='./plot_tools/Multiscale_2'
!!!  point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0

!    nombre='./geometries/msh_files/A380_Final.msh'
!    filename='./plot_tools/A380_Final'
!!!  point inside to check Gauss integral
!    x0=4.0d0
!    y0=0.0d0
!    z0=0.0d0




  ! load in the msh file
  call readgeometry(Geometry1, nombre, norder_skel, &
      norder_smooth)

  ifflatten = 1
  if (ifflatten .eq. 1) then
    call cisurf_quad2flat(Geometry1)
  end if
  
  
  ! plot the skeleton mesh
  plot_name = 'skeleton.vtk'
  call plotskeletonvtk(Geometry1, plot_name)

  
  !call refineskeleton(Geometry1, nrefine)
  
  ! dump out discretization points on the skeleton mesh
  call funcion_skeleton(Geometry1)
  call funcion_normal_vert(Geometry1)

  call start_Feval_tree(Feval_stuff_1, Geometry1, rlam)
  call funcion_Base_Points(Geometry1)


  !! Esto es para el modo sin FMM
  if (fmm_flag .eq. 0) then
    print *, 'do not run with fmm_flag = 0 !!!'
    stop
    call start_Feval_local(Feval_stuff_1,Geometry1)
  endif


  !
  ! before finding smooth surface, compute the Gauss integral on the
  ! skeleton mesh
  !

!  print *
!  print *, '. . . checking gauss identity on skeleton'
!  call check_gauss_skeleton(Geometry1, x0, y0, z0,err_skel)

  
  call find_smooth_surface(Geometry1, Feval_stuff_1, adapt_flag)

  !print *
  !name_aux=trim(filename)// '_r00.gov'
  !print *, '. . . saving *.gov file: ', trim(name_aux)
  !call record_Geometry(Geometry1,name_aux)

  print *
  print *, '. . . checking gauss identity on smooth surface'
  call check_Gauss(Geometry1,x0,y0,z0,error_report(1))

  
  !
  ! plot the smoothed surface
  !
  ifplot = 1
  if (ifplot .eq. 1) then
    print *
    print *
    print *, '. . . plotting vtk smoothed geometry'

    plot_name = 'smoothed.vtk'
    call plotsmoothgeometryvtk(Geometry1, plot_name)
    print *, '. . . finished plotting vtk smoothed geometry'
  end if

!  write (*,*) 'Empezando la parte critica de refinar'
!  read (*,*)
  
  stop


  
  !
  ! refinement not working properly, must rewrite
  !
  

  !
  ! do some refinement and experiment
  !
  
   do count=1,nrefine
!     write (*,*) 'Refinement num: ',count
!     read (*,*)

     call cpu_time(t1)
!$    t1 = omp_get_wtime()     
     call refine_geometry_smart(Geometry1)
     call cpu_time(t2)
!$    t2 = omp_get_wtime()    

     call prin2("Refine geometry time=*",t2-t1,1)

     call cpu_time(t1)
!$    t1 = omp_get_wtime()     
     call funcion_Base_Points(Geometry1)
     call cpu_time(t2)
!$    t2 = omp_get_wtime()    
     call prin2("funcion base time=*",t2-t1,1)


     call cpu_time(t1)
!$    t1 = omp_get_wtime()     
     call find_smooth_surface(Geometry1,Feval_stuff_1,adapt_flag)
     call cpu_time(t2)
!$    t2 = omp_get_wtime()    
     call prin2("find smooth surface time=*",t2-t1,1)
     !write (*,*) 'SAVING .GOV FILE'
     !write(istr1,"(I2.2)") count
     !name_aux = trim(filename)// '_r'//trim(istr1)//'.gov'
     !call record_Geometry(Geometry1,name_aux)
     
     plot_name = 'smoothed1.vtk'
     call plotsmoothgeometryvtk(Geometry1, plot_name)
     call check_Gauss(Geometry1,x0,y0,z0,error_report(count+1))
     write (*,*) 'error_report: ',error_report(count+1)
   enddo

   call prin2('error_report=*',error_report,nrefine+1)
   write (*,*) 'error_report final: ',error_report


    write (*,*) 'FINAL REPORT'
    do count=0,nrefine
        write (*,*) 'Refinement nº: ',int(count,4), '  Error: ', &
        &real(error_report(count+1),4)!, '  Time: ',real(time_report(count+1),4),'sec'
    enddo


  ! write (*,*) 'FINAL REPORT'
  ! do count=0,nrefine
  !   write (*,*) 'Refinement num: ',int(count), '  Error: ', &
  !       &error_report(count+1)
  ! enddo

end program smoother





subroutine check_gauss_skeleton(Geometry1, x0, y0, z0,err_rel)
  use ModType_Smooth_Surface
  implicit none

  !List of calling arguments
  type (Geometry) :: Geometry1
  double precision :: x0,y0,z0

  !List of local variables
  integer :: umio,count1,count2,flag,n_order_sf
  double precision :: F,Ex,Ey,Ez,R,x,y,z,pi,w,nx,ny,nz, done
  double precision :: err_rel

  done = 1
  pi = 4*atan(done)

  F = 0
  write (*,*) 'Num Smooth points',Geometry1%n_Sk_points
  read (*,*)
  do count1=1,Geometry1%n_Sk_points

    x=Geometry1%skeleton_Points(1,count1)
    y=Geometry1%skeleton_Points(2,count1)
    z=Geometry1%skeleton_Points(3,count1)
    w=Geometry1%skeleton_w(count1)

    nx=Geometry1%skeleton_N(1,count1)
    ny=Geometry1%skeleton_N(2,count1)
    nz=Geometry1%skeleton_N(3,count1)


    R=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)

    Ex=(x-x0)/(4*pi*R**3)
    Ey=(y-y0)/(4*pi*R**3)
    Ez=(z-z0)/(4*pi*R**3)
    F = F + (Ex*nx+Ey*ny+Ez*nz)*w
  enddo

  err_rel=abs(F-1)
  call prin2('value of integral = *', F, 1)
  call prin2('relative error = *', err_rel, 1)
    
  return
end subroutine check_gauss_skeleton








