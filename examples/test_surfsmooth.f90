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

  integer :: N, count,nrefine
  integer :: adapt_flag, ifflatten
  integer :: interp_flag,fmm_flag
  integer :: norder_skel, norder_smooth

  character (len=100) :: nombre, filename,plot_name,name_aux
  character (len=8) :: istr1
  double precision :: x0,y0,z0
  double precision, allocatable :: time_report(:), error_report(:)


  call prini(6,13)
  
  ! order with which to discretize the skeleton patches (pick
  ! something high-order)
  norder_skel = 12
  
  ! order with which to discretize the smooth patches, choose
  ! something reasonable: 4, 6, 8, 10, etc.
  norder_smooth = 6

  ! Specify the numnber of refinements to do starting from 0
  ! nrefine=1  

  ! this is to enable adaptativity (otherwise sigma is constant)
  ! adapt_flag = 0  ->  no adaptivity, mean triangle size
  ! adapt_flag = 1  ->  some adaptivity, alpha form
  ! adapt_flag = 2  ->  full recursive definition, slightly slower
  adapt_flag = 2    

  ! this is to enable FMM (if =1) otherwise ( =0) iterates with stokes
  ! identity (local surface integral + contour integral)
  fmm_flag=1      

  call prinf('. . . printing flags and options*', norder_skel, 0)
  call prinf('norder_skel = *', norder_skel, 1)
  call prinf('norder_smooth = *', norder_smooth, 1)
  ! call prinf('nrefine = *', nrefine, 1)
  call prinf('adapt_flag = *', adapt_flag, 1)

  

  allocate(Geometry1)
  allocate(error_report(nrefine+100))

  !
  ! specify the msh file to read in
  !

  !nombre='./geometries/sphere.msh'
  nombre='./geometries/torus_384.gidmsh'
  !filename='./plot_files/high_genus'

  ! point inside to check Gauss integral
  x0 = 2
  y0 = 0
  z0 = 0


  ! load in the msh file
  call readgeometry(Geometry1, nombre, norder_skel, &
      norder_smooth)

  ifflatten = 0
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
  
  call start_Feval_tree(Feval_stuff_1, Geometry1)
  call funcion_Base_Points(Geometry1)


  !! Esto es para el modo sin FMM
  if (fmm_flag .eq. 0) then
    print *, 'do not run with fmm_flag = 0 !!!'
    stop
    call start_Feval_local(Feval_stuff_1,Geometry1)
  endif

  
  call find_smooth_surface(Geometry1, Feval_stuff_1, adapt_flag)

  print *
  !name_aux=trim(filename)// '_r00.gov'
  !print *, '. . . saving *.gov file: ', trim(name_aux)
  !call record_Geometry(Geometry1,name_aux)

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
  
  ! do count=1,nrefine
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
  ! do count=0,nrefine
  !   write (*,*) 'Refinement num: ',int(count), '  Error: ', &
  !       &error_report(count+1)
  ! enddo

end program smoother








