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
