program Smooth_Surface_Quadratic
use some_types
implicit none


interface
subroutine setup_tree_sigma_geometry(Main_box,Geometry1)
    use some_types
        type (Geometry), intent(inout) :: Geometry1
        type (Box), pointer :: Main_box
    end subroutine
end interface

interface
    subroutine find_smooth_surface(Geometry1,alpha,Main_box)
        use some_types
        type (Geometry), intent(inout) :: Geometry1
        real ( kind = 8 ), intent(in) :: alpha
        type (Box), pointer :: Main_box
    end subroutine
end interface

integer,parameter :: seed = 86456   !This is for random stuff
type (Box), pointer :: Main_box
type (Geometry), allocatable :: Geometry1
integer ( kind = 8 ) N,n_order_sk,n_order_sf, count,n_refinement,adaptive_flag
character(len=30) nombre,filename
real ( kind = 8 ) alpha,sgma,Gx,Gy,Gz,X,Y,Z,F,grad_F(3),x0,y0,z0,a_nada,b_nada,c_nada
    n_order_sk=78
    n_order_sf=45
    allocate(Geometry1)

    n_refinement=0
    adaptive_flag=1

!Uncomment one of the following options

!    nombre='Round_1.msh'
!    filename='Round_1.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='Round_2.msh'
!    filename='Round_2.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='Cube_substraction.msh'
!    filename='Cube_substraction.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0


!    nombre='Multiscale_1.msh'
!    filename='Multiscale_1.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0

!    nombre='Multiscale_2.msh'
!    filename='Multiscale_2.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.0d0


!    nombre='Round_genus_2.msh'
!    filename='Round_genus_2.gov'
!    point inside to check Gauss integral
!    x0=3.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='Multires_cavidad_esfera.msh'
!    filename='Multires_cavidad_esfera.gov'
!    point inside to check Gauss integral
!    x0=-1.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='sphere_union.msh'
!    filename='sphere_union.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0

!    nombre='sphere_substraction.msh'
!    filename='sphere_substraction.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='substraction2_v2.msh'
!    filename='substraction2_v2.gov'
!    point inside to check Gauss integral
!    x0=0.7d0
!    y0=0.3d0
!    z0=1.0d0


!    nombre='cunya1.msh'
!    filename='cunya1.gov'
!    point inside to check Gauss integral
!    x0=-0.3d0
!    y0=-0.3d0
!    z0=1.0d0


!    nombre='cunya_local.msh'
!    filename='cunya_local.gov'
!    point inside to check Gauss integral
!    x0=0.6d0
!    y0=0.2d0
!    z0=1.0d0


    nombre='../geometries/cunya/cunya_2_1.msh'
    filename='cunya_2_1.gov'
!    point inside to check Gauss integral
    x0=-0.2d0
    y0=0.2d0
    z0=1.0d0

!    nombre='torus_box.msh'
!    filename='torus_box.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=1.0d0
!    z0=0.0d0



!    nombre='open_cavity_30deg_v2.msh'
!    filename='open_cavity_30deg_v2.gov'
!    point inside to check Gauss integral
!    x0=1.5d0
!    y0=0.0d0
!    z0=0.0d0


!    nombre='cubo_esfera_2.msh'
!    filename='cubo_esfera_2_r1.gov'
!    point inside to check Gauss integral
!    x0=-.6d0
!    y0=-.6d0
!    z0=0.6d0


!    nombre='esfera_esfera.msh'
!    filename='esfera_esfera_r1.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.5d0


!    nombre='capsule_multiscale.msh'
!    filename='capsule_multiscale.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.1d0
!    z0=.2d0

!    nombre='sci_fi_3.msh'
!    filename='sci_fi_3.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.01d0
!    z0=.02d0

!    nombre='parabolic_antenna.msh'
!    filename='parabolic_antenna.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.75d0
!    z0=.00d0


!    nombre='mi_barco_simple_7.msh'
!    filename='mi_barco_simple_7.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=-.50d0

!    nombre='../geometries/big_genus/huge_genus_4.msh'
!    filename='huge_genus_4.gov'
!    point inside to check Gauss integral
!    x0=1.5d0
!    y0=1.5d0
!    z0=1.5d0



!    nombre='simplest_cube_quadratic.msh'
!    filename='simplest_cube_quadratic.gov'
!    point inside to check Gauss integral
!    x0=0.0d0
!    y0=0.0d0
!    z0=1.5d0

!    nombre='jet_fighter.msh'
!    filename='jet_fighter.gov'
!    point inside to check Gauss integral
!    x0=3.0d0
!    y0=0.0d0
!    z0=0.0d0


!    point inside to check Gauss integral

    call leemsh(Geometry1,nombre,n_order_sk,n_order_sf)


    call funcion_skeleton(Geometry1,n_order_sk)
!    do count=1,Geometry1%npoints
!        write (*,*) Geometry1%skeleton_Points(1,count),Geometry1%skeleton_Points(2,count)
!        write (*,*) Geometry1%skeleton_Points(3,count)
!        read (*,*)
!    enddo

    call funcion_normal_vert(Geometry1)
    call find_centroids_sigmas(Geometry1)
    alpha=1.0d0/maxval(Geometry1%sgmas)**2/5.0d0
    write (*,*) 'Alpha: ',alpha

    do count=1,n_refinement
!        write (*,*) 'ntri antes: ', Geometry1%ntri,Geometry1%n_Sk_points
        call refine_geometry(Geometry1,n_order_sf)
!        write (*,*) 'ntri despues: ', Geometry1%ntri,Geometry1%n_Sk_points
!        read (*,*)
    enddo

    call funcion_Base_Points(Geometry1,n_order_sf)

    write (*,*) 'Geometry information '
    write (*,*) 'File name           : ', nombre
    write (*,*) 'Number of triangles : ', Geometry1%ntri
    write (*,*) 'Number of points    : ', Geometry1%npoints
    write (*,*) 'Number of points on the smooth surface: ',Geometry1%n_Sf_points
    write (*,*) 'Number of points on the skeleton: ',Geometry1%n_Sk_points


    call setup_tree_sigma_geometry(Main_box,Geometry1)
    write (*,*) 'No llega nunca'
    call find_smooth_surface(Geometry1,alpha,Main_box)


!    do count=1,Geometry1%n_Sf_points
!        write (*,*) 'Points :', Geometry1%S_smooth(count,1),Geometry1%S_smooth(count,2),Geometry1%S_smooth(count,3)
!        write (*,*) 'Normals :', Geometry1%N_smooth(count,1),Geometry1%N_smooth(count,2),Geometry1%N_smooth(count,3)
!        write (*,*) 'U vect :', Geometry1%ru_smooth(count,1),Geometry1%ru_smooth(count,2),Geometry1%ru_smooth(count,3)
!        write (*,*) 'V vect :', Geometry1%rv_smooth(count,1),Geometry1%rv_smooth(count,2),Geometry1%rv_smooth(count,3)
!        write (*,*) 'w:', Geometry1%w_smooth(count)
!        read (*,*)
!    enddo

    call record_Geometry(Geometry1,filename)
    call check_Gauss(Geometry1,x0,y0,z0)
stop
end program




subroutine setup_tree_sigma_geometry(Main_box,Geometry1)
use some_types
implicit none

interface
    subroutine Setup_tree_add_all(Current_box,Pts,N,radius,max_depth,current_depth,Box_limits)
        use some_types
        integer (kind = 8 ), intent(in) :: N
        type (Box), pointer :: Current_box
        real ( kind = 8 ), intent(in) :: Pts(3,N),Box_limits(3,2),radius
        integer ( kind = 8 ), intent(in) :: max_depth,current_depth
    end subroutine
end interface

interface
    subroutine allocate_tree_points(Current_box)
        use some_types
        type (Box), pointer :: Current_box
    end subroutine
end interface

interface
    subroutine Locate_all_on_tree(Current_box,Pts,N,sgma_v,radius,Box_limits)
        use some_types
        integer ( kind = 8 ), intent(in) :: N
        type (Box), pointer :: Current_box
        real ( kind = 8 ), intent(in) :: Pts(3,N),Box_limits(3,2),sgma_v(N),radius
    end subroutine
end interface

interface
    subroutine Setup_Colleagues(Current_box,current_depth)
        use some_types
        type (Box), pointer :: Current_box
        integer ( kind = 8 ), intent(in) :: current_depth
    end subroutine
end interface


!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
type (Box), pointer :: Main_box

!List of local variables
integer ( kind = 8 ) max_depth, current_depth,count
real ( kind = 8 ) alpha, radius,multiplier,epsil,Box_limits(3,2)

    allocate(Main_box)
    Main_box%n_points=0
    current_depth=0
    max_depth=10
    epsil=1.0d-13
    Box_limits(1,1)=minval(Geometry1%skeleton_Points(1,:))
    Box_limits(1,2)=maxval(Geometry1%skeleton_Points(1,:))
    Box_limits(2,1)=minval(Geometry1%skeleton_Points(2,:))
    Box_limits(2,2)=maxval(Geometry1%skeleton_Points(2,:))
    Box_limits(3,1)=minval(Geometry1%skeleton_Points(3,:))
    Box_limits(3,2)=maxval(Geometry1%skeleton_Points(3,:))
    alpha=1.0d0/maxval(Geometry1%sgmas)**2/5.0d0
    radius=sqrt(log(epsil)/(-alpha))
    radius=radius*10.0d0
!    radius=1000000000000000000000000.0d0
!    write (*,*) 'Radius: ', radius

!    do count=1,Geometry1%ntri_sk
!        write (*,*) Geometry1%sgmas(count)
!    enddo

call Setup_tree_add_all(Main_box,Geometry1%Centroids,Geometry1%ntri_sk,radius,max_depth,current_depth,Box_limits)
    call allocate_tree_points(Main_box)
!    do count=1,Geometry1%n_Sk_points
!        write (*,*) 'sgmas en tree: ',Geometry1%sgmas(count)
!    enddo
call Locate_all_on_tree(Main_box,Geometry1%Centroids,Geometry1%ntri_sk,Geometry1%sgmas,radius,Box_limits)
    call Setup_Colleagues(Main_box,current_depth)


return
end





subroutine leemsh(Geometry1,nombre,n_order_sk,n_order_sf)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
character(len=30), intent(in) :: nombre
integer ( kind = 8 ), intent(in) :: n_order_sk,n_order_sf

!List of local variables

integer ( kind = 8 ) umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7
integer ( kind = 8 ) aux8,aux9,aux10,aux11,aux12,aux13,aux14
integer :: ierror

    open(UNIT=8, FILE=nombre, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
    read(8,*) aux1,aux2,aux3,m, N
    Geometry1%npoints=m
    Geometry1%ntri=N
    Geometry1%ntri_sk=N
    Geometry1%n_Sf_points=N*n_order_sf
    Geometry1%n_Sk_points=N*n_order_sk
    allocate(Geometry1%S_smooth(3,N*n_order_sf))
    allocate(Geometry1%N_smooth(3,N*n_order_sf))
    allocate(Geometry1%skeleton_Points(3,N*n_order_sk))
    allocate(Geometry1%skeleton_w(N*n_order_sk))
    allocate(Geometry1%skeleton_N(3,N*n_order_sk))
    allocate(Geometry1%ru_smooth(3,N*n_order_sf))
    allocate(Geometry1%rv_smooth(3,N*n_order_sf))
    allocate(Geometry1%Base_Points(3,N*n_order_sf))
    allocate(Geometry1%Base_Points_N(3,N*n_order_sf))
    allocate(Geometry1%Base_Points_U(3,N*n_order_sf))
    allocate(Geometry1%Base_Points_V(3,N*n_order_sf))
    allocate(Geometry1%w_smooth(N*n_order_sf))
    allocate(Geometry1%Centroids(3,N))
    allocate(Geometry1%sgmas(N))
    allocate(Geometry1%Points(3,m))
    allocate(Geometry1%Tri(6,N))
    allocate(Geometry1%Normal_Vert(3,m))
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
end


subroutine fun_roots_derivative(h,Geometry1,alpha,F,dF,Main_box,max_step)
use some_types
implicit none

interface
    subroutine eval_density_grad_FMM(Geometry1,targets,alpha,F,grad_F,Main_box,max_step)
    use some_types
        !List of calling arguments
        type (Geometry), intent(in) :: Geometry1
        real ( kind = 8 ) , intent(in) :: targets(3,Geometry1%n_Sf_points),alpha
        real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), grad_F(3,Geometry1%n_Sf_points)
        real ( kind = 8 ), intent(out) :: max_step(Geometry1%n_Sf_points)
        type (Box), pointer :: Main_box
    end subroutine
end interface

!List of calling arguments
type (Geometry), intent(in) :: Geometry1
real ( kind = 8 ), intent(in) :: alpha,h(Geometry1%n_Sf_points)
real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), dF(Geometry1%n_Sf_points)
real ( kind = 8 ), intent(out) :: max_step(Geometry1%n_Sf_points)
type (Box), pointer :: Main_box

!List of local variables
real ( kind = 8 ) r_t(3,Geometry1%n_Sf_points), Norm_N, grad_F(3,Geometry1%n_Sf_points)
integer (kind = 8 ) count

    do count=1,Geometry1%n_Sf_points
        Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+&
        &Geometry1%Base_Points_N(2,count)**2+Geometry1%Base_Points_N(3,count)**2)
r_t(:,count)=Geometry1%Base_Points(:,count)+Geometry1%Base_Points_N(:,count)*h(count)/Norm_N
    enddo
!call eval_density_grad(Geometry1,r_t(1,:),r_t(2,:),r_t(3,:),alpha,F,grad_F)
    call eval_density_grad_FMM(Geometry1,r_t,alpha,F,grad_F,Main_box,max_step)
    do count=1,Geometry1%n_Sf_points
        Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+&
        &Geometry1%Base_Points_N(2,count)**2+Geometry1%Base_Points_N(3,count)**2)
        dF(count)=(grad_F(1,count)*Geometry1%Base_Points_N(1,count)+grad_F(2,count)*Geometry1%Base_Points_N(2,count)&
        &+grad_F(3,count)*Geometry1%Base_Points_N(3,count))/Norm_N
    enddo
return
end




subroutine eval_density_grad_FMM(Geometry1,targets,alpha,F,grad_F,Main_box,max_step)
use some_types
implicit none

interface
    subroutine fast_gaussian_global(Current_box,targ_v,alpha,pot_v,N,Gx,Gy,Gz)
        use some_types
        type (Box), pointer :: Current_box
        integer ( kind = 8 ), intent(in) :: N
        real ( kind = 8 ), intent(in) :: targ_v(3,N),alpha
        real ( kind = 8 ), intent(out) :: pot_v(N),Gx(N),Gy(N),Gz(N)
    end subroutine
end interface

!List of calling arguments
type (Geometry), intent(in) :: Geometry1
real ( kind = 8 ) , intent(in) :: targets(3,Geometry1%n_Sf_points),alpha
real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), grad_F(3,Geometry1%n_Sf_points)
real ( kind = 8 ), intent(out) :: max_step(Geometry1%n_Sf_points)
type (Box), pointer :: Main_box
!! OJO CON LA SALIDA; REAL O COMPLEJA!!!

!List of local variables
real ( kind = 8 ) R, sgma(Geometry1%n_Sf_points), sgma_grad(3,Geometry1%n_Sf_points)
real ( kind = 8 ) H, HH, HHH, prod, w,trads(Geometry1%n_Sf_points)
complex ( kind = 8 ) sigma(Geometry1%n_Sk_points), mu(Geometry1%n_Sk_points)
integer ( kind = 8 ) ier, iprec,ifcharge,ifdipole,ifpottarg,iffldtarg
integer ( kind = 8 ) count,count2
complex ( kind = 8 ) pottarg(Geometry1%n_Sf_points),fldtarg(3,Geometry1%n_Sf_points)
integer ( kind = 4 ) n_sources,n_targets
real ( kind = 8 ) tfmm
    do count2=1,Geometry1%n_Sk_points
        sigma(count2)=0.0d0
        mu(count2)=1.0d0
    enddo
!    write (*,*) 'modificado'
!    do count2=1,Geometry1%n_Sf_points
!        call my_interp_grad(targets(1,count2),targets(2,count2),targets(3,count2),Geometry1,alpha,sgma(count2),&
!        &sgma_grad(1,count2),sgma_grad(2,count2),sgma_grad(3,count2))
!        sgma(count2)=0.12d0
!        sgma_grad(1,count2)=0.0d0
!        sgma_grad(2,count2)=0.0d0
!        sgma_grad(3,count2)=0.0d0
!        trads(count2)=6.0d0*sgma(count2)
!        write (*,*) sgma(count2),sgma_grad(1,count2),sgma_grad(2,count2),sgma_grad(3,count2)
!        read (*,*)
!    enddo
!    write (*,*) 'modificado2'

!!!!!Esta es la parte nueva
!write (*,*) 'No llega nunca2'
call fast_gaussian_global(Main_box,targets,alpha,sgma,Geometry1%n_Sf_points,sgma_grad(1,:),sgma_grad(2,:),sgma_grad(3,:))
!write (*,*) 'numero de puntos skeleton',Geometry1%n_Sk_points
!    sgma=sgma/5.0d0
!    sgma_grad=sgma_grad/5.0d0
    trads=6.0d0*sgma
!    do count2=1,Geometry1%n_Sf_points
!        trads(count2)=10000000000.0d0
!        write (*,*) 'targets: ', targets(:,count2)
!    enddo
    ifcharge=0
    ifdipole=1
    ifpottarg=1
    iffldtarg=1
    iprec=3
    ier=0
    n_sources=Geometry1%n_Sk_points
    n_targets=Geometry1%n_Sf_points
!    write (*,*) 'justo antes de fmm',n_sources,n_targets
!    read (*,*)
    call tfmm3dwrap(ier,iprec,Geometry1%skeleton_Points,Geometry1%skeleton_N,n_sources,&
    &Geometry1%skeleton_w,ifcharge,sigma,ifdipole,mu,targets,&
    &n_targets,trads,sgma,sgma_grad,ifpottarg,pottarg,iffldtarg,fldtarg,tfmm)
!    write (*,*) 'out of fmm'
    do count2=1,Geometry1%n_Sf_points
        pottarg(count2)=pottarg(count2)-0.5d0
        max_step(count2)=sgma(count2)*.5
    enddo
    F=real(pottarg)
    grad_F=-1.0d0*real(fldtarg)
!    do count=1,Geometry1%n_Sf_points
!        write (*,*) pottarg(count),grad_F(1,count),grad_F(2,count),grad_F(3,count)
!        read (*,*)
!    enddo
return
end


subroutine eval_density_grad(Geometry1,X_t,Y_t,Z_t,alpha,F,grad_F)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(in) :: Geometry1
real ( kind = 8 ) , intent(in) :: X_t(Geometry1%n_Sf_points),Y_t(Geometry1%n_Sf_points),Z_t(Geometry1%n_Sf_points),alpha
real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), grad_F(3,Geometry1%n_Sf_points)

!List of local variables
real ( kind = 8 ) R, sgma, Gx, Gy, Gz, X_s, Y_s, Z_s, N_x, N_y, N_z
real ( kind = 8 ) H, HH, HHH, prod, w
integer ( kind = 8 ) count,count2
    do count2=1,Geometry1%n_Sf_points
        call my_interp_grad(X_t(count2),Y_t(count2),Z_t(count2),Geometry1,alpha,sgma,Gx,Gy,Gz)
        F(count2)=-0.5d0
        grad_F(1,count2)=0.0d0
        grad_F(2,count2)=0.0d0
        grad_F(3,count2)=0.0d0
        do count=1,Geometry1%n_Sk_points
            X_s=Geometry1%skeleton_Points(1,count)
            Y_s=Geometry1%skeleton_Points(2,count)
            Z_s=Geometry1%skeleton_Points(3,count)
            N_x=Geometry1%skeleton_N(1,count)
            N_y=Geometry1%skeleton_N(2,count)
            N_z=Geometry1%skeleton_N(3,count)
            w=Geometry1%skeleton_w(count)
            R=sqrt((X_s-X_t(count2))**2+(Y_s-Y_t(count2))**2+(Z_s-Z_t(count2))**2)
            call my_erf_like_v44(R,sgma,H,HH,HHH)
!            write (*,*) 'Hs 1: ',R,sgma,H,HH,HHH
!            call my_erf_like_v3(R,sgma,H,HH,HHH)
!            write (*,*) 'Hs 2: ',R,sgma,H,HH,HHH
!            read (*,*)
            prod=(X_s-X_t(count2))*N_x+(Y_s-Y_t(count2))*N_y+(Z_s-Z_t(count2))*N_z
            F(count2)=F(count2)+H*prod*w
            grad_F(1,count2)=grad_F(1,count2)+((HH*(X_t(count2)-X_s)*prod+HHH*prod*Gx-H*N_x)*w)
            grad_F(2,count2)=grad_F(2,count2)+((HH*(Y_t(count2)-Y_s)*prod+HHH*prod*Gy-H*N_y)*w)
            grad_F(3,count2)=grad_F(3,count2)+((HH*(Z_t(count2)-Z_s)*prod+HHH*prod*Gz-H*N_z)*w)
        enddo
    enddo
return
end



subroutine my_erf_like_v44(r,sgma,H,HH,HHH)
implicit none

!List of calling arguments
real ( kind = 8 ) , intent(in) :: r,sgma
real ( kind = 8 ), intent(out) ::  H, HH, HHH

!List of local variables
real ( kind = 8 ) pi, my_exp, denom,x(15),w(15)
integer ( kind = 8) count

    pi=3.141592653589793238462643383d0
    my_exp=exp(-r**2/(2*sgma**2))
    HHH=-1/(4*pi)*sqrt(2/pi)*(1/sgma**4)*my_exp

    if (r==0.0d0) then
        H=sqrt(2.0d0)/(12*sgma**3*sqrt(pi**3))
        HH=-1.0d0*sqrt(2.0d0)/(20.0d0*sgma**5.0d0*sqrt(pi**3))
    else if (r/sgma<0.9d0) then
        x(1)=6.003740989757256e-03
        x(2)=3.136330379964697e-02
        x(3)=7.589670829478640e-02
        x(4)=1.377911343199150e-01
        x(5)=2.145139136957306e-01
        x(6)=3.029243264612183e-01
        x(7)=3.994029530012828e-01
        x(8)=5.000000000000000e-01
        x(9)=6.005970469987173e-01
        x(10)=6.970756735387817e-01
        x(11)=7.854860863042694e-01
        x(12)=8.622088656800850e-01
        x(13)=9.241032917052137e-01
        x(14)=9.686366962003530e-01
        x(15)=9.939962590102427e-01
        w(1)=1.537662099805856e-02
        w(2)=3.518302374405402e-02
        w(3)=5.357961023358594e-02
        w(4)=6.978533896307718e-02
        w(5)=8.313460290849692e-02
        w(6)=9.308050000778104e-02
        w(7)=9.921574266355579e-02
        w(8)=1.012891209627806e-01
        w(9)=9.921574266355579e-02
        w(10)=9.308050000778104e-02
        w(11)=8.313460290849692e-02
        w(12)=6.978533896307718e-02
        w(13)=5.357961023358594e-02
        w(14)=3.518302374405402e-02
        w(15)=1.537662099805856e-02
        x=x*r/sgma/sqrt(2.0d0)
        w=w*r/sgma/sqrt(2.0d0)
        H=0.0d0
        HH=0.0d0
        do count=1,15
!            x(count)=x(count)*r/sgma/sqrt(2.0d0)
!            w(count)=w(count)*r/sgma/sqrt(2.0d0)
            H=H+(x(count)**2)*exp(-x(count)**2)/(sqrt(pi**3))*w(count)/(r**3);
            HH=HH+(x(count)**4)*exp(-x(count)**2)*-2.0d0*w(count)/(r**5*sqrt(pi**3));
        enddo
    else
        H=(erf((sqrt(2.0d0)*r)/(2*sgma))/(4*r**2*pi) - (sqrt(2.0d0)*my_exp)/(4*sgma*r*sqrt(pi**3)))/r
        denom=(4.*r**5*sgma**3*sqrt(pi**3))
        HH=(sqrt(2.0d0)*r**3*my_exp - 3*sgma**3*sqrt(pi)*erf((sqrt(2.0d0)*r)/(2*sgma)) + 3*sqrt(2.0d0)*r*sgma**2*my_exp)/denom
    end if
return
end



subroutine find_centroids_sigmas(Geometry1)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1

!List of local variables
integer ( kind = 8 ) count
real ( kind = 8 ) lambda,d1,d2,d3,dmax,sgma_mean
real ( kind = 8 ) P1(3),P2(3),P3(3),Pc(3),Pd1(3),Pd2(3),Pd3(3)
    lambda=0.1d0   !!Mucho ojo con esto tío, antes era lambda=0.03d0
    sgma_mean=0.0d0
    do count=1,Geometry1%ntri
        P1=Geometry1%Points(:,Geometry1%Tri(1,count))
        P2=Geometry1%Points(:,Geometry1%Tri(2,count))
        P3=Geometry1%Points(:,Geometry1%Tri(3,count))
        Pc=(P1+P2+P3)/3.0d0
        Pd1=P2-P1
        Pd2=P3-P2
        Pd3=P3-P1
        d1=sqrt(Pd1(1)**2+Pd1(2)**2+Pd1(3)**2)
        d2=sqrt(Pd2(1)**2+Pd2(2)**2+Pd2(3)**2)
        d3=sqrt(Pd3(1)**2+Pd3(2)**2+Pd3(3)**2)
        dmax=max(d1,d2,d3)
        Geometry1%Centroids(:,count)=Pc
        Geometry1%sgmas(count)=lambda*dmax
    enddo

!Esta parte es por si quiero homogeneizar. Mucho ojo, esto elimina la adaptatividad!!!
!    sgma_mean=sum(Geometry1%sgmas)/Geometry1%ntri
!    do count=1,Geometry1%ntri
!        Geometry1%sgmas(count)=sgma_mean
!    enddo
!Hasta aquí elimina la adaptatividad


return
end

subroutine my_interp_grad(X,Y,Z,Geometry1,alpha,sgma,Gx,Gy,Gz)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(in) :: Geometry1
real ( kind = 8 ), intent(in) :: alpha, X, Y, Z
real ( kind = 8 ), intent(out) :: sgma,Gx,Gy,Gz

!List of local variables
real ( kind = 8 ) F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,di,my_exp
integer ( kind = 8 ) count
    F=0.0d0
    D=0.0d0
    dF_x=0.0d0
    dF_y=0.0d0
    dF_z=0.0d0
    dD_x=0.0d0
    dD_y=0.0d0
    dD_z=0.0d0
    do count=1,Geometry1%ntri
        di=sqrt((Geometry1%Centroids(1,count)-X)**2+(Geometry1%Centroids(2,count)-Y)**2+(Geometry1%Centroids(3,count)-Z)**2)
        my_exp=exp(-alpha*di**2)
        F=F+Geometry1%sgmas(count)*my_exp
        dF_x=dF_x+(-2*alpha)*Geometry1%sgmas(count)*my_exp*(X-Geometry1%Centroids(1,count))
        dF_y=dF_y+(-2*alpha)*Geometry1%sgmas(count)*my_exp*(Y-Geometry1%Centroids(2,count))
        dF_z=dF_z+(-2*alpha)*Geometry1%sgmas(count)*my_exp*(Z-Geometry1%Centroids(3,count))
        dD_x=dD_x+(-2*alpha)*my_exp*(X-Geometry1%Centroids(1,count))
        dD_y=dD_y+(-2*alpha)*my_exp*(Y-Geometry1%Centroids(2,count))
        dD_z=dD_z+(-2*alpha)*my_exp*(Z-Geometry1%Centroids(3,count))
        D=D+my_exp
    enddo
    Gx=(dF_x*D-F*dD_x)/D**2
    Gy=(dF_y*D-F*dD_y)/D**2
    Gz=(dF_z*D-F*dD_z)/D**2
    sgma=F/D
return
end



subroutine funcion_normal_vert(Geometry1)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1

!List of local variables
type (My_cell) My_cell1
integer ( kind = 8 ) count, n_order
real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
real ( kind = 8 ) U(3),V(3)
real ( kind = 8 ) F_x(3),F_y(3),F_z(3),dS(3)
real ( kind = 8 ) nP_x(3),nP_y(3),nP_z(3)
real ( kind = 8 ) current_Normal_vert(3)

    U= (/0.0d0, 1.0d0, 0.0d0/)
    V= (/0.0d0, 0.0d0, 1.0d0/)
    n_order=3
    My_cell1%n_Cell=Geometry1%npoints
    allocate(My_cell1%Var_Mat(Geometry1%npoints))
    do count=1,Geometry1%npoints
        My_cell1%Var_Mat(count)%n_Mat=0
        My_cell1%Var_Mat(count)%current_n_Mat=1
    enddo
    do count=1,Geometry1%ntri
        My_cell1%Var_Mat(Geometry1%Tri(1,count))%n_Mat=My_cell1%Var_Mat(Geometry1%Tri(1,count))%n_Mat+1
        My_cell1%Var_Mat(Geometry1%Tri(2,count))%n_Mat=My_cell1%Var_Mat(Geometry1%Tri(2,count))%n_Mat+1
        My_cell1%Var_Mat(Geometry1%Tri(3,count))%n_Mat=My_cell1%Var_Mat(Geometry1%Tri(3,count))%n_Mat+1
    enddo
    do count=1,Geometry1%npoints
        allocate(My_cell1%Var_Mat(count)%Mat(3,My_cell1%Var_Mat(count)%n_Mat))
    enddo
    do count=1,Geometry1%ntri
        P1=Geometry1%Points(:,Geometry1%Tri(1,count))
        P2=Geometry1%Points(:,Geometry1%Tri(2,count))
        P3=Geometry1%Points(:,Geometry1%Tri(3,count))
        P4=Geometry1%Points(:,Geometry1%Tri(4,count))
        P5=Geometry1%Points(:,Geometry1%Tri(5,count))
        P6=Geometry1%Points(:,Geometry1%Tri(6,count))
        call eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,n_order)
My_cell1%Var_Mat(Geometry1%Tri(1,count))%Mat(:,My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat)=(/nP_x(1), nP_y(1), nP_z(1)/)
My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat=My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat+1
My_cell1%Var_Mat(Geometry1%Tri(2,count))%Mat(:,My_cell1%Var_Mat(Geometry1%Tri(2,count))%current_n_Mat)=(/nP_x(1), nP_y(1), nP_z(1)/)
My_cell1%Var_Mat(Geometry1%Tri(2,count))%current_n_Mat=My_cell1%Var_Mat(Geometry1%Tri(2,count))%current_n_Mat+1
My_cell1%Var_Mat(Geometry1%Tri(3,count))%Mat(:,My_cell1%Var_Mat(Geometry1%Tri(3,count))%current_n_Mat)=(/nP_x(1), nP_y(1), nP_z(1)/)
My_cell1%Var_Mat(Geometry1%Tri(3,count))%current_n_Mat=My_cell1%Var_Mat(Geometry1%Tri(3,count))%current_n_Mat+1
    enddo
    do count=1,Geometry1%npoints
        if (My_cell1%Var_Mat(count)%n_Mat>0) then
            call mimean(My_cell1%Var_Mat(count)%Mat,My_cell1%Var_Mat(count)%n_Mat,current_Normal_vert)
            Geometry1%Normal_Vert(:,count)=current_Normal_vert
        else
            Geometry1%Normal_Vert(:,count)=(/0.d0, 0.d0, 0.d0/)
        endif
    enddo
return
end


subroutine mimean(All_Normals,n_Normals,Current_Normal)
use some_types
implicit none

!List of calling arguments
integer ( kind = 8 ), intent(in) :: n_Normals
real ( kind = 8 ), intent(in) :: All_Normals(3,n_Normals)
real ( kind = 8 ), intent(out) :: current_Normal(3)

!List of local variables
integer ( kind = 8 ) count
real ( kind = 8 ) coef(n_Normals),sum_coef

    current_Normal=(/ 0.0d0, 0.0d0, 0.0d0 /)
    sum_coef=0
    do count=1,n_Normals
        call find_cos(All_Normals(:,count),All_Normals,coef(count),n_Normals)
        sum_coef=sum_coef+coef(count)
    enddo
    do count=1,n_Normals
        current_Normal=current_Normal+All_Normals(:,count)*coef(count)
    enddo
    current_Normal=current_normal/sum_coef
return
end


subroutine find_cos(v,All_Normals,coef,n_Normals)
use some_types
implicit none

!List of calling arguments
integer ( kind = 8 ), intent(in) :: n_Normals
real ( kind = 8 ), intent(out) :: coef
real ( kind = 8 ), intent(in) :: All_Normals(3,n_Normals), v(3)

!List of local variables
integer ( kind = 8 ) count
real ( kind = 8 ) n_times, tol, my_dot, norm1, norm2
    n_times=0
    tol=1.0d-10
    norm1=sqrt(v(1)**2+v(2)**2+v(3)**2)
    do count=1,n_Normals
        norm2=sqrt(All_Normals(1,count)**2+All_Normals(2,count)**2+All_Normals(3,count)**2)
        my_dot=dot_product(v,All_Normals(:,count))/(norm1*norm2)
        if (my_dot<0) then
            my_dot=0
        endif
        n_times=n_times+my_dot**2
    enddo
    coef=1/n_times
return
end




subroutine funcion_skeleton(Geometry1,n_order_sk)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
integer ( kind = 8 ), intent(in) :: n_order_sk


!List of local variables
real ( kind = 8 ) U(n_order_sk),V(n_order_sk),w(n_order_sk)
real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
real ( kind = 8 ) F_x(n_order_sk),F_y(n_order_sk),F_z(n_order_sk),dS(n_order_sk)
real ( kind = 8 ) nP_x(n_order_sk),nP_y(n_order_sk),nP_z(n_order_sk)

integer ( kind = 8 ) count
    if (n_order_sk==45) then
        call GaussTri45(U,V,w)
    else if (n_order_sk==78) then
        call GaussTri78(U,V,w)
    end if
    do count=1,Geometry1%ntri_sk
        P1=Geometry1%Points(:,Geometry1%Tri(1,count))
        P2=Geometry1%Points(:,Geometry1%Tri(2,count))
        P3=Geometry1%Points(:,Geometry1%Tri(3,count))
        P4=Geometry1%Points(:,Geometry1%Tri(4,count))
        P5=Geometry1%Points(:,Geometry1%Tri(5,count))
        P6=Geometry1%Points(:,Geometry1%Tri(6,count))
        call eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,n_order_sk)
        Geometry1%skeleton_Points(1,(count-1)*n_order_sk+1:(count)*n_order_sk)=F_x
        Geometry1%skeleton_Points(2,(count-1)*n_order_sk+1:(count)*n_order_sk)=F_y
        Geometry1%skeleton_Points(3,(count-1)*n_order_sk+1:(count)*n_order_sk)=F_z
        Geometry1%skeleton_w((count-1)*n_order_sk+1:(count)*n_order_sk)=w*dS
        Geometry1%skeleton_N(1,(count-1)*n_order_sk+1:(count)*n_order_sk)=nP_x
        Geometry1%skeleton_N(2,(count-1)*n_order_sk+1:(count)*n_order_sk)=nP_y
        Geometry1%skeleton_N(3,(count-1)*n_order_sk+1:(count)*n_order_sk)=nP_z
    enddo
return
end


subroutine funcion_Base_Points(Geometry1,n_order_sf)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
integer ( kind = 8 ), intent(in) :: n_order_sf


!List of local variables
real ( kind = 8 ) U(n_order_sf),V(n_order_sf),w(n_order_sf)
real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),N1(3),N2(3),N3(3)
real ( kind = 8 ) F_x(n_order_sf),F_y(n_order_sf),F_z(n_order_sf),dS(n_order_sf)
real ( kind = 8 ) nP_x(n_order_sf),nP_y(n_order_sf),nP_z(n_order_sf)
real ( kind = 8 ) U_x(n_order_sf),U_y(n_order_sf),U_z(n_order_sf)
real ( kind = 8 ) V_x(n_order_sf),V_y(n_order_sf),V_z(n_order_sf)

integer ( kind = 8 ) count
    if (n_order_sf==45) then
        call GaussTri45(U,V,w)
    else if (n_order_sf==78) then
        call GaussTri78(U,V,w)
    end if

    do count=1,Geometry1%ntri
        P1=Geometry1%Points(:,Geometry1%Tri(1,count))
        P2=Geometry1%Points(:,Geometry1%Tri(2,count))
        P3=Geometry1%Points(:,Geometry1%Tri(3,count))
        P4=Geometry1%Points(:,Geometry1%Tri(4,count))
        P5=Geometry1%Points(:,Geometry1%Tri(5,count))
        P6=Geometry1%Points(:,Geometry1%Tri(6,count))
        N1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
        N2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
        N3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
    call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_sf)
        nP_x=N1(1)*(1.0d0-U-V)+N2(1)*U+N3(1)*V
        nP_y=N1(2)*(1.0d0-U-V)+N2(2)*U+N3(2)*V
        nP_z=N1(3)*(1.0d0-U-V)+N2(3)*U+N3(3)*V
        Geometry1%Base_Points(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_x
        Geometry1%Base_Points(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_y
        Geometry1%Base_Points(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=F_z
        Geometry1%Base_Points_N(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_x
        Geometry1%Base_Points_N(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_y
        Geometry1%Base_Points_N(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=nP_z
        Geometry1%Base_Points_U(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_x
        Geometry1%Base_Points_U(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_y
        Geometry1%Base_Points_U(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=U_z

        Geometry1%Base_Points_V(1,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_x
        Geometry1%Base_Points_V(2,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_y
        Geometry1%Base_Points_V(3,(count-1)*n_order_sf+1:(count)*n_order_sf)=V_z

        Geometry1%w_smooth((count-1)*n_order_sf+1:(count)*n_order_sf)=w
    enddo
return
end






subroutine eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,n_order)
use some_types
implicit none

!List of calling arguments
integer ( kind = 8 ), intent(in) :: n_order
real ( kind = 8 ), intent(in) :: P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
real ( kind = 8 ), intent(in) :: U(n_order),V(n_order)
real ( kind = 8 ), intent(out) :: F_x(n_order),F_y(n_order),F_z(n_order),dS(n_order)
real ( kind = 8 ), intent(out) :: nP_x(n_order),nP_y(n_order),nP_z(n_order)

!List of local variables
real ( kind = 8 ) coef_x(6),coef_y(6),coef_z(6),U_x,U_y,U_z,V_x,V_y,V_z
integer ( kind = 8 ) count
    coef_x(1)=P1(1)
    coef_x(2)=-3*P1(1)-P2(1)+4*P4(1)
    coef_x(3)=-3*P1(1)-P3(1)+4*P6(1)
    coef_x(4)=2*P1(1)+2*P2(1)-4*P4(1)
    coef_x(5)=2*P1(1)+2*P3(1)-4*P6(1)
    coef_x(6)=4*P1(1)-4*P4(1)+4*P5(1)-4*P6(1)

    coef_y(1)=P1(2)
    coef_y(2)=-3*P1(2)-P2(2)+4*P4(2)
    coef_y(3)=-3*P1(2)-P3(2)+4*P6(2)
    coef_y(4)=2*P1(2)+2*P2(2)-4*P4(2)
    coef_y(5)=2*P1(2)+2*P3(2)-4*P6(2)
    coef_y(6)=4*P1(2)-4*P4(2)+4*P5(2)-4*P6(2)

    coef_z(1)=P1(3)
    coef_z(2)=-3*P1(3)-P2(3)+4*P4(3)
    coef_z(3)=-3*P1(3)-P3(3)+4*P6(3)
    coef_z(4)=2*P1(3)+2*P2(3)-4*P4(3)
    coef_z(5)=2*P1(3)+2*P3(3)-4*P6(3)
    coef_z(6)=4*P1(3)-4*P4(3)+4*P5(3)-4*P6(3)

    do count=1,n_order
F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
        U_x=coef_x(2)+2*coef_x(4)*U(count)+coef_x(6)*V(count)
        U_y=coef_y(2)+2*coef_y(4)*U(count)+coef_y(6)*V(count)
        U_z=coef_z(2)+2*coef_z(4)*U(count)+coef_z(6)*V(count)
        V_x=coef_x(3)+2*coef_x(5)*V(count)+coef_x(6)*U(count)
        V_y=coef_y(3)+2*coef_y(5)*V(count)+coef_y(6)*U(count)
        V_z=coef_z(3)+2*coef_z(5)*V(count)+coef_z(6)*U(count)
        nP_x(count)=U_y*V_z-U_z*V_y;
        nP_y(count)=U_z*V_x-U_x*V_z;
        nP_z(count)=U_x*V_y-U_y*V_x;
        dS(count)=sqrt(nP_x(count)**2+nP_y(count)**2+nP_z(count)**2);
        nP_x(count)=nP_x(count)/dS(count)
        nP_y(count)=nP_y(count)/dS(count)
        nP_z(count)=nP_z(count)/dS(count)
    enddo
return
end



subroutine eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order)
use some_types
implicit none

!List of calling arguments
integer ( kind = 8 ), intent(in) :: n_order
real ( kind = 8 ), intent(in) :: P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
real ( kind = 8 ), intent(in) :: U(n_order),V(n_order)
real ( kind = 8 ), intent(out) :: F_x(n_order),F_y(n_order),F_z(n_order),dS(n_order)
real ( kind = 8 ), intent(out) :: nP_x(n_order),nP_y(n_order),nP_z(n_order)
real ( kind = 8 ), intent(out) :: U_x(n_order),U_y(n_order),U_z(n_order),V_x(n_order),V_y(n_order),V_z(n_order)

!List of local variables
real ( kind = 8 ) coef_x(6),coef_y(6),coef_z(6)
integer ( kind = 8 ) count
    coef_x(1)=P1(1)
    coef_x(2)=-3*P1(1)-P2(1)+4*P4(1)
    coef_x(3)=-3*P1(1)-P3(1)+4*P6(1)
    coef_x(4)=2*P1(1)+2*P2(1)-4*P4(1)
    coef_x(5)=2*P1(1)+2*P3(1)-4*P6(1)
    coef_x(6)=4*P1(1)-4*P4(1)+4*P5(1)-4*P6(1)

    coef_y(1)=P1(2)
    coef_y(2)=-3*P1(2)-P2(2)+4*P4(2)
    coef_y(3)=-3*P1(2)-P3(2)+4*P6(2)
    coef_y(4)=2*P1(2)+2*P2(2)-4*P4(2)
    coef_y(5)=2*P1(2)+2*P3(2)-4*P6(2)
    coef_y(6)=4*P1(2)-4*P4(2)+4*P5(2)-4*P6(2)

    coef_z(1)=P1(3)
    coef_z(2)=-3*P1(3)-P2(3)+4*P4(3)
    coef_z(3)=-3*P1(3)-P3(3)+4*P6(3)
    coef_z(4)=2*P1(3)+2*P2(3)-4*P4(3)
    coef_z(5)=2*P1(3)+2*P3(3)-4*P6(3)
    coef_z(6)=4*P1(3)-4*P4(3)+4*P5(3)-4*P6(3)

    do count=1,n_order
F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
        U_x(count)=coef_x(2)+2*coef_x(4)*U(count)+coef_x(6)*V(count)
        U_y(count)=coef_y(2)+2*coef_y(4)*U(count)+coef_y(6)*V(count)
        U_z(count)=coef_z(2)+2*coef_z(4)*U(count)+coef_z(6)*V(count)
        V_x(count)=coef_x(3)+2*coef_x(5)*V(count)+coef_x(6)*U(count)
        V_y(count)=coef_y(3)+2*coef_y(5)*V(count)+coef_y(6)*U(count)
        V_z(count)=coef_z(3)+2*coef_z(5)*V(count)+coef_z(6)*U(count)
        nP_x(count)=U_y(count)*V_z(count)-U_z(count)*V_y(count);
        nP_y(count)=U_z(count)*V_x(count)-U_x(count)*V_z(count);
        nP_z(count)=U_x(count)*V_y(count)-U_y(count)*V_x(count);
        dS(count)=sqrt(nP_x(count)**2+nP_y(count)**2+nP_z(count)**2);
        nP_x(count)=nP_x(count)/dS(count)
        nP_y(count)=nP_y(count)/dS(count)
        nP_z(count)=nP_z(count)/dS(count)
    enddo
return
end

subroutine find_smooth_surface(Geometry1,alpha,Main_box)
use some_types
implicit none

interface
    subroutine My_Newton(x,tol,maxiter,Geometry1,alpha,flag,Main_box) !!Hay que modificar esto y vectorizarlo
    use some_types
        !List of calling arguments
        type (Geometry), intent(in) :: Geometry1
        integer ( kind = 8 ), intent(in) :: maxiter
        real ( kind = 8 ), intent(inout) :: x(Geometry1%n_Sf_points)
        real ( kind = 8 ), intent(in) :: tol,alpha
        integer ( kind = 8 ), intent(out) :: flag
        type (Box), pointer :: Main_box
    end subroutine
end interface

interface
    subroutine eval_density_grad_FMM(Geometry1,targets,alpha,F,grad_F,Main_box,max_step)
    use some_types
        !List of calling arguments
        type (Geometry), intent(in) :: Geometry1
        real ( kind = 8 ) , intent(in) :: targets(3,Geometry1%n_Sf_points),alpha
        real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), grad_F(3,Geometry1%n_Sf_points)
        real ( kind = 8 ), intent(out) ::  max_step(Geometry1%n_Sf_points)
        type (Box), pointer :: Main_box
    end subroutine
end interface


!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
real ( kind = 8 ), intent(in) :: alpha
type (Box), pointer :: Main_box

!List of local variables
integer ( kind = 8 ) flag,maxiter,ipoint,count,n_order_sf,itri
real ( kind = 8 ) h(Geometry1%n_Sf_points),tol,Norm_N
real ( kind = 8) r_t(3,Geometry1%n_Sf_points),F(Geometry1%n_Sf_points),grad_F(3,Geometry1%n_Sf_points)
real ( kind = 8 ) dVdu(3),dVdv(3),dhdu,dhdv,h_var
real ( kind = 8 ) max_step(Geometry1%n_Sf_points)
    do count=1,Geometry1%n_Sf_points
        h(count)=0.0d0
    enddo
    tol=1.0d-10
    maxiter=40
    flag=0
    write (*,*) 'No llega nunca1'
    call My_Newton(h,tol,maxiter,Geometry1,alpha,flag,Main_box)

    if (flag==1) then
        write (*,*) 'ERROR DE CONVERGENCIA NEWTON'
    end if
    do count=1,Geometry1%n_Sf_points
        Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+Geometry1%Base_Points_N(2,count)&
        &**2+Geometry1%Base_Points_N(3,count)**2)
r_t(:,count)=Geometry1%Base_Points(:,count)+Geometry1%Base_Points_N(:,count)*h(count)/Norm_N
        Geometry1%S_smooth(:,count)=r_t(:,count)
    enddo

    call eval_density_grad_FMM(Geometry1,r_t,alpha,F,grad_F,Main_box,max_step)

!    call eval_density_grad(Geometry1,r_t(1,:),r_t(2,:),r_t(3,:),alpha,F,grad_F)
    n_order_sf=Geometry1%n_Sf_points/Geometry1%ntri
    do count=1,Geometry1%n_Sf_points
        Geometry1%N_smooth(:,count)=-1.0d0*grad_F(:,count)/(sqrt(grad_F(1,count)**2+grad_F(2,count)**2+grad_F(3,count)**2))
        itri=(count-1)/n_order_sf+1
        Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+Geometry1%Base_Points_N(2,count)&
        &**2+Geometry1%Base_Points_N(3,count)**2)
        h_var=h(count)/Norm_N
        dVdu=Geometry1%Normal_Vert(:,Geometry1%Tri(2,itri))-Geometry1%Normal_Vert(:,Geometry1%Tri(1,itri))
        dVdv=Geometry1%Normal_Vert(:,Geometry1%Tri(3,itri))-Geometry1%Normal_Vert(:,Geometry1%Tri(1,itri))
        dhdu=-1.0d0*dot_product(Geometry1%N_smooth(:,count),dVdu*h_var+Geometry1%Base_Points_U(:,count))&
        &/dot_product(Geometry1%N_smooth(:,count),Geometry1%Base_Points_N(:,count))
        dhdv=-1.0d0*dot_product(Geometry1%N_smooth(:,count),dVdv*h_var+Geometry1%Base_Points_V(:,count))&
        &/dot_product(Geometry1%N_smooth(:,count),Geometry1%Base_Points_N(:,count))
        Geometry1%ru_smooth(:,count)=dVdu*h_var+Geometry1%Base_Points_N(:,count)*dhdu+Geometry1%Base_Points_U(:,count)
        Geometry1%rv_smooth(:,count)=dVdv*h_var+Geometry1%Base_Points_N(:,count)*dhdv+Geometry1%Base_Points_V(:,count)
        Geometry1%w_smooth(count)=Geometry1%w_smooth(count)*norm2(my_cross(Geometry1%ru_smooth(:,count)&
        &,Geometry1%rv_smooth(:,count)))
        Geometry1%ru_smooth(:,count)=Geometry1%ru_smooth(:,count)/norm2(Geometry1%ru_smooth(:,count))
        Geometry1%rv_smooth(:,count)=my_cross(Geometry1%N_smooth(:,count),Geometry1%ru_smooth(:,count))
    enddo

return
end


subroutine My_Newton(x,tol,maxiter,Geometry1,alpha,flag,Main_box) !!Hay que modificar esto y vectorizarlo
use some_types
implicit none

interface
    subroutine fun_roots_derivative(h,Geometry1,alpha,F,dF,Main_box,max_step)
        use some_types
        !List of calling arguments
        type (Geometry), intent(in) :: Geometry1
        real ( kind = 8 ), intent(in) :: alpha,h(Geometry1%n_Sf_points)
        real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), dF(Geometry1%n_Sf_points)
        real ( kind = 8 ), intent(out) :: max_step(Geometry1%n_Sf_points)
        type (Box), pointer :: Main_box
    end subroutine
end interface


!List of calling arguments
type (Geometry), intent(in) :: Geometry1
integer ( kind = 8 ), intent(in) :: maxiter
real ( kind = 8 ), intent(inout) :: x(Geometry1%n_Sf_points)
real ( kind = 8 ), intent(in) :: tol,alpha
integer ( kind = 8 ), intent(out) :: flag
type (Box), pointer :: Main_box

!List of local variables
integer ( kind = 8 ) count,count2
real ( kind = 8 ) F(Geometry1%n_Sf_points),dF(Geometry1%n_Sf_points),err(Geometry1%n_Sf_points)
real ( kind = 8 ) max_step(Geometry1%n_Sf_points)


    do count2=1,Geometry1%n_Sf_points
        err(count2)=tol+1.0d0
    enddo
    count=1
    flag=0
    do while ((maxval(err)>tol).and.(count<maxiter))
        call fun_roots_derivative(x,Geometry1,alpha,F,dF,Main_box,max_step)
        count=count+1
        do count2=1,Geometry1%n_Sf_points
            err(count2)=F(count2)/dF(count2)
            write (*,*) max_step(count2)
            if (abs(err(count2))>max_step(count2)) then
                err(count2)=err(count2)/abs(err(count2))*max_step(count2)
            endif
            x(count2)=x(count2)-err(count2)
            err(count2)=abs(err(count2))
        enddo
!        do count2=1,Geometry1%n_Sf_points
!            write (*,*) 'Newton iteration all: ', err(count2)
!        enddo
    write (*,*) 'Newton iteration: ', count, maxval(err)
    end do
    if (maxval(err)>tol) then
        flag=1
!        write (*,*) 'contador máximo Newton', count, err
    endif
return
end

subroutine GaussTri78(U,V,w)
use some_types
implicit none

!List of calling arguments
real ( kind = 8 ), intent(out) :: U(78),V(78),w(78)

    U(1)=2.506716198593823d-01
    U(2)=2.454819973818980d-01
    U(3)=1.414775223467523d-01
    U(4)=6.246218145439392d-02
    U(5)=1.351673429820326d-01
    U(6)=5.762583474888089d-02
    U(7)=1.157285643227313d-02
    U(8)=5.727409232524566d-02
    U(9)=1.477135713273572d-01
    U(10)=1.086565104434900d-02
    U(11)=1.145731953841411d-02
    U(12)=6.314190530269786d-02
    U(13)=1.262977548009919d-02
    U(14)=5.355756368497966d-02
    U(15)=9.131247716681390d-03
    U(16)=1.165344533744900d-02
    U(17)=1.414775223467523d-01
    U(18)=1.351673429820326d-01
    U(19)=5.762583474888089d-02
    U(20)=1.157285643227313d-02
    U(21)=5.727409232524566d-02
    U(22)=1.086565104434900d-02
    U(23)=1.145731953841411d-02
    U(24)=6.314190530269786d-02
    U(25)=1.262977548009919d-02
    U(26)=9.131247716681390d-03
    U(27)=3.746641900703088d-01
    U(28)=2.454819973818980d-01
    U(29)=3.657625546357003d-01
    U(30)=4.687689092728030d-01
    U(31)=2.505575451033266d-01
    U(32)=3.746561894876218d-01
    U(33)=4.392371829061964d-01
    U(34)=2.500294235413509d-01
    U(35)=1.477135713273572d-01
    U(36)=3.205405197076403d-01
    U(37)=2.059297451605091d-01
    U(38)=1.370596870191102d-01
    U(39)=1.185917096025605d-01
    U(40)=5.355756368497955d-02
    U(41)=5.517190361060625d-02
    U(42)=1.165344533744889d-02
    U(43)=4.927599230175472d-01
    U(44)=6.142751119146407d-01
    U(45)=5.677179757634971d-01
    U(46)=5.491899606615303d-01
    U(47)=6.926964841334032d-01
    U(48)=6.685938292480105d-01
    U(49)=7.826129353010766d-01
    U(50)=7.997984076781918d-01
    U(51)=8.687785149173401d-01
    U(52)=9.356968486727122d-01
    U(53)=3.746641900703089d-01
    U(54)=5.090360052362038d-01
    U(55)=4.927599230175473d-01
    U(56)=4.687689092728031d-01
    U(57)=6.142751119146407d-01
    U(58)=5.677179757634973d-01
    U(59)=5.491899606615305d-01
    U(60)=6.926964841334033d-01
    U(61)=7.045728573452854d-01
    U(62)=6.685938292480106d-01
    U(63)=7.826129353010766d-01
    U(64)=7.997984076781918d-01
    U(65)=8.687785149173401d-01
    U(66)=8.928848726300407d-01
    U(67)=9.356968486727122d-01
    U(68)=9.766931093251019d-01
    U(69)=3.657625546357005d-01
    U(70)=2.505575451033268d-01
    U(71)=3.746561894876220d-01
    U(72)=4.392371829061966d-01
    U(73)=2.500294235413512d-01
    U(74)=3.205405197076406d-01
    U(75)=2.059297451605094d-01
    U(76)=1.370596870191105d-01
    U(77)=1.185917096025609d-01
    U(78)=5.517190361060659d-02

    V(1)=3.746641900703088d-01
    V(2)=5.090360052362038d-01
    V(3)=4.927599230175473d-01
    V(4)=4.687689092728030d-01
    V(5)=6.142751119146408d-01
    V(6)=5.677179757634971d-01
    V(7)=5.491899606615303d-01
    V(8)=6.926964841334033d-01
    V(9)=7.045728573452855d-01
    V(10)=6.685938292480105d-01
    V(11)=7.826129353010767d-01
    V(12)=7.997984076781919d-01
    V(13)=8.687785149173401d-01
    V(14)=8.928848726300405d-01
    V(15)=9.356968486727122d-01
    V(16)=9.766931093251019d-01
    V(17)=3.657625546357003d-01
    V(18)=2.505575451033266d-01
    V(19)=3.746561894876219d-01
    V(20)=4.392371829061965d-01
    V(21)=2.500294235413509d-01
    V(22)=3.205405197076404d-01
    V(23)=2.059297451605092d-01
    V(24)=1.370596870191103d-01
    V(25)=1.185917096025606d-01
    V(26)=5.517190361060631d-02
    V(27)=2.506716198593823d-01
    V(28)=2.454819973818981d-01
    V(29)=1.414775223467523d-01
    V(30)=6.246218145439392d-02
    V(31)=1.351673429820326d-01
    V(32)=5.762583474888089d-02
    V(33)=1.157285643227318d-02
    V(34)=5.727409232524572d-02
    V(35)=1.477135713273573d-01
    V(36)=1.086565104434900d-02
    V(37)=1.145731953841417d-02
    V(38)=6.314190530269792d-02
    V(39)=1.262977548009930d-02
    V(40)=5.355756368497977d-02
    V(41)=9.131247716681501d-03
    V(42)=1.165344533744911d-02
    V(43)=1.414775223467523d-01
    V(44)=1.351673429820325d-01
    V(45)=5.762583474888089d-02
    V(46)=1.157285643227313d-02
    V(47)=5.727409232524566d-02
    V(48)=1.086565104434895d-02
    V(49)=1.145731953841406d-02
    V(50)=6.314190530269781d-02
    V(51)=1.262977548009914d-02
    V(52)=9.131247716681334d-03
    V(53)=3.746641900703088d-01
    V(54)=2.454819973818980d-01
    V(55)=3.657625546357003d-01
    V(56)=4.687689092728029d-01
    V(57)=2.505575451033266d-01
    V(58)=3.746561894876217d-01
    V(59)=4.392371829061963d-01
    V(60)=2.500294235413509d-01
    V(61)=1.477135713273571d-01
    V(62)=3.205405197076402d-01
    V(63)=2.059297451605090d-01
    V(64)=1.370596870191101d-01
    V(65)=1.185917096025604d-01
    V(66)=5.355756368497949d-02
    V(67)=5.517190361060620d-02
    V(68)=1.165344533744878d-02
    V(69)=4.927599230175472d-01
    V(70)=6.142751119146407d-01
    V(71)=5.677179757634970d-01
    V(72)=5.491899606615301d-01
    V(73)=6.926964841334032d-01
    V(74)=6.685938292480105d-01
    V(75)=7.826129353010767d-01
    V(76)=7.997984076781918d-01
    V(77)=8.687785149173401d-01
    V(78)=9.356968486727122d-01


    w(1)=1.595446593388906d-02
    w(2)=1.479574329197294d-02
    w(3)=1.205412183409247d-02
    w(4)=5.389662484702490d-03
    w(5)=1.032579481921157d-02
    w(6)=7.683338871370425d-03
    w(7)=3.330907634877638d-03
    w(8)=7.790759996149548d-03
    w(9)=9.233537320239574d-03
    w(10)=3.386799957729980d-03
    w(11)=3.009102424412810d-03
    w(12)=6.931582989269429d-03
    w(13)=2.346864909958631d-03
    w(14)=4.059587091313051d-03
    w(15)=1.332341144974900d-03
    w(16)=8.504413804547386d-04
    w(17)=1.205412183409247d-02
    w(18)=1.032579481921157d-02
    w(19)=7.683338871370425d-03
    w(20)=3.330907634877638d-03
    w(21)=7.790759996149548d-03
    w(22)=3.386799957729980d-03
    w(23)=3.009102424412810d-03
    w(24)=6.931582989269429d-03
    w(25)=2.346864909958631d-03
    w(26)=1.332341144974900d-03
    w(27)=1.595446593388906d-02
    w(28)=1.479574329197294d-02
    w(29)=1.205412183409247d-02
    w(30)=5.389662484702490d-03
    w(31)=1.032579481921157d-02
    w(32)=7.683338871370425d-03
    w(33)=3.330907634877638d-03
    w(34)=7.790759996149548d-03
    w(35)=9.233537320239574d-03
    w(36)=3.386799957729980d-03
    w(37)=3.009102424412810d-03
    w(38)=6.931582989269429d-03
    w(39)=2.346864909958631d-03
    w(40)=4.059587091313051d-03
    w(41)=1.332341144974900d-03
    w(42)=8.504413804547386d-04
    w(43)=1.205412183409247d-02
    w(44)=1.032579481921157d-02
    w(45)=7.683338871370425d-03
    w(46)=3.330907634877638d-03
    w(47)=7.790759996149548d-03
    w(48)=3.386799957729980d-03
    w(49)=3.009102424412810d-03
    w(50)=6.931582989269429d-03
    w(51)=2.346864909958631d-03
    w(52)=1.332341144974900d-03
    w(53)=1.595446593388906d-02
    w(54)=1.479574329197294d-02
    w(55)=1.205412183409247d-02
    w(56)=5.389662484702490d-03
    w(57)=1.032579481921157d-02
    w(58)=7.683338871370425d-03
    w(59)=3.330907634877638d-03
    w(60)=7.790759996149548d-03
    w(61)=9.233537320239574d-03
    w(62)=3.386799957729980d-03
    w(63)=3.009102424412810d-03
    w(64)=6.931582989269429d-03
    w(65)=2.346864909958631d-03
    w(66)=4.059587091313051d-03
    w(67)=1.332341144974900d-03
    w(68)=8.504413804547386d-04
    w(69)=1.205412183409247d-02
    w(70)=1.032579481921157d-02
    w(71)=7.683338871370425d-03
    w(72)=3.330907634877638d-03
    w(73)=7.790759996149548d-03
    w(74)=3.386799957729980d-03
    w(75)=3.009102424412810d-03
    w(76)=6.931582989269429d-03
    w(77)=2.346864909958631d-03
    w(78)=1.332341144974900d-03

return
end



subroutine GaussTri45(U,V,w)
use some_types
implicit none

!List of calling arguments
real ( kind = 8 ), intent(out) :: U(45),V(45),w(45)

    U(1)= 2.208782312184409d-01
    U(2)= 2.279434110906403d-01
    U(3)= 9.768561674285928d-02
    U(4)= 1.993517255300858d-02
    U(5)= 9.242062436313472d-02
    U(6)= 1.826986667642783d-02
    U(7)= 1.843714987046580d-02
    U(8)= 9.924857696792372d-02
    U(9)= 1.942704386325927d-02
    U(10)= 1.366715331364010d-02
    U(11)= 9.768561674285928d-02
    U(12)= 9.242062436313472d-02
    U(13)= 1.826986667642783d-02
    U(14)= 1.843714987046580d-02
    U(15)= 1.942704386325927d-02
    U(16)= 2.208782312184409d-01
    U(17)= 3.860282944546798d-01
    U(18)= 3.691524881881350d-01
    U(19)= 4.900324137234956d-01
    U(20)= 2.213846197875926d-01
    U(21)= 3.483600913607328d-01
    U(22)= 1.999439759389136d-01
    U(23)= 9.924857696792364d-02
    U(24)= 8.094245831218494d-02
    U(25)= 1.366715331363999d-02
    U(26)= 5.331618950690056d-01
    U(27)= 6.861947558492726d-01
    U(28)= 6.333700419628391d-01
    U(29)= 7.816188741906205d-01
    U(30)= 8.996304978245557d-01
    U(31)= 5.582435375631182d-01
    U(32)= 3.860282944546799d-01
    U(33)= 5.331618950690056d-01
    U(34)= 4.900324137234958d-01
    U(35)= 6.861947558492725d-01
    U(36)= 6.333700419628392d-01
    U(37)= 7.816188741906205d-01
    U(38)= 8.015028460641525d-01
    U(39)= 8.996304978245555d-01
    U(40)= 9.726656933727196d-01
    U(41)= 3.691524881881352d-01
    U(42)= 2.213846197875928d-01
    U(43)= 3.483600913607331d-01
    U(44)= 1.999439759389139d-01
    U(45)= 8.094245831218522d-02

    V(1)= 5.582435375631182d-01
    V(2)= 3.860282944546798d-01
    V(3)= 5.331618950690056d-01
    V(4)= 4.900324137234957d-01
    V(5)= 6.861947558492726d-01
    V(6)= 6.333700419628392d-01
    V(7)= 7.816188741906205d-01
    V(8)= 8.015028460641525d-01
    V(9)= 8.996304978245557d-01
    V(10)= 9.726656933727196d-01
    V(11)= 3.691524881881351d-01
    V(12)= 2.213846197875926d-01
    V(13)= 3.483600913607329d-01
    V(14)= 1.999439759389137d-01
    V(15)= 8.094245831218494d-02
    V(16)= 2.208782312184409d-01
    V(17)= 2.279434110906403d-01
    V(18)= 9.768561674285930d-02
    V(19)= 1.993517255300864d-02
    V(20)= 9.242062436313481d-02
    V(21)= 1.826986667642788d-02
    V(22)= 1.843714987046591d-02
    V(23)= 9.924857696792375d-02
    V(24)= 1.942704386325939d-02
    V(25)= 1.366715331364021d-02
    V(26)= 9.768561674285925d-02
    V(27)= 9.242062436313464d-02
    V(28)= 1.826986667642788d-02
    V(29)= 1.843714987046574d-02
    V(30)= 1.942704386325922d-02
    V(31)= 2.208782312184408d-01
    V(32)= 3.860282944546798d-01
    V(33)= 3.691524881881350d-01
    V(34)= 4.900324137234956d-01
    V(35)= 2.213846197875925d-01
    V(36)= 3.483600913607328d-01
    V(37)= 1.999439759389134d-01
    V(38)= 9.924857696792355d-02
    V(39)= 8.094245831218488d-02
    V(40)= 1.366715331363988d-02
    V(41)= 5.331618950690056d-01
    V(42)= 6.861947558492726d-01
    V(43)= 6.333700419628390d-01
    V(44)= 7.816188741906205d-01
    V(45)= 8.996304978245557d-01

    w(1)= 2.268433084372136d-02
    w(2)= 2.596237134029915d-02
    w(3)= 1.712942075066845d-02
    w(4)= 6.888599436262718d-03
    w(5)= 1.391190101815005d-02
    w(6)= 7.041187557933284d-03
    w(7)= 6.511662181692097d-03
    w(8)= 1.113176401316034d-02
    w(9)= 4.726400776735904d-03
    w(10)= 1.358456462863544d-03
    w(11)= 1.712942075066845d-02
    w(12)= 1.391190101815005d-02
    w(13)= 7.041187557933284d-03
    w(14)= 6.511662181692097d-03
    w(15)= 4.726400776735904d-03
    w(16)= 2.268433084372136d-02
    w(17)= 2.596237134029915d-02
    w(18)= 1.712942075066845d-02
    w(19)= 6.888599436262718d-03
    w(20)= 1.391190101815005d-02
    w(21)= 7.041187557933284d-03
    w(22)= 6.511662181692097d-03
    w(23)= 1.113176401316034d-02
    w(24)= 4.726400776735904d-03
    w(25)= 1.358456462863544d-03
    w(26)= 1.712942075066845d-02
    w(27)= 1.391190101815005d-02
    w(28)= 7.041187557933284d-03
    w(29)= 6.511662181692097d-03
    w(30)= 4.726400776735904d-03
    w(31)= 2.268433084372136d-02
    w(32)= 2.596237134029915d-02
    w(33)= 1.712942075066845d-02
    w(34)= 6.888599436262718d-03
    w(35)= 1.391190101815005d-02
    w(36)= 7.041187557933284d-03
    w(37)= 6.511662181692097d-03
    w(38)= 1.113176401316034d-02
    w(39)= 4.726400776735904d-03
    w(40)= 1.358456462863544d-03
    w(41)= 1.712942075066845d-02
    w(42)= 1.391190101815005d-02
    w(43)= 7.041187557933284d-03
    w(44)= 6.511662181692097d-03
    w(45)= 4.726400776735904d-03

return
end




subroutine record_Geometry(Geometry1,filename)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(in) :: Geometry1
character (len=30) filename

!List of local variables
integer ( kind = 8 ) umio,count1,count2,flag,n_order_sf
integer :: ierror
    open(8, FILE=filename,STATUS='REPLACE')
    n_order_sf=Geometry1%n_Sf_points/Geometry1%ntri
    write(8,*) n_order_sf
    write(8,*) Geometry1%ntri
    write(8,*) Geometry1%n_Sf_points
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
        write(8,*) Geometry1%ru_smooth(1,count1)
    enddo
    do count1=1,Geometry1%n_Sf_points
        write(8,*) Geometry1%ru_smooth(2,count1)
    enddo
    do count1=1,Geometry1%n_Sf_points
        write(8,*) Geometry1%ru_smooth(3,count1)
    enddo

    do count1=1,Geometry1%n_Sf_points
        write(8,*) Geometry1%rv_smooth(1,count1)
    enddo
    do count1=1,Geometry1%n_Sf_points
        write(8,*) Geometry1%rv_smooth(2,count1)
    enddo
    do count1=1,Geometry1%n_Sf_points
        write(8,*) Geometry1%rv_smooth(3,count1)
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

    do count1=1,Geometry1%n_Sf_points
        write(8,*) Geometry1%w_smooth(count1)
    enddo
    close (8)
return
end



subroutine check_Gauss(Geometry1,x0,y0,z0)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
real ( kind = 8 ), intent(in) :: x0,y0,z0
!List of local variables
integer ( kind = 8 ) umio,count1,count2,flag,n_order_sf
real ( kind = 8 )  F,Ex,Ey,Ez,R,x,y,z,pi,w,nx,ny,nz,err_rel
    pi=3.141592653589793238462643383d0
    F=0.0d0

    do count1=1,Geometry1%n_Sf_points
        x=Geometry1%S_smooth(1,count1)
        y=Geometry1%S_smooth(2,count1)
        z=Geometry1%S_smooth(3,count1)
        w=Geometry1%w_smooth(count1)
        nx=Geometry1%N_smooth(1,count1)
        ny=Geometry1%N_smooth(2,count1)
        nz=Geometry1%N_smooth(3,count1)
        R=sqrt((x-x0)**2+(y-y0)**2+(z-z0)**2)
        Ex=(x-x0)/(4*pi*R**3)
        Ey=(y-y0)/(4*pi*R**3)
        Ez=(z-z0)/(4*pi*R**3)
        F=F+(Ex*nx+Ey*ny+Ez*nz)*w
    enddo
    err_rel=abs(F-1)
    write (*,*) 'Relative Error: ',err_rel,'Value: ', F
return
end



subroutine refine_geometry(Geometry1,n_order_sf)
use some_types
implicit none

!List of calling arguments
type (Geometry), intent(inout) :: Geometry1
integer (kind = 8 ), intent(in) :: n_order_sf

!List of local variables
integer ( kind = 8 ) count,contador_indices
real ( kind = 8 ), allocatable :: Points(:,:), Normal_Vert(:,:)
integer ( kind = 8 ), allocatable :: Tri(:,:)
real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
real ( kind = 8 ) Pa(3),Pb(3),Pc(3),Pd(3),Pe(3),Pf(3),Pg(3),Ph(3),Pi(3)
real ( kind = 8 ) Nor1(3),Nor2(3),Nor3(3),Nor4(3),Nor5(3),Nor6(3)
real ( kind = 8 ) Nor_a(3),Nor_b(3),Nor_c(3),Nor_d(3),Nor_e(3),Nor_f(3),Nor_g(3),Nor_h(3),Nor_i(3)
real ( kind = 8 ) U(9),V(9)
real ( kind = 8 ) F_x(9),F_y(9),F_z(9),dS(9)
real ( kind = 8 ) nP_x(9),nP_y(9),nP_z(9)
real ( kind = 8 ) U_x(9),U_y(9),U_z(9),V_x(9),V_y(9),V_z(9)
integer ( kind = 8 ) m,N,n_order_aux

    allocate(Points(3,Geometry1%ntri*15))
    allocate(Normal_Vert(3,Geometry1%ntri*15))
    allocate(Tri(6,Geometry1%ntri*4))
    contador_indices=1
    n_order_aux=9
U=[0.2500d0,0.7500d0,0d0,0.2500d0,0.5000d0,0.7500d0,0.2500d0,0d0,0.2500d0]
V=[0d0,0d0,0.2500d0,0.2500d0,0.2500d0,0.2500d0,0.5000d0,0.7500d0,0.7500d0]
    do count=1,Geometry1%ntri
        P1=Geometry1%Points(:,Geometry1%Tri(1,count))
        P2=Geometry1%Points(:,Geometry1%Tri(2,count))
        P3=Geometry1%Points(:,Geometry1%Tri(3,count))
        P4=Geometry1%Points(:,Geometry1%Tri(4,count))
        P5=Geometry1%Points(:,Geometry1%Tri(5,count))
        P6=Geometry1%Points(:,Geometry1%Tri(6,count))
        call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order_aux)
        Pa=[F_x(1),F_y(1),F_z(1)]
        Pb=[F_x(2),F_y(2),F_z(2)]
        Pc=[F_x(3),F_y(3),F_z(3)]
        Pd=[F_x(4),F_y(4),F_z(4)]
        Pe=[F_x(5),F_y(5),F_z(5)]
        Pf=[F_x(6),F_y(6),F_z(6)]
        Pg=[F_x(7),F_y(7),F_z(7)]
        Ph=[F_x(8),F_y(8),F_z(8)]
        Pi=[F_x(9),F_y(9),F_z(9)]
        Nor1=Geometry1%Normal_Vert(:,Geometry1%Tri(1,count))
        Nor2=Geometry1%Normal_Vert(:,Geometry1%Tri(2,count))
        Nor3=Geometry1%Normal_Vert(:,Geometry1%Tri(3,count))
        Nor4=(Nor1+Nor2)/2.0d0
        Nor5=(Nor2+Nor3)/2.0d0
        Nor6=(Nor3+Nor1)/2.0d0
        Nor_a=(Nor1+Nor4)/2.0d0
        Nor_b=(Nor4+Nor2)/2.0d0
        Nor_c=(Nor1+Nor6)/2.0d0
        Nor_d=(Nor4+Nor6)/2.0d0
        Nor_e=(Nor4+Nor5)/2.0d0
        Nor_f=(Nor5+Nor2)/2.0d0
        Nor_g=(Nor5+Nor6)/2.0d0
        Nor_h=(Nor3+Nor6)/2.0d0
        Nor_i=(Nor3+Nor5)/2.0d0
        Points(:,contador_indices)=P1
        Points(:,contador_indices+1)=Pa
        Points(:,contador_indices+2)=P4
        Points(:,contador_indices+3)=Pb
        Points(:,contador_indices+4)=P2
        Points(:,contador_indices+5)=Pc
        Points(:,contador_indices+6)=Pd
        Points(:,contador_indices+7)=Pe
        Points(:,contador_indices+8)=Pf
        Points(:,contador_indices+9)=P6
        Points(:,contador_indices+10)=Pg
        Points(:,contador_indices+11)=P5
        Points(:,contador_indices+12)=Ph
        Points(:,contador_indices+13)=Pi
        Points(:,contador_indices+14)=P3

        Normal_Vert(:,contador_indices)=Nor1
        Normal_Vert(:,contador_indices+1)=Nor_a
        Normal_Vert(:,contador_indices+2)=Nor4
        Normal_Vert(:,contador_indices+3)=Nor_b
        Normal_Vert(:,contador_indices+4)=Nor2
        Normal_Vert(:,contador_indices+5)=Nor_c
        Normal_Vert(:,contador_indices+6)=Nor_d
        Normal_Vert(:,contador_indices+7)=Nor_e
        Normal_Vert(:,contador_indices+8)=Nor_f
        Normal_Vert(:,contador_indices+9)=Nor6
        Normal_Vert(:,contador_indices+10)=Nor_g
        Normal_Vert(:,contador_indices+11)=Nor5
        Normal_Vert(:,contador_indices+12)=Nor_h
        Normal_Vert(:,contador_indices+13)=Nor_i
        Normal_Vert(:,contador_indices+14)=Nor3
        Tri(:,(count-1)*4+1)=contador_indices*[1, 1, 1, 1, 1, 1]+[0, 2, 9, 1, 6, 5]
        Tri(:,(count-1)*4+2)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 11, 9, 7, 10, 6]
        Tri(:,(count-1)*4+3)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 4, 11, 3, 8, 7]
        Tri(:,(count-1)*4+4)=contador_indices*[1, 1, 1, 1, 1, 1]+[9, 11, 14, 10, 13, 12]
        contador_indices=contador_indices+15
    enddo
    Geometry1%npoints=Geometry1%ntri*15
    Geometry1%ntri=Geometry1%ntri*4
    m=Geometry1%npoints
    N=Geometry1%ntri
    Geometry1%n_Sf_points=N*n_order_sf

    deallocate(Geometry1%S_smooth)
    deallocate(Geometry1%N_smooth)
    deallocate(Geometry1%ru_smooth)
    deallocate(Geometry1%rv_smooth)
    deallocate(Geometry1%Base_Points)
    deallocate(Geometry1%Base_Points_N)
    deallocate(Geometry1%Base_Points_U)
    deallocate(Geometry1%Base_Points_V)
    deallocate(Geometry1%w_smooth)
    deallocate(Geometry1%Points)
    deallocate(Geometry1%Tri)
    deallocate(Geometry1%Normal_Vert)

    allocate(Geometry1%S_smooth(3,N*n_order_sf))
    allocate(Geometry1%N_smooth(3,N*n_order_sf))
    allocate(Geometry1%ru_smooth(3,N*n_order_sf))
    allocate(Geometry1%rv_smooth(3,N*n_order_sf))
    allocate(Geometry1%Base_Points(3,N*n_order_sf))
    allocate(Geometry1%Base_Points_N(3,N*n_order_sf))
    allocate(Geometry1%Base_Points_U(3,N*n_order_sf))
    allocate(Geometry1%Base_Points_V(3,N*n_order_sf))
    allocate(Geometry1%w_smooth(N*n_order_sf))

    allocate(Geometry1%Points(3,m))
    allocate(Geometry1%Tri(6,N))
    allocate(Geometry1%Normal_Vert(3,m))

    Geometry1%Points=Points
    Geometry1%Tri=Tri
    Geometry1%Normal_Vert=Normal_Vert

!    do count=1,Geometry1%npoints
!        write (*,*) 'puntos: ', Geometry1%Points(:,count)
!    enddo
!    read (*,*)
return
end

