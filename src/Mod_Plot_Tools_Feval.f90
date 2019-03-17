Module Mod_Plot_Tools_Feval
use Mod_TreeLRD
use ModType_Smooth_Surface
use Mod_Fast_Sigma
use Mod_Feval
implicit none

public :: plot_Feval

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF SUBROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

subroutine plot_Feval(Feval_stuff_1, Geometry1)
implicit none

!List of calling arguments
type ( Feval_stuff ), pointer :: Feval_stuff_1
type ( Geometry ), pointer :: Geometry1


!List of local variables
type ( TreeLRD ), pointer :: TreeLRD_1
type ( Fast_Sigma_stuff ), pointer ::  FSS_1
character ( len=100 ) nombre,filename,plot_name
real ( kind = 8 ) x_min,x_max,y_min,y_max,z_min,z_max,Lx,Ly,Lz
real ( kind = 8 ) xp_min,xp_max,yp_min,yp_max,zp_min,zp_max
real ( kind = 8 ), allocatable :: F_plot(:,:),targ_vect(:,:),F(:),F_grad(:,:),sgma_y(:),sgma_z(:)
integer ( kind = 8 ) N_plot,M_plot,count,count1,count2,icount,adapt_flag,speed_flag,n_targ


    TreeLRD_1 => Feval_stuff_1%FSS_1%TreeLRD_1

!    do count=1,Geometry1%ntri_sk
!        write (*,*) TreeLRD_1%W_Pts_mem(1,count),TreeLRD_1%W_Pts_mem(2,count),TreeLRD_1%W_Pts_mem(3,count)
!        write (*,*) TreeLRD_1%W_sgmas_mem(count)
!    enddo
    plot_name='plot_points'
    call plot_curve_3D(Geometry1%skeleton_Points(1,:),Geometry1%skeleton_Points(2,:),&
     &Geometry1%skeleton_Points(3,:),Geometry1%n_Sk_points,plot_name)
    write (*,*) 'mucho ojo nen1'

!    call plot_curve_3D(Geometry1%skeleton_N(1,:),Geometry1%skeleton_N(2,:),&
!     &Geometry1%skeleton_N(3,:),Geometry1%n_Sk_points,plot_name)

    write (*,*) 'mucho ojo nen2'

!call plot_curve_3D(TreeLRD_1%W_Pts_mem(1,:),TreeLRD_1%W_Pts_mem(2,:),TreeLRD_1%W_Pts_mem(3,:),Geometry1%ntri_sk,plot_name)

    write(*,*) 'inside2'

    M_plot=100
    N_plot=100

    allocate(F_plot(M_plot,N_plot))
    allocate(targ_vect(3,M_plot*N_plot))
    allocate(F(M_plot*N_plot))
    allocate(F_grad(3,M_plot*N_plot))

    x_min=minval(TreeLRD_1%W_Pts_mem(1,:))
    x_max=maxval(TreeLRD_1%W_Pts_mem(1,:))
    y_min=minval(TreeLRD_1%W_Pts_mem(2,:))
    y_max=maxval(TreeLRD_1%W_Pts_mem(2,:))
    z_min=minval(TreeLRD_1%W_Pts_mem(3,:))
    z_max=maxval(TreeLRD_1%W_Pts_mem(3,:))
    Lx=x_max-x_min
    Ly=y_max-y_min
    Lz=z_max-z_min

    xp_max=x_max+Lx/4.0d0
    xp_min=x_min-Lx/4.0d0
    yp_max=y_max+Ly/4.0d0
    yp_min=y_min-Ly/4.0d0
    zp_max=z_max+Lz/4.0d0
    zp_min=z_min-Lz/4.0d0


    icount=1
    do count2=1,M_plot
        do count1=1,N_plot
!!! Plane ZY
            targ_vect(3,icount)=zp_min+(zp_max-zp_min)*(count1-1.0d0)/M_plot
            targ_vect(2,icount)=yp_min+(yp_max-yp_min)*(count2-1.0d0)/N_plot
            targ_vect(1,icount)=(x_max+x_min)/2.0d0+(x_max-x_min)/4.0d0*0.0d0
!            call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,icount),sgma(icount))
            icount=icount+1
        enddo
    enddo
    adapt_flag=1
    speed_flag=1
    n_targ=N_plot*M_plot
    write (*,*) 'inicio primer FMM'
    call eval_density_grad_FMM(Geometry1,targ_vect,n_targ,F,F_grad,Feval_stuff_1,adapt_flag)
    write (*,*) 'fin primer FMM'
!    call function_eval_sigma(FSS_1,targ_vect,n_targ,sgma,sgma_x,sgma_y,sgma_z,adapt_flag,speed_flag)
    plot_name='./plot_tools/plot_sigma_1'
    write (*,*) 'inicio primer plot'
    call plot2D_v3(targ_vect(1,:),targ_vect(2,:),targ_vect(3,:),F,N_plot,M_plot,plot_name)

    write (*,*) 'fin primer plot'

    icount=1
    do count2=1,M_plot
        do count1=1,N_plot
!!! Plane XY
            targ_vect(1,icount)=xp_min+(xp_max-xp_min)*(count1-1.0d0)/M_plot
            targ_vect(2,icount)=yp_min+(yp_max-yp_min)*(count2-1.0d0)/N_plot
            targ_vect(3,icount)=(z_max+z_min)/2.0d0+(z_max-z_min)/4.0d0*0.0d0
!            call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,icount),sgma(icount))
            icount=icount+1
        enddo
    enddo
    n_targ=N_plot*M_plot
    write (*,*) 'inicio segundo FMM'
    call eval_density_grad_FMM(Geometry1,targ_vect,n_targ,F,F_grad,Feval_stuff_1,adapt_flag)
    write (*,*) 'fin segundo FMM'

!    call function_eval_sigma(FSS_1,targ_vect,n_targ,sgma,sgma_x,sgma_y,sgma_z,adapt_flag,speed_flag)

    plot_name='./plot_tools/plot_sigma_2'
    call plot2D_v3(targ_vect(1,:),targ_vect(2,:),targ_vect(3,:),F,N_plot,M_plot,plot_name)

    icount=1
    do count2=1,M_plot
        do count1=1,N_plot
!!! Plane XZ
            targ_vect(1,icount)=xp_min+(xp_max-xp_min)*(count1-1.0d0)/M_plot
            targ_vect(3,icount)=zp_min+(zp_max-zp_min)*(count2-1.0d0)/N_plot
            targ_vect(2,icount)=(y_max+y_min)/2.0d0
!            call initial_guess_sgma(TreeLRD_1%Main_box,targ_vect(:,icount),sgma(icount))
            icount=icount+1
        enddo
    enddo
    n_targ=N_plot*M_plot

    call eval_density_grad_FMM(Geometry1,targ_vect,n_targ,F,F_grad,Feval_stuff_1,adapt_flag)
!    call function_eval_sigma(FSS_1,targ_vect,n_targ,sgma,sgma_x,sgma_y,sgma_z,adapt_flag,speed_flag)

    plot_name='./plot_tools/plot_sigma_3'
    call plot2D_v3(targ_vect(1,:),targ_vect(2,:),targ_vect(3,:),F,N_plot,M_plot,plot_name)

    deallocate(F_plot)
    deallocate(targ_vect)
    deallocate(F)
    deallocate(F_grad)


    call plot_tree_tool(TreeLRD_1)

return
end


subroutine plot_tree_tool(TreeLRD_1)
implicit none

!List of calling arguments
type ( TreeLRD ), pointer :: TreeLRD_1

!List of local variables
character ( len=100 ) plot_name
real ( kind = 8 ), allocatable :: W_boxes(:)
integer ( kind = 8 ) icount,tot_num_box

    plot_name='plot_tree'
    tot_num_box=total_number_leaf_boxes(TreeLRD_1%Main_box)
    allocate(W_boxes(tot_num_box*4))
    icount=1
    call pile_leaf_boxes(TreeLRD_1%Main_box,W_boxes,icount,tot_num_box)
    call plot_tree(W_boxes,tot_num_box,TreeLRD_1%W_Pts_mem,TreeLRD_1%ntot_W_Pts_mem,plot_name)
    deallocate(W_boxes)

return
end


recursive subroutine pile_boxes(Current_box,W_boxes,icount,tot_num_box)
implicit none

!List of calling arguments
type ( Box ), pointer :: Current_box
integer ( kind = 8 ), intent(in) :: tot_num_box!List of local variables
integer ( kind = 8 ), intent(inout) :: icount!List of local variables
real ( kind = 8 ), intent(inout) :: W_boxes(4*tot_num_box)

!List of local variables
integer ( kind = 8 ) count1

    if (Current_box%is_leaf) then
        W_boxes(icount)=Current_box%Box_center(1)
        W_boxes(icount+1)=Current_box%Box_center(2)
        W_boxes(icount+2)=Current_box%Box_center(3)
        W_boxes(icount+3)=Current_box%Box_size
        icount=icount+4
    else
        do count1=1,8
            call pile_boxes(Current_box%Children(count1)%BP,W_boxes,icount,tot_num_box)
        enddo
        W_boxes(icount)=Current_box%Box_center(1)
        W_boxes(icount+1)=Current_box%Box_center(2)
        W_boxes(icount+2)=Current_box%Box_center(3)
        W_boxes(icount+3)=Current_box%Box_size
        icount=icount+4
    endif

return
end


recursive subroutine pile_leaf_boxes(Current_box,W_boxes,icount,tot_num_box)
implicit none

!List of calling arguments
type ( Box ), pointer :: Current_box
integer ( kind = 8 ), intent(in) :: tot_num_box!List of local variables
integer ( kind = 8 ), intent(inout) :: icount!List of local variables
real ( kind = 8 ), intent(inout) :: W_boxes(4*tot_num_box)

!List of local variables
integer ( kind = 8 ) count1

    if (Current_box%is_leaf) then
        W_boxes(icount)=Current_box%Box_center(1)
        W_boxes(icount+1)=Current_box%Box_center(2)
        W_boxes(icount+2)=Current_box%Box_center(3)
        W_boxes(icount+3)=Current_box%Box_size
        icount=icount+4
    else
        do count1=1,8
            call pile_leaf_boxes(Current_box%Children(count1)%BP,W_boxes,icount,tot_num_box)
        enddo
    endif

return
end


end module Mod_Plot_Tools_Feval
