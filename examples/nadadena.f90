Module some_types


  implicit none



  type Box_array
     type (Box), pointer :: BP => null ()
  end type Box_array





  
  type Box
     real ( kind = 8 ) Box_center(3)
     real ( kind = 8 ), allocatable :: Points(:,:)
     real ( kind = 8 ), allocatable :: sgmas(:)
     integer ( kind = 8 ) n_points,current_point
     type (Box_array) Children(8)
     type (Box_array) Colleague(27)
     type (Box), pointer :: Parent => null ()
  end type Box





type Geometry
    real ( kind = 8 ), allocatable :: S_smooth(:,:)             !Points on the real smooth surface
    real ( kind = 8 ), allocatable :: N_smooth(:,:)             !Normals on the real smooth surface
    real ( kind = 8 ), allocatable :: skeleton_Points(:,:)      !Integration nodes on the skeleton to compute F (F=0.5 is the surface)
    real ( kind = 8 ), allocatable :: skeleton_w(:)             !Integration weights on the skeleton to compute F (F=0.5 is the surface)
    real ( kind = 8 ), allocatable :: skeleton_N(:,:)           !Normal Vector on the skeleton at the nodes (to compute the double layer)
    real ( kind = 8 ), allocatable :: ru_smooth(:,:)            !u vectro on the real smooth surface
    real ( kind = 8 ), allocatable :: rv_smooth(:,:)            !v vectro on the real smooth surface
    real ( kind = 8 ), allocatable :: w_smooth(:)               !Integration weigths on the real smooth surface
    real ( kind = 8 ), allocatable :: Centroids(:,:)            !Centroids of each triangle of the skeleton (to compute sgma(x))
    real ( kind = 8 ), allocatable :: sgmas(:)                  !Values of sgma on each centroid, proportional to the side of the triangle
    real ( kind = 8 ), allocatable :: Points(:,:)               !Points that define the msh file (Each quadratic triangle has 6 points)
    real ( kind = 8 ), allocatable :: Normal_Vert(:,:)          !Pseudo-normals defined on each vertex of each triangle of the skeleton
    real ( kind = 8 ), allocatable :: Base_Points(:,:)          !base base points of eachs mooth point on the skeleton
    real ( kind = 8 ), allocatable :: Base_Points_N(:,:)        !Pseudo-normals defined on each base point of the surface
    real ( kind = 8 ), allocatable :: Base_Points_U(:,:)        !U vector defined on each base point of the smooth surface
    real ( kind = 8 ), allocatable :: Base_Points_V(:,:)        !U vector defined on each base point of the smooth surface
    integer ( kind = 8 ), allocatable :: Tri(:,:)               !Triangles of the skeleton (each triangle 6 points)
    integer ( kind = 8 ) npoints                                !Total number of points in the skeleton
    integer ( kind = 8 ) n_Sf_points                            !total number of points on the real smooth surface
    integer ( kind = 8 ) n_Sk_points                            !total number of integration nodes on the skeleton
    integer ( kind = 8 ) ntri                                   !Total number of triangles on the smooth surface
    integer ( kind = 8 ) ntri_sk                                   !Total number of triangles on the skeleton

end type Geometry



type Variable_Matrix
    integer ( kind = 8 ) n_Mat, current_n_Mat
    real ( kind = 8 ), allocatable :: Mat(:,:)
end type Variable_Matrix
type My_cell
    type (Variable_Matrix), allocatable :: Var_Mat(:)
    integer ( kind = 8 ) n_Cell
end type My_cell
contains
FUNCTION my_cross(a, b)
real ( kind = 8 ) my_cross(3)
real ( kind = 8 ), intent(in) :: a(3),b(3)
my_cross(1) = a(2) * b(3) - a(3) * b(2)
my_cross(2) = a(3) * b(1) - a(1) * b(3)
my_cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION my_cross
end module some_types


recursive subroutine allocate_tree_points(Current_box)
use some_types
implicit none

!List of calling arguments
type (Box), pointer :: Current_box

    if (Current_box%n_points>0) then
        allocate(Current_box%Points(3,Current_box%n_points))
        allocate(Current_box%sgmas(Current_box%n_points))
        Current_box%current_point=1
    endif
    if (associated(Current_box%Children(1)%BP)) then
        call allocate_tree_points(Current_box%Children(1)%BP)
    endif
    if (associated(Current_box%Children(2)%BP)) then
        call allocate_tree_points(Current_box%Children(2)%BP)
    endif
    if (associated(Current_box%Children(3)%BP)) then
        call allocate_tree_points(Current_box%Children(3)%BP)
    endif
    if (associated(Current_box%Children(4)%BP)) then
        call allocate_tree_points(Current_box%Children(4)%BP)
    endif
    if (associated(Current_box%Children(5)%BP)) then
        call allocate_tree_points(Current_box%Children(5)%BP)
    endif
    if (associated(Current_box%Children(6)%BP)) then
        call allocate_tree_points(Current_box%Children(6)%BP)
    endif
    if (associated(Current_box%Children(7)%BP)) then
        call allocate_tree_points(Current_box%Children(7)%BP)
    endif
    if (associated(Current_box%Children(8)%BP)) then
        call allocate_tree_points(Current_box%Children(8)%BP)
    endif
return
end


recursive subroutine Setup_tree_add_all(Current_box,Pts,N,radius,max_depth,current_depth,Box_limits)
use some_types
implicit none

interface
    subroutine Setup_tree_add_point(Current_box,Pt,radius,max_depth,current_depth,Box_limits)
        use some_types
        type (Box), pointer :: Current_box
        real ( kind = 8 ), intent(in) :: Pt(3),Box_limits(3,2),radius
        integer ( kind = 8 ), intent(in) :: max_depth,current_depth
    end subroutine
end interface


!List of calling arguments
integer ( kind = 8 ), intent(in) :: N
type (Box), pointer :: Current_box
real ( kind = 8 ), intent(in) :: Pts(3,N),Box_limits(3,2),radius
integer ( kind = 8 ), intent(in) :: max_depth,current_depth

!List of local variables
integer (kind = 8 ) count1

    do count1=1,N
call Setup_tree_add_point(Current_box,Pts(:,count1),radius,max_depth,current_depth,Box_limits)
    enddo
return
end


recursive subroutine Setup_tree_add_point(Current_box,Pt,radius,max_depth,current_depth,Box_limits)
use some_types
implicit none

!List of calling arguments
type (Box), pointer :: Current_box
real ( kind = 8 ), intent(in) :: Pt(3),Box_limits(3,2),radius
integer ( kind = 8 ), intent(in) :: max_depth,current_depth

!List of local variables

real ( kind = 8 ) x_center,y_center,z_center
real ( kind = 8 ) Box_limits_children(3,2)
real (kind = 8 ) dmax

    dmax=(Box_limits(1,2)-Box_limits(1,1))

    x_center=(Box_limits(1,2)+Box_limits(1,1))/2.0d0
    y_center=(Box_limits(2,2)+Box_limits(2,1))/2.0d0
    z_center=(Box_limits(3,2)+Box_limits(3,1))/2.0d0
    Current_box%Box_center(1)=x_center
    Current_box%Box_center(2)=y_center
    Current_box%Box_center(3)=z_center
    if ((max_depth>current_depth).and.(dmax>radius)) then
        if (Pt(1)<=x_center .and. Pt(2)<=y_center .and. Pt(3)<=z_center) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
            if (associated(Current_box%Children(1)%BP)) then
                call Setup_tree_add_point(Current_box%Children(1)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(1)%BP)
                Current_box%Children(1)%BP%Parent=>Current_box
                Current_box%Children(1)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(1)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        else if (Pt(1)>x_center .and. Pt(2)<=y_center .and. Pt(3)<=z_center) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
            if (associated(Current_box%Children(2)%BP)) then
                call Setup_tree_add_point(Current_box%Children(2)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(2)%BP)
                Current_box%Children(2)%BP%Parent=>Current_box
                Current_box%Children(2)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(2)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        else if (Pt(1)<=x_center .and. Pt(2)>y_center .and. Pt(3)<=z_center) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
            if (associated(Current_box%Children(3)%BP)) then
                call Setup_tree_add_point(Current_box%Children(3)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(3)%BP)
                Current_box%Children(3)%BP%Parent=>Current_box
                Current_box%Children(3)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(3)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        else if (Pt(1)>x_center .and. Pt(2)>y_center .and. Pt(3)<=z_center) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
            if (associated(Current_box%Children(4)%BP)) then
                call Setup_tree_add_point(Current_box%Children(4)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(4)%BP)
                Current_box%Children(4)%BP%Parent=>Current_box
                Current_box%Children(4)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(4)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        else if (Pt(1)<=x_center .and. Pt(2)<=y_center .and. Pt(3)>z_center) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
            if (associated(Current_box%Children(5)%BP)) then
                call Setup_tree_add_point(Current_box%Children(5)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(5)%BP)
                Current_box%Children(5)%BP%Parent=>Current_box
                Current_box%Children(5)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(5)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        else if (Pt(1)>x_center .and. Pt(2)<=y_center .and. Pt(3)>z_center) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
            if (associated(Current_box%Children(6)%BP)) then
                call Setup_tree_add_point(Current_box%Children(6)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(6)%BP)
                Current_box%Children(6)%BP%Parent=>Current_box
                Current_box%Children(6)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(6)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        else if (Pt(1)<=x_center .and. Pt(2)>y_center .and. Pt(3)>z_center) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
            if (associated(Current_box%Children(7)%BP)) then
                call Setup_tree_add_point(Current_box%Children(7)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(7)%BP)
                Current_box%Children(7)%BP%Parent=>Current_box
                Current_box%Children(7)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(7)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif

        else if (Pt(1)>x_center .and. Pt(2)>y_center .and. Pt(3)>z_center) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
            if (associated(Current_box%Children(8)%BP)) then
                call Setup_tree_add_point(Current_box%Children(8)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            else
                allocate(Current_box%Children(8)%BP)
                Current_box%Children(8)%BP%Parent=>Current_box
                Current_box%Children(8)%BP%n_points=0
                call Setup_tree_add_point(Current_box%Children(8)%BP,Pt,radius,max_depth,current_depth+1,Box_limits_children)
            endif
        endif
    else
        Current_box%n_points=Current_box%n_points+1
!        print*, 'N particles: ',Current_box%n_points
    endif


return
end




subroutine Locate_all_on_tree(Current_box,Pts,N,sgma_v,radius,Box_limits)
use some_types
implicit none

interface
    subroutine Locate_points_on_tree(Current_box,Pt,sgma,radius,Box_limits)
        use some_types
        type (Box), pointer :: Current_box
        real ( kind = 8 ), intent(in) :: Pt(3),Box_limits(3,2),sgma,radius
    end subroutine
end interface


!List of calling arguments
integer ( kind = 8 ), intent(in) :: N
type (Box), pointer :: Current_box
real ( kind = 8 ), intent(in) :: Pts(3,N),Box_limits(3,2),sgma_v(N),radius

!List of local variables
integer ( kind = 8 ) count1

    do count1=1,N
        call Locate_points_on_tree(Current_box,Pts(:,count1),sgma_v(count1),radius,Box_limits)
    enddo

return
end




recursive subroutine Locate_points_on_tree(Current_box,Pt,sgma,radius,Box_limits)
use some_types
implicit none

!List of calling arguments
type (Box), pointer :: Current_box
real ( kind = 8 ), intent(in) :: Pt(3),Box_limits(3,2),sgma,radius

!List of local variables

real ( kind = 8 ) x_center,y_center,z_center,dmax
real ( kind = 8 ) Box_limits_children(3,2)


    dmax=(Box_limits(1,2)-Box_limits(1,1))

    x_center=(Box_limits(1,2)+Box_limits(1,1))/2.0d0
    y_center=(Box_limits(2,2)+Box_limits(2,1))/2.0d0
    z_center=(Box_limits(3,2)+Box_limits(3,1))/2.0d0
    if (dmax.le.radius) then
        Current_box%Points(1,Current_box%current_point)=Pt(1)
        Current_box%Points(2,Current_box%current_point)=Pt(2)
        Current_box%Points(3,Current_box%current_point)=Pt(3)
        Current_box%sgmas(Current_box%current_point)=sgma
        Current_box%current_point=Current_box%current_point+1
    else
        if (Pt(1)<=x_center .and. Pt(2)<=y_center .and. Pt(3)<=z_center) then
            if (associated(Current_box%Children(1)%BP)) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
                call Locate_points_on_tree(Current_box%Children(1)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        elseif (Pt(1)>x_center .and. Pt(2)<=y_center .and. Pt(3)<=z_center) then
            if (associated(Current_box%Children(2)%BP)) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
                call Locate_points_on_tree(Current_box%Children(2)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        elseif (Pt(1)<=x_center .and. Pt(2)>y_center .and. Pt(3)<=z_center) then
            if (associated(Current_box%Children(3)%BP)) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
                call Locate_points_on_tree(Current_box%Children(3)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        elseif (Pt(1)>x_center .and. Pt(2)>y_center .and. Pt(3)<=z_center) then
            if (associated(Current_box%Children(4)%BP)) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=Box_limits(3,1)
                Box_limits_children(3,2)=z_center
                call Locate_points_on_tree(Current_box%Children(4)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        elseif (Pt(1)<=x_center .and. Pt(2)<=y_center .and. Pt(3)>z_center) then
            if (associated(Current_box%Children(5)%BP)) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
                call Locate_points_on_tree(Current_box%Children(5)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        elseif (Pt(1)>x_center .and. Pt(2)<=y_center .and. Pt(3)>z_center) then
            if (associated(Current_box%Children(6)%BP)) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=Box_limits(2,1)
                Box_limits_children(2,2)=y_center
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
                call Locate_points_on_tree(Current_box%Children(6)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        elseif (Pt(1)<=x_center .and. Pt(2)>y_center .and. Pt(3)>z_center) then
            if (associated(Current_box%Children(7)%BP)) then
                Box_limits_children(1,1)=Box_limits(1,1)
                Box_limits_children(1,2)=x_center
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
                call Locate_points_on_tree(Current_box%Children(7)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif

        elseif (Pt(1)>x_center .and. Pt(2)>y_center .and. Pt(3)>z_center) then
            if (associated(Current_box%Children(8)%BP)) then
                Box_limits_children(1,1)=x_center
                Box_limits_children(1,2)=Box_limits(1,2)
                Box_limits_children(2,1)=y_center
                Box_limits_children(2,2)=Box_limits(2,2)
                Box_limits_children(3,1)=z_center
                Box_limits_children(3,2)=Box_limits(3,2)
                call Locate_points_on_tree(Current_box%Children(8)%BP,Pt,sgma,radius,Box_limits_children)
            else
                Current_box%Points(1,Current_box%current_point)=Pt(1)
                Current_box%Points(2,Current_box%current_point)=Pt(2)
                Current_box%Points(3,Current_box%current_point)=Pt(3)
                Current_box%sgmas(Current_box%current_point)=sgma
                Current_box%current_point=Current_box%current_point+1
            endif
        endif

    endif
return
end


recursive subroutine Setup_Colleagues(Current_box,current_depth)
use some_types

interface
    subroutine are_colleagues(box1,box2,current_depth,answer)
        use some_types
        type (Box), pointer :: box1,box2
        integer ( kind = 8 ), intent(in) :: current_depth
        logical, intent(out) :: answer
    end subroutine

end interface


!List of calling arguments
type (Box), pointer :: Current_box
integer ( kind = 8 ), intent(in) :: current_depth

!List of local variables
integer (kind = 8 ) count1,count2,ipointer
logical answer

    ipointer=1
    if (associated(Current_box%Parent)) then
        do count1=1,27
            if (associated(Current_box%Parent%Colleague(count1)%BP)) then
                do count2=1,8
                  if (associated(Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP)) then
            call are_colleagues(Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP,Current_box,current_depth,answer)
                    if (answer) then
                        Current_box%Colleague(ipointer)%BP=>Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP
                        ipointer=ipointer+1
                    endif
                  endif
                enddo
            endif
        enddo
    else
        Current_box%Colleague(ipointer)%BP=>Current_box
    endif

    do count2=1,8
        if (associated(Current_box%Children(count2)%BP)) then
            call Setup_Colleagues(Current_box%Children(count2)%BP,current_depth+1)
        endif
    enddo

return
end

subroutine are_colleagues(box1,box2,current_depth,answer)
use some_types

!List of calling arguments
type (Box), pointer :: box1,box2
integer ( kind = 8 ), intent(in) :: current_depth
logical, intent(out) :: answer

!List of local variables
real (kind = 8 ) d1,d2,d3,length

    d1=abs(box1%Box_center(1)-box2%Box_center(1))
    d2=abs(box1%Box_center(2)-box2%Box_center(2))
    d3=abs(box1%Box_center(3)-box2%Box_center(3))
    length=1.0d0/(2**(current_depth))

    if (max(d1,d2,d3).le.(length+1.0d-14)) then
        answer=.true.
    else
        answer=.false.
    endif

return
end


subroutine brute_force_gaussian_all(Pts,sgma_v,alpha,N,targ_v,pot_v)
!List of calling arguments
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: Pts(3,N),sgma_v(N),targ_v(3,N),alpha
real ( kind = 8 ), intent(out) :: pot_v(N)

!List of local variables
integer (kind = 8 ) count

    do count=1,N
        call brute_force_gaussian(Pts,sgma_v,alpha,N,targ_v(:,count),pot_v(count))
    enddo


return
end



subroutine brute_force_gaussian(Pts,sgma_v,alpha,N,targ,pot)
!List of calling arguments
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: Pts(3,N),sgma_v(N),targ(3),alpha
real ( kind = 8 ), intent(out) :: pot

!List of local variables
integer (kind = 8 ) count
real ( kind = 8 ) d2, numerator, denominator

    pot=0.0d0
    numerator=0.0d0
    denominator=0.0d0
    do count=1,N
        d2=(Pts(1,count)-targ(1))**2+(Pts(2,count)-targ(2))**2+(Pts(3,count)-targ(3))**2
        numerator=numerator+sgma_v(count)*exp(-alpha*d2)
        denominator=denominator+exp(-alpha*d2)
    enddo
    pot=numerator/denominator
return
end

subroutine fast_gaussian_global(Current_box,targ_v,alpha,pot_v,N,Gx,Gy,Gz)
use some_types

interface
    subroutine fast_gaussian(Current_box,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
        use some_types
        type (Box), pointer :: Current_box
        real ( kind = 8 ), intent(in) :: targ(3),alpha
        real ( kind = 8 ), intent(inout) :: F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z
    end subroutine

end interface

!List of calling arguments
type (Box), pointer :: Current_box
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: targ_v(3,N),alpha
real ( kind = 8 ), intent(out) :: pot_v(N),Gx(N),Gy(N),Gz(N)

!List of local variables
integer ( kind = 8 ) count
real ( kind = 8 ) F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z
!!!!!$omp parallel do private (count1,numerator,denominator) shared(Current_box,targ_v,pot_v,alpha)
    do count=1,N
call fast_gaussian(Current_box,targ_v(:,count),alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
        pot_v(count)=F/D
        Gx(count)=(dF_x*D-F*dD_x)/D**2
        Gy(count)=(dF_y*D-F*dD_y)/D**2
        Gz(count)=(dF_z*D-F*dD_z)/D**2
    enddo
!!!!!$omp end parallel do
return
end


recursive subroutine fast_gaussian(Current_box,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
use some_types

interface
    subroutine fast_gaussian_box(Current_box,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
        use some_types
        type (Box), pointer :: Current_box
        real ( kind = 8 ), intent(in) :: targ(3),alpha
        real ( kind = 8 ), intent(inout) :: F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z
    end subroutine
end interface


!List of calling arguments
type (Box), pointer :: Current_box
real ( kind = 8 ), intent(in) :: targ(3),alpha
real ( kind = 8 ), intent(inout) :: F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z

!List of local variables
integer (kind = 8 ) count1

!write (*,*) 'control point',Current_box%Box_center(1),Current_box%Box_center(2),Current_box%Box_center(3)
!write (*,*) 'target',targ(1),targ(2),targ(3)

    if (.not. (associated(Current_box%Parent))) then
        F=0.0d0
        D=0.0d0
        dF_x=0.0d0
        dF_y=0.0d0
        dF_z=0.0d0
        dD_x=0.0d0
        dD_y=0.0d0
        dD_z=0.0d0
    endif
    do count1=1,27
        if (associated(Current_box%Colleague(count1)%BP)) then
            call fast_gaussian_box(Current_box%Colleague(count1)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
        endif
    enddo
    if (associated(Current_box%Children(1)%BP)) then
    if ((targ(1)<=Current_box%Box_center(1)).and.(targ(2)<=Current_box%Box_center(2)).and.(targ(3)<=Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(1)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(2)%BP)) then
    if ((targ(1)>Current_box%Box_center(1)).and.(targ(2)<=Current_box%Box_center(2)).and.(targ(3)<=Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(2)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(3)%BP)) then
    if ((targ(1)<=Current_box%Box_center(1)).and.(targ(2)>Current_box%Box_center(2)).and.(targ(3)<=Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(3)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(4)%BP)) then
    if ((targ(1)>Current_box%Box_center(1)).and.(targ(2)>Current_box%Box_center(2)).and.(targ(3)<=Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(4)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(5)%BP)) then
    if ((targ(1)<=Current_box%Box_center(1)).and.(targ(2)<=Current_box%Box_center(2)).and.(targ(3)>Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(5)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(6)%BP)) then
    if ((targ(1)>Current_box%Box_center(1)).and.(targ(2)<=Current_box%Box_center(2)).and.(targ(3)>Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(6)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(7)%BP)) then
    if ((targ(1)<=Current_box%Box_center(1)).and.(targ(2)>Current_box%Box_center(2)).and.(targ(3)>Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(7)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
    if (associated(Current_box%Children(8)%BP)) then
    if ((targ(1)>Current_box%Box_center(1)).and.(targ(2)>Current_box%Box_center(2)).and.(targ(3)>Current_box%Box_center(3))) then
            call fast_gaussian(Current_box%Children(8)%BP,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
    endif
    endif
return
end

subroutine fast_gaussian_box(Current_box,targ,alpha,F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z)
use some_types

!List of calling arguments
type (Box), pointer :: Current_box
real ( kind = 8 ), intent(in) :: targ(3),alpha
real ( kind = 8 ), intent(inout) :: F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z

!List of local variables
integer (kind = 8 ) count
real ( kind = 8 ) d2,my_exp
!    write (*,*) 'n_puntos_box: ',Current_box%n_points
!    read (*,*)
    do count=1,Current_box%n_points
!        write (*,*) 'control point',pot
d2=(Current_box%Points(1,count)-targ(1))**2+(Current_box%Points(2,count)-targ(2))**2+(Current_box%Points(3,count)-targ(3))**2
        my_exp=exp(-alpha*d2)
        F=F+Current_box%sgmas(count)*my_exp
        D=D+my_exp
        dF_x=dF_x+(-2.0d0*alpha)*Current_box%sgmas(count)*my_exp*(targ(1)-Current_box%Points(1,count))
        dF_y=dF_y+(-2.0d0*alpha)*Current_box%sgmas(count)*my_exp*(targ(2)-Current_box%Points(2,count))
        dF_z=dF_z+(-2.0d0*alpha)*Current_box%sgmas(count)*my_exp*(targ(3)-Current_box%Points(3,count))
        dD_x=dD_x+(-2.0d0*alpha)*my_exp*(targ(1)-Current_box%Points(1,count))
        dD_y=dD_y+(-2.0d0*alpha)*my_exp*(targ(2)-Current_box%Points(2,count))
        dD_z=dD_z+(-2.0d0*alpha)*my_exp*(targ(3)-Current_box%Points(3,count))
!        write (*,*) 'dentro_box:',F,D,Current_box%sgmas(count),my_exp
    enddo
return
end





subroutine brute_force_gaussian_all_derivatives(Pts,sgma_v,alpha,N,targ_v,pot_v,Gx,Gy,Gz)
implicit none
!List of calling arguments
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: Pts(3,N),sgma_v(N),targ_v(3,N),alpha
real ( kind = 8 ), intent(out) :: pot_v(N),Gx(N),Gy(N),Gz(N)

!List of local variables
integer (kind = 8 ) count

    do count=1,N
call brute_force_gaussian_derivatives(Pts,sgma_v,alpha,N,targ_v(:,count),pot_v(count),Gx(count),Gy(count),Gz(count))
    enddo

return
end



subroutine brute_force_gaussian_derivatives(Pts,sgma_v,alpha,N,targ,pot,Gx,Gy,Gz)
implicit none
!List of calling arguments
integer ( kind = 8 ), intent(in) :: N
real ( kind = 8 ), intent(in) :: Pts(3,N),sgma_v(N),targ(3),alpha
real ( kind = 8 ), intent(out) :: pot,Gx,Gy,Gz

!List of local variables
integer (kind = 8 ) count
real ( kind = 8 ) d2
real ( kind = 8 ) F,D,dF_x,dF_y,dF_z,dD_x,dD_y,dD_z,my_exp
    F=0.0d0
    D=0.0d0
    dF_x=0.0d0
    dF_y=0.0d0
    dF_z=0.0d0
    dD_x=0.0d0
    dD_y=0.0d0
    dD_z=0.0d0
    pot=0.0d0
    do count=1,N
        d2=(Pts(1,count)-targ(1))**2+(Pts(2,count)-targ(2))**2+(Pts(3,count)-targ(3))**2
        my_exp=exp(-alpha*d2)
        F=F+sgma_v(count)*my_exp
        dF_x=dF_x+(-2*alpha)*sgma_v(count)*my_exp*(targ(1)-Pts(1,count))
        dF_y=dF_y+(-2*alpha)*sgma_v(count)*my_exp*(targ(2)-Pts(2,count))
        dF_z=dF_z+(-2*alpha)*sgma_v(count)*my_exp*(targ(3)-Pts(3,count))
        dD_x=dD_x+(-2*alpha)*my_exp*(targ(1)-Pts(1,count))
        dD_y=dD_y+(-2*alpha)*my_exp*(targ(2)-Pts(2,count))
        dD_z=dD_z+(-2*alpha)*my_exp*(targ(3)-Pts(3,count))
        D=D+my_exp
    enddo
    Gx=(dF_x*D-F*dD_x)/D**2
    Gy=(dF_y*D-F*dD_y)/D**2
    Gz=(dF_z*D-F*dD_z)/D**2
    pot=F/D
!    write (*,*) 'num y denom brute force: ', F,D
!    call sleep(9)
return
end


