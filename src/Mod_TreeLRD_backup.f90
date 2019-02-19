!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                    MODULE MAIN INFORMATION                          !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! DESCRIPTION:
!!!
!!! This module contains functions to make a fully adaptive tree and to turn this tree
!!! into a level restricted tree and defragment memory
!!!
!!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!!
!!! Public:
!!!
!!!  1-start_tree                            !! generate a fully adaptive tree
!!!  2-make_level_restricted                 !! transform a fully adaptive tree into a level restricted tree
!!!  3-setup_nearest                         !! Set the pointers to the nearest neig. among leaf nodes
!!!  3-defrag_tree_Points                    !! this defragments the tree. Now Points and sgmas of each box are
!!!                                          !! located in a particular order in memory
!!!                                          !! that is, each box at each level will have its points
!!!                                          !! at contiguous locations in memory
!!!
!!! Private:
!!!
!!!  4-subdivide_box
!!!  5-setup_colleagues
!!!  6-is_touching
!!!  7-subdivide_box_once
!!!  8-flag_primary_violators
!!!  9-flag_primary_violators_pair
!!! 10-is_primary_violating
!!! 11-fix_primary_violators
!!! 12-flag_secondary_violators
!!! 13-flag_secondary_violators_pair
!!! 14-is_secondary_violating
!!! 15-fix_secondary_violators
!!! 16-review_secondary
!!!
!!! INDEX OF DERIVED TYPES:
!!!
!!! Public:
!!!
!!! 1-TreeLRD       !! Main type. Contains the full tree, including the main box of the tree, and all points and sigmas
!!! 2-Box           !! The tree is a 'linked list' in the form of boxes pointing to other boxes (children, parent, Colleagues)
!!!
!!! Private:
!!!
!!! 3-Box_array     !! Array of boxes, allows to group all the children or Colleagues of a box in one variable
!!!
!!! MODULES REQUIRED: (none)
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module Mod_TreeLRD
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

type, public :: TreeLRD
    type ( Box ), pointer :: Main_box => null ()        !! Pointer to the main box of the tree
    real ( kind = 8 ), pointer :: W_Pts_mem(:,:) => null () !! Pointer to a big array that contains all the points (centroids)
                                                            !! of the tree in the usual contiguous order
    real ( kind = 8 ), pointer :: W_sgmas_mem(:) => null () !! Pointer to a big array that contains all the values of sigma
                                                            !! associated to each centroid
    integer ( kind = 8 ) :: ntot_W_Pts_mem                  !! Total number of points (centroids) on the tree
    integer ( kind = 8 ) :: ntot_W_sgmas_mem                !! Total number of sigmas (same as centroids...)

end type TreeLRD

type, private :: Box_array
    type ( Box ), pointer :: BP => null ()
end type Box_array

type, public :: Box
    !! This type is a box, that contains information about the box and pointers to other
    !! boxes
    real ( kind = 8 ) Box_center(3)                       !! Center of the box
    real ( kind = 8 ) Box_size                            !! Size of the box (x_max-x_min)
    real ( kind = 8 ), pointer :: Points(:,:) => null ()  !! Array to points inside the box, size will be (3,n_points)
    real ( kind = 8 ), pointer :: sgmas(:) => null ()     !! Values of sigma associated to each point
    logical :: is_leaf                                    !! This flag indicates if the box is a leaf node
    logical :: primary_violator                           !! This flag is used to build a level restricted tree
    logical :: secondary_violator                         !! This flag is used to build a level restricted tree
    integer ( kind = 8 ) n_points                         !! This number indicates the number of points inside the box
    type ( Box_array ) Children(8)                        !! This concains pointers to child boxes
    type ( Box_array ) Colleague(27)                      !! This contains pointers to Colleague boxes
    type ( Box_array ) Nearest(56)                        !! This contains pointers to Colleague boxes
    type ( Box ), pointer :: Parent => null ()            !! This contains a pointer to the parent's box
end type Box


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    END DEFINITION OF TYPES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF SUBROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

public :: start_tree
public :: make_level_restricted
public :: defrag_tree_Points

private :: subdivide_box
private :: setup_colleagues
private :: is_touching
private :: subdivide_box_once
private :: flag_primary_violators
private :: flag_primary_violators_pair
private :: is_primary_violating
private :: fix_primary_violators
private :: flag_secondary_violators
private :: flag_secondary_violators_pair
private :: is_secondary_violating
private :: fix_secondary_violators
private :: review_secondary

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! These functions are for the purpose of making a fully adaptive tree !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!
!! 1-subroutine start_tree          !! This is the only function that the user will use
                                    !! example >> call start_tree(Main_box,n_max_leaf,Pts,sgmas,n)
!! 2-subroutine subdivide_box
!! 3-subroutine setup_colleagues
!! 4-function is_touching
!! 5-subroutine subdivide_box_once
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine start_tree(TreeLRD_1,n_max_leaf,Pts,sgmas,n)
implicit none
!! This subroutine starts a tree with n input points Pts and sgma values on each point

    !List of calling arguments
    type ( TreeLRD ), pointer :: TreeLRD_1          !! pointer to null where the tree will be allocated
    integer ( kind = 8 ), intent(in) :: n           !! number of points to allocate in the tree
    integer ( kind = 8 ), intent(in) :: n_max_leaf  !! maximum number of points allowed in a box to stop splitting
    real ( kind = 8 ), intent(in) :: Pts(3,n)       !! real matrix with the points that will be allocated in the tree
    real ( kind = 8 ), intent(in) :: sgmas(n)       !! real matrix with the values of sigma associated to each source

    !List of local variables
    real ( kind = 8 ) d_aux
    type ( Box ), pointer :: Main_box
    real ( kind = 8 ) d_all_dim(3),max_x,min_x,max_y,min_y,max_z,min_z

        allocate(TreeLRD_1)                         !! Allocates the tree
        allocate(TreeLRD_1%Main_box)                !! Allocates the main box and then initializes all the values
        Main_box => TreeLRD_1%Main_box

        max_x=maxval(Pts(1,:))
        min_x=minval(Pts(1,:))
        max_y=maxval(Pts(2,:))
        min_y=minval(Pts(2,:))
        max_z=maxval(Pts(3,:))
        min_z=minval(Pts(3,:))
        d_all_dim(1)=max_x-min_x
        d_all_dim(2)=max_y-min_y
        d_all_dim(3)=max_z-min_z
        Main_box%Box_size=maxval(d_all_dim)         !! Initizalyze the size of the box
                                                    !! Initialyze points and sgmas
        allocate(Main_box%Points(3,n))
        allocate(Main_box%sgmas(n))
        Main_box%Points=Pts
        Main_box%sgmas=sgmas
        Main_box%n_points=n
                                                    !! Initialyze flags related to Level Restricted Tree
        Main_box%is_leaf=(.true.)
        Main_box%primary_violator=(.false.)
        Main_box%secondary_violator=(.false.)
                                                    !! Initialyze the center of the box
        Main_box%Box_center(1)=(max_x+min_x)/2.0d0
        Main_box%Box_center(2)=(max_y+min_y)/2.0d0
        Main_box%Box_center(3)=(max_z+min_z)/2.0d0
        call subdivide_box(Main_box,n_max_leaf)     !! After initialyzing the values of the main box, splits recursively

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine subdivide_box(Current_box,n_max_leaf)
implicit none
!! This subroutine splits the current box in 8 until the number of points per
!! box is lower than n_max_leaf. The resulting tree is completely adaptive and with
!! no restrictions. Other functions have to be used afterwards to transform it in a Level
!! Restricted tree. The subroutine acts recursively until all the tree is built

    !List of calling arguments
    type ( Box ), pointer :: Current_box              !! pointer to the box that we want to split
    integer ( kind = 8 ), intent(in) :: n_max_leaf    !! maximum number of points allowed in a box to stop splitting

    !List of local variables
    integer ( kind = 8 ) count1

        if (Current_box%n_points>n_max_leaf) then
            call subdivide_box_once(Current_box)
            do count1=1,8
                call subdivide_box(Current_box%Children(count1)%BP,n_max_leaf)
            enddo
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine setup_colleagues(Current_box)
implicit none
!! This subroutine setup the colleagues in the subtree starting from Current Box (all decendents)

    !List of calling arguments
    type ( Box ), pointer :: Current_box          !! Pointer to the subtree that we want to set up

    !List of local variables
    integer ( kind = 8 ) count1,count2,ipointer,num_col
    logical answer

        ipointer=1
        !! Loop through all the childrens of the colleagues of my parent and check if they are touching
        if (associated(Current_box%Parent)) then
            do count1=1,27
                if (associated(Current_box%Parent%Colleague(count1)%BP)) then
                    do count2=1,8
                        if (associated(Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP)) then
                            if (is_touching(Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP,Current_box)) then
                                num_col=number_colleague(Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP,Current_box)
                                Current_box%Colleague(num_col)%BP=>Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP
                                Current_box%Parent%Colleague(count1)%BP%Children(count2)%BP%Colleague(28-num_col)%BP=>Current_box
                            endif
                        endif
                    enddo
                endif
            enddo
        else
            Current_box%Colleague(ipointer)%BP=>Current_box
        endif
        !! Do a recursive call to all my children
        do count2=1,8
            if (associated(Current_box%Children(count2)%BP)) then
                call setup_colleagues(Current_box%Children(count2)%BP)
            endif
        enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine setup_nearest(Current_box)
implicit none
!! This subroutine setup the nearest neighbours of the leaf boxes in the full tree. Procede by recursion from Main_Box

    !List of calling arguments
    type ( Box ), pointer :: Current_box          !! Pointer to the subtree that we want to set up

    !List of local variables
    integer ( kind = 8 ) count1

        if (Current_box%is_leaf) then
            call setup_nearest_leaf(Current_box)
        else
            do count1=1,8
                call setup_nearest(Current_box%Children(count1)%BP)
            enddo
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine setup_nearest_leaf(Current_box)
implicit none
!! This subroutine setup the nearest neighbours of the current box (which is, must be, a leaf box)

    !List of calling arguments
    type ( Box ), pointer :: Current_box          !! Pointer to the subtree that we want to set up

    !List of local variables
    integer ( kind = 8 ) count1,count2,ipointer
    type ( Box ), pointer :: aux_box
    logical answer

        ipointer=1
        do count1=1,27
            if (associated(Current_box%Colleague(count1)%BP)) then
                if (.not.(associated(Current_box,Current_box%Colleague(count1)%BP))) then
                    if (Current_box%Colleague(count1)%BP%is_leaf) then
                        Current_box%Nearest(ipointer)%BP => Current_box%Colleague(count1)%BP
                        ipointer=ipointer+1
                    else
                        do count2=1,8
                            if (is_touching(Current_box%Colleague(count1)%BP%Children(count2)%BP,Current_box)) then
                                Current_box%Nearest(ipointer)%BP => Current_box%Colleague(count1)%BP%Children(count2)%BP
                                ipointer=ipointer+1
                            endif
                        enddo
                    endif
                endif
            endif
        enddo
        if (associated(Current_box%Parent)) then
            aux_box => Current_box%Parent
            do count1=1,27
                if (associated(aux_box%Colleague(count1)%BP)) then
                    if (aux_box%Colleague(count1)%BP%is_leaf) then
                        if (is_touching(aux_box%Colleague(count1)%BP,Current_box)) then
                            Current_box%Nearest(ipointer)%BP => aux_box%Colleague(count1)%BP
                            ipointer=ipointer+1
                        endif
                    endif
                endif
            enddo
        endif
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function is_touching(Box_1,Box_2)
implicit none
!! This subroutine evals if two boxes are touching

    !List of calling arguments
    type ( Box ), pointer :: Box_1,Box_2  !! Pointers to two boxes under consideration
    logical :: is_touching              !! Output result

    !List of local variables
    real ( kind = 8 ) Box_2_size

        if ((associated(Box_1)).and.(associated(Box_2))) then
            Box_2_size=Box_2%Box_size+1.0d-14
            if (((Box_1%Box_center(1)+Box_1%Box_size/2.0d0) >= (Box_2%Box_center(1)-Box_2_size/2.0d0)) .and. &
               &((Box_1%Box_center(1)-Box_1%Box_size/2.0d0) <= (Box_2%Box_center(1)+Box_2_size/2.0d0)) .and. &
               &((Box_1%Box_center(2)+Box_1%Box_size/2.0d0) >= (Box_2%Box_center(2)-Box_2_size/2.0d0)) .and. &
               &((Box_1%Box_center(2)-Box_1%Box_size/2.0d0) <= (Box_2%Box_center(2)+Box_2_size/2.0d0)) .and. &
               &((Box_1%Box_center(3)+Box_1%Box_size/2.0d0) >= (Box_2%Box_center(3)-Box_2_size/2.0d0)) .and. &
               &((Box_1%Box_center(3)-Box_1%Box_size/2.0d0) <= (Box_2%Box_center(3)+Box_2_size/2.0d0))) then
                is_touching=(.true.)
            else
                is_touching=(.false.)
            endif
        else
            is_touching=(.false.)
        endif

end function is_touching

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine subdivide_box_once(Current_box)
implicit none
!! This subroutine splits the current box in 8 one time.
!! This will create empty boxes if necessary

    !List of calling arguments
    type ( Box ), pointer :: Current_box  !! pointer to the box that we want to split

    !List of local variables
    integer ( kind = 8 ) count1,count2,count3,count_aux
    integer ( kind = 8 ) n_part_array(2,2,2)    !! This will contain the number of points on each child
    real ( kind = 8 ) r_aux(3),d_aux

        Current_box%is_leaf=(.false.)       !! Declare the current box as not leaf node
        n_part_array( : , : , : ) = 0       !! Initialize the number of points on each child
        do count1=1,Current_box%n_points    !! Loop through all points and see its destination
            if (Current_box%Points(1,count1)>Current_box%Box_center(1)) then
                if (Current_box%Points(2,count1)>Current_box%Box_center(2)) then
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        n_part_array(1,1,1)=n_part_array(1,1,1)+1
                    else
                        n_part_array(1,1,2)=n_part_array(1,1,2)+1
                    endif
                else
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        n_part_array(1,2,1)=n_part_array(1,2,1)+1
                    else
                        n_part_array(1,2,2)=n_part_array(1,2,2)+1
                    endif
                endif
            else
                if (Current_box%Points(2,count1)>Current_box%Box_center(2)) then
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        n_part_array(2,1,1)=n_part_array(2,1,1)+1
                    else
                        n_part_array(2,1,2)=n_part_array(2,1,2)+1
                    endif
                else
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        n_part_array(2,2,1)=n_part_array(2,2,1)+1
                    else
                        n_part_array(2,2,2)=n_part_array(2,2,2)+1
                    endif
                endif
            endif
        enddo

        d_aux=Current_box%Box_size/2.0d0    !! set the sizee of the new boxes created
        do count1=1,2
            do count2=1,2
                do count3=1,2
                    allocate(Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP) !! Allocate all children and set up initial values
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%Parent=>Current_box
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%n_points=n_part_array(count1,count2,count3)
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%is_leaf=.true.
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%primary_violator=.false.
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%secondary_violator=.false.
                    r_aux(1)=Current_box%Box_center(1)-(count1-1.5d0)*d_aux
                    r_aux(2)=Current_box%Box_center(2)-(count2-1.5d0)*d_aux
                    r_aux(3)=Current_box%Box_center(3)-(count3-1.5d0)*d_aux
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%Box_center=r_aux
                    Current_box%Children(count1+2*(count2-1)+4*(count3-1))%BP%Box_size=d_aux
                    if (n_part_array(count1,count2,count3)>0) then
                        allocate(Current_box%Children(count1+2*(count2-1)+4*(count3-1))&
                            &%BP%Points(3,n_part_array(count1,count2,count3)))
                        allocate(Current_box%Children(count1+2*(count2-1)+4*(count3-1))&
                            &%BP%sgmas(n_part_array(count1,count2,count3)))
                    endif
                enddo
            enddo
        enddo
        do count1=1,Current_box%n_points
            if (Current_box%Points(1,count1)>Current_box%Box_center(1)) then
                if (Current_box%Points(2,count1)>Current_box%Box_center(2)) then
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        Current_box%Children(1)%BP%Points(:,n_part_array(1,1,1))=Current_box%Points(:,count1)
                        Current_box%Children(1)%BP%sgmas(n_part_array(1,1,1))=Current_box%sgmas(count1)
                        n_part_array(1,1,1)=n_part_array(1,1,1)-1
                    else
                        Current_box%Children(5)%BP%Points(:,n_part_array(1,1,2))=Current_box%Points(:,count1)
                        Current_box%Children(5)%BP%sgmas(n_part_array(1,1,2))=Current_box%sgmas(count1)
                        n_part_array(1,1,2)=n_part_array(1,1,2)-1
                    endif
                else
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        Current_box%Children(3)%BP%Points(:,n_part_array(1,2,1))=Current_box%Points(:,count1)
                        Current_box%Children(3)%BP%sgmas(n_part_array(1,2,1))=Current_box%sgmas(count1)
                        n_part_array(1,2,1)=n_part_array(1,2,1)-1
                    else
                        Current_box%Children(7)%BP%Points(:,n_part_array(1,2,2))=Current_box%Points(:,count1)
                        Current_box%Children(7)%BP%sgmas(n_part_array(1,2,2))=Current_box%sgmas(count1)
                        n_part_array(1,2,2)=n_part_array(1,2,2)-1
                    endif
                endif
            else
                if (Current_box%Points(2,count1)>Current_box%Box_center(2)) then
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        Current_box%Children(2)%BP%Points(:,n_part_array(2,1,1))=Current_box%Points(:,count1)
                        Current_box%Children(2)%BP%sgmas(n_part_array(2,1,1))=Current_box%sgmas(count1)
                        n_part_array(2,1,1)=n_part_array(2,1,1)-1
                    else
                        Current_box%Children(6)%BP%Points(:,n_part_array(2,1,2))=Current_box%Points(:,count1)
                        Current_box%Children(6)%BP%sgmas(n_part_array(2,1,2))=Current_box%sgmas(count1)
                        n_part_array(2,1,2)=n_part_array(2,1,2)-1
                    endif
                else
                    if (Current_box%Points(3,count1)>Current_box%Box_center(3)) then
                        Current_box%Children(4)%BP%Points(:,n_part_array(2,2,1))=Current_box%Points(:,count1)
                        Current_box%Children(4)%BP%sgmas(n_part_array(2,2,1))=Current_box%sgmas(count1)
                        n_part_array(2,2,1)=n_part_array(2,2,1)-1
                    else
                        Current_box%Children(8)%BP%Points(:,n_part_array(2,2,2))=Current_box%Points(:,count1)
                        Current_box%Children(8)%BP%sgmas(n_part_array(2,2,2))=Current_box%sgmas(count1)
                        n_part_array(2,2,2)=n_part_array(2,2,2)-1
                    endif
                endif
            endif
        enddo
        Current_box%n_points=0
        if (associated(Current_box%Points)) then
            deallocate(Current_box%Points)
        endif
        if (associated(Current_box%sgmas)) then
            deallocate(Current_box%sgmas)
        endif
        call setup_colleagues(Current_box)
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! These functions are for the purpose of making a level restricted tree !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!
!! 1-subroutine make_level_restricted       !! this is the only function that the user will have to use
                                            !! the rest of functions are called from that one
                                            !! example: >> call make_level_restricted(TreeLRD_1)
!! 2-subroutine flag_primary_violators
!! 3-subroutine flag_primary_violators_pair
!! 4-function is_primary_violating
!! 5-subroutine fix_primary_violators
!! 6-subroutine flag_secondary_violators
!! 7-subroutine flag_secondary_violators_pair
!! 8-function is_secondary_violating
!! 9-subroutine fix_secondary_violators
!! 10-subroutine review_secondary
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine make_level_restricted(TreeLRD_1)
implicit none
!! This subroutine transform a fully adaptive tree into a level restricted tree

    !List of calling arguments
    type ( TreeLRD ), pointer :: TreeLRD_1

        call flag_primary_violators(TreeLRD_1%Main_box)       !! Flag and fix whatever is needed
        call flag_secondary_violators(TreeLRD_1%Main_box)
        call fix_primary_violators(TreeLRD_1%Main_box)
        call fix_secondary_violators(TreeLRD_1%Main_box)
        call clean_flags(TreeLRD_1%Main_box)                  !! Clean the flags

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine clean_flags(Current_box)
implicit none
!! This subroutine clean flags related to primary and secondary violators in all the tree
!! it is done recursively through all the tree

    !List of calling arguments
    type ( Box ), pointer :: Current_box    !! Pointer to the current box

    !List of local variables
    integer ( kind = 8 ) count1

        Current_box%primary_violator=.false.
        Current_box%secondary_violator=.false.
        if (.not.(Current_box%is_leaf)) then
            do count1=1,8
                call clean_flags(Current_box%Children(count1)%BP)
            enddo
        endif
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine flag_primary_violators(Current_box)
implicit none
!! This subroutine flags all primary violators in the tree recursively starting from Main_Box

    !List of calling arguments
    type ( Box ), pointer :: Current_box    !! Pointer to the current box

    !List of local variables
    integer ( kind = 8 ) count1,count2

        if (Current_box%is_leaf) then
            do count1=1,27
                if (associated(Current_box%Colleague(count1)%BP)) then
                    call flag_primary_violators_pair(Current_box%Colleague(count1)%BP,Current_box)
                endif
            enddo
        else
            do count2=1,8
                call flag_primary_violators(Current_box%Children(count2)%BP)
            enddo
        endif
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine flag_primary_violators_pair(Box_1,Box_2)
implicit none
!! This subroutine finds if Box_2 is primary violator due to Box_1 or any descendent of Box_1
!! it is done recursively

    !List of calling arguments
    type ( Box ), pointer :: Box_1,Box_2

    !List of local variables
    integer ( kind = 8 ) count1

        if (is_touching(Box_1,Box_2)) then
            if (is_primary_violating(Box_1,Box_2)) then
                Box_2%primary_violator=.true.
            else
                if (.not.(Box_1%is_leaf)) then
                    do count1=1,8
                        call flag_primary_violators_pair(Box_1%Children(count1)%BP,Box_2)
                    enddo
                endif
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function is_primary_violating(Box_1,Box_2)
implicit none
!! This function finds if Box_2 is primary violator due to Box_1 due to size discrepancy
!! They must be touching also, but this is checked before entering to this function

    !List of calling arguments
    type ( Box ), pointer :: Box_1,Box_2
    logical :: is_primary_violating

    !List of local variables
    real ( kind = 8 ) ratio_sizes

        ratio_sizes=Box_1%Box_size/Box_2%Box_size
        is_primary_violating=((ratio_sizes>2.1d0).or.(ratio_sizes<0.49d0))

end function is_primary_violating

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine fix_primary_violators(Current_box)
implicit none
!! This subroutine fix a primary violator by refining as much as necessary
!! It does it for the full tree recursively starting from Current_box

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,count2

        if (Current_box%primary_violator) then      !! First check if there is something to fix in the current box.
                                                    !! if not it will call the recurson to all children
            Current_box%primary_violator=.false.    !! Fix the flag and subdivide, and flag children and fix children
            call subdivide_box_once(Current_box)
            do count1=1,8
                call flag_primary_violators(Current_box%Children(count1)%BP)
            enddo
            do count1=1,8
                call fix_primary_violators(Current_box%Children(count1)%BP)
            enddo
        else
            if (.not.(Current_box%is_leaf)) then
                do count2=1,8
                    call fix_primary_violators(Current_box%Children(count2)%BP)
                enddo
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine flag_secondary_violators(Current_box)
implicit none
!! This subroutine flags all secondary violators in the tree recursively starting from Main_Box

    !List of calling arguments
    type ( Box ), pointer :: Current_box

    !List of local variables
    integer ( kind = 8 ) count1,count2

        if ((Current_box%is_leaf).and.(.not.(Current_box%primary_violator))) then
            do count1=1,27
                if (associated(Current_box%Colleague(count1)%BP)) then
                    call flag_secondary_violators_pair(Current_box%Colleague(count1)%BP,Current_box)
                endif
            enddo
        else
            if (.not.(Current_box%is_leaf)) then
                do count2=1,8
                    call flag_secondary_violators(Current_box%Children(count2)%BP)
                enddo
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine flag_secondary_violators_pair(Box_1,Box_2)
implicit none
!! This function finds if Box_2 is secondary violator due to Box_1 or any descendent of Box_1

    !List of calling arguments
    type ( Box ), pointer :: Box_1,Box_2

    !List of local variables
    integer ( kind = 8 ) count1

        if (is_touching(Box_1,Box_2)) then
            if (is_secondary_violating(Box_1,Box_2)) then
                Box_2%secondary_violator=.true.
            else
                if ((.not.(Box_1%is_leaf)).and.(dabs(Box_1%Box_size/Box_2%Box_size-1.0d0)<1.0d-9)) then
                    do count1=1,8
                        call flag_secondary_violators_pair(Box_1%Children(count1)%BP,Box_2)
                    enddo
                endif
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function is_secondary_violating(Box_1,Box_2)
implicit none
!! This function finds if Box_2 is secondary violator due to Box_1

    !List of calling arguments
    type ( Box ), pointer :: Box_1,Box_2
    logical :: is_secondary_violating

    !List of local variables
    real ( kind = 8 ) ratio_sizes

        is_secondary_violating=((dabs(Box_2%Box_size/Box_1%Box_size-2.0d0)<1.0d-9).and.(Box_1%primary_violator))

end function is_secondary_violating

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine fix_secondary_violators(Current_box)
implicit none
!! This subroutine fixes secondary violators through all the tree starting from Main Box

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,count2

        if (Current_box%secondary_violator) then
            Current_box%secondary_violator=.false.
            call subdivide_box_once(Current_box)
            if (associated(Current_box%Parent)) then
                do count1=1,27
                    if (associated(Current_box%Parent%Colleague(count1)%BP)) then
                        call review_secondary(Current_box%Parent%Colleague(count1)%BP)
                    endif
                enddo
            endif
        else
            if (.not.(Current_box%is_leaf)) then
                do count2=1,8
                    call fix_secondary_violators(Current_box%Children(count2)%BP)
                enddo
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine review_secondary(Current_box)
implicit none
!! This function reviews possible level restrictions due to secondary violators refinement

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1

        if (Current_box%is_leaf) then
            call flag_primary_violators(Current_box)
            if (Current_box%primary_violator) then
                Current_box%primary_violator=.false.
                call subdivide_box_once(Current_box)
                if (associated(Current_box%Parent)) then
                    do count1=1,27
                        if (associated(Current_box%Parent%Colleague(count1)%BP)) then
                            call review_secondary(Current_box%Parent%Colleague(count1)%BP)
                        endif
                    enddo
                endif
            endif
        endif

return
end

function number_colleague(Box_1,Box_2)
implicit none

    !List of calling arguments
    type (Box), pointer :: Box_1,Box_2
    integer ( kind = 8 ) number_colleague

    !List of local variables
    integer ( kind = 8 ) num_points,count1,count2,ipointer
    real ( kind = 8 ) size_aux


    if ((is_touching(Box_1,Box_2)).and.((Box_1%Box_size/Box_2%Box_size-1.0d0)<1.0d-8)) then
        size_aux=Box_1%Box_size/2.0d0
        if (Box_1%Box_center(1)<(Box_2%Box_center(1)-size_aux)) then
                if (Box_1%Box_center(2)<(Box_2%Box_center(2)-size_aux)) then
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=1
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=19
                    else
                        number_colleague=10
                    endif
                elseif (Box_1%Box_center(2)>(Box_2%Box_center(2)+size_aux)) then
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=7
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=25
                    else
                        number_colleague=16
                    endif
                else
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=4
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=22
                    else
                        number_colleague=13
                    endif
                endif
        elseif (Box_1%Box_center(1)>(Box_2%Box_center(1)+size_aux)) then
              if (Box_1%Box_center(2)<(Box_2%Box_center(2)-size_aux)) then
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=3
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=21
                    else
                        number_colleague=12
                    endif
                elseif (Box_1%Box_center(2)>(Box_2%Box_center(2)+size_aux)) then
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=9
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=27
                    else
                        number_colleague=18
                    endif
                else
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=6
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=24
                    else
                        number_colleague=15
                    endif
                endif
        else
              if (Box_1%Box_center(2)<(Box_2%Box_center(2)-size_aux)) then
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=2
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=20
                    else
                        number_colleague=11
                    endif
                elseif (Box_1%Box_center(2)>(Box_2%Box_center(2)+size_aux)) then
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=8
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=26
                    else
                        number_colleague=17
                    endif
                else
                    if (Box_1%Box_center(3)<(Box_2%Box_center(3)-size_aux)) then
                        number_colleague=5
                    elseif (Box_1%Box_center(3)>(Box_2%Box_center(3)+size_aux)) then
                        number_colleague=23
                    else
                        number_colleague=14
                    endif
                endif
        endif
    else
        number_colleague=0
    endif

end function number_colleague



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! These functions are for the purpose of making a defragmented tree    !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!
!! 1-call defrag_tree_Points(TreeLRD_1)  !! this is the only function that the user will have to use
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine defrag_tree_Points(TreeLRD_1)
implicit none
!! This function defragments all the tree placing Points in a particular order in memory

    !List of calling arguments
    type ( TreeLRD ), pointer :: TreeLRD_1

    !List of local variables
    integer ( kind = 8 ) num_points,count1,count2,ipointer

        num_points=total_number_Points(TreeLRD_1%Main_box)
        TreeLRD_1%ntot_W_Pts_mem=num_points
        if (associated(TreeLRD_1%W_Pts_mem)) then
            deallocate(TreeLRD_1%W_Pts_mem)
        endif
        if (associated(TreeLRD_1%W_sgmas_mem)) then
            deallocate(TreeLRD_1%W_sgmas_mem)
        endif
        allocate(TreeLRD_1%W_Pts_mem(3,TreeLRD_1%ntot_W_Pts_mem))
        allocate(TreeLRD_1%W_sgmas_mem(TreeLRD_1%ntot_W_Pts_mem))
        ipointer=1
        call defrag_box_Points(TreeLRD_1,TreeLRD_1%Main_box,ipointer)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine defrag_box_Points(TreeLRD_1, Current_box,ipointer)
implicit none
!! This function defragments all the boxes in a given level
!! (Jexp of each box) placing it in a contiguous array

    !List of calling arguments
    type ( TreeLRD ), pointer :: TreeLRD_1
    type (Box), pointer :: Current_box
    integer ( kind = 8 ), intent(inout) :: ipointer

    !List of local variables
    integer ( kind = 8 ) count1,count2

        if (.not.(Current_box%is_leaf)) then
            do count1=1,8
                call defrag_box_Points(TreeLRD_1, Current_box%Children(count1)%BP,ipointer)
            enddo
            if (Current_box%n_points>0) then
                TreeLRD_1%W_Pts_mem(:,ipointer:ipointer+Current_box%n_points-1)=Current_box%Points
                deallocate(Current_box%Points)
                Current_box%Points => TreeLRD_1%W_Pts_mem(:,ipointer:ipointer+Current_box%n_points-1)
                ipointer=ipointer+Current_box%n_points
            endif
        else
            if (Current_box%n_points>0) then
                TreeLRD_1%W_Pts_mem(:,ipointer:ipointer+Current_box%n_points-1)=Current_box%Points
                deallocate(Current_box%Points)
                Current_box%Points => TreeLRD_1%W_Pts_mem(:,ipointer:ipointer+Current_box%n_points-1)
                TreeLRD_1%W_sgmas_mem(ipointer:ipointer+Current_box%n_points-1)=Current_box%sgmas
                deallocate(Current_box%sgmas)
                Current_box%sgmas => TreeLRD_1%W_sgmas_mem(ipointer:ipointer+Current_box%n_points-1)
                ipointer=ipointer+Current_box%n_points
            endif
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function total_number_Points(Current_box) result(num_box_aux)
implicit none
!! This will tell the total number of boxes in the full tree

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,num_box_aux

        if (Current_box%is_leaf) then
            num_box_aux=Current_box%n_points
        else
            num_box_aux=Current_box%n_points
            do count1=1,8
                num_box_aux=num_box_aux+total_number_Points(Current_box%Children(count1)%BP)
            enddo
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function find_max_depth(Current_box) result(max_depth)
implicit none
!! This function will find the maximum depth of the tree

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,max_depth,max_depth_vect(8)

        if (Current_box%is_leaf) then
            max_depth=1
        else
            do count1=1,8
                max_depth_vect(count1)=find_max_depth(Current_box%Children(count1)%BP)
            enddo
            max_depth=1+maxval(max_depth_vect)
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function number_boxes_at_level(Current_box,n_level) result(num_box_aux)
implicit none
!! This will tell the number of boxes at a particular level n_level

    !List of calling arguments
    type (Box), pointer :: Current_box
    integer ( kind = 8 ), intent(in) :: n_level

    !List of local variables
    integer (kind = 8 ) count1,num_box_aux
        
	 if (n_level.eq.1) then
            num_box_aux=1
        elseif ((Current_box%is_leaf).and.(n_level>1)) then
            num_box_aux=0
        else
            num_box_aux=0
            do count1=1,8
                num_box_aux=num_box_aux+number_boxes_at_level(Current_box%Children(count1)%BP,n_level-1)
            enddo
        endif
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function total_number_boxes(Current_box) result(num_box_aux)
implicit none
!! This will tell the total number of boxes in the full tree

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,num_box_aux
        
        if (Current_box%is_leaf) then
            num_box_aux=1
        else
            num_box_aux=1
            do count1=1,8
                num_box_aux=num_box_aux+total_number_boxes(Current_box%Children(count1)%BP)
            enddo
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive function total_number_leaf_boxes(Current_box) result(num_box_aux)
implicit none
!! This will tell the total number of boxes in the full tree

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,num_box_aux
        if (Current_box%is_leaf) then
            num_box_aux=1
        else
            num_box_aux=0
            do count1=1,8
                num_box_aux=num_box_aux+total_number_leaf_boxes(Current_box%Children(count1)%BP)
            enddo
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine deallocate_tree(TreeLRD_1)
implicit none
!! This is to deallocate the full tree

    !List of calling arguments
    type (TreeLRD), pointer :: TreeLRD_1

    if (associated(TreeLRD_1%W_sgmas_mem)) then
        deallocate(TreeLRD_1%W_sgmas_mem)
    endif
    if (associated(TreeLRD_1%W_Pts_mem)) then
        deallocate(TreeLRD_1%W_Pts_mem)
    endif

    call deallocate_box(TreeLRD_1%Main_box)

    TreeLRD_1%ntot_W_Pts_mem=0
    TreeLRD_1%ntot_W_sgmas_mem=0

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

recursive subroutine deallocate_box(Current_box)
implicit none
!! This is to deallocate the full tree

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer ( kind = 8 ) count1

    if (.not.(Current_box%is_leaf)) then
        do count1=1,8
            call deallocate_box(Current_box%Children(count1)%BP)
        enddo
    endif
!    if (associated(Current_box%Points)) then
!        if (allocated(Current_box%Points)) then
!            deallocate(Current_box%Points)
!        endif
!    endif
!    if (associated(Current_box%sgmas)) then
!        if (allocated(Current_box%sgmas)) then
!            deallocate(Current_box%sgmas)
!        endif
!    endif
    deallocate(Current_box)

return
end

subroutine all_centers_sgmas(TreeLRD_1,centers,sgmas_2,n_boxes)
implicit none

    !List of calling arguments
    type ( TreeLRD ), pointer :: TreeLRD_1
    integer ( kind =  8 ), intent(in) :: n_boxes
    real ( kind = 8 ), intent(inout) :: centers(3,n_boxes),sgmas_2(n_boxes)

    !List of local variables
    integer ( kind = 8 ) ipointer

    ipointer=1
    call all_centers_sgmas_box(TreeLRD_1%Main_box,centers,sgmas_2,n_boxes,ipointer)

return
end


recursive subroutine all_centers_sgmas_box(Current_box,centers,sgmas_2,n_boxes,ipointer)
implicit none

    !List of calling arguments
    type ( Box ), pointer :: Current_box
    integer ( kind =  8 ), intent(in) :: n_boxes
    real ( kind = 8 ), intent(inout) :: centers(3,n_boxes),sgmas_2(n_boxes)
    integer ( kind = 8 ), intent(inout) :: ipointer

    !List of local variables
    integer ( kind = 8 ) count1

    if (Current_box%is_leaf) then
        centers(:,ipointer)=Current_box%Box_center
        sgmas_2(ipointer)=Current_box%Box_size/10.0d0
        ipointer=ipointer+1

    else
        do count1=1,8
            call all_centers_sgmas_box(Current_box%Children(count1)%BP,centers,sgmas_2,n_boxes,ipointer)
        enddo
    endif

return
end



!!!!!!!!!!! This function is for the purpose of debugg !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
recursive subroutine print_all_tree(Current_box)
implicit none
!! This is for debug. It will print flags and sizes of boxes

    !List of calling arguments
    type (Box), pointer :: Current_box

    !List of local variables
    integer (kind = 8 ) count1,count2


    write (*,*) 'Box size:', Current_box%Box_size
!!,Current_box%primary_violator,&
!    & Current_box%secondary_violator,Current_box%Box_size
!    if (Current_box%secondary_violator) then
!        write (*,*) Current_box%Box_size
!    endif
        if (.not.(Current_box%is_leaf)) then
            do count1=1,8
                call print_all_tree(Current_box%Children(count1)%BP)
            enddo
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Mod_TreeLRD
