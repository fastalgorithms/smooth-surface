!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                    MODULE MAIN INFORMATION                          !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!
!!! This module contains functions to open a msh file, put all the information into Geometry
!!! and record the resulting smooth surface into a .gov file
!!!
!!!
!!! INDEX OF SUBROUTINES AND FUNCTIONS:
!!!
!!! Public:
!!!
!!! 1- readmsh          !! Read a geometry froma  msh file (read the skeleton and record it into
!!!                     !! a data of type Geometry)
!!! 2- record_Geometry  !! Save the final smooth geometry in a file together with u,v,n vectors and weights w.
!!!
!!! Private: (none)
!!!
!!! MODULES REQUIRED:
!!!
!!! 1-ModType_Smooth_Surface
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Module Mod_InOut
use ModType_Smooth_Surface
implicit none

public :: readmsh
public :: record_Geometry

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    DEFINITION OF SUBROUTINES     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains


subroutine readgeometry(Geometry1,filename,n_order_sk,n_order_sf)
implicit none
!! This subroutine open a msh file and load the information in a variable of type Geometry

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
    character(len=100), intent(in) :: filename         !! name of the msh file
    integer ( kind = 8 ), intent(in) :: n_order_sk,n_order_sf   !! order of the skeleton and surface

    !List of local variables

        if (INDEX(filename,'.msh')>0) then
            call readmsh(Geometry1,filename,n_order_sk,n_order_sf)
        elseif (INDEX(filename,'.tri')>0) then
            call readtri(Geometry1,filename,n_order_sk,n_order_sf)
        else
            write (*,*) 'Geometry type not recognized'
            stop
        endif


return
end


subroutine readtri(Geometry1,filename,n_order_sk,n_order_sf)
implicit none
!! This subroutine open a msh file and load the information in a variable of type Geometry

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
    character(len=100), intent(in) :: filename         !! name of the msh file
    integer ( kind = 8 ), intent(in) :: n_order_sk,n_order_sf   !! order of the skeleton and surface

    !List of local variables
    integer ( kind = 8 ) umio,i,m,N,j,aux1,aux2,aux3,ipointer
    integer :: ierror

        open(UNIT=8, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
        read(8,*) m, N
        write (*,*) 'npoints: ',m,'ntri: ',N,n_order_sf,n_order_sk
        Geometry1%n_order_sf=n_order_sf
        Geometry1%npoints=m+N*3
        Geometry1%ntri=N
        Geometry1%ntri_sk=N
        Geometry1%n_Sf_points=N*n_order_sf
        Geometry1%n_Sk_points=N*n_order_sk
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
end



subroutine readmsh(Geometry1,filename,n_order_sk,n_order_sf)
implicit none
!! This subroutine open a msh file and load the information in a variable of type Geometry

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1     !! where the geometry will be loaded
    character(len=100), intent(in) :: filename         !! name of the msh file
    integer ( kind = 8 ), intent(in) :: n_order_sk,n_order_sf   !! order of the skeleton and surface

    !List of local variables
    integer ( kind = 8 ) umio,i,m,N,j,aux1,aux2,aux3,aux4,aux5,aux6,aux7,aux8
    integer :: ierror
        open(UNIT=8, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=ierror)
        read(8,*) aux1,aux2,aux3,m, N
        write (*,*) 'npoints: ',m,'ntri: ',N,n_order_sf,n_order_sk
        Geometry1%n_order_sf=n_order_sf
        Geometry1%npoints=m
        Geometry1%ntri=N
        Geometry1%ntri_sk=N
        Geometry1%n_Sf_points=N*n_order_sf
        Geometry1%n_Sk_points=N*n_order_sk
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
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine record_Geometry(Geometry1,filename)
implicit none

    !List of calling arguments
    type (Geometry), intent(in) :: Geometry1
    character (len=100) filename

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!    END DEFINITION OF SUBROUTINES    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module Mod_InOut
