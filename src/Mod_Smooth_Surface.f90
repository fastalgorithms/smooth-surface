Module Mod_Smooth_Surface
use Mod_Tri_Tools
!use Mod_InOut
use Mod_Plot_Tools
use Mod_Feval
use ModType_Smooth_Surface


implicit none

contains

function my_cross(a, b)

    real ( kind = 8 ) my_cross(3)
    real ( kind = 8 ), intent(in) :: a(3),b(3)

        my_cross(1) = a(2) * b(3) - a(3) * b(2)
        my_cross(2) = a(3) * b(1) - a(1) * b(3)
        my_cross(3) = a(1) * b(2) - a(2) * b(1)

end function my_cross


function find_angle(u_x, u_y, u_z, v_x, v_y, v_z)

    real ( kind = 8 ), intent(in) :: u_x, u_y, u_z, v_x, v_y, v_z
    real ( kind = 8 ) find_angle
    real ( kind = 8 ) num,denom

        num=u_x*v_x+u_y*v_y+u_z*v_z
        denom=sqrt(u_x**2+u_y**2+u_z**2)*sqrt(v_x**2+v_y**2+v_z**2)
        find_angle=dacos(num/denom)
!        write (*,*) u_x, u_y, u_z
!        write (*,*) v_x, v_y, v_z
!        write (*,*) find_angle
!        read (*,*)

end function find_angle



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine funcion_skeleton(Geometry1,n_order_sk)
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
        if (allocated(Geometry1%skeleton_Points)) then
            deallocate(Geometry1%skeleton_Points)
        endif
        if (allocated(Geometry1%skeleton_w)) then
            deallocate(Geometry1%skeleton_w)
        endif
        if (allocated(Geometry1%skeleton_N)) then
            deallocate(Geometry1%skeleton_N)
        endif
        allocate(Geometry1%skeleton_Points(3,Geometry1%n_Sk_points))
        allocate(Geometry1%skeleton_w(Geometry1%n_Sk_points))
        allocate(Geometry1%skeleton_N(3,Geometry1%n_Sk_points))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,n_order)
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
            F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)&
             &*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
            F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)&
             &*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
            F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)&
             &*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order)
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
            F_x(count)=coef_x(1)+coef_x(2)*U(count)+coef_x(3)*V(count)+coef_x(4)&
             &*U(count)**2+coef_x(5)*V(count)**2+coef_x(6)*U(count)*V(count)
            F_y(count)=coef_y(1)+coef_y(2)*U(count)+coef_y(3)*V(count)+coef_y(4)&
             &*U(count)**2+coef_y(5)*V(count)**2+coef_y(6)*U(count)*V(count)
            F_z(count)=coef_z(1)+coef_z(2)*U(count)+coef_z(3)*V(count)+coef_z(4)&
             &*U(count)**2+coef_z(5)*V(count)**2+coef_z(6)*U(count)*V(count)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine funcion_normal_vert(Geometry1)
implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    type (My_cell) My_cell1
    integer ( kind = 8 ) count, n_order
    real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    real ( kind = 8 ) U(3),V(3)
    real ( kind = 8 ) F_x(3),F_y(3),F_z(3),dS(3)
    real ( kind = 8 ) U_x(3),U_y(3),U_z(3),V_x(3),V_y(3),V_z(3)
    real ( kind = 8 ) angle_vertex(3)
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
            if (My_cell1%Var_Mat(count)%n_Mat>0) then
                allocate(My_cell1%Var_Mat(count)%Mat(4,My_cell1%Var_Mat(count)%n_Mat))
            endif
        enddo
        do count=1,Geometry1%ntri
            P1=Geometry1%Points(:,Geometry1%Tri(1,count))
            P2=Geometry1%Points(:,Geometry1%Tri(2,count))
            P3=Geometry1%Points(:,Geometry1%Tri(3,count))
            P4=Geometry1%Points(:,Geometry1%Tri(4,count))
            P5=Geometry1%Points(:,Geometry1%Tri(5,count))
            P6=Geometry1%Points(:,Geometry1%Tri(6,count))
!            call eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,nP_x,nP_y,nP_z,dS,n_order)

         call eval_quadratic_patch_UV(P1,P2,P3,P4,P5,P6,U,V,F_x,F_y,F_z,U_x,U_y,U_z,V_x,V_y,V_z,nP_x,nP_y,nP_z,dS,n_order)


            My_cell1%Var_Mat(Geometry1%Tri(1,count))%Mat(1:3,My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat)&
            &=(/nP_x(1), nP_y(1), nP_z(1)/)
            angle_vertex(1)=find_angle(U_x(1),U_y(1),U_z(1),V_x(1),V_y(1),V_z(1))
            My_cell1%Var_Mat(Geometry1%Tri(1,count))%Mat(4,My_cell1%Var_Mat&
             &(Geometry1%Tri(1,count))%current_n_Mat)=angle_vertex(1)
            My_cell1%Var_Mat(Geometry1%Tri(1,count))%current_n_Mat=My_cell1%Var_Mat&
             &(Geometry1%Tri(1,count))%current_n_Mat+1

!            write (*,*) Geometry1%Points(:,Geometry1%Tri(1,count))
!            write (*,*) U_x(1),U_y(1),U_z(1)
!            write (*,*) V_x(1),V_y(1),V_z(1)
!            write (*,*) angle_vertex(1)/3.14159268d0

            My_cell1%Var_Mat(Geometry1%Tri(2,count))%Mat(1:3,My_cell1%Var_Mat&
             &(Geometry1%Tri(2,count))%current_n_Mat)=(/nP_x(2), nP_y(2), nP_z(2)/)
            angle_vertex(2)=find_angle(-U_x(2),-U_y(2),-U_z(2),V_x(2)-U_x(2),V_y(2)-U_y(2),V_z(2)-U_z(2))
            My_cell1%Var_Mat(Geometry1%Tri(2,count))%Mat(4,My_cell1%Var_Mat&
             &(Geometry1%Tri(2,count))%current_n_Mat)=angle_vertex(2)
            My_cell1%Var_Mat(Geometry1%Tri(2,count))%current_n_Mat=My_cell1%Var_Mat&
             &(Geometry1%Tri(2,count))%current_n_Mat+1

!            write (*,*) Geometry1%Points(:,Geometry1%Tri(2,count))
!            write (*,*) -U_x(2),-U_y(2),-U_z(2)
!            write (*,*) V_x(2)-U_x(2),V_y(2)-U_y(2),V_z(2)-U_z(2)
!            write (*,*) angle_vertex(2)/3.14159268d0


            My_cell1%Var_Mat(Geometry1%Tri(3,count))%Mat(1:3,My_cell1%Var_Mat&
             &(Geometry1%Tri(3,count))%current_n_Mat)=(/nP_x(3), nP_y(3), nP_z(3)/)
            angle_vertex(3)=find_angle(-V_x(3)+U_x(3),-V_y(3)+U_y(3),-V_z(3)+U_z(3),-V_x(3),-V_y(3),-V_z(3))
            My_cell1%Var_Mat(Geometry1%Tri(3,count))%Mat(4,My_cell1%Var_Mat&
             &(Geometry1%Tri(3,count))%current_n_Mat)=angle_vertex(3)
            My_cell1%Var_Mat(Geometry1%Tri(3,count))%current_n_Mat=My_cell1%Var_Mat&
             &(Geometry1%Tri(3,count))%current_n_Mat+1

!            write (*,*) Geometry1%Points(:,Geometry1%Tri(3,count))
!            write (*,*) -V_x(3)+U_x(3),-V_y(3)+U_y(3),-V_z(3)+U_z(3)
!            write (*,*) -V_x(3),-V_y(3),-V_z(3)
!            write (*,*) angle_vertex(3)/3.14159268d0

!            read (*,*)

        enddo
        if (allocated(Geometry1%Normal_Vert)) then
            deallocate(Geometry1%Normal_Vert)
        endif
        allocate(Geometry1%Normal_Vert(3,Geometry1%npoints))
        do count=1,Geometry1%npoints
            if (My_cell1%Var_Mat(count)%n_Mat>0) then
!                write (*,*) count
!                write (*,*) My_cell1%Var_Mat(count)%Mat(1:3,:)
!                write (*,*) count
!                write (*,*) My_cell1%Var_Mat(count)%Mat(4,:)
!                write (*,*) 'n point: ',count
!                write (*,*) sum(My_cell1%Var_Mat(count)%Mat(4,:))/3.14159268d0
!                write (*,*)
!                write (*,*) Geometry1%Points(:,count)
!                read (*,*)
!                write (*,*) 'n_mat: ',My_cell1%Var_Mat(count)%n_Mat
!                write (*,*) My_cell1%Var_Mat(count)%Mat(4,:)
!                write (*,*) My_cell1%Var_Mat(count)%Mat(1:3,:)

                call mimean(My_cell1%Var_Mat(count)%Mat,My_cell1%Var_Mat(count)%n_Mat,current_Normal_vert)
                Geometry1%Normal_Vert(:,count)=current_Normal_vert
!                write (*,*) 'aqui: ', current_Normal_vert
!                read (*,*)
            else
                Geometry1%Normal_Vert(:,count)=(/0.d0, 0.d0, 0.d0/)
            endif
        enddo
!        write (*,*) Geometry1%Normal_Vert

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine mimean(All_Normals,n_Normals,Current_Normal)
implicit none

    !List of calling arguments
    integer ( kind = 8 ), intent(in) :: n_Normals
    real ( kind = 8 ), intent(in) :: All_Normals(4,n_Normals)
    real ( kind = 8 ), intent(out) :: Current_Normal(3)

    !List of local variables
    integer ( kind = 8 ) count
    real ( kind = 8 ) coef(n_Normals),sum_coef

        current_Normal=(/ 0.0d0, 0.0d0, 0.0d0 /)
        sum_coef=0
        do count=1,n_Normals
            sum_coef=sum_coef+All_Normals(4,count)
        enddo
        do count=1,n_Normals
            Current_Normal=Current_Normal+All_Normals(1:3,count)*All_Normals(4,count)
        enddo
        Current_Normal=Current_Normal/sum_coef
!        write (*,*) 'all normals: ', All_Normals(1:3,:)
!        write (*,*) 'all coefs: ', All_Normals(4,:)
!        write (*,*) Current_Normal,sum_coef
!        read (*,*)
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine find_cos(v,All_Normals,coef,n_Normals)
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine refine_geometry(Geometry1,n_order_sf)
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

        if (allocated(Geometry1%Points)) then
            deallocate(Geometry1%Points)
        endif
        allocate(Geometry1%Points(3,m))
        if (allocated(Geometry1%Tri)) then
            deallocate(Geometry1%Tri)
        endif
        allocate(Geometry1%Tri(6,N))
        if (allocated(Geometry1%Normal_Vert)) then
            deallocate(Geometry1%Normal_Vert)
        endif
        allocate(Geometry1%Normal_Vert(3,m))

        Geometry1%Points=Points
        Geometry1%Tri=Tri
        Geometry1%Normal_Vert=Normal_Vert

        deallocate(Points)
        deallocate(Normal_Vert)
        deallocate(Tri)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine funcion_Base_Points(Geometry1)
implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    real ( kind = 8 ) U(Geometry1%n_order_sf),V(Geometry1%n_order_sf),w(Geometry1%n_order_sf)
    real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3),N1(3),N2(3),N3(3)
    real ( kind = 8 ) F_x(Geometry1%n_order_sf),F_y(Geometry1%n_order_sf)
    real ( kind = 8 ) F_z(Geometry1%n_order_sf),dS(Geometry1%n_order_sf)
    real ( kind = 8 ) nP_x(Geometry1%n_order_sf),nP_y(Geometry1%n_order_sf),nP_z(Geometry1%n_order_sf)
    real ( kind = 8 ) U_x(Geometry1%n_order_sf),U_y(Geometry1%n_order_sf),U_z(Geometry1%n_order_sf)
    real ( kind = 8 ) V_x(Geometry1%n_order_sf),V_y(Geometry1%n_order_sf),V_z(Geometry1%n_order_sf)
    integer ( kind = 8 ) count,n_order_sf

        n_order_sf=Geometry1%n_order_sf

        if (n_order_sf==45) then
            call GaussTri45(U,V,w)
        else if (n_order_sf==78) then
            call GaussTri78(U,V,w)
        end if
        if (allocated(Geometry1%Base_Points)) then
            deallocate(Geometry1%Base_Points)
        endif
        if (allocated(Geometry1%Base_Points_N)) then
            deallocate(Geometry1%Base_Points_N)
        endif
        if (allocated(Geometry1%Base_Points_U)) then
            deallocate(Geometry1%Base_Points_U)
        endif
        if (allocated(Geometry1%Base_Points_V)) then
            deallocate(Geometry1%Base_Points_V)
        endif
        if (allocated(Geometry1%w_smooth)) then
            deallocate(Geometry1%w_smooth)
        endif
        allocate(Geometry1%Base_Points(3,Geometry1%n_Sf_points))
        allocate(Geometry1%Base_Points_N(3,Geometry1%n_Sf_points))
        allocate(Geometry1%Base_Points_U(3,Geometry1%n_Sf_points))
        allocate(Geometry1%Base_Points_V(3,Geometry1%n_Sf_points))
        allocate(Geometry1%w_smooth(Geometry1%n_Sf_points))
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine find_smooth_surface(Geometry1,Feval_stuff_1,adapt_flag)
implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    type ( Feval_stuff ), pointer :: Feval_stuff_1      !! data type that
    integer ( kind = 8 ), intent(in) :: adapt_flag

    !List of local variables
    integer ( kind = 8 ) flag,maxiter,ipoint,count,n_order_sf,itri
    real ( kind = 8 ) h(Geometry1%n_Sf_points),tol,Norm_N
    real ( kind = 8) r_t(3,Geometry1%n_Sf_points),F(Geometry1%n_Sf_points),grad_F(3,Geometry1%n_Sf_points)
    real ( kind = 8 ) dVdu(3),dVdv(3),dhdu,dhdv,h_var
write (*,*) 'telita2'


        if (allocated(Geometry1%S_smooth)) then
            deallocate(Geometry1%S_smooth)
        endif
        allocate(Geometry1%S_smooth(3,Geometry1%n_Sf_points))
        if (allocated(Geometry1%N_smooth)) then
            deallocate(Geometry1%N_smooth)
        endif
        allocate(Geometry1%N_smooth(3,Geometry1%n_Sf_points))
        if (allocated(Geometry1%ru_smooth)) then
            deallocate(Geometry1%ru_smooth)
        endif
        allocate(Geometry1%ru_smooth(3,Geometry1%n_Sf_points))
        if (allocated(Geometry1%rv_smooth)) then
            deallocate(Geometry1%rv_smooth)
        endif
        allocate(Geometry1%rv_smooth(3,Geometry1%n_Sf_points))
        if (.not.allocated(Geometry1%height)) then
            write (*,*) 'not allocated'
!            read (*,*)
            allocate (Geometry1%height(Geometry1%n_Sf_points))
            do count=1,Geometry1%n_Sf_points
                Geometry1%height(count)=0.0d0
                h(count)=0.0d0
            enddo
        else

!            write (*,*) 'not allocated'
!            read (*,*)
!            deallocate (Geometry1%height)
!            allocate (Geometry1%height(Geometry1%n_Sf_points))
!            do count=1,Geometry1%n_Sf_points
!                Geometry1%height(count)=0.0d0
!                h(count)=0.0d0
!            enddo




            write (*,*) 'yes allocated'
!            read (*,*)
            do count=1,Geometry1%n_Sf_points
                h(count)=Geometry1%height(count)
            enddo
        endif
        tol=1.0d-14
        maxiter=14
        flag=0
        write (*,*) 'No llega nunca1'
        call My_Newton(h,tol,maxiter,Geometry1,flag,Feval_stuff_1,adapt_flag,grad_F,r_t)

        if (flag==1) then
            write (*,*) 'ERROR DE CONVERGENCIA NEWTON'
        end if
        do count=1,Geometry1%n_Sf_points
!            Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+Geometry1%Base_Points_N(2,count)&
!            &**2+Geometry1%Base_Points_N(3,count)**2)


            !!!r_t(:,count)=Geometry1%Base_Points(:,count)+Geometry1%Base_Points_N(:,count)*h(count)/Norm_N
!            r_t(:,count)=Geometry1%Base_Points(:,count)+Geometry1%Base_Points_N(:,count)*h(count)

            Geometry1%S_smooth(:,count)=r_t(:,count)
        enddo
!        call eval_density_grad_FMM(Geometry1,r_t,Geometry1%n_Sf_points,F,grad_F,Feval_stuff_1,adapt_flag)
!       call eval_density_grad(Geometry1,r_t(1,:),r_t(2,:),r_t(3,:),alpha,F,grad_F)
        n_order_sf=Geometry1%n_Sf_points/Geometry1%ntri
        do count=1,Geometry1%n_Sf_points
            Geometry1%N_smooth(:,count)=-1.0d0*grad_F(:,count)/(sqrt(grad_F(1,count)**2+grad_F(2,count)**2+grad_F(3,count)**2))
            itri=(count-1)/n_order_sf+1
            Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+Geometry1%Base_Points_N(2,count)&
            &**2+Geometry1%Base_Points_N(3,count)**2)

            !!!!h_var=h(count)/Norm_N
            h_var=h(count)

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

        do count=1,Geometry1%n_Sf_points
            Geometry1%height(count)=h(count)
        enddo

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine My_Newton(x,tol,maxiter,Geometry1,flag,Feval_stuff_1,adapt_flag,grad_F,r_t) !!Hay que modificar esto y vectorizarlo
implicit none

    !List of calling arguments
    type (Geometry), intent(in) :: Geometry1
    integer ( kind = 8 ), intent(in) :: maxiter
    real ( kind = 8 ), intent(inout) :: x(Geometry1%n_Sf_points)
    real ( kind = 8 ), intent(in) :: tol
    integer ( kind = 8 ), intent(out) :: flag
    type ( Feval_stuff ), pointer :: Feval_stuff_1      !! data type that
    integer ( kind = 8 ), intent(in) :: adapt_flag
    real ( kind = 8), intent(inout) :: r_t(3,Geometry1%n_Sf_points),grad_F(3,Geometry1%n_Sf_points)


    !List of local variables
    integer ( kind = 8 ) count,count2
    real ( kind = 8 ) F(Geometry1%n_Sf_points),dF(Geometry1%n_Sf_points),err(Geometry1%n_Sf_points)
    integer (kind = 8 ) flag_con(Geometry1%n_Sf_points)

        do count2=1,Geometry1%n_Sf_points
            err(count2)=tol+1.0d0
            flag_con(count2)=0
        enddo
        count=1
        flag=0
        do while ((maxval(err)>tol).and.(count<maxiter))
            write (*,*) 'Number of living targets: ', Geometry1%n_Sf_points-sum(flag_con),sum(flag_con)
            write (*,*) 'Ratio of living targets: ',(real(Geometry1%n_Sf_points-sum(flag_con),8))/(real(Geometry1%n_Sf_points,8))
            call fun_roots_derivative(x,Geometry1,F,dF,Feval_stuff_1,adapt_flag,flag_con,grad_F,r_t)
            count=count+1
            do count2=1,Geometry1%n_Sf_points
                if (flag_con(count2)==0) then
                    err(count2)=F(count2)/dF(count2)
                    x(count2)=x(count2)-err(count2)
                    err(count2)=abs(err(count2))
                    if (err(count2)<tol) then
                        flag_con(count2)=1
                    endif
                endif
            enddo
    !       do count2=1,Geometry1%n_Sf_points
    !           write (*,*) 'Newton iteration all: ', err(count2)
    !       enddo
        write (*,*) 'Newton iteration: ', count, maxval(err)
        end do
        if (maxval(err)>tol) then
            flag=1
    !       write (*,*) 'contador máximo Newton', count, err
        endif

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine check_Gauss(Geometry1,x0,y0,z0,err_rel)
implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1
    real ( kind = 8 ), intent(in) :: x0,y0,z0
    real ( kind = 8 ), intent(out) :: err_rel

    !List of local variables
    integer ( kind = 8 ) umio,count1,count2,flag,n_order_sf
    real ( kind = 8 )  F,Ex,Ey,Ez,R,x,y,z,pi,w,nx,ny,nz

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fun_roots_derivative(h,Geometry1,F,dF,Feval_stuff_1,adapt_flag,flag_con,grad_F,r_t)
implicit none

    !List of calling arguments
    type (Geometry), intent(in) :: Geometry1
    real ( kind = 8 ), intent(in) :: h(Geometry1%n_Sf_points)
    real ( kind = 8 ), intent(out) ::  F(Geometry1%n_Sf_points), dF(Geometry1%n_Sf_points)
    type ( Feval_stuff ), pointer :: Feval_stuff_1      !! data type that
    integer ( kind = 8 ), intent(in) :: adapt_flag
    integer ( kind = 8 ), intent(in) :: flag_con(Geometry1%n_Sf_points)
    real ( kind = 8), intent(inout) :: r_t(3,Geometry1%n_Sf_points),grad_F(3,Geometry1%n_Sf_points)



    !List of local variables
    real ( kind = 8 ) Norm_N
    integer (kind = 8 ) count,ipointer
real ( kind = 8 ), allocatable :: r_t2(:,:),v_norm(:,:),grad_F2(:,:),F2(:)

        allocate(r_t2(3,Geometry1%n_Sf_points-sum(flag_con)))
        allocate(v_norm(3,Geometry1%n_Sf_points-sum(flag_con)))
        allocate(grad_F2(3,Geometry1%n_Sf_points-sum(flag_con)))
        allocate(F2(Geometry1%n_Sf_points-sum(flag_con)))

        write (*,*) 'a punto'
        ipointer=1
        do count=1,Geometry1%n_Sf_points
            Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+&
            &Geometry1%Base_Points_N(2,count)**2+Geometry1%Base_Points_N(3,count)**2)

            !!!!r_t(:,count)=Geometry1%Base_Points(:,count)+Geometry1%Base_Points_N(:,count)*h(count)/Norm_N
            r_t(:,count)=Geometry1%Base_Points(:,count)+Geometry1%Base_Points_N(:,count)*h(count)

            if (flag_con(count).eq.0) then
                r_t2(:,ipointer)=r_t(:,count)
                v_norm(:,ipointer)=Geometry1%Base_Points_N(:,count)
                ipointer=ipointer+1
            endif
        enddo
!        write (*,*) Geometry1%n_Sf_points-sum(flag_con),size(r_t2)/3,size(F2),size(grad_F2)/3
!        read (*,*)

        call eval_density_grad_FMM(Geometry1,r_t2,v_norm,Geometry1%n_Sf_points-sum(flag_con),F2,grad_F2,Feval_stuff_1,adapt_flag)
        ipointer=1
        do count=1,Geometry1%n_Sf_points

            if (flag_con(count).eq.0) then
                Norm_N=sqrt(Geometry1%Base_Points_N(1,count)**2+&
                &Geometry1%Base_Points_N(2,count)**2+Geometry1%Base_Points_N(3,count)**2)

                !!!!dF(count)=(grad_F(1,count)*Geometry1%Base_Points_N(1,count)+grad_F(2,count)*Geometry1%Base_Points_N(2,count)&
                !!!!&+grad_F(3,count)*Geometry1%Base_Points_N(3,count))/Norm_N
                dF(count)=(grad_F2(1,ipointer)*Geometry1%Base_Points_N(1,count)+grad_F2(2,ipointer)*&
                &Geometry1%Base_Points_N(2,count)+grad_F2(3,ipointer)*Geometry1%Base_Points_N(3,count))
                F(count)=F2(ipointer)
                r_t(:,count)=r_t2(:,ipointer)
                grad_F(:,count)=grad_F2(:,ipointer)
                ipointer=ipointer+1
            endif

        enddo
        deallocate(r_t2)
        deallocate(v_norm)
        deallocate(grad_F2)
        deallocate(F2)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine refine_geometry_smart(Geometry1)
implicit none

    !List of calling arguments
    type (Geometry), intent(inout) :: Geometry1

    !List of local variables
    character ( len=100 ) plot_name
    integer ( kind = 8 ) count,contador_indices
    real ( kind = 8 ), allocatable :: Points(:,:), Normal_Vert(:,:)
    integer ( kind = 8 ), allocatable :: Tri(:,:)
    real ( kind = 8 ), allocatable :: h_new(:)
    real ( kind = 8 ) h_tri(Geometry1%n_order_sf)
    real ( kind = 8 ) h_1(Geometry1%n_order_sf),h_2(Geometry1%n_order_sf)
    real ( kind = 8 ) h_3(Geometry1%n_order_sf),h_4(Geometry1%n_order_sf)
    real ( kind = 8 ) P1(3),P2(3),P3(3),P4(3),P5(3),P6(3)
    real ( kind = 8 ) Pa(3),Pb(3),Pc(3),Pd(3),Pe(3),Pf(3),Pg(3),Ph(3),Pi(3)
    real ( kind = 8 ) Nor1(3),Nor2(3),Nor3(3),Nor4(3),Nor5(3),Nor6(3)
    real ( kind = 8 ) Nor_a(3),Nor_b(3),Nor_c(3),Nor_d(3),Nor_e(3),Nor_f(3),Nor_g(3),Nor_h(3),Nor_i(3)
    real ( kind = 8 ) U(9),V(9)
    real ( kind = 8 ) F_x(9),F_y(9),F_z(9),dS(9)
    real ( kind = 8 ) nP_x(9),nP_y(9),nP_z(9)
    real ( kind = 8 ) U_x(9),U_y(9),U_z(9),V_x(9),V_y(9),V_z(9)
    integer ( kind = 8 ) m,N,n_order_aux,n_order_sf
    real( kind = 8 ) U45(45),V45(45),w45(45)
    real( kind = 8 ) coef_h(45)

        n_order_sf=Geometry1%n_order_sf
        call GaussTri45(U45,V45,w45)
        allocate(Points(3,Geometry1%ntri*15))
        allocate(Normal_Vert(3,Geometry1%ntri*15))
        allocate(Tri(6,Geometry1%ntri*4))
        allocate(h_new(Geometry1%n_Sf_points*4))
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
            Tri(:,(count-1)*4+2)=contador_indices*[1, 1, 1, 1, 1, 1]+[11, 9, 2, 10, 6, 7]![2, 11, 9, 7, 10, 6]
            Tri(:,(count-1)*4+3)=contador_indices*[1, 1, 1, 1, 1, 1]+[2, 4, 11, 3, 8, 7]
            Tri(:,(count-1)*4+4)=contador_indices*[1, 1, 1, 1, 1, 1]+[9, 11, 14, 10, 13, 12]

            h_tri=Geometry1%height((count-1)*n_order_sf+1:(count)*n_order_sf)
 !           write (*,*) h_tri
 !           read (*,*)
!            get coefs from h_tri
!evaluate the new height on each of the 4 triangles


!!!            call fast_Matvec(M_inv,h_tri,coef_h,45)
!            write (*,*) coef_h
!            read (*,*)

!!!            call pol_val_2D_45_fast2(U45/2.0d0,V45/2.0d0,45,coef_h,h_1)
!!!            call pol_val_2D_45_fast2(0.5d0-U45/2.0d0,0.5d0-V45/2.0d0,45,coef_h,h_2)
!!!            call pol_val_2D_45_fast2(0.5d0+U45/2.0d0,V45/2.0d0,45,coef_h,h_3)
!!!            call pol_val_2D_45_fast2(U45/2.0d0,0.5d0+V45/2.0d0,45,coef_h,h_4)


        call refine_tri45(h_tri,h_1,h_2,h_3,h_4)

            h_new((count-1)*n_order_sf*4+1:(count-1)*n_order_sf*4+n_order_sf)=h_1
            h_new((count-1)*n_order_sf*4+n_order_sf+1:(count-1)*n_order_sf*4+2*n_order_sf)=h_2
            h_new((count-1)*n_order_sf*4+2*n_order_sf+1:(count-1)*n_order_sf*4+3*n_order_sf)=h_3
            h_new((count-1)*n_order_sf*4+3*n_order_sf+1:(count-1)*n_order_sf*4+4*n_order_sf)=h_4

            contador_indices=contador_indices+15

            !plot_name='./plot_tools/h_tri'
            !call plot_curve_3D(U45,V45,h_tri,n_order_sf,plot_name)
            !plot_name='./plot_tools/h_1'
            !call plot_curve_3D(U45/2.0d0,V45/2.0d0,h_1,n_order_sf,plot_name)
            !plot_name='./plot_tools/h_2'
            !call plot_curve_3D(0.5d0-U45/2.0d0,0.5d0-V45/2.0d0,h_2,n_order_sf,plot_name)
            !plot_name='./plot_tools/h_3'
            !call plot_curve_3D(0.5d0+U45/2.0d0,V45/2.0d0,h_3,n_order_sf,plot_name)
            !plot_name='./plot_tools/h_4'
            !call plot_curve_3D(U45/2.0d0,0.5d0+V45/2.0d0,h_4,n_order_sf,plot_name)
            !write (*,*) 'DONE'
            !read (*,*)


        enddo
        Geometry1%npoints=Geometry1%ntri*15
        Geometry1%ntri=Geometry1%ntri*4
        m=Geometry1%npoints
        N=Geometry1%ntri
        Geometry1%n_Sf_points=N*n_order_sf

        if (allocated(Geometry1%Points)) then
            deallocate(Geometry1%Points)
        endif
        allocate(Geometry1%Points(3,m))
        if (allocated(Geometry1%Tri)) then
            deallocate(Geometry1%Tri)
        endif
        allocate(Geometry1%Tri(6,N))
        if (allocated(Geometry1%Normal_Vert)) then
            deallocate(Geometry1%Normal_Vert)
        endif
        allocate(Geometry1%Normal_Vert(3,m))

        Geometry1%Points=Points
        Geometry1%Tri=Tri
        Geometry1%Normal_Vert=Normal_Vert

        deallocate(Points)
        deallocate(Normal_Vert)
        deallocate(Tri)


        deallocate(Geometry1%height)

        allocate(Geometry1%height(Geometry1%n_Sf_points))
        Geometry1%height=h_new
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end module Mod_Smooth_Surface

