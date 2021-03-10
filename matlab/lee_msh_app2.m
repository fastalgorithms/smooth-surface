function lee_msh_app2(filename,ha,N,flag_col,flag_rot)%(filename)
    fid = fopen(filename);
%    N=40;
%    fid = fopen('prueba_quadratic.msh');


    F = fscanf(fid, '%g %g', [1 inf]);
    fclose(fid);
    
    n_points=F(4);
    n_tri=F(5);
    F=F(6:end);
    X=F(9:11:9+11*(n_points-1));
    Y=F(10:11:10+11*(n_points-1));
    Z=F(11:11:11+11*(n_points-1));
    F=F(11+11*(n_points-1)+2:end);
    P1=F(9:14:9+14*(n_tri-1));
    P2=F(10:14:10+14*(n_tri-1));
    P3=F(11:14:11+14*(n_tri-1));
    P4=F(12:14:12+14*(n_tri-1));
    P5=F(13:14:13+14*(n_tri-1));
    P6=F(14:14:14+14*(n_tri-1));
    X=X';
    Y=Y';
    Z=Z';
    P1=P1';
    P2=P2';
    P3=P3';
    P4=P4';
    P5=P5';
    P6=P6';
    Points=[X,Y,Z];
    Tri=[P1 P2 P3 P4 P5 P6];
    plot_quadratic_mesh(Points,Tri,n_tri,ha,N,flag_col,flag_rot)
end



function plot_quadratic_mesh(Points,Tri,n_tri,ha,N,flag_col,flag_rot)
    Mm1=[    1     0     0     0     0     0;
    -3    -1     0     4     0     0;
    -3     0    -1     0     0     4;
     2     2     0    -4     0     0;
     2     0     2     0     0    -4;
     4     0     0    -4     4    -4];
    u=linspace(0,1,N+1);
    [U,V]=meshgrid(u);  
    icount1=1;
    icount2=1;
    for count1=0:(N-1)
        for count2=0:(N-1-count1)
            vert_flat(1,icount1)=U(count1+1,count2+1);
            vert_flat(2,icount1)=V(count1+1,count2+1);
            vert_flat(1,icount1+1)=U(count1+1+1,count2+1);
            vert_flat(2,icount1+1)=V(count1+1+1,count2+1);
            vert_flat(1,icount1+2)=U(count1+1,count2+1+1);
            vert_flat(2,icount1+2)=V(count1+1,count2+1+1);  
            tri_flat(1,icount2)=icount1;
            tri_flat(2,icount2)=icount1+1;
            tri_flat(3,icount2)=icount1+2;
            icount1=icount1+3;
            icount2=icount2+1;
        end
    end

    for count1=1:(N-1)
        for count2=1:(N-1-count1+1)
            vert_flat(1,icount1)=U(count1+1,count2+1);
            vert_flat(2,icount1)=V(count1+1,count2+1);
            vert_flat(1,icount1+1)=U(count1+1-1,count2+1);
            vert_flat(2,icount1+1)=V(count1+1-1,count2+1);
            vert_flat(1,icount1+2)=U(count1+1,count2+1-1);
            vert_flat(2,icount1+2)=V(count1+1,count2+1-1);   
            tri_flat(1,icount2)=icount1;
            tri_flat(2,icount2)=icount1+1;
            tri_flat(3,icount2)=icount1+2;
            icount1=icount1+3;
            icount2=icount2+1;
        end
    end
    tri_flat=tri_flat';
%    trisurf(tri_flat',vert_flat(1,:),vert_flat(2,:),vert_flat(1,:)*0)
    
    vert_flat_tot=[];
    tri_flat_tot=[];
    icount=0;
%    size(tri_flat)
%    N^2
%    max(max(tri_flat))
%    return;
    for count1=1:n_tri
%        vert_flat_tot=[vert_flat_tot;vert_flat];
        tri_flat_tot=[tri_flat_tot;tri_flat+icount*(max(max(tri_flat)))];
        icount=icount+1;
    end
 U=vert_flat(1,:);
 V=vert_flat(2,:);
 U=U(:);
 V=V(:);
 
        P1=Points(Tri(:,1),:);
        P2=Points(Tri(:,2),:);
        P3=Points(Tri(:,3),:);
        P4=Points(Tri(:,4),:);
        P5=Points(Tri(:,5),:);
        P6=Points(Tri(:,6),:);
        B_x=[P1(:,1) P2(:,1) P3(:,1) P4(:,1) P5(:,1) P6(:,1)]';
        B_y=[P1(:,2) P2(:,2) P3(:,2) P4(:,2) P5(:,2) P6(:,2)]';
        B_z=[P1(:,3) P2(:,3) P3(:,3) P4(:,3) P5(:,3) P6(:,3)]';
        coef_x=Mm1*B_x;
        coef_y=Mm1*B_y;
        coef_z=Mm1*B_z;
        F_x=eval_pol_tri2(coef_x,U,V);
        F_y=eval_pol_tri2(coef_y,U,V);
        F_z=eval_pol_tri2(coef_z,U,V);
%        F_x=F_x';
%        F_y=F_y';
%        F_z=F_z';
        
%         plot3(F_x(:),F_y(:),F_z(:),'*')
%  surf(F_x,F_y,F_z,'edgealpha', 0.2,'Parent',ha)
 vert_flat=[F_x(:) F_y(:) F_z(:)];
 
% trisurf(tri_flat_tot(1:10,:),vert_flat(1,:),vert_flat(2,:),vert_flat(3,:))
col=rand(1,n_tri);
aux=ones(N^2*3,1);
COL=aux*col;
COL=COL(:);


    if strcmp(flag_rot,'z>z (default)')    
    elseif strcmp(flag_rot,'z>x')
        aux=F_x;
        F_x=F_z;
        F_z=-aux;
    elseif strcmp(flag_rot,'z>y')
        aux=F_y;
        F_y=F_z;
        F_z=-aux;
    end

    if strcmp(flag_col,'uniform')
        trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),COL*0,'Parent',ha)
    elseif strcmp(flag_col,'tri')
        trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),COL,'Parent',ha)
    elseif strcmp(flag_col,'z')
        trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),F_z(:),'Parent',ha)
    elseif strcmp(flag_col,'x')
        trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),F_x(:),'Parent',ha)
    elseif strcmp(flag_col,'y')
        trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),F_y(:),'Parent',ha)
    end

end


function F=eval_pol_tri2(coef,U,V)
    F=coef(1,:)+U*coef(2,:)+V*coef(3,:)+(U.^2)*coef(4,:)+(V.^2)*coef(5,:)+(U.*V)*coef(6,:);
end

function F=eval_pol_tri(coef,U,V)
    F=coef(1)+coef(2)*U+coef(3)*V+coef(4)*U.^2+coef(5)*V.^2+coef(6)*U.*V;
end


function [F_x,F_y,F_z,nP_x,nP_y,nP_z,dS]=eval_quadratic_patch(P1,P2,P3,P4,P5,P6,U,V)
    B_x=[P1(1) P2(1) P3(1) P4(1) P5(1) P6(1)]';
    B_y=[P1(2) P2(2) P3(2) P4(2) P5(2) P6(2)]';
    B_z=[P1(3) P2(3) P3(3) P4(3) P5(3) P6(3)]';
    coef_x=Mm1*B_x;
    coef_y=Mm1*B_y;
    coef_z=Mm1*B_z;
    F_x=coef_x(1)+coef_x(2)*U+coef_x(3)*V+coef_x(4)*U.^2+coef_x(5)*V.^2+coef_x(6)*U.*V;
    F_y=coef_y(1)+coef_y(2)*U+coef_y(3)*V+coef_y(4)*U.^2+coef_y(5)*V.^2+coef_y(6)*U.*V;
    F_z=coef_z(1)+coef_z(2)*U+coef_z(3)*V+coef_z(4)*U.^2+coef_z(5)*V.^2+coef_z(6)*U.*V;
    U_x=coef_x(2)+2*coef_x(4)*U+coef_x(6)*V;
    U_y=coef_y(2)+2*coef_y(4)*U+coef_y(6)*V;
    U_z=coef_z(2)+2*coef_z(4)*U+coef_z(6)*V;
    V_x=coef_x(3)+2*coef_x(5)*V+coef_x(6)*U;
    V_y=coef_y(3)+2*coef_y(5)*V+coef_y(6)*U;
    V_z=coef_z(3)+2*coef_z(5)*V+coef_z(6)*U;
    nP_x=U_y.*V_z-U_z.*V_y;
    nP_y=U_z.*V_x+U_x.*V_z;
    nP_z=U_x.*V_y-U_y.*V_x;
    dS=sqrt(nP_x.^2+nP_y.^2+nP_z.^2);
end