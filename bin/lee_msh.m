function [Geometry]=lee_msh(filename)
    fid = fopen(filename);

%   fid = fopen('prueba_quadratic.msh');


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
    
     Points=[X,Y,Z]
     Tri=[P1 P2 P3 P4 P5 P6]

     return;

     plot_quadratic_mesh(Points,Tri,n_tri,n_points)
    
    
    Geometry.npoints=n_points;
    Geometry.ntri=n_tri;
    Geometry.Points=Points;
    Geometry.Tri=Tri;

    
%    Geometry.npoints=F(1)
%    Geometry.ntri=F(2)
%    Geometry.Points=[F(3:3:Geometry.npoints*3+2).' F(4:3:Geometry.npoints*3+2).' F(5:3:Geometry.npoints*3+2).'];
%    Geometry.Tri=reshape(F(Geometry.npoints*3+3:end),3,Geometry.ntri)';
%    Geometry.S_smooth{1}=[];
%    Geometry.S_smooth{2}=[];
    
    
end



function plot_quadratic_mesh(Points,Tri,n_tri,n_points)
    Mm1=[    1     0     0     0     0     0;
    -3    -1     0     4     0     0;
    -3     0    -1     0     0     4;
     2     2     0    -4     0     0;
     2     0     2     0     0    -4;
     4     0     0    -4     4    -4];
    u=linspace(0,1,10);
    [U,V]=meshgrid(u);
    V=V.*(1-U);
%     surf(U,V,U*0)
%     B=[0 0 0 0 0 1]';
%     coef=Mm1*B;
%     F=eval_pol_tri(coef,U,V);
    
    for count=1:n_tri
        P1=Points(Tri(count,1),:);
        P2=Points(Tri(count,2),:);
        P3=Points(Tri(count,3),:);
        P4=Points(Tri(count,4),:);
        P5=Points(Tri(count,5),:);
        P6=Points(Tri(count,6),:);
        B_x=[P1(1) P2(1) P3(1) P4(1) P5(1) P6(1)]';
        B_y=[P1(2) P2(2) P3(2) P4(2) P5(2) P6(2)]';
        B_z=[P1(3) P2(3) P3(3) P4(3) P5(3) P6(3)]';
        coef_x=Mm1*B_x;
        coef_y=Mm1*B_y;
        coef_z=Mm1*B_z;
        F_x=eval_pol_tri(coef_x,U,V);
        F_y=eval_pol_tri(coef_y,U,V);
        F_z=eval_pol_tri(coef_z,U,V);
        surf(F_x,F_y,F_z,'edgealpha', 0.2)
        hold on
    end
%    shading interp
    axis equal
    grid
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