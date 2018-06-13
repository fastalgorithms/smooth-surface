function [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz]=plot_gov(filename,x0,y0,z0)

%x0,y0,z0 Must be a point inside the surface to test the Gauss theorem

%    filename='fortran_plot'
    fid = fopen(filename);

%   fid = fopen('prueba_quadratic.msh');


    F = fscanf(fid, '%g %g', [1 inf]);
    fclose(fid);
    n_Sf_points=F(3);
    order=F(1);
    F=F(4:end);
    [a,b]=size(F);
    3*n_Sf_points*4+n_Sf_points;
    
    pts=F(1:3*n_Sf_points);
    u=F(3*n_Sf_points+1:3*n_Sf_points*2);
    v=F(3*n_Sf_points*2+1:3*n_Sf_points*3);
    n=F(3*n_Sf_points*3+1:3*n_Sf_points*4);
    W=F(3*n_Sf_points*4+1:end);
    PTS=[pts(1:end/3)',pts(end/3+1:end/3*2)',pts(end/3*2+1:end)'];
    U=[u(1:end/3)',u(end/3+1:end/3*2)',u(end/3*2+1:end)'];
    V=[v(1:end/3)',v(end/3+1:end/3*2)',v(end/3*2+1:end)'];
    N=[n(1:end/3)',n(end/3+1:end/3*2)',n(end/3*2+1:end)'];

    
%    plot3(PTS(:,1),PTS(:,2),PTS(:,3),'.')
    
    axis equal
    
    X=PTS(:,1);
    Y=PTS(:,2);
    Z=PTS(:,3);
    
%     R=sqrt((X-x0).^2+(Y-y0).^2+(Z-z0).^2);
%     Ex=(X-x0)./(4*pi*R.^3);
%     Ey=(Y-y0)./(4*pi*R.^3);
%     Ez=(Z-z0)./(4*pi*R.^3);
%     Integral=sum((Ex.*N(:,1)+Ey.*N(:,2)+Ez.*N(:,3)).*W')
%     Err_rel=abs(Integral-1)

    
    nSx=N(:,1)';
    nSy=N(:,2)';
    nSz=N(:,3)';
    Ux=U(:,1)';
    Uy=U(:,2)';
    Uz=U(:,3)';
    Vx=V(:,1)';
    Vy=V(:,2)';
    Vz=V(:,3)';
    X=X';
    Y=Y';
    Z=Z';
    
    
    
        grid
    axis equal
    plot_dots_colors(X,Y,Z,order)
    figure
%    plot3(X,Y,Z,'.')

    hold on
    quiver3(X,Y,Z,nSx,nSy,nSz)
    quiver3(X,Y,Z,Ux,Uy,Uz)
    quiver3(X,Y,Z,Vx,Vy,Vz)
    axis equal
    grid

    
    
end



function plot_dots_colors(X,Y,Z,order)
    Ntri=length(X)/order;
    figure
    for count=0:(Ntri-1)
        X_local=X(count*order+1:(count+1)*order);
        Y_local=Y(count*order+1:(count+1)*order);
        Z_local=Z(count*order+1:(count+1)*order);
        plot3(X_local,Y_local,Z_local,'.','MarkerSize',10)
        hold on
    end
    axis equal
    hold off
end