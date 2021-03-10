function [Points,Tri,n_tri]=open_msh(filename)
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
%    plot_quadratic_mesh(Points,Tri,n_tri,ha,N,flag_col,flag_rot)
end