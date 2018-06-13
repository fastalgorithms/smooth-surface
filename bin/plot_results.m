function plot_results
%      filename1='Multiscale_1_n45_r0.gov'
%      filename2='Multiscale_1.msh';
%      x0=0;
%      y0=0;
%      z0=-.5;
%      [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz]=plot_gov(filename1,x0,y0,z0);
%      Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%      [Geometry]=lee_msh(filename2);
% 
% 
%     filename1='mi_barco_simple_13.gov'
%     filename2='mi_barco_simple_13.msh';
%     x0=0;
%     y0=.0;
%     z0=-.5;
%     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%     Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%     figure
%     plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%     title('45 Gaussian points per triangle')
%     camproj('perspective')
%     figure
%     plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
%     title('Normal, U and V vectors (zoom)')
%     camproj('perspective')
%    figure
%    [Geometry]=lee_msh(filename2);



     filename1='open_cavity_30deg_v2.gov'
     filename2='open_cavity_30deg_v2.msh';
     x0=1.5;
     y0=1.5;
     z0=-4;
     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
     Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
     figure
     plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
     title('45 Gaussian points per triangle')
     camproj('perspective')
     figure
     plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
     title('Normal, U and V vectors (zoom)')
     camproj('perspective')
     figure
    [Geometry]=lee_msh(filename2);



%     filename1='pico_1.gov'
%     filename2='pico_1.msh';
%     x0=0;
%     y0=0;
%     z0=1;
%     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz]=plot_gov(filename1,x0,y0,z0);
%     Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%     [Geometry]=lee_msh(filename2);
    
%     filename1='pico_2.gov'
%     filename2='pico_2.msh';
%     x0=0;
%     y0=0;
%     z0=2;
%     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz]=plot_gov(filename1,x0,y0,z0);
%     Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%     [Geometry]=lee_msh(filename2);

%     filename1='sci_fi_3.gov'
%     filename2='sci_fi_3.msh';
%     x0=0;
%     y0=0;
%     z0=1;
%     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%     length(X)/order
%     subplot(1,2,1)
%     plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%     title('45 Gaussian points per triangle')
%     camproj('perspective')
%     subplot(1,2,2)
%     plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
%     title('Normal, U and V vectors (zoom)')
%     camproj('perspective')
%     Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    figure
%   [Geometry]=lee_msh(filename2);

  
%     filename1='parabolic_antenna.gov'
%     filename2='parabolic_antenna.msh';
%     x0=0;
%     y0=0.75;
%     z0=0;
%     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%     length(X)/order
%     subplot(1,2,1)
%     plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%     title('45 Gaussian points per triangle')
%     camproj('perspective')
%     subplot(1,2,2)
%     plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
%     title('Normal, U and V vectors (zoom)')
%     camproj('perspective')
%Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)

%     figure
%     [Geometry]=lee_msh(filename2);


%     filename1='two_cavity_filter.gov'
%     filename2='two_cavity_filter.msh';
%     x0=0;
%     y0=0;
%     z0=-0.5;
%     [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%     length(X)/order
% figure
%     plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%     title('45 Gaussian points per triangle')
%     camproj('perspective')
%     figure
%     plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
%     title('Normal, U and V vectors (zoom)')
%     camproj('perspective')
%Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    figure
%    [Geometry]=lee_msh(filename2);

%    filename1='huge_genus_4.gov'
%    filename2='huge_genus_4.msh';
%    x0=0;
%    y0=0;
%    z0=-0.5;
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz]=plot_gov(filename1,x0,y0,z0);
%    Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    [Geometry]=lee_msh(filename2);


%{
    x0=0;
    y0=0;
    z0=0;
    filename1='Round_1_n45_r0.gov'
    filename2='Round_1_n45_r1.gov'
    filename3='Round_1_n45_r2.gov'
    filename_msh='Round_1.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,3,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')
    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors')
    camproj('perspective')
    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,3,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')
    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
    subplot(2,3,3)
    
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 2nd refinement')
    camproj('perspective')
    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2,Err_rel_45_3];
    
    
    
    filename1='Round_1_n78_r0.gov'
    filename2='Round_1_n78_r1.gov'
    filename3='Round_1_n78_r2.gov'
    filename_msh='Round_1.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2,Err_rel_78_3];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
    figure
    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')


%}







%{
    x0=0;
    y0=0;
    z0=0;
    filename1='Round_2_n45_r0.gov'
    filename2='Round_2_n45_r1.gov'
    filename3='Round_2_n45_r2.gov'
    filename_msh='Round_2.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,3,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')
    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors')
    camproj('perspective')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,3,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')
    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
    subplot(2,3,3)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 2nd refinement')
    camproj('perspective')
    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2,Err_rel_45_3];
    
    
    
    filename1='Round_2_n78_r0.gov'
    filename2='Round_2_n78_r1.gov'
    filename3='Round_2_n78_r2.gov'
    filename_msh='Round_2.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2,Err_rel_78_3];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
    figure
    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')


%    [Geometry]=lee_msh(filename_msh);

%}










%{
    x0=0;
    y0=0;
    z0=1;
    filename1='cube_substraction_n45_r0.gov'
    filename2='cube_substraction_n45_r1.gov'
    filename3='cube_substraction_n45_r2.gov'
    filename_msh='cube_substraction_2.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,3,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')

    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors')
    camproj('perspective')
    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,3,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')
    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
    subplot(2,3,3)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 2nd refinement')
    camproj('perspective')
    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2,Err_rel_45_3];
    
    
    
    filename1='cube_substraction_n78_r0.gov'
    filename2='cube_substraction_n78_r1.gov'
    filename3='cube_substraction_n78_r2.gov'
    filename_msh='cube_substraction.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2,Err_rel_78_3];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
    figure
    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')

%}









%{
    x0=0;
    y0=0;
    z0=1;
    filename1='Multiscale_1_n45_r0.gov'
    filename2='Multiscale_1_n45_r1.gov'
    filename3='Multiscale_1_n45_r2.gov'
    filename_msh='Multiscale_1.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,3,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')
    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors')
    camproj('perspective')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,3,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')

    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
    subplot(2,3,3)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 2nd refinement')
    camproj('perspective')

    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2,Err_rel_45_3];
    
    
    
    filename1='Multiscale_1_n78_r0.gov'
    filename2='Multiscale_1_n78_r1.gov'
    filename3='Multiscale_1_n78_r2.gov'
    filename_msh='Multiscale_1.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2, n_tri_3];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2,Err_rel_78_3];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
    figure
    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')

%}






%{
    x0=0;
    y0=0;
    z0=2;
    filename1='pico_2_n45_r0.gov'
    filename2='pico_2_n45_r1.gov'
%    filename3='pico_2_n45_r2.gov'
    filename_msh='pico_2.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,2,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);    
    title('45 Gaussian points per triangle')
    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,2,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    grid
%    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2];
    
    
    
    filename1='pico_2_n78_r0.gov'
    filename2='pico_2_n78_r1.gov'
    filename3='pico_2_n78_r2.gov'
    filename_msh='pico_2.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
    figure
    [Geometry]=lee_msh(filename_msh);

%}
















%{
    x0=0.5;
    y0=0.5;
    z0=0.5;
    filename1='high_genus_3_n45_r0.gov'
    filename2='high_genus_3_n45_r1.gov'
%    filename3='pico_2_n45_r2.gov'
    filename_msh='high_genus_3.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,2,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')

    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors (zoom)')
    camproj('perspective')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,2,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')


    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;

%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    grid
%    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2];
    
    
    
    filename1='high_genus_3_n78_r0.gov'
    filename2='high_genus_3_n78_r1.gov'
    filename3='high_genus_3_n78_r2.gov'
    filename_msh='high_genus_3.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
figure

    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')

%}













%{
    x0=0;
    y0=0.1;
    z0=0.2;
    filename1='capsule_multiscale_n45_r0.gov'
    filename2='capsule_multiscale_n45_r1.gov'
%    filename3='pico_2_n45_r2.gov'
    filename_msh='capsule_multiscale.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,2,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')

    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors (zoom)')
    camproj('perspective')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,2,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')


    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;

%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    grid
%    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2];
    
    
    
    filename1='capsule_multiscale_n78_r0.gov'
    filename2='capsule_multiscale_n78_r1.gov'
    filename_msh='capsule_multiscale.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
figure

    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')

%}














%{
    x0=0;
    y0=0;
    z0=-1.5;
    filename1='esfera_esfera_n45_r0.gov'
    filename2='esfera_esfera_n45_r1.gov'
%    filename3='pico_2_n45_r2.gov'
    filename_msh='esfera_esfera.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,2,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')

    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors (zoom)')
    camproj('perspective')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,2,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')


    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;

%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    grid
%    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2];
    
    
    
    filename1='esfera_esfera_n78_r0.gov'
    filename2='esfera_esfera_n78_r1.gov'
    filename_msh='esfera_esfera.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
figure

    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')

%}

















%{
    x0=-0.6;
    y0=-0.6;
    z0=-0.6;
    filename1='cubo_esfera_multires_n45_r0.gov'
    filename2='cubo_esfera_multires_n45_r1.gov'
%    filename3='pico_2_n45_r2.gov'
    filename_msh='cubo_esfera_multires.msh';
    figure
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
    subplot(2,2,1)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle')
    camproj('perspective')

    subplot(2,2,3)
    plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    title('Normal, U and V vectors (zoom)')
    camproj('perspective')

    Err_rel_45_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
    subplot(2,2,2)
    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    title('45 Gaussian points per triangle 1st refinement')
    camproj('perspective')


    Err_rel_45_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;

%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    grid
%    Err_rel_45_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_45=[Err_rel_45_1,Err_rel_45_2];
    
    
    
    filename1='cubo_esfera_multires_n78_r0.gov'
    filename2='cubo_esfera_multires_n78_r1.gov'
    filename_msh='cubo_esfera_multires.msh';

    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename1);
%    subplot(1,3,1)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_1=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_1=length(X)/order;
    
    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename2);
%    subplot(1,3,2)
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
    Err_rel_78_2=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
    n_tri_2=length(X)/order;
    
%    [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename3);
%    subplot(1,3,3)   
%    plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order);
%    Err_rel_78_3=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)
%    n_tri_3=length(X)/order
    
    n_tri_vect=[n_tri_1, n_tri_2];
    Err_rel_78=[Err_rel_78_1,Err_rel_78_2];

    
    subplot(2,2,4)
    loglog(n_tri_vect,Err_rel_45,n_tri_vect,Err_rel_78)
    xlabel('Number of triangles')
    ylabel('Error in Gauss integral identity')
    legend('45 points per triangle','78 points per triangle')
    grid
figure

    [Geometry]=lee_msh(filename_msh);
    camproj('perspective')

%}


end




function Err_rel=test_Gauss_integral(X,Y,Z,W,nSx,nSy,nSz,x0,y0,z0)


    R=sqrt((X-x0).^2+(Y-y0).^2+(Z-z0).^2);
    Ex=(X-x0)./(4*pi*R.^3);
    Ey=(Y-y0)./(4*pi*R.^3);
    Ez=(Z-z0)./(4*pi*R.^3);
%       Ex=Ex*0;
%       Ey=Ey*0;
%       Ez=Ez*0;
    Integral=sum((Ex.*nSx+Ey.*nSy+Ez.*nSz).*W);
    Err_rel=abs(Integral-1);

end


function [X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order]=load_gov(filename)
    fid = fopen(filename);
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

    
    X=PTS(:,1);
    Y=PTS(:,2);
    Z=PTS(:,3);
    
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
    
end

function plot_gov_dots(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    Ntri=length(X)/order;
    for count=0:(Ntri-1)
        X_local=X(count*order+1:(count+1)*order);
        Y_local=Y(count*order+1:(count+1)*order);
        Z_local=Z(count*order+1:(count+1)*order);
        plot3(X_local,Y_local,Z_local,'.','MarkerSize',2)
        hold on
    end
    axis equal
    hold off
    grid
end


function plot_gov_arrows(X,Y,Z,W,nSx,nSy,nSz,Ux,Uy,Uz,Vx,Vy,Vz,order)
    quiver3(X,Y,Z,nSx,nSy,nSz)
    hold on
    quiver3(X,Y,Z,Ux,Uy,Uz)
    quiver3(X,Y,Z,Vx,Vy,Vz)
    axis equal
    grid
    hold off
end

