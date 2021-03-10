function [ntri,norder,npoints,srcvals]=open_gov3(filename)
    fid = fopen(filename);
    F = fscanf(fid, '%g %g', [1 inf]);
    fclose(fid);
    norder=F(1);
    ntri=F(2);
    npols=(norder+1)*(norder+2)/2;
    npoints=npols*ntri;
    F=F(3:end);
    srcvals=zeros(12,npoints);
    srcvals(1,:)=F(1:npoints);
    srcvals(2,:)=F(npoints+1:2*npoints);
    srcvals(3,:)=F(2*npoints+1:3*npoints);
    srcvals(4,:)=F(3*npoints+1:4*npoints);
    srcvals(5,:)=F(4*npoints+1:5*npoints);
    srcvals(6,:)=F(5*npoints+1:6*npoints);
    srcvals(7,:)=F(6*npoints+1:7*npoints);
    srcvals(8,:)=F(7*npoints+1:8*npoints);
    srcvals(9,:)=F(8*npoints+1:9*npoints);
    srcvals(10,:)=F(9*npoints+1:10*npoints);
    srcvals(11,:)=F(10*npoints+1:11*npoints);
    srcvals(12,:)=F(11*npoints+1:12*npoints);
end

