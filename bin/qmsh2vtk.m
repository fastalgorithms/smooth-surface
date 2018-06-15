function qmsh2vtk(filename)

  % converts a quadratic msh, as output by GiD, into a vtk file
  % for plotting in paraview, for example

  fid = fopen(filename);
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
  Points = Points';
  
  Tri=[P1 P2 P3 P4 P5 P6];
  Tri = Tri - 1;
  Tri = Tri';
  
  fileout = strrep(filename, '.msh', '.vtk')
  fid = fopen(fileout, 'wt')
  
  fprintf(fid, '# vtk DataFile Version 3.0\n');
  fprintf(fid, 'MSH file from GiD\n');
  fprintf(fid, 'ASCII\n');
  fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
  fprintf(fid, 'POINTS %d float\n', n_points);
  
  fprintf(fid, '%12.6e %12.6e %12.6e\n', Points);
  
  fprintf(fid, '\n');
  fprintf(fid, 'CELLS %d %d\n', n_tri, 7*n_tri);
  
  fprintf(fid, '6 %d %d %d %d %d %d\n', Tri);
  
  fprintf(fid, '\n');
  fprintf(fid, 'CELL_TYPES %d\n', n_tri);
  
  for i=1:n_tri
      fprintf(fid, '22\n');
  end

  fprintf(fid, '\n');
  fprintf(fid, 'POINT_DATA %d\n', n_points);
  fprintf(fid, 'SCALARS height float 1\n');
  fprintf(fid, 'LOOKUP_TABLE default\n');
  fprintf(fid, '%12.6e\n', Z);
  
  fclose(fid)
  
  return;
  
  
end
