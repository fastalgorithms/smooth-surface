function plot_go3(filename,ha,N,flag_col,flag_rot)
[ntri,norder,npoints,srcvals]=open_gov3(filename);
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
[a,b]=size(tri_flat);
tri_flat_tot=zeros(a*ntri,b);
dd=max(max(tri_flat));
    for count1=1:ntri
%        vert_flat_tot=[vert_flat_tot;vert_flat];
         tri_flat_tot(a*(count1-1)+1:a*count1,:)=tri_flat+(count1-1)*dd;
%        tri_flat_tot=[tri_flat_tot;tri_flat+icount*(max(max(tri_flat)))];
%        icount=icount+1;
    end


%     for count1=1:ntri
% %        vert_flat_tot=[vert_flat_tot;vert_flat];
%         tri_flat_tot=[tri_flat_tot;tri_flat+icount*(max(max(tri_flat)))];
%         icount=icount+1;
%     end
    U=vert_flat(1,:);
    V=vert_flat(2,:);
    U=U(:);
    V=V(:);
    [uvs,wts]=get_vioreanu_nodes(norder);
    npols=(norder+1)*(norder+2)/2;
    M1=koorn_vals2coefs(norder, npols, uvs);
    
    UVS=[U'; V'];
    sizeUVS=size(UVS)
    M2=koorn_coefs2vals(norder, npols, UVS);
    
    MM=M2*M1;
    
%     
%     x=uvs(1,:);
%     y=uvs(2,:);
%     z=x.*y;
%     plot3(x,y,z,'*')
%     
%     
%     u=linspace(0,1,N+1);
%     [U,V]=meshgrid(u);  
% 
%     index=find(U+V<=1);
%     U=U(index);
%     V=V(index);
%     U=U(:);
%     V=V(:);
%     Z=U.*V;
%     figure
%     plot3(U,V,Z,'*')
%     figure
%     size(z)
%     Z2=MM*(z');
%     plot3(U,V,Z2,'*')
%     
%     return
    
    Px=srcvals(1,:)';
    Py=srcvals(2,:)';
    Pz=srcvals(3,:)';
    Bx=reshape(Px,npols,ntri);
    By=reshape(Py,npols,ntri);
    Bz=reshape(Pz,npols,ntri);

    PPx=MM*Bx;
    PPy=MM*By;
    PPz=MM*Bz;
    size(PPx(:))
%    tri_flat_tot
%    vert_flat_tot=[PPx(:),PPy(:),PPz(:)];

%    trisurf(tri_flat_tot,PPx(:),PPy(:),PPz(:),'Parent',ha)
%    axis equal
%    shading interp
%    camlight headlight
    

    F_x=PPx;
    F_y=PPy;
    F_z=PPz;
    
    
    col=rand(1,ntri);
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


function [uvs, umatr, vmatr, wts]=vioreanu_simplex_quad(norder, npols)
% !
% !  This subroutine extracts the precomputed vioreanu quadrature
% !  or discretization nodes, the corresponding weights, the coeffs
% !  to values matrix, and the values to coeffs matrix
% !  
% !  Input arguments
% !    
% !    - norder: integer
% !        input order. Must be between 0 and 20
% !    - npols: integer
% !        npols = (norder+1)*(norder+2)/2
% !
% !  Output arguments
% !
% !    - uvs: double precision(2,npols)
% !        vioreanu quadrature nodes of order=norder+1
% !    - umatr: double precision (npols,npols)
% !        vals to coeffs matrix
% !    - vmatr: double precision (npols,npols)
% !        coefs to vals matrix
% !    - wts: double precision(npols)
% !        vioreanu quadrature weights of order=norder+1
% !
  [uvs,wts]=get_vioreanu_nodes(norder);
  [umatr,vmatr]=koorn_vals2coefs_coefs2vals(norder,npols);


end




function [umatr,vmatr]=koorn_vals2coefs_coefs2vals(korder,kpols)
%   !
%   !!   compute vals2coefs and coefs2vals matrix 
%   !
%   !    input
%   !    korder     in: integer
%   !                 order of rokhlin vioreanu (rv) nodes at which function is
%   !                 sampled
%   !    kpols      in: integer
%   !                 number of nodes corresponding to korder rv nodes
%   !
%   !    output
%   !    umatr      out: real *8 (kpols,kpols)
%   !               vals2coefs matrix
%   ! 
%   !    vmatr      out: real *8 (kpols,kpols)
%   !               coefs2vals matrix
  [xys,wts]=get_vioreanu_nodes(korder);
  vmatr=koorn_coefs2vals(korder,kpols,xys);
  umatr=inv(vmatr);
end
  
function amat=koorn_vals2coefs(nmax, npols, uvs)
%   !
%   ! This routine returns a square matrix that maps point values of a
%   ! function to coefficients in an orthonormal polynomial expansion on
%   ! the simplex (0,0), (1,0), (0,1). The point locations can be
%   ! arbitrary, but note that they WILL affect the conditioning of this
%   ! matrix.
%   !
%   ! Input:
%   !   nmax, npols - order of expansion, it should be the case that
%   !       npols = (nmax+1)(nmax+2)/2
%   !   uvs - point locations, npols of them
%   !
%   ! Output:
%   !   amat - matrix such that coefs = amat*vals
%   !
%   !

  bmat=koorn_coefs2vals(nmax, npols, uvs);
  amat=inv(bmat);
end

function amat=koorn_coefs2vals(nmax, npols, uvs)
%   !
%   ! This routine returns a square matrix that maps coefficients in an
%   ! orthonormal polynomial expansion on the simplex (0,0), (1,0),
%   ! (0,1) to function values at the points uvs. The point locations
%   ! can be arbitrary, but note that they WILL affect the conditioning
%   ! of this matrix.
%   !
%   ! Input:
%   !   nmax, npols - order of expansion, it should be the case that
%   !       npols = (nmax+1)(nmax+2)/2
%   !   uvs - point locations, npols of them
%   !
%   ! Output:
%   !   amat - matrix such that vals = amat*coefs
%   !
%   !
[a,b]=size(uvs);

  for i=1:b
      [npols2, pols]=koorn_pols(uvs(:,i), nmax);
%    call koorn_pols(uvs(1,i), nmax, npols2, pols)
    for j=1:npols
      amat(i,j) = pols(j);
    end
  end
end
  


function [npols, pols]=koorn_pols(uv, nmax)

%   !
%   ! This subroutine evalutes a bunch of orthogonal polynomials on the
%   ! simplex with vertices (0,0), (1,0), (0,1). The polynomials
%   ! computed by this routine are normalized and rotated Koornwinder
%   ! polynomials, which are classically defined on the triangle with
%   ! vertices (0,0), (1,0), (1,1), and given analytically as:
%   !
%   !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
%   !
%   ! After mapping to the uv-simplex via the change of variables
%   !
%   !      u = y,  v = 1-x
%   !
%   ! these polynomials are given as
%   !
%   !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
%   !        P_k((2u+v-1)/(1-v))
%   !
%   ! See Koornwinder 1975 for details, or the NIST handbook.
%   !
%   ! Input:
%   !   uv - a point in the simplex (0,0), (1,0), (0,1)
%   !   nmax - maximum degree polynomials to compute, a total of
%   !       (nmax+1)(nmax+2)/2 values are returned
%   !
%   ! Output:
%   !   npols - number of pols = (nmax+1)*(nmax+2)/2
%   !   pols - values of all the polynomials, ordered in the following
%   !       (n,k) manner:
%   !
%   !            (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), etc.
%   !
%   !
% 
%  implicit real *8 (a-h,o-z)
%  real *8 :: uv(2), pols(*)

%  real *8 :: legpols(0:100), jacpols(0:100,0:100)
legpols=zeros(1,101);
jacpols=zeros(101,101);

  done = 1;
  


  u = uv(1);
  v = uv(2);
  z = 2*u+v-1;
  y = 1-v;

  legpols(1+0) = 1;
  legpols(1+1) = z;

  for k=1:nmax
    legpols(1+k+1) = ((2*k+1)*z*legpols(1+k) - k*legpols(1+k-1)*y*y)/(k+1);
  end

  x = 1-2*v;
  
  for k = 0:nmax
    beta = 2*k+1    ;
    jacpols(1+0,1+k) = 1;
    jacpols(1+1,1+k) = (-beta + (2+beta)*x)/2;
    
    for n = 1:nmax-k-1
      an = (2*n+beta+1)*(2*n+beta+2)/2/(n+1)/(n+beta+1);
      bn = (-beta^2)*(2*n+beta+1)/2/(n+1)/(n+beta+1)/(2*n+beta);
      cn = n*(n+beta)*(2*n+beta+2)/(n+1)/(n+beta+1)/(2*n+beta);
      jacpols(1+n+1,1+k) = (an*x + bn)*jacpols(1+n,1+k) - cn*jacpols(1+n-1,1+k);
    end

  end



  iii = 0;
  for n = 0:nmax
    for k = 0:n
      sc = sqrt(done/(2*k+1)/(2*n+2));
      iii = iii + 1;
      pols(iii) = legpols(1+k)*jacpols(1+n-k,1+k)/sc;
    end
  end

  npols = iii;
  
end 




function val=koorn_evalexp(nmax, npols, uv, coefs)

%   !
%   ! Evaluate the orthgonal polynomial expansion with given
%   ! coefficients at the points uv
%   !
%   ! Input:
%   !   nmax, npols - number of terms in expansion,
%   !       npols = (nmax+1)*(nmax+2)/2
%   !   uv - point in the simplex at which to evaluate
%   !   coefs - coefs of expansion, ordered the same as koorn_pols
%   !
%   ! Output:
%   !   val - the value of the expansion at uv
%   !
%   !

  [npols2, pols]=koorn_pols(uv, nmax);
  val = 0;
  for i = 1:npols
    val = val + coefs(i)*pols(i);
  end

end 












