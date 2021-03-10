function [tri_flat_tot,F_x,F_y,F_z,COL,Err_L2,Err_L1,Err_Lsup]=load_flattri_smooth_geometry_cell(ntri,norder,npoints,srcvals,N)
    n_cell=length(ntri);
    for count=1:n_cell
        [tri_flat_tot{count},F_x{count},F_y{count},F_z{count},COL{count},Err_L2{count},Err_L1{count},Err_Lsup{count}]=load_flattri_smooth_geometry(ntri{count},norder{count},npoints{count},srcvals{count},N)
    end
end