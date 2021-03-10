function refresh_geometry_cell(tri_flat_tot,F_x,F_y,F_z,COL,Err_L2,Err_L1,Err_Lsup,ha,flag_col,flag_rot,flag_err, err_act,alph)
    n_cell=length(tri_flat_tot);
    for count=1:n_cell
       refresh_geometry(tri_flat_tot{count},F_x{count},F_y{count},F_z{count},COL{count},Err_L2{count},Err_L1{count},Err_Lsup{count},ha,flag_col,flag_rot,flag_err,err_act,alph{count})
       hold on
    end
    hold off
end
