function refresh_geometry(tri_flat_tot,F_x,F_y,F_z,COL,Err_L2,Err_L1,Err_Lsup,ha,flag_col,flag_rot,flag_err, err_act,alph)

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

    if err_act
        if strcmp(flag_err,'L2 Tail')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),Err_L2,'Parent',ha,'FaceAlpha',alph)
        elseif strcmp(flag_err,'L1 Tail')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),Err_L1,'Parent',ha,'FaceAlpha',alph)
        elseif strcmp(flag_err,'Lsup Tail')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),Err_Lsup,'Parent',ha,'FaceAlpha',alph)
        end
        colorbar(ha);
    else
        if strcmp(flag_col,'uniform')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),COL*0,'Parent',ha,'FaceAlpha',alph)
        elseif strcmp(flag_col,'tri')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),COL,'Parent',ha,'FaceAlpha',alph)
        elseif strcmp(flag_col,'z')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),F_z(:),'Parent',ha,'FaceAlpha',alph)
        elseif strcmp(flag_col,'x')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),F_x(:),'Parent',ha,'FaceAlpha',alph)
        elseif strcmp(flag_col,'y')
            hhh=trisurf(tri_flat_tot,F_x(:),F_y(:),F_z(:),F_y(:),'Parent',ha,'FaceAlpha',alph)
        end
    end
end