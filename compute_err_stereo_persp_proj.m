function res = compute_err_stereo_persp_proj(x, pts_obj_ini, pts_cam_ini, pts_cam2_ini, cam1_T_cam2_ini)

persistent pts_obj pts_cam pts_cam2 N cam2_R_cam1 cam2_t_cam1

if nargin>1
	pts_obj = pts_obj_ini;
	pts_cam = pts_cam_ini;
	pts_cam2 = pts_cam2_ini;
	N = size(pts_obj_ini,2);
	res = 0;
    cam2_R_cam1 = (cam1_T_cam2_ini(1:3,1:3))';
    cam2_t_cam1 = -cam2_R_cam1 * cam1_T_cam2_ini(1:3,4);

else

    R = expm(skew(x(1:3)));
    t = x(4:6);
    t = t(:);

    P1 = R * pts_obj + t;                   
    P2 = cam2_R_cam1 * P1 + cam2_t_cam1;     
    e1 = P1(1:2,:) ./ P1(3,:) - pts_cam(1:2,:);
    e2 = P2(1:2,:) ./ P2(3,:) - pts_cam2(1:2,:);
	res = [e1(:); e2(:)];                  
end