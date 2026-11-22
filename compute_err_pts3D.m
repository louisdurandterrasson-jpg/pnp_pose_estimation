function err_pts_3D = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam)

% To complete - Question 2.1
% Compute the 3D reconstruction error on 3D points expressed in the camera
% frame
pts_est = cam_R_obj_est * object + cam_t_obj_est;  
err_pts_3D = mean(sqrt(sum((pts_est - list_pts_in_cam).^2, 1)));  
end


