function [cam_R_obj, cam_t_obj] = pose_3D3D(object, pts_in_cam)

c_obj = mean(object, 2);              
c_cam = mean(pts_in_cam, 2);

H = (object - c_obj) * (pts_in_cam - c_cam)';
[U, ~, V] = svd(H);

cam_R_obj = V * diag([1, 1, det(V * U')]) * U';   % det = +1
cam_t_obj = c_cam - cam_R_obj * c_obj;
end