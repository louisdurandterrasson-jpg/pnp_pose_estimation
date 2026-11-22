function res = compute_err_rep_in_pix(Rcam_obj_est, tcam_obj_est, object, pts_im, K)

% Question 2.2
nb_pts = size(object, 2);
err_sum = 0;
for i = 1:nb_pts
    P_cam = Rcam_obj_est * object(:, i) + tcam_obj_est;
    p_proj = K * (P_cam / P_cam(3));
    err_sum = err_sum + norm(p_proj(1:2) - pts_im(:, i));
end
res = err_sum / nb_pts;

end
        

