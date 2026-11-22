function [cam_R_obj_est, tcam_obj_est, Rnon_ortho] = Linear_PnP(object, pts_im, K)

nb_pts = size(object, 2);

% Linear solving method using translation and rotation matrix coefficients as unknowns

pts_norm = K \ [pts_im; ones(1, nb_pts)];

A = zeros(2*nb_pts, 12);
for i = 1:nb_pts
    P = [object(:, i)', 1];   
    x = pts_norm(1, i);
    y = pts_norm(2, i);
    A(2*i-1, :) = [P, zeros(1,4), -x*P];  
    A(2*i,   :) = [zeros(1,4), P, -y*P]; 
end

[~, ~, V] = svd(A);
M = reshape(V(:, end), 4, 3)';  

% Scaling solution so that det=1 (note that it does not impose orthogonality
% constraints)
scale = nthroot(det(M(:, 1:3)), 3);  

tcam_obj_est = M(:, 4) / scale;
Rnon_ortho = M(:, 1:3) / scale; % R before enforcing orthogonality constraints

% Enforcing orthogonality constraints
% Projection sur SO(3) : matrice de rotation la plus proche de Rnon_ortho
[U, ~, Vr] = svd(Rnon_ortho);

cam_R_obj_est = U * diag([1, 1, det(U * Vr')]) * Vr';

end