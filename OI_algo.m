function [cam_R_obj_est, cam_t_obj_est] = OI_algo(object, pts_im, K)

% Pose estimation using the iterative method "Orthogonal Iteration"
% proposed by Lu et al, 1998

n = size(object, 2);

v = K \ [pts_im; ones(1, n)];

V = zeros(3, 3, n);
for i = 1:n
    V(:,:,i) = (v(:,i) * v(:,i)') / (v(:,i)' * v(:,i));
end

Vsum = sum(V, 3);
F = inv(eye(3) - Vsum / n) / n;

p_mean = mean(object, 2);
p_c = object - p_mean;

q = v;                                
R = abs_orientation(p_c, q - mean(q, 2));

err_old = inf;
for it = 1:1000
    t = t_opt(R, object, V, F);

    q = zeros(3, n);
    err = 0;
    for i = 1:n
        Pc = R * object(:,i) + t;
        q(:,i) = V(:,:,i) * Pc;
        err = err + sum((Pc - q(:,i)).^2);  
    end

    R = abs_orientation(p_c, q - mean(q, 2));

    if abs(err_old - err) < 1e-12 * max(err_old, 1e-20) || err < 1e-20
        break;
    end
    err_old = err;
end

cam_R_obj_est = R;
cam_t_obj_est = t_opt(R, object, V, F);

end

% translation optimale
function t = t_opt(R, object, V, F)
n = size(object, 2);
s = zeros(3, 1);
for i = 1:n
    s = s + (V(:,:,i) - eye(3)) * R * object(:,i);
end
t = F * s;
end

% rotation qui minimise sum
function R = abs_orientation(p_c, q_c)
M = q_c * p_c';                          
[U, ~, W] = svd(M);
R = U * diag([1, 1, det(U * W')]) * W';  
end