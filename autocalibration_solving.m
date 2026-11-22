function [cam_t_obj_est, cam_R_obj_est, K_est] = autocalibration_solving(object, pts_im)

% Solve for pose and intrinsic parameters in persp. proj. from a single image
%   object : 3 x N (repere objet), pts_im : 2 x N (pixels)

N = size(object, 2);

% première étape : solution linéaire
A = zeros(2*N, 12);
for i = 1:N
    P = [object(:, i)', 1];                  
    u = pts_im(1, i);
    v = pts_im(2, i);
    A(2*i-1, :) = [P, zeros(1,4), -u*P];      
    A(2*i,   :) = [zeros(1,4), P, -v*P];   
end
[~, ~, V] = svd(A);
M = reshape(V(:, end), 4, 3)';               


if mean(M(3, :) * [object; ones(1, N)]) < 0
    M = -M;
end

% deuxième étape:decomposition RQ 
B  = M(:, 1:3);
Pf = flipud(eye(3));                          
[Q, U] = qr((Pf * B)');
K_raw = Pf * U' * Pf;                        
R     = Pf * Q';                             

D = diag(sign(diag(K_raw)));                  
K_raw = K_raw * D;
R     = D * R;

[Ur, ~, Vr] = svd(R);
R = Ur * diag([1, 1, det(Ur * Vr')]) * Vr';

t     = K_raw \ M(:, 4);                     
K_est = K_raw / K_raw(3,3);                  

%  troisième étape R
w0 = r2w(R);
x0 = [K_est(1,1); K_est(2,2); K_est(1,3); K_est(2,3); w0; t];

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display', 'off');
x = lsqnonlin(@(p) err_autocalib(p, object, pts_im), x0, [], [], options);

K_est         = [x(1), 0, x(3); 0, x(2), x(4); 0, 0, 1];
cam_R_obj_est = expm(skew(x(5:7)));
cam_t_obj_est = x(8:10);

end

% fonction minimisé par skew
function e = err_autocalib(p, object, pts_im)
K = [p(1), 0, p(3); 0, p(2), p(4); 0, 0, 1];
R = expm(skew(p(5:7)));
t = p(8:10);
t = t(:);
P = R * object + t;
q = K * (P ./ P(3,:));
e = q(1:2,:) - pts_im;
e = e(:);
end

% pour convertir la matrice de rotation en vecteur
function w = r2w(R)
c = max(min((trace(R) - 1) / 2, 1), -1);
theta = acos(c);
if theta < 1e-10
    w = [0; 0; 0];
else
    w = theta / (2*sin(theta)) * [R(3,2)-R(2,3); R(1,3)-R(3,1); R(2,1)-R(1,2)];
end
end