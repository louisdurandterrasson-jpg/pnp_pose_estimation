function res = compute_err_persp_proj(x, pts_obj_ini, pts_im_ini)

persistent pts_obj pts_im_m N

if nargin>1
    % Initialization of the parameters that are kept constant during optimization
	pts_obj = pts_obj_ini;
	pts_im_m = pts_im_ini;
	N = size(pts_obj,2);
	res = 0;

else

    % Computation of the error for the current estimated pose x
	% Note that the error is here computed in the normlized image plane for
    % limiting computations
    % After initialization, pts_obj and pts_im_m contains the 3D coordinates of the 
    % object points and the image points in the normalized plane

    R = expm(skew(x(1:3)));        
    t = x(4:6);
    t = t(:);
    P = R * pts_obj + t;            
    p = P(1:2,:) ./ P(3,:);        
    e = p - pts_im_m(1:2,:);       
    res = e(:);                    
end