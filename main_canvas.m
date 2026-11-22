% Canvas for pose estimation tutorial

clear all; close all;

addpath EPnP;

% Loading data
load(fullfile('data', 'data_tutorial_PnP.mat'));
pause;

nb_config_test = numel(list_T);
nb_pts = size(object,2);

trans_err = [];
rot_err = [];
err_pts_3D = [];

list_tested_methods = [1,2,3,4,5,6,7]; % list that can be modified depending on the methods you want to test
% 1 : analytic PnP method (provided)
% 2 : linear method BEFORE orthogonality constraints enforcement
% 3 : linear method AFTER orthogonality constraints enforcement
% 4 : numerical optimization
% 5 : iterative algorithm (OI)
% 6 : autocalibration
% 7 : stereoscopy (triangulation (+ optimization))

noise_ampl = [0, 0.1, 0.5, 1, 2]; % different levels of image noise amplitudes
nb_noise = length(noise_ampl);

% Question 1.1
% Compute perfect image points (without noise) for each configurationto test
for kk=1:nb_config_test
    
    % To complete
    pts_proj = K * (list_pts_in_cam{kk} ./ list_pts_in_cam{kk}(3,:));
    list_pts_im_perf{kk} = pts_proj(1:2,:);
end

% Question 1.2
% To complete
for ii=1:nb_noise 
	
	for kk=1:nb_config_test
		
		ang_rand = 2*pi*rand(1, nb_pts);
		pts_im_noisy = list_pts_im_perf{kk} + noise_ampl(ii) * [cos(ang_rand); sin(ang_rand)];
        
        mat_pts_im{ii}{kk} = pts_im_noisy;
        
		
    end	
end

% Question 2.3
for ii=1:nb_noise 	

    tic
    for kk=1:nb_config_test
    
        % Exact values (= ground truth) for the tested configuration
        cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};	
        
        pts_im = mat_pts_im{ii}{kk};
        
		% Results from PnP method
        [cam_R_obj_est,cam_t_obj_est]=efficient_pnp([object', ones(nb_pts,1)], [pts_im', ones(nb_pts,1)],K);
		%[cam_R_obj_est,cam_t_obj_est] =efficient_pnp_gauss([object', ones(nb_pts,1)], [pts_im', ones(nb_pts,1)],K)
       
        
        pos_err{1}(ii,kk,1) = cam_t_obj_ref(1) - cam_t_obj_est(1);
		pos_err{1}(ii,kk,2) = cam_t_obj_ref(2) - cam_t_obj_est(2);
		pos_err{1}(ii,kk,3) = cam_t_obj_ref(3) - cam_t_obj_est(3);
		trans_err{1}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
		
        [theta,u] = r2thetau(cam_R_obj_ref* cam_R_obj_est');
		rot_err{1}(ii,kk) = abs(theta);
		
        err_pts_3D{1}(ii,kk) = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam{kk});
        err_rep{1}(ii,kk) = compute_err_rep_in_pix(cam_R_obj_est, cam_t_obj_est, object, pts_im, K);
        	
    end
    
     analytic_duration(ii) = toc;
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), err_pts_3D{1}(ii,list_ind_no_rot))
     strtitle = sprintf('Position errors of 3D points for analytical method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), trans_err{1}(ii,list_ind_no_rot))
     strtitle = sprintf('translation error for analytical method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), err_rep{1}(ii,list_ind_no_rot))
     strtitle = sprintf('image reprojection error for analytical method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
      
     % Comparison of translation errors along individual axes
     figure;
     plot(pos_err{1}(ii,:,1), 'r');
     hold on;
     plot(pos_err{1}(ii,:,2), 'g');
     plot(pos_err{1}(ii,:,3), 'b');
     legend('x', 'y', 'z');
     strtitle = sprintf('translation errors for analytical method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
     
end
analytic_duration 

% Statistics (median and max errors) per noise level
for jj=1:nb_noise	
	trans_err_cons = trans_err{1}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{1}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{1}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{1}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{1}(jj) = median(trans_err_filt,2); 
	%rot_err_med{1}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{1}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{ii}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{1}(jj) = median(points_err_filt,2);
    points_err_max{1}(jj) = max(points_err_filt,[],2);
        
    im_err_med{1}(jj) = median(im_err_filt,2);
    im_err_max{1}(jj) = max(im_err_filt,[],2);
        
end

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %% Question 3
 for ii=1:nb_noise 
   tic
   for kk=1:nb_config_test
        cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};
	
        pts_im = mat_pts_im{ii}{kk};
       
        % Resolution lineaire sur R, t 
        [cam_R_obj_est, cam_t_obj_est, Rnon_ortho] = Linear_PnP(object, pts_im, K);
    	
        trans_err{2}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
    	trans_err{3}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
    	
        [theta,u] = r2thetau(cam_R_obj_ref* cam_R_obj_est');
    	rot_err{2}(ii,kk) = abs(theta);
       
        err_rep{2}(ii,kk) = compute_err_rep_in_pix(cam_R_obj_est, cam_t_obj_est, object, pts_im, K);
        err_rep{3}(ii,kk) = compute_err_rep_in_pix(Rnon_ortho, cam_t_obj_est, object, pts_im, K);
       
        err_pts_3D{2}(ii,kk) = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam{kk});
        
        err_pts_3D{3}(ii,kk) = compute_err_pts3D(Rnon_ortho, cam_t_obj_est, object, list_pts_in_cam{kk});
       
        %pause;
    end
  	linear_duration(ii) = toc;
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), err_pts_3D{2}(ii,list_ind_no_rot))
    strtitle = sprintf('Position errors of 3D points for linear method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), trans_err{2}(ii,list_ind_no_rot))
    strtitle = sprintf('translation error for linear method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), err_rep{2}(ii,list_ind_no_rot))
    strtitle = sprintf('image reprojection error for linear method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
     
   
 end
 
 linear_duration
 
% Statistics (median and mex errors) per noise level

for jj=1:nb_noise	
	trans_err_cons = trans_err{2}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{2}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{2}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{2}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{2}(jj) = median(trans_err_filt,2); 
	%rot_err_med{2}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{2}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{2}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{2}(jj) = median(points_err_filt,2);
    points_err_max{2}(jj) = max(points_err_filt,[],2);
        
    im_err_med{2}(jj) = median(im_err_filt,2);
    im_err_max{2}(jj) = max(im_err_filt,[],2);
        
end

for jj=1:nb_noise	
	trans_err_cons = trans_err{3}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{2}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{3}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{3}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{3}(jj) = median(trans_err_filt,2); 
	%rot_err_med{2}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{3}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{2}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{3}(jj) = median(points_err_filt,2);
    points_err_max{3}(jj) = max(points_err_filt,[],2);
        
    im_err_med{3}(jj) = median(im_err_filt,2);
    im_err_max{3}(jj) = max(im_err_filt,[],2);
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 4
for ii=1:nb_noise   
    tic
    for kk=1:nb_config_test
        cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};
	
        % To complete : compute image points in the normalized image plane
        pts_im = mat_pts_im{ii}{kk};
        pts_normed_im_plane = K \ [pts_im; ones(1, nb_pts)];
        
        % Numerical minimization
        % Initialization of objective function to minimize
		compute_err_persp_proj(0, object, pts_normed_im_plane);
	
        % Initialization of pose
        % to complete
        x_ini = [0; 0; 0; 0; 0; 60];

        fprintf('initial error\n');
		compute_err_persp_proj(x_ini)
		
        % Optimization using fminunc or lsqnonlin function
        % to complete
        % x = fminunc(@compute_err_persp_proj, x_ini);
	   
        % or
        options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt')
        x = lsqnonlin(@compute_err_persp_proj, x_ini, [], [], options);
		
    
        % Extraction of cam_R_obj_est and cam_t_obj
        cam_R_obj_est = expm(skew(x(1:3)));
		cam_t_obj_est = x(4:6);
    
        trans_err{4}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
		err_rep{4}(ii,kk) = compute_err_rep_in_pix(cam_R_obj_est, cam_t_obj_est, object, pts_im, K)
     
        err_pts_3D{4}(ii,kk) = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam{kk})
       
        
    
    optim_duration(ii) = toc;

    end
   %pause;
 
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), err_pts_3D{4}(ii,list_ind_no_rot))
     strtitle = sprintf('Position errors of 3D points for optimization method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), trans_err{4}(ii,list_ind_no_rot))
     strtitle = sprintf('translation error for optimization method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), err_rep{4}(ii,list_ind_no_rot))
     strtitle = sprintf('image reprojection error for optimization method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
      
 end
optim_duration

% Statistics (median and mex errors) per noise level

for jj=1:nb_noise	
	trans_err_cons = trans_err{4}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{4}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{4}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{4}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{4}(jj) = median(trans_err_filt,2); 
	%rot_err_med{4}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{4}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{4}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{4}(jj) = median(points_err_filt,2);
    points_err_max{4}(jj) = max(points_err_filt,[],2);
        
    im_err_med{4}(jj) = median(im_err_filt,2);
    im_err_max{4}(jj) = max(im_err_filt,[],2);
        
end
     
%%%%%%%%%%%%%%%%%%%%%
%% Question 5
%% Iterative method
for ii=1:nb_noise   
    tic
    for kk=1:nb_config_test
        cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};
	
        %%%%%%%%
        % Try implementing the Orthogonal Iteration method
		pts_im = mat_pts_im{ii}{kk};   
        [cam_R_obj_est, cam_t_obj_est] = OI_algo(object, pts_im, K);
    
        trans_err{5}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
		err_rep{5}(ii,kk) = compute_err_rep_in_pix(cam_R_obj_est, cam_t_obj_est, object, pts_im, K)
     
        err_pts_3D{5}(ii,kk) = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam{kk})
       
    
        OI_duration(ii) = toc;

    end
   %pause;
 
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), err_pts_3D{5}(ii,list_ind_no_rot))
    strtitle = sprintf('Position errors of 3D points for OI method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), trans_err{5}(ii,list_ind_no_rot))
    strtitle = sprintf('translation error for OI method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), err_rep{5}(ii,list_ind_no_rot))
    strtitle = sprintf('image reprojection error for OI method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
      
end
OI_duration

% Statistics (median and mex errors) per noise level

for jj=1:nb_noise	
	trans_err_cons = trans_err{5}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{5}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{5}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{5}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{5}(jj) = median(trans_err_filt,2); 
	%rot_err_med{5}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{5}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{5}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{5}(jj) = median(points_err_filt,2);
    points_err_max{5}(jj) = max(points_err_filt,[],2);
        
    im_err_med{5}(jj) = median(im_err_filt,2);
    im_err_max{5}(jj) = max(im_err_filt,[],2);
        
end
     

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
 %% Question 6 -
 %% Auto-calibration
 
for ii=1:nb_noise 
   tic
   for kk=1:nb_config_test
        cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};
	
        pts_im = mat_pts_im{ii}{kk};
       
        % Solving with autocalibration
        [cam_t_obj_est, cam_R_obj_est, K_est] = autocalibration_solving(object, pts_im);

        if kk == list_ind_no_rot(25)
            list_K_est{ii} = K_est;
            t_z_est(ii) = cam_t_obj_est(3);
        end
           	
        trans_err{6}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
    	
        [theta,u] = r2thetau(cam_R_obj_ref*cam_R_obj_est');
    	rot_err{6}(ii,kk) = abs(theta);
       
        err_rep{6}(ii,kk) = compute_err_rep_in_pix(cam_R_obj_est, cam_t_obj_est, object, pts_im, K);
        
        err_pts_3D{6}(ii,kk) = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam{kk});
        
       
    end
  	autocalib_duration(ii) = toc;
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), err_pts_3D{6}(ii,list_ind_no_rot))
     strtitle = sprintf('Position errors of 3D points for autocalibration method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
    
     figure;
     plot3Derrors( list_trans(:,list_ind_no_rot), trans_err{6}(ii,list_ind_no_rot))
    strtitle = sprintf('translation error for autocalibration method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);
    
     figure;
      plot3Derrors( list_trans(:,list_ind_no_rot), err_rep{6}(ii,list_ind_no_rot))
        strtitle = sprintf('image reprojection error for autocalibration method for %.2f pix. noise', noise_ampl(ii)); 
     title(strtitle);

    
     
   
 end
 
 autocalib_duration

 kk_ref = list_ind_no_rot(25);
K
fprintf('vrai rapport f/tz = %.2f\n', K(1,1) / list_cam_t_obj{kk_ref}(3));
for ii = 1:nb_noise
    fprintf('bruit %.1f px : f = %.0f, tz = %.1f, f/tz = %.2f\n', noise_ampl(ii), ...
            list_K_est{ii}(1,1), t_z_est(ii), list_K_est{ii}(1,1)/t_z_est(ii));
    disp(list_K_est{ii});
end
 
% Statistics (median and mex errors) per noise level

for jj=1:nb_noise	
	trans_err_cons = trans_err{6}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{6}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{6}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{6}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{6}(jj) = median(trans_err_filt,2); 
	%rot_err_med{6}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{6}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{6}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{6}(jj) = median(points_err_filt,2);
    points_err_max{6}(jj) = max(points_err_filt,[],2);
        
    im_err_med{6}(jj) = median(im_err_filt,2);
    im_err_max{6}(jj) = max(im_err_filt,[],2);
        
end

fprintf('mediane translation (mm) : analytique %s | autocalib %s\n', ...
        mat2str(trans_err_med{1},3), mat2str(trans_err_med{6},3));
fprintf('max translation (mm)     : analytique %s | autocalib %s\n', ...
        mat2str(trans_err_max{1},3), mat2str(trans_err_max{6},3));
     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%% Stereoscopy

% cameras configuration  
cam2_t_cam1 = [10;0;0];
cam2_R_cam1 = eye(3);

cam2_T_cam1 = [cam2_R_cam1 cam2_t_cam1; 0 0 0 1];
cam1_T_cam2 = inv(cam2_T_cam1);                               


% Generation of second camera data
% Compute perfect image points (without noise) for each configurationto test
for kk=1:nb_config_test
    
    % To complete
    P2 = cam2_R_cam1 * list_pts_in_cam{kk} + cam2_t_cam1;  
    list_pts_im_perf_cam2{kk} = K(1:2,:) * (P2 ./ P2(3,:));
end

% To complete
for ii=1:nb_noise 
	
	for kk=1:nb_config_test
		cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};
	
	
		ang_rand = 2*pi*rand(1, nb_pts);
		pts_im_noisy_cam2 = list_pts_im_perf_cam2{kk} + noise_ampl(ii) * [cos(ang_rand); sin(ang_rand)];
        
        mat_pts_im_cam2{ii}{kk} = pts_im_noisy_cam2;                
        
		
    end	
end




for ii=1:nb_noise 
   tic
   for kk=1:nb_config_test
        cam_t_obj_ref = list_cam_t_obj{kk}; 
		cam_R_obj_ref = list_cam_R_obj{kk};
	
        pts_im = mat_pts_im{ii}{kk};
       
        % Solving with basic triangulation
        for jj=1:nb_pts
            % Estimate 3D position of the point in camera frame
            m1 = K \ [pts_im(:,jj); 1];                         
            m2 = K \ [mat_pts_im_cam2{ii}{kk}(:,jj); 1];         
            Z = [cam2_R_cam1 * m1, -m2] \ (-cam2_t_cam1);       
            pos3D_est(1:3,jj) = Z(1) * m1;
        end
        % Estimate pose from 3D positions in camera frame
        [cam_R_obj_est, cam_t_obj_est] = pose_3D3D(object, pos3D_est);
        trans_err_triang(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
    
        % Question 7.3
        % Numerical optimization
		
        % compute image points in the normalized image plane for both
        % cameras
        pts_normed_im_plane = K \ [pts_im; ones(1, nb_pts)];
        pts_normed_im_plane_cam2 = K \ [mat_pts_im_cam2{ii}{kk}; ones(1, nb_pts)];
       
         % Initialization of objective function to minimize
		compute_err_stereo_persp_proj(0, object, pts_normed_im_plane, pts_normed_im_plane_cam2, cam1_T_cam2);
	
        % Initialization of pose
        % to complete
        [theta, u] = r2thetau(cam_R_obj_est);      
        x_ini = [theta * u; cam_t_obj_est];

        
        % Optimization using fminunc or lsqnonlin function
        % to complete
        % x = fminunc(@compute_err_stereo_persp_proj, x_ini);
	   
        % or
        options = optimoptions('lsqnonlin', 'Algorithm','levenberg-marquardt', 'Display','off'); 
        x = lsqnonlin(@compute_err_stereo_persp_proj, x_ini, [], [], options);
		
    
        % Extraction of cam_R_obj_est and cam_t_obj
        cam_R_obj_est = expm(skew(x(1:3)));
		cam_t_obj_est = x(4:6);
    
       

        trans_err{7}(ii,kk) = norm(cam_t_obj_ref - cam_t_obj_est);
    	
        [theta,u] = r2thetau(cam_R_obj_ref*cam_R_obj_est');
    	rot_err{7}(ii,kk) = abs(theta);
       
        err_rep{7}(ii,kk) = compute_err_rep_in_pix(cam_R_obj_est, cam_t_obj_est, object, pts_im, K);
        
        err_pts_3D{7}(ii,kk) = compute_err_pts3D(cam_R_obj_est, cam_t_obj_est, object, list_pts_in_cam{kk});
        
       
    end
  	stereo_duration(ii) = toc;
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), err_pts_3D{7}(ii,list_ind_no_rot))
    strtitle = sprintf('Position errors of 3D points for stereo method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), trans_err{7}(ii,list_ind_no_rot))
    strtitle = sprintf('translation error for stereo method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);
    
    figure;
    plot3Derrors( list_trans(:,list_ind_no_rot), err_rep{7}(ii,list_ind_no_rot))
    strtitle = sprintf('image reprojection error for stereo method for %.2f pix. noise', noise_ampl(ii)); 
    title(strtitle);

    
     
   
 end
 
 stereo_duration
 
% Statistics (median and mex errors) per noise level

for jj=1:nb_noise	
	trans_err_cons = trans_err{7}(jj,:);
	trans_err_filt = trans_err_cons(find(trans_err_cons ~= inf));	
				
	%rot_err_cons = rot_err{7}(jj,:);
	%rot_err_filt = rot_err_cons(find(rot_err_cons ~= inf));	
		
    points_err_cons = err_pts_3D{7}(jj,:);
    points_err_filt = points_err_cons(find(points_err_cons ~= inf));
        
    im_err_cons = err_rep{7}(jj,:);
    im_err_filt = im_err_cons(find(im_err_cons ~= inf));
        
        
    trans_err_med{7}(jj) = median(trans_err_filt,2); 
	%rot_err_med{7}(jj) = median(rot_err_filt,2); 
	
	trans_err_max{7}(jj) = max(trans_err_filt,[],2); 
	%rot_err_max{7}(jj) = max(rot_err_filt,[],2); 
	
    points_err_med{7}(jj) = median(points_err_filt,2);
    points_err_max{7}(jj) = max(points_err_filt,[],2);
        
    im_err_med{7}(jj) = median(im_err_filt,2);
    im_err_max{7}(jj) = max(im_err_filt,[],2);
        
end
     
 

col = ['r', 'g', 'b', 'k', 'm', 'c', 'y'];
figure;
hold on;
P = [];
for kk=list_tested_methods
	color = sprintf('%s*-', col(kk));
	P = [P plot(noise_ampl, trans_err_med{kk}, color)];
end
legend('Analytic', 'linear', 'linear wo norm', 'optimisation', 'OI', 'autocalib', 'stereo', 'Location', 'NorthWest');
xlabel('amplitude of image noise');
ylabel('translation error (mm)')
title('Median translation error');
set(P, 'linewidth', 2);




figure;
hold on;
P = [];
for kk=list_tested_methods
	color = sprintf('%s*-', col(kk));
	P = [P plot(noise_ampl, trans_err_max{kk}, color)];
end
legend('Analytic', 'linear', 'linear wo norm', 'optimisation', 'OI', 'autocalib', 'stereo', 'Location', 'NorthWest');
xlabel('amplitude of image noise');
ylabel('translation error (mm)')
title('Max translation error');
set(P, 'linewidth', 2);

 

figure;
hold on;
P = [];
for kk=list_tested_methods
	color = sprintf('%s*-', col(kk));
	P = [P plot(noise_ampl, points_err_med{kk}, color)];
end
legend('Analytic', 'linear', 'linear wo norm', 'optimisation', 'OI', 'autocalib', 'stereo', 'Location', 'NorthWest');
xlabel('amplitude of image noise');
ylabel('3D points error (mm)')
title('Median error on 3D points');
set(P, 'linewidth', 2);

figure;
hold on;
P = [];
for kk=list_tested_methods
	color = sprintf('%s*-', col(kk));
	P = [P plot(noise_ampl, points_err_max{kk}, color)];
end
legend('Analytic', 'linear', 'linear wo norm', 'optimisation', 'OI', 'autocalib', 'stereo', 'Location', 'NorthWest');
xlabel('amplitude of image noise');
ylabel('3D points error (mm)')
title('Max error on 3D points');
set(P, 'linewidth', 2);



figure;
hold on;
P = [];
for kk=list_tested_methods
	color = sprintf('%s*-', col(kk));
	P = [P plot(noise_ampl, im_err_med{kk}, color)];
end
legend('Analytic', 'linear', 'linear wo norm', 'optimisation', 'OI', 'autocalib', 'stereo', 'Location', 'NorthWest');
xlabel('amplitude of image noise');
ylabel('image points error (pix)')
title('Median error in image');
set(P, 'linewidth', 2);
figure;
hold on;
P = [];
for kk=list_tested_methods
	color = sprintf('%s*-', col(kk));
	P = [P plot(noise_ampl, im_err_max{kk}, color)];
end
legend('Analytic', 'linear', 'linear wo norm', 'optimisation', 'OI', 'autocalib', 'stereo', 'Location', 'NorthWest');
xlabel('amplitude of image noise');
ylabel('image points error (pix)')
title('Max error in image');
set(P, 'linewidth', 2);
