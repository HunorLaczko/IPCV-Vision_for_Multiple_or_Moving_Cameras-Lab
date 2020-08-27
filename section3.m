% include ACT_lite path
ACT_path = 'ACT_lite';
addpath(genpath(ACT_path));
% include extra funs
extra_funs_path = 'extra_funs';
addpath(genpath(extra_funs_path));

clear all

load('section1.mat');
load('section2.mat');

warning off
disp('************************************* START')

q_data = n_view_matching(points, features, ima, 0.45, params.Metric, 30);
q_data = homogenize_coords(q_data);

ncam = size(q_data,3);
npoints = size(q_data,2);

% ------------------------------------------------------------------------
% 2. Compute the fundamental matrix using the first and last cameras
% of the camera set (N cameras)
% ------------------------------------------------------------------------

q2_cams = zeros (3, npoints);
q2_cams(:,:,1) = q_data(:,:,1);
q2_cams(:,:,2) = q_data(:,:,ncam);

[F, P_2cam,Q_2cam,q_2cam_est] = MatFunProjectiveCalib(q2_cams);

disp(['Residual reprojection error. 8 point algorithm. First and alst camera   = ' num2str( ErrorRetroproy(q2_cams,P_2cam,Q_2cam)/2 )]);
draw_reproj_error(q2_cams,P_2cam,Q_2cam);

% ------------------------------------------------------------------------
% 3. Resection. Obtain the projection matrices of the rest of cameras using the PDLT_NA function 
% ------------------------------------------------------------------------

P_rep = zeros(3,4,ncam);
P_rep(:,:,[1 ncam]) = P_2cam;
for i = 2:ncam-1
    P_rep(:,:,i) =  PDLT_NA(q_data(:,:,i), Q_2cam);
end

disp(['Reprojection error for initial projective reconstruction using all cameras = ' num2str( ErrorRetroproy(q_data,P_rep,Q_2cam)/2 )]);
draw_reproj_error(q_data,P_rep,Q_2cam);

% % auxiliary matrix that indicates that all points are visible in all the cameras
vp = ones(npoints,ncam);

[P_ba,Q_ba] = BAProjectiveCalib(q_data,P_rep,Q_2cam,vp);

disp(['Reprojection error after bundle adjustment = ' num2str( ErrorRetroproy(q_data,P_ba,Q_ba)/2 )]);
draw_reproj_error(q_data,P_ba,Q_ba);

disp(['Reprojection error of views used for euclidean reconstruction = ' num2str( ErrorRetroproy(q_data(:,:,1:2),P_ba(:,:,1:2),Q_ba)/2 )]);
draw_reproj_error(q_data(:,:,1:2),P_ba(:,:,1:2),Q_ba);

% 4. Re-compute the Fundamental matrix between two of the cameras, 
% using the projection matrices obtained after the Projective Bundle Adjustment step

F = vgg_F_from_P(P_ba(:,:,1), P_ba(:,:,2));
E = K' * F * K;
[R_est,T_est] = factorize_E(E);

Rcam = zeros(3,3,2,4);
Tcam = zeros(3,2,4);

Rcam(:,:,2,1) = R_est(:,:,1);
Rcam(:,:,2,2) = R_est(:,:,1);
Rcam(:,:,2,3) = R_est(:,:,2);
Rcam(:,:,2,4) = R_est(:,:,2);

Tcam(:,2,1) = T_est;
Tcam(:,2,2) = -T_est;
Tcam(:,2,3) = T_est;
Tcam(:,2,4) = -T_est;

Rcam(:,:,1,1) = eye(3,3);
Rcam(:,:,1,2) = eye(3,3);
Rcam(:,:,1,3) = eye(3,3);
Rcam(:,:,1,4) = eye(3,3);

Q_euc = zeros(4,npoints); % Variable for recontructed points
P_euc = zeros(3,4,2);       % Variable for projection matrices
figNo=figure;
K2 = zeros(3,3,2);
K2(:,:,1) = K;
K2(:,:,2) = K;
Q_rep = project_points(P_ba, Q_ba);

figure();
for sol=1:4
    % Euclidean triangulation to obtain the 3D points (use TriangEuc)
    Q_euc = TriangEuc(Rcam(:,:,2,sol),Tcam(:,2,sol),K2,Q_rep);
       
    % visualize 3D reconstruction
    subplot(2,2,sol);
    draw_scene(Q_euc, K2, Rcam(:,:,:,sol), Tcam(:,:,sol));
    title(sprintf('Solution %d', sol));
     
%     % Compute the projection matrices from K, Rcam, Tcam
%     for k=1:2
%        P_euc(:,:,k) = K * [Rcam(:,:,k,sol), -Rcam(:,:,k,sol) * Tcam(:,k,sol)];
%     end
%     
%     % Obtain the re-projected points q_rep
%     q_rep = zeros(3, npoints, 2);
%     q_rep(:,:,1) = P_euc(:,:,1) * Q_euc;
%     q_rep(:,:,2) = P_euc(:,:,2) * Q_euc;
%     
%     % Visualize reprojectd points to check that all solutions correspond to
%     % the projected images
%     q_rep = un_homogenize_coords(q_rep);
%     for k=1:2
%       figure(figNo); subplot(4,2,2*(sol-1)+k); scatter(q_rep(1,:,k),q_rep(2,:,k),30,[1,0,0]);
%       title(sprintf('Reprojection %d, image %d', sol, k));
%       daspect([1, 1, 1]);
%       pbaspect([1, 1, 1]);
%       axis([-1000, 1000, -1000, 1000]);
%     end
end

