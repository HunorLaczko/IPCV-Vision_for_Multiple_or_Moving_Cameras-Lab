function [F,P,X3d,xc] = MatFunProjectiveCalib(x)
%  MatFunProjectiveCalib Computes a projective calibration from points in both images
%
% Input:
%  - x(3,npoints,2): hom. coords of the points in the images
%
% Output:
%  - F(3,3): fundamental matrix 
%  - P(3,4,2): projection matrices
%  - X3d(4,npoints): 3D points in hom. coordinates
%  - xc(3,npoints,ncam): hom. coords of the reprojected points

[dim,npoints,ncam] = size(x);
P = zeros(3,4,ncam);

% ------------------------------------------------------------------------
% Fundamental matrix using FDLT_Norm
% ------------------------------------------------------------------------
% ...
[F,cost]=FDLT_Norm(x(:,:,1),x(:,:,2));
% [F,inliersIndex] = estimateFundamentalMatrix(x(1:2,:,1)',x(1:2,:,2)', 'Method','RANSAC',...
%     'NumTrials',2000,'DistanceThreshold',2);
% ------------------------------------------------------------------------
% Compute the projection matrices from F using the definition in Hartley.
% You can use the function NumKernel to compute the right-side kernel of 
% matrix F (remember that the left-side kernel of F is equal the right-side
% kernel of transposed F)
% Warning! Be careful with the order of the cameras to be compatible with
% the computed F. See Hartley p. 256 
% ------------------------------------------------------------------------
P = zeros(3,4,2); % P(:,:,1) must contain proj. matrix of camera 1. P(:,:,2) must contain proj. matrix of camera 2 
% ...
P(:,1:3,1) = eye(3);

[Y,~] = NumKernel(F');
M = Cross2Matrix(Y);
v = [0;0;0]; %some vector V
lambda = 1; %some vector lambda

P(:,1:3,2) = (M*F + Y*v');
P(:,4,2) = lambda * Y;
% ------------------------------------------------------------------------
% Compute the 3D points with the function linear_triangulation
% ------------------------------------------------------------------------
% ...
X3d = linear_triangulation(x,P);
% ------------------------------------------------------------------------
% Compute the reprojected points 
% ------------------------------------------------------------------------
% ...
xc(:,:,1) = P(:,:,1)*X3d;
xc(:,:,2) = P(:,:,2)*X3d;

end
