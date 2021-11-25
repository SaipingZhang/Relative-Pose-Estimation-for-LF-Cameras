clear all;
clc;

% Please change to your own paths
run('D:\vlfeat-0.9.21\toolbox\vl_setup.m');
folder = 'D:\Relative Pose Estimation for LF Cameras\data\monkeyking\';

imageFiles={'IMG_1__Decoded', 'IMG_2__Decoded'};

% The number of iterations in Levenberg-Marquardt-on-manifold optimization
iter=10;

LF = cell(length(imageFiles),1);
pointsC = cell(length(imageFiles),1);
featuresC = cell(length(imageFiles),1);
vptsC = cell(length(imageFiles),1);
DescriptorsC = cell(length(imageFiles),1);

% Index of matched pairs
indexPairs = cell(length(imageFiles)-1,1);
matchedPoints1 = cell(length(imageFiles)-1,1);
matchedPoints2 = cell(length(imageFiles)-1,1);

% Disparities of matched feature points
zi1 = cell(length(imageFiles)-1,1);
zi2 = cell(length(imageFiles)-1,1);

% LF-Points: [x;y;lambda;1]
FeaturePoints1 = cell(length(imageFiles)-1,1);
FeaturePoints2 = cell(length(imageFiles)-1,1);

% index of selected feature points after ransac
Index = cell(length(imageFiles)-1,1);

% Poses of LF cameras
R = cell(length(imageFiles),1);
T = cell(length(imageFiles),1);
% When the first LF camera coordinate system is regarded as the world coordinate system
% R{2} and T{2} are absolute rotation and translation of the second LF camera, but also
% the relative rotation and translation between the first LF camera and the second LF camera
R{1} = eye(3);
T{1} = [0;0;0];

% Refined estimated relative pose after Levenberg-Marquardt-on-manifold optimization 
RefinedR = cell(length(imageFiles),1);
RefinedT = cell(length(imageFiles),1);
RefinedR{1} = eye(3);
RefinedT{1} = [0;0;0];

P = cell(length(imageFiles),1);

% disparity map
D = cell(length(imageFiles),1);

% projection matrix acquired by our calibration method
load('HMatrix.mat');

% the diameter of each micro-lens image in pixel
pixelsize=14.1;

% intrinsic parameter
fx = 1/HMatrix(3,3)/pixelsize;
fy = 1/HMatrix(4,4)/pixelsize;
cx = -HMatrix(3,5)/HMatrix(3,3)/pixelsize;
cy = -HMatrix(4,5)/HMatrix(4,4)/pixelsize;
k1 = HMatrix(3,1)/HMatrix(3,3)/pixelsize;
k2 = HMatrix(1,1)/HMatrix(3,3)/pixelsize;
KMatrix = [fx, 0, cx, 0;
    0, fy, cy, 0;
    0,  0,-k1,-k2;
    0,  0,  1, 0];

for j=1:length(imageFiles) %for j=1:2
    
    LF{j} = load([folder imageFiles{j} '.mat']);
    LF{j} = uint8(LF{j}.LF(:,:,:,:,1:3));
    fprintf('%s\n', imageFiles{j})
    
    D{j} = load([folder 'depth/E4/' imageFiles{j} '.mat']);
    D{j} = 0.02*D{j}.E4;

    %SIFT
    [featuresC{j},DescriptorsC{j}] = lf_DecExt_Siftfeatures(LF{j});
    % depth of integer pixel only
    ix = 1:size(D{j},2);
    iy = 1:size(D{j},1);
    iz = D{j};
    % interpolation for the depth of feature points which is at sub-pixel
    xi = featuresC{j}(1,:);
    yi = featuresC{j}(2,:);
    zi = interp2(ix,iy,iz,xi,yi,'linear');
    %featuresC{j}: [x;y;lambda;scale;orientation];
    featuresC{j} = [featuresC{j};zi];
    featuresC{j} = featuresC{j}([1 2 5 3 4],:);
end

for i = 1:length(imageFiles)-1 %for i = 1:1

    % Matching Features
    indexPairs{i} = vl_ubcmatch(DescriptorsC{i},DescriptorsC{i+1},5)';
    matchedPoints1{i} = featuresC{i}(1:2,indexPairs{i}(:, 1));
    matchedPoints2{i} = featuresC{i+1}(1:2,indexPairs{i}(:, 2));
    
    % disparity interpolation
    ix1 = 1:size(D{i},2);
    iy1 = 1:size(D{i},1);
    iz1 = D{i};
    ix2 = 1:size(D{i+1},2);
    iy2 = 1:size(D{i+1},1);
    iz2 = D{i+1};
   
    xi1 = matchedPoints1{i}(1,:);
    yi1 = matchedPoints1{i}(2,:);
    zi1{i} = interp2(ix1,iy1,iz1,xi1,yi1,'linear');
    
    xi2 = matchedPoints2{i}(1,:);
    yi2 = matchedPoints2{i}(2,:);
    zi2{i} = interp2(ix2,iy2,iz2,xi2,yi2,'linear');
    
    xy1 = matchedPoints1{i}(1:2,:);
    lambda1 = zi1{i};
    FeaturePoints1{i} = [xy1;lambda1;ones(1,size(xy1,2))];
    
    xy2 = matchedPoints2{i}(1:2,:);
    lambda2 = zi2{i};
    FeaturePoints2{i} = [xy2;lambda2;ones(1,size(xy2,2))];
    
    % Relative Pose Estimation
    [R{i+1}, T{i+1}] = estimation_RT_LinearMethod(FeaturePoints1{i},FeaturePoints2{i},KMatrix);
   
    % Refinement: lsq-on-manifold¡ª¡ªmargin effect if initial solution if good enough
    [invR,invT] = invert_Rt(R{i+1},T{i+1});
    RotationVector = rodriguesMatrixToVector(invR);%invR means the rotation of 2 relative to 1
    phi0 = [invT;RotationVector];
    f = @(phi)OptimizationFunction(phi,KMatrix,FeaturePoints1,FeaturePoints2);
    opts = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','display','iter','MaxFunEvals',inf,'MaxIter',100,'TolFun',eps,'TolX',eps);
    for loop = 1:iter
        [phi,fval] = lsqnonlin(f,phi0,[], [], opts);
        phi0 = logse3(expse3(phi-phi0) * expse3(phi0));
    end
    [RefinedR{i+1}, RefinedT{i+1}]  = invert_Rt(rodriguesVectorToMatrix(phi(4:6)), phi(1:3)); 
end

% Calculate pointclouds
TM = [1 0 0 0;0 1 0 0;0 0 1 0;0,0,0,1];

for i=1:length(imageFiles)
    xs = 1:size(D{i},2);
    ys = 1:size(D{i},1);
    [xs,ys] = meshgrid(xs,ys);
    xs = reshape(xs,[1,size(D{i},1)*size(D{i},2)]);
    ys = reshape(ys,[1,size(D{i},1)*size(D{i},2)]);
    zs = reshape(D{i},[1,size(D{i},1)*size(D{i},2)]);
    CP = [xs;ys;zs;ones(1,size(xs,2))];
    C=squeeze(LF{i}((size(LF{i},1)+1)/2,(size(LF{i},2)+1)/2,:,:,:));
    C_R = C(:,:,1);C_G= C(:,:,2);C_B = C(:,:,3);
    WP = KMatrix \ CP;
    WP = WP./WP(end,:);
    [invR,invT] = invert_Rt(RefinedR{i},RefinedT{i});
    TM = TM * [invR,invT;0,0,0,1];
    P{i} = TM(1:3,:) * WP;
    ptCloud(i) = pointCloud(P{i}','Color',[C_R(:),C_G(:),C_B(:)]);  
end

% show pointcloulds
xyzPoints=[];
c=[];
for i=1:length(imageFiles)
    tmp_xyz = ptCloud(i).Location;
    tmp_c = ptCloud(i).Color;  
    xyzPoints = [xyzPoints;tmp_xyz];
    c = [c;tmp_c];
end
xyzPoints(:,1) = -xyzPoints(:,1);
ptCloud_merge = pointCloud(xyzPoints,'Color',c);
figure;
pcshow(ptCloud_merge, 'parent', gca);


% Functions of REFINEMENT PART lsq-on-manifold
function FinalMeanError = OptimizationFunction(phi,KMatrix,FeaturePoints1,FeaturePoints2)

ScenePoints2 = KMatrix \ FeaturePoints2{1};
ScenePoints2 = ScenePoints2./ScenePoints2(4,:);
RotationMatrix = rotationVectorToMatrix(phi(4:6));
RotationMatrix = RotationMatrix';
TMatrix = [RotationMatrix,phi(1:3);0,0,0,1];
lfPoints1 = KMatrix * TMatrix * ScenePoints2;
lfPoints1 = lfPoints1./lfPoints1(4,:);
diff1 = (lfPoints1 - FeaturePoints1{1})./lfPoints1;

FinalMeanError = double(mean(vecnorm(diff1)));
end

function SE3 = expse3(xi)
if size(xi, 1) == 1
    xi = xi';
end
rho = xi(1:3,1);
phi = xi(4:6,1);
phi_wedge = vec2skew(phi);
J = 0;
for i = 0:10
    J = J + 1/factorial(i+1)*phi_wedge^i;
end
C = rotationVectorToMatrix(phi)';
SE3 = [C, J*rho;
    zeros(1,3),1];
end

function xi = logse3(T)
C = T(1:3,1:3);
phi = rotationMatrixToVector(C');
phi_wedge = vec2skew(phi);
J = 0;
for i = 0:10
    J = J + 1/factorial(i+1)*phi_wedge^i;
end
rho = J\T(1:3,4);
xi = [rho',phi]';
end

function skew = vec2skew(vec)
skew = [0, -vec(3), vec(2);
    vec(3), 0, -vec(1);
    -vec(2), vec(1), 0];
end

