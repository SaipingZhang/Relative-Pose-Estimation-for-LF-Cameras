function [R, T]=estimation_RT_LinearMethod(FeaturePoints1,FeaturePoints2,KMatrix)

N = size(FeaturePoints1,2);

Ar = zeros(3*N,9);
At = zeros(3*N,4);

fx = KMatrix(1,1);
fy = KMatrix(2,2);
cx = KMatrix(1,3);
cy = KMatrix(2,3);
k1 = -KMatrix(3,3);
k2 = -KMatrix(3,4);

n=1;
for i=1:N
    lfp1=FeaturePoints1(:,i)';
    lfp2=FeaturePoints2(:,i)';
    
    x1=lfp1(1,1);
    y1=lfp1(1,2);
    lambda1=lfp1(1,3);
    x2=lfp2(1,1); 
    y2=lfp2(1,2);
    lambda2=lfp2(1,3);
    
    Ar(n,1) = x1-cx;
    Ar(n,2) = (y1*fx)/fy-(cy*fx)/fy;
    Ar(n,3) = fx;
    Ar(n,7) = (x1*cx)/fx-cx^2/fx-x2*(x1/fx-cx/fx);
    Ar(n,8) = (y1*cx)/fy-x2*(y1/fy-cy/fy)-(cx*cy)/fy;
    Ar(n,9) = cx-x2;
    At(n,1)=- (fx*k1)/k2 - (fx*lambda1)/k2;
    At(n,3)=x2*(k1/k2 + lambda1/k2) - (cx*k1)/k2 - (cx*lambda1)/k2;
    
    Ar(n+1,4)=(x1*fy)/fx - (cx*fy)/fx;
    Ar(n+1,5)=y1 - cy;
    Ar(n+1,6)=fy;
    Ar(n+1,7) = (x1*cy)/fx-y2*(x1/fx-cx/fx)-(cx*cy)/fx;
    Ar(n+1,8) = (y1*cy)/fy - cy^2/fy - y2*(y1/fy - cy/fy);
    Ar(n+1,9) = cy - y2;
    At(n+1,2)= - (fy*k1)/k2 - (fy*lambda1)/k2;
    At(n+1,3)=y2*(k1/k2 + lambda1/k2) - (cy*k1)/k2 - (cy*lambda1)/k2;
    
    Ar(n+2,7)=(cx*k1)/fx - (x1*k1)/fx - lambda2*(x1/fx - cx/fx);
    Ar(n+2,8)=(cy*k1)/fy - (y1*k1)/fy - lambda2*(y1/fy - cy/fy);
    Ar(n+2,9)=- k1 - lambda2;
    At(n+2,3)=lambda2*(k1/k2 + lambda1/k2) + k1^2/k2 + (k1*lambda1)/k2;
    At(n+2,4)= k1 + lambda1; 
    n=n+3;
    
end

Matrix = (At * pinv(At) - eye(3*N,3*N)) * Ar;

[~,~,V] = svd( Matrix,'econ');
Rsolve = V(:,9);

Rest = reshape( Rsolve, 3,3 )';
R = project_to_rotation( Rest );
R1=[R(1,1);R(1,2);R(1,3);R(2,1);R(2,2);R(2,3);R(3,1);R(3,2);R(3,3)];

T = -pinv(At) * Ar * R1;
T1= T./T(end);
T = T1(1:3);

