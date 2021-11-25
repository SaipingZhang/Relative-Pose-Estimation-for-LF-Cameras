function [rotationVector, dvdM] = rodriguesMatrixToVector(rotationMatrix)
% rodriguesMatrixToVector Convert a 3D rotation matrix into a rotation
% vector.
%
% rotationVector = rodriguesMatrixToVector(rotationMatrix) computes a 
% rotation vector (axis-angle representation) corresponding to a 3D 
% rotation matrix using the Rodrigues formula.
%
% rotationMatrix is a 3x3 3D rotation matrix.
%
% rotationVector is a 3-element rotation vector corresponding to the 
% rotationMatrix. The vector represents the axis of rotation in 3D, and its 
% magnitude is the rotation angle in radians.
%
% [..., dvdM] = rodriguesMatrixToVector(rotationMatrix) additionally
% returns the 3-by-9 derivatives of rotationVector w.r.t rotationMatrix.
%
% See also vision.internal.calibration.rodriguesVectorToMatrix

% References:
% [1] R. Hartley, A. Zisserman, "Multiple View Geometry in Computer
%     Vision," Cambridge University Press, 2003.
% 
% [2] E. Trucco, A. Verri. "Introductory Techniques for 3-D Computer
%     Vision," Prentice Hall, 1998.

% Copyright 2012 MathWorks, Inc.

%#codegen

% get the rotation matrix that is the closest approximation to the input
[U, ~, V] = svd(rotationMatrix);
rotationMatrix = U * V';

t = trace(rotationMatrix);
% t is the sum of the eigenvalues of the rotationMatrix.
% The eigenvalues are 1, cos(theta) + i sin(theta), cos(theta) - i sin(theta)
% t = 1 + 2 cos(theta), -1 <= t <= 3

tr = (t - 1) / 2;
theta = real(acos(tr));

r = [rotationMatrix(3,2) - rotationMatrix(2,3); ...
     rotationMatrix(1,3) - rotationMatrix(3,1); ...
     rotationMatrix(2,1) - rotationMatrix(1,2)];

needJacobian = nargout > 1;
outType = class(rotationMatrix);

if needJacobian
    dtrdM = [1 0 0 0 1 0 0 0 1] / 2;
    dtrdM = cast(dtrdM, outType);
    
    drdM = [0 0 0 0 0 1 0 -1 0;...
            0 0 -1 0 0 0 1 0 0;...
            0 1 0 -1 0 0 0 0 0];
    drdM = cast(drdM, outType);
end

threshold = cast(1e-4, outType); 

if sin(theta) >= threshold
    % theta is not close to 0 or pi
    vth = 1 / (2*sin(theta));
    v = r * vth;
    rotationVector = theta * v;

    if needJacobian        
        dthetadtr = -1/sqrt(1-tr^2);
        dthetadM = dthetadtr * dtrdM;
        
        % var1 = [vth; theta]
        dvthdtheta = -vth*cos(theta)/sin(theta);
        dvar1dtheta = [dvthdtheta; 1];
        dvar1dM =  dvar1dtheta * dthetadM;

        % var = [r;vth;theta]
        dvardM = [drdM;dvar1dM];

        % var2 = [r;theta]
        dvdvar = [vth*eye(3,outType) r zeros(3,1,outType)];
        dthetadvar = cast([0 0 0 0 1], outType);
        dvar2dvar = [dvdvar; dthetadvar];

        dVdvar2 = [theta*eye(3,outType) v];

        dvdM = dVdvar2 * dvar2dvar * dvardM;
    end
elseif t-1 > 0
    % theta is close to 0
    rotationVector = (.5 - (t - 3) / 12) * r;
    if needJacobian   
        dvdM = (-1 / 12) * r * dtrdM + (.5 - (t - 3) / 12) * drdM;
    end
else
    % theta is close to pi
    rotationVector = ...
        computeRotationVectorForAnglesCloseToPi(rotationMatrix, theta);
    if needJacobian
        % No definition for this case
        dvdM = zeros(3, 9, outType);
    end
end

function rotationVector = ...
    computeRotationVectorForAnglesCloseToPi(rotationMatrix, theta)
% r = theta * v / |v|, where (w, v) is a unit quaternion.
% This formulation is derived by going from rotation matrix to unit
% quaternion to axis-angle

% choose the largest diagonal element to avoid a square root of a negative
% number
[~, a] = max(diag(rotationMatrix));
a = a(1);
b = mod(a, 3) + 1;
c = mod(a+1, 3) + 1;

% compute the axis vector
s = sqrt(rotationMatrix(a, a) - rotationMatrix(b, b) - rotationMatrix(c, c) + 1);
v = zeros(3, 1, 'like', rotationMatrix);
v(a) = s / 2;
v(b) = (rotationMatrix(b, a) + rotationMatrix(a, b)) / (2 * s);
v(c) = (rotationMatrix(c, a) + rotationMatrix(a, c)) / (2 * s);

rotationVector = theta * v / norm(v);
