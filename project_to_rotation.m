%% Find nearest rotation matrix to a given matrix
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [R] = project_to_rotation( M )
    [V D] = eig( M'*M );
    isqrtMtM = 1/sqrt(D(1,1)) * V(:,1) * V(:,1)' + 1/sqrt(D(2,2)) * V(:,2) * V(:,2)' + 1/sqrt(D(3,3)) * V(:,3) * V(:,3)';
    R = M * isqrtMtM;
    R = R / det(R);
end
