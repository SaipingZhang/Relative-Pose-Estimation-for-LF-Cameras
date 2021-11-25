%% Compute the inverse of a rotation/translation combo [R t]
%% (first R, then t)
%
% (c) Bastian Goldluecke 11/2014, University of Konstanz
% License: Creative Commons BY-SA 4.0
%          Please cite our paper if you use our code in your work.
%
function [Ri ti] = invert_Rt( R,t )
    Ri = R';
    ti = -Ri * t;
end
