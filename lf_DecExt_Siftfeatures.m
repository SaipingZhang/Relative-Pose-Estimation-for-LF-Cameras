function [ featuresC,DescriptorsC] = lf_DecExt_Siftfeatures( LF )
T = size( LF,1);
S = size( LF,2);
%center view index
Sc = (S+1)/2;
Tc = (T+1)/2;

% detect and extract Sift features for center view sub-aperture image
for t=1:T
    for s =1:S
        if s == Sc && t == Tc
            I = im2single(rgb2gray(squeeze(LF(t,s,:,:,:))));
            [featuresC,DescriptorsC] = vl_sift(I);
        end
    end
end