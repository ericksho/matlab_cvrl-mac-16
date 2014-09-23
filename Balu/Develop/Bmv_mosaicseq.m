function mosaic = Bmv_mosaicseq(f, seq, show)
% SIFT_MOSAIC Demonstrates matching two images using SIFT and RANSAC
%
%   SIFT_MOSAIC demonstrates matching two images based on SIFT
%   features and RANSAC and computing their mosaic.
%
%   SIFT_MOSAIC by itself runs the algorithm on two standard test
%   images. Use SIFT_MOSAIC(IM1,IM2) to compute the mosaic of two
%   custom images IM1 and IM2.
% AUTORIGHTS
%
% Adapted from function sift_mosaic written by A. Vedaldi
% http://www.vlfeat.org/applications/sift-mosaic-code.html



a = seq(1);
Ia = Bio_loadimg(f,a);

for b=seq(2:end)
    Ib = Bio_loadimg(f,b);
    mosaic = Bmv_mosaic2(Ia,Ib,0);
    Ia = double(mosaic);
    if show
        imshow(Ia,[])
        enterpause(0.1)
    end
end