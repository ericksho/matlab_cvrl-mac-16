% J = Bim_rgb2hcm(RGB)
%
% Toolbox: Balu
%    Conversion RGB to high contrast image.
%
%    RGB: color image
%    J  : hcm image
%
%  See details in:
%  Mery, D.; Pedreschi, F. (2005): Segmentation of Colour Food Images using 
%  a Robust Algorithm. Journal of Food Engineering 66(3): 353-360.
%
%  Example:
%     I = imread('testimg2.jpg');
%     J = Bim_rgb2hcm(I);
%     figure(1)
%     imshow(I); title('control image')
%     figure(2)
%     imshow(J); title('high contrast image')
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl


function J = Bim_rgb2hci(RGB)

RGB = double(RGB);
if (size(RGB,3)==1)
    I = RGB;
else
    RGB64 = imresize(RGB,[64 64]);
    k = fminsearch(@Bstd2,[1 1],[],RGB64);
    I = k(1)*RGB(:,:,1) + k(2)*RGB(:,:,2) + RGB(:,:,3);
end
J = I - min(I(:));
J = J/max(J(:));
n = fix(size(J,1)/4);
if (mean2(J(1:n,1:n)) > 0.4)
    J = 1 - J;
end
end


% s = Bstd2(k,RGB)
%
% Toolbox: Balu
%    Standard deviation of normalized image I, where
%    I = k(1)*R+k(2)*G+B (R = RGB(:,:,1), G=RGB(:,:,2), B = RGB(:,:,3)
%    s: Standard deviation
%
%    This function is called by Brgb2hcm, Bsegbalu.
%
% See also Brgb2hcm, Bsegbalu.
%
% D.Mery, PUC-DCC, Apr. 2008
% http://dmery.ing.puc.cl
%

function s = Bstd2(k,RGB)
% I1 = Bim_maxmin(RGB(:,:,1));
% I2 = Bim_maxmin(RGB(:,:,2));
% I3 = Bim_maxmin(RGB(:,:,3));
% I1n = I1(:)-graythresh(I1);
% I2n = I2(:)-graythresh(I2);
% I3n = I3(:)-graythresh(I3);
% 
% max1 = 1-g1;
% min1 = -g1;
% max2 = 1-g2;
% min2 = -g2;
% max3 = 1-g3;
% min3 = -g3;
% 
% 
% I1s = zeros(size(I1));
% I2s = I1s;
% I3s = I1s;
% 
% I1s(I1n>=0) = I1n(I1n>=0)/max1;
% I1s(I1n<0)  = I1n(I1n<0)/min1;
% 
% I2s(I2n>=0) = I1n(I2n>=0)/max1;
% I2s(I2n<0)  = I1n(I2n<0)/min1;
% 
% I3s(I3n>=0) = I1n(I3n>=0)/max1;
% I3s(I3n<0)  = I1n(I3n<0)/min1;
% 


I = k(1)*RGB(:,:,1) + k(2)*RGB(:,:,2) + RGB(:,:,3);
I = Bim_maxmin(I);
g = graythresh(I);
i1 = find(I>g);
i2 = find(I<=g);
m1 = mean(I(i1));
m2 = mean(I(i2));
v1 = var(I(i1));
v2 = var(I(i2));
s = -(m1-m2)/(v1+v2);

end