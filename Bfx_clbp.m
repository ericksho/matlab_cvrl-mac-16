% [X,Xn,options] = Bfx_clbp(I,R,options)
% [X,Xn,options] = Bfx_clbp(I,options)
% [X,Xn] = Bfx_clbp(I,R,options)
% [X,Xn] = Bfx_clbp(I,options)
%
% Toolbox: Balu
%    Completed Local Binary Patterns features
%
%    X is the features vector, Xn is the list of feature names (see Example
%    to see how it works).
%
%    It calculates the CLBP over the a regular grid of patches. The function
%    uses Heikkila & Ahonen (see http://www.cse.oulu.fi/MVG/Research/LBP).
%
%    It returns a matrix of uniform lbp82 descriptors for I, made by
%    concatenating histograms of each grid cell in the image.
%    Grid size is options.hdiv * options.vdiv
%
%    R is a binary image or empty. If R is given the lbp will be computed
%    the corresponding pixles R==0 in image I will be set to 0.
%
%     Output:
%     X is a matrix of size ((hdiv*vdiv) x 59), each row has a
%         histogram corresponding to a grid cell. We use 59 bins.
%     options.x of size hdiv*vdiv is the x coordinates of center of ith grid cell
%     options.y of size hdiv*vdiv is the y coordinates of center of ith grid cell
%     Both coordinates are calculated as if image was a square of side length 1.
%
%     References:
%
%     CLBP Reference mising
%
%     Ojala, T.; Pietikainen, M. & Maenpaa, T. Multiresolution gray-scale
%     and rotation invariant texture classification with local binary
%     patterns. IEEE Transactions on Pattern Analysis and Machine
%     Intelligence, 2002, 24, 971-987.
%
%     Mu, Y. et al (2008): Discriminative Local Binary Patterns for Human
%     Detection in Personal Album. CVPR-2008.
%
%     Example 1:
%      options.vdiv = 1;                  % one vertical divition
%      options.hdiv = 1;                  % one horizontal divition
%      options.samples  = 8;              % number of neighbor samples
%      options.mappingtype = 'u2';        % uniform LBP
%      I = imread('testimg1.jpg');        % input image
%      figure(1);imshow(I,[])             % image to be analyzed
%      [X,Xn] = Bfx_lbp(I,[],options);    % LBP features
%      figure(2);bar(X)                   % histogram
%
%     Example 2:
%      options.vdiv = 2;                  % one vertical divition
%      options.hdiv = 2;                  % one horizontal divition
%      options.samples  = 8;              % number of neighbor samples
%      options.mappingtype = 'ri';        % rotation-invariant LBP
%      I = imread('testimg1.jpg');        % input image
%      J = I(120:219,120:239,2);          % region of interest (green)
%      figure(1);imshow(J,[])             % image to be analyzed
%      [X,Xn] = Bfx_lbp(J,[],options);    % LBP features
%      bar(X)                             % histogram
%
%   See also Bfx_lbp, Bfx_bsif, Bfx_gabor, Bfx_clp, Bfx_fourier, Bfx_dct.
%
% (c) GRIMA-DCCUC, 2011
% http://grima.ing.puc.cl
%
function [X,Xn,options] = Bfx_clbp(I,R,options)

if nargin==2;
    options = R;
    R = ones(size(I));
end

vdiv = options.vdiv;
hdiv = options.hdiv;

if ~isfield(options,'show')
    options.show = 0;
end

if options.show == 1
    disp('--- extracting completed local binary patern...');
end

if ~isempty(R);
    I(R==0) = 0;
end

if ~isfield(options, 'radius')
    options.radius = 1;
end

mapping=getmapping(options.samples,options.mappingtype); 

[N,M] = size(I);
w = N/vdiv;
h = M/hdiv;

X = [];
Xn = [];

for i = 1:vdiv
    for j = 1:hdiv
        part = I(i*w-w+1:w*i,j*h-h+1:h*j);
        [CLBP_SH,CLBP_MH] = clbp(I,option.radius,options.samples,mapping,'h');
        X = [X CLBP_SH];
        Xn = [Xn CLBP_MH];
    end
end

end



