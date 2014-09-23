clt
I = imread('X1.png');
ix = 315:394; jx=60:139; 
m1 = [mean(ix) mean(jx) 1]';
Iq = I(ix,jx);
It = imread('X2.png');
figure(1);
imshow(I,[]);
title('Query image: blue box')
hold on;
plot(m1(2),m1(1),'rx')
plot([min(jx) min(jx) max(jx) max(jx) min(jx)],[max(ix) min(ix) min(ix) max(ix) max(ix)])
figure(2);
imshow(It,[]);
title('Target image')
enterpause
disp('Searching matching region without epiolar restriction...')
op.q=80; op.d=5; op.show=1;op.nkp=3;op.fast=0;
D1 = Bmv_tqsift(Iq,It,op);
figure(3)
Bio_edgeview(It,bwperim(D1>10))
title('without epiplar line')
enterpause
disp('Searching matching region with epiolar restriction...')
close all
F  = Bmv_fundamentalSIFT(I,It);        % F matrix estimation
ell = F*m1;                            % epipolar line
R = Bmv_line2img(ell,size(It));
op.roi = imdilate(R,ones(20,20));
D2 = Bmv_tqsift(Iq,It,op);
figure(3)
Bio_edgeview(It,bwperim(D2>10))
title('with epiplar line')
hold on
Bmv_epiplot(F,m1);
