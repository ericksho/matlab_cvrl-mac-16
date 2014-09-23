% D = Bmv_tqlark(Iq,It,options)
%
% Toolbox: Balu
%
%    Search of image query (Iq) in target image (It) using LARK.
%
%    options.q    : sliding windows's size in pixels
%    options.d    : sliding step in pixels
%    options.nkp  : minimal number of matching keypoints
%    options.fast : '1' computes all SIFT keypoints of It at once
%                   '0' computes the SIFT keypoints for each sliding window
%    options.show : display results
%    options.roi  : region of interest where the matching in target image
%                   will be searched. If roi is not given, it will be
%                   considered that the search are is the whole image It.
%
%    D is the detection map.
%
%
%    Example 1:
%       I = imread('X1.png');
%       Iq = I(165:194,80:109);
%       It = imread('X2.png');
%       op.q=30; op.d=5; op.show=1;op.nkp=2;op.fast=1;
%       D = Bmv_tqlark(Iq,It,op);
%       figure
%       Bio_edgeview(It,bwperim(D>18))
%
%    Example 2:
%       I = imread('X1.png');
%       ix = 315:394; jx=60:139; 
%       m1 = [mean(ix) mean(jx) 1]';
%       Iq = I(ix,jx);
%       It = imread('X2.png');
%       figure(1);
%       imshow(I,[]);
%       title('Query image: blue box')
%       hold on;
%       plot(m1(2),m1(1),'rx')
%       plot([min(jx) min(jx) max(jx) max(jx) min(jx)],[max(ix) min(ix) min(ix) max(ix) max(ix)])
%       figure(2);
%       imshow(It,[]);
%       title('Target image')
%       enterpause
%       disp('Searching matching region without epiolar restriction...')
%       op.q=80; op.d=5; op.show=1;op.nkp=3;op.fast=0;
%       D1 = Bmv_tqsift(Iq,It,op);
%       figure(3)
%       Bio_edgeview(It,bwperim(D1>10))
%       title('without epiplar line')
%       enterpause
%       disp('Searching matching region with epiolar restriction...')
%       close all
%       F  = Bmv_fundamentalSIFT(I,It);        % F matrix estimation
%       ell = F*m1;                            % epipolar line
%       R = Bmv_line2img(ell,size(It));
%       op.roi = imdilate(R,ones(20,20));
%       D2 = Bmv_tqsift(Iq,It,op);
%       figure(3)
%       Bio_edgeview(It,bwperim(D2>10))
%       title('with epiplar line')
%       hold on
%       Bmv_epiplot(F,m1);
%
% (c) D.Mery, PUC-DCC, 2011
% http://dmery.ing.puc.cl

function [RM,options] = Bmv_tqlark(Iq,It,options)

show   = options.show;
energy = options.energy;
q      = options.q;
t      = options.t;
di     = options.di;

twoD   = size(Iq,3)==1;

Iq = double(Iq);
It = double(It);

if max2(Iq)>10
    Iq = Iq/256;
end
if max2(It)>10
    It = It/256;
end

if show
    if twoD
        close all
        figure(1)
        subplot(1,2,1);imshow(Iq,[]);title('Query Image')
        subplot(1,2,2);imshow(It,[]);title('Target Image')
    end
end


if isfield(options,'Sq')&&isfield(options,'Wq')&&...
        isfield(options,'indq')
    Sq   = options.Sq;
    Wq   = options.Wq;
    indq = options.indq;
    if show
        disp('LARK descriptors for query image already computed...');
    end
else
    [Sq,Wq,indq] = Bfx_lark(Iq,options);
end

if isfield(options,'St')&&isfield(options,'Wt')&&...
        isfield(options,'indt')
    St   = options.St;
    Wt   = options.Wt;
    indt = options.indt;
    if show
        disp('LARK descriptors for target image already computed...');
    end
else
    [St,Wt,indt] = Bfx_lark(It,options);
end


[Fq,lambda,A,Ws,mx,B] = Bft_pca(Wq',energy);
Fq = Fq';
N = size(Wt,2);
mx = mean(Wt,2)';
MX = ones(N,1)*mx;
W0 = Wt' - MX;
Ft = W0*B;
Ft = Ft';

nFq = norm(Fq(:));
Fqn = Fq(:)'/nFq;
mi = indq(1);
mj = indq(2);
i2 = round(mi/2)-1;
j2 = round(mj/2)-1;

if twoD==1 % 2D
    
    RMo = zeros(indt(1),indt(2));
    
    for i=1:di:indt(1)-mi+1
        iss = i:i+mi-1;
        ion = ones(mj,1)*iss; ion = ion(:);
        for j=1:di:indt(2)-mj+1
            jss = j:j+mj-1;
            jon = jss'*ones(1,mi); jon = jon(:);
            Fti = Ft(:,(ion-1)*indt(2)+jon);
            rho = Fqn*Fti(:)/norm(Fti(:)); % = trace((Fq'*Fti)/(nFq*norm(Fti)));
            RMo(i+i2,j+j2) = rho^2/(1-rho^2);
        end
    end
    RM1 = imresize(RMo,indt*q);
    [N,M] = size(It);
    RM  = zeros(N,M);
    [N1,M1] = size(RM1);
    RM(1+t:N1+t,1+t:M1+t) = RM1;
    %RM = imresize(RMo,[N M]);
    
    if show
        if isfield(options,'thRM')
            figure
            imshow((RM/max2(RM)*3+It)*16,jet)
            %figure
            %imshow(It,[])
            hold on
            Do = Bim_nonmaxsup(RMo,options.thRM,max([mi mj])/2);
            D  = zeros(N,M);
            warning off
            D1 = imresize(Do,indt*q);
            warning on
            D(1+t:N1+t,1+t:M1+t) = D1;
            D = imdilate(D,ones(mi*q,mj*q));
            options.D = D;
%            [ii,jj] = find(D==1);
%            ii = ii-1;
%            jj = jj-1;
%            nn = length(ii);
%            for i=1:nn
%                x = [jj(i)-j2    jj(i)+j2 jj(i)+j2 jj(i)-j2 jj(i)-j2];
%                y = [ii(i)+i2 ii(i)+i2 ii(i)-i2    ii(i)-i2 ii(i)+i2];
%                plot((x+0.5)*q+1+t,(y+0.5)*q+1+t,'r')
%            end
       end
    end
    
    
else %3D
    mk = indq(3);
    RMo = zeros(indt(1),indt(2),indt(3));
    k2 = round(indq(3)/2)-1;
    
    for i=1:di:indt(1)-mi+1
        iss = i:i+mi-1;
        ion = ones(mj*mk,1)*iss; ion = ion(:);
        for j=1:di:indt(2)-mj+1
            jss = j:j+mj-1;
            jon = jss'*ones(1,mk);jon=jon';jon = jon(:);jon=jon*ones(1,mi);jon=jon(:);
            for k=1:di:indt(3)-indq(3)+1
                kss = k:k+indq(3)-1;
                kon = kss'*ones(1,mi*mj);kon = kon(:);
                Fti = Ft(:,(ion-1)*indt(2)*indt(3)+(jon-1)*indt(3)+kon);
                rho = Fqn*Fti(:)/norm(Fti(:)); % = trace((Fq'*Fti)/(nFq*norm(Fti)));
                RMo(i+i2,j+j2,k+k2) = rho^2/(1-rho^2);
            end
        end
    end
    [N,M,P] = size(It);
    sp = fix((P - indt(3))/2);
    RM = zeros(N,M,P);
    RM(:,:,sp:sp+indt(3)-1) = imresize(RMo,[N M]);
    
    
    if show
        close all
        if isfield(options,'thRM')
            %             for k=1:P
            %                 figure(k)
            %                 clf
            %                 imshow((RM(:,:,k)/max2(RM)*3+It(:,:,k))*16,jet)
            %                 hold on
            %             end
            
            D = Bim_nonmaxsup(RMo,options.thRM,max([mi mj mk])/2);
            
            DM = zeros(N,M,P);
            DM(:,:,sp:sp+indt(3)-1) = imdilate(imresize(D,[N M]),ones(11,11));
            for k=1:P
                figure(k)
                clf
%                subplot(1,2,1);imshow((DM(:,:,k)*3+It(:,:,k))*16,jet)
                subplot(1,2,1);imshow(It(:,:,k),[])
                subplot(1,2,2);imshow(DM(:,:,k),[])
                enterpause
            end
            
            
            
            options.D    = DM;
            %             [ii,jj,kk] = find(D==1);
            %             nn = length(ii);
            %             for i=1:nn
            %                 figure(kk(i))
            %                 x = [jj(i)-j2    jj(i)+j2 jj(i)+j2 jj(i)-j2 jj(i)-j2];
            %                 y = [ii(i)+i2 ii(i)+i2 ii(i)-i2    ii(i)-i2 ii(i)+i2];
            %                 plot((x+0.5)*q(2)+1+t,(y+0.5)*q(1)+1+t,'r')
            %                 enterpause
            %             end
        end
    end
end
options.RMo  = RMo;
options.Sq   = Sq;
options.Wq   = Wq;
options.Fq   = Fq;
options.St   = St;
options.Wt   = Wt;
options.Ft   = Ft;
options.indq = indq;
options.indt = indt;
