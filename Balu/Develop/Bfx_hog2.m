
function [H,w,J] = Bfx_hog2(I,options)

nj    = options.nj;       % number of HOG windows per bound box
ni    = options.ni;       % in i (vertical) and j (horizaontal) direction
B     = options.B;        % number of histogram bins
show  = options.show;     % show histograms

[N,M] = size(I);          % N num of lines ; M num of columns
H     = zeros(nj*ni*B,1); % column vector with zeros
I     = double(I);
dj    = floor(M/(nj+1));
di    = floor(N/(ni+1));
t     = 0;
hj    = [-1,0,1];
hi    = -hj';
Gj    = imfilter(I,hj);
Gi    = imfilter(I,hi);
A     = atan2(Gi,Gj);
A     = mod(A,pi);
G     = ((Gi.^2)+(Gj.^2)).^.5;
K     = 4*di*dj;
ang   = zeros(K,1);
mag   = zeros(K,1);
if show
    J    = zeros(N,M);
    ss   = round(min([N/ni M/nj])*0.45);
    bs   = 2*ss+1;
    bim1 = zeros(bs, bs);
    bim1(:,ss) = 1;
    bim = zeros([size(bim1) B]);
    bim(:,:,1) = bim1;
    for i = 2:B
        bim(:,:,i) = imrotate(bim1, -(i-1)*180/B, 'crop');
    end
end
w = zeros(ni,nj,B);
for i = 1:ni
    ii = (i-1)*di+1:(i+1)*di;
    i0 = round(mean(ii));
    for j = 1:nj
        jj     = (j-1)*dj+1:(j+1)*dj;
        j0     = round(mean(jj));
        t      = t+1;
        ang(:) = A(ii,jj);
        mag(:) = G(ii,jj);
        %assembling the histogram with B bins
        H2     = zeros(B,1);
        for b = 1:B
            q      = find(ang<=pi*b/B);
            H2(b)  = H2(b)+sum(mag(q));
            ang(q) = 5;
        end

        H2                = H2/(norm(H2)+0.01);
        H(indices(t,B),1) = H2;
        w(i,j,:)          = H2;
        if show
            for b=1:B
                %ang_lim = pi*b/B+pi/2;
                %for q=-ss:ss
                %    qi = round(i0+q*cos(ang_lim));
                %    qj = round(j0+q*sin(ang_lim));
                %    J(qi,qj) = J(qi,qj)+H2(b);
                %end
            J(i0-ss:i0+ss,j0-ss:j0+ss) = J(i0-ss:i0+ss,j0-ss:j0+ss) + bim(:,:,b)*H2(b);
                
            end
            % J(round(i0),round(j0)) = 1;
            %                imshow(J,[])
            %                enterpause(0)
        end
    end
end
if show
    imshow(J,[])
end