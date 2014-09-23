%Image descriptor based on Histogram of Orientated Gradients for gray-level images. This code 
%was developed for the work: O. Ludwig, D. Delgado, V. Goncalves, and U. Nunes, 'Trainable 
%Classifier-Fusion Schemes: An Application To Pedestrian Detection,' In: 12th International IEEE 
%Conference On Intelligent Transportation Systems, 2009, St. Louis, 2009. V. 1. P. 432-437. In 
%case of publication with this code, please cite the paper above.

function H = Bfx_hog(I,options)

nj    = options.nj;       % number of HOG windows per bound box
ni    = options.ni;       % in i (vertical) and j (horizaontal) direction
B     = options.B;        % number of histogram bins

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
G     = ((Gi.^2)+(Gj.^2)).^.5;
K     = 4*di*dj;
ang   = zeros(K,1);
mag   = zeros(K,1);
for n=0:ni-1
    ii = n*di+1:(n+2)*di;
    for m=0:nj-1
        jj     = m*dj+1:(m+2)*dj;
        t      = t+1;
        ang(:) = A(ii,jj); 
        mag(:) = G(ii,jj);
        %assembling the histogram with 9 bins (range of 20 degrees per bin)
        b      = 0;
        H2     = zeros(B,1);
        for ang_lim=-pi+2*pi/B:2*pi/B:pi
            b = b+1;
            for k=1:K
                if ang(k)<ang_lim
                    ang(k) = 100;
                    H2(b)  = H2(b)+mag(k);
                end
            end
        end
        H((t-1)*B+1:t*B,1)=H2/(norm(H2)+0.01);
    end
end
