% Neighbors using proximity search algorithm:
%
% Chavez, Figueroa, Navarro (2011): Effective Proximity Retrieval by
% Ordering Permutations. PAMI 30(9):1647-1658.
%
% U   : Training data (one row per instance)
% P   : Permutants
% Pix : Permutations
%
% op.k=15; op.k2= 15; op.f=20;tic;ds = Bcl_knn_prox(X1,d1,X2(1,:),op);toc;Bev_performance(ds,d2(1))

function ds = Bcl_knn_prox(X,d,Xt,options)

K   = options.k;  % number of neighbors
k   = options.k2; % variable k of the paper
f   = options.f;
U   = X';
n   = size(U,2);
% m   = size(U,1);
if f>1
    N = f;
else
    N   = f*n;
end
t   = rand(n,1);[i,j]=sort(t);
P   = U(:,j(1:k));

% Indexing
kd  = vl_kdtreebuild(P);
Pix = vl_kdtreequery(kd,P,U,'NumNeighbors',k);
toc
% Searching
% j         = (1:k)';
jj        = repmat((1:k)',1,n);
Nt        = size(Xt,1);
ds        = zeros(Nt,1);
PIq       = vl_kdtreequery(kd,P,Xt','NumNeighbors',k);
[PIs,PJs] = sort(single(PIq));
toc
for t = 1:Nt
    Piq1    = PJs(:,t);
    D1      = Piq1(Pix)-jj;
    [ai,aj] = sort(sum(D1.*D1));
    D2      = U(:,aj(1:N))-repmat(Xt(t,:)',1,N);
    [id,jd] = sort(sum(D2.*D2));
    ds(t)   = mode(d(aj(jd(1:K))));
end