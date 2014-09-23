% ds      = Bcl_knn_prox(X,d,Xt,options)  Training & Testing together
% options = Bcl_knn_prox(X,d,options)     Training only
% ds      = Bcl_knn_prox(Xt,options)      Testing only
%
% Toolbox: Balu
%    KNN (k-nearest neighbors) classifier using proximity search algorithm:
%
%    Chavez, Figueroa, Navarro (2008): Effective Proximity Retrieval by
%    Ordering Permutations. PAMI 30(9):1647-1658.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.k is the number of neighbors
%       options.k2 is the number of permutants (variable k in the paper)
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%
%    Example: Training & Test together:
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.k  = 15;
%       op.k2 = 15;
%       op.f  = 20;
%       ds = Bcl_knn_prox(X,d,Xt,op);   % knn with 15 neighbors
%       p = Bev_performance(ds,dt) % performance on test data
%
%    Example: Training only
%       load datagauss             % simulated data (2 classes, 2 features)
%       Bio_plotfeatures(X,d)      % plot feature space
%       op.k  = 15;
%       op.k2 = 15;
%       op.f  = 20;
%       op = Bcl_knn_prox(X,d,op);   % knn with 15 neighbors
%
%    Example: Testing only (after training only example):
%       ds = Bcl_knn_prox(Xt,op);       
%       p = Bev_performance(ds,dt) % performance on test data
%
%
% D.Mery, PUC-DCC, 2013
% http://dmery.ing.puc.cl


function ds = Bcl_knn_prox(varargin)

% U   : Training data (one row per instance)
% P   : Permutants
% Pix : Permutations

[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});
options.string = sprintf('knnp,%2d ',options.k);
K   = options.k;  % number of neighbors
k   = options.k2; % variable k of the paper
f   = options.f;
if train
    U         = X';
    options.U = U;
    options.d = d;
    n         = size(U,2);
    options.n = n;
    if f>1
        N = f;
    else
        N   = f*n;
    end
    t         = rand(n,1);[i,j]=sort(t);
    P         = U(:,j(1:k));
    options.P = P;
    options.n = n;
    options.N = N;

    % Indexing
    options.kd  = vl_kdtreebuild(P);
    options.Pix = vl_kdtreequery(options.kd,P,U,'NumNeighbors',k);
    ds = options;
end
if test
    U         = options.U;
    N         = options.N;
    d         = options.d;
    jj        = repmat((1:k)',1,options.n);
    Nt        = size(Xt,1);
    ds        = zeros(Nt,1);
    PIq       = vl_kdtreequery(options.kd,options.P,Xt','NumNeighbors',k);
    [PIs,PJs] = sort(single(PIq));
    for t = 1:Nt
        Piq1    = PJs(:,t);
        D1      = Piq1(options.Pix)-jj;
        [ai,aj] = sort(sum(D1.*D1));
        D2      = U(:,aj(1:N))-repmat(Xt(t,:)',1,N);
        [id,jd] = sort(sum(D2.*D2));
        ds(t)   = mode(d(aj(jd(1:K))));
    end

end
