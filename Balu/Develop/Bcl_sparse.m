% ds      = Bcl_sparse(X,d,Xt,options)  Training & Testing together
% options = Bcl_sparse(X,d,options)     Training only
% ds      = Bcl_sparse(Xt,options)      Testing only
%
% Toolbox: Balu
%    Classifications using sparse representations.
%
%    Design data:
%       X is a matrix with features (columns)
%       d is the ideal classification for X
%       options.type : '1' for one dictionary per class.
%       options.T    : sparsity (number of atoms used for the representation)
%       options.K    : number of atoms of the dictionary
%       options.iternum : number of iterations of KSVD algorithm
%
%    Test data:
%       Xt is a matrix with features (columns)
%
%    Output:
%       ds is the classification on test data
%       options.D dictionaries.
%
%    Example: Training & Test together:
%       load datareal                 % real data (2 classe, 279 features) 
%       [X,d,Xt,dt] = Bds_stratify(f,d,0.50);
%       op.type = 1;                  % one dictionary per class
%       op.T = 2;                     % sparsity (number of atoms used for the representation)
%       op.K = 40;                    % number of atoms of the dictionary
%       op.iternum = 30;              % number of iterations of KSV algorithm
%       ds = Bcl_sparse(X,d,Xt,op);   % sparse classifier
%       p = Bev_performance(ds,dt)    % performance on test data
%
%    Example: Training only
%       load datareal                 % real data (2 classe, 279 features) 
%       [X,d,Xt,dt] = Bds_stratify(f,d,0.75);
%       op.type = 1;                  % one dictionary per class
%       op.T = 2;                     % sparsity (number of atoms used for the representation)
%       op.K = 40;                    % number of atoms of the dictionary
%       op.iternum = 30;              % number of iterations of KSV algorithm
%       op = Bcl_sparse(X,d,op);      % sparse - training
%
%    Example: Testing only (after training only example):
%       ds = Bcl_sparse(Xt,op);       % sparse - testing
%       p = Bev_performance(ds,dt)    % performance on test data
%
%
% (c) GRIMA-DCCUC, 2012
% http://grima.ing.puc.cl


function [ds,options] = Bcl_sparse(varargin)
[train,test,X,d,Xt,options] = Bcl_construct(varargin{:});


options.string = 'sparse  ';
K    = options.K;
if train
    dmin = min(d);
    dmax = max(d);
    d    = d-dmin+1;
    C = dmax-dmin+1; % number of classes
    if C>2
        error('Bcl_sparse works only for two classes');
    end
    T    = options.T;
    i1   = find(d==1);
    i2   = find(d==2);
    Y1   = X(i1,:);
    Y2   = X(i2,:);

    params.Tdata = T;
    params.dictsize = K;
    params.iternum = options.iternum;
    params.memusage = 'high';

    params.data = Y1';
    [D1,g1] = ksvd(params,'');

    params.data = Y2';
    [D2,g2] = ksvd(params,'');


    options.D1    = D1;
    options.D2    = D2;
    options.g1    = g1;
    options.g2    = g2;

    
    ds = options;




end

if test
    D1 = options.D1;
    D2 = options.D2;
    X1 = full(options.g1)';
    X2 = full(options.g2)';
    Yt = Xt;
    Xt1 = omp(D1'*Yt',D1'*D1,K)'; % Sparsity-constrained Orthogonal Matching Pursuit
    Xt1 = full(Xt1); % pasa de sparse representation a normal

    Xt2 = omp(D2'*Yt',D2'*D2,K)'; % Sparsity-constrained Orthogonal Matching Pursuit
    Xt2 = full(Xt2); % pasa de sparse representation a normal
    switch options.type
        case 1 % minimal reconstruction error
            R1 = (Yt'-D1*Xt1')';
            n1 = sqrt(sum(R1.*R1,2));

            R2 = (Yt'-D2*Xt2')';
            n2 = sqrt(sum(R2.*R2,2));

            ds = (n1>n2)+1;
        case 2 % nearest neigbor to the sparse representation

            k = 1;

            kdtree1    = vl_kdtreebuild(X1');
            [i1,dist1] = vl_kdtreequery(kdtree1,X1',Xt1','NumNeighbors',k);

            kdtree2    = vl_kdtreebuild(X2');
            [i2,dist2] = vl_kdtreequery(kdtree2,X2',Xt2','NumNeighbors',k);

            ds = (dist1'>dist2')+1;

    end
end

