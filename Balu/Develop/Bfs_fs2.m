% Xs = Bfs_fs2(X,d,options)
%
% Toolbox: Balu
%    Feature selection of selected features.
%
%    options(1).fs defines the first feature selection algorithm, and
%    options(2).fs the second. In Xs are the selected features that
%    algorithm 2 selects from the features selected by algorithm 1.
%
%
% Example: 
%    load datareal
%    op.show = 1;                 % display results
%    op.b.name = 'fisher';        % SFS with Fisher
%    op.m0 = 40;                  % 40 features will be selected with SFS
%    op.m = 10;                   % from them 10 will be selected with LSEF
%    s = Bfs_sfslsef(f,d,op);     % index of selected features
%    Xn = fn(s,:)                 % list of feature names
%    Xs = f(:,s);                 % selected features
%    op_lda.p = [];
%    ds = Bcl_lda(Xs,d,Xs,op_lda);% LDA classifier
%    p = Bev_performance(d,ds)    % performance with sfs + lsef 
%
% (c) D.Mery, PUC-DCC, Jul. 2011
% http://dmery.ing.puc.cl

function selec = Bfs_fs2(X,d,options)


fs

op    = options;
op.m  = min([size(X,2) op.m0+20]);
s0    = Bfs_balu(X,d,op);
X0    = X(:,s0);
op.m  = options.m0;
s1    = Bfs_lsef(X0,d,op);
Xs    = X(:,s1);
op.m  = options.m;
s2    = Bfs_sfs(Xs,d,op);
selec = s0(s1(s2));
