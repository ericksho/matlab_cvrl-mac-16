BaluSetup                  % simulated data (2 classes, 2 features)
op.p = [];
ds = Bcl_lda(X,d,Xt,op);   % LDA classifier
p = Bev_performance(ds,dt) % performance on test data