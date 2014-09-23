I = imread('testimg7.bmp');
J = imread('testimg8.bmp');

Itrain = I(:,1:250);
Jtrain = J(:,1:250);
Itest  = I(:,251:500);
Jtest  = J(:,251:500);

Itest = I;
Jtest = J;


[ftrain,dtrain] = vl_sift(single(Itrain));
i  = round(ftrain(2,:));
j  = round(ftrain(1,:));
k  = sub2ind(size(Itrain),i,j);
d  = Jtrain(k);
i1 = find(d==1);
i0 = find(d==0);

f0 = ftrain(:,i0);
f1 = ftrain(:,i1);
d0 = dtrain(:,i0);
d1 = dtrain(:,i1);

V = 50;
c0 = single(vl_kmeans(single(d0),V,'Algorithm','Elkan'));
c1 = single(vl_kmeans(single(d1),V,'Algorithm','Elkan'));


k0 = vl_kdtreebuild(c0);
k1 = vl_kdtreebuild(c1);

[ftest,dtest] = vl_sift(single(Itest));

[j0,D0]  = vl_kdtreequery(k0,c0,single(dtest),'NumNeighbors',1);
[j1,D1]  = vl_kdtreequery(k1,c1,single(dtest),'NumNeighbors',1);

ds = D0>D1;


i  = round(ftest(2,ds));
j  = round(ftest(1,ds));
k  = sub2ind(size(Itest),i,j);
Jtest = zeros(size(Itest));
Jtest(k) = 1;
figure(1)
imshow(Itest)
hold on
plot(j,i,'r.')
