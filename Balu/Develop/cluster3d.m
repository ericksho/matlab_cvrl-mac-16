X = Ms(1:3,:)';


m = size(X,2);


N = size(X,1);
d = zeros(N,1);
c = 1;

d(1) = c;


while(sum(d>0)<N)
    i0 = find(d==0);
    ic = find(d==c);
    X0 = X(i0,:);
    n0 = length(i0);
    Xc = ones(n0,1)*mean(X(ic,:),1);
    Xd = Xc-X0;
    M2 = sqrt(sum(Xd.*Xd,2));
    i2 = find(M2<th);
    if ~isempty(i2)
        d(i0(i2)) = c;
    else
        c = c+1;
        d(i0(1))=c;
    end
end
Xc = zeros(c,m);
for i=1:c
    Xc(i,:) = mean(X(d==i,:),1);
end