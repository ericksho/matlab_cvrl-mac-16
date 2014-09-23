function D = Bim_nonmaxsup(S,th,r)

[N,M,P] = size(S);

So = S(:);
Sm = So>th;

r2 = r^2;
if P==1 % 2D
    
    NM = [N M];
    D = zeros(N,M);
    
    while sum(Sm) > 1
        
        [smx,g] = max(So);
        [i1,j1] = ind2sub(NM,g);
        Sm(g) = 0;
        So(g) = 0;
        D(i1,j1) = 1;
        ii = find(Sm==1);
        [i2,j2] = ind2sub(NM,ii);
        n2 = length(i2);
        d = ones(n2,1)*[i1 j1]-[i2 j2];
        dm = sum(d.*d,2);
        ii = find(dm<r2);
        if not(isempty(ii))>0
            g = sub2ind(NM,i2(ii),j2(ii));
            Sm(g) = 0;
            So(g) = 0;
        end
        
    end
    if sum(Sm)==1
        Do = zeros(N,M);
        Do(:) = Sm;
        D = or(D,Do);
    end
else % 3D
    
    NMP = [N M P];
    D = zeros(N,M,P);
    
    while sum(Sm) > 1
        
        [smx,g] = max(So);
        [i1,j1,k1] = ind2sub(NMP,g);
        Sm(g) = 0;
        So(g) = 0;
        D(i1,j1,k1) = 1;
        ii = find(Sm==1);
        [i2,j2,k2] = ind2sub(NMP,ii);
        n2 = length(i2);
        d = ones(n2,1)*[i1 j1 k1]-[i2 j2 k2];
        dm = sum(d.*d,2);
        ii = find(dm<r2);
        if not(isempty(ii))>0
            g = sub2ind(NMP,i2(ii),j2(ii),k2(ii));
            Sm(g) = 0;
            So(g) = 0;
        end
        
    end
    if sum(Sm)==1
        Do = zeros(N,M,P);
        Do(:) = Sm;
        D = or(D,Do);
    end
    
end