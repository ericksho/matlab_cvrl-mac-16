function [S,W,ind,indr] = Bfx_lark(I,options)

q     = options.q; % mask window
t     = options.t; % mask for covariance matrix
alpha = options.alpha;
show  = options.show;
s     = options.s;
hd2   = (options.hd)^2;

[N,M,P] = size(I);
S       = zeros(size(I));



if P==1 %2D
    if show
        figure
        disp('Bfx_lark: computing 2D-LARK descriptors...');
    end
    
    % gradient
    Gx = conv2(I,s','same');
    Gy = conv2(I,s ,'same');
    
    % covariance matrices
    C = zeros(2,2,N,M);
    e     = options.e;
    tau   = options.tau;
    i_im  = sqrt(-1);
    for i=2+t:N-t-1
        for j=2+t:M-t-1
            gx      = Gx(i-t:i+t,j-t:j+t);
            gy      = Gy(i-t:i+t,j-t:j+t);
            J       = [gx(:) gy(:)];
            c       = J'*J;
            lambda1 = ((c(1,1)+c(2,2))+sqrt((c(1,1)-c(2,2))^2+4*c(1,2)*c(2,1)))/2;
            lambda2 = ((c(1,1)+c(2,2))-sqrt((c(1,1)-c(2,2))^2+4*c(1,2)*c(2,1)))/2;
            r       = (c(2,2)+c(1,2)-lambda1)-i_im*(c(1,1)+c(2,1)-lambda1);
            %theta   = atan(-(c(1,1)+c(2,1)-lambda1)/(c(2,2)+c(1,2)-lambda1));
            theta = angle(r);
            u1      = [cos(theta) sin(theta)];
            u2      = [-u1(2) u1(1)]; % = [-sin(theta) cos(theta)];
            s1      = sqrt(lambda1);
            s2      = sqrt(lambda2);
            creg    = ((s1*s2+e)^alpha)*((s1+tau)/(s2+tau)*(u1'*u1)+(s2+tau)/(s1+tau)*(u2'*u2));
            
            x = sum(isnan(creg(:))==1);
            if x>1
                C(:,:,i,j) = c;
            else
                C(:,:,i,j) = creg;
            end
        end
    end
    
    
    
    % LARK
    r0 = ones(2,1)*fix(q/2)+1;
    W = zeros(q*q,length(2+t:q:N-q-1)*length(2+t:q:M-q-1));
    
    g = 0;
    i0 = 0;
    indr = [];
    for i=2+t:q:N-q-1
        j0 = 0;
        i0 = i0+1;
        for j=2+t:q:M-q-1
            j0 = j0+1;
            % y = I(i:i+q-1,j:j+q-1);
            c = C(:,:,i:i+q-1,j:j+q-1);
            s = zeros(q,q); % LARK similarity
            for a=1:q
                for b=1:q
                    dx = [a b]' - r0;
                    cl = c(:,:,a,b);
                    % s(a,b) = exp(-(dx'*dx)/2/hd2);
                    s(a,b) = exp(-(dx'*cl*dx)/2/hd2);
                    % s(a,b) = sqrt(det(cl))/2/pi/hd2*exp(-(dx'*cl*dx)/2/hd2);
                    % p(a,b) = sqrt(det(cl));
                end
            end
            S(i:i+q-1,j:j+q-1) = s;% /sum2(s);
            g = g+1;
            W(:,g) = s(:);
            indr = [indr; g i0 j0];
        end
        if show
            imshow(S*64,jet)
            pause(0)
        end
    end
    %imax = i0;
    %jmax = j0;
    ind = [i0 j0];
    if show
        imshow(S*64,jet)
        pause(0)
    end
    
    
else %3D
    la1 = options.la1;
    la2 = options.la2;
    if show
        figure
        disp('Bfx_lark: computing 3D-LARK sequence descriptors...');
    end
    
    if length(q)==1
        qi = q;
        qj = q;
        qk = q;
    else
        qi = q(1);
        qj = q(2);
        qk = q(3);
    end
    
    
    
    % gradient
    
    Gx = zeros(size(I));
    Gy = zeros(size(I));
    Gt = zeros(size(I));
    
    
    for i=1:N-1
        for j=1:M-1
            for k=1:P-1
                Gx(i,j,k) = I(i,j,k)-I(i+1,j,k);
                Gy(i,j,k) = I(i,j,k)-I(i,j+1,k);
                Gt(i,j,k) = I(i,j,k)-I(i,j,k+1);
            end
        end
    end
    % covariance matrices
    C = zeros(3,3,N,M,P);
    p   = (2*t+1)^3;
    rho = zeros(3,1);
    for k=2+t:P-t-1
        k
        for i=2+t:N-t-1
            for j=2+t:M-t-1
                gx      = Gx(i-t:i+t,j-t:j+t,k-t:k+t);
                gy      = Gy(i-t:i+t,j-t:j+t,k-t:k+t);
                gt      = Gt(i-t:i+t,j-t:j+t,k-t:k+t);
                J       = [gx(:) gy(:) gt(:)];
                c       = J'*J;
                [UU,SS,VV] = svd(J);
                rho(1) = (SS(1,1)+la1)/((SS(2,2)*SS(3,3))+la1);
                rho(2) = (SS(2,2)+la1)/((SS(1,1)*SS(3,3))+la1);
                rho(3) = (SS(3,3)+la1)/((SS(1,1)*SS(2,2))+la1);
                ga     = ((SS(1,1)*SS(2,2)*SS(3,3)+la2)/p)^alpha;
                creg   = p*ga*(rho(1)*VV(:,1)*VV(:,1)'+rho(2)*VV(:,2)*VV(:,2)'+rho(3)*VV(:,3)*VV(:,3)');
                x = sum(isnan(creg(:))==1);
                if x>1
                    C(:,:,i,j,k) = c;
                else
                    C(:,:,i,j,k) = creg;
                end
            end
        end
    end
    
    
    
    % LARK
    r0 = fix([qi qj qk]'/2)+1;
    W = zeros(qi*qj*qk,length(2+t:qi:N-qi-1)*length(2+t:qj:M-qj-1)*length(2+t:qk:P-qk-1));
    indr = [];
    hd23 = sqrt((2*pi*hd2)^3);
    g = 0;
    i0 = 0;
    for i=2+t:qi:N-qi-1
        i
        i0 = i0+1;
        
        j0 = 0;
        for j=2+t:qj:M-qj-1
            j0 = j0+1;
            
            k0 = 0;
            for k=2+t:qk:P-qk-1
                k0 = k0+1;
                
                % y = I(i:i+q-1,j:j+q-1);
                creg = C(:,:,i:i+qi-1,j:j+qj-1,k:k+qk-1);
                s = zeros(qi,qj,qk); % LARK similarity
                for a=1:qi
                    for b=1:qj
                        for c=1:qk
                            dx = [a b c]' - r0;
                            cl = creg(:,:,a,b,c);
                            % cl = eye(3,3);
                            % s(a,b,c) = exp(-(dx'*dx)/2/hd2);
                            % s(a,b,c) = exp(-(dx'*cl*dx)/2/hd2);
                            s(a,b,c) = sqrt(det(cl))/hd23*exp(-(dx'*cl*dx)/2/hd2);
                            % p(a,b,c) = sqrt(det(cl));
                        end
                    end
                end
                S(i:i+qi-1,j:j+qj-1,k:k+qk-1) = s/sum(abs(s(:)));
                g = g+1;
                W(:,g) = s(:);
            indr = [indr; g i0 j0 k0];
            end
        end
    end
        ind = [i0 j0 k0];

    
    %[N,M] = size(S);
    %S0 = zeros(size(I));
    %S0(1+t:t+N,1+t:t+M) = S;
    %S = S0;
    
    
    if show
        for k=1:P
            imshow(S(:,:,k)*64,jet);
            title(num2str(k));
            pause(0.5)
        end
    end
end