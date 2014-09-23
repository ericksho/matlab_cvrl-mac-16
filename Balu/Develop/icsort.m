% It returns the k smallest values of x using 
% Paredes & Navarro. "Optimal incremental sorting.", Proc. 8th Workshop on 
% Algorithm Engineering and Experiments and 3rd Workshop on Analytic 
% Algorithmics and Combinatorics (ALENEX-ANALCO?06). 2006.

function s = icsort(x,k)

ok = 0;

s = [];

while not(ok)
    p = x(1);
    i0 = find(x<p);
    if (length(s)+length(i0))<=k
        s = [s x(x<=p)];
        if length(s)>=k
            ok = 1;
        end
        x1 = x(x>p);
        if ~isempty(x1)
            x = x1;
        else
            x = x(x<=p);
         end
    else
        x1 = x(x<p);
        if ~isempty(x1)
            x = x1;
        else
            x = min(x);
        end
    end
end

%j = sort(s);
length(s)
s = s(1:k);
