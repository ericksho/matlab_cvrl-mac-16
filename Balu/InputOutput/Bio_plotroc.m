% Bio_edgeview(B,E,c,g)
%
% Toolbox: Balu
%   Plot ROC curve and fit to an exponetial curve
%
% Example:
%
% th = 3;x = 0:0.05:1; y = 1-exp(-3*x)+randn(1,21)*0.05;
% Bio_plotroc(x,y)
%
% D.Mery, PUC-DCC, Apr. 2013
% http://dmery.ing.puc.cl
%

function Bio_plotroc(x,y)


clf
plot(x,y,'r*')

ths = fminsearch(@thest,1,[],x,y);
xs = 0:0.05:1;
ys = 1-exp(-ths*xs);

hold on
plot(xs,ys)

end



function err = thest(th,x,y)
ys = 1-exp(-th*x);
err = norm(y-ys);
end




