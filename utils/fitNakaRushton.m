function [f,c50,rmax,n,offset] = fitNakaRushton(c,r, makePlot)

if ~exist('makePlot', 'var') || isempty(makePlot)
    makePlot = false;
end

% find contrast that evokes closest to half-maximal response
rMid = ((max(r)-min(r))/2) + min(r);
[~,rMidIndex] = min(abs(r-rMid));
initC50 = c(rMidIndex(1));

% Set starting points and bounds:

% c50 n offset rmax
sp = [initC50 2 0 rMid];
lb = [0 1 0 -inf];
ub = [max(c) 5 0 max(r)];
% Note that we fix the offset to be zero, forcing the function to go
% through 0.

% Define function and algorithm
searchopts = optimoptions(@lsqcurvefit,'Algorithm','levenberg-marquardt', 'Display', 'off');
fun = @(x,c)x(4)* ((c.^x(2))./((c.^x(2))+x(1).^x(2)))+x(3);

% Fit function
x = lsqcurvefit(fun,sp,c,r,lb,ub,searchopts);
 
% Save out parameters
c50 = x(1);
n = x(2);
offset = x(3);
rmax = x(4);
f = fun;

% Plot
if makePlot
    c1 = 0:max(c)/1000:max(c)*2;
    y1 = fun(x,c1);
    figure;hold on
    plot(c,r, c1, y1);%pause(1);%close;
    line([c50 c50], [0 rmax],'LineStyle', ':', 'Color', 'k');
end



end

%%% from  http://gru.stanford.edu/svn/matlab/fitNakaRushton.m
% % parmaeters
%              %Rmax          c50     n     offsets x 5
% initParams = [max(r)        initC50 2  min(r)];
% minParams =  [0             0       1  -inf];
% maxParams =  [inf           1       5  inf];

% % OLD (using different fitting function, yielding similar results)
% formula_to_fit = 'rmax * ((x.^n)./((x.^n)+c50.^n))+o';
% f = fit(c, r, formula_to_fit,  'StartPoint', sp, 'Lower', lb, 'Upper', ub);
% c50 = f.c50;
% rmax = f.rmax;
% n = f.n;
% offset = f.o;

