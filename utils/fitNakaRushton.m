function [f,c50,rmax,n,offset] = fitNakaRushton(c,r, makePlot)

if ~exist('makePlot', 'var') || isempty(makePlot)
    makePlot = false;
end

% find contrast that evokes closest to half-maximal response
rMid = ((max(r)-min(r))/2) + min(r);
[~,rMidIndex] = min(abs(r-rMid));
initC50 = c(rMidIndex(1));

formula_to_fit = 'rmax * ((x.^n)./((x.^n)+c50.^n))+o';

% c50 n offset rmax
sp = [initC50 2 min(r) max(r)];
lb = [0 1 min(r)-0.1*min(r) 0];
ub = [max(c) 3 rMid max(r)+0.1*max(r)];

f = fit(c, r, formula_to_fit,  'StartPoint', sp, 'Lower', lb, 'Upper', ub);

if makePlot
    x1 = 0:max(c)/1000:max(c);
    y1 = f(x1);
    figure;hold on
    plot(c,r, x1, y1);pause(1);close;
end

c50 = f.c50;
rmax = f.rmax;
n = f.n;
offset = f.o;

end

%%% from  http://gru.stanford.edu/svn/matlab/fitNakaRushton.m
% % parmaeters
%              %Rmax          c50     n     offsets x 5
% initParams = [max(r)        initC50 2  min(r)];
% minParams =  [0             0       1  -inf];
% maxParams =  [inf           1       5  inf];
