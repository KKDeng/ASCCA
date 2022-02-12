function perf(T,logplot,a)
%PERF    Performace profiles
%
% PERF(T,logplot)-- produces a performace profile as described in
%   Benchmarking optimization software with performance profiles,
%   E.D. Dolan and J.J. More',
%   Mathematical Programming, 91 (2002), 201--213.
% Each column of the matrix T defines the performance data for a solver.
% Failures on a given problem are represented by a NaN.
% The optional argument logplot is used to produce a
% log (base 2) performance plot.
%
% This function is based on the perl script of Liz Dolan.
%
% Jorge J. More', June 2004

if (nargin < 3) logplot = 0; end

colors  = ['m' 'b' 'r' 'c' 'k' 'g' 'y'];
lines   = [':' '-' '-.' '--','-','--',':'];
markers = ['x' '*' 's' 'd' 'v' '^' 'yo'];
lines   = {'k-' 'b--' 'r-.' 'c:' 'g-*' 'm-' 'r-'};

[np,ns] = size(T);

% Minimal performance per solver

minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
    r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log(r)/log(a); end

max_ratio = max(max(r));
if max_ratio < 1; max_ratio = 1; end

% Replace all NaN's with twice the max_ratio and sort.

r(find(isnan(r))) = 2*max_ratio;
r = sort(r);

xlim =  max_ratio*1.1;
if xlim > 10; xlim = 10; end
if xlim < 1e-6;  xlim = 1e-6; end
% Plot stair graphs with markers.

clf;
% for s = 1: ns
%  [xs,ys] = stairs(r(:,s),[1:np]/np);
%  option = ['-' colors(s) markers(s)];
% %  plot([xs;max_ratio],[ys;1],option,'MarkerSize',3);
%  plot([xs;max_ratio],[ys;1],lines{s},'LineWidth',2,'MarkerSize',3);
% %  plot([r(:,s);max_ratio],[[1:np]'/np;1],option,'MarkerSize',3);
%  hold on;
% end

for s = 1: ns
    [xs,ys] = stairs(r(:,s),[1:np]/np);
    option = ['-' colors(s) markers(s)];
    %  plot([xs;max_ratio],[ys;1],option,'MarkerSize',3);
    plot([xs;max_ratio],[ys;1],lines{s},'LineWidth',2,'MarkerSize',3);
%     plot([xs;max_ratio],[ys;1], option,'LineWidth',2,'MarkerSize',3);
%  plot([r(:,s);max_ratio],[[1:np]'/np;1],option,'MarkerSize',3);
    hold on;
end

% for s = 5: ns
%     [xs,ys] = stairs(r(:,s),[1:np]/np);
%     option = ['-' colors(s) markers(s)];
%     %  plot([xs;max_ratio],[ys;1],option,'MarkerSize',3);
%     plot([xs;max_ratio],[ys;1],lines{s},'LineWidth',1,'MarkerSize',3);
%     %  plot([r(:,s);max_ratio],[[1:np]'/np;1],option,'MarkerSize',3);
%     hold on;
% end

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

if logplot
    axis([0 xlim 0 1.05 ]);
else
    axis([ 1 xlim 0 1.05 ]);
end

% Legends and title should be added.



