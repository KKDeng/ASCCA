function perf(T,logplot)
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

if (nargin < 2) logplot = 0; end

colors  = ['k' 'b' 'r' 'g' 'c' 'k' 'y'];
lines   = ['-' '--' '-.' ':'];
markers = ['x' '*' 's' 'd' 'v' '^' 'o'];

[np,ns] = size(T);

% Minimal performance per solver

minperf = min(T,[],2);

% Compute ratios and divide by smallest element in each row.

r = zeros(np,ns);
for p = 1: np
  r(p,:) = T(p,:)/minperf(p);
end

if (logplot) r = log2(r); end

max_ratio = max(max(r));

% Replace all NaN's with twice the max_ratio and sort.

r(find(isnan(r))) = 2*max_ratio;
r = sort(r);

% Plot stair graphs with markers.

clf;
for s = 1: ns
 [xs(:,s),ys(:,s)] = stairs(r(:,s),[1:np]/np);

% [xs , ys] = stairs(r(:,s),[1:np]/np);

% option = [lines(s) colors(s) markers(s)];
% plot(xs,ys,option,'MarkerSize',2);
% hold on;
end

%plot(xs(:,1), ys(:,1), 'r-', xs(:,2), ys(:,2), 'b-.', xs(:,3), ys(:,3), 'k--');

plot(xs(:,1), ys(:,1), 'r-', xs(:,2), ys(:,2), 'b-.');

% Axis properties are set so that failures are not shown,
% but with the max_ratio data points shown. This highlights
% the "flatline" effect.

axis([ 0 1.1*max_ratio 0 1 ]);
% Legends and title should be added.

legend('ipoptpen', 'ipoptfilt', 'pennlp');
%legend('ipoptpen', 'ipoptfilt');
title('Iterations');
xlabel('iteration');

%axis([ 0 1000 0.4 1 ]);
axis([ 0 10 0.4 1 ]);


