function y = linop_crosprod( p1, p2, SXroot, SYroot, A, mode )

% disp(size(A));
% disp(size(SXroot));
% disp(size(SYroot));

switch mode
    case 0 
        y = { [p1, p2], [p1, p2] };
    case 1
        y = (SXroot * A) * SYroot;
    case 2
        y = (SXroot * A) * SYroot;
end


% SCCALab v0.1 by Zongming Ma.
% All rights reserved.
