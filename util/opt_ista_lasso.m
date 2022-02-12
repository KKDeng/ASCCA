
function x = opt_ista_lasso(A,b,lama,x0)
%dim=size(A);
%x0=zeros(dim(2),1);
evs = eig(A);
L = max(evs);t0=1;y=x0;
while(1)
    z=y-1/L*(A*y-b);
    x=subplus(abs(z)-lama/L) .* sign(z);
    t=(1+sqrt(1+4*t0^2))/2;
    y=x+(t0-1)/t*(x-x0);
    if(norm(x-x0)<1e-3)
        break;
    end
    t0=t;
    x0=x;
end