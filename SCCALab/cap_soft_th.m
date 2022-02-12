function t = cap_soft_th(d, r, tol)
amin = 0;
amax = max(d);

while 1
    amid = (amin + amax) / 2;
    t = min(1, max(0, d-amid));
    if sum(t) <= r
        amax = amid;
    else
        amin = amid;
    end
    if amax - amin < tol^2
        break
    end
end