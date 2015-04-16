function X=randlowrankmat(sz, rr)

[U,R]=qr(randn(sz(1),rr),0);
[V,R]=qr(randn(sz(2),rr),0);

X=U*randn(rr)*V';

