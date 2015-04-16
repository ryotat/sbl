n=100;
r=5;

while 1
  M=rand_mask_regular(n,2*r);
  if diffdim(M, r)>=r*(2*n-r) && rank(stressmatgen(M,r))==n-r;
    break;
  end
end

ind=find(M);
[I,J]=ind2sub([n,n],ind);

% generate truth
X0=randlowrankmat([100,100],5);

% convex
X_cvx=matrix_adm(zeros(n),{I,J}, X0(ind), 0);

err_cvx = norm(X_cvx(:)-X0(:))^2;

% SBL
[X_sbl,Z,Psi,fval]=sblmatcomp(n, n, ind, X0(ind), 0.001);

err_sbl = norm(X_sbl(:)-X0(:))^2;

figure
subplot(1,2,1);
plot(X_cvx(:), X0(:), 'x');
grid on;
title(sprintf('Convex (err=%g)', err_cvx));

subplot(1,2,2);
plot(X_sbl(:), X0(:), 'x');
grid on;
title(sprintf('SBl (err=%g)', err_sbl));



figure, plot(fval);