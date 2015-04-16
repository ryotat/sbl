function [X,Z,Psi,fval]=sblmatcomp(m, n, ind, yy, lambda)

maxiter=50;

[I,J]=ind2sub([m,n],ind);

N=length(ind);

Oind=zeros(N,N);
for ii=1:N
  for jj=1:N
    if J(ii)==J(jj)
      Oind(ii,jj)=I(ii)+(I(jj)-1)*m;
    else
      Oind(ii,jj)=nan;
    end
  end
end
ival=~isnan(Oind);

X=zeros(m,n);
Psi=eye(m,m);
M=zeros(m,n);
Psi1=zeros(N,N);
for kk=1:maxiter
  % update X
  Psi1(ival)=Psi(Oind(ival));
  if ~isequal(Psi1,Psi1')
    keyboard;
  end
  T=inv(Psi1+lambda*eye(N));
  M(ind)=T*yy;
  X=Psi*M;
  if ~isreal(X)
    keyboard;
  end
  
  
  % update Z
  Z=zeros(m,m);
  for jj=1:m
    ix=find(J==jj);
    Ijj=I(ix);
    Z(Ijj,Ijj)=Z(Ijj,Ijj)+T(ix,ix);
  end
  if ~isreal(Z)
    keyboard;
  end

  % update Psi
  Zh=sqrtm(Z);
  ZhX=Zh*X;
  Psi=Zh\(sqrtm(ZhX*ZhX')/Zh);
  Psi=(Psi+Psi')/2;
  Psi=real(Psi);
  Psi1(ival)=Psi(Oind(ival));

  % report
  err=norm(yy-X(ind))^2;
  ee=eig(Psi1+lambda*eye(N));
  fval(kk)=err/lambda+trace(X'*(Psi\X))+sum(log(ee));
  ss=svd(X);
  nc=sum(ss>=0.01*max(ss));
  fprintf('kk=%d fval=%g err=%g rank=%d\n', kk, fval(kk), err, nc);
end