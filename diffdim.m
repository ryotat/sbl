function dimcot = diffdim(A,r,method)
%diffdim         calculates the dimension of the differentials of the mask A
%                 in rank r
%
%usage
%  dimcot = diffdim(A,r)
%
%input
%  A              (m,n)-adjacency matrix of the bipartite graph
%  r              the rank of matrix A
%
%output 
%  dimcot         the dimension of the differential
%
%author
%  franz.j.kiraly@tu-berlin.de

A=A~=0;

[m,n]=size(A);

randU=randn(m,r);
randV=randn(r,n);

diff = zeros(sum(A(:)),r*(m+n));

[I,J]=find(A);
for kk=1:length(I)
  i=I(kk);
  j=J(kk);
  Iu=(i-1)*r+(1:r);
  Iv=m*r+(j-1)*r+(1:r);
  diff(kk,Iu) = randV(:,j)';
  diff(kk,Iv) = randU(i,:);
end

if ~exist('method','var') || isequal(method,'svd')
  dimcot = rank(diff,1e-12);
elseif isequal(method,'qr')
  [Q,R]=qr(diff,0);
  ind=find(abs(diag(R))>1e-12);
  dimcot=max(ind)+rank(R(max(ind)+1:end,:));
else
  error('Unknown method [%s]\n', method);
end

