function stressbasis = stresses(M,r,U,V)
%stresses   gives a basis for the stress space of the mask M
%
%usage
%  stressbasis = stresses(M,r)
%
%input
%  M              (m,n)-adjacency matrix of the bipartite graph
%  r              target rank
%optional
%  U              (m,r)-matrix for left singular vectors
%  V              (n,r)-matrix of right singular vectors
%              (U and V are initialized random if not supplied)
%
%output 
%  stressbasis   basis for the space of stresses of the mask M
%
%author
%  theran@math.fu-berlin.de
%  franz.j.kiraly@tu-berlin.de

M = M ~= 0;

[m,n] = size(M);


if nargin < 3
    
    randU=randn(m,r);
    randV=randn(n,r);
    
else
    
    randU = U;
    randV = V;
    
end


diff = zeros(sum(M(:)),r*(m+n));

[I,J]=find(M);
for kk=1:length(I)
  i=I(kk);
  j=J(kk);
  Iu=(i-1)*r+(1:r);
  Iv=m*r+(j-1)*r+(1:r);
  diff(kk,Iu) = randV(j,:);
  diff(kk,Iv) = randU(i,:);
end

stressbasis = null(diff');

end % end of function