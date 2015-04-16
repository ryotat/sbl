function genstressmatrix = stressmatgen(M,r,U,V)
%stresses   outputs a generic stress matrix for the mask M
%
%usage
%  genstressmatrix = stressmatgen(M,r)
%
%input
%  M                 (m,n)-adjacency matrix of the bipartite graph
%  r                 the true rank
%optional
%  U              (m,r)-matrix for left singular vectors
%  V              (n,r)-matrix of right singular vectors
%              (U and V are initialized random if not supplied)
%
%
%output 
%  genstressmatrix   generic stress matrix for the mask A in rank r
%
%author
%  theran@math.fu-berlin.de
%  franz.j.kiraly@tu-berlin.de

M = M ~= 0;
if nargin < 3
    stressvecs = stresses(M,r);
else
    stressvecs = stresses(M,r,U,V);
end

[m,n] = size(M);
[~,numstresses] = size(stressvecs);

IJ = find(M);

randlin = randn(numstresses,1);

stressvec = stressvecs*randlin;


omega = zeros(m,n);
omega(IJ) = stressvec;
genstressmatrix = omega;

end % end of function
