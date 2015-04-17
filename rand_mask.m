function A = rand_mask(m,n,p)
%rand_mask     outputs a random m x n mask with p known entries
%               uniformly drawn from all such masks
%
%usage
%  A = rand_mask(m,n,p)
%
%input
%  m,n            number of rows, columns             
%  p              number of known entries
%
%output 
%  A              the m x n mask with p known entries
%
%author
%  franz.j.kiraly@tu-berlin.de

mnswitched = 0;
pswitched = 0;

if m > n
    nt = n;
    n = m;
    m = nt;
    mnswitched = 1;
end

if p > m*n/2
   p = m*n - p;
   pswitched = 1; 
end

if p >= m*n
  A= ones(m,n);
else
freetot = m*n;
freerow = n*(1:m);
A = zeros(m,n);

for k=1:p
    indi=randi(freetot);
    indj=randi(freetot);
    i = find(freerow >= indi,1);
    freeinrow = find(A(i,:) == 0);
    j = freeinrow(randi(length(freeinrow)));
    freerow(i:m)=freerow(i:m)-1;
    freetot = freetot -1;
    A(i,j)=1;
end
end

if pswitched == 1
   A = ones(m,n) - A; 
end

if mnswitched == 1
A = A';    
end
end