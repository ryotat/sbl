function M=rand_mask_regular(n,d)
%rand_mask_regular     generate a n x n mask with d entries per row and
%column
%
%usage
%  mask = rank_mask_regular(n,d)
%
%input
%  n           matrix size 
%  d           number of ones per row and column
%
%output
%  M          d-regular mask
%
%reference
%  Kim-Vu algorithm
%
%author
%  Ryota Tomioka <tomioka@mist.i.u-tokyo.ac.jp>

M=sparse([],[],[],n,n,n*d);
while ~isequal(unique(sum(M)),d) || ~isequal(unique(sum(M')),d)
  fprintf('Sampling...\n');
  ind1=repmat(1:n, [1,d]);
  ind2=ind1; ind2=ind2(randperm(n*d));
  
  while length(ind1)>0
    ix1=randi(length(ind1),1); ii=ind1(ix1);
    ix2=randi(length(ind2),1); jj=ind2(ix2);

    if ~M(ii,jj)
      M(ii,jj)=1;
      ind1(ix1)=[];
      ind2(ix2)=[];
    else
      if all(all(M(ind1,ind2)))
        fprintf('Giving up at length=%g\n',length(ind1));
        M=sparse([],[],[],n,n,n*d);
        break;
      end
    end
  end

end

% Original code
% ind1=repmat(1:n, [1,d]);

%check=0;
%while ~check
%  ind2=ind1; ind2=ind2(randperm(n*d));
%  M=sparse(ind1, ind2, ones(1,n*d), n, n);
%  
%  check=max(M(:))==1;
%end
