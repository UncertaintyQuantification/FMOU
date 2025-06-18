function S=subtractionmatrix(n,i,d)
%SubtractionMatrix    S=subtractionmatrix(numstats,refstat,dimension)
%
%Makes a matrix that when multiplied by a displacement vector subtracts the reference
%station from all the others.  Optional 'dimension' defaults to 3.

if nargin < 3
   d = 3;
end

S=speye((n-1)*d,n*d);
j=(i-1)*d+1;
S(:,j+d:end)=S(:,j:end-d);
S(:,j:j+d-1)=repmat(-eye(d),n-1,1);
