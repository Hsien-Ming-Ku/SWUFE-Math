function [b,Lambda]=matrix1s(x)

global MyLAMBDA MyMATRIX

  if isempty(x), [b,Lambda]=makematrix(999); return, end

  % b=MyMATRIX*x; 
  b=MyLAMBDA.*x;

return
%====================================================
function [n,Lambda]=makematrix(n)

global MyLAMBDA MyMATRIX

  Lambda=[1,1000];
  h=(Lambda(2)-Lambda(1))/n; 

  shift=0;
  shift=2.001; shift=Lambda(1)-10*eps;
  % shift=8.993537413297;
  MyLAMBDA=(Lambda(1):h:Lambda(2))'-shift;
  % MyLAMBDA=(sqrt(MyLAMBDA))- shift;
  n=length(MyLAMBDA);
 
  % MyMATRIX=diag(MyLAMBDA);
  % MyMATRIX=spdiags(MyLAMBDA,0,n,n);
 
return
