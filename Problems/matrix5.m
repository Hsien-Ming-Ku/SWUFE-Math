function [b,Lambda]=matrix5(x)

global MyLAMBDA MyMATRIX

  if isempty(x), [b,Lambda]=makematrix(100); return, end

  b=MyMATRIX*x; 

return
%====================================================
function [n,Lambda]=makematrix(n)

global MyLAMBDA MyMATRIX

  A=2*rand(n,n)-1; A=A'+A;
  Lambda=eig(A);
  MyMATRIX=A;

 
return
