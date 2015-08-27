function [b,Lambda]=matrix4(x)

global MyLAMBDA MyMATRIX

  if isempty(x), [b,Lambda]=makematrix(50); return, end

  b=MyMATRIX*x; 

return
%====================================================
function [n,Lambda]=makematrix(n)

global MyLAMBDA MyMATRIX

  A=2*rand(n,n)-1; % A=A'+A;
  Lambda=eig(A);
  MyMATRIX=A;

 
return
