function [b,Lambda]=matrix1(x,flag)

global MyLAMBDA MyMATRIX

  if isempty(x), [b,Lambda]=makematrix(1000); return, end

  % b=MyMATRIX*x; 
  b=MyLAMBDA.*x;

return
%====================================================
function [n,Lambda]=makematrix(n)

global MyLAMBDA MyMATRIX
 
  Lambda=[1,1000];
  h=(Lambda(2)-Lambda(1))/n; 

  MyLAMBDA=(Lambda(1):h:Lambda(2))';
  MyLAMBDA=(sqrt(Lambda(1):h:Lambda(2)))';
  n=length(MyLAMBDA);
 
  % MyMATRIX=diag(MyLAMBDA);
  % MyMATRIX=spdiags(MyLAMBDA,0,n,n);
 
return
