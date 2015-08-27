function [b,Lambda]=matrix2(x)

global MyLAMBDA MyMATRIX V_HOUSEHOLDER

  if isempty(x), [b,Lambda]=makematrix(1000); return, end

  V=V_HOUSEHOLDER;
  x=x-V*(V'*(2*x));
  b=MyMATRIX*x; 
  b=b-V*(V'*(2*b));

return
%====================================================
function [n,Lambda]=makematrix(n)

global MyLAMBDA MyMATRIX V_HOUSEHOLDER

  Lambda=[1,1000];
  h=(Lambda(2)-Lambda(1))/n; 

  

  MyLAMBDA=(Lambda(1):h:Lambda(2))';
  n=length(MyLAMBDA);
 
  shift=MyLAMBDA(3,1)-1.0e-10;
  MyMATRIX=spdiags(MyLAMBDA-shift,0,n,n);

  %%% diagonal matrices may give the wrong impression 
  %%% of the eefects of rounding errors
  V=rand(n,1); [V,R]=qr(V,0);
  V_HOUSEHOLDER=V;

return
