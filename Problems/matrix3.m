function [b,B]=matrix3(x,b)

global MyMATRIX

  if isempty(x), 
    makematrix(b); return, 
  end

  b=MyMATRIX*x;
 
return
%====================================================
function makematrix(n)

global MyMATRIX

   e = ones(n,1);
   MyMATRIX = spdiags([e -2*e e], -1:1, n, n);
   % full(MyMATRIX)


return 
