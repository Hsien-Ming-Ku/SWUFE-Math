function [b,B,hbtype]=matrixhbo(x,str)

global MyMATRIX

  if isempty(x), 
    if nargin==1, str=''; end, 
    [b,B,hbtype]=makematrix(str); return, 
  end

  
  if nargin==2, b=MyMATRIX'*x; return; end
  
  b=MyMATRIX*x;
 
return
%====================================================
function [n,b,hbtype]=makematrix(str)

global MyMATRIX

matrix='meier01';
if (~isempty(str)), matrix=str; end

[MyMATRIX,xy,b,hbtype,rhstype]=hbo(matrix);
[n,m]=size(MyMATRIX);

return 
