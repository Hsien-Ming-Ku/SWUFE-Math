function [A,b,hbtype]=myhbo(matrix,mdir,dirLOCAL)
%MakeHB(MATRIX)
%  transfers the matrix /data/harwell-boeing/MATRIX.rua.gz
%  in the harwell-boeing collection to MATLAB format,
%  puts it in /tmp/MATRIX.mat, and extends the MATLABPATH
%  with /tmp and /data/harwell-boeing/
%
%MakeHB(MATRIX,MDIR)
%  transfers the matrix /data/harwell-boeing/MDIR/MATRIX.rua.gz
%  in the harwell-boeing collection to MATLAB format 
%  and puts it in /tmp/MATRIX.mat and extends the MATLABPATH
%  with /tmp and /data/harwell-boeing
%
%MakeHB(MATRIX,MDIR,dirLOCAL)
%  transfers the matrix /data/harwell-boeing/MDIR/MATRIX.rua.gz
%  in the harwell-boeing collection to MATLAB format 
%  and puts it in /dirLOCAL/MATRIX.mat and extends the MATLABPATH
%  with /dirLOCAL and /data/harwell-boeing

%   Copyright (c) 98
%   Gerard Sleijpen.


if ~exist([matrix,'.mat'])

dirHB='/data/harwell_boeing/'; 
if nargin==3,dirLOCAL=[dirLOCAL,'/']; else, dirLOCAL='/tmp/'; end
eval(sprintf('addpath %s %s',dirLOCAL,dirHB))


dirHB1=dirHB;
if nargin>1 & ~isempty(mdir), dirHB1=[dirHB1,mdir,'/']; end


if ~exist([matrix,'.mat'])

  ext=extension(dirHB1,matrix);
  if isempty(ext)
    dirHB1='/scratch/Sleijpen/Harwell_Boeing/';
    if nargin>1 & ~isempty(mdir), dirHB1=[dirHB1,mdir,'/']; end
    ext=extension(dirHB1,matrix);
  end

  if ~isempty(ext)
    eval(sprintf('! cp %s%s.%s.gz %s.',dirHB1,matrix,ext,dirLOCAL))
    eval(sprintf('! gunzip %s%s.%s.gz',dirLOCAL,matrix,ext))
    if ~strcmp(ext,'mat')
      eval(sprintf('! cd %s; %shbo2mat %s.%s',dirLOCAL,dirHB,matrix,ext))
      eval(sprintf('! rm %s%s.%s',dirLOCAL,matrix,ext))
    end
  else
     msg=['  can not find the matrix   ',matrix];
     error(msg)
  end

end
end


eval('[A,b,hbtype]=hbo(matrix);')

return
%--------------------------------------------------------------------------
function ext=extension(direc,matrix)

  ext_options=['pse';'rsa';'rza';'rua';'rra';'csa';'psa';'pua';'pra';'mat'];

  ext=[];
  for j=1:size(ext_options,1)
    ext0=ext_options(j,:);
    check_matrix=[direc,matrix,'.',ext0,'.gz'];
    if exist(check_matrix), ext=ext0; end
  end

return


% PSE  - Pattern symmetric unassembled
% RSA  - Real symmetric
% RZA  - Real skew symmetric
% RUA  - Real unsymmetric
% RRA  - Real rectangular
% CSA  - Complex symmetric
% PSA  - Pattern symmetric
% PUA  - Pattern unsymmetric
% PRA  - Pattern rectangular
