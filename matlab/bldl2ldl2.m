function [L,D]=bldl2ldl2(BL,BD)
% 1) [L,D]=bldl2ldl2(BL,BiD), [L,D]=bldl2ldl2(P)
%
% 1) transform block BL BiD^{-1} BL^T decomposition into scalar LDL^T factorization


%    Authors:
%
%	Matthias Bollhoefer, TU Braunschweig
%
%    Date:
%
%	December 25, 2016. JANUS Block ILU R1.0.  
%
%    Notice:
%
%	Copyright (c) 2016 by TU Braunschweig.  All Rights Reserved.
%
%	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
%	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
%
%    Availability:
%
%	This file is located at
%
%	http://bilu.tu-bs.de/

if nargin==1
   BD=BL.BiD;
   BL=BL.BL;
end    
   n=max(BL{end}.J);
   L=speye(n);
   D=speye(n);
   for k=1:length(BL)
       I=BL{k}.I;
       J=BL{k}.J;
       BD_loc=BD{k}.D;
       BD_loc=inv(tril(BD_loc)+tril(BD_loc,-1)');
       LL=chol(BD_loc,'lower'); DD=diag(diag(LL)); LL=LL/DD; DD=DD*DD;
       
       D(J,J)=DD;
       L(J,J)=LL;
       L(I,J)=BL{k}.L*LL;
   end % for k
