function [L,D,U]=bldu2ldu2(BL,BD,BUT)
% 1) [L,D,U]=bldu2ldu2(BL,BiD,BUT), [L,D,U]=bldu2ldu2(P)
%
% 1) transform block BL BiD^{-1} BUT^T decomposition into scalar LDU factorization


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
   BUT=BL.BUT;
   BD =BL.BiD;
   BL =BL.BL;
end    
   n=max(BL{end}.J);
   L=speye(n);
   D=speye(n);
   U=speye(n);
   for k=1:length(BL)
       Ic=BL{k}.I;
       Ir=BUT{k}.I;
       J=BL{k}.J;
       BD_loc=inv(BD{k}.D);
       [LL,UU]=lu(BD_loc,'lower'); DD=diag(diag(UU)); LL=LL/DD; DD=DD*DD;
       
       D(J,J)=DD;
       L(J,J)=LL;
       L(Ic,J)=BL{k}.L*LL;
       U(J,Ir)=UU'*BUT{k}.L.';
   end % for k
