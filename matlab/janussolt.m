function z=janussolt(PREC,y)
% z=janussolt(PREC,y)
%
% solve system with given block ILU , (SL^{-1} Q LDU P^T SR^{-1})^T z = y
%
% Input
% -----
% PREC    structure with factorization obtained from janus
% y       right hand side(s)
%
% Output
% ------
% z       computed solution


%    Authors:
%
%	Matthias Bollhoefer, TU Braunschweig
%
%    Date:
%
%	August 28, 2017. JANUS Block ILU R1.0.  
%
%    Notice:
%
%	Copyright (c) 2017 by TU Braunschweig.  All Rights Reserved.
%
%	THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
%	EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
%
%    Availability:
%
%	This file is located at
%
%	http://bilu.tu-bs.de/
% $Id: janussolt.m 4135 2018-04-18 15:34:46Z bolle $

if isfield(PREC,'pivots')
   pivots=PREC.pivots;
else
   pivots=[];
end   
 
% Cholesky does not need pivoting 
if PREC.isdefinite && PREC.ishermitian & PREC.invert_blocks==0
   pivots=-(PREC.n+1);
end
 
if PREC.isreal
   % real symmetric case
   if PREC.issymmetric
      z=PREC.SL*(PREC.P*bildlsol(PREC.BL,PREC.BiD,pivots,PREC.P'*(PREC.SL*y)));
      % z=bildlsol(PREC.BL,PREC.BiD,y);
   % general real case
   else
      z=PREC.SR*(PREC.P*bilusolt(PREC.BL,PREC.BiD,PREC.BUT,pivots,PREC.Q'*(PREC.SL*y)));
      % z=bilusol(PREC.BL,PREC.BiD,PREC.BUT,y);
   end
else % non-real case
   % non-real complex-symmetric case
   if PREC.issymmetric
      z=PREC.SL*(PREC.P*bildlsols(PREC.BL,PREC.BiD,pivots,PREC.P'*(PREC.SL*y)));
      % z=bildlsols(PREC.BL,PREC.BiD,y);
   % non-real hermitian case
   elseif PREC.hermitian
      z=conj(PREC.SL*(PREC.P*bildlsol(PREC.BL,PREC.BiD,pivots,PREC.P'*(PREC.SL*conj(y)))));
      % z=bildlsol(PREC.BL,PREC.BiD,y);
   % general non-real case
   else
      z=PREC.SR*(PREC.P*bilusolt(PREC.BL,PREC.BiD,PREC.BUT,pivots,PREC.Q'*(PREC.SL*y)));
      % z=bilusol(PREC.BL,PREC.BiD,PREC.BUT,y);
   end
end
