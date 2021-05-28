function z=janussol(PREC,y)
% z=janussol(PREC,y)
%
% solve system with given block ILU , (SL^{-1} Q LDU P^T SR^{-1}) z = y
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
% $Id: janussol.m 4154 2018-04-23 19:54:44Z bolle $

if isfield(PREC,'pivots')
   pivots=PREC.pivots;
else
   pivots=[];
end   
 
% Cholesky does not need pivoting 
if PREC.isdefinite && PREC.ishermitian & PREC.invert_blocks==0
   % dummy pivot vector
   if ~isfield(PREC,'pivots')
      pivots=-(PREC.n+1);
   end
end

if PREC.isreal
   % real symmetric case
   if PREC.issymmetric
      z=PREC.SL*(PREC.P*bildlsol(PREC.BL,PREC.BiD,pivots,PREC.P'*(PREC.SL*y)));
      % z=bildlsol(PREC.BL,PREC.BiD,y);
   % general real case
   else
      z=PREC.SR*(PREC.P*bilusol(PREC.BL,PREC.BiD,PREC.BUT,pivots,PREC.Q'*(PREC.SL*y)));
      % z=bilusol(PREC.BL,PREC.BiD,PREC.BUT,y);
   end
else % non-real case
   % non-real complex-symmetric case
   if PREC.issymmetric
      z=PREC.SL*(PREC.P*bildlsols(PREC.BL,PREC.BiD,pivots,PREC.P'*(PREC.SL*y)));
      % z=bildlsols(PREC.BL,PREC.BiD,y);
   % non-real hermitian case
   elseif PREC.hermitian
      z=PREC.SL*(PREC.P*bildlsol(PREC.BL,PREC.BiD,pivots,PREC.P'*(PREC.SL*y)));
      % z=bildlsol(PREC.BL,PREC.BiD,y);
   % general non-real case
   else
      z=PREC.SR*(PREC.P*bilusol(PREC.BL,PREC.BiD,PREC.BUT,pivots,PREC.Q'*(PREC.SL*y)));
      % z=bilusol(PREC.BL,PREC.BiD,PREC.BUT,y);
   end
end

