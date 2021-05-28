function nz=janusnnz(PREC)
% nz=bilunnz(PREC)
%
% number of nonzeros of the block ILU
%
% Input
% -----
% PREC  preconditioner structure
% PREC.BL   lower block triangular matrix in cell array format, BL~L
% PREC.BiD  block diagonal matrix in cell array format BiD~D^{-1}
% PREC.BUT  optionally lower block triangular matrix in cell array format, BUT^T~U
%
% Output
% ------
% nz      number of nonzeros


%    Authors:
%
%	Matthias Bollhoefer, TU Braunschweig
%
%    Date:
%
%	August 29, 2017. JANUS Block ILU R1.0.  
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
% $Id: janusnnz.m 3621 2017-09-09 15:20:49Z bolle $

if ~isfield(PREC,'BUT')
   nz=bilunnz(PREC.BL,PREC.BiD,PREC.BL);
else
   nz=bilunnz(PREC.BL,PREC.BiD,PREC.BUT);
end


