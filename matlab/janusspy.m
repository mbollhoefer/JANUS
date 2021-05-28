function janusspy(A,PREC)
% janusspy(PREC)
%
% display multilevel preconditioner PREC
%
% janusspy(A,PREC)
%
% display remapped original matrix A associated with the sequence of 
% reorderings given by PREC

% $Id: janusspy.m 4006 2018-02-14 20:38:39Z bolle $


if nargin==1
   PREC=A;

   n=PREC.n;
   nB=length(PREC.BiD);
   
   L=sparse(n,n);
   D=sparse(n,n);
   if isfield(PREC,'BUT')
      UT=sparse(n,n);
   end

   nlev=length(PREC);
   sumnB=0;
   for i=1:nB
       
       % DD=PREC.BiD{i}.D;
       J=PREC.BiD{i}.J;
       nJ=length(J);
       % LL=PREC.BL{i}.L;
       IL=PREC.BL{i}.I;
       nIL=length(IL);
       % UUT=PREC.BUT{i}.L;

       if isfield(PREC,'BUT')
	  IUT=PREC.BUT{i}.I;
	  nIUT=length(IUT);
	  % UT(IUT,J)=ones(nIUT,nJ);
	  UT(IUT,J)=1;
       else
	  IUT=IL;
	  nIUT=nIL;
       end
       
       % D(J,J)=ones(nJ,nJ);
       % L(IL,J)=ones(nIL,nJ);
       
       D(J,J)=1;
       L(IL,J)=1;
       
   end % for

   spy(L,'g'); hold on
   if isfield(PREC,'BUT')
      spy(UT','b');
   else
      spy(L','b');
   end
   spy(D,'r'); 
 
    
   title(['JANUS block ILU (' num2str(nB) ' blocks)'])
   xlabel(['nz=' num2str(janusnnz(PREC))])
   hold off;
    
    
    
else % two input arguments, display A
  
   if isfield(PREC,'Q')
      spy(PREC.Q'*A*PREC.P);
   else
      spy(PREC.P'*A*PREC.P);
   end
end % if-else nargin==1
