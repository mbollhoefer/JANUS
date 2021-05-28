function [nz,mxblock,avgblock,stddev]=blockstatistics(BL,BD,BUT)
% [nz,mxblock,avgblock,stddev]=blockstatistics(PREC)
% [nz,mxblock,avgblock,stddev]=blockstatistics(BL)
% [nz,mxblock,avgblock,stddev]=blockstatistics(BL,BD)
% [nz,mxblock,avgblock,stddev]=blockstatistics(BL,BD,BUT)
%
% statistics about fill and diagonal block sizes
%  
% Input
% -----  
% PREC       structured JANUS block ILU PREC
%  
% BL         structure with block lower triangular matrix
% BD         structure with block diagonal matrix
% BUT        structure with transposed block upper triangular matrix
%
% Output
% ------
% nz         total number of nonzeros
% mxblock    maximum diagonal block size
% avgblock   average diagonal block size
% stddev     standard deviation w.r.t. the block diagonal sizes
  
  

% $Id: blockstatistics.m 6123 2020-04-01 12:26:30Z bolle $

if nargin==1
   if isfield(BL,'BUT')
      BUT=BL.BUT;
   end
   if isfield(BL,'BiD')
      BD =BL.BiD;
   end
   if isfield(BL,'BL')
      BL =BL.BL;
   end
end



nblocks=length(BL);
n=max(BL{nblocks}.J);

d=zeros(nblocks,1);
for k=1:nblocks 
    d(k)=length(BL{k}.J); 
end % for k 
if exist('BUT')
   du=zeros(nblocks,1);
   for k=1:nblocks 
       du(k)=length(BUT{k}.J); 
   end % for k 
end

if ~exist('BUT')
   if exist('BD')
      nz=bilunnz(BL,BD,BL);
   else
      nz=0;
      for k=1:nblocks 
	  nz=nz+nnz(BL{k}.D)+nnz(BL{k}.L);
      end % for k 
   end
else
   nz=bilunnz(BL,BD,BUT);
end



mb=max(d);
b=zeros(1,mb);
for k=1:nblocks 
    % blocksize
    bs=d(k);
    % count number of blocks
    b(bs)=b(bs)+1;
end % for lk
if exist('BUT')
   mbu=max(du);
   bu=zeros(1,mbu);
   for k=1:nblocks 
       % blocksize
       bs=du(k);
       % count number of blocks
       bu(bs)=bu(bs)+1;
   end % for lk
end

% mean value
mean=sum(d)/length(d);
if exist('BUT')
   meanu=sum(du)/length(du);
end

% standard deviation
stddev=sqrt(sum((d-mean).^2)/(nblocks-1));
if exist('BUT')
   stddevu=sqrt(sum((du-meanu).^2)/(nblocks-1));
end

% median
d=sort(d);
l=length(d)/2;
if l==ceil(l)
   median=d(l);
else
   median=(d(floor(l))+d(ceil(l)))/2;
end   
if exist('BUT')
   du=sort(du);
   l=length(du)/2;
   if l==ceil(l)
      medianu=du(l);
   else
      medianu=(du(floor(l))+du(ceil(l)))/2;
   end   
end   


figure;
% plot(d,'-b')
% title('block size sorted in increasing order','FontSize',16,'FontWeight','bold')
% 
% figure(2);
% bar(100*b./sum(b))
bar(100*b.*(1:length(b))./n)
title('relative block distribution L','FontSize',16,'FontWeight','bold')

str=sprintf('block size L, mean: %5.1f, std. dev. %5.1f, median: %5.1f, max: %3d',mean,stddev,median,mb);
mxblock=mb;
avgblock=mean;
fprintf('%s\n',str);
xlabel(str,'FontSize',16,'FontWeight','bold');
ylabel('percentage of blocks','FontSize',16,'FontWeight','bold')
if exist('BUT')
   figure;
   % bar(100*bu./sum(bu))
   bar(100* bu.*(1:length(bu))./n)
   title('relative block distribution U','FontSize',16,'FontWeight','bold')

   str=sprintf('block size U, mean: %5.1f, std. dev. %5.1f, median: %5.1f, max: %3d',meanu,stddevu,medianu,mbu);
   mxblock =round((mxblock+mbu)/2);
   avgblock=(avgblock+meanu)/2;
   stddev  =(stddev+stddevu)/2;
   fprintf('%s\n',str);
   xlabel(str,'FontSize',16,'FontWeight','bold');
   ylabel('percentage of blocks','FontSize',16,'FontWeight','bold')
end
 
