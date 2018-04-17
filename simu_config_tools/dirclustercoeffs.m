function [cc cccyc ccmid ccin ccout nf]=dirclustercoeffs(A,weighted,normalized)
% DIRCLUSTERCOEFFS Compute clustering coefficients for a directed graph
%
% cc=dirclustercoeffs(A) returns the directed clustering coefficients
% (which generalize the clustering coefficients of an undirected graph, 
% and so calling this function on an undirected graph will produce the same
% answer as clustercoeffs, but less efficiently.)
%
% This function implements the algorithm from Fagiolo, Phys Rev. E. 76
% 026107 (doi:10:1103/PhysRevE.76.026107).  
%
% [cc,cccyc,ccmid,ccin,ccout,nf]=dirclusteringcoeffs(A) returns different 
% components of the clustering coefficients corresponding to cycles,
% middles, in triangles, and out triangles.  See the manuscript for a 
% description of the various types of triangles counted in the above
% metrics.
%
% See also CLUSTERCOEFFS
%
% Example:
%   load_gaimc_graph('celegans'); % load the C elegans nervous system network
%   cc=dirclustercoeffs(A);
%   [maxval maxind]=max(cc)
%   labels(maxind) % most clustered vertex in the nervous system

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-22: Initial coding
% 2009-05-15: Documentation and example

if ~exist('normalized','var') || isempty(normalized), normalized=true; end
if ~exist('weighted','var') || isempty(weighted), weighted=true; end
donorm=1; usew=1;
if ~normalized, donorm=0; end
if ~weighted, usew=0; end

if isstruct(A)
    rp=A.rp; ci=A.ci; %ofs=A.offset;
    cp=A.cp; ri=A.ri; % get 
    if usew, ai=A.ai; ati=A.ati; end
else
    if usew, [rp ci ai]=sparse_to_csr(A); [cp ri ati]=sparse_to_csr(A'); 
    else [rp ci]=sparse_to_csr(A);  [cp ri]=sparse_to_csr(A'); 
    end
    if any(ai)<0, error('gaimc:clustercoeffs',...
            ['only positive edge weights allowed\n' ...
             'try dirclustercoeffs(A,0) for an unweighted comptuation']); 
    end
end
n=length(rp)-1;

% initialize all the variables
cc=zeros(n,1); ind=false(n,1); cache=zeros(n,1); degs=zeros(n,1);
if nargout>1, cccyc=zeros(n,1); end
if nargout>2, ccmid=zeros(n,1); end
if nargout>3, ccin=zeros(n,1); end
if nargout>4, ccout=zeros(n,1); end
if nargout>5, nf=zeros(n,1); end
% precompute degrees
for v=1:n, 
    for rpi=rp(v):rp(v+1)-1, 
        w=ci(rpi); 
        if v==w, continue; else degs(w)=degs(w)+1; degs(v)=degs(v)+1; end
    end
end
ew=1; ew2=1;
for v=1:n
    % setup counts for the different cycle types
    bilatedges=0; curcccyc=0; curccmid=0; curccin=0; curccout=0;
    % 1.  
    % find triangles with out links as last step, so precompute the inlinks
    % back to node v
    for cpi=cp(v):cp(v+1)-1
        w=ri(cpi); if usew, ew=ati(cpi); end
        if v~=w, ind(w)=1; cache(w)=ew^(1/3); end
    end
    % count cycles (cycles are out->out->out)
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); if v==w, continue; end % discount self-loop
        if usew, ew=ai(rpi)^(1/3); end
        for rpi2=rp(w):rp(w+1)-1
            x=ci(rpi2); if x==w, continue; end
            if x==v, bilatedges=bilatedges+1; continue; end
            if ind(x)
                if usew, ew2=ai(rpi2); end
                curcccyc=curcccyc+ew*ew2^(1/3)*cache(x);
            end
        end
    end
    % count middle-man circuits (out->in->out)
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); if v==w, continue; end % discount self-loop
        if usew, ew=ai(rpi)^(1/3); end
        for cpi=cp(w):cp(w+1)-1
            x=ri(cpi); if x==w, continue; end
            if ind(x)
                if usew, ew2=ati(cpi); end
                curccmid=curccmid+ew*ew2^(1/3)*cache(x);
            end
        end
    end
    % count in-link circuits (in->out->out)
    for cpi=cp(v):cp(v+1)-1
        w=ri(cpi); if v==w, continue; end % discount self-loop
        if usew, ew=ati(cpi)^(1/3); end
        for rpi2=rp(w):rp(w+1)-1
            x=ci(rpi2); if x==w, continue; end
            if ind(x)
                if usew, ew2=ai(rpi2); end
                curccin=curccin+ew*ew2^(1/3)*cache(x);
            end
        end
    end
    % reset and reinit the cache for outlinks
    for cpi=cp(v):cp(v+1)-1, w=ri(cpi); ind(w)=0; end % reset indicator
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); if usew, ew=ai(rpi); end
        if v~=w, ind(w)=1; cache(w)=ew^(1/3); end
    end
    % count out-link circuits (out->out->in)
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); if v==w, continue; end % discount self-loop
        if usew, ew=ai(rpi)^(1/3); end
        for rpi2=rp(w):rp(w+1)-1
            x=ci(rpi2); if x==w, continue; end
            if ind(x)
                if usew, ew2=ai(rpi2); end
                curccout=curccout+ew*ew2^(1/3)*cache(x);
            end
        end
    end
    for rpi=rp(v):rp(v+1)-1, w=ci(rpi); ind(w)=0; end % reset indicator
    % store the values
    nf=degs(v)*(degs(v)-1) - 2*bilatedges;
    curcc=curcccyc+curccmid+curccin+curccout;
    if nf>0 && donorm, curcc=curcc/nf; end
    cc(v)=curcc;
    if nargout>1, cccyc(v)=curcccyc; end
    if nargout>2, ccmid(v)=curccmid; end
    if nargout>3, ccin(v)=curccin; end
    if nargout>4, ccout(v)=curccout; end
    if nargout>5, nf(v)=nf; end
end
end



%% Helper function
function [rp ci ai ncol]=sparse_to_csr(A,varargin)
% SPARSE_TO_CSR Convert a sparse matrix into compressed row storage arrays
% 
% [rp ci ai] = sparse_to_csr(A) returns the row pointer (rp), column index
% (ci) and value index (ai) arrays of a compressed sparse representation of
% the matrix A.
%
% [rp ci ai] = sparse_to_csr(i,j,v,n) returns a csr representation of the
% index sets i,j,v with n rows.
%
% [rp ci ai] = sparse_to_csr(struct('rp',rp,'ci',ci)); returns the same
% precomputed answers 
%
% Example:
%   A=sparse(6,6); A(1,1)=5; A(1,5)=2; A(2,3)=-1; A(4,1)=1; A(5,6)=1; 
%   [rp ci ai]=sparse_to_csr(A)
%
% See also CSR_TO_SPARSE, SPARSE

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-07: Initial version
% 2008-04-24: Added triple array input
% 2009-05-01: Added ncol output
% 2009-05-15: Fixed triplet input
% 2010-01-29: Added precomputed struct input

error(nargchk(1, 5, nargin, 'struct'))
retc = nargout>1; reta = nargout>2;

if nargin>1
    if nargin>4, ncol = varargin{4}; end
    nzi = A; nzj = varargin{1};
    if reta && length(varargin) > 2, nzv = varargin{2}; end    
    if nargin<4, n=max(nzi); else n=varargin{3}; end
    nz = length(A);
    if length(nzi) ~= length(nzj), error('gaimc:invalidInput',...
            'length of nzi (%i) not equal to length of nzj (%i)', nz, ...
            length(nzj)); 
    end
    if reta && length(varargin) < 3, error('gaimc:invalidInput',...
            'no value array passed for triplet input, see usage'); 
    end
    if ~isscalar(n), error('gaimc:invalidInput',...
            ['the 4th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
    if nargin < 5, ncol = max(nzj); 
    elseif ~isscalar(ncol), error('gaimc:invalidInput',...
            ['the 5th input to sparse_to_csr with triple input was not ' ...
             'a scalar']); 
    end
else
    if isstruct(A), % handle the precomputed case
        rp=A.rp; ci=A.ci;
        if nargout>2, ai=A.ai; end
        if nargout>3, ncol=A.ncol; end
        return;
    end
    n = size(A,1); nz = nnz(A); ncol = size(A,2);
    retc = nargout>1; reta = nargout>2;
    if reta,     [nzi nzj nzv] = find(A); 
    else         [nzi nzj] = find(A);
    end
end
if retc, ci = zeros(nz,1); end
if reta, ai = zeros(nz,1); end
rp = zeros(n+1,1);
for i=1:nz
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp);
if ~retc && ~reta, rp=rp+1; return; end
for i=1:nz
    if reta, ai(rp(nzi(i))+1)=nzv(i); end
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;
end

