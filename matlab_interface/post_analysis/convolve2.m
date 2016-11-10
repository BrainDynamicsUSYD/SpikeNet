function y = convolve2(x, m, shape, tol)
%CONVOLVE2 Two dimensional convolution.
%   Y = CONVOLVE2(X, M) performs the 2-D convolution of matrices X and
%   M. If [mx,nx] = size(X) and [mm,nm] = size(M), then size(Y) =
%   [mx+mm-1,nx+nm-1]. Values near the boundaries of the output array are
%   calculated as if X was surrounded by a border of zero values.
%
%   Y = CONVOLVE2(X, M, SHAPE) where SHAPE is a string returns a
%   subsection of the 2-D convolution with size specified by SHAPE:
%
%       'full'      - (default) returns the full 2-D convolution
% 
%       'valid'     - returns only those parts of the convolution
%                     that can be computed without padding; size(Y) =
%                     [mx-mm+1,nx-nm+1] when size(X) > size(M)
% 
%       'same'      - returns the central part of the convolution
%                     that is the same size as X using zero padding
% 
%       'wrap' or
%       'circular'  - as for 'same' except that instead of using
%                     zero-padding the input X is taken to wrap round as
%                     on a toroid
% 
%       'reflect' or
%       'symmetric' - as for 'same' except that instead of using
%                     zero-padding the input X is taken to be reflected at
%                     its boundaries
% 
%       'replicate' - as for 'same' except that instead of using
%                     zero-padding the rows at the array boundary are
%                     replicated
%
%   CONVOLVE2 is fastest when mx > mm and nx > nm - i.e. the first
%   argument is the input and the second is the mask.
%
%   If the rank of the mask M is low, CONVOLVE2 will decompose it into a
%   sum of outer product masks, each of which is applied efficiently as
%   convolution with a row vector and a column vector, by calling CONV2.
%   The function will often be faster than CONV2 or FILTER2 (in some
%   cases much faster) and will produce the same results as CONV2 to
%   within a small tolerance.
%
%   Y = CONVOLVE2(... , TOL) where TOL is a number in the range 0.0 to
%   1.0 computes the convolution using a reduced-rank approximation to
%   M, provided this will speed up the computation. TOL limits the
%   relative sum-squared error in the effective mask; that is, if the
%   effective mask is E, the error is controlled such that
%
%       sum(sum( (M-E) .* (M-E) ))
%       --------------------------    <=  TOL
%            sum(sum( M .* M ))
%
%   See also CONV2, FILTER2, EXINDEX

% Copyright David Young, Feb 2002, revised Jan 2005, Jan 2009, Apr 2011,
% Feb 2014

% Deal with optional arguments
narginchk(2,4);
if nargin < 3
    shape = 'full';    % shape default as for CONV2
    tol = 0;
elseif nargin < 4
    if isnumeric(shape)
        tol = shape;
        shape = 'full';
    else
        tol = 0;
    end
end;

% Set up to do the wrap & reflect operations, not handled by conv2
if ismember(shape, {'wrap' 'circular' 'reflect' 'symmetric' 'replicate'})
    x = extendarr(x, m, shape);
    shape = 'valid';
end

% do the convolution itself
y = doconv(x, m, shape, tol);
end

%-----------------------------------------------------------------------

function y = doconv(x, m, shape, tol)
% Carry out convolution
[mx, nx] = size(x);
[mm, nm] = size(m);

% If the mask is bigger than the input, or it is 1-D already,
% just let CONV2 handle it.
if mm > mx || nm > nx || mm == 1 || nm == 1
    y = conv2(x, m, shape);
else
    % Get svd of mask
    if mm < nm; m = m'; end        % svd(..,0) wants m > n
    [u,s,v] = svd(m, 0);
    s = diag(s);
    rank = trank(m, s, tol);
    if rank*(mm+nm) < mm*nm         % take advantage of low rank
        if mm < nm;  t = u; u = v; v = t; end  % reverse earlier transpose
        vp = v';
        % For some reason, CONV2(H,C,X) is very slow, so use the normal call
        y = conv2(conv2(x, u(:,1)*s(1), shape), vp(1,:), shape);
        for r = 2:rank
            y = y + conv2(conv2(x, u(:,r)*s(r), shape), vp(r,:), shape);
        end
    else
        if mm < nm; m = m'; end     % reverse earlier transpose
        y = conv2(x, m, shape);
    end
end
end

%-----------------------------------------------------------------------

function r = trank(m, s, tol)
% Approximate rank function - returns rank of matrix that fits given
% matrix to within given relative rms error. Expects original matrix
% and vector of singular values.
if tol < 0 || tol > 1
    error('Tolerance must be in range 0 to 1');
end
if tol == 0             % return estimate of actual rank
    tol = length(m) * max(s) * eps;
    r = sum(s > tol);
else
    ss = s .* s;
    t = (1 - tol) * sum(ss);
    r = 0;
    sm = 0;
    while sm < t
        r = r + 1;
        sm = sm + ss(r);
    end
end
end

%-----------------------------------------------------------------------

function y = extendarr(x, m, shape)
% Extend x so as to wrap around on both axes, sufficient to allow a
% "valid" convolution with m to return a result the same size as x.
% We assume mask origin near centre of mask for compatibility with
% "same" option.

[mx, nx] = size(x);
[mm, nm] = size(m);

mo = floor((1+mm)/2); no = floor((1+nm)/2);  % reflected mask origin
ml = mo-1;            nl = no-1;             % mask left/above origin
mr = mm-mo;           nr = nm-no;            % mask right/below origin

% deal with shape option terminology - was inconsistent with exindex
switch shape
    case 'wrap'
        shape = 'circular';
    case 'reflect'
        shape = 'symmetric';
end
y = exindex(x, 1-ml:mx+mr, 1-nl:nx+nr, shape);

end

function arr = exindex(arr, varargin)
%EXINDEX extended array indexing
%   ARROUT = EXINDEX(ARRIN, S1, S2, ...) indexes a virtual array made by
%   extending ARRIN with zeros in all directions, using subscripts S1, S2
%   etc.
%
%   ARROUT = EXINDEX(ARRIN, S1, R1, S2, R2, ...) extends ARRIN using rule
%   R1 on the first dimension, R2 on the second dimension etc.
%
%   ARROUT = EXINDEX(ARRIN, S1, S2, ..., R) extends ARRIN using rule R on
%   every dimension.
%
%   Subscripts
%   ----------
%
%   Broadly, if V is the virtual extended array, ARROUT = V(S1, S2, ...)
%
%   The elements of the subscript arguments S1, S2 etc must be integers.
%   They need not be positive and are not restricted in any way by the size
%   of ARRIN. Logical indexing and linear indexing are not supported.
%
%   There must be at least one subscript argument for each dimension of
%   ARRIN as reported by NDIMS, except that row and column vectors may have
%   1 or 2 subscripts. A single subscript is taken to refer to the
%   dimension along which the vector lies, as in normal vector indexing.
%   Scalars require 2 subscripts. If there are more subscripts than
%   dimensions, ARRIN is taken to have trailing singleton dimensions, as in
%   normal array indexing.
%
%   The number of dimensions of ARROUT will be the number of subscript
%   arguments, though trailing singleton dimensions will, as usual, be
%   suppressed. The size of ARROUT is given by the normal Matlab rules for
%   the result of indexing into ARRIN: that is
%
%       size(ARROUT) = size( ARRIN(ones(size(S1)), ones(size(S2)), ...) )
%
%   A subscript argument may be the string ':'. This behaves like a colon
%   in ordinary subscripting: a colon for the K'th subscript stands for
%   1:size(ARRIN, K). The 'end' keyword is not supported.
%
%   Rules
%   -----
%
%   Each rule may be one of the following:
%
%   A scalar cell: ARRIN is padded with elements equal to the contents of
%   the cell. The class of the cell contents must be compatible with the
%   class of ARRIN.
%
%       If different constants are used on different dimensions, padding is
%       done in the order of the subscripts. For example, a 2D array is
%       extended first in the row index direction and then in the column
%       index direction. For all other cases, the order in which dimensions
%       are extended has no effect.
%
%   'circular': ARRIN is extended with copies of itself; i.e. V is tiled
%   with ARRIN.
%
%   'symmetric': ARRIN is extended with copies of itself with reflection at
%   its boundaries; i.e. V is tiled with [ARRIN fliplr(ARRIN);
%   flipud(ARRIN) fliplr(flipud(ARRIN))].
%
%   'replicate': ARRIN is extended by copying its border elements; i.e. an
%   element of V is equal to the nearest element of ARRIN.
%
%   If no rule is given, padding is with zeros.
%
%   Examples
%   --------
%
%   Pad a 2D matrix with K extra rows and columns with reflection on both
%   axes:
%
%       b = exindex(a, 1-k:size(a,1)+k, 1-k:size(a,2)+k, 'symmetric');
%
%   Circularly shift a 2D matrix by R rows downwards and C columns
%   rightwards:
%
%       b = exindex(a, 1-r:size(a,1)-r, 1-c:size(a,2)-c, 'circular');
%
%   Force a row or column vector to be 1024 elements long, trimming or
%   padding with zeros as necessary:
%
%       u = exindex(v, 1:1024);
%
%   The same, with a non-zero padding value:
%
%       u = exindex(v, 1:1024, {-1});   % note constant in cell
%
%   Truncate or extend all the rows of a matrix to 1024 columns:
%
%       b = exindex(a, ':', 1:1024);
%
%   Extend a 2-D array into the third dimension by copying it:
%
%       b = exindex(a, ':', ':', 1:3, 'replicate');
%
%   Pad a 1-D cell array with cells containing the empty matrix:
%
%       cellout = exindex(cellin, 0:10, {{[]}}); 
%
%   See also: padarray, circshift, repmat

% Copyright David Young 2010

% Sort out arguments
[exindices, rules, nd, sz] = getinputs(arr, varargin{:});
consts = cellfun(@iscell, rules);  % Check for constants, as can be
constused = any(consts);           % more efficient if there are none

% Setup for constant padding
if constused
    tofill = cell(1, nd);
end

% Main loop over subscript arguments, transforming them into valid
% subscripts into arr using the rule for each dimension
if constused
    for i = 1:nd
        [exindices{i}, tofill{i}] = extend(exindices{i}, rules{i}, sz(i));
    end
else % no need for information for doing constants
    for i = 1:nd
        exindices{i} = extend(exindices{i}, rules{i}, sz(i));
    end
end

% Create the new array by indexing into arr. If there are no constants,
% this does the whole job
arr = arr(exindices{:});

% Fill areas that need constants
if constused
    % Get full range of output array indices
    ranges = arrayfun(@(x) {1:x}, size(arr));
    for i = nd:-1:1    % order matters
        if consts(i)
            ranges{i} = tofill{i};      % don't overwrite original
            c = rules{i};               % get constant and fill ...
            arr(ranges{:}) = c{1};      % we've checked c is scalar
            ranges{i} = ~tofill{i};     % don't overwrite
        end
    end
end

end

% -------------------------------------------------------------------------

function [exindices, rules, nd, sz] = getinputs(arr, varargin)
% Sort out and check arguments. Inputs are as given in the help comments
% for exindex. Outputs are cell arrays; each element of exindices is a
% set of integer extended indices which has been checked for validity; each
% element of rules is a rule which has not been checked for validity.

% Use index/rules arguments only to establish no. dimensions - ndims(arr)
% is no use, as trailing singleton dimensions truncated and vectors can be
% 2D or 1D
nd = length(varargin);
if nd == 0
    error('exindex:missingargs', 'Not enough arguments');
elseif nd == 1
    exindices = varargin;
    rules = {{0}};
elseif ~(isnumeric(varargin{2}) || strcmp(varargin{2}, ':'))
    % have alternating indices and rule
    nd = nd/2;
    if round(nd) ~= nd
        error('exindex:badnumargs', ...
            'Odd number of arguments after initial index/rule pair');
    end
    exindices = varargin(1:2:end);
    rules = varargin(2:2:end);
elseif nd > 2 && ~(isnumeric(varargin{end}) || strcmp(varargin{end}, ':'))
    % have a general rule at end
    nd = nd - 1;
    exindices = varargin(1:nd);
    [rules{1:nd}] = deal(varargin{end});
else
    % no rule is specified
    exindices = varargin;
    [rules{1:nd}] = deal({0});
end

% Sort out mismatch of apparent array size and number of dimensions
% indexed
sz = size(arr);
ndarr = ndims(arr);
if nd < ndarr
    if nd == 1 && ndarr == 2
        % Matlab allows vectors to be indexed with a single subscript and
        % to retain their shape. In all other cases (including scalars) a
        % single subscript causes the output to take the same shape as the
        % subscript array - we can't deal with this.
        if sz(1) == 1 && sz(2) > 1
            % have a row vector
            exindices = [{1} exindices {1}];
            rules = [rules rules];  % 1st rule doesn't matter
        elseif sz(2) == 1 && sz(1) > 1
            % have a column vector
            exindices = [exindices {1}];
            rules = [rules rules];  % 2nd rule doesn't matter
        else
            error('exindex:wantvector', ...
                'Only one index but array is not a vector');
        end
    else
        error('exindex:toofewindices', ...
            'Array has more dimensions than there are index arguments');
    end
    nd = 2;
elseif nd > ndarr
    % Effective array size
    sz = [sz ones(1, nd-ndarr)];
end

% Expand any colons now to simplify checking.
% It's tempting to allow the 'end' keyword here: easy to substitute the
% size of the dimension. However, to be worthwhile it would be necessary to
% use evalin('caller',...) so that expressions using end could be given as
% in normal indexing. This would mean moving the code up to exindex itself,
% and evalin makes for inefficiency and fragility, so this hasn't been
% done.
colons = strcmp(exindices, ':');
if any(colons)  % saves a little time
    exindices(colons) = arrayfun(@(x) {1:x}, sz(colons));
end

% Check the indices (rules are checked as required in extend)
checkindex = @(ind) validateattributes(ind, {'numeric'}, ...
    {'integer'}, 'exindex', 'index');
cellfun(checkindex, exindices);

end

% -------------------------------------------------------------------------

function [ind, tofill] = extend(ind, rule, s)
% The core function: maps extended array subscripts into valid input array
% subscripts.

if ischar(rule)    % pad with rule
    
    tofill = [];  % never used
    switch rule
        case 'replicate'
            ind = min( max(1,ind), s );
        case 'circular'
            ind = mod(ind-1, s) + 1;
        case 'symmetric'
            ind = mod(ind-1, 2*s) + 1;
            ott = ind > s;
            ind(ott) = 2*s + 1 - ind(ott);
        otherwise
            error('exindex:badopt', 'Unknown option');
    end
    
elseif iscell(rule) && isscalar(rule)     % pad with constant
    
    % The main messiness is due to constant padding. This can't be done
    % with indexing into the original array, but we want the indexing
    % structure to be preserved, so for now we index to element 1 on each
    % dimension, and record the indices of the regions that need to be
    % fixed.
    
    tofill = ind < 1 | ind > s;
    ind(tofill) = 1;
    
else
    
    error('exindex:badconst', 'Expecting string or scalar cell');
    
end

end

