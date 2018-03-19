function imagesc_nan(varargin)


if nargin == 1
    data = varargin{1};
    data(isinf(data)) = NaN;
    [nr,nc] = size(data);
    pcolor([data nan(nr,1); nan(1,nc+1)]);
elseif nargin == 3
    x0= varargin{1};
    y0 = varargin{2};
    x = linspace(min(x0),max(x0),length(y0));
    y = linspace(min(y0),max(y0),length(x0));
    data = varargin{3};
    data(isinf(data)) = NaN;
    [nr,nc] = size(data);
    pcolor([x(:)' x(end)],[y(:)' y(end)],[data nan(nr,1); nan(1,nc+1)]);
end
shading flat;

end
