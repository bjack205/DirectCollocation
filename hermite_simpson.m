function [weights, funcs] = hermite_simpson(nSegments)

% Hermite Simpson (separated)
nGrid = 2*nSegments+1;
weights = ones(nGrid,1)*2/3;
weights(2:2:end) = 4/3;
weights([1 end]) = 1/3;

funcs.constraints = @constraints;
funcs.interp_ctrl = @pwPoly2;
funcs.interp_state = @pwPoly3;
end

function [c, ceq] = constraints(t, x, u, dynFcn)
    nGrid = size(x,2);
    n = size(x,1);
    dt = (t(end) - t(1))/(nGrid-1);
    
    iLow = 1:2:nGrid-1;
    iMid = iLow + 1;
    iUpp = iMid + 1;
    
    xLow = x(:,iLow);
    xMid = x(:, iMid);
    xUpp = x(:, iUpp);
    
    f = dynFcn(t, x, u);
    fLow = f(:, iLow);
    fMid = f(:, iMid);
    fUpp = f(:, iUpp);
    
    midpoints = xMid - (xUpp+xLow)/2 - dt*(fLow-fUpp)/4;
    collocation = xUpp - xLow - dt*(fLow + 4*fMid + fUpp)/3;
    
    defects = zeros(n,nGrid-1);
    defects(:,iLow) = collocation;
    defects(:,iMid) = midpoints;

    c = [];
    ceq = reshape(defects, n*(nGrid-1), 1);
end


function u = pwPoly2(t,x,tq)
    % Collocation times
    t_colloc = t(1:2:end);
    
    n = floor((length(t)-1)/2);
    m = size(x,1);
    k = length(tq);
    u = zeros(m,k);
    
    % Bin the query times
    edges = [-inf, t_colloc, inf];
    [~,bin] = histc(tq, edges);
    
    % Loop over each quadratic segment
    for i=1:n
        idx = bin==(i+1);
        if sum(idx) > 0
            gridIdx = 2*(i-1) + [1,2,3];
            u(:,idx) = quadInterp(t(gridIdx),x(:,gridIdx),tq(idx));
        end
    end

    % Replace any out-of-bounds queries with NaN
    outOfBounds = bin==1 | bin==(n+2);
    u(:,outOfBounds) = nan;

    % Check for any points that are exactly on the upper grid point:
    if sum(tq==t(end))>0
        u(:,tq==t(end)) = x(:,end);
    end
    
    
end


function x = quadInterp(tGrid,xGrid,t)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   tGrid = [1, 3] = time grid
%   xGrid = [m, 3] = function grid
%   t = [1, p] = query times, spanned by tGrid
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

% Rescale the query points to be on the domain [-1,1]
t = 2*(t-tGrid(1))/(tGrid(3)-tGrid(1)) - 1;

% Compute the coefficients:
a = 0.5*(xGrid(:,3) + xGrid(:,1)) - xGrid(:,2);
b = 0.5*(xGrid(:,3)-xGrid(:,1));
c = xGrid(:,2);

% Evaluate the polynomial for each dimension of the function:
p = length(t);
m = size(xGrid,1);
x = zeros(m,p);
tt = t.^2;
for i=1:m
    x(i,:) = a(i)*tt + b(i)*t + c(i);
end

end

function x = pwPoly3(tGrid,xGrid,fGrid,t)
% x = pwPoly3(tGrid,xGrid,fGrid,t)
%
% This function does piece-wise quadratic interpolation of a set of data,
% given the function value at the edges and midpoint of the interval of
% interest.
%
% INPUTS:
%   tGrid = [1, 2*n-1] = time grid, knot idx = 1:2:end
%   xGrid = [m, 2*n-1] = function at each grid point in time
%   fGrid = [m, 2*n-1] = derivative at each grid point in time
%   t = [1, k] = vector of query times (must be contained within tGrid)
%
% OUTPUTS:
%   x = [m, k] = function value at each query time
%
% NOTES:
%   If t is out of bounds, then all corresponding values for x are replaced
%   with NaN
%

nGrid = length(tGrid);
if mod(nGrid-1,2)~=0 || nGrid < 3
    error('The number of grid-points must be odd and at least 3');
end

% Figure out sizes
n = floor((length(tGrid)-1)/2);
m = size(xGrid,1);
k = length(t);
x = zeros(m, k);

% Figure out which segment each value of t should be on
edges = [-inf, tGrid(1:2:end), inf];
[~, bin] = histc(t,edges);

% Loop over each quadratic segment
for i=1:n
    idx = bin==(i+1);
    if sum(idx) > 0
        kLow = 2*(i-1) + 1;
        kMid = kLow + 1;
        kUpp = kLow + 2;
        h = tGrid(kUpp)-tGrid(kLow);
        xLow = xGrid(:,kLow);
        fLow = fGrid(:,kLow);
        fMid = fGrid(:,kMid);
        fUpp = fGrid(:,kUpp);
        alpha = t(idx) - tGrid(kLow);
        x(:,idx) = cubicInterp(h,xLow, fLow, fMid, fUpp,alpha);
    end
end

% Replace any out-of-bounds queries with NaN
outOfBounds = bin==1 | bin==(n+2);
x(:,outOfBounds) = nan;

% Check for any points that are exactly on the upper grid point:
pnts = sum(t==tGrid(end));
if pnts > 0
    if pnts > 1
        a = 1;
    end
    x(:,t==tGrid(end)) = repmat(xGrid(:,end),1,pnts);
end

end


function x = cubicInterp(h,xLow, fLow, fMid, fUpp,del)
%
% This function computes the interpolant over a single interval
%
% INPUTS:
%   h = time step (tUpp-tLow)
%   xLow = function value at tLow
%   fLow = derivative at tLow
%   fMid = derivative at tMid
%   fUpp = derivative at tUpp
%   del = query points on domain [0, h]
%
% OUTPUTS:
%   x = [m, p] = function at query times
%

%%% Fix matrix dimensions for vectorized calculations
nx = length(xLow);
nt = length(del);
xLow = xLow*ones(1,nt);
fLow = fLow*ones(1,nt);
fMid = fMid*ones(1,nt);
fUpp = fUpp*ones(1,nt);
del = ones(nx,1)*del;

a = (2.*(fLow - 2.*fMid + fUpp))./(3.*h.^2);
b = -(3.*fLow - 4.*fMid + fUpp)./(2.*h);
c = fLow;
d = xLow;

x = d + del.*(c + del.*(b + del.*a));

end
