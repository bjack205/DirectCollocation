function soln = direct_collocation(costFcn, dynFcn, b, nSegments, method, guess, verbose)
if nargin < 7
    verbose = true;
end


%% Setup
% Useful constants
n = length(b.xf.u);
m = length(b.u.u);
tf = b.tf.l;

% Collocation Method
if strcmp(method, 'hermite-simpson')
    [weights, funcs] = hermite_simpson(nSegments);
elseif strcmp(method, 'trapezoid')
    [weights, funcs] = trapezoid(nSegments);
end
nGrid = length(weights);

% Initial Guess
ti = linspace(guess.time(1), guess.time(end), nGrid);
xi = interp1(guess.time, guess.state', ti)';
ui = interp1(guess.time, guess.control, ti);
[zi, pk] = pack(ti([1,end]), xi, ui);

% BCs
tLow = [0; b.tf.l];
tUpp = [0; b.tf.u];
xLow = [b.x0.l, ones(n, nGrid-2).*b.x.l, b.xf.l];
xUpp = [b.x0.u, ones(n, nGrid-2).*b.x.u, b.xf.u];
uLow = ones(m, nGrid).*b.u.l;
uUpp = ones(m, nGrid).*b.u.u;

zLow = pack(tLow, xLow, uLow);
zHigh = pack(tUpp, xUpp, uUpp);

%% Pass to NLP
% Formulate Problem
objFcn = @(z) objective(z, pk, costFcn, weights);
conFcn = @(z) nlcon(z, pk, funcs.constraints, dynFcn);

% Solve the Problem
if verbose
    display = 'iter';
else
    display = 'none';
end
nlpOpt = optimset(...
    'Display',display,...
    'MaxFunEvals',1e5);

tic;
[zSoln, objVal,exitFlag,output] = fmincon(objFcn,zi,[],[],[],[],zLow,zHigh,conFcn,nlpOpt);
t_soln = toc;
[tS, xS, uS] = unpack(zSoln, pk);

%% Post-processing
% Store the results
soln.grid.time = tS;
soln.grid.state = xS;
soln.grid.control = uS;

% Iterpolation consistent with the method
fS = dynFcn(tS,xS,uS);
soln.interp.state = @(t)(funcs.interp_state(tS,xS,fS,t) );
soln.interp.control = @(t)(funcs.interp_ctrl(tS,uS,t));
soln.interp.state_derivative = @(t)(funcs.interp_ctrl(tS,fS,t));

% Extra info
soln.info = output;
soln.info.exitFlag = exitFlag;
soln.info.objVal = objVal;
soln.info.t_computation = t_soln;

% Error metrics
soln.interp.collCst = @(t)( ...
    dynFcn(t, soln.interp.state(t), soln.interp.control(t))...
    - soln.interp.state_derivative(t) );

% Use multi-segment simpson quadrature to estimate the absolute local error
% along the trajectory.
absColErr = @(t)(abs(soln.interp.collCst(t)));
nSegment = nGrid-1;
nState = size(xS,1);
quadTol = 1e-12;   %Compute quadrature to this tolerance  
soln.info.error = zeros(nState,nSegment);
for i=1:nSegment
    soln.info.error(:,i) = rombergQuadrature(absColErr,tS([i,i+1]),quadTol);
    if verbose
        fprintf('Finished Error on Segment %i of %i\n', i, nSegment)
    end
end
soln.info.maxError = max(max(soln.info.error));


end

function [z, pk] = pack(t,x,u)
    N = numel(t) + numel(x) + numel(u);
    nGrid = size(x,2);
    n = size(x,1);
    m = size(u,1);
    
    indz = reshape(2+(1:N-2), n+m, nGrid);
    indx = indz(1:n, :);
    indu = indz(n+(1:m), :);
    
    z = zeros(N,1);
    
    z(1:2) = t;
    z(indx(:)) = x(:);
    z(indu(:)) = u(:);
    
    pk.indz = indz;
    pk.n = n;
end

function [t,x,u] = unpack(z, pk)
N = length(z);
n = pk.n;
nGrid = size(pk.indz,2);
m = size(pk.indz,1)-n;

t = linspace(z(1), z(2), nGrid);
x = z(pk.indz(1:n, :));
u = z(pk.indz(n+(1:m), :));
x = reshape(x, n, nGrid);
u = reshape(u, m, nGrid);
end

function cost = objective(z, pk, costFcn, weights)
    [t,x,u] = unpack(z, pk);
    nGrid = size(x,2);
    dt = (t(2) - t(1))/(nGrid-1);
    
    cost = costFcn(t,x,u);
    cost = dt*cost*weights;
end

function [c, ceq] = nlcon(z, pk, confun, dynfun)
[t,x,u] = unpack(z, pk);
[c, ceq] = confun(t,x,u,dynfun);
end
