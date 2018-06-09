addpath('..\OptimTraj\demo\cartPole')
addpath('..\OptimTraj')

%% Problem Declaration
% Declare Constants
p.M = 2;  % kg
p.m = 0.5;  % kg
p.ell = 0.5;  % m
p.I = 0.006;  % kg m^2
p.b = 0.1;  % friction N/(m/s)
p.g = 9.81;  % gravity m/s/s

% Initial State
x0 = [0;0;0;0];
b.x0.u = x0;
b.x0.l = x0;

% Goal State
xf = [0.5; pi; 0; 0];
tf = 2;
b.xf.u = xf;
b.xf.l = xf;
b.tf.u = tf;
b.tf.l = tf;

% Constraints
ulim = 30; % Bound on force (N)
xlim = 0.6; % Bound on x (m)
b.u.u = ulim;  
b.u.l = -ulim;
b.x.u = [xlim; 2*pi; inf; inf];
b.x.l = [-xlim; -2*pi; -inf; -inf];

% Cost function
costFcn = @(t,x,u) u.^2;

% Dynamics
p.m1 = p.M;
p.m2 = p.m;
p.l = p.ell;
dynFcn = @(t,x,u) cartPoleDynamics(x,u,p);
% dynFcn = @(x,u) invertedPendulumDynamics(x,u,p);

% Setup Params
nSegments = 20;
method = 'trapezoid';
method = 'hermite-simpson';

% Initial Guess
guess.time = [0, tf];
guess.state = [b.x0.l, b.xf.l];
guess.control = [0;0];

%% Solve the problem
soln = direct_collocation(costFcn, dynFcn, b, nSegments, method, guess);

%%%% Unpack the simulation
t = linspace(soln.grid.time(1), soln.grid.time(end), 60);
z = soln.interp.state(t);
u = soln.interp.control(t);

%% Plots:

%%%% Draw Trajectory:
[p1,p2] = cartPoleKinematics(z,p);

figure(2); clf;
nFrame = 20;  %Number of frames to draw
drawCartPoleTraj(t,p1,p2,nFrame);

% Then we can plot an estimate of the error along the trajectory
figure(5); clf;

cc = soln.interp.collCst(t);

subplot(2,2,1);
plot(t,cc(1,:))
title('Collocation Error:   dx/dt - f(t,x,u)')
ylabel('d/dt cart position')

subplot(2,2,3);
plot(t,cc(2,:))
xlabel('time')
ylabel('d/dt pole angle')

idx = 1:length(soln.info.error);
subplot(2,2,2); hold on;
plot(idx,soln.info.error(1,:),'ko');
title('State Error')
ylabel('cart position')

subplot(2,2,4); hold on;
plot(idx,soln.info.error(2,:),'ko');
xlabel('segment index')
ylabel('pole angle');

[ax,l1,l2] = plotyy(t,z',t,cc(1:2,:));
legend({'position','velocity','\theta','\omega','pos-error','\theta-error'})
grid on
xlabel('time')
ylabel(ax(1),'state values')
ylabel(ax(2),'error values')
l2(1).LineStyle = '--';
l2(2).LineStyle = ':';
l2(1).Color = 'r';
l2(2).Color = 'r';



%% Gif
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'testAnimated.gif';
fig = figure(50);
set(gcf,'color','w')
[xLow, xUpp, yLow, yUpp] = getBounds(p1,p2);
for i = 1:length(p1)
    % Draw plot for y = x.^n
    drawCartPoleAnim(fig,[p1(:,i);p2(:,i)],t,xLow, xUpp, yLow, yUpp);
    drawnow
    frame = getframe(fig);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    
    if u == 1
        images = zeros(size(A,1),size(A,2),1,length(p1));
    end
    images(:,:,1,i) = A;
%     
%     if i == 1
%         gif('mygif.gif','frame',fig,'DelayTime',.00001)
%     else 
%         gif('frame',fig,'DelayTime',.00001);
%     end 
end
imwrite(images,map,'testgif.gif','DelayTime',t(2),'LoopCount',Inf);

function [xLow, xUpp, yLow, yUpp] = getBounds(p1,p2)
%
% Returns the upper and lower bound on the data in val
%

val = [p1,p2];
xLow = min(val(1,:));
xUpp = max(val(1,:));
yLow = min(val(2,:));
yUpp = max(val(2,:));

end








