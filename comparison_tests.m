main;

methods = {'trapezoid','hermite-simpson'};
nSegmentsList = [2,5,10,15,20,30,40,50,75,100];

solutions = cell(2,length(nSegmentsList));
for m = 1:2
    method = methods{m};
    for n = 1:length(nSegmentsList)
        nSegments = nSegmentsList(n);
        soln = direct_collocation(costFcn, dynFcn, b, nSegments, method, guess, false);
        solutions{m,n} = soln;
    end
end

%%
solutions_guess_end = cell(2,length(nSegmentsList));
progressbar('Method','nSegments')
for m = 1:2
    method = methods(m);
    for n = 1:length(nSegmentsList)-1
        nSegments = nSegmentsList(end);
        guess_soln = solutions{m,n};
        soln = direct_collocation(costFcn, dynFcn, b, nSegments, method, guess_soln.grid, false);
        solutions_guess_end{m,n} = soln;
        progressbar((m-1)/2, (n-1)/length(nSegmentsList))
    end
end

%% Use init from low segment collocation
solutions_guess = cell(2,length(nSegmentsList));
guess_soln = solutions{1,2};
guess.t = guess_soln.grid.time;
guess.x = guess_soln.grid.state;
guess.u = guess_soln.grid.control;
progressbar('Method','nSegments')
for m = 1:2
    for n = 1:length(nSegmentsList)
        nSegments = nSegmentsList(n);
        if isempty(solutions_guess{m,n})
            soln = direct_collocation(costFcn, dynFcn, b, nSegments, method, guess, false);
            solutions_guess{m,n} = soln;
        end
        progressbar(m/2, n/length(nSegmentsList))
    end
end

%% Plots
load comparison_data.mat
% Original Guess
t = linspace(0,tf,100);
for m = 1:2
    method = methods{m};
    for n = 1:length(nSegmentsList)
        nSegments = nSegmentsList(n);
        soln = solutions{m,n};
        runtime{m}(n) = soln.info.t_computation;
        maxerror{m}(n) = soln.info.maxError;
        objVal{m}(n) = soln.info.objVal;
        status{m}(n) = soln.info.exitFlag;
        traj{m,n} = soln.interp.state(t);
    end
end

% Runtime Comparison
figure(1); subplot(3,1,1); hold on; grid on;
plot(nSegmentsList, runtime{1}, 'linewidth', 2)
plot(nSegmentsList, runtime{2}, 'linewidth', 2)
xlabel('Number of Segments')
ylabel('Runtime (sec)')
legend(methods,'Location','northwest')

% Error Comparison
subplot(3,1,2); 
semilogy(nSegmentsList, maxerror{1}, 'linewidth', 2);  hold on;  grid on
semilogy(nSegmentsList, maxerror{2}, 'linewidth', 2)
xlabel('Number of Segments')
ylabel('Max error')
legend(methods)

% Error Comparison
subplot(3,1,3); 
semilogy(nSegmentsList, objVal{1}, 'linewidth', 2);  hold on;  grid on
semilogy(nSegmentsList, objVal{2}, 'linewidth', 2)
xlabel('Number of Segments')
ylabel('Objective Value')
legend(methods)

% Trajectories
figure; hold on
lbl = {'position', 'velocity', 'angle (rad)', 'ang vel (rad/s)'};
for m = 1:2
    for n = 1:length(nSegmentsList)
        for i = 1:4
            ind = (i-1)*2+m;
            ax(m,i) = subplot(4,2,ind); hold on; grid on
            plot(t,traj{m,n}(i,:))
            xlabel('time (sec)')
            ylabel(lbl{i})
        end
    end
end
for i = 1:4
    ylim1 = ax(1,i).YLim;
    ylim2 = ax(2,i).YLim;
    if diff(ylim1) > diff(ylim2)
        ylim(ax(2,i),ylim1);
    else
        ylim(ax(1,i),ylim2);
    end
end
title(ax(1,1),'Trapezoidal')
title(ax(2,1),'Hermite-Simpson')

%% Guess 2
for m = 1:2
    method = methods{m};
    for n = 1:length(nSegmentsList)
        nSegments = nSegmentsList(n);
        soln = solutions_guess{m,n};
        runtime{2,m}(n) = soln.info.t_computation;
        maxerror{2,m}(n) = soln.info.maxError;
        objVal{2,m}(n) = soln.info.objVal;
        status{2,m}(n) = soln.info.exitFlag;
    end
end

% Runtime Comparison
figure(1); hold on; grid on;
plot(nSegmentsList, runtime{2,1}, 'linewidth', 2)
plot(nSegmentsList, runtime{2,2}/2, 'linewidth', 2)
xlabel('Number of Segments')
ylabel('Runtime (sec)')
legend(methods)

%% Guess End
for m = 1:2
    method = methods{m};
    for n = 1:length(nSegmentsList)-1
        nSegments = nSegmentsList(n);
        soln = solutions_guess_end{m,n};
        runtime{3,m}(n) = soln.info.t_computation;
        maxerror{3,m}{m}(n) = soln.info.maxError;
        objVal{3,m}(n) = soln.info.objVal;
        status{3,m}(n) = soln.info.exitFlag;
    end
end

% Runtime Comparison
figure(2); hold on; grid on;
plot(nSegmentsList(1:end-1), runtime{1,2}(end) - runtime{3,1}, 'linewidth', 2)
plot(nSegmentsList(1:end-1), runtime{1,2}(end) - runtime{3,2}, 'linewidth', 2)
xlabel('Number of Segments')
ylabel('Runtime Savings (sec)')
legend(methods)
title('Initialing with Results of Smaller Problem')
