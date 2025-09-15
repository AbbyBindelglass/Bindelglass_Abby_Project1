% AMS 595 MATLAB Project 1 
% Abby Bindelglass

% 1. 

clc; % clear command window text
clear; % removes all variables from workspace
close all; % closes open figure windows

% uniformly select random points (x,y) in [0,1]^2

% a point lies in the quarter-unit circle if x^2 + y^2 <= 1

% probability of being inside equals area(quarter circle)/area(square) =
% (pi/4)/1, pi = 4 * (inside_count / total_samples)

% use for loop and log the running estimate at checkpoints
% plot pi_hat vs N and |error| vs N
% benchmark independent runs over range of N to show error vs time and
% error vs N

real_pi = pi; % store real value of pi

N_total = 5e5; % Total number of random points
checkpoint = 500; % log every 500 iterations

num_logs = floor(N_total / checkpoint); % how many times logged during trajectory
hits = 0; % points inside quarter circle

N_log = zeros(num_logs,1); % N values to record at
pi_log = zeros(num_logs,1); % running estimates of pi at N
err_log = zeros(num_logs,1); % |pi_true - pi_est| at N

t0 = tic; % timer
idx = 0; % initialize index to track how many times logs happen
for n = 1:N_total % iterate N_total times, simulate one point per iteration
    x = rand; % random x coordinate from U(0,1)
    y = rand; % random y coordinate from U(0,1)
    if x*x + y*y <= 1 % Check if (x,y) falls inside quarter circle
        hits = hits + 1; % if point is inside quarter circle, add one to inside_count
    end
    % every 'checkpoint' iterations, record current N, pi estimate, and error
    if mod(n,checkpoint) == 0
        idx = idx + 1; % add 1 to log index to add new entry
        N_log(idx) = n; % store sample size where logging is happening
        pi_hat = 4 * (hits / n); % monte carlo estimator for pi at N
        pi_log(idx) = pi_hat; % store running estimate
        err_log(idx) = abs(real_pi - pi_hat); % compute absolute error against true pi
    end
end

t_for = toc(t0); % total time

% plot of running estimate of pi vs N
figure; % new figure window
plot(N_log, pi_log, 'LineWidth', 1.2); % plot recorded N values (x-axis) against recorded pi estimates (y-axis)
hold on; % hold on so something else can be added to same plot
yline(real_pi, '--'); % dashed pine at real value of pi with label
xlabel('N'); % x axis label for number of samples used
ylabel('\pi estimate'); % y axis label for estimate of pi
title(['\pi vs N (N_{total} = ' num2str(N_total) ... % title, including N_total and elapsed time
    ', time = ' num2str(t_for,'%.2f') ' s)']);
grid on; % add grid

% plot of absolute error vs N
figure; % new figure
plot(N_log,err_log,'LineWidth',1.2); % plot absolute error at each checkpoint as function of N
set(gca,'XScale','log'); % x axis to log scale
xlabel('N (log)'); % label x axis
ylabel('| \pi - \pi_{est} |'); % label y axis
title('Error vs N'); % title for error plot
grid on; % grid lines for readability

fprintf('Task 1: N_total = %d, pi_hat = %.6f |error|=%.6f, time=%.2fs\n', ...
    N_total, pi_log(end), err_log(end), t_for); % summary line for trajectory run

N_list = [10, 50, 100, 500, 1e3, 5e3, 1e4, 5e4, 1e5]; % independent sample sizes (small --> large) for accuracy/runtime tradeoff 

% preallocate arrays to store runtime, estimate, and error for each N in N_list
time_list = zeros(size(N_list)); % elapsed seconds for each run
err_list = zeros(size(N_list)); % abs(error) for each run

for k = 1:numel(N_list) % loop over each independent sample size
    Nk = N_list(k); % read current total sample size
    c_in = 0; % reset inside counter

    t1 = tic; % start timing just this run
    for j = 1:Nk % draw Nk points and count hits
        x = rand; y = rand; % draw (x,y) pair from [0,1]^2
        if x*x + y*y <= 1 % if point lies inside quarter circle
            c_in = c_in + 1; % add 1 to hit counter
        end
    end
    time_list(k) = toc(t1); % capture elapsed seconds
    pi_k = 4 * (c_in / Nk); % compute monte carlo estimator
    err_list(k) = abs(real_pi - pi_k); % store absolute error

    % print summary line for each run to see scaling across different values of N
    fprintf('N = %8d | pi_hat = %.6f | |err| = %.6f | t = %.3fs\n', ...
        Nk, pi_k, err_list(k),time_list(k));
end

% Plot error vs time
figure; % new figure
plot(time_list,err_list,'o-', 'LineWidth', 1.2); % elapsed time on x axis, absolute error on y, points connected with lines
xlabel('Time (s)'); % label x axis as time in seconds
ylabel('| \pi - \pi_{est} |'); % label y axis as absolute error
title('Precision vs time'); % title showing plot visualizes tradeoff between precision and compute time
grid on; % add grid

% Plot error vs N with log axes
figure; % new figure
plot(N_list, err_list, 'o-', 'LineWidth', 1.2); % error (y) vs N (x) as line with markers
set(gca, 'XScale', 'log', 'YScale', 'log'); % log scaling on both axes
xlabel('N (log)'); ylabel(' | \pi - \pi_{est} | (log)'); % label x axis for N and y axis for error
title('Error vs N (independent runs)'); % title of plot
grid on; % add grid

% 2.
% Monte Carlo estimate of pi using a While Loop 
% run monte carlo batches until estimate stabilizes to s sig figs
% do not compare true pi in stopping rule

clc; clear; % clear environment 
sig_list = [2,3,4]; % target sig figs to achieve
batch_size = 1e4; % number of points to simulate per while-iteration
M = 5; % set M (larger M = more confidence, but takes longer)

% safety caps to prevent long runs
max_batches = 5e4; % max number of while-iterations
max_points = 1e9; % max total points

for s = sig_list % for each target s.f. in sig_list, run stability loop
    % reset counters
    total_in = 0; % total inside
    total_n = 0; % total points
    same_cnt = 0; % # consecutive times rounded value stayed same
    last_r = NaN; % previous rounded-to-s estimate

    t2 = tic; % time stability run

    while true % while loop that goes until stability is achieved or safety caps hit

        % Generate vectorized batch of random x and y values
        x = rand(batch_size,1);
        y = rand(batch_size,1);

        total_in = total_in + sum(x.^2 + y.^2 <= 1); % cumulative number of inside points
        total_n = total_n + batch_size; % cumulative number of generated points

        % Current Monte Carlo estimate
        pi_hat = 4 * (total_in / total_n);

        % Round to s significant figures
        r_now = round_sig(pi_hat,s);

        if isequal(r_now,last_r) % if rounded value matches previous
            same_cnt = same_cnt + 1; % increment stability streak
        else
            same_cnt = 1; % reset consecutive count to 1
            last_r = r_now; % update last rounded value
        end

        if same_cnt >= M % if reached M consecutive matches
            break; % exit while loop
        end
    end
    t2 = toc(t2); % record time it took to achieve stability

    % Print summary
    fprintf('Task 2: s=%d s.f., pi_hat=%.*g, N=%d, time=%.2fs\n', s, s, pi_hat, total_n, t2);
    fprintf(' |error| vs MATLAB pi: %.6g\n', abs(pi-pi_hat));
end

% 3.
% function with live visualization
% call task 3 function with target of 3 sig figs

pi3 = mc_pi_visual(3); % calls function defined below, pi3 is final estimate

% Local functions:

function y = round_sig(x,s)

    y = zeros(size(x),'like',x); % preallocate output
    keep = (x~=0) & isfinite(x); % only nonzero, finite values need rounding math
    if any(keep(:)) 
        e = floor(log10(abs(x(keep)))); % decimal exponent
        y(keep) = round(x(keep).*10.^(-e+(s-1))).*10.^(e-(s-1)); % scale, then round, then unscale
    end
    y(isinf(x)) = x(isinf(x)); % keep +/- Inf unchanged
    y(isnan(x)) = x(isnan(x)); % keep NaNs unchanged
end

function pi_est = mc_pi_visual(s,batch_size,M,show_limit)
    % defaults if caller doesn't include
    if nargin<2 || isempty(batch_size), batch_size=1e4;
    end
    if nargin<3 || isempty(M), M=5;
    end
    if nargin<4 || isempty(show_limit), show_limit=5e4;
    end

    total_in = 0; total_n = 0; % counts for hits, total points
    same_cnt = 0; last_r = NaN; % stability status

    % buffers for points drawn
    xin=[];yin=[];xout=[];yout=[];

    % figure setup
    figure('Color','w'); % new figure with white background
    hold on; % so more can be added to the same figure

    % quarter-circle boundary in unit square
    tt = linspace(0, pi/2, 200); % 200 points from 0 to pi/2
    plot(cos(tt),sin(tt),'k-'); % plot boundary arc as thin black line

    % graphics objects for inside/outside
    hIn = plot(NaN, NaN, '.', 'DisplayName','Inside'); % inside pts
    hOut = plot(NaN, NaN, '.', 'DisplayName','Outside'); % outside pts

    axis([0 1 0 1]); axis square; % fix axes to unit square, make aspect ratio 1:1
    grid on; % add grid
    xlabel('x'); ylabel('y'); % label x and y axes
    title('Monte Carlo \pi (live)'); % add title
    legend('Location','southoutside'); % legend below plot

    % text object for live status updates
    txt = text(0.02, 0.98, '', 'Units','normalized',...
    'VerticalAlignment','top');

    while true
        % generate vectorized batch of points uniformly in [0,1]^2
        x = rand(batch_size,1);
        y = rand(batch_size,1);
        in = (x.^2 + y.^2) <= 1; % find out which points are in the quarter circle

        total_in = total_in + sum(in); % update cumulative inside count
        total_n = total_n + batch_size; % update cumulative total points

        % monte carlo estimate
        pi_hat = 4 * (total_in / total_n);

        r_now = round_sig(pi_hat,s); % round to s significant figures to test for stability
        if isequal(r_now, last_r) % if rounded estimate = previous rounded value
            same_cnt = same_cnt + 1; % increase stability streak
        else
            same_cnt = 1; % reset streak
            last_r = r_now; % track new rounded value
        end

        remain = max(0,show_limit - (numel(xin) + numel(xout))); % add only subset of batch to plot buffers
        if remain > 0 % if there is capacity to display more points
            take = min(remain,numel(in)); % determine how many points from batch allowed to display without exceeding cap
            x_b = x(1:take); y_b = y(1:take);
            in_b = in(1:take);
            xin = [xin; x_b(in_b)]; % append inside pts
            yin = [yin; y_b(in_b)]; 
            xout = [xout;x_b(~in_b)]; % append outside pts
            yout = [yout;y_b(~in_b)];

        end

        % push buffers into plotted objects
        set(hIn, 'XData', xin, 'YData', yin);
        set(hOut, 'XData', xout, 'YData', yout);

        % update status text
        set(txt, 'String', sprintf('N = %d\n\\pi\\approx %.*g (to %d s.f.)\nStable %d/%d',... % update on figure status text
            total_n, s, pi_hat, s, min(same_cnt,M),M));

        drawnow limitrate; % reduce overhead

        % stop when stable for M consecutive batches
        if same_cnt >= M
            break;
        end
    end

    pi_est = pi_hat; % save final estimate

    % summary to command window
    fprintf('Task 3: s=%d s.f., pi_hat=%.*g, N=%d\n', s, s, pi_est, total_n);
end