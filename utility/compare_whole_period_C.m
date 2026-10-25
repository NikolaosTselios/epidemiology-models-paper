%% compare_whole_period_C.m
% Standalone script - does NOT modify any existing file. Fits SIR, SEIR
% and SEIRS ONCE EACH over the ENTIRE study period (2020-04-22 to
% 2022-07-10) instead of the usual 27 disjoint 30-day windows, using the
% existing, unmodified fit_SIR_model_C / fit_SEIR_model_C /
% fit_SEIRS_model_C (they already accept any start_date/end_date, so no
% change to them is needed).
%
% WHY: quantifies, with one real number per model, how much worse a
% single fixed set of parameters does over the whole 2+ year period
% compared to the sum of AICc over the 27 independently-fitted windows.
% This is expected to fail badly -- a single SIR/SEIR/SEIRS curve cannot
% reproduce several distinct epidemic waves (original strain, Delta,
% Omicron) with one fixed (beta,gamma,...) -- and this script gives you
% that failure as a concrete number and a picture, instead of only a
% qualitative argument.
%
% Note on AICc here: with T ~ 810 days (n = 2T ~ 1620), the small-sample
% correction term becomes negligible (n >> k), so AIC and AICc end up
% almost equal in this run -- unlike the 30-day windows, where the
% correction mattered.

clear; close all; clc;
rng(1,'twister');

%% ---- Settings ----
csv_file   = fullfile('Data','Covid_19_Greece.csv');
country    = 'Greece';
N          = 1e7;
weights    = [0.9 0.1];
normalize  = true;
num_trials = 800;   % only 3 fits total here (not 27x3), so full 800 is cheap

start_date = datetime(2020,4,22);
end_date   = datetime(2022,7,10);

out_dir = fullfile('results_whole_period');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

% Reference: sum of AICc over the 27 independently-fitted 30-day windows,
% from your existing Table 1 (run_all_C.m output). Update these three
% numbers if you re-run the main pipeline and they change.
ref_AICc_SIR   = 22682.79;
ref_AICc_SEIR  = 22347.08;
ref_AICc_SEIRS = 22481.93;

dataTable = readtable(csv_file);

fprintf('\n===================================================\n');
fprintf(' SINGLE FIT OVER THE WHOLE PERIOD (%s to %s)\n', ...
    datestr(start_date,'yyyy-mm-dd'), datestr(end_date,'yyyy-mm-dd'));
fprintf(' country=%s | num_trials=%d\n', country, num_trials);
fprintf('===================================================\n');

%% ---- Fit each model once, over the whole period ----
fprintf('\nFitting SIR ...\n');
resSIR = fit_SIR_model_C(dataTable, start_date, end_date, N, weights, normalize, num_trials);

fprintf('\nFitting SEIR ...\n');
resSEIR = fit_SEIR_model_C(dataTable, start_date, end_date, N, weights, normalize, num_trials);

fprintf('\nFitting SEIRS ...\n');
resSEIRS = fit_SEIRS_model_C(dataTable, start_date, end_date, N, weights, normalize, num_trials);

results  = {resSIR, resSEIR, resSEIRS};
models   = {'SIR','SEIR','SEIRS'};
ref_AICc = [ref_AICc_SIR, ref_AICc_SEIR, ref_AICc_SEIRS];

%% ---- Readable summary ----
fprintf('\n===================================================\n');
fprintf(' RESULT: single whole-period fit vs. sum of 27 windows\n');
fprintf('===================================================\n');
fprintf('  %-6s  %10s  %10s  %14s  %14s  %9s  %8s\n', ...
    'Model','RMSE_I','RMSE_R','AICc(whole)','AICc(27 wins)','Ratio','R0');

for m = 1:3
    r = results{m};
    rmse_I = sqrt(r.RSS_I / r.T);
    rmse_R = sqrt(r.RSS_R / r.T);
    ratio  = r.AICc / ref_AICc(m);
    fprintf('  %-6s  %10.1f  %10.1f  %14.1f  %14.2f  %8.1fx  %8.2f\n', ...
        models{m}, rmse_I, rmse_R, r.AICc, ref_AICc(m), ratio, r.R0);
end

fprintf('\n"Ratio" = AICc(single whole-period fit) / AICc(sum of 27 windows).\n');
fprintf('A ratio much greater than 1 means a single fixed parameter set is\n');
fprintf('far worse than letting the parameters change every 30 days, as the\n');
fprintf('main pipeline does.\n');
fprintf('===================================================\n');

%% ---- Picture: how badly the whole-period fit misses the real curve ----
fig = figure('Visible','off','Color','w','Position',[100 100 1500 700]);

for m = 1:3
    r = results{m};

    subplot(2,3,m);
    plot(r.dates, r.infected_real, 'k-', 'LineWidth', 1.8); hold on;
    plot(r.dates, r.infected_model, '--', 'Color',[0.85 0 0], 'LineWidth', 1.6);
    title(models{m}, 'FontWeight','bold');
    if m==1, ylabel('Infected'); end
    grid on; box off;

    subplot(2,3,3+m);
    plot(r.dates, r.removed_real, 'k-', 'LineWidth', 1.8); hold on;
    plot(r.dates, r.removed_model, '--', 'Color',[0 0.35 0.9], 'LineWidth', 1.6);
    if m==1, ylabel('Removed'); end
    xlabel('Date');
    grid on; box off;
end

try
    sgtitle('Single fit over the whole study period (real = black, model = dashed)', ...
        'FontWeight','bold');
catch
end

fig_path = fullfile(out_dir, sprintf('whole_period_fit_%s.png', country));
print(fig, fig_path, '-dpng', '-r200');
close(fig);

fprintf('\nFigure saved -> %s\n', fig_path);
