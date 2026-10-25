%% run_all_C.m  —  MASTER pipeline, one clean run from scratch
% =====================================================================
% Runs the whole chain end-to-end and gathers every paper deliverable
% into ONE timestamped folder, without modifying any validated function.
%
% Chain:
%   1. Fitting (SIR, SEIR, SEIRS) over all 30-day windows
%   2. Headline model comparison  (stored AICc totals -> Table 1)
%   3. Consolidated per-window metrics + candidate windows
%   4. Residual diagnostics (box plots, heatmaps, ranking, per-window I/O)
%   5. The 4 representative trajectory plots (windows 5, 13, 14, 20)
%   5b. ALL 27 trajectory plots (QA only, not for the paper)
%   6. A printed manifest of everything created
%
% Output folder:  results_FINAL/<timestamp>/
%
% Run from the PROJECT ROOT (the folder that contains Data/, results/,
% and has scripts/ models/ evaluation/ utility/ on the path via startup).
% =====================================================================

clear; close all; clc;
rng(1,'twister');          % fixed seed -> this final run is reproducible

%% ---- Configuration (identical to run_fitting_models_C) ----
N          = 1e7;
weights    = [0.9 0.1];
normalize  = true;
num_trials = 800;          % <-- set to ~50 for a quick pipeline test, 800 for the real run
step_days  = 30;
country    = 'Greece';
csv_file   = fullfile('Data','Covid_19_Greece.csv');
plotW      = [5 13 14 20]; % the 4 representative windows we chose

%% ---- One clean, timestamped output folder ----
stamp   = datestr(now,'yyyy-mm-dd_HHMM');
outRoot = fullfile('results_FINAL', stamp);
dir_fit      = fullfile(outRoot,'01_fitting');
dir_metr     = fullfile(outRoot,'02_metrics');
dir_diag     = fullfile(outRoot,'03_diagnostics');
dir_traj     = fullfile(outRoot,'04_trajectories');
dir_traj_all = fullfile(outRoot,'05_trajectories_all');   % QA: all windows
mkdir_if(dir_fit); mkdir_if(dir_metr); mkdir_if(dir_diag); mkdir_if(dir_traj);
mkdir_if(dir_traj_all);

fprintf('\n================ MASTER RUN ================\n');
fprintf('Output -> %s\n', outRoot);
fprintf('num_trials = %d  (this may take a while)\n', num_trials);
fprintf('===========================================\n');

%% ===== 1. Fitting =====
% compare_*_periods_C save canonical .mat under results/Results_<M>/Greece/
fprintf('\n[1/6] Fitting SIR / SEIR / SEIRS ...\n');
results_SIR   = compare_SIR_periods_C(  csv_file, country, N, weights, normalize, num_trials, step_days);
results_SEIR  = compare_SEIR_periods_C( csv_file, country, N, weights, normalize, num_trials, step_days);
results_SEIRS = compare_SEIRS_periods_C(csv_file, country, N, weights, normalize, num_trials, step_days);

assert(isequal(results_SIR.period_idx, results_SEIR.period_idx, results_SEIRS.period_idx), ...
    'Models were not evaluated on identical windows.');

% copy the fitting outputs (.mat + .csv) into the clean folder
copy_glob(dir_fit, fullfile('results','Results_SIR',  country,'SIR_all_periods_results1_C.*'));
copy_glob(dir_fit, fullfile('results','Results_SEIR', country,'SEIR_all_periods_results1_C.*'));
copy_glob(dir_fit, fullfile('results','Results_SEIRS',country,'SEIRS_all_periods_results1_C.*'));

%% ===== 2. Headline model comparison (Table 1) =====
% Uses the STORED per-window AICc summed inside each compare_* (summary.AICc_tot),
% exactly as run_fitting_models_C does. This is the paper's Table 1.
fprintf('\n[2/6] Building Table 1 (model comparison) ...\n');
models   = {'SIR';'SEIR';'SEIRS'};
AIC_tot  = [results_SIR.AIC_tot;  results_SEIR.AIC_tot;  results_SEIRS.AIC_tot];
AICc_tot = [results_SIR.AICc_tot; results_SEIR.AICc_tot; results_SEIRS.AICc_tot];
dAICc    = AICc_tot - min(AICc_tot);
wA       = exp(-0.5*dAICc);  wA = wA/sum(wA);
Table1   = table(models, AIC_tot, AICc_tot, dAICc, wA, ...
    'VariableNames', {'Model','AIC','AICc','DeltaAICc','AkaikeWeight'});
Table1   = sortrows(Table1,'AICc','ascend');
disp(Table1);
writetable(Table1, fullfile(dir_metr,'Table1_model_comparison.csv'));

%% ===== 3. Consolidated per-window metrics + candidate windows =====
% collect_all_metrics reloads the canonical .mat, writes Metrics_long.csv
% and Window_summary.csv to the current folder, and prints candidate windows.
fprintf('\n[3/6] Consolidated metrics (collect_all_metrics) ...\n');
collect_all_metrics();
move_glob(dir_metr, 'Metrics_long.csv');
move_glob(dir_metr, 'Window_summary.csv');

%% ===== 4. Residual diagnostics =====
% Box plots + heatmaps + rank-based ranking (all 27 windows) and per-window
% ACF/QQ/residuals-vs-fitted for the 4 chosen windows -> all into dir_diag.
fprintf('\n[4/6] Residual diagnostics ...\n');
diag_SIR   = compute_residual_diagnostics_C(results_SIR,   'PlotWindows', plotW, 'OutputDir', dir_diag);
diag_SEIR  = compute_residual_diagnostics_C(results_SEIR,  'PlotWindows', plotW, 'OutputDir', dir_diag);
diag_SEIRS = compute_residual_diagnostics_C(results_SEIRS, 'PlotWindows', plotW, 'OutputDir', dir_diag);
comparison = compare_residual_diagnostics_C({diag_SIR,diag_SEIR,diag_SEIRS}, {'SIR','SEIR','SEIRS'}, ...
    'SavePlots', true, 'ShowPlots', false, 'OutputDir', dir_diag);

%% ===== 5. The 4 representative trajectory plots =====
fprintf('\n[5/6] Trajectory plots for windows %s ...\n', mat2str(plotW));
plot_windows_trajectories(results_SIR, results_SEIR, results_SEIRS, plotW, dir_traj);

%% ===== 5b. ALL windows, for full-picture QA (not for the paper) =====
n_windows = numel(results_SIR.all_results);
allW = 1:n_windows;
fprintf('\n[5b/6] Trajectory plots for ALL %d windows (QA) ...\n', n_windows);
plot_windows_trajectories(results_SIR, results_SEIR, results_SEIRS, allW, dir_traj_all);

%% ===== 6. Manifest =====
fprintf('\n[6/6] Manifest\n');
print_manifest(outRoot);
fprintf('\n================ DONE ================\n');
fprintf('Everything is in: %s\n', outRoot);
fprintf('=====================================\n');


%% =====================================================================
%%  Local helper functions
%% =====================================================================

function mkdir_if(d)
if ~exist(d,'dir'), mkdir(d); end
end

function copy_glob(destDir, pattern)
% copy every file matching pattern into destDir (warn if none)
L = dir(pattern);
if isempty(L)
    warning('run_all:copy','No file matched: %s', pattern); return;
end
for i = 1:numel(L)
    copyfile(fullfile(L(i).folder, L(i).name), destDir);
end
end

function move_glob(destDir, pattern)
% move every file matching pattern into destDir (warn if none)
L = dir(pattern);
if isempty(L)
    warning('run_all:move','No file matched: %s', pattern); return;
end
for i = 1:numel(L)
    movefile(fullfile(L(i).folder, L(i).name), destDir);
end
end

function print_manifest(root)
L = dir(fullfile(root,'**','*'));
L = L(~[L.isdir]);
fprintf('\nFiles created under %s :\n', root);
for i = 1:numel(L)
    rel = strrep(fullfile(L(i).folder, L(i).name), [root filesep], '');
    fprintf('   %-55s  %8.1f KB\n', rel, L(i).bytes/1024);
end
end


function plot_windows_trajectories(rSIR, rSEIR, rSEIRS, wins, destDir)
% 2x3 comparison figure per window: Infected (top) / Removed (bottom),
% one column per model. Style matches plot_model_comparisons1.
R = {rSIR.all_results, rSEIR.all_results, rSEIRS.all_results};
mdl = {'SIR','SEIR','SEIRS'};

for w = wins(:)'
    a = R{1}{w}; b = R{2}{w}; c = R{3}{w};   % SIR, SEIR, SEIRS results for this window
    res = {a,b,c};

    % shared axis limits across the three models (fair comparison)
    Iall = [a.infected_real(:); a.infected_model(:); b.infected_model(:); c.infected_model(:)];
    Rall = [a.removed_real(:);  a.removed_model(:);  b.removed_model(:);  c.removed_model(:)];
    yI = [min(Iall) max(Iall)];  yR = [min(Rall) max(Rall)];

    fig = figure('Visible','off','Color','w','Position',[100 100 1500 850]);

    for m = 1:3
        r = res{m};
        % --- Infected (top row) ---
        subplot(2,3,m);
        plot(r.dates, r.infected_real,  'k-',  'LineWidth', 2.6); hold on;
        plot(r.dates, r.infected_model, '--', 'Color',[0.85 0 0], 'LineWidth', 1.8);
        title(mdl{m},'FontWeight','bold','FontSize',14);
        if m==1, ylabel('Infected'); end
        ylim(yI); grid on; box off; set(gca,'GridAlpha',0.08);

        % --- Removed (bottom row) ---
        subplot(2,3,3+m);
        plot(r.dates, r.removed_real,  'k-',  'LineWidth', 2.6); hold on;
        plot(r.dates, r.removed_model, '--', 'Color',[0 0.35 0.9], 'LineWidth', 1.8);
        title(mdl{m},'FontWeight','bold','FontSize',14);
        if m==1, ylabel('Removed'); end
        xlabel('Date');
        ylim(yR); grid on; box off; set(gca,'GridAlpha',0.08);
    end

    try
        sgtitle(sprintf('Model Comparison — Window %d (%s to %s)', w, ...
            datestr(a.dates(1),'yyyy-mm-dd'), datestr(a.dates(end),'yyyy-mm-dd')), ...
            'FontWeight','bold','FontSize',15);
    catch
    end

    fname = fullfile(destDir, sprintf('Trajectory_Window_%02d.png', w));
    print(fig, fname, '-dpng', '-r300');
    close(fig);
    fprintf('   saved %s\n', fname);
end
end