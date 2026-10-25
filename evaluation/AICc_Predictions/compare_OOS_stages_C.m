%% compare_OOS_stages_C.m
% Standalone script - does NOT modify any existing file. Runs the
% out-of-sample validation (see fit_*_model_OOS_C.m) across three
% training-set sizes ("stages"), all in ONE run, and prints a single,
% easy-to-read summary table - no CSV needed to see the result.
%
% Stages (all applied to the same 27 disjoint 30-day windows):
%   Early  -> train_days = 10  (test ~20 days)
%   Middle -> train_days = 15  (test ~15 days)
%   Late   -> train_days = 20  (test ~10 days)
%
% Same spirit as the early/middle/late staged validation in Yang et al.
% (2021, Epidemics, "Rational evaluation of various epidemic models
% based on the COVID-19 data of China") - adapted to this study's
% structure: instead of growing the training portion within a single
% continuous wave, we grow it within each of the 27 independent 30-day
% windows and aggregate the result across all of them.
%
% REQUIRES (unmodified, already on path): fit_SIR_model_OOS_C,
% fit_SEIR_model_OOS_C, fit_SEIRS_model_OOS_C, generate_fitting_periods.

clear; close all; clc;
rng(1,'twister');

%% ---- Settings ----
csv_file   = fullfile('Data','Covid_19_Greece.csv');
country    = 'Greece';
N          = 1e7;
weights    = [0.9 0.1];
normalize  = true;
num_trials = 800;      % reduced from 800: this script fits 3x as many
                        % times as a single-stage run. Rankings/patterns
                        % are stable well below 800 (same reasoning as
                        % weight_sensitivity.m). Raise if you have time.
step_days  = 30;

stage_names      = {'Early','Middle','Late'};
stage_train_days = [10, 15, 20];

start_date = "2020-04-22";
end_date   = "2022-07-10";

out_dir = fullfile('results_OOS_stages');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

fprintf('\n===================================================\n');
fprintf(' OUT-OF-SAMPLE VALIDATION ACROSS 3 TRAINING STAGES\n');
fprintf(' country=%s | num_trials=%d\n', country, num_trials);
fprintf(' Early(train=10d) / Middle(train=15d) / Late(train=20d)\n');
fprintf('===================================================\n');

periods = generate_fitting_periods(start_date, end_date, step_days);
nP = size(periods,1);
dataTable = readtable(csv_file);

models  = {'SIR','SEIR','SEIRS'};
nStages = numel(stage_train_days);
nModels = numel(models);

all_rows = cell(1, nStages);   % raw per-window results, one cell per stage

%% ---- Run all three stages ----
for s = 1:nStages
    train_days = stage_train_days(s);
    fprintf('\n--- Stage %d/%d: %s (train_days=%d) ---\n', s, nStages, stage_names{s}, train_days);

    rows = {};
    for i = 1:nP
        sd = datetime(periods{i,1});
        ed = datetime(periods{i,2});

        rS  = fit_SIR_model_OOS_C(  dataTable, sd, ed, N, weights, normalize, num_trials, train_days);
        rE  = fit_SEIR_model_OOS_C( dataTable, sd, ed, N, weights, normalize, num_trials, train_days);
        rES = fit_SEIRS_model_OOS_C(dataTable, sd, ed, N, weights, normalize, num_trials, train_days);

        if isempty(rS) || isempty(rE) || isempty(rES)
            fprintf('  window %2d/%d: skipped (too short for this split)\n', i, nP);
            continue;
        end

        aicc_train = [rS.AICc_train, rE.AICc_train, rES.AICc_train];
        rmse_test  = [rS.RMSE_test_I, rE.RMSE_test_I, rES.RMSE_test_I];

        [~, i_aicc] = min(aicc_train);
        [~, i_rmse] = min(rmse_test);

        rows(end+1,:) = {i, aicc_train(1), aicc_train(2), aicc_train(3), ...
                          rmse_test(1), rmse_test(2), rmse_test(3), ...
                          i_aicc, i_rmse}; %#ok<SAGROW>

        fprintf('  window %2d/%d done\n', i, nP);
    end

    all_rows{s} = rows;
end

%% ---- Build and print the summary: one readable block per stage ----
fprintf('\n\n');
fprintf('=================================================================================\n');
fprintf(' SUMMARY  (country=%s, num_trials=%d, %d windows total)\n', country, num_trials, nP);
fprintf('=================================================================================\n');

summary_rows = {};   % kept for the optional CSV export at the end

for s = 1:nStages
    rows  = all_rows{s};
    nUsed = size(rows,1);

    aicc_sum  = zeros(1,nModels);
    rmse_all  = cell(1,nModels);
    win_count = zeros(1,nModels);
    match_count = 0;

    for r = 1:nUsed
        aicc_train = [rows{r,2}, rows{r,3}, rows{r,4}];
        rmse_test  = [rows{r,5}, rows{r,6}, rows{r,7}];
        i_aicc = rows{r,8};
        i_rmse = rows{r,9};

        aicc_sum = aicc_sum + aicc_train;
        for m = 1:nModels
            rmse_all{m}(end+1) = rmse_test(m); %#ok<AGROW>
        end
        win_count(i_rmse) = win_count(i_rmse) + 1;
        if i_aicc == i_rmse
            match_count = match_count + 1;
        end
    end

    rmse_median = cellfun(@median, rmse_all);
    win_pct     = 100 * win_count / max(nUsed,1);
    match_pct   = 100 * match_count / max(nUsed,1);

    fprintf('\n--- %s stage (train_days=%d) | windows used: %d/%d | AICc-vs-RMSE match: %.0f%% ---\n', ...
        stage_names{s}, stage_train_days(s), nUsed, nP, match_pct);
    fprintf('  %-6s  %14s  %14s  %10s\n', 'Model', 'AICc(train)', 'RMSE(test)', 'Win rate');
    fprintf('  %-6s  %14s  %14s  %10s\n', '', 'sum', 'median', '');
    for m = 1:nModels
        fprintf('  %-6s  %14.1f  %14.1f  %9.0f%%\n', models{m}, aicc_sum(m), rmse_median(m), win_pct(m));
    end

    for m = 1:nModels
        summary_rows(end+1,:) = {stage_names{s}, stage_train_days(s), nUsed, models{m}, ...
            aicc_sum(m), rmse_median(m), win_count(m), win_pct(m), match_pct}; %#ok<SAGROW>
    end
end

fprintf('\n=================================================================================\n');
fprintf(' HOW TO READ THIS TABLE:\n');
fprintf('  AICc(train) sum    -> lower = better fit to the days the model actually saw.\n');
fprintf('  RMSE(test) median  -> lower = better prediction on days it never saw. Median is\n');
fprintf('                        used (not mean) so a couple of huge-scale windows do not\n');
fprintf('                        drown out the typical-window picture.\n');
fprintf('  Win rate           -> %% of windows where this model had the single lowest\n');
fprintf('                        RMSE(test), i.e. was the actual best forecaster.\n');
fprintf('  Match              -> %% of windows where the model AICc(train) picked as best\n');
fprintf('                        was ALSO the model with the lowest RMSE(test) - i.e. how\n');
fprintf('                        often AICc correctly predicted the best forecaster.\n');
fprintf('=================================================================================\n');

%% ---- Same numbers as CSV, kept only as backup / appendix material ----
Tsum = cell2table(summary_rows, 'VariableNames', ...
    {'Stage','TrainDays','WindowsUsed','Model','AICc_train_sum','RMSE_test_median', ...
     'WinCount','WinRate_pct','MatchRate_pct'});
outfile = fullfile(out_dir, sprintf('oos_stages_summary_%s.csv', country));
writetable(Tsum, outfile);
fprintf('\n(Same numbers also saved to %s - only needed if you want an appendix table.)\n', outfile);
