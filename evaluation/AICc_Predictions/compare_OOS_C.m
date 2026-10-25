%% compare_OOS_C.m
% Standalone out-of-sample validation script - does NOT touch any
% existing file. Calls the fit_*_model_OOS_C functions for every disjoint
% 30-day window, with the first train_days days as training and the
% remainder as held-out test. Same spirit as the train/test split in
% Yang et al. (2021, Epidemics): AICc on training data, RMSE on test data.
%
% WHAT THIS CHECKS: whether the model with the lowest AICc (on training
% data) actually has the lowest RMSE on test data - i.e. whether AICc
% correctly predicts which model will generalize better to unseen data,
% not just which fits best in-sample.
%
% REQUIRES: fit_SIR_model_OOS_C, fit_SEIR_model_OOS_C,
% fit_SEIRS_model_OOS_C, and the existing generate_fitting_periods must
% be on the MATLAB path.

clear; close all; clc;
rng(1,'twister');

%% ---- Settings ----
csv_file   = fullfile('Data','Covid_19_Greece.csv');
country    = 'Greece';
N          = 1e7;
weights    = [0.9 0.1];
normalize  = true;
num_trials = 800;      % <-- lower (e.g. 100-200) for a quick test run
step_days  = 30;
train_days = 20;       % first 20 days = training; the rest (~10) = test

start_date = "2020-04-22";
end_date   = "2022-07-10";

out_dir = fullfile('results_OOS');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

fprintf('\n===================================================\n');
fprintf(' OUT-OF-SAMPLE VALIDATION (AICc train vs RMSE test)\n');
fprintf(' country=%s | train_days=%d | num_trials=%d\n', country, train_days, num_trials);
fprintf('===================================================\n');

%% ---- Windows (same generate_fitting_periods as the main pipeline) ----
periods = generate_fitting_periods(start_date, end_date, step_days);
nP = size(periods,1);
fprintf('Total windows: %d\n', nP);

dataTable = readtable(csv_file);

models = {'SIR','SEIR','SEIRS'};
rows = {};

%% ---- Window loop ----
for i = 1:nP

    sd = datetime(periods{i,1});
    ed = datetime(periods{i,2});

    fprintf('\n=== Window %d/%d: %s to %s ===\n', i, nP, ...
        datestr(sd,'yyyy-mm-dd'), datestr(ed,'yyyy-mm-dd'));

    rS  = fit_SIR_model_OOS_C(  dataTable, sd, ed, N, weights, normalize, num_trials, train_days);
    rE  = fit_SEIR_model_OOS_C( dataTable, sd, ed, N, weights, normalize, num_trials, train_days);
    rES = fit_SEIRS_model_OOS_C(dataTable, sd, ed, N, weights, normalize, num_trials, train_days);

    if isempty(rS) || isempty(rE) || isempty(rES)
        fprintf('  -> Skipped (window too short for this split).\n');
        continue;
    end

    aicc_train = [rS.AICc_train, rE.AICc_train, rES.AICc_train];
    rmse_test  = [rS.RMSE_test_I, rE.RMSE_test_I, rES.RMSE_test_I];

    [~, i_aicc] = min(aicc_train);
    [~, i_rmse] = min(rmse_test);
    match = (i_aicc == i_rmse);

    fprintf('  AICc_train:  SIR %.1f | SEIR %.1f | SEIRS %.1f  -> best: %s\n', ...
        aicc_train(1), aicc_train(2), aicc_train(3), models{i_aicc});
    fprintf('  RMSE_test_I: SIR %.1f | SEIR %.1f | SEIRS %.1f  -> best: %s\n', ...
        rmse_test(1), rmse_test(2), rmse_test(3), models{i_rmse});
    fprintf('  Match: %s\n', mat2str(match));

    rows(end+1,:) = { i, datestr(sd,'yyyy-mm-dd'), datestr(ed,'yyyy-mm-dd'), ...
        rS.T, rS.T_train, rS.T_test, ...
        aicc_train(1), aicc_train(2), aicc_train(3), models{i_aicc}, ...
        rmse_test(1), rmse_test(2), rmse_test(3), models{i_rmse}, match }; %#ok<SAGROW>
end

%% ---- Results table ----
Tres = cell2table(rows, 'VariableNames', ...
    {'Window','StartDate','EndDate','T','T_train','T_test', ...
     'AICc_train_SIR','AICc_train_SEIR','AICc_train_SEIRS','BestByAICc_train', ...
     'RMSE_test_SIR','RMSE_test_SEIR','RMSE_test_SEIRS','BestByRMSE_test', ...
     'Match'});

outfile = fullfile(out_dir, sprintf('oos_validation_%s.csv', country));
writetable(Tres, outfile);

nWin   = height(Tres);
nMatch = sum(Tres.Match);

fprintf('\n===================================================\n');
fprintf('OOS validation complete: %s\n', country);
fprintf('AICc(train)-best == RMSE(test)-best: %d/%d windows (%.0f%%)\n', ...
    nMatch, nWin, 100*nMatch/nWin);
fprintf('===================================================\n');

fprintf('\nMean RMSE_test_I per model (lower = better generalization):\n');
fprintf('  SIR:   %.1f\n', mean(Tres.RMSE_test_SIR));
fprintf('  SEIR:  %.1f\n', mean(Tres.RMSE_test_SEIR));
fprintf('  SEIRS: %.1f\n', mean(Tres.RMSE_test_SEIRS));

fprintf('\nHow many times each model wins (RMSE_test_I):\n');
[~, winIdx] = min([Tres.RMSE_test_SIR, Tres.RMSE_test_SEIR, Tres.RMSE_test_SEIRS], [], 2);
for m = 1:3
    fprintf('  %-6s: %d/%d windows\n', models{m}, sum(winIdx==m), nWin);
end

fprintf('\nSaved -> %s\n', outfile);
