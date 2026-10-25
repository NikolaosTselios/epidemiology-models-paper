function T = weight_sensitivity(csv_file, country)
% WEIGHT_SENSITIVITY
% Re-fits SIR/SEIR/SEIRS under several I/R weightings and reports how the
% AICc ranking changes. Confirms whether the model-selection conclusion is
% robust to the (subjective) 0.9/0.1 choice, or driven by it.
%
% Usage:
%   weight_sensitivity('Covid_19_Greece.csv', 'Greece')
%   weight_sensitivity('Covid_19_Finland.csv',  'Finland')
%
% NOTE: this re-runs the fitting, so it OVERWRITES results/Results_*/<country>/*.mat
% with reduced-trial versions. Your canonical deliverables are already saved in
% results_FINAL/<timestamp>/, so nothing important is lost; but if you want the
% canonical .mat back, re-run the master afterwards.

N         = 5.6e6;
normalize = true;
num_trials = 800;      % reduced for speed; rankings are stable well below 800
step_days = 30;

% weight settings to test: [w_I  w_R]
Wset = {[0.9 0.1], [0.7 0.3], [0.5 0.5], [1 0]};
models = {'SIR','SEIR','SEIRS'};

rows = cell(numel(Wset), 9);
fprintf('\n===== Weight sensitivity: %s =====\n', country);

for w = 1:numel(Wset)
    W = Wset{w};
    fprintf('\n--- weights [%.2f %.2f] ---\n', W(1), W(2));

    rS  = compare_SIR_periods_C(  csv_file, country, N, W, normalize, num_trials, step_days);
    rE  = compare_SEIR_periods_C( csv_file, country, N, W, normalize, num_trials, step_days);
    rES = compare_SEIRS_periods_C(csv_file, country, N, W, normalize, num_trials, step_days);

    aicc = [rS.AICc_tot, rE.AICc_tot, rES.AICc_tot];   % [SIR SEIR SEIRS]
    d    = aicc - min(aicc);
    wk   = exp(-0.5*d);  wk = wk / sum(wk);            % Akaike weights
    [~, ord] = sort(aicc, 'ascend');
    ranking  = strjoin(models(ord), ' < ');

    rows(w,:) = { mat2str(W), aicc(1), aicc(2), aicc(3), ...
                  ranking, wk(1), wk(2), wk(3), models{ord(1)} };

    fprintf('AICc: SIR %.1f | SEIR %.1f | SEIRS %.1f\n', aicc(1), aicc(2), aicc(3));
    fprintf('Akaike w: SIR %.3g | SEIR %.3g | SEIRS %.3g\n', wk(1), wk(2), wk(3));
    fprintf('Ranking: %s   (winner: %s)\n', ranking, models{ord(1)});
end

T = cell2table(rows, 'VariableNames', ...
    {'Weights','AICc_SIR','AICc_SEIR','AICc_SEIRS', ...
     'Ranking','w_SIR','w_SEIR','w_SEIRS','Winner'});

outfile = sprintf('weight_sensitivity_7_9_%s.csv', country);
writetable(T, outfile);

fprintf('\n===== Summary (%s) =====\n', country);
disp(T(:, {'Weights','Ranking','w_SIR','w_SEIR','w_SEIRS'}));
fprintf('Saved -> %s\n', outfile);

% the key question, answered in one line:
allWin = string(T.Winner);
if all(allWin == "SEIRS")
    fprintf('==> SEIRS wins under EVERY weighting: conclusion robust to the weight choice.\n');
else
    fprintf('==> Winner changes across weightings: report this — the conclusion is weight-sensitive.\n');
end
end