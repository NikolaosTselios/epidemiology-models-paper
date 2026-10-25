%% run_residual_diagnostics_C.m
% Τρέχει τα διαγνωστικά residuals και για τα τρία μοντέλα και παράγει
% συγκριτικό πίνακα ανά κριτήριο.
%
% Προϋποθέτει ότι τα results_SIR / results_SEIR / results_SEIRS υπάρχουν
% ήδη στο workspace (από το run_fitting_models_C) ή φορτώνονται από .mat.

%% === Διαγνωστικά ανά μοντέλο ===
plotW = [5 13 14 20];

diag_SIR   = compute_residual_diagnostics_C(results_SIR,   'PlotWindows', plotW);
diag_SEIR  = compute_residual_diagnostics_C(results_SEIR,  'PlotWindows', plotW);
diag_SEIRS = compute_residual_diagnostics_C(results_SEIRS, 'PlotWindows', plotW);

diagnostics = {diag_SIR, diag_SEIR, diag_SEIRS};
model_names = {'SIR','SEIR','SEIRS'};

%% === Συγκριτικά σχήματα ===
comparison = compare_residual_diagnostics_C(diagnostics, model_names, ...
    'SavePlots', true, 'ShowPlots', false);

%% === Πίνακας ανά κριτήριο (αυτό μπαίνει στο paper) ===
% ΔΕΝ αθροίζουμε ασύμβατες ποσότητες σε ενιαίο score. Κάθε κριτήριο
% αναφέρεται χωριστά· ο αναγνώστης βλέπει τα trade-offs.

nM = numel(diagnostics);
rows = cell(nM,1);
for m = 1:nM
    s = diagnostics{m}.summary;
    rows{m} = { model_names{m}, ...
        s.meanRMSE_I, s.meanRMSE_R, ...
        s.meanMAE_I,  s.meanMAE_R, ...
        s.relBias_I,  s.relBias_R, ...
        s.meanDW_I,   s.meanDW_R, ...
        s.normalityRate_I,    s.normalityRate_R, ...
        s.independenceRate_I, s.independenceRate_R };
end
crit = cell2table(vertcat(rows{:}), 'VariableNames', ...
    {'Model','RMSE_I','RMSE_R','MAE_I','MAE_R', ...
     'relBias_I','relBias_R','DW_I','DW_R', ...
     'NormRate_I','NormRate_R','IndepRate_I','IndepRate_R'});

fprintf('\n===== Residual diagnostics ανά κριτήριο =====\n');
disp(crit);
writetable(crit, 'residual_diagnostics_by_criterion_C_5_13_14_20.csv');

%% === Προαιρετικά: rank-based composite (κλιμακοανεξάρτητο) ===
% Αν χρειάζεσαι οπωσδήποτε ενιαία κατάταξη, ΜΗΝ αθροίζεις τις τιμές.
% ’θροισε τις ΘΕΣΕΙΣ (1 = καλύτερο) ώστε καμία κλίμακα να μην κυριαρχεί.

lowerBetter = {'RMSE_I','RMSE_R','MAE_I','MAE_R','relBias_I','relBias_R'};
higherBetter= {'NormRate_I','NormRate_R','IndepRate_I','IndepRate_R'};

R = zeros(nM,0);
for c = lowerBetter
    [~,ord] = sort(crit.(c{1}),'ascend');   [~,rk] = sort(ord);   R(:,end+1) = rk; %#ok<SAGROW>
end
for c = higherBetter
    [~,ord] = sort(crit.(c{1}),'descend');  [~,rk] = sort(ord);   R(:,end+1) = rk; %#ok<SAGROW>
end
% DW: κοντά στο 2 είναι το καλύτερο
for c = {'DW_I','DW_R'}
    [~,ord] = sort(abs(crit.(c{1})-2),'ascend'); [~,rk] = sort(ord); R(:,end+1) = rk; %#ok<SAGROW>
end

rankTable = table(model_names(:), sum(R,2), mean(R,2), ...
    'VariableNames', {'Model','RankSum','MeanRank'});
rankTable = sortrows(rankTable,'RankSum','ascend');

fprintf('\n===== Rank-based composite (μικρότερο = καλύτερο) =====\n');
disp(rankTable);
writetable(rankTable,'residual_ranking_rankbased_C_5_13_14_20.csv');
