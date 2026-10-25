function plot_model_comparisons1()
%% plot_model_comparisons.m
% Δημιουργεί συγκριτικά γραφήματα για SIR, SEIR και SEIRS μοντέλα ανά περίοδο.

clear; close all; clc;

% === Όνομα χώρας ===
countryName = 'Greece';

%% === Φόρτωση αποτελεσμάτων ===
load(fullfile('results', 'Results_SIR', countryName, 'SIR_all_periods_results1_C.mat'));    
all_results_SIR = all_results;

load(fullfile('results', 'Results_SEIR', countryName, 'SEIR_all_periods_results1_C.mat'));  
all_results_SEIR = all_results;

load(fullfile('results', 'Results_SEIRS', countryName, 'SEIRS_all_periods_results1_C.mat'));
all_results_SEIRS = all_results;

%% === Έλεγχος συμβατότητας ===
nPeriods = length(all_results_SIR);
assert(length(all_results_SEIR) == nPeriods, 'Mismatch SIR vs SEIR');
assert(length(all_results_SEIRS) == nPeriods, 'Mismatch SIR vs SEIRS');

%% === Δημιουργία φακέλου για plots ===
plot_dir = fullfile('plots', 'Compare_Models_fixed_axis_10_prettier_plots_10_days_GR_07-9-2026', countryName);
if ~exist(plot_dir, 'dir')
    mkdir(plot_dir);
end

%% === Loop εμφάνισης και αποθήκευσης για κάθε περίοδο ===
for i = 1:nPeriods
    fprintf('\n=== Δημιουργία plot για περίοδο %d/%d ===\n', i, nPeriods);

    res_SIR   = all_results_SIR{i};
    res_SEIR  = all_results_SEIR{i};
    res_SEIRS = all_results_SEIRS{i};

    %% === Global axis limits for fair comparison ===

infected_all = [
    res_SIR.infected_real(:);
    res_SIR.infected_model(:);
    res_SEIR.infected_model(:);
    res_SEIRS.infected_model(:)
];

removed_all = [
    res_SIR.removed_real(:);
    res_SIR.removed_model(:);
    res_SEIR.removed_model(:);
    res_SEIRS.removed_model(:)
];

infected_min = min(infected_all);
infected_max = max(infected_all);

removed_min = min(removed_all);
removed_max = max(removed_all);

%% === Δημιουργία Figure ===

figure('Color','w', ...
       'Name', sprintf('Comparison - Period %d', i), ...
       'Position', [100 100 1700 980]);

%% =========================
%% Infected Subplots
%% =========================

subplot(2,3,1);

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',11,...
    'LineWidth',1);

plot(res_SIR.dates,res_SIR.infected_real,...
    'Color',[0 0 0],...
    'LineWidth',2.8);
hold on;

plot(res_SIR.dates, res_SIR.infected_model, ...
    '--', ...
    'Color', [0.85 0 0], ...
    'LineWidth', 2);

title('SIR',...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'FontWeight','bold');

ylabel('Infected',...
    'FontSize',13,...
    'FontName','Times New Roman')

ylim([infected_min infected_max]);

grid on
box off

ax = gca;
ax.GridAlpha = 0.08;

%% -------------------------

subplot(2,3,2);

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',11,...
    'LineWidth',1);

plot(res_SEIR.dates, res_SEIR.infected_real, ...
    'Color', [0 0 0], ...
    'LineWidth', 2.8);
hold on;

plot(res_SEIR.dates, res_SEIR.infected_model, ...
    '--', ...
    'Color', [0.85 0 0], ...
    'LineWidth', 2);

title('SEIR',...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'FontWeight','bold');

ylim([infected_min infected_max]);

grid on
box off

ax = gca;
ax.GridAlpha = 0.08;

%% -------------------------

subplot(2,3,3);

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',11,...
    'LineWidth',1);

plot(res_SEIRS.dates, res_SEIRS.infected_real, ...
    'Color', [0 0 0], ...
    'LineWidth', 2.8);
hold on;

plot(res_SEIRS.dates, res_SEIRS.infected_model, ...
    '--', ...
    'Color', [0.85 0 0], ...
    'LineWidth', 2);

title('SEIRS',...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'FontWeight','bold');

ylim([infected_min infected_max]);

grid on
box off

ax = gca;
ax.GridAlpha = 0.08;

%% =========================
%% Removed Subplots
%% =========================

subplot(2,3,4);

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',11,...
    'LineWidth',1);

plot(res_SIR.dates, res_SIR.removed_real, ...
    'Color', [0 0 0], ...
    'LineWidth', 2.8);
hold on;

plot(res_SIR.dates, res_SIR.removed_model, ...
    '--', ...
    'Color', [0 0.35 0.9], ...
    'LineWidth', 2);

title('SIR',...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'FontWeight','bold');

ylabel('Removed',...
    'FontSize',13,...
    'FontName','Times New Roman')
xlabel('Date');

ylim([removed_min removed_max]);

grid on
box off

ax = gca;
ax.GridAlpha = 0.08;

%% -------------------------

subplot(2,3,5);

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',11,...
    'LineWidth',1);

plot(res_SEIR.dates, res_SEIR.removed_real, ...
    'Color', [0 0 0], ...
    'LineWidth', 2.8);
hold on;

plot(res_SEIR.dates, res_SEIR.removed_model, ...
    '--', ...
    'Color', [0 0.35 0.9], ...
    'LineWidth', 2);

title('SEIR',...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'FontWeight','bold');

xlabel('Date');

ylim([removed_min removed_max]);

grid on
box off

ax = gca;
ax.GridAlpha = 0.08;
%% -------------------------

subplot(2,3,6);

set(gca,...
    'FontName','Times New Roman',...
    'FontSize',11,...
    'LineWidth',1);

plot(res_SEIRS.dates, res_SEIRS.removed_real, ...
    'Color', [0 0 0], ...
    'LineWidth', 2.8);
hold on;

plot(res_SEIRS.dates, res_SEIRS.removed_model, ...
    '--', ...
    'Color', [0 0.35 0.9], ...
    'LineWidth', 2);

title('SEIRS',...
    'FontName','Times New Roman',...
    'FontSize',15,...
    'FontWeight','bold');

xlabel('Date');

ylim([removed_min removed_max]);

grid on
box off

ax = gca;
ax.GridAlpha = 0.08;
%% === Common title ===

sgtitle('Model Comparison',...
    'FontName','Times New Roman',...
    'FontSize',17,...
    'FontWeight','bold');

    % === Αποθήκευση ===
    saveas(gcf, fullfile(plot_dir, sprintf('Compare_Models_Period_%02d.png', i)));
    close(gcf);
    
    fprintf('SIR removed max: %.2f\n', max(res_SIR.removed_real));
    fprintf('SEIR removed max: %.2f\n', max(res_SEIR.removed_real));
    fprintf('SEIRS removed max: %.2f\n', max(res_SEIRS.removed_real));
end

%% === Plot R0 over time ===
fprintf('\n=== Δημιουργία R0 over time ===\n');

R0_SIR   = zeros(nPeriods,1);
R0_SEIR  = zeros(nPeriods,1);
R0_SEIRS = zeros(nPeriods,1);

for i = 1:nPeriods
    R0_SIR(i)   = all_results_SIR{i}.R0;
    R0_SEIR(i)  = all_results_SEIR{i}.R0;
    R0_SEIRS(i) = all_results_SEIRS{i}.R0;
end

figure('Name','R0 over time','Position',[100,100,800,600]);
plot(1:nPeriods, R0_SIR, 'r-o', 'LineWidth', 2); hold on;
plot(1:nPeriods, R0_SEIR, 'g-o', 'LineWidth', 2);
plot(1:nPeriods, R0_SEIRS, 'b-o', 'LineWidth', 2);

grid on;
xlabel('Period');
ylabel('R0');
legend({'SIR','SEIR','SEIRS'}, 'Location','best');
title('R0 over time per model');

saveas(gcf, fullfile(plot_dir, 'Compare_R0_over_time_FN.png'));
close(gcf);

%% === Δημιουργία πίνακα ===
Period = (1:nPeriods)';
R0_table = table(Period, R0_SIR, R0_SEIR, R0_SEIRS);

writetable(R0_table, fullfile(plot_dir, 'Summary_R0_Table_FN.csv'));

fprintf('\n=== Ολοκληρώθηκε η δημιουργία όλων των γραφημάτων και του πίνακα R0 ===\n');

end
