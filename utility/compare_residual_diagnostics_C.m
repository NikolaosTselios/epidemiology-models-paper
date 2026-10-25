function comparison = compare_residual_diagnostics_C(diagnostics, model_names, varargin)
%% compare_residual_diagnostics_C   (συμβατό με MATLAB R2019a)
%
% Συγκριτική παρουσίαση διαγνωστικών residuals για SIR / SEIR / SEIRS.
%
% INPUT
%   diagnostics : cell array { sir_diag, seir_diag, seirs_diag }
%                 (έξοδοι του compute_residual_diagnostics_C)
%   model_names : cell array {'SIR','SEIR','SEIRS'}
%
% Optional
%   'SavePlots' (true), 'ShowPlots' (false), 'OutputDir'
%
% OUTPUT
%   comparison.models   : ανά μοντέλο
%   comparison.summary  : συγκεντρωτικά
%   comparison.criteria : πίνακας ανά κριτήριο (αυτός μπαίνει στο paper)
%   comparison.ranking  : RANK-BASED κατάταξη
%
% ΣΗΜΑΝΤΙΚΟ - ΓΙΑΤΙ Η ΚΑΤΑΤΑΞΗ ΕΙΝΑΙ RANK-BASED
% Το άθροισμα RMSE + MAE + Bias + DW - rates ΔΕΝ είναι έγκυρο: το RMSE_R
% είναι της τάξης 10^4-10^5 ενώ το DW και τα ποσοστά είναι της τάξης 1.
% Το σύνθετο score θα ήταν στην πράξη σκέτο RMSE_R και τα διαγνωστικά θα
% συνεισέφεραν ~0.001%. Εδώ αθροίζονται οι ΘΕΣΕΙΣ κατάταξης (1 = καλύτερο)
% ανά κριτήριο, ώστε καμία κλίμακα να μην κυριαρχεί.
%
% ΣΥΜΒΑΤΟΤΗΤΑ: μόνο subplot / print (όχι tiledlayout, nexttile,
% exportgraphics). Το boxplot έχει fallback χωρίς Statistics Toolbox.

%% ===================== Parse inputs =====================
p = inputParser;
addParameter(p,'SavePlots',true,@islogical);
addParameter(p,'ShowPlots',false,@islogical);
addParameter(p,'OutputDir',fullfile('plots','Residual_Comparison_free_g_a_d_5_13_14_20'), ...
    @(x)ischar(x)||isstring(x));
parse(p,varargin{:});
opt = p.Results;

assert(iscell(diagnostics),'diagnostics must be a cell array.');
assert(iscell(model_names),'model_names must be a cell array.');
assert(numel(diagnostics)==numel(model_names), ...
    'diagnostics and model_names must have same length.');

nModels = numel(diagnostics);

fprintf('\n=========================================================\n');
fprintf('Residual Diagnostics Comparison  (%d models)\n', nModels);
fprintf('=========================================================\n');

if opt.SavePlots && ~exist(opt.OutputDir,'dir'), mkdir(opt.OutputDir); end

%% ===================== Extract =====================
comparison = struct();
comparison.models = struct([]);

for m = 1:nModels
    d = diagnostics{m};
    comparison.models(m).name     = model_names{m};
    comparison.models(m).table    = d.table;
    comparison.models(m).summary  = d.summary;
    comparison.models(m).nWindows = height(d.table);
end

nW = height(diagnostics{1}.table);
Z  = @() zeros(nW,nModels);

RMSE_I=Z(); RMSE_R=Z(); MAE_I=Z(); MAE_R=Z();
DW_I=Z();   DW_R=Z();   BIAS_I=Z(); BIAS_R=Z();
JB_I=Z();   JB_R=Z();   LB_I=Z();   LB_R=Z();
SKEW_I=Z(); SKEW_R=Z(); KURT_I=Z(); KURT_R=Z();
RELB_I=Z(); RELB_R=Z();

for m = 1:nModels
    t = diagnostics{m}.table;
    RMSE_I(:,m)=t.RMSE_I;   RMSE_R(:,m)=t.RMSE_R;
    MAE_I(:,m) =t.MAE_I;    MAE_R(:,m) =t.MAE_R;
    DW_I(:,m)  =t.DW_I;     DW_R(:,m)  =t.DW_R;
    BIAS_I(:,m)=abs(t.Mean_I);  BIAS_R(:,m)=abs(t.Mean_R);
    JB_I(:,m)  =t.JB_p_I;   JB_R(:,m)  =t.JB_p_R;
    LB_I(:,m)  =t.LB_p_I;   LB_R(:,m)  =t.LB_p_R;
    SKEW_I(:,m)=t.Skew_I;   SKEW_R(:,m)=t.Skew_R;
    KURT_I(:,m)=t.Kurtosis_I; KURT_R(:,m)=t.Kurtosis_R;
    % σχετική μεροληψία: κλιμακοανεξάρτητη, συγκρίσιμη I vs R
    RELB_I(:,m)=abs(t.Mean_I)./max(t.RMSE_I,eps);
    RELB_R(:,m)=abs(t.Mean_R)./max(t.RMSE_R,eps);
end

labels = cell(1,2*nModels);
for i=1:nModels, labels{i}          = sprintf('%s-I',model_names{i}); end
for i=1:nModels, labels{nModels+i}  = sprintf('%s-R',model_names{i}); end

%% ===================== Summary =====================
s = struct();
s.Model              = model_names(:);
s.meanRMSE_I         = mean(RMSE_I)';   s.meanRMSE_R        = mean(RMSE_R)';
s.meanMAE_I          = mean(MAE_I)';    s.meanMAE_R         = mean(MAE_R)';
s.meanDW_I           = mean(DW_I)';     s.meanDW_R          = mean(DW_R)';
s.meanRelBias_I      = mean(RELB_I)';   s.meanRelBias_R     = mean(RELB_R)';
s.meanSkew_I         = mean(SKEW_I)';   s.meanSkew_R        = mean(SKEW_R)';
s.meanKurtosis_I     = mean(KURT_I)';   s.meanKurtosis_R    = mean(KURT_R)';
s.normalityRate_I    = 100*sum(JB_I>0.05)'/nW;
s.normalityRate_R    = 100*sum(JB_R>0.05)'/nW;
s.independenceRate_I = 100*sum(LB_I>0.05)'/nW;
s.independenceRate_R = 100*sum(LB_R>0.05)'/nW;
comparison.summary = s;

%% ===================== Πίνακας ανά κριτήριο =====================
criteria = table(model_names(:), ...
    s.meanRMSE_I, s.meanRMSE_R, s.meanMAE_I, s.meanMAE_R, ...
    s.meanRelBias_I, s.meanRelBias_R, s.meanDW_I, s.meanDW_R, ...
    s.normalityRate_I, s.normalityRate_R, ...
    s.independenceRate_I, s.independenceRate_R, ...
    'VariableNames',{'Model','RMSE_I','RMSE_R','MAE_I','MAE_R', ...
    'relBias_I','relBias_R','DW_I','DW_R', ...
    'NormRate_I','NormRate_R','IndepRate_I','IndepRate_R'});

comparison.criteria = criteria;
fprintf('\n--- Diagnostics by criterion ---\n');
disp(criteria);

%% ===================== Σχήμα 1: boxplots =====================
fprintf('\nGenerating comparison figures...\n');
fig = local_newfig(opt,[50 50 1750 980]);

subplot(3,4,1);  local_box(RMSE_I, model_names); ylabel('RMSE'); title('RMSE (Infected)');  grid on
subplot(3,4,2);  local_box(RMSE_R, model_names); ylabel('RMSE'); title('RMSE (Removed)');   grid on
subplot(3,4,3);  local_box(MAE_I,  model_names); ylabel('MAE');  title('MAE (Infected)');   grid on
subplot(3,4,4);  local_box(MAE_R,  model_names); ylabel('MAE');  title('MAE (Removed)');    grid on

subplot(3,4,5);  local_box([DW_I DW_R], labels); hold on
                 local_href(2,'r--'); ylabel('DW'); title('Durbin-Watson (ideal 2)');
                 xtickangle(45); grid on

subplot(3,4,6);  local_box([RELB_I RELB_R], labels);
                 ylabel('|mean| / RMSE'); title('Relative residual bias');
                 xtickangle(45); grid on

subplot(3,4,7);  local_box([JB_I JB_R], labels); hold on
                 local_href(0.05,'r--'); ylabel('p-value'); title('Jarque-Bera');
                 xtickangle(45); grid on

subplot(3,4,8);  local_box([LB_I LB_R], labels); hold on
                 local_href(0.05,'r--'); ylabel('p-value'); title('Ljung-Box');
                 xtickangle(45); grid on

subplot(3,4,9);  local_box([SKEW_I SKEW_R], labels); hold on
                 local_href(0,'r--'); ylabel('Skewness'); title('Residual skewness');
                 xtickangle(45); grid on

subplot(3,4,10); local_box([KURT_I KURT_R], labels); hold on
                 local_href(3,'r--'); ylabel('Kurtosis'); title('Residual kurtosis');
                 xtickangle(45); grid on

subplot(3,4,[11 12]);
Y = [100*sum(JB_I>0.05)/nW, 100*sum(JB_R>0.05)/nW;
     100*sum(LB_I>0.05)/nW, 100*sum(LB_R>0.05)/nW]';
bar(Y);
set(gca,'XTick',1:numel(labels),'XTickLabel',labels);
xtickangle(45); ylabel('Passing windows (%)'); ylim([0 100]); grid on
legend({'Normality (JB)','Independence (LB)'},'Location','northoutside');
title('Diagnostic success rate');

local_sgtitle(fig,'Comparison of Residual Diagnostics');
if opt.SavePlots
    local_save(fig, fullfile(opt.OutputDir,'Residual_Diagnostics_Comparison'));
end
if ~opt.ShowPlots, close(fig); end

%% ===================== Σχήμα 2: heatmaps =====================
fprintf('Generating heatmaps...\n');
fig = local_newfig(opt,[100 100 1700 900]);

subplot(2,3,1); imagesc(RMSE_I); colorbar
set(gca,'XTick',1:nModels,'XTickLabel',model_names);
xlabel('Model'); ylabel('Window'); title('RMSE (Infected)');

subplot(2,3,2); imagesc(RMSE_R); colorbar
set(gca,'XTick',1:nModels,'XTickLabel',model_names);
xlabel('Model'); ylabel('Window'); title('RMSE (Removed)');

subplot(2,3,3); imagesc([DW_I DW_R]); colorbar
set(gca,'XTick',1:numel(labels),'XTickLabel',labels); xtickangle(45);
xlabel('Model'); ylabel('Window'); title('Durbin-Watson');

subplot(2,3,4); imagesc([JB_I JB_R]); colorbar; caxis([0 1]);
set(gca,'XTick',1:numel(labels),'XTickLabel',labels); xtickangle(45);
xlabel('Model'); ylabel('Window'); title('Jarque-Bera p-values');

subplot(2,3,5); imagesc([LB_I LB_R]); colorbar; caxis([0 1]);
set(gca,'XTick',1:numel(labels),'XTickLabel',labels); xtickangle(45);
xlabel('Model'); ylabel('Window'); title('Ljung-Box p-values');

subplot(2,3,6); imagesc([RELB_I RELB_R]); colorbar
set(gca,'XTick',1:numel(labels),'XTickLabel',labels); xtickangle(45);
xlabel('Model'); ylabel('Window'); title('Relative residual bias');

local_sgtitle(fig,'Residual Diagnostics by Fitting Window');
if opt.SavePlots
    local_save(fig, fullfile(opt.OutputDir,'Residual_Diagnostics_Heatmaps'));
end
if ~opt.ShowPlots, close(fig); end

%% ===================== Rank-based κατάταξη =====================
fprintf('Computing rank-based ranking...\n');

lowerBetter  = {'RMSE_I','RMSE_R','MAE_I','MAE_R','relBias_I','relBias_R'};
higherBetter = {'NormRate_I','NormRate_R','IndepRate_I','IndepRate_R'};
dwCols       = {'DW_I','DW_R'};   % το ιδανικό είναι το 2, όχι το ελάχιστο

Rk = [];
for c = lowerBetter
    Rk(:,end+1) = local_rank(criteria.(c{1}), 'ascend'); %#ok<AGROW>
end
for c = higherBetter
    Rk(:,end+1) = local_rank(criteria.(c{1}), 'descend'); %#ok<AGROW>
end
for c = dwCols
    Rk(:,end+1) = local_rank(abs(criteria.(c{1})-2), 'ascend'); %#ok<AGROW>
end

ranking = table(model_names(:), sum(Rk,2), mean(Rk,2), ...
    'VariableNames',{'Model','RankSum','MeanRank'});
ranking = sortrows(ranking,'RankSum','ascend');
comparison.ranking = ranking;

fprintf('\n--- Rank-based ranking (lower = better) ---\n');
disp(ranking);

fig = local_newfig(opt,[200 200 900 460]);
bar(ranking.RankSum);
set(gca,'XTick',1:height(ranking),'XTickLabel',ranking.Model);
ylabel('Sum of criterion ranks'); title('Overall residual ranking'); grid on
for k = 1:height(ranking)
    text(k, ranking.RankSum(k), sprintf('#%d',k), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom', ...
        'FontWeight','bold');
end
if opt.SavePlots
    local_save(fig, fullfile(opt.OutputDir,'Residual_Ranking'));
end
if ~opt.ShowPlots, close(fig); end

%% ===================== Αποθήκευση πινάκων =====================
if opt.SavePlots
    writetable(criteria, fullfile(opt.OutputDir,'diagnostics_by_criterion_5_13_14_20.csv'));
    writetable(ranking,  fullfile(opt.OutputDir,'residual_ranking_rankbased_5_13_14_20.csv'));
end

fprintf('\nComparison completed.\n');
end % ============================ main ============================


%% ======================================================================
%% Βοηθητικές συναρτήσεις
%% ======================================================================

function fig = local_newfig(opt, pos)
vis = 'off'; if opt.ShowPlots, vis = 'on'; end
fig = figure('Visible',vis,'Color','w','Position',pos);
end


function local_href(y, style)
% οριζόντια γραμμή αναφοράς - χωρίς yline (ασφαλέστερο σε παλιές εκδόσεις)
xl = xlim;
plot(xl, [y y], style);
xlim(xl);
end


function local_box(M, labs)
% boxplot με fallback αν λείπει το Statistics Toolbox
if exist('boxplot','file') == 2
    try
        boxplot(M, 'Labels', labs);
        return
    catch
    end
end
local_boxfallback(M, labs);
end


function local_boxfallback(M, labs)
% Απλό boxplot από βασικά γραφικά: κουτί Q1-Q3, διάμεσος, whiskers 1.5*IQR
n = size(M,2);
hold on
for j = 1:n
    x = M(:,j); x = x(~isnan(x));
    if isempty(x), continue; end
    q  = local_quantile(x,[0.25 0.5 0.75]);
    iqr_ = q(3)-q(1);
    lo = min(x(x >= q(1)-1.5*iqr_));
    hi = max(x(x <= q(3)+1.5*iqr_));
    w  = 0.3;
    % κουτί
    patch([j-w j+w j+w j-w],[q(1) q(1) q(3) q(3)], [0.85 0.9 0.97], ...
          'EdgeColor',[0.2 0.2 0.2]);
    % διάμεσος
    plot([j-w j+w],[q(2) q(2)],'r-','LineWidth',1.6);
    % whiskers
    plot([j j],[q(3) hi],'k-'); plot([j j],[lo q(1)],'k-');
    plot([j-w/2 j+w/2],[hi hi],'k-'); plot([j-w/2 j+w/2],[lo lo],'k-');
    % outliers
    out = x(x<lo | x>hi);
    if ~isempty(out), plot(j*ones(size(out)), out, 'k+'); end
end
xlim([0.5 n+0.5]);
set(gca,'XTick',1:n,'XTickLabel',labs);
hold off
end


function q = local_quantile(x, pr)
% quantile χωρίς Statistics Toolbox (γραμμική παρεμβολή, τύπος MATLAB)
x = sort(x(:));
n = numel(x);
if n == 1, q = repmat(x,size(pr)); return; end
pos = n*pr(:) + 0.5;
pos = min(max(pos,1),n);
lo  = floor(pos); hi = ceil(pos); fr = pos - lo;
q   = x(lo).*(1-fr) + x(hi).*fr;
q   = reshape(q, size(pr));
end


function rk = local_rank(v, direction)
% θέσεις κατάταξης: 1 = καλύτερο
[~,ord] = sort(v(:), direction);
rk = zeros(numel(v),1);
rk(ord) = 1:numel(v);
end


function local_sgtitle(fig, txt)
if exist('sgtitle','file') || exist('sgtitle','builtin')
    try
        sgtitle(txt,'FontWeight','bold','FontSize',15);
        return
    catch
    end
end
annotation(fig,'textbox',[0 0.95 1 0.05],'String',txt, ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontWeight','bold','FontSize',15);
end

function local_save(fig, basepath)

if exist('exportgraphics','file')
    exportgraphics(fig,[basepath '.png'],'Resolution',300);
    exportgraphics(fig,[basepath '.pdf'],'ContentType','vector');
else
    set(fig,'PaperPositionMode','auto','InvertHardcopy','off','Color','w');
    print(fig, [basepath '.png'], '-dpng', '-r300');
    print(fig, [basepath '.pdf'], '-dpdf', '-painters', '-bestfit');
end
end
