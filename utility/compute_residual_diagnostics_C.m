function diag_out = compute_residual_diagnostics_C(summary, varargin)
%% compute_residual_diagnostics_C   (συμβατό με MATLAB R2019a)
%
% Υπολογίζει διαγνωστικά residuals ανά χρονικό παράθυρο για ένα μοντέλο
% (SIR / SEIR / SEIRS) και επιστρέφει δομή συμβατή με το
% compare_residual_diagnostics_C.
%
% INPUT
%   summary : η δομή που επιστρέφει το compare_<MODEL>_periods_C
%             (χρειάζεται το πεδίο summary.all_results)
%
% Optional name-value
%   'Lags'        : lags για Ljung-Box (default: min(10, floor(T/3)))
%   'AdjustDoF'   : true (default) -> dof = lags - p, όπως αρμόζει σε
%                   residuals από προσαρμοσμένο μοντέλο
%   'PlotWindows' : δείκτες παραθύρων για αναλυτικά σχήματα. [] = κανένα
%   'SavePlots'   : true (default)
%   'ShowPlots'   : false (default)
%   'OutputDir'   : φάκελος αποθήκευσης σχημάτων
%
% OUTPUT
%   diag_out.table   : πίνακας ανά παράθυρο
%   diag_out.summary : συγκεντρωτικά ανά μοντέλο
%   diag_out.model   : όνομα μοντέλου
%
% ΣΗΜΕΙΩΣΗ ΕΡΜΗΝΕΙΑΣ
% Το fitting γίνεται σε εξομαλυμένα (7-day moving average) δεδομένα, άρα
% τα residuals είναι εκ κατασκευής αυτοσυσχετισμένα και το Ljung-Box
% αναμένεται να απορρίπτει την ανεξαρτησία στα περισσότερα παράθυρα.
% Αυτό υποστηρίζει την ερμηνεία του AIC/AICc ως ΣΧΕΤΙΚΟΥ δείκτη επιλογής.
%
% Το Jarque-Bera είναι ασυμπτωτικό: με T ~ 30 η p-value είναι ενδεικτική.
%
% ΣΥΜΒΑΤΟΤΗΤΑ: χρησιμοποιεί μόνο subplot / print (όχι tiledlayout,
% nexttile, exportgraphics) ώστε να τρέχει σε R2019a και παλαιότερες.

%% ===================== Parse inputs =====================
p_in = inputParser;
addParameter(p_in,'Lags',      [],    @(x)isnumeric(x));
addParameter(p_in,'AdjustDoF', true,  @islogical);
addParameter(p_in,'PlotWindows',[],   @(x)isnumeric(x));
addParameter(p_in,'SavePlots', true,  @islogical);
addParameter(p_in,'ShowPlots', false, @islogical);
addParameter(p_in,'OutputDir', fullfile('plots','Residual_Diagnostics_free_g_a_d_5_13_14_20'), ...
                               @(x)ischar(x)||isstring(x));
parse(p_in,varargin{:});
opt = p_in.Results;

assert(isfield(summary,'all_results'), ...
    'Το summary πρέπει να περιέχει το πεδίο all_results.');

ar = summary.all_results;
ar = ar(~cellfun(@isempty, ar));          % κράτα μόνο τα επιτυχή παράθυρα
W  = numel(ar);
assert(W > 0, 'Δεν βρέθηκαν αποτελέσματα παραθύρων.');

model_name = ar{1}.model;

fprintf('\n=== Residual diagnostics: %s (%d windows) ===\n', model_name, W);

if opt.SavePlots && ~isempty(opt.PlotWindows) && ~exist(opt.OutputDir,'dir')
    mkdir(opt.OutputDir);
end

%% ===================== Preallocate =====================
Period   = strings(W,1);
Tw       = zeros(W,1);
p_model  = zeros(W,1);

[RSS_I,RSS_R,RMSE_I,RMSE_R,MAE_I,MAE_R]         = deal(zeros(W,1));
[Mean_I,Mean_R,Std_I,Std_R]                      = deal(zeros(W,1));
[DW_I,DW_R,Skew_I,Skew_R,Kurtosis_I,Kurtosis_R]  = deal(zeros(W,1));
[JB_p_I,JB_p_R,LB_p_I,LB_p_R,ACF1_I,ACF1_R]      = deal(zeros(W,1));

%% ===================== Loop παραθύρων =====================
for w = 1:W

    r = ar{w};

    % --- residuals (raw: χωρίς βάρη, χωρίς normalization) ---
    eI = r.infected_real(:) - r.infected_model(:);
    eR = r.removed_real(:)  - r.removed_model(:);

    T  = numel(eI);
    Tw(w)      = T;
    p_model(w) = r.p;

    d = r.dates;
    Period(w) = sprintf('%s to %s', ...
        datestr(d(1),'yyyy-mm-dd'), datestr(d(end),'yyyy-mm-dd'));

    % --- lags για Ljung-Box ---
    if isempty(opt.Lags)
        L = min(10, max(1, floor(T/3)));
    else
        L = min(opt.Lags, T-2);
    end

    % --- dof correction: residuals από προσαρμοσμένο μοντέλο ---
    if opt.AdjustDoF
        dofI = max(1, L - r.p);
        dofR = max(1, L - r.p);
    else
        dofI = L;  dofR = L;
    end

    % --- βασικές μετρικές ---
    RSS_I(w)  = sum(eI.^2);        RSS_R(w)  = sum(eR.^2);
    RMSE_I(w) = sqrt(mean(eI.^2)); RMSE_R(w) = sqrt(mean(eR.^2));
    MAE_I(w)  = mean(abs(eI));     MAE_R(w)  = mean(abs(eR));
    Mean_I(w) = mean(eI);          Mean_R(w) = mean(eR);
    Std_I(w)  = std(eI);           Std_R(w)  = std(eR);

    % --- διαγνωστικά ---
    DW_I(w) = local_dw(eI);        DW_R(w) = local_dw(eR);

    [Skew_I(w), Kurtosis_I(w)] = local_moments(eI);
    [Skew_R(w), Kurtosis_R(w)] = local_moments(eR);

    JB_p_I(w) = local_jb(eI);      JB_p_R(w) = local_jb(eR);

    [LB_p_I(w), ACF1_I(w)] = local_lb(eI, L, dofI);
    [LB_p_R(w), ACF1_R(w)] = local_lb(eR, L, dofR);

    % --- αναλυτικά σχήματα για επιλεγμένα παράθυρα ---
    if ismember(w, opt.PlotWindows)
        local_window_figure(eI, eR, r, model_name, w, Period(w), opt);
    end
end

%% ===================== Πίνακας ανά παράθυρο =====================
diag_out = struct();
diag_out.model = model_name;

diag_out.table = table(Period, Tw, p_model, ...
    RSS_I, RSS_R, RMSE_I, RMSE_R, MAE_I, MAE_R, ...
    Mean_I, Mean_R, Std_I, Std_R, ...
    DW_I, DW_R, ACF1_I, ACF1_R, ...
    JB_p_I, JB_p_R, LB_p_I, LB_p_R, ...
    Skew_I, Skew_R, Kurtosis_I, Kurtosis_R, ...
    'VariableNames', {'Period','T','p', ...
        'RSS_I','RSS_R','RMSE_I','RMSE_R','MAE_I','MAE_R', ...
        'Mean_I','Mean_R','Std_I','Std_R', ...
        'DW_I','DW_R','ACF1_I','ACF1_R', ...
        'JB_p_I','JB_p_R','LB_p_I','LB_p_R', ...
        'Skew_I','Skew_R','Kurtosis_I','Kurtosis_R'});

%% ===================== Συγκεντρωτικά =====================
s = struct();
s.model            = model_name;
s.nWindows         = W;

s.meanRMSE_I       = mean(RMSE_I);       s.meanRMSE_R       = mean(RMSE_R);
s.meanMAE_I        = mean(MAE_I);        s.meanMAE_R        = mean(MAE_R);
s.meanDW_I         = mean(DW_I);         s.meanDW_R         = mean(DW_R);
s.meanAbsBias_I    = mean(abs(Mean_I));  s.meanAbsBias_R    = mean(abs(Mean_R));
s.meanSkew_I       = mean(Skew_I);       s.meanSkew_R       = mean(Skew_R);
s.meanKurtosis_I   = mean(Kurtosis_I);   s.meanKurtosis_R   = mean(Kurtosis_R);
s.meanACF1_I       = mean(ACF1_I);       s.meanACF1_R       = mean(ACF1_R);

% πλήθος παραθύρων που ΔΕΝ απορρίπτουν την υπόθεση (α = 0.05)
s.normalityPass_I    = sum(JB_p_I > 0.05);   s.normalityPass_R    = sum(JB_p_R > 0.05);
s.independencePass_I = sum(LB_p_I > 0.05);   s.independencePass_R = sum(LB_p_R > 0.05);
s.normalityRate_I    = 100*s.normalityPass_I/W;
s.normalityRate_R    = 100*s.normalityPass_R/W;
s.independenceRate_I = 100*s.independencePass_I/W;
s.independenceRate_R = 100*s.independencePass_R/W;

% σχετική μεροληψία: |mean residual| / RMSE  (κλιμακοανεξάρτητη)
s.relBias_I = mean(abs(Mean_I)./max(RMSE_I,eps));
s.relBias_R = mean(abs(Mean_R)./max(RMSE_R,eps));

diag_out.summary = s;

%% ===================== Εκτύπωση =====================
fprintf('  Normality   (JB p>0.05): I %2d/%2d (%3.0f%%)  |  R %2d/%2d (%3.0f%%)\n', ...
    s.normalityPass_I, W, s.normalityRate_I, s.normalityPass_R, W, s.normalityRate_R);
fprintf('  Independence(LB p>0.05): I %2d/%2d (%3.0f%%)  |  R %2d/%2d (%3.0f%%)\n', ...
    s.independencePass_I, W, s.independenceRate_I, s.independencePass_R, W, s.independenceRate_R);
fprintf('  Durbin-Watson (ideal 2): I %.2f  |  R %.2f\n', s.meanDW_I, s.meanDW_R);
fprintf('  Lag-1 ACF     (ideal 0): I %.2f  |  R %.2f\n', s.meanACF1_I, s.meanACF1_R);
fprintf('  Relative bias |mean|/RMSE: I %.2f  |  R %.2f\n', s.relBias_I, s.relBias_R);

if s.independenceRate_I < 50 || s.independenceRate_R < 50
    fprintf(['  [EXPECTED] Low independence rate - consequence of the 7-day\n' ...
             '  smoothing. Report explicitly; AIC/AICc is interpreted as a\n' ...
             '  relative selection index.\n']);
end

end % ============================ main ============================


%% ======================================================================
%% Βοηθητικές συναρτήσεις (χωρίς εξάρτηση από toolboxes)
%% ======================================================================

function dw = local_dw(e)
% Durbin-Watson: ~2 = απουσία autocorrelation, <2 = θετική, >2 = αρνητική
den = sum(e.^2);
if den < eps, dw = NaN; return; end
dw = sum(diff(e).^2) / den;
end


function [S,K] = local_moments(e)
% Ασυμμετρία και κύρτωση (population, μη διορθωμένη κύρτωση: 3 = Gaussian)
n  = numel(e);
ec = e - mean(e);
s2 = sum(ec.^2)/n;
if s2 < eps, S = NaN; K = NaN; return; end
S = (sum(ec.^3)/n) / s2^1.5;
K = (sum(ec.^4)/n) / s2^2;
end


function pval = local_jb(e)
% Jarque-Bera. Η στατιστική ακολουθεί ασυμπτωτικά chi2 με 2 β.ε., για την
% οποία η ουρά είναι ακριβώς exp(-x/2) - δεν χρειάζεται chi2cdf.
n = numel(e);
if n < 8, pval = NaN; return; end
[S,K] = local_moments(e);
if isnan(S), pval = NaN; return; end
JB = n/6 * (S^2 + (K-3)^2/4);
pval = min(max(exp(-JB/2),0),1);
end


function [pval, acf1] = local_lb(e, L, dof)
% Ljung-Box Q με διόρθωση βαθμών ελευθερίας για residuals μοντέλου.
% Η ουρά της chi2 υπολογίζεται με gammainc (base MATLAB).
n  = numel(e);
ec = e - mean(e);
den = sum(ec.^2);
if den < eps || n < 5, pval = NaN; acf1 = NaN; return; end

L  = min(L, n-2);
rk = zeros(L,1);
for k = 1:L
    rk(k) = sum(ec(1:end-k).*ec(k+1:end)) / den;
end
acf1 = rk(1);

Q    = n*(n+2) * sum( rk.^2 ./ (n - (1:L)') );
dof  = max(1, dof);
pval = gammainc(Q/2, dof/2, 'upper');   % = 1 - chi2cdf(Q, dof)
pval = min(max(pval,0),1);
end


function local_window_figure(eI, eR, r, model_name, w, period_str, opt)
% Αναλυτικό σχήμα ενός παραθύρου: ACF, QQ, residuals vs fitted
% R2019a-safe: subplot αντί για tiledlayout/nexttile

vis = 'off'; if opt.ShowPlots, vis = 'on'; end
fig = figure('Visible',vis,'Color','w','Position',[100 100 1400 680]);

data = {eI, eR};
fit  = {r.infected_model(:), r.removed_model(:)};
lab  = {'Infected','Removed'};

for j = 1:2
    e  = data{j};
    n  = numel(e);
    ec = e - mean(e);
    L  = min(10, max(1,floor(n/3)));
    base = (j-1)*3;                  % γραμμή 1 -> 1..3, γραμμή 2 -> 4..6

    % --- ACF ---
    subplot(2,3,base+1)
    rk = zeros(L,1);
    for k = 1:L
        rk(k) = sum(ec(1:end-k).*ec(k+1:end)) / sum(ec.^2);
    end
    stem(1:L, rk, 'filled'); hold on
    ci = 1.96/sqrt(n);
    plot([0 L+1],[ ci  ci],'r--');
    plot([0 L+1],[-ci -ci],'r--');
    plot([0 L+1],[  0   0],'k-');
    xlim([0 L+1]); ylim([-1 1]);
    xlabel('Lag'); ylabel('ACF');
    title(sprintf('ACF (%s)', lab{j})); grid on

    % --- QQ plot (χωρίς Statistics Toolbox) ---
    subplot(2,3,base+2)
    es = sort(e);
    pp = ((1:n)' - 0.5)/n;
    q  = sqrt(2)*erfinv(2*pp - 1);            % θεωρητικά κανονικά quantiles
    plot(q, es, 'o', 'MarkerFaceColor',[0 0.45 0.74], ...
                     'MarkerEdgeColor','none','MarkerSize',5); hold on
    mu = mean(e); sg = std(e);
    plot(q, mu + sg*q, 'r-','LineWidth',1.2);
    xlabel('Theoretical quantiles'); ylabel('Sample quantiles');
    title(sprintf('QQ plot (%s)', lab{j})); grid on

    % --- residuals vs fitted ---
    subplot(2,3,base+3)
    plot(fit{j}, e, 'o', 'MarkerFaceColor',[0.85 0.33 0.10], ...
                         'MarkerEdgeColor','none','MarkerSize',5); hold on
    xl = xlim; plot(xl,[0 0],'r--'); xlim(xl);
    xlabel('Fitted'); ylabel('Residual');
    title(sprintf('Residuals vs fitted (%s)', lab{j})); grid on
end

local_sgtitle(fig, sprintf('%s - window %d: %s', model_name, w, period_str));

if opt.SavePlots
    fname = sprintf('%s_window%02d_residuals', model_name, w);
    local_save(fig, fullfile(opt.OutputDir, fname));
end
if ~opt.ShowPlots, close(fig); end
end


function local_sgtitle(fig, txt)
% sgtitle υπάρχει από R2018b, με fallback σε annotation
if exist('sgtitle','file') || exist('sgtitle','builtin')
    try
        sgtitle(txt,'FontWeight','bold','FontSize',14);
        return
    catch
    end
end
annotation(fig,'textbox',[0 0.94 1 0.06],'String',txt, ...
    'HorizontalAlignment','center','EdgeColor','none', ...
    'FontWeight','bold','FontSize',14);
end


function local_save(fig, basepath)

if exist('exportgraphics','file')
    exportgraphics(fig,[basepath '.png'],'Resolution',200);
    exportgraphics(fig,[basepath '.pdf'],'ContentType','vector');
else
    set(fig,'PaperPositionMode','auto','InvertHardcopy','off','Color','w');
    print(fig, [basepath '.png'], '-dpng', '-r200');
    print(fig, [basepath '.pdf'], '-dpdf', '-painters', '-bestfit');
end
end
