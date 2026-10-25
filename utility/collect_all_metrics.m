function collect_all_metrics()
% COLLECT_ALL_METRICS
% Builds a consolidated record of every metric of interest across all
% fitting windows and all three models (SIR / SEIR / SEIRS), then surfaces
% candidate representative windows so you don't have to scan 27x3 rows by eye.
%
% OUTPUTS:
%   Metrics_long.csv     one row per (window x model) -- the full record
%   Window_summary.csv   one row per window -- cross-model comparison
%   Console printout      ranked candidate windows + aggregate totals (Table 1)
%
% ASSUMPTIONS about each *_all_periods_results1_C.mat:
%   contains a cell array `all_results`, and each all_results{i} has at least
%     .dates .infected_real .infected_model .removed_real .removed_model .R0
%   (these are the fields used by your plot_model_comparisons script).
%
% Fit metrics and residual diagnostics are RECOMPUTED from the real/model
% arrays, so the script does not depend on stored metric field names and
% doubles as a verification of your existing numbers. Model parameters are
% auto-detected (named fields or a params vector); missing ones stay NaN.
%
% Diagnostics use only base MATLAB (gammainc for the chi-square CDF), so no
% Statistics/Econometrics toolbox is required.

clc;
countryName = 'Greece';

%% ---- Load the three result sets ----
S_SIR   = load(fullfile('results','Results_SIR',   countryName, 'SIR_all_periods_results1_C.mat'));
S_SEIR  = load(fullfile('results','Results_SEIR',  countryName, 'SEIR_all_periods_results1_C.mat'));
S_SEIRS = load(fullfile('results','Results_SEIRS', countryName, 'SEIRS_all_periods_results1_C.mat'));

R.SIR   = S_SIR.all_results;
R.SEIR  = S_SEIR.all_results;
R.SEIRS = S_SEIRS.all_results;

models  = {'SIR','SEIR','SEIRS'};
nParams = struct('SIR',2,'SEIR',3,'SEIRS',4);   % beta,gamma[,alpha[,delta]]

nP = numel(R.SIR);
assert(numel(R.SEIR)==nP && numel(R.SEIRS)==nP, 'Period count mismatch between models');

% Show available fields once so you can adjust parameter detection if needed
fprintf('Fields available in one result struct:\n');

SIR_results = R.SIR;
fields = fieldnames(SIR_results{1});

fprintf('%s ', fields{:});
fprintf('\n\n');

%% ---- LONG table: one row per (window x model) ----
longRows = cell(nP*numel(models), 30);
row = 0;
for m = 1:numel(models)
    mdl = models{m};
    model_results = R.(mdl);

    for i = 1:nP
        res = model_results{i};

        Ir = res.infected_real(:);  Im = res.infected_model(:);
        Rr = res.removed_real(:);   Rm = res.removed_model(:);
        T  = numel(Ir);

        rI = Ir - Im;   % infected residuals
        rR = Rr - Rm;   % removed  residuals

        % --- fit metrics ---
        RMSE_I = sqrt(mean(rI.^2));  MAE_I = mean(abs(rI));
        RMSE_R = sqrt(mean(rR.^2));  MAE_R = mean(abs(rR));

        % --- AICc (diagonal Gaussian, separate variances, k=p+2, n=2T) ---
        AICc = local_aicc(rI, rR, nParams.(mdl));

        % --- residual diagnostics ---
        [DW_I,JBp_I,LBp_I,RB_I,SK_I,KU_I] = local_diag(rI);
        [DW_R,JBp_R,LBp_R,RB_R,SK_R,KU_R] = local_diag(rR);

        % --- parameters (auto-detected; NaN if absent) ---
        pv = get_params_vec(res, mdl);   % [beta gamma alpha delta]
        R0 = get_field(res, {'R0','r0'}, NaN);

        % --- window shape indicators (from real infected) ---
        [rangeRatio, nTurns, interior] = local_shape(Ir);

        row = row + 1;
        longRows(row,:) = { i, mdl, ...
            datestr(res.dates(1),'yyyy-mm-dd'), datestr(res.dates(end),'yyyy-mm-dd'), T, ...
            pv(1), pv(2), pv(3), pv(4), R0, ...
            RMSE_I, MAE_I, RMSE_R, MAE_R, AICc, ...
            DW_I, JBp_I, LBp_I, RB_I, SK_I, KU_I, ...
            DW_R, JBp_R, LBp_R, RB_R, SK_R, KU_R, ...
            rangeRatio, nTurns, interior };
    end
end

longVars = {'Window','Model','StartDate','EndDate','T', ...
    'beta','gamma','alpha','delta','R0', ...
    'RMSE_I','MAE_I','RMSE_R','MAE_R','AICc', ...
    'DW_I','JBp_I','LBp_I','RelBias_I','Skew_I','Kurt_I', ...
    'DW_R','JBp_R','LBp_R','RelBias_R','Skew_R','Kurt_R', ...
    'RangeRatio_I','Turns_I','InteriorExtremum_I'};
Tlong = cell2table(longRows, 'VariableNames', longVars);
writetable(Tlong, 'Metrics_long_Greece.csv');
fprintf('Wrote Metrics_long.csv  (%d rows)\n', height(Tlong));

%% ---- WIDE table: one row per window (cross-model comparison) ----
rmseI = @(res) sqrt(mean((res.infected_real(:) - res.infected_model(:)).^2));

wideRows = cell(nP, 14);
for i = 1:nP
    resS = R.SIR{i};
    T = numel(resS.infected_real);

    rI_SIR   = rmseI(R.SIR{i});
    rI_SEIR  = rmseI(R.SEIR{i});
    rI_SEIRS = rmseI(R.SEIRS{i});

    impr1 = 100*(rI_SIR  - rI_SEIR ) / rI_SIR;    % SIR  -> SEIR  improvement (%): exposed helps
    impr2 = 100*(rI_SEIR - rI_SEIRS) / rI_SEIR;   % SEIR -> SEIRS improvement (%): immunity helps

        SEIRS_results = R.SEIRS;
    pvS = get_params_vec(SEIRS_results{i}, 'SEIRS');
    delta_SEIRS = pvS(4);

    [rangeRatio, nTurns, interior] = local_shape(resS.infected_real(:));

    [~, bi] = min([rI_SIR rI_SEIR rI_SEIRS]);
    bestModel = models{bi};

    wideRows(i,:) = { i, ...
        datestr(resS.dates(1),'yyyy-mm-dd'), datestr(resS.dates(end),'yyyy-mm-dd'), T, ...
        rI_SIR, rI_SEIR, rI_SEIRS, impr1, impr2, delta_SEIRS, bestModel, ...
        rangeRatio, nTurns, interior };
end

wideVars = {'Window','StartDate','EndDate','T', ...
    'RMSE_I_SIR','RMSE_I_SEIR','RMSE_I_SEIRS', ...
    'Impr_SIR_SEIR_pct','Impr_SEIR_SEIRS_pct','delta_SEIRS','BestByRMSE_I', ...
    'RangeRatio_I','Turns_I','InteriorExtremum_I'};
Twide = cell2table(wideRows, 'VariableNames', wideVars);
writetable(Twide, 'Window_summary_Greece.csv');
fprintf('Wrote Window_summary.csv  (%d rows)\n', height(Twide));

%% ---- Candidate windows (printed to console) ----
k = min(6, nP);

fprintf('\n=== Exposed compartment helps most (rank by SIR->SEIR improvement) ===\n');
[~, ord] = sort(Twide.Impr_SIR_SEIR_pct, 'descend');
disp(Twide(ord(1:k), {'Window','StartDate','EndDate','Impr_SIR_SEIR_pct','RMSE_I_SIR','RMSE_I_SEIR'}));

fprintf('\n=== Immunity helps most (rank by SEIR->SEIRS improvement) ===\n');
[~, ord] = sort(Twide.Impr_SEIR_SEIRS_pct, 'descend');
disp(Twide(ord(1:k), {'Window','StartDate','EndDate','Impr_SEIR_SEIRS_pct','delta_SEIRS'}));

fprintf('\n=== Hardest windows (all models struggle: highest RMSE_I_SEIRS) ===\n');
[~, ord] = sort(Twide.RMSE_I_SEIRS, 'descend');
disp(Twide(ord(1:k), {'Window','StartDate','EndDate','RMSE_I_SEIRS','InteriorExtremum_I','Turns_I'}));

fprintf('\n=== Windows where delta is active (> 1e-4) ===\n');
act = Twide(Twide.delta_SEIRS > 1e-4, {'Window','StartDate','EndDate','delta_SEIRS','Impr_SEIR_SEIRS_pct'});
if isempty(act), fprintf('   (none)\n'); else, disp(act); end

%% ---- Aggregate totals (these feed the summary Table 1) ----
fprintf('\n=== Aggregate totals across all windows (for Table 1) ===\n');
fprintf('%-6s  %12s  %12s  %12s\n','Model','sumRMSE_I','sumMAE_I','sumAICc');
for m = 1:numel(models)
    mdl = models{m};
    idx = strcmp(Tlong.Model, mdl);
    fprintf('%-6s  %12.1f  %12.1f  %12.1f\n', mdl, ...
        sum(Tlong.RMSE_I(idx)), sum(Tlong.MAE_I(idx)), sum(Tlong.AICc(idx),'omitnan'));
end
fprintf('\nDone.\n');
end


%% ======================= local helper functions =======================

function AICc = local_aicc(rI, rR, p)
% AICc for joint I/R fit: diagonal Gaussian likelihood with separate
% MLE variances for I and R, k = p + 2 (params + 2 variances), n = 2T.
T   = numel(rI);
s2I = mean(rI.^2);  if s2I <= 0, s2I = eps; end
s2R = mean(rR.^2);  if s2R <= 0, s2R = eps; end
logL = -0.5*( T*log(2*pi*s2I) + T + T*log(2*pi*s2R) + T );
k = p + 2;
n = 2*T;
AIC = 2*k - 2*logL;
denom = n - k - 1;
if denom > 0
    AICc = AIC + 2*k*(k+1)/denom;
else
    AICc = NaN;   % too few points for the AICc correction
end
end


function [DW, JBp, LBp, RB, SK, KU] = local_diag(r)
% Residual diagnostics computed from scratch (base MATLAB only).
r = r(:);
n = numel(r);

DW = sum(diff(r).^2) / sum(r.^2);          % Durbin-Watson (ideal 2)
RB = abs(mean(r)) / sqrt(mean(r.^2));      % |mean| / RMSE  (relative bias)

m = mean(r);
s = std(r, 1);                             % population std
if s <= 0, s = eps; end
SK = mean(((r - m)/s).^3);                 % skewness
KU = mean(((r - m)/s).^4);                 % kurtosis (normal = 3)

% Jarque-Bera p-value (chi-square, 2 dof)
JB  = n/6 * (SK^2 + (KU - 3)^2 / 4);
JBp = 1 - chi2cdf_local(JB, 2);

% Ljung-Box p-value over h lags
h   = min(10, n - 1);
rho = local_acf(r, h);
kk  = (1:h)';
Q   = n*(n+2) * sum( (rho.^2) ./ (n - kk) );
LBp = 1 - chi2cdf_local(Q, h);
end


function rho = local_acf(r, h)
% Sample autocorrelation at lags 1..h.
r = r(:) - mean(r);
denom = sum(r.^2);
n = numel(r);
rho = zeros(h,1);
if denom <= 0, return; end
for lag = 1:h
    rho(lag) = sum(r(1+lag:n) .* r(1:n-lag)) / denom;
end
end


function p = chi2cdf_local(x, k)
% Chi-square CDF via the regularized lower incomplete gamma (base MATLAB).
p = gammainc(x/2, k/2);
end


function [rangeRatio, nTurns, interior] = local_shape(x)
% Shape indicators of the observed infected curve within a window.
x = x(:);
rangeRatio = (max(x) - min(x)) / mean(x);          % how much it moves

dx  = diff(x);
sgn = sign(dx);
sgn(sgn == 0) = [];
nTurns = sum(diff(sgn) ~= 0);                      % slope sign changes

[~, imax] = max(x);
[~, imin] = min(x);
n    = numel(x);
edge = max(2, round(0.15*n));
interior = double( (imax > edge && imax < n-edge) || (imin > edge && imin < n-edge) );
end


function pv = get_params_vec(res, mdl)
% Returns [beta gamma alpha delta] with NaN where not applicable / not found.
pv    = nan(1,4);
named = {'beta','gamma','alpha','delta'};
got   = false;
for j = 1:4
    if isfield(res, named{j}) && ~isempty(res.(named{j}))
        pv(j) = res.(named{j})(1);
        got = true;
    end
end
if ~got
    % Fallback: a single parameter vector stored under one of these names
    for vf = {'params','theta','p','x','popt','par'}
        if isfield(res, vf{1}) && isnumeric(res.(vf{1})) && ~isempty(res.(vf{1}))
            v = res.(vf{1})(:)';
            switch mdl
                case 'SIR',   np = 2;
                case 'SEIR',  np = 3;
                case 'SEIRS', np = 4;
            end
            np = min(np, numel(v));
            pv(1:np) = v(1:np);
            break;
        end
    end
end
end


function v = get_field(res, names, default)
% First matching numeric field, else default.
v = default;
for j = 1:numel(names)
    if isfield(res, names{j})
        val = res.(names{j});
        if isnumeric(val) && ~isempty(val)
            v = val(1);
            return;
        end
    end
end
end
