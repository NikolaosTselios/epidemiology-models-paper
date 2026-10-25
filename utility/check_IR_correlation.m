%% check_IR_correlation.m
clc;
country = 'Greece';
models  = {'SIR','SEIR','SEIRS'};
np      = [2 3 4];                       % ODE params per model
sumAICc_diag = zeros(1,3); sumAICc_full = zeros(1,3);

for m = 1:3
    S = load(fullfile('results',['Results_' models{m}],country, ...
                      [models{m} '_all_periods_results1_C.mat']));
    A = S.all_results; nW = numel(A);
    rho = nan(nW,1); aicc_d = nan(nW,1); aicc_f = nan(nW,1);

    for i = 1:nW
        rI = A{i}.infected_real(:) - A{i}.infected_model(:);
        rR = A{i}.removed_real(:)  - A{i}.removed_model(:);
        T  = numel(rI);

        sII = mean(rI.^2); sRR = mean(rR.^2); sIR = mean(rI.*rR);
        rho(i) = sIR / sqrt(sII*sRR);        % residual correlation
        p = np(m); n = 2*T;

        % diagonal (your current AIC): independent I,R
        logL_d = -T*log(2*pi) - (T/2)*log(sII*sRR) - T;
        k_d = p + 2;
        aicc_d(i) = 2*k_d - 2*logL_d + 2*k_d*(k_d+1)/max(n-k_d-1,1);

        % full 2x2 covariance: accounts for I-R correlation
        detS   = max(sII*sRR - sIR^2, eps);
        logL_f = -T*log(2*pi) - (T/2)*log(detS) - T;
        k_f = p + 3;
        aicc_f(i) = 2*k_f - 2*logL_f + 2*k_f*(k_f+1)/max(n-k_f-1,1);
    end

    sumAICc_diag(m) = sum(aicc_d,'omitnan');
    sumAICc_full(m) = sum(aicc_f,'omitnan');
    fprintf('%-6s  I-R residual corr: median %+.2f  range [%+.2f, %+.2f]\n', ...
        models{m}, median(rho,'omitnan'), min(rho), max(rho));
end

fprintf('\n%-6s  %12s  %12s\n','Model','AICc_diag','AICc_full');
for m=1:3, fprintf('%-6s  %12.1f  %12.1f\n', models{m}, sumAICc_diag(m), sumAICc_full(m)); end

[~,od]=sort(sumAICc_diag); [~,of]=sort(sumAICc_full);
fprintf('\nRanking (diagonal): %s < %s < %s\n', models{od});
fprintf('Ranking (full 2x2): %s < %s < %s\n', models{of});
fprintf('==> Ranking %s\n', ternary(isequal(od,of),'UNCHANGED — assumption harmless to your conclusion','CHANGED — must report'));

function s = ternary(c,a,b), if c, s=a; else, s=b; end, end