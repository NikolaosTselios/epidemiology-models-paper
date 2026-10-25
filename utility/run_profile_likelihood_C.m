%% run_profile_likelihood_C.m
% Αυτόνομο script: profile likelihood για το γ του SIR, σε επιλεγμένα παράθυρα.
%
% ΤΙ ΚΑΝΕΙ
%   Για κάθε επιλεγμένο 30-ήμερο παράθυρο:
%     1. Κλειδώνει το γ σε ένα πλέγμα τιμών.
%     2. Για κάθε τιμή, ξανα-βελτιστοποιεί ΜΟΝΟ το β (profiling).
%     3. Υπολογίζει την profile NLL και βρίσκει το 95% CI (κατώφλι 1.92).
%     4. Σχεδιάζει την καμπύλη και αποθηκεύει σχήμα + CSV.
%
% Συμβατότητα: MATLAB R2019a

clear;
close all;
clc;

%% ===================== Ρυθμίσεις =====================

N          = 1e7;
step_days  = 30;

csv_file   = fullfile('data','Covid_19_Greece.csv');

start_date = "2020-04-22";
end_date   = "2022-07-10";

windows_to_profile = [6 14 22];

%% ==========================================================
%% ΝΕΟ adaptive gamma grid
%% ==========================================================

gamma_grid = unique([ ...
    linspace(0.04,0.30,25), ...
    linspace(0.05,0.10,40) ...
]);

beta_lb = 0.01;
beta_ub = 1.00;

n_starts  = 20;

smooth_win = 7;

out_dir = fullfile('plots','Profile_Likelihood');

if ~exist(out_dir,'dir')
    mkdir(out_dir);
end

thr = 1.92;

%% ===================== Φόρτωση δεδομένων =====================

dataTable = readtable(csv_file);

dates_all = datetime(dataTable.Date);

cum_all = dataTable.Total_Cases;

w = 10;

infected_all = movsum(dataTable.Confirmed_Cases,[w-1 0]);
infected_all = infected_all(:);
infected_all(infected_all<0)=0;

removed_all = cum_all - infected_all;
removed_all = removed_all(:);
removed_all(removed_all<0)=0;

%% ===================== Παράθυρα =====================

periods = generate_fitting_periods(start_date,end_date,step_days);

nP = size(periods,1);

fprintf('Σύνολο παραθύρων: %d\n',nP);
fprintf('Profiling στα παράθυρα: %s\n',mat2str(windows_to_profile));

summary_rows = {};

%% ==========================================================
%% LOOP
%% ==========================================================

for wi = windows_to_profile

    if wi<1 || wi>nP
        warning('Window %d εκτός ορίων.',wi);
        continue
    end

    sd = datetime(periods{wi,1});
    ed = datetime(periods{wi,2});

    idx = dates_all>=sd & dates_all<=ed;

    dates_fit = dates_all(idx);

    if smooth_win>0

        I_fit = smoothdata(infected_all(idx),'movmean',smooth_win);
        R_fit = smoothdata(removed_all(idx),'movmean',smooth_win);

    else

        I_fit = infected_all(idx);
        R_fit = removed_all(idx);

    end

    I_fit = I_fit(:);
    R_fit = R_fit(:);

    T = numel(I_fit);

    I0 = max(I_fit(1),1);
    R0 = R_fit(1);
    S0 = N-I0-R0;

    if S0<=0
        warning('Window %d: S0<=0.',wi);
        continue
    end

    fprintf('\n=========================================\n');
    fprintf('Window %d\n',wi);
    fprintf('%s -> %s\n',datestr(sd),datestr(ed));
    fprintf('=========================================\n');

    %% =====================================================
    %% Ελεύθερο fit
    %% =====================================================

    opts = optimoptions('fmincon',...
        'Display','off',...
        'Algorithm','sqp');

    free_obj = @(p) nll_SIR( ...
        p,...
        S0,I0,R0,N,T,...
        I_fit,R_fit);

    best = inf;
    theta_hat = [0.10 0.10];

    for k=1:40

        p0 = [ ...
            beta_lb + (beta_ub-beta_lb)*rand(), ...
            gamma_grid(1) + (gamma_grid(end)-gamma_grid(1))*rand() ...
            ];

        [p,f] = fmincon( ...
            free_obj,...
            p0,...
            [],[],[],[],...
            [beta_lb gamma_grid(1)],...
            [beta_ub gamma_grid(end)],...
            [],opts);

        if f<best
            best=f;
            theta_hat=p;
        end

    end

    beta_hat  = theta_hat(1);
    gamma_hat = theta_hat(2);

    fprintf('beta_hat  = %.5f\n',beta_hat);
    fprintf('gamma_hat = %.5f\n',gamma_hat);

    %% =====================================================
    %% PROFILE LIKELIHOOD
    %% =====================================================

    nll = nan(size(gamma_grid));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Warm start
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    b_prev = beta_hat;

    for j=1:numel(gamma_grid)

        g = gamma_grid(j);

        obj = @(beta) nll_SIR( ...
            [beta g],...
            S0,I0,R0,N,T,...
            I_fit,R_fit);

        bbest = inf;
                %==========================================================
        % Warm-start + multi-start
        %==========================================================

        b0_list = [ ...
            b_prev, ...
            beta_lb + (beta_ub-beta_lb)*rand(1,n_starts-1) ...
            ];

        for t=1:numel(b0_list)

            [b,f] = fmincon( ...
                obj,...
                b0_list(t),...
                [],[],[],[],...
                beta_lb,...
                beta_ub,...
                [],opts);

            if f<bbest
                bbest = f;
                b_prev = b;
            end

        end

        nll(j)=bbest;

    end

    %% =====================================================
    %% Profile likelihood
    %% =====================================================

    dnll = nll-min(nll);

    [~,imin] = min(nll);

    %% =====================================================
    %% Interpolated 95% CI
    %% =====================================================

    il = find(dnll(1:imin)>thr,1,'last');

    if ~isempty(il)

        ci_lo = interp1( ...
            dnll([il il+1]),...
            gamma_grid([il il+1]),...
            thr);

    else

        ci_lo = gamma_grid(1);

    end

    ir = imin-1 + find(dnll(imin:end)>thr,1,'first');

    if ~isempty(ir)

        ci_hi = interp1( ...
            dnll([ir-1 ir]),...
            gamma_grid([ir-1 ir]),...
            thr);

    else

        ci_hi = gamma_grid(end);

    end

    flat = max(dnll)<thr;

    touches_lo = (ci_lo==gamma_grid(1));
    touches_hi = (ci_hi==gamma_grid(end));

    fprintf('\n');
    fprintf('gamma_hat = %.5f\n',gamma_hat);
    fprintf('95%% CI    = [%.5f , %.5f]\n',ci_lo,ci_hi);

    if flat
        fprintf('WARNING: Flat profile likelihood (parameter may not be identifiable)\n');
    end

    if touches_lo
        fprintf('WARNING: CI touches LOWER boundary\n');
    end

    if touches_hi
        fprintf('WARNING: CI touches UPPER boundary\n');
    end

    %% =====================================================
    %% Plot
    %% =====================================================

    fig = figure( ...
        'Visible','off',...
        'Color','w',...
        'Position',[100 100 720 520]);

    plot( ...
        gamma_grid,...
        dnll,...
        'b-o',...
        'LineWidth',1.4,...
        'MarkerFaceColor',[0.2 0.4 0.8]);

    hold on

    plot( ...
        [gamma_grid(1) gamma_grid(end)],...
        [thr thr],...
        'r--',...
        'LineWidth',1.2);

    plot( ...
        gamma_hat,...
        0,...
        'g*',...
        'MarkerSize',13,...
        'LineWidth',1.5);

    plot([ci_lo ci_lo],[0 thr],'k:');
    plot([ci_hi ci_hi],[0 thr],'k:');

    xlabel('\gamma');
    ylabel('Profile NLL - min');

    title(sprintf( ...
        'SIR profile likelihood | window %d',wi));

    grid on

    base = fullfile( ...
        out_dir,...
        sprintf('profile_gamma_SIR_window%02d',wi));

    save_fig(fig,base);

    close(fig);

    summary_rows(end+1,:) = { ...
        wi,...
        sprintf('%s to %s', ...
        datestr(sd,'yyyy-mm-dd'),...
        datestr(ed,'yyyy-mm-dd')),...
        gamma_hat,...
        ci_lo,...
        ci_hi,...
        flat,...
        (touches_lo||touches_hi)};

end

%% =====================================================
%% Summary
%% =====================================================

if ~isempty(summary_rows)

    Tsum = cell2table( ...
        summary_rows,...
        'VariableNames',{ ...
        'Window',...
        'Period',...
        'gamma_hat',...
        'CI_low',...
        'CI_high',...
        'Flat_unidentifiable',...
        'CI_touches_bound'});

    disp(Tsum)

    writetable( ...
        Tsum,...
        fullfile( ...
        out_dir,...
        'profile_gamma_summary.csv'));

end

fprintf('\n');
fprintf('=============================================\n');
fprintf('Profile likelihood ολοκληρώθηκε.\n');
fprintf('=============================================\n');

%% =====================================================
%% Local functions
%% =====================================================

function c = nll_SIR(params,S0,I0,R0,N,T,I_data,R_data)

Isim = sim_SIR(params,S0,I0,R0,N,T,'I');
Rsim = sim_SIR(params,S0,I0,R0,N,T,'R');

eI = I_data(:)-Isim(:);
eR = R_data(:)-Rsim(:);

RSS_I = sum(eI.^2);
RSS_R = sum(eR.^2);

if RSS_I<=0 || RSS_R<=0
    c = 1e12;
    return
end

c = 0.5*T*log(RSS_I/T) + ...
    0.5*T*log(RSS_R/T);

end

function out = sim_SIR(params,S0,I0,R0,N,T,which)

beta  = params(1);
gamma = params(2);

S=zeros(1,T);
I=zeros(1,T);
R=zeros(1,T);

S(1)=S0;
I(1)=I0;
R(1)=R0;

for k=1:T-1

    dS = -beta*S(k)*I(k)/N;
    dI =  beta*S(k)*I(k)/N - gamma*I(k);
    dR =  gamma*I(k);

    S(k+1)=S(k)+dS;
    I(k+1)=I(k)+dI;
    R(k+1)=R(k)+dR;

end

switch which

    case 'I'
        out = I(:);

    case 'R'
        out = R(:);

end

end

function save_fig(fig,basepath)

if exist('exportgraphics','file')

    exportgraphics( ...
        fig,...
        [basepath '.png'],...
        'Resolution',200);

    exportgraphics( ...
        fig,...
        [basepath '.pdf'],...
        'ContentType','vector');

else

    set(fig,...
        'PaperPositionMode','auto',...
        'InvertHardcopy','off',...
        'Color','w');

    print(fig,...
        [basepath '.png'],...
        '-dpng',...
        '-r200');

    print(fig,...
        [basepath '.pdf'],...
        '-dpdf',...
        '-painters',...
        '-bestfit');

end

end
