function summary = compare_SEIRS_periods_C(csv_filename, country_name, N, weights, normalize, num_trials, step_days)
%% compare_SEIRS_periods.m
% Εκτελεί fitting του μοντέλου SEIRS σε διαδοχικά (disjoint) χρονικά παράθυρα.
%
% Για κάθε χρονικό παράθυρο πραγματοποιείται προσαρμογή (fitting) του
% μοντέλου SEIRS και αποθηκεύονται οι αντίστοιχες εκτιμήσεις των παραμέτρων,
% καθώς και οι μετρικές RSS, MSE, MAE, AIC και AICc.
%
% Το συνολικό AIC_tot και AICc_tot υπολογίζονται ως άθροισμα των
% αντίστοιχων τιμών όλων των παραθύρων.
%
% Το summary που επιστρέφεται περιέχει:
%    τα αποτελέσματα κάθε παραθύρου,
%    τις συνολικές μετρικές,
%    τον πίνακα αποτελεσμάτων,
%    τους δείκτες των περιόδων που χρησιμοποιήθηκαν.
%
% Η αποθήκευση των period_idx είναι απαραίτητη ώστε τα μοντέλα
% SIR, SEIR και SEIRS να συγκρίνονται στα ίδια ακριβώς χρονικά
% παράθυρα κατά τον υπολογισμό των συνολικών AIC/AICc.

fprintf('\n=== Running compare_SEIRS_periods for %s ===\n', country_name);

model_name = 'SEIRS';

%% === Φόρτωση δεδομένων ===

dataTable = readtable(csv_filename);

%% === Δημιουργία χρονικών παραθύρων προσαρμογής ===

start_date = "2020-04-22";
end_date   = "2022-07-10";

periods = generate_fitting_periods(start_date,end_date,step_days);
nP = size(periods,1);

fprintf('Δημιουργήθηκαν %d χρονικά παράθυρα (%d ημέρες το καθένα).\n', ...
    nP,step_days);

%% === Δημιουργία φακέλων αποθήκευσης ===

plots_dir   = fullfile('plots','Plots_SEIRS',country_name);
results_dir = fullfile('results','Results_SEIRS',country_name);

if ~exist(plots_dir,'dir')
    mkdir(plots_dir);
end

if ~exist(results_dir,'dir')
    mkdir(results_dir);
end

%% === Αρχικοποίηση μεταβλητών ===

all_results = cell(nP,1);

% Δείκτες των παραθύρων που χρησιμοποιήθηκαν επιτυχώς
period_idx = [];

Period = strings(0,1);

Beta  = [];
Gamma = [];
Alpha = [];
Delta = [];
R0    = [];
Tw    = [];

RSS_I = [];
RSS_R = [];

MSE_I = [];
MAE_I = [];

MSE_R = [];
MAE_R = [];

AICw  = [];
AICcw = [];

%% === Προσαρμογή του μοντέλου σε κάθε χρονικό παράθυρο ===

for i = 1:nP

    fprintf('\n=== Fitting period %d/%d ===\n',i,nP);

    start_date_dt = datetime(periods{i,1});
    end_date_dt   = datetime(periods{i,2});

    res = fit_SEIRS_model_C( ...
            dataTable,...
            start_date_dt,...
            end_date_dt,...
            N,...
            weights,...
            normalize,...
            num_trials);

    % Παράλειψη παραθύρου εάν δεν επιστράφηκε αποτέλεσμα
    if isempty(res)

        fprintf('Skipping period %d.\n',i);
        continue;

    end

    all_results{i} = res;

    % Αποθήκευση του δείκτη του παραθύρου
    period_idx(end+1,1) = i; %#ok<AGROW>

    % Χρονικό διάστημα του παραθύρου
    Period(end+1,1) = sprintf('%s to %s', ...
        datestr(start_date_dt,'yyyy-mm-dd'), ...
        datestr(end_date_dt,'yyyy-mm-dd')); %#ok<AGROW>

    % Εκτιμημένες παράμετροι
    Beta(end+1,1)  = res.beta;
    Gamma(end+1,1) = res.gamma;
    Alpha(end+1,1) = res.alpha;
    Delta(end+1,1) = res.delta;
    R0(end+1,1)    = res.R0;
    Tw(end+1,1)    = res.T;

    % RSS
    RSS_I(end+1,1) = res.RSS_I;
    RSS_R(end+1,1) = res.RSS_R;

    % MSE / MAE
    MSE_I(end+1,1) = res.mse_I;
    MAE_I(end+1,1) = res.mae_I;

    MSE_R(end+1,1) = res.mse_R;
    MAE_R(end+1,1) = res.mae_R;

    % Κριτήρια επιλογής μοντέλου
    AICw(end+1,1)  = res.AIC;
    AICcw(end+1,1) = res.AICc;

    %% === Αποθήκευση γραφήματος ===

    %fig_name = sprintf( ...
        %'SEIRS_fit_%s_to_%s_beta%.3f_gamma%.3f_alpha%.3f_delta%.3f_R0%.2f.png', ...
        %datestr(start_date_dt,'yyyy-mm-dd'), ...
        %datestr(end_date_dt,'yyyy-mm-dd'), ...
        %res.beta,...
        %res.gamma,...
        %res.alpha,...
        %res.delta,...
        %res.R0);

    %saveas(gcf, fullfile(plots_dir,fig_name));

    %close(gcf);

end
%% === Πίνακας αποτελεσμάτων ανά χρονικό παράθυρο ===
%
% Ο πίνακας αυτός περιέχει όλες τις εκτιμημένες παραμέτρους και τις
% μετρικές αξιολόγησης για κάθε χρονικό παράθυρο ξεχωριστά.
%
% Χρησιμοποιείται τόσο για αποθήκευση σε αρχείο CSV όσο και για
% μεταγενέστερη ανάλυση ή σύγκριση μεταξύ διαφορετικών μοντέλων.

results_table = table( ...
    Period,...
    Beta,...
    Gamma,...
    Alpha,...
    Delta,...
    R0,...
    Tw,...
    RSS_I,...
    RSS_R,...
    MSE_I,...
    MAE_I,...
    MSE_R,...
    MAE_R,...
    AICw,...
    AICcw,...
    'VariableNames', ...
    {'Period','Beta','Gamma','Alpha','Delta','R0','T', ...
     'RSS_I','RSS_R', ...
     'MSE_I','MAE_I', ...
     'MSE_R','MAE_R', ...
     'AIC','AICc'});

disp(results_table);

%% === Συνολικές μετρικές όλων των χρονικών παραθύρων ===
%
% Σύμφωνα με τη θεωρία του AIC, όταν οι περίοδοι είναι ανεξάρτητες,
% το συνολικό AIC προκύπτει ως άθροισμα των επιμέρους τιμών.
%
% Το ίδιο ισχύει και για το AICc.

AIC_tot  = sum(AICw);
AICc_tot = sum(AICcw);

fprintf( ...
    '\n[%s] windows=%d | RSS_I_tot=%.3g | RSS_R_tot=%.3g | AIC_tot=%.2f | AICc_tot=%.2f\n', ...
    model_name,...
    numel(period_idx),...
    sum(RSS_I),...
    sum(RSS_R),...
    AIC_tot,...
    AICc_tot);

%% === Δημιουργία δομής εξόδου ===
%
% Η δομή summary χρησιμοποιείται από το run_fitting_models ώστε
% να πραγματοποιηθεί η τελική σύγκριση μεταξύ SIR, SEIR και SEIRS
% χρησιμοποιώντας τα ίδια χρονικά παράθυρα.

summary = struct();

summary.model       = model_name;

% Δείκτες των παραθύρων που χρησιμοποιήθηκαν επιτυχώς
summary.period_idx  = period_idx;

% Συνολικός αριθμός παραθύρων
summary.W           = numel(period_idx);

% Μετρικές ανά χρονικό παράθυρο
summary.AIC_w       = AICw;
summary.AICc_w      = AICcw;

summary.RSS_I_w     = RSS_I;
summary.RSS_R_w     = RSS_R;

% Συνολικές μετρικές
summary.AIC_tot     = AIC_tot;
summary.AICc_tot    = AICc_tot;

% Αναλυτικός πίνακας αποτελεσμάτων
summary.table       = results_table;

% Όλα τα αποτελέσματα των επιμέρους fittings
summary.all_results = all_results;

%% === Αποθήκευση αποτελεσμάτων ===

save( ...
    fullfile(results_dir,'SEIRS_all_periods_results1_C.mat'), ...
    'summary', ...
    'all_results', ...
    'results_table');

writetable( ...
    results_table,...
    fullfile(results_dir,'SEIRS_all_periods_results1_C.csv'));

fprintf('\n=== Ολοκληρώθηκε το SEIRS fitting για τη χώρα %s ===\n', ...
        country_name);

end