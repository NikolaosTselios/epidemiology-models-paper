function results = fit_SEIRS_model_C(dataTable, start_date, end_date, N, weights, normalize, num_trials)
% fit_SEIRS_model - Εκτελεί fitting του μοντέλου SEIRS σε ένα χρονικό παράθυρο.
%
% Επιστρέφει, εκτός από τις παραμέτρους, τα per-window RSS_I, RSS_R και το
% AIC/AICc του παραθύρου. Το συνολικό AIC_tot = sum(AIC ανά window) και το
% AICc_tot = sum(AICc ανά window) υπολογίζονται έξω από τη συνάρτηση, στον
% driver που καλεί το fit_SEIRS_model μία φορά ανά disjoint window.
%
% Στατιστική υπόθεση: diagonal Gaussian, με ξεχωριστές διακυμάνσεις για I και R.
% k = p + 2 (p = αριθμός παραμέτρων του ODE μοντέλου, +2 για σ_I^2 και σ_R^2).
% n = 2*T (οι χρονοσειρές I και R στοιβαγμένες).
% Το AIC δίνεται μέχρι μία κοινή προσθετική σταθερά για όλα τα μοντέλα και
% χρησιμοποιείται ως σχετικός δείκτης επιλογής (relative ranking), όχι ως
% απόλυτη τιμή AIC.

% === Προεπεξεργασία δεδομένων ===
dates     = datetime(dataTable.Date);
confirmed = dataTable.Confirmed_Cases_corrected;
%deaths    = dataTable.Total_Deaths; %#ok<NASGU>  % Δεν χρησιμοποιείται εδώ.

% Συσσωρευτικά κρούσματα
cumulative_cases = dataTable.Total_Cases;

% Ενεργά κρούσματα = άθροισμα των επιβεβαιωμένων κρουσμάτων των τελευταίων 10 ημερών
window   = 10;
infected = movsum(confirmed,[window-1 0]);      % Για κάθε ημέρα: sum(t-13:t)
infected = infected(:);
infected(infected<0)=0;

removed = cumulative_cases - infected;
removed = removed(:);
removed(removed<0)=0;

% === Επιλογή χρονικού παραθύρου προσαρμογής (fitting) ===
idx          = (dates>=start_date & dates<=end_date);
dates_fit    = dates(idx);
infected_fit = smoothdata(infected(idx),'movmean',7);   % Εξομαλυμένη χρονοσειρά I
removed_fit  = smoothdata(removed(idx),'movmean',7);    % Εξομαλυμένη χρονοσειρά R

T = length(dates_fit);

% === Αρχικές συνθήκες του μοντέλου SEIRS (σταθερές - δεν γίνεται fitting) ===
I0      = max(infected_fit(1),1);     % Ελάχιστη τιμή I για αποφυγή μηδενικής αρχικής κατάστασης
R0_init = removed_fit(1);             % Αρχική τιμή Removed
E0      = I0;                         % Υπόθεση: αρχικά εκτεθειμένοι = αρχικά μολυσμένοι
S0      = N - I0 - R0_init - E0;      % Αρχικός αριθμός ευπαθών

assert(S0>0,'S0 <= 0: Ελέγξτε το N ή τις αρχικές συνθήκες.');      % === ΝΕΟ (guard) ===

% === Ορισμός της αντικειμενικής συνάρτησης σφάλματος ===
wI = weights(1);
wR = weights(2);

objective = @(params) compute_error(params,S0,E0,I0,R0_init,N,T,...
                                    infected_fit,removed_fit,...
                                    wI,wR,normalize);

% === Fitting με multi-start fmincon ===
opts = optimoptions('fmincon','Display','off','Algorithm','sqp');

bounds = [ ...
    0.01  1.00;
    0.05  0.50;
    0.05  0.5;
    0.0 0.012];

min_error = inf;
best_params = [];

for trial = 1:num_trials

    % Lower και upper bounds
    lo = max(bounds(:,1),1e-6);
    ub = bounds(:,2);

    % Log-uniform αρχικοποίηση
    init = exp( ...
        log(lo) + ...
        rand(size(lo)).*(log(ub)-log(lo)));

    init = init';

    [params,fval] = fmincon( ...
        objective,...
        init,...
        [],[],[],[],...
        bounds(:,1),...
        bounds(:,2),...
        [],opts);

    if fval < min_error

        min_error = fval;
        best_params = params;

    end

end

assert(~isempty(best_params), ...
    'Το fmincon δεν επέστρεψε λύση σε κανένα trial.');        % === ΝΕΟ ===

% === Έλεγχος ενεργών περιορισμών (η ερμηνεία ως MLE απαιτεί εσωτερικό ελάχιστο) ===
tol = 1e-4;                                                   % === ΝΕΟ ===
bp = best_params(:);

if any(abs(bp-bounds(:,1))<tol | abs(bp-bounds(:,2))<tol)

    warning('fit_SEIRS_model:binding',...
        'Window %s%s: παράμετρος βρίσκεται στο όριο (binding constraint)  ερμηνεύστε το AIC με προσοχή.',...
        datestr(start_date,'yyyy-mm-dd'),...
        datestr(end_date,'yyyy-mm-dd'));

end

% === Τελικές βέλτιστες παράμετροι ===
beta_fit  = best_params(1);
gamma_fit = best_params(2);
alpha_fit = best_params(3);
delta_fit = best_params(4);

R0_est = beta_fit/gamma_fit;

% === Προσομοίωση του μοντέλου SEIRS ===
[I_model,R_model] = simulate(best_params,S0,E0,I0,R0_init,N,T);

% === Υπολογισμός residuals (χωρίς normalization και χωρίς βάρη) ως προς τις εξομαλυμένες χρονοσειρές ===
resI = infected_fit(:)-I_model(:);
resR = removed_fit(:)-R_model(:);

% === RSS (άθροισμα τετραγώνων σφαλμάτων) και βασικές μετρικές για τα δύο compartments ===
RSS_I = sum(resI.^2);
mse_I = mean(resI.^2);
mae_I = mean(abs(resI));

RSS_R = sum(resR.^2);
mse_R = mean(resR.^2);
mae_R = mean(abs(resR));

% === Υπολογισμός AIC / AICc ανά παράθυρο (diagonal Gaussian με δύο διακυμάνσεις) ===
p = 4;                 % SEIRS: beta, gamma, alpha, delta
k = p + 2;             % + σ_I² και σ_R²
n = 2*T;               % Οι χρονοσειρές I και R στοιβαγμένες

assert(RSS_I>0 && RSS_R>0,...
    'RSS = 0 -> log(-Inf). Ελέγξτε το συγκεκριμένο παράθυρο.');

assert(n-k-1>0,...
    'Το παράθυρο είναι πολύ μικρό για τον υπολογισμό του AICc (n-k-1 <= 0).');

% AIC μέχρι κοινή προσθετική σταθερά (έγκυρο για σχετική σύγκριση μοντέλων)
AIC  = T*log(RSS_I/T) + T*log(RSS_R/T) + 2*k;
AICc = AIC + 2*k*(k+1)/(n-k-1);

% === Δημιουργία δομής εξόδου ===
results = struct(...
    'model','SEIRS',...
    'beta',beta_fit,...
    'gamma',gamma_fit,...
    'alpha',alpha_fit,...
    'delta',delta_fit,...
    'R0',R0_est,...
    'T',T,...
    'p',p,'k',k,'n',n,...
    'RSS_I',RSS_I,...
    'RSS_R',RSS_R,...
    'objective', min_error, ...
    'mse_I',mse_I,...
    'mae_I',mae_I,...
    'mse_R',mse_R,...
    'mae_R',mae_R,...
    'AIC',AIC,...
    'AICc',AICc,...
    'dates',dates_fit,...
    'infected_real',infected_fit,...
    'removed_real',removed_fit,...
    'infected_model',I_model,...
    'removed_model',R_model);
% === Γράφημα αποτελεσμάτων ===
%figure('Name','SEIRS Fitting');
%ModelName = 'SEIRS';
%
%subplot(2,1,1);
%plot(dates_fit, infected_fit, 'r-', 'LineWidth', 2); hold on;
%plot(dates_fit, I_model, 'k--', 'LineWidth', 1.5);
%legend({'Infected (Real)', 'Infected (Model)'}, 'Location', 'best');
%ylabel('Infected');
%title(sprintf('%s Fitting | %s to %s | R0 = %.2f | AICc = %.1f', ...
%    ModelName, ...
%    datestr(dates_fit(1),'yyyy-mm-dd'), ...
%    datestr(dates_fit(end),'yyyy-mm-dd'), ...
%    R0_est, AICc));
%grid on;
%
%subplot(2,1,2);
%plot(dates_fit, removed_fit, 'b-', 'LineWidth', 2); hold on;
%plot(dates_fit, R_model, 'k--', 'LineWidth', 1.5);
%legend({'Removed (Real)', 'Removed (Model)'}, 'Location', 'best');
%xlabel('Ημερομηνία');
%ylabel('Removed');
%grid on;

%% === Εσωτερικές συναρτήσεις ===

function err = compute_error(params, S0, E0, I0, R0_0, N, T, ...
                             I_data, R_data, wI, wR, use_norm)

    [I_sim, R_sim] = simulate(params, S0, E0, I0, R0_0, N, T);

    I_data = I_data(:);
    R_data = R_data(:);

    I_sim = I_sim(:);
    R_sim = R_sim(:);

    if use_norm

        sI = max(eps, std(I_data));
        sR = max(eps, std(R_data));

        % Κανονικοποίηση των residuals με την τυπική απόκλιση
        rI = (I_data - I_sim) / sI;
        rR = (R_data - R_sim) / sR;

        % Σταθμισμένο Mean Squared Error
        err = wI * mean(rI.^2) + ...
              wR * mean(rR.^2);

    else

        % Σταθμισμένο Mean Squared Error χωρίς κανονικοποίηση
        err = wI * mean((I_data - I_sim).^2) + ...
              wR * mean((R_data - R_sim).^2);

    end

end

function [I_out, R_out] = simulate(params, S0, E0, I0, R0_0, N, T)

    % Εκτελεί προσομοίωση του μοντέλου SEIRS και επιστρέφει
    % τις χρονοσειρές των Infectious (I) και Removed (R)
    % για χρονικό ορίζοντα T ημερών.

    beta  = params(1);
    gamma = params(2);
    alpha = params(3);
    delta = params(4);

    r = 1;

    S = zeros(1,T);
    E = zeros(1,T);
    I = zeros(1,T);
    R = zeros(1,T);

    % Αρχικές συνθήκες
    S(1) = S0;
    E(1) = E0;
    I(1) = I0;
    R(1) = R0_0;

    % Προσομοίωση με διακριτή ολοκλήρωση (Euler, dt = 1 ημέρα)
    for t = 1:T-1

        dS = -r * beta * S(t) * I(t) / N + delta * R(t);
        dE =  r * beta * S(t) * I(t) / N - alpha * E(t);
        dI =  alpha * E(t) - gamma * I(t);
        dR =  gamma * I(t) - delta * R(t);

        S(t+1) = S(t) + dS;
        E(t+1) = E(t) + dE;
        I(t+1) = I(t) + dI;
        R(t+1) = R(t) + dR;

    end

    % Επιστροφή των χρονοσειρών Infectious και Removed
    I_out = I(:);
    R_out = R(:);

end

end