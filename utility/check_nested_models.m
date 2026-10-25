function check_nested_models(countryName)

if nargin < 1
    countryName = 'Greece';
end

%% ===============================================================
% Load Results
%% ===============================================================

load(fullfile('results','Results_SEIR',countryName,...
    'SEIR_all_periods_results1_C.mat'));
SEIR = all_results;

load(fullfile('results','Results_SEIRS',countryName,...
    'SEIRS_all_periods_results1_C.mat'));
SEIRS = all_results;

assert(numel(SEIR)==numel(SEIRS),...
    'Different number of fitting windows.');

nWindows = numel(SEIR);

fprintf('\n=============================================\n');
fprintf('      NESTED MODEL CONSISTENCY CHECK\n');
fprintf('=============================================\n\n');

tol = 1e-8;

viol_obj = [];
viol_rss = [];

for i = 1:nWindows

    %% ---------- objective -----------------

    if SEIRS{i}.objective > SEIR{i}.objective + tol
        viol_obj(end+1)=i;
    end

    %% ---------- RSS -----------------------

    rss_seir  = SEIR{i}.RSS_I  + SEIR{i}.RSS_R;
    rss_seirs = SEIRS{i}.RSS_I + SEIRS{i}.RSS_R;

    if rss_seirs > rss_seir + tol
        viol_rss(end+1)=i;
    end

end

%% ===============================================================
% Report
%% ===============================================================

fprintf('Windows checked : %d\n\n',nWindows);

fprintf('Objective violations : %d\n',numel(viol_obj));

if isempty(viol_obj)

    fprintf('   ? No objective violations.\n');

else

    fprintf('   Windows : ');
    fprintf('%d ',viol_obj);
    fprintf('\n');

end

fprintf('\n');

fprintf('RSS violations : %d\n',numel(viol_rss));

if isempty(viol_rss)

    fprintf('   ? No RSS violations.\n');

else

    fprintf('   Windows : ');
    fprintf('%d ',viol_rss);
    fprintf('\n');

end

fprintf('\n=============================================\n');

end