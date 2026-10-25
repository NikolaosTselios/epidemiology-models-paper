function check_pipeline_reproducibility()

%% Parameters

country_name = 'Greece';

csv_filename = fullfile( ...
    'Covid_Modeling','data','COVID_19_Greece.csv');

N = 1e7;

weights = [0.9 0.1];

normalize = true;

num_trials = 800;

step_days = 30;

fprintf('\n===============================================\n');
fprintf(' PIPELINE REPRODUCIBILITY CHECK\n');
fprintf('===============================================\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(1);

sir1 = compare_SIR_periods_C(...
    csv_filename,country_name,N,...
    weights,normalize,num_trials,step_days);

seir1 = compare_SEIR_periods_C(...
    csv_filename,country_name,N,...
    weights,normalize,num_trials,step_days);

seirs1 = compare_SEIRS_periods_C(...
    csv_filename,country_name,N,...
    weights,normalize,num_trials,step_days);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RUN 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng(2);

sir2 = compare_SIR_periods_C_SUMMARY_CHECK(...
    csv_filename,country_name,N,...
    weights,normalize,num_trials,step_days);

seir2 = compare_SEIR_periods_C_SUMMARY_CHECK(...
    csv_filename,country_name,N,...
    weights,normalize,num_trials,step_days);

seirs2 = compare_SEIRS_periods_C_SUMMARY_CHECK(...
    csv_filename,country_name,N,...
    weights,normalize,num_trials,step_days);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SUMMARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

models = {
    sir1,sir2;
    seir1,seir2;
    seirs1,seirs2};

names = {
    'SIR'
    'SEIR'
    'SEIRS'};

fprintf('\n');

for i=1:3

    r1=models{i,1};
    r2=models{i,2};

    relAIC  = abs(r1.AIC_tot-r2.AIC_tot)/abs(r1.AIC_tot);

    relAICc = abs(r1.AICc_tot-r2.AICc_tot)/abs(r1.AICc_tot);

   RSS1 = sum(r1.RSS_I_w) + sum(r1.RSS_R_w);
   RSS2 = sum(r2.RSS_I_w) + sum(r2.RSS_R_w);

    relRSS=abs(RSS1-RSS2)/RSS1;

    fprintf('%s\n',names{i});

    fprintf('   RSS difference   : %.4f %%\n',100*relRSS);

    fprintf('   AIC difference   : %.4f %%\n',100*relAIC);

    fprintf('   AICc difference  : %.4f %%\n\n',100*relAICc);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RANKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AICc_run1=[sir1.AICc_tot seir1.AICc_tot seirs1.AICc_tot];

AICc_run2=[sir2.AICc_tot seir2.AICc_tot seirs2.AICc_tot];

[~,idx1]=sort(AICc_run1);

[~,idx2]=sort(AICc_run2);

fprintf('Ranking Run 1 : %s > %s > %s\n',...
    names{idx1(1)},...
    names{idx1(2)},...
    names{idx1(3)});

fprintf('Ranking Run 2 : %s > %s > %s\n',...
    names{idx2(1)},...
    names{idx2(2)},...
    names{idx2(3)});

if isequal(idx1,idx2)

    fprintf('\n? Ranking identical.\n');

else

    fprintf('\n? Ranking changed.\n');

end

fprintf('\n===============================================\n');

end