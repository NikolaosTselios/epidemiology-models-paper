function T = generate_country_csv_from_hopkins_Remake(country, reportsDir, varargin)
% GENERATE_COUNTRY_CSV_FROM_HOPKINS  (R2019a-compatible)
% Build a clean daily CSV for one country from the Johns Hopkins (JHU CSSE)
% daily-report files, ready for the fitting pipeline.
%
% Output columns (same as the Greek CSV the pipeline expects):
%   Date            calendar day
%   Confirmed_Cases daily NEW confirmed cases
%   Total_Cases     cumulative confirmed cases
%
% Handles: provinces (summed), schema changes across 2020-2023, missing
% report days (filled by the daily grid), and weekly-only reporting
% (linear interpolation of the cumulative series -> smooth daily new cases,
% no zeros/spikes; works for gaps of any length).
%
% Usage:
%   generate_country_csv_from_hopkins('Finland', reportsDir)
%   generate_country_csv_from_hopkins('Greece',  reportsDir, ...
%        'StartDate','2020-04-09','EndDate','2022-07-10')

%% ---- options ----
p = inputParser;
p.addRequired('country');
p.addRequired('reportsDir');
p.addParameter('StartDate', '');
p.addParameter('EndDate',   '');
p.addParameter('OutFile',   '');
p.addParameter('EnforceMonotone', true, @islogical);
p.parse(country, reportsDir, varargin{:});
opt = p.Results;
country = strtrim(string(country));

files = dir(fullfile(reportsDir, '*.csv'));
assert(~isempty(files), 'No CSV files found in %s', reportsDir);

% silence the "modified variable names" warning that R2019a prints per file
wstate = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
cleaner = onCleanup(@() warning(wstate));

dates    = NaT(numel(files),1);
cumConf  = nan(numel(files),1);
readCount = 0;   % how many report files were actually read

%% ---- read every daily report ----
for i = 1:numel(files)
    [~, name] = fileparts(files(i).name);

    % date from filename MM-DD-YYYY.csv (skip non-date files, e.g. README)
    try
        d = datetime(name, 'InputFormat','MM-dd-yyyy', 'Locale','en_US');
    catch
        continue;
    end
    if isnat(d), continue; end

    % read the report (R2019a default naming; no VariableNamingRule)
    try
        Tr = readtable(fullfile(files(i).folder, files(i).name));
    catch
        continue;
    end
    readCount = readCount + 1;

    vn = string(Tr.Properties.VariableNames);
    cCol = vn(contains(lower(vn),'country'));
    kCol = vn(strcmpi(vn,'Confirmed'));
    if isempty(cCol) || isempty(kCol), continue; end
    cCol = cCol(1); kCol = kCol(1);

    ctry = strtrim(string(Tr.(cCol)));
    conf = Tr.(kCol);
    if ~isnumeric(conf), conf = str2double(string(conf)); end

    mask = strcmpi(ctry, country);
    if any(mask)
        dates(i)   = d;
        cumConf(i) = sum(conf(mask), 'omitnan');   % sum over provinces
    end
end

% distinguish the two failure modes clearly
assert(readCount > 0, ...
    'Could not read any report file. Check the folder path and MATLAB version.');

%% ---- assemble and clean ----
ok = ~isnat(dates) & ~isnan(cumConf);
assert(any(ok), ...
    'Country "%s" not found in %d report files. Check the exact JHU spelling (e.g. "US", "Korea, South").', ...
    country, readCount);

dates = dates(ok);  cumConf = cumConf(ok);
[dates, ix] = sort(dates);  cumConf = cumConf(ix);

[udates, ~, g] = unique(dates);
cumConf = accumarray(g, cumConf, [], @max);
dates   = udates;

if opt.EnforceMonotone
    cumConf = cummax(cumConf);
end

%% ---- daily grid + linear interpolation of the cumulative ----
grid = (dates(1):caldays(1):dates(end))';
isUpd = [true; diff(cumConf) > 0];
Total_Cases = interp1(datenum(dates(isUpd)), cumConf(isUpd), datenum(grid), 'linear');

keep = ~isnan(Total_Cases);
grid = grid(keep);  Total_Cases = Total_Cases(keep);

Confirmed_Cases = [Total_Cases(1); diff(Total_Cases)];
Confirmed_Cases = max(Confirmed_Cases, 0);

T = table(grid, Confirmed_Cases, Total_Cases, ...
    'VariableNames', {'Date','Confirmed_Cases','Total_Cases'});

%% ---- optional clipping ----
if ~isempty(opt.StartDate), T = T(T.Date >= datetime(opt.StartDate),:); end
if ~isempty(opt.EndDate),   T = T(T.Date <= datetime(opt.EndDate),:);   end

%% ---- write (Date like the Greek file: dd-mmm-yyyy 00:00:00) ----
if isempty(opt.OutFile)
    opt.OutFile = sprintf('Covid_19_%s_Hopkins_rmk.csv', country);
end
Tout = T;
Tout.Date = string(datestr(T.Date, 'dd-mmm-yyyy HH:MM:SS'));
writetable(Tout, opt.OutFile);

fprintf('%s: %d days, %s to %s  ->  %s\n', country, height(T), ...
    datestr(T.Date(1),'yyyy-mm-dd'), datestr(T.Date(end),'yyyy-mm-dd'), opt.OutFile);
end