[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

datatable = array2table(d_primes');
group_info = results(1).group_info;
within_table = table(categorical(group_info(:,1)), categorical((group_info(:,2)>4) + 1), 'VariableNames', {'NoiseType', 'PosRelNoise'});
rm = fitrm(datatable,'Var1-Var16~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'NoiseType*PosRelNoise'); % C is used to handle epsilon correction for violations of sphericity
report_main_effect(rm, C{2}, ranovatbl, 'NoiseType')
report_main_effect(rm, C{3}, ranovatbl, 'PosRelNoise')
report_interaction(rm, C{4}, ranovatbl, 'NoiseType', 'PosRelNoise')