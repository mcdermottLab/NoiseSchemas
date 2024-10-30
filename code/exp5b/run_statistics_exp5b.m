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

prev_exp_name = strcat(exp_name(1:end-1), 'a');
prev_results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', prev_exp_name)); % get path to results
load(prev_results_file)

d_primes_prev = nan(length(results(1).d_primes),length(results));
keep_subs_prev = nan(1,length(results));
for i = 1:length(results)    
    d_primes_prev(:,i) = results(i).d_primes;
    keep_subs_prev(i) = results(i).pass;
end

d_primes_prev = d_primes_prev(:, logical(keep_subs_prev));
n_subs_prev = size(d_primes_prev,2);

datatable = array2table(horzcat(d_primes([[5:8, 13:16]],:), d_primes_prev([5:8, 13:16], :))'); % only comparing timepoints AFTER interruption
datatable.InterrupterDuration = categorical(vertcat(ones(n_subs,1), 2*ones(n_subs_prev,1)));
group_info = results(1).group_info([5:8, 13:16],:);
within_table = table(categorical(group_info(:,1)), categorical(group_info(:,2)-4), 'VariableNames', {'NoiseType', 'Time'});
rm = fitrm(datatable,'Var1-Var8~InterrupterDuration','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'NoiseType*Time'); % C is used to handle epsilon correction for violations of sphericity

varname = 'InterrupterDuration'; % between subjects only
row_index = find(strcmp(ranovatbl.Properties.RowNames, sprintf('%s', varname)));
F = ranovatbl(row_index,:).F;
df1 = ranovatbl(row_index,:).DF;
df2 = ranovatbl(row_index+1,:).DF;
pvalue = ranovatbl(row_index,:).pValue;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
partial_eta_squared = ranovatbl(row_index,:).SumSq / (ranovatbl(row_index,:).SumSq + ranovatbl(row_index+1,:).SumSq);
    
fprintf('Main Effect of %s\n', varname)
fprintf('F(%.2f,%.2f)=%.2f, %s\n', df1, df2, F, pvalue_str)
fprintf('Partial Eta Squared: %.2f\n\n', partial_eta_squared)
