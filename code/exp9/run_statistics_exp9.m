[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
timing_analysis_file = fullfile(file_path, sprintf('../../results/%s/timing_analysis.mat', exp_name)); % get path to results
load(timing_analysis_file)

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
within_table = table(categorical(group_info(:,1)), categorical(group_info(:,2)), 'VariableNames', {'SNR', 'Time'});
rm = fitrm(datatable,'Var1-Var20~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'SNR*Time'); % C is used to handle epsilon correction for violations of sphericity
report_main_effect(rm, C{2}, ranovatbl, 'SNR')
report_main_effect(rm, C{3}, ranovatbl, 'Time')
report_interaction(rm, C{4}, ranovatbl, 'SNR', 'Time')


% statistics comparing exp9 to exp7
d_primes = (d_primes(1:10,:) + d_primes(11:20,:)) / 2;
comparison_exp = 'exp7';
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp)); % get path to results
load(results_file)

d_primes_comparison = nan(size(results(1).d_primes,1),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_comparison(:, i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes_comparison = d_primes_comparison(:, logical(keep_subs));
n_subs_comparison = size(d_primes_comparison,2);

datatable = array2table(horzcat(d_primes, d_primes_comparison)'); 
datatable.BGType = categorical(vertcat(ones(n_subs,1), 2*ones(n_subs_comparison,1)));
group_info = (1:10)'; % foreground onset times
within_table = table(categorical(group_info), 'VariableNames', {'Time'});
rm = fitrm(datatable,'Var1-Var10~BGType','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'Time'); % C is used to handle epsilon correction for violations of sphericity

fprintf('Statistics for %s vs %s\n', exp_name, comparison_exp)
varname = 'BGType'; % between subjects only
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

varname = 'BGType:Time'; % between-within interaction
row_index = find(strcmp(ranovatbl.Properties.RowNames, sprintf('%s', varname)));
F = ranovatbl(row_index,:).F;
df1 = ranovatbl(row_index,:).DF;
df2 = ranovatbl(row_index+1,:).DF;
pvalue = ranovatbl(row_index,:).pValue;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.4f', pvalue);
end
partial_eta_squared = ranovatbl(row_index,:).SumSq / (ranovatbl(row_index,:).SumSq + ranovatbl(row_index+1,:).SumSq);    
fprintf('Interaction Effect of %s\n', varname)
fprintf('F(%.2f,%.2f)=%.2f, %s\n', df1, df2, F, pvalue_str)
fprintf('Partial Eta Squared: %.2f\n\n', partial_eta_squared)

permutation_file = fullfile(file_path, sprintf('../../results/%s/permutation_test.mat', exp_name)); % get path to results
load(permutation_file)
p_timescale = mean(abs(timescale_nullstats) >= abs(timescale_statistic));
p_magnitude = mean(abs(magnitude_nullstats) >= abs(magnitude_statistic));

pvalue = p_timescale;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
fprintf('Timescale Permutation Test: %s\n', pvalue_str)
fprintf('Timescale difference: %.0f ms\n\n', timescale_statistic)

pvalue = p_magnitude;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
fprintf('Magnitude Permutation Test: %s\n', pvalue_str)
fprintf('Magnitude difference: %.2f\n\n', magnitude_statistic)

disp('Matched statistics')

permutation_file = fullfile(file_path, sprintf('../../results/%s/matched_permutation_test.mat', exp_name)); % get path to results
load(permutation_file)
p_timescale = mean(abs(timescale_nullstats) >= abs(timescale_statistic));
p_magnitude = mean(abs(magnitude_nullstats) >= abs(magnitude_statistic));

pvalue = p_timescale;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
fprintf('Timescale Permutation Test: %s\n', pvalue_str)
fprintf('Timescale difference: %.0f ms\n\n', timescale_statistic)

pvalue = p_magnitude;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
fprintf('Magnitude Permutation Test: %s\n', pvalue_str)
fprintf('Magnitude difference: %.2f\n\n', magnitude_statistic)

x = linspace(250,2500,10);
y = mean(d_primes,2)';
est_params_avg = fit_elbow_func(x,y);
actual_elbow_point = est_params_avg(3);
elbow_points = est_params(:,3);

fprintf('Elbow point: %.0f ms, 95%% CI: [%.0f, %.0f] ms\n', actual_elbow_point, prctile(elbow_points, 2.5), prctile(elbow_points, 97.5))