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

% statistics for exp7 effect of time only
datatable = array2table(d_primes');
group_info = results(1).group_info;
within_table = table(categorical(group_info(:,2)), 'VariableNames', {'Time'});
rm = fitrm(datatable,'Var1-Var10~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'Time'); % C is used to handle epsilon correction for violations of sphericity

fprintf('Statistics for %s only\n', exp_name)
report_main_effect(rm, C{2}, ranovatbl, 'Time')

% statistics comparing exp7 to exp1
comparison_exp = 'exp1';
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', comparison_exp)); % get path to results
load(results_file)

d_primes_comparison = nan(size(results(1).d_primes,1),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes_comparison(:, i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes_comparison = d_primes_comparison(:, logical(keep_subs));
d_primes_comparison = (d_primes_comparison(1:10,:) + d_primes_comparison(11:20,:)) / 2;
n_subs_comparison = size(d_primes_comparison,2);

datatable = array2table(horzcat(d_primes, d_primes_comparison)'); 
datatable.BGFGPairing = categorical(vertcat(ones(n_subs,1), 2*ones(n_subs_comparison,1)));
group_info = (1:10)'; % foreground onset times
within_table = table(categorical(group_info), 'VariableNames', {'Time'});
rm = fitrm(datatable,'Var1-Var10~BGFGPairing','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'Time'); % C is used to handle epsilon correction for violations of sphericity

fprintf('Statistics for %s vs %s\n', exp_name, comparison_exp)
varname = 'BGFGPairing'; % between subjects only
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

varname = 'BGFGPairing:Time'; % between-within interaction
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
fprintf('Interaction Effect of %s\n', varname)
fprintf('F(%.2f,%.2f)=%.2f, %s\n', df1, df2, F, pvalue_str)
fprintf('Partial Eta Squared: %.2f\n\n', partial_eta_squared)

[r, pvalue] = corr(mean(d_primes, 2), mean(d_primes_comparison,2));
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
fprintf('Exp 7 vs Exp 1 Correlation\n')
fprintf('R^2=%.2f, %s\n\n', r^2, pvalue_str)

permutation_file = fullfile(file_path, sprintf('../../results/%s/permutation_test_%s_vs_%s.mat', exp_name, exp_name, comparison_exp)); % get path to results
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

% statistics comparing exp7 to exp6
comparison_exp = 'exp6';
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
datatable.Repetition = categorical(vertcat(ones(n_subs,1), 2*ones(n_subs_comparison,1)));
group_info = (1:10)'; % foreground onset times
within_table = table(categorical(group_info), 'VariableNames', {'Time'});
rm = fitrm(datatable,'Var1-Var10~Repetition','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'Time'); % C is used to handle epsilon correction for violations of sphericity

fprintf('Statistics for %s vs %s\n', exp_name, comparison_exp)
varname = 'Time'; 
row_index = find(strcmp(ranovatbl.Properties.RowNames, sprintf('(Intercept):%s', varname)));
F = ranovatbl(row_index,:).F;
df1 = ranovatbl(row_index,:).DF;
df2 = ranovatbl(row_index+2,:).DF;
pvalue = ranovatbl(row_index,:).pValue;
if pvalue < 0.001
    pvalue_str = 'p<0.001';
else
    pvalue_str = sprintf('p=%.2f', pvalue);
end
partial_eta_squared = ranovatbl(row_index,:).SumSq / (ranovatbl(row_index,:).SumSq + ranovatbl(row_index+2,:).SumSq);
    
fprintf('Main Effect of %s\n', varname)
fprintf('F(%.2f,%.2f)=%.2f, %s\n', df1, df2, F, pvalue_str)
fprintf('Partial Eta Squared: %.2f\n\n', partial_eta_squared)

varname = 'Repetition'; % between subjects only
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

varname = 'Repetition:Time'; % between subjects only
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
    
fprintf('Interaction Effect of %s\n', varname)
fprintf('F(%.2f,%.2f)=%.2f, %s\n', df1, df2, F, pvalue_str)
fprintf('Partial Eta Squared: %.2f\n\n', partial_eta_squared)

permutation_file = fullfile(file_path, sprintf('../../results/%s/permutation_test_%s_vs_%s.mat', exp_name, exp_name, comparison_exp)); % get path to results
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

fprintf('Matched statistics for %s vs %s\n', exp_name, comparison_exp)
permutation_file = fullfile(file_path, sprintf('../../results/%s/matched_permutation_test_%s_vs_%s.mat', exp_name, exp_name, comparison_exp)); % get path to results
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

pairings_file = fullfile(file_path, sprintf('../../results/%s/spectrotemporal_split_pairings.mat', exp_name)); % get path to timing analysis
load(pairings_file)
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)

low_pairs(:,2) = goodinds(low_pairs(:,2));
low_pairs(:, 3:4) = [d_spec_low, d_spectemp_low];
high_pairs(:,2) = goodinds(high_pairs(:,2));
high_pairs(:, 3:4) = [d_spec_high, d_spectemp_high];

all_groups = [];
for i = 1:3
    for j = 1:10
        all_groups(end+1, :) = [i j];
    end
end

d_primes = nan(30, length(results));
distances = nan(4, length(results));
remove_subs = [];
for i = 1:length(results)
    if results(i).pass
        responses = results(i).responses;
        on_trials = find(responses(:,5)==2);
        low_trials = on_trials(ismember(responses(on_trials,1:2), low_pairs(:,1:2), 'rows'));
        responses(low_trials,3) = 3;
        high_trials = on_trials(ismember(responses(on_trials,1:2), high_pairs(:,1:2), 'rows'));
        responses(high_trials,3) = 2;
        responses_to_use = responses(:, 2:end);
        responses_to_use(sum(responses_to_use,2)==0,:) = []; % remove trials with no data                           
        [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
        group_info = [ID1 ID2];
        hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
        mean_false_alarm_rate = compute_false_alarm_rate(responses_to_use(:, 4), responses_to_use(:, 5)); % compute false alarms over all "conditions" since they only apply with foregrounds present
        inds = find(ismember(all_groups, [ID1 ID2], 'rows'));
        d_primes(inds, i) = norminv(hit_rates) - norminv(mean_false_alarm_rate);
        
        low_inds = find(ismember(low_pairs(:, 1:2), responses(low_trials,1:2), 'rows'));
        high_inds = find(ismember(high_pairs(:, 1:2), responses(high_trials,1:2), 'rows'));
        distances(1, i) = mean(low_pairs(low_inds, 3));        
        distances(2, i) = mean(low_pairs(low_inds, 4));
        distances(3, i) = mean(high_pairs(high_inds, 3));
        distances(4, i) = mean(high_pairs(high_inds, 4));
    else
        remove_subs(end+1) = i;
    end
end

d_primes = d_primes(11:30, :);
d_primes(:, remove_subs) = [];
n_subs = size(d_primes,2);

datatable = array2table(d_primes');
group_info = vertcat(results(1).group_info,results(1).group_info);
group_info(11:20,1) = 2;
within_table = table(categorical(group_info(:,1)), categorical(group_info(:,2)), 'VariableNames', {'STDist', 'Time'});
rm = fitrm(datatable,'Var1-Var20~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'STDist*Time'); % C is used to handle epsilon correction for violations of sphericity
report_main_effect(rm, C{2}, ranovatbl, 'STDist')
report_main_effect(rm, C{3}, ranovatbl, 'Time')
report_interaction(rm, C{4}, ranovatbl, 'STDist', 'Time')

permutation_file = fullfile(file_path, sprintf('../../results/%s/spectrotemporal_split_permutation_test.mat', exp_name)); % get path to results
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