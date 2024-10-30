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

datatable = array2table(d_primes');
group_info = results(1).group_info;
within_table = table(categorical(group_info(:,2)), 'VariableNames', {'Time'});
rm = fitrm(datatable,'Var1-Var10~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'Time'); % C is used to handle epsilon correction for violations of sphericity

report_main_effect(rm, C{2}, ranovatbl, 'Time')

% re-analyze data by splitting into first half versus second half of block
temp_block_position = repmat(1:2, 20, 1);
block_position = repmat(temp_block_position(:), 8, 1);
d_primes = nan(20,length(results));
for i = 1:length(results)
    responses = results(i).responses;
    responses_to_use = responses(:, 2:end);    
    responses_to_use(:, 2) = block_position;
    responses_to_use(sum(responses,2)==0,:) = []; % remove trials with no data                            

    [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
    group_info = [ID1 ID2];
    hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
    mean_false_alarm_rates = splitapply(@compute_false_alarm_rate, responses_to_use(:, 4), responses_to_use(:, 5), responses_to_use(:, 2)); % compute false alarms over all "conditions" since they only apply with foregrounds present
    temp_false_alarm_rates = repmat(mean_false_alarm_rates, 1, 10)';
    d_primes(:, i) = norminv(hit_rates) - norminv(temp_false_alarm_rates(:));
end

d_primes = d_primes(:, logical(keep_subs));

datatable = array2table(d_primes');
within_table = table(categorical(group_info(:,1)), categorical(group_info(:,2)), 'VariableNames', {'BlockPosition', 'Time'});
rm = fitrm(datatable,'Var1-Var20~1','WithinDesign',within_table);
[ranovatbl, ~, C, ~] = ranova(rm, 'WithinModel', 'BlockPosition*Time'); % C is used to handle epsilon correction for violations of sphericity

report_main_effect(rm, C{2}, ranovatbl, 'BlockPosition')