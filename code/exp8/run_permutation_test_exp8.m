[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
permutation_test_file = fullfile(file_path, sprintf('../../results/%s/permutation_test.mat', exp_name)); % get path to save permutation test results

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

d_primes_rep = d_primes(1:10, :);
d_primes_nonrep = d_primes(11:20, :);

% fit elbow function to actual data
x = linspace(250,2500,10);
x_test = 250:2500;

y = mean(d_primes_rep,2)';
est_params = fit_elbow_func(x,y);
y_pred = elbow_func(x_test, est_params);

d_primes_comparison = d_primes_nonrep;
y_comparison = mean(d_primes_comparison,2)';
est_params_comparison = fit_elbow_func(x,y_comparison);
y_pred_comparison = elbow_func(x_test, est_params_comparison);

magnitude_statistic = (y_pred(end)-y_pred(1)) - (y_pred_comparison(end)-y_pred_comparison(1));

rng(1)
niters = 10000;
magnitude_nullstats = nan(niters, 1);

distcomp.feature( 'LocalUseMpiexec', false ); % added to help the pool work...
num_workers = 20;
parallel_pool = parpool(num_workers);

n_trials = 320; % exp has 320 trials
n_conditions = 10;
repcode = repmat(1:2, 1, n_trials/2)'; % 1==repeated bg and 2==non-repeated bg

tic
parfor n = 1:niters
    d_primes_perm = nan(length(results(1).d_primes),length(results));
    for i = 1:length(results)
        responses = results(i).responses;
        responses_to_use = responses(:, [1 3:6]); % get relevant stimulus/condition info   
        for j = 1:n_conditions % number of onset times
            trial_inds = find(responses_to_use(:,3)==j & responses_to_use(:,4)==2);
            temp_repcode = repcode(trial_inds);            
            responses_to_use(trial_inds, 2) = temp_repcode(randperm(length(temp_repcode)));
        end
        trial_inds = find(responses_to_use(:,3)==1 & responses_to_use(:,4)==1);
        temp_repcode = repcode(trial_inds);    
        responses_to_use(trial_inds, 2) = temp_repcode(randperm(length(temp_repcode)));
        responses_to_use(sum(responses,2)==0,:) = []; % remove trials with no data                                            
        [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
        group_info = [ID1 ID2];
        hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
        mean_false_alarm_rates = splitapply(@compute_false_alarm_rate, responses_to_use(:, 4), responses_to_use(:, 5), responses_to_use(:, 2)); % compute false alarms over all "conditions" since they only apply with foregrounds present
        temp_false_alarm_rates = repmat(mean_false_alarm_rates, 1, 10)';
        d_primes_perm(:, i) = norminv(hit_rates) - norminv(temp_false_alarm_rates(:));
    end
    d_primes_perm = d_primes_perm(:, logical(keep_subs));

    perm_dp = d_primes_perm(1:10, :);
    perm_y = mean(perm_dp,2)';
    perm_est_params = fit_elbow_func(x,perm_y);
    perm_y_pred = elbow_func(x_test, perm_est_params);

    perm_dp_comparison = d_primes_perm(11:20, :);
    perm_y_comparison = mean(perm_dp_comparison,2)';
    perm_est_params_comparison = fit_elbow_func(x,perm_y_comparison);
    perm_y_pred_comparison = elbow_func(x_test, perm_est_params_comparison);

    magnitude_nullstats(n) = (perm_y_pred(end)-perm_y_pred(1)) - (perm_y_pred_comparison(end)-perm_y_pred_comparison(1));
end
toc

save(permutation_test_file, 'magnitude_nullstats', 'magnitude_statistic')
delete(parallel_pool)