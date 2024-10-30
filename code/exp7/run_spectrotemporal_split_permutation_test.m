[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)

permutation_test_file = fullfile(file_path, sprintf('../../results/%s/spectrotemporal_split_permutation_test.mat', exp_name)); % get path to save permutation test results
pairings_file = fullfile(file_path, sprintf('../../results/%s/spectrotemporal_split_pairings.mat', exp_name)); % get path to save permutation test results
load(pairings_file)

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

low_pairs(:,2) = goodinds(low_pairs(:,2));
low_pairs(:, 3:4) = [d_spec_low, d_spectemp_low];
high_pairs(:,2) = goodinds(high_pairs(:,2));
high_pairs(:, 3:4) = [d_spec_high, d_spectemp_high];
% low pairs contain fg-bg pairs with low spectrotemporal distance
% high pairs contain fg-bg pairs with high spectrotemporal distance
% columns are: [bg index, fg index, spectral distance, spectrotemporal distance]

all_groups = [];
for i = 1:3 % 1==trials not included here, 2==high spectemp distance trials, 3==low spectemp distance trials
    for j = 1:10 % time condition
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

d_primes = d_primes(11:30, :); % remove data for pairings not considered here
d_primes(:, remove_subs) = [];
n_subs = size(d_primes,2);

d_primes_high = d_primes(1:10, :);
d_primes_low = d_primes(11:20, :);

% fit elbow function to actual data
x = linspace(250,2500,10);
x_test = 250:2500;

y = mean(d_primes_high,2,'omitnan')';
est_params = fit_elbow_func(x,y);
y_pred = elbow_func(x_test, est_params);

d_primes_comparison = d_primes_low;
y_comparison = mean(d_primes_comparison,2,'omitnan')';
est_params_comparison = fit_elbow_func(x,y_comparison);
y_pred_comparison = elbow_func(x_test, est_params_comparison);

timescale_statistic = est_params(3) - est_params_comparison(3);
magnitude_statistic = (y_pred(end)-y_pred(1)) - (y_pred_comparison(end)-y_pred_comparison(1));

rng(1)
niters = 10000;
magnitude_nullstats = nan(niters, 1);
timescale_nullstats = nan(niters, 1);

distcomp.feature( 'LocalUseMpiexec', false ); % added to help the pool work...
num_workers = 20;
parallel_pool = parpool(num_workers);

tic
parfor n = 1:niters
    rng(n)
    if mod(n,100)==0
        fprintf('Starting iteration %d\n', n)
    end
    d_primes_perm = nan(length(results(1).d_primes),length(results));
    d_primes_perm = nan(30, length(results));
    remove_subs = [];
    for i = 1:length(results)
        if results(i).pass
            responses = results(i).responses;
            on_trials = find(responses(:,5)==2);
            low_trials = on_trials(ismember(responses(on_trials,1:2), low_pairs(:,1:2), 'rows'));
            high_trials = on_trials(ismember(responses(on_trials,1:2), high_pairs(:,1:2), 'rows'));
            all_trials = vertcat(low_trials,high_trials);

            rp = randperm(length(all_trials));
            random_low_trials = all_trials(rp(1:length(low_trials)));
            random_high_trials = all_trials(rp((1+length(low_trials):end)));
            responses(random_low_trials,3) = 3;
            responses(random_high_trials,3) = 2;
    
            responses_to_use = responses(:, 2:end);
            responses_to_use(sum(responses_to_use,2)==0,:) = []; % remove trials with no data  
    
            [G, ID1, ID2] = findgroups(responses_to_use(:,2),responses_to_use(:,3)); % get the group info
            group_info = [ID1 ID2];
            hit_rates = splitapply(@compute_hit_rate, responses_to_use(:, 4), responses_to_use(:, 5), G);  % get hit rates and false alarm rates for all conditions (SNR by onset time)      
            mean_false_alarm_rate = compute_false_alarm_rate(responses_to_use(:, 4), responses_to_use(:, 5)); % compute false alarms over all "conditions" since they only apply with foregrounds present
            inds = find(ismember(all_groups, [ID1 ID2], 'rows'));
            d_primes_perm(inds, i) = norminv(hit_rates) - norminv(mean_false_alarm_rate);
        else
            remove_subs(end+1) = i;
        end
    end

    d_primes_perm = d_primes_perm(11:30, :);
    d_primes_perm(:, remove_subs) = [];

    perm_dp = d_primes_perm(1:10, :);
    perm_y = mean(perm_dp,2,'omitnan')';
    perm_est_params = fit_elbow_func(x,perm_y);
    perm_y_pred = elbow_func(x_test, perm_est_params);

    perm_dp_comparison = d_primes_perm(11:20, :);
    perm_y_comparison = mean(perm_dp_comparison,2,'omitnan')';
    perm_est_params_comparison = fit_elbow_func(x,perm_y_comparison);
    perm_y_pred_comparison = elbow_func(x_test, perm_est_params_comparison);
    
    timescale_nullstats(n) = perm_est_params(3) - perm_est_params_comparison(3);
    magnitude_nullstats(n) = (perm_y_pred(end)-perm_y_pred(1)) - (perm_y_pred_comparison(end)-perm_y_pred_comparison(1));
end
toc

save(permutation_test_file, 'magnitude_nullstats', 'magnitude_statistic', 'timescale_nullstats', 'timescale_statistic')
delete(parallel_pool)
