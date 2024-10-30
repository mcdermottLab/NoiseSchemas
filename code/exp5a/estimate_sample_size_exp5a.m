[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
data_path = fullfile(file_path, sprintf('../../results/%s/raw_pilot/', exp_name)); % get path to raw pilot data

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

files = dir(fullfile(data_path, 'Batch*.csv')); % get list of files
n_trials = 320; % exp1 has 320 trials
threshold = 0.8; % performance threshold for average d'
results = struct; % create an empty struct to store data in
k = 0; % counter for number of subjects with usable data
for f = 1:length(files)
    batch_file_name = fullfile(data_path, files(f).name);   
    batch_info = readtable(batch_file_name,'Delimiter',',', 'HeaderLines',0); % reads in the batch csv from mturk
    for sub = 1:size(batch_info,1) % loop through all subjects
        if isnan(batch_info.Answer_BonusAmount(sub)) % if they don't have any bonus, then they failed the headphone check and didn't do the experiment OR didn't finish it
        elseif batch_info.Answer_BonusAmount(sub)==0 % if they don't have any bonus, then they failed the headphone check and didn't do the experiment OR didn't finish it
        elseif strcmp(batch_info.Answer_Q6hearingloss{sub}, 'yesHL') % do not include if they self-report hearing loss
        else
            k = k + 1; % increase counter since subject has usable data
            responses = zeros(n_trials,5);  % preallocate an array
            for trial = 1:n_trials % loop through all the trials
                    stimID_var = sprintf('Answer_stimID_for_expt_1trial_%d',trial); % store the variable name to index the table for the given trial
                    responses(trial,1:4) = str2num(strrep(batch_info{sub,stimID_var}{:}, '_', ' ')); % extract 4 stimulus parameters (see below)
                    resp_var = sprintf('Answer_Resp_for_expt_1trial_%d',trial); % store the variable name to index the table for the given trial                    
                    responses(trial,5) = batch_info{sub,resp_var}; % extract subject response
            end
            results(k).responses = responses; % store stimulus info and subject responses organized along columns as: | bg/fg paring | SNR | fg onset time | fg on/off | subject response |             
            results(k).demographics = {batch_info.Answer_Q1Gender{sub}, batch_info.Answer_Q2age(sub), batch_info.Answer_Q7music{sub}}; % store demographic info
            results(k).overall_percent_correct = mean(responses(:,4) == responses(:,5)); % store overall performance

            [G, ID1, ID2] = findgroups(responses(:,2),responses(:,3)); % get the group info
            group_info = [ID1 ID2];
            hit_rates = splitapply(@compute_hit_rate, responses(:, 4), responses(:, 5), G);  % get hit rates and false alarm rates for all conditions (interrupter type by onset time)      
            mean_false_alarm_rate = splitapply(@compute_false_alarm_rate, responses(:, 4), responses(:, 5), responses(:,2)); % compute false alarms over interrupter types

            results(k).hit_rates = hit_rates; % store hit rates
            results(k).mean_false_alarm_rate = mean_false_alarm_rate; % store mean false alarm rate
            results(k).d_primes = vertcat(norminv(hit_rates(1:8)) - norminv(mean_false_alarm_rate(1)), ...
                                          norminv(hit_rates(9:16)) - norminv(mean_false_alarm_rate(2))); % store d'
            results(k).overall_performance = mean(results(k).d_primes); % store average d'
            if results(k).overall_performance >= threshold
                results(k).pass = true;
            else
                results(k).pass = false;
            end
            results(k).group_info = group_info; % store grouping info
        end
    end 
end

d_primes = nan(length(results(1).d_primes),length(results));
keep_subs = nan(1,length(results));
for i = 1:length(results)    
    d_primes(:,i) = results(i).d_primes;
    keep_subs(i) = results(i).pass;
end

d_primes = d_primes(:, logical(keep_subs));
n_subs = size(d_primes,2);

seed = 1;
rng(seed)
n_iters = 10000;
step_size = 2;
min_subs = 10;
max_subs = floor(n_subs/step_size) * step_size;
sample_sizes = min_subs:step_size:max_subs;
shr_means = nan(length(sample_sizes), 1);
for s = 1:length(sample_sizes)
    split_size = round(sample_sizes(s)/2);    
    split_half_reliabilites = nan(n_iters,1);
    for n = 1:n_iters
        random_permutation = randperm(n_subs);
        split1 = random_permutation(1:split_size);
        split2 = random_permutation((1+split_size):(2*split_size));
        mean_d_primes_split1 = mean(d_primes(:, split1),2);
        mean_d_primes_split2 = mean(d_primes(:, split2),2);
        r = corr(mean_d_primes_split1([1:4, 9:12]), mean_d_primes_split2([1:4, 9:12]));
        split_half_reliabilites(n) = spearman_brown(r);
    end
    shr_means(s) = tanh(mean(atanh(split_half_reliabilites)));
end

target_sample_size = extrapolate_reliability(sample_sizes', ...
                                              shr_means, ...
                                              0.90, ...
                                              100, ... 
                                              [-10 10 1] ...
                                              );

fprintf('number of subjects: %d\n', n_subs)
fprintf('target sample size: %d\n', target_sample_size)
fprintf('min reliability: %.2f (n=%d)\n', shr_means(1), min_subs)
fprintf('max reliability: %.2f (n=%d)\n', shr_means(end), max_subs)