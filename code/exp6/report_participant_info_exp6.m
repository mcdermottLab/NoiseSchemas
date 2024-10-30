[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
data_path = fullfile(file_path, sprintf('../../results/%s/raw/', exp_name)); % get path to raw data
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)

files = dir(fullfile(data_path, 'sub*.csv')); % get list of files
total_participants = length(files);

num_pass = 0;
ages = [];
sex = {};
for i = 1:length(results)
    pass = results(i).pass;
    num_pass = num_pass + pass;
    if pass
        if iscell(results(i).demographics{2})
            ages(end+1) = NaN;
        else
            ages(end+1) = results(i).demographics{2};
        end
        sex{end+1} = results(i).demographics{1};
    end
end

n_female = sum(ismember(sex, 'Female'));
n_male = sum(ismember(sex, 'Male'));
n_nonbinary = sum(ismember(sex, 'Nonbinary'));

ages(ages > 99) = NaN;

clc
fprintf('=== %s ===\n', exp_name)
fprintf('total: %d\n', total_participants)
fprintf('HC/HL excluded: %d\n', total_participants - length(results))
fprintf('passed HC: %d\n', length(results))
fprintf('failed thresh: %d\n', length(results) - num_pass)
fprintf('passed thresh: %d\n', num_pass)
fprintf('females: %d\n', n_female)
fprintf('males: %d\n', n_male)
fprintf('nonbinary: %d\n', n_nonbinary)
fprintf('unknown: %d\n', num_pass - (n_female + n_male + n_nonbinary))
fprintf('mean age: %.1f\n', mean(ages, 'omitnan'))
fprintf('median age: %.1f\n', median(ages, 'omitnan'))
fprintf('stdev age: %.1f\n', std(ages, 'omitnan'))
