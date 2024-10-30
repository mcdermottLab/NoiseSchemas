[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
data_path = fullfile(file_path, sprintf('../../results/%s/raw/', exp_name)); % get path to raw data
save_path = fullfile(file_path, sprintf('../../results/%s/', exp_name)); % get path to save results

files = dir(fullfile(data_path, 'v*_subject_*.csv')); % get list of files
results = struct; % create an empty struct to store data in
azimuths = linspace(-90, 90, 19);
rows = ['ABCDEFG'];
elevations = linspace(40, -20, 7);
n_trials = 160;
threshold = 30; % performance threshold for average azimuth error
demographic_info = readtable(fullfile(data_path, 'HiN_localization_demographic_info.csv'));

for f = 1:length(files)
    experiment_data = readtable(fullfile(files(f).folder, files(f).name));
    azimuth_error = nan(n_trials, 1);
    elevation_error = nan(n_trials, 1);
    for t = 1:n_trials
        correct_azimuth = azimuths(str2num(experiment_data.correct_response{t}(2:end)));
        subject_azimuth = azimuths(str2num(experiment_data.subject_response{t}(2:end)));
        azimuth_error(t) = abs(correct_azimuth - subject_azimuth);        
        correct_elevation = elevations(experiment_data.correct_response{t}(1) == rows);
        subject_elevation = elevations(experiment_data.subject_response{t}(1) == rows);
        elevation_error(t) = abs(correct_elevation - subject_elevation);          
    end
    results(f).responses = experiment_data.subject_response(1:n_trials);
    results(f).correct_responses = experiment_data.correct_response(1:n_trials);
    results(f).onset_times = experiment_data.onset_time(1:n_trials);
    results(f).azimuth_errors = azimuth_error;
    results(f).elevation_errors = elevation_error;
    results(f).angular_distances = experiment_data.angular_distance(1:n_trials);
    results(f).mean_azimuth_error = splitapply(@mean, azimuth_error, experiment_data.onset_time(1:n_trials));
    results(f).mean_elevation_error = splitapply(@mean, elevation_error, experiment_data.onset_time(1:n_trials));
    results(f).mean_angular_distance = splitapply(@mean, experiment_data.angular_distance(1:n_trials), experiment_data.onset_time(1:n_trials));
    
    [~, subjectID] = fileparts(files(f).name);
    subject_ind = find(ismember(demographic_info.SubjectID, subjectID));
    if strcmp(demographic_info.HaveYouEverKnownHowToPlayAMusicalInstrument_{subject_ind}, 'Yes')
        musician = 'musician';
    else
        musician = 'nonmusician';
    end
    results(f).demographics = {demographic_info.WhatIsYourGender_{subject_ind}, demographic_info.WhatIsYourAge_(subject_ind), musician}; % store demographic info
    if mean(results(f).mean_azimuth_error) >= threshold
        results(f).pass = false;
    else
        results(f).pass = true;
    end
end

save_name = fullfile(save_path, 'results.mat');
save(save_name, 'results')