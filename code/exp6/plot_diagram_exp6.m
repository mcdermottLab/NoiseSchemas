%% initialization and settings
[file_path, file_name] = fileparts(mfilename('fullpath')); % get full path to this file
path_parts = strsplit(file_path,filesep); % split path into parts
exp_name = path_parts{end}; % extract experiment name
wav_path = fullfile(file_path, sprintf('../../figures/%s/wavs/', exp_name)); % get path to raw waveforms
results_file = fullfile(file_path, sprintf('../../results/%s/results.mat', exp_name)); % get path to results
load(results_file)
save_path = fullfile(file_path, sprintf('../../figures/%s/', exp_name)); % get path to save figure
if ~exist(save_path, 'dir'); mkdir(save_path); end

paths = [pathsep, path, pathsep];
utils_path = fullfile(file_path, '../utils/');
on_path  = contains(paths, [pathsep, utils_path, pathsep], 'IgnoreCase', ispc);
if ~ on_path; addpath(utils_path); end

% settings for waveforms...note must match saved waveforms in wav_path
repeated_bgs = [2, 29, 26];
stims_fg = [1 2 3 1 1, ...
            25 27 1 1 7, ...
            1 13 19 14 11];
stims_bg = [2 2 2 2 2, ...
            29 29 29 29 29, ...
            26 26 26 26 26];
outliers = [1 2 2 1 1, ...
            2 2 1 1 2, ...
            1 2 2 2 2];
times = [NaN 7 1 NaN NaN, ...
         8 6 NaN NaN 3, ...
         NaN 5 2  9 4]; % NaNs where no outlier
times(isnan(times)) = 1;
trials = [1:3 39 40, ...
          [1:3 39 40] + 40, ...
          [1:3 39 40] + 280]; % trial numbers
n_stims = 15;
n_stims_per_col = 5;
curr_reps = zeros(max(stims_bg), 1);
blocks = [1 2 8];

% stimulus settings (taken from stimulus creation script)
seglen = 500; % of the foreground (in MS)
bg_length = 3.25; % length of background in seconds
audio_sr = 20000; % sample rate of all the audio
SNR = -8; % set of possible SNRs
onset_times = linspace(0.25,2.50,10); % set of possible positions of the outlier

% waveform plot settings
n_positions = 6;
skip_inds = 4;
clip_thresh = 0.25; % value at which to clip waveforms to make them look a bit nicer
y_positions = flip((1:(n_positions)))*clip_thresh*2.4; % vector of y-positions of waveforms
skip_position = y_positions(skip_inds); % positions that are skipped and replaced with ellipsis
y_positions(skip_inds) = []; % actually remove the skipped positions
y_positions(skip_inds(1):end) = y_positions(skip_inds(1):end) + clip_thresh*2.4/2;
skip_position = skip_position + clip_thresh*2.4/4;
rep_colors = [152 9 15; 217 47 49; 192 9 19]/255;

%% plot diagram

% initialize diagram figure
figure('units','inches','innerposition',[0 0 11 5])
hold on

% loop over trials
for n = 1:n_stims
    % determine which column to plot in
    col = ceil(n/n_stims_per_col);
    if col==3
        x_spacing = 1.6*audio_sr; 
    else
        x_spacing = 1.5*audio_sr; 
    end
    x_start = ((bg_length*audio_sr) + x_spacing)*(col-1);
    y_pos_ind = n - n_stims_per_col*(col-1);

    % load background waveform
    bg = stims_bg(n);
    rep = curr_reps(bg) + 1;
    curr_reps(bg) = rep;
    background_file = sprintf('background_%d_rep_%d.wav', bg, rep);
    background = audioread(fullfile(wav_path, background_file));
    background(abs(background) > clip_thresh) = clip_thresh;

    % load foreground waveform
    fg = stims_fg(n);
    foreground_file = sprintf('foreground_%d.wav', fg);
    foreground = audioread(fullfile(wav_path, foreground_file));

    % get proper onset_time
    t = times(n);
    ind1 = onset_times(t)*audio_sr; % convert the ms timing to a sample

    % set SNR of foreground relative to background
    clip1 = background(ind1:ind1+(audio_sr*seglen/1000)-1); % extract corresponding segment of background            
    setsnr = SNR; % get SNR level
    c = (rms(clip1)/rms(foreground))*10^(setsnr/20); % compute the scaling factor to get the proper SNR
    outlier_switch = outliers(n);
    if outlier_switch == 1 % then no outlier
        c = 0; % will kill the outlier
    else % otherwise keep the right SNR level
    end
    foreground = c*foreground; % actually set proper SNR
    foreground(abs(foreground) > clip_thresh) = clip_thresh; % apply clip threshold
    
    % select background color
    if any(bg == repeated_bgs)
        which_bg = find(bg == repeated_bgs); 
        bg_color = rep_colors(which_bg, :);
    else
        bg_color = 0.6*ones(3,1);
    end

    % plot stimulus
    bg_timevec = (1:length(background)) + x_start;
    plot(bg_timevec, background + y_positions(y_pos_ind), 'LineWidth', 1.5, 'Color', bg_color)
    if outlier_switch==2 % then also plot the foreground
        temp_timevec = (1:length(background)) + x_start;
        fg_timevec = temp_timevec(ind1:ind1+(audio_sr*seglen/1000)-1);
        plot(fg_timevec, foreground + y_positions(y_pos_ind), 'LineWidth', 1.5, 'Color', 'k')
    end
   
    % add trial number label
    text(-0.4*audio_sr + x_start, y_positions(y_pos_ind), sprintf('%d', trials(n)), 'HorizontalAlignment','center', 'FontSize', 18, 'VerticalAlignment','middle')
    
    if mod(n, n_stims_per_col)==0 % if done plotting everything in column then add extra details
        % add ellipsis
        text(-0.4*audio_sr + 0.025*audio_sr + x_start, skip_position(1), '. . .', 'HorizontalAlignment','center', 'FontSize', 22, 'Rotation', 90, 'VerticalAlignment','baseline')
        text(length(background)/2 + x_start, skip_position(1), '. . .', 'HorizontalAlignment','center', 'FontSize', 22, 'Rotation', 90, 'VerticalAlignment','baseline')                
        
        % add block text
        block_text = text(x_start - 0.9*audio_sr, skip_position(1), sprintf('block %d', blocks(col)), ...
                    'HorizontalAlignment','center', 'FontSize', 20, 'Rotation', 90, 'VerticalAlignment','middle', ...
                    'Color',0.6*ones(3,1));
        
        % plot lines denoting block
        plot([x_start - 0.9*audio_sr + 0.01*audio_sr, x_start - 0.9*audio_sr + 0.01*audio_sr], ...
             [y_positions(end) - clip_thresh*1.05, block_text.Extent(2) - clip_thresh*0.75], ...
             'LineWidth', 1.5, 'Color', 0.6*ones(3,1))
        plot([x_start - 0.9*audio_sr + 0.01*audio_sr, x_start - 0.9*audio_sr + 0.1*audio_sr], ...
             [y_positions(end) - clip_thresh*1.05, y_positions(end) - clip_thresh*1.05], ...
             'LineWidth', 1.5, 'Color', 0.6*ones(3,1))
        plot([x_start - 0.9*audio_sr + 0.01*audio_sr, x_start - 0.9*audio_sr + 0.01*audio_sr], ...
             [block_text.Extent(2) + block_text.Extent(4) + clip_thresh*0.75, y_positions(1) + clip_thresh*1.05], ...
             'LineWidth', 1.5, 'Color', 0.6*ones(3,1))
        plot([x_start - 0.9*audio_sr + 0.01*audio_sr, x_start - 0.9*audio_sr + 0.1*audio_sr], ...
             [y_positions(1) + clip_thresh*1.05, y_positions(1) + clip_thresh*1.05], ...
             'LineWidth', 1.5, 'Color', 0.6*ones(3,1))

        if col==3 % if last column
            % add block ellipsis
            text(x_start - 1.45*audio_sr, skip_position(1), '. . .', 'HorizontalAlignment','center', 'FontSize', 22, 'VerticalAlignment','baseline')
        end
    end

end

% set axes limits
x_end = x_start + (bg_length*audio_sr);
xlim([-1*audio_sr, x_end + 0.3*audio_sr])
ylim([0.5, 5.1])
axis off

% add more text to figure
text(-0.4*audio_sr, 4.15, 'trial', 'HorizontalAlignment','center', 'FontSize', 20, 'VerticalAlignment','middle')
text(length(background)/2, 4.15, 'stimulus', 'HorizontalAlignment','center', 'FontSize', 20, 'VerticalAlignment','middle')

% save diagram figure
set(gca,'LooseInset',get(gca,'TightInset')); % remove space around edges
ax = gca;
ax.Position = [0.01 0.01 1 1.1];
save_name = fullfile(save_path, sprintf('%s_diagram.pdf', exp_name));
set(gcf, 'PaperPositionMode', 'auto');
fig = gcf;
fig.Renderer = 'Painters'; % to override forced bitmap due to so many graphics
pos = fig.Position;
set(fig, 'PaperSize', [pos(3), pos(4)]);
print(gcf, save_name, '-dpdf');