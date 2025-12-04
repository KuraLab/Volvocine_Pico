function plot_phase_evolution_all_springs(base_dir, n_seconds_to_cut, plot_duration, apply_filter, filter_window_size, do_save_figure, n_sync, m_sync, sample_window, close_after_plot)
% Run plot_phase_evolution_omega245 for every data folder inside Spring* directories under base_dir.
%
% Usage:
%   plot_phase_evolution_all_springs()
%   plot_phase_evolution_all_springs('EstimateF', ...)
%
% Parameters mirror plot_phase_evolution_omega245; defaults match that function.
% Additional option close_after_plot controls whether each figure is closed after processing.

    if nargin < 1 || isempty(base_dir)
        base_dir = 'EstimateF';
    end
    if nargin < 2
        n_seconds_to_cut = [];
    end
    if nargin < 3
        plot_duration = [];
    end
    if nargin < 4
        apply_filter = [];
    end
    if nargin < 5
        filter_window_size = [];
    end
    if nargin < 6
        do_save_figure = [];
    end
    if nargin < 7
        n_sync = [];
    end
    if nargin < 8
        m_sync = [];
    end
    if nargin < 9
        sample_window = [];
    end
    if nargin < 10 || isempty(close_after_plot)
        close_after_plot = false;
    end

    if ~isfolder(base_dir)
        error('Base directory not found: %s', base_dir);
    end

    entries = dir(fullfile(base_dir, 'Spring*'));
    entries = entries([entries.isdir]);
    entries = entries(~ismember({entries.name},{'.','..'}));

    if isempty(entries)
        warning('No Spring* folders found under %s.', base_dir);
        return;
    end

    freq_vals = [];
    mean_vals = [];
    stderr_vals = [];
    spring_idx = [];
    spring_labels = {entries.name};

    for i = 1:numel(entries)
        spring_path = fullfile(base_dir, entries(i).name);
        spring_subdirs = dir(spring_path);
        spring_subdirs = spring_subdirs([spring_subdirs.isdir]);
        spring_subdirs = spring_subdirs(~ismember({spring_subdirs.name},{'.','..'}));

        fprintf('[INFO] Processing %s (%d subfolders).\n', spring_path, numel(spring_subdirs));

        for j = 1:numel(spring_subdirs)
            data_path = fullfile(spring_path, spring_subdirs(j).name);
            csvs = dir(fullfile(data_path, '*.csv'));
            if isempty(csvs)
                continue;
            end

            fprintf('  [INFO] Plotting %s ...\n', data_path);
            cluster_info = [];
            try
                cluster_info = plot_phase_evolution_omega245(data_path, n_seconds_to_cut, plot_duration, ...
                    apply_filter, filter_window_size, do_save_figure, n_sync, m_sync, sample_window);
                drawnow;
                if close_after_plot
                    fig_handle = gcf;
                    if ishghandle(fig_handle)
                        close(fig_handle);
                    end
                end
            catch ME
                warning('  [WARN] Failed for %s: %s', data_path, ME.message);
            end

            freq_val = str2double(spring_subdirs(j).name);
            if isnan(freq_val)
                continue;
            end

            if isempty(cluster_info)
                continue;
            end

            k = max(1, cluster_info(1).n_clusters);
            if k == 1
                offsets = 0;
            else
                offsets = linspace(-0.12, 0.12, k);
            end
            for c = 1:numel(cluster_info)
                cid = cluster_info(c).cluster_id;
                offset = offsets(min(max(cid,1), numel(offsets)));
                freq_vals(end+1) = freq_val + offset; %#ok<AGROW>
                mean_vals(end+1) = cluster_info(c).mean; %#ok<AGROW>
                stderr_vals(end+1) = cluster_info(c).stderr; %#ok<AGROW>
                spring_idx(end+1) = i; %#ok<AGROW>
            end
        end
    end

    if isempty(freq_vals)
        warning('No numeric Spring frequencies produced phase statistics; skipping combined plot.');
        return;
    end

    figure;
    hold on;
    colors = lines(max(1, numel(entries)));
    marker_size = 40;
    legend_handles = gobjects(0);
    legend_labels = {};

    for i = 1:numel(entries)
        mask = spring_idx == i;
        if ~any(mask)
            continue;
        end
        xvals = mean_vals(mask);
        yvals = freq_vals(mask);
        h = scatter(xvals, yvals, marker_size, 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', colors(i,:), 'DisplayName', spring_labels{i});
        legend_handles(end+1) = h; %#ok<AGROW>
        legend_labels{end+1} = spring_labels{i}; %#ok<AGROW>
    end

    if ~isempty(legend_handles)
        legend(legend_handles, legend_labels, 'Location', 'best');
    end
    ylabel('Frequency (folder name)');
    xlabel('Mean phase difference (rad)');
    grid on;
    xlim([-pi, pi]);
    xticks([-pi, -pi/2, 0, pi/2, pi]);
    xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
    hold off;

    tuneFigure;
end
