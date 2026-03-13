function out = plot_exported_gamma(source_path, psi_grid)
% plot_exported_gamma Load and display exported Gamma(psi) results.
%
% Usage:
%   plot_exported_gamma(fullfile('EstimateQ', 'Spring1', '255', 'gamma_exports'))
%   plot_exported_gamma(fullfile('EstimateQ', 'Spring1', '255'))
%   plot_exported_gamma(fullfile('EstimateQ', 'Spring1', '255', 'gamma_exports', 'gamma_export_latest.mat'))
%   plot_exported_gamma(fullfile('EstimateQ', 'Spring1', '255', 'gamma_exports', 'gamma_curve_agent2_derived.csv'))
%
% Inputs:
%   source_path : gamma_exports folder, its parent folder, exported MAT file,
%                 or exported CSV file. If omitted, the function first looks
%                 for a gamma_exports folder under the current directory and
%                 then searches recursively below it.
%   psi_grid    : optional evaluation grid used when plotting MAT exports.
%                 CSV exports are plotted on their stored grid.
%
% Output:
%   out: struct describing the loaded source, plotted entries, and figure.

    existing_figures = findall(0, 'Type', 'figure');
    if ~isempty(existing_figures)
        close(existing_figures);
        drawnow;
    end

    % ===== User setting =====
    % Single source: set a folder/MAT/CSV path, or [] for auto-discovery.
    configured_source_path = [];

    % Multi-spring overlay: list gamma_exports folders for each spring.
    % When this cell array has 2+ entries, true_gamma curves are overlaid
    % on one figure.  Set to {} to disable.
    configured_overlay_paths = {
        fullfile('EstimateQ', 'Spring1', '255', 'gamma_exports'), ...
        fullfile('EstimateQ', 'Spring2', '255', 'gamma_exports'), ...
        fullfile('EstimateQ', 'Spring3', '255', 'gamma_exports'), ...
        fullfile('EstimateQ', 'Spring4', '255', 'gamma_exports'), ...
        fullfile('EstimateQ', 'Spring5', '255', 'gamma_exports'), ...
    };
    % ========================

    if nargin < 2 || isempty(psi_grid)
        psi_grid = [];
    else
        validateattributes(psi_grid, {'numeric'}, {'vector', 'finite'}, mfilename, 'psi_grid');
        psi_grid = psi_grid(:);
    end

    if nargin < 1 || isempty(source_path)
        % Multi-spring overlay takes priority when paths are configured.
        if ~isempty(configured_overlay_paths)
            out = plot_true_gamma_overlay(configured_overlay_paths, psi_grid);
            return;
        end
        source_path = resolve_initial_source_path(configured_source_path);
    end

    source_info = resolve_gamma_source_path(source_path);
    switch source_info.source_type
        case 'bundle'
            out = plot_gamma_bundle(source_info, psi_grid);
        case 'csv'
            out = plot_single_gamma_csv(source_info.file_path);
        case 'csv-folder'
            out = plot_gamma_csv_folder(source_info);
        otherwise
            error('Unsupported source type: %s', source_info.source_type);
    end
end

function out = plot_true_gamma_overlay(source_paths, psi_grid)
% true_gamma: all springs overlaid on one figure.
% derived / a2 gamma: one figure per spring directory.
    spring_colors = get(groot, 'DefaultAxesColorOrder');
    agent_styles  = {'-', '--', ':', '-.'};  % line style per agent index

    empty_lines = struct('x', {}, 'y', {}, 'name', {}, 'line_style', {}, 'color', {}, 'line_width', {});
    true_lines = empty_lines;
    figures    = gobjects(1, 0);

    for idx = 1:numel(source_paths)
        folder   = source_paths{idx};
        mat_path = fullfile(folder, 'gamma_export_latest.mat');
        if exist(mat_path, 'file') ~= 2
            warning('plot_exported_gamma:missingFile', 'gamma_export_latest.mat not found: %s', folder);
            continue;
        end
        loaded = load(mat_path);
        if ~isfield(loaded, 'gamma_export')
            warning('plot_exported_gamma:noExport', 'No gamma_export in: %s', mat_path);
            continue;
        end
        gamma_export = loaded.gamma_export;
        color        = spring_colors(mod(idx - 1, size(spring_colors, 1)) + 1, :);
        spring_label = extract_spring_label(folder);

        % --- true_gamma: collect for shared overlay + individual figure ---
        if isfield(gamma_export, 'true_gamma') && isstruct(gamma_export.true_gamma) && ...
                isfield(gamma_export.true_gamma, 'available') && gamma_export.true_gamma.available
            try
                [psi_v, gam_v] = reconstruct_exported_gamma(gamma_export.true_gamma, psi_grid, 'true');
                if ~isempty(psi_v)
                    true_lines(end + 1) = make_line(psi_v, gam_v, spring_label, '-', color, 1.8); %#ok<AGROW>
                    individual_line = make_line(psi_v, gam_v, 'true gamma', '-', color, 1.8);
                    figures(end + 1) = make_overlay_figure(individual_line, ...
                        sprintf('%s: true gamma', spring_label)); %#ok<AGROW>
                end
            catch err
                warning('plot_exported_gamma:reconstructFailed', '%s true_gamma: %s', spring_label, err.message);
            end
        end

        % --- derived / a2: one figure per spring ---
        if isfield(gamma_export, 'agents') && ~isempty(gamma_export.agents)
            derived_lines = empty_lines;
            a2_lines      = empty_lines;
            for agent_idx = 1:numel(gamma_export.agents)
                agent_export = gamma_export.agents(agent_idx);
                agent_id     = agent_export.agent_id;
                line_style   = agent_styles{mod(agent_idx - 1, numel(agent_styles)) + 1};
                line_name    = sprintf('agent %d', agent_id);

                derived_lines = append_overlay_line(derived_lines, agent_export.derived_gamma, ...
                    psi_grid, line_name, line_style, color);
                a2_lines = append_overlay_line(a2_lines, agent_export.a2_gamma, ...
                    psi_grid, line_name, line_style, color);
            end
            if ~isempty(derived_lines)
                figures(end + 1) = make_overlay_figure(derived_lines, ...
                    sprintf('%s: derived gamma', spring_label)); %#ok<AGROW>
            end
            if ~isempty(a2_lines)
                figures(end + 1) = make_overlay_figure(a2_lines, ...
                    sprintf('%s: a2 gamma', spring_label)); %#ok<AGROW>
            end
        end
    end

    % --- true_gamma overlay figure (all springs) ---
    if ~isempty(true_lines)
        fig = make_overlay_figure(true_lines, ...
            sprintf('true gamma overlay (%d springs)', numel(source_paths)));
        figures = [fig, figures];
    end

    if isempty(figures)
        error('No gamma data could be loaded from any of the specified paths.');
    end

    out = struct();
    out.source_type  = 'overlay';
    out.source_paths = source_paths;
    out.figure       = figures(1);
    out.figures      = figures;
end

function lines = append_overlay_line(lines, gamma_source, psi_grid, line_name, line_style, color)
    gamma_definition = extract_gamma_definition(gamma_source);
    if ~is_resonant_gamma_plottable(gamma_definition)
        return;
    end
    selected_component = get_gamma_definition_component(gamma_definition);
    if strcmp(selected_component, 'antisymmetric') || strcmp(selected_component, 'symmetric')
        plot_component = 'full';
    else
        plot_component = 'selected';
    end
    try
        [psi_v, gam_v] = reconstruct_exported_gamma(gamma_source, psi_grid, plot_component);
    catch
        return;
    end
    if isempty(psi_v) || isempty(gam_v)
        return;
    end
    lines(end + 1) = make_line(psi_v, gam_v, line_name, line_style, color, 1.4); %#ok<AGROW>
end

function line = make_line(psi_v, gam_v, name, line_style, color, line_width)
    line = struct('x', psi_v(:), 'y', gam_v(:), 'name', name, ...
        'line_style', line_style, 'color', color, 'line_width', line_width);
end

function fig = make_overlay_figure(lines, fig_name)
    plot_spec = struct('title', '', 'lines', lines, 'info_text', '');
    fig = figure('Color', 'w', 'Name', fig_name);
    ax  = axes('Parent', fig);
    draw_gamma_plot_spec(ax, plot_spec);
    finalize_gamma_figure(fig);
end

function spring_label = extract_spring_label(folder)
    parts = strsplit(folder, filesep);
    spring_label = 'Spring';
    for p = 1:numel(parts)
        if ~isempty(regexp(parts{p}, '^Spring\d+$', 'once'))
            spring_label = parts{p};
            return;
        end
    end
end

function source_path = resolve_initial_source_path(configured_source_path)
    if nargin >= 1 && ~isempty(configured_source_path)
        if isstring(configured_source_path)
            configured_source_path = char(configured_source_path);
        end
        if ischar(configured_source_path) && (isfolder(configured_source_path) || exist(configured_source_path, 'file') == 2)
            source_path = configured_source_path;
            return;
        end
    end

    source_path = get_default_gamma_source_path();
end

function source_path = get_default_gamma_source_path()
    candidate = fullfile(pwd, 'gamma_exports');
    if isfolder(candidate)
        source_path = candidate;
        return;
    end

    candidate = discover_gamma_export_dir(pwd);
    if ~isempty(candidate)
        source_path = candidate;
    else
        source_path = pwd;
    end
end

function source_info = resolve_gamma_source_path(source_path)
    if isstring(source_path)
        source_path = char(source_path);
    end
    if ~ischar(source_path)
        error('source_path must be a character vector or string scalar.');
    end

    source_info = struct('source_type', '', 'source_path', source_path, 'folder_path', '', 'file_path', '', 'csv_files', {{}});

    if isfolder(source_path)
        export_dir = source_path;
        nested_export_dir = fullfile(source_path, 'gamma_exports');
        if isfolder(nested_export_dir)
            export_dir = nested_export_dir;
        end

        latest_mat = fullfile(export_dir, 'gamma_export_latest.mat');
        if exist(latest_mat, 'file') == 2
            source_info.source_type = 'bundle';
            source_info.folder_path = export_dir;
            source_info.file_path = latest_mat;
            return;
        end

        csv_listing = dir(fullfile(export_dir, '*.csv'));
        if ~isempty(csv_listing)
            [~, sort_idx] = sort({csv_listing.name});
            csv_listing = csv_listing(sort_idx);
            source_info.source_type = 'csv-folder';
            source_info.folder_path = export_dir;
            source_info.csv_files = arrayfun(@(item) fullfile(item.folder, item.name), csv_listing, 'UniformOutput', false);
            return;
        end

        discovered_export_dir = discover_gamma_export_dir(source_path);
        if ~isempty(discovered_export_dir) && ~strcmp(discovered_export_dir, export_dir)
            source_info = resolve_gamma_source_path(discovered_export_dir);
            source_info.source_path = source_path;
            return;
        end

        error('No gamma export files were found under %s.', export_dir);
    end

    if exist(source_path, 'file') ~= 2
        error('Source not found: %s', source_path);
    end

    [parent_dir, ~, ext] = fileparts(source_path);
    switch lower(ext)
        case '.mat'
            source_info.source_type = 'bundle';
            source_info.folder_path = parent_dir;
            source_info.file_path = source_path;
        case '.csv'
            source_info.source_type = 'csv';
            source_info.folder_path = parent_dir;
            source_info.file_path = source_path;
        otherwise
            error('Unsupported source extension: %s', ext);
    end
end

function export_dir = discover_gamma_export_dir(root_dir)
    export_dir = '';
    if ~isfolder(root_dir)
        return;
    end

    latest_bundle_listing = dir(fullfile(root_dir, '**', 'gamma_export_latest.mat'));
    latest_bundle_listing = latest_bundle_listing(~[latest_bundle_listing.isdir]);
    if ~isempty(latest_bundle_listing)
        [~, newest_idx] = max([latest_bundle_listing.datenum]);
        export_dir = latest_bundle_listing(newest_idx).folder;
        return;
    end

    csv_listing = dir(fullfile(root_dir, '**', 'gamma_*.csv'));
    csv_listing = csv_listing(~[csv_listing.isdir]);
    if isempty(csv_listing)
        return;
    end

    valid_mask = false(size(csv_listing));
    for idx = 1:numel(csv_listing)
        valid_mask(idx) = strcmpi(get_last_path_part(csv_listing(idx).folder), 'gamma_exports');
    end
    csv_listing = csv_listing(valid_mask);
    if isempty(csv_listing)
        return;
    end

    [~, newest_idx] = max([csv_listing.datenum]);
    export_dir = csv_listing(newest_idx).folder;
end

function name = get_last_path_part(dir_path)
    [~, name] = fileparts(dir_path);
end

function out = plot_gamma_bundle(source_info, psi_grid)
    loaded = load(source_info.file_path);
    if ~isfield(loaded, 'gamma_export')
        error('The MAT file does not contain a gamma_export variable: %s', source_info.file_path);
    end

    gamma_export = loaded.gamma_export;
    plot_specs = collect_gamma_bundle_plot_specs(gamma_export, psi_grid);
    if isempty(plot_specs)
        error('No plottable gamma entries were found in %s.', source_info.file_path);
    end

    n_plots = numel(plot_specs);
    figures = gobjects(1, n_plots);
    for idx = 1:numel(plot_specs)
        fig = figure('Color', 'w', 'Name', sprintf('Exported gamma [%d/%d]: %s', idx, n_plots, plot_specs(idx).title));
        ax = axes('Parent', fig);
        draw_gamma_plot_spec(ax, plot_specs(idx));
        finalize_gamma_figure(fig);
        figures(idx) = fig;
    end

    out = struct();
    out.source_type = 'bundle';
    out.source_path = source_info.file_path;
    out.figure = figures(1);
    out.figures = figures;
    out.plot_specs = plot_specs;
    out.gamma_export = gamma_export;
    out.bundle_title = build_bundle_title(gamma_export);
end

function plot_specs = collect_gamma_bundle_plot_specs(gamma_export, psi_grid)
    plot_specs = struct('title', {}, 'lines', {}, 'info_text', {});

    if isfield(gamma_export, 'agents') && ~isempty(gamma_export.agents)
        empty_lines = struct('x', {}, 'y', {}, 'name', {}, 'line_style', {}, 'color', {}, 'line_width', {});
        derived_lines = empty_lines;
        a2_lines      = empty_lines;

        for idx = 1:numel(gamma_export.agents)
            agent_export = gamma_export.agents(idx);
            agent_color  = get_agent_color(idx);
            derived_lines = append_agent_gamma_line(derived_lines, agent_export.derived_gamma, psi_grid, agent_export.agent_id, agent_color);
            a2_lines      = append_agent_gamma_line(a2_lines,      agent_export.a2_gamma,     psi_grid, agent_export.agent_id, agent_color);
        end

        if ~isempty(derived_lines)
            plot_specs(end + 1) = struct('title', 'derived gamma', 'lines', derived_lines, 'info_text', ''); %#ok<AGROW>
        end
        if ~isempty(a2_lines)
            plot_specs(end + 1) = struct('title', 'a2 gamma', 'lines', a2_lines, 'info_text', ''); %#ok<AGROW>
        end
    end

    if isfield(gamma_export, 'true_gamma') && isstruct(gamma_export.true_gamma) && ...
            isfield(gamma_export.true_gamma, 'available') && gamma_export.true_gamma.available
        true_spec = build_true_gamma_plot_spec(gamma_export.true_gamma, psi_grid);
        if ~isempty(true_spec)
            plot_specs(end + 1) = true_spec; %#ok<AGROW>
        end
    end
end

function lines = append_agent_gamma_line(lines, gamma_source, psi_grid, agent_id, agent_color)
    gamma_definition = extract_gamma_definition(gamma_source);
    if ~is_resonant_gamma_plottable(gamma_definition)
        return;
    end
    selected_component = get_gamma_definition_component(gamma_definition);
    if strcmp(selected_component, 'antisymmetric') || strcmp(selected_component, 'symmetric')
        plot_component = 'full';
    else
        plot_component = 'selected';
    end
    try
        [psi_values, gamma_values] = reconstruct_exported_gamma(gamma_source, psi_grid, plot_component);
    catch
        return;
    end
    if isempty(psi_values) || isempty(gamma_values)
        return;
    end
    lines(end + 1) = struct( ...
        'x', psi_values(:), ...
        'y', gamma_values(:), ...
        'name', sprintf('agent %d', agent_id), ...
        'line_style', '-', ...
        'color', agent_color, ...
        'line_width', 1.8); %#ok<AGROW>
end

function agent_color = get_agent_color(agent_order_idx)
    palette = get(groot, 'DefaultAxesColorOrder');
    agent_color = palette(mod(agent_order_idx - 1, size(palette, 1)) + 1, :);
end

function spec = build_true_gamma_plot_spec(true_gamma, psi_grid)
    spec = struct([]);
    if ~isstruct(true_gamma) || ~isfield(true_gamma, 'available') || ~true_gamma.available
        return;
    end

    line_specs = {
        'true', 'true gamma', '-', [0.85, 0.33, 0.10], 2.0};

    lines = struct('x', {}, 'y', {}, 'name', {}, 'line_style', {}, 'color', {}, 'line_width', {});
    for idx = 1:size(line_specs, 1)
        try
            [psi_values, gamma_values] = reconstruct_exported_gamma(true_gamma, psi_grid, line_specs{idx, 1});
        catch
            continue;
        end

        if isempty(psi_values) || isempty(gamma_values)
            continue;
        end

        lines(end + 1) = struct( ...
            'x', psi_values(:), ...
            'y', gamma_values(:), ...
            'name', line_specs{idx, 2}, ...
            'line_style', line_specs{idx, 3}, ...
            'color', line_specs{idx, 4}, ...
            'line_width', line_specs{idx, 5}); %#ok<AGROW>
    end

    if isempty(lines)
        return;
    end

    spec = struct();
    spec.title = sprintf('True gamma: agent %d to %d', true_gamma.agent_id_1, true_gamma.agent_id_2);
    spec.lines = lines;
    spec.info_text = sprintf('gamma_true = agent %d minus agent %d', true_gamma.agent_id_2, true_gamma.agent_id_1);
end

function out = plot_single_gamma_csv(csv_path)
    table_data = readtable(csv_path);
    plot_spec = build_csv_plot_spec(table_data, get_file_title(csv_path));

    fig = figure('Color', 'w', 'Name', sprintf('Exported gamma CSV: %s', csv_path));
    ax = axes('Parent', fig);
    draw_gamma_plot_spec(ax, plot_spec);
    finalize_gamma_figure(fig);

    out = struct();
    out.source_type = 'csv';
    out.source_path = csv_path;
    out.figure = fig;
    out.table = table_data;
    out.plot_specs = plot_spec;
end

function out = plot_gamma_csv_folder(source_info)
    n_files = numel(source_info.csv_files);
    if n_files < 1
        error('No CSV files were found under %s.', source_info.folder_path);
    end

    plot_specs = repmat(struct('title', '', 'lines', struct([]), 'info_text', ''), 1, n_files);
    tables = cell(1, n_files);
    figures = gobjects(1, n_files);

    for idx = 1:n_files
        tables{idx} = readtable(source_info.csv_files{idx});
        plot_specs(idx) = build_csv_plot_spec(tables{idx}, get_file_title(source_info.csv_files{idx}));
        fig = figure('Color', 'w', 'Name', sprintf('Exported gamma CSV [%d/%d]: %s', idx, n_files, plot_specs(idx).title));
        ax = axes('Parent', fig);
        draw_gamma_plot_spec(ax, plot_specs(idx));
        finalize_gamma_figure(fig);
        figures(idx) = fig;
    end

    out = struct();
    out.source_type = 'csv-folder';
    out.source_path = source_info.folder_path;
    out.figure = figures(1);
    out.figures = figures;
    out.tables = tables;
    out.plot_specs = plot_specs;
end

function plot_spec = build_csv_plot_spec(table_data, plot_title)
    if ~istable(table_data) || width(table_data) < 2
        error('CSV data must contain at least two columns.');
    end
    if ~strcmp(table_data.Properties.VariableNames{1}, 'psi')
        error('The first CSV column must be psi.');
    end

    variable_names = table_data.Properties.VariableNames;
    has_true_gamma = any(strcmpi(variable_names, 'gamma_true'));
    psi_values = table_data.psi(:);
    lines = struct('x', {}, 'y', {}, 'name', {}, 'line_style', {}, 'color', {}, 'line_width', {});
    for idx = 2:width(table_data)
        variable_name = variable_names{idx};
        if strcmpi(variable_name, 'gamma_symmetric') || strcmpi(variable_name, 'gamma_antisymmetric')
            continue;
        end
        if has_true_gamma && (strcmpi(variable_name, 'gamma_agent_1') || strcmpi(variable_name, 'gamma_agent_2'))
            continue;
        end
        style = get_csv_variable_style(variable_name);
        lines(end + 1) = struct( ...
            'x', psi_values, ...
            'y', table_data.(variable_name)(:), ...
            'name', strrep(variable_name, '_', ' '), ...
            'line_style', style.line_style, ...
            'color', style.color, ...
            'line_width', style.line_width); %#ok<AGROW>
    end

    plot_spec = struct();
    plot_spec.title = plot_title;
    plot_spec.lines = lines;
    plot_spec.info_text = sprintf('%d samples', numel(psi_values));
end

function draw_gamma_plot_spec(ax, plot_spec)
    hold(ax, 'on');
    for idx = 1:numel(plot_spec.lines)
        plot(ax, plot_spec.lines(idx).x, plot_spec.lines(idx).y, ...
            'LineStyle', plot_spec.lines(idx).line_style, ...
            'LineWidth', plot_spec.lines(idx).line_width, ...
            'Color', plot_spec.lines(idx).color, ...
            'DisplayName', plot_spec.lines(idx).name);
    end

    plot(ax, [-pi, pi], [0, 0], ':', 'LineWidth', 0.8, 'Color', [0.6, 0.6, 0.6], 'HandleVisibility', 'off');
    xlabel(ax, '$$\psi$$', 'Interpreter', 'latex');
    ylabel(ax, '$$\Gamma(\psi)$$', 'Interpreter', 'latex');

    grid(ax, 'on');
    box(ax, 'on');
    xlim(ax, [-pi, pi]);
    xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax, {'$$-\pi$$', '$$-\pi/2$$', '0', '$$\pi/2$$', '$$\pi$$'});
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    if numel(plot_spec.lines) > 1
        legend(ax, 'Location', 'best');
    end
end

function title_text = build_bundle_title(gamma_export)
    title_text = 'Exported gamma bundle';
    if ~isstruct(gamma_export)
        return;
    end

    if isfield(gamma_export, 'source_dirpath') && isfield(gamma_export, 'gamma_ratio')
        ratio_text = mat2str(gamma_export.gamma_ratio);
        title_text = sprintf('Exported gamma bundle: %s, gamma ratio %s', gamma_export.source_dirpath, ratio_text);
    elseif isfield(gamma_export, 'source_dirpath')
        title_text = sprintf('Exported gamma bundle: %s', gamma_export.source_dirpath);
    end
end

function gamma_definition = extract_gamma_definition(gamma_source)
    gamma_definition = gamma_source;
    if isstruct(gamma_source) && isfield(gamma_source, 'gamma_resonance')
        gamma_definition = gamma_source.gamma_resonance;
    end
end

function tf = is_resonant_gamma_plottable(gamma_definition)
    tf = false;
    if ~isstruct(gamma_definition)
        return;
    end

    if isfield(gamma_definition, 'psi_grid_centered') && ~isempty(gamma_definition.psi_grid_centered)
        tf = true;
        return;
    end
    if isfield(gamma_definition, 'psi_grid') && ~isempty(gamma_definition.psi_grid)
        tf = true;
        return;
    end
    if isfield(gamma_definition, 'harmonic_index') && ~isempty(gamma_definition.harmonic_index)
        tf = true;
    end
end

function component_name = get_gamma_definition_component(gamma_definition)
    component_name = 'full';
    if isstruct(gamma_definition) && isfield(gamma_definition, 'component') && ~isempty(gamma_definition.component)
        component_name = lower(strtrim(char(gamma_definition.component)));
    end
end

function style = get_csv_variable_style(variable_name)
    switch lower(variable_name)
        case 'gamma_selected'
            style = struct('line_style', '-', 'color', [0.00, 0.45, 0.74], 'line_width', 2.0);
        case 'gamma_full'
            style = struct('line_style', '--', 'color', [0.20, 0.20, 0.20], 'line_width', 1.4);
        case 'gamma_antisymmetric'
            style = struct('line_style', '-.', 'color', [0.00, 0.45, 0.74], 'line_width', 1.2);
        case 'gamma_true'
            style = struct('line_style', '-', 'color', [0.85, 0.33, 0.10], 'line_width', 2.0);
        case 'gamma_agent_1'
            style = struct('line_style', ':', 'color', [0.45, 0.45, 0.45], 'line_width', 1.2);
        case 'gamma_agent_2'
            style = struct('line_style', '--', 'color', [0.00, 0.45, 0.74], 'line_width', 1.2);
        otherwise
            style = struct('line_style', '-', 'color', [0.30, 0.30, 0.30], 'line_width', 1.2);
    end
end

function title_text = get_file_title(file_path)
    [~, title_text, ext] = fileparts(file_path);
    title_text = [title_text ext];
end

function finalize_gamma_figure(fig)
    if exist('tuneFigure', 'file') == 2
        tuneFigure;
    else
        set(fig, 'Renderer', 'painters');
    end
end