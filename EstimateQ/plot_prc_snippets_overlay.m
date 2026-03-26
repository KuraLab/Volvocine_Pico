function out = plot_prc_snippets_overlay(input_path, n_samples)
% Reconstruct and overlay PRCs from snippet txt files.
%
% Usage:
%   plot_prc_snippets_overlay()
%   plot_prc_snippets_overlay(fullfile('Spring1','255','gamma_exports'))
%   plot_prc_snippets_overlay(fullfile('Spring1','255','gamma_exports','prc_snippet_ref_w1.txt'))
%
% If input_path is a folder, all prc_snippet_*.txt files are plotted.
% If input_path is a file, only that snippet is plotted.

    if nargin < 1 || isempty(input_path)
        input_path = resolve_default_snippet_dir();
    else
        input_path = resolve_input_path(input_path);
    end
    if nargin < 2 || isempty(n_samples)
        n_samples = 1024;
    end

    txt_files = collect_snippet_files(input_path);
    if isempty(txt_files)
        error('No snippet txt files found for: %s\nRun from EstimateQ or pass an explicit path (e.g., fullfile(''Spring1'',''255'',''gamma_exports'')).', input_path);
    end

    psi = linspace(-pi, pi, n_samples).';
    n_files = numel(txt_files);

    labels = cell(n_files, 1);
    prc_harmonics = zeros(n_files, 1);
    z_matrix = zeros(n_samples, n_files);

    for k = 1:n_files
        file_path = txt_files{k};
        [prc_a, prc_b, label_text] = parse_prc_snippet_file(file_path);
        labels{k} = label_text;

        n_h = min(numel(prc_a), numel(prc_b)) - 1;
        prc_harmonics(k) = n_h;

        z = zeros(size(psi));
        for n = 1:n_h
            z = z + prc_a(n + 1) * cos(n * psi) + prc_b(n + 1) * sin(n * psi);
        end
        z_matrix(:, k) = z;
    end

    fig = figure('Color', 'w', 'Name', 'Reconstructed PRC snippets');
    ax = axes('Parent', fig);
    hold(ax, 'on');

    cmap = lines(max(n_files, 3));
    for k = 1:n_files
        plot(ax, psi, z_matrix(:, k), 'LineWidth', 1.8, 'Color', cmap(k, :), 'DisplayName', labels{k});
    end
    plot(ax, psi, zeros(size(psi)), ':', 'LineWidth', 1.0, 'Color', [0.6, 0.6, 0.6], 'DisplayName', '0 line');

    xlabel(ax, '$$\psi$$', 'Interpreter', 'latex');
    ylabel(ax, '$$z(\psi)$$', 'Interpreter', 'latex');
    title(ax, 'Reconstructed PRCs from snippet txt', 'Interpreter', 'none');

    xlim(ax, [-pi, pi]);
    xticks(ax, [-pi, -pi/2, 0, pi/2, pi]);
    xticklabels(ax, {'$$-\pi$$', '$$-\pi/2$$', '0', '$$\pi/2$$', '$$\pi$$'});
    grid(ax, 'on');
    box(ax, 'on');

    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';

    legend(ax, 'Location', 'best', 'Interpreter', 'none');

    out = struct();
    out.files = txt_files;
    out.labels = labels;
    out.prc_harmonics = prc_harmonics;
    out.psi = psi;
    out.z_matrix = z_matrix;
    out.figure = fig;

    fprintf('[INFO] Reconstructed and plotted %d snippet file(s).\n', n_files);
end

function txt_files = collect_snippet_files(input_path)
    txt_files = {};

    if isfolder(input_path)
        files = dir(fullfile(input_path, 'prc_snippet_*.txt'));
        if isempty(files)
            return;
        end

        [~, order] = sort({files.name});
        files = files(order);

        txt_files = cell(numel(files), 1);
        for i = 1:numel(files)
            txt_files{i} = fullfile(files(i).folder, files(i).name);
        end
        return;
    end

    if isfile(input_path)
        txt_files = {input_path};
    end
end

function resolved_path = resolve_default_snippet_dir()
    candidate_pwd = fullfile(pwd, 'Spring1', '255', 'gamma_exports');
    if isfolder(candidate_pwd)
        resolved_path = candidate_pwd;
        return;
    end

    base_dir = fileparts(mfilename('fullpath'));
    resolved_path = fullfile(base_dir, 'Spring1', '255', 'gamma_exports');
end

function resolved_path = resolve_input_path(input_path)
    if isfolder(input_path) || isfile(input_path)
        resolved_path = input_path;
        return;
    end

    candidate_pwd = fullfile(pwd, input_path);
    if isfolder(candidate_pwd) || isfile(candidate_pwd)
        resolved_path = candidate_pwd;
        return;
    end

    base_dir = fileparts(mfilename('fullpath'));
    resolved_path = fullfile(base_dir, input_path);
end

function [prc_a, prc_b, label_text] = parse_prc_snippet_file(file_path)
    raw_text = fileread(file_path);
    lines_txt = splitlines(string(raw_text));

    prc_harmonics = [];
    label_text = '';

    for i = 1:numel(lines_txt)
        line = char(lines_txt(i));

        tok_label = regexp(line, '^\s*#\s*Auto-generated from\s*(.*)\s*$', 'tokens', 'once');
        if ~isempty(tok_label)
            label_text = strtrim(tok_label{1});
        end

        tok_h = regexp(line, '^\s*prc_harmonics\s*=\s*(\d+)\s*$', 'tokens', 'once');
        if ~isempty(tok_h)
            prc_harmonics = str2double(tok_h{1});
        end
    end

    if isempty(label_text)
        [~, name, ext] = fileparts(file_path);
        label_text = [name ext];
    end

    if isempty(prc_harmonics) || ~isfinite(prc_harmonics) || prc_harmonics < 1
        error('Could not parse prc_harmonics from %s', file_path);
    end

    prc_harmonics = floor(prc_harmonics);
    prc_a = zeros(prc_harmonics + 1, 1);
    prc_b = zeros(prc_harmonics + 1, 1);

    num_pattern = '([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)';
    pat_a = ['^\s*prc_a\[(\d+)\]\s*=\s*' num_pattern '\s*$'];
    pat_b = ['^\s*prc_b\[(\d+)\]\s*=\s*' num_pattern '\s*$'];

    for i = 1:numel(lines_txt)
        line = char(lines_txt(i));

        tok_a = regexp(line, pat_a, 'tokens', 'once');
        if ~isempty(tok_a)
            idx = str2double(tok_a{1});
            val = str2double(tok_a{2});
            if isfinite(idx) && isfinite(val) && idx >= 1 && idx <= prc_harmonics
                prc_a(idx + 1) = val;
            end
            continue;
        end

        tok_b = regexp(line, pat_b, 'tokens', 'once');
        if ~isempty(tok_b)
            idx = str2double(tok_b{1});
            val = str2double(tok_b{2});
            if isfinite(idx) && isfinite(val) && idx >= 1 && idx <= prc_harmonics
                prc_b(idx + 1) = val;
            end
        end
    end
end
