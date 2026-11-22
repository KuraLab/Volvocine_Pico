function [prc_a, prc_b, label_text] = parse_prc_snippet_file(file_path)
% Parse the Fourier coefficients and label in a firmware PRC snippet.
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

    if isempty(prc_harmonics) || ~isfinite(prc_harmonics) || prc_harmonics < 0
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
            if isfinite(idx) && isfinite(val) && idx >= 0 && idx <= prc_harmonics
                prc_a(idx + 1) = val;
            end
            continue;
        end

        tok_b = regexp(line, pat_b, 'tokens', 'once');
        if ~isempty(tok_b)
            idx = str2double(tok_b{1});
            val = str2double(tok_b{2});
            if isfinite(idx) && isfinite(val) && idx >= 0 && idx <= prc_harmonics
                prc_b(idx + 1) = val;
            end
        end
    end
end
