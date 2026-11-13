function tune_analysis_figure(fig)
% Apply the user's tuneFigure style while preserving mathematical labels.
% tuneFigure selects TeX for labels without double-dollar delimiters. Keep
% each original interpreter so existing single-dollar LaTeX remains valid.
% Preserve the size of multi-panel figures, which otherwise shrink to the
% single-panel default and become crowded at the tuned font sizes.

    assert(exist('tuneFigure', 'file') == 2, ...
        'tuneFigure.m must be available on the MATLAB path.');
    text_objects = findall(fig, '-property', 'Interpreter');
    interpreters = cell(numel(text_objects), 1);
    for k = 1:numel(text_objects)
        interpreters{k} = text_objects(k).Interpreter;
    end
    position = fig.Position;
    has_multiple_axes = numel(findall(fig, 'Type', 'axes')) > 1;

    tuneFigure(fig);

    if has_multiple_axes, fig.Position = position; end
    for k = 1:numel(text_objects)
        if isgraphics(text_objects(k))
            text_objects(k).Interpreter = interpreters{k};
        end
    end
    drawnow;
end
