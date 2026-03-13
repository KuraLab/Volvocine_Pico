function [psi_grid, gamma_values, meta] = reconstruct_exported_gamma(gamma_source, psi_grid, component)
% reconstruct_exported_gamma Re-evaluate exported Gamma(psi) data.
%
% Usage:
%   S = load(fullfile(dirpath, 'gamma_exports', 'gamma_export_latest.mat'));
%   [psi, gamma_true] = reconstruct_exported_gamma(S.gamma_export.true_gamma);
%   [psi, gamma_full] = reconstruct_exported_gamma(S.gamma_export.agents(1).derived_gamma, [], 'full');
%   [psi, gamma_selected] = reconstruct_exported_gamma(S.gamma_export.agents(1).derived_gamma.gamma_resonance);
%
% Inputs:
%   gamma_source : exported gamma struct. This can be either
%                  - gamma_export.true_gamma
%                  - gamma_export.agents(k).derived_gamma
%                  - gamma_export.agents(k).derived_gamma.gamma_resonance
%                  - gamma_export.agents(k).a2_gamma
%   psi_grid     : optional column vector of psi values. If omitted, the
%                  stored sampling grid is returned, or a default
%                  linspace(-pi, pi, 512) grid is used for harmonic data.
%   component    : optional component selector.
%                  For resonant gamma: 'selected', 'full', 'symmetric',
%                  'antisymmetric'.
%                  For true gamma: 'true', 'agent1', 'agent2'.
%
% Outputs:
%   psi_grid     : evaluation grid.
%   gamma_values : Gamma evaluated on psi_grid.
%   meta         : metadata describing the resolved source and component.

    if nargin < 2
        psi_grid = [];
    end
    if nargin < 3
        component = '';
    end

    if nargin == 2 && (ischar(psi_grid) || (isstring(psi_grid) && isscalar(psi_grid)))
        component = psi_grid;
        psi_grid = [];
    end

    if isstring(component)
        component = char(component);
    end
    if isstring(psi_grid)
        psi_grid = char(psi_grid);
    end

    gamma_definition = unwrap_exported_gamma_source(gamma_source);
    if ~isstruct(gamma_definition)
        error('gamma_source must be an exported gamma struct.');
    end

    if isfield(gamma_definition, 'gamma_true') && isfield(gamma_definition, 'psi_grid')
        [psi_grid, gamma_values, meta] = reconstruct_sampled_true_gamma(gamma_definition, psi_grid, component);
        return;
    end

    if isfield(gamma_definition, 'harmonic_index') && isfield(gamma_definition, 'gamma_cos') && isfield(gamma_definition, 'gamma_sin')
        [psi_grid, gamma_values, meta] = reconstruct_harmonic_gamma(gamma_definition, psi_grid, component);
        return;
    end

    if isfield(gamma_definition, 'psi_grid_centered') || isfield(gamma_definition, 'psi_grid')
        [psi_grid, gamma_values, meta] = reconstruct_sampled_gamma(gamma_definition, psi_grid, component);
        return;
    end

    error('Unsupported gamma_source structure.');
end

function gamma_definition = unwrap_exported_gamma_source(gamma_source)
    gamma_definition = gamma_source;
    if isstruct(gamma_source) && isfield(gamma_source, 'gamma_resonance')
        gamma_definition = gamma_source.gamma_resonance;
    end
end

function [psi_grid, gamma_values, meta] = reconstruct_harmonic_gamma(gamma_definition, psi_grid, component)
    harmonic_index = gamma_definition.harmonic_index(:);
    gamma_cos = gamma_definition.gamma_cos(:);
    gamma_sin = gamma_definition.gamma_sin(:);
    selected_component = get_selected_component_name(gamma_definition);
    resolved_component = normalize_resonant_component(component, selected_component);

    if isempty(psi_grid)
        psi_grid = get_default_psi_grid(gamma_definition, 512);
    else
        validateattributes(psi_grid, {'numeric'}, {'vector', 'finite'}, mfilename, 'psi_grid');
        psi_grid = psi_grid(:);
    end

    gamma_symmetric = zeros(size(psi_grid));
    gamma_antisymmetric = zeros(size(psi_grid));
    for idx = 1:numel(harmonic_index)
        k = harmonic_index(idx);
        gamma_symmetric = gamma_symmetric + gamma_cos(idx) * cos(k * psi_grid);
        gamma_antisymmetric = gamma_antisymmetric + gamma_sin(idx) * sin(k * psi_grid);
    end
    gamma_full = gamma_symmetric + gamma_antisymmetric;

    switch resolved_component
        case 'full'
            gamma_values = gamma_full;
        case 'symmetric'
            gamma_values = gamma_symmetric;
        case 'antisymmetric'
            gamma_values = gamma_antisymmetric;
        otherwise
            error('Unsupported resonant component: %s', resolved_component);
    end

    meta = struct();
    meta.source_type = 'harmonic';
    meta.component = resolved_component;
    meta.selected_component = selected_component;
    meta.psi_label = get_optional_field(gamma_definition, 'psi_label', '');
    meta.harmonic_index = harmonic_index;
    meta.gamma_cos = gamma_cos;
    meta.gamma_sin = gamma_sin;
    meta.available = get_optional_field(gamma_definition, 'enabled', true);
end

function [psi_grid, gamma_values, meta] = reconstruct_sampled_gamma(gamma_definition, psi_grid, component)
    selected_component = get_selected_component_name(gamma_definition);
    resolved_component = normalize_resonant_component(component, selected_component);
    sample_psi = get_default_psi_grid(gamma_definition, 512);
    sample_values = get_sampled_resonant_values(gamma_definition, resolved_component);

    if isempty(sample_psi) || isempty(sample_values)
        error('The exported gamma structure does not contain sampled values for component %s.', resolved_component);
    end

    if isempty(psi_grid)
        psi_grid = sample_psi(:);
        gamma_values = sample_values(:);
    else
        validateattributes(psi_grid, {'numeric'}, {'vector', 'finite'}, mfilename, 'psi_grid');
        psi_grid = psi_grid(:);
        gamma_values = interp1(sample_psi(:), sample_values(:), psi_grid, 'linear', 'extrap');
    end

    meta = struct();
    meta.source_type = 'sampled-resonant';
    meta.component = resolved_component;
    meta.selected_component = selected_component;
    meta.psi_label = get_optional_field(gamma_definition, 'psi_label', '');
    meta.available = get_optional_field(gamma_definition, 'enabled', true);
end

function [psi_grid, gamma_values, meta] = reconstruct_sampled_true_gamma(gamma_definition, psi_grid, component)
    if isempty(component)
        resolved_component = 'true';
    else
        resolved_component = lower(strtrim(char(component)));
    end

    sample_psi = gamma_definition.psi_grid(:);
    if isempty(sample_psi)
        error('The exported true gamma structure does not contain psi_grid.');
    end

    switch resolved_component
        case {'true', 'selected', 'full'}
            sample_values = gamma_definition.gamma_true(:);
            resolved_component = 'true';
        case 'agent1'
            sample_values = gamma_definition.gamma_agent_1(:);
        case 'agent2'
            sample_values = gamma_definition.gamma_agent_2(:);
        otherwise
            error('Unsupported true gamma component: %s', resolved_component);
    end

    if isempty(psi_grid)
        psi_grid = sample_psi;
        gamma_values = sample_values;
    else
        validateattributes(psi_grid, {'numeric'}, {'vector', 'finite'}, mfilename, 'psi_grid');
        psi_grid = psi_grid(:);
        gamma_values = interp1(sample_psi, sample_values, psi_grid, 'linear', 'extrap');
    end

    meta = struct();
    meta.source_type = 'sampled-true';
    meta.component = resolved_component;
    meta.selected_component = 'true';
    meta.psi_label = '\psi';
    meta.available = get_optional_field(gamma_definition, 'available', true);
    meta.agent_id_1 = get_optional_field(gamma_definition, 'agent_id_1', []);
    meta.agent_id_2 = get_optional_field(gamma_definition, 'agent_id_2', []);
end

function psi_grid = get_default_psi_grid(gamma_definition, default_count)
    if isfield(gamma_definition, 'psi_grid_centered') && ~isempty(gamma_definition.psi_grid_centered)
        psi_grid = gamma_definition.psi_grid_centered(:);
        return;
    end
    if isfield(gamma_definition, 'psi_grid') && ~isempty(gamma_definition.psi_grid)
        psi_grid = gamma_definition.psi_grid(:);
        return;
    end
    psi_grid = linspace(-pi, pi, default_count).';
end

function component_name = get_selected_component_name(gamma_definition)
    component_name = 'full';
    if isfield(gamma_definition, 'component') && ~isempty(gamma_definition.component)
        component_name = lower(strtrim(char(gamma_definition.component)));
    end
end

function resolved_component = normalize_resonant_component(component, selected_component)
    if isempty(component)
        resolved_component = selected_component;
        return;
    end

    resolved_component = lower(strtrim(char(component)));
    if strcmp(resolved_component, 'selected')
        resolved_component = selected_component;
    end
end

function sample_values = get_sampled_resonant_values(gamma_definition, component)
    switch component
        case 'full'
            sample_values = get_optional_field(gamma_definition, 'gamma_values_full_centered', []);
            if isempty(sample_values)
                sample_values = get_optional_field(gamma_definition, 'gamma_values_full', []);
            end
        case 'symmetric'
            sample_values = get_optional_field(gamma_definition, 'gamma_values_symmetric_centered', []);
            if isempty(sample_values)
                sample_values = get_optional_field(gamma_definition, 'gamma_values_symmetric', []);
            end
        case 'antisymmetric'
            sample_values = get_optional_field(gamma_definition, 'gamma_values_antisymmetric_centered', []);
            if isempty(sample_values)
                sample_values = get_optional_field(gamma_definition, 'gamma_values_antisymmetric', []);
            end
        otherwise
            error('Unsupported resonant component: %s', component);
    end

    if isempty(sample_values)
        sample_values = get_optional_field(gamma_definition, 'gamma_values_centered', []);
        if isempty(sample_values)
            sample_values = get_optional_field(gamma_definition, 'gamma_values', []);
        end
    end
end

function value = get_optional_field(source_struct, field_name, default_value)
    value = default_value;
    if isstruct(source_struct) && isfield(source_struct, field_name)
        value = source_struct.(field_name);
    end
end
