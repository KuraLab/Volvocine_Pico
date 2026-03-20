S = load('EstimateQ/Spring1/255/gamma_exports/gamma_export_latest.mat', 'gamma_export');
G = S.gamma_export;
phase_ids = G.phase_agent_ids(:).';
zid = G.z_agent_id;
fprintf('phase_agent_ids: [%d %d]\n', phase_ids(1), phase_ids(2));
fprintf('z_agent_id: %d\n', zid);
if zid == phase_ids(1)
    fprintf('selected_label: s_1\n');
elseif zid == phase_ids(2)
    fprintf('selected_label: s_2\n');
else
    fprintf('selected_label: s_{agent_id=%d} (not phase index 1 or 2)\n', zid);
end
