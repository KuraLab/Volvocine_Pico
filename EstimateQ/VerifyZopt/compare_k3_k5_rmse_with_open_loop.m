function out = compare_k3_k5_rmse_with_open_loop(varargin)
% COMPARE_K3_K5_RMSE_WITH_OPEN_LOOP Compare high-gain s2 RMSE vs detuning.
%   out = compare_k3_k5_rmse_with_open_loop('ShowFigures', false);
%
% Uses all numeric frequency folders in Spring3/k3w2 and Spring3/k5w2.
% Same 240_0531 open-loop model, full-file p10/p90 normalization and clipping,
% and feedback-time 0--115 s interval as the k=1 analysis. Delta omega is
% omega1-omega2. Results are saved separately from k=1 and prior subset runs.
% Six Module 1/2 zoom overlays (0--5 and 60--65 s) are also saved for each
% gain at Delta omega=0.2*pi,0.3*pi,0.5*pi, under closed_loop_timeseries/overlaid_zoom
% in separate k{gain}_delta_omega_{0p2,0p3,0p5}pi folders (PNG and PDF).
% To export from the existing analysis without refitting:
%   out = compare_k3_k5_rmse_with_open_loop('ReplotOnly',true,'ShowFigures',false);
% SaveZoomPanels=false disables these additional exports.
% The dotphi12 panels compare the logged phase slopes with rates from predicted s.
% ShowReconstructedRate=false hides the rates computed from predicted s.
% Other name/value options are forwarded to the common comparison function.

    here = fileparts(mfilename('fullpath'));
    out = compare_delta0_rmse_with_open_loop('Conditions', {'k3w2','k5w2'}, ...
        'HighGainDeltaOmegaOverPi', [], ...
        'OutputDir', fullfile(here, 'k3_k5_all_detunings_rmse_240'), varargin{:});
end
