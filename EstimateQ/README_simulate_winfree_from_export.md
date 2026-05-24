# simulate_winfree_from_export

This adds a Winfree-type two-phase simulator that reconstructs both:

- `s_2(phi1, phi2)` from `gamma_export_latest.mat`
- `z_opt(theta)` from the exported `s_2` basis by maximizing
  `R(psi_plus) = int_0^{2pi} W(theta; psi_plus)^2 dtheta`, with
  `W(theta) = s_2(theta-psi_plus, theta) - s_2(theta+psi_plus, theta)`

The simulated model is:

- `dphi1/dt = omega1 + sigma * z(phi1) * s(phi1, phi2)`
- `dphi2/dt = omega2 + sigma * z(phi2) * s(phi2, phi1)`

## Files

- `simulate_winfree_from_export.m`: main function.
- `run_simulate_winfree_from_export.m`: quick runner with plots.

## Quick start

From MATLAB, in `EstimateQ`:

```matlab
run_simulate_winfree_from_export
```

## Main options

Configure through `opts` in `simulate_winfree_from_export`:

- `mat_path`: optional export file path. Empty means auto-pick latest.
- `signal_role`: `'a2'` or `'derived'`.
- `agent_mode`: `'phase_id_2'` (default), `'phase_id_1'`, or `'agent_id'`.
- `agent_id`: used only when `agent_mode='agent_id'`.
- `omega`: `[omega1, omega2]`.
- `sigma`: coupling strength.
- `tspan`: integration range.
- `phi0`: initial phase values.

## Notes

- The reconstructed `s_2` includes exported `z_mean` if present.
- `z_opt` is optionally power-normalized to match unit sine-wave power.
- Interpolation for `z_opt(theta)` is periodic on `[0, 2*pi]`.
