% phi の範囲を -pi から pi まで取得
phi = linspace(-pi, pi, 1000);

% 関数値
domega = 0.01*pi*5;
y = domega + 0.01*pi*(-10*sin(phi) - 3*sin(2*phi)+0.01);

% プロット
figure;
plot(phi, y, 'LineWidth', 1.5);
xlabel('$$\psi$$','Interpreter','latex');
ylabel('$$\dot{\psi}$$','Interpreter','latex');
xlim([-pi pi]);
xticks([-pi -pi/2 0 pi/2 pi]);
xticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
grid on;
% set y-axis ticks only at 0 and domega
ymin = min(y);
ymax = max(y);
% ensure domega is available as a tick
ticks = unique([0, domega]);
% make sure ticks fall within the plotted range; if not, expand limits
new_ymin = min([ymin, ticks]);
new_ymax = max([ymax, ticks]);
yticks(ticks);
% label 0 and domega using LaTeX
ytl = cell(size(ticks));
for ii = 1:numel(ticks)
    if abs(ticks(ii)) < 1e-12
        ytl{ii} = '0';
    else
        ytl{ii} = '$$\Delta\omega$$';
    end
end
yticklabels(ytl);
set(gca,'TickLabelInterpreter','latex');
% add padding to y-limits so ticks/markers don't touch the edges
yrange = new_ymax - new_ymin;
pad = max(0.02 * pi, 0.05 * yrange);
ylim([new_ymin - pad, new_ymax + pad]);
% plot horizontal y=0 line and mark zeros colored by derivative sign
hold on;
plot(phi, zeros(size(phi)), 'k-', 'LineWidth', 1);

% numerical derivative
dy = gradient(y, phi);

% find sign changes (zero crossings) between consecutive samples
cross_idx = find(y(1:end-1) .* y(2:end) <= 0);
root_phi = [];
root_d = [];
for ii = 1:numel(cross_idx)
	i = cross_idx(ii);
	% avoid degenerate case where both samples are zero
	if y(i) == 0 && y(i+1) == 0
		rp = phi(i);
	else
		% linear interpolation for root between phi(i) and phi(i+1)
		rp = phi(i) - y(i) * (phi(i+1) - phi(i)) / (y(i+1) - y(i));
	end
	root_phi(end+1) = rp; %#ok<AGROW>
	% interpolate derivative at root
	rd = interp1(phi, dy, rp, 'linear');
	root_d(end+1) = rd; %#ok<AGROW>
end

% plot markers: black filled when derivative > 0, white filled when < 0
for ii = 1:numel(root_phi)
	if root_d(ii) > 0
		plot(root_phi(ii), 0, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 0.8);
	else
		plot(root_phi(ii), 0, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 10, 'LineWidth', 0.8);
	end
end

% place arrows on y=0: one per interval between roots (N roots -> N+1 intervals)
edges = [-pi, sort(root_phi), pi];
% handle case with no roots
if isempty(root_phi)
	edges = [-pi, pi];
end
mid_phi = (edges(1:end-1) + edges(2:end)) / 2;
% arrow length (in x units)
arrow_len = 0.12 * pi;
for j = 1:numel(mid_phi)
	% evaluate function at interval midpoint
	val = interp1(phi, y, mid_phi(j), 'linear');
	if isempty(val) || isnan(val)
		val = 0;
	end
	% compute interval width and limit arrow length to slightly less than that
	interval_width = edges(j+1) - edges(j);
	max_allowed = 0.1 * interval_width;
	mag = min(arrow_len, max_allowed);
	% if interval is extremely small, fall back to a tiny visible arrow
	if mag <= 0
		mag = 0.02 * pi;
	end
	if val >= 0
		u = mag; % right
	else
		u = -mag; % left
	end
	% draw arrow on y=0 using fixed-size triangular head (data units)
	% fixed head size (in data x-units); adjust these constants as needed
	head_len_fixed = 0.04 * pi; % arrowhead length
	head_width_fixed = 0.02 * pi; % full arrowhead width

	dir = sign(u);
	head_tip = mid_phi(j) + u; % tip position
	head_base = head_tip - dir * head_len_fixed; % base of triangle
	% draw shaft from midpoint to base (avoid overlapping the head)
	plot([mid_phi(j), head_base], [0, 0], 'k-', 'LineWidth', 10);

	% triangle vertices (tip, base upper, base lower)
	tri_x = [head_tip, head_base, head_base];
	tri_y = [0,  head_width_fixed/2, -head_width_fixed/2];
	if dir < 0
		% for left-pointing arrow the tip is left of the base (tri still valid)
		% no further change needed other than ensuring correct orientation
	end
	patch(tri_x, tri_y, 'k', 'EdgeColor', 'k');
end

tuneFigure;
%saveFigure;