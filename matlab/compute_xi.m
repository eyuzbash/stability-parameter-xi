% compute_xi.m
% Computes the stability parameter ξ and optimal T from α²F(ω)

clear; clc;

% Load alpha^2 F(omega) data
data = load('');  % <-- change filename as needed
omega = data(:,1);
a2F = data(:,2);

% Remove duplicates by averaging y-values
[unique_omega, ~, idx] = unique(omega);
unique_a2F = accumarray(idx, a2F, [], @mean);

omega = unique_omega;
a2F = unique_a2F;

omegaMin = min(omega);
omegaMax = max(omega);

% Interpolate α²F(ω)
a2F_interp = @(w) interp1(omega, a2F, w, 'linear', 0);

% Define g(x)
g = @(x) 6.*x.^2 + 12.*x.^3 .*imag(hurwitzZeta(2, 1i.*x)) + 6.*x.^4 .* real(-2*hurwitzZeta(3, 1i.*x));

% Tabulate g(x)
x_tab = linspace(1e-6, 20, 2000); % <- density of the interpolation and integration mesh
g_tab = g(x_tab);
g_interp = @(x) interp1(x_tab, g_tab, x, 'pchip', 0);

% Define integrand
integrand = @(omega, T) g_interp(omega ./ (2 * pi * T)) .* ...
                        (2 .* a2F_interp(omega) ./ omega);

% Log-linear adaptive ω splitting
omegaSplit = unique([logspace(log10(omegaMin), log10(0.003), 10), ...
                     linspace(0.003, omegaMax, 10)]);

% Define xi(T) via integration
xi_integral = @(T) sum(arrayfun(@(i) ...
    integral(@(w) integrand(w, T), omegaSplit(i), omegaSplit(i+1), ...
    'RelTol',1e-4,'AbsTol',1e-10), 1:length(omegaSplit)-1));

% Maximize ξ(T)
Tmin = 1e-5;
Tmax = omegaMax;
obj = @(T) -xi_integral(T);
[Topt, negXiMax] = fminbnd(obj, Tmin, Tmax);
xiMax = -negXiMax;

% Output results
fprintf('Optimal temperature (frequency units): %.6f\n', Topt);
fprintf('Stability parameter ξ: %.6f\n', xiMax);

beep;
