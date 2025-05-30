% compute_xi.m
% Computes the stability parameter ξ and optimal T from α²F(ω)

clear; clc;

% Load alpha^2 F(omega) data
data = load('a2F_Hg_0GPa.txt');  % <-- change filename as needed
omega = data(:,1);
a2F = data(:,2);

omegaMin = min(omega);
omegaMax = max(omega);

% Interpolate α²F(ω)
a2F_interp = @(w) interp1(omega, a2F, w, 'linear', 0);

% Define g(x)
g = @(x) 6.*x.^2 + ...
         12.*x.^3 .* imag(psi(1, 1i.*x)) + ...
         6.*x.^4 .* real(psi(2, 1i.*x));

% Tabulate g(x)
x_tab = linspace(1e-6, 20, 20000);
g_tab = g(x_tab);
g_interp = @(x) interp1(x_tab, g_tab, x, 'pchip', 0);

% Define integrand
integrand = @(omega, T) g_interp(omega ./ (2 * pi * T)) .* ...
                        (2 .* a2F_interp(omega) ./ omega);

% Log-linear adaptive ω splitting
omegaSplit = unique([logspace(log10(omegaMin), log10(0.003), 5), ...
                     linspace(0.003, omegaMax, 4)]);

% Define xi(T) via integration
xi_integral = @(T) sum(arrayfun(@(i) ...
    integral(@(w) integrand(w, T), omegaSplit(i), omegaSplit(i+1), ...
    'RelTol',1e-4,'AbsTol',1e-10), 1:length(omegaSplit)-1));

% Maximize ξ(T)
Tmin = 1e-5;
Tmax = omegaMax / 2;
obj = @(T) -xi_integral(T);
[Topt, negXiMax] = fminbnd(obj, Tmin, Tmax);
xiMax = -negXiMax;

% Output results
fprintf('Optimal temperature (frequency units): %.6f\n', Topt);
fprintf('Stability parameter ξ: %.6f\n', xiMax);
