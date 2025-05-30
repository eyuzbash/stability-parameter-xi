# compute_xi.py
# Computes the stability parameter ξ from an Eliashberg function a2F(ω)

import numpy as np
from scipy import interpolate, integrate, optimize
import mpmath

# Set mpmath precision for complex polygamma
mpmath.mp.dps = 25

# Load α²F(ω) data
filename = "data/a2F_Hg_0GPa.txt"  # <-- update path if needed
data = np.loadtxt(filename)
omega_vals = data[:, 0]
a2F_vals = data[:, 1]

# Interpolate α²F(ω)
a2F_interp = interpolate.interp1d(omega_vals, a2F_vals, kind='linear',
                                  bounds_error=False, fill_value=0.0)
omega_min = omega_vals.min()
omega_max = omega_vals.max()

# Tabulate g(x) using mpmath for complex-valued polygamma
x_tab = np.linspace(1e-6, 20, 20000)
g_tab = np.array([
    float(6 * x**2 +
          12 * x**3 * mpmath.im(mpmath.polygamma(1, 1j * x)) +
          6 * x**4 * mpmath.re(mpmath.polygamma(2, 1j * x)))
    for x in x_tab
])
g_interp = interpolate.interp1d(x_tab, g_tab, kind='cubic',
                                bounds_error=False, fill_value=0.0)

# Define the integrand
def integrand(omega, T):
    if T <= 0:
        return 0.0
    x = omega / (2 * np.pi * T)
    return g_interp(x) * (2 * a2F_interp(omega) / omega)

# Log-linear adaptive split of ω domain
log_part = np.logspace(np.log10(omega_min), np.log10(0.003), 5)
lin_part = np.linspace(0.003, omega_max, 4)
split_points = np.unique(np.concatenate([log_part, lin_part]))

# Define ξ(T)
def xi_integral(T):
    if T <= 0:
        return 0.0
    result = 0.0
    for i in range(len(split_points) - 1):
        a = split_points[i]
        b = split_points[i + 1]
        subres, _ = integrate.quad(integrand, a, b, args=(T,), limit=200)
        result += subres
    return result

# Maximize ξ(T)
Tmin = 1e-5
Tmax = omega_max / 2
opt_result = optimize.minimize_scalar(lambda T: -xi_integral(T),
                                      bounds=(Tmin, Tmax), method='bounded')
T_opt = opt_result.x
xi_max = -opt_result.fun

# Output results
print(f"Optimal temperature (frequency units): {T_opt:.6f}")
print(f"Stability parameter ξ: {xi_max:.6f}")
