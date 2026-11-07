# H1g (C22, exploratory): does the stress-kernel SHAPE decide the optimal schedule shape?
# s(t) = sum_j u_j [G(t - t_j) - G(t - t_{j+1})], G = unit-step response normalised so constant rate gives s(T) = 1.
# N = int_0^{T+post} (exp(Pi s) - 1) dt (Dieterich small-t limit, as H1f).  Volume: mean(u) = 1, 0 <= u <= umax.
import numpy as np, json, sys
from scipy.optimize import minimize
def kernels(theta):
    return {
      'exponential (first-order, saturating)': lambda t: np.where(t > 0, 1 - np.exp(-np.maximum(t, 0) / theta), 0.0),
      'radial diffusive (log, infinite-acting)': lambda t: np.where(t > 0, np.log1p(np.maximum(t, 0) / (0.02 * theta)), 0.0),
      'linear diffusive (sqrt, 1-D)': lambda t: np.where(t > 0, np.sqrt(np.maximum(t, 0)), 0.0),
    }
def solve(G, Pi, nc=12, post=0.5, umax=3.0):
    tj = np.linspace(0, 1, nc + 1); t = np.linspace(0, 1 + post, 1501)[1:]
    A = np.stack([G(t - tj[j]) - G(t - tj[j + 1]) for j in range(nc)], 1)
    A /= (A.sum(1)[np.argmin(abs(t - 1))])                       # constant rate -> s(T) = 1
    dt = t[1] - t[0]
    f = lambda u: np.sum(np.expm1(Pi * (A @ u))) * dt
    g = lambda u: Pi * dt * (A.T @ np.exp(Pi * (A @ u)))
    u0 = np.ones(nc)
    r = minimize(f, u0, jac=g, method='SLSQP', bounds=[(0, umax)] * nc,
                 constraints=[{'type': 'eq', 'fun': lambda u: u.mean() - 1, 'jac': lambda u: np.ones(nc) / nc}],
                 options={'maxiter': 500, 'ftol': 1e-14})
    s = A @ r.x
    return dict(gain=1 - r.fun / f(u0), u=np.round(r.x, 3).tolist(), t_peak_const=float(t[np.argmax(A @ u0)]))
out = {}
for Pi in (3.0, 10.0):
    for name, G in kernels(0.3).items():
        rr = solve(G, Pi); out[f'{name} | Pi {Pi}'] = rr
        print(f'Pi {Pi:4.1f} {name:42s} gain {100*rr["gain"]:6.2f} %  u = {rr["u"]}', flush=True)
json.dump(out, open('h1g_results.json', 'w'), indent=1)
