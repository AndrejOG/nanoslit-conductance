import numpy as np
from dataclasses import dataclass
from scipy.integrate import solve_bvp, trapezoid

CHARGE = 1.602176634e-19
K_B = 1.380649e-23
N_A = 6.02214076e23
EPSILON_0 = 8.88541878128e-12


@dataclass(frozen=True)
class Electrolyte:
    c_molar: float
    T: float = 293.15 
    eps_r: float = 80.4 
    D_plus: float = 1.96e-9
    D_minus: float = 2.03e-9

    @property
    def n_0(self): 
        return 1000 * self.c_molar * N_A

    @property 
    def debye_length(self):
        return np.sqrt(self.eps_r*EPSILON_0*K_B*self.T/(2*(CHARGE**2)*self.n_0))

    @property 
    def bjerrum_length(self):
       return (CHARGE**2)/(4*np.pi*self.eps_r*EPSILON_0*K_B*self.T)

    def dukhin_length(self, sigma):
        return abs(sigma) / (CHARGE * self.n_0)

    @property
    def mu_plus(self):
        return CHARGE * self.D_plus / (K_B * self.T)

    @property
    def mu_minus(self):
            return CHARGE * self.D_minus / (K_B * self.T)
    

def wall_field(sigma, electrolyte): 
    return (4*np.pi*electrolyte.bjerrum_length*electrolyte.debye_length*sigma)/CHARGE

def gouy_chapman(z, s):
    return 4 * np.arctanh(np.tanh((2*np.arcsinh(s/2))/4)*np.exp(-z))

def grid(a, n, p=3):
    t = np.linspace(0.0, 1.0, n)
    return a * (1 - (1 - t)**p)

@dataclass
class Result:
    z: np.ndarray
    psi: np.ndarray
    n_plus: np.ndarray
    n_minus: np.ndarray
    h: float
    sigma: float

def solver(h, sigma, electrolyte, tol=1e-6, max_nodes=200_000, n=2000):
    lambda_D = electrolyte.debye_length
    a = h / (2 * lambda_D)
    s = wall_field(sigma, electrolyte)

    def fun(Z, y):
        return np.vstack([y[1], np.sinh(y[0])])

    def make_bc(s):
        return lambda ya, yb: np.array([ya[1], yb[1] - s])

    Z0 = grid(a, 500)
    if a < 2:
        guess = np.arcsinh(s / a) + (Z0**2 * s)/(2*a)
    else:
        guess = gouy_chapman(a - Z0, s) + gouy_chapman(a + Z0, s)
    y_guess = np.vstack([guess, np.gradient(guess, Z0)])
    sol = solve_bvp(fun, make_bc(s), Z0, y_guess, tol=tol, max_nodes=max_nodes)

    if not sol.success:
        x, y = Z0, np.zeros((2, Z0.size))
        for s_i in np.linspace(0.05, 1, 20) * s:
            sol = solve_bvp(fun, make_bc(s_i), x, y, tol=tol, max_nodes=max_nodes)
            x, y = sol.x, sol.y
        if not sol.success:
            raise RuntimeError(sol.message)

    Z = grid(a, n)
    psi = sol.sol(Z)[0]
    return Result(z=Z*lambda_D, psi=psi, n_plus=electrolyte.n_0 * np.exp(-psi), n_minus=electrolyte.n_0 * np.exp(psi), h=h, sigma=sigma)

def conductance(h, sigma, electrolyte, w, L):
    pb = solver(h, sigma, electrolyte)
    j = CHARGE * (electrolyte.mu_plus*pb.n_plus + electrolyte.mu_minus*pb.n_minus)
    return 2*(w/L)*trapezoid(j, pb.z) 

def conductance_model(h, sigma, electrolyte, w, L): 
    mu_counter = electrolyte.mu_plus if sigma < 0 else electrolyte.mu_minus
    G_vol = w / L * CHARGE * (electrolyte.mu_plus + electrolyte.mu_minus) * electrolyte.n_0 * h
    G_surf = w / L * 2 * mu_counter * abs(sigma)
    return G_vol, G_surf

def crossover_concentration(h, sigma):
    return abs(sigma) / (CHARGE * h) / (1e3 * N_A)


        


    




