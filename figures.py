import pathlib

import matplotlib.pyplot as plt
import numpy as np

from nanoconductance import (Electrolyte, conductance, conductance_model, crossover_concentration, gouy_chapman, solver, wall_field)

pathlib.Path("figures").mkdir(exist_ok=True)
SIGMA = -0.05     
W= 50e-6
L = 5e-3

def fig_double_layer():
    electrolyte = Electrolyte(1e-3)
    lambda_D = electrolyte.debye_length
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), constrained_layout=True)

    for H in [0.5, 2, 10, 50]:
        pb = solver(H * lambda_D, -0.02, electrolyte)
        x = pb.z / (pb.h / 2)
        line, = ax1.semilogy(x, pb.n_plus / electrolyte.n_0, label=f"$h/\\lambda_D$ = {H}")
        ax1.semilogy(x, pb.n_minus / electrolyte.n_0, "--", color=line.get_color())
    ax1.set(xlabel="$2z/h$ (0 = centre, 1 = paroi)", ylabel=r"$n_\pm/n_0$", title = "Counter-ions (solid line) and co-ions (dashed line)")
    ax1.legend()

    h = 60 * lambda_D
    pb = solver(h, SIGMA, electrolyte)
    x = (h / 2 - pb.z) / lambda_D
    ax2.plot(x, pb.psi, lw=3, alpha=0.6, label="Numerical Poisson-Boltzmann")
    ax2.plot(x, gouy_chapman(x, wall_field(SIGMA, electrolyte)), "k--",
             label="Gouy-Chapman (exact)")
    ax2.set(xlim=(0, 6), xlabel=r"distance à la paroi / $\lambda_D$",
            ylabel=r"$e\phi/k_BT$", title="Validation : isolated wall")
    ax2.legend()
    fig.savefig("figures/fig1_double_layer.png", dpi=200)



def fig_conductance():
    cs = np.logspace(-5, 0, 26)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2), constrained_layout=True)

    for h in [50e-9, 200e-9, 1000e-9]:
        G = np.array([conductance(h, SIGMA, Electrolyte(c), W, L) for c in cs])
        G_vol, G_surf = np.array(
            [conductance_model(h, SIGMA, Electrolyte(c), W, L) for c in cs]).T
        line, = ax1.loglog(cs, 1e9 * G, "o", ms=4, label=f"h = {h*1e9:g} nm")
        color = line.get_color()
        ax1.loglog(cs, 1e9 * (G_vol + G_surf), "-", color=color, lw=1)
        ax1.axvline(crossover_concentration(h, SIGMA), color=color, ls=":", lw=1)

        l_du = np.array([Electrolyte(c).dukhin_length(SIGMA) for c in cs])
        ax2.loglog(h / l_du, G / G_surf, "o", ms=4, color=color)

    ax1.axhline(1e9 * G_surf[0], color="k", ls="--", lw=1,
                label=r"plateau $2\mu_+|\Sigma|\,w/L$")
    ax1.set(xlabel="c (mol/L)", ylabel="G (nS)",
            title="Dots : Poisson-Boltzmann ; lines : $G_{vol}+G_{surf}$")
    ax1.legend()

    el = Electrolyte(1e-3)
    k = (el.mu_plus + el.mu_minus) / (2 * el.mu_plus)
    x = np.logspace(-3.5, 3.5, 200)
    ax2.loglog(x, 1 + k * x, "k-", lw=1, label=r"$1 + \frac{\mu_+ + \mu_-}{2\mu_+}\,h/\ell_{Du}$")
    ax2.set(xlabel=r"$h/\ell_{Du}$", ylabel=r"$G/G_{surf}$",
            title="Master curve")
    ax2.legend()
    fig.savefig("figures/fig2_conductance.png", dpi=200)


fig_double_layer()
fig_conductance()