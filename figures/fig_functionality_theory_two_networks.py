try:
    from multifrailty import theory, simtools
except:
    import theory, simtools
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb
from multiprocessing import Pool
from functools import partial

sb.set("paper", font_scale=1.9)

k_vec = (4, 8)
grgd_max = 4
grgd_steps = 200
grgd_min = 1e-3
alpha = 1
fig_w, fig_h = (8.5, 8.5)
def get_mu_from_g(g_):
    return np.vectorize(lambda p: (1 - g_(1 - p)))


if __name__ == "__main__":
    grgd = np.linspace(grgd_min, grgd_max, grgd_steps)
    to_do = partial(theory.pinf_two_networks, k=k_vec, solpoints=200)
    inputs = zip(grgd, grgd)

    # sols = list(map(to_do, inputs))
    with Pool() as p:
        sols = p.map(to_do, inputs)
    
    split_sols = theory.pack_two_network_solutions_into_nan_arrays(sols)

    g = map(theory.u_factory, k_vec)
    mu = list(map(get_mu_from_g, g))

    plt.figure(figsize=(fig_w, fig_h))
    plt.xlabel(r"Repair Rate $(\gamma_{r,0}/\gamma_d)$")
    plt.ylabel("Functional Fraction")
    colors = ["C1", "C2"]
    for mu_, sols_, k_, c_ in zip(mu, split_sols, k_vec, colors):
        plt.plot(grgd, mu_(sols_[:, 0]), color=c_, lw=3, label=f"k = {k_}")
        plt.plot(grgd, mu_(sols_[:, 1]), color=c_, lw=3, ls="--")
        plt.plot(grgd, mu_(sols_[:, 2]), color=c_, lw=3)

    plt.legend(loc="best")
    plt.tight_layout()
    plt.savefig(f"coupled_different_k{k_vec}.png", transparent=False)
    plt.savefig(f"coupled_different_k{k_vec}.pdf", transparent=False)
