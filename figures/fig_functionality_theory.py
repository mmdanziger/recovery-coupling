try:
    from multifrailty import theory, simtools
except:
    import theory, simtools
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sb

sb.set("paper", font_scale=1.9)

k = 8
grgd_max = 3
grgd_steps = 1000
grgd_min = 1e-3
alpha = 1
fig_w, fig_h = (8.5, 8.5)
grgd_vec = np.linspace(grgd_min, grgd_max, grgd_steps)
sols_1 = theory.pinf_solutions_vec(alpha_i=alpha, k=k, grgd_vec=grgd_vec)
sols_0 = theory.pinf_vec(grgd_vec)
g = theory.u_factory(k)
mu = np.vectorize(lambda p: (1 - g(1 - p)))
plt.figure(figsize=(fig_w, fig_h))

plt.plot(grgd_vec, mu(sols_0), color="blue", lw=3, label="uncoupled")
plt.xlabel(r"Repair Rate $(\gamma_{r,0}/\gamma_d)$")
plt.ylabel("Functional Fraction")
plt.tight_layout()
plt.savefig(f"uncoupled_{k}.png", transparent=False)

plt.plot(grgd_vec, mu(sols_1[:, 0]), color="orange", lw=3, label="coupled")
plt.plot(grgd_vec, mu(sols_1[:, 1]), color="orange", ls="--", lw=3)
plt.plot(grgd_vec, mu(sols_1[:, 2]), color="orange", lw=3)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(f"coupled_uncoupled_comparison_{k}.png", transparent=False)
plt.savefig(f"coupled_uncoupled_comparison_{k}.pdf", transparent=False)
