import matplotlib.pyplot as plt
try:
    from multifrailty import theory
except:
    import theory
import numpy as np
import seaborn as sb

sb.set("paper", font_scale=2)

k = 10
grgd_max = 20
grgd_steps = 1000
grgd_min = 1e-3
alpha = 1
fig_w, fig_h = (5.5, 5.5)

plt.figure(figsize=(fig_w, fig_h))
gr = 0.5
x = np.linspace(1e-3, 1, 1000)
plt.plot(x, gr * x, lw=3, color="blue", label="uncoupled (elastic)")
g = theory.u_factory(k)
alpha = 1
plt.xlabel("Damaged Fraction")
plt.ylabel("Recovery Fraction / unit time")
plt.tight_layout()
plt.savefig(f"stress_strain_uncoupled_k{k}.png")
plt.savefig(f"stress_strain_uncoupled_k{k}.pdf")

mu = np.vectorize(lambda p: (1 - alpha * g(1 - p)))
grx = lambda gr, x: gr * x * mu(x)
plt.plot(x, grx(gr, x), lw=3, color="orange", label=f"coupled (inelastic)")
plt.axis("auto")

if False:
    for alpha in [1, .9, .6, .3, .01]:
        g = theory.u_factory(k)
        mu = np.vectorize(lambda p: (1 - alpha * g(1 - p)))
        grx = lambda gr, x: gr * x * mu(x)
        plt.plot(x,
                 grx(gr, x),
                 lw=3,
                 label=f"coupled (inelastic) $\\alpha={alpha}$")
        #plt.plot(x, grx(gr,x),lw=3,color="orange",label="coupled (inelastic)")
        plt.axis([0, 0.3, 0, 0.3])
plt.xlabel("Damaged Fraction")
plt.ylabel("Recovery Fraction / unit time")
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(f"stress_strain_theory_k{k}.png")
plt.savefig(f"stress_strain_theory_k{k}.pdf")
plt.close()
