"""Module for calculating theoretical solutions for recovery coupled networks."""
from typing import Callable, List, Tuple
import numpy as np
from scipy.optimize import fsolve, root
from scipy.special import zeta
from tqdm import tqdm
import warnings
import pychebfun
from functools import partial
from itertools import product as iproduct
from operator import itemgetter

methods = [
    'hybr', 'lm', 'broyden1', 'broyden2', 'anderson', 'linearmixing',
    'excitingmixing', 'krylov', 'df-sane'
]

avoid_chebfun = True


def P_inf_factory(k):
    '''
    ER network percolation function returns the function, not the value
    :param k: average degree
    :return: percolation profile function
    '''
    def f(p):
        if np.isnan(p):
            return np.nan
        return fsolve(lambda x: x - p * (1 - np.exp(-k * x)),
                      0.5,
                      full_output=False)

    return f


def u_factory(k: float) -> Callable[[float], float]:
    """A factory to return u(p) functions for ER networks.

    Args:
        k (float): the average degree

    Returns:
        Callable[[float], float]: a function that is the solution of the
            ER self-consistency condition
    """
    def f(p: float) -> float:
        def ER_u_func(x: float) -> float:
            return (1 - p) + p * (np.exp(k * (x - 1))) - x

        x, info, ierr, message = fsolve(ER_u_func, 0.5, full_output=True)
        if ierr == 1:
            return x
        else:
            return np.nan

    return f


def u_factory_sf(alpha):
    from mpmath.fp import polylog

    z = zeta(alpha - 1)

    def f(p):
        return fsolve(lambda x: (1 - p) + p * polylog(alpha - 1, x) /
                      (x * z) - x,
                      0.5,
                      full_output=False)

    return f


def generate_p_sequence(gr=0.5, gd=0.01, k=3, tmax=100, init=None):
    g = P_inf_factory(k)
    plist = []
    if not init:
        init = gd
    plist.append(init)

    while len(plist) < tmax:
        # plist.append(plist[-1] * (1 - gr))
        # plist.append(plist[-1] + (1 - plist[-1]) * gd)
        plist.append(plist[-1] + (1 - plist[-1]) * gd - plist[-1] * gr)
    return np.array(plist)


pinf = np.vectorize(lambda gr, gd: 1 / (1 - gr * (1 - 1 / gd)))

pinf_vec = np.vectorize(lambda grgd: 1 / (1 + grgd))

# @np.vectorize
# def pinf_coupled(gr, gd, k=3, solpoints=10, eps=1e-5, method="hybr"):
#     g = u_factory(k)
#     Pinf = P_inf_factory(k)
#     mu = lambda p: (1 - g(1 - p))
#     # to_solve = lambda x: 1 / (1 - gr * mu(x) * (1 - 1 / gd)) - x
#     to_solve = lambda x: 1 / (1 + (gr * mu(x) / gd)) - x
#     return get_all_sols(to_solve)


@np.vectorize
def pinf_coupled(grgd: float,
                 k: float = 3,
                 alpha_i: float = 1,
                 solpoints: int = 10,
                 eps: float = 1e-5,
                 method: str = "hybr"):
    """Solve the case where two symmetric ER networks are recovery coupled to each other.

    Args:
        grgd (float): ratio of gamma_r to gamma_d
        k (float, optional): degree of network. Defaults to 3.
        alpha_i (float, optional): coupling strength. Defaults to 1.
        solpoints (int, optional): number of guess points to generate for solver. Defaults to 10.
        eps (float, optional): epsilon within which two solutions will be considered identical. Defaults to 1e-5.
        method (str, optional): method to use in `root`. Defaults to "hybr".

    Returns:
        np.ndarray: [description]
    """
    g = u_factory(k)
    mu = lambda p: (1 - alpha_i * g(1 - p))
    to_solve = lambda x: 1 / (1 + (grgd * mu(x))) - x
    return get_all_sols(
        to_solve,
        eps=eps,
        method=method,
        solpoints=solpoints,
    )


def pinf_two_networks(grgd: Tuple[float, float],
                      k: Tuple[float, float] = (3, 3),
                      alpha_i: Tuple[float, float] = (1, 1),
                      solpoints: int = 10,
                      eps: float = 1e-5,
                      method: str = "hybr"):
    """Find the fixed points for two recovery coupled ER networks (not-symmetric)

    Args:
        grgd (Tuple[float, float]): gamma_r / gamma_d ratio in each network
        k (Tuple[float, float], optional): avg degree in each network. Defaults to (3, 3).
        alpha_i (Tuple[float, float], optional): coupling strength in each network. Defaults to (1, 1).
        solpoints (int, optional): number of guesses to feed solver. Defaults to 10.
        eps (float, optional): precision of solution. Defaults to 1e-5.
        method (str, optional): method to pass to `root`. Defaults to "hybr".

    Returns:
        List[np.ndarray]: a list of all solutions found
    """
    g = list(map(u_factory, k))
    mu = lambda p: (1 - alpha_i[0] * g[0]
                    (1 - p)), lambda p: (1 - alpha_i[1] * g[1](1 - p))

    def two_networks_self_consistent(f1f2):
        cond1 = 1 / (1 + (grgd[0] * mu[1](f1f2[1]))) - f1f2[0]
        cond2 = 1 / (1 + (grgd[1] * mu[0](f1f2[0]))) - f1f2[1]
        return np.array([cond1, cond2], dtype=float).squeeze()

    return get_all_sols_two_networks(
        two_networks_self_consistent,
        eps=eps,
        method=method,
        solpoints=solpoints,
    )


def get_all_sols_two_networks(to_solve: Callable[[Tuple[float, float]],
                                                 Tuple[float, float]],
                              solpoints: int = 10,
                              eps: float = 1e-5,
                              method: str = "hybr",
                              bounds: list = [0, 1]) -> np.ndarray:
    """Solve `to_solve` returning all unique solutions within `bounds`.

    We need to find all unique solutions in the bounds of interest. By default, 
    `root` will only find a single solution. This way we scan the whole space and supply
    many guesses, find many solutions which we then limit to the unique ones up to a given
    precision.

    Args:
        to_solve (Callable): the function to solve
        solpoints (int, optional): the number of points to guess. Defaults to 10.
        eps (float, optional): the precision of solutions. Defaults to 1e-5.
        method (str, optional): method to pass to `root`. Defaults to "hybr".
        bounds (list, optional): bounds within which to search. Defaults to [0, 1].

    Returns:
        np.ndarray: array of floats representing solutions, expected to be of length 1 or 3
    """

    minsolguess, maxsolguess = bounds[0] + eps, bounds[1] - eps
    solpoints_1d = 20
    guesses = get_list_of_guesses_for_solver(minsolguess,
                                             maxsolguess,
                                             solpoints=solpoints_1d)
    gen = np.random.default_rng()
    all_guesses_2d = list(iproduct(guesses, guesses))
    if solpoints is None:
        guesses_2d = all_guesses_2d
    else:
        guesses_2d = gen.choice(all_guesses_2d, size=solpoints, replace=False)
    sols = get_solutions_using_root(to_solve, method, guesses_2d)

    return limit_to_unique_sols(sols)


def limit_to_unique_sols(sols):
    def stringify(arr, sep=","):
        return sep.join(["%.5f" % i for i in arr])

    def destringify(arrstring, sep=","):
        return np.array([float(i) for i in arrstring.split(sep)], dtype=float)

    unique_sols = list(map(destringify, set(map(stringify, sols))))
    return unique_sols


def pinf_solutions_vec(alpha_i: float, k: float,
                       grgd_vec: np.ndarray) -> np.ndarray:
    """Generate a vector of recovery coupled network fixed points for a symmetric system of ER networks.

    The solutions are in units of fraction of nodes disabled. So a solution
    of 1 means total collapse. Every value of grgd (ratio of gamma_r to gamma_d)
    has three solutions in the returned array. If one of them is not real the value is NaN.

    Args:
        alpha_i (float): coupling strength alpha
        k (float): ER network avg degree
        grgd_vec (np.ndarray): values of grgd to check 

    Returns:
        np.ndarray: solutions of shape (len(grgd_vec), 3)
    """
    partial_pinf = partial(pinf_coupled, alpha_i=alpha_i, k=k)
    solvec = list(map(partial_pinf, tqdm(grgd_vec)))
    sols = pack_solutions_into_nan_array(solvec)
    return sols


def pack_solutions_into_nan_array(solvec: List[List[float]]) -> np.ndarray:
    """Pack solutions into sorted array that has nans where no solutions are found.

    We have up to three solutions so the first column should be the first solution
    the second and third columns, the others if they exist. The first solution is 
    the always true 1 solution (total collapse)

    Args:
        solvec (List[List[float]]): list of solutions for each point. Each solution
           is itself a list of up to three solutions   

    Returns:
        np.ndarray: a sorted array of solutions as described above
    """

    sols = np.nan * np.ones((len(solvec), 3))
    for i, r in enumerate(solvec):
        for j, v in enumerate(sorted(r, reverse=True)):
            sols[i, j] = v
    return sols


def pack_two_network_solutions_into_nan_arrays(
        sols: List[List[np.ndarray]]) -> Tuple[np.ndarray, np.ndarray]:
    """Pack a list of two network solutions into two arrays with NaNs where no solution found.

    Returns two outputs that follow the same format as `pack_solutions_into_nan_array`.

    Args:
        sols (List[List[np.ndarray]]): the solutions of the two network case

    Returns:
        List[np.ndarray]: two arrays each with three columns and len(sols) rows
    """
    output = []
    for i in [0, 1]:

        def extract_element_i(sol):
            return sorted(map(itemgetter(i), sol), reverse=True)

        this_sols = list(map(extract_element_i, sols))
        output.append(pack_solutions_into_nan_array(this_sols))

    return output[0], output[1]


def get_grgd_c(k: float = 3):
    from scipy.misc import derivative
    g = u_factory(k)
    Pinf = P_inf_factory(k)
    mu = lambda p: (1 - g(1 - p))
    utag = lambda p: derivative(g, p, dx=1e-5)
    # to_solve = lambda x: 1 / (1 - gr * mu(x) * (1 - 1 / gd)) - x
    to_solve1 = lambda grgd, x: 1 / (1 + (grgd * mu(x))) - x
    to_solve2 = lambda grgd, x: grgd * utag(1 - x) / (1 +
                                                      (grgd * mu(x)))**2 + 1
    to_solve2d = lambda grgd_x: np.array(
        [to_solve1(grgd_x[0], grgd_x[1]),
         to_solve2(grgd_x[0], grgd_x[1])]).reshape(-1)
    return fsolve(to_solve2d, np.array([0.5, 0.5]).reshape(-1))


@np.vectorize
def pinf_repair_capacity(gr, gd, beta):
    beta /= 2
    Nc = 1 / beta

    #p = (-gd - gr + np.sqrt(4 * beta * gd * gr ** 2 + (gd + gr) ** 2)) / (2 * beta * gr ** 2)
    #p = ((-Nc) * gd + Nc * gr + gd * gr - np.sqrt(4 * Nc * gd ** 2 * gr + (Nc*gd-Nc*gr-gd*gr)**2))/2 * gd * gr)
    b = Nc / gr + Nc / gd - 1
    p = (-b + np.sqrt(b * b + 4 * Nc / gr)) / 2
    return p


def pinf_repair_capacity_coupled(gr0,
                                 gd,
                                 beta,
                                 k=3,
                                 solpoints=10,
                                 eps=1e-5,
                                 method="hybr"):
    g = u_factory(k)
    Pinf = P_inf_factory(k)
    mu = lambda p: (1 - g(1 - p))
    # to_solve = lambda x: 1 / (1 - gr * mu(x) * (1 - 1 / gd)) - x
    gr = lambda x: gr0 / (1 + gr0 * x * beta)
    to_solve = lambda x: 1 / (1 + (gr(x) * mu(x) / gd)) - x
    return get_all_sols(to_solve)


def get_all_sols(to_solve: Callable,
                 solpoints: int = 10,
                 eps: float = 1e-5,
                 method: str = "hybr",
                 bounds: list = [0, 1]) -> np.ndarray:
    """Solve `to_solve` returning all unique solutions within `bounds`.

    We need to find all unique solutions in the bounds of interest. By default, 
    `root` will only find a single solution. This way we scan the whole space and supply
    many guesses, find many solutions which we then limit to the unique ones up to a given
    precision.

    Args:
        to_solve (Callable): the function to solve
        solpoints (int, optional): the number of points to guess. Defaults to 10.
        eps (float, optional): the precision of solutions. Defaults to 1e-5.
        method (str, optional): method to pass to `root`. Defaults to "hybr".
        bounds (list, optional): bounds within which to search. Defaults to [0, 1].

    Returns:
        np.ndarray: array of floats representing solutions, expected to be of length 1 or 3
    """

    minsolguess, maxsolguess = find_solution_bounds(to_solve, method)
    if maxsolguess - minsolguess < eps:
        return np.array([0.5 * (minsolguess + maxsolguess)])

    sols = get_solutions(to_solve, method, minsolguess, maxsolguess, solpoints)
    sols = [
        float(i) for i in np.unique(["%.6f" % j for j in sols])
        if bounds[0] <= float(i) <= bounds[1]
    ]

    return np.array(sols).astype(float).reshape(-1)


def find_solution_bounds(to_solve, method):
    minsolguess = 0.00001
    maxsolguess = 0.99999
    sol0 = root(to_solve, minsolguess, method=method)
    sol1 = root(to_solve, maxsolguess, method=method)
    if sol0.success and sol1.success:
        minsolguess = sol0.x
        maxsolguess = sol1.x
    return minsolguess, maxsolguess


def get_solutions(to_solve: Callable[[float],
                                     float], method: str, minsolguess: float,
                  maxsolguess: float, solpoints: int) -> List[float]:
    """Get all solutions for a given callable by generating many guesses and checking them all.

    Args:
        to_solve (Callable): the function to
        method (str): the method to pass to `root`
        minsolguess (float): the minimum guess
        maxsolguess (float): the maximum guess
        solpoints (int): the number of guesses (will end up with 3x this number)

    Returns:
        List[float]: solutions to `to_solve`
    """

    guesses = get_list_of_guesses_for_solver(minsolguess, maxsolguess,
                                             solpoints)
    sols = get_solutions_using_root(to_solve, method, guesses)
    return sols


def get_solutions_using_pychebfun(
    to_solve,
    method,
    Nch,
):
    raise DeprecationWarning("This method does not work well. Caution.")
    with warnings.catch_warnings():
        warnings.filterwarnings('error')
        while len(sols) in [0, 2]:
            approx_roots = pychebfun.chebfun(np.vectorize(to_solve),
                                             N=Nch).roots()
            sols = [root(to_solve, x0, method=method).x for x0 in approx_roots]
            if len(sols) > 3:
                raise ValueError("Too Many Sols")
            Nch += 400
            if Nch > 1000:
                raise ValueError("Not enough sols")
    return sols


def get_solutions_using_root(to_solve: Callable, method: str,
                             guesses: List) -> List:
    """Solve function using root with multiple initial guesses between minsolguess and maxsolguess


    Args:
        to_solve (Callable): function to solve
        method (str): solver method string

    Returns:
        list: solutions found of length up to `len(guesses)` maybe less if solver fails
    """
    sols = []
    for guess in guesses:
        sol = root(to_solve, guess, method=method)
        if sol.success:
            sols.append(sol.x)
    return sols


def get_list_of_guesses_for_solver(minsolguess: float, maxsolguess: float,
                                   solpoints: int) -> list:
    """Get a list of guesses to pass to solver.
    
    This is done by logspacing both directions and linspacing and then removing duplicates.

    Args:
        minsolguess (float): minimum guess value
        maxsolguess (float): maximum guess value
        solpoints (int): number of guess points

    Returns:
        list: a list of length up to 3*`solpoints` (less if they overlap)
    """
    asc_log_spaced = np.logspace(np.log10(minsolguess), np.log10(maxsolguess),
                                 solpoints)

    desc_log_spaced = maxsolguess + minsolguess - np.logspace(
        np.log10(maxsolguess), np.log10(minsolguess), solpoints)[-1::-1]

    linear_spaced = np.linspace(minsolguess, maxsolguess, solpoints)
    guesses = np.hstack((asc_log_spaced, desc_log_spaced, linear_spaced))

    guessorder = sorted(np.unique(guesses))

    return guessorder


def generate_p_coupled_sequence(gr=0.5, gd=0.01, k=3, tmax=100, init=None):
    u = u_factory(k)
    mu = lambda p: 1 - u(1 - p)
    plist = []
    if not init:
        init = gd
    plist.append(init)

    while len(plist) < tmax:
        # plist.append(plist[-1] * (1 - gr * mu(plist[-1])))
        plist.append(plist[-1] + (1 - plist[-1]) * gd - gr * mu(plist[-1]))
    return np.array(plist)


def generate_p_coupled_repair_capacity_sequence(gr0=0.5,
                                                gd=0.01,
                                                beta=25,
                                                k=3,
                                                tmax=100,
                                                init=None):
    u = u_factory(k)
    mu = lambda p: 1 - u(1 - p)
    plist = []
    if not init:
        init = gd
    plist.append(init)
    gr = lambda x: gr0 / (1 + gr0 * x * beta)
    bounded = lambda x: max(min(x, 1), 0)
    while len(plist) < tmax:
        # plist.append(plist[-1] * (1 - gr * mu(plist[-1])))
        plist.append(
            bounded(plist[-1] + (1 - plist[-1]) * gd -
                    gr(plist[-1]) * mu(plist[-1]) * plist[-1]))
    return np.array(plist)


def make_recurrence_plot(gr=0.5, gd=0.01, k=3, coupled=True):
    import matplotlib.pyplot as plt
    g = u_factory(k)
    mu = lambda p: (1 - g(1 - p)) if coupled else 1
    curve = np.vectorize(lambda x: gd + x * (1 - gr * mu(x)) * (1 - gd))
    x = np.linspace(0, 1, 1000)
    plt.figure()
    plt.plot(x, x)
    plt.plot(x, curve(x))
    plt.tight_layout()
    plt.show()
