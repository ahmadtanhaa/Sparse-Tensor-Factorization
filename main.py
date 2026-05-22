import math
import itertools
from collections import defaultdict, deque
import matplotlib.pyplot as plt
import numpy as np

# ============================================================
# IEEE-style plotting setup
# ============================================================

plt.rcParams.update({
    "figure.figsize": (3.5, 2.6),     # single-column friendly
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "lines.linewidth": 1.0,
    "axes.linewidth": 0.6,
    "grid.linewidth": 0.4,
    "savefig.dpi": 600,
    "ps.fonttype": 42,
    "pdf.fonttype": 42,
    "text.usetex": False,             # keep False unless your local LaTeX is configured
})

# ============================================================
# Basic helpers
# ============================================================

def ceil_div(a: int, b: int) -> int:
    return (a + b - 1) // b


def nCk(n: int, k: int) -> int:
    return math.comb(n, k) if 0 <= k <= n else 0


# ============================================================
# Windows on active values {2,...,P}
# ============================================================

def build_windows(P: int, Lambda: int):
    active_vals = list(range(2, P + 1))
    windows = [tuple(active_vals[i:i + Lambda]) for i in range(0, len(active_vals), Lambda)]

    value_to_window = {}
    for w_idx, vals in enumerate(windows):
        for v in vals:
            value_to_window[v] = w_idx

    return windows, value_to_window


# ============================================================
# Admissible tuples
# ============================================================

def generate_admissible_tuples(L: int, P: int, Gamma: int):
    tuples_list = []
    supports_list = []

    for d in range(Gamma + 1):
        for supp in itertools.combinations(range(L), d):
            for vals in itertools.product(range(2, P + 1), repeat=d):
                p = [1] * L
                for i, coord in enumerate(supp):
                    p[coord] = vals[i]
                tuples_list.append(tuple(p))
                supports_list.append(tuple(supp))

    return tuples_list, supports_list


# ============================================================
# Tile construction
# ============================================================

def enumerate_gamma_tiles(L: int, P: int, Gamma: int, Lambda: int):
    windows, value_to_window = build_windows(P, Lambda)

    tiles = []
    tile_info = {}
    tid = 0

    for Q in itertools.combinations(range(L), Gamma):
        for win_tuple in itertools.product(range(len(windows)), repeat=Gamma):
            closure_capacity = 1
            for w_idx in win_tuple:
                closure_capacity *= (1 + len(windows[w_idx]))

            tiles.append(tid)
            tile_info[tid] = {
                "Q": tuple(Q),
                "win_tuple": tuple(win_tuple),
                "closure_capacity": closure_capacity,
            }
            tid += 1

    return tiles, tile_info, windows, value_to_window


def feasible_tiles_for_tuple(p, supp, Gamma: int, tile_info, value_to_window):
    supp_set = set(supp)
    if len(supp) > Gamma:
        return []

    fixed_window = {ell: value_to_window[p[ell]] for ell in supp}
    feasible = []

    for tid, info in tile_info.items():
        Q = info["Q"]
        if not supp_set.issubset(Q):
            continue

        ok = True
        for ell in supp:
            pos = Q.index(ell)
            if info["win_tuple"][pos] != fixed_window[ell]:
                ok = False
                break

        if ok:
            feasible.append(tid)

    return feasible


def build_graph(L, P, Gamma, Lambda):
    tuples_list, supports_list = generate_admissible_tuples(L, P, Gamma)
    tiles, tile_info, windows, value_to_window = enumerate_gamma_tiles(L, P, Gamma, Lambda)

    adjacency = defaultdict(list)
    reverse_adj = defaultdict(list)

    for u, (p, supp) in enumerate(zip(tuples_list, supports_list)):
        feas = feasible_tiles_for_tuple(p, supp, Gamma, tile_info, value_to_window)
        adjacency[u] = feas
        for t in feas:
            reverse_adj[t].append(u)

    return tuples_list, tiles, tile_info, adjacency, reverse_adj


# ============================================================
# Max flow (Dinic)
# ============================================================

class Edge:
    def __init__(self, to, rev, cap):
        self.to = to
        self.rev = rev
        self.cap = cap


class Dinic:
    def __init__(self, n):
        self.g = [[] for _ in range(n)]

    def add_edge(self, u, v, cap):
        fwd = Edge(v, len(self.g[v]), cap)
        bwd = Edge(u, len(self.g[u]), 0)
        self.g[u].append(fwd)
        self.g[v].append(bwd)

    def max_flow(self, s, t):
        flow = 0

        while True:
            level = [-1] * len(self.g)
            q = deque([s])
            level[s] = 0

            while q:
                u = q.popleft()
                for e in self.g[u]:
                    if e.cap > 0 and level[e.to] < 0:
                        level[e.to] = level[u] + 1
                        q.append(e.to)

            if level[t] < 0:
                return flow

            it = [0] * len(self.g)

            def dfs(u, f):
                if u == t:
                    return f
                while it[u] < len(self.g[u]):
                    e = self.g[u][it[u]]
                    if e.cap > 0 and level[e.to] == level[u] + 1:
                        pushed = dfs(e.to, min(f, e.cap))
                        if pushed > 0:
                            e.cap -= pushed
                            self.g[e.to][e.rev].cap += pushed
                            return pushed
                    it[u] += 1
                return 0

            while True:
                pushed = dfs(s, 10**18)
                if pushed == 0:
                    break
                flow += pushed


# ============================================================
# Assignment algorithm
# ============================================================

def achievable_per_block(L, P, Gamma, Lambda, Delta, verbose=False):
    tuples_list, tiles, tile_info, adjacency, reverse_adj = build_graph(L, P, Gamma, Lambda)

    remaining = set(range(len(tuples_list)))
    capacity = {t: tile_info[t]["closure_capacity"] for t in tiles}
    loads = {t: 0 for t in tiles}

    while remaining:
        rem = list(remaining)
        tlist = [t for t in tiles if capacity[t] > 0]

        s = 0
        left0 = 1
        right0 = left0 + len(rem)
        sink = right0 + len(tlist)

        mf = Dinic(sink + 1)

        pos_u = {u: i for i, u in enumerate(rem)}
        pos_t = {t: i for i, t in enumerate(tlist)}

        for u in rem:
            mf.add_edge(s, left0 + pos_u[u], 1)

        for u in rem:
            for t in adjacency[u]:
                if capacity[t] > 0:
                    mf.add_edge(left0 + pos_u[u], right0 + pos_t[t], 1)

        for t in tlist:
            mf.add_edge(right0 + pos_t[t], sink, capacity[t])

        flow = mf.max_flow(s, sink)

        if flow < len(rem):
            u = next(iter(remaining))
            feasible_nonfull = [t for t in adjacency[u] if capacity[t] > 0]
            if not feasible_nonfull:
                raise RuntimeError(f"No feasible tile with capacity remains for tuple {u}.")
            t = feasible_nonfull[0]
            remaining.remove(u)
            capacity[t] -= 1
            loads[t] += 1
        else:
            for u in rem:
                u_node = left0 + pos_u[u]
                matched_tile = None

                for e in mf.g[u_node]:
                    if right0 <= e.to < sink and e.cap == 0:
                        t_idx = e.to - right0
                        matched_tile = tlist[t_idx]
                        break

                if matched_tile is None:
                    raise RuntimeError(f"Flow succeeded but no assigned tile found for tuple {u}.")

                capacity[matched_tile] -= 1
                loads[matched_tile] += 1

            remaining.clear()

    total_servers_block = sum(min(Delta, x) for x in loads.values() if x > 0)
    return total_servers_block, loads


def achievable_total(K, L, P, Gamma, Lambda, Delta, verbose=False):
    per_block, loads = achievable_per_block(L, P, Gamma, Lambda, Delta, verbose=verbose)
    total = ceil_div(K, Delta) * per_block
    return total, per_block, loads


def default_tensor(K, L, P, Gamma, Lambda, Delta):
    return ceil_div(K, Delta) * nCk(L, Gamma) * (ceil_div(P, Lambda) ** Gamma) * min(Delta, Lambda ** Gamma)


def tdc_baseline(K, L, P, Gamma, Delta):
    L_prime = P ** L
    return ceil_div(K, Delta) * (L_prime / Gamma) * min(Delta, Gamma)


# ============================================================
# Histogram plot with FULL font control
# ============================================================

def run_histogram_plot(K, L, P, Lambda, Delta, verbose=False, save_eps=True):
    import numpy as np
    import matplotlib.pyplot as plt

    # ---- FULL FONT CONTROL ----
    plt.rcParams.update({
        'font.size': 5,            # base size
        'axes.labelsize': 5,       # axis labels (Γ, N)
        'axes.titlesize': 5,
        'xtick.labelsize': 4,      # x tick labels
        'ytick.labelsize': 4,      # y tick labels
        'legend.fontsize': 4,      # legend text
    })

    gammas = list(range(1, L + 1))
    ach = []
    deff = []
    tdc = []

    print(f"\nSetting: K={K}, L={L}, P={P}, Lambda={Lambda}, Delta={Delta}")

    for g in gammas:
        a, per_block, loads = achievable_total(K, L, P, g, Lambda, Delta, verbose=verbose)
        d = default_tensor(K, L, P, g, Lambda, Delta)
        t = tdc_baseline(K, L, P, g, Delta)

        ach.append(a)
        deff.append(d)
        tdc.append(t)

        filled_tiles = sum(1 for x in loads.values() if x > 0)
        print(
            f"Gamma={g} | Achievable={a} | Default={d} | TDC={t:.2f} "
            f"| per-block={per_block} | filled tiles={filled_tiles}"
        )

    x = np.arange(len(gammas))
    width = 0.16

    # ---- Figure ----
    fig, ax = plt.subplots(figsize=(3.5, 2.3), dpi=300)

    # ---- Bars ----
    ax.bar(x - width, ach, width=width, label='Algorithm')
    ax.bar(x,         deff, width=width, label='Default')
    ax.bar(x + width, tdc,  width=width, label='TDC')

    # ---- Labels ----
    ax.set_xlabel(r'$\Gamma$', fontsize=5, labelpad=1)
    ax.set_ylabel(r'$N$', fontsize=5, labelpad=1)

    # ---- Ticks ----
    ax.set_xticks(x)
    ax.set_xticklabels([str(g) for g in gammas], fontsize=4)
    ax.tick_params(axis='y', labelsize=4)
    ax.tick_params(axis='both', width=0.4, length=2.5)

    # ---- Grid ----
    ax.grid(axis='y', linestyle=':', linewidth=0.35, alpha=0.7)

    # ---- Spines ----
    for spine in ax.spines.values():
        spine.set_linewidth(0.4)

    # ---- Legend ----
    ax.legend(
        loc='upper right',
        fontsize=5,
        frameon=True,
        borderpad=0.2,
        handlelength=0.8,
        handletextpad=0.3,
        labelspacing=0.2,
        borderaxespad=0.2
    )

    # ---- Layout ----
    fig.tight_layout(pad=0.15)

    # ---- Save ----
    pdf_name = f"comparison_L{L}_small.pdf"
    fig.savefig(pdf_name, format='pdf', bbox_inches='tight', pad_inches=0.005)
    print(f"Saved: {pdf_name}")

    if save_eps:
        eps_name = f"comparison_L{L}_small.eps"
        fig.savefig(eps_name, format='eps', bbox_inches='tight', pad_inches=0.005)
        print(f"Saved: {eps_name}")

    plt.show()
    plt.close(fig)


# ============================================================
# Run separate IEEE-style histogram plots
# ============================================================

if __name__ == "__main__":
    # Case 1: L = 3
    run_histogram_plot(K=6, L=3, P=6, Lambda=3, Delta=6, verbose=False, save_eps=True)

    # Case 2: L = 5
    run_histogram_plot(K=8, L=5, P=5, Lambda=2, Delta=4, verbose=False, save_eps=True)

    # Optional:
    # run_histogram_plot(K=8, L=6, P=5, Lambda=2, Delta=4, verbose=False, save_eps=True)
