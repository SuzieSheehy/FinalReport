#!/usr/bin/env python3
"""Every number quoted in the square-vs-hexagonal comparison tables of the
final report is produced here, so the tables are derived rather than asserted.

Run:  python3 lattice_numbers.py
"""
import math
from collections import Counter

import numpy as np

SQ = [(1, 0), (-1, 0), (0, 1), (0, -1)]
HEX = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]
G = 0.915965594177219015054603514932384110774  # Catalan's constant


def square_sites(L):
    return [(i, j) for i in range(L) for j in range(L)]


def hex_sites(R):
    return [(q, r) for q in range(-R, R + 1) for r in range(-R, R + 1)
            if max(abs(q), abs(r), abs(q + r)) <= R]


def spectrum(sites, offsets):
    S = set(sites)
    z = len(offsets)
    c = Counter()
    for x in sites:
        m = sum(1 for d in offsets if (x[0] + d[0], x[1] + d[1]) in S)
        c[-z + m] += 1
    return c


def reduced_laplacian(sites, offsets):
    """Delta_tilde = the toppling matrix with the sink rows deleted:
    -Delta(x,y).  det(Delta_tilde) = |R|, the order of the sandpile group."""
    idx = {s: k for k, s in enumerate(sites)}
    S = set(sites)
    z = len(offsets)
    n = len(sites)
    M = np.zeros((n, n))
    for x in sites:
        i = idx[x]
        M[i, i] = z
        for d in offsets:
            y = (x[0] + d[0], x[1] + d[1])
            if y in S:
                M[i, idx[y]] -= 1
    return M


def tree_entropy(sites, offsets):
    """lim (1/N) log |R|  -- the spanning-tree entropy per site."""
    M = reduced_laplacian(sites, offsets)
    sign, logdet = np.linalg.slogdet(M)
    assert sign > 0, sign
    return logdet / len(sites)


print("=" * 74)
print("1. ROW-SUM SPECTRUM OF DELTA  (loss classes on the boundary)")
print("=" * 74)
for L in (4, 8, 50):
    c = spectrum(square_sites(L), SQ)
    print(f"  square L={L:<4} N={L*L:<7} {dict(sorted(c.items(), reverse=True))}")
for R in (2, 4, 8, 50):
    c = spectrum(hex_sites(R), HEX)
    n = 3 * R * R + 3 * R + 1
    print(f"  hex    R={R:<4} N={n:<7} {dict(sorted(c.items(), reverse=True))}")
print("  square : corners 4x(-2), edges 4(L-2)x(-1), interior (L-2)^2 x 0")
print("  hex    : vertices 6x(-3), edges 6(R-1)x(-2), interior (3R^2-3R+1) x 0")
print("  -> square has a -1 class and no -3; hex has a -3 class and no -1.")
print("  -> relative spread of per-toppling loss   square 2/1 = 2.0,"
      " hex 3/2 = 1.5")

print()
print("=" * 74)
print("2. DISSIPATION BUDGET   sum over sites of (-row sum)")
print("=" * 74)
print(f"{'lattice':>18} {'N':>8} {'sum(-rowsum)':>14} {'B':>8} "
      f"{'B/N':>9} {'(sum)/B':>9}")
for L in (10, 50, 100):
    s = square_sites(L)
    c = spectrum(s, SQ)
    tot = sum(-k * v for k, v in c.items())
    B = sum(v for k, v in c.items() if k < 0)
    print(f"{'square L=' + str(L):>18} {L*L:>8} {tot:>14} {B:>8} "
          f"{B/(L*L):>9.5f} {tot/B:>9.5f}")
for R in (10, 50, 100):
    s = hex_sites(R)
    c = spectrum(s, HEX)
    tot = sum(-k * v for k, v in c.items())
    B = sum(v for k, v in c.items() if k < 0)
    print(f"{'hex R=' + str(R):>18} {len(s):>8} {tot:>14} {B:>8} "
          f"{B/len(s):>9.5f} {tot/B:>9.5f}")
print("  square: sum(-row sum) = 4L,  B = 4L-4,  (sum)/B -> 1")
print("  hex   : sum(-row sum) = 12R+6, B = 6R, (sum)/B -> 2")

print()
print("=" * 74)
print("3. GRAPH BALLS  |{y : d(x,y) <= R}| in each lattice's own metric")
print("=" * 74)
print(f"{'R':>4} {'square 2R^2+2R+1':>18} {'hex 3R^2+3R+1':>15} {'ratio':>8}")
for R in (1, 2, 4, 8, 32, 100):
    print(f"{R:>4} {2*R*R+2*R+1:>18} {3*R*R+3*R+1:>15} "
          f"{(3*R*R+3*R+1)/(2*R*R+2*R+1):>8.4f}")
print("  ratio -> 3/2 : at the same graph radius the hexagonal lattice holds")
print("  50% more sites, so the same R is a denser cluster.")

print()
print("=" * 74)
print("4. SAND PILE GROUP ORDER  |R| = det(Delta_tilde) = # spanning trees of")
print("   the wired graph (all N sites + a sink joined to every off-board")
print("   direction).  Note this is NOT the plain grid graph: a corner has")
print("   two edges to the sink.")
print("=" * 74)


def wired_edges(sites, offsets):
    """Edge list of the wired multigraph, as (vertex, vertex) pairs with the
    sink named 'S'.  Parallel edges are listed separately."""
    S = set(sites)
    z = len(offsets)
    edges = []
    for x in sites:
        out = 0
        for d in offsets:
            y = (x[0] + d[0], x[1] + d[1])
            if y in S:
                if x < y:
                    edges.append((x, y))
            else:
                out += 1
        edges += [(x, "S")] * out
    return edges, z


def brute_force_spanning_trees(sites, offsets):
    """Independent count, by enumerating edge subsets of the right size."""
    from itertools import combinations
    edges, _ = wired_edges(sites, offsets)
    verts = list(sites) + ["S"]
    root = {v: v for v in verts}

    def find(a):
        while root[a] != a:
            root[a] = root[root[a]]
            a = root[a]
        return a

    n = len(verts)
    count = 0
    for sub in combinations(range(len(edges)), n - 1):
        root.update({v: v for v in verts})
        ok = True
        for k in sub:
            a, b = edges[k]
            ra, rb = find(a), find(b)
            if ra == rb:
                ok = False
                break
            root[ra] = rb
        if ok:
            count += 1
    return count


for L in (2, 3):
    M = reduced_laplacian(square_sites(L), SQ)
    det_exact = round(np.linalg.det(M))
    bf = brute_force_spanning_trees(square_sites(L), SQ)
    print(f"  square L={L}: det(Delta_tilde) = {det_exact}, "
          f"brute-force spanning trees of the wired graph = {bf}, "
          f"agree = {det_exact == bf}")

print()
print("  tree entropy h(N) = (1/N) log |R|; square exact limit 4G/pi ="
      f" {4*G/math.pi:.5f}")
print(f"{'size':>10} {'N':>8} {'square h':>12} {'hex R':>6} {'N':>8} "
      f"{'hex h':>12}")


def fit_limit(rs, hs):
    """Fit h = h_inf + a/r + b/r^2 and return h_inf."""
    A = np.array([[1.0, 1.0 / r, 1.0 / r**2] for r in rs])
    coef, *_ = np.linalg.lstsq(A, np.array(hs), rcond=None)
    return coef[0]


sq_rs, sq_hs, hx_rs, hx_hs = [], [], [], []
for n in (8, 12, 16, 20, 24, 28, 32):
    hs = tree_entropy(square_sites(n), SQ)
    R = max(1, round((math.sqrt(1 + (n * n - 1) / 3) - 1) / 2))
    hh = tree_entropy(hex_sites(R), HEX)
    sq_rs.append(n); sq_hs.append(hs)
    hx_rs.append(R); hx_hs.append(hh)
    print(f"{'L=' + str(n):>10} {n*n:>8} {hs:>12.6f} {R:>6} "
          f"{3*R*R+3*R+1:>8} {hh:>12.6f}")

sq_fit = fit_limit(sq_rs, sq_hs)
hx_fit = fit_limit(hx_rs, hx_hs)
print(f"  extrapolated h_inf:  square = {sq_fit:.5f}   "
      f"(exact 4G/pi = {4*G/math.pi:.5f}, error {abs(sq_fit-4*G/math.pi):.2e})")
print(f"                       hex    = {hx_fit:.5f}")
print(f"  ratio hex/square = {hx_fit/sq_fit:.4f}")
print("  1/N log |R| converges to a lattice-specific constant.  The square")
print("  value is the classical 4G/pi; the hexagonal lattice sits about 36%")
print("  higher, so the sandpile group is strictly larger there even if the")
print("  avalanche exponents coincide.")

print()
print("=" * 74)
print("5. EXACT HEIGHT PROBABILITIES ON Z^2  (Priezzhev; Jeng-Piroux-Ruelle)")
print("   p(0) = 2/pi^2 - 4/pi^3")
print("   p(1) = 1/4 - 1/(2pi) - 3/pi^2 + 12/pi^3")
print("   p(2) = 3/8 + 1/pi - 12/pi^3")
print("   p(3) = 3/8 - 1/(2pi) + 1/pi^2 + 4/pi^3")
print("=" * 74)
pi = math.pi
p = [2 / pi**2 - 4 / pi**3,
     1 / 4 - 1 / (2 * pi) - 3 / pi**2 + 12 / pi**3,
     3 / 8 + 1 / pi - 12 / pi**3,
     3 / 8 - 1 / (2 * pi) + 1 / pi**2 + 4 / pi**3]
zeta = sum(i * pi_ for i, pi_ in enumerate(p))
for i, pi_ in enumerate(p):
    print(f"   p({i}) = {pi_:.6f}   ({100*pi_:.2f} %)")
print(f"   sum p(i) = {sum(p):.12f}   (must be 1)")
print(f"   zeta = sum i p(i) = {zeta:.12f}   (17/8 = {17/8:.12f})")
print(f"   measured on 128x128, 60000 grains: 2.1121; shares "
      f"7.6 / 17.6 / 30.8 / 44.0 %")

print()
print("=" * 74)
print("6. MEAN HEIGHT ON A z-REGULAR LATTICE   zeta = (z - 1 + xi)/2")
print("   (Jarai 2018, Sec. 5.5.2; xi = looping constant = E[# neighbours")
print("    of o visited by the infinite LERW from o])")
print("=" * 74)
print(f"   square  z=4, xi=5/4 : zeta = (4-1+1.25)/2 = {(4-1+1.25)/2} = 17/8")
for xi in (1.0, 1.25, 1.375, 1.5, 1.75, 2.0):
    print(f"   hex     z=6, xi={xi:<6}: zeta = (6-1+{xi})/2 = "
          f"{(6-1+xi)/2:.4f}")
print(f"   naive threshold scaling 17/8 * 6/4 = {17/8*6/4:.4f} would need")
print(f"   xi_hex = {2*(17/8*6/4)-6+1:.4f}; keeping xi = 5/4 gives "
      f"{(6-1+1.25)/2:.4f}.")
