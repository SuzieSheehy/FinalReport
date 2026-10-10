#!/usr/bin/env python3
"""Figures for the final report, in both language versions.

  lattice.pdf / lattice-cn.pdf   (a) square region, (b) hexagonal region, both
                                 coloured by the row sum of the toppling matrix
  scaling.pdf / scaling-cn.pdf   (a) sites within graph distance R, (b) the two
                                 graph balls at the same radius, (c) the share of
                                 the loss carried by the deepest boundary class

Run:  python3 make_figures.py            # both languages
      python3 make_figures.py en         # English only
      python3 make_figures.py cn         # Chinese only

The hexagonal tessellation is checked rather than eyeballed: `check_tiling`
counts how many cells meet at each vertex, and a correct tiling must have
three cells at every interior vertex.
"""
import math
import sys
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SQ = [(1, 0), (-1, 0), (0, 1), (0, -1)]
HEX = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, -1), (-1, 1)]

# row sum -> colour: 0 (conserving) ... 3 (deepest loss)
LOSS_COLORS = ["#f2f2ef", "#9dc3e6", "#3a6ea5", "#9e2a2b"]
LOSS_INDEX = {0: 0, -1: 1, -2: 2, -3: 3}

TEXT = {
 "en": {
  "sq_title": r"(a) square region, $z=4$, $k=4$",
  "sq_legend": [r"interior $0$", r"edge $-1$", r"corner $-2$"],
  "hx_title": r"(b) hexagonal region, $z=6$, $k=6$",
  "hx_legend": [r"interior $0$", r"side $-2$", r"vertex $-3$"],
  "a_title": r"(a) balls of radius $R$",
  "a_x": r"graph radius $R$", "a_y": r"sites within distance $R$",
  "a_l1": r"square, $2R^2+2R+1$", "a_l2": r"hexagonal, $3R^2+3R+1$",
  "a_ratio": r"ratio $\to 3/2$",
  "b_title": r"(b) the same $R$ is a diamond or a hexagon",
  "b_hex": r"hexagon: hexagonal lattice, $N=127$",
  "b_dia": r"diamond: square $L^1$ ball, $N=85$",
  "b_circ": r"circle of radius $R$",
  "c_title": r"(c) the deepest class keeps a $1/\sqrt{N}$ share",
  "c_y": r"deepest class share of loss  [%]",
  "c_l1": r"square $2/\sqrt{N}$", "c_l2": r"hexagon $\approx 2.6/\sqrt{N}$",
 },
 "cn": {
  "sq_title": r"(a) 方形区域，$z=4$，$k=4$",
  "sq_legend": [r"内点 $0$", r"边点 $-1$", r"角点 $-2$"],
  "hx_title": r"(b) 六边形区域，$z=6$，$k=6$",
  "hx_legend": [r"内点 $0$", r"边 $-2$", r"顶点 $-3$"],
  "a_title": r"(a) 半径 $R$ 的球",
  "a_x": r"图距离 $R$", "a_y": r"距离 $R$ 内的格点数",
  "a_l1": r"方形，$2R^2+2R+1$", "a_l2": r"六边形，$3R^2+3R+1$",
  "a_ratio": r"比值 $\to 3/2$",
  "b_title": r"(b) 同样的 $R$：菱形与六边形",
  "b_hex": r"六边形：六边形晶格，$N=127$",
  "b_dia": r"菱形：方形 $L^1$ 球，$N=85$",
  "b_circ": r"半径 $R$ 的圆",
  "c_title": r"(c) 最深一类保持 $1/\sqrt{N}$ 占比",
  "c_y": r"最深一类损失占比  [%]",
  "c_l1": r"方形 $2/\sqrt{N}$", "c_l2": r"六边形 $\approx 2.6/\sqrt{N}$",
 },
}


def use_style(lang):
    """Serif text matching the document; usetex for EN, CJK font for CN."""
    if lang == "en":
        try:
            matplotlib.rcParams["text.usetex"] = True
            matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
        except Exception:
            pass
        matplotlib.rcParams["font.family"] = "serif"
    else:
        matplotlib.rcParams["text.usetex"] = False
        matplotlib.rcParams["font.family"] = ["Songti SC", "STSong",
                                              "PingFang SC", "Arial Unicode MS"]
        matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["mathtext.fontset"] = "cm"
    matplotlib.rcParams["font.size"] = 9
    matplotlib.rcParams["axes.linewidth"] = 0.6


def square_sites(L):
    return [(i, j) for i in range(L) for j in range(L)]


def hex_sites(R):
    return [(q, r) for q in range(-R, R + 1) for r in range(-R, R + 1)
            if max(abs(q), abs(r), abs(q + r)) <= R]


def row_sums(sites, offsets):
    S = set(sites)
    z = len(offsets)
    return {x: -z + sum(1 for d in offsets
                        if (x[0] + d[0], x[1] + d[1]) in S) for x in sites}


HEX_R = 1 / math.sqrt(3)          # circumradius of the Voronoi cell


def hex_vertices(cx, cy, radius=HEX_R):
    """Vertices of the Voronoi cell of a site at (cx, cy).

    The neighbour directions are 0, 60, ..., 300 degrees, so the cell's EDGES
    face those directions and its VERTICES lie at 30 + 60k degrees.  Getting
    this 30 degrees wrong makes the tiles touch corner-to-corner instead of
    edge-to-edge, which is not a hexagonal tiling at all.
    """
    return [(cx + radius * math.cos(math.pi / 6 + k * math.pi / 3),
             cy + radius * math.sin(math.pi / 6 + k * math.pi / 3))
            for k in range(6)]


def check_tiling(sites):
    """How many cells meet at each vertex.  A correct hexagonal tiling puts
    exactly three cells at every interior vertex."""
    V = Counter()
    for x in sites:
        cx, cy = x[0] + x[1] / 2.0, x[1] * math.sqrt(3) / 2.0
        for vx, vy in hex_vertices(cx, cy):
            V[(round(vx, 9), round(vy, 9))] += 1
    return dict(sorted(Counter(V.values()).items()))


def draw_square(ax, L, T):
    sites = square_sites(L)
    rs = row_sums(sites, SQ)
    for x in sites:
        ax.add_patch(plt.Rectangle((x[0] - 0.5, x[1] - 0.5), 1, 1,
                                   facecolor=LOSS_COLORS[LOSS_INDEX[rs[x]]],
                                   edgecolor="white", linewidth=0.7))
    o = ((L - 1) / 2, (L - 1) / 2)
    for d in SQ:
        ax.annotate("", xy=(o[0] + d[0] * 0.95, o[1] + d[1] * 0.95),
                    xytext=o, arrowprops=dict(arrowstyle="-|>", lw=1.1,
                                              color="#222222",
                                              shrinkA=0, shrinkB=0))
    ax.set_xlim(-1.3, L + 0.3)
    ax.set_ylim(-1.3, L + 0.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(T["sq_title"], fontsize=9)
    for i, t in enumerate(T["sq_legend"]):
        ax.text(-1.2, -0.75 - 0.55 * i, t, fontsize=7.2, va="center")
        ax.add_patch(plt.Rectangle((-2.6, -0.85 - 0.55 * i), 0.28, 0.28,
                                   facecolor=LOSS_COLORS[i],
                                   edgecolor="white", linewidth=0.6,
                                   clip_on=False))


def draw_hex(ax, R, T):
    sites = hex_sites(R)
    rs = row_sums(sites, HEX)
    for x in sites:
        cx = x[0] + x[1] / 2.0
        cy = x[1] * math.sqrt(3) / 2.0
        ax.add_patch(plt.Polygon(hex_vertices(cx, cy), closed=True,
                                 facecolor=LOSS_COLORS[LOSS_INDEX[rs[x]]],
                                 edgecolor="white", linewidth=0.7))
    for d in HEX:
        ax.annotate("", xy=(d[0] + d[1] / 2.0, d[1] * math.sqrt(3) / 2.0),
                    xytext=(0, 0),
                    arrowprops=dict(arrowstyle="-|>", lw=1.1, color="#222222",
                                    shrinkA=0, shrinkB=0))
    ax.set_xlim(-R - 1.4, R + 1.4)
    ax.set_ylim(-R * math.sqrt(3) / 2 - 1.4, R * math.sqrt(3) / 2 + 1.4)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(T["hx_title"], fontsize=9)
    x0 = -R - 1.35
    y0 = -R * math.sqrt(3) / 2 - 1.05
    for i, t in enumerate(T["hx_legend"]):
        ax.text(x0 + 0.45, y0 + 0.82 - 0.55 * i, t, fontsize=7.2, va="center")
        ax.add_patch(plt.Polygon(hex_vertices(x0, y0 + 0.82 - 0.55 * i,
                                              radius=0.30),
                                 closed=True, facecolor=LOSS_COLORS[i],
                                 edgecolor="white", linewidth=0.6,
                                 clip_on=False))


def make_lattice(lang, T):
    fig = plt.figure(figsize=(6.4, 2.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.0], wspace=0.16)
    L, R = 8, 4                       # N = 64 and N = 61: matched sizes
    draw_square(fig.add_subplot(gs[0, 0]), L, T)
    draw_hex(fig.add_subplot(gs[0, 1]), R, T)
    suffix = "" if lang == "en" else "-cn"
    fig.savefig(f"lattice{suffix}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote lattice{suffix}.pdf  (square N=64, hexagon N=61)")


def make_scaling(lang, T):
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(7.0, 2.7))

    Rs = list(range(1, 13))
    sq_n = [2 * r * r + 2 * r + 1 for r in Rs]
    hx_n = [3 * r * r + 3 * r + 1 for r in Rs]
    ax1.plot(Rs, sq_n, "o-", ms=3.4, lw=1.1, color="#3a6ea5", label=T["a_l1"])
    ax1.plot(Rs, hx_n, "s-", ms=3.4, lw=1.1, color="#9e2a2b", label=T["a_l2"])
    ax1.set_xlabel(T["a_x"])
    ax1.set_ylabel(T["a_y"])
    ax1.set_title(T["a_title"], fontsize=8.6)
    ax1.legend(fontsize=6.6, frameon=False, loc="upper left")
    ax1.grid(alpha=0.25, lw=0.5)
    ax1.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax1.spines[s].set_visible(False)
    ax1.annotate(T["a_ratio"], xy=(9.6, 0.5 * (sq_n[8] + hx_n[8])),
                 xytext=(7.4, 108), fontsize=7,
                 arrowprops=dict(arrowstyle="-", lw=0.6, color="#444444"))

    Rb = 6
    ax2.add_patch(plt.Polygon([(Rb, 0), (0, Rb), (-Rb, 0), (0, -Rb)],
                              closed=True, facecolor="#3a6ea5", alpha=0.22,
                              ec="#3a6ea5", lw=1.0))
    ax2.add_patch(plt.Polygon([(Rb * math.cos(k * math.pi / 3),
                                Rb * math.sin(k * math.pi / 3))
                               for k in range(6)], closed=True,
                              facecolor="#9e2a2b", alpha=0.22, ec="#9e2a2b",
                              lw=1.0))
    th = [k * math.pi / 180 for k in range(0, 361, 2)]
    ax2.plot([Rb * math.cos(t) for t in th], [Rb * math.sin(t) for t in th],
             lw=1.1, ls=":", color="#333333")
    ax2.text(0, -Rb - 1.75, T["b_hex"], fontsize=7, color="#7a1f20",
             ha="center", va="top")
    ax2.text(0, -Rb - 2.80, T["b_dia"], fontsize=7, color="#27486b",
             ha="center", va="top")
    ax2.text(0, -Rb - 0.45, T["b_circ"], fontsize=7, color="#333333",
             ha="center", va="top")
    ax2.set_xlim(-Rb - 1.2, Rb + 1.2)
    ax2.set_ylim(-Rb - 3.70, Rb + 1.0)
    ax2.set_aspect("equal")
    ax2.axis("off")
    ax2.set_title(T["b_title"], fontsize=8.6)

    Ns = [10 ** (k / 8) for k in range(24, 33)]
    sq_share = [200.0 / math.sqrt(n) for n in Ns]
    hx_share = [100.0 * 3 / (2 * math.sqrt(n / 3) + 1) for n in Ns]
    ax3.loglog(Ns, sq_share, "-", lw=1.1, color="#3a6ea5", label=T["c_l1"])
    ax3.loglog(Ns, hx_share, "-", lw=1.1, color="#9e2a2b", label=T["c_l2"])
    ax3.set_xlabel(r"$N$")
    ax3.set_ylabel(T["c_y"])
    ax3.set_title(T["c_title"], fontsize=8.6)
    ax3.legend(fontsize=6.6, frameon=False)
    ax3.grid(alpha=0.25, lw=0.5, which="both")
    ax3.tick_params(labelsize=7)
    for s in ("top", "right"):
        ax3.spines[s].set_visible(False)

    fig.tight_layout()
    suffix = "" if lang == "en" else "-cn"
    fig.savefig(f"scaling{suffix}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote scaling{suffix}.pdf")


if __name__ == "__main__":
    want = sys.argv[1] if len(sys.argv) > 1 else "both"
    langs = ["en", "cn"] if want == "both" else [want]

    # the tessellation is verified once, before anything is drawn
    share = check_tiling(hex_sites(4))
    assert share.get(3, 0) >= 90, f"hexagonal tiling broken: {share}"
    print("hexagonal tiling check (cells meeting at a vertex):", share)

    for lg in langs:
        use_style(lg)
        make_lattice(lg, TEXT[lg])
        make_scaling(lg, TEXT[lg])
    print("  (b) at R=6 the square L1 ball has 85 sites, the hexagonal ball 127")
    print("  (c) square 2/sqrt(N), hexagon 3*sqrt(3)/(2*sqrt(N)) = 2.5981/sqrt(N),"
          " ratio 1.2990")
