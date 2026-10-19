"""Generator form of the transversal diagonal-gate feasibility check.

Same decision as bockstein_check.analyze, but the logicality condition is
imposed through products of generators instead of by enumerating all of
rowspan(H_X), so it scales to codes with hundreds of qubits.

Constraints
-----------
Let h run over a basis of C_2 = rowspan(H_X) and k_1, k_2, ... over a basis of
C_1 = ker(H_Z).  The integer lift of an XOR sum is

    lift(x_1 xor ... xor x_r) = sum_{S nonempty} (-2)^{|S|-1} prod_{i in S} x_i,

so the phase theta . lift(y) is constant on every coset y + C_2 inside C_1 iff

    theta . (h * k_1 * ... * k_{t-1}) == 0  (mod 2^{m-t+1}),   t = 1, ..., m,

i.e. the level-m matrix H^(m) has rows h, 2 h*k, 4 h*k*k', ... (mod 2^m).
(Expand theta . lift(y) as a multilinear polynomial in the coordinates of y
in a basis of C_1 that extends a basis of C_2; constancy on the cosets means
that every monomial containing a C_2 coordinate has vanishing coefficient.)

The top-level products (t = m) are constrained only mod 2, so they are
reduced over GF(2) as they are generated and contribute at most n rows.
The remaining rows are fed one by one into an incremental kernel computation
over Z_{2^m}.

Logical phases
--------------
On a solution theta, the logical class of y in C_1 gets phase theta . y
(mod 2^m).  Its parity is linear in y, so a genuinely level-m gate exists iff
some solution has an odd phase on some logical-X generator; only the k
generators are needed for the decision.  Per-class phase subgroups are listed
for all 2^k classes when k <= MAX_K_CLASSES, otherwise for the generators only.

Vectors are Python ints used as bitmasks; vectors over Z_{2^m} are stored as m
bit planes.  Pure Python, no external dependencies.

Run:  python generator_form.py HX.txt HZ.txt [m]
(each file: one row per line, entries 0/1, separators optional)
"""

import random
import sys
from itertools import combinations_with_replacement
from math import gcd

from bockstein_check import gf2_kernel_basis

MAX_K_CLASSES = 10

# ----------------------------------------------------------------------
# GF(2) vectors as int bitmasks
# ----------------------------------------------------------------------


def to_mask(v):
    return sum(1 << j for j, x in enumerate(v) if x)


def from_mask(x, n):
    return [x >> j & 1 for j in range(n)]


def gf2_insert(basis, x):
    """Reduce x against basis {pivot bit: row} and insert it if independent.
    Returns True iff x was independent."""
    while x:
        p = x.bit_length() - 1
        b = basis.get(p)
        if b is None:
            basis[p] = x
            return True
        x ^= b
    return False


def gf2_reduce_mask(basis, x):
    """Fully reduce x against basis: the result has no pivot bits set, so it is
    a canonical representative of x modulo the span."""
    for p in sorted(basis, reverse=True):
        if x >> p & 1:
            x ^= basis[p]
    return x


# ----------------------------------------------------------------------
# Vectors over Z_{2^m} as m bit planes (plane i = bit i of every entry)
# ----------------------------------------------------------------------


def _add(a, b):
    out = []
    carry = 0
    for x, y in zip(a, b):
        out.append(x ^ y ^ carry)
        carry = (x & y) | (carry & (x ^ y))
    return out


def _neg(a, full):
    ones = [full] + [0] * (len(a) - 1)
    return _add([x ^ full for x in a], ones)


def _shift(a, e):
    """Multiply by 2^e (mod 2^m)."""
    m = len(a)
    return [0] * min(e, m) + a[:max(m - e, 0)]


def _scale(a, c):
    """Multiply by the integer c (mod 2^m)."""
    out = [0] * len(a)
    e = 0
    while c and e < len(a):
        if c & 1:
            out = _add(out, _shift(a, e))
        c >>= 1
        e += 1
    return out


def _dot(a, w, q):
    """a . w (mod q) for a over Z_q (planes) and a 0/1 vector w (mask)."""
    return sum((x & w).bit_count() << i for i, x in enumerate(a)) % q


def _to_list(a, n):
    return [sum((x >> j & 1) << i for i, x in enumerate(a)) for j in range(n)]


def _valuation(v):
    return (v & -v).bit_length() - 1


# ----------------------------------------------------------------------
# Constraints and solver
# ----------------------------------------------------------------------


def code_bases(HX, HZ):
    """Bitmask bases: (stabilizer echelon {pivot: row}, C_2 basis, logical-X
    generators, C_1 basis = C_2 basis + logical generators)."""
    n = len(HX[0])
    xs = [to_mask(row) for row in HX]
    for z in map(to_mask, HZ):
        assert all((x & z).bit_count() % 2 == 0 for x in xs), "H_X H_Z^T != 0"
    stab = {}
    for row in HX:
        gf2_insert(stab, to_mask(row))
    c2 = list(stab.values())
    ext = dict(stab)
    logicals = []
    for v in gf2_kernel_basis(HZ, n):
        x = to_mask(v)
        if gf2_insert(ext, x):
            logicals.append(gf2_reduce_mask(stab, x))
    return stab, c2, logicals, c2 + logicals


def constraint_levels(c2, c1, m, n):
    """Constraint rows of H^(m), grouped by level t = 1..m.

    Level t < m: list of distinct masks w = h * k_1 * ... * k_{t-1}, each
    meaning theta . w == 0 (mod 2^{m-t+1}).  Level m: a GF(2) basis of the
    span of all such products (they are constrained only mod 2)."""
    full = (1 << n) - 1
    levels = []
    for t in range(1, m + 1):
        seen = set()
        top = {}
        for h in c2:
            for ks in combinations_with_replacement(range(len(c1)), t - 1):
                w = h
                for i in ks:
                    w &= c1[i]
                if not w or w in seen:
                    continue
                seen.add(w)
                if t == m:
                    gf2_insert(top, w)
                    if len(top) == n:
                        break
            if t == m and len(top) == n:
                break
        levels.append(list(top.values()) if t == m else sorted(seen))
    assert all(w & ~full == 0 for lev in levels for w in lev)
    return levels


def solve_levels(levels, n, m):
    """Generators of {theta in Z_{2^m}^n : all level constraints hold}.

    Incremental kernel: keep generators G of the solution module of the
    constraints so far; for a new constraint with values c(g), take g* of
    minimal 2-adic valuation v*, clear the other values with g -> g - a g*,
    and replace g* by 2^{m-v*} g*."""
    q = 1 << m
    full = (1 << n) - 1
    gens = [[1 << j] + [0] * (m - 1) for j in range(n)]
    # top level first: it is the cheapest and usually shrinks the module most
    order = [m] + list(range(1, m))
    for t in order:
        coef = 1 << (t - 1)
        for w in levels[t - 1]:
            vals = [(coef * _dot(g, w, q)) % q for g in gens]
            nz = [i for i, v in enumerate(vals) if v]
            if not nz:
                continue
            p = min(nz, key=lambda i: _valuation(vals[i]))
            vp = _valuation(vals[p])
            uinv = pow(vals[p] >> vp, -1, q)
            gp = gens[p]
            for i in nz:
                if i != p:
                    a = ((vals[i] >> vp) * uinv) % q
                    gens[i] = _add(gens[i], _neg(_scale(gp, a), full))
            gens[p] = _shift(gp, m - vp)
            gens = [g for g in gens if any(g)]
    return gens


def analyze_generators(HX, HZ, m):
    """Generator-form counterpart of bockstein_check.analyze (same keys, plus
    'rows': number of constraint rows per level after reduction)."""
    n = len(HX[0])
    q = 1 << m
    stab, c2, logicals, c1 = code_bases(HX, HZ)
    levels = constraint_levels(c2, c1, m, n)
    gens = solve_levels(levels, n, m)

    # sanity: every generator satisfies every constraint
    for g in gens:
        for t, lev in enumerate(levels, start=1):
            for w in lev:
                assert ((1 << (t - 1)) * _dot(g, w, q)) % q == 0, "solver bug"

    if len(logicals) <= MAX_K_CLASSES:
        reps = [0]
        for x in logicals:
            reps = reps + [r ^ x for r in reps]
        reps = reps[1:]
    else:
        reps = logicals

    classes = []
    example = None
    example_class = None
    for l in reps:
        vals = [_dot(g, l, q) for g in gens]
        d = q
        for v in vals:
            d = gcd(d, v)
        classes.append((from_mask(l, n), d))
        if d % 2 == 1 and example is None:
            i = next(i for i, v in enumerate(vals) if v % 2 == 1)
            example = [(x * pow(vals[i], -1, q)) % q
                       for x in _to_list(gens[i], n)]
            example_class = from_mask(l, n)
    # parity is linear in the class, so the generators decide feasibility
    feasible = any(_dot(g, l, q) % 2 for g in gens for l in logicals)
    return {
        "gens": [_to_list(g, n) for g in gens],
        "classes": classes,
        "feasible": feasible,
        "example": example,
        "example_class": example_class,
        "n": n,
        "m": m,
        "rows": [len(lev) for lev in levels],
    }


def verify_theta_sampled(HX, HZ, m, theta, samples=2000, seed=0):
    """Randomized check that U(theta) preserves the code space, independent of
    the constraint derivation: for random y in ker(H_Z) and s in rowspan(H_X),
    compare theta . lift(y xor s) with theta . lift(y) mod 2^m.  Returns the
    logical phases of the logical-X generators; raises AssertionError on a
    violation."""
    n = len(HX[0])
    q = 1 << m
    _, c2, logicals, c1 = code_bases(HX, HZ)
    rng = random.Random(seed)

    def phase(x):
        return sum(theta[j] for j in range(n) if x >> j & 1) % q

    def rand_comb(basis):
        x = 0
        for b in basis:
            if rng.random() < 0.5:
                x ^= b
        return x

    for _ in range(samples):
        y = rand_comb(c1)
        s = rand_comb(c2)
        assert phase(y ^ s) == phase(y), "theta does not preserve the code space"
    return [phase(l) for l in logicals]


# ----------------------------------------------------------------------
# Command line
# ----------------------------------------------------------------------


def load_matrix(path):
    rows = []
    with open(path) as f:
        for line in f:
            bits = [int(c) for c in line if c in "01"]
            if bits:
                rows.append(bits)
    return rows


def main(argv):
    if len(argv) not in (3, 4):
        print(__doc__.strip().splitlines()[-2])
        return 2
    HX, HZ = load_matrix(argv[1]), load_matrix(argv[2])
    ms = [int(argv[3])] if len(argv) == 4 else [1, 2, 3]
    for x in HX:
        for z in HZ:
            assert sum(a & b for a, b in zip(x, z)) % 2 == 0, "H_X H_Z^T != 0"
    for m in ms:
        res = analyze_generators(HX, HZ, m)
        q = 1 << m
        d = min((dd for _, dd in res["classes"]), default=q)
        sub = "{0}" if d == q else "%d*Z_%d" % (d, q)
        print("m=%d  feasible=%s  phase subgroup=%s  rows/level=%s" % (
            m, "YES" if res["feasible"] else "NO", sub, res["rows"]))
        if res["example"]:
            print("      example theta = %s" % res["example"])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
