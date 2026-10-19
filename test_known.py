"""Known-fact tests for the transversal diagonal-gate feasibility checker.

Run:  python test_known.py
Prints a result table and PASS/FAIL lines (ASCII only), exits nonzero on FAIL.

Expected (known) facts:
  Steane  [[7,1,3]] : transversal S  YES   (arXiv:2602.14499, Cor V.1)
  Steane  [[7,1,3]] : transversal T  NO    (arXiv:2602.14499, Cor V.1; holds
                                            even for non-uniform angles)
  RM      [[15,1,3]]: transversal T  YES   (classic result, T^x15 = logical T^dag)
  Surface [[9,1,3]] : transversal T  NO    (known)
  Surface [[9,1,3]] : transversal S  NO    (verified here: achievable
                                            phases are 2*Z_4, logical Z only)
"""

import random
import sys
import time

from bockstein_check import (analyze, verify_theta, constraint_rows,
                             gf2_kernel_basis, gf2_rank)
from codes import (ALL_CODES, check_code, steane, reed_muller_15,
                   rotated_surface, reed_muller)
from generator_form import (analyze_generators, verify_theta_sampled,
                            code_bases, constraint_levels, gf2_reduce_mask,
                            to_mask)

GATE_NAME = {1: "Z", 2: "S", 3: "T"}

failures = []


def same_solutions(HX, HZ, m):
    """Generator form and enumeration have the same solution set (each
    module's generators satisfy the other's constraints) and the same phase
    subgroup on every logical class."""
    n = len(HX[0])
    q = 2 ** m
    a = analyze(HX, HZ, m)
    b = analyze_generators(HX, HZ, m)
    rows, _, _ = constraint_rows(HX, HZ, n, m)
    for g in b["gens"]:
        if any(sum(c * x for c, x in zip(row, g)) % q for row in rows):
            return False
    stab, c2, _, c1 = code_bases(HX, HZ)
    levels = constraint_levels(c2, c1, m, n)
    for g in a["gens"]:
        for t, lev in enumerate(levels, start=1):
            for w in lev:
                if (2 ** (t - 1) * sum(g[j] for j in range(n) if w >> j & 1)) % q:
                    return False
    da = {gf2_reduce_mask(stab, to_mask(l)): d for l, d in a["classes"]}
    db = {gf2_reduce_mask(stab, to_mask(l)): d for l, d in b["classes"]}
    return da == db and a["feasible"] == b["feasible"]


def random_css_code(rng):
    """Random CSS code with n <= 10 and k >= 1 (H_X spans a random subspace
    of ker H_Z), or None if the draw is degenerate."""
    n = rng.randint(3, 10)
    HZ = [[rng.randint(0, 1) for _ in range(n)] for _ in range(rng.randint(1, n - 1))]
    C1 = gf2_kernel_basis(HZ, n)
    HX = []
    for _ in range(rng.randint(1, max(len(C1) - 1, 1))):
        v = [0] * n
        for b in C1:
            if rng.random() < 0.5:
                v = [x ^ y for x, y in zip(v, b)]
        if any(v):
            HX.append(v)
    if not HX or n - gf2_rank(HX) - gf2_rank(HZ) < 1:
        return None
    return HX, HZ


def random_rm15_subcode(rng):
    """Random subcode of the [[15,1,3]] RM code: H_X spans a random subspace
    of its X stabilizers, H_Z a random subset of its Z checks.  Returns None
    if the draw has k > 5 (too many classes for the enumeration)."""
    HX0, HZ0 = reed_muller_15()
    HX = []
    for _ in range(rng.randint(1, 4)):
        v = [0] * 15
        for b in HX0:
            if rng.random() < 0.5:
                v = [x ^ y for x, y in zip(v, b)]
        if any(v):
            HX.append(v)
    HZ = [row for row in HZ0 if rng.random() < 0.8]
    if not HX or not HZ or 15 - gf2_rank(HX) - gf2_rank(HZ) > 5:
        return None
    return HX, HZ


def report(label, ok, detail=""):
    status = "PASS" if ok else "FAIL"
    print("%s  %s%s" % (status, label, ("  [" + detail + "]") if detail else ""))
    if not ok:
        failures.append(label)


def main():
    t0 = time.perf_counter()

    # ------------------------------------------------------------------
    # Code self-checks
    # ------------------------------------------------------------------
    print("== Code self-checks ==")
    codes = []
    for name, ctor, (n, k, d) in ALL_CODES:
        HX, HZ = ctor()
        try:
            check_code(name, HX, HZ, n, k, d)
            report("%s is a valid [[%d,%d,%d]] CSS code" % (name, n, k, d), True)
        except AssertionError as e:
            report(name + " self-check", False, str(e))
        codes.append((name, HX, HZ))

    # ------------------------------------------------------------------
    # Feasibility analysis: each code at levels m = 1, 2, 3
    # ------------------------------------------------------------------
    print()
    print("== Feasibility (genuinely level-m transversal diagonal logical gate) ==")
    results = {}
    for name, HX, HZ in codes:
        for m in (1, 2, 3):
            res = analyze(HX, HZ, m)
            results[(name, m)] = res
            # cross-check any feasible example by brute-force phase enumeration
            if res["feasible"]:
                phases = verify_theta(HX, HZ, m, res["example"])
                ph = phases[tuple(res["example_class"])]
                assert ph % 2 == 1, "example does not realize an odd phase"

    header = "%-26s %-6s %-9s %-16s %s" % ("code", "gate", "feasible",
                                           "phase subgroup", "example theta")
    print(header)
    print("-" * len(header))
    for name, HX, HZ in codes:
        for m in (1, 2, 3):
            res = results[(name, m)]
            q = 2 ** m
            d = min((dd for _, dd in res["classes"]), default=q)
            sub = "{0}" if d == q else "%d*Z_%d" % (d, q)
            ex = str(res["example"]) if res["example"] else "-"
            print("%-26s %-6s %-9s %-16s %s" % (
                name, GATE_NAME[m] + " (m=%d)" % m,
                "YES" if res["feasible"] else "NO", sub, ex))

    # ------------------------------------------------------------------
    # Known-fact assertions
    # ------------------------------------------------------------------
    print()
    print("== Known-fact tests ==")
    st, rm, sf = codes[0][0], codes[1][0], codes[2][0]

    for name in (st, rm, sf):
        report("%s: transversal logical Z (m=1) feasible (sanity)" % name,
               results[(name, 1)]["feasible"])

    report("Steane: transversal S feasible (arXiv:2602.14499 Cor V.1)",
           results[(st, 2)]["feasible"])
    report("Steane: transversal T infeasible even non-uniform (Cor V.1)",
           not results[(st, 3)]["feasible"])
    report("RM15: transversal T feasible (known)",
           results[(rm, 3)]["feasible"])
    report("Surface: transversal T infeasible (known)",
           not results[(sf, 3)]["feasible"])
    report("Surface: transversal S infeasible (verified here)",
           not results[(sf, 2)]["feasible"])

    # ------------------------------------------------------------------
    # Direct reproduction of the canonical gate implementations
    # ------------------------------------------------------------------
    print()
    print("== Canonical implementations (independent brute-force check) ==")

    # Steane: U(3 * 1_7) at m=2, i.e. (S^dag)^x7, is a logical S (phase i on |1bar>)
    HX, HZ = steane()
    try:
        phases = verify_theta(HX, HZ, 2, [3] * 7)
        ph = [p for c, p in phases.items() if any(c)][0]
        report("Steane: theta = 3*1_7 (S^dag each) gives odd logical phase %d mod 4" % ph,
               ph % 2 == 1)
    except AssertionError as e:
        report("Steane: theta = 3*1_7 preserves code space", False, str(e))

    # RM15: U(1_15) at m=3, i.e. T^x15, is a logical T^dag (odd phase on |1bar>)
    HX, HZ = reed_muller_15()
    try:
        phases = verify_theta(HX, HZ, 3, [1] * 15)
        ph = [p for c, p in phases.items() if any(c)][0]
        report("RM15: theta = 1_15 (T each) gives odd logical phase %d mod 8" % ph,
               ph % 2 == 1)
    except AssertionError as e:
        report("RM15: theta = 1_15 preserves code space", False, str(e))

    # ------------------------------------------------------------------
    # Generator form vs. enumeration
    # ------------------------------------------------------------------
    print()
    print("== Generator form (generator_form.py) vs. enumeration ==")
    for name, HX, HZ in codes:
        ok = all(same_solutions(HX, HZ, m) for m in (1, 2, 3, 4))
        report("%s: same solution module and phase subgroups, m=1..4" % name, ok)

    rng = random.Random(2026)
    n_codes = 0
    n_ok = 0
    while n_codes < 300:
        code = random_css_code(rng)
        if code is None:
            continue
        n_codes += 1
        n_ok += all(same_solutions(code[0], code[1], m) for m in (1, 2, 3, 4))
    report("random CSS codes, n <= 10, m=1..4: %d/%d agree" % (n_ok, n_codes),
           n_ok == n_codes)

    # subcodes of RM15 (T feasible on some classes, k up to 5)
    n_codes = 0
    n_ok = 0
    n_t = 0
    while n_codes < 50:
        code = random_rm15_subcode(rng)
        if code is None:
            continue
        n_codes += 1
        n_ok += all(same_solutions(code[0], code[1], m) for m in (1, 2, 3))
        n_t += analyze(code[0], code[1], 3)["feasible"]
    report("random RM15 subcodes, m=1..3: %d/%d agree (%d with T feasible)"
           % (n_ok, n_codes, n_t), n_ok == n_codes and n_t > 0)

    # ------------------------------------------------------------------
    # Generator form at sizes out of reach of the enumeration
    # ------------------------------------------------------------------
    print()
    print("== Generator form at larger sizes ==")
    HX, HZ = rotated_surface(5)
    try:
        check_code("Rotated surface d=5", HX, HZ, 25, 1, 5)
        report("rotated_surface(5) is a valid [[25,1,5]] CSS code", True)
    except AssertionError as e:
        report("rotated_surface(5) self-check", False, str(e))

    for d in (9, 13):
        HX, HZ = rotated_surface(d)
        t1 = time.perf_counter()
        r2 = analyze_generators(HX, HZ, 2)
        r3 = analyze_generators(HX, HZ, 3)
        dt = time.perf_counter() - t1
        report("Surface d=%d (n=%d, r_X=%d): S, T infeasible, phases 2*Z_4, 4*Z_8"
               % (d, d * d, (d * d - 1) // 2),
               not r2["feasible"] and not r3["feasible"]
               and [dd for _, dd in r2["classes"]] == [2]
               and [dd for _, dd in r3["classes"]] == [4],
               "%.2f s, level-3 rows %s" % (dt, r3["rows"]))

    HX, HZ = reed_muller(8)
    t1 = time.perf_counter()
    res = analyze_generators(HX, HZ, 3)
    dt = time.perf_counter() - t1
    ok = res["feasible"]
    if ok:
        try:
            verify_theta_sampled(HX, HZ, 3, res["example"])
        except AssertionError:
            ok = False
    report("RM [[255,1,3]]: transversal T feasible (example passes sampled check)",
           ok, "%.2f s" % dt)

    dt = time.perf_counter() - t0
    print()
    print("Total time: %.3f s" % dt)
    if failures:
        print("RESULT: FAIL (%d failed)" % len(failures))
        return 1
    print("RESULT: ALL PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
