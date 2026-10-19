# bockstein-gate-check

[日本語版 README](README_ja.md)

A minimal prototype that decides whether a CSS code admits a transversal diagonal gate.
It is a minimal reference implementation of the decision principle of arXiv:2602.14499 (J. Haruna, "Homological origin of transversal implementability of logical diagonal gates in quantum CSS codes"), and checks the decision numerically on well-known small codes.

Pure Python, no external dependencies.

## Decision principle

For a CSS code (X-check matrix $`H_X`$, Z-check matrix $`H_Z`$, $`H_X H_Z^{\mathsf{T}} = 0`$), consider a level-$`m`$ transversal diagonal gate

$$U(\theta) = \bigotimes_j \mathrm{diag}\bigl(1,\ \omega^{\theta_j}\bigr),\qquad \omega = e^{2\pi i/2^m},\qquad \theta \in (\mathbb{Z}_{2^m})^n$$

where the angles $`\theta_j`$ may be non-uniform.
Codewords are $`|\bar{x}\rangle = \sum_{s\in\mathrm{rowspan}(H_X)} |l\oplus s\rangle`$, with $`l \in \ker H_Z`$ a logical-class representative.
Since $`U`$ is diagonal, $`U`$ preserves the code space and acts as a diagonal logical gate if and only if, on every logical class, the phase $`\sum_j \theta_j y_j \bmod 2^m`$ is constant (here $`y`$ is the integer lift of a GF(2) vector).
Using the integer lift $`(l\oplus s)_j = l_j + s_j - 2 l_j s_j`$, this constancy reduces to the system of linear congruences over $`\mathbb{Z}_{2^m}`$

$$\sum_j \theta_j\, s_j\, (1 - 2 l_j) \equiv 0 \pmod{2^m}$$

for all logical-class representatives $`l`$ and all $`s \in \mathrm{rowspan}(H_X)`$. The map $`s \mapsto`$ (integer lift) is nonlinear over GF(2), so the constraints are imposed by **enumerating all $`2^{r_X}`$ elements of $`\mathrm{rowspan}(H_X)`$** (the generators alone are not enough). All $`2^k`$ logical classes are enumerated as well.

The solution set is computed exactly by integer diagonalization of the constraint matrix (entries $`\in \lbrace -1,0,+1\rbrace`$), using Smith-normal-form-style row and column operations; there is no brute-force search over $`\theta`$. For a solution $`\theta`$, the logical phase is $`\varphi(l) = \sum_j \theta_j l_j \bmod 2^m`$ (on solutions it does not depend on the choice of representative).

**Decision**: "a genuinely level-$`m`$ diagonal logical gate exists (an odd power of logical $`S`$ for $`m=2`$, of logical $`T`$ for $`m=3`$)" ⟺ there is a solution $`\theta`$ giving an odd $`\varphi`$ on some nontrivial logical class. This is decided by whether the gcd of the values of $`\varphi`$ on the generators of the solution module (together with $`2^m`$) equals 1. If the gcd is $`2^\nu`$, the logical phases achievable at level $`m`$ are only multiples of $`2^\nu`$ (that is, only gates of effective level $`m-\nu`$ can be realized).

Relation to arXiv:2602.14499: the paper characterizes this feasibility as the vanishing of the obstruction map (a Bockstein-type homomorphism) of a lifting problem (Thm IV.1 / Cor IV.2). This prototype implements the same feasibility question as a direct check, solving the logicality condition directly as a system of linear congruences over $`\mathbb{Z}_{2^m}`$. It is a minimal implementation of the paper's claim that the criterion at each level can be decided concretely by Gaussian elimination / Smith normal form.

## Files

- `bockstein_check.py` — the decision procedure (GF(2) linear algebra, integer diagonalization, homogeneous solver mod $`2^m`$, `analyze`, brute-force verification `verify_theta`)
- `codes.py` — definitions of the three codes and their self-checks ($`H_X H_Z^{\mathsf{T}} = 0`$, $`[[n,k,d]]`$ verified by full enumeration)
- `test_known.py` — tests against known facts (prints PASS/FAIL, exits 1 on failure)

## Usage

```
python codes.py        # code self-checks only
python test_known.py   # decision + known-fact tests (prints the result table)
```

## Results (run on 2026-08-11, all tests PASS, total time 0.02 s)

| Code | Gate | Feasible | Achievable logical-phase subgroup | Example $`\theta`$ |
|---|---|---|---|---|
| Steane $`[[7,1,3]]`$ | $`Z`$ ($`m=1`$) | **Yes** | all of $`\mathbb{Z}_2`$ | 1 on the support of logical $`X`$ |
| Steane $`[[7,1,3]]`$ | $`S`$ ($`m=2`$) | **Yes** | all of $`\mathbb{Z}_4`$ | $`[1,1,3,1,3,3,1]`$ |
| Steane $`[[7,1,3]]`$ | $`T`$ ($`m=3`$) | **No** | $`2\cdot\mathbb{Z}_8`$ (even only) | — |
| Reed–Muller $`[[15,1,3]]`$ | $`S`$ ($`m=2`$) | **Yes** | all of $`\mathbb{Z}_4`$ | $`[1,1,3,3,1,1,3,0,\dots,0]`$ |
| Reed–Muller $`[[15,1,3]]`$ | $`T`$ ($`m=3`$) | **Yes** | all of $`\mathbb{Z}_8`$ | $`[1,1,7,1,7,7,1,1,7,7,1,7,1,1,7]`$ |
| Rotated surface code $`[[9,1,3]]`$ | $`S`$ ($`m=2`$) | **No** | $`2\cdot\mathbb{Z}_4`$ (logical $`Z`$ only) | — |
| Rotated surface code $`[[9,1,3]]`$ | $`T`$ ($`m=3`$) | **No** | $`4\cdot\mathbb{Z}_8`$ (logical $`Z`$ only) | — |

Comparison with known facts (all reproduced):

- **Steane: transversal $`S`$ feasible, transversal $`T`$ infeasible** (reproduces arXiv:2602.14499 Cor V.1). The infeasibility of $`T`$ holds **even with non-uniform angles** (all logical phases of level-3 solutions are even, i.e. $`2\cdot\mathbb{Z}_8`$). Moreover, brute-force verification confirms that $`\theta = 3\cdot 1_7`$ ($`S^\dagger`$ on every qubit) is a logical $`S`$ (matching the explicit solution in Cor V.1).
- **RM $`[[15,1,3]]`$: transversal $`T`$ feasible**. Brute-force verification confirms that $`\theta = 1_{15}`$ ($`T`$ on every qubit) gives logical phase $`7 \bmod 8`$ (= logical $`T^\dagger`$), matching the classic result.
- **Rotated surface code: transversal $`T`$ infeasible** (known). We also confirm that **transversal $`S`$ is infeasible** (the achievable phases are $`2\cdot\mathbb{Z}_4`$, i.e. up to logical $`Z`$).

Side observation: the level-3 phase subgroup of Steane being $`2\cdot\mathbb{Z}_8`$ means that "level-3 transversal gates can reach at most logical $`S`$ ($`\omega_8^2 = i`$)", consistent with the non-vanishing of the level $`2\to 3`$ lifting obstruction (the paper's $`\beta`$ map).

## Implementation notes (interpretive assumptions)

1. "Transversal diagonal gate" is restricted to tensor products of single-qubit diagonal gates $`U(\theta) = \bigotimes_j \mathrm{diag}(1, \omega^{\theta_j})`$ (the same setting as the paper; multi-qubit diagonal gates such as CZ are out of scope).
2. The criterion for "a level-$`m`$ gate exists" is taken to be "there is a solution whose logical phase $`\varphi`$ is odd on some nontrivial logical class" (if $`\varphi`$ is only ever even, $`\omega^{\varphi}`$ can be realized at a coarser level). For the three $`k=1`$ codes there is a single nontrivial class.
3. The decision does not depend on the choice of logical-class representative $`l`$ (under the constancy condition, $`\varphi(l\oplus s) = \varphi(l)`$; this is checked by the solution-verification asserts).
4. The constraints are imposed over all elements of $`\mathrm{rowspan}(H_X)`$ and all logical classes, not just generators (because of the nonlinearity of the lift, $`a\oplus b = a + b - 2ab`$). Since $`2^{r_X} \le 16`$, full enumeration is exact.
5. Whenever the answer is feasible, the returned example $`\theta`$ is always independently verified by brute-force phase enumeration that uses no linear algebra (`verify_theta`).
