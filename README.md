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

for all logical-class representatives $`l`$ and all $`s \in \mathrm{rowspan}(H_X)`$. The map $`s \mapsto`$ (integer lift) is nonlinear over GF(2), so the constraints are imposed by **enumerating all $`2^{r_X}`$ elements of $`\mathrm{rowspan}(H_X)`$** (the generators alone are not enough; adding componentwise products of generators is, see the generator form below). All $`2^k`$ logical classes are enumerated as well.

The solution set is computed exactly by integer diagonalization of the constraint matrix (entries $`\in \lbrace -1,0,+1\rbrace`$), using Smith-normal-form-style row and column operations; there is no brute-force search over $`\theta`$. For a solution $`\theta`$, the logical phase is $`\varphi(l) = \sum_j \theta_j l_j \bmod 2^m`$ (on solutions it does not depend on the choice of representative).

**Decision**: "a genuinely level-$`m`$ diagonal logical gate exists (an odd power of logical $`S`$ for $`m=2`$, of logical $`T`$ for $`m=3`$)" ⟺ there is a solution $`\theta`$ giving an odd $`\varphi`$ on some nontrivial logical class. This is decided by whether the gcd of the values of $`\varphi`$ on the generators of the solution module (together with $`2^m`$) equals 1. If the gcd is $`2^\nu`$, the logical phases achievable at level $`m`$ are only multiples of $`2^\nu`$ (that is, only gates of effective level $`m-\nu`$ can be realized).

Relation to arXiv:2602.14499: the paper characterizes this feasibility as the vanishing of the obstruction map (a Bockstein-type homomorphism) of a lifting problem (Thm IV.1 / Cor IV.2). This prototype implements the same feasibility question as a direct check, solving the logicality condition directly as a system of linear congruences over $`\mathbb{Z}_{2^m}`$. It is a minimal implementation of the paper's claim that the criterion at each level can be decided concretely by Gaussian elimination / Smith normal form.

## Generator form (for larger codes)

Enumerating $`\mathrm{rowspan}(H_X)`$ costs $`2^{r_X}`$ and is only usable for small codes. `generator_form.py` imposes the same condition through products of generators. Let $`h`$ run over a basis of $`C_2 = \mathrm{rowspan}(H_X)`$ and $`k_1, k_2, \dots`$ over a basis of $`C_1 = \ker H_Z`$. Expanding the integer lift $`\mathrm{lift}(x_1\oplus\cdots\oplus x_r) = \sum_{S\neq\emptyset} (-2)^{|S|-1} \prod_{i\in S} x_i`$ shows that the phase is constant on every coset of $`C_2`$ if and only if

$$\theta\cdot(h * k_1 * \cdots * k_{t-1}) \equiv 0 \pmod{2^{m-t+1}},\qquad t = 1,\dots,m,$$

that is, $`H^{(m)}\theta \equiv 0 \pmod{2^m}`$ for the matrix $`H^{(m)}`$ with rows $`h,\ 2\,h*k,\ 4\,h*k*k',\ \dots`$ ($`*`$ is the componentwise product). The top-level products are constrained only mod 2, so they are reduced over GF(2) while they are generated, and they contribute at most $`n`$ rows. The remaining rows are fed into an incremental kernel computation over $`\mathbb{Z}_{2^m}`$. These lower-level rows (for $`m=3`$, the rows $`2\,h*k`$, constrained mod 4) must not be reduced by XOR of 0/1 vectors: the integer lift of a sum brings in products of four generators, which are not among the constraints. They can be reduced only as integer vectors mod $`2^{m-t+1}`$; the implementation does not pre-reduce them and lets the kernel computation, which keeps at most $`n`$ generators, absorb them. GF(2) vectors are int bitmasks, and $`\mathbb{Z}_{2^m}`$ vectors are stored as $`m`$ bit planes. The parity of the logical phase is linear in the class, so the decision needs only the $`k`$ logical-X generators.

The tests check that the two forms give the same solution module and the same phase subgroups: on the three codes and on 300 random CSS codes with $`n \le 10`$, for $`m = 1,\dots,4`$, and on 50 random subcodes of the $`[[15,1,3]]`$ code (most with $`T`$ feasible), for $`m = 1,2,3`$.

## Files

- `bockstein_check.py` — the decision procedure (GF(2) linear algebra, integer diagonalization, homogeneous solver mod $`2^m`$, `analyze`, brute-force verification `verify_theta`)
- `generator_form.py` — the generator form (`analyze_generators`, same output as `analyze`), a sampled verification `verify_theta_sampled` for codes too large for `verify_theta`, and a command line for codes given as files
- `codes.py` — definitions of the three codes and their self-checks ($`H_X H_Z^{\mathsf{T}} = 0`$, $`[[n,k,d]]`$ verified by full enumeration), plus the families `rotated_surface(d)` and `reed_muller(r)` for the scaling tests
- `test_known.py` — tests against known facts and agreement of the two forms (prints PASS/FAIL, exits 1 on failure)

## Usage

```
python codes.py                            # code self-checks only
python test_known.py                       # decision + known-fact tests (prints the result table)
python generator_form.py HX.txt HZ.txt [m] # check your own code (one 0/1 row per line; m = 1, 2, 3 if omitted)
```

## Results (run on 2026-10-20, all tests PASS, total time 3.4 s)

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

Generator form at larger sizes (2026-10-20, pure Python on a desktop PC):

| Code | $`n`$ | $`r_X`$ | Result | Rows of $`H^{(3)}`$ (levels 1, 2, 3; level 2 before reduction) | Time |
|---|---|---|---|---|---|
| Rotated surface code, $`d=13`$ | 169 | 84 | $`S`$, $`T`$ infeasible ($`2\cdot\mathbb{Z}_4`$, $`4\cdot\mathbb{Z}_8`$) | 84, 234, 157 | 0.3 s |
| Reed–Muller $`[[255,1,3]]`$ | 255 | 8 | $`T`$ feasible (example passes the sampled check) | 8, 44, 92 | 0.2 s |
| Random CSS code, $`k=1`$ | 185 | 90 | $`T`$ infeasible ($`4\cdot\mathbb{Z}_8`$) | 90, 4185, 185 | 1.2 s |

The random code is a dense case for the number of products, and its time is from a separate run that is not part of `test_known.py`. There the GF(2) span of the triple products fills $`\mathbb{F}_2^{n}`$ early, and generation stops. A code with an odd-phase solution never fills it. Generating and reducing all 376,740 triple products of the same code without stopping early takes 0.7 s.

Side observation: the level-3 phase subgroup of Steane being $`2\cdot\mathbb{Z}_8`$ means that "level-3 transversal gates can reach at most logical $`S`$ ($`\omega_8^2 = i`$)", consistent with the non-vanishing of the level $`2\to 3`$ lifting obstruction (the paper's $`\beta`$ map).

## Implementation notes (interpretive assumptions)

1. "Transversal diagonal gate" is restricted to tensor products of single-qubit diagonal gates $`U(\theta) = \bigotimes_j \mathrm{diag}(1, \omega^{\theta_j})`$ (the same setting as the paper; multi-qubit diagonal gates such as CZ are out of scope).
2. The criterion for "a level-$`m`$ gate exists" is taken to be "there is a solution whose logical phase $`\varphi`$ is odd on some nontrivial logical class" (if $`\varphi`$ is only ever even, $`\omega^{\varphi}`$ can be realized at a coarser level). For the three $`k=1`$ codes there is a single nontrivial class.
3. The decision does not depend on the choice of logical-class representative $`l`$ (under the constancy condition, $`\varphi(l\oplus s) = \varphi(l)`$; this is checked by the solution-verification asserts).
4. The constraints are imposed over all elements of $`\mathrm{rowspan}(H_X)`$ and all logical classes, not just generators (because of the nonlinearity of the lift, $`a\oplus b = a + b - 2ab`$). Since $`2^{r_X} \le 16`$, full enumeration is exact. The generator form imposes the same condition with $`O(r_X\,(\dim C_1)^{m-1})`$ products instead of $`2^{r_X+k}`$ rows.
5. Whenever the answer is feasible, the returned example $`\theta`$ is always independently verified by brute-force phase enumeration that uses no linear algebra (`verify_theta`). For codes too large for that, `verify_theta_sampled` checks the phase on random pairs $`y \in \ker H_Z`$, $`s \in \mathrm{rowspan}(H_X)`$ instead.
