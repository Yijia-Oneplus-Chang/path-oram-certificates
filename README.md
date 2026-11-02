# Exact Finite Certificates for Path ORAM

This package verifies the finite arithmetic premises used in **The Minimum
Bucket Capacity for Reliable Path ORAM**.
It recomputes both certificates and compares every stored entry with the
supplied table. The paper proves why these computations give valid bounds
and how the finite premises imply the general results.

## Run the Verification

From the directory containing `verify_certificates.py` and `certificates.json`,
run:

```sh
python3 verify_certificates.py
```

Only Python 3's standard library is needed. No LaTeX installation, network
access, or third-party Python package is required. The calculations use exact
integers and `fractions.Fraction`, not floating-point estimates or simulations.
Default verification is read-only. Success returns exit status zero. A failed
arithmetic check or a mismatch with the supplied table returns a nonzero status.

The standalone archive includes `expected-results.txt`, which records the
successful terminal output. `MANUSCRIPT.json` identifies the accompanying
manuscript by its title, PDF SHA256 hash, and the hashes of its active TeX and
bibliography files. It also records the Python version and elapsed verification
time measured for this release. `SHA256SUMS` records the hashes of the files in
the certificate package.

In the manuscript source directory, this document is named
`README-certificates.md`. In the standalone certificate archive, it is the
main `README.md`.

## What the Certificates Establish

In the independent geometric-load tree, each node receives an independent
load `G` with `P[G = k] = 2^(-k-1)` for integers `k >= 0`. Let `Q_h^(Z)`
be the number of blocks passed above the root after upward packing in a tree
with `h` bucket levels and capacity `Z`. These auxiliary trees satisfy

```text
Q_0^(Z) = 0
Q_(h+1)^(Z) = max(0, G + Q_h,left^(Z) + Q_h,right^(Z) - Z),
```

where the three inputs are independent and the child overflows have the law
of `Q_h^(Z)`. Here `h` counts bucket levels. It is not the physical Path ORAM
depth `L`, which counts edges.

For capacity four, `M_h = E[(7/5)^(Q_h^(4))]`. The term `B_h` appears in the
exact recurrence `M_(h+1) = a M_h^2 + B_h`, where `a = 3125/7203`.
The term `B_h` is defined by

```text
u = 7/5
c(k) = sum over g >= 0 of 2^(-g-1) * max(0, 1 - u^(g+k-4))
B_h = E[c(Q_h,left^(4) + Q_h,right^(4))].
```

The sum defining `c(k)` has only finitely many nonzero terms. With
`r = 23/20` and `b = 57621/100000`, the third capacity-four check is

```text
r - a*r*r - b = 19687/720300000 > 0.
```

| Calculation | Exact finite conclusion | Role in the paper |
|---|---|---|
| Capacity four, height 12 | `M_12 < 111/100 < 23/20`, `B_12 < 57621/100000`, and a positive closure margin | Lemma 4.21 verifies the three conditions of Lemma 4.20, which bounds the moment at every height. |
| Capacity three, height 14 | A computed lower mean is greater than `11/2`, hence `E[min(32, Q_14^(3))] > 11/2`. | Lemma 5.5 supplies the local surplus used in the lower proof. |
| Lower-bound error arithmetic | `11/2 - 1/8192 - 3/64 = 44671/8192 > 5` | Together with the error bounds proved in Lemma 5.4 and Appendix C, the local surplus survives transfer to a real scan region. |


The computation establishes finite premises, not a bound on every height by
enumeration. Lemma 4.20 proves the extension to all heights. On the lower side,
the paper proves the distributional lower bound, controls the transfer losses,
and aggregates the real regions' overflow.

## How to Read certificates.json

The top-level objects are `z4` and `z3`. Fractions are stored as exact strings
such as `23/20`. Integer mass numerators are stored as JSON integers.

### Capacity Four: Upper Bounds

`z4.grid = 2^80`, `z4.height = 12`, and `z4.u = 7/5`.
The fields `a`, `r`, and `margin` are the constants defined above.
For each row with height `h`:

| Field | Meaning |
|---|---|
| `probability_upper[k]` | The stored upper bound `P_h(k)` on `P[Q_h^(4) = k]`. |
| `M_upper` | The separately propagated upper bound on `M_h`. |
| `B_upper` | An upper bound on `B_h`, computed from that row's child-sum coefficients. |

At height `h`, the retained indices run from zero through
`K_h = 4*(12-h)+3`. The 13 rows contain 364 probability upper bounds.
The paper proves that omitted indices cannot affect the retained outputs.
Probability and moment updates are rounded upward to the `2^-80` grid.
The coefficients are **upper bounds, not a normalized probability
distribution**. Normalizing them could invalidate the inequalities.

The final moment bound is

```text
1334003097123611627266701 / 1208925819614629174706176 < 111/100.
```

The final `B_upper` is

```text
1010971287715732434488932565388985262516425950508827 /
1754532715615748953303523641675897765096947017842688 < 57621/100000.
```

Lemma 4.21 proves the validity of the updates and their rounding direction.
Appendix B, Algorithm 1, specifies the complete calculation.

### Capacity Three: A Distribution Giving a Lower Bound

`z3.grid = 2^48`, `z3.cap = 32`, and `z3.height = 14`.
Row `h` contains 33 nonnegative integer mass numerators. Entry `k`, divided by
`2^48`, is the probability of value `k` under the computed distribution
`mu_h`. Each row sums to `2^48`. The 15 rows contain 495 integers.

Positive masses are rounded downward and the discarded mass is moved to zero.
The paper proves that every positive tail of `mu_h` is at most the corresponding
tail of `min(32, Q_h^(3))`. The individual masses need not be pointwise lower
bounds. At every update, the program checks agreement between two formulas:
one sums individual outcome probabilities, and the other uses differences
of survival probabilities.

The final field `mean_lower` is the mean of `mu_14`:

```text
775231701030647 / 140737488355328 > 11/2.
```

It is a lower bound on `E[min(32, Q_14^(3))]`, not that expectation's exact
value. Section 5.4 proves this comparison. Appendix B, Algorithm 2, specifies
the finite calculation.

The fields `time_loss_upper = 1/8192` and `coupling_loss_upper = 3/64`
record the bounds used in Appendix C. The program checks their arithmetic
combination. The validity of those error bounds for the stochastic processes
is proved in the paper.

## Optional Commands

To check a certificate at a different path:

```sh
python3 verify_certificates.py --certificate /path/to/certificates.json
```

To regenerate the exact table at a new path:

```sh
python3 verify_certificates.py --write-output /path/to/new-certificates.json
```

Generation checks the arithmetic premises but does not compare against an
existing table. The program refuses to overwrite an existing output file.
Use the default verification command to check the supplied certificate.

## Separate Diagnostic

The manuscript source also contains `check_small_cases.py`. It enumerates two
small instances of the actual Path ORAM process. This diagnostic is not needed
to verify either certificate and supplies no premise of the asymptotic
theorems. It is excluded from the minimal certificate archive. If used, run it
in a disposable source copy, because it overwrites `small-cases.json`.
