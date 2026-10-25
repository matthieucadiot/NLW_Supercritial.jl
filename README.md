# Computer-assisted proofs of discretely self-similar blow-up for the defocusing NLW

This repository contains the Julia codes and approximate profiles for the computer-assisted proofs in

**Finite-time blow-up for the real-valued defocusing energy-supercritical NLW**,  
M. Cadiot, A. D. Ionescu and S. Palasek.  
[Preprint](https://arxiv.org/abs/2610.06678)

The packages cover the three cases

$$
(d,p)=(9,3),\qquad (7,5),\qquad (6,7)
$$

of the real-valued defocusing wave equation

$$
\partial_t^2u-\Delta u+u^p=0.
$$

The profile is periodic in similarity time and real-analytic across the light cone. It gives a solution

$$
u(t,x)=(T-t)^{-\alpha}W^*\left(\omega^*\log\frac1{T-t},\frac{|x|}{T-t}\right),
\qquad \alpha=\frac2{p-1}.
$$

The analytic reduction from a validated profile to blow-up of smooth, compactly supported radial data is Proposition 2.4 of the paper. The numerical profile validation is Theorem 5.11. References below follow the manuscript accompanying these three proof packages.

## Proof packages

Download the archive for the desired case. Each package is self-contained.

| $(d,p)$ | Archive | Approximation $N_W$ | Finite inverse $N$ | Tail cutoff $K$ |
|---|---|---|---|---|
| $(9,3)$ | [NLW_proof_d9_p3.zip](NLW_proof_d9_p3.zip) | $(103,1290)$ | $(65,1400)$ | $(165,10500)$ |
| $(7,5)$ | [NLW_proof_d7_p5.zip](NLW_proof_d7_p5.zip) | $(101,600)$ | $(65,900)$ | $(203,7000)$ |
| $(6,7)$ | [NLW_proof_d6_p7.zip](NLW_proof_d6_p7.zip) | $(145,600)$ | $(95,1000)$ | $(205,8000)$ |

Each extracted folder contains four files:

- `main_proof.jl`: the complete computer-assisted proof for that case.
- `approximate_solution_d9_p3.jld2`, `approximate_solution_d7_p5.jld2`, or `approximate_solution_d6_p7.jld2`: the corresponding approximate profile and frequency, together with `d` and `p`.
- `Project.toml`: the Julia environment.
- `Manifest.toml`: the pinned dependency versions.

The saved profile uses Fourier modes in $\theta$ and Chebyshev polynomials in $x=2\rho^2-1$. Only odd Fourier modes occur, so $W(\theta+\pi,\rho)=-W(\theta,\rho)$. The profile and frequency are stored as BigFloat intervals. The nonlinear products and residual are evaluated with 192-bit arithmetic.

## Running a proof

Use Julia 1.12; the supplied manifests were generated with Julia 1.12.4. Run one of the following commands from the directory containing the downloaded archive:

```bash
unzip NLW_proof_d9_p3.zip && julia --threads=32 --project=NLW_proof_d9_p3 NLW_proof_d9_p3/main_proof.jl
```

```bash
unzip NLW_proof_d7_p5.zip && julia --threads=32 --project=NLW_proof_d7_p5 NLW_proof_d7_p5/main_proof.jl
```

```bash
unzip NLW_proof_d6_p7.zip && julia --threads=32 --project=NLW_proof_d6_p7 NLW_proof_d6_p7/main_proof.jl
```

The script activates its own environment and installs the pinned dependencies on the first run. Subsequent runs need only the `julia` command. Set `--threads` to the number of CPU cores available to the computation.

These are full proofs at the parameters of Tables 2–3 and require a high-memory machine or cluster allocation. The floating-point inverse alone occupies approximately the memory listed below; peak memory is higher because the interval computations also require working arrays.

| $(d,p)$ | Finite matrix dimension | Dense inverse alone | Reported time on 32 cores |
|---|---:|---:|---:|
| $(9,3)$ | 92,467 | 68.4 GB | 9.4 hours |
| $(7,5)$ | 59,467 | 28.3 GB | 4.2 hours |
| $(6,7)$ | 96,097 | 73.9 GB | 14.5 hours |

The timings are those reported in Remark 5.12 of the paper; they depend on the hardware and software environment.

## What the code verifies

The unknown is $U=(\omega,W)$, and the profile equation, with its phase condition, is

$$
F(U)=\bigl(\ell(W),\Lambda(\omega)W+W^p\bigr)=0.
$$

Starting from the saved approximation $\bar U$, the code follows Sections 4–5 of the paper:

1. Compute the full residual bound $Y\geq\|F(\bar U)\|$, including every generated nonlinear mode: Lemma 5.1.
2. Validate the finite inverse columns and the infinite-tail bounds $\kappa,\kappa',\kappa''$: Section 5.3.1 and Lemmas 5.7 and 5.9. This includes hypotheses (5.11)–(5.12).
3. Construct the finite approximate inverse and enclose the defect $I-DF(\bar U)A$, yielding $Z_1$: Section 4.4 and Lemma 5.2.
4. Bound the nonlinear variation by $Z_2(r)$: Lemma 5.8.
5. Verify the radii-polynomial conditions (4.15) and apply Theorem 4.5.

The correction is written as $U^*=\bar U+A\Psi$. The validated radius bounds $\|\Psi\|$; the frequency and profile errors are obtained by multiplying it by $\|\pi_\omega A\|$ and $\|\pi_W A\|$.

A successful run prints `EXISTENCE PROVEN`, the bounds $Y,Z_1,Z_2$, a validated radius, and explicit frequency and profile error bounds. Failed validation conditions raise an error.

For comparison, Table 2 of the paper gives the following error bounds:

| $(d,p)$ | $\lvert\omega^*-\bar\omega\rvert$ | $\lVert W^*-\bar W\rVert_{\ell^1}$ |
|---|---:|---:|
| $(9,3)$ | $1.39\times10^{-16}$ | $2.30\times10^{-14}$ |
| $(7,5)$ | $4.93\times10^{-13}$ | $3.71\times10^{-12}$ |
| $(6,7)$ | $3.12\times10^{-11}$ | $1.25\times10^{-10}$ |

Finally, the code verifies the uniform profile-size bounds of Remark 5.13. It evaluates the traces at $\rho=0,1$ on 2001 time points and uses a Lipschitz bound to cover all $\theta$. Together with the validated profile error, this supplies the upper and positive lower bounds used for the self-similar blow-up rate.

## Dependencies

The pinned environment includes:

- [RadiiPolynomial.jl](https://github.com/OlivierHnt/RadiiPolynomial.jl) 0.9.11, for sequence spaces, operators and interval FFTs.
- [IntervalArithmetic.jl](https://github.com/JuliaIntervals/IntervalArithmetic.jl) 1.0.10, for rigorous enclosures.
- [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) 0.6.7, for loading the approximate profiles.
- Julia's `LinearAlgebra` standard library.

## License and citation

This code is distributed under the [MIT License](LICENSE).

Please cite the paper above and this repository when using the proofs or data:

```bibtex
@software{CadiotIonescuPalasek2026NLW,
  author = {Cadiot, Matthieu and Ionescu, Alexandru D. and Palasek, Stan},
  title  = {NLW_Supercritial.jl},
  year   = {2026},
  url    = {https://github.com/matthieucadiot/NLW_Supercritial.jl}
}
```

## Contact

[matthieu.cadiot@polytechnique.edu](mailto:matthieu.cadiot@polytechnique.edu)
