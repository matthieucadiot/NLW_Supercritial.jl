# Computer-assisted proofs of real-valued discretely self-similar blow-up profiles in the defocusing energy-supercritical wave equation.


Table of contents:


* [Introduction](#introduction)
* [The defocusing wave equation](#the-defocusing-wave-equation)
   * [Proof of existence of a discretely self-similar profile](#proof-of-existence-of-a-discretely-self-similar-profile)
   * [A posteriori bounds on the profile](#a-posteriori-bounds-on-the-profile)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)


# Introduction

This Julia code is a complement to the article

#### [[1]](https://arxiv.org/abs/2610.06678) : "Finite-time blow-up for the real-valued defocusing energy-supercritical NLW", M. Cadiot, A. D. Ionescu and S. Palasek, [ArXiv Link](https://arxiv.org/abs/2610.06678)

as it provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl).


# The defocusing wave equation

The defocusing wave equation


```math
\partial_t^2 u - \Delta u + u^p = 0
```


for real-valued functions on $`\mathbb{R}^{1+d}`$ is energy-supercritical for the three cases


```math
(d,p)=(9,3),\quad (7,5),\quad (6,7).
```


In [[1]](https://arxiv.org/abs/2610.06678), we construct a discretely self-similar solution in a backward light cone, of the form


```math
u(t,x) = (T-t)^{-\alpha}\, W\left(\omega \log\frac{1}{T-t}, \frac{|x|}{T-t}\right),\qquad \alpha=\frac{2}{p-1},
```


where the profile $`W = W(\theta,\rho)`$ is $`2\pi`$-periodic in its first variable. In particular, cutting off the initial data provides smooth, compactly supported, radial data whose solution blows up in finite time (see Theorem 1.1 in [[1]](https://arxiv.org/abs/2610.06678)). The frequency $`\omega`$ of the profile is an unknown of the problem, and the couple $`(\omega,W)`$ solves


```math
\Lambda(\omega)W + W^p = 0,
```


where, using the variable $`x = 2\rho^2-1 \in [-1,1]`$,


```math
\begin{aligned}
\Lambda(\omega)={}&\omega^2\partial_\theta^2+(2\alpha+1)\omega\partial_\theta+\alpha(\alpha+1)\\
&-4(1-x^2)\partial_x^2-\big(4d-(4\alpha+6)(1+x)\big)\partial_x\\
&+4\omega(1+x)\partial_\theta\partial_x.
\end{aligned}
```


We refer to Section 2 and Section 4 in [[1]](https://arxiv.org/abs/2610.06678) for the derivation of this equation and for a complete description of the problem. The profile $`W`$ is represented by a Fourier series in $`\theta`$ and by a Chebyshev series in $`x`$. In particular, $`W`$ satisfies $`W(\theta+\pi,\rho) = -W(\theta,\rho)`$, meaning that only the odd Fourier modes are used.

## Proof of existence of a discretely self-similar profile

We provide three compressed files, one for each of the cases $`(d,p)=(9,3)`$, $`(7,5)`$ and $`(6,7)`$. Each archive contains its own `main_proof.jl` and the approximate solution for that case. The three scripts follow the same proof, using the analysis of [[1]](https://arxiv.org/abs/2610.06678), with the corresponding values of $`d`$ and $`p`$.

For each case, the unknowns are the frequency $`\omega`$ and the profile $`W`$. The script proves that the following system has a solution near the saved approximate solution :


```math
F(\omega,W)=
\begin{pmatrix}
\ell(W)\\[2pt]
\Lambda(\omega)W+W^p
\end{pmatrix}
=0.
```


The second component is the profile equation. The first component, $`\ell(W)=0`$, is a phase condition fixing the freedom to replace $`W(\theta,\rho)`$ by a translate $`W(\theta+\theta_0,\rho)`$.

We provide a candidate solution for each proof, which is given in the file `approximate_solution_d9_p3.jld2`, `approximate_solution_d7_p5.jld2` or `approximate_solution_d6_p7.jld2`, respectively. Each file contains the approximate frequency $`\omega`$ and the Fourier-Chebyshev coefficients of the approximate profile $`W`$, both stored in interval arithmetic with BigFloat precision, as well as the values of $`d`$ and $`p`$. These correspond to the approximate solution $`\bar{U} = (\bar{\omega},\bar{W})`$ in Section 5.5. In particular, $`\bar{W}`$ is real and possesses the Fourier modes $`|n_1| \leq N_{W,1}`$ and the Chebyshev modes $`n_2 \leq N_{W,2}`$, with $`N_W=(103,1290)`$, $`(101,600)`$ and $`(145,600)`$ for the three cases, respectively. The approximate profiles were obtained by numerical continuation from rotating self-similar profiles of the complex-valued equation (see Section 3). This step is not rigorous and is not part of the present code.

Given this approximate solution, the code computes rigorously the bounds $`Y`$, $`Z_1`$ and $`Z_2`$ of the Newton-Kantorovich approach. The code follows the structure and the notations of Section 5 in [[1]](https://arxiv.org/abs/2610.06678) :

 - the bound $`Y`$ is computed following Section 5.2 (Lemma 5.1). The nonlinear terms $`pW^{p-1}`$ and $`W^p`$, as well as the residual, are computed with BigFloat precision (192 bits).
 - the constants $`\kappa`$, $`\kappa'`$ and $`\kappa''`$, controlling the norm of $`\Lambda^{-1}\pi^{>N}`$, $`\partial_\omega\Lambda \Lambda^{-1}\pi^{>N}`$ and $`\mathcal{N}^2\Lambda^{-1}\pi^{>N}`$, are computed following Section 5.3.1 (Lemmas 5.3 to 5.7 and Lemma 5.9). In particular, the assumptions (5.11) and (5.12) are verified.
 - the finite dimensional part of $`DF(\omega,W)`$ and its approximate inverse $`B_N`$ are constructed following Section 4.4. The product $`DF(\omega,W) A_N`$ is achieved by blocks without constructing the whole matrix $`DF(\omega,W)`$, but rather using its action. This helps reduce the memory cost of the code.
 - the bound $`Z_1 = \max\{Z_{1,0} + Z_{1,1}, Z_{1,2}\}`$ is computed following Section 5.3 (Lemma 5.2).
 - the bound $`Z_2`$ is computed following Section 5.4 (Lemma 5.8).

Finally, the code verifies the conditions (4.15) and validates the computer-assisted proof (cf. Theorem 5.11). If the computer-assisted proof succeeds, the values of $`Y`$, $`Z_1`$, $`Z_2`$ and of the radius $`r_0`$ are displayed. In particular, we obtain that there exists a zero $`(\omega^\ast,W^\ast)`$ of $`F`$ in a neighborhood of the approximate solution, and the code displays explicit bounds for $`|\omega^\ast - \bar\omega|`$ and $`\lVert W^\ast - \bar W\rVert`$, making the proof constructive. With the values of the parameters given in the code, we obtain the bounds presented in Tables 2 and 3 of [[1]](https://arxiv.org/abs/2610.06678). The code stops with an error if one of the required conditions is not verified.

## A posteriori bounds on the profile

If the proof of the profile is achieved, the code computes a posteriori bounds on $`W^\ast`$, following Remark 5.13 in [[1]](https://arxiv.org/abs/2610.06678). More precisely, we obtain a rigorous upper bound for $`|W^\ast(\theta,\rho)|`$, valid for all $`\theta`$ and all $`\rho \in [0,1]`$, and a rigorous lower bound for


```math
\max\left(|W^\ast(\theta,0)|, |W^\ast(\theta,1)|\right)
```


which is valid for all $`\theta`$. The lower bound is obtained by evaluating $`W`$ on a uniform grid of 2001 points of $`[0,\pi]`$ and using a Lipschitz constant in $`\theta`$. In particular, we validate that


```math
\begin{array}{ll}
27.55 \leq \displaystyle\max_{0 \leq \rho \leq 1}|W^\ast(\theta,\rho)| \leq 336.9, & (d,p)=(9,3),\\[4pt]
1.244 \leq \displaystyle\max_{0 \leq \rho \leq 1}|W^\ast(\theta,\rho)| \leq 23.07, & (d,p)=(7,5),\\[4pt]
1.008 \leq \displaystyle\max_{0 \leq \rho \leq 1}|W^\ast(\theta,\rho)| \leq 11.26, & (d,p)=(6,7),
\end{array}
```


for all $`\theta`$, which provides the blow-up rate stated in Remark 1.2 in [[1]](https://arxiv.org/abs/2610.06678).


# Utilisation and References

The interested user needs to download and extract the compressed file corresponding to the desired proof :

 - [NLW_proof_d9_p3.zip](NLW_proof_d9_p3.zip) for $`(d,p)=(9,3)`$,
 - [NLW_proof_d7_p5.zip](NLW_proof_d7_p5.zip) for $`(d,p)=(7,5)`$,
 - [NLW_proof_d6_p7.zip](NLW_proof_d6_p7.zip) for $`(d,p)=(6,7)`$.

The four files `main_proof.jl`, the corresponding approximate solution, `Project.toml` and `Manifest.toml` must be kept in the same folder. The code `main_proof.jl` can then be run directly : it activates the environment given by `Project.toml` and `Manifest.toml`, and installs the required packages if needed. For instance, after extracting the first archive, one can run `julia --threads=32 NLW_proof_d9_p3/main_proof.jl`.

The values of the parameters used for the proof are given at the beginning of the section "Stored approximation and parameters" of the code. In particular, $`N=(65,1400)`$, $`(65,900)`$ and $`(95,1000)`$ are the Fourier and Chebyshev cutoffs of the finite dimensional part, and $`K=(165,10500)`$, $`(203,7000)`$ and $`(205,8000)`$ are the cutoffs used for the control of $`\Lambda^{-1}`$, for the three cases, respectively. The weights of the norm are $`\nu=1.1`$ and $`\nu_\theta=1.01`$ in all three proofs.

The proofs require a high-memory machine. The floating point inverse alone occupies approximately 68.4 GB, 28.3 GB and 73.9 GB for the three cases, respectively, and the peak memory is higher because of the interval computations. The proofs took 9.4, 4.2 and 14.5 hours on 32 cores (see Remark 5.12 in [[1]](https://arxiv.org/abs/2610.06678)).

The code is built using the following packages :

 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) (version 0.9.11)
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl) (version 1.0.10)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [JLD2](https://github.com/JuliaIO/JLD2.jl) (version 0.6.7)

The code can be run with Julia 1.12. The provided manifests were generated with Julia 1.12.4.


# License and Citation

This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).

If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{NLW_Supercritial.jl,
  author = {Matthieu Cadiot and Alexandru D. Ionescu and Stan Palasek},
  title  = {NLW_Supercritial.jl},
  url    = {https://github.com/matthieucadiot/NLW_Supercritial.jl},
  year   = {2026}
}
```

# Contact

You can contact me at :

matthieu.cadiot@polytechnique.edu
