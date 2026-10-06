# Computer-assisted proof of a real-valued discretely self-similar blow-up profile in the defocusing energy-supercritical quintic wave equation.



Table of contents:


* [Introduction](#introduction)
* [The defocusing quintic wave equation](#the-defocusing-quintic-wave-equation)
   * [Proof of existence of a discretely self-similar profile](#proof-of-existence-of-a-discretely-self-similar-profile)
   * [A posteriori bounds on the profile](#a-posteriori-bounds-on-the-profile)
* [Utilisation and References](#utilisation-and-references)
* [License and Citation](#license-and-citation)
* [Contact](#contact)



# Introduction

This Julia code is a complement to the article 

#### [[1]](https://arxiv.org/abs/2610.06678) : "Finite-time blow-up for the real-valued defocusing energy-supercritical quintic NLW", M. Cadiot, A. D. Ionescu and S. Palasek, [ArXiv Link](https://arxiv.org/abs/2610.06678)

as it provides the necessary rigorous computations of the bounds presented along the paper. The computations are performed using the package [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl). The mathematical objects (spaces, sequences, operators,...) are built using the package [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl). 


# The defocusing quintic wave equation

The defocusing quintic wave equation
$$\partial_t^2 u - \Delta u + u^5 = 0$$
for real-valued functions on $\mathbb{R}^{1+10}$ is energy-supercritical. In [[1]](https://arxiv.org/abs/2610.06678), we construct a discretely self-similar solution in a backward light cone, of the form
$$u(t,x) = (T-t)^{-1/2}\, W\left(\omega \log\frac{1}{T-t}, \frac{|x|}{T-t}\right),$$
where the profile $W = W(\theta,\rho)$ is $2\pi$-periodic in its first variable. In particular, cutting off the initial data provides smooth, compactly supported, radial data whose solution blows up in finite time (see Theorem 1.1 in [[1]](https://arxiv.org/abs/2610.06678)). The frequency $\omega$ of the profile is an unknown of the problem, and the couple $(\omega,W)$ solves
$$\Lambda(\omega)W + W^5 = 0,$$
where, using the variable $x = 2\rho^2-1 \in [-1,1]$,
$$\Lambda(\omega) = \omega^2\partial_\theta^2 + 2\omega\partial_\theta + \frac{3}{4} - 4(1-x^2)\partial_x^2 - (32-8x)\partial_x + 4\omega(1+x)\partial_\theta\partial_x.$$
We refer to Section 2 and Section 4 in [[1]](https://arxiv.org/abs/2610.06678) for the derivation of this equation and for a complete description of the problem. The profile $W$ is represented by a Fourier series in $\theta$ and by a Chebyshev series in $x$. In particular, $W$ satisfies $W(\theta+\pi,\rho) = -W(\theta,\rho)$, meaning that only the odd Fourier modes are used.

## Proof of existence of a discretely self-similar profile

The code main_proof.jl provides the computer-assisted component for the constructive proof of existence of a discretely self-similar profile, using the analysis of [[1]](https://arxiv.org/abs/2610.06678). Specifically, we look for a zero of
$$F(\omega,W) = \left(\ell(W),~ \Lambda(\omega)W + W^5\right),$$
where $\ell$ is a phase condition removing the invariance by translation in $\theta$. 

We provide a candidate solution for the proof, which is given in the file approximate_solution.jld2. It contains the approximate frequency $\omega$ and the Fourier-Chebyshev coefficients of the approximate profile $W$, both stored in interval arithmetic with BigFloat (256 bits) precision. These correspond to the approximate solution $\bar{U} = (\bar{\omega},\bar{W})$ in Section 5.5. In particular, $\bar{W}$ is real and possesses the Fourier modes $|n_1| \leq 65$ and the Chebyshev modes $n_2 \leq 100$. The approximate profile was obtained by numerical continuation from rotating self-similar profiles of the complex-valued equation (see Section 3). This step is not rigorous and is not part of the present code.

Given this approximate solution, the code computes rigorously the bounds $Y$, $Z_1$ and $Z_2$ of the Newton-Kantorovich approach. The code follows the structure and the notations of Section 5 in [[1]](https://arxiv.org/abs/2610.06678) :

 - the bound $Y$ is computed following Section 5.2 (Lemma 5.1). The nonlinear terms $5W^4$ and $W^5$, as well as the residual, are computed with BigFloat precision.
 - the constants $\kappa$, $\kappa'$ and $\kappa''$, controlling the norm of $\Lambda^{-1}\pi^{>N}$, $\partial_\omega\Lambda \Lambda^{-1}\pi^{>N}$ and $\mathcal{N}^2\Lambda^{-1}\pi^{>N}$, are computed following Section 5.3.1 (Lemmas 5.3 to 5.7 and Lemma 5.9). In particular, the assumptions (5.12) and (5.13) are verified.
 - the finite dimensional part of $DF(\omega,W)$ and its approximate inverse $B_N$ are constructed following Section 4.4. Note that the product $DF(\omega,W) B_N$ is achieved column by column without constructing $DF(\omega,W)$, but rather using its action. This helps reduce the memory cost of the code.
 - the bound $Z_1 = \max\{Z_{1,0} + Z_{1,1}, Z_{1,2}\}$ is computed following Section 5.3 (Lemma 5.2).
 - the bound $Z_2$ is computed following Section 5.4 (Lemma 5.8).

Finally, the code verifies the conditions (4.15) and validates (or not) the computer-assisted proof (cf. Theorem 5.11). If the computer-assisted proof succeeds, the values of $Y$, $Z_1$, $Z_2$ and of the radius $r_0$ are displayed. In particular, we obtain that there exists a zero $(\omega^\ast,W^\ast)$ of $F$ in a neighborhood of the approximate solution, and the code displays explicit bounds for $|\omega^\ast - \omega|$ and $\|W^\ast - W\|$, making the proof constructive. With the values of the parameters given in the code, we obtain $r_0 = 7.216 \cdot 10^{-11}$, $|\omega^\ast - \omega| \leq 4.05 \cdot 10^{-10}$ and $\|W^\ast - W\| \leq 8.42 \cdot 10^{-10}$.

## A posteriori bounds on the profile

If the proof of the profile is achieved, the code computes a posteriori bounds on $W^\ast$, following Remark 5.12 in [[1]](https://arxiv.org/abs/2610.06678). More precisely, we obtain a rigorous upper bound for $|W^\ast(\theta,\rho)|$, valid for all $\theta$ and all $\rho \in [0,1]$, and a rigorous lower bound for 
$$\max\left(|W^\ast(\theta,0)|, |W^\ast(\theta,1)|\right)$$ 
which is valid for all $\theta$. The lower bound is obtained by evaluating $W$ on a uniform grid of 2001 points of $[0,\pi]$ and using a Lipschitz constant in $\theta$. In particular, we validate that
$$2.26 \leq \max_{0 \leq \rho \leq 1} |W^\ast(\theta,\rho)| \leq 8.25$$
for all $\theta$, which provides the blow-up rate stated in Remark 1.2 in [[1]](https://arxiv.org/abs/2610.06678).


 # Utilisation and References

 The interested user needs to download all files in the same folder. The code main_proof.jl can then be run directly : it activates the environment given by Project.toml and Manifest.toml, and installs the required packages if needed. 
 
 The values of the parameters used for the proof are given at the beginning of the section "Computer-assisted proof of the profile" of the code. In particular, $N_1 = 65$ and $N_2 = 200$ are the numbers of Fourier and Chebyshev modes of the finite dimensional part, $K_1 = 99$ and $K_2 = 1500$ are the sizes used for the control of $\Lambda^{-1}$, and $\nu = 1.3$, $\nu_\theta = 1.01$ are the weights of the norm. With these values, the proof requires around 7GB at peak memory and can therefore be run on a laptop. A second set of values is given in comment for a quick test of the code. Note that the proof is not conclusive with these values.

 The code is build using the following packages :
 - [RadiiPolynomial](https://github.com/OlivierHnt/RadiiPolynomial.jl) (version 0.9.11)
 - [IntervalArithmetic](https://github.com/JuliaIntervals/IntervalArithmetic.jl) (version 1.0.10)
 - [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/)
 - [JLD2](https://github.com/JuliaIO/JLD2.jl)
 
 The code has been written for Julia 1.10.
 
 
 # License and Citation
 
  This code is available as open source under the terms of the [MIT License](http://opensource.org/licenses/MIT).
  
If you wish to use this code in your publication, research, teaching, or other activities, please cite it using the following BibTeX template:

```
@software{NLW_Supercritial.jl,
  author = {Matthieu Cadiot and Alexandru D. Ionescu and Stan Palasek},
  title  = {NLW_Supercritial.jl},
  url    = {https://github.com/matthieucadiot/NLW_Supercritial.jl},
  note = {\url{ https://github.com/matthieucadiot/NLW_Supercritial.jl},
  year   = {2026}
}
```

# Contact

You can contact me at :

matthieu.cadiot@polytechnique.edu
