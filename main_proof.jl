using Pkg ; Pkg.activate(@__DIR__)
isfile(joinpath(@__DIR__, "Manifest.toml")) || Pkg.Registry.update()
Pkg.resolve() ; Pkg.instantiate()
using RadiiPolynomial, IntervalArithmetic, LinearAlgebra, JLD2
using Logging ; disable_logging(Logging.Info)


################ Some basic functions ##################
function big2float(U)
    #### rigorous enclosure in Float64 of a sequence U with complex coefficients in Interval{BigFloat}
    Ur = interval.(Float64.(inf.(real.(U)),RoundDown),Float64.(sup.(real.(U)),RoundUp))
    Ui = interval.(Float64.(inf.(imag.(U)),RoundDown),Float64.(sup.(imag.(U)),RoundUp))
    return complex(Ur,Ui)
end

function weight_nu(ν,N)
    #### weights of the ℓ¹ norm with geometric weight ν for the Chebyshev coefficients 0,...,N. The factor 2 comes from the convention a_0 + 2∑ a_k T_k of RadiiPolynomial.
    w = interval.(zeros(N+1))
    w[1] = interval(1)
    for i=1:N
        w[i+1] = interval(2)*ν^interval(i)
    end
    return w
end

function norm_cols(B,w1,w2)
    #### norms of the columns of the matrix B, with the weights w1 for the rows and w2 for the columns : the maximum is the operator norm of B for the weighted ℓ¹ norms
    #### the entries are scaled by the weights before taking the absolute values : the weights are of order 1e171 for 1500 Chebyshev modes, and the absolute value of a complex number uses squares, which would overflow
    return vec(sum(abs.(w1.*B./w2'), dims=1))
end

function norm_col(B,k,ν)
    #### ℓ¹_ν operator norm of the column k ≥ 1 of an operator B defined on Chebyshev sequences, that is ‖B π_k‖ (only used in main_proof_old.jl)
    return norm(Sequence(codomain(B), B[:,k]), ℓ¹(GeometricWeight(ν)))/(interval(2)*ν^interval(k))
end

#################### Computer-assisted proof of the profile ######################

t_start = time()

N1 = 65 ; N2 = 200 ; K1 = 99 ; K2 = 1500  ; ν = I"1.3" ; νθ = I"1.01"   #### values for the proof.
# N1 = 5 ; N2 = 30 ; K1 = 99 ; K2 = 150   #### values for a quick test of the code (the proof is not conclusive with these values)
NW1 = min(N1,65) ; NW2 = min(N2,100) ; NV1 = N1+10 ; NV2 = N2+10 ; ν0 = I"1.8"
setprecision(256)
if iseven(N1)||iseven(K1)||iseven(NW1)
    error("N1, K1 and NW1 have to be odd : only the odd Fourier modes are used")
end
if (K1 <= N1)||(K2 <= N2)||(NW1 > N1)||(NW2 > N2)
    error("The sizes are not compatible : K1 > N1, K2 > N2, NW1 ≤ N1 and NW2 ≤ N2 are needed")
end

### spaces of sequences : Fourier in θ and Chebyshev in x = 2ρ²-1. W is in SW, V in SV and S is the space of the finite part
SW = Fourier(NW1,interval(1))⊗Chebyshev(NW2)
S = Fourier(N1,interval(1))⊗Chebyshev(N2)
SV = Fourier(NV1,interval(1))⊗Chebyshev(NV2)
### 5W^4 and W^5 are computed for all their Fourier modes and for the Chebyshev modes n2 ≤ NV2 and n2 ≤ N2
S4 = Fourier(4*NW1,interval(1))⊗Chebyshev(NV2)
S5 = Fourier(max(5*NW1,N1),interval(1))⊗Chebyshev(N2)
### spaces for the Chebyshev coefficients of one Fourier mode
SN = Chebyshev(N2) ; SK = Chebyshev(K2)
### norms : ℓν for Chebyshev sequences, ℓθν for Fourier-Chebyshev sequences, ℓθν0 for the a priori bounds
ℓν = ℓ¹(GeometricWeight(ν)) ; ℓθν = ℓ¹(GeometricWeight(νθ),GeometricWeight(ν)) ; ℓθν0 = ℓ¹(GeometricWeight(νθ),GeometricWeight(ν0))

### The file approximate_solution.jld2 contains the approximate solution (ω,W) in interval arithmetic : W is a sequence of Fourier(65,interval(1))⊗Chebyshev(300) with coefficients in Complex{Interval{BigFloat}}, with the convention a_0 + 2∑ a_k T_k of RadiiPolynomial, and ω is an Interval{BigFloat}.
### W only contains the modes n1 ≤ NW1 and n2 ≤ NW2
sol = load(joinpath(@__DIR__, "approximate_solution.jld2"))
W_big = project(sol["W"], SW) ; ω_big = sol["ω"]
#### the approximate solution is saved in bigfloats, we save a copy in float for computations where the bigfloats are not needed
W = big2float(W_big) ; ω = interval(Float64(inf(ω_big),RoundDown),Float64(sup(ω_big),RoundUp))

println("
" * "="^74)
println("CONSTRUCTIVE EXISTENCE PROOF of a discretely self-similar profile for the quintic NLW")
println("-"^74)


############ Section 5.2 : computation of the Y bound (Lemma 5.1). Nonlinear terms and residual, computed in bigfloats ############

println("  Computing 5W^4, W^5 and the residual with bigfloats ($(precision(BigFloat)) bits) ...")

### values of W on a grid, the grid is fine enough to recover exactly the coefficients of W^4 and W^5
W_grid = fft(W_big, fft_size(Fourier(5*NW1,interval(1))⊗Chebyshev(5*NW2)))
V_big = interval(5)*ifft!(W_grid.^4, S4)  ### 5W^4 for n2 ≤ NV2, so that DN(W) is the multiplication operator by 5W^4
W5_big = ifft!(W_grid.^5, S5)  ### W^5 for n2 ≤ N2

#### F(ω,W) = (ℓ(W), Λ(ω)W + W^5), where Λ(ω) = ω²∂θ² + 2ω∂θ + 3/4 - 4(1-x²)∂x² - (32-8x)∂x + 4ω(1+x)∂θ∂x. ΛW_W5_big contains Λ(ω)W + π^{≤(∞,N2)}W^5, as in Lemma 5.1
x = Sequence(Fourier(0,interval(1))⊗Chebyshev(1), [interval(0) ; I"0.5"])  ### the function x
∂xW = differentiate(W_big,(0,1)) ; ∂xxW = differentiate(∂xW,(0,1)) ; ∂θW = differentiate(W_big,(1,0)) ; ∂θθW = differentiate(W_big,(2,0)) ; ∂θxW = differentiate(W_big,(1,1))
ΛW_W5_big = ω_big^2*∂θθW + interval(2)*ω_big*∂θW + I"0.75"*W_big - interval(4)*(interval(1) - x*x)*∂xxW - (interval(32) - interval(8)*x)*∂xW + interval(4)*ω_big*(interval(1) + x)*∂θxW + W5_big
#### ∂ωΛ(ω)W = 2ω∂θ²W + 2∂θW + 4(1+x)∂θ∂xW and 𝒩²W = -∂θ²W
dΛW_big = interval(2)*ω_big*∂θθW + interval(2)*∂θW + interval(4)*(interval(1) + x)*∂θxW
#### phase condition ℓ(W) = Im ∑ w_{1,n2}, that is the imaginary part of the first Fourier mode at x = 1
ℓW = imag(Sequence(Chebyshev(NW2), W_big[(1,0:NW2)])(interval(1))) ; ℓW = interval(Float64(inf(ℓW),RoundDown),Float64(sup(ℓW),RoundUp))

### copies in float. From now on, the bigfloats are not used anymore : we free the memory of the grid
ΛW_W5 = big2float(ΛW_W5_big) ; V = big2float(V_big) ; dΛW = project(big2float(dΛW_big),S) ; 𝒩²W = big2float(-∂θθW)
W_grid = nothing ; GC.gc()

#### A priori bounds (5.7) for the Chebyshev modes which are not computed : for every sequence g, ‖g - π^{≤(∞,L)}g‖ ≤ (ν/ν0)^(L+1) ‖g‖_0, and ‖W^k‖_0 ≤ ‖W‖_0^k (Banach algebra), where ‖·‖_0 is the norm with the weight ν0
#### tail_W5 ≥ ‖W^5 - π^{≤(∞,N2)}W^5‖  and  tail_5W4 ≥ ‖5W^4 - π^{≤(∞,NV2)}5W^4‖
norm_W = norm(W,ℓθν) ; norm_W_0 = norm(W,ℓθν0)
tail_W5 = (ν/ν0)^interval(N2+1)*norm_W_0^5 ; tail_5W4 = interval(5)*(ν/ν0)^interval(NV2+1)*norm_W_0^4
norm_5W4 = norm(V,ℓθν) + tail_5W4 ; norm_𝒩²W = norm(𝒩²W,ℓθν)
#### In the bound Z1, 5W^4 is replaced by its truncation V of size NV and the error ‖5W^4 - V‖ is added
norm_5W4_V = norm(tail(V,(NV1,NV2)),ℓθν) + tail_5W4 ; V = project(V,SV)

println("‖W‖ = $norm_W ,  ‖5W^4‖ ≤ $norm_5W4 ,  ‖𝒩²W‖ = $norm_𝒩²W")
println("‖5W^4 - V‖ ≤ $norm_5W4_V ,  a priori bounds : ‖5W^4 - π^{≤(∞,NV2)}5W^4‖ ≤ $tail_5W4 ,  ‖W^5 - π^{≤(∞,N2)}W^5‖ ≤ $tail_W5")

#### Lemma 5.1 : Y = |ℓ(W)| + ‖Λ(ω)W + π^{≤(∞,N2)}W^5‖ + (ν/ν0)^(N2+1) ‖W‖_0^5 ≥ ‖F(ω,W)‖
Y = abs(ℓW) + norm(ΛW_W5,ℓθν) + tail_W5
println("Bound Y : Y = $Y   (|ℓ(W)| ≤ $(abs(ℓW)) , ‖Λ(ω)W + π^{≤(∞,N2)}W^5‖ ≤ $(norm(ΛW_W5,ℓθν)))")


############ Section 5.3.1 : computation of the norm of Λ^{-1}π^{>N} (Lemmas 5.3 to 5.7), and of ∂ωΛ Λ^{-1}π^{>N} and 𝒩²Λ^{-1}π^{>N} (Lemma 5.9) : constants κ, κ′, κ′′ ############

#### constants β (5.11), θ1, θ2, θ3, θ3′ (5.22), η, η′ (5.21), and verification of the assumptions (5.12) and (5.13)
γp(n) = sqrt(interval(64) + interval(n)^2*ω^2)
β(n) = interval(4)*(interval(18) + (ν-interval(1))*γp(n))/(ν-interval(1))^2
θ1 = (interval(1)+ν)^2/(ν^2*(ν-interval(1))*interval(K2+1))*(γp(K1) + (interval(1)+ν)*(interval(4)*interval(K2)^2 + interval(K1)^2*ω^2)/(interval(16)*interval(K2-1)))
θ2 = (interval(1)+ν)^2*(interval(4)*interval(K2+1)^2 + interval(K1)^2*ω^2)/(interval(8)*ν*(ν-interval(1))*interval(K2)*interval(K2+1))
θ3 = interval(1)/(interval(4)*interval(K2-1)^2)*(interval(1) + β(K1)/interval(K2+1))
θ3′ = interval(1)/interval(K2+1)*(interval(1) + interval(4)/(ν-interval(1)) + β(K1)/interval(K2-1))
η = interval(1)/ω^2*(interval(1) + β(K1+2)/(interval(K1+2)*ω))
η′ = interval(1)/ω*(interval(2) + interval(4)/(ν-interval(1)) + interval(2)*β(K1+2)/(interval(K1+2)*ω))

sK = interval(K1+2)*ω   #### s_{K1+2} in (5.13)
H1 = interval(36)/((ν-interval(1))*interval(K2-1)) + sqrt(interval(1) + interval(32)/interval(K2-1) + interval(288)/interval(K2-1)^2 + interval(6)/interval(K2-1)^3)
H2 = interval(36)/((ν-interval(1))*sK) + sqrt(interval(1) + I"15.25"/sK + interval(288)/sK^2 + I"15.6"/sK^3)
H = (sup(H1) <= inf((interval(1)+ν)/interval(2)))&&(sup(H2) <= inf((interval(1)+ν)/interval(2)))&&(K2 > N2)&&(K1 > N1)&&(NW1 <= N1)&&(NW2 <= N2)&&(inf(ν0) > sup(ν))
if H
    println("The assumptions (5.12) and (5.13) are satisfied")
else
    println("The assumptions (5.12) and (5.13) are not satisfied : the proof will not be conclusive. For ν = 1.3, K1 ≥ 83 and K2 ≥ 908 are needed")
end

#### Λ_n(ω) = L0 + i n ω L1 - n²ω² and ∂ωΛ_n(ω) = i n L1 - 2n²ω on the Chebyshev modes 0,...,K2, where L0 = -4(1-x²)∂x² - (32-8x)∂x + 3/4 and L1 = 4(1+x)∂x + 2. These matrices are upper triangular and exact.
Id = LinearOperator(SK, SK, interval.(1.0*I[1:K2+1, 1:K2+1]))
xc = Sequence(Chebyshev(1), [interval(0) ; I"0.5"])
D = project(Derivative(1), SK, SK, Interval{Float64}) ; D2 = D*D
Mx = project(Multiplication(xc), SK, SK) ; Mx2 = project(Multiplication(xc*xc), SK, SK)
L0 = -interval(4)*(D2 - Mx2*D2) - interval(32)*D + interval(8)*Mx*D + I"0.75"*Id
L1 = interval(4)*(D + Mx*D) + interval(2)*Id

println("  Computing the operators Λ_n^{-1}, n = 1,3,...,$K1 ...")
#### For n = 1,3,...,K1 we construct an approximate inverse X of Λ_n and control Λ_n^{-1} with a Neumann series. The operators for -n are the complex conjugates.
#### κ, κ′ and κ′′ (Lemma 5.7 and Lemma 5.9) are bounds for the norms of the columns of Λ^{-1}, ∂ωΛ Λ^{-1} and 𝒩²Λ^{-1} which are not in π^{≤N}. We initialize with the bounds for the modes |n| > K1 (Lemma 5.5)
wK = weight_nu(ν,K2)
κ = η/interval(K1+2)^2 ; κ′ = η′ ; κ′′ = η ; ε_max = interval(0)
for n = 1:2:K1
    global κ, κ′, κ′′, ε_max
    Λ = L0 + (interval(im)*interval(n)*ω)*L1 - (interval(n)*ω)^2*Id
    dΛ = (interval(im)*interval(n))*L1 - (interval(2)*interval(n)^2*ω)*Id
    ### approximate inverse and its defect : Λ^{-1} = X(Id - E)^{-1} on the modes 0,...,K2
    X = interval.(inv(mid.(Λ)))
    E = Id - Λ*X ; ε = opnorm(E,ℓν)
    ε_max = maximum([ε_max ε])
    ### λ and λ′ : the bounds (5.8) and (5.9) for the norms of the columns 0,...,K2 of Λ^{-1} and ∂ωΛ Λ^{-1}. The second term is the error of the Neumann series
    λ = norm_cols(coefficients(X), wK, wK) ; λ = λ .+ maximum(λ)*ε/(interval(1)-ε)
    λ′ = norm_cols(coefficients(dΛ*X), wK, wK) ; λ′ = λ′ .+ maximum(λ′)*ε/(interval(1)-ε)
    ### we only keep the columns which are not in π^{≤N} (n2 > N2 if n ≤ N1), and we add the bound of Lemma 5.6 for the columns n2 > K2, which uses the columns K2-1 and K2
    k0 = n <= N1 ? N2+2 : 1
    λ = [λ[k0:end] ; θ1*λ[K2]+θ2*λ[K2+1]+θ3] ; λ′ = [λ′[k0:end] ; θ1*λ′[K2]+θ2*λ′[K2+1]+interval(n)*θ3′]
    κ = maximum([κ ; λ]) ; κ′ = maximum([κ′ ; λ′]) ; κ′′ = maximum([κ′′ ; interval(n)^2*λ])
end
println("Defect of the approximate inverses of Λ_n : ε ≤ $ε_max")
println("κ = $κ ≥ ‖Λ^{-1}π^{>N}‖ ,  κ′ = $(κ′) ≥ ‖∂ωΛ Λ^{-1}π^{>N}‖ ,  κ′′ = $(κ′′) ≥ ‖𝒩²Λ^{-1}π^{>N}‖")


############ Section 4.4 : finite dimensional part, DF and its approximate inverse B_N ############

println("  Computing DF and its approximate inverse B_N (dimension $((N1+1)*(N2+1)+1)) ...")
#### The vectors of ℂ × π^{≤N}ℓ¹ are ordered as (ω, (W_{n,n2})_{0≤n2≤N2} for n = -N1,-N1+2,...,N1). ind(n) gives the indices of the Fourier mode n
d1 = N1+1 ; d2 = N2+1 ; d = d1*d2
ind(n) = div(n+N1,2)*d2 .+ (2:d2+1)
#### weights for the norm on X
w = [interval(1) ; vec([νθ^interval(abs(n))*wk for wk = weight_nu(ν,N2), n = -N1:2:N1])]

#### MV[k+NV1+1] is the multiplication by the Fourier mode k of V, it maps the Chebyshev modes 0,...,N2 to 0,...,N2+NV2
MV = [project(Multiplication(Sequence(Chebyshev(NV2), V[(k,0:NV2)])), SN, Chebyshev(N2+NV2)) for k = -NV1:NV1]
#### Λ_n(ω) and ∂ωΛ_n(ω) on the Chebyshev modes 0,...,N2
L0N = coefficients(project(L0, SN, SN)) ; L1N = coefficients(project(L1, SN, SN)) ; IdN = interval.(1.0*I[1:d2, 1:d2])
ΛN(n) = L0N + (interval(im)*interval(n)*ω)*L1N - (interval(n)*ω)^2*IdN
dΛN(n) = (interval(im)*interval(n))*L1N - (interval(2)*interval(n)^2*ω)*IdN
ξ = project(Evaluation(interval(1)), SN, ParameterSpace(), Interval{Float64})  ### evaluation at x = 1

#### DN is the matrix of π^{≤N}DF(ω,W)π^{≤N}, where 5W^4 is replaced by V : Λ + multiplication by V for W, ∂ωΛW for ω, and the phase condition (first row).
#### DN is not stored (it would need 13 GB) : DN_rows(R) computes its rows R, where R is a block of indices of blk (2 Fourier modes, and the component ω for the first block)
blk = [(n0 == -N1 ? 1 : ind(n0)[1]):ind(min(n0+2,N1))[end] for n0 = -N1:4:N1]
function DN_rows(R)
    B = zeros(Complex{Interval{Float64}}, length(R), d+1)
    if R[1] == 1
        B[1,ind(1)] = -I"0.5"*interval(im)*vec(coefficients(ξ)) ; B[1,ind(-1)] = I"0.5"*interval(im)*vec(coefficients(ξ))
    end
    for n = filter(n -> issubset(ind(n),R), -N1:2:N1)
        r = ind(n) .- (R[1]-1)
        for m = filter(isodd, max(-N1,n-NV1):min(N1,n+NV1))
            B[r,ind(m)] = coefficients(MV[n-m+NV1+1])[1:d2,:]
        end
        B[r,ind(n)] += ΛN(n)
        B[r,1] = dΛW[(n,0:N2)]
    end
    return B
end

##### approximate inverse B_N, stored in floats (BN). BN first contains the midpoint of DN, scaled with the weights of the norm (the entries of DN are very unbalanced, and the scaled matrix is well conditioned), and it is replaced by its inverse
BN = zeros(ComplexF64, d+1, d+1)
for R = blk
    BN[R,:] = mid.(DN_rows(R))
end
wm = mid.(w)
BN .= wm.*BN./wm'
BN = LinearAlgebra.inv!(lu!(BN))   ### inverse in place : inv(BN) would need a copy of the matrix
BN .= BN.*wm'./wm
GC.gc()

println("-"^74)
println("  Computing the Newton–Kantorovich bounds (Z₁, Z₂) ...   ($(round((time()-t_start)/60, digits=1)) min)")
println("-"^74)


############ Section 5.3 : computation of the Z1 bound (Lemma 5.2) ############

#### Z1_0 ≥ ‖π^{≤N} - π^{≤N}DF(ω,W)B_N π^{≤N}‖ (Lemma 5.2), where DN is the matrix of π^{≤N}DF(ω,W)π^{≤N} with V instead of 5W^4. The product is computed by blocks, for the rows and for the columns : BNJ contains the columns J of B_N, converted to intervals, and s the sums of the columns
#### GC.gc() frees the memory used by each product
Z1_0 = interval(0)
for J = blk
    global Z1_0
    BNJ = interval.(BN[:,J])
    s = interval.(zeros(length(J)))
    for R = blk
        B = DN_rows(R)*BNJ
        if R == J
            B = B - interval.(1.0*I[1:length(J), 1:length(J)])
        end
        s += vec(sum(abs.(w[R].*B), dims=1))
        GC.gc()
    end
    Z1_0 = maximum([Z1_0 ; s./w[J]])
    println("    columns $(J[1]) to $(J[end]) of $(d+1) done   ($(round((time()-t_start)/60, digits=1)) min)")
end
println("Bound Z_{1,0} : Z1_0 = $Z1_0")

#### 𝒯 contains the norms of the columns of π^{>N}(V * π^{≤N}), that is the part of V*h which is not in π^{≤N} : Chebyshev modes n2 > N2 if the Fourier mode n = m+k is in π^{≤N}, all the Chebyshev modes otherwise
wN = weight_nu(ν,N2) ; wNV = weight_nu(ν,N2+NV2)
𝒯 = interval.(zeros(d+1))
for k = filter(iseven, -NV1:NV1)
    Mk = coefficients(MV[k+NV1+1])
    c_out = norm_cols(Mk[d2+1:end,:], wNV[d2+1:end], wN) ; c_all = norm_cols(Mk, wNV, wN)
    for m = -N1:2:N1
        𝒯[ind(m)] += νθ^interval(abs(m+k)-abs(m))*(abs(m+k) <= N1 ? c_out : c_all)
    end
end

#### norms of A (5.28) : norm_πW_A ≥ ‖π_W A‖, norm_πω_A ≥ ‖π_ω A‖, norm_∂ωΛπW_A ≥ ‖∂ωΛ π_W A‖, norm_𝒩²πW_A ≥ ‖𝒩² π_W A‖, and Z1_1 ≥ ‖π^{>N} M_V π_W B_N‖ + ‖5W^4 - V‖ ‖π_W B_N‖ (Lemma 5.2), where the first term is bounded by ‖𝒯 π_W B_N‖. We start with the columns of B_N, which are converted to intervals
wW = [interval(0) ; w[2:end]]   ### weights for the component W
wq = wW.*[interval(0) ; vec([interval(n)^2 for k = 0:N2, n = -N1:2:N1])]   ### weights for 𝒩² π_W
BN1 = interval.(BN[:,1])
norm_πW_A = sum(abs.(wW.*BN1)) ; norm_πω_A = abs(BN1[1]) ; norm_𝒩²πW_A = sum(abs.(wq.*BN1)) ; Z1_1 = sum(abs.(w.*𝒯.*BN1))
for n = -N1:2:N1
    global norm_πW_A, norm_πω_A, norm_𝒩²πW_A, Z1_1
    BNn = interval.(BN[:,ind(n)]) ; wn = w[ind(n)]
    norm_πW_A = maximum([norm_πW_A ; norm_cols(BNn, wW, wn)])
    norm_πω_A = maximum([norm_πω_A ; abs.(BNn[1,:])./wn])
    norm_𝒩²πW_A = maximum([norm_𝒩²πW_A ; norm_cols(BNn, wq, wn)])
    Z1_1 = maximum([Z1_1 ; norm_cols(BNn, w.*𝒯, wn)])
end
#### ∂ωΛ π_W B_N is computed for each Fourier mode m of the rows, sp contains the sums of the columns
sp = interval.(zeros(d+1))
for m = -N1:2:N1
    global sp
    sp += vec(sum(abs.(w[ind(m)].*(dΛN(m)*interval.(BN[ind(m),:]))), dims=1))
    GC.gc()
end
norm_∂ωΛπW_A = maximum(sp./w)
#### the last term of Z1_1 is ‖5W^4 - V‖ ‖π_W B_N‖
Z1_1 = Z1_1 + norm_πW_A*norm_5W4_V
println("Bound Z_{1,1} : Z1_1 = $Z1_1")
#### columns of A which are not in π^{≤N} : A = Λ^{-1} on these columns
println("‖π_W B_N‖ ≤ $norm_πW_A ,  ‖π_ω B_N‖ ≤ $norm_πω_A ,  ‖∂ωΛ π_W B_N‖ ≤ $norm_∂ωΛπW_A ,  ‖𝒩² π_W B_N‖ ≤ $norm_𝒩²πW_A")
norm_πW_A = maximum([norm_πW_A κ]) ; norm_∂ωΛπW_A = maximum([norm_∂ωΛπW_A κ′]) ; norm_𝒩²πW_A = maximum([norm_𝒩²πW_A κ′′])

#### Z1_2 = (‖V‖ + ‖5W^4 - V‖ + 1/(2νθ)) ‖Λ^{-1}π^{>N}‖ (Lemma 5.2), with ‖V‖ + ‖5W^4 - V‖ = ‖5W^4‖ and ‖Λ^{-1}π^{>N}‖ ≤ κ
Z1_2 = (norm_5W4 + I"0.5"/νθ)*κ
println("Bound Z_{1,2} : Z1_2 = $Z1_2")

#### Z1 = max{Z_{1,0} + Z_{1,1}, Z_{1,2}} (Lemma 5.2)
Z1 = maximum([Z1_0+Z1_1 Z1_2])
println("Bound Z1 : Z1 = $Z1")

############ Section 5.4 : computation of the Z2 bound (Lemma 5.8) ############

### Z2 is the bound Z_2(r) of Lemma 5.8 for r = r_Z2 = 1e-6 : it is valid for all r ≤ r_Z2, and r0 < r_Z2 is verified a posteriori.
r_Z2 = interval(1)/interval(1000000)

Z2 = interval(20)*norm_πW_A^2*(norm_W + norm_πW_A*r_Z2)^3 + interval(2)*norm_πω_A*norm_∂ωΛπW_A + norm_πω_A^2*(interval(2)*norm_𝒩²W + interval(3)*norm_𝒩²πW_A*r_Z2)
println("Bound Z2 : Z2 = $Z2")

############ Section 5.5 : radii polynomial conditions (4.15) and conclusion (Theorem 5.12) ############

if (sup(Z1) < 1)&&(inf((interval(1)-Z1)^2) > sup(interval(2)*Z2*Y))

    r0 = (1 - sup(big(Z1)) - sqrt( (1-sup(big(Z1)))^2 - 2*sup(big(Z2))*sup(Y) ) )/inf(Z2)
    r0 = interval(r0)

    p1 = I"0.5"*Z2*r0^2 - (I"1" - Z1)*r0 + Y
    p2 = Z1 + Z2*r0

    #### the last two conditions are the assumptions (5.12), (5.13), and ε < 1 for the approximate inverses of Λ_n (Lemmas 5.7 and 5.9)
    if (sup(p1) < 0)&&(sup(p2) < 1)&&(sup(r0)<inf(r_Z2))&&H&&(sup(ε_max) < 1)
        sup_r = Float64(sup(r0),RoundUp)
        println("-"^74)
        println("  ✓ EXISTENCE PROVEN : a real zero (ω*,W*) of F exists near the approximate solution.")
        println("    Newton–Kantorovich bounds :  Y  = $(sup(Y))")
        println("                                 Z₁ = $(sup(Z1))")
        println("                                 Z₂ = $(sup(Z2))")
        println("    Validated radius          :  r₀ = $sup_r")
        println("    ⇒ |ω* - ω| ≤ $(Float64(sup(norm_πω_A*r0),RoundUp))  and  ‖W* - W‖ ≤ $(Float64(sup(norm_πW_A*r0),RoundUp)).")
    else
        println("  ✗ second condition (p₁ < 0, p₂ < 1, r₀ < r_Z2, (5.12), (5.13) and ε < 1) not verified — proof inconclusive.")
    end
else
    println("  ✗ first condition (Z₁ < 1 and (1-Z₁)² > 2 Z₂ Y) not verified — proof inconclusive.")
end

println("Elapsed time : $(round((time()-t_start)/60, digits=1)) minutes")




########## A posteriori computations on the norm of W (blow-up rate) ##############
M = 2000
####  w[j,k+1] is the coefficient of exp(i(2j-1)θ) T_k(2ρ²-1). RadiiPolynomial uses the convention a_0 + 2∑ a_k T_k, hence w_k = 2a_k for k ≥ 1. The coefficients for the negative Fourier modes are the complex conjugates
w = [interval(Float64(inf(real(c)),RoundDown),Float64(sup(real(c)),RoundUp)) + interval(im)*interval(Float64(inf(imag(c)),RoundDown),Float64(sup(imag(c)),RoundUp)) for c = [interval(k == 0 ? 1 : 2)*W[(2*j-1,k)] for j = 1:div(NW1+1,2), k = 0:NW2]]

#### W(θ,ρ) = 2 Re ∑_{n > 0} V_n(ρ) exp(inθ), where V_n(ρ) = ∑_k w_{n,k} T_k(2ρ²-1). Since T_k(1) = 1 and T_k(-1) = (-1)^k, V1 and V0 contain the values of V_n at ρ = 1 and at ρ = 0
V1 = [sum(w[j,:]) for j = 1:div(NW1+1,2)]
V0 = [sum(w[j,1:2:end]) - sum(w[j,2:2:end]) for j = 1:div(NW1+1,2)]

#### S ≥ |W(θ,ρ)| and L ≥ |∂θW(θ,ρ)| for all θ and all ρ ∈ [0,1], since |T_k| ≤ 1 on [-1,1]
S = interval(2)*sum(abs.(w))
L = interval(2)*sum(interval.(1:2:NW1).*vec(sum(abs.(w), dims=2)))

#### values of W(θ,0) and W(θ,1) on the grid θ = iπ/M, i = 0,...,M. Since W(θ+π,ρ) = -W(θ,ρ), it is enough to consider θ ∈ [0,π]
Wθ(V,θ) = interval(2)*sum(real(V[j])*cos(interval(2*j-1)*θ) - imag(V[j])*sin(interval(2*j-1)*θ) for j = 1:div(NW1+1,2))
m = S
for i = 0:M
    global m
    θ = interval(i)*interval(π)/interval(M)
    m = min(m, max(abs(Wθ(V0,θ)), abs(Wθ(V1,θ))))
end
#### every θ ∈ [0,π] is at distance at most π/(2M) of the grid
m = m - L*interval(π)/interval(2*M)

println("Lipschitz constant in θ : L = $L")
δ = norm_πW_A*r0
println("For all θ :  max(|W*(θ,0)|, |W*(θ,1)|) ≥ $(Float64(inf(m - δ),RoundDown))  and  |W*(θ,ρ)| ≤ $(Float64(sup(S + δ),RoundUp)) for all ρ ∈ [0,1]")