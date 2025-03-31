## <a name="xf527981b70fb60f622623f96ff30c1d9863b944"></a>**Advantages of Ultraweak dPG for Heston-SLV**
The ultraweak dPG formulation offers several advantages for the Heston-SLV option pricing model:

1. **Robustness in Degenerate Regions**: The method is stable for the full range of model parameters, including cases where V≈0, where the PDE becomes degenerate and standard methods may fail.
1. **Handling of Convection-Dominated Regimes**: For high values of κ or r−q, where convection dominates diffusion, the dPG method remains stable without requiring additional stabilization terms.
1. **Built-in Adaptivity**: The residual-based error estimator allows for efficient adaptive refinement, focusing computational resources where they are most needed, such as near the strike price or in regions of high volatility gradients.
1. **Natural Upwinding**: The optimal test functions automatically introduce appropriate upwinding in convection-dominated regions, eliminating the need for ad-hoc stabilization parameters.
1. **Flexibility for Low-Regularity Solutions**: The discontinuous trial space can better capture solutions with low regularity or sharp features, which can occur in option pricing problems with non-smooth payoffs or barriers.
1. **Element-Local Computations**: The computation of optimal test functions is local to each element, making the method highly parallelizable and efficient for large-scale problems.
1. **High-Order Accuracy**: The method naturally accommodates high-order approximations, allowing for spectral convergence rates for smooth solutions.
1. **Robust Treatment of Mixed Derivatives**: The mixed derivative term ∂2u∂S∂V, which is characteristic of the Heston-SLV model, is handled naturally through the auxiliary variable σSV without introducing additional complexity.
## <a name="comparison-with-standard-weak-form"></a>**Comparison with Standard Weak Form**
Compared to the standard weak form presented in the original document, the ultraweak dPG formulation differs in several key aspects:

1. **Trial Space**:
   - Standard weak form: Uses H1-conforming elements (continuous across element boundaries)
   - Ultraweak dPG: Uses discontinuous L2 elements (no continuity requirements)
1. **Test Space**:
   - Standard weak form: Uses the same space for test and trial functions
   - Ultraweak dPG: Uses enriched test spaces with higher regularity
1. **Derivatives**:
   - Standard weak form: Retains first derivatives on trial functions (e.g., ∂u∂S)
   - Ultraweak dPG: Shifts all derivatives to test functions, with trial functions appearing without any derivatives
1. **Numerical Fluxes**:
   - Standard weak form: No explicit numerical fluxes (continuity enforced strongly)
   - Ultraweak dPG: Requires careful handling of interface terms through numerical fluxes for all variables
1. **Stability**:
   - Standard weak form: May require additional stabilization terms for convection-dominated problems
   - Ultraweak dPG: Achieves stability through problem-dependent optimal test functions
1. **Error Estimation**:
   - Standard weak form: No built-in error estimator
   - Ultraweak dPG: Natural residual-based error estimator for adaptive refinement
1. **System Size**:
   - Standard weak form: Smaller system with fewer unknowns
   - Ultraweak dPG: Larger system due to auxiliary variables, but can be reduced through static condensation

These differences make the ultraweak dPG formulation particularly well-suited for challenging option pricing problems where the standard weak form might struggle, such as low-volatility regimes, non-smooth payoffs, or exotic options with barriers or early exercise features.
## <a name="implementation-guidelines"></a>**Implementation Guidelines**
For practical implementation of the ultraweak dPG method for the Heston-SLV model, the following steps are recommended:

1. **Mesh Generation**:
   - Create a suitable mesh for the computational domain 0,Smax×0,Vmax
   - Consider non-uniform meshes with refinement near S=K (strike price) and V=0
1. **Basis Functions**:
   - Choose appropriate basis functions for the trial space (e.g., discontinuous polynomials)
   - Define enriched basis functions for the test space (typically one degree higher)
1. **Optimal Test Functions**:
   - Implement the local adjoint problems to compute optimal test functions for each trial basis function
   - Store these functions or their coefficients for assembly
1. **Assembly**:
   - Assemble the global system matrices and load vectors
   - Include all interface terms to enforce weak continuity across elements
1. **Boundary Conditions**:
   - Incorporate boundary conditions through the boundary integrals
   - Pay special attention to the degenerate boundary at V=0
1. **Static Condensation**:
   - Implement element-local static condensation to reduce system size
   - Solve the condensed system for the primary variable u
   - Recover auxiliary variables as needed
1. **Time Stepping**:
   - Choose an appropriate time integration scheme (Backward Euler, Crank-Nicolson, or BDF2)
   - Implement the time loop with suitable initial conditions
1. **Adaptivity**:
   - Compute error estimators for each element
   - Implement a marking strategy and mesh refinement algorithm
   - Iterate until desired accuracy is achieved
1. **Validation**:
   - Verify the implementation against benchmark cases with known analytical solutions
   - Compare with the standard weak form for simple test cases
1. **Optimization**:
   - Optimize performance through parallelization of local solves
   - Consider matrix-free implementations for large-scale problems

Following these guidelines should result in a robust and efficient implementation of the ultraweak dPG method for solving the Heston-SLV option pricing model. ### Adaptivity and Error Estimation

A significant advantage of the dPG method is its built-in error estimator, which can be used to guide adaptive refinement of the mesh. The error estimator for each element K is defined as:

ηK=∥rK∥VK′

where rK is the element residual and VK′ is the dual of the test space.
#### <a name="x1d695b5633b552d9415516b887436fb72fb6c88"></a>*Practical Computation of the Error Estimator*
To compute this error estimator in practice, we follow these steps:

1. Compute the residual rK for element K:

   For the trial solution uh,σS,h,σV,h,σSS,h,σSV,h,σVV,h, the residual functional rK is defined as:

   ⟨rK,v⟩=bKuh,v−lKv ∀v∈VK

   where bK⋅,⋅ is the bilinear form restricted to element K, and lK⋅ is the load functional.
1. Represent the dual norm ∥rK∥VK′ using the Riesz representation theorem:

   Let eK∈VK be the solution to:

   eK,vVK=⟨rK,v⟩ ∀v∈VK

   Then:

   ∥rK∥VK′=∥eK∥VK
1. In practice, solve for eK using a discrete basis for VK (typically enriched beyond the test space):

   eK=j​cjφj

   where {φj} is a basis for VK, leading to the linear system:

   Gc=r

   with Gij=φi,φjVK and ri=⟨rK,φi⟩.
1. The error estimator is then computed as:

   ηK=cTGc
#### <a name="adaptive-refinement-strategy"></a>*Adaptive Refinement Strategy*
Using the error estimator, we can implement an adaptive refinement strategy:

1. Solve the dPG problem on the current mesh.
1. Compute the error estimator ηK for each element K.
1. Mark elements for refinement using a criterion such as:
   - Maximum strategy: Refine elements with ηK>αmaxJηJ (typically α=0.5).
   - Bulk/Dörfler marking: Refine a minimal subset of elements whose combined error exceeds a fraction of the total error.
1. Refine the marked elements and repeat until a desired accuracy is achieved.

This adaptive process allows the method to automatically concentrate computational resources in regions with high error, such as near singularities or sharp features in the solution.### Time Discretization

For the time-dependent Heston-SLV PDE, the ultraweak formulation can be combined with various time-stepping schemes. The time derivative term ∂u∂t appears in the weak form as:

Ω​∂u∂tv dS dV

We can apply different discretization strategies:

1. Backward Euler (fully implicit, first-order accurate):

   un+1−unΔt+Lun+1=0

   This gives: Ω​un+1−unΔtv dS dV+Ω​Lun+1v dS dV=0

   Where Lun+1 represents all spatial terms evaluated at the new time level, including: r−qSσSn+1+κθ−VσVn+1+12L2VS2σSSn+1+ρξLVSσSVn+1+12ξ2VσVVn+1−run+1
1. Crank-Nicolson (second-order accurate):

   un+1−unΔt+12Lun+1+Lun=0

   This gives: Ω​un+1−unΔtv dS dV+12Ω​Lun+1+Lunv dS dV=0

   With spatial terms evaluated at both time levels.
1. BDF2 (second-order accurate, better stability):

   3un+1−4un+un−12Δt+Lun+1=0

   This gives: Ω​3un+1−4un+un−12Δtv dS dV+Ω​Lun+1v dS dV=0

Each time step requires solving the full ultraweak system for the new solution values un+1,σSn+1,σVn+1,σSSn+1,σSVn+1,σVVn+1.# Ultraweak Discontinuous Petrov-Galerkin Form of Heston-SLV PDE
## <a name="x6c6a93db70259fc209b4b914f29b78f39caf457"></a>**Review of Heston-SLV Model and Strong Form**
The Heston-SLV model combines stochastic volatility with local volatility, governed by the following SDEs:

dSt=r−qStdt+LSt,tVtStdWt1

dVt=κθ−Vtdt+ξVtdWt2

where EdWt1dWt2=ρdt.

The standard PDE for option price uS,V,t is:

∂u∂t+r−qS∂u∂S+κθ−V∂u∂V+12L2S,tVS2∂2u∂S2+ρξLS,tVS∂2u∂S∂V+12ξ2V∂2u∂V2−ru=0
## <a name="first-order-system-reformulation"></a>**First-Order System Reformulation**
The ultraweak formulation begins by recasting the second-order PDE as a first-order system through the introduction of auxiliary variables:

σS=∂u∂S

σV=∂u∂V

σSS=∂σS∂S=∂2u∂S2

σSV=∂σS∂V=∂2u∂S∂V

σVV=∂σV∂V=∂2u∂V2

This gives us the following system:

1. Primary PDE:

∂u∂t+r−qSσS+κθ−VσV+12L2S,tVS2σSS+ρξLS,tVSσSV+12ξ2VσVV−ru=0

1. Auxiliary equations:

σS−∂u∂S=0

σV−∂u∂V=0

σSS−∂σS∂S=0

σSV−∂σS∂V=0

σVV−∂σV∂V=0
## <a name="xa8e81182a3473b0fd4202b5d6eb009ec44fc58c"></a>**Function Spaces for Ultraweak Formulation**
For the ultraweak formulation, we define the function spaces with appropriate regularity requirements:
### <a name="trial-spaces"></a>**Trial Spaces**
All trial functions are in L2Ω (square-integrable functions, no continuity requirements): - u∈L2Ω: The primary variable (option price) - σS,σV∈L2Ω: First derivatives - σSS,σSV,σVV∈L2Ω: Second derivatives
### <a name="test-spaces"></a>**Test Spaces**
Test functions require different regularities depending on their role: - v∈H2Ω: Test function for primary PDE (requires second derivatives) - wS,wV∈H1Ω: Test functions for first derivative equations (require first derivatives) - wSS,wSV,wVV∈H1Ω: Test functions for second derivative equations

The higher regularity of test functions is essential to absorb the derivatives transferred from trial functions during integration by parts. In a practical dPG implementation, these test spaces are typically enriched (using higher polynomial degree) to ensure inf-sup stability.

The domain is Ω=0,Smax×0,Vmax.
## <a name="ultraweak-variational-formulation"></a>**Ultraweak Variational Formulation**
### <a name="multiplication-by-test-functions"></a>**Multiplication by Test Functions**
We multiply each equation by appropriate test functions and integrate over Ω:

1. Primary PDE with test function v:

Ω​∂u∂t+r−qSσS+κθ−VσV+12L2VS2σSS+ρξLVSσSV+12ξ2VσVV−ruv dS dV=0

1. Auxiliary equations with test functions wS,wV,wSS,wSV,wVV:

Ω​σS−∂u∂SwS dS dV=0

Ω​σV−∂u∂VwV dS dV=0

Ω​σSS−∂σS∂SwSS dS dV=0

Ω​σSV−∂σS∂VwSV dS dV=0

Ω​σVV−∂σV∂VwVV dS dV=0
### <a name="integration-by-parts-for-ultraweak-form"></a>**Integration by Parts for Ultraweak Form**
The key distinction of the ultraweak formulation is that we transfer all derivatives from trial functions to test functions through integration by parts.
#### <a name="primary-pde-term"></a>*Primary PDE Term*
The time derivative term remains unchanged:

Ω​∂u∂tv dS dV

For the remaining terms in the primary PDE, no integration by parts is needed as they contain no derivatives of trial functions.
#### <a name="auxiliary-equations"></a>*Auxiliary Equations*
For each auxiliary equation, we apply integration by parts to transfer derivatives from trial to test functions:

1. S-derivative equation:

Ω​σS−∂u∂SwS dS dV=0

Apply integration by parts to the term with ∂u∂S:

Ω​∂u∂SwS dS dV=−Ω​u∂wS∂S dS dV+∂Ω​uwSnS dΓ

Where nS is the S-component of the outward normal vector to the boundary ∂Ω.

The equation becomes:

Ω​σSwS dS dV+Ω​u∂wS∂S dS dV−∂Ω​uwSnS dΓ=0

1. V-derivative equation:

Ω​σV−∂u∂VwV dS dV=0

Apply integration by parts:

Ω​∂u∂VwV dS dV=−Ω​u∂wV∂V dS dV+∂Ω​uwVnV dΓ

The equation becomes:

Ω​σVwV dS dV+Ω​u∂wV∂V dS dV−∂Ω​uwVnV dΓ=0

1. SS-derivative equation:

Ω​σSS−∂σS∂SwSS dS dV=0

Apply integration by parts:

Ω​∂σS∂SwSS dS dV=−Ω​σS∂wSS∂S dS dV+∂Ω​σSwSSnS dΓ

The equation becomes:

Ω​σSSwSS dS dV+Ω​σS∂wSS∂S dS dV−∂Ω​σSwSSnS dΓ=0

1. SV-derivative equation:

Ω​σSV−∂σS∂VwSV dS dV=0

Apply integration by parts:

Ω​∂σS∂VwSV dS dV=−Ω​σS∂wSV∂V dS dV+∂Ω​σSwSVnV dΓ

The equation becomes:

Ω​σSVwSV dS dV+Ω​σS∂wSV∂V dS dV−∂Ω​σSwSVnV dΓ=0

1. VV-derivative equation:

Ω​σVV−∂σV∂VwVV dS dV=0

Apply integration by parts:

Ω​∂σV∂VwVV dS dV=−Ω​σV∂wVV∂V dS dV+∂Ω​σVwVVnV dΓ

The equation becomes:

Ω​σVVwVV dS dV+Ω​σV∂wVV∂V dS dV−∂Ω​σVwVVnV dΓ=0
### <a name="complete-ultraweak-formulation"></a>**Complete Ultraweak Formulation**
Combining all equations, the complete ultraweak formulation is:

Ω​∂u∂tv dS dV+Ω​r−qSσSv dS dV+Ω​κθ−VσVv dS dV+Ω​12L2VS2σSSv dS dV+Ω​ρξLVSσSVv dS dV+Ω​12ξ2VσVVv dS dV−Ω​ruv dS dV=0

Ω​σSwS dS dV+Ω​u∂wS∂S dS dV−∂Ω​uwSnS dΓ=0

Ω​σVwV dS dV+Ω​u∂wV∂V dS dV−∂Ω​uwVnV dΓ=0

Ω​σSSwSS dS dV+Ω​σS∂wSS∂S dS dV−∂Ω​σSwSSnS dΓ=0

Ω​σSVwSV dS dV+Ω​σS∂wSV∂V dS dV−∂Ω​σSwSVnV dΓ=0

Ω​σVVwVV dS dV+Ω​σV∂wVV∂V dS dV−∂Ω​σVwVVnV dΓ=0
## <a name="xd78d05f557c68bbd0f7dd75c58e013a8f771632"></a>**Discontinuous Aspects and Inter-Element Coupling**
For a domain discretized into elements (mesh) Th={K}, we need to account for discontinuities in the trial functions across element boundaries. Let Eh be the set of all interior facets (edges in 2D) of the mesh.
### <a name="element-local-formulation"></a>**Element-Local Formulation**
In the discontinuous Petrov-Galerkin method, the trial functions are allowed to be completely discontinuous across element boundaries. The test functions, however, are chosen from a conforming space with appropriate regularity for each variable.

For each element K∈Th, we have a local contribution to the global system:

Primary PDE (on element K):

K​∂u∂tv dS dV+K​r−qSσSv dS dV+K​κθ−VσVv dS dV+K​12L2VS2σSSv dS dV+K​ρξLVSσSVv dS dV+K​12ξ2VσVVv dS dV−K​ruv dS dV=0

Auxiliary equations (on element K):

K​σSwS dS dV+K​u∂wS∂S dS dV−∂K​uwSnS ds=0

K​σVwV dS dV+K​u∂wV∂V dS dV−∂K​uwVnV ds=0

K​σSSwSS dS dV+K​σS∂wSS∂S dS dV−∂K​σSwSSnS ds=0

K​σSVwSV dS dV+K​σS∂wSV∂V dS dV−∂K​σSwSVnV ds=0

K​σVVwVV dS dV+K​σV∂wVV∂V dS dV−∂K​σVwVVnV ds=0
### <a name="inter-element-coupling"></a>**Inter-Element Coupling**
To account for discontinuities across element boundaries, we introduce numerical fluxes and traces. For any function ϕ and interior edge e∈Eh shared by elements K+ and K−, we define:

- Jump: ​ϕ​=ϕ+n++ϕ−n−, where n+ and n− are the outward normal vectors from K+ and K−.
- Average: {​{ϕ}​}=12ϕ++ϕ−

For each interior edge e∈Eh, we add the following terms to enforce weak continuity for all variables:

1. Primary variable u:

e​​u​⋅{​{wSnS+wVnV}​} ds+e​{​{u}​}⋅​wSnS+wVnV​ ds

1. First derivatives σS and σV:

e​​σS​⋅{​{wSSnS+wSVnV}​} ds+e​{​{σS}​}⋅​wSSnS+wSVnV​ ds

e​​σV​⋅{​{wVVnV}​} ds+e​{​{σV}​}⋅​wVVnV​ ds

1. Second derivatives σSS, σSV, and σVV:

e​​σSS​⋅{​{wSS}​} ds+e​{​{σSS}​}⋅​wSS​ ds

e​​σSV​⋅{​{wSV}​} ds+e​{​{σSV}​}⋅​wSV​ ds

e​​σVV​⋅{​{wVV}​} ds+e​{​{σVV}​}⋅​wVV​ ds

These flux terms are critical for stability, especially for the terms involving mixed derivatives σSV which are characteristic of the Heston-SLV model.
## <a name="optimal-test-functions"></a>**Optimal Test Functions**
A defining feature of the DPG method is the computation of problem-dependent optimal test functions that guarantee stability. In the operator form, for a given trial space basis function ϕi, the corresponding optimal test function ψi is found by solving:

ψi,μV=bϕi,μ=B\*μ,ϕiU ∀μ∈V

where: - B is the PDE operator in ultraweak form - B\* is the adjoint operator - V is the test space (enriched space with higher polynomial degree) - U is the trial space - b⋅,⋅ is the bilinear form of the ultraweak formulation
### <a name="x2993add3c2460e8cb5ab501b6ab1bc7205e50f3"></a>**Rigorous Formulation of the Adjoint Problem**
For each trial basis function ϕi, the optimal test function ψi solves:

ψi,μV=bϕi,μ ∀μ∈VK

where the inner product ⋅,⋅V is typically defined as:

ψ,μV=ψ,μL2K+∇ψ,∇μL2K+higher-order terms as needed

This inner product must be sufficiently strong to control all derivatives appearing in the test space.
### <a name="x195198e848adda3e66969b703e40de41c3810a2"></a>**Element-Local Computation of Optimal Test Functions**
For each element K∈Th, we compute the optimal test functions by solving local problems. Let’s express this more precisely for each trial variable:

1. For a basis function ϕui associated with u, we find ψui=vui,wS,ui,wV,ui,wSS,ui,wSV,ui,wVV,ui by solving:

ψui,μV=K​−rϕuiμ+ϕui∂μS∂S+ϕui∂μV∂V dS dV+boundary terms ∀μ∈VK

1. For a basis function ϕσSi associated with σS, we find ψσSi by solving:

ψσSi,μV=K​r−qSϕσSiμ+ϕσSiμS+ϕσSi∂μSS∂S+ϕσSi∂μSV∂V dS dV+boundary terms ∀μ∈VK

1. For a basis function ϕσVi associated with σV, we find ψσVi by solving:

ψσVi,μV=K​κθ−VϕσViμ+ϕσViμV+ϕσVi∂μVV∂V dS dV+boundary terms ∀μ∈VK

1. For a basis function ϕσSSi associated with σSS, we find ψσSSi by solving:

ψσSSi,μV=K​12L2VS2ϕσSSiμ+ϕσSSiμSS dS dV+boundary terms ∀μ∈VK

1. For a basis function ϕσSVi associated with σSV, we find ψσSVi by solving:

ψσSVi,μV=K​ρξLVSϕσSViμ+ϕσSViμSV dS dV+boundary terms ∀μ∈VK

1. For a basis function ϕσVVi associated with σVV, we find ψσVVi by solving:

ψσVVi,μV=K​12ξ2VϕσVViμ+ϕσVViμVV dS dV+boundary terms ∀μ∈VK
### <a name="enriched-test-space-implementation"></a>**Enriched Test Space Implementation**
In practice, we implement the optimal test functions using an enriched polynomial space: - If the trial space uses polynomials of degree p, then the test space typically uses polynomials of degree p+Δp where Δp≥1. - This enrichment ensures the inf-sup stability condition is satisfied.

The local problems are solved on each element independently, making the approach highly parallelizable.
## <a name="boundary-conditions"></a>**Boundary Conditions**
At the domain boundaries S=0, S=Smax, V=0, and V=Vmax, appropriate boundary conditions must be imposed. These typically include:

1. At S=0: Option price u=0 for call options, or specified values for other option types.
1. At S=Smax: Far-field conditions, e.g., u≈S−Ke−rT−t for call options.
1. At V=0: Feller condition, which is automatically satisfied by the scheme.
1. At V=Vmax: Far-field conditions, often Neumann-type ∂u∂V=0.
### <a name="x12ad2695b09b1550935375976db2a2029476fb1"></a>**Imposing Boundary Conditions in Ultraweak Form**
In the ultraweak formulation, boundary conditions are incorporated in a specific manner:
#### <a name="dirichlet-boundary-conditions"></a>*Dirichlet Boundary Conditions*
For Dirichlet conditions (e.g., at S=0 where u=0 for call options):

1. Replace u with the prescribed value gD in boundary integrals:

∂ΩD​gDwSnS dΓ

Instead of:

∂ΩD​uwSnS dΓ

1. Similarly, for derivative variables at boundaries where they are known:

∂ΩD​gσSwSSnS dΓ

Where gσS is the known value of σS on the Dirichlet boundary.
#### <a name="neumann-boundary-conditions"></a>*Neumann Boundary Conditions*
For Neumann conditions (e.g., at V=Vmax where ∂u∂V=0):

1. For the σV variable (which represents ∂u∂V):

∂ΩN​0⋅wVVnV dΓ=0

1. The boundary term in the equation for u becomes:

∂ΩN​uwVnV dΓ

This term remains in the formulation as an unknown.
### <a name="implementation-in-weak-form"></a>**Implementation in Weak Form**
The complete set of boundary integrals in our formulation includes:

∂Ω​uwSnS dΓ,∂Ω​uwVnV dΓ,∂Ω​σSwSSnS dΓ,∂Ω​σSwSVnV dΓ,∂Ω​σVwVVnV dΓ

For specific boundaries:

1. At S=0 (Dirichlet for call option):

ΓS=0​0⋅wSnS dΓ

1. At S=Smax (Far-field condition):

ΓS=Smax​S−Ke−rT−twSnS dΓ

1. At V=Vmax (Neumann condition):

ΓV=Vmax​σVwVVnV dΓ=0

This explicit treatment of boundary conditions ensures that the ultraweak formulation properly enforces the physical constraints of the Heston-SLV model.
## <a name="matrix-structure"></a>**Matrix Structure**
The ultraweak dPG formulation leads to a block system of the form:

AuuAuσSAuσVAuσSSAuσSVAuσVVAσSuAσSσS0000AσVu0AσVσV0000AσSSσS0AσSSσSS000AσSVσS00AσSVσSV000AσVVσV00AσVVσVVuσSσVσSSσSVσVV=fufσSfσVfσSSfσSVfσVV

where each block matrix is derived from the corresponding bilinear form. For example:

Auuij=−Ω​rϕiuψju dS dV+jump terms

AuσSij=Ω​r−qSϕiσSψju dS dV

AσSuij=Ω​ϕiu∂ψjσS∂S dS dV−boundary terms

and so on, where ϕiu and ψju are the trial and optimal test basis functions, respectively.

The zero blocks in the matrix reflect the decoupling between certain variables in the auxiliary equations. For instance, σS is only directly coupled with u, σSS, and σSV, but not with σV or σVV.
### <a name="static-condensation"></a>**Static Condensation**
To improve computational efficiency, static condensation can be applied by eliminating the auxiliary variables (σS, σV, σSS, σSV, σVV) at the element level, resulting in a reduced system involving only the primary variable u.

The static condensation process works as follows:

1. We start with the block system: AuuAuσAσuAσσuσ=fufσ

   Where σ represents the collective auxiliary variables.
1. From the second row, we express σ in terms of u: σ=Aσσ−1fσ−Aσuu
1. Substitute back into the first row: Auu−AuσAσσ−1Aσuu=fu−AuσAσσ−1fσ
1. This gives us the condensed system: Au=f

   Where: A=Auu−AuσAσσ−1Aσu f=fu−AuσAσσ−1fσ
1. After solving for u, we can recover the auxiliary variables: σ=Aσσ−1fσ−Aσuu

The key advantage is that the condensed system for u is much smaller than the full system, particularly for the Heston-SLV model where we have five auxiliary variables. This significantly reduces the computational cost while maintaining the accuracy of the solution.

For the Heston-SLV problem, the static condensation can be performed element-by-element before assembly, taking advantage of the local nature of the DPG formulation.
