"""property_docs.py -- the IQARIS property dictionary: one machine-readable entry per
database property (name, symbol, unit, type, method, prose description and formal equation).

This is the SINGLE SOURCE OF TRUTH for "what does this column mean?".  It is a faithful,
machine-readable mirror of the manuscript's property dictionary (the paper's Supporting
Information: the complete property dictionary and its master property table), keyed to the
ACTUAL public database column names.

It feeds four consumers, all from this one file:

  1. ``registry.ColumnSpec.description`` -- the per-column blurb.
  2. the CLI / API ``info`` command  (``ie.info("q_QTAIM")`` /
     ``python -m iqaris_export info q_QTAIM``) -- prints name, symbol, type, method, unit,
     description and the equation in plain-text / unicode.
  3. the interactive viewer popup (``viz3d_qt``) -- the same information, with the equation
     typeset by matplotlib and shown next to the property selector.
  4. ``glossary()`` -- a rendered Markdown / HTML reference document of every property.

Each entry carries the equation twice: ``eq_text`` (a unicode one-liner for the terminal and
as a rendering fallback) and ``eq_latex`` (a single-line, matplotlib-mathtext-safe LaTeX
string for the typeset popup / a copy-paste source).  Units for a *specific* column are taken
from the live :mod:`registry` when a database connection is available (they are exact there);
the ``unit`` field below is the family-level unit used by the standalone glossary.

Component families (vectors, tensors: dipole X/Y/Z, Hessian/stress eigenvalues 1/2/3, ...) are
stored as indexed component columns in the database, so one entry's ``columns`` lists the whole
family (exact names and/or ``fnmatch`` patterns).  Looking up any real column returns its entry.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict
from fnmatch import fnmatchcase

# Level-of-theory labels (shared). QTAIM density topology is computed at both M06-2X and PBE0
# (a *_QTAIM column is queryable at --level m062x and --level pbe0); ab initio IQA is PBE0 only.
M_DFT = "M06-2X/def2-TZVP and PBE0/def2-TZVP (QTAIM density topology, gas phase)"
M_SQM = "PM7 / PM6 / PM6-D3H4 / PM6-D3H4X / PM6-ORG (IQA-SQM)"
M_IQADFT = "PBE0/def2-TZVP (ab initio IQA, gas phase)"


@dataclass(frozen=True)
class PropDoc:
    key: str                     # canonical public key (the manuscript "Key")
    name: str                    # human-readable property name
    symbol: str = ""             # unicode symbol, e.g. "q(A)", "delta(A,B)"
    unit: str = ""               # family-level unit (glossary); live unit comes from registry
    type: str = ""               # Molecular | Atomic | Interatomic | CP | BCP
    method: str = M_DFT          # level(s) of theory
    text: str = ""               # prose description
    eq_text: str = ""            # unicode / plain-text equation (terminal + fallback)
    eq_latex: str = ""           # single-line, mathtext-safe LaTeX (typeset popup)
    bonus: bool = False          # a "bonus" property (stored + documented, not central)
    columns: tuple = ()          # real DB column names / fnmatch patterns this entry covers
    see_also: tuple = ()         # related keys

    def to_dict(self):
        d = asdict(self)
        return d


#  The dictionary.  Organized by family, mirroring the Supporting Information.
DOCS: list = [

    #  QTAIM (M06-2X, PBE0) -- molecular energy contributions
    PropDoc(
        key="E_SCF", name="SCF total energy (IQA reference)", symbol="E_SCF", unit="Eh",
        type="Molecular", method=M_IQADFT,
        text="Self-consistent-field total electronic energy of the PBE0/def2-TZVP wavefunction -- the "
             "reference energy the ab initio IQA decomposition reconstructs (closure residual "
             "DE_IQA_WFN = E_IQA_MOL - E_SCF). Lives in iqa_pbe0_molecular. NOT provided at M06-2X: the "
             "additive total-energy decomposition is IQA's job (PBE0 + IQA-SQM); M06-2X carries the QTAIM "
             "density topology and the density-derived kinetic/potential fields only.",
        eq_text="E = ⟨ Ψ | Ĥ | Ψ ⟩",
        eq_latex=r"E = \langle \Psi | \hat{H} | \Psi \rangle",
        columns=("E_SCF",), see_also=("E_IQA_MOL", "DE_IQA_WFN")),

    #  QTAIM (M06-2X, PBE0) -- atomic energy terms
    PropDoc(
        key="L_QTAIM", name="Atomic Lagrangian", symbol="L(A)", unit="Eh", type="Atomic",
        text="One fourth of the negative integral of the Laplacian of the electron density over "
             "the atomic basin. For a properly integrated QTAIM basin L(A) is essentially zero "
             "(the zero-flux condition), so |L(A)| is a standard integration-quality check.",
        eq_text="L(A) = -¼ ∫_ΩA ∇²ρ(r) dr",
        eq_latex=r"L(A) = -\frac{1}{4}\int_{\Omega_A}\nabla^2\rho(\mathbf{r})\,d\mathbf{r}",
        columns=("L_QTAIM",), see_also=("G_QTAIM", "K_QTAIM")),
    PropDoc(
        key="G_QTAIM", name="Atomic Lagrangian kinetic energy", symbol="G(A)", unit="Eh",
        type="Atomic",
        text="Positive-definite (Lagrangian) form of the kinetic energy density integrated over "
             "the atomic basin, from the mixed second derivatives of the one-particle density matrix.",
        eq_text="G(A) = ½ ∫_ΩA ∇·∇' ρ(r,r')|_{r'→r} dr",
        eq_latex=r"G(A) = \frac{1}{2}\int_{\Omega_A}\nabla\!\cdot\!\nabla'\rho(\mathbf{r},\mathbf{r}')\,d\mathbf{r}",
        columns=("G_QTAIM",), see_also=("K_QTAIM", "L_QTAIM")),
    PropDoc(
        key="K_QTAIM", name="Atomic Hamiltonian kinetic energy", symbol="K(A)", unit="Eh",
        type="Atomic",
        text="Hamiltonian form of the atomic kinetic energy. Related to G and L by K(A) = G(A) + L(A); "
             "since L(A)≈0 for a well-integrated basin, K(A) ≈ G(A).",
        eq_text="K(A) = G(A) + L(A)   (≈ G(A) since L≈0)",
        eq_latex=r"K(A) = G(A) + L(A)",
        columns=("K_QTAIM",), see_also=("G_QTAIM", "L_QTAIM")),
    PropDoc(
        key="T_QTAIM", name="Atomic kinetic energy", symbol="T(A)", unit="Eh", type="Atomic",
        text="Electronic kinetic energy of atom A in Hamiltonian form, the local kinetic "
             "contribution of basin ΩA.",
        eq_text="T(A) = ½ ∫_ΩA ∇Ψ* · ∇Ψ dr",
        eq_latex=r"T(A) = \frac{1}{2}\int_{\Omega_A}\nabla\Psi^{*}\!\cdot\!\nabla\Psi\,d\mathbf{r}",
        columns=("T_QTAIM",)),
    PropDoc(
        key="EE_QTAIM", name="Electronic energy contribution of atom", symbol="E_e(A)", unit="Eh",
        type="Atomic",
        text="Contribution of atom A to the total electronic energy, E_e(A) = -T(A). Summed over "
             "atoms it reconstructs the molecular energy when the electronic virial theorem holds.",
        eq_text="E_e(A) = -T(A)",
        eq_latex=r"E_e(A) = -T(A)",
        columns=("EE_QTAIM",)),
    PropDoc(
        key="VEN_AA_QTAIM", name="Electron–nucleus attraction (own nucleus)",
        symbol="V_en(A)", unit="Eh", type="Atomic",
        text="Attraction between the electron density of atom A and its OWN nucleus (Z_A at R_A); "
             "the intra-atomic attractive component (AIMALL VenO).",
        eq_text="V_en(A) = -∫_ΩA [Z_A / |r - R_A|] ρ(r) dr",
        eq_latex=r"V_{en}(A) = -\int_{\Omega_A}\frac{Z_A}{|\mathbf{r}-\mathbf{R}_A|}\rho(\mathbf{r})\,d\mathbf{r}",
        columns=("VEN_AA_QTAIM",), see_also=("VEN_AM_QTAIM",)),
    PropDoc(
        key="VEN_AM_QTAIM", name="Electron–nuclei attraction (all nuclei)",
        symbol="V_en^mol(A)", unit="Eh", type="Atomic",
        text="Attraction between the electron density of atom A and ALL nuclei in the molecule "
             "(B runs over every atom).",
        eq_text="V_en^mol(A) = -Σ_B ∫_ΩA [Z_B / |r - R_B|] ρ(r) dr",
        eq_latex=r"V_{en}^{mol}(A) = -\sum_{B}\int_{\Omega_A}\frac{Z_B}{|\mathbf{r}-\mathbf{R}_B|}\rho(\mathbf{r})\,d\mathbf{r}",
        columns=("VEN_AM_QTAIM",), see_also=("VEN_AA_QTAIM",)),
    PropDoc(
        key="VNN_AM_QTAIM", name="Atomic nuclear repulsion (half)", symbol="V_nn(A,Mol)/2",
        unit="Eh", type="Atomic",
        text="Per-atom half-contribution of atom A to the total nuclear–nuclear repulsion "
             "(AIMALL Vnn(A,Mol)/2). Summed over atoms it gives the molecular nuclear repulsion below.",
        eq_text="V_nn^tot = ½ Σ_A Σ_{B≠A} Z_A Z_B / R_AB",
        eq_latex=r"V_{nn}^{tot} = \frac{1}{2}\sum_{A}\sum_{B\neq A}\frac{Z_A Z_B}{R_{AB}}",
        columns=("VNN_AM_QTAIM",)),

    #  QTAIM (M06-2X, PBE0) -- electron metrics
    PropDoc(
        key="q_QTAIM", name="Atomic charge", symbol="q(A)", unit="e", type="Atomic",
        text="Net electronic charge of an atom: its nuclear charge minus the electron density "
             "integrated over the QTAIM atomic basin ΩA.",
        eq_text="q(A) = Z_A - ∫_ΩA ρ(r) dr",
        eq_latex=r"q_A = Z_A - \int_{\Omega_A}\rho(\mathbf{r})\,d\mathbf{r}",
        columns=("q_QTAIM",), see_also=("N_QTAIM",)),
    PropDoc(
        key="N_QTAIM", name="Electron population", symbol="N(A)", unit="e", type="Atomic",
        text="Total electron count of an atom, the integral of ρ over its basin. "
             "Equivalently N_A = Z_A - q_A.",
        eq_text="N(A) = ∫_ΩA ρ(r) dr  =  Z_A - q_A",
        eq_latex=r"N_A = \int_{\Omega_A}\rho(\mathbf{r})\,d\mathbf{r} = Z_A - q_A",
        columns=("N_QTAIM",), see_also=("q_QTAIM",)),
    PropDoc(
        key="LI_QTAIM", name="Localization index", symbol="λ(A)", unit="e", type="Atomic",
        text="Number of electrons localized exclusively within basin ΩA (not shared with "
             "neighbors): the population minus its variance, equivalently minus the integral of the "
             "exchange-correlation pair density over ΩA×ΩA.",
        eq_text="λ_A = N_A - σ²(N_A)",
        eq_latex=r"\lambda_A = N_A - \sigma^2_{N_A}",
        columns=("LI_QTAIM",), see_also=("DI_QTAIM", "N_QTAIM")),
    PropDoc(
        key="MU_INTRA_QTAIM", name="Atomic dipole moment (intra-atomic / polarization)",
        symbol="μ_Intra(A)", unit="e·bohr", type="Atomic",
        text="Asymmetry (polarization) of the electron density within basin ΩA relative to its "
             "nucleus; a 3-vector. Stored as X/Y/Z components (magnitude MU_INTRA_MAG_QTAIM).",
        eq_text="μ_Intra(A) = -∫_ΩA (r - R_A) ρ(r) dr",
        eq_latex=r"\boldsymbol{\mu}_{Intra}(A) = -\int_{\Omega_A}(\mathbf{r}-\mathbf{R}_A)\rho(\mathbf{r})\,d\mathbf{r}",
        columns=("MU_INTRA_X_QTAIM", "MU_INTRA_Y_QTAIM", "MU_INTRA_Z_QTAIM"),
        see_also=("MU_INTER_QTAIM", "MU_QTAIM", "MU_INTRA_MAG_QTAIM")),
    PropDoc(
        key="MU_INTER_QTAIM", name="Atomic dipole moment (inter-atomic / charge-transfer)",
        symbol="μ_Bond(A)", unit="e·bohr", type="Atomic",
        text="Charge-transfer contribution to the atomic dipole: the effect of neighboring atomic "
             "charges, referenced to the bond critical points. Stored as X/Y/Z components.",
        eq_text="μ_Bond(A) = Σ_{B≠A} q_B [R_BCP(A|B) - R_A]",
        eq_latex=r"\boldsymbol{\mu}_{Bond}(A) = \sum_{B\neq A} q_B\,[\mathbf{R}_{BCP}(A|B) - \mathbf{R}_A]",
        columns=("MU_INTER_X_QTAIM", "MU_INTER_Y_QTAIM", "MU_INTER_Z_QTAIM"),
        see_also=("MU_INTRA_QTAIM", "MU_QTAIM")),
    PropDoc(
        key="MU_QTAIM", name="Atomic dipole moment (total)", symbol="μ(A)", unit="e·bohr",
        type="Atomic",
        text="Total contribution of atom A to the molecular dipole: the sum of its intra-atomic "
             "(polarization) and inter-atomic (charge-transfer) dipoles. Σ_A μ(A) is the "
             "molecular dipole. Stored as X/Y/Z components.",
        eq_text="μ(A) = μ_Intra(A) + μ_Bond(A)",
        eq_latex=r"\boldsymbol{\mu}(A) = \boldsymbol{\mu}_{Intra}(A) + \boldsymbol{\mu}_{Bond}(A)",
        columns=("MU_X_QTAIM", "MU_Y_QTAIM", "MU_Z_QTAIM"),
        see_also=("MU_INTRA_QTAIM", "MU_INTER_QTAIM")),
    PropDoc(
        key="MU4_QTAIM", name="Atomic traceless quadrupole tensor", symbol="Q(A)",
        unit="e·bohr²", type="Atomic",
        text="Second-rank electrostatic moment describing the anisotropy of the electronic charge "
             "of atom A about its nucleus. Stored as the 6 independent components XX, XY, XZ, YY, YZ, ZZ.",
        eq_text="Q(A) = -∫_ΩA ρ(r) [3 r_A r_A - r_A² I] dr",
        eq_latex=r"\mathbf{Q}(A) = -\int_{\Omega_A}\rho(\mathbf{r})\,[3\,\mathbf{r}_A\mathbf{r}_A - r_A^2\mathbf{I}]\,d\mathbf{r}",
        columns=("MU4_XX_QTAIM", "MU4_XY_QTAIM", "MU4_XZ_QTAIM",
                 "MU4_YY_QTAIM", "MU4_YZ_QTAIM", "MU4_ZZ_QTAIM")),

    #  QTAIM (M06-2X, PBE0) -- delocalization index (pair)
    PropDoc(
        key="DI_QTAIM", name="Delocalization index", symbol="δ(A,B)", unit="e (pairs)",
        type="Interatomic",
        text="Number of electron pairs shared between basins ΩA and ΩB: minus twice the "
             "covariance of their populations, equivalently minus twice the exchange-correlation "
             "pair density integrated over ΩA×ΩB. The QTAIM measure of bond order. "
             "With LI it reconstructs the population: N_A = λ_A + ½ Σ_{B≠A} δ(A,B).",
        eq_text="δ(A,B) = -2 σ(N_A,N_B) = -2 ∫∫_{ΩA×ΩB} ρ_xc(r1,r2) dr1 dr2",
        eq_latex=r"\delta_{A,B} = -2\,\sigma_{N_A,N_B} = -2\int_{\Omega_A}\!\int_{\Omega_B}\rho^{xc}_2(\mathbf{r}_1,\mathbf{r}_2)\,d\mathbf{r}_1 d\mathbf{r}_2",
        columns=("DI_QTAIM",), see_also=("LI_QTAIM",)),

    #  QTAIM (M06-2X, PBE0) -- molecular graph / CP counts / Poincare-Hopf
    PropDoc(
        key="N_BCP_QTAIM", name="Critical-point counts (molecular graph)",
        symbol="n_NACP, n_NNACP, n_BCP, n_RCP, n_CCP", unit="count", type="Molecular",
        text="Number of each type of critical point of ρ in the molecular graph: nuclear (NACP), "
             "non-nuclear attractor (NNACP), bond (BCP), ring (RCP) and cage (CCP) critical points.",
        columns=("N_BCP_QTAIM", "N_RCP_QTAIM", "N_CCP_QTAIM", "N_NACP_QTAIM", "N_NNACP_QTAIM"),
        see_also=("PH_QTAIM", "CP_TYPE_QTAIM")),
    PropDoc(
        key="PH_QTAIM", name="Poincaré–Hopf relation", symbol="PH_r", unit="bool",
        type="Molecular",
        text="Topological-consistency flag for the molecular graph via the Euler characteristic of "
             "the density field: PH_r should equal 1 for a finite molecule. A graph satisfying it is "
             "topologically complete (all gradient paths and basins consistently connected).",
        eq_text="PH_r = n_NACP - n_BCP + n_RCP - n_CCP   (= 1)",
        eq_latex=r"PH_r = n_{NACP} - n_{BCP} + n_{RCP} - n_{CCP}",
        columns=("PH_QTAIM",), see_also=("N_BCP_QTAIM",)),

    #  QTAIM (M06-2X, PBE0) -- critical-point fields
    PropDoc(
        key="CP_POS_QTAIM", name="Critical-point coordinates", symbol="R_CP", unit="Angstrom",
        type="CP",
        text="Position of the critical point in space, locating it within the molecular graph. "
             "Stored as CP_X_ANG / CP_Y_ANG / CP_Z_ANG (converted from bohr to Å).",
        columns=("CP_X_ANG", "CP_Y_ANG", "CP_Z_ANG")),
    PropDoc(
        key="CP_TYPE_QTAIM", name="Critical-point type", symbol="(r,s)", unit="label", type="CP",
        text="Classification of the CP by rank r (nonzero Hessian eigenvalues) and signature s "
             "(sum of eigenvalue signs): (3,-3) nuclear (NCP/NNA), (3,-1) bond (BCP), (3,+1) ring "
             "(RCP), (3,+3) cage (CCP). CP_ATOMS lists the atoms it connects; CP_RANK / CP_SIG give r, s.",
        columns=("CP_TYPE_QTAIM", "CP_RANK_QTAIM", "CP_SIG_QTAIM", "CP_ATOMS_QTAIM")),
    PropDoc(
        key="CP_RHO_QTAIM", name="Electron density at the CP", symbol="ρ_CP", unit="e/bohr³",
        type="CP",
        text="Value of the electron density at the critical point. At a BCP it correlates with the "
             "strength of the interaction within a comparable series of atom pairs.",
        eq_text="ρ_CP = ρ(R_CP)",
        eq_latex=r"\rho_{CP} = \rho(\mathbf{R}_{CP})",
        columns=("CP_RHO_QTAIM",)),
    PropDoc(
        key="CP_GRADRHO_QTAIM", name="Density gradient at the CP", symbol="∇ρ_CP",
        unit="e/bohr⁴", type="CP",
        text="Gradient vector of the electron density at the CP (ideally zero at a true CP; the "
             "stored value gauges how well the point is converged). Stored as X/Y/Z + magnitude.",
        eq_text="∇ρ_CP = ∇ρ(R_CP)",
        eq_latex=r"\nabla\rho_{CP} = \nabla\rho(\mathbf{R}_{CP})",
        columns=("CP_GRADRHO_X_QTAIM", "CP_GRADRHO_Y_QTAIM", "CP_GRADRHO_Z_QTAIM",
                 "CP_GRADRHO_MAG_QTAIM")),
    PropDoc(
        key="CP_DELSQRHO_QTAIM", name="Laplacian of the density at the CP",
        symbol="∇²ρ_CP", unit="e/bohr⁵", type="CP",
        text="Sum of the Hessian eigenvalues of ρ at the CP. Negative values indicate local "
             "charge concentration (shared/covalent interactions), positive values depletion "
             "(closed-shell / ionic interactions).",
        eq_text="∇²ρ_CP = α₁ + α₂ + α₃",
        eq_latex=r"\nabla^2\rho_{CP} = \sum_{i=1}^{3}\alpha_i",
        columns=("CP_DELSQRHO_QTAIM",), see_also=("CP_HESSEIG_QTAIM",)),
    PropDoc(
        key="CP_HESSEIG_QTAIM", name="Hessian eigenvalues of the density at the CP",
        symbol="α₁ ≤ α₂ ≤ α₃", unit="e/bohr⁵", type="CP",
        text="The three eigenvalues (principal curvatures) of the density Hessian at the CP, in "
             "ascending order. Their signs define the CP type; their trace is the Laplacian. "
             "Stored as CP_HESSEIG_1/2/3.",
        columns=("CP_HESSEIG_1_QTAIM", "CP_HESSEIG_2_QTAIM", "CP_HESSEIG_3_QTAIM"),
        see_also=("CP_ELLIP_QTAIM", "CP_DELSQRHO_QTAIM")),
    PropDoc(
        key="CP_ELLIP_QTAIM", name="Bond ellipticity", symbol="ε", unit="dimensionless",
        type="BCP",
        text="Departure of the density along a bond path from cylindrical symmetry, from the two "
             "negative Hessian curvatures at the BCP. ε gauges π-character / bond strain "
             "(0 for a cylindrically symmetric bond).",
        eq_text="ε = α₁/α₂ - 1",
        eq_latex=r"\varepsilon = \frac{\alpha_1}{\alpha_2} - 1",
        columns=("CP_ELLIP_QTAIM",), see_also=("CP_HESSEIG_QTAIM",)),
    PropDoc(
        key="CP_V_QTAIM", name="Virial field at the CP", symbol="V", unit="Eh/bohr³", type="CP",
        text="Potential-energy density at the CP, equal to the trace of the quantum stress tensor; "
             "characterizes the local balance of forces.",
        eq_text="V = Tr(σ)",
        eq_latex=r"V = \mathrm{Tr}(\sigma)",
        columns=("CP_V_QTAIM",)),
    PropDoc(
        key="CP_G_QTAIM", name="Kinetic-energy densities at the CP", symbol="G, K, L",
        unit="Eh/bohr³", type="CP",
        text="Local kinetic-energy densities evaluated at the CP: G (Lagrangian, positive-definite), "
             "K (Hamiltonian), and L = K - G = -¼∇²ρ. The G/ρ and |V|/G ratios "
             "at a BCP diagnose the interaction type.",
        eq_text="L = K - G = -¼ ∇²ρ",
        eq_latex=r"L = K - G = -\tfrac{1}{4}\nabla^2\rho",
        columns=("CP_G_QTAIM", "CP_K_QTAIM", "CP_L_QTAIM")),
    PropDoc(
        key="CP_VNUC_QTAIM", name="Potentials at the CP (nuclear / e-n / repulsive)",
        symbol="V_nuc, V_en, V_rep", unit="Eh/e (V_nuc); Eh/bohr³ (V_en, V_rep)", type="CP",
        text="Potential components at the CP: V_nuc (electrostatic potential from the nuclei), "
             "V_en (electron–nuclear contribution to the virial field) and V_rep = V - V_en "
             "(repulsive remainder).",
        eq_text="V_rep = V - V_en",
        eq_latex=r"V_{rep} = V - V_{en}",
        columns=("CP_VNUC_QTAIM", "CP_VEN_QTAIM", "CP_VREP_QTAIM")),
    PropDoc(
        key="CP_DELSQV_QTAIM", name="Laplacians of the CP potentials / kinetic densities",
        symbol="∇²V, ∇²V_en, ∇²V_rep, ∇²G, ∇²K",
        unit="a.u.", type="CP",
        text="Local curvatures (Laplacians) of the potential-energy and kinetic-energy densities at "
             "the CP.",
        columns=("CP_DELSQV_QTAIM", "CP_DELSQVEN_QTAIM", "CP_DELSQVREP_QTAIM",
                 "CP_DELSQG_QTAIM", "CP_DELSQK_QTAIM")),
    PropDoc(
        key="CP_STRESSEIG_QTAIM", name="Stress-tensor eigenvalues at the CP",
        symbol="σ₁, σ₂, σ₃", unit="a.u.", type="CP",
        text="Eigenvalues of the quantum-mechanical stress tensor σ at the CP (its trace is the "
             "virial field V). Stored as CP_STRESSEIG_1/2/3.",
        columns=("CP_STRESSEIG_1_QTAIM", "CP_STRESSEIG_2_QTAIM", "CP_STRESSEIG_3_QTAIM"),
        see_also=("CP_EHRENFEST_QTAIM", "CP_V_QTAIM")),
    PropDoc(
        key="CP_EHRENFEST_QTAIM", name="Ehrenfest force density at the CP", symbol="-∇·σ",
        unit="a.u.", type="CP",
        text="Ehrenfest force density at the CP, minus the divergence of the stress tensor. "
             "Stored as X/Y/Z components.",
        eq_text="F_Ehrenfest = -∇·σ",
        eq_latex=r"\mathbf{F} = -\nabla\!\cdot\!\sigma",
        columns=("CP_EHRENFEST_X_QTAIM", "CP_EHRENFEST_Y_QTAIM", "CP_EHRENFEST_Z_QTAIM")),
    PropDoc(
        key="CP_ESP_QTAIM", name="Electrostatic potential at the CP", symbol="ESP, ESP_e, ESP_n",
        unit="Eh/e", type="CP",
        text="Molecular electrostatic potential at the CP and its electronic (ESP_e) and nuclear "
             "(ESP_n) components: ESP = ESP_e + ESP_n.",
        eq_text="ESP = ESP_e + ESP_n",
        eq_latex=r"\mathrm{ESP} = \mathrm{ESP}_e + \mathrm{ESP}_n",
        columns=("CP_ESP_QTAIM", "CP_ESPE_QTAIM", "CP_ESPN_QTAIM")),

    #  QTAIM (M06-2X, PBE0) -- bond-path properties (BCP only)
    PropDoc(
        key="BP_QTAIM", name="Bond path (bonded atom pair)", symbol="BP", unit="—",
        type="Interatomic",
        text="The bond path is the gradient line of ρ linking two nuclear attractors through "
             "their BCP; its existence is the QTAIM signature of a bonded interaction. CP_BCP_I / "
             "CP_BCP_J give the indices of the two bonded atoms.",
        columns=("BP_QTAIM", "CP_BCP_I", "CP_BCP_J")),
    PropDoc(
        key="BPL_QTAIM", name="Bond path length", symbol="BPL", unit="Angstrom", type="Interatomic",
        text="Length measured ALONG the (generally curved) bond path between the two attractors. "
             "BPL ≥ the straight-line internuclear distance; the excess measures bond-path curvature.",
        columns=("BPL_QTAIM",), see_also=("GBL_I_QTAIM",)),
    PropDoc(
        key="GBL_I_QTAIM", name="Geometric bond lengths (I–IV)",
        symbol="GBL_I … GBL_IV", unit="Angstrom", type="Interatomic",
        text="Four straight-line bond-length measures for a BCP: I distance between nuclear "
             "attractors; II distance between nuclei; III sum of BCP→attractor distances; "
             "IV sum of BCP→nucleus distances.",
        columns=("GBL_I_QTAIM", "GBL_II_QTAIM", "GBL_III_QTAIM", "GBL_IV_QTAIM"),
        see_also=("BPL_QTAIM",)),

    #  QTAIM (M06-2X, PBE0) -- interatomic-surface (IAS) integrals
    PropDoc(
        key="N_IAS_QTAIM", name="Density integral over the interatomic surface",
        symbol="N_IAS(A|B)", unit="e/bohr", type="Interatomic",
        text="Electron density integrated over the zero-flux interatomic surface separating basins "
             "A and B — the amount of density crossing the A|B interface.",
        eq_text="N_IAS(A|B) = ∫_IAS(A|B) ρ(r) dS",
        eq_latex=r"N_{IAS}(A|B) = \int_{IAS(A|B)}\rho(\mathbf{r})\,dS",
        columns=("N_IAS_QTAIM",)),
    PropDoc(
        key="G_IAS_QTAIM", name="Kinetic / Laplacian / virial integrals over the IAS",
        # Heterogeneous family: G/K/L/V_IAS are Eh/bohr, DELSQRHO_IAS is e/bohr^3.  Use the ';'-joined
        # composite convention so registry._resolve_unit falls back to the per-column _unit_for and
        # `info`/`list-props`/exports show each column's EXACT unit, rather than a single parenthetical
        # form that would lead with the wrong Eh/bohr for DELSQRHO_IAS_QTAIM.
        symbol="G_IAS, K_IAS, L_IAS, ∇²ρ_IAS, V_IAS", unit="Eh/bohr; e/bohr³ (∇²ρ_IAS)",
        type="Interatomic",
        text="Local scalar fields integrated over the interatomic surface separating basins A and B: "
             "the Lagrangian (G) and Hamiltonian (K) kinetic-energy densities, L = K - G, the density "
             "Laplacian, and the virial field V. Together they characterize the A|B interface.",
        eq_text="X_IAS(A|B) = ∫_IAS(A|B) X(r) dS   (X = G, K, L, ∇²ρ, V)",
        eq_latex=r"X_{IAS}(A|B) = \int_{IAS(A|B)} X(\mathbf{r})\,dS",
        columns=("G_IAS_QTAIM", "K_IAS_QTAIM", "L_IAS_QTAIM", "DELSQRHO_IAS_QTAIM", "V_IAS_QTAIM"),
        see_also=("N_IAS_QTAIM",)),

    #  QTAIM (M06-2X, PBE0) -- bond-path / interatomic angles
    PropDoc(
        key="BPA_QTAIM", name="Bond-path and geometric bond angles",
        symbol="BPA, GBA_I–IV", unit="degree", type="Interatomic",
        text="Angle subtended at a central atom by two bond paths (BPA), and the corresponding "
             "straight-line geometric bond angles GBA_I–IV (mirroring the four GBL definitions).",
        columns=("BPA_QTAIM", "GBA_I_QTAIM", "GBA_II_QTAIM", "GBA_III_QTAIM", "GBA_IV_QTAIM")),
    PropDoc(
        key="ELEMENT_IJK_QTAIM", name="Angle-triple element labels",
        symbol="El(I,J,K)", unit="label", type="Interatomic",
        text="Element symbols of the three atoms (I=terminal, J=vertex, K=terminal) that define a "
             "bond angle in the angles family. VARCHAR labels, not a numeric measurement.",
        columns=("ELEMENT_I", "ELEMENT_J", "ELEMENT_K")),

    #  QTAIM bonus properties (stored + documented; not central to the paper)
    PropDoc(
        key="MU_INTRA_MAG_QTAIM", name="Intra-atomic dipole magnitude", symbol="|μ_Intra(A)|",
        unit="e·bohr", type="Atomic", bonus=True,
        text="Magnitude of the intra-atomic (polarization) dipole vector μ_Intra(A).",
        columns=("MU_INTRA_MAG_QTAIM",), see_also=("MU_INTRA_QTAIM",)),
    PropDoc(
        key="LOC_PCT_QTAIM", name="Localization percentage", symbol="%Loc(A)", unit="percent",
        type="Atomic", bonus=True,
        text="Percentage of the atomic population that is localized within the basin, 100·λ_A/N_A.",
        columns=("LOC_PCT_QTAIM",), see_also=("LI_QTAIM",)),
    PropDoc(
        key="DI_ATOM_QTAIM", name="Atomic delocalization (total / bonded / non-bonded)",
        symbol="DI(A)/2", unit="e", type="Atomic", bonus=True,
        text="Half the sum of an atom's delocalization indices — the electrons it shares with the "
             "rest of the molecule — split into bonded and non-bonded partner contributions.",
        columns=("DI_ATOM_QTAIM", "DI_BOND_QTAIM", "DI_NONBOND_QTAIM"), see_also=("DI_QTAIM",)),
    PropDoc(
        key="D2_AA_QTAIM", name="Atomic electron-pair counts", symbol="D2(A,A), D2(A,A'), D2(A,Mol)",
        unit="pairs", type="Atomic", bonus=True,
        text="Numbers of electron pairs: within atom A, between A and the other atoms, and between A "
             "and the whole molecule.",
        columns=("D2_AA_QTAIM", "D2_AAP_QTAIM", "D2_AMOL_QTAIM")),
    PropDoc(
        key="AREA_QTAIM", name="Atomic surface area / volume", symbol="Area(A), Vol(A), NVol(A)",
        unit="bohr² (Area); bohr³ (Vol); e (NVol); e/bohr³ (NVol density); percent (NVol pct)",
        type="Atomic", bonus=True,
        text="Size of the atomic basin bounded by the 0.001 a.u. isodensity envelope: the surface "
             "area, the enclosed volume, and the volume weighted by density (with derived density / "
             "percentage columns).",
        columns=("AREA_QTAIM", "VOL_QTAIM", "NVOL_QTAIM", "NVOL_DENS_QTAIM", "NVOL_PCT_QTAIM")),
    PropDoc(
        key="R_QTAIM", name="Atomic radial moments", symbol="R_{-2..+2}(A)", unit="a.u.",
        type="Atomic", bonus=True,
        text="Radial moments of the electron density within the basin, ∫_ΩA r^k ρ dr for "
             "k = -2, -1, 0, +1, +2 (referenced to the nucleus).",
        columns=("R_M2_QTAIM", "R_M1_QTAIM", "R_0_QTAIM", "R_P1_QTAIM", "R_P2_QTAIM")),
    PropDoc(
        key="GR_QTAIM", name="Atomic radial distortion moments", symbol="GR_{-2..+2}(A)",
        unit="a.u.", type="Atomic", bonus=True,
        text="Radial moments of the density weighted by the kinetic-energy-density distortion "
             "(AIMALL GR_k), for k = -2, -1, 0, +1, +2.",
        columns=("GR_M2_QTAIM", "GR_M1_QTAIM", "GR_0_QTAIM", "GR_P1_QTAIM", "GR_P2_QTAIM")),
    PropDoc(
        key="N_ALPHA_QTAIM", name="Atomic spin populations", symbol="N_α, N_β, N_spin",
        unit="e", type="Atomic", bonus=True,
        text="Alpha and beta electron populations of the basin and their difference (the atomic spin "
             "population). Nonzero only for open-shell references; IQARIS is closed-shell, so these are ~0.",
        columns=("N_ALPHA_QTAIM", "N_BETA_QTAIM", "N_SPIN_QTAIM")),
    PropDoc(
        key="Q_BOND_RECON_QTAIM", name="Reconstructed atomic charge", symbol="q_recon(A)", unit="e",
        type="Atomic", bonus=True,
        text="Atomic charge reconstructed from the bonding (LI/DI) decomposition, used as an internal "
             "consistency check against q_QTAIM.",
        columns=("Q_BOND_RECON_QTAIM",), see_also=("q_QTAIM",)),
    PropDoc(
        key="ESP_IDS_QTAIM", name="Atomic ESP descriptors on the isodensity surface",
        symbol="ESP…_IDS",
        unit="Eh/e (value stats: max/min/avg/mad); Eh·bohr²/e (integrals); "
             "Eh²·bohr²/e² (2nd-moment integrals); (Eh/e)² (variances); bohr² (areas)",
        type="Atomic", bonus=True,
        text="Statistical descriptors of the molecular electrostatic potential sampled over the atom's "
             "0.001 a.u. isodensity surface (IDS): integrals, averages, extrema, variances and mean "
             "absolute deviations, split into positive / negative regions, plus the positive/negative "
             "surface areas. Widely used as reactivity / non-covalent-interaction descriptors.",
        columns=("ESP*_IDS_QTAIM", "ABSESP_INT_IDS_QTAIM", "AREA_ESP_*_IDS_QTAIM")),
    PropDoc(
        key="AREA_IDS_QTAIM", name="Atomic isodensity-envelope area",
        symbol="A_IDS(A)", unit="bohr²", type="Atomic", bonus=True,
        text="Area of atom A's molecular isodensity envelope at the 0.0004 / 0.001 / 0.002 a.u. "
             "density contours (the three AREA_IDS_0004/001/002 columns) -- a plain surface area, "
             "not an ESP statistic.",
        columns=("AREA_IDS_0004_QTAIM", "AREA_IDS_001_QTAIM", "AREA_IDS_002_QTAIM"),
        see_also=("ESP_IDS_QTAIM",)),
    PropDoc(
        key="TWOD2_QTAIM", name="Electron pairs between atoms", symbol="2 D2(A,B)", unit="pairs",
        type="Interatomic", bonus=True,
        text="Total number of electron pairs shared between atoms A and B (twice the pair-density "
             "integral); closely related to the delocalization index.",
        columns=("TWOD2_QTAIM",), see_also=("DI_QTAIM",)),
    PropDoc(
        key="DELOC_AB_QTAIM", name="Shared-electron percentages", symbol="%Deloc(A,B), %Deloc(B,A)",
        unit="percent", type="Interatomic", bonus=True,
        text="Percentage of each atom's delocalized electrons that is shared with the specific "
             "partner (asymmetric: A→B and B→A).",
        columns=("DELOC_AB_QTAIM", "DELOC_BA_QTAIM")),
    PropDoc(
        key="BONDED_QTAIM", name="Bonded flag", symbol="bonded", unit="bool", type="Interatomic",
        bonus=True,
        text="True when a bond path (and BCP) links the atom pair — the QTAIM criterion for a "
             "bonded interaction.",
        columns=("BONDED_QTAIM",), see_also=("BP_QTAIM",)),
    PropDoc(
        key="AREA_IAS_QTAIM", name="Interatomic-surface areas", symbol="Area_IAS", unit="bohr²",
        type="Interatomic", bonus=True,
        text="Area of the interatomic (zero-flux) surface between the atom pair, measured at the "
             "0.0004 / 0.001 / 0.002 a.u. isodensity envelopes.",
        columns=("AREA_IAS_0004_QTAIM", "AREA_IAS_001_QTAIM", "AREA_IAS_002_QTAIM")),
    PropDoc(
        key="MU_BOND_IJ_QTAIM", name="Directed bond dipole contribution",
        symbol="μ_bond(i→j)", unit="e·bohr", type="Interatomic", bonus=True,
        text="Per-pair charge-transfer dipole contribution referenced to the BCP, in both directions "
             "(i→j and j→i). Summed over partners it builds MU_INTER(A).",
        columns=("MU_BOND_IJ_QTAIM", "MU_BOND_JI_QTAIM",   # SI family-key aliases (info <key> resolves)
                 "MU_BOND_IJ_X_QTAIM", "MU_BOND_IJ_Y_QTAIM", "MU_BOND_IJ_Z_QTAIM",
                 "MU_BOND_JI_X_QTAIM", "MU_BOND_JI_Y_QTAIM", "MU_BOND_JI_Z_QTAIM"),
        see_also=("MU_INTER_QTAIM",)),
    PropDoc(
        key="Q_CONTRIB_IJ_QTAIM", name="Directed charge-transfer contribution",
        symbol="q(i→j)", unit="e", type="Interatomic", bonus=True,
        text="Per-pair contribution to the bonding dipole from the partner's net charge, in both "
             "directions.",
        columns=("Q_CONTRIB_IJ_QTAIM", "Q_CONTRIB_JI_QTAIM")),
    PropDoc(
        key="CP_HESSEVEC_QTAIM", name="Hessian eigenvectors at the CP", symbol="e_1, e_2, e_3",
        unit="unit vectors", type="CP", bonus=True,
        text="The three (unit) eigenvectors of the density Hessian at the CP, giving the principal "
             "curvature directions. Stored as CP_HESSEVEC_{1,2,3}_{X,Y,Z}.",
        columns=("CP_HESSEVEC_*",), see_also=("CP_HESSEIG_QTAIM",)),
    PropDoc(
        key="CP_STRESSEVEC_QTAIM", name="Stress-tensor eigenvectors at the CP",
        symbol="s_1, s_2, s_3", unit="unit vectors", type="CP", bonus=True,
        text="The three (unit) eigenvectors of the quantum stress tensor at the CP. Stored as "
             "CP_STRESSEVEC_{1,2,3}_{X,Y,Z}.",
        columns=("CP_STRESSEVEC_*",), see_also=("CP_STRESSEIG_QTAIM",)),
    PropDoc(
        key="CP_RHO_NUC_QTAIM", name="Density and ESP at the nucleus (NACP)",
        symbol="ρ_nuc, ESP_nuc", unit="a.u.", type="CP", bonus=True,
        text="Values evaluated exactly at the nuclear critical point: the electron density at the "
             "nucleus and the electrostatic potential (and its e/n components) at the nucleus.",
        columns=("CP_RHO_NUC_QTAIM", "CP_ESP_NUC_QTAIM", "CP_ESPE_NUC_QTAIM", "CP_ESPN_NUC_QTAIM")),
    PropDoc(
        key="N_AIMINT_QTAIM", name="Integration-quality checks (AIMInt sum vs analytic)",
        symbol="X_AIMInt / X_Analytic", unit="a.u.", type="Molecular", bonus=True,
        text="Molecular totals of N, G, K, L and V_en obtained two ways — by summing the numerically "
             "integrated atomic basins (AIMInt) and from the analytic wavefunction — whose agreement "
             "is a global integration-quality diagnostic.",
        columns=("*_AIMINT_QTAIM", "*_ANALYTIC_QTAIM", "VENT_AIMINT_QTAIM", "VENT_ANALYTIC_QTAIM")),
    PropDoc(
        key="NELEC_MOL_QTAIM", name="Molecular electron counts", symbol="N, N_α, N_β",
        unit="e", type="Molecular", bonus=True,
        text="Total number of electrons of the molecule and the alpha/beta counts (from the "
             "wavefunction).",
        columns=("NELEC_MOL_QTAIM", "NALPHA_MOL_QTAIM", "NBETA_MOL_QTAIM")),

    #  IQA-SQM properties (MOPAC / ENPART), available at 5 semiempirical levels
    PropDoc(
        key="q_SQM", name="Atomic charge (semiempirical)", symbol="q(A)", unit="e", type="Atomic",
        method=M_SQM,
        text="Mulliken-type net atomic charge from the semiempirical density, q_A = Z_A - N_A, with "
             "N_A the gross AO population summed over the orbitals centered on A.",
        eq_text="q(A) = Z_A - N_A",
        eq_latex=r"q_A = Z_A - N_A",
        columns=("q_SQM",), see_also=("q_QTAIM", "N_SQM")),
    PropDoc(
        key="Etot_SQM", name="Total molecular energy (IQA-SQM)", symbol="E_tot", unit="Eh",
        type="Molecular", method=M_SQM,
        text="Total energy of the system as the sum of all intra-atomic and interatomic IQA-SQM "
             "contributions.",
        eq_text="E_tot = E_tot^intra + E_tot^inter",
        eq_latex=r"E_{tot} = E_{tot}^{intra} + E_{tot}^{inter}",
        columns=("Etot_SQM",), see_also=("Etot_intra_SQM", "Etot_inter_SQM")),
    PropDoc(
        key="Etot_intra_SQM", name="Total intra-atomic energy", symbol="E_tot^intra", unit="Eh",
        type="Molecular", method=M_SQM,
        text="Sum of the intra-atomic energies of all atoms.",
        eq_text="E_tot^intra = Σ_A E_intra(A)",
        eq_latex=r"E_{tot}^{intra} = \sum_A E_{intra}(A)",
        columns=("Etot_intra_SQM",), see_also=("Eintra_SQM",)),
    PropDoc(
        key="Etot_inter_SQM", name="Total interatomic energy", symbol="E_tot^inter", unit="Eh",
        type="Molecular", method=M_SQM,
        text="Sum of all pairwise interaction energies between distinct atoms.",
        eq_text="E_tot^inter = Σ_{A<B} E_inter(A,B)",
        eq_latex=r"E_{tot}^{inter} = \sum_{A<B} E_{inter}(A,B)",
        columns=("Etot_inter_SQM",), see_also=("Einter_SQM",)),
    PropDoc(
        key="Eintra_SQM", name="Intra-atomic energy per atom", symbol="E_intra(A)", unit="Eh",
        type="Atomic", method=M_SQM,
        text="Energy of atom A in its molecular environment: its one-center electron–electron "
             "repulsion plus its one-center electron–nuclear + kinetic term.",
        eq_text="E_intra(A) = E_intra^Vee(A) + E_intra^(Ven+T)(A)",
        eq_latex=r"E_{intra}(A) = E_{intra}^{V_{ee}}(A) + E_{intra}^{V_{en}+T}(A)",
        columns=("Eintra_SQM",), see_also=("Eintra_ee_SQM", "Eintra_en_SQM")),
    PropDoc(
        key="Eintra_ee_SQM", name="Intra-atomic electron–electron repulsion",
        symbol="E_intra^Vee(A)", unit="Eh", type="Atomic", method=M_SQM,
        text="One-center electron–electron repulsion within atom A's own basis (the E–E "
             "one-center term).",
        columns=("Eintra_ee_SQM",)),
    PropDoc(
        key="Eintra_en_SQM", name="Intra-atomic electron–nucleus + kinetic energy",
        symbol="E_intra^(Ven+T)(A)", unit="Eh", type="Atomic", method=M_SQM,
        text="One-center, one-electron term of atom A: electron–nuclear attraction together with "
             "electronic kinetic energy (folded into the semiempirical core integrals; the E–N "
             "one-center term).",
        columns=("Eintra_en_SQM",)),
    PropDoc(
        key="Einter_SQM", name="Interatomic interaction energy", symbol="E_inter(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_SQM,
        text="Total interaction energy of atoms A and B, partitioned (Salazar-Lozas et al.) into the "
             "classical electrostatic (V_cl), exchange (V_xc) and semiempirical resonance (V_R) "
             "channels. V_xc + V_R is the covalent channel. Dispersion is a separate post-SCF correction.",
        eq_text="E_inter(A,B) = V_cl(A,B) + V_xc(A,B) + V_R(A,B)",
        eq_latex=r"E_{inter}(A,B) = V_{cl}(A,B) + V_{xc}(A,B) + V_R(A,B)",
        columns=("Einter_SQM",),
        see_also=("Eelstat_SQM", "Eexchange_SQM", "Eresonance_SQM", "Edisp_SQM")),
    PropDoc(
        key="Eelstat_SQM", name="Classical electrostatic (Coulombic) interaction",
        symbol="V_cl(A,B)", unit="Eh", type="Interatomic (pair)", method=M_SQM,
        text="Total classical electrostatic interaction between A and B: the sum of the "
             "electron–electron, electron–nucleus and nucleus–nucleus terms (MOPAC's C column).",
        eq_text="V_cl(A,B) = V_ee(A,B) + V_en(A,B) + V_nn(A,B)",
        eq_latex=r"V_{cl}(A,B) = V_{ee}(A,B) + V_{en}(A,B) + V_{nn}(A,B)",
        columns=("Eelstat_SQM",),
        see_also=("Eelstat_ee_SQM", "Eelstat_en_SQM", "Eelstat_nn_SQM")),
    PropDoc(
        key="Eelstat_ee_SQM", name="Electrostatic components (e–e, e–n, n–n)",
        symbol="V_ee, V_en, V_nn", unit="Eh", type="Interatomic (pair)", method=M_SQM,
        text="The three classical electrostatic sub-terms between atoms A and B: electron–electron "
             "repulsion, electron–nucleus attraction and nucleus–nucleus repulsion. They sum to V_cl.",
        columns=("Eelstat_ee_SQM", "Eelstat_en_SQM", "Eelstat_nn_SQM"), see_also=("Eelstat_SQM",)),
    PropDoc(
        key="Eexchange_SQM", name="Exchange interaction", symbol="V_xc(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_SQM,
        text="Two-center exchange contribution from wavefunction antisymmetry — the quantum, "
             "covalent part of the electron–electron interaction. In a single-determinant "
             "semiempirical method this channel is purely exchange (no separate correlation).",
        columns=("Eexchange_SQM",), see_also=("Einter_SQM", "Eresonance_SQM")),
    PropDoc(
        key="Eresonance_SQM", name="Resonance energy", symbol="V_R(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_SQM,
        text="Two-center one-electron (bond-formation) energy from the off-diagonal core-Hamiltonian "
             "elements H_μν = β_μν S_μν. The dominant covalency descriptor "
             "in semiempirical IQA; it has no direct ab-initio IQA counterpart.",
        columns=("Eresonance_SQM",), see_also=("Einter_SQM", "Eexchange_SQM")),
    # Molecular-total counterparts (the `*_mol_SQM` columns): the whole-molecule sums of the
    # corresponding two-center channel, one value per structure.  Distinct columns/tables from the
    # per-pair terms above, so they get their own entries.
    PropDoc(
        key="Eresonance_mol_SQM", name="Molecular resonance energy (total)", symbol="V_R(mol)",
        unit="Eh", type="Molecular", method=M_SQM,
        text="Whole-molecule sum of the two-center resonance energy V_R(A,B) over all atom pairs "
             "(the molecular-total counterpart of the per-pair Eresonance_SQM).",
        columns=("Eresonance_mol_SQM",), see_also=("Eresonance_SQM", "Exr_SQM")),
    PropDoc(
        key="Eexchange_mol_SQM", name="Molecular exchange energy (total)", symbol="V_xc(mol)",
        unit="Eh", type="Molecular", method=M_SQM,
        text="Whole-molecule sum of the two-center exchange energy V_xc(A,B) over all atom pairs "
             "(the molecular-total counterpart of the per-pair Eexchange_SQM).",
        columns=("Eexchange_mol_SQM",), see_also=("Eexchange_SQM", "Exr_SQM")),
    PropDoc(
        key="Eelstat_mol_SQM", name="Molecular electrostatic interaction (total)", symbol="V_cl(mol)",
        unit="Eh", type="Molecular", method=M_SQM,
        text="Whole-molecule sum of the two-center classical electrostatic interaction V_cl(A,B) "
             "over all atom pairs (the molecular-total counterpart of the per-pair Eelstat_SQM).",
        columns=("Eelstat_mol_SQM",), see_also=("Eelstat_SQM",)),
    PropDoc(
        key="Eadd_SQM", name="Additive atomic energy", symbol="E_add(A)", unit="Eh",
        type="Atomic", method=M_SQM,
        text="Effective energy of atom A: its intra-atomic energy plus half of its interatomic "
             "interactions. Σ_A E_add(A) = E_tot.",
        eq_text="E_add(A) = E_intra(A) + ½ Σ_{B≠A} E_inter(A,B)",
        eq_latex=r"E_{add}(A) = E_{intra}(A) + \tfrac{1}{2}\sum_{B\neq A} E_{inter}(A,B)",
        columns=("Eadd_SQM",)),
    PropDoc(
        key="Edisp_SQM", name="Dispersion correction (pairwise)", symbol="E_disp(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_SQM,
        text="Empirical van-der-Waals correction between atoms A and B (Grimme D3/D2 style), a "
             "post-SCF term outside the IQA partition. Zero for PM6; for PM6-ORG only the molecular "
             "total (Edisp_mol) is stored.",
        eq_text="E_disp(A,B) = -C6^AB / R_AB⁶ · f_damp(R_AB)",
        eq_latex=r"E_{disp}(A,B) = -\frac{C_6^{AB}}{R_{AB}^{6}}\,f_{damp}(R_{AB})",
        columns=("Edisp_SQM",), see_also=("Edisp_mol_SQM",)),

    #  IQA-SQM bonus properties
    PropDoc(
        key="N_SQM", name="Electron population (semiempirical)", symbol="N(A)", unit="e",
        type="Atomic", method=M_SQM, bonus=True,
        text="Gross atomic electron population from the semiempirical density (Σ P·S over the "
             "orbitals on A); N_A = Z_A - q_A.",
        columns=("N_SQM",), see_also=("q_SQM",)),
    PropDoc(
        key="s_pop_SQM", name="Valence s/p/d populations", symbol="N_A^{s,p,d}", unit="e",
        type="Atomic", method=M_SQM, bonus=True,
        text="Partition of the atomic valence population into s, p and d shell occupations.",
        columns=("s_pop_SQM", "p_pop_SQM", "d_pop_SQM")),
    PropDoc(
        key="Eone_center_SQM", name="One-center total energy", symbol="E_1c", unit="Eh",
        type="Molecular", method=M_SQM, bonus=True,
        text="Grand total of the one-center (E–E + E–N) energy terms over all atoms.",
        columns=("Eone_center_SQM",)),
    PropDoc(
        key="Etwo_center_SQM", name="Two-center grand total", symbol="E_2c", unit="Eh",
        type="Molecular", method=M_SQM, bonus=True,
        text="Grand total of all two-center interaction energies.",
        columns=("Etwo_center_SQM",)),
    PropDoc(
        key="Exr_SQM", name="Exchange + resonance energy (total)", symbol="E_xr", unit="Eh",
        type="Molecular", method=M_SQM, bonus=True,
        text="Molecular sum of the exchange and resonance channels — the total covalent energy.",
        columns=("Exr_SQM",), see_also=("Eexchange_SQM", "Eresonance_SQM")),
    PropDoc(
        key="Etot_reconstructed_SQM", name="Reconstructed total energy", symbol="E_tot(recon)",
        unit="Eh", type="Molecular", method=M_SQM, bonus=True,
        text="Total energy rebuilt from the summed IQA components, used as an internal consistency "
             "check against Etot_SQM (the reconstruction residual gates curation).",
        columns=("Etot_reconstructed_SQM",), see_also=("Etot_SQM",)),
    PropDoc(
        key="Edisp_mol_SQM", name="Molecular dispersion energy", symbol="E_disp(mol)", unit="Eh",
        type="Molecular", method=M_SQM, bonus=True,
        text="Total molecular empirical dispersion correction (also the pairwise sum "
             "Edisp_sum_pairwise where a per-pair partition exists).",
        columns=("Edisp_mol_SQM", "Edisp_sum_pairwise_SQM"), see_also=("Edisp_SQM",)),
    PropDoc(
        key="Ehbond_SQM", name="Hydrogen-bond energy", symbol="E_Hbond", unit="Eh",
        type="Molecular", method=M_SQM, bonus=True,
        text="Empirical hydrogen-bond correction energy (H4-type; n_hbonds counts the detected bonds).",
        columns=("Ehbond_SQM",), see_also=("n_hbonds",)),
    PropDoc(
        key="Eelec_SQM", name="Electronic energy", symbol="E_elec", unit="Eh", type="Molecular",
        method=M_SQM, bonus=True,
        text="Total electronic energy from the SCF (before adding the core–core repulsion).",
        columns=("Eelec_SQM",), see_also=("Ecore_SQM",)),
    PropDoc(
        key="Ecore_SQM", name="Core–core repulsion", symbol="E_core", unit="Eh",
        type="Molecular", method=M_SQM, bonus=True,
        text="Total core–core (nuclear) repulsion energy of the semiempirical model.",
        columns=("Ecore_SQM",), see_also=("Eelec_SQM",)),
    PropDoc(
        key="HOF_SQM", name="Heat of formation", symbol="ΔH_f", unit="Eh", type="Molecular",
        method=M_SQM, bonus=True,
        text="Standard semiempirical heat of formation of the molecule (converted to Hartree from "
             "MOPAC's kcal/mol).",
        columns=("HOF_SQM",)),
    PropDoc(
        key="HOMO_SQM", name="Frontier orbital energies and gap", symbol="E_HOMO, E_LUMO, gap",
        unit="Eh", type="Molecular", method=M_SQM, bonus=True,
        text="Highest-occupied and lowest-unoccupied molecular-orbital energies and the HOMO–LUMO "
             "gap.",
        columns=("HOMO_SQM", "LUMO_SQM", "gap_SQM")),
    PropDoc(
        key="IP_SQM", name="Ionization potential", symbol="IP", unit="Eh", type="Molecular",
        method=M_SQM, bonus=True,
        text="Vertical ionization potential (Koopmans, −E_HOMO) reported by MOPAC.",
        columns=("IP_SQM",)),
    PropDoc(
        key="dipole_SQM", name="Molecular dipole moment (semiempirical)", symbol="μ_mol",
        unit="e·bohr", type="Molecular", method=M_SQM, bonus=True,
        text="Total molecular dipole moment and its X/Y/Z components from the semiempirical density.",
        columns=("dipole_SQM", "dipole_x_SQM", "dipole_y_SQM", "dipole_z_SQM")),
    PropDoc(
        key="cosmo_area_SQM", name="COSMO area / volume", symbol="A_COSMO, V_COSMO",
        unit="Angstrom² (area); Angstrom³ (volume)", type="Molecular", method=M_SQM, bonus=True,
        text="Solvent-accessible surface area and enclosed volume from the COSMO continuum-solvation "
             "cavity.",
        columns=("cosmo_area_SQM", "cosmo_vol_SQM")),
    PropDoc(
        key="MW_SQM", name="Molecular weight", symbol="MW", unit="g/mol", type="Molecular",
        method=M_SQM, bonus=True,
        text="Molecular weight of the structure.",
        columns=("MW_SQM",)),
    PropDoc(
        key="point_group_SQM", name="Point group", symbol="—", unit="label", type="Molecular",
        method=M_SQM, bonus=True,
        text="Molecular point-group symmetry label detected by MOPAC.",
        columns=("point_group_SQM",)),
    PropDoc(
        key="n_filled_levels", name="Number of filled MO levels", symbol="—", unit="count",
        type="Molecular", method=M_SQM, bonus=True,
        text="Number of doubly occupied molecular-orbital levels.",
        columns=("n_filled_levels",)),
    PropDoc(
        key="n_hbonds", name="Number of hydrogen bonds", symbol="—", unit="count",
        type="Molecular", method=M_SQM, bonus=True,
        text="Count of hydrogen bonds detected by the semiempirical H-bond model.",
        columns=("n_hbonds",), see_also=("Ehbond_SQM",)),
    PropDoc(
        key="element_pair", name="Pair element labels", symbol="—", unit="label",
        type="Interatomic",
        text="Chemical-element symbols of the two atoms of an interatomic pair row.",
        columns=("element_i", "element_j")),

    #  IQA properties (ab initio Interacting Quantum Atoms), PBE0/def2-TZVP.  From AIMAll
    #  -encomp=4 on the PBE0 wavefunction; a wholly ab initio energy partition, the DFT
    #  counterpart of the semiempirical IQA-SQM block above.  All energies in Hartree.
    #
    #  Two AIMAll conventions the reader must keep straight:
    #   * FULL vs HALVED.  The per-PAIR (i,j) table stores the FULL two-atom quantity
    #     (e.g. VNE_IJ_IQA = V_ne(i,j) between the two atoms).  A handful of PER-ATOM columns
    #     tagged *_HALF are AIMAll's atom-partitioned HALVES (its "/2" tables): the atom-A share
    #     of a two-atom quantity, so that summing the halves over atoms reproduces the pair sum.
    #     Halves are stored for convenience/round-tripping and flagged `bonus`.
    #   * V_ne DIRECTION.  V_ne(i,j) is nucleus(i)-with-electrons(j); V_en(i,j) is
    #     electrons(i)-with-nucleus(j).  The two are NOT equal for i!=j (only their sum, V_neen,
    #     is symmetric).  IQARIS stores both with this direction convention (validated against
    #     AIMAll's halved directional tables to ~1e-10 Ha).

    # ---- IQA (PBE0), molecular (diagnostics / roll-ups; all bonus) -------------------------
    PropDoc(
        key="E_IQA_MOL", name="IQA-recovered molecular energy", symbol="E_IQA(mol)", unit="Eh",
        type="Molecular", method=M_IQADFT, bonus=True,
        text="Total molecular energy reconstructed from the ab initio IQA partition: the sum of "
             "every atom's intra-atomic self-energy plus every distinct interatomic interaction. It "
             "closes to the PBE0 SCF energy to within the AIMAll integration error (DE_IQA_WFN), a "
             "few kJ/mol for these systems.",
        eq_text="E_IQA(mol) = Σ_A E_intra(A) + Σ_{A<B} E_int(A,B)",
        eq_latex=r"E_{IQA}(\mathrm{mol}) = \sum_A E_{intra}(A) + \sum_{A<B} E_{int}(A,B)",
        columns=("E_IQA_MOL",),
        see_also=("E_IQA", "EINT_IQA", "DE_IQA_WFN", "E_PBE0_Eh")),
    PropDoc(
        key="EINTRA_SUM_IQA", name="Total intra-atomic energy", symbol="Σ E_intra", unit="Eh",
        type="Molecular", method=M_IQADFT, bonus=True,
        text="Sum over all atoms of the intra-atomic (self) energy E_intra(A).",
        eq_text="Σ_A E_intra(A)",
        eq_latex=r"\sum_A E_{intra}(A)",
        columns=("EINTRA_SUM_IQA",), see_also=("EINTRA_IQA", "E_IQA_MOL")),
    PropDoc(
        key="EINTER_SUM_IQA", name="Total interatomic energy (atom-sum form)",
        symbol="Σ E_inter", unit="Eh", type="Molecular", method=M_IQADFT, bonus=True,
        text="Sum over all atoms of the interatomic energy share E_inter(A) (= half of each atom's "
             "interactions). Equals the pair-sum form EINT_PAIRSUM_IQA to the per-atom closure.",
        eq_text="Σ_A E_inter(A) = Σ_{A<B} E_int(A,B)",
        eq_latex=r"\sum_A E_{inter}(A) = \sum_{A<B} E_{int}(A,B)",
        columns=("EINTER_SUM_IQA",), see_also=("EINTER_IQA", "EINT_PAIRSUM_IQA")),
    PropDoc(
        key="EINT_PAIRSUM_IQA", name="Total interatomic energy (pair-sum form)",
        symbol="Σ E_int", unit="Eh", type="Molecular", method=M_IQADFT, bonus=True,
        text="Sum over all distinct atom pairs of the full interaction energy E_int(i,j); the "
             "pair-table counterpart of EINTER_SUM_IQA.",
        eq_text="Σ_{i<j} E_int(i,j)",
        eq_latex=r"\sum_{i<j} E_{int}(i,j)",
        columns=("EINT_PAIRSUM_IQA",), see_also=("EINTER_SUM_IQA", "EINT_IQA")),
    PropDoc(
        key="VCL_SUM_IQA", name="Total classical interaction energy", symbol="Σ V_cl", unit="Eh",
        type="Molecular", method=M_IQADFT, bonus=True,
        text="Sum over all distinct atom pairs of the classical (Coulombic) interaction V_cl(i,j).",
        eq_text="Σ_{i<j} V_cl(i,j)",
        eq_latex=r"\sum_{i<j} V_{cl}(i,j)",
        columns=("VCL_SUM_IQA",), see_also=("VCL_IJ_IQA", "VXC_SUM_IQA")),
    PropDoc(
        key="VXC_SUM_IQA", name="Total exchange-correlation interaction energy",
        symbol="Σ V_xc", unit="Eh", type="Molecular", method=M_IQADFT, bonus=True,
        text="Sum over all distinct atom pairs of the exchange-correlation interaction V_xc(i,j) — "
             "the molecule's total ab initio covalent (quantum) interaction energy.",
        eq_text="Σ_{i<j} V_xc(i,j)",
        eq_latex=r"\sum_{i<j} V_{xc}(i,j)",
        columns=("VXC_SUM_IQA",), see_also=("VXC_IJ_IQA", "VCL_SUM_IQA")),
    PropDoc(
        key="DE_IQA_WFN", name="IQA closure residual vs SCF", symbol="ΔE_IQA", unit="Eh",
        type="Molecular", method=M_IQADFT, bonus=True,
        text="Difference between the IQA-recovered molecular energy and the PBE0 SCF energy — the "
             "AIMAll integration/reconstruction error. Its magnitude gates the IQA curation "
             "(|ΔE_IQA| ≤ max(2.0, 0.5·N_atoms) kJ/mol); it is computed in double precision at parse "
             "time and is the empirical bound quoted for the RIJCOSX approximation.",
        eq_text="ΔE_IQA = E_IQA(mol) − E_SCF(PBE0)",
        eq_latex=r"\Delta E_{IQA} = E_{IQA}(\mathrm{mol}) - E_{SCF}(\mathrm{PBE0})",
        columns=("DE_IQA_WFN",), see_also=("E_IQA_MOL", "DE_IQA_RECON")),
    PropDoc(
        key="DE_IQA_RECON", name="IQA atom-vs-pair reconstruction residual",
        symbol="ΔE_recon", unit="Eh", type="Molecular", method=M_IQADFT, bonus=True,
        text="AIMAll's own consistency residual between the atom-summed and pair-summed forms of the "
             "interatomic energy — an internal integration check, distinct from the SCF closure "
             "DE_IQA_WFN.",
        columns=("DE_IQA_RECON",), see_also=("DE_IQA_WFN",)),

    # ---- IQA (PBE0), atomic (headline) ------------------------------------------------------
    PropDoc(
        key="E_IQA", name="IQA atomic energy", symbol="E_IQA(A)", unit="Eh", type="Atomic",
        method=M_IQADFT,
        text="Additive ab initio energy of atom A in the molecule: its intra-atomic self-energy plus "
             "its (halved) share of every interatomic interaction. The atomic energies sum to the "
             "molecular energy — the ab initio analogue of the semiempirical E_add(A).",
        eq_text="E_IQA(A) = E_intra(A) + E_inter(A),   Σ_A E_IQA(A) = E_IQA(mol)",
        eq_latex=r"E_{IQA}(A) = E_{intra}(A) + E_{inter}(A)",
        columns=("E_IQA",), see_also=("EINTRA_IQA", "EINTER_IQA", "Eadd_SQM", "E_IQA_MOL")),
    PropDoc(
        key="EINTRA_IQA", name="Intra-atomic (self) energy", symbol="E_intra(A)", unit="Eh",
        type="Atomic", method=M_IQADFT,
        text="Self-energy of atom A: its electronic kinetic energy plus the electron-nucleus "
             "attraction and electron-electron repulsion contained entirely within its own basin.",
        eq_text="E_intra(A) = T(A) + V_ne(A,A) + V_ee(A,A)",
        eq_latex=r"E_{intra}(A) = T(A) + V_{ne}(A,A) + V_{ee}(A,A)",
        columns=("EINTRA_IQA",), see_also=("VNE_AA_IQA", "VEE_AA_IQA", "Eintra_SQM")),
    PropDoc(
        key="EINTER_IQA", name="Interatomic energy share", symbol="E_inter(A)", unit="Eh",
        type="Atomic", method=M_IQADFT,
        text="Atom A's additive share of all its interatomic interactions — half of the sum of "
             "E_int(A,B) over every partner B (the other half is charged to B).",
        eq_text="E_inter(A) = ½ Σ_{B≠A} E_int(A,B)",
        eq_latex=r"E_{inter}(A) = \tfrac{1}{2}\sum_{B\neq A} E_{int}(A,B)",
        columns=("EINTER_IQA",), see_also=("E_IQA", "EINT_IQA")),
    PropDoc(
        key="V_IQA", name="Atomic potential energy", symbol="V(A)", unit="Eh", type="Atomic",
        method=M_IQADFT,
        text="Total potential-energy contribution of atom A (all electron-nucleus, electron-electron "
             "and nucleus-nucleus terms charged to A), split into classical and exchange-correlation "
             "parts.",
        eq_text="V(A) = V_C(A) + V_X(A)",
        eq_latex=r"V(A) = V_{C}(A) + V_{X}(A)",
        columns=("V_IQA",), see_also=("VC_IQA", "VX_IQA")),
    PropDoc(
        key="VC_IQA", name="Atomic classical potential", symbol="V_C(A)", unit="Eh", type="Atomic",
        method=M_IQADFT,
        text="Classical (Coulombic) part of atom A's potential energy — its electron-nucleus, "
             "Coulomb electron-electron and nucleus-nucleus contributions.",
        columns=("VC_IQA",), see_also=("V_IQA", "VX_IQA")),
    PropDoc(
        key="VX_IQA", name="Atomic exchange-correlation potential", symbol="V_X(A)", unit="Eh",
        type="Atomic", method=M_IQADFT,
        text="Exchange-correlation part of atom A's potential energy — the quantum (covalent) "
             "channel charged to A.",
        columns=("VX_IQA",), see_also=("V_IQA", "VC_IQA")),
    PropDoc(
        key="V_INTRA_IQA", name="Intra-atomic potential (total / classical)",
        symbol="V_intra(A)", unit="Eh", type="Atomic", method=M_IQADFT,
        text="The intra-atomic (within-basin) part of atom A's potential energy, and its classical "
             "component V_C^intra(A).",
        columns=("V_INTRA_IQA", "VC_INTRA_IQA"), see_also=("V_IQA", "EINTRA_IQA")),
    PropDoc(
        key="V_INTER_IQA", name="Interatomic potential (total / classical / xc)",
        symbol="V_inter(A)", unit="Eh", type="Atomic", method=M_IQADFT,
        text="The interatomic part of atom A's potential energy (its share of the between-atom "
             "potential), with its classical V_C^inter(A) and exchange-correlation V_X^inter(A) "
             "components.",
        columns=("V_INTER_IQA", "VC_INTER_IQA", "VX_INTER_IQA"),
        see_also=("V_IQA", "EINTER_IQA")),
    PropDoc(
        key="VNEEN_IQA", name="Atomic electron-nucleus attraction (total)",
        symbol="V_neen(A)", unit="Eh", type="Atomic", method=M_IQADFT,
        text="Total electron-nucleus attraction charged to atom A, summing the nucleus-with-electrons "
             "and electrons-with-nucleus contributions over the whole molecule.",
        eq_text="V_neen(A) = V_ne(A,mol) + V_en(A,mol)",
        eq_latex=r"V_{neen}(A) = V_{ne}(A,\mathrm{mol}) + V_{en}(A,\mathrm{mol})",
        columns=("VNEEN_IQA",), see_also=("VNE_AM_HALF_IQA", "VEN_AM_HALF_IQA")),

    # ---- IQA (PBE0), atomic (bonus: intra-atomic components + AIMAll halves) ----------------
    PropDoc(
        key="VEE_A_IQA", name="Atomic electron-electron repulsion (total / C / X)",
        symbol="V_ee(A)", unit="Eh", type="Atomic", method=M_IQADFT, bonus=True,
        text="Total electron-electron repulsion charged to atom A and its Coulomb (V_ee^C) and "
             "exchange-correlation (V_ee^X) parts.",
        columns=("VEE_A_IQA", "VEEC_A_IQA", "VEEX_A_IQA"), see_also=("VX_IQA",)),
    PropDoc(
        key="VNE_AA_IQA", name="Intra-atomic electron-nucleus attraction",
        symbol="V_ne(A,A)", unit="Eh", type="Atomic", method=M_IQADFT, bonus=True,
        text="Attraction between the electrons and the nucleus of atom A, both within A's own basin "
             "(a component of E_intra(A)).",
        columns=("VNE_AA_IQA",), see_also=("EINTRA_IQA",)),
    PropDoc(
        key="VEE_AA_IQA", name="Intra-atomic electron-electron repulsion (total / C / X)",
        symbol="V_ee(A,A)", unit="Eh", type="Atomic", method=M_IQADFT, bonus=True,
        text="Electron-electron repulsion within atom A's own basin and its Coulomb (V_ee^C(A,A)) "
             "and exchange-correlation (V_ee^X(A,A)) parts. V_ee^X(A,A) is atom A's self-exchange, "
             "identical to the atomic exchange V_X^intra(A).",
        columns=("VEE_AA_IQA", "VEEC_AA_IQA", "VEEX_AA_IQA"), see_also=("EINTRA_IQA", "VX_IQA")),
    PropDoc(
        key="VNE_AM_HALF_IQA", name="Atom-molecule e-n attraction, halved",
        symbol="V_ne(A,mol)/2, V_en(A,mol)/2", unit="Eh", type="Atomic", method=M_IQADFT, bonus=True,
        text="AIMAll's atom-partitioned halves of the electron-nucleus attraction between atom A and "
             "the whole molecule: V_ne (nucleus A with all electrons) and V_en (electrons of A with "
             "all nuclei), each divided by two so the atomic shares sum to the molecular total.",
        columns=("VNE_AM_HALF_IQA", "VEN_AM_HALF_IQA"),
        see_also=("VNEEN_IQA", "VNEEN_INTER_HALF_IQA")),
    PropDoc(
        key="VNEEN_INTER_HALF_IQA", name="Interatomic e-n attraction share, halved",
        symbol="V_neen^inter(A)/2", unit="Eh", type="Atomic", method=M_IQADFT, bonus=True,
        text="Atom A's halved share of the interatomic electron-nucleus attraction (AIMAll '/2' "
             "convention).",
        columns=("VNEEN_INTER_HALF_IQA",), see_also=("VNEEN_IQA",)),
    PropDoc(
        key="VNE_AAP_HALF_IQA", name="Interatomic pair components, halved (per atom)",
        symbol="V_·(A,A')/2", unit="Eh", type="Atomic", method=M_IQADFT, bonus=True,
        text="AIMAll's per-atom halves of the six interatomic pair components summed over A's "
             "partners A': electron-nucleus (both directions), electron-electron (total, Coulomb, "
             "exchange-correlation) and nucleus-nucleus. Each is half of the corresponding pair-table "
             "quantity; provided for round-tripping the atomic and pair tables.",
        columns=("VNE_AAP_HALF_IQA", "VEN_AAP_HALF_IQA", "VEE_AAP_HALF_IQA",
                 "VNN_AAP_HALF_IQA", "VEEC_AAP_HALF_IQA", "VEEX_AAP_HALF_IQA"),
        see_also=("VNE_IJ_IQA", "VEE_IJ_IQA", "VNN_IJ_IQA")),

    # ---- IQA (PBE0), pairs (headline; FULL two-atom quantities) -----------------------------
    PropDoc(
        key="EINT_IQA", name="IQA interaction energy", symbol="E_int(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_IQADFT,
        text="Total ab initio interaction energy between atoms A and B, exactly split into a "
             "classical (Coulombic) and an exchange-correlation channel — the DFT counterpart of the "
             "semiempirical E_inter(A,B). V_xc is the covalent channel; there is no separate "
             "resonance term and no empirical dispersion.",
        eq_text="E_int(A,B) = V_cl(A,B) + V_xc(A,B)",
        eq_latex=r"E_{int}(A,B) = V_{cl}(A,B) + V_{xc}(A,B)",
        columns=("EINT_IQA",), see_also=("VCL_IJ_IQA", "VXC_IJ_IQA", "Einter_SQM")),
    PropDoc(
        key="VCL_IJ_IQA", name="Classical interaction (Coulombic)", symbol="V_cl(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_IQADFT,
        text="Classical electrostatic interaction between the two atoms: the electron-nucleus "
             "(both directions), Coulomb electron-electron and nucleus-nucleus terms.",
        eq_text="V_cl(A,B) = V_ne(A,B) + V_en(A,B) + V_ee^C(A,B) + V_nn(A,B)",
        eq_latex=r"V_{cl}(A,B) = V_{ne}(A,B) + V_{en}(A,B) + V_{ee}^{C}(A,B) + V_{nn}(A,B)",
        columns=("VCL_IJ_IQA",),
        see_also=("VNE_IJ_IQA", "VEN_IJ_IQA", "VEEC_IJ_IQA", "VNN_IJ_IQA", "Eelstat_SQM")),
    PropDoc(
        key="VXC_IJ_IQA", name="Exchange-correlation interaction", symbol="V_xc(A,B)", unit="Eh",
        type="Interatomic (pair)", method=M_IQADFT,
        text="Exchange-correlation interaction between atoms A and B — the quantum, covalent part of "
             "E_int, equal to the interatomic exchange-correlation electron-electron term V_ee^X(A,B). "
             "The ab initio counterpart of the semiempirical exchange V_xc(A,B), and the strongest "
             "cross-fidelity correlate of the QTAIM delocalization index.",
        eq_text="V_xc(A,B) = V_ee^X(A,B)",
        eq_latex=r"V_{xc}(A,B) = V_{ee}^{X}(A,B)",
        columns=("VXC_IJ_IQA",), see_also=("EINT_IQA", "Eexchange_SQM", "DI_QTAIM")),
    PropDoc(
        key="VNE_IJ_IQA", name="Electron-nucleus attraction (directional)",
        symbol="V_ne(A,B), V_en(A,B)", unit="Eh", type="Interatomic (pair)", method=M_IQADFT,
        text="The two directional electron-nucleus attractions of the pair: V_ne is nucleus A with "
             "the electrons of B, V_en is the electrons of A with nucleus B. They differ for A≠B "
             "(only their sum is symmetric).",
        eq_text="V_ne(A,B) = nucleus(A)·electrons(B);  V_en(A,B) = electrons(A)·nucleus(B)",
        eq_latex=r"V_{ne}(A,B)\ [\mathrm{nuc}_A\!\cdot\!e_B],\quad V_{en}(A,B)\ [e_A\!\cdot\!\mathrm{nuc}_B]",
        columns=("VNE_IJ_IQA", "VEN_IJ_IQA"), see_also=("VNEEN_IJ_IQA", "VCL_IJ_IQA")),
    PropDoc(
        key="VEE_IJ_IQA", name="Interatomic electron-electron repulsion (total / Coulomb)",
        symbol="V_ee(A,B)", unit="Eh", type="Interatomic (pair)", method=M_IQADFT,
        text="Total electron-electron repulsion between the two atoms' basins and its Coulomb part "
             "V_ee^C(A,B). The exchange-correlation part is V_xc(A,B); V_ee = V_ee^C + V_xc.",
        eq_text="V_ee(A,B) = V_ee^C(A,B) + V_ee^X(A,B)",
        eq_latex=r"V_{ee}(A,B) = V_{ee}^{C}(A,B) + V_{ee}^{X}(A,B)",
        columns=("VEE_IJ_IQA", "VEEC_IJ_IQA"), see_also=("VXC_IJ_IQA", "VCL_IJ_IQA")),
    PropDoc(
        key="VNN_IJ_IQA", name="Interatomic nucleus-nucleus repulsion", symbol="V_nn(A,B)",
        unit="Eh", type="Interatomic (pair)", method=M_IQADFT,
        text="Nucleus-nucleus Coulomb repulsion between atoms A and B, Z_A·Z_B/R_AB. Summed over all "
             "pairs it reproduces the molecular nuclear-repulsion energy (a parser-integrity check).",
        eq_text="V_nn(A,B) = Z_A Z_B / R_AB",
        eq_latex=r"V_{nn}(A,B) = Z_A Z_B / R_{AB}",
        columns=("VNN_IJ_IQA",), see_also=("VCL_IJ_IQA",)),

    # ---- IQA (PBE0), pairs (bonus roll-up) --------------------------------------------------
    PropDoc(
        key="VNEEN_IJ_IQA", name="Interatomic electron-nucleus attraction (sum)",
        symbol="V_neen(A,B)", unit="Eh", type="Interatomic (pair)", method=M_IQADFT, bonus=True,
        text="Sum of the two directional electron-nucleus attractions of the pair, V_ne(A,B)+V_en(A,B) "
             "— the symmetric total that enters V_cl.",
        columns=("VNEEN_IJ_IQA",), see_also=("VNE_IJ_IQA", "VCL_IJ_IQA")),
]


#  Lookup + formatting
_EXACT: dict = {}
_PATTERNS: list = []          # (pattern, PropDoc) for wildcard columns
_BY_KEY: dict = {}

for _d in DOCS:
    _BY_KEY[_d.key] = _d
    for _c in (_d.key,) + tuple(_d.columns):
        if "*" in _c or "?" in _c:
            _PATTERNS.append((_c, _d))
        else:
            _EXACT.setdefault(_c, _d)


def _strip_method_suffix(column: str) -> str:
    """`Etot_SQM__PM7` -> `Etot_SQM` (the all-sqm exporter appends `__<METHOD>`)."""
    return column.split("__", 1)[0]


def lookup(column: str):
    """Return the :class:`PropDoc` documenting a database column (or its public Key), or None.

    Handles component columns (dipole X/Y/Z, Hessian 1/2/3, ...) and the `__METHOD` suffix
    that the all-levels exporter adds to SQM columns.
    """
    if not column:
        return None
    for cand in (column, _strip_method_suffix(column)):
        if cand in _EXACT:
            return _EXACT[cand]
        if cand in _BY_KEY:
            return _BY_KEY[cand]
    base = _strip_method_suffix(column)
    for pat, doc in _PATTERNS:
        if fnmatchcase(column, pat) or fnmatchcase(base, pat):
            return doc
    return None


import re as _re
# component tokens of a vector/tensor family: Cartesian axes (incl. lowercase dipole_x),
# quadrupole pairs, and Hessian/stress eigen-index.
_AXIS_RE = _re.compile(r"(?:^|_)(XX|XY|XZ|YY|YZ|ZZ|[XYZxyz])(?:_|$)")
_IDX_RE = _re.compile(r"(?:HESSEIG|STRESSEIG|HESSEVEC|STRESSEVEC|EIG|EVEC)_?([123])")


def _component_note(column: str, doc) -> str:
    """A short component label ('X', 'XX', '#2', ...) when `column` is a specific member of a
    multi-column vector/tensor family, else None -- so each component's own text can say which
    piece of the vector/tensor it is instead of repeating the parent/magnitude blurb."""
    if not doc or len(doc.columns) < 2:
        return None
    base = _strip_method_suffix(column)
    m = _IDX_RE.search(base)
    if m:
        return f"#{m.group(1)}"
    m = _AXIS_RE.search(base)
    if m:
        return m.group(1).upper()
    return None


def format_text(column: str, unit: str = None, level_aware: bool = None, width: int = 78) -> str:
    """Human-readable multi-line block for a column (terminal / `info`).  `unit` overrides the
    family unit with the exact per-column unit from the live registry when provided."""
    doc = lookup(column)
    lines = []
    rule = "─" * min(width, 60)
    if doc is None:
        lines.append(column)
        lines.append(rule)
        lines.append("(no dictionary entry; see `list-props` for its family and unit)")
        if unit:
            lines.append(f"Unit   : {unit}")
        return "\n".join(lines)
    comp = _component_note(column, doc)
    head = f"{column}"
    if doc.symbol:
        head += f"   {doc.symbol}" + (f",{comp.lower()}" if comp else "")
    lines.append(head)
    lines.append(rule)
    name = doc.name + (f" — {comp} component" if comp else "")
    lines.append(name + ("  [bonus]" if doc.bonus else ""))
    lines.append("")
    meta = [("Unit", unit or doc.unit), ("Type", doc.type), ("Method", doc.method)]
    for label, val in meta:
        if val:
            lines.append(f"{label:7}: {val}")
    lines.append("")
    lines.extend(_wrap(doc.text, width))
    if doc.eq_text:
        lines.append("")
        lines.append("    " + doc.eq_text)
    if doc.see_also:
        lines.append("")
        lines.append("See also: " + ", ".join(doc.see_also))
    return "\n".join(lines)


def _wrap(text: str, width: int):
    import textwrap
    return textwrap.wrap(text, width=width) or [""]


def _method_group(d):
    """Three-way level bucket for the glossary: QTAIM (M06-2X) / ab initio IQA (PBE0) /
    IQA-SQM.  Uses object identity of the shared method sentinels."""
    if d.method is M_SQM:
        return "sqm"
    if d.method is M_IQADFT:
        return "iqa"
    return "qtaim"


def as_records():
    """Every entry as a plain dict (for JSON export of the whole dictionary)."""
    return [d.to_dict() for d in DOCS]


def _glossary_unit(d) -> str:
    """Unit text for a glossary entry.  A single unit is shown verbatim; a family PropDoc whose `unit`
    is a composite ';'/','-joined description (AREA/VOL, CP_V*, ESP…_IDS, COSMO, IAS) is expanded to
    the exact per-column units (DB-free, via registry._resolve_unit) so the document never shows the
    family blob as if it were one column's unit."""
    u = d.unit or ""
    if u and (";" in u or "," in u) and d.columns:
        from .registry import _resolve_unit          # local import: registry imports this module
        per = {}
        for c in d.columns:
            per.setdefault(_resolve_unit(c), []).append(c)
        return "; ".join(f"{unit} ({', '.join(cols)})" for unit, cols in per.items())
    return u


def glossary(fmt: str = "markdown") -> str:
    """Render the whole dictionary as a standalone reference document.

    fmt: 'markdown' (default) or 'html'.
    """
    # Group by the REAL three-way level distinction, not a binary SQM/not test: the ab
    # initio IQA (PBE0) entries carry M_IQADFT and must render in their own section, not merged
    # into the QTAIM (M06-2X) topology one.
    fams = {}
    for d in DOCS:
        fams.setdefault((_method_group(d), d.bonus), []).append(d)
    order = [("qtaim", False), ("qtaim", True), ("iqa", False), ("iqa", True),
             ("sqm", False), ("sqm", True)]
    titles = {
        ("qtaim", False): "QTAIM properties",       ("qtaim", True): "QTAIM bonus properties",
        ("iqa", False): "ab initio IQA properties", ("iqa", True): "ab initio IQA bonus properties",
        ("sqm", False): "IQA-SQM properties",       ("sqm", True): "IQA-SQM bonus properties",
    }
    if fmt == "html":
        return _glossary_html(order, titles, fams)
    out = ["# IQARIS property dictionary\n",
           "One entry per database property: name, symbol, unit, type, level of theory, "
           "description and formal equation. Generated from `property_docs.py`.\n"]
    for grp in order:
        items = fams.get(grp)
        if not items:
            continue
        out.append(f"\n## {titles[grp]}\n")
        for d in items:
            out.append(f"### {d.name} — `{d.key}`")
            bits = [b for b in (f"**Symbol:** {d.symbol}" if d.symbol else "",
                                f"**Unit:** {_glossary_unit(d)}" if d.unit else "",
                                f"**Type:** {d.type}" if d.type else "",
                                f"**Method:** {d.method}") if b]
            out.append("  \n".join(bits))
            out.append("")
            out.append(d.text)
            if d.eq_latex:
                out.append(f"\n$$ {d.eq_latex} $$")
            if d.columns:
                out.append(f"\nColumns: {', '.join('`%s`' % c for c in d.columns)}")
            out.append("")
    return "\n".join(out)


def _glossary_html(order, titles, fams):
    esc = lambda s: (s or "").replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    out = ["<h1>IQARIS property dictionary</h1>",
           "<p>One entry per database property: name, symbol, unit, type, level of theory, "
           "description and formal equation.</p>"]
    for grp in order:
        items = fams.get(grp)
        if not items:
            continue
        out.append(f"<h2>{titles[grp]}</h2>")
        for d in items:
            out.append(f"<h3>{esc(d.name)} <code>{esc(d.key)}</code></h3>")
            out.append("<p>" + " &middot; ".join(
                x for x in (f"<b>{esc(d.symbol)}</b>" if d.symbol else "",
                            f"Unit: {esc(_glossary_unit(d))}" if d.unit else "",
                            f"Type: {esc(d.type)}" if d.type else "",
                            f"Method: {esc(d.method)}") if x) + "</p>")
            out.append(f"<p>{esc(d.text)}</p>")
            if d.eq_text:
                out.append(f"<pre>{esc(d.eq_text)}</pre>")
            if d.columns:                       # the DB-column mapping, as in markdown
                out.append("<p>Columns: " + ", ".join(
                    f"<code>{esc(c)}</code>" for c in d.columns) + "</p>")
    return "\n".join(out)
