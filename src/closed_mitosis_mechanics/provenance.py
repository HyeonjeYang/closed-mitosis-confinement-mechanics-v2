"""Evidence-quality status vocabulary and primary-source references.

Every physical parameter this package reads carries a `status` string from
this module, so its evidential status is explicit in the input metadata.
"""

from __future__ import annotations


class Status:
    """Evidence-quality tiers. A value's status never changes based on
    whether it is convenient for a particular calculation."""

    PRIMARY_SOURCE_VERIFIED = "PRIMARY_SOURCE_VERIFIED"
    """Read directly from a primary source statement or equation. Any related
    convention cross-check is described explicitly in the parameter metadata."""

    DERIVED = "DERIVED"
    """An exact arithmetic/algebraic consequence of PRIMARY_SOURCE_VERIFIED values."""

    MODEL_CALIBRATED = "MODEL_CALIBRATED"
    """A baseline scenario value adopted for illustration, not a directly
    measured universal constant (e.g. the 31-pN spindle buckling threshold)."""

    UNRESOLVED_PROVENANCE = "UNRESOLVED_PROVENANCE"
    """A value cited to a source but not reproducible from that source's own
    stated equations under any interpretation identified. Never a confirmed
    error without author clarification; never used in main calculations."""

    CONTEXT_ONLY = "CONTEXT_ONLY"
    """A value from a different biological system, included only to bound a
    plausible pressure scale, never as a direct measurement of the system
    under study."""


# DOIs for the primary sources this package's parameter provenance rests on.
LIM_2007_DOI = "10.1371/journal.pone.0000948"
"""Lim GHW, Huber G, Torii Y, Hirata A, Miller J, Sazer S (2007). Vesicle-Like
Biomechanics Governs Important Aspects of Nuclear Geometry in Fission Yeast.
PLoS ONE 2(9):e948. Fits the membrane tension-to-bending ratio of 65 um^-2 and adopts a typical
B=4.00e-19 J attributed to Powers et al. 2002. The absolute B is an adopted
value, not an absolute modulus inferred by the shape fit."""

LIM_HUBER_2009_DOI = "10.1016/j.bpj.2008.10.074"
"""Lim GHW, Huber G (2009). The Tethered Infinitesimal Tori and Spheres
Algorithm. Biophysical Journal 96(6):2064-2081. Eq. 76 gives the whole-NE
energy F_NE=4*B1*int(H^2)dA+2*sigma_NE*A with B1=200 pN*nm (per bilayer),
sigma_NE=0.013 pN/nm (per bilayer). This is convention-consistent with the
Lim 2007 fitted ratio and adopted absolute modulus; it is not an independent
fit. The single-surface mapping gives kappa_eff=2*B1=400 pN*nm and
sigma_eff=2*sigma_NE=0.026 pN/nm."""

BEGLEY_2026_DOI = "10.1002/cm.22035"
"""Begley et al. (2026). Mechanical Coupling With the Nuclear Envelope Shapes
the Schizosaccharomyces pombe Mitotic Spindle. Cytoskeleton 83(1):40-52.
Table 1 reports kappa=40 pN*nm, sigma=0.013 pN/nm (status=UNRESOLVED_PROVENANCE
here; reproduces their own reported f0=6.4 pN exactly, but neither cited Lim
source contains a 40 pN*nm value under any interpretation identified)."""

WARD_2014_DOI = "10.7554/eLife.03398"
"""Ward JJ, Roque H, Antony C, Nedelec F (2014). Mechanical design principles
of a mitotic spindle. eLife 3:e03398. Figure 4E wild-type fit F_c=A*Ls^-4;
Begley et al. 2026 extracted A=1.94e4 pN*um^4 from this figure (not published
as a number in Ward 2014's own text)."""

SORCE_2015_DOI = "10.1038/ncomms9872"
"""Sorce B et al. (2015). Mitotic cell rounding forces stress the cell cortex.
Nature Communications 6:8872. Source of the 170-890 Pa animal-epithelial
contextual pressure range (status=CONTEXT_ONLY; not an S. pombe measurement)."""
