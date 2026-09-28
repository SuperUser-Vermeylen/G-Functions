from __future__ import annotations

import numpy as np
import pandas as pd
import sympy as sym


COLS = [
    "Temperature",
    "alphaBeta_L",
    "alphaBeta_R",
    "alphaLiquid_L",
    "alphaLiquid_R",
    "betaLiquid_L",
    "betaLiquid_R",
]


def build_phase_diagram() -> pd.DataFrame:
    """Build the temperature-dependent phase diagram used by the plotting logic."""
    TE, XEL, XEFCC, XEBCT, XS, XL, XS1, XS2 = sym.symbols("TE, XEL, XEFCC, XEBCT, XS, XL, XS1, XS2")

    G_Pb_fccL = 4810 - 8.017 * TE
    G_Sn_bctL = 7179 - 14.216 * TE
    G_Pb_fccbct = 489 + 3.52 * TE
    G_Sn_bctfcc = 5510 - 8.46 * TE

    G_Pb_fccL_log = 8.314 * TE * sym.log((1 - XEL) / (1 - XEFCC))
    G_Sn_bctL_log = 8.314 * TE * sym.log((XEL) / (XEBCT))
    G_Pb_fccbct_log = 8.314 * TE * sym.log((1 - XEBCT) / (1 - XEFCC))
    G_Sn_bctfcc_log = 8.314 * TE * sym.log((XEFCC) / (XEBCT))

    G_Pb_fccL_RK = 4700 * XEL ** 2 - 5200 * XEFCC ** 2
    G_Sn_bctL_RK = 4700 * (1 - XEL) ** 2 - 12000 * (1 - XEBCT) ** 2
    G_Pb_fccbct_RK = 12000 * XEBCT ** 2 - 5200 * XEFCC ** 2
    G_Sn_bctfcc_RK = 5200 * (1 - XEFCC) ** 2 - 12000 * (1 - XEBCT) ** 2

    eqn1 = G_Pb_fccL + G_Pb_fccL_log + G_Pb_fccL_RK
    eqn2 = G_Sn_bctL + G_Sn_bctL_log + G_Sn_bctL_RK
    eqn3 = G_Pb_fccbct + G_Pb_fccbct_log + G_Pb_fccbct_RK
    eqn4 = G_Sn_bctfcc + G_Sn_bctfcc_log + G_Sn_bctfcc_RK

    eut_T, eut_XFCC, eut_XL, eut_XBCT = sym.nsolve(
        (eqn1, eqn2, eqn3, eqn4),
        (TE, XEFCC, XEL, XEBCT),
        (724, 0.1, 0.5, 0.7),
    )

    phase_diagram = pd.DataFrame(
        [[eut_T, eut_XFCC, eut_XBCT, np.nan, eut_XL, np.nan, np.nan]],
        columns=COLS,
    )
    phase_diagram = phase_diagram.set_index("Temperature")

    for temp in np.arange(300, 445, 1):
        GFCCBCTPb = 489 + 3.52 * temp
        GBCTFCCSn = 5510 - 8.468 * temp

        logFCCBCTPb = 8.314 * temp * sym.log((1 - XS2) / (1 - XS1))
        logFCCBCTSn = 8.314 * temp * sym.log(XS2 / XS1)

        RKFCCBCTPb = 12000 * XS2 ** 2 - 5200 * XS1 ** 2
        RKFCCBCTSn = 12000 * (1 - XS2) ** 2 - 5200 * (1 - XS1) ** 2

        eqn1 = GFCCBCTPb + logFCCBCTPb + RKFCCBCTPb
        eqn2 = -1 * GBCTFCCSn + logFCCBCTSn + RKFCCBCTSn

        x1, x2 = sym.nsolve((eqn1, eqn2), (XS1, XS2), (eut_XFCC, eut_XBCT))
        row = pd.DataFrame(
            [[temp, x1, x2, np.nan, np.nan, np.nan, np.nan]],
            columns=COLS,
        )
        row = row.set_index("Temperature")
        phase_diagram = pd.concat([phase_diagram, row])

    for temp in np.arange(446, 600, 1):
        GFCCLPb = 4810 - 8.017 * temp
        GBCTLSn = 7179 - 14.126 * temp
        GBCTFCCSn = 5510 - 8.468 * temp

        GFCCLSn = GBCTLSn - GBCTFCCSn

        logFCCLPb = 8.314 * temp * sym.log((1 - XL) / (1 - XS))
        logFCCLSn = 8.314 * temp * sym.log(XL / XS)

        RKFCCLPb = 4700 * XL ** 2 - 5200 * XS ** 2
        RKFCCLSn = 4700 * (1 - XL) ** 2 - 5200 * (1 - XS) ** 2

        eqn1 = GFCCLPb + logFCCLPb + RKFCCLPb
        eqn2 = GFCCLSn + logFCCLSn + RKFCCLSn

        x1, x2 = sym.nsolve((eqn1, eqn2), (XS, XL), (0.1, 0.5))
        row = pd.DataFrame(
            [[temp, np.nan, np.nan, x1, x2, np.nan, np.nan]],
            columns=COLS,
        )
        row = row.set_index("Temperature")
        phase_diagram = pd.concat([phase_diagram, row])

    for temp in np.arange(446, 508, 1):
        GFCCLPb = 4810 - 8.017 * temp
        GFCCBCTPb = 489 + 3.52 * temp

        GBCTLSn = 7179 - 14.126 * temp
        GBCTLPb = -1 * GFCCBCTPb + GFCCLPb

        logBCTLPb = 8.314 * temp * sym.log((1 - XL) / (1 - XS))
        logBCTLSn = 8.314 * temp * sym.log(XL / XS)

        RKBCTLPb = 4700 * XL ** 2 - 12000 * XS ** 2
        RKBCTLSn = 4700 * (1 - XL) ** 2 - 12000 * (1 - XS) ** 2

        eqn1 = GBCTLPb + logBCTLPb + RKBCTLPb
        eqn2 = GBCTLSn + logBCTLSn + RKBCTLSn

        x2, x1 = sym.nsolve((eqn1, eqn2), (XS, XL), (0.9, 0.5))
        phase_diagram.loc[temp, ["betaLiquid_L", "betaLiquid_R"]] = [x1, x2]

    phase_diagram = phase_diagram.sort_index()
    phase_diagram = phase_diagram.reset_index()
    phase_diagram["Temperature"] = phase_diagram["Temperature"].astype(int)
    phase_diagram = phase_diagram.set_index("Temperature")
    return phase_diagram


phaseDiagram = build_phase_diagram()
