from __future__ import annotations

import numpy as np

from .phase_diagram import build_phase_diagram

PRECISION = 3
DELTA = 10 ** -PRECISION
X = np.arange(DELTA, 1, DELTA)


def chem_potential(temperature: float, G_Total_FCC, G_Total_BCT, G_Total_L, phase_diagram=None):
    """Compute the chemical potential curves for the current temperature."""
    if phase_diagram is None:
        phase_diagram = build_phase_diagram()

    if temperature <= 445:
        x1 = np.around(phase_diagram.loc[[temperature]]["alphaBeta_L"].values.astype(np.double), PRECISION)
        x2 = np.around(phase_diagram.loc[[temperature]]["alphaBeta_R"].values.astype(np.double), PRECISION)

        index1 = np.where(np.isclose(X, x1))
        index2 = np.where(np.isclose(X, x2))

        Mu1 = (G_Total_FCC[index1] - G_Total_BCT[index2]) / (X[index1] - X[index2]) * (X - X[index1]) + G_Total_FCC[index1]
        Mu2 = np.empty((len(X),))
        Mu2[:] = np.nan

        return [Mu1, Mu2, index1, index2, np.nan, np.nan, np.nan, np.nan]

    if 445 < temperature <= 507:
        x1 = np.around(phase_diagram.loc[[temperature]]["alphaLiquid_L"].values.astype(np.double), PRECISION)
        x2 = np.around(phase_diagram.loc[[temperature]]["alphaLiquid_R"].values.astype(np.double), PRECISION)
        x3 = np.around(phase_diagram.loc[[temperature]]["betaLiquid_L"].values.astype(np.double), PRECISION)
        x4 = np.around(phase_diagram.loc[[temperature]]["betaLiquid_R"].values.astype(np.double), PRECISION)

        index1 = np.where(np.isclose(X, x1))
        index2 = np.where(np.isclose(X, x2))
        index3 = np.where(np.isclose(X, x3))
        index4 = np.where(np.isclose(X, x4))

        Mu1 = (G_Total_FCC[index1] - G_Total_L[index2]) / (X[index1] - X[index2]) * (X - X[index1]) + G_Total_FCC[index1]
        Mu2 = (G_Total_BCT[index4] - G_Total_L[index3]) / (X[index4] - X[index3]) * (X - X[index3]) + G_Total_L[index3]

        return [Mu1, Mu2, index1, index2, index1, index2, index3, index4]

    x1 = np.around(phase_diagram.loc[[temperature]]["alphaLiquid_L"].values.astype(np.double), PRECISION)
    x2 = np.around(phase_diagram.loc[[temperature]]["alphaLiquid_R"].values.astype(np.double), PRECISION)

    index1 = np.where(np.isclose(X, x1))
    index2 = np.where(np.isclose(X, x2))

    Mu1 = (G_Total_FCC[index1] - G_Total_L[index2]) / (X[index1] - X[index2]) * (X - X[index1]) + G_Total_FCC[index1]
    Mu2 = np.empty(len(X))
    Mu2[:] = np.nan

    return [Mu1, Mu2, np.nan, np.nan, index1, index2, np.nan, np.nan]


def _phase_coefficients(temperature: float):
    if 300 <= temperature <= 445:
        return {
            "G_A_fcc": 0,
            "G_A_L": 4810 - 8.017 * temperature,
            "G_A_bct": 489 + 3.52 * temperature,
            "G_B_bct": 0,
            "G_B_L": 7179 - 14.216 * temperature,
            "G_B_fcc": 5510 - 8.46 * temperature,
        }

    if 445 < temperature <= 505:
        return {
            "G_A_fcc": 0,
            "G_A_L": 4810 - 8.017 * temperature,
            "G_A_bct": 489 + 3.52 * temperature,
            "G_B_bct": 0,
            "G_B_L": 7179 - 14.216 * temperature,
            "G_B_fcc": 5510 - 8.46 * temperature,
        }

    if 505 < temperature <= 599:
        return {
            "G_A_fcc": 0,
            "G_A_L": 4810 - 8.017 * temperature,
            "G_A_bct": 489 + 3.52 * temperature,
            "G_B_bct": -7179 + 14.216 * temperature,
            "G_B_L": 0,
            "G_B_fcc": -1669 + 5.756 * temperature,
        }

    raise ValueError("temperature must be between 300 and 599 K")


def gfunctions(temperature: float, phase_diagram=None):
    """Compute Gibbs energy curves and associated chemical potential data."""
    if phase_diagram is None:
        phase_diagram = build_phase_diagram()

    coeffs = _phase_coefficients(temperature)

    G_Total_FCC = (
        (1 - X) * coeffs["G_A_fcc"]
        + X * coeffs["G_B_fcc"]
        + 8.314 * temperature * ((1 - X) * np.log(1 - X) + X * np.log(X))
        + 5200 * (1 - X) * X
    )
    G_Total_BCT = (
        (1 - X) * coeffs["G_A_bct"]
        + X * coeffs["G_B_bct"]
        + 8.314 * temperature * ((1 - X) * np.log(1 - X) + X * np.log(X))
        + 12000 * (1 - X) * X
    )
    G_Total_L = (
        (1 - X) * coeffs["G_A_L"]
        + X * coeffs["G_B_L"]
        + 8.314 * temperature * ((1 - X) * np.log(1 - X) + X * np.log(X))
        + 4700 * (1 - X) * X
    )

    mu = chem_potential(temperature, G_Total_FCC, G_Total_BCT, G_Total_L, phase_diagram=phase_diagram)

    return [
        G_Total_FCC,
        G_Total_BCT,
        G_Total_L,
        mu[0],
        mu[1],
        mu[2],
        mu[3],
        mu[4],
        mu[5],
        mu[6],
        mu[7],
    ]
