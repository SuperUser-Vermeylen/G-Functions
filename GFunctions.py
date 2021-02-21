import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import numpy as np
from PhaseDiagram import phaseDiagram

precision = 3
delta = 10 ** -precision
X = np.arange(delta, 1, delta)


def chem_potential(temperature, G_Total_FCC, G_Total_BCT, G_Total_L):

    if temperature <= 445:

        x1 = np.around(phaseDiagram.loc[[temperature]]["alphaBeta_L"].values.astype(np.double), precision)
        x2 = np.around(phaseDiagram.loc[[temperature]]["alphaBeta_R"].values.astype(np.double), precision)

        index1 = np.where(np.isclose(X, x1))
        index2 = np.where(np.isclose(X, x2))

        Mu1 = (G_Total_FCC[index1] - G_Total_BCT[index2])/(X[index1] - X[index2]) * (X - X[index1]) + G_Total_FCC[index1]
        Mu2 = np.empty((len(X),))
        Mu2[:] = np.NaN

        alphaBeta_L = index1
        alphaBeta_R = index2
        alphaLiquid_L = np.nan
        alphaLiquid_R = np.nan
        betaLiquid_L = np.nan
        betaLiquid_R = np.nan

        return [Mu1, Mu2, alphaBeta_L, alphaBeta_R, alphaLiquid_L, alphaLiquid_R, betaLiquid_L, betaLiquid_R]

    elif 445 < temperature <= 507:

        x1 = np.around(phaseDiagram.loc[[temperature]]["alphaLiquid_L"].values.astype(np.double), precision)
        x2 = np.around(phaseDiagram.loc[[temperature]]["alphaLiquid_R"].values.astype(np.double), precision)
        x3 = np.around(phaseDiagram.loc[[temperature]]["betaLiquid_L"].values.astype(np.double), precision)
        x4 = np.around(phaseDiagram.loc[[temperature]]["betaLiquid_R"].values.astype(np.double), precision)

        index1 = np.where(np.isclose(X, x1))
        index2 = np.where(np.isclose(X, x2))
        index3 = np.where(np.isclose(X, x3))
        index4 = np.where(np.isclose(X, x4))

        Mu1 = (G_Total_FCC[index1] - G_Total_L[index2])/(X[index1] - X[index2]) * (X - X[index1]) + G_Total_FCC[index1]
        Mu2 = (G_Total_BCT[index4] - G_Total_L[index3])/(X[index4] - X[index3]) * (X - X[index3]) + G_Total_L[index3]

        alphaBeta_L = np.nan
        alphaBeta_R = np.nan
        alphaLiquid_L = index1
        alphaLiquid_R = index2
        betaLiquid_L = index3
        betaLiquid_R = index4

        return [Mu1, Mu2, alphaBeta_L, alphaBeta_R, alphaLiquid_L, alphaLiquid_R, betaLiquid_L, betaLiquid_R]

    elif 507 < temperature:

        x1 = np.around(phaseDiagram.loc[[temperature]]["alphaLiquid_L"].values.astype(np.double), precision)
        x2 = np.around(phaseDiagram.loc[[temperature]]["alphaLiquid_R"].values.astype(np.double), precision)

        index1 = np.where(np.isclose(X, x1))
        index2 = np.where(np.isclose(X, x2))

        Mu1 = (G_Total_FCC[index1] - G_Total_L[index2])/(X[index1] - X[index2]) * (X - X[index1]) + G_Total_FCC[index1]
        Mu2 = np.empty(len(X))
        Mu2[:] = np.NaN

        alphaBeta_L = np.nan
        alphaBeta_R = np.nan
        alphaLiquid_L = index1
        alphaLiquid_R = index2
        betaLiquid_L = np.nan
        betaLiquid_R = np.nan

        return [Mu1, Mu2, alphaBeta_L, alphaBeta_R, alphaLiquid_L, alphaLiquid_R, betaLiquid_L, betaLiquid_R]


def gfunctions(temperature):
    flag = 1
    while flag:
        if 300 <= temperature <= 445:

            G_A_fcc = 0
            G_A_L = 4810 - 8.017 * temperature
            G_A_bct = 489 + 3.52 * temperature
            G_B_bct = 0
            G_B_L = 7179 - 14.216 * temperature
            G_B_fcc = 5510 - 8.46 * temperature
            flag = 0

        elif 445 < temperature <= 505:

            G_A_fcc = 0
            G_A_L = 4810 - 8.017 * temperature
            G_A_bct = 489 + 3.52 * temperature
            G_B_bct = 0
            G_B_L = 7179 - 14.216 * temperature
            G_B_fcc = 5510 - 8.46 * temperature
            flag = 0

        elif 505 < temperature <= 599:

            G_A_fcc = 0
            G_A_L = 4810 - 8.017 * temperature
            G_A_bct = 489 + 3.52 * temperature
            G_B_bct = -7179 + 14.216 * temperature
            G_B_L = 0
            G_B_fcc = -1669 + 5.756 * temperature
            flag = 0

    G_Total_FCC = (1 - X) * G_A_fcc + X * G_B_fcc + 8.314 * temperature * ((1 - X) * np.log(1 - X) + X * np.log(X)) + 5200 * (1 - X) * X
    G_Total_BCT = (1 - X) * G_A_bct + X * G_B_bct + 8.314 * temperature * ((1 - X) * np.log(1 - X) + X * np.log(X)) + 12000 * (1 - X) * X
    G_Total_L = (1 - X) * G_A_L + X * G_B_L + 8.314 * temperature * ((1 - X) * np.log(1 - X) + X * np.log(X)) + 4700 * (1 - X) * X

    mu = chem_potential(temperature, G_Total_FCC, G_Total_BCT, G_Total_L)

    return [G_Total_FCC, G_Total_BCT, G_Total_L, mu[0], mu[1], mu[2], mu[3], mu[4], mu[5], mu[6], mu[7]]


def update(val):
    newFunctions = gfunctions(slider_temperature.val)
    A.set_ydata(newFunctions[0])
    B.set_ydata(newFunctions[1])
    L.set_ydata(newFunctions[2])
    mu1.set_ydata(newFunctions[3])
    mu2.set_ydata(newFunctions[4])

    if val <= 445:
        alphaBeta_L.set_xdata([X[newFunctions[5]]] * 2)
        alphaBeta_L.set_ydata([-3000, newFunctions[3][newFunctions[5]]])
        alphaBeta_R.set_xdata([X[newFunctions[6]]] * 2)
        alphaBeta_R.set_ydata([-3000, newFunctions[3][newFunctions[6]]])
        alphaLiquid_L.set_xdata([0] * 2)
        alphaLiquid_L.set_ydata([0] * 2)
        alphaLiquid_R.set_xdata([0] * 2)
        alphaLiquid_R.set_ydata([0] * 2)
        betaLiquid_L.set_xdata([0] * 2)
        betaLiquid_L.set_ydata([0] * 2)
        betaLiquid_R.set_xdata([0] * 2)
        betaLiquid_R.set_ydata([0] * 2)

        texts1[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (np.around(newFunctions[3][0], 0),),
                                         (r'$\mu_{\beta}^{A}=%.2f$' % (np.around(newFunctions[3][-1], 0),)))))
        texts2[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{B}=%.2f$' % (np.around(newFunctions[3][0], 0),),
                                         (r'$\mu_{\beta}^{B}=%.2f$' % (np.around(newFunctions[3][-1], 0),)))))

    elif 445 < val <= 507:

        alphaBeta_L.set_xdata([0] * 2)
        alphaBeta_L.set_ydata([0] * 2)
        alphaBeta_R.set_xdata([0] * 2)
        alphaBeta_R.set_ydata([0] * 2)
        alphaLiquid_L.set_xdata([X[newFunctions[7]]] * 2)
        alphaLiquid_L.set_ydata([-3000, newFunctions[3][newFunctions[7]]])
        alphaLiquid_R.set_xdata([X[newFunctions[8]]] * 2)
        alphaLiquid_R.set_ydata([-3000, newFunctions[3][newFunctions[8]]])
        betaLiquid_L.set_xdata([X[newFunctions[9]]] * 2)
        betaLiquid_L.set_ydata([-3000, newFunctions[4][newFunctions[9]]])
        betaLiquid_R.set_xdata([X[newFunctions[10]]] * 2)
        betaLiquid_R.set_ydata([-3000, newFunctions[4][newFunctions[10]]])

        texts1[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (np.around(newFunctions[3][0], 0),),
                                         (r'$\mu_{\ell}^{A}=%.2f$' % (np.around(newFunctions[3][-1], 0),)))))
        texts2[0][0].set_text('\n'.join((r'$\mu_{\ell}^{B}=%.2f$' % (np.around(newFunctions[4][0], 0),),
                                         (r'$\mu_{\beta}^{B}=%.2f$' % (np.around(newFunctions[4][-1], 0),)))))

    elif 507 < val:

        alphaBeta_L.set_xdata([0] * 2)
        alphaBeta_L.set_ydata([0] * 2)
        alphaBeta_R.set_xdata([0] * 2)
        alphaBeta_R.set_ydata([0] * 2)
        alphaLiquid_L.set_xdata([X[newFunctions[7]]] * 2)
        alphaLiquid_L.set_ydata([-3000, newFunctions[3][newFunctions[7]]])
        alphaLiquid_R.set_xdata([X[newFunctions[8]]] * 2)
        alphaLiquid_R.set_ydata([-3000, newFunctions[3][newFunctions[8]]])
        betaLiquid_L.set_xdata([0] * 2)
        betaLiquid_L.set_ydata([0] * 2)
        betaLiquid_R.set_xdata([0] * 2)
        betaLiquid_R.set_ydata([0] * 2)

        texts1[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (np.around(newFunctions[3][0], 0),),
                                         (r'$\mu_{\ell}^{A}=%.2f$' % (np.around(newFunctions[3][-1], 0),)))))
        texts2[0][0].set_text('\n'.join((r'$\mu_{\ell}^{B}=%.2f$' % (np.around(newFunctions[4][0], 0),),
                                         (r'$\mu_{\beta}^{B}=%.2f$' % (np.around(newFunctions[4][-1], 0),)))))

    df_index = int(val - 300)
    phase_Data = diagramData.iloc[0:df_index, :]
    df_Temp = phase_Data.Temperature
    df_alphaBeta_L = phase_Data.alphaBeta_L
    df_alphaBeta_R = phase_Data.alphaBeta_R
    df_alphaLiquid_L = phase_Data.alphaLiquid_L
    df_alphaLiquid_R = phase_Data.alphaLiquid_R
    df_betaLiquid_L = phase_Data.betaLiquid_L
    df_betaLiquid_R = phase_Data.betaLiquid_R

    s1.set_xdata(df_alphaBeta_L)
    s1.set_ydata(df_Temp)
    s2.set_xdata(df_alphaBeta_R)
    s2.set_ydata(df_Temp)
    s3.set_xdata(df_alphaLiquid_L)
    s3.set_ydata(df_Temp)
    s4.set_xdata(df_alphaLiquid_R)
    s4.set_ydata(df_Temp)
    s5.set_xdata(df_betaLiquid_L)
    s5.set_ydata(df_Temp)
    s6.set_xdata(df_betaLiquid_R)
    s6.set_ydata(df_Temp)

    fig1.canvas.draw_idle()
    fig2.canvas.draw_idle()


fig1, ax1 = plt.subplots()
ax1.set_ylim([-3000, 3000])
fig2, ax2 = plt.subplots()
ax2.set_ylim([300, 600])
ax2.set_xlim([0, 1])

plt.sca(ax2)
diagramData = phaseDiagram.reset_index()
phase_Data = diagramData.iloc[0:145, :]
s1, = ax2.plot(phase_Data.alphaBeta_L, phase_Data.Temperature, 'k*')
s2, = ax2.plot(phase_Data.alphaBeta_R, phase_Data.Temperature, 'k*')
s3, = ax2.plot(phase_Data.alphaLiquid_L, phase_Data.Temperature, 'k*')
s4, = ax2.plot(phase_Data.alphaLiquid_R, phase_Data.Temperature, 'k*')
s5, = ax2.plot(phase_Data.betaLiquid_L, phase_Data.Temperature, 'k*')
s6, = ax2.plot(phase_Data.betaLiquid_R, phase_Data.Temperature, 'k*')

plt.title('Phase Diagram of the Binary A-B System')
plt.xlabel('$X_{B}$')

plt.sca(ax1)
Y = gfunctions(445)
A, = plt.plot(X, Y[0], lw=2)
B, = plt.plot(X, Y[1], lw=2)
L, = plt.plot(X, Y[2], lw=2)
mu1, = plt.plot(X, Y[3], lw=2)
mu2, = plt.plot(X, Y[4], lw=2)

alphaBeta_L, = plt.plot([X[Y[5]]]*2, [-3000, Y[3][Y[5]]], 'k--', lw=2)
alphaBeta_R, = plt.plot([X[Y[6]]]*2, [-3000, Y[3][Y[6]]], 'k--', lw=2)
alphaLiquid_L, = plt.plot([X[Y[5]]]*2, [-3000, Y[3][Y[5]]], 'k--', lw=2)
alphaLiquid_R, = plt.plot([X[Y[6]]]*2, [-3000, Y[3][Y[6]]], 'k--', lw=2)
betaLiquid_L, = plt.plot([X[Y[5]]]*2, [-3000, Y[3][Y[5]]], 'k--', lw=2)
betaLiquid_R, = plt.plot([X[Y[6]]]*2, [-3000, Y[3][Y[6]]], 'k--', lw=2)

plt.subplots_adjust(left=0.05, right=0.8, bottom=0.1, top=0.9)
plt.axhline(y=0, color='black', linestyle='-')
plt.title('Gibbs Energy of Phases in the Binary A-B System')
plt.xlabel('$X_{B}$')
plt.yticks([0.0])

textstr1 = '\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (np.nan,), (r'$\mu_{\beta}^{A}=%.2f$' % (np.nan,))))
textstr2 = '\n'.join((r'$\mu_{\alpha}^{B}=%.2f$' % (np.nan,), (r'$\mu_{\beta}^{B}=%.2f$' % (np.nan,))))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
texts1 = [[plt.text(0.83, 0.50, textstr1, fontsize=14, transform=plt.gcf().transFigure, bbox=props)]]
texts2 = [[plt.text(0.83, 0.25, textstr2, fontsize=14, transform=plt.gcf().transFigure, bbox=props)]]


ax1.margins(x=0)
plt.legend([A, B, L, mu1, mu2], [r'${\alpha}$', r'${\beta}$', 'L', r'$\mu^A_{\alpha,\beta}$', r'$\mu^B_{\alpha,\beta}$']
           , bbox_to_anchor=(1.04, 1), loc="upper left")

ax_color = 'lightgreen'
temp_ax = plt.axes([0.02, 0.1, 0.01, 0.8], facecolor=ax_color)
slider_temperature = Slider(temp_ax, 'T/K', 300, 599, valinit=445, valstep=1, orientation='vertical', valfmt='%0.0f')
slider_temperature.on_changed(update)

plt.show()
