from __future__ import annotations

from pathlib import Path

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

from .g_functions import X, gfunctions
from .phase_diagram import build_phase_diagram


def create_interactive_plot() -> None:
    """Create the interactive Gibbs-energy and phase-diagram plot."""
    phase_diagram = build_phase_diagram()
    diagram_data = phase_diagram.reset_index()

    fig1, ax1 = plt.subplots()
    ax1.set_ylim([-3000, 3000])
    fig2, ax2 = plt.subplots()
    ax2.set_ylim([300, 600])
    ax2.set_xlim([0, 1])

    plt.sca(ax2)
    phase_data = diagram_data.iloc[0:145, :]
    s1, = ax2.plot(phase_data.alphaBeta_L, phase_data.Temperature, 'k*')
    s2, = ax2.plot(phase_data.alphaBeta_R, phase_data.Temperature, 'k*')
    s3, = ax2.plot(phase_data.alphaLiquid_L, phase_data.Temperature, 'k*')
    s4, = ax2.plot(phase_data.alphaLiquid_R, phase_data.Temperature, 'k*')
    s5, = ax2.plot(phase_data.betaLiquid_L, phase_data.Temperature, 'k*')
    s6, = ax2.plot(phase_data.betaLiquid_R, phase_data.Temperature, 'k*')

    plt.title('Phase Diagram of the Binary A-B System')
    plt.xlabel('$X_{B}$')

    plt.sca(ax1)
    Y = gfunctions(445, phase_diagram)
    A, = plt.plot(X, Y[0], lw=2)
    B, = plt.plot(X, Y[1], lw=2)
    L, = plt.plot(X, Y[2], lw=2)
    mu1, = plt.plot(X, Y[3], lw=2)
    mu2, = plt.plot(X, Y[4], lw=2)

    alphaBeta_L, = plt.plot([X[Y[5]]] * 2, [-3000, Y[3][Y[5]]], 'k--', lw=2)
    alphaBeta_R, = plt.plot([X[Y[6]]] * 2, [-3000, Y[3][Y[6]]], 'k--', lw=2)
    alphaLiquid_L, = plt.plot([X[Y[5]]] * 2, [-3000, Y[3][Y[5]]], 'k--', lw=2)
    alphaLiquid_R, = plt.plot([X[Y[6]]] * 2, [-3000, Y[3][Y[6]]], 'k--', lw=2)
    betaLiquid_L, = plt.plot([X[Y[5]]] * 2, [-3000, Y[3][Y[5]]], 'k--', lw=2)
    betaLiquid_R, = plt.plot([X[Y[6]]] * 2, [-3000, Y[3][Y[6]]], 'k--', lw=2)

    plt.subplots_adjust(left=0.05, right=0.8, bottom=0.1, top=0.9)
    plt.axhline(y=0, color='black', linestyle='-')
    plt.title('Gibbs Energy of Phases in the Binary A-B System')
    plt.xlabel('$X_{B}$')
    plt.yticks([0.0])

    textstr1 = '\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (float('nan'),), (r'$\mu_{\beta}^{A}=%.2f$' % (float('nan'),))))
    textstr2 = '\n'.join((r'$\mu_{\alpha}^{B}=%.2f$' % (float('nan'),), (r'$\mu_{\beta}^{B}=%.2f$' % (float('nan'),))))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    texts1 = [[plt.text(0.83, 0.50, textstr1, fontsize=14, transform=plt.gcf().transFigure, bbox=props)]]
    texts2 = [[plt.text(0.83, 0.25, textstr2, fontsize=14, transform=plt.gcf().transFigure, bbox=props)]]

    ax1.margins(x=0)
    plt.legend(
        [A, B, L, mu1, mu2],
        [r'${\alpha}$', r'${\beta}$', 'L', r'$\mu^A_{\alpha,\beta}$', r'$\mu^B_{\alpha,\beta}$'],
        bbox_to_anchor=(1.04, 1),
        loc='upper left',
    )

    def update(val):
        new_functions = gfunctions(slider_temperature.val, phase_diagram)
        A.set_ydata(new_functions[0])
        B.set_ydata(new_functions[1])
        L.set_ydata(new_functions[2])
        mu1.set_ydata(new_functions[3])
        mu2.set_ydata(new_functions[4])

        if val <= 445:
            alphaBeta_L.set_xdata([X[new_functions[5]]] * 2)
            alphaBeta_L.set_ydata([-3000, new_functions[3][new_functions[5]]])
            alphaBeta_R.set_xdata([X[new_functions[6]]] * 2)
            alphaBeta_R.set_ydata([-3000, new_functions[3][new_functions[6]]])
            alphaLiquid_L.set_xdata([0] * 2)
            alphaLiquid_L.set_ydata([0] * 2)
            alphaLiquid_R.set_xdata([0] * 2)
            alphaLiquid_R.set_ydata([0] * 2)
            betaLiquid_L.set_xdata([0] * 2)
            betaLiquid_L.set_ydata([0] * 2)
            betaLiquid_R.set_xdata([0] * 2)
            betaLiquid_R.set_ydata([0] * 2)

            texts1[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (float(np.around(new_functions[3][0], 0)),), (r'$\mu_{\beta}^{A}=%.2f$' % (float(np.around(new_functions[3][-1], 0)),))))
            texts2[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{B}=%.2f$' % (float(np.around(new_functions[3][0], 0)),), (r'$\mu_{\beta}^{B}=%.2f$' % (float(np.around(new_functions[3][-1], 0)),))))
        elif 445 < val <= 507:
            alphaBeta_L.set_xdata([0] * 2)
            alphaBeta_L.set_ydata([0] * 2)
            alphaBeta_R.set_xdata([0] * 2)
            alphaBeta_R.set_ydata([0] * 2)
            alphaLiquid_L.set_xdata([X[new_functions[7]]] * 2)
            alphaLiquid_L.set_ydata([-3000, new_functions[3][new_functions[7]]])
            alphaLiquid_R.set_xdata([X[new_functions[8]]] * 2)
            alphaLiquid_R.set_ydata([-3000, new_functions[3][new_functions[8]]])
            betaLiquid_L.set_xdata([X[new_functions[9]]] * 2)
            betaLiquid_L.set_ydata([-3000, new_functions[4][new_functions[9]]])
            betaLiquid_R.set_xdata([X[new_functions[10]]] * 2)
            betaLiquid_R.set_ydata([-3000, new_functions[4][new_functions[10]]])

            texts1[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (float(np.around(new_functions[3][0], 0)),), (r'$\mu_{\ell}^{A}=%.2f$' % (float(np.around(new_functions[3][-1], 0)),))))
            texts2[0][0].set_text('\n'.join((r'$\mu_{\ell}^{B}=%.2f$' % (float(np.around(new_functions[4][0], 0)),), (r'$\mu_{\beta}^{B}=%.2f$' % (float(np.around(new_functions[4][-1], 0)),))))
        elif 507 < val:
            alphaBeta_L.set_xdata([0] * 2)
            alphaBeta_L.set_ydata([0] * 2)
            alphaBeta_R.set_xdata([0] * 2)
            alphaBeta_R.set_ydata([0] * 2)
            alphaLiquid_L.set_xdata([X[new_functions[7]]] * 2)
            alphaLiquid_L.set_ydata([-3000, new_functions[3][new_functions[7]]])
            alphaLiquid_R.set_xdata([X[new_functions[8]]] * 2)
            alphaLiquid_R.set_ydata([-3000, new_functions[3][new_functions[8]]])
            betaLiquid_L.set_xdata([0] * 2)
            betaLiquid_L.set_ydata([0] * 2)
            betaLiquid_R.set_xdata([0] * 2)
            betaLiquid_R.set_ydata([0] * 2)

            texts1[0][0].set_text('\n'.join((r'$\mu_{\alpha}^{A}=%.2f$' % (float(np.around(new_functions[3][0], 0)),), (r'$\mu_{\ell}^{A}=%.2f$' % (float(np.around(new_functions[3][-1], 0)),))))
            texts2[0][0].set_text('\n'.join((r'$\mu_{\ell}^{B}=%.2f$' % (float(np.around(new_functions[4][0], 0)),), (r'$\mu_{\beta}^{B}=%.2f$' % (float(np.around(new_functions[4][-1], 0)),))))

        phase_slice = diagram_data.iloc[0:int(val - 300), :]
        df_temp = phase_slice.Temperature
        df_alphaBeta_L = phase_slice.alphaBeta_L
        df_alphaBeta_R = phase_slice.alphaBeta_R
        df_alphaLiquid_L = phase_slice.alphaLiquid_L
        df_alphaLiquid_R = phase_slice.alphaLiquid_R
        df_betaLiquid_L = phase_slice.betaLiquid_L
        df_betaLiquid_R = phase_slice.betaLiquid_R

        s1.set_xdata(df_alphaBeta_L)
        s1.set_ydata(df_temp)
        s2.set_xdata(df_alphaBeta_R)
        s2.set_ydata(df_temp)
        s3.set_xdata(df_alphaLiquid_L)
        s3.set_ydata(df_temp)
        s4.set_xdata(df_alphaLiquid_R)
        s4.set_ydata(df_temp)
        s5.set_xdata(df_betaLiquid_L)
        s5.set_ydata(df_temp)
        s6.set_xdata(df_betaLiquid_R)
        s6.set_ydata(df_temp)

        fig1.canvas.draw_idle()
        fig2.canvas.draw_idle()

    ax_color = 'lightgreen'
    temp_ax = plt.axes([0.02, 0.1, 0.01, 0.8], facecolor=ax_color)
    slider_temperature = Slider(temp_ax, 'T/K', 300, 599, valinit=445, valstep=1, orientation='vertical', valfmt='%0.0f')
    slider_temperature.on_changed(update)

    plt.show()


def main() -> None:
    create_interactive_plot()


if __name__ == "__main__":
    main()
