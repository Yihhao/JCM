import matplotlib.pyplot as plt
import numpy as np
from numpy import pi


def plot_energy_b(delta, evals_mat, g, wc_f, plot_range, plot_line,
                  fig=None, ax=None, figsize=(8, 12)):
    color = ['b', 'r', 'k']
    if not fig and not ax:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    for n in np.arange(plot_line):
        if n == 0:
            c = color[2]
        else:
            c = color[n % 2]
        ax.plot(delta / g, (evals_mat[:, n] - evals_mat[:, 0]) / (2 * pi), c)

    ax.set_xlabel('$\Delta$ /g ', fontsize=16)
    plt.xlim(-1 * plot_range, plot_range)

    y_ticks = []
    y_position = []
    for i in np.arange(plot_line//2 + 1):
        y_position.append(i * wc_f)
        if i == 0:
            y_ticks.append(0)
        else:
            y_ticks.append(r'$%s\omega_c$' % i)
    plt.yticks(y_position, y_ticks, fontsize=16)
    return fig, ax


def plot_energy_d(evals_mat, delta, wc, g, plot_range,
                  fig=None, ax=None, figsize=(6, 6)):
    if not fig and not ax:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    text = ['first rung', 'second rung', 'third rung']
    style = ['-', '--', ':']
    for n in [1, 2, 3]:
        for i in [2 * n - 1, 2 * n]:
            c = 'k'
            s = style[n-1]
            y = (evals_mat[:, i] - evals_mat[:, 0] - n * wc) / g
            if i % 2 == 0:
                ax.plot(delta / g, y, c + s, label=text[n - 1])
            else:
                ax.plot(delta / g, y, c + s)
    plt.ylim(-4, 4)
    plt.xlim(-1 * plot_range, plot_range)
    plt.ylabel(r'$\omega / g$')
    plt.xlabel(r'$\Delta / g$')

    return fig, ax


def plot_energy_e(evals_mat, walist, wc, g, delta, plot_range, plot_line,
                  fig=None, ax=None, figsize=(6, 6)):
    if not fig and not ax:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    rung = np.zeros([plot_line, len(walist)])
    for n in np.arange(plot_line):
        rung[n, :] = (evals_mat[:, n] - evals_mat[:, 0])
    c = ['g', 'b']
    s = ['-', '--', ':']
    for n in [1, 2, 3]:
        for i in [2 * n - 1, 2 * n]:
            style = s[n - 1]
            if n == 1:
                style = c[i // 2] + style
                y = (rung[i] - rung[0] - wc) / g
                ax.plot(delta / g, y, style, linewidth=3, label='0 -> %s' % i)
                # ax.plot(delta / g, y, label='0 -> %s' % i)
            else:
                pre_n = (n - 1)
                y1 = (rung[i] - rung[2 * pre_n - 1] - wc) / g
                ax.plot(delta / g, y1, 'g' + style, label='%s -> %s' % (2 * pre_n - 1, i))
                y2 = (rung[i] - rung[2 * pre_n] - wc) / g
                ax.plot(delta / g, y2, 'b' + style, label='%s -> %s' % (2 * pre_n, i))


    plt.xlim(-1 * plot_range, plot_range)
    plt.ylim(-4, 4)
    plt.xlabel(r'$\Delta /g$')
    plt.ylabel(r'$\omega / g$')
    return fig, ax