import numpy as np
import matplotlib.pyplot as plt
from itertools import product

time = np.arange(0, 40 * np.pi) / 10

# Number of bits per QAM symbol
n = 4
# 'grid' shows IQ waveforms in grid, 'rows' shows IQ | Sum in rows
mode = 'rows'


def map_iq(i_bits: list[int], q_bits: list[int]) -> complex:
    ''' takes in-phase and quadrature bits to generate a complex IQ symbol '''
    i = 0
    for bit in i_bits:
        i = (i << 1) | bit
    q = 0
    for bit in q_bits:
        q = (q << 1) | bit

    i -= (2**len(i_bits)-1)/2
    q -= (2**len(q_bits)-1)/2

    return i*2 + q*2j


# iterate over all combinations of I and Q bits to generate symbols
# n//2 bits are assigned to Q and the rest to I (equal exept in case of odd n)
iq_symbols: list[complex] = []
for i_bits in product((0, 1), repeat=n-n//2):
    for q_bits in product((0, 1), repeat=n//2):
        iq_symbols.append(map_iq(i_bits, q_bits))


# turn IQ symbols into sinusiods, real = cos, imag = sin
iq_waves: list[np.complex128] = []
for iq in iq_symbols:
    iq_waves.append(
        iq.real * np.cos(time) + iq.imag * np.sin(time) * 1j
    )

if mode == 'rows':
    # set up IQ and Sum plot axis for side-by-side plotting
    num_rows = 2**n
    ax_iq: list[plt.Axes] = [plt.subplot2grid((num_rows, 2), (0, 0))]
    ax_sum: list[plt.Axes] = [plt.subplot2grid((num_rows, 2), (0, 1))]
    for x in range(1, num_rows):
        ax_iq.append(plt.subplot2grid((num_rows, 2), (x, 0)))
        ax_sum.append(plt.subplot2grid((num_rows, 2), (x, 1)))
        ax_iq[-1].sharey(ax_iq[0])
        ax_sum[-1].sharey(ax_sum[0])
elif mode == 'grid':
    # set up IQ axis for grid plotting
    num_rows = 2**(n // 2)
    num_cols = 2**(n - (n // 2))
    ax_iq: list[plt.Axes] = [plt.subplot2grid((num_rows, num_cols), (0, 0))]
    ax_sum = None
    for x in range(num_cols):
        for y in range(num_rows):
            if x == 0 and y == 0:
                continue
            ax_iq.append(plt.subplot2grid((num_rows, num_cols), (y, x)))
            ax_iq[-1].sharey(ax_iq[0])
else:
    raise ValueError("Error, invalid value for mode. Must be either 'rows' or 'grid'")


# plot all I and Q waves
if mode == 'rows':
    for ax1, ax2, iq in zip(ax_iq, ax_sum, iq_waves):
        ax1.plot(time, iq.real)
        ax1.plot(time, iq.imag)
        ax2.plot(time, iq.real+iq.imag)
elif mode == 'grid':
    for ax1, iq in zip(ax_iq, iq_waves):
        ax1.plot(time, iq.real)
        ax1.plot(time, iq.imag)
else:
    raise ValueError("Error, invalid value for mode. Must be either 'rows' or 'grid'")


plt.figure()
# setup constellation diagram axis
ax_const = plt.subplot()

# plot constellation diagram
ax_const.scatter([iq.real for iq in iq_symbols], [iq.imag for iq in iq_symbols])
ax_const.axhline(y=0, color='k', linewidth=.5)
ax_const.axvline(x=0, color='k', linewidth=.5)
ax_const.set_aspect('equal', 'box')

ax_const.set_title(f"QAM - {n} bits per symbol")


plt.show()
