import numpy as np
import matplotlib.pyplot as plt
from itertools import product

time = np.arange(0, 100 * np.pi)

# Number of bits per QAM symbol
n = 4
mode = 'grid'

def split_bits(bits: list[int]) -> tuple[list[int]]:
    ''' split list of bits in half (I-Q) '''
    return bits[len(bits)//2:], bits[:len(bits)//2]

def map_iq(i_bits: list[int], q_bits: list[int]) -> tuple[float]:
    ''' take in-phase and quadrature bits and generate I-Q symbol '''
    i = 0
    for bit in i_bits:
        i = (i << 1) | bit
    q = 0
    for bit in q_bits:
        q = (q << 1) | bit

    i -= (2**len(i_bits)-1)/2
    q -= (2**len(q_bits)-1)/2

    return i*2, q*2


# loop over all possible bit combinations i.e. (0,0),(0,1),(1,0),(1,1)
# and generate lists of IQ symbols and IQ waveforms
i_symbols, q_symbols = [], []
i_waves, q_waves = [], []
for a in product((0, 1), repeat=n):
    I, Q = map_iq(*split_bits(a))
    i_symbols.append(I)
    q_symbols.append(Q)
    
    i_waves.append(I * np.cos(time/10))
    q_waves.append(Q * np.sin(time/10))


# set up IQ and Sum plot axis for side-by-side plotting
if mode != 'grid':
    num_plot_rows = 2**n
    ax_iq = [plt.subplot2grid((num_plot_rows, 2), (0, 0))]
    ax_sum = [plt.subplot2grid((num_plot_rows, 2), (0, 1))]
    for x in range(1, num_plot_rows):
        ax_iq.append(plt.subplot2grid((num_plot_rows, 2), (x, 0)))
        ax_sum.append(plt.subplot2grid((num_plot_rows, 2), (x, 1)))
        ax_iq[-1].sharey(ax_iq[0])
        ax_sum[-1].sharey(ax_sum[0])
else:
    # set up IQ axis for grid plotting
    grid_sz = 2**(n//2)
    ax_iq = [plt.subplot2grid((grid_sz, grid_sz), (0, 0))]
    ax_sum = None
    for x in range(grid_sz):
        for y in range(grid_sz):
            if x == 0 and y == 0:
                continue
            ax_iq.append(plt.subplot2grid((grid_sz, grid_sz), (y, x)))
            ax_iq[-1].sharey(ax_iq[0])


# plot all I and Q waves
if ax_sum is not None:
    for ax1, ax2, i, q in zip(ax_iq, ax_sum, i_waves, q_waves):
        ax1.plot(time, i)
        ax1.plot(time, q)
        ax2.plot(time, i+q)
else:
    for ax1, i, q in zip(ax_iq, i_waves, q_waves):
        ax1.plot(time, i)
        ax1.plot(time, q)


plt.figure()
# setup constellation diagram axis
ax_const = plt.subplot()

# plot constellation diagram
ax_const.scatter(i_symbols, q_symbols)
ax_const.axhline(y=0, color='k', linewidth=.5)
ax_const.axvline(x=0, color='k', linewidth=.5)
ax_const.set_aspect('equal', 'box')

ax_const.set_title(f"QAM - {n} bits per symbol")


plt.show()
