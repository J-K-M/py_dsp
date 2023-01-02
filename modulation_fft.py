import numpy as np
import matplotlib.pyplot as plt


bits_per_symbol = 2

carrier_freq = 10
symbol_rate = 5
num_symbols = 50
resolution = 1000


samp_freq = resolution*carrier_freq
samp_period = 1/samp_freq
symbol_duration = 1/symbol_rate


def generate_psk_symbol(val):
    phase = (2*np.pi / 2**bits_per_symbol) * val
    t = np.arange(0, symbol_duration, samp_period)
    return np.sin(carrier_freq * 2*np.pi * t + phase)


def generate_ask_symbol(val):
    t = np.arange(0, symbol_duration, samp_period)
    return val * np.sin(carrier_freq * 2*np.pi * t)


def generate_fsk_symbol(val):
    t = np.arange(0, symbol_duration, samp_period)
    freq_step = 2
    freq_shift = freq_step * (val - ((2**bits_per_symbol)/2))
    return np.sin((carrier_freq + freq_shift) * 2*np.pi * t)


symbols = []
for i in range(num_symbols):
    mod = np.random.randint(0, 2**bits_per_symbol)
    symbols.append(generate_fsk_symbol(mod))


sig = np.concatenate(symbols)
time = np.arange(0, sig.shape[-1]*samp_period, samp_period)

fft = np.fft.fft(sig)
freq = np.fft.fftfreq(sig.shape[-1], d=samp_period)


fig, ax = plt.subplot_mosaic([["top"], ['bottom']])
ax['top'].plot(time, sig)
ax['top'].grid(True)

ax['bottom'].plot(freq, np.abs(fft))
plt.xlim([0, 10*carrier_freq])

plt.show()
