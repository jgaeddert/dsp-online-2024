#!/usr/bin/env python
'''script to generate interpolation figures'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('-output',  default=None,          help='save output file instead of plotting')
p.add_argument('-plotsyms',action='store_true',   help='enable plotting symbols')
p.add_argument('-nstd',    default=0, type=float, help='noise standard deviation')
p.add_argument('-fc',      default=0, type=float, help='noise standard deviation')
args = p.parse_args()

# initialization
plt.style.use('seaborn-v0_8-darkgrid')
modmap = np.array((1,-1))
rng    = np.random.default_rng(12345)

# design interpolator from prototype
M, m, As, num_symbols = 4, 6, 60., 120
interp = dsp.firinterp(M, m, As)

# generate random symbols and interpolate
symbols = rng.choice(modmap,num_symbols).astype(np.csingle)
samples = interp.execute(np.concatenate((symbols,np.zeros(2*m))))

num_samples = (num_symbols+2*m)*M

t0 = np.arange(num_symbols)
t1 = np.arange(num_samples)/M - m

samples *= np.exp(2j*np.pi*args.fc*t1)

noise = rng.normal(0,args.nstd,2*num_samples).astype(np.single).view(np.csingle)
samples += noise

# plot impulse and spectral responses
fig, ax = plt.subplots(1,figsize=(12,4))
ax.plot(t1, np.real(samples))
if not args.plotsyms:
    ax.plot(t1, np.imag(samples))
ax.set_xlabel('Time [samples]')
ax.set_ylabel('Signal')
ax.grid(True, zorder=5)

if args.plotsyms:
    ax.plot(t0, np.real(symbols), 'o')

if args.output is not None:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
else:
    plt.show()

