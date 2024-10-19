#!/usr/bin/env python
'''script to generate partitioned figures'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('-output',  default=None,          help='save output file instead of plotting')
p.add_argument('-plotcomp',action='store_true',   help='enable plotting composite sequenc')
p.add_argument('-plotsyms',action='store_true',   help='enable plotting symbols')
p.add_argument('-nstd',    default=0, type=float, help='noise standard deviation')
p.add_argument('-fc',      default=0, type=float, help='noise standard deviation')
p.add_argument('-plotcos', action='store_true',   help='enable plotting cosine of carrier offset')
args = p.parse_args()

# initialization
plt.style.use('seaborn-v0_8-darkgrid')
modmap = np.array((1,-1))
rng    = np.random.default_rng(12345)

# design interpolator from prototype
M, m, As, L, P = 8, 5, 60., 15, 8
num_symbols = L * P
interp = dsp.firinterp(M, m, As)

# generate random symbols and interpolate
symbols = rng.choice(modmap,num_symbols).astype(np.csingle)

# group symbols into discrete partitions
symbol_partitions = symbols.reshape((P,L))

# interpolate each partition
sample_partitions = np.empty((P,(L+2*m)*M), dtype=np.csingle )
for p in range(P):
    interp.reset()
    sample_partitions[p,:] = interp.execute(np.concatenate((symbol_partitions[p,:],np.zeros(2*m))))

# plot impulse and spectral responses
fig, ax = plt.subplots(1,figsize=(12,4))
for p in range(P):
    tp = np.arange( (L+2*m)*M )/M - m + p*L
    ax.plot(tp, np.real(sample_partitions[p,:]))
ax.set_xlabel('Time [symbols]')
ax.set_ylabel('Signal')
ax.grid(True, zorder=5)
ax.set(xlim=(-m,num_symbols+m),ylim=(-1.80,1.80))

if args.plotsyms:
    t0 = np.arange(num_symbols)
    ax.plot(t0, np.real(symbols), 'o', markersize=3, color='black')

if args.plotcomp:
    '''plot composite signal by adding individual partitions together'''
    samples = np.zeros(((num_symbols+2*m)*M,), dtype=np.csingle)
    t1 = np.arange((num_symbols+2*m)*M)/M - m
    for p in range(P):
        idx = p*L*M
        num = (L+2*m)*M
        samples[idx:(idx+num)] += sample_partitions[p,:]
    ax.plot(t1, np.real(samples), '-', linewidth=0.7, color='black')

if args.output is not None:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
else:
    plt.show()

