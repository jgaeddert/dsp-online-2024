#!/usr/bin/env python
'''script to generate partitioned figures'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('-output',  default=None,          help='save output file instead of plotting')
p.add_argument('-N',       default=240, type=int, help='number of symbols')
p.add_argument('-P',       default=8,   type=int, help='number of partitions')
p.add_argument('-plotcomp',action='store_true',   help='enable plotting composite sequence')
p.add_argument('-plotsyms',action='store_true',   help='enable plotting symbols')
p.add_argument('-plotimag',action='store_true',   help='enable plotting imaginary component')
p.add_argument('-fc',      default=0, type=float, help='noise standard deviation')
p.add_argument('-fcapprox',action='store_true',   help='enable setting approximate fc offset for each partition')
p.add_argument('-plotcos', action='store_true',   help='enable plotting cosine of carrier offset')
p.add_argument('-plotcor', action='store_true',   help='enable plotting cross-correlation')
args = p.parse_args()

# initialization
plt.style.use('seaborn-v0_8-darkgrid')
modmap = np.array((1,-1))
rng    = np.random.default_rng(12345)

# design interpolator from prototype
M, m, As, num_symbols, P = 8, 5, 60., args.N, args.P
L = num_symbols // P
num_samples = (num_symbols + 2*m)*M
if num_symbols != L*P:
    raise BaseException(f'number of partitions must evenly divide number of symbols ({num_symbols})')
interp = dsp.firinterp(M, m, As)

# time vector
t0 = np.arange(num_symbols)
tp = np.arange( (L+2*m)*M )/M - m
t1 = np.arange(num_samples)/M - m

# generate random symbols and interpolate
symbols = rng.choice(modmap,num_symbols).astype(np.csingle)

# group symbols into discrete partitions
symbol_partitions = symbols.reshape((P,L))

# interpolate each partition
signal_partitions = np.empty((P,(L+2*m)*M), dtype=np.csingle )
for p in range(P):
    interp.reset()
    signal_partitions[p,:] = interp.execute(np.concatenate((symbol_partitions[p,:],np.zeros(2*m))))

# calculate composite signal
signal = np.zeros(((num_symbols+2*m)*M,), dtype=np.csingle)
for p in range(P):
    idx = p*L*M
    num = (L+2*m)*M
    signal[idx:(idx+num)] += signal_partitions[p,:]

# apply frequency offset
phasor = np.exp(2j*np.pi*args.fc*t1)
samples = signal * phasor
sample_partitions = np.copy(signal_partitions)
for p in range(P):
    idx = p*L*M
    num = (L+2*m)*M
    sample_partitions[p,:] *= phasor[idx+num//2] if args.fcapprox else phasor[idx:(idx+num)]

# plot time series
if args.plotcor:
    fig, (ax,ax2) = plt.subplots(2,1,figsize=(12,8))
else:
    fig, ax = plt.subplots(1,figsize=(12,4))
#plt.xticks(np.arange(P+1)*L)
for p in range(P):
    ax.plot(tp + p*L, np.real(sample_partitions[p,:]))
if args.plotimag:
    ax.set_prop_cycle(None) # reset color cycler
    for p in range(P):
        ax.plot(tp + p*L, np.imag(sample_partitions[p,:]), '-', linewidth=0.5)

ax.set_xlabel('Time [symbols]')
ax.set_ylabel('Signal')
ax.grid(True, zorder=5)
ax.set(xlim=(-m,num_symbols+m),ylim=(-1.80,1.80))

if args.plotsyms:
    ax.plot(t0, np.real(symbols), 'o', markersize=3, color='black')

if args.plotcomp:
    '''plot composite signal by adding individual partitions together'''
    ax.plot(t1, np.real(signal), '-', linewidth=0.7, color='black')

if args.plotcos:
    ax.plot(t1, np.real(phasor), ':', linewidth=0.5, color='black')

if args.plotcor:
    '''plot cross correlation vs. time lag'''
    txy = np.arange(num_samples)/M - m
    rxy = np.empty((P,num_samples),dtype=np.csingle)
    for p in range(P):
        partition = signal_partitions[p,:]
        rxy[p,:] = np.fft.ifft(np.fft.fft(samples) * np.conj(np.fft.fft(partition,num_samples)))
        #print('rxy',max(abs(rxy)),'n',len(samples), 'E{s^2}', np.sum(np.abs(signal)**2),)
        rxy[p,:] /= np.sum(np.abs(partition)**2)
        rxy[p,:] = np.roll(rxy[p,:], m*M)
        ax2.plot(txy, np.real(rxy[p,:]))
    ax2.set_prop_cycle(None) # reset color cycler
    for p in range(P):
        ax2.plot(txy, np.imag(rxy[p,:]),'-', linewidth=0.7)
    for p in range(P):
        ax2.plot(txy, np.abs(rxy[p,:]), ':', linewidth=0.5, color='black')
    ax2.set_xlabel('Lag [symbols]')
    ax2.set_ylabel('Cross Correlation')
    ax2.set(xlim=(-m,num_symbols+m),ylim=(-1.1,1.1))

if args.output is not None:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
else:
    plt.show()

