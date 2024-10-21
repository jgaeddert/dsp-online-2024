#!/usr/bin/env python
'''script to generate interpolation figures'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('-output',  default=None,          help='save output file instead of plotting')
p.add_argument('-N',       default=240, type=int, help='number of symbols')
p.add_argument('-R',       default=None,type=int, help='number of symbols that aren not just zeros')
p.add_argument('-xticks',  default=None,type=int, help='xticks spacing')
p.add_argument('-plotsyms',action='store_true',   help='enable plotting symbols')
p.add_argument('-nstd',    default=0, type=float, help='noise standard deviation')
p.add_argument('-fc',      default=0, type=float, help='noise standard deviation')
p.add_argument('-plotcos', action='store_true',   help='enable plotting cosine of carrier offset')
p.add_argument('-plotcor', action='store_true',   help='enable plotting cross-correlation')
args = p.parse_args()

# initialization
plt.style.use('seaborn-v0_8-darkgrid')
modmap = np.array((1,-1))
rng    = np.random.default_rng(12345)

# design interpolator from prototype
M, m, As, num_symbols = 8, 5, 60., args.N
interp = dsp.firinterp(M, m, As)

# generate random symbols and interpolate
symbols = rng.choice(modmap,num_symbols).astype(np.csingle)
if args.R is not None:
    symbols[args.R:] = 0
signal  = interp.execute(np.concatenate((symbols,np.zeros(2*m))))

num_samples = (num_symbols+2*m)*M

t0 = np.arange(num_symbols)
t1 = np.arange(num_samples)/M - m

phasor = np.exp(2j*np.pi*args.fc*t1)
samples = signal * phasor

noise = rng.normal(0,args.nstd,2*num_samples).astype(np.single).view(np.csingle)
samples += noise

# plot time series
if args.plotcor:
    fig, (ax,ax2) = plt.subplots(2,1,figsize=(12,8))
else:
    fig, ax = plt.subplots(1,figsize=(12,4))
ax.plot(t1, np.real(samples))
if not args.plotsyms:
    ax.plot(t1, np.imag(samples))
ax.set_xlabel('Time [symbols]')
ax.set_ylabel('Signal')
ax.grid(True, zorder=5)
ax.set(xlim=(-m,num_symbols+m),ylim=(-1.80,1.80))
if args.xticks is not None:
    P = num_symbols // args.xticks
    ax.set_xticks(np.arange(P+1)*args.xticks)

if args.plotsyms:
    ax.plot(t0, np.real(symbols), 'o', markersize=3, color='black')

if args.plotcos:
    ax.plot(t1, np.real(phasor), ':', linewidth=0.5, color='black')

if args.plotcor:
    '''plot cross correlation vs. time lag'''
    rxy = np.fft.fftshift( np.fft.ifft( np.fft.fft(samples) * np.conj(np.fft.fft(signal)) ))
    #print('rxy',max(abs(rxy)),'n',len(samples), 'E{s^2}', np.sum(np.abs(signal)**2),)
    rxy /= np.sum(np.abs(signal)**2)
    txy = np.arange(num_samples)/M - m - num_symbols/2
    ax2.plot(txy, np.real(rxy), txy, np.imag(rxy))
    ax2.plot(txy, np.abs(rxy), ':', linewidth=0.5, color='black')
    ax2.set_xlabel('Lag [symbols]')
    ax2.set_ylabel('Cross Correlation')
    ax2.set(xlim=(-num_symbols/2-m,num_symbols/2+m),ylim=(-1.1,1.1))
    if args.xticks is not None:
        P = num_symbols // args.xticks
        ax2.set_xticks(np.arange(P+1)*args.xticks - P*args.xticks//2)

if args.output is not None:
    fig.savefig(args.output, dpi=200, bbox_inches='tight')
else:
    plt.show()

