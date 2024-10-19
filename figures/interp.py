#!/usr/bin/env python
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt

# initialization
plt.style.use('seaborn-v0_8-darkgrid')
modmap = np.array((1+1j,1-1j,-1+1j,-1-1j))
rng    = np.random.default_rng(12345)

# design interpolator from prototype
M, m, As, num_symbols = 4, 6, 60., 120
interp = dsp.firinterp(M, m, As)

# generate random symbols and interpolate
symbols = rng.choice(modmap,num_symbols).astype(np.csingle)
samples = interp.execute(np.concatenate((symbols,np.zeros(2*m))))

# plot results
t0 = np.arange(num_symbols)
t1 = np.arange((num_symbols+2*m)*M)/M - m

# plot impulse and spectral responses
fig, (axi, axq) = plt.subplots(2,figsize=(12,8))
for ax,func,label in ((axi,np.real,'Real'),(axq,np.imag,'Imag')):
    ax.plot(t1, func(samples))
    ax.plot(t0, func(symbols), 'o')
    ax.set_xlabel('Time [samples]')
    ax.set_ylabel(label)
    ax.grid(True, zorder=5)

fig.savefig("interp.png", dpi=200, bbox_inches='tight')

#plt.show()

