#!/usr/bin/env python
'''partition-based detector'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.tri as tri
import matplotlib.patches as patches

class qpart:
    def __init__(self, num_symbols=240, partitions=12, interp=2, m=5, As=60., seed=12345):
        self.num_symbols = num_symbols
        self.partitions  = partitions
        self.M = interp
        self.m = m
        self.L = self.num_symbols // self.partitions # symbols per partition

        # generate sequence
        modmap = np.array((1,-1))
        rng    = np.random.default_rng(seed)
        interp = dsp.firinterp(self.M, m, As)

        # generate random symbols and interpolate
        self.symbols = rng.choice(modmap,self.num_symbols).astype(np.csingle)

        # interpolate full signal
        self.signal = interp.execute(np.concatenate((self.symbols,np.zeros(2*self.m))))

        # group symbols into discrete partitions
        self.symbol_partitions = self.symbols.reshape((self.partitions,self.L))

        # interpolate each partition
        self.signal_partitions = np.empty((self.partitions,(self.L+2*self.m)*self.M), dtype=np.csingle )
        for p in range(self.partitions):
            interp.reset()
            self.signal_partitions[p,:] = interp.execute(np.concatenate((self.symbol_partitions[p,:],np.zeros(2*self.m))))

        # transform across partition
        self.nfft_0 = max( 2*self.L*self.M, (self.L+2*self.m)*self.M )
        self.R = np.empty((self.partitions, self.nfft_0),dtype=np.csingle)
        for p in range(self.partitions):
            part = self.signal_partitions[p,:]
            self.R[p,:] = np.fft.fft(part, self.nfft_0) / np.sum(np.abs(part)**2)

        # input buffers, transformed
        self.buf = np.zeros((self.input_len,), dtype=np.csingle) # input buffer, time series
        self.B = np.zeros(self.R.shape, dtype=np.csingle)

        # correlation output
        self.rxy = np.zeros((self.partitions,self.nfft_0), dtype=np.csingle)

        # grid
        self.nfft_1 = 2*self.partitions
        self.grid = np.zeros((self.nfft_1, self.nfft_0), dtype=np.csingle)

    def __repr__(self,):
        '''object representation'''
        return f'qpart(num_symbols={self.num_symbols},partitions={self.partitions},interp={self.M},m={self.m})'

    @property
    def sequence(self,):
        '''get a copy of the clean sequence'''
        return np.copy(self.signal)

    @property
    def input_len(self,):
        '''get the expected number of input samples'''
        return self.L * self.M

    def plot_partitions(self,):
        '''plot partitions'''
        fig, (ax,ax2) = plt.subplots(2,figsize=(12,8))
        tp = np.arange( (self.L+2*self.m)*self.M )/self.M - self.m
        for p in range(self.partitions):
            ax.plot(tp + p*self.L, np.real(self.signal_partitions[p,:]))
        for p in range(self.partitions):
            nfft = self.signal_partitions.shape[1]
            f = np.arange(nfft)/nfft - 0.5
            P = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.signal_partitions[p],nfft))))
            ax2.plot(f,P)
        plt.show()

    def plot_rxy(self,):
        '''plot correlation output'''
        '''
        fig, ax = plt.subplots(1,figsize=(12,4))
        txy = np.arange(self.nfft_0)/self.M - self.m
        for p in range(self.partitions):
            r = np.roll(self.rxy[p,:],self.m * self.M)
            ax.plot(txy + p*self.L, np.abs(r))
        #ax.set_prop_cycle(None) # reset color cycler
        #for p in range(self.partitions):
        #    r = np.roll(rxy[p,:],self.m * self.M) / 215
        #    ax.plot(txy + p*self.L, np.imag(r))
        plt.show()
        '''
        plot_box(v=self.rxy.T/np.max(np.abs(self.rxy)), labels=False)

    def plot_grid(self,):
        '''plot full grid'''
        plot_box(v=self.grid/np.max(np.abs(self.grid)), labels=False)

    def execute(self, buf: np.ndarray, plot=False):
        '''push block of samples'''
        # validate input
        if buf.shape != (self.input_len,):
            raise BaseException(f'expected input shape to be ({self.input_len},)')
        # shift transform buffers down
        self.B = np.roll(self.B, -1, axis=0)
        # shift samples into time-domain buffer
        self.buf = np.concatenate((self.buf[-self.input_len:], buf))
        # compute transform of input buffer and append to end of frequency buffer
        self.B[-1,:] = np.fft.fft(self.buf, self.nfft_0)
        # compute partitioned correlation output
        self.rxy = np.fft.ifft(self.B * np.conj(self.R), axis=1)
        #print('rxy',rxy.shape)
        self.grid = np.fft.fft(self.rxy, self.nfft_1, axis=0)
        if plot:
            self.plot_rxy()
            self.plot_grid()
        print('grid',self.grid.shape)
        return 0

def plot_box(v:np.ndarray, x=None, y=None, output=None, labels=True,
             marker=lambda x: np.abs(x)>0.9):
    '''plot boxes'''
    if x is None:
        x = np.arange(v.shape[0])
    if y is None:
        y = np.arange(v.shape[1])
    # scores as boxes
    fig1, ax1 = plt.subplots(figsize=(10,10))
    my_cmap = plt.cm.get_cmap('tab10')
    my_cmap.set_under('black')
    im = ax1.imshow(np.abs(v).T, vmin=0, vmax=1) # cmap=my_cmap) #, vmin=0.1)
    # add all tick marks...
    ax1.set_xticks(np.arange(len(x)))
    ax1.set_yticks(np.arange(len(y)))
    # ... and label them with the respective list entries
    ax1.set_xticklabels(x)
    ax1.set_yticklabels(y)
    # rotate the tick labels and set their alignment.
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right",rotation_mode="anchor")
    # loop over data dimensions and create text annotations.
    for i in range(len(y)):
        for j in range(len(x)):
            if labels:
                text = ax1.text(j, i, '%.2f\n%d' % (np.abs(v[j,i]), np.angle(v[j,i])*180/np.pi),ha="center", va="center", color="w")
            # add marker if converged
            #if marker(v[i,j]):
            #    rect = patches.Rectangle((j-0.45,i-0.45),0.9,0.9,fill=False,edgecolor='gold',linewidth=1)
            #    ax1.add_artist(rect)
    plt.show()

if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-output',  default=None,          help='save output file instead of plotting')
    p.add_argument('-N',       default=240, type=int, help='number of symbols')
    p.add_argument('-P',       default=24,   type=int, help='number of partitions')
    p.add_argument('-plotcomp',action='store_true',   help='enable plotting composite sequence')
    p.add_argument('-plotsyms',action='store_true',   help='enable plotting symbols')
    p.add_argument('-plotimag',action='store_true',   help='enable plotting imaginary component')
    p.add_argument('-fc',      default=0, type=float, help='noise standard deviation')
    p.add_argument('-fcapprox',action='store_true',   help='enable setting approximate fc offset for each partition')
    p.add_argument('-plotcos', action='store_true',   help='enable plotting cosine of carrier offset')
    p.add_argument('-plotcor', action='store_true',   help='enable plotting cross-correlation')
    args = p.parse_args()

    # create detector object
    det = qpart(num_symbols=args.N, partitions=args.P)
    print(det)
    det.plot_partitions()

    # get clean signal and apply offsets
    s = det.sequence
    # extend length
    s = np.concatenate((s,np.zeros(300,dtype=np.csingle)))
    n = len(s)
    s *= np.exp(0)

    #s = np.arange(n)

    # operate in blocks...
    num_blocks = len(s) // det.input_len
    print('len(s)',len(s), 'input_len', det.input_len, 'num_blocks', num_blocks)
    for i in range(num_blocks):
        num = det.input_len
        print(i)
        rxy_max = det.execute(s[(i*num):((i+1)*det.input_len)], plot=(i==args.P))
        #print(i, rxy_max)

