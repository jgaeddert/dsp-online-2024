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
        if self.L * self.partitions != self.num_symbols:
            raise BaseException(f'number of partitions must evenly divide number of symbols ({self.num_symbols})')

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
        #self.nfft_1 = 2*self.partitions
        self.nfft_1 = self.nfft_0 + 2 # artifically increase transform size for better grid
        self.grid = np.zeros((self.nfft_1, self.nfft_0), dtype=np.csingle)
        self.lag  = (np.arange(self.nfft_0) - self.nfft_0/2)/self.M # time offsets
        self.df   = (np.arange(self.nfft_1) - self.nfft_1/2)/(self.nfft_1*self.M*self.L) # frequency offsets

        # estimates
        self.rxy_max = 0
        self.dt_hat  = 0
        self.df_hat  = 0

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

    def execute(self, buf: np.ndarray):
        '''push block of samples and process'''
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
        # compute grid as transforms across time (scaled)
        self.grid = np.fft.fft(self.rxy, self.nfft_1, axis=0) / self.partitions
        # apply shift for convenience (plotting and extracting indices)
        self.grid = np.fft.fftshift(self.grid)
        # compute estimates
        row,col = np.unravel_index(np.argmax(np.abs(self.grid),axis=None), self.grid.shape)
        self.rxy_max = np.max(np.abs(self.grid))
        self.df_hat  = self.df[row]
        self.dt_hat  = self.lag[col]
        #print('max grid', self.rxy_max, 'at dt =', self.dt_hat, ', df =', self.df_hat)
        return self.rxy_max, self.dt_hat, self.df_hat

    # the remainder of this class is purely for plotting

    def plot_partitions(self,output=None):
        '''plot partitions'''
        fig, ax = plt.subplots(1,figsize=(12,4))
        tp = np.arange( (self.L+2*self.m)*self.M )/self.M - self.m
        for p in range(self.partitions):
            ax.plot(tp + p*self.L, np.real(self.signal_partitions[p,:]))
        #for p in range(self.partitions):
        #    nfft = self.signal_partitions.shape[1]
        #    f = np.arange(nfft)/nfft - 0.5
        #    P = 20*np.log10(np.abs(np.fft.fftshift(np.fft.fft(self.signal_partitions[p],nfft))))
        #    ax2.plot(f,P)
        # set figure properties
        ax.set_xlabel('Time [symbols]')
        ax.set_ylabel('Signal')
        ax.set_xticks(np.arange(self.partitions+1)*self.L)
        ax.set(xlim=(-self.m,self.num_symbols+self.m),ylim=(-1.8,1.8))
        if output is not None:
            fig.savefig(output, dpi=200, bbox_inches='tight')
        else:
            plt.show()

    def plot_rxy(self,output=None):
        '''plot correlation output'''
        fig, ax = plt.subplots(1,figsize=(12,4))
        txy = np.arange(self.nfft_0)/self.M - self.m
        # plot real
        for p in range(self.partitions):
            r = np.roll(self.rxy[p,:],self.m * self.M)
            ax.plot(txy + p*self.L, np.real(r))
        # plot imag + abs
        ax.set_prop_cycle(None) # reset color cycler
        for p in range(self.partitions):
            r = np.roll(self.rxy[p,:],self.m * self.M)
            ax.plot(txy + p*self.L, np.imag(r), '-', linewidth=0.5)
        # set figure properties
        ax.set_xlabel('Lag [symbols]')
        ax.set_ylabel('Cross Correlation')
        ax.set_xticks(np.arange(self.partitions+1)*self.L)
        ax.set(xlim=(-self.m,self.num_symbols+self.m),ylim=(-1.1,1.1))
        if output is not None:
            fig.savefig(output, dpi=200, bbox_inches='tight')
        else:
            plt.show()

    def plot_rxy_stacked(self,output=None):
        '''plot correlation output with stacked figures'''
        fig, ax = plt.subplots(self.partitions,figsize=(12,12))
        txy = np.arange(self.nfft_0)/self.M - self.m
        # plot each partition
        for p in range(self.partitions):
            for _ in range(p+1): # get appropriate color for plot
                c = ax[p]._get_lines.get_next_color()
            r = np.roll(self.rxy[p,:],self.m * self.M)
            ax[p].plot(txy, np.real(r), '-', color=c)
            ax[p].plot(txy, np.imag(r), '-', color=c, linewidth=0.5)
            ax[p].set(xlim=(-self.m,self.L+self.m),ylim=(-1,1))
            ax[p].set_ylabel(f'Correlation ({p})')
            ax[p].set_yticks((-1,-.5,0,.5,1))
        # set figure properties
        ax[-1].set_xlabel('Lag [symbols]')
        #ax[-1].set_xticks(np.arange(self.L+2*self.m+1)-self.m)
        if output is not None:
            fig.savefig(output, dpi=200, bbox_inches='tight')
        else:
            plt.show()

    def plot_grid(self,output=None):
        '''plot full grid'''
        fig,ax = plt.subplots(1,figsize=(8,8))
        my_cmap = plt.get_cmap('summer')
        #my_cmap.set_under('black')
        ax.pcolormesh(self.lag,self.df,np.abs(self.grid),shading='auto',vmin=0,vmax=1,cmap=my_cmap)
        ax.set_xlabel('Lag [symbols]')
        ax.set_ylabel('Frequency Offset')
        ax.grid(True, which='minor')
        # add marker for maximum
        ax.plot(self.dt_hat, self.df_hat, '.', color='black')
        if output is not None:
            fig.savefig(output, dpi=200, bbox_inches='tight')
        else:
            plt.show()

if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-output',  default=None,          help='save output file instead of plotting')
    p.add_argument('-N',       default=240, type=int, help='number of symbols')
    p.add_argument('-P',       default=8,   type=int, help='number of partitions')
    p.add_argument('-interp',  default=2,   type=int, help='interpolation rate')
    p.add_argument('-dt',      default=0,   type=int, help='timing offset [samples]')
    p.add_argument('-fc',      default=0, type=float, help='carrier offset [f/Fs]')
    args = p.parse_args()

    # set plot style
    plt.style.use('seaborn-v0_8-darkgrid')

    # create detector object
    det = qpart(num_symbols=args.N, partitions=args.P, interp=args.interp)
    print(det)
    #det.plot_partitions()

    # get clean signal and apply offsets
    s = det.sequence
    # extend length
    s = np.concatenate((s,np.zeros(700,dtype=np.csingle)))
    n = len(s)

    # add time delay
    s = np.roll(s, args.dt)

    # add carrier offset
    s *= np.exp(2j*np.pi*args.fc*np.arange(n))

    #s = np.arange(n)

    # operate in blocks...
    num_blocks = len(s) // det.input_len
    print('len(s)',len(s), 'input_len', det.input_len, 'num_blocks', num_blocks)
    for i in range(num_blocks):
        num = det.input_len
        rxy_max, dt_hat, df_hat = det.execute(s[(i*num):((i+1)*det.input_len)])
        if i==args.P:
            print('max grid', rxy_max, 'at dt =', dt_hat, ', df =', df_hat)
            det.plot_rxy_stacked()
            det.plot_grid() #'grid.png')

