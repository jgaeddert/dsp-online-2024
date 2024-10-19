#!/usr/bin/env python
'''partition-based detector'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt

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
            self.R[p,:] = np.fft.fft(self.signal_partitions[p,:], self.nfft_0)

        # input buffers, transformed
        self.buf = np.zeros((10,), dtype=np.csingle) # input buffer, time series
        self.B = np.zeros(self.R.shape, dtype=np.csingle)
        #print(self.B.shape)

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
        fig, ax = plt.subplots(1,figsize=(12,4))
        tp = np.arange( (self.L+2*self.m)*self.M )/self.M - self.m
        for p in range(self.partitions):
            ax.plot(tp + p*self.L, np.real(self.signal_partitions[p,:]))
        plt.show()

    def execute(self, buf: np.ndarray):
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
        # compute partitioned...
        rxy = np.fft.ifft(self.B * np.conj(self.R), axis=0)
        print(rxy.shape)
        tmp = np.mean( np.abs(rxy), axis=1 )
        print(tmp.shape)
        #print(np.abs(rxy))
        return np.max(tmp)

if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-output',  default=None,          help='save output file instead of plotting')
    p.add_argument('-N',       default=120, type=int, help='number of symbols')
    p.add_argument('-P',       default=10,  type=int, help='number of partitions')
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
    s = np.concatenate((s,np.zeros(100,dtype=np.csingle)))
    n = len(s)
    s *= np.exp(0)

    # operate in blocks...
    num_blocks = len(s) // det.input_len
    print('len(s)',len(s), 'input_len', det.input_len, 'num_blocks', num_blocks)
    for i in range(num_blocks):
        rxy_max = det.execute(s[i:(i+det.input_len)])
        print(i, rxy_max)

