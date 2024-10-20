#!/usr/bin/env python
'''partition-based detector'''
import argparse
import liquid as dsp
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.tri as tri
import matplotlib.patches as patches
import qpart

if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-output',  default=None,          help='save output file instead of plotting')
    p.add_argument('-N',       default=240, type=int, help='number of symbols')
    p.add_argument('-P',       default=8,   type=int, help='number of partitions')
    p.add_argument('-interp',  default=8,   type=int, help='interpolation rate')
    p.add_argument('-dt',      default=0,   type=int, help='timing offset [samples]')
    p.add_argument('-fc',      default=0, type=float, help='carrier offset [f/Fs]')
    args = p.parse_args()

    # set plot style
    plt.style.use('seaborn-v0_8-darkgrid')

    # create detector object
    det = qpart.qpart(num_symbols=args.N, partitions=args.P, interp=args.interp)
    print(det)
    #det.plot_partitions()

    # get clean signal and apply offsets
    s = det.sequence
    # extend length
    num_extend = max(0, det.input_len*(args.P+1) - len(s))
    s = np.concatenate((s,np.zeros(num_extend,dtype=np.csingle)))
    n = len(s)

    # add time delay
    s = np.roll(s, args.dt)

    # add carrier offset
    s *= np.exp(2j*np.pi*args.fc*np.arange(n))

    # add noise
    #s = np.arange(n)

    # operate in blocks...
    num_blocks = len(s) // det.input_len
    print('len(s)',len(s), 'input_len', det.input_len, 'num_blocks', num_blocks)
    for i in range(num_blocks):
        num = det.input_len
        rxy_max, dt_hat, df_hat = det.execute(s[(i*num):((i+1)*det.input_len)])
        if i==args.P:
            print('max grid', rxy_max, 'at dt =', dt_hat, ', df =', df_hat)
            break

    # save output figures
    det.plot_rxy(output = args.output + '.png')
    det.plot_rxy_stacked(output = args.output + '_rxy_stacked.png')
    det.plot_grid(output = args.output + '_grid.png')

    # plot correlator output across partitions
    lag  = args.dt
    corr = det.rxy[:,lag]
    tc   = np.arange(args.P)
    fig, ax = plt.subplots(1,figsize=(12,4))
    ax.plot(tc, np.real(corr), '-', linewidth=1.2, color='black')
    ax.plot(tc, np.imag(corr), '-', linewidth=0.5, color='black')
    for p in range(args.P):
        ax.plot([tc[p], tc[p]], [np.real(corr[p]), np.imag(corr[p])], 'o')
    ax.set_xlabel('Partition Index')
    ax.set_ylabel(f'Correlation at Lag Index {lag}')
    ax.set(xlim=(-0.25,args.P-0.75),ylim=(-1.2,1.2))
    fig.savefig(args.output + '_corr.png', dpi=200, bbox_inches='tight')

