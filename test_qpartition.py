#!/usr/bin/env python
'''simple qpartition tests'''
import argparse
import liquid as dsp
import numpy as np
import qpart
import pytest

def harness_qpartition(num_symbols: int, partitions: int, interp: int, fc: float = 0,
                       dt: int = 0):
    '''basic qpartition testing'''
    # create detector object
    det = qpart.qpart(num_symbols, partitions, interp)

    # get clean signal and apply offsets
    s = det.sequence
    # extend length
    num_extend = max(0, det.input_len*(partitions+1) - len(s))
    s = np.concatenate((s,np.zeros(num_extend,dtype=np.csingle)))
    n = len(s)

    # add time delay
    s = np.roll(s, dt)

    # add carrier offset
    s *= np.exp(2j*np.pi*fc*np.arange(n))

    # operate in blocks
    num_blocks = len(s) // det.input_len
    #print('len(s)',len(s), 'input_len', det.input_len, 'num_blocks', num_blocks)
    for i in range(num_blocks):
        num = det.input_len
        rxy_max, dt_hat, fc_hat = det.execute(s[(i*num):((i+1)*num)])
        dt_hat *= interp # scale by interpolation length
        if i==partitions:
            break

    print('max grid', rxy_max, 'at dt =', dt_hat, ', df =', fc_hat)
    assert np.abs(rxy_max -  1) < 0.05
    assert np.abs(fc_hat  - fc) < 1e-5
    assert np.abs(dt_hat  - dt) < 2

def test_qpartition():
    '''run tests for qpartition object'''
    harness_qpartition(240, 12, 4, fc=     0, dt= 3)

    harness_qpartition(240, 12, 4, fc= 0.001, dt= 0)
    harness_qpartition(240, 12, 4, fc=-0.001, dt= 0)

    harness_qpartition(240, 12, 4, fc=-0.001, dt= 0)
    harness_qpartition(240, 12, 4, fc=-0.001, dt= 1)
    harness_qpartition(240, 12, 4, fc=-0.001, dt= 8)
    harness_qpartition(240, 12, 4, fc=-0.001, dt=-7)

if __name__=='__main__':
    test_qpartition()
