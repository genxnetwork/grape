import scipy
import scipy.optimize
import random
import math
import gzip
import logging
from collections import namedtuple
import statsmodels.api as sm
lowess = sm.nonparametric.lowess
import re
import operator as op
from functools import reduce


import numpy as np

w_rho_vals = []

def compute_mafs(vcf_input,max_window_size):
        global w_rho_vals
       # w_rho_vals.append(1)
        f = gzip.open(vcf_input,'rt')
        entries_started = False

        expected_vals = []

        site_counter = 0
        mafs = []
        for line in f:
                if ('#CHROM' in line):
                        entries_started = True
                elif (entries_started):
                        _values = line.split()
                        _pos = _values[1]
                        alt = _values[4]
                        if (len(alt.split(',')) > 1):
                           continue
                        i = 2
                        while(i < len(_values) and _values[i] != 'GT'):
                           i += 1
                        i += 1
                        tags = _values[7]
                        i = 9
                        num_ones = 0
                        num_zeros = 0
                        while (i < len(_values)):
                                        vals = re.split('\||/',_values[i])
                                        if (vals[0] == '1'):
                                                num_ones = num_ones + 1
                                        elif (vals[0] == '0'):
                                                num_zeros = num_zeros + 1
                                        if (vals[1] == '1'):
                                                num_ones = num_ones + 1
                                        elif (vals[1] == '0'):
                                                num_zeros = num_zeros + 1
                                        i = i + 1
                        v =  min(num_zeros,num_ones)/float(num_zeros+num_ones)
                        mafs.append(v)
        f.close()
        x_vals = [0]
        y_vals = [0]
        for o in range(1,max_window_size):
                window_size = o
                expected_vals = []
                for i in range(0,len(mafs),window_size):
                        expected_maf = 0
                        _sum = 0.0
                        for j in range(0,window_size):
                                if (i + j  >= len(mafs)):
                                        continue
                                _sum += mafs[i+j]
                        for j in range(0,window_size):
                                if (i + j >= len(mafs)):
                                        continue
                                if (_sum != 0):
                                        expected_maf += (mafs[i+j] * mafs[i+j]/_sum)
                        expected_vals.append(expected_maf)

                import numpy
                _a = numpy.array(expected_vals)

                pa = numpy.percentile(_a,1)
                rho_v = pa*pa + (1-pa)*(1-pa)
                x_vals.append(o)
                y_vals.append(rho_v)

        w_rho_vals = lowess(y_vals, x_vals,frac=0.1,return_sorted=False)
        #for k in range(0,len(x_vals)):
        #        print (x_vals[k],w_rho_vals[k])
        #       for j in w_rho_vals:
        #       print l,rho_v

def ncr_old(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer / denom


def w_rho_func(_w):

    if (_w >= len(w_rho_vals)):
            return w_rho_vals[len(w_rho_vals)-1]
    else:
     #   w = int(round(_w,0))
        return w_rho_vals[int(_w)]


def ncr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

def fp(e,N,L,rho,r,c,w):
    sum = 0.0
    #rho = w_rho_func(w)
    for i in range(0,c):
        sum += (ncr(r,i)* ((rho**(L/w))**i) * ( (1-rho**(L/w))**(r-i) ))
    return 1 - sum

def tp(er,N,L,rho,r,c,w):
    sum = 0.0
    for i in range(0,c):
        sum += (ncr(r,i) * (math.exp(-(er*L)/w)**i) * ((1-math.exp(-(L*er)/w))**(r-i)) )
    return 1 - sum


def compute_w(error_rate,N,L,rho,r,c,max_w=300):
	fun = lambda w: 0.5* N*(N-1)*fp(error_rate,N,L,w_rho_func(w),r,c,w) / tp(error_rate,N,L,w_rho_func(int(w)),r,c,w)

    #bnds = ((0, None))
    #x0 = 40
	lambda_val = 0.5* N*(N-1)
	w_min = -1
	_started = False
	for w in range(1,max_w):
		_tp = tp(error_rate,N,L,w_rho_func(int(w)),r,c,w)
		_fp = fp(error_rate,N,L,w_rho_func(w),r,c,w)
		if (round(_tp,2) -  0.5* N*(N-1)* round(_fp,2)) == 1:
			if (_started):
				continue
			else:
				w_min = w
				_started = True
		elif _started:
			print (w_min, w)
			return w_min, w
	return w_min, max_w


	#cons = ({'type': 'ineq', 'fun': lambda w:  0.5*N*(N-1)*fp(error_rate,N,L,w_rho_func(w),r,c,w) - tp(error_rate,N,L,w_rho_func(int(w)),r,c,w)})
    #res = scipy.optimize.minimize(fun,[x0], method='COBYLA', tol=1e-1,bounds=bnds)#,constraints=cons)


    #print (res)



if __name__ == '__main__':

    try:
        snakemake
    except NameError:
        test_dir = 'test_data/rapid'
        Snakemake = namedtuple('Snakemake', ['input', 'output', 'params', 'log'])
        snakemake = Snakemake(
            input={'vcf': f'test_data/rapid/phased/chr1.phased.vcf.gz',
                   'num_haplotypes': f'{2504*2}'},

            output=[f'{test_dir}/rapid_params'],
            params={'error_rate': '0.0025',
                    'min_snps': '600',
                    'num_runs': '10',
                    'num_success': '2'},
            log=[f'{test_dir}/estimate_rapid_params.log']
        )

    logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG, format='%(levelname)s:%(asctime)s %(message)s')

    vcf_input = snakemake.input['vcf']
    error_rate = float(snakemake.params['error_rate'])
    num_haplotypes = int(snakemake.input['num_haplotypes'])
    min_snps = int(snakemake.params['min_snps'])
    
    num_runs = int(snakemake.params['num_runs'])
    num_success = int(snakemake.params['num_success'])

    max_window_size = 300
    #error_rate = 0.0025
    #min_length_SNPs = 12000#14000
    rho_initial = 0.9
    compute_mafs(vcf_input, max_window_size)
    print(f'MAFS were computed successfully')
    w_min, w_max = compute_w(error_rate, num_haplotypes, min_snps, rho_initial, num_runs, num_success, max_window_size)
    print(f'w_min is {w_min}, w_max is {w_max}')
    with open(snakemake.output[0], 'w') as out:
        out.write(f'W={w_min}')