#!/usr/bin/env python
## Swift Biosciences 16S snapp workflow
## Author Benli Chai & Sukhinder Sandhu 20200502

from scipy.optimize import minimize
from math import log
import numpy as np
import pandas as pd

#to prepare and minimize matrix with constants and variables
def minimize_var(df, Sums):
    rowNames = df.index.values
    columnNames = df.columns.values
    Array = pd.DataFrame(df).to_numpy()
    print '\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'

    #get the mask position of the read to be minimized
    R = Sums.index(filter(lambda x: x != -1, Sums)[0])
    Sum = Sums[R]
    print 'ASV read', R, columnNames[R], Sum, 'read count to be allocated to:'
    print rowNames
    print 'Pre-minimization'
    print df.transpose().shape

    #the fractions of total read count to each reference
    frac = Array.sum(axis = 1)/Array.sum()
    print '\n', 'Mask =', Sums
    preFun = 'fun = lambda x:'#start the objetive function definition code
    preCons = {} #' assume multiple constraints keyed by the column number
    preBnds = 'bnds = [' #set bounds
    preX0 = 'x0 = [' #start x0

    for i in range(len(Array)): #iterate over references
        row = Array[i]
        sect = 'log(np.var(['
        values = ''
        for j in range(len(row)):
            value = Array[i][j] #all values
            if j == R:#variable setup
                value = 'x[%s]'%i
                preX0 += '%s,'%(Sums[j]*float(frac[i]))
                if not j in preCons:
                    preCons[j] = 'con%s = lambda x:'%j
                preCons[j] += 'x[%s] +'%i
                preBnds += '[0.01, %s],'%Sums[j]
            values += ',%s'%value
        sect += values.strip(',') + '])) +'
        preFun += sect
    preFun = preFun.strip('+') #finishing up the objective function

    preBnds = preBnds.strip(',') + ']' # finishing up the bounds
    preX0 = preX0.strip(',') + ']' #finishing test statement 
    print 'Initial test values:', preX0
    preCons = preCons[R].strip('+') + '- %s'%Sum
    preConAll = "cons = {'type':'eq','fun':con%s}"%R

    try:
        exec("""%s"""%preFun)
        exec("""%s"""%preCons)
        exec("""%s"""%preBnds)
        exec("""%s"""%preConAll)
        exec("""%s"""%preX0)
        exec("""sol = minimize(fun, x0, method='SLSQP', bounds=bnds, constraints=cons)""")
        print sol
        Array[:, R] = sol.x #insert the minimized values into the column to update the array
    except ValueError:# rare cases of overflow
        preX0 = [float(i) for i in preX0.split('[')[1].split(']')[0].split(',')] #obtain the values from string
        Array[:, R] = preX0

    print 'Post-minimization'
    Array = np.around(Array, 2)
    df = pd.DataFrame(Array, index=rowNames, columns = columnNames)
    print df.T.shape
    return df #return the minimized DataFrame
