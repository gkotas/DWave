from sympy import *
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

from util import makeCoeffMatrix

# Makes printing sympy pretty
init_printing()

# Enable printing of
DEBUG_PRINTS = True

# Values of each edge
S_vals = [1, 2, 3, 4, 5]

# Sigma looping symbols
i = symbols('i', cls=Idx)
# Subscript symbols
S, x = symbols('S x', cls=IndexedBase)
# Other symbols
g = symbols('g')

# DWave parameters
chainstrength = 1000
numruns = 40
gamma = 100

# A dictionary of substitutions
# Key: symbol
# Value: symbol's value
subs = {}
subs[g] = gamma
for j in range(len(S_vals)):
    subs[S[j]] = S_vals[j]

# Constraint: Only use three numbers
const1 = expand((summation(x[i], (i, 0, 4)) - 3)**2)
# Constraint: Sum of numbers is 8
const2 = expand((summation(S[i]*x[i], (i, 0, 4)) - 8)**2)

# Make the polynomial of constraints and list of x variables
poly = Poly(g*(const1 + const2), [x[j] for j in range(len(S_vals))])

# The matrix of symbols
mat = makeCoeffMatrix(poly, len(S_vals))
if DEBUG_PRINTS:
    pprint(mat)
# substitute in the values
mat_subs = mat.subs(subs)
if DEBUG_PRINTS:
    pprint(mat_subs)

# Q Matix from the sympy matrix
Q = {}
for j in range(len(S_vals)):
    for k in range(len(S_vals)):
        Q[(j, k)] = mat_subs[j, k]

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)

## ------- Return results to user -------
R = iter(response)
E = iter(response.data())
for line in response:
    sample = next(R)
    S1 = [S_vals[i] for i in sample if sample[i] > 0]
    S0 = [S_vals[i] for i in sample if sample[i] < 1]
    print("S1 Sum: ", sum(S1), "\t", S1)
