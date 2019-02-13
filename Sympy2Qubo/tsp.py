from sympy import *
from dwave.system.samplers import DWaveSampler
from dwave.system.composites import EmbeddingComposite

from util import makeCoeffMatrix

# Makes printing sympy pretty
init_printing()

# Enable printing of
DEBUG_PRINTS = False

# Values of each edge
d_vals = [6, 10, 4, 2, 5, 2, 1, 8, 8, 7]  # 5 cities: Min = 20
# Number of edges
E = len(d_vals)

# Sigma looping symbols
i = symbols('i', cls=Idx)
# Subscript symbols
d, x = symbols('d x', cls=IndexedBase)
# Other symbols
g = symbols('g')

# DWave parameters
chainstrength = 100
numruns = 100
gamma = 50

# A dictionary of substitutions
# Key: symbol
# Value: symbol's value
subs = {}
subs[g] = gamma
for j in range(len(d_vals)):
    subs[d[j]] = d_vals[j]

# Constraints: Sum of all edges with "a" must equal 2
# E.g. (ab + ac + ad + ad - 1)**2
# Repeat for all cities a, b, c, d, and e
const1 = expand((x[0] + x[1] + x[2] + x[3] - 2)**2)
const2 = expand((x[0] + x[4] + x[5] + x[6] - 2)**2)
const3 = expand((x[1] + x[4] + x[7] + x[8] - 2)**2)
const4 = expand((x[2] + x[5] + x[7] + x[9] - 2)**2)
const5 = expand((x[3] + x[6] + x[8] + x[9] - 2)**2)

# Objective Funtion: Minimize sum of edges in path
obj = expand(summation(d[i]*x[i], (i, 0, E - 1)))

# Make the polynomial of constraints plus objective function, and list of
# x variables
poly = Poly(g*(const1 + const2 + const3 + const4 + const5) + obj,
            [x[j] for j in range(E)])

# The matrix of symbols
mat = makeCoeffMatrix(poly, E)
if DEBUG_PRINTS:
    pprint(mat)
# substitute in the values
mat_subs = mat.subs(subs)
if DEBUG_PRINTS:
    pprint(mat_subs)

# Q Matix from the sympy matrix
Q = {}
for j in range(E):
    for k in range(E):
        Q[(j, k)] = mat_subs[j, k]

sampler = EmbeddingComposite(DWaveSampler())
response = sampler.sample_qubo(Q, chain_strength=chainstrength, num_reads=numruns)

for index in range(len(response.record.sample)):
    sample = response.record.sample[index]
    s = sum([a*b for a,b in zip(sample, d_vals)])
    print("Sample:", sample, "\tEnergy:", response.record.energy[index],
          "\tCost:", s)
