import sympy as sp
from sympy.core.sympify import sympify

import dimod
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector

#DWAVE_API_TOKEN=DEV-754a213785424e0b988b54a1614c70f03ed5d88f

num = int(input("number: "))

# converts int into array of base 2 bits
def convert_to_binary(num):
  binaryForm = []

  while (num > 0):
    binaryForm.append(num % 2)
    num = num // 2

  return binaryForm

n = convert_to_binary(num)

print(n)
print("bits: " + str(len(n)))

# TODO: see if you can implement better method
# for guessing bit length of factors
# ex: 91 -> 7*13 (3b * 4b), but this method
# returns 4b * 4b
#p_bits = 4
p_bits = ((len(n) + 1) // 2)
q_bits = ((len(n) + 1) // 2) + ((len(n) + 1) % 2)
#q_bits = 3

# binary variables of form 'x_i'
x = sp.symbols('x_1:%d' % ((p_bits - 1) + (q_bits - 1) + 1))

p = "1"
for i in range(p_bits - 1):
  p += " + " + str(2**(i+1)) + "*" + str(x[i])

# forms polynomial for bits in p
p = sp.sympify(p)

q = "1"
for i in range(q_bits - 1):
  q += " + " + str(2**(i+1)) + "*" + str(x[i + (p_bits - 1)])

# forms polynomial for bits in q
q = sp.sympify(q)

cost = sp.expand((num - (p * q))**2)

# substitues x_i**2 with x_i
# since every bit value squared is the same 
# as its original value
for i in range(len(x)):
  cost = cost.subs(sympify(str(x[i]) + "**2"), sympify(str(x[i])))

#ORDER REDUCTION:
# cxyz
# c < 0 -> cw(x + y + z -2)
# c > 0 -> c{w(x+y+z−1) + (xy+yz+zx) − (x+y+z) + 1}
# loops twice to cover quartic terms
for i in range(2):
  for term in cost.as_ordered_terms():
    # if higerorder term of degree 3
    if(len(term.free_symbols) == 3):
      # adds auxilary variable
      x = sp.symbols('x_1:%d' % (len(x) + 2))

      # variables in term
      vars = list(term.free_symbols)

      # coefficient of term
      vars_str = ""
      for var in vars:
        vars_str += str(var) + " * "
      vars_str = vars_str[:-3]
      coeff = term.as_coefficient(sympify(vars_str))

      # c > 0 -> c{w(x+y+z−1) + (xy+yz+zx) − (x+y+z) + 1}
      if(coeff > 0):
        substitution = sp.sympify(
          coeff * ((x[len(x) - 1]) * (vars[0] + vars[1] + vars[2] - 1) + \
          (vars[0] * vars[1] + vars[1] * vars[2] + vars[2] * vars[0]) - \
          (vars[0] + vars[1] + vars[2]) + 1)
        )

        cost = sp.expand(cost.subs(term, substitution))
      # c < 0 -> cw(x + y + z -2)
      elif (coeff < 0):
        substitution = coeff * (x[len(x) - 1]) * (vars[0] + vars[1] + vars[2] - 2)
      cost = sp.expand(cost.subs(term, substitution))


print ("COST:")
print (cost)

Q = {}
offset = 0
max_bias = 0

for term in cost.as_ordered_terms():
  vars = list(term.free_symbols)

  coeff = term
  if (len(vars) >= 1):
    for var in vars:
      coeff = coeff.coeff(var)

  if (len(vars) == 2):
    Q[((x.index(vars[0]) + 1), (x.index(vars[1]) + 1))] = coeff
    if(abs(coeff) > max_bias):
      max_bias = abs(coeff)
  elif (len(vars) == 1):
    Q[((x.index(vars[0]) + 1), (x.index(vars[0]) + 1))] = coeff
  elif (len(vars) == 0):
    offset = coeff

print (Q)
print(offset)
print(f"chain strength: {max_bias}")

bqm = dimod.binary.BinaryQuadraticModel.from_qubo(Q, offset)

sampler = EmbeddingComposite(DWaveSampler(token="DEV-754a213785424e0b988b54a1614c70f03ed5d88f"))
sample_set = sampler.sample(bqm, chain_strength=(int(max_bias)), num_reads=2000)
print("Using DWaveSampler()")
#print(sample_set)
print(f"p: {p}")
print(f"q: {q}")


real = p_bits + q_bits - 2

vars = (sample_set.variables)[:real]

print(dimod.keep_variables(sample_set, vars).aggregate())

dwave.inspector.show(sample_set)
