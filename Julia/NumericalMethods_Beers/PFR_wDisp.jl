#
#  PFR with dispersion and rxn
#  Rxn: A + B -> C, B + C - > D
#


L = 10 # length
vz = 1 # velocity
D = 1.0e-4 # dispersion coefficient

# Reaction moddel parameters
k1 = 1.0
k2 = 1.0

# inlet concentrations
ca0 = 1.0
cb0 = 1.0
cc0 = 0.0
cd0 = 0.0

D_uL = D/(vz * L)

function disprxn(t, c, dc)

  D_uL = 1.0e-5
  nc = 4
  ngrid = 10
  c0 = [1.0;1.0;0.0;0.0]

  dz = 1/(ngrid+1)

  Peloc = dz/D_uL

  alpha_lo = 1.0/dz + D_uL/dz^2
  alpha_mid = -1.0/dz - 2.0 * D_uL/dz^2
  alpha_hi = D_uL/dz^2

  dconc = zeros(ngrid, nc)

  k = 0
  for i in 1:ngrid
    for j in 1:nc
      k = k + 1
      conc[i,j] = c[k]
    end
  end

  for i in 2:(ngrid-1)
    for j in 1:nc
      dconc[i,j] = alpha_lo * conc[i-1,j] + alpha_mid * conc[i,j] + alpha_hi * conc[i+1,j]
    end
  end

  for j in 1:nc
    dconc[1,j] = (alpha_mid + alpha_lo/(1.0+Peloc)) * conc[1,j] + alpha_lo * Peloc * c0[j]/(1 + Peloc) + alpha_hi * conc[j, 2]
  end

  for j in 1:nc
    dconc[ngrid,j] = alpha_lo * conc[ngrid - 1,j] + (alpha_mid + alpha_hi) * conc[j, ngrid]
  end

  k = 0
  for i in 1:ngrid
    for j in 1:nc
      k = k + 1
      dc[k] = dconc[i,j]
    end
  end


end

using DifferentialEquations

cinit = zeros(4*10)

prob = ODEProblem(disprxn, cinit, (0.0,3.0))
csol = solve(prob, CVODE_BDF())

typeof(disprxn)
