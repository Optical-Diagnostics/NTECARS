using Revise
using NTECARS

T      = 1200.0
df     = CO2.MultiTemperatureDistribution(T_12 = T, T_3 = T, T_rot = T)
Q_calc = df.Q


states = CO2.database_of_allowed_states(;v₁_max = 3, l₂_max=6, v₃_max= 3, J_max = 90)
Q = sum([df(s) for s in states])