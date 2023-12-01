import math as m

r = 20 # m
rho = 1.225 # kg/m3
w = 2 # rad/s
Vo = 10 # m/s
c = 1 #m
profile = "FFA-W3-241.txt"

# Q1: The velocity in the wake is measured to u1=6m/s. Calculate the axial induction factor

# u1 = (1 - 2a) * Vo

u1 = 6 # m/s

a = (1 - u1/Vo)/2

print('Induction factor:', a)

# Q2: Neglecting the tangential induction factor compute the flow angle ϕ

# tan_flowangle = ((1-a)*Vo)/(w*r)

tan_flowangle = ((1-a)*Vo)/(w*r)

flowangle = m.atan(tan_flowangle)

print('Flow angle:', flowangle)

# Q3: The normal loading is measured to pn=1000N/m. Estimate the lift when neglecting the drag

# pn = l * cos(flowangle) + d * sin(flowangle)

pn = 1000 # N/m

l = pn / m.cos(flowangle)

print('Lift:', l)

# Q4: Again, neglecting the tangential induction factor calculate Vrel, Cl and the angle of attack α

# Vrel * sin(flowangle) = Vo * (1 - a)

# l = (rho * Vrel**2 * c * Cl) / 2

Vrel = (Vo * (1-a)) / m.sin(flowangle)

print('Vrel:', Vrel)

Cl = 2*l / (rho * Vrel**2 * c)

print('Cl:', Cl)

aoa = 47 # from table

print('Angle of attack:', aoa)
