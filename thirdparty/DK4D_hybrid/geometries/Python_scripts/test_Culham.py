a = 1.
R0 = 5.
N = 100
B0 = 1.
q0 = 0.8
qa = 0.7
alpha = 1
p0 = 10**5
Nomega = 50
Nr = 10
pa = p0/10
h = a/N
r = h*np.arange((N+1))
q = q0 + (qa - q0)*(r/a)**2
p = (p0 - pa)*(1-(r/a)**2)**alpha + pa
#Delta, E, T, P = Culham_equilibrum(q, p, R0, B0, 0.25, 0.1, h, Nomega, Nr)
Diff_Delta, Delta, E, Diff_E, T, Diff_T, P, f, g, grad_r, gr_cdot_gom, grad_omega = Culham_equilibrum(q, p, R0, B0, 0.25, 0.1, h, Nomega, Nr)
print(Delta)
