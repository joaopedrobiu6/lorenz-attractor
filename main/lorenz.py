import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

WIDTH, HEIGHT, DPI = 1000, 750, 100

sigma, beta, rho = 10, 2.667, 28
x0, y0, z0 = 0, 1, 1.05

tmax, n = 100, 10000


def lorenz(t, X, sigma, beta, rho):
    x, y, z = X
    xp = -sigma * (x - y)
    yp = rho * x - y - x * z
    zp = -beta * z + x * y
    return xp, yp, zp


solved = solve_ivp(lorenz, (0, tmax), (x0, y0, z0),
                   args=(sigma, beta, rho),
                   dense_output=True)

t = np.linspace(0, tmax, n)

x, y, z = solved.sol(t)

print(f"({x}, {y}, {z}")

fig = plt.figure(facecolor='k', figsize=(WIDTH / DPI, HEIGHT / DPI))
ax = fig.gca(projection='3d')
ax.set_facecolor('k')
fig.subplots_adjust(left=0, right=1, bottom=0, top=1)

s = 10
cmap = plt.cm.winter
for i in range(0, n - s, s):
    ax.plot(x[i:i + s + 1],
            y[i:i + s + 1],
            z[i:i + s + 1],
            color=cmap(i / n),
            alpha=0.4)

#ax.set_axis_off()

plt.savefig('lorenz.png', dpi=DPI)
plt.show()