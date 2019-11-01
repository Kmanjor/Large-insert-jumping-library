import numpy as np
import math
import matplotlib.pyplot as plt


def lognormfun(x, mu, sigma):
    pdf = np.exp(-((np.log(x) - mu)**2)/(2 * sigma**2)) / (x * sigma * np.sqrt(2 * np.pi))
    return pdf

def normfun(x, mu, sigma):
    pdf = np.exp(-((x - mu)**2)/(2*sigma**2)) / (sigma * np.sqrt(2*np.pi))
    return pdf

E = 5500
V = 2000**2

mu = np.log((E**2)/math.sqrt(V+E**2))
sigma = math.sqrt(np.log(V/(E**2) + 1))

x = abs(np.round(np.random.lognormal(mu, sigma, 10000)))
print(x.mean(), x.std())
print(x.min(), x.max())
plt.hist(x, bins=100, density=True)
mm = np.arange(x.min(), x.max(), 0.1)
yy = lognormfun(mm, mu, sigma)
yy2 = normfun(mm, E, np.sqrt(V))
plt.plot(mm, yy2)
plt.plot(mm, yy)
plt.show()
