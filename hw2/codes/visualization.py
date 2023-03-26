import matplotlib.pyplot as plt
import numpy as np

X = []
Y = []
U = []
with open("sol.txt", "r") as file:
    for line in file:
        floats = [float(i) for i in line.split(' ')]
        X+=[floats[0]]
        Y+=[floats[1]]
        U+=[floats[2]]

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
ax.plot_trisurf(np.array(X), np.array(Y), np.array(U), cmap=plt.cm.coolwarm)

plt.savefig('solution_plot.svg')

plt.show()

