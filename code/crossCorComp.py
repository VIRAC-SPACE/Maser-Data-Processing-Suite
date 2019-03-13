import numpy as np
import matplotlib.pyplot as plt

t = [0, 1, 2, 4]
tfull = [0, 1, 2, 3, 4, 5, 6, 7]

x1 = [1, -3, 0.5, 4]
x2 = [-6, 1,5, 4, 1]

y_same = np.correlate(x1, x2, "same")
y_full = np.correlate(x1, x2, "full")

print ("y_same", y_same)
print("y_full", y_full)

plt.plot(t, x1)
plt.plot(tfull, y_full)
plt.show()