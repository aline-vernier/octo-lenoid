import numpy as np
from matplotlib import pyplot as plt

w_dir = "G:\GPT\Salle Noire\Axially Magnetized Solenoid\\test"
w_file = "\\fieldmap.TXT"

data = []
with open(w_dir + w_file) as my_file:
    for line in my_file:
        data.append(
            line.split()
        )


data = np.transpose(np.array([map(float, x) for x in data[1:-1]]))
numR = np.size(np.unique(data[0]))
numP = np.size(data[0])


print np.argmax(data[3])

z = [data[1, x] for x in range(0, numP, numP/numR)]
Bz = [data[3, x] for x in range(0, numP, numP/numR)]


plt.plot(z, Bz)
plt.show()


