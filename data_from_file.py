import matplotlib.pyplot as plt
import csv

x = []
y = []

with open('out.csv','r') as csvfile:
    plots = csv.reader(csvfile, delimiter=',')
    for row in plots:
        x.append(float(row[0]))
        y.append(float(row[1]))

plt.plot(x,y, label='Data loaded from file: out.cvs')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Hydroinformatics II (Olaf Kolditz)\nEX09: Finite-Difference-Method (implicit)')
plt.legend()
plt.savefig("EX08-fdm-implicit.png")
plt.show()
