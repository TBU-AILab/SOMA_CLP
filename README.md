# SOMA_CLP

Implementation of the ESP-SOMA in Python

Full paper: ESP-SOMA Solving CEC 2019 100-Digit Challenge

Full paper: ESP-SOMA solving Constrained Technological Design Optimization Problem

Usage
Sample use can also be seen at the end of the file main.py.

dim = 10 #dimension size
NP = 10 #population size
maxFEs = 10000 * dim #maximum number of objective function evaluations
gap = 7
adaptivePRT = 1
pathLength = 3.0
step = 0.11

sphere = Sphere(dim) #defined test function
soma = ESP_SOMA(dim, maxFEs, sphere, NP, gap, pathLength, step, adaptivePRT)
resp = soma.run()
print(resp)
Output resp then includes optimized values features and value of objective function ofv. Also, the id of particle is included.

File descriptions
main.py
The main file contains the main class ESP_SOMA and one sample test function class Sphere.
