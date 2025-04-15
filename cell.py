import numpy as np
import matplotlib.pyplot as plt

class Cell:
    def __init__(self, T, dt, c0):
        self.T = T # ms
        self.dt = dt # ms
        self.N = int(self.T / self.dt) # number of steps
        self.c0 = c0 # initial concentration in the cell [num of ions / length ** 3]
        self.C = np.zeros(shape=(self.N, 1)) # concentration of ions at a given time point [num of ions / length ** 3]
        self.t = 0

    def set_j_prod(self, j_prod):
        self.j_prod = j_prod

    def set_j_deg(self, j_deg):
        self.j_deg = j_deg

    def set_j_in(self, j_in):
        self.j_in = j_in        

    def set_j_out(self, j_out):
        self.j_out = j_out
        
    def run(self):
        print('jjajaja')
        self.C[0] = self.c0 # assign the initial value
        print('2 jjajaja')

        for i in range(0, self.N-1):
           print('i: ', i)
           self.C[i+1] = self.C[i] + self.dt * (self.j_prod(self) - self.j_deg(self) + self.j_in(self) - self.j_out(self))
           self.t += self.dt