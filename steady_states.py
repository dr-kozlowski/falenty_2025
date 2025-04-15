import numpy as np
import matplotlib.pyplot as plt

class SteadyStates():

    def __init__(self, V0, N_gates, N, T) -> None:
        self.T = T # [deg C]
        self.Vs = np.linspace(-V0, V0, N)
        self.alphas = np.zeros((N_gates, N))
        self.betas = np.zeros((N_gates, N))
        # self.taus = np.zeros((3, N))
        self.taus_inf = np.zeros((N_gates, N))
        # self.gates = np.zeros((3, N))
        self.gates_inf = np.zeros((N_gates, N))
        self.k = 3 ** (0.1 * (self.T - 6.3))

    def calculate_gates_infs(self, i, V):
        self.gates_inf[:, i] = 1 / (1 + np.exp(V - self.V_gates) / self.as)

    def calculate_rate_equations(self, i, V):

        



        alpha_m = 1e3 * (2.5 - 100 * V) / (np.exp(2.5 - 100 * V) - 1)
        alpha_n = 1e3 * (0.1 - 10 * V) / (np.exp(1 - 100 * V) - 1)
        alpha_h = 70 * np.exp(-50 * V)

        self.alphas[:, i] = np.array([alpha_m, alpha_n, alpha_h])

        beta_m = 4e3 * np.exp(-500 * V / 9)
        beta_n = 125 * np.exp(-25 * V / 2)
        beta_h = 1e3 / (np.exp(3 - 100 * V) + 1)
        self.betas[:, i] = np.array([beta_m, beta_n, beta_h])



    def run_simulation(self):


        for i, V in enumerate(self.Vs):

            self.calculate_rate_equations(i, V)
            self.taus_inf[:, i] = 1 / (self.k * (self.alphas[:, i] + self.betas[:, i]))
            self.gates_inf[:, i] = self.alphas[:, i] / (self.alphas[:, i] + self.betas[:, i])
