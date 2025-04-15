import numpy as np
import matplotlib.pyplot as plt
from steady_states import SteadyStates

class HH_neuron():
    def __init__(self, T, t_end, dt, g_na, g_k, g_l,
                 V_na, V_k, V_l, V_rest, Cm, I_stimulus: np.array) -> None:
        self.T = T # [C], + 273.15 # [K], temperature with converstion from grad C to Kelvins
        self.calculate_temperature_correction()
        
        self.g_na = g_na # [S]
        self.g_k = g_k # [S]
        self.g_l = g_l # [S]

        self.V_na = V_na # [V]
        self.V_k = V_k # [V]
        self.V_l = V_l # [V]
        self.V_rest = V_rest # [V]
        self.V0 = V_rest # [V]

        self.Cm = Cm # [F]
        
        self.t = t_end # duration [s]
        self.dt = dt # [s]
        self.N = int(self.t / self.dt) # number of steps
        self.ts = np.linspace(0, self.t, self.N)
        self.I_stimulus = I_stimulus # [A]

# ------------------------------ delayed rectifier
        self.E_k = E_k # [V]
        
        self.gate_d_inf = np.zeros(self.N)
        self.gate_d_relaxation = np.zeros(self.N)
        self.V_d = V_n # half-maximum for delayed rectifier [V]
        self.s_n = s_n # step width for delayed rectifier [V]
        self.V_kn = V_kn # half-maximum for delayed rectifier [V]
        self.s_kn = s_kn # step width for delayed rectifier [V]
        

# ----------------------------- calcium activated outward current
        self.z_Ca = z_Ca # valence of calcium
        self.F = F # Faraday constant [C / mol]
        self.cell_vol = cell_vol # [L] or m**3
        self.c_iCa = self.calculate_c_iCa() # [M / C]


        self.Ca_conc = np.zeros(self.N)
        self.Ca_0 = Ca_0 # [M]
        self.k_Ca = k_Ca # [Hz]

# ----------------------------- transient A gate
        self.V_X = V_X # half-maximum potential for voltage-dependent weighting factor [V]
        self.s_X = s_X # [V]

        self.V_A = V_A # [V]
        self.gate_A = np.zeros(self.N)
        self.gate_A_inf = np.zeros(self.N)
        self.s_A = s_A # step width for transient A gate activation [V]

        self.gate_b_A1 = np.zeros(self.N)
        self.V_bA1 = V_bA1 # [V]        
        self.b_A1_inf = np.zeros(self.N)
        self.s_bA1 = s_bA1 # step width for transient A gate inactivation 1 [V]

        self.gate_b_A2 = np.zeros(self.N)
        self.b_A2_inf = np.zeros(self.N)
        self.V_bA2 = V_bA2 # [V]
        self.s_bA2 = s_bA2 # step width for transient A gate inactivation 2 [V]

        self.k_b_A1 = k_b_A1 # relaxation rate for transient A gate inactivation 1

        self.c_A = c_A 
        self.V_kA = V_kA # [V]
        self.s_kA = s_kA # step width for transient A gate relaxation A2 [V]

    def calculate_temperature_correction(self):
        self.k = 3 ** (0.1 * (self.T - 6.3)) # 6.3 is in Celsius

    def calculate_c_iCa(self):
        return 1 / (self.z_Ca * self.F * self.cell_vol)

    def steady_states(self, i, V):
    # MUST BE RUN AT THE BEGINNING OF EACH LOOP ENTRY
        # delayed rectifier activation and relaxation gates
        self.gate_d_inf[i] = 1 / (1 + np.exp((V - self.V_d) / self.s_n))
        self.gate_d_relaxation[i] = self.c_n / ((1 + np.exp(V - self.V_kn) / self.s_kn))

        # calcium activated outward current
        self.a_oCa_inf[i] = (1 / ((V - self.V_ao1 + 
                               self.f * self.Ca_conc[i]) / self.s_ao1)) * (1 / 1 + 
                                                                           ((V - self.V_ao2 
                                                                             + self.f * self.Ca_conc) / self.s_ao2)) * (self.Ca_conc[i] / (self.c1 + self.Ca_conc[i]))
        self.b_oCa_inf[i] = self.c2 / (self.c3 + self.Ca_conc[i])

        # transient A gate
        self.gate_A_inf[i] = 1 / (1 + np.exp((V - self.V_A) / self.s_A))
        self.gate_b_A1_inf[i] = 1 / (1 + np.exp((V - self.V_bA1) / self.s_bA1))
        self.gate_b_A2_inf[i] = 1 / (1 + np.exp((V - self.V_bA2) / self.s_bA2))
        self.gate_A_relaxation[i] = self.c_A / ((1 + np.exp(V - self.V_kA) / self.s_kA))

    def calculate_delayed_rectifier_a(self, i):
        self.gate_d[i+1] = self.gate_d[i] + self.dt * (self.gate_d_inf[i] - self.gate_d[i]) * self.gate_d_relaxation[i]

    def calculate_Ca_activated_outward_gate(self, i, V):
        self.gate_a_oCa[i+1] = self.gate_a_oCa[i] + self.dt * (self.a_oCa_inf[i] - self.gate_a_oCa[i])
        self.gate_b_oCa[i+1] = self.gate_b_oCa[i] + self.dt * (self.b_oCa_inf[i] - self.gate_b_oCa[i])

    def calculate_Ca_concentration(self, i, V):
        self.Ca_conct[i+1] = self.Ca_conc[i] + self.dt * (-self.c_iCa * self.i_Ca[i] - 
                                                          self.k_Ca * self.Ca_conc[i] + self.k_Ca * self.Ca_0)

    def calculate_transient_A_gates(self, i, V):
        self.gate_A[i+1] = self.gate_A[i] + self.dt * (self.gate_A_inf[i] - self.gate_A[i]) * self.gate_A_relaxation[i]
        self.gate_b_A1[i+1] = self.gate_b_A1[i] + self.dt * (self.b_A1_inf[i] - self.gate_b_A1[i]) * self.k_b_A1
        self.gate_b_A2[i+1] = self.gate_b_A2[i] + self.dt * (self.b_A2_inf[i] - self.gate_b_A2[i]) * self.gate_A_relaxation[i]


    def calculate_ionic_currents(self, i, V):
        self.i_d[i] = self.g_d * (self.gate_d[i] ** 4) * (V - self.E_k) # delayed rectifier current

        self.i_oCa[i] = self.g_oCa * self.a_o * self.b_o * (V - self.E_K) # calcium activated outward current

        # transient A-like current
        X = lambda V: 1 / (1 + np.exp((V - self.V_X) / self.s_X)) # voltage-dependent weighting factor
        self.i_A[i] = self.g_A * (self.gate_A[i] ** 3) * (self.X(V) * self.gate_b_A1[i] + (1 - self.X(V) * self.gate_b_A2)) * (V - self.E_K)


    def calculate_A(self):

        self.Am = - self.k * (self.alpha_m + self.beta_m)
        self.Ah = - self.k * (self.alpha_h + self.beta_h)
        self.An = - self.k * (self.alpha_n + self.beta_n)

    def calculate_B(self):

        self.Bm = self.alpha_m * self.k
        self.Bh = self.alpha_h * self.k
        self.Bn = self.alpha_n * self.k


    def calculate_rate_equations(self, V):
        self.alpha_m = 1e3 * (2.5 - 100 * V) / (np.exp(2.5 - 100 * V) - 1)
        self.alpha_n = 1e3 * (0.1 - 10 * V) / (np.exp(1 - 100 * V) - 1)
        self.alpha_h = 70 * np.exp(-50 * V)

        self.beta_m = 4e3 * np.exp(-500 * V / 9)
        self.beta_n = 125 * np.exp(-25 * V / 2)
        self.beta_h = 1e3 / (np.exp(3 - 100 * V) + 1)

        self.k = 3 ** (0.1 * (self.T - 6.3))

        self.calculate_A()
        self.calculate_B()

    def ode(self, x, A, B):
        return x * np.exp(A * self.dt) + (B / A) * (np.exp(A * self.dt) - 1)

    def calculate_gates(self, i):
        # self.calculate_rate_equations(np.Vs[i])

        # ode = lambda x, A, B: x * np.exp(A * self.dt) + (B / A) * (np.exp(A * self.dt) - 1)
        return np.array([self.ode(self.gates[0, i], self.Am, self.Bm), 
                         self.ode(self.gates[1, i], self.An, self.Bn), 
                         self.ode(self.gates[2, i], self.Ah, self.Bh)]).T

    def calculate_ionic_currents(self, i, V):
        m, n, h = self.gates[:, i]
        I_na = self.g_na * (m ** 3) * h * (V - self.V_na)
        I_k = self.g_k * (n ** 4) * (V - self.V_k)
        I_l = self.g_l * (V - self.V_l)
        return np.array([I_na, I_k, I_l]).T

    def run_HH_model(self):
        # Potential vector 1xlength
        self.Vs = np.zeros(self.N)

        # all matrices are 3xlenght, with the rows always being m - n - h
        self.gates = np.zeros((3, self.N))
        self.I_ions = np.zeros((3, self.N))
        
    	# initialize first voltage value and first gates
        self.Vs[0] = self.V0

        self.calculate_rate_equations(self.V0)

        steady_states = SteadyStates(V0=self.V0, N=1, T=self.T)
        steady_states.run_simulation()
        gates_init = steady_states.gates_inf
        self.gates[:, 0] = np.reshape(gates_init, (gates_init.shape[0], ))

    	## iterative calculation of the membrane potential
        for i in range(0, self.N-1):

    		# calculate ionic currents for current timestep
            self.I_ions[:, i] = self.calculate_ionic_currents(i, self.Vs[i])
	    
		# calculate membrane potential of a future timestep
            self.Vs[i+1] = self.Vs[i] + (1 / self.Cm) * (-np.sum(self.I_ions[:, i]) + self.I_stimulus[i]) * self.dt
	    
		# calculate gating variables of a future timestep
            self.calculate_rate_equations(self.Vs[i])
            self.gates[:, i+1] = self.calculate_gates(i)
