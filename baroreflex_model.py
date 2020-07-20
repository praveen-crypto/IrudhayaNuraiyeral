import os

import numpy as np
from julia import Main
from scipy.integrate import solve_ivp

import arterymodel as am

Main.include("baroreflex_model.jl")

def calc(HR):
    try:
        path = os.path.join('Datas', 'RLCtru.txt')
        rlctru = np.genfromtxt(path, delimiter=',')
        Rs = rlctru[:, 0]
        L = rlctru[:, 1]
        C = rlctru[:, 2]
        Rp = rlctru[:, 3]
        f = open(path, 'w')
        f.truncate()
        f.close()
        HR = int(HR)
        PF = 0
        dt = 0.001
        st = 0
        et = 10
        clock = np.arange(st, et, dt)
        T = 60 / HR
        t1 = clock / T
        t2 = t1 - np.floor(t1)
        t3 = np.multiply(T, t2)
        x1 = PF * (np.square(np.sin(3.14 * t3 / 0.3)))
        x2 = np.floor(t3 + 0.7)
        R1 = 0.11
        L = 0.011
        R2 = 1.11
        C = 0.91
        it = np.multiply((1 - x2), x1)
        # GENERATION OF PULSE
        # initial value
        #pulse_generator = lambda t, x: am.pulsegen(t, x, R1, L, R2, C, clock, it)  # input for solver function
        temp = Main.main_baro(HR)
        system_initial = np.zeros(256)
        system_finder = lambda t, x: am.integrated_ode(t, x, temp, clock)
        sol = solve_ivp(system_finder, [st, et], system_initial, method='Radau', t_eval=np.arange(st, et, dt))
        to = sol.t
        xo = sol.y
        t = to[:7500]
        x = xo[:, 2500:]
        return t, x
    except:
        return -1, -100000

if __name__ == "__main__":
    import time
    import Stenosis
    import matplotlib.pyplot as plt
    from julia import Main
    Main.include("baroreflex_model.jl")
    try:
        start = time.time()
        stn_dat = {'0': None, '1': 3, '7': 4, '13':6, '3': 7, '11': 20, '10': 0, '51': 0, '46': 0,
                '74': 0, '56': None, '70': 0, '62': 0, '63': 0, '108': 0, '109': 0, '102': 0, '107': 0, '96': 0, '92': 55}
        Stenosis.steno(0.04, 1.05, 0.6, **stn_dat)
        t, x = calc(75)
        #print("\n",x)
        end = time.time()
        print("Compilation time is: ", (end-start))
        plt.plot(t, x[132, :])
        plt.show()
    except Exception as e:
        print(str(e))
