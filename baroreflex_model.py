import matplotlib.pyplot as plt
import julia
import time

st = time.time()
from julia import Main
Main.include('baroreflex_model.jl')


MyResult = Main.main_baro()
et = time.time()
print(et-st)
plt.plot(MyResult)
plt.show()
