from julia import Main
import matplotlib.pyplot as plt
import time

st = time.time()
Main.include("cardiac.jl")
res = Main.cardiac()
et = time.time()

print("total time:", (et-st))
plt.plot(res)
plt.show()



