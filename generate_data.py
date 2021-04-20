import math
import numpy as np
import random

if __name__ == "__main__":
   grid_size = 500
   nqty = 5000

   xgrid = np.linspace(0, 50, num=grid_size)
   with open("x_grid.dat", 'w') as f:
       for i in xgrid:
           f.write(str(i)+"\n")

   with open("fx.dat", 'w') as f:
       for i in range(nqty):
           for j in range(grid_size):
               f.write(str(random.uniform(-500,500))+"\n")
               
