from matplotlib import pyplot as plt
from matplotlib import animation
import numpy as np
import csv
import pandas as pd
# read in the observables data
data = pd.read_csv("data.csv")
# now read in the data from the lattice files
latt_file = open("latt.txt", 'r')
lines = latt_file.readlines()

# create a list to store each 2d array in
lattices = []
curr = []
# go through the file one by one
for line in lines:
    # if its a newline, we're done with the current lattice and can start another one
    if line == "\n":
        lattices.append(curr)
        curr = []
    # otherwise add the line to the lattice after processing
    else:
        # split the line into a list of integers
        # this maps the fn that converts from string to int to each string
        str_temp = line.split(',')
        # get rid of the last element
        str_temp = str_temp[:-1]
        curr.append(list(map(int, str_temp)))

# plotting observables
plt.subplot(1,2,1)
plt.scatter(data["Temperature"], data["Average Energy"])
plt.subplot(1,2,2)
plt.scatter(data["Temperature"], data["Average Magnetization"])


plt.show()

fig = plt.figure()

dim = len(lattices[0][0])
a = lattices[0]
im = plt.imshow(a, interpolation='none', cmap = plt.cm.seismic)

def init():
    im.set_data(lattices[0])
    return [im]

def animate(i):
    im.set_array(lattices[i])
    return [im]

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=len(lattices), interval=100, blit=True)
plt.show()
