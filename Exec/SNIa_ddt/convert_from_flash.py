# This routine converts a FLASH-prepared 1D SN progenitor model
# for the DDT problem into a CASTRO-compatible format.

import numpy as np

flash_filename = "flash.dat"

# Import the data as a numpy array using multiple columns

data = np.loadtxt(flash_filename, skiprows=2)

# Get the number of points from the first non-comment line

data_file = open(flash_filename, 'r')

first_line = data_file.readline()
npts = int(data_file.readline())

data_file.close()

# Now store radius, density, and temperature.

radius = data[:,0]
density = data[:,1]
temperature = data[:,2]

# We don't get pressure, but we don't need it anyway,
# so just set it to zero.

pressure = np.zeros(npts)

# Store composition

c12  = data[:,3]
ne22 = data[:,4]

# Oxygen is just 1 - C - Ne

o16 = 1.0 - c12 - ne22

output_filename = "wd_hse.dat"

output_file = open(output_filename, 'w')

output_file.write("# npts = " + str(npts) + "\n")
output_file.write("# num of variables = " + "6" + "\n")
output_file.write("# density\n")
output_file.write("# temperature\n")
output_file.write("# pressure\n")
output_file.write("# carbon-12\n")
output_file.write("# oxygen-16\n")
output_file.write("# neon-22\n")

for n in range(0,int(npts)):

    R    = radius[n]
    rho  = density[n]
    P    = pressure[n]
    T    = temperature[n]
    XC12 = c12[n]
    XO16 = o16[n]
    XNe22 = ne22[n]
    output_file.write(str(R) + " " + str(rho) + " " + str(T) + " " + str(P) + " " + \
                      str(XC12) + " " + str(XO16) + " " + str(XNe22) + "\n")

output_file.close()
