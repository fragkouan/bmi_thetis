"""
Script to run the Thetis model
"""
from bmithetis import BmiThetis
from thetis import *
import time
from datetime import timedelta
from math import ceil
import numpy as np

# Start time
t1 = time.time()

# Initialize Thetis
thetis = BmiThetis()
thetis.initialize('./inputs/thetis_01_configuration.json')

# Calculate the times of the Thetis simulation. The
times = np.arange(
    thetis.t_start + thetis.dt,
    thetis.t_end + thetis.dt,
    thetis.dt
)

# If you want to run Thetis "in segments", uncomment
# Loop through the Thetis timesteps
# for t in times:
#     thetis.update_until(t)

# Advance Thetis till the specified time
thetis.update_until(thetis.t_end)



# Run the model until the specified end time
# thetis.update_until(thetis.t_end)

# Finalise model
thetis.finalize()
# End time
t2 = time.time()

# Print time duration
td_str =str(timedelta(seconds=t2-t1))
x = td_str.split(':')
x[2] = ceil(float(x[2]))
print("Simulation duration : ",x[0],"Hours",x[1],"Minutes",x[2],"Seconds")

