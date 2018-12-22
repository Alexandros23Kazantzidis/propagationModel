from sgp4.earth_gravity import wgs72
from sgp4.io import twoline2rv
from cowell import propagator
import numpy as np
import datetime


line1 = ('1 36508U 10013A   18248.50316386  .00000076  00000-0  17991-4 0  9996')
line2 = ('2 36508  92.0318 310.1069 0004900 265.0518  95.0130 14.52173178445782')


satellite=twoline2rv(line1,line2,wgs72)
position,velocity=satellite.propagate(
     2018,9,5,12,5,33)

print(satellite.epoch)    #nonzeroonerror

print(satellite.error_message)

print(position)

print(velocity)

y = np.array([4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03])
t0, tf = 0, 60.00
final = np.zeros((200, 6))
final_2 = np.zeros((200, 6))
diffs = []

propagatorObj = propagator()
data = propagatorObj.read_geopotential_coeffs("restruct_EGM2008.gfc", False)
propagatorObj.createNumpyArrayForCoeffs(data, 10, 10)

final[0, :] = propagatorObj.rk4(y, t0, tf, 1, 2, False, 0.5, False, False, 2.1, 1, propagatorObj.C, propagatorObj.S, 720, "2012-11-20", "17:40:00", datetime.datetime.now())

print(final)
