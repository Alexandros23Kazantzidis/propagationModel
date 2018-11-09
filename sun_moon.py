from jplephem.spk import SPK
import julian
import datetime
import numpy as np
import math as mp


def position_sun_moon(bspFileName):

	jd = julian.to_jd(datetime.datetime.now(), fmt='jd')
	jd = 2458355.500000

	kernel = SPK.open(bspFileName)
	positionMoon = kernel[3, 301].compute(jd)

	positionSun = kernel[0, 10].compute(jd)
	positionSun -= kernel[0, 3].compute(jd)
	positionSun -= kernel[3, 399].compute(jd)

	positionSun = positionSun * 1000
	positionMoon = positionMoon * 1000

	return positionSun, positionMoon


def find_sun_moon_acc(y, positionSun, positionMoon):

	G = 6.674e-11
	sunM = 1.989e+30
	moonM = 7.34767309e+22

	accelerationSun = - ((G * sunM) / positionSun ** 3) * y[0:3]
	accelerationMoon = - ((G * moonM) / positionMoon ** 3) * y[0:3]

	return accelerationSun, accelerationMoon


if __name__ == "__main__":

	y = np.array([4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03])
	bspFileName = 'de430.bsp'
	positionSun, positionMoon = position_sun_moon(bspFileName)
	accelerationSun, accelerationMoon = find_sun_moon_acc(y, positionSun, positionMoon)
	print(accelerationSun, accelerationMoon)

