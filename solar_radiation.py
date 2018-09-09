import numpy as np
import math as mp
from sun_moon import position_sun_moon


def get_solar_radiation_petrubation(sunr, satr, epsilon, mSat, dimensions):

	earthR = 6378137
	sunR = 695508000
	# s = satr
	# a = mp.asin(sunR / np.linalg.norm(sunr - satr))
	# b = mp.asin(earthR / np.linalg.norm(s))
	# c = mp.acos((np.dot(-np.transpose(s), (sunr - satr))) / (np.linalg.norm(s) * np.linalg.norm(sunr - satr)))
	#
	# x = (c**2 + a**2 - b**2) / 2*c
	# y = mp.sqrt(a**2 - x**2)
	# A = a**2 * mp.acos(x, a) + b**2 * mp.acos(c - x, b) - c*y

	##Simple formula
	A = 2 * (dimensions[0] * dimensions[1] + dimensions[0] * dimensions[2] + dimensions[1] * dimensions[2])
	Po = 4.56e-06
	Cr = 1 + epsilon
	Au = 149597870700

	accelerationRadiation = -Po * Cr * (A / mSat) * (sunr / sunr**3) * Au**2

	return accelerationRadiation


if __name__ == "__main__":

	y = np.array([4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03])
	bspFileName = 'de430.bsp'
	positionSun, positionMoon = position_sun_moon(bspFileName)
	accelerationRadiation = get_solar_radiation_petrubation(positionSun, y[0:3], 0.5, 720, np.array([4.6, 2.34, 2.20]))


