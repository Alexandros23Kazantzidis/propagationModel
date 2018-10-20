import numpy as np
import math as mp
import re


def find_drag_petrubation(Cd, mSat, dimensions, y):

	we = 7.292115e-5
	v_atm = we * np.array([-y[1], y[0], 0])
	satv = y[3:6]
	vrel = satv - v_atm
	rForDrag = (mp.sqrt(y[0]**2 + y[1]**2 + y[2]**2) - 6378137) / 1000
	A = 2 * (dimensions[0] * dimensions[1] + dimensions[0] * dimensions[2] + dimensions[1] * dimensions[2])

	## Find the density p of the atmosphere based on Harris-Priester atmospheric density model
	density = dict()
	density[300] = 35.26
	density[320] = 25.11
	density[340] = 18.19
	density[360] = 13.37
	density[380] = 9.955
	density[400] = 7.492
	density[420] = 5.684
	density[440] = 4.355
	density[460] = 3.362
	density[480] = 2.612
	density[500] = 2.042
	density[520] = 1.605
	density[540] = 1.267
	density[560] = 1.005
	density[580] = 0.7997
	density[600] = 0.6390
	density[620] = 0.5123
	density[640] = 0.4121
	density[660] = 0.3325
	density[680] = 0.2691
	density[700] = 0.2185
	density[720] = 0.1779
	density[740] = 0.1452
	density[760] = 0.1190
	density[780] = 0.09776
	density[800] = 0.08059
	density[840] = 0.05741
	density[880] = 0.04210
	density[920] = 0.03130
	density[960] = 0.02360
	density[1000] = 0.01810

	try:
		index = 0
		keepKey = []
		keepValue = []
		for key, value in density.items():
			keepValue.append(value)
			keepKey.append(key)
			if rForDrag < key:
				indexMin = index - 1
				indexMax = index
				break

			index += 1

		p = keepValue[indexMin] + (rForDrag - keepKey[indexMin]) * (keepValue[indexMax] - keepValue[indexMin]) / (keepKey[indexMax] - keepKey[indexMin])
	except:
		if rForDrag < 300:
			p = 35.26
		else:
			p = 0.018

	vrel = vrel / 1000
	A = A / 1000
	p = p / 1000
	accelerationDrag = -1/2 * Cd * (A / mSat) * p * vrel ** 2 * (vrel / np.linalg.norm(vrel))
	accelerationDrag = accelerationDrag * 1000

	return accelerationDrag


def drag(s):

	mu = 398600.4418  # gravitational parameter mu
	J2 = 1.08262668e-3  # J2 coefficient
	Re = 6378.137  # equatorial radius of the Earth
	we = 7.292115e-5  # rotation rate of the Earth in rad/s
	ee = 0.08181819  # eccentricity of the Earth's shape
	"""Returns the drag acceleration for a given state.
	   Args:
		   s(1x6 numpy array): the state vector [rx,ry,rz,vx,vy,vz]
	   Returns:
		   1x3 numpy array: the drag acceleration [ax,ay,az]
	"""

	r = np.linalg.norm(s[0:3])
	v_atm = we*np.array([-s[1],s[0],0])   # calculate velocity of atmosphere
	v_rel = s[3:6] - v_atm

	rs = Re*(1-(ee*s[2]/r)**2)   # calculate radius of surface
	h = r-rs
	p = 0.6*np.exp(-(h-175)*(29.4-0.012*h)/915) # in kg/km^3
	coeff = 3.36131e-9     # in km^2/kg
	acc = -p*coeff*np.linalg.norm(v_rel)*v_rel

	return acc


def doorbnos_drag(time, date, y, goceFileName, Cd):

	handler = open(goceFileName,"r")

	data = handler.read()

	dateIndex = [m.start() for m in re.finditer(date, data)]

	yearData = data[dateIndex[0]:dateIndex[-1]]

	timeIndex = [m.start() for m in re.finditer(time, yearData)][0]
	endOfRowIndex = [m.start() for m in re.finditer("\n", yearData[timeIndex:timeIndex+200])][0]

	timeRow = yearData[timeIndex:timeIndex+endOfRowIndex]

	firstSplit = timeRow.split(" ")

	properSplit = []
	for i in range(0, len(firstSplit)):
		if firstSplit[i] == "":
			pass
		else:
			properSplit.append(firstSplit[i])

	density = float(properSplit[7])

	dimensions = np.array([2, 2, 2])
	mSat = 1077
	we = 7.292115e-5
	v_atm = we * np.array([-y[1], y[0], 0])
	satv = y[3:6]
	vrel = satv - v_atm
	rForDrag = (mp.sqrt(y[0]**2 + y[1]**2 + y[2]**2) - 6378137) / 1000
	A = 2 * (dimensions[0] * dimensions[1] + dimensions[0] * dimensions[2] + dimensions[1] * dimensions[2])

	vrel = vrel
	A = A
	density = density
	accelerationDrag = -1/2 * Cd * (A / mSat) * density * vrel ** 2 * (vrel / np.linalg.norm(vrel))
	accelerationDrag = accelerationDrag * 1000

	print(accelerationDrag)

	return accelerationDrag


if __name__ == "__main__":

	# y = np.array([4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03])
	y = np.array([+6465819, -1429000, +14827, -206.5, -865.7, +7707.7])
	# accelerationDrag = find_drag_petrubation(2.1, 720, np.array([4.6, 2.34, 2.20]), y)
	accelerationDrag = find_drag_petrubation(2.1, 1077, np.array([2, 2, 2]), y)
	print(accelerationDrag)
	# y = y / 1000
	# acc = drag(y)
	# print(acc)

	y = np.array([+6465819, -1429000, +14827, -206.5, -865.7, +7707.7])
	date = "2012-11-20"
	time = "17:40:00.000"
	goceFileName = "goce_denswind_v1_5_2012-11.txt"
	Cd = 2.1
	accDrag = doorbnos_drag(time, date, y, goceFileName, Cd)
	print(accDrag)
