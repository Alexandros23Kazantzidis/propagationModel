import numpy as np
import matplotlib.pyplot as plt
import re
import pickle
from scipy import special as sp
from tabulate import tabulate
import math as mp
from sun_moon import find_sun_moon_acc


def restruct_geopotential_file(gravityModelFileName):
	"""Creates a new restructured file holding the geopotential coefficients Cnm,Snm of a gravity model
		Args:
			gravityModelFileName (string): the name of the file holding the coeffs
	"""

	handler = open(gravityModelFileName, 'r')
	writeHandler = open("restruct_" + gravityModelFileName, "a")
	for line in handler:
		reStructLine = re.sub("\s+", ",", line.strip())
		writeHandler.write(reStructLine)
		writeHandler.write("\n")

	handler.close()
	writeHandler.close()


def read_geopotential_coeffs(restructGravityFileName,firstRead):
	"""Returns an array holding the geopotential coefficients Cnm,Snm of a gravity model
		Args:
			gravityModelFileName (string): the name of the file holding the coeffs
			firstRead (boolean): if you are using this method the first time or the .p file has been created already
		Returns:
			nxm numpy array: coeffs in proper format to be used by ydot_geopotential
	"""

	if firstRead is True:

		data = np.genfromtxt(restructGravityFileName, delimiter=",")
		data = data[1:, 1:]

		pickle.dump(data, open("data.p", "wb"))
	elif firstRead is False:

		data = pickle.load(open("data.p", "rb"))

	return data


def cart2geodetic(x, y, z):
	"""cartesian to geodetic coordinates transformation ellipsoid used WGS84
		Args:
			x (float)
			y (float)
			z (float)
		Returns:
			r (float)
			phi (float)
			l (float)
	"""

	a = 6378137
	f = 1 / 298.257223563
	b = a - f*a
	e = mp.sqrt(2*f - f**2)
	epsilon = mp.sqrt((a**2 - b**2) / b**2)
	p = mp.sqrt(x**2 + y**2)
	theta = mp.atan2((z*a), (p*b))

	phi = mp.atan2((z + (epsilon**2)*b*mp.sin(theta)**3), (p - (e**2)*a*mp.cos(theta)**3))
	l = mp.atan2(y, x)
	r = p / mp.cos(phi) - (a / mp.sqrt(1-(e**2)*mp.sin(phi)**2))

	return r, phi, l


def createNumpyArrayForCoeffs(gravityCoeff, n, m):
	"""
		Return Cnm and Snm in separate numpy arrays for certain degree n and order m
	"""
	C = np.zeros((n, m))
	S = np.zeros((n, m))

	for i in range(0, len(gravityCoeff)):
		nIndex = gravityCoeff[i, 0]
		mIndex = gravityCoeff[i, 1]

		if (nIndex > n-1) or (mIndex > m-1):
			continue
		else:
			C[int(nIndex), int(mIndex)] = gravityCoeff[i, 2]
			S[int(nIndex), int(mIndex)] = gravityCoeff[i, 3]

	C[0, 0] = 1
	return C, S


def ydot_geopotential(y, C, S, n, m):
	"""Returns the time derivative of only the geopotential effect a given state.
		Args:
			y(1x6 numpy array): the state vector [rx,ry,rz,vx,vy,vz]
			gravityCoeff (nxm numpy array): the array where you have stored the geopotential coefficients Cnm,Snm
			n (integer): degree of coeffs you want to use for the calculation
			m (integer): order of coeffs you want to use for the calculation
		Returns:
			1x6 numpy array: the time derivative of y [vx,vy,vz,ax,ay,az]
	"""

	r, phi, l = cart2geodetic(y[0], y[1], y[2])

	earthR = 6378137
	mu = 398600.4405e+09

	r = r + earthR

	legendreP = sp.lpmn(n+1, m+1, mp.sin(phi))
	legendreP = legendreP[0]
	legendreP = np.transpose(legendreP)
	V = np.zeros((n+1, m+1))
	W = np.zeros((n+1, m+1))

	for i in range(0, n+1):
		for j in range(0, i+1):

			V[i, j] = ((earthR/r) ** (i+1)) * (legendreP[i, j] * mp.cos(j * l))
			W[i, j] = ((earthR/r) ** (i+1)) * (legendreP[i, j] * mp.sin(j * l))

	ax = 0
	ay = 0
	az = 0
	for i in range(0, n-2):
		for j in range(0, i+1):

			if j == 0:
				ax = ax + (mu / earthR**2) * (-C[i, 0] * V[i+1, 1])

				ay = ay + (mu / earthR ** 2) * (-C[i, 0] * W[i + 1, 1])
			else:
				ax = ax + (mu / earthR ** 2) * 0.5 * (-C[i, j] * V[i + 1, j + 1] - S[i, j] * W[i + 1, j + 1]) + \
					 (mp.factorial(i - j + 2) / mp.factorial(i - j)) * (
					 C[i, j] * V[i + 1, j - 1] + S[i, j] * W[i + 1, j - 1])

				ay = ay + (mu / earthR ** 2) * 0.5 * (-C[i, j] * W[i + 1, j + 1] + S[i, j] * V[i + 1, j + 1]) + \
				     (mp.factorial(i - j + 2) / mp.factorial(i - j)) * (
					 -C[i, j] * W[i + 1, j - 1] + S[i, j] * V[i + 1, j - 1])

			az = az + (mu / earthR ** 2) * ((i-j+1) * (-C[i, j] * V[i+1, j] - S[i, j] * W[i+1, j]))

	ax = -ax
	ay = -ay
	# print(ax, ay, az)
	return np.array([ax, ay, az])


def ydot(y, C, S, n, m):
	"""Returns the time derivative of a given state.
		Args:
			y(1x6 numpy array): the state vector [rx,ry,rz,vx,vy,vz]
		Returns:
			1x6 numpy array: the time derivative of y [vx,vy,vz,ax,ay,az]
	"""

	mu = 398600.4405e+09
	r = np.linalg.norm(y[0:3])
	# a = -mu/(r**3)*y[0:3]

	a = ydot_geopotential(y, C, S, n, m)

	accelerationSun, accelerationMoon = find_sun_moon_acc(y)

	a = a + accelerationSun + accelerationMoon

	# p_j2 = j2_pert(y)
	# p_drag = drag(y)
	#
	# a = a+p_j2+p_drag
	# print(a)
	return np.array([*y[3:6], *a])


def rk4(y, t0, tf, h=5):
	"""Runge-Kutta 4th Order Numerical Integrator
		Args:
			y(1x6 numpy array): the state vector [rx,ry,rz,vx,vy,vz]
			t0(float)  : initial time
			tf(float)  : final time
			h(float)   : step-size
		Returns:
			1x6 numpy array: the state at time tf
	"""

	t = t0

	if tf < t0:
		h = -h

	while(abs(tf-t) > 0.00001):
		if (abs(tf-t) < abs(h)):
			h = tf-t

		# k1 = h*ydot(y)
		# k2 = h*ydot(y+k1/2)
		# k3 = h*ydot(y+k2/2)
		# k4 = h*ydot(y+k3)

		k1 = h * ydot(y, C, S, 10, 10)
		k2 = h * ydot(y + k1 / 2, C, S, 10, 10)
		k3 = h * ydot(y + k2 / 2, C, S, 10, 10)
		k4 = h * ydot(y + k3, C, S, 10, 10)

		y = y+(k1+2*k2+2*k3+k4)/6
		t = t+h

	return y


if __name__ == "__main__":
	y = np.array([4.57158479e+06, -5.42842773e+06, 1.49451936e+04, -2.11034321e+02, -1.61886788e+02, 7.48942330e+03])
	t0, tf = 0, 100.00
	final = np.zeros((100, 6))

	# restruct_geopotential_file("EGM2008.gfc")
	data = read_geopotential_coeffs("restruct_EGM2008.gfc", False)

	C, S = createNumpyArrayForCoeffs(data, 10, 10)

	# ydot_geopotential(y, C, S, 10, 10)
	# ydot(y)

	for i in range(0, 100):
		final[i, :] = rk4(y, t0, tf, 2)
		t0 = tf
		tf = tf + 100
		y = final[i, :]
		print(y)

	plt.plot(final[:, 1])
	plt.show()

