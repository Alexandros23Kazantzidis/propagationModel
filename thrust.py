import numpy as np


def compute_thrush_acc(system, m):

	if system == 1:
		mFuel = 1.3
		ve = 3000
		Iesp = 300
	elif system == 2:
		mFuel = 0.130
		ve = 3500
		Iesp = 350
	elif system == 3:
		mFuel = 0.003
		ve = 3500
		Iesp = 350
	elif system == 4:
		mFuel = 0.0000008
		ve = 25000
		Iesp = 2500

	accel = ((mFuel / m) * ve)

	return accel


if __name__ == "__main__":

	accel = compute_thrush_acc(1, 1000)
	print(accel)
	accel = compute_thrush_acc(2, 1000)
	print(accel)
	accel = compute_thrush_acc(3, 1000)
	print(accel)
	accel = compute_thrush_acc(4, 1000)
	print(accel)