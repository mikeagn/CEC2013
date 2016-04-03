#!/usr/bin/env python
###############################################################################
# Version: 1.1
# Last modified on: 3 April, 2016 
# Developers: Michael G. Epitropakis
#      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
###############################################################################
from cec2013.cec2013 import *
import numpy as np

def main():
	print 70*"="
	# Demonstration of all functions
	for i in range(1,21):
		# Create function
		f = CEC2013(i)

		# Create position vectors
		x = np.ones( f.get_dimension() )

		# Evaluate :-)
		value = f.evaluate(x)
		print "f", i, "(", x, ") = ", f.evaluate(x)

	print 70*"="
	# Demonstration of using how_many_goptima function
	for i in range(1,21):
		# Create function
		f = CEC2013(i)
		dim = f.get_dimension()

		# Create population of position vectors
		pop_size = 10
		X = np.zeros( (pop_size, dim) )
		ub =np.zeros( dim )
		lb =np.zeros( dim )
		# Get lower, upper bounds
		for k in range(dim):
			ub[k] = f.get_ubound(k)
			lb[k] = f.get_lbound(k)

		# Create population within bounds
		fitness = np.zeros( pop_size )
		for j in range(pop_size):
			X[j] = lb + (ub - lb) * np.random.rand( 1 , dim )
			fitness[j] = f.evaluate(X[j])

		# Calculate how many global optima are in the population
		accuracy = 0.001
		count, seeds = how_many_goptima(X, f, accuracy)
		print "In the current population there exist", count, "global optimizers."
		print "Global optimizers:", seeds

	print 70*"="

if __name__ == "__main__":
	main()
