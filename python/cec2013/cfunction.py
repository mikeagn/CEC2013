###############################################################################
# Version: 1.1
# Last modified on: 3 April, 2016 
# Developers: Michael G. Epitropakis
#      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
###############################################################################
import numpy as np
# UNCOMMENT APPROPRIATELY
#MINMAX = 1		# Minimization 
MINMAX = -1		# Maximization

class CFunction(object):
	__dim_ = -1
	__nofunc_ = -1
	__C_ = 2000.0
	__lambda_ = None
	__sigma_ = None
	__bias_ = None
	__O_ = None
	__M_ = None
	__weight_ = None
	__lbound_ = None
	__ubound_ = None
	__fi_ = None
	__z_ = None
	__f_bias_ = 0
	__fmaxi_ = None
	__tmpx_ = None
	__function_ = None

	def __init__(self, dim, nofunc):
		self.__dim_ = dim
		self.__nofunc_ = nofunc

	def evaluate(self, x): 
		pass
	
	def get_lbound(self, ivar):
		assert (ivar >= 0 and ivar < self.__dim_), ["ivar is not in valid variable range: %d not in [0,%d]" % ivar,self.__dim_]
		return self.__lbound_[ivar];

	def get_ubound(self, ivar):
		assert (ivar >= 0 and ivar < self.__dim_), ["ivar is not in valid variable range: %d not in [0,%d]" % ivar,self.__dim_]
		return self.__ubound_[ivar];

	def __evaluate_inner_(self, x):
		if self.__function_ == None:
			raise NameError('Composition functions\' dict is uninitialized')
		self.__fi_ = np.zeros(self.__nofunc_)

		self.__calculate_weights(x)
		for i in range(self.__nofunc_):
			self.__transform_to_z(x, i)
			self.__fi_[i] = self.__function_[i](self.__z_)

		tmpsum = np.zeros(self.__nofunc_)
		for i in range(self.__nofunc_):
			tmpsum[i] = self.__weight_[i] * (self.__C_ * self.__fi_[i] / self.__fmaxi_[i] + self.__bias_[i])
		return sum(tmpsum) * MINMAX + self.__f_bias_

	def __calculate_weights(self, x):
		self.__weight_ = np.zeros(self.__nofunc_)
		for i in range(self.__nofunc_):
			mysum = sum( (x-self.__O_[i])**2 )
			self.__weight_[i] = np.exp( -mysum/(2.0 * self.__dim_ * self.__sigma_[i] * self.__sigma_[i]) )
		maxw = np.max(self.__weight_)
		#maxi = self.__weight_.argmax(axis=0)
	
		maxw10 = maxw**10
		for i in range(self.__nofunc_):
			if self.__weight_[i] != maxw:
			#if i != maxi:
				self.__weight_[i] = self.__weight_[i] * (1.0 - maxw10)

		mysum = np.sum( self.__weight_ )
		for i in range(self.__nofunc_):
			if mysum == 0.0:
				self.__weight_[i] = 1.0 / (1.0 * self.__nofunc_ )
			else:
				self.__weight_[i] = self.__weight_[i] / mysum

	def __calculate_fmaxi(self):
		self.__fmaxi_ = np.zeros(self.__nofunc_)
		if self.__function_ == None:
			raise NameError('Composition functions\' dict is uninitialized')

		x5 = 5 * np.ones(self.__dim_)

		for i in range(self.__nofunc_):
			self.__transform_to_z_noshift(x5, i)
			self.__fmaxi_[i] = self.__function_[i](self.__z_)

	def __transform_to_z_noshift(self, x, index):
		# z_i = (x)/\lambda_i
		tmpx = np.divide(x, self.__lambda_[index])
		# Multiply z_i * M_i
		self.__z_ = np.dot( tmpx, self.__M_[index] )

	def __transform_to_z(self, x, index):
		# Calculate z_i = (x - o_i)/\lambda_i
		tmpx = np.divide((x - self.__O_[index]), self.__lambda_[index])
		# Multiply z_i * M_i
		self.__z_ = np.dot( tmpx, self.__M_[index] )

	def __load_rotmat(self, fname):
		self.__M_ = []

		with open(fname, 'r') as f:
			tmp = np.zeros( (self.__dim_,self.__dim_) )
			cline = 0
			ctmp = 0
			for line in f:
				line = line.split()
				if line:
					line = [float(i) for i in line]
					# re initialize array when reached dim
					if (ctmp%self.__dim_ == 0):
						tmp = np.zeros( (self.__dim_,self.__dim_) )
						ctmp = 0

					# add line to tmp
					tmp[ctmp] = line[:self.__dim_]
					# if we loaded self.__nofunc_ * self.__dim_-1 lines break
					if cline >= self.__nofunc_ * self.__dim_-1:
						break
					# add array to M_ when it is fully created
					if (cline%self.__dim_ == 0):
						self.__M_.append(tmp)
					ctmp = ctmp + 1
					cline = cline + 1


# Sphere function
def FSphere(x):
    return (x**2).sum()

# Rastrigin's function
def FRastrigin(x):
    return np.sum( x**2-10.*np.cos( 2.*np.pi*x )+10)

# Griewank's function
def FGrienwank(x):
    i = np.sqrt(np.arange(x.shape[0])+1.0)
    return np.sum(x**2)/4000.0 - np.prod(np.cos(x/i)) + 1.0

# Weierstrass's function
def FWeierstrass(x):
    alpha = 0.5
    beta = 3.0
    kmax = 20
    D = x.shape[0]
    exprf = 0.0

    c1 = alpha**np.arange(kmax+1)
    c2 = 2.0*np.pi*beta**np.arange(kmax+1)
    f = 0
    c = -D*np.sum(c1*np.cos(c2*0.5))
    
    for i in range(D):
        f += np.sum(c1*np.cos(c2*(x[i]+0.5)))
    return f + c

def F8F2(x):
    f2 = 100.0 * ( x[0]**2 - x[1] )**2 + (1.0 - x[0])**2
    return 1.0 + (f2**2)/4000.0 - np.cos(f2)

# FEF8F2 function
def FEF8F2(x):
    D = x.shape[0]
    f = 0
    for i in range(D-1):
        f += F8F2( x[ [i,i+1] ] + 1 )
    f += F8F2( x[ [D-1,0] ] + 1 )
    return f

