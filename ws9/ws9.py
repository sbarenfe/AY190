import numpy as np
import matplotlib.pyplot as plt

def two():
	for i in range(1,6):
		A,b=read(i)
		for j in range(len(b)):
			if A[j][j]==0:
				print "Problem with:"
				print A[j]
		x=GENP(A,b)
		print "x equals:"
		print x

def three():
	for i in range(1,6):
		A,b=read(i)
		for j in range(len(b)):
			if A[j][j]==0:
				print "Problem with:"
				print A[j]
		x=np.linalg.solve(A,b)
		print "x equals:"
		print x

def read(n):
	Afile='LSE%i_m.dat' %n
	bfile='LSE%i_bvec.dat' %n
	A=np.loadtxt(Afile)
	b=np.loadtxt(bfile)
	print np.shape(A)
	print len(b)
	print np.linalg.slogdet(A)
	return A,b

def GENP(A, b):
    '''
    Gaussian elimination with no pivoting.
    % input: A is an n x n nonsingular matrix
    %        b is an n x 1 vector
    % output: x is the solution of Ax=b.
    % post-condition: A and b have been modified. 
    '''
    A=np.float64(A)
    b=np.float64(b)
    n =  len(A)
    if b.size != n:
        raise ValueError("Invalid argument: incompatible sizes between A & b.", b.size, n)
    for pivot_row in xrange(n-1):
        for row in xrange(pivot_row+1, n):
            multiplier = A[row][pivot_row]/A[pivot_row][pivot_row]
            #the only one in this column since the rest are zero
            A[row][pivot_row] = 0 #my change to 0
            for col in xrange(pivot_row + 1, n):
                A[row][col] = A[row][col] - multiplier*A[pivot_row][col]
            #Equation solution column
            b[row] = b[row] - multiplier*b[pivot_row]
    x = np.zeros(n)
    k = n-1
    x[k] = b[k]/A[k,k]
    k=k-1
    while k >= 0:
        x[k] = (b[k] - np.dot(A[k,k+1:],x[k+1:]))/A[k,k]
        k = k-1
    return x
