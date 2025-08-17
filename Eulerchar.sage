import time
from itertools import combinations_with_replacement
from itertools import combinations
from sympy.combinatorics import Permutation
from sage.combinat.combinat import eulerian_polynomial
#ground set should be numbers rather than letters
#otherwise there is an issue with sorting, flats needs to be ordered correctly to be coerced into the Chow ring
#additionally, the ground set should contain 0

#
#
#functions to create classes in A(M)
#
#

#in the chow ring of some matroid, creates the element corresponding to a flat F
def flatToChow(F, A):
	s = 'A'
	flatlist = list(F)
	flatlist.sort()
	for i in flatlist:
		s = s + str(i)
	return A(s)

#creates the divisor alpha
def alpha(A, M):
	alpha = A(0)
	for i in range(1, M.rank()):
		for F in M.flats(i):
			if (0 in F):
				alpha = alpha + flatToChow(F, A)
	return alpha


#takes a matroid N and a matroid M, together with chow ring A
#computes the Chern roots of S_N restricted to A(M)
#Even though there is a dual, this does indeed compute the Chern roots of S_N
def chernRootSN(M, A, N):
	N.dual()
	chernlist = []
	a = alpha(A, M)
	for i in range(N.rank()):
		chern = a #case when S = E
		for j in range(1, M.rank()):
			for F in M.flats(j):
				if (N.rank(F) > i):
					chern = chern - flatToChow(F, A)
		chernlist.append(chern)
	return chernlist

#takes a matroid N
#returns the first Chern class of the line bundle corresponding to P(N), in A(M)
def matroidToChow(M, A, N):
	Ndual = N.dual()
	a = alpha(A, M)
	c1 = -(Ndual.rank())*a
	for i in range(1, M.rank()):
		for F in M.flats(i):
			c1 += Ndual.rank(F)*flatToChow(F, A)
	return c1



#takes a matroid N and a matroid M, together with a chow ring A
#computes Chern roots of Q_N restrict to A(M)
def chernRootQN(M, A, N):
	trunclasses = []
	Ndual = N.dual()
	for i in range(N.size() - N.rank() + 1):
		trunclasses.append(matroidToChow(M, A, Ndual))
		Ndual = Ndual.truncation()
	chernlist = []
	for i in range(N.size() - N.rank()):
		chernlist.append(trunclasses[i] - trunclasses[i+1])
	return chernlist


#takes a list of Chern roots
#returns the Chern classes
def chernClasses(L):
	chernclass = [1]
	for i in range(1, len(L) + 1):
		ck = 0
		for term in combinations(L, i):
			prod = 1
			for i in term:
				prod = prod*i
			ck += prod
		chernclass.append(ck)
	return chernclass



#
#
#
#functions which compute things related to polynomials
#
#
#
R = PolynomialRing(QQ, 't')
var('t')

#takes a polynomial P in variable t
#returns the h^* polynomial, i.e., the polynomial satisfying \sum_{i \ge 0}P(i) t^i = h^*(t)/(1 -t)^(d + 1)
#d = deg(P)
def ehrhartToHstar(P):
	P = R(P)
	coef = P.list()
	d = P.degree()
	ans = (1-t)^d*coef[0]
	for i in range(1, len(coef)):
		ans += coef[i]*t*(1-t)^(d-i)*eulerian_polynomial(i)
	return ans.list()

#takes the h-vector of a (d-1)-dimensional complex (may have trailing zeros)
#returns the f-vector (f_{-1}, f_0, f_1, ..., f_{d-1})
def htofvector(L, d):
	for i in range(d + 1 - len(L)):
		L.append(0)
	fvec = [1]
	for i in range(1, d+1):
		termi = 0
		for k in range(0, i+1):
			termi += binomial(d-k, i-k)*L[k]
		fvec.append(termi)
	return fvec

#takes the h-vector
#returns the f-vector (may have trailing zeros)
def ftohvector(L):
	hvec = []
	d = len(L) - 1
	for k in range(len(L)):
		hk = 0
		for i in range(k + 1):
			hk += (-1)^(k-i)*binomial(d - i, k-i)*L[i]
		hvec.append(hk)
	return hvec


#
#
#functions to compute Euler characteristics
#
#




#computes \chi(M, [P(N)]^a)
#uses the HRR-type formula to write this as deg_M((1 + alpha + ...)c(S_N^{\vee})^a)
def euler(M, N,num):
	A = M.chow_ring(QQ)
	N = N.dual()
	chernSList = chernClasses(chernRootSN(M, A, N))
	totalclass = 0
	for i in range(len(chernSList)):
		totalclass += (-1)^i*chernSList[i]
	#creates 1 + alpha + ..
	a = alpha(A, M)
	totalalpha = 1
	for i in range(1, M.rank()):
		totalalpha += a^i
	deg1 = a^(M.rank() - 1)
	prod = totalclass^num*totalalpha
	#adds something to make sure there is something in top degree
	prodplus = prod + 1/2*deg1
	return QQ(prodplus.lc()/deg1.lc()) - 1/2

#Compute the h^* polynomial of P(N) in A(M)
def hstarMN(M, N):
	A = M.chow_ring(QQ)
	points = [[0, 1]]
	for i in range(1, M.rank()):
		points.append([i, euler(M, N, i)])
	ehr = R.lagrange_polynomial(points)
	degans = ehr.degree()
	ans = ehrhartToHstar(ehr)
	for i in range(degans + 1 - len(ans)):
		ans.append(0)
	return ans

#Compute the Ehrhart polynomial of P(N) on M
def ehrhartPoly(M, N):
	A = M.chow_ring(QQ)
	points = [[0, 1]]
	for i in range(1, M.rank()):
		points.append([i, euler(M, N, A, i)])
	ehr = R.lagrange_polynomial(points)
	return ehr




#takes 2 matroids, N, P
#and also M
#computes \chi(M, P(N) \otimes P(P)^{-1})
def eulerMatroidPolytopeDiff(M, N, P):
	A = M.chow_ring(QQ)
	a = alpha(A, M)
	totalalpha = 1
	for i in range(1, M.rank()):
		totalalpha += a^i
	deg1 = a^(M.rank() - 1)
	Ndual = N.dual()
	Pdual = P.dual()
	chernSList = chernClasses(chernRootSN(M, A, Ndual))
	chernQList = chernClasses(chernRootQN(M, A, Pdual))
	totalclassS = 0
	for i in range(len(chernSList)):
		totalclassS += (-1)^i*chernSList[i]
	totalclassQ = 0
	for i in range(len(chernQList)):
		totalclassQ += (-1)^i*chernQList[i]
	prod = totalclassS*totalclassQ*totalalpha
	#print(prod)
	#need to add something 
	prodplus = prod + 1/2*deg1
	return QQ(prodplus.lc()/deg1.lc()) - 1/2





