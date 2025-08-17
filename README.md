This folder, which was written by Matt Larson, contains Sage code to do certain Euler characteristic computations in K-rings of wonderful varieties or matroids. 

Associated to a linear subspace L of k^n which is not contained in any coordinate hyperplane, there is a smooth projective variety W_L, obtained by blowing up PL at the coordinate subspaces in increasing order of dimension. The Grothendieck group of coherent sheaves on W_L, together with its Euler characteristic map, only depends on the (loopless) matroid associated to L. In the paper "K-rings of wonderful varieties and matroids," by Matt Larson, Shiyue Li, Sam Payne, and Nick Proudfoot, the authors introduce a combinatorially-defined analogue of the K-ring for an arbitrary loopless. 

The code in this folder computes the Euler characteristics of certain interesting classes. It uses the Matroid class in sage, and it computes them by doing a computation in the Chow ring of a matroid. 

WARNING: In order to make this computations work correctly, the ground set of the matroids involved must contain 0. 


For every matroid N on {0, ..., n}, the matroid polytope P(N) defines a line bundle on the n-dimensional permutohedral toric variety. Because the wonderful variety W_L is embedded in the permutohedral toric variety, this defines a line bundle on W_L. For example, if N = U_{n, n+1}, then the first Chern class of this line bundle is the class usually denoted by \beta in the Chow ring of a matroid. The construction of this line bundle can be extended to define a K-class for non-realizable matroids. 

The code can be used to compute Euler characteristics of tensor powers of these line bundles, using the function euler(M, N, num) to compute \chi(M, P(N)^{\otimes a}). The function which sends an integer a to \chi(M, P(N)^{\otimes a}) is a polynomial. The code is also able to compute this polynomial (using ehrhartPoly(M, N)) and the h^* polynomial of P(N) on M, which is a certain transformation of this polynomial. See Theorem 1.5 in the paper "K-theoretic positivity for matroids" by Chris Eur and Matt Larson for the definition.

For instance, Example 5.7 in "K-theoretic positivity for matroids" states that the h^*(M, U_{n, n+1}) is the h-vector of the broken circuit complex of M. The following code verified this for the Fano matroid.

A = Matrix(GF(2), [[1,0,0,1,1,0,1],[0,1,0,1,0,1,1],[0,0,1,0,1,1,1]])

M = Matroid(A)

N = matroids.Uniform(6,7)

hstarMN(M, N)

M.broken_circuit_complex().h_vector()


The code is also able to compute the Euler characteristic of a difference of two line bundles corresponding to matroid polytopes, i.e., if M, N, and P are matroids on the same ground set, then it can compute \chi(M, P(N) \otimes P(P)^{-1}). These Euler characteristics empirically display interesting vanishing behavior. 
