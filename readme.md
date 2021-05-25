2D Uniform Black Box FMM
========================

This file is best read on a markdown reader; preferably the app named Mou on Mac OS X; Atom or GitBook or other such markdown editors/readers on Linux.

This is a Black box FMM on a uniform tree, mainly for testing the speed of the FMM. You would need to have g++ (preferably g++-6, else change it in Makefile2D.mk) and the linear algebra library `Eigen` available @ <http://eigen.tuxfamily.org> installed to run this code.

To run the code, you would first need to make the Makefile. Currently, the code is optimized for the kernels 1/R and log(R)

	make -f Makefile2D.mk

Then we need to run the executable, i.e.,

	./testFMM2D nLevels nNodes L

where `nLevels` is a positive integer intdicating the number of levels in the tree, `nNodes` is a positive integer indicating the number of Chebyshev nodes along one dimension and `L` indicates the semi-length of the simulation box which is from [-L,L]^2. The number of particles would be `4^(nLevels)*nNodes^2`.

For instance, below is a sample execution of the code and its output.

It is always good to clean using make clean before running the code, i.e.,
	
	make -f Makefile2D.mk clean

Then make the file

	make -f Makefile2D.mk

Run the generated executable as, for instance,

	./testFMM2D 10 6 1

The generated output will be as follows

	Number of particles is: 37748736

	Time taken to create the tree is: 0.556001

	Time taken to assemble the operators is: 0.00377488

	Time taken to assemble the charges is: 0.52678

	Time taken for multipole to multipole is: 0.189195

	Time taken for multipole to local is: 4.61744

	Time taken for local to local is: 0.189695

	Time taken for self and neighbors at leaf is: 1.1813

	Total time taken is: 7.26796

	Apply time taken is: 6.17763

	Total Speed in particles per second is: 5.19386e+06

	Apply Speed in particles per second is: 6.11055e+06

	Number of particles is: 37748736

	Performing Error check...

	Box number is: 360061

	Box center is: (0.799805, -0.768555);

	Error is: 2.00299e-07