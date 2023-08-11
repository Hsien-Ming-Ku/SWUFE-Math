# SWUFE-Math Test Matrices Library
The library named "SWUFE-Math Test Matrices Library" will be established as another popular test matrix collection, which is different from some previous test matrix collections, such as Matrix Market (http://math.nist.gov/MatrixMarket/matrices.html) and The SuiteSparse Matrix Collection (formerly known as the University of Florida Sparse Matrix Collection, http://www.cise.ufl.edu/research/sparse/matrices/). Our new test matrix library will have following new features:

(1) Almost all of test matrices (including the corresponding right-hand side vectors) are from practical simulation applications; In particular, the size, physical/discretization parameters of various test problems can be flexibly modified on users' requirement;

(2) Collect more complex-valued test matrices, which are common in many engineering applications, such as Computational Electromagnetism and Structural Dynamics;

(3) Dense matrices (including some structured matrices) are particularly included in our test matrix library, the matrix data and the subrontinue for running the dense matrix-vector product are provided for the interested users;

(4) Many the matrix generator files are included and provide for users, then they can modified for their specific research;

(5) Some subrontinues for using classical preconditioning techniques (such as Jacobi, SSOR, ILU, etc) and recent Krylov subspace solvers 
(e.g., BiCOR, CORS, BiCORSTAB, GPBiCOR, (QMR-)COCR, etc) are also provied for interested users and students;

(6) Test matrices for novel class of matrix computations (such as eigenvalue problems, matrix function, shifted linear systems, linear systems with multiple right-hand sides, etc) are especially added;

At last but not least, we will also plan to give a systematic framework for making a relible numerical comparsion among different algorithms, and then users can conclude some useful information/suggestions for practical applications. We hope that our test matrix library can provide test problems for algorithm researchers, they can provide some robust algorithms for scientific and engineering 
applications, finally, researchers who work on real applications can provide us more new and difficult test problems.
 
              ----------------------         ---------------             --------------------
             | Numerical simultions | ----> | Test problems |   ------> | Algorithm research | --------->
              ----------------------         ---------------             --------------------           |
                       |<----------------------------------------------------------------------------<--
                       
Scientific and engineering applications:

1. Computational Electromagnetic, collaborate with researches from ANSYS Inc., UESTC, Curtin Univ., Univ. of Texas at Austin, INRIA, etc;

2. Computational fluid dynamics, collaborate with researches from Univ. of Groningen, TU Delft, TU Denmark, Univ. of British Columbia,  etc;

3. Electronic structure computations, collaborate with researches from Univ. of Tsukuba, Peking Univ., Univ. of Wisconsin Madison, etc;

4. Image processing, collaborate with researches from Xi'an Jiaotong Univ., HKBU, UCLA, TU Denmark, Univ. of Insubria, etc;

5. Numerical solutions of PDEs/FDEs, collaborate with Russian Academy of Sciences, Purdue Univ., UC Irvine, Univ. of Macau, Southeast Univ., etc.
                       
------------------------------------------------------------------------------------------------------------------------------------
!! Work/update notes: (please also ensure that the source for academic research only, not used for any commercial purposes)
 
1) In September 2016, our "UESTC-Math Matrix Library" has some sub-collections, i.e., "Non-Hermitian matrices", "Complex symmetric matrices", "Hermitian matrices", "Structured matrices" and "dense matrices". Then a User's guide will be given and updated. 

2) 4 June, 2016, I release some simple preconditioning techniques, who named "scaling2.m", "SSOR.m" and "indefiniteILU.m";

3) 10 July, 2016, I upload two file "Ren.tex"  and "How to prepare your test matrices.pdf" for users, who wants 
to provide their interesting test matrices for us;

4）At the end of 2017, we fail to release the first version of the User's Guide and test matrix data with ten sub-collections of test matrices;

5）Before the end of 2018, we will try our best to release the first version of the User's Guide and test matrix data with ten sub-collections of test matrices; (in progress)

!! If you intend to use the software, please do not forget to place the proper acknowledgements. Or maybe cite them as

X.-M. Gu, Y. Liu, B. Carpentieri, Y.-L. Zhao, SWUFE-Math Test Matrix Library based on MATLAB: Linear System and Eigenvalue Problem, Technical Report, School of Mathematics, Southwestern University of Finance and Economics, July 26, 2017. Available online at \url{https://github.com/Hsien-Ming-Ku/SWUFE-Math}.
