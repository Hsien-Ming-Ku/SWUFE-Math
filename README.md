# Test_matrices
Test matrices for (shifted) Krylov subspace methods.

Certainly, these matrices also can be used for investigating the performance 
of Krylov subspace method for A*x = b and the corresponding preconditioning 
techniques.

%------------------------------------------------------------------

!!! NOTE:

Some test matrices are used for investigating the performance of shifted Krylov subspace methods

(A - sigma_i *I) = b, sigma_i \in \mathbb{C}, i = 1,2,...,t.

These test matrices A are from the Numerical Linear Algebra (NLA) Course of Prof. Gerard L.G. Sleijpen 
(Mathematical Institute, Utrecht University, The Netherlands, Homepage: http://www.staff.science.uu.nl/~sleij101/, but recently the link of his NLA Course is often badly available). Moreover, some matrices come from two popular test matrix collections:

1. Matrix Market: http://math.nist.gov/MatrixMarket/matrices.html
2. The University of Florida Sparse Matrix Collection: http://www.cise.ufl.edu/research/sparse/matrices/

!!!! So please also ensure that the source for academic research only, not used for any commercial purposes.

!!!  I put them here just for the convenience of my research program.

% Good news: 
Dear everyone, recently I plan to post some frequently-used test matrices in my repositories. We named this collection as 
"UESTC-Math Matrix Library". I will devide them into "Non-Hermitian matrices", "Complex symmetric matrices", "Hermitian 
matrices", "Structured matrices" and "dense matrices". Then a User's guide will be given and updated. 

!!!!! If you intend to use the software, please do not forget to place the proper acknowledgements. Or maybe cite them as

X.-M. Gu, T.-Z. Huang, L. Li, Y. Zhang, H.-B. Li, B. Carpentieri, UESTC-Math Matrix Library based on MATLAB: Linear System and Eigenvalue Problem, Technical Report, version 1.0, School of Mathematical Sciences, University of Electronic Science and Technology of China, April, 2016. Available online at \url{https://github.com/Hsien-Ming-Ku/Test_matrices}.
