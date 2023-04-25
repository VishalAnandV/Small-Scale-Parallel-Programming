## Small Scale Parallel Programming

This project describes the design and development of a sparse matrixvector product kernel. The implemented kernel is capable of computing y ← Ax, where A is a sparse matrix stored in CSR and ELLPACK storage formats. In order to utilise the computing capabilities of Crescent (Super Computer @ Cranfield University), the kernel is parallelized in both OpenMP and CUDA versions and the code for both formats is developed in C language. The correctness of the results is tested against a serial implementation reference. Necessary pre-processing of the matrix market data is done and converted into desired formats. The performance for given set of matrices is computed, and analysed. Also, the performance data of OpenMP and CUDA version are compared and analysed.
 

