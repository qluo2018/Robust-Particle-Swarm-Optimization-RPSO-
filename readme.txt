This is a simple Matlab algorithm for RPSO written by Qiang Luo, PhD student in Department of Mathematics and Systems Science, National University of Defense Technology, China. It is free for any academic users, but be aware of that there is no guarantee of bug free. Any discussion about the theory or the application of this algorithm are warmly welcomed. If you use this code, please kindly cite the following paper:

Qiang Luo, Dongyun Yi, A Co-evolving Framework for Robust Particle Swarm Optimization£¬Applied Mathematics and Computation, 2008, 199(2):611-622.


--------------------------------------------------------------
How to use?

The core algorithms are coded in PSOed200.m which minimizes a given objective function.
Four algorihms (SPSO, LPSO, FPSP, RPSO) have been developed in this programe which can be called by specifying different AlgPara.
Many objective functions have been implemented in myFun.m. If you want to include your own objective function in the programe, please simply try to add the code of your function in myFun.m as another 'case' of 'switch'.
The main algorithm is in runPSO.m which shows how to call the function of PSOed200. Simply call the PSOed200.m with different parameters.

---------------------------------------------------------------

Corresponding email: luoqiang@nudt.edu.cn


