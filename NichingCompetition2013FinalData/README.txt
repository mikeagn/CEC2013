For each algorithm there are two files with the PR,SR values produced for the
competition:

- Each <algorithm name>_PR.dat contains the Average Peak Ratio values for each
problem and each accuracy level.
- Each <algorithm name>_SR.dat contains the Average Success Rate values for each
problem and each accuracy level.

The structure of the files contain a tab-separated matrix (20 x 5, 20 problems
x 5 accuracy levels) e.g.:

for all problems and each accuracy level: IMAGINE_PR.dat

1.0    1.0    1.0    1.0    1.0
1.0    0.999    0.998    0.997    0.996
...
...


Similarly for the IMAGINE_SR.dat:

1.0    1.0    1.0    1.0    1.0
1.0    0.99    0.99    0.99    0.99
...
... 

If you are using this material please acknowledge our help and cite the
competition's technical report:

X. Li, A. Engelbrecht, and M.G. Epitropakis, ``Benchmark Functions for
CEC'2013 Special Session and Competition on Niching Methods for Multimodal
Function Optimization'', Technical Report, Evolutionary Computation and
Machine Learning Group, RMIT University, Australia, 2013.

In addition, please properly cite the algorithms:

For the dADE/nrand/1,2 please properly cite: 

M. G. Epitropakis, Li, X., and Burke, E. K., “A Dynamic Archive Niching
Differential Evolution Algorithm for Multimodal Optimization”, IEEE Congress
on Evolutionary Computation, 2013. CEC 2013. Cancun, Mexico, pp. 79-86, 2013

For the DE/nrand/1,2 please properly cite:

M. G. Epitropakis, V. P. Plagianakos, and M. N. Vrahatis, "Finding multiple
global optima exploiting differential evolution’s niching capability," in 2011
IEEE Symposium on Differential Evolution (SDE), April 2011, pp. 1-8.

For the NEA1,2 please properly cite:

 M. Preuss. "Niching the CMA-ES via nearest-better clustering." In
Proceedings of the 12th annual conference companion on Genetic and evolutionary
computation (GECCO ’10). ACM, New York, NY, USA, pp. 1711-1718, 2010.

For the CrowdingDE (CDE) please properly cite:

R. Thomsen, "Multimodal optimization using crowding-based differential
evolution," In the IEEE Congress on Evolutionary Computation, 2004. CEC2004,
vol.2, pp. 1382-1389, 19-23 June, 2004

For the DECG, DELG, DELS_ajitter please properly cite:

J. Ronkkonen, Continuous Multimodal Global Optimization with Differential
Evolution-Based Methods, Ph.D. thesis, Lappeenranta University of Technology,
2009

For the CMA-ES please properly cite:

N. Hansen and A. Ostermeier (2001). Completely Derandomized
Self-Adaptation in Evolution Strategies. Evolutionary Computation, 9(2), pp.
159-195;

For the iPOP-CMA-ES please properly cite:

A. Auger and N. Hansen, "A restart CMA evolution strategy with increasing
population size," In the 2005 IEEE Congress on Evolutionary Computation,
2005. vol.2, pp.1769-1776, 2-5 Sept. 2005

For the ANSGAII and PNANSGAII please properly cite:

A parameterless-niching-assisted bi-objective approach to multimodal
optimization, http://dx.doi.org/10.1109/CEC.2013.6557558

For the NVMO (Molina) please properly cite: 

Variable mesh optimization for the 2013 CEC Special Session Niching Methods
for Multimodal Optimization, http://dx.doi.org/10.1109/CEC.2013.6557557
