This repository contains the Matlab implementations of the Poisson multi-Bernoulli mixture (PMBM) tracker proposed in 

Granstršm, K., Svensson, L., Xia, Y., Williams, J., & Garc’a-Fem‡ndez, ç. F. (2018, July). Poisson multi-Bernoulli mixture trackers: continuity through random finite sets of trajectories. In 2018 21st International Conference on Information Fusion (FUSION) IEEE.

Full texts can be found in https://arxiv.org/abs/1812.05131


The filters are evaluated using the generalised optimal subpattern-assignment (GOSPA)

[C] A. S. Rahmathullah, A. F. García-Fernández, and L. Svensson, “Generalized optimal sub-pattern assignment metric,” in 20th International
Conference on Information Fusion, 2017.

[D] Video on GOSPA: https://www.youtube.com/watch?v=M79GTTytvCM


- Simulation.m runs the PMBM tracker

- density class handle: 1. @trajectoryGaussianFilteringForm: only Kalman filtering is performed; 2. @trajectoryAccumulatedStateDensity: performs smoothing-while-filtering.

- data association algorithm: 1. kBest2DAssign: Murty's algorithm, copied from https://github.com/USNavalResearchLaboratory/TrackerComponentLibrary; 2. assign2DByGibbs: Gibb's sampling, adapted from http://ba-tuong.vo-au.com/rfs_tracking_toolbox_new.7z, gibbswrap_jointpredupdt_custom.m

