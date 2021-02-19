# 3D-2D_Kibble_Lazarides_Shafi_Domain_Wall_in_Polar_Distorted_BPhase

These codes were used to calculate the equlibrium state configurations of 2D(3D in translation symmetry case) Kiibble-Lazarides-Shafi Sring Wall and the Spin vector soliton connecting on this composite objects.

The numerical results coming from these codes were publised in two mf my papers:

Phys. Rev. Research 2, 043356, link : https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.2.043356

Nat. Commun. 10, 237 (2019). https://doi.org/10.1038/s41467-018-08204-8

The algorithm is non-linear BFGS optimazation, and was implemented by myself in the matlab envriment. These codes run with parallrization calculation acceleration. The shell .sh script was used to control the runing state on the super computer.

C++ language was used to implement the labriry, which calculate the "A-Matrix". It is one of the most important entity when the configuraion and free energy are discretized with finite element strategy. To compile the library and the interface, run the .m generation script.

( I can not upload the .sh script now, may be githut error.)
