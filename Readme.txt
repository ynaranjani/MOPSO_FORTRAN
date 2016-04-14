This code is the Fortran 90 implementation of the Multi-Objective Particle Swarm
 Optimization (MOPSO). To learn more about it you can refer to:
	C. A. C. Coello, G. T. Pulido and M. S. Lechuga, "Handling multiple 
		objectives with particle swarm optimization" in IEEE Transactions on 
		Evolutionary Computation, vol. 8, no. 3, pp. 256-279, June 2004.
	Margarita Reyes-sierra and Carlos A. Coello Coello, "Multi-Objective 
		particle swarm optimizers: A survey of the state-of-the-art" in 
		INTERNATIONAL JOURNAL OF COMPUTATIONAL INTELLIGENCE RESEARCH,
		vol. 2, no. 3, pp. 287-308, 2006
To run the code you simply need to do 3 things:
1) modify or create an INCLUDE file in the folder called /in. look at one of the
	samples and make one for your own problem. 
2) make sure the solver you defined in the include file (line~83) exists in the
	ProblemBank.f90. The objective functions and constraints must be defined there
3) change the INCLUDE file in the MOPSO.f90. to the one you created/modified.

 The sample problems are taken for the following book:
	"Evolutionary Algorithms for Solving Multi-Objective Problems" by 
	Carlos Coello Coello, Gary B. Lamont and David A. van Veldhuizen
	take a look at sample problems on page Chapter 4, page 183-...

to run simply:
$ make
$ ./main

The results will be written in the out/ProbName folder (where ProbName is the 
name you chose for your problem at the problem INCLUDE file line: 45). 'PS.dat' and 'PF.dat'
 contain the Pareto set (points on design space) and Pareto front 
 (corresponding objective values). The other output might -depending on 
 your problem setup- be a 'log.dat' file which inlcudes some report of the
 events in the simulation time including timing of different parts of the
 code (helpful for debugging) and intermediate solutions 'progress.dat'
 that includes the archive contents as the generations evolve.

Please feel free to contact me if you have any questions.
	Yousef Naranjani
	yo.na1118@gmail.com
	April 2016
