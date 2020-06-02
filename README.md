# khcore

This code can be compiled with g++. If the compilation is not successful, please use the later version of g++ and c++11 or higher standard libraries.

Compilation: 

	g++ "-std=c++11" -O3 -o run KHCore.cpp

The first parameter is the algorithm that the program runs. Please select one of {txt-to-bin, decompose, compare, scal, casestudy}

	txt-to-bin: format the txt file to the bin file.
	decompose: run the decomposition algorithms.
	compare: compute the precision of the sampling algrithm.
	scal: the scalability testing.
	casestudy: case studies testing on dblp dataset.

Running examples: 

	./run txt-to-bin data/example.txt
	./run decompose data/example.bin 2
	./run decompose data/example.bin 2 4 0.1 0.01
	./run decompose data/example.bin 2 4 0.1 0.01 data/ core-number.bin
	./run compare path file1 file2

If running 'scal'

	usage: ./run scal scal-rate vary infile h t sampling-rate error-rate outfile-path outfile
	scal-rate: the range is (0,100].
	vary: '1' for sampling vertices, or '2' for sampling edges.

If running 'casestudy'

	usage:./run compare file1 file2 file3 a-professor
	file1: the binary graph file.
	file2: the binary core number file.
	file3: the txt file with names of each id.
	a-professor: the name of a professor. Replace ' ' between the first name and the last name with '_'. 
