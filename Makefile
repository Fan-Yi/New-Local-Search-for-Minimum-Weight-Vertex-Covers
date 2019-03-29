all: myMinWVC

myMinWVC: myMinWVC.cpp graph.h localSearch.h constants.h \
					myBijection.h weightBuckets.h hugeInt.h operandSets.h coverHash.h config.h
	g++ -std=gnu++0x -g -O3 -static myMinWVC.cpp -o myMinWVC

clean: rm -f *~
