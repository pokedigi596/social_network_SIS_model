all:
	g++ -O3 *.cpp -I/usr/local/include/igraph -L/usr/local/lib -ligraph -Wno-deprecated-declarations -o igraph_test
clean:
	rm -rf *.exe