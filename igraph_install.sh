git clone https://github.com/igraph/igraph.git
cd igraph
mkdir build
cd build
cmake ..
make -j$(nproc)
sudo make install