Install NEURON 7.4 with Python enabled
======================================

1. Install build dependencies
   ```
   sudo apt install g++5 gcc-5 build-essential libx11-dev libxext-dev libpython2.7-dev ncurses-dev python-scipy python-matplotlib
   ```

2. Decide where to put NEURON, for example `$HOME/nrn/7.4/`, which we
   will call `$NRNPREFIX`. We can set this in the terminal:
   ```
   NRNPREFIX=$HOME/nrn/7.4/
   ```

3. Download
   [iv-1.9.tar.gz](https://neuron.yale.edu/ftp/neuron/versions/v7.4/iv-1.9.tar.gz)
   from https://neuron.yale.edu/ftp/neuron/versions/v7.4/

4. Unpack, configure and compile
   ```
   tar zxvf iv-1.9.tar.gz
   cd iv
   CC=gcc-5 CXX=g++5 ./configure --prefix=$NRNPREFIX/iv
   make
   make install
   cd ..
   ```

5. Download
   [nrn-7.4.tar.gz](https://neuron.yale.edu/ftp/neuron/versions/v7.4/nrn-7.4.tar.gz)
   from https://neuron.yale.edu/ftp/neuron/versions/v7.4/

6. Unpack, configure and compile
   ```
   tar zxvf nrn-7.4.tar.gz
   cd nrn
   CC=gcc-5 CXX=g++5 ./configure --prefix=$NRNPREFIX/iv --with-iv=$NRNPREFIX/iv --with-nrnpython=pyton2.7
   make
   make install
   cd ..
   ```

7. Set the `PYTHONPATH`:
   ```
   export PYTHONPATH=$HOME/nrn/7.4/nrn/lib/python/:$PYTHONPATH
   ```
   This should be added to your shell configuration, e.g. `~/.profile`
   or `~/.bash_profile`.
