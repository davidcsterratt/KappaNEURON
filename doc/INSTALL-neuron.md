Install NEURON 7.4 with Python enabled
======================================

1. Install build dependencies
   ```
   sudo apt install g++5 gcc-5 libx11-dev libxext-dev libpython2.7-dev ncurses-dev python-scipy python-matplotlib
   ```

2. Decide where to put NEURON, for example `$HOME/nrn/7.4/`, which we
   will call `$NRNPREFIX`. We can set this in the terminal:
   ```
   NRNPREFIX=$HOME/nrn/7.4/
   ```

3. Download
   [iv-19.tar.gz](https://neuron.yale.edu/ftp/neuron/versions/v7.4/iv-19.tar.gz)
   from https://neuron.yale.edu/ftp/neuron/versions/v7.4/

4. Unpack, configure and compile
   ```
   tar zxvf iv-19.tar.gz
   cd iv-19
   CC=gcc-5 CXX=g++-5 ./configure --prefix=$NRNPREFIX/iv
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
   cd nrn-7.4
   CC=gcc-5 CXX=g++-5 ./configure --prefix=$NRNPREFIX/nrn --with-iv=$NRNPREFIX/iv --with-nrnpython=python2.7
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
   
   8. Test that it works by typing:
   ```
   python2.7 -m neuron
   ```
   The output should be:
   ```
   NEURON -- Release 7.4 (1370:16a7055d4a86) 2015-11-09
   Duke, Yale, and the BlueBrain Project -- Copyright 1984-2015
   See http://www.neuron.yale.edu/neuron/credits

   /usr/bin/python: No module named neuron.__main__; 'neuron' is a package and cannot be directly executed
   ```
