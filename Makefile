all:
	javac -cp SpatialKappa-v2.1.1.jar SpatialKappaSim.java 
	java -cp './*' SpatialKappaSim main

pid = $(shell /usr/sbin/ss -lp  '( dport = :25333 or sport = :25333 )' | tail -n +2 |  perl -p -e 's/.*java\",([0-9]+),.*/\1/;')

CLASSFILES = SpatialKappaSim.class SpatialKappaSimEntryPoint.class
JAVAFILES := $(CLASSFILES:.class=.java)

$(CLASSFILES): $(JAVAFILES)
	javac -cp /home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar $(JAVAFILES)

.PHONY test: $(CLASSFILES)
	make -C test test

py4j: 
	javac -cp /home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar SpatialKappaSim.java SpatialKappaSimEntryPoint.java
  ## Make sure no gateway servers are running
	echo $(pid)
	if [ "x$(pid)" != "x" ]; then kill -9 $(pid) ; fi
	CLASSPATH=.:/home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../SpatialKappa/*':../py4j/py4j0.7.jar	java SpatialKappaSimEntryPoint main &
	sleep 1
	@ CLASSPATH=.:/home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../SpatialKappa/*':../py4j/py4j0.7.jar python2.7 -i hh_autapse_ca_rxd.py

test/x86_64/special: test/capulse.mod test/caPump.mod
	cd test; /disk/scratch/sterratt/x86_64/bin/nrnivmodl

test-bla: 
	javac -cp /home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar SpatialKappaSim.java SpatialKappaSimEntryPoint.java
  ## Make sure no gateway servers are running
	echo $(pid)
	if [ "x$(pid)" != "x" ]; then kill -9 $(pid) ; fi
	CLASSPATH=.:/home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../SpatialKappa/*':../py4j/py4j0.7.jar	java SpatialKappaSimEntryPoint main &
	sleep 1
	cd test && make test

