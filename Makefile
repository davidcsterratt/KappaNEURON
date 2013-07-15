all:
	javac -cp SpatialKappa-v2.1.1.jar SpatialKappaSim.java 
	java -cp './*' SpatialKappaSim main

pid = $(shell /usr/sbin/ss -lp  '( dport = :25333 or sport = :25333 )' | tail -n +2 |  perl -p -e 's/.*java\",([0-9]+),.*/\1/;')

py4j: 
	javac -cp SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar SpatialKappaSim.java SpatialKappaSimEntryPoint.java
  ## Make sure no gateway servers are running
	echo $(pid)
	if [ "x$(pid)" != "x" ]; then kill -9 $(pid) ; fi
	CLASSPATH='./*':../py4j/py4j0.7.jar	java SpatialKappaSimEntryPoint main &
	sleep 4
	@ CLASSPATH='./*':../py4j/py4j0.7.jar python spatialKappa.py

spatialKappaNeuron: 
	javac -cp SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar SpatialKappaSim.java SpatialKappaSimEntryPoint.java
  ## Make sure no gateway servers are running
	echo $(pid)
	if [ "x$(pid)" != "x" ]; then kill -9 $(pid) ; fi
	CLASSPATH='./*':../py4j/py4j0.7.jar	java SpatialKappaSimEntryPoint main &
	sleep 4
	@ CLASSPATH='./*':../py4j/py4j0.7.jar python spatialKappaNeuron.py

