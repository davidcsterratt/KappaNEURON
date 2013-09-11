CLASSFILES = SpatialKappaSim.class SpatialKappaSimEntryPoint.class
JAVAFILES := $(CLASSFILES:.class=.java)

$(CLASSFILES): $(JAVAFILES)
	javac -cp /home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar $(JAVAFILES)

pid = $(shell /usr/sbin/ss -lp  '( dport = :25333 or sport = :25333 )' | tail -n +2 |  perl -p -e 's/.*java\",([0-9]+),.*/\1/;')

gateway: $(CLASSFILES)
	echo $(pid)
	if [ "x$(pid)" != "x" ]; then kill -9 $(pid) ; fi
	CLASSPATH=.:/home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../SpatialKappa/*':../py4j/py4j0.7.jar	java SpatialKappaSimEntryPoint main &
	sleep 1

.PHONY test: gateway
	make -C test test

py4j: gateway
	@ CLASSPATH=.:/home/sterratt/projects/mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:'../SpatialKappa/*':../py4j/py4j0.7.jar python2.7 -i hh_autapse_ca_rxd.py

