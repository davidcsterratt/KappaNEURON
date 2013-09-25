CLASSFILES = SpatialKappaSim.class SpatialKappaSimEntryPoint.class
JAVAFILES := $(CLASSFILES:.class=.java)

$(CLASSFILES): $(JAVAFILES)
	javac -cp ../SpatialKappa/SpatialKappa-v2.1.1.jar:../py4j/py4j0.7.jar $(JAVAFILES)

all: $(CLASSFILES)

.PHONY test: all
	make -C test test

# py4j: gateway
#	@ CLASSPATH=.:../SpatialKappa/SpatialKappa-v2.1.1.jar:'../SpatialKappa/*':../py4j/py4j0.7.jar python2.7 -i hh_autapse_ca_rxd.py

