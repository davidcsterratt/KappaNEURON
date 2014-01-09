CLASSFILES = SpatialKappaSim.class SpatialKappaSimEntryPoint.class
JAVAFILES := $(CLASSFILES:.class=.java)

$(CLASSFILES): $(JAVAFILES)
	javac -cp ../mlm/SpatialKappa/SpatialKappa-v2.1.1.jar:../mlm/py4j/py4j0.7.jar $(JAVAFILES)

all: $(CLASSFILES)

.PHONY test: all
	make -C test test

