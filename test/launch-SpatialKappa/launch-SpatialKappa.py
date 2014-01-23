from SpatialKappa import *

sk = SpatialKappa(classpath=os.pathsep.join(('../..', '../../../mlm/SpatialKappa/SpatialKappa-v2.1.1.jar', '../../../mlm/SpatialKappa/ant-antlr-3.2.jar')))
kappa_sim = sk.kappa_sim("ms", True)
kappa_sim.loadFile("caPump.ka")
kappa_sim.runToTime(0.01)
print kappa_sim.getObservation("P")
print kappa_sim.getObservation("ca")
