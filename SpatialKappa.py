import py4j.java_gateway as jg

class SpatialKappa:
    """Runs gateway to SpatialKappa"""

    port = None
    gateway = None
    gateway_client = None

    def __init__(self):
        self.port = jg.launch_gateway(classpath='.:../mlm/SpatialKappa/SpatialKappa-v2.1.1.jar', die_on_exit=True)
        self.gateway_client = jg.GatewayClient(port=self.port)
        self.gateway = jg.JavaGateway(self.gateway_client)        
        jg.java_import(self.gateway.jvm, 'SpatialKappaSim.*')
        print self.gateway

    def kappa_sim(self, time_units, verbose):
        return self.gateway.jvm.SpatialKappaSim(time_units, verbose)


