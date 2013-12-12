import py4j.GatewayServer;

public class SpatialKappaSimEntryPoint {

    public SpatialKappaSimEntryPoint() {
    }

    public SpatialKappaSim newSpatialKappaSim(String timeUnits, Boolean verbose) {
        SpatialKappaSim sks = new SpatialKappaSim(timeUnits, verbose);
        return sks;
    }

    public SpatialKappaSim newSpatialKappaSim() {
        return this.newSpatialKappaSim("ms", false);
    }

    public static void main(String[] args) {
        GatewayServer gatewayServer = new GatewayServer(new SpatialKappaSimEntryPoint());
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }

}