import py4j.GatewayServer;

public class SpatialKappaSimEntryPoint {

    public SpatialKappaSimEntryPoint() {
    }

    public SpatialKappaSim newSpatialKappaSim(String timeUnits) {
        SpatialKappaSim sks = new SpatialKappaSim(timeUnits);
        return sks;
    }

    public SpatialKappaSim newSpatialKappaSim() {
        return this.newSpatialKappaSim("ms");
    }

    public static void main(String[] args) {
        GatewayServer gatewayServer = new GatewayServer(new SpatialKappaSimEntryPoint());
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }

}