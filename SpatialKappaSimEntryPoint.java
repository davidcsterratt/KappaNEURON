
import py4j.GatewayServer;

public class SpatialKappaSimEntryPoint {

    public SpatialKappaSimEntryPoint() {
    }

    public SpatialKappaSim newSpatialKappaSim() {
        SpatialKappaSim sks = new SpatialKappaSim();
        return sks;
    }

    public static void main(String[] args) {
        GatewayServer gatewayServer = new GatewayServer(new SpatialKappaSimEntryPoint());
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }

}