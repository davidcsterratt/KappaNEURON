
import py4j.GatewayServer;

public class SpatialKappaSimEntryPoint {

    private SpatialKappaSim sks;

    public SpatialKappaSimEntryPoint() {
        sks = new SpatialKappaSim();
    }

    public SpatialKappaSim getSpatialKappaSim() {
        return sks;
    }
    
    public void runByTime(int steps, int stepSize) {
        sks.runByTime(steps, stepSize);
    }

    public static void main(String[] args) {
        GatewayServer gatewayServer = new GatewayServer(new SpatialKappaSimEntryPoint());
        gatewayServer.start();
        System.out.println("Gateway Server Started");
    }

}