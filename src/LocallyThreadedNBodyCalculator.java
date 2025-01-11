import data.*;

import java.util.List;

public class LocallyThreadedNBodyCalculator {

    private final int threadCount;
    private final List<Point> points;
    private Quadrant rootQuadrant;

    private final int partitionSize;

    public LocallyThreadedNBodyCalculator() {
        this.threadCount = Runtime.getRuntime().availableProcessors();
        this.points = Generator.points(CommonCore.bodyCount);
        partitionSize = CommonCore.bodyCount / threadCount;

        System.out.println("Initialization:");
        System.out.println("Bodies Count = " + CommonCore.bodyCount);
        System.out.println("Iteration Count = " + CommonCore.iterationCount);
        System.out.println("Thread Count = " + threadCount);
        System.out.println("Particles per thread = " + partitionSize);

    }

    public void run() {
        rootQuadrant = generateQuadrants(CommonCore.maxWidth, CommonCore.maxHeight);

        rootQuadrant = new Quadrant();
        rootQuadrant.setDimensions(new Dimension(CommonCore.maxWidth, CommonCore.maxHeight));
        points.forEach(point -> {
           rootQuadrant.insert(point);
        });
    }

    private Quadrant generateQuadrants(double maxWidth, double maxHeight) {

        var quadrant = new Quadrant();
        quadrant.bottomLeftCorner = new Position(0, 0)
        quadrant.size = new Dimension(spaceDimension.dx, spaceDimension.dy)
        quadrant.subquadrants = []
        quadrant.subquadrantCount = 0
        quadrant.centerOfMass = new Position(0, 0)
        quadrant.overallMass = new Mass(0)
        return quadrant
    }

    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
//        simulation.run();
    }
}