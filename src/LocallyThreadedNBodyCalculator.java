import data.*;
import data.Quadrant.RelativePosition;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

import static data.Quadrant.RelativePosition.*;

public class LocallyThreadedNBodyCalculator {

    private final int threadCount;
    private final CyclicBarrier barrier;
    private final ExecutorService executor;
    private final AtomicInteger activeWorkers;
    private final List<Point> points;
    private Quadrant rootQuadrant;

    private final int partitionSize;

    private final Random random = new Random();

    public LocallyThreadedNBodyCalculator() {
        this.threadCount = Runtime.getRuntime().availableProcessors();
        this.barrier = new CyclicBarrier(threadCount);
        this.executor = Executors.newFixedThreadPool(threadCount);
        this.activeWorkers = new AtomicInteger(threadCount);
        this.points = new ArrayList<>();

        partitionSize = CommonCore.bodyCount / threadCount;

        System.out.println("Initialization:");
        System.out.println("Bodies Count = " + CommonCore.bodyCount);
        System.out.println("Iteration Count = " + CommonCore.iterationCount);
        System.out.println("Thread Count = " + threadCount);
        System.out.println("Particles per thread = " + partitionSize);

        initializeArrays();
        initializeSpace();
    }


    private void initializeArrays() {
        for (int i = 0; i < CommonCore.bodyCount; i++) {
            points.add(new Point(
                    new Position(),
                    0,
                    new Mass(),
                    new Velocity()
            ));
        }
    }

    private class Worker implements Runnable {
        private final int threadId;
        private final int startIndex;
        private final int endIndex;

        Worker(int threadId) {
            this.threadId = threadId;
            this.startIndex = threadId * partitionSize;
            this.endIndex = startIndex + partitionSize;
        }

        @Override
        public void run() {
            try {
                for (int iter = 0; iter < CommonCore.iterationCount; iter++) {
                    for (var i = startIndex; i < endIndex; i++) {
                        computeForceForParticle(i);
                    }
                    barrier.await();

                    for (var i = startIndex; i < endIndex; i++) {
                        updateVelocityAndPosition(i);
                    }
                    barrier.await();

                    if (threadId == 0) {
                        writePositions(iter);
                        if ((iter + 1) % 50 == 0) {
                            System.out.printf("Completed iteration %d of %d%n", iter + 1, CommonCore.iterationCount);
                        }
                    }
                    barrier.await();
                }
            } catch (Exception e) {
                System.err.println("Error in worker " + threadId + ": " + e.getMessage());
            } finally {
                activeWorkers.decrementAndGet();
            }
        }
    }

    private double generateRand() {
        return random.nextDouble();
    }

    private double generateRandEx() {
        return 2 * generateRand() - 1;
    }

    private void initializeSpace() {
        double innerXBound = CommonCore.maxWidth - CommonCore.maxBodyRadius;
        double innerYBound = CommonCore.maxHeight - CommonCore.maxBodyRadius;

        for (int i = 0; i < CommonCore.bodyCount; i++) {
            points.set(i, new Point(
                    new Position(
                            generateRand() * innerXBound,
                            generateRand() * innerYBound),
                    CommonCore.maxBodyRadius * generateRand(),
                    new Mass(CommonCore.defaultMass * generateRand()),
                    new Velocity(generateRandEx(), generateRandEx())
            ));
        }
    }

    private double computeDistance(Position a, Position b) {
        return Math.sqrt(Math.pow(a.horizontal() - b.horizontal(), 2.0) +
                Math.pow(a.vertical() - b.vertical(), 2.0));
    }

    private void computeForceForParticle(int index) {
        double forceX = 0.0, forceY = 0.0;

        for (int jindex = 0; jindex < CommonCore.bodyCount; jindex++) {
            if (jindex == index) {
                continue;
            }

            final var distance = computeDistance(points.get(index).position(), points.get(jindex).position());
            if (distance == 0) {
                continue;
            }
            final var multiplier = (
                    CommonCore.GravitationalConstant * (points.get(index).mass().value() * points.get(jindex).mass().value())
                            / (Math.pow(distance, 2.0))) / distance;

            forceX += multiplier * ((points.get(jindex).position().horizontal() - points.get(index).position().horizontal()));
            forceY += multiplier * ((points.get(jindex).position().vertical() - points.get(index).position().vertical()));

        }
        points.get(index).setForce(new Force(forceX, forceY));
    }

    private void updateVelocityAndPosition(int index) {
        final var currentPosition = points.get(index).position();
        var currentVelocity = points.get(index).velocity();
        final var currentForce = points.get(index).force();
        final var currentMass = points.get(index).mass();
        currentVelocity = new Velocity(
                currentVelocity.horizontal() + (currentForce.horizontal() / currentMass.value()) * CommonCore.timeIncrement,
                currentVelocity.vertical() + (currentForce.vertical() / currentMass.value()) * CommonCore.timeIncrement);
        points.get(index).setVelocity(currentVelocity);
        points.get(index).setPosition(new Position(
                currentPosition.horizontal() + currentVelocity.horizontal() * CommonCore.timeIncrement,
                currentPosition.vertical() + currentVelocity.vertical() * CommonCore.timeIncrement));
        handleBoundaryConditions(index);
    }

    private void handleBoundaryConditions(int i) {
        var current = points.get(i);
        var pos = current.position();
        var x = pos.horizontal();
        var y = pos.vertical();
        var radius = current.radius();
        if ((x + radius) >= CommonCore.maxWidth || ((x - radius) <= 0)) {
            x *= -1;
        }
        if ((y + radius >= CommonCore.maxHeight) || ((y - radius) <= 0)) {
            y *= -1;
        }
        current.setPosition(new Position(x, y));
    }

    private void writePositions(int iteration) {
        String filename = String.format("thread_positions_iter_%d.csv", iteration);
        try (final var writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)))) {
            writer.println("particle_id,x,y,z");
            for (var index = 0; index < CommonCore.bodyCount; index++) {
                writer.printf("%d,%.6f,%.6f%n",
                        index,
                        points.get(index).position().horizontal(),
                        points.get(index).position().vertical());
            }
        } catch (IOException e) {
            System.err.println("Error writing to output file: " + e.getMessage());
        }
    }

    public void run() {
        writePositions(0);

        for (int i = 0; i < threadCount; i++) {
            executor.submit(new Worker(i));
        }
        executor.close();
    }

    private void generateQuadTree() {
        rootQuadrant = new Quadrant(
                -1,
                new Dimension(CommonCore.maxWidth,
                        CommonCore.maxHeight),
                new Position(0, 0));
        rootQuadrant.setIndex(0);

        for (var i = 1; i < CommonCore.bodyCount; i++) {
            var currentCell = rootQuadrant;
            while (currentCell.innerQuadrantCount() != 0) {
                final var subcell = locateQuadrant(currentCell, i);
                currentCell = currentCell.innerQuadrants().get(subcell.ordinal());
            }
            addToCell(currentCell, i);
        }
    }

    private void setLocationOfSubcells(Quadrant quadrant, double width, double height) {
        quadrant.innerQuadrants().getFirst().setBottomLeftCorner(quadrant.bottomLeftCorner());
        quadrant.innerQuadrants().get(1).setBottomLeftCorner(
                new Position(quadrant.bottomLeftCorner().horizontal() + width, quadrant.bottomLeftCorner().vertical()));
        quadrant.innerQuadrants().get(2).setBottomLeftCorner(
                new Position(quadrant.bottomLeftCorner().horizontal(), quadrant.bottomLeftCorner().vertical() + height));
        quadrant.innerQuadrants().get(3).setBottomLeftCorner(
                new Position(quadrant.bottomLeftCorner().horizontal() + width, quadrant.bottomLeftCorner().vertical() + height));
    }

    private void generateSubcells(Quadrant quadrant) {
        double width = quadrant.dimensions().horizontal() / 2.0;
        double height = quadrant.dimensions().vertical() / 2.0;

        quadrant.setInnerQuadrantCount(4);

        for (int i = 0; i < quadrant.innerQuadrantCount(); i++) {
            quadrant.innerQuadrants().set(i, new Quadrant(
                    -1,
                    new Dimension(width, height),
                    new Position(0, 0)));
        }

        setLocationOfSubcells(quadrant, width, height);
    }

    private RelativePosition locateQuadrant(Quadrant quadrant, int index) {
        final var bottomLeft = quadrant.innerQuadrants().getFirst().bottomLeftCorner();
        final var bottomLeftWidth = bottomLeft.horizontal();
        final var bottomLeftHeight = bottomLeft.vertical();
        final var width = points.get(index).position().horizontal();
        final var height = points.get(index).position().vertical();
        if (width == bottomLeftWidth && height == bottomLeftHeight) {
            return BOTTOM_LEFT;
        }
        if (width == bottomLeftWidth && height != bottomLeftHeight) {
            return TOP_LEFT;
        }
        if (height == bottomLeftHeight) {
            return TOP_RIGHT;
        }
        return BOTTOM_LEFT;
    }

    private void addToCell(Quadrant quadrant, int index) {
        if (quadrant.index() == -1) {
            quadrant.setIndex(index);
            return;
        }

        generateSubcells(quadrant);

        final var firstSubcellId = locateQuadrant(quadrant, quadrant.index());
        quadrant.innerQuadrants().get(firstSubcellId.ordinal()).setIndex(quadrant.index());

        final var secondSubcellId = locateQuadrant(quadrant, index);

        if (firstSubcellId == secondSubcellId) {
            addToCell(quadrant.innerQuadrants().get(firstSubcellId.ordinal()), index);
        } else {
            quadrant.innerQuadrants().get(secondSubcellId.ordinal()).setIndex(index);
        }
    }

    private void computeForceFromCell(Quadrant quadrant, int index) {
        final var distance = computeDistance(points.get(index).position(), points.get(quadrant.index()).position());
        final var multiplier = (CommonCore.GravitationalConstant * (points.get(index).mass().value() * points.get(quadrant.index()).mass().value()) / (Math.pow(distance, 2.0)));
        final var force = points.get(index).force();
        var fx = force.horizontal();
        var fy = force.vertical();
//        fx += ;
//        fy += ;
//        forces.get(index).horizontal() += multiplier * ((positions.get(cell.index).horizontal() - positions.get(index).horizontal()) / distance);
//        forces.get(index).vertical() += multiplier * ((positions.get(cell.index).vertical() - positions.get(index).vertical()) / distance);
    }

    private void computeForceFromQuadtree(Quadrant quadrant, int index) {
        if (quadrant.innerQuadrantCount() == 0) {
            if (quadrant.index() != -1 && quadrant.index() != index) {
                computeForceFromCell(quadrant, index);
            }
        } else {
            double d = computeDistance(points.get(index).position(), points.get(quadrant.index()).position());

            if (CommonCore.angle > (quadrant.bottomLeftCorner().horizontal() / d)) {
                computeForceFromCell(quadrant, index);
            } else {
                for (int i = 0; i < quadrant.innerQuadrantCount(); i++) {
                    computeForceFromQuadtree(quadrant.innerQuadrants().get(i), index);
                }
            }
        }
    }

    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
        simulation.run();
    }
}