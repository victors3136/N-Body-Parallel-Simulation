import data.*;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

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
                    new Velocity(),
                    new Force()
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
                    new Velocity(generateRandEx(), generateRandEx()),
                    new Force()
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
        velocities.set(index, new Velocity(
                velocities.get(index).horizontal() + (forces.get(index).horizontal() / masses.get(index).value()) * CommonCore.timeIncrement,
                velocities.get(index).vertical() + (forces.get(index).vertical() / masses.get(index).value()) * CommonCore.timeIncrement,
                velocities.get(index).depth() + (forces.get(index).depth() / masses.get(index).value()) * CommonCore.timeIncrement
        ));
        positions.set(index, new Position(
                positions.get(index).horizontal() + velocities.get(index).horizontal() * CommonCore.timeIncrement,
                positions.get(index).vertical() + velocities.get(index).vertical() * CommonCore.timeIncrement,
                positions.get(index).depth() + velocities.get(index).depth() * CommonCore.timeIncrement));
        handleBoundaryConditions(index);
    }

    private void handleBoundaryConditions(int i) {
        var x = positions.get(i).horizontal();
        var y = positions.get(i).vertical();
        var z = positions.get(i).depth();
        if ((x + radiuses.get(i)) >= CommonCore.maxWidth || ((x - radiuses.get(i)) <= 0)) {
            x *= -1;
        }
        if ((y + radiuses.get(i) >= CommonCore.maxHeight) || ((y - radiuses.get(i)) <= 0)) {
            y *= -1;
        }
        if ((z >= CommonCore.maxDepth) || ((z - radiuses.get(i)) <= 0)) {
            z *= -1;
        }
        positions.set(i, new Position(x, y, z));
    }

    private void writePositions(int iteration) {
        String filename = String.format("thread_positions_iter_%d.csv", iteration);
        try (final var writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)))) {
            writer.println("particle_id,x,y,z");
            for (var index = 0; index < CommonCore.bodyCount; index++) {
                writer.printf("%d,%.6f,%.6f,%.6f%n",
                        index,
                        positions.get(index).horizontal(),
                        positions.get(index).vertical(),
                        positions.get(index).depth());
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

    private void generateOcttree() {
        rootQuadrant = new Quadrant(
                new Dimension(
                        CommonCore.maxWidth,
                        CommonCore.maxHeight,
                        CommonCore.maxDepth),
                new Position(0, 0, 0));
        rootQuadrant.index = 0;

        for (var i = 1; i < CommonCore.bodyCount; i++) {
            var currentCell = rootQuadrant;
            while (currentCell.innerQuadrantCount != 0) {
                final var subcell = locateSubcell(currentCell, i);
                currentCell = currentCell.innerQuadrants.get(subcell);
            }
            addToCell(currentCell, i);
        }
    }

    private void setLocationOfSubcells(Quadrant quadrant, double width, double height, double depth) {
        quadrant.innerQuadrants.getFirst().setBottomLeftCorner(quadrant.bottomLeftCorner);
        quadrant.innerQuadrants.get(1).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal() + width, quadrant.bottomLeftCorner.vertical(), quadrant.bottomLeftCorner.depth()));
        quadrant.innerQuadrants.get(2).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal() + width, quadrant.bottomLeftCorner.vertical(), quadrant.bottomLeftCorner.depth() + depth));
        quadrant.innerQuadrants.get(3).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal(), quadrant.bottomLeftCorner.vertical(), quadrant.bottomLeftCorner.depth() + depth));
        quadrant.innerQuadrants.get(4).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal(), quadrant.bottomLeftCorner.vertical() + height, quadrant.bottomLeftCorner.depth()));
        quadrant.innerQuadrants.get(5).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal() + height, quadrant.bottomLeftCorner.vertical() + height, quadrant.bottomLeftCorner.depth()));
        quadrant.innerQuadrants.get(6).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal() + height, quadrant.bottomLeftCorner.vertical() + height, quadrant.bottomLeftCorner.depth() + depth));
        quadrant.innerQuadrants.get(7).setBottomLeftCorner(new Position(quadrant.bottomLeftCorner.horizontal(), quadrant.bottomLeftCorner.vertical() + height, quadrant.bottomLeftCorner.depth() + depth));
    }

    private void generateSubcells(Quadrant quadrant) {
        double width = quadrant.dimensions.horizontal() / 2.0;
        double height = quadrant.dimensions.vertical() / 2.0;
        double depth = quadrant.dimensions.depth() / 2.0;

        quadrant.innerQuadrantCount = 8;

        for (int i = 0; i < quadrant.innerQuadrantCount; i++) {
            quadrant.innerQuadrants.set(i, new Quadrant(new Dimension(width, height, depth), new Position(0, 0, 0)));
        }

        setLocationOfSubcells(quadrant, width, height, depth);
    }

    private int locateSubcell(Quadrant quadrant, int index) {
        final var sixth = quadrant.innerQuadrants.get(6).bottomLeftCorner;
        final var sixthX = sixth.horizontal();
        final var sixthY = sixth.vertical();
        final var sixthZ = sixth.depth();
        final var x = positions.get(index).horizontal();
        final var y = positions.get(index).vertical();
        final var z = positions.get(index).depth();
        if (x > sixthX) {
            if (y > sixthY) {
                return z > sixthZ ? 6 : 5;
            } else {
                return z > sixthZ ? 2 : 1;
            }
        } else {
            if (y > sixthY) {
                return z > sixthZ ? 7 : 4;
            } else {
                return z > sixthZ ? 3 : 0;
            }
        }
    }

    private void addToCell(Quadrant quadrant, int index) {
        if (quadrant.index == -1) {
            quadrant.index = index;
            return;
        }

        generateSubcells(quadrant);

        int firstSubcellId = locateSubcell(quadrant, quadrant.index);
        quadrant.innerQuadrants.get(firstSubcellId).index = quadrant.index;

        int secondSubcellId = locateSubcell(quadrant, index);

        if (firstSubcellId == secondSubcellId) {
            addToCell(quadrant.innerQuadrants.get(firstSubcellId), index);
        } else {
            quadrant.innerQuadrants.get(secondSubcellId).index = index;
        }
    }

    private void computeForceFromCell(Quadrant quadrant, int index) {
        double d = computeDistance(positions.get(index), positions.get(quadrant.index));
        double f = (CommonCore.GravitationalConstant * (masses.get(index).value() * masses.get(quadrant.index).value()) / (Math.pow(d, 2.0)));

//        forces.get(index).horizontal() += f * ((positions.get(cell.index).horizontal() - positions.get(index).horizontal()) / d);
//        forces.get(index).vertical() += f * ((positions.get(cell.index).vertical() - positions.get(index).vertical()) / d);
//        forces.get(index).depth() += f * ((positions.get(cell.index).depth() - positions.get(index).depth()) / d);
    }

    private void computeForceFromOcttree(Quadrant quadrant, int index) {
        if (quadrant.innerQuadrantCount == 0) {
            if (quadrant.index != -1 && quadrant.index != index) {
                computeForceFromCell(quadrant, index);
            }
        } else {
            double d = computeDistance(positions.get(index), positions.get(quadrant.index));

            if (CommonCore.angle > (quadrant.bottomLeftCorner.horizontal() / d)) {
                computeForceFromCell(quadrant, index);
            } else {
                for (int i = 0; i < quadrant.innerQuadrantCount; i++) {
                    computeForceFromOcttree(quadrant.innerQuadrants.get(i), index);
                }
            }
        }
    }

    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
        simulation.run();
    }
}