import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

public class LocallyThreadedNBodyCalculator {

    // Thread-specific variables
    private final int threadCount;
    private final CyclicBarrier barrier;
    private final ExecutorService executor;
    private final AtomicInteger activeWorkers;

    private List<Position> positions;
    private List<Velocity> velocities;
    private List<Mass> masses;
    private List<Double> radiuses;
    private List<Force> forces;
    private Cell rootCell;

    private final int partitionSize;

    private final Random random = new Random();

    public LocallyThreadedNBodyCalculator() {
        this.threadCount = Runtime.getRuntime().availableProcessors();
        this.barrier = new CyclicBarrier(threadCount);
        this.executor = Executors.newFixedThreadPool(threadCount);
        this.activeWorkers = new AtomicInteger(threadCount);

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
        masses = List.of(new Mass[CommonCore.bodyCount]);
        radiuses = List.of(new Double[CommonCore.bodyCount]);
        positions = List.of(new Position[CommonCore.bodyCount]);
        velocities = List.of(new Velocity[CommonCore.bodyCount]);
        forces = List.of(new Force[CommonCore.bodyCount]);

        for (int i = 0; i < CommonCore.bodyCount; i++) {
            positions.set(i, new Position(0, 0, 0));
            velocities.set(i, new Velocity(0, 0, 0));
            forces.set(i, new Force(0, 0, 0));
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
        double ixbound = CommonCore.maxWidth - CommonCore.maxBodyRadius;
        double iybound = CommonCore.maxHeight - CommonCore.maxBodyRadius;
        double izbound = CommonCore.maxDepth - CommonCore.maxBodyRadius;

        for (int i = 0; i < CommonCore.bodyCount; i++) {
            masses.set(i, new Mass(CommonCore.defaultMass * generateRand()));
            radiuses.set(i, CommonCore.maxBodyRadius * generateRand());
            positions.set(i, new Position(
                    generateRand() * ixbound,
                    generateRand() * iybound,
                    generateRand() * izbound));
            velocities.set(i, new Velocity(generateRandEx(), generateRandEx(), generateRandEx()));
        }
    }

    private double computeDistance(Position a, Position b) {
        return Math.sqrt(Math.pow(a.horizontal() - b.horizontal(), 2.0) +
                Math.pow(a.vertical() - b.vertical(), 2.0) +
                Math.pow(a.depth() - b.depth(), 2.0));
    }

    private void computeForceForParticle(int index) {
        double forceX = 0.0, forceY = 0.0, forceZ = 0.0;

        for (int jindex = 0; jindex < CommonCore.bodyCount; jindex++) {
            if (jindex == index) {
                continue;
            }

            final var distance = computeDistance(positions.get(index), positions.get(jindex));
            if (distance == 0) {
                continue;
            }
            final var multiplier = (
                    CommonCore.GravitationalConstant * (masses.get(index).value() * masses.get(jindex).value())
                            / (Math.pow(distance, 2.0))) / distance;

            forceX += multiplier * ((positions.get(jindex).horizontal() - positions.get(index).horizontal()));
            forceY += multiplier * ((positions.get(jindex).vertical() - positions.get(index).vertical()));
            forceZ += multiplier * ((positions.get(jindex).depth() - positions.get(index).depth()));
        }
        forces.set(index, new Force(forceX, forceY, forceZ));
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
        rootCell = new Cell(
                new Dimension(
                        CommonCore.maxWidth,
                        CommonCore.maxHeight,
                        CommonCore.maxDepth),
                new Position(0, 0, 0));
        rootCell.index = 0;

        for (var i = 1; i < CommonCore.bodyCount; i++) {
            var currentCell = rootCell;
            while (currentCell.subcellCount != 0) {
                final var subcell = locateSubcell(currentCell, i);
                currentCell = currentCell.subcells.get(subcell);
            }
            addToCell(currentCell, i);
        }
    }

    private void setLocationOfSubcells(Cell cell, double width, double height, double depth) {
        cell.subcells.getFirst().setPosition(cell.position);
        cell.subcells.get(1).setPosition(new Position(cell.position.horizontal() + width, cell.position.vertical(), cell.position.depth()));
        cell.subcells.get(2).setPosition(new Position(cell.position.horizontal() + width, cell.position.vertical(), cell.position.depth() + depth));
        cell.subcells.get(3).setPosition(new Position(cell.position.horizontal(), cell.position.vertical(), cell.position.depth() + depth));
        cell.subcells.get(4).setPosition(new Position(cell.position.horizontal(), cell.position.vertical() + height, cell.position.depth()));
        cell.subcells.get(5).setPosition(new Position(cell.position.horizontal() + height, cell.position.vertical() + height, cell.position.depth()));
        cell.subcells.get(6).setPosition(new Position(cell.position.horizontal() + height, cell.position.vertical() + height, cell.position.depth() + depth));
        cell.subcells.get(7).setPosition(new Position(cell.position.horizontal(), cell.position.vertical() + height, cell.position.depth() + depth));
    }

    private void generateSubcells(Cell cell) {
        double width = cell.size.horizontal() / 2.0;
        double height = cell.size.vertical() / 2.0;
        double depth = cell.size.depth() / 2.0;

        cell.subcellCount = 8;

        for (int i = 0; i < cell.subcellCount; i++) {
            cell.subcells.set(i, new Cell(new Dimension(width, height, depth), new Position(0, 0, 0)));
        }

        setLocationOfSubcells(cell, width, height, depth);
    }

    private int locateSubcell(Cell cell, int index) {
        final var sixth = cell.subcells.get(6).position;
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

    private void addToCell(Cell cell, int index) {
        if (cell.index == -1) {
            cell.index = index;
            return;
        }

        generateSubcells(cell);

        int firstSubcellId = locateSubcell(cell, cell.index);
        cell.subcells.get(firstSubcellId).index = cell.index;

        int secondSubcellId = locateSubcell(cell, index);

        if (firstSubcellId == secondSubcellId) {
            addToCell(cell.subcells.get(firstSubcellId), index);
        } else {
            cell.subcells.get(secondSubcellId).index = index;
        }
    }

    private void computeForceFromCell(Cell cell, int index) {
        double d = computeDistance(positions.get(index), positions.get(cell.index));
        double f = (CommonCore.GravitationalConstant * (masses.get(index).value() * masses.get(cell.index).value()) / (Math.pow(d, 2.0)));

//        forces.get(index).horizontal() += f * ((positions.get(cell.index).horizontal() - positions.get(index).horizontal()) / d);
//        forces.get(index).vertical() += f * ((positions.get(cell.index).vertical() - positions.get(index).vertical()) / d);
//        forces.get(index).depth() += f * ((positions.get(cell.index).depth() - positions.get(index).depth()) / d);
    }

    private void computeForceFromOcttree(Cell cell, int index) {
        if (cell.subcellCount == 0) {
            if (cell.index != -1 && cell.index != index) {
                computeForceFromCell(cell, index);
            }
        } else {
            double d = computeDistance(positions.get(index), positions.get(cell.index));

            if (CommonCore.angle > (cell.position.horizontal() / d)) {
                computeForceFromCell(cell, index);
            } else {
                for (int i = 0; i < cell.subcellCount; i++) {
                    computeForceFromOcttree(cell.subcells.get(i), index);
                }
            }
        }
    }

    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
        simulation.run();
    }
}