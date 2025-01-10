import java.io.*;
import java.util.Random;
import java.util.concurrent.*;
import java.util.concurrent.atomic.AtomicInteger;

public class LocallyThreadedNBodyCalculator {

    // Thread-specific variables
    private final int numThreads;
    private final CyclicBarrier barrier;
    private final ExecutorService executor;
    private final Thread[] workers;
    private final AtomicInteger activeWorkers;

    // Instance variables
    private Position[] position;   // Current positions for all particles
    private Velocity[] velocity;   // Velocity of particles
    private double[] mass;         // Mass of each particle
    private double[] radius;       // Radius of each particle
    private Force[] force;         // Force experienced by all particles
    private Cell rootCell;         // Root of BH octtree

    private int partSize;       // Number of particles each thread is responsible for

    // Inner classes for data structures
    static class Position {
        double px, py, pz;
    }

    static class Velocity {
        double vx, vy, vz;
    }

    static class Force {
        double fx, fy, fz;
    }

    static class Cell {
        int index;
        int noSubcells;
        double mass;
        double x, y, z;
        double cx, cy, cz;
        double width, height, depth;
        Cell[] subcells;

        Cell(double width, double height, double depth) {
            this.mass = 0;
            this.noSubcells = 0;
            this.index = -1;
            this.cx = 0;
            this.cy = 0;
            this.cz = 0;
            this.width = width;
            this.height = height;
            this.depth = depth;
            this.subcells = new Cell[8];
        }
    }

    private Random random = new Random();

    // Constructor
    public LocallyThreadedNBodyCalculator() {
        // Use number of available processors for threads
        this.numThreads = Runtime.getRuntime().availableProcessors();
        this.barrier = new CyclicBarrier(numThreads);
        this.executor = Executors.newFixedThreadPool(numThreads);
        this.workers = new Thread[numThreads];
        this.activeWorkers = new AtomicInteger(numThreads);


        // Calculate partition size
        partSize = CommonCore.N / numThreads;

        System.out.println("Initialization:");
        System.out.println("CommonCore.N = " + CommonCore.N);
        System.out.println("CommonCore.TIME = " + CommonCore.TIME);
        System.out.println("Threads = " + numThreads);
        System.out.println("Particles per thread = " + partSize);

        initializeArrays();
        initializeSpace();
    }

    private void initializeArrays() {
        mass = new double[CommonCore.N];
        radius = new double[CommonCore.N];
        position = new Position[CommonCore.N];
        velocity = new Velocity[CommonCore.N];
        force = new Force[CommonCore.N];

        // Initialize all objects
        for (int i = 0; i < CommonCore.N; i++) {
            position[i] = new Position();
            velocity[i] = new Velocity();
            force[i] = new Force();
        }
    }

    private class Worker implements Runnable {
        private final int threadId;
        private final int startIndex;
        private final int endIndex;

        Worker(int threadId) {
            this.threadId = threadId;
            this.startIndex = threadId * partSize;
            this.endIndex = startIndex + partSize;
        }

        @Override
        public void run() {
            try {
                for (int iter = 0; iter < CommonCore.TIME; iter++) {
                    // Compute forces for this partition
                    for (int i = startIndex; i < endIndex; i++) {
                        computeForceForParticle(i);
                    }
                    barrier.await();

                    // Update velocities and positions
                    for (int i = startIndex; i < endIndex; i++) {
                        updateVelocityAndPosition(i);
                    }
                    barrier.await();

                    // Let thread 0 handle output
                    if (threadId == 0) {
                        writePositions(iter);
                        if ((iter + 1) % 50 == 0) {
                            System.out.printf("Completed iteration %d of %d%n", iter + 1, CommonCore.TIME);
                        }
                    }
                    barrier.await();
                }
            } catch (Exception e) {
                System.err.println("Error in worker " + threadId + ": " + e.getMessage());
                e.printStackTrace();
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
        double ixbound = CommonCore.XBOUND - CommonCore.RBOUND;
        double iybound = CommonCore.YBOUND - CommonCore.RBOUND;
        double izbound = CommonCore.ZBOUND - CommonCore.RBOUND;

        for (int i = 0; i < CommonCore.N; i++) {
            mass[i] = CommonCore.MASS_OF_UNKNOWN * generateRand();
            radius[i] = CommonCore.RBOUND * generateRand();
            position[i].px = generateRand() * ixbound;
            position[i].py = generateRand() * iybound;
            position[i].pz = generateRand() * izbound;
            velocity[i].vx = generateRandEx();
            velocity[i].vy = generateRandEx();
            velocity[i].vz = generateRandEx();
        }
    }

    private double computeDistance(Position a, Position b) {
        return Math.sqrt(Math.pow(a.px - b.px, 2.0) +
                Math.pow(a.py - b.py, 2.0) +
                Math.pow(a.pz - b.pz, 2.0));
    }

    private void computeForceForParticle(int index) {
        force[index].fx = 0.0;
        force[index].fy = 0.0;
        force[index].fz = 0.0;

        for (int j = 0; j < CommonCore.N; j++) {
            if (j == index) continue;

            double d = computeDistance(position[index], position[j]);
            double f = (CommonCore.G * (mass[index] * mass[j]) / (Math.pow(d, 2.0)));

            force[index].fx += f * ((position[j].px - position[index].px) / d);
            force[index].fy += f * ((position[j].py - position[index].py) / d);
            force[index].fz += f * ((position[j].pz - position[index].pz) / d);
        }
    }

    private void updateVelocityAndPosition(int i) {
        // Update velocity
        velocity[i].vx += (force[i].fx / mass[i]) * CommonCore.DELTAT;
        velocity[i].vy += (force[i].fy / mass[i]) * CommonCore.DELTAT;
        velocity[i].vz += (force[i].fz / mass[i]) * CommonCore.DELTAT;

        // Update position
        position[i].px += velocity[i].vx * CommonCore.DELTAT;
        position[i].py += velocity[i].vy * CommonCore.DELTAT;
        position[i].pz += velocity[i].vz * CommonCore.DELTAT;

        // Handle boundary conditions
        handleBoundaryConditions(i);
    }

    private void handleBoundaryConditions(int i) {
        if ((position[i].px + radius[i]) >= CommonCore.XBOUND ||
                (position[i].px - radius[i]) <= 0) {
            velocity[i].vx *= -1;
        }
        if ((position[i].py + radius[i] >= CommonCore.YBOUND) ||
                (position[i].py - radius[i]) <= 0) {
            velocity[i].vy *= -1;
        }
        if ((position[i].pz + radius[i]) >= CommonCore.ZBOUND ||
                (position[i].pz - radius[i]) <= 0) {
            velocity[i].vz *= -1;
        }
    }

    private void writePositions(int iteration) {
        String filename = String.format("thread_positions_iter_%d.csv", iteration);
        try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)))) {
            writer.println("particle_id,x,y,z");
            for (int i = 0; i < CommonCore.N; i++) {
                writer.printf("%d,%.6f,%.6f,%.6f%n",
                        i,
                        position[i].px,
                        position[i].py,
                        position[i].pz);
            }
        } catch (IOException e) {
            System.err.println("Error writing to output file: " + e.getMessage());
        }
    }

    public void run() {
        // Write initial positions
        writePositions(0);

        // Create and start worker threads
        for (int i = 0; i < numThreads; i++) {
            workers[i] = new Thread(new Worker(i));
            workers[i].start();
        }

        // Wait for all workers to complete
        try {
            while (activeWorkers.get() > 0) {
                Thread.sleep(100);
            }
        } catch (InterruptedException e) {
            System.err.println("Main thread interrupted: " + e.getMessage());
        } finally {
            executor.shutdown();
        }
    }

    // Barnes-Hut methods
    private void generateOcttree() {
        rootCell = new Cell(CommonCore.XBOUND, CommonCore.YBOUND, CommonCore.ZBOUND);
        rootCell.index = 0;
        rootCell.x = 0;
        rootCell.y = 0;
        rootCell.z = 0;

        for (int i = 1; i < CommonCore.N; i++) {
            Cell cell = rootCell;
            while (cell.noSubcells != 0) {
                int sc = locateSubcell(cell, i);
                cell = cell.subcells[sc];
            }
            addToCell(cell, i);
        }
    }

    private void setLocationOfSubcells(Cell cell, double width, double height, double depth) {
        cell.subcells[0].x = cell.x;
        cell.subcells[0].y = cell.y;
        cell.subcells[0].z = cell.z;

        cell.subcells[1].x = cell.x + width;
        cell.subcells[1].y = cell.y;
        cell.subcells[1].z = cell.z;

        cell.subcells[2].x = cell.x + width;
        cell.subcells[2].y = cell.y;
        cell.subcells[2].z = cell.z + depth;

        cell.subcells[3].x = cell.x;
        cell.subcells[3].y = cell.y;
        cell.subcells[3].z = cell.z + depth;

        cell.subcells[4].x = cell.x;
        cell.subcells[4].y = cell.y + height;
        cell.subcells[4].z = cell.z;

        cell.subcells[5].x = cell.x + width;
        cell.subcells[5].y = cell.y + height;
        cell.subcells[5].z = cell.z;

        cell.subcells[6].x = cell.x + width;
        cell.subcells[6].y = cell.y + height;
        cell.subcells[6].z = cell.z + depth;

        cell.subcells[7].x = cell.x;
        cell.subcells[7].y = cell.y + height;
        cell.subcells[7].z = cell.z + depth;
    }

    private void generateSubcells(Cell cell) {
        double width = cell.width / 2.0;
        double height = cell.height / 2.0;
        double depth = cell.depth / 2.0;

        cell.noSubcells = 8;

        for (int i = 0; i < cell.noSubcells; i++) {
            cell.subcells[i] = new Cell(width, height, depth);
        }

        setLocationOfSubcells(cell, width, height, depth);
    }

    private int locateSubcell(Cell cell, int index) {
        if (position[index].px > cell.subcells[6].x) {
            if (position[index].py > cell.subcells[6].y) {
                return position[index].pz > cell.subcells[6].z ? 6 : 5;
            } else {
                return position[index].pz > cell.subcells[6].z ? 2 : 1;
            }
        } else {
            if (position[index].py > cell.subcells[6].y) {
                return position[index].pz > cell.subcells[6].z ? 7 : 4;
            } else {
                return position[index].pz > cell.subcells[6].z ? 3 : 0;
            }
        }
    }

    private void addToCell(Cell cell, int index) {
        if (cell.index == -1) {
            cell.index = index;
            return;
        }

        generateSubcells(cell);

        int sc1 = locateSubcell(cell, cell.index);
        cell.subcells[sc1].index = cell.index;

        int sc2 = locateSubcell(cell, index);

        if (sc1 == sc2) {
            addToCell(cell.subcells[sc1], index);
        } else {
            cell.subcells[sc2].index = index;
        }
    }

    private Cell computeCellProperties(Cell cell) {
        if (cell.noSubcells == 0) {
            if (cell.index != -1) {
                cell.mass = mass[cell.index];
                return cell;
            }
        } else {
            double tx = 0, ty = 0, tz = 0;
            for (int i = 0; i < cell.noSubcells; i++) {
                Cell temp = computeCellProperties(cell.subcells[i]);
                if (temp != null) {
                    cell.mass += temp.mass;
                    tx += position[temp.index].px * temp.mass;
                    ty += position[temp.index].py * temp.mass;
                    tz += position[temp.index].pz * temp.mass;
                }
            }

            if (cell.mass > 0) {
                cell.cx = tx / cell.mass;
                cell.cy = ty / cell.mass;
                cell.cz = tz / cell.mass;
                return cell;
            }
        }
        return null;
    }

    private void computeForceFromCell(Cell cell, int index) {
        double d = computeDistance(position[index], position[cell.index]);
        double f = (CommonCore.G * (mass[index] * mass[cell.index]) / (Math.pow(d, 2.0)));

        force[index].fx += f * ((position[cell.index].px - position[index].px) / d);
        force[index].fy += f * ((position[cell.index].py - position[index].py) / d);
        force[index].fz += f * ((position[cell.index].pz - position[index].pz) / d);
    }

    private void computeForceFromOcttree(Cell cell, int index) {
        if (cell.noSubcells == 0) {
            if (cell.index != -1 && cell.index != index) {
                computeForceFromCell(cell, index);
            }
        } else {
            double d = computeDistance(position[index], position[cell.index]);

            if (CommonCore.THETA > (cell.width / d)) {
                computeForceFromCell(cell, index);
            } else {
                for (int i = 0; i < cell.noSubcells; i++) {
                    computeForceFromOcttree(cell.subcells[i], index);
                }
            }
        }
    }

    public static void main(String[] args) {
        LocallyThreadedNBodyCalculator simulation = new LocallyThreadedNBodyCalculator();
        simulation.run();
    }
}