import mpi.*;
import java.io.*;
import java.util.Random;

public class DistributedNBodyCalculator {
    private Position[] position;   // Current positions for all particles
    private Velocity[] ivelocity;  // Initial velocity for all particles
    private Velocity[] velocity;   // Velocity of particles in current processor
    private double[] mass;         // Mass of each particle
    private double[] radius;       // Radius of each particle
    private Force[] force;         // Force experienced by all particles
    private Cell rootCell;         // Root of BH octtree

    private int rank;           // Rank of process
    private int size;           // Number of processes in the group
    private int partSize;       // Number of particles each processor is responsible for
    private int pindex;         // Index pointer for current processor's data

    // Inner classes for data structures
    static class Position implements Serializable {
        private static final long serialVersionUID = 1L;
        double px, py, pz;
    }

    static class Velocity implements Serializable {
        private static final long serialVersionUID = 1L;
        double vx, vy, vz;
    }

    static class Force implements Serializable {
        private static final long serialVersionUID = 1L;
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
    public DistributedNBodyCalculator(String[] args) throws MPIException {
        try {
            MPI.Init(args);

            rank = MPI.COMM_WORLD.Rank();
            size = MPI.COMM_WORLD.Size();

            if (size < 1) {
                throw new MPIException("Need at least one process");
            }


            if (rank == 0) {
                System.out.println("Using hardcoded values:");
                System.out.println("N = " + CommonCore.bodyCount);
                System.out.println("TIME = " + CommonCore.iterationCount);
            }

            // Validate that CommonCore.N is divisible by the number of processes
            if (CommonCore.bodyCount % size != 0) {
                if (rank == 0) {
                    System.out.println("Warning: Adjusting CommonCore.N to be divisible by number of processes");
                }
            }

            partSize = CommonCore.bodyCount / size;
            pindex = rank * partSize;

            if (rank == 0) {
                System.out.println("Initialization complete:");
                System.out.println("CommonCore.N = " + CommonCore.bodyCount);
                System.out.println("CommonCore.TIME = " + CommonCore.iterationCount);
                System.out.println("Processes = " + size);
                System.out.println("Particles per process = " + partSize);
            }

            initializeArrays();
        } catch (Exception e) {
            System.err.println("Initialization error: " + e.getMessage());
            e.printStackTrace();
            MPI.COMM_WORLD.Abort(1);
            throw e;
        }
    }

    private void initializeArrays() {
        mass = new double[CommonCore.bodyCount];
        radius = new double[CommonCore.bodyCount];
        position = new Position[CommonCore.bodyCount];
        ivelocity = new Velocity[CommonCore.bodyCount];
        velocity = new Velocity[partSize];
        force = new Force[partSize];

        // Initialize Position and Velocity objects
        for (int i = 0; i < CommonCore.bodyCount; i++) {
            position[i] = new Position();
            ivelocity[i] = new Velocity();
        }
        for (int i = 0; i < partSize; i++) {
            velocity[i] = new Velocity();
            force[i] = new Force();
        }
    }

    // Main simulation methods would follow here...
    // I'll continue with more implementation in subsequent parts

    public static void main(String[] args) {
        try {
            DistributedNBodyCalculator simulation = new DistributedNBodyCalculator(args);
            simulation.run();
        } catch (Exception e) {
            System.err.println("Fatal error: " + e.getMessage());
            System.exit(1);
        } finally {
            try {
                MPI.Finalize();
            } catch (MPIException e) {
                System.err.println("Error during MPI finalization: " + e.getMessage());
                System.exit(1);
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
        // Inner bounds to prevent generating particles outside space boundaries
        double ixbound = CommonCore.maxWidth - CommonCore.maxBodyRadius;
        double iybound = CommonCore.maxHeight - CommonCore.maxBodyRadius;
        double izbound = CommonCore.maxDepth - CommonCore.maxBodyRadius;

        for (int i = 0; i < CommonCore.bodyCount; i++) {
            mass[i] = CommonCore.defaultMass * generateRand();
            radius[i] = CommonCore.maxBodyRadius * generateRand();
            position[i].px = generateRand() * ixbound;
            position[i].py = generateRand() * iybound;
            position[i].pz = generateRand() * izbound;
            ivelocity[i].vx = generateRandEx();
            ivelocity[i].vy = generateRandEx();
            ivelocity[i].vz = generateRandEx();
        }
    }

    private boolean checkCollision(int index1, int index2) {
        return Math.pow(position[index1].px - position[index2].px, 2.0) +
                Math.pow(position[index1].py - position[index2].py, 2.0) +
                Math.pow(position[index1].pz - position[index2].pz, 2.0) <
                Math.pow(radius[index1] + radius[index2], 2.0);
    }

    private double computeDistance(Position a, Position b) {
        return Math.sqrt(Math.pow(a.px - b.px, 2.0) +
                Math.pow(a.py - b.py, 2.0) +
                Math.pow(a.pz - b.pz, 2.0));
    }

    private void reinitializeRadius() {
        for (int i = 0; i < CommonCore.bodyCount; i++) {
            for (int j = i + 1; j < CommonCore.bodyCount; j++) {
                if (checkCollision(i, j)) {
                    double d = computeDistance(position[i], position[j]);
                    radius[i] = radius[j] = d/2.0;
                }
            }
        }
    }

    private void computeForce() {
        for (int i = 0; i < partSize; i++) {
            force[i].fx = 0.0;
            force[i].fy = 0.0;
            force[i].fz = 0.0;

            for (int j = 0; j < CommonCore.bodyCount; j++) {
                if (j == (i + pindex)) continue;

                double d = computeDistance(position[i + pindex], position[j]);

                // Compute gravitational force according to Newton's law
                double f = (CommonCore.GravitationalConstant * (mass[i + pindex] * mass[j]) / (Math.pow(d, 2.0)));

                // Resolve forces in each direction
                force[i].fx += f * ((position[j].px - position[i + pindex].px) / d);
                force[i].fy += f * ((position[j].py - position[i + pindex].py) / d);
                force[i].fz += f * ((position[j].pz - position[i + pindex].pz) / d);
            }
        }
    }

    private void computeVelocity() {
        for (int i = 0; i < partSize; i++) {
            velocity[i].vx += (force[i].fx / mass[i + pindex]) * CommonCore.timeIncrement;
            velocity[i].vy += (force[i].fy / mass[i + pindex]) * CommonCore.timeIncrement;
            velocity[i].vz += (force[i].fz / mass[i + pindex]) * CommonCore.timeIncrement;
        }
    }

    private void computePositions() {
        for (int i = 0; i < partSize; i++) {
            position[i + pindex].px += velocity[i].vx * CommonCore.timeIncrement;
            position[i + pindex].py += velocity[i].vy * CommonCore.timeIncrement;
            position[i + pindex].pz += velocity[i].vz * CommonCore.timeIncrement;

            // Check if particles attempt to cross boundary
            if ((position[i + pindex].px + radius[i + pindex]) >= CommonCore.maxWidth ||
                    (position[i + pindex].px - radius[i + pindex]) <= 0) {
                velocity[i].vx *= -1;
            }
            else if ((position[i + pindex].py + radius[i + pindex] >= CommonCore.maxHeight) ||
                    (position[i + pindex].py - radius[i + pindex]) <= 0) {
                velocity[i].vy *= -1;
            }
            else if ((position[i + pindex].pz + radius[i + pindex]) >= CommonCore.maxDepth ||
                    (position[i + pindex].pz - radius[i + pindex]) <= 0) {
                velocity[i].vz *= -1;
            }
        }
    }

    private void writePositions(int iteration) {
        if (rank == 0) {
            String filename = String.format("positions_iter_%d.csv", iteration);
            try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)))) {
                // Write CSV header
                writer.println("particle_id,x,y,z");

                // Write each particle's position
                for (int i = 0; i < CommonCore.bodyCount; i++) {
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
    }

    public void run() throws MPIException {
        try {
            if (rank == 0) {
                System.out.printf("\nRunning simulation for %d bodies with %d iterations, and DELTAT = %f..%n%n",
                        CommonCore.bodyCount, CommonCore.iterationCount, CommonCore.timeIncrement);
                initializeSpace();

                // Write initial positions
                writePositions(0);
            }

            // Broadcast mass and position arrays
            MPI.COMM_WORLD.Bcast(mass, 0, CommonCore.bodyCount, MPI.DOUBLE, 0);
            MPI.COMM_WORLD.Bcast(position, 0, CommonCore.bodyCount, MPI.OBJECT, 0);

            // Scatter initial velocities
            MPI.COMM_WORLD.Scatter(
                    ivelocity, 0, partSize, MPI.OBJECT,
                    velocity, 0, partSize, MPI.OBJECT, 0
            );

            // Main simulation loop
            for (int i = 0; i < CommonCore.iterationCount; i++) {
                generateOcttree();
                computeCellProperties(rootCell);
                computeForceFromOcttree();
                deleteOcttree(rootCell);

                computeVelocity();
                computePositions();

                // Gather updated positions
                MPI.COMM_WORLD.Allgather(
                        position, rank * partSize, partSize, MPI.OBJECT,
                        position, 0, partSize, MPI.OBJECT
                );

                // Write positions for this iteration
                writePositions(i + 1);

                // Optional: Print progress
                if (rank == 0 && (i + 1) % 50 == 0) {
                    System.out.printf("Completed iteration %d of %d%n", i + 1, CommonCore.iterationCount);
                }
            }

        } catch (MPIException e) {
            System.err.println("MPI Error during simulation: " + e.getMessage());
            MPI.COMM_WORLD.Abort(1);
            throw e;
        }
    }

    private void generateOcttree() {
        // Initialize root of octtree
        rootCell = new Cell(CommonCore.maxWidth, CommonCore.maxHeight, CommonCore.maxDepth);
        rootCell.index = 0;
        rootCell.x = 0;
        rootCell.y = 0;
        rootCell.z = 0;

        // Add remaining particles to the tree
        for (int i = 1; i < CommonCore.bodyCount; i++) {
            Cell cell = rootCell;

            // Find which node to add the body to
            while (cell.noSubcells != 0) {
                int sc = locateSubcell(cell, i);
                cell = cell.subcells[sc];
            }

            addToCell(cell, i);
        }
    }

    private void setLocationOfSubcells(Cell cell, double width, double height, double depth) {
        // Set location of new cells
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
        // Calculate subcell dimensions
        double width = cell.width / 2.0;
        double height = cell.height / 2.0;
        double depth = cell.depth / 2.0;

        // Cell no longer a leaf
        cell.noSubcells = 8;

        // Create and initialize new subcells
        for (int i = 0; i < cell.noSubcells; i++) {
            cell.subcells[i] = new Cell(width, height, depth);
        }

        setLocationOfSubcells(cell, width, height, depth);
    }

    private int locateSubcell(Cell cell, int index) {
        // Determine which subcell to add the body to
        if (position[index].px > cell.subcells[6].x) {
            if (position[index].py > cell.subcells[6].y) {
                if (position[index].pz > cell.subcells[6].z)
                    return 6;
                else
                    return 5;
            } else {
                if (position[index].pz > cell.subcells[6].z)
                    return 2;
                else
                    return 1;
            }
        } else {
            if (position[index].py > cell.subcells[6].y) {
                if (position[index].pz > cell.subcells[6].z)
                    return 7;
                else
                    return 4;
            } else {
                if (position[index].pz > cell.subcells[6].z)
                    return 3;
                else
                    return 0;
            }
        }
    }

    private void addToCell(Cell cell, int index) {
        if (cell.index == -1) {
            cell.index = index;
            return;
        }

        generateSubcells(cell);

        // The current cell's body must now be re-added to one of its subcells
        int sc1 = locateSubcell(cell, cell.index);
        cell.subcells[sc1].index = cell.index;

        // Locate subcell for new body
        int sc2 = locateSubcell(cell, index);

        if (sc1 == sc2)
            addToCell(cell.subcells[sc1], index);
        else
            cell.subcells[sc2].index = index;
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

            // Compute center of mass
            cell.cx = tx / cell.mass;
            cell.cy = ty / cell.mass;
            cell.cz = tz / cell.mass;

            return cell;
        }
        return null;
    }

    private void computeForceFromCell(Cell cell, int index) {
        double d = computeDistance(position[index], position[cell.index]);

        // Compute gravitational force according to Newton's law
        double f = (CommonCore.GravitationalConstant * (mass[index] * mass[cell.index]) / (Math.pow(d, 2.0)));

        // Resolve forces in each direction
        force[index - pindex].fx += f * ((position[cell.index].px - position[index].px) / d);
        force[index - pindex].fy += f * ((position[cell.index].py - position[index].py) / d);
        force[index - pindex].fz += f * ((position[cell.index].pz - position[index].pz) / d);
    }

    private void computeForceFromOcttree(Cell cell, int index) {
        if (cell.noSubcells == 0) {
            if (cell.index != -1 && cell.index != index) {
                computeForceFromCell(cell, index);
            }
        } else {
            double d = computeDistance(position[index], position[cell.index]);

            if (CommonCore.angle > (cell.width / d)) {
                // Use approximation
                computeForceFromCell(cell, index);
            } else {
                for (int i = 0; i < cell.noSubcells; i++) {
                    computeForceFromOcttree(cell.subcells[i], index);
                }
            }
        }
    }

    private void computeForceFromOcttree() {
        for (int i = 0; i < partSize; i++) {
            force[i].fx = 0.0;
            force[i].fy = 0.0;
            force[i].fz = 0.0;

            computeForceFromOcttree(rootCell, i + pindex);
        }
    }

    private void deleteOcttree(Cell cell) {
        if (cell.noSubcells == 0) {
            return;
        }

        for (int i = 0; i < cell.noSubcells; i++) {
            deleteOcttree(cell.subcells[i]);
        }
    }

    private void synchronizeProcesses() throws MPIException {
        MPI.COMM_WORLD.Barrier();
    }
}