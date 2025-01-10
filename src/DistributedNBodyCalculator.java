import data.*;
import mpi.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class DistributedNBodyCalculator {
    private List<Point> points;
    private Quadrant rootCell;

    private int rank;
    private int size;
    private int partSize;
    private int pindex;

    private Random random = new Random();

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
        points = new ArrayList<>(CommonCore.bodyCount);
        for (int i = 0; i < CommonCore.bodyCount; i++) {
            points.add(new Point(
                    new Position(0, 0),
                    0,
                    new Mass(0),
                    new Velocity(0, 0)
            ));
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

        for (int i = 0; i < CommonCore.bodyCount; i++) {
            double mass = Math.max(CommonCore.defaultMass * generateRand(), CommonCore.defaultMass * 0.1);
            double radius = Math.max(CommonCore.maxBodyRadius * generateRand(), CommonCore.maxBodyRadius * 0.1);
            
            Position position = new Position(
                    CommonCore.maxBodyRadius + (generateRand() * (ixbound - 2 * CommonCore.maxBodyRadius)),
                    CommonCore.maxBodyRadius + (generateRand() * (iybound - 2 * CommonCore.maxBodyRadius))
            );
            
            double maxInitialVelocity = 100.0;
            Velocity velocity = new Velocity(
                    generateRandEx() * maxInitialVelocity,
                    generateRandEx() * maxInitialVelocity
            );

            points.set(i, new Point(position, radius, new Mass(mass), velocity));
        }
    }

    private boolean checkCollision(Point p1, Point p2) {
        return Math.pow(p1.position().horizontal() - p2.position().horizontal(), 2.0) +
                Math.pow(p1.position().vertical() - p2.position().vertical(), 2.0) <
                Math.pow(p1.radius() + p2.radius(), 2.0);
    }

    private double computeDistance(Position a, Position b) {
        return Math.sqrt(Math.pow(a.horizontal() - b.horizontal(), 2.0) +
                Math.pow(a.vertical() - b.vertical(), 2.0));
    }

    private void computeVelocity() {
        for (int i = 0; i < partSize; i++) {
            Point point = points.get(i + pindex);

            if (point.mass().value() <= 0) {
                continue;
            }

            double vx = point.velocity().horizontal() +
                    (point.force().horizontal() / point.mass().value()) *
                            CommonCore.timeIncrement;
            double vy = point.velocity().vertical() +
                    (point.force().vertical() / point.mass().value()) *
                            CommonCore.timeIncrement;

            if (Double.isNaN(vx) || Double.isNaN(vy)) {
                continue;
            }
            points.set(i + pindex, new Point(
                    point.position(),
                    point.radius(),
                    point.mass(),
                    new Velocity(vx, vy),
                    point.force()
            ));
        }
    }

    private void computePositions() {
        for (int i = 0; i < partSize; i++) {
            Point point = points.get(i + pindex);
            double px = point.position().horizontal() +
                    point.velocity().horizontal() * CommonCore.timeIncrement;
            double py = point.position().vertical() +
                    point.velocity().vertical() * CommonCore.timeIncrement;

            Velocity newVelocity = point.velocity();

            if ((px + point.radius()) >= CommonCore.maxWidth ||
                    (px - point.radius()) <= 0) {
                newVelocity = new Velocity(-point.velocity().horizontal(),
                        point.velocity().vertical());
            }
            if ((py + point.radius() >= CommonCore.maxHeight) ||
                    (py - point.radius()) <= 0) {
                newVelocity = new Velocity(point.velocity().horizontal(),
                        -point.velocity().vertical());
            }

            points.set(i + pindex, new Point(
                    new Position(px, py),
                    point.radius(),
                    point.mass(),
                    newVelocity,
                    point.force()
            ));
        }
    }

    private void writePositions(int iteration) {
        if (rank == 0) {
            String filename = String.format("src/visualise/csv/positions_iter_%d.csv", iteration);
            try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(filename)))) {
                writer.println("particle_id,x,y");
                for (int i = 0; i < CommonCore.bodyCount; i++) {
                    Point point = points.get(i);
                    writer.printf("%d,%.6f,%.6f%n",
                            i,
                            point.position().horizontal(),
                            point.position().vertical());
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
                writePositions(0);
            }

            MPI.COMM_WORLD.Bcast(points.toArray(), 0, CommonCore.bodyCount, MPI.OBJECT, 0);

            for (int i = 0; i < CommonCore.iterationCount; i++) {
                generateQuadtree();
                computeQuadrantProperties(rootCell);
                computeForceFromQuadtree();
                computeVelocity();
                computePositions();

                MPI.COMM_WORLD.Allgather(
                        points.subList(rank * partSize, (rank + 1) * partSize).toArray(),
                        0, partSize, MPI.OBJECT,
                        points.toArray(), 0, partSize, MPI.OBJECT
                );

                writePositions(i + 1);

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

    private void generateQuadtree() {
        rootCell = new Quadrant(
                0,
                new Dimension(CommonCore.maxWidth, CommonCore.maxHeight),
                new Position(0, 0)
        );

        for (int i = 1; i < CommonCore.bodyCount; i++) {
            Quadrant quadrant = rootCell;

            while (quadrant.innerQuadrantCount() != 0) {
                int subQuadrant = locateSubQuadrant(quadrant, i);
                quadrant = quadrant.innerQuadrants().get(subQuadrant);
            }

            addToQuadrant(quadrant, i);
        }
    }

    private int locateSubQuadrant(Quadrant quadrant, int index) {
        Point point = points.get(index);
        double midX = quadrant.bottomLeftCorner().horizontal() + quadrant.dimensions().horizontal() / 2;
        double midY = quadrant.bottomLeftCorner().vertical() + quadrant.dimensions().vertical() / 2;

        if (point.position().horizontal() > midX) {
            return point.position().vertical() > midY ? 1 : 0;
        } else {
            return point.position().vertical() > midY ? 2 : 3;
        }
    }

    private void addToQuadrant(Quadrant quadrant, int index) {
        if (quadrant.index() == -1) {
            quadrant.setIndex(index);
            return;
        }

        if (quadrant.innerQuadrantCount() == 0) {
            double width = quadrant.dimensions().horizontal() / 2;
            double height = quadrant.dimensions().vertical() / 2;

            if (width < CommonCore.maxBodyRadius || height < CommonCore.maxBodyRadius) {
                return;
            }

            Position bottomLeft = quadrant.bottomLeftCorner();

            List<Quadrant> subQuadrants = new ArrayList<>();
            subQuadrants.add(new Quadrant(-1, new Dimension(width, height),
                    new Position(bottomLeft.horizontal() + width, bottomLeft.vertical())));
            subQuadrants.add(new Quadrant(-1, new Dimension(width, height),
                    new Position(bottomLeft.horizontal() + width, bottomLeft.vertical() + height)));
            subQuadrants.add(new Quadrant(-1, new Dimension(width, height),
                    new Position(bottomLeft.horizontal(), bottomLeft.vertical() + height)));
            subQuadrants.add(new Quadrant(-1, new Dimension(width, height), bottomLeft));

            quadrant.setInnerQuadrants(subQuadrants);
            quadrant.setInnerQuadrantCount(4);

            int sq1 = locateSubQuadrant(quadrant, quadrant.index());
            quadrant.innerQuadrants().get(sq1).setIndex(quadrant.index());

            int sq2 = locateSubQuadrant(quadrant, index);
            if (sq1 == sq2) {
                addToQuadrant(quadrant.innerQuadrants().get(sq1), index);
            } else {
                quadrant.innerQuadrants().get(sq2).setIndex(index);
            }
        } else {
            int sq = locateSubQuadrant(quadrant, index);
            addToQuadrant(quadrant.innerQuadrants().get(sq), index);
        }
    }

    private Quadrant computeQuadrantProperties(Quadrant quadrant) {
        if (quadrant.innerQuadrantCount() == 0) {
            if (quadrant.index() != -1) {
                return new Quadrant(
                        quadrant.index(),
                        0,
                        points.get(quadrant.index()).mass(),
                        quadrant.bottomLeftCorner(),
                        points.get(quadrant.index()).position(),
                        quadrant.dimensions(),
                        quadrant.innerQuadrants()
                );
            }
        } else {
            double totalMass = 0;
            double weightedX = 0;
            double weightedY = 0;

            for (Quadrant subQuadrant : quadrant.innerQuadrants()) {
                Quadrant computed = computeQuadrantProperties(subQuadrant);
                if (computed != null) {
                    totalMass += computed.mass().value();
                    weightedX += points.get(computed.index()).position().horizontal() * computed.mass().value();
                    weightedY += points.get(computed.index()).position().vertical() * computed.mass().value();
                }
            }

            return new Quadrant(
                    quadrant.index(),
                    quadrant.innerQuadrantCount(),
                    new Mass(totalMass),
                    quadrant.bottomLeftCorner(),
                    new Position(weightedX / totalMass, weightedY / totalMass),
                    quadrant.dimensions(),
                    quadrant.innerQuadrants()
            );
        }
        return null;
    }

    private void computeForceFromQuadrant(Quadrant quadrant, int index) {
        Point point = points.get(index);
        Point otherPoint = points.get(quadrant.index());

        double d = computeDistance(point.position(), otherPoint.position());
        
        if (d < CommonCore.maxBodyRadius || d <= 0) {
            return;
        }

        double f = (CommonCore.GravitationalConstant *
                (point.mass().value() * otherPoint.mass().value()) /
                (Math.pow(d, 2.0)));

        double dx = otherPoint.position().horizontal() - point.position().horizontal();
        double dy = otherPoint.position().vertical() - point.position().vertical();
        
        double fx = f * (dx / d);
        double fy = f * (dy / d);

        if (Double.isNaN(fx) || Double.isNaN(fy)) {
            return;
        }

        points.set(index, new Point(
                point.position(),
                point.radius(),
                point.mass(),
                point.velocity(),
                new Force(
                        point.force().horizontal() + fx,
                        point.force().vertical() + fy
                )
        ));
    }

    private void computeForceFromQuadtree(Quadrant quadrant, int index) {
        if (quadrant == null || quadrant.index() == -1) {
            return;
        }

        if (quadrant.innerQuadrantCount() == 0) {
            if (quadrant.index() != -1 && quadrant.index() != index) {
                computeForceFromQuadrant(quadrant, index);
            }
        } else {
            double d = computeDistance(
                    points.get(index).position(),
                    quadrant.centerOfMass()
            );

            if (d <= 0) {
                return;
            }

            if (CommonCore.angle > (quadrant.dimensions().horizontal() / d)) {
                computeForceFromQuadrant(quadrant, index);
            } else {
                for (Quadrant subQuadrant : quadrant.innerQuadrants()) {
                    computeForceFromQuadtree(subQuadrant, index);
                }
            }
        }
    }

    private void computeForceFromQuadtree() {
        for (int i = 0; i < partSize; i++) {
            Point point = points.get(i + pindex);
            points.set(i + pindex, new Point(
                    point.position(),
                    point.radius(),
                    point.mass(),
                    point.velocity(),
                    new Force(0, 0)
            ));

            if (rootCell != null) {
                computeForceFromQuadtree(rootCell, i + pindex);
            }
        }
    }

    private void deleteQuadtree(Quadrant quadrant) {
        if (quadrant.innerQuadrantCount() == 0) {
            return;
        }

        for (Quadrant subQuadrant : quadrant.innerQuadrants()) {
            deleteQuadtree(subQuadrant);
        }
    }

    private void synchronizeProcesses() throws MPIException {
        MPI.COMM_WORLD.Barrier();
    }
}