import data.*;
import mpi.*;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DistributedNBodyCalculator {
    private List<Point> points;
    private Quadrant rootQuadrant;

    private final int rank;
    private final int partSize;
    private final int pindex;

    public DistributedNBodyCalculator(String[] args) throws MPIException {
        try {
            MPI.Init(args);

            rank = MPI.COMM_WORLD.Rank();
            final var size = MPI.COMM_WORLD.Size();

            if (size < 1) {
                throw new MPIException("Need at least one process");
            }

            partSize = CommonCore.bodyCount / size;
            pindex = rank * partSize;

            if (rank == 0) {
                System.out.println("Initialization complete:");
                System.out.println("Bodies = " + CommonCore.bodyCount);
                System.out.println("Iterations = " + CommonCore.iterationCount);
                System.out.println("Processes = " + size);
                System.out.println("Particles per process = " + partSize);
            }

            initializeArrays();
        } catch (Exception e) {
            System.err.println("Initialization error: " + e.getMessage());
            MPI.COMM_WORLD.Abort(1);
            throw e;
        }
    }

    private void initializeArrays() {
        points = new ArrayList<>(CommonCore.bodyCount);
        for (int i = 0; i < CommonCore.bodyCount; i++) {
            points.add(new Point());
        }
    }

    private void initializeSpace() {
        for (int i = 0; i < CommonCore.bodyCount; i++) {
            points.set(i, RandomGenerator.nextPoint());
        }
    }

    private void computeVelocity() {
        for (var i = 0; i < partSize; i++) {
            final var point = points.get(i + pindex);

            if (point.mass().value() <= 0) {
                continue;
            }

            final var horizontalVelocity = point.velocity().horizontal() +
                    (point.force().horizontal() / point.mass().value()) *
                            CommonCore.timeIncrement;
            final var verticalVelocity = point.velocity().vertical() +
                    (point.force().vertical() / point.mass().value()) *
                            CommonCore.timeIncrement;

            double boundedVX = Math.min(Math.max(horizontalVelocity, -CommonCore.maxVelocity), 
                    CommonCore.maxVelocity);
            double boundedVY = Math.min(Math.max(verticalVelocity, -CommonCore.maxVelocity), 
                    CommonCore.maxVelocity);

            if (Double.isNaN(boundedVX) || Double.isNaN(boundedVY)) {
                System.out.println("NaN velocity detected for particle " + (i + pindex));
                continue;
            }

            point.setVelocity(new Velocity(boundedVX, boundedVY));
        }
    }

    private void computePositions() {
        IntStream.range(0, partSize).map(i -> i + pindex).mapToObj(i -> points.get(i)).forEach(
                point -> {
                    final var newPX = point.position().horizontal() +
                            point.velocity().horizontal() * CommonCore.timeIncrement;
                    final var newPY = point.position().vertical() +
                            point.velocity().vertical() * CommonCore.timeIncrement;

                    final var newVX = point.velocity().horizontal() * (
                            ((newPX + point.radius().value()) >= (CommonCore.maxWidth - CommonCore.boundingBoxTolerance) ||
                                    (newPX - point.radius().value()) <= 0) ? -1 : 1);
                    final var newVY = point.velocity().vertical() * (
                            ((newPY + point.radius().value()) >= (CommonCore.maxHeight - CommonCore.boundingBoxTolerance) ||
                                    (newPY - point.radius().value()) <= 0) ? -1 : 1);

                    double boundedPX = Math.min(Math.max(newPX, point.radius().value()), 
                            CommonCore.maxWidth - point.radius().value());
                    double boundedPY = Math.min(Math.max(newPY, point.radius().value()), 
                            CommonCore.maxHeight - point.radius().value());

                    point.setPosition(new Position(boundedPX, boundedPY));
                    point.setVelocity(new Velocity(newVX, newVY));
                }
        );
    }

    private void writePositions(int iteration) {
        if (rank != 0) {
            return;
        }
        CommonCore.write(iteration, points, CommonCore.Mode.DISTRIBUTED);
    }

    public void run() throws MPIException {
        try {
            if (rank == 0) {
                System.out.printf("\nRunning simulation for %d bodies with %d iterations, and time increment = %f..%n%n",
                        CommonCore.bodyCount, CommonCore.iterationCount, CommonCore.timeIncrement);
                initializeSpace();
                writePositions(0);
            }

            Point[] pointsArray = points.toArray(new Point[0]);
            MPI.COMM_WORLD.Bcast(pointsArray, 0, CommonCore.bodyCount, MPI.OBJECT, 0);
            for (int i = 0; i < CommonCore.bodyCount; i++) {
                points.set(i, pointsArray[i]);
            }

            for (var i = 0; i < CommonCore.iterationCount; i++) {
                generateQuadtree();
                computeForceFromQuadtree();
                updatePositions();
                ensureUniquePositions();

                Point[] localPoints = points.subList(rank * partSize, (rank + 1) * partSize)
                        .toArray(new Point[0]);
                Point[] gatheredPoints = new Point[CommonCore.bodyCount];
                
                MPI.COMM_WORLD.Allgather(
                        localPoints, 0, partSize, MPI.OBJECT,
                        gatheredPoints, 0, partSize, MPI.OBJECT
                );

                for (int j = 0; j < CommonCore.bodyCount; j++) {
                    points.set(j, gatheredPoints[j]);
                }

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
        rootQuadrant = new Quadrant(
                new Dimension(CommonCore.maxWidth, CommonCore.maxHeight),
                new Position()
        );
        // Insert all points into the quadtree
        points.forEach(point -> rootQuadrant.insert(point));
    }

    private void computeForceFromQuadtree() {
        if (rootQuadrant == null) {
            return;
        }
        
        // Only compute forces for this process's partition
        IntStream.range(0, partSize)
                .mapToObj(i -> points.get(i + pindex))
                .forEach(point -> rootQuadrant.addForceActingOn(point));
    }

    private void updatePositions() {
        // Update positions only for this process's partition
        IntStream.range(0, partSize)
                .mapToObj(i -> points.get(i + pindex))
                .forEach(Point::updatePosition);
    }

    private void ensureUniquePositions() {
        // Only check this process's partition
        for (var i = pindex; i < pindex + partSize; i++) {
            for (var j = i + 1; j < CommonCore.bodyCount; j++) {
                if (points.get(i).position().equals(points.get(j).position())) {
                    final var position = points.get(j).position();
                    points.get(j).setPosition(new Position(
                            position.horizontal() + 2 * RandomGenerator.nextDouble() - 1,
                            position.vertical() + 2 * RandomGenerator.nextDouble() - 1
                    ));
                }
            }
        }
    }
}