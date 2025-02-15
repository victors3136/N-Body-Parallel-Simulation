import data.*;
import mpi.*;

import java.util.ArrayList;
import java.util.List;
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
            for (var i = 0; i < CommonCore.bodyCount; i++) {
                points.set(i, pointsArray[i]);
            }

            for (var iteration = 0; iteration < CommonCore.iterationCount; iteration++) {
                if(iteration % 50  == 0 || iteration == CommonCore.iterationCount - 1) {
                    CommonCore.displayProgress(iteration);
                }
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

                writePositions(iteration + 1);

                if (rank == 0 && (iteration + 1) % 50 == 0) {
                    System.out.printf("Completed iteration %d of %d%n", iteration + 1, CommonCore.iterationCount);
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
        points.forEach(point -> rootQuadrant.insert(point));
    }

    private void computeForceFromQuadtree() {
        if (rootQuadrant == null) {
            return;
        }

        IntStream.range(0, partSize)
                .mapToObj(i -> points.get(i + pindex))
                .forEach(point -> rootQuadrant.addForceActingOn(point));
    }

    private void updatePositions() {
        IntStream.range(0, partSize)
                .mapToObj(i -> points.get(i + pindex))
                .forEach(Point::updatePosition);
    }

    private void ensureUniquePositions() {
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