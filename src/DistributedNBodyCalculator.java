import data.*;
import mpi.*;

import java.io.*;
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

            if (rank == 0) {
                System.out.println("Bodies = " + CommonCore.bodyCount);
                System.out.println("Iterations = " + CommonCore.iterationCount);
            }

            if (CommonCore.bodyCount % size != 0) {
                if (rank == 0) {
                    System.out.println("Warning: Adjust body count to be divisible by number of processes");
                }
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
        points = IntStream.range(0, CommonCore.bodyCount).mapToObj(_ -> new Point()).toList();
    }

    private void initializeSpace() {
        for (final var point : points) {
            point.setMass(RandomGenerator.nextMass());
            point.setRadius(RandomGenerator.nextRadius());
            point.setPosition(RandomGenerator.nextPosition());
            point.setVelocity(RandomGenerator.nextVelocity());
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

            if (Double.isNaN(horizontalVelocity) || Double.isNaN(verticalVelocity)) {
                System.out.println("Got a nan?!");
                continue;
            }
            point.setVelocity(new Velocity(horizontalVelocity, verticalVelocity));
        }
    }

    private void computePositions() {
        IntStream.range(0, partSize).map(i -> i + pindex).mapToObj(i -> points.get(i)).forEach(
                point -> {
                    final var newPX = point.position().horizontal() +
                            point.velocity().horizontal() * CommonCore.timeIncrement;
                    final var newPY = point.position().vertical() +
                            point.velocity().vertical() * CommonCore.timeIncrement;
                    final var newVY = point.velocity().horizontal() * (
                            ((newPX + point.radius().value()) >= CommonCore.maxWidth ||
                                    (newPX - point.radius().value()) <= 0) ? -1 : 1);
                    final var newVX = point.velocity().vertical() * (
                            ((newPY + point.radius().value() >= CommonCore.maxHeight) ||
                                    (newPY - point.radius().value()) <= 0) ? -1 : 1);
                    point.setPosition(new Position(newPX, newPY));
                    point.setVelocity(new Velocity(newVY, newVX));

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

            MPI.COMM_WORLD.Bcast(points.toArray(), 0, CommonCore.bodyCount, MPI.OBJECT, 0);

            for (var i = 0; i < CommonCore.iterationCount; i++) {
                generateQuadtree();
                computeQuadrantProperties(rootQuadrant);
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
        rootQuadrant = new Quadrant(
                0,
                new Dimension(CommonCore.maxWidth, CommonCore.maxHeight),
                new Position()
        );
        for (var i = 1; i < CommonCore.bodyCount; i++) {
            var currentQuadrant = rootQuadrant;
            while (currentQuadrant.innerQuadrantCount() != 0) {
                final var subQuadrant = locateSubQuadrant(currentQuadrant, i);
                currentQuadrant = currentQuadrant.innerQuadrants().get(subQuadrant);
            }

            addToQuadrant(currentQuadrant, i);
        }
    }

    private int locateSubQuadrant(Quadrant quadrant, int index) {
        final var point = points.get(index);
        final var midX = quadrant.bottomLeftCorner().horizontal() + quadrant.dimensions().horizontal() / 2;
        final var midY = quadrant.bottomLeftCorner().vertical() + quadrant.dimensions().vertical() / 2;

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

            final var sq1 = locateSubQuadrant(quadrant, quadrant.index());
            quadrant.innerQuadrants().get(sq1).setIndex(quadrant.index());

            final var sq2 = locateSubQuadrant(quadrant, index);
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
            if (quadrant.index() == -1) {
                return null;
            }
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
        var totalMass = 0.0;
        var centerX = 0.0;
        var centerY = 0.0;

        for (Quadrant subQuadrant : quadrant.innerQuadrants()) {
            Quadrant computed = computeQuadrantProperties(subQuadrant);
            if (computed != null) {
                totalMass += computed.mass().value();
                centerX += points.get(computed.index()).position().horizontal() * computed.mass().value();
                centerY += points.get(computed.index()).position().vertical() * computed.mass().value();
            }
        }

        return new Quadrant(
                quadrant.index(),
                quadrant.innerQuadrantCount(),
                new Mass(totalMass),
                quadrant.bottomLeftCorner(),
                new Position(centerX / totalMass, centerY / totalMass),
                quadrant.dimensions(),
                quadrant.innerQuadrants()
        );

    }

    private void computeForceFromQuadtree() {
        if (rootQuadrant == null) {
            return;
        }
        IntStream.range(0, partSize)
                .map(i -> i + pindex)
                .mapToObj(i -> points.get(i)).forEach(
                        point -> rootQuadrant.addForceActingOn(point)
                );
    }
}