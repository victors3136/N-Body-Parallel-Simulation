import mpi.MPI;
import mpi.MPIException;

public class Main {
    public static void main(String[] args) {
            DistributedNBodyCalculator simulation = new DistributedNBodyCalculator(args);
            simulation.run();
            MPI.Finalize();
    }
}