import mpi.MPI;
import mpi.MPIException;

public class Main {
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
}