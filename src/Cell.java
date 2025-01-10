import java.util.List;

public class Cell {
    int index;
    int subcellCount;
    Mass mass;
    Position position;
    Position centerOfMass;
    Dimension size;
    List<Cell> subcells;

    Cell(Dimension size, Position position) {
        this.mass = new Mass(0.0);
        this.subcellCount = 0;
        this.index = -1;
        this.centerOfMass = new Position(0.0, 0.0, 0.0);
        this.size = size;
        this.position = position;
        this.subcells = List.of(new Cell[8]);
    }
    public void setPosition(Position position) {
        this.position = position;
    }
}
