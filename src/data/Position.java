package data;

public record Position(
        double horizontal,
        double vertical
) {
    Position(){
        this(0, 0);
    }
}