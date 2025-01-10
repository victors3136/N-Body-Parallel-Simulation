package data;

public record Position(
        double horizontal,
        double vertical
) {
    public Position(){
        this(0, 0);
    }
}