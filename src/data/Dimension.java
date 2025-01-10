package data;

public record Dimension (
        double horizontal,
        double vertical
){
    Dimension(){
        this(0, 0);
    }
}
