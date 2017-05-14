package integrationMethod;

/**
 * Created by amounier on 4/21/17.
 */
public class ParticleState {

    double x;
    double y;

    double velX;
    double velY;

    public ParticleState(double x, double y, double velX, double velY) {
        this.x = x;
        this.y = y;
        this.velX = velX;
        this.velY = velY;
    }

    @Override
    public String toString() {
        return x + "";
    }
}
