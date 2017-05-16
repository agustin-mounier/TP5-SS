package IntegrationMethod;

/**
 * Created by amounier on 4/21/17.
 */
public class Beeman {

    private GranularForce f;
    private double dT;
    private Particle particle;

    private double prevAcelX = 0;
    private double prevAcelY = 0;

    public Beeman(GranularForce f, double dT) {
        this.f = f;
        this.dT = dT;
        this.particle = f.getParticle();
    }

    public Particle moveParticle() {
        double aX = f.getXForce() / particle.mass;
        double aY = f.getYForce() / particle.mass;

        double newX = particle.x + (dT * particle.velX) + (2.0 / 3) * aX * (dT * dT) - (1.0 / 6) * prevAcelX * (dT * dT);
        double newY = particle.y + (dT * particle.velY) + (2.0 / 3) * aY * (dT * dT) - (1.0 / 6) * prevAcelY * (dT * dT);

        double predictedVelX = particle.velX + (3.0 / 2) * aX * dT - (1.0 / 2) * prevAcelX * dT;
        double predictedVelY = particle.velY + (3.0 / 2) * aY * dT - (1.0 / 2) * prevAcelY * dT;

        double currVelX = particle.velX;
        double currVelY = particle.velY;

        particle.x = newX;
        particle.y = newY;
        particle.velX = predictedVelX;
        particle.velY = predictedVelY;

        double correctedVelX = currVelX + (1.0 / 3) * (f.getXForce() / particle.mass) * dT + (5.0 / 6) * aX * dT - (1.0 / 6) * prevAcelX * dT;
        double correctedVelY = currVelY + (1.0 / 3) * (f.getYForce() / particle.mass) * dT + (5.0 / 6) * aY * dT - (1.0 / 6) * prevAcelY * dT;

        particle.velX = correctedVelX;
        particle.velY = correctedVelY;

        prevAcelX = aX;
        prevAcelY = aY;

        return particle;
    }

    public GranularForce getForce() {
        return f;
    }

    public void reset() {
        prevAcelX = 0;
        prevAcelY = 0;
        f.reset();
    }
}
