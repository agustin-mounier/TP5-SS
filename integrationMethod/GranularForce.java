package IntegrationMethod;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by amounier on 5/13/17.
 */
public class GranularForce {

    Particle p;
    private static final double kn = 1E5;
    private static final double kt = 2 * kn;
    private static final double G = -10;
    private Set<Particle> neighbours;

    public GranularForce(Particle p, Set<Particle> neighbours) {
        super();
        this.p = p;
        this.neighbours = neighbours;
    }

    public double getXForce() {
        double fXSum = 0;
        for (Particle neighbour : neighbours) {
            fXSum += getFN(p, neighbour) * getENX(p, neighbour) + getFT(p, neighbour) * (-(getENY(p, neighbour)));
        }

        return fXSum;
    }

    public double getYForce() {
        double fYSum = p.mass * G;
        for (Particle neighbour : neighbours) {
            fYSum += getFN(p, neighbour) * getENY(p, neighbour) + getFT(p, neighbour) * getENX(p, neighbour);
        }
        return fYSum;
    }

    private double getENY(Particle p1, Particle p2) {
        return (p2.y - p1.y) / getDistance(p2.x, p2.y, p1.x, p1.y);
    }

    private double getENX(Particle p1, Particle p2) {
        return (p2.x - p1.x) / getDistance(p2.x, p2.y, p1.x, p1.y);
    }

    private double getFN(Particle p1, Particle p2) {
        return getEpsilon(p1, p2) * -kn;
    }

    private double getFT(Particle p1, Particle p2) {
        return -kt * getEpsilon(p1, p2) * (((p1.velX - p2.velX) * (-getENY(p1, p2)))
                + ((p1.velY - p2.velY) * (getENX(p1, p2))));
    }

    private double getEpsilon(Particle p1, Particle p2) {
        double ep = p1.radius + p2.radius - (getDistance(p1.x, p1.y, p2.x, p2.y));
        if (ep > p1.radius || ep > p2.radius) {
            ep = Math.min(p1.radius, p2.radius) / 2.0;
        }
        return ep;
    }

    private double getDistance(double x0, double y0, double x1, double y1) {
        return Math.sqrt(Math.pow(x0 - x1, 2) + Math.pow(y0 - y1, 2));
    }

    public void setNeighbours(Set<Particle> neighbours) {
        Set<Particle> filteredNeighbours = new HashSet<>();
        for (Particle neighbour : neighbours) {
            if (getDistance(p.getX(), p.getY(), neighbour.getX(), neighbour.getY()) <= p.getRadius() + neighbour.getRadius()) {
                filteredNeighbours.add(neighbour);
            }
        }
        this.neighbours = filteredNeighbours;
    }

    public Particle getParticle() {
        return p;
    }

    public void reset() {
        neighbours = null;
    }
}
