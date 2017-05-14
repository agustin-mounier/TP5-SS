package integrationMethod;

import java.util.DoubleSummaryStatistics;
import java.util.HashSet;
import java.util.Set;

/**
 * Created by amounier on 5/13/17.
 */
public class GranularForce {

    Particle p;
    double kn = 1E5;
    double kt = 2 * kn;
    private static final double G = -10;
    private Set<Particle> neighbours;

    public GranularForce(Particle p, Set<Particle> neighbours) {
        super();
        this.p = p;
        this.neighbours = neighbours;
    }

    public double getXForce() {
        double fX = 0;
        for (Particle particle : neighbours) {
            fX += getFN(p, particle) * getENX(p, particle) + getFT(p, particle) * (-(getENY(p, particle)));
        }

        return fX;
    }

    public double getYForce() {
        double fY = p.mass * G;
        for (Particle particle : neighbours) {
            fY += getFN(p, particle) * getENY(p, particle) + getFT(p, particle) * getENX(p, particle);
        }
        return fY;
    }

    private double getENY(Particle p, Particle particle) {
        return (particle.y - p.y) / getDistance(particle.x, particle.y, p.x, p.y);
    }

    private double getENX(Particle p, Particle particle) {
        return (particle.x - p.x) / getDistance(particle.x, particle.y, p.x, p.y);
    }

    private double getFN(Particle p, Particle other) {
        return -kn * getEpsilon(p, other);
    }

    private double getFT(Particle p, Particle other) {
        return -kt * getEpsilon(p, other) * (((p.velX - other.velX) * (-getENY(p, other)))
                + ((p.velY - other.velY) * (getENX(p, other))));
    }

    private double getEpsilon(Particle p, Particle other) {
        double ep = p.radius + other.radius - (getDistance(p.x, p.y, other.x, other.y));
        if (ep > p.radius || ep > other.radius) {
            ep = Math.min(p.radius, other.radius) / 2.0;
        }
        return ep;
    }

    private double getDistance(double x0, double y0, double x1, double y1) {
        return Math.sqrt(Math.pow(x0 - x1, 2) + Math.pow(y0 - y1, 2));
    }

    public void setNeighbours(Set<Particle> neighbours) {
        Set<Particle> filteredNeighbours = new HashSet<>();
        for(Particle neighbour : neighbours) {
            if(getDistance(p.getX(), p.getY(), neighbour.getX(), neighbour.getY()) <= p.getRadius() + neighbour.getRadius()) {
                filteredNeighbours.add(neighbour);
            }
        }
        this.neighbours = filteredNeighbours;
    }

    public Set<Particle> getNeighbours() {
        return neighbours;
    }

    public Particle getParticle() {
        return p;
    }

    public void reset() {
        neighbours = null;
    }
}
