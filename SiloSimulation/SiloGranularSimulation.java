package SiloSimulation;

import cellindexmethod.CellIndexMethod;
import integrationMethod.Beeman;
import integrationMethod.GranularForce;
import integrationMethod.Particle;

import java.util.*;

/**
 * Created by amounier on 5/13/17.
 */
public class SiloGranularSimulation {

    double L;
    double W;
    double D;

    int N;

    static final double MASS = 0.01;
    static final double MAX_TRIES = 100;

    static final double BELOW_SPACE = 1;

    List<Particle> particles;

    int id = 1;
    double tf;
    double dt = 0.000001;
    int fpsModule = (int) (1 / dt) / 60;

    boolean open = true;
    Map<Particle, Beeman> integrationMethod;

    public SiloGranularSimulation(double L, double W, double D, int maxParticles, double tf) {
        this.L = L;
        this.W = W;
        this.D = D;
        this.tf = tf;
        integrationMethod = new HashMap<>();
        locateParticles(maxParticles);
    }

    public void run() {
        double currentTime = 0;

        CellIndexMethod CIM = new CellIndexMethod(L + BELOW_SPACE, D / 5, particles, false);
        int iteration = 0;
        while (currentTime < tf) {
            if (iteration % fpsModule == 0) {
                printStateForOvito(iteration / fpsModule, particles);
            }
            CIM.reloadMatrix(particles);

            Map<Particle, Set<Particle>> neighbors = CIM.findNeighbors(particles);

            List<Particle> nextGen = new ArrayList<>();

            addWallParticles(neighbors);

            for (Particle particle : particles) {
                Beeman beeman = integrationMethod.get(particle);
                beeman.getForce().setNeighbours(neighbors.get(particle));
                Particle p = beeman.moveParticle();
                if (p.getY() < BELOW_SPACE / 2) {
                    relocateParticle(p);
                }
                nextGen.add(p);
            }

            this.particles = nextGen;
            currentTime += dt;
            iteration++;
        }


    }

    private void locateParticles(int size) {
        this.particles = new ArrayList<>();

        boolean flag = true;

        while (flag && particles.size() < size) {
            double diameter = Math.random() * (this.D / 5 - this.D / 7) + (this.D / 7);
            double r = diameter / 2;
            double x, y;
            int tries = 0;
            do {
                x = Math.random() * (this.W - 2 * r) + r;
                y = BELOW_SPACE + Math.random() * (this.L - 2 * r) + r;
                tries++;
                if (tries == MAX_TRIES) {
                    flag = false;
                    break;
                }
            } while (overlap(x, y, diameter / 2));
            if (flag) {
                Particle p = new Particle(id++, x, y, 0, 0, diameter / 2, MASS);
                GranularForce gf = new GranularForce(p, null);
                Beeman beeman = new Beeman(gf, dt);
                integrationMethod.put(p, beeman);
                this.particles.add(p);
            }
        }
        this.N = this.particles.size();
    }

    private void relocateParticle(Particle p) {
        double x, y;
        int tries = 0;
        do {
            x = Math.random() * (this.W - 2 * p.getRadius()) + p.getRadius();
            y = BELOW_SPACE + L - (2 * p.getRadius()) + ((tries/10) * 2 * p.getRadius());
            tries++;
        } while (overlap(x, y, p.getRadius()));
        p.setX(x);
        p.setY(y);
        p.setVelX(0);
        p.setVelY(0);
        GranularForce gf = new GranularForce(p, null);
        Beeman beeman = new Beeman(gf, dt);
        integrationMethod.put(p, beeman);
    }


    private boolean overlap(double x, double y, double r) {
        for (Particle p : particles) {
            if (getDistance(p.getX(), p.getY(), x, y) < (p.getRadius() + r))
                return true;
        }
        return false;
    }

    private double getDistance(double x0, double y0, double x1, double y1) {
        return Math.sqrt(Math.pow(x0 - x1, 2) + Math.pow(y0 - y1, 2));
    }

    public void addWallParticles(Map<Particle, Set<Particle>> neighbors) {
        for (Particle p : neighbors.keySet()) {
            Particle wallParticle = getWallParticle(p);
            if (wallParticle != null) {
                wallParticle.setWall(true);
                neighbors.get(p).add(wallParticle);
            }
        }
    }

    public Particle getWallParticle(Particle p) {
        if (p.getX() + p.getRadius() > W) {
            return new Particle(id++, W + p.getRadius(), p.getY(), 0, 0, p.getRadius(), MASS);
        } else if (p.getX() - p.getRadius() < 0) {
            return new Particle(id++, -p.getRadius(), p.getY(), 0, 0, p.getRadius(), MASS);
        } else if (p.getY() - p.getRadius() < BELOW_SPACE) {
            if (open) {
                if (Math.abs(p.getY() - (BELOW_SPACE)) < p.getRadius()) {
                    if (p.getX() <= (W / 2 - D / 2) || p.getX() >= (W / 2 + D / 2)) {
                        return new Particle(id++,p.getX(), (BELOW_SPACE) - p.getRadius(), 0,0,p.getRadius(), p.getMass());
                    } else if (p.getX() - p.getRadius() <= (W / 2 - D / 2) && getDistance(p.getX(), p.getY(),
                            (W / 2 - D / 2), (BELOW_SPACE)) < p.getRadius()) {
                         return new Particle(id++,(W / 2 - D / 2), (BELOW_SPACE ),0, 0,p.getRadius(),p.getMass());
                    } else if (p.getX() + p.getRadius() >= (W / 2 + D / 2) && getDistance(p.getX(), p.getY(),
                            (W / 2 + D / 2), (BELOW_SPACE )) < p.getRadius()) {
                        return new Particle(id++, (W / 2 + D / 2),(BELOW_SPACE ),0, 0,p.getRadius(),p.getMass());
                    }
                }
            } else {
                return new Particle(id++, p.getX(), BELOW_SPACE - p.getRadius(), 0, 0, p.getRadius(), MASS);
            }
        }
        return null;
    }

    public void printStateForOvito(int iteration, List<Particle> particles) {

        System.out.println(N);
        System.out.println(iteration);
        for (Particle p : particles) {
            System.out.println(p);
        }
    }

    public static void main(String[] args) {
        SiloGranularSimulation SGM = new SiloGranularSimulation(5, 2, 1, 100, 5);

        SGM.run();
    }
}
