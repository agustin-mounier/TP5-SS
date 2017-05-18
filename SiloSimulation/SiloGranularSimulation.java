package SiloSimulation;

import CellIndexMethod.CellIndexMethod;
import IntegrationMethod.Beeman;
import IntegrationMethod.GranularForce;
import IntegrationMethod.Particle;
import javafx.util.Pair;

import java.io.IOException;
import java.io.PrintWriter;
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
    // Cant particles that went out vs time it went out
    private static Map<Integer, Double> flowint;

    private static List<Pair<Double, Double>> systemEnergy;

    int id = 1;
    double tf;
    double dt = 0.000001;
    int fpsModule = (int) (1 / dt) / 60;

    int secondsModule = (int) (1 / dt);

    boolean open = true;
    Map<Particle, Beeman> integrationMethod;

    int particlesFlowCount;
    int particlesFlowCountGlobal;

    public SiloGranularSimulation(double L, double W, double D, int maxParticles, double tf) {
        this.L = L;
        this.W = W;
        this.D = D;
        this.tf = tf;
        integrationMethod = new HashMap<>();
        systemEnergy = new ArrayList<>();
        flowint = new HashMap<>();
        locateParticles(maxParticles);
    }

    public void run() {
        double currentTime = 0;

        CellIndexMethod CIM = new CellIndexMethod(L + BELOW_SPACE, D / 5, particles, false);
        int iteration = 0;
        while (currentTime < tf) {
            if (iteration % fpsModule == 0) {
                printMeanParticlesEnergy(currentTime);
                printState(iteration / fpsModule, particles);
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
                if (p.getY() < BELOW_SPACE && !p.isOut()) {
                    p.setOut(true);

                    particlesFlowCount++;
                    particlesFlowCountGlobal++;
                    flowint.put(particlesFlowCountGlobal, currentTime);
                }

                nextGen.add(p);
            }
            if (iteration % secondsModule == 0) {
//                System.out.println("Flow: " + particlesFlowCount + " particles/s");
                particlesFlowCount = 0;
            }
            this.particles = nextGen;
            currentTime += dt;
            iteration++;
        }

//        System.out.println("Global Flow: " + particlesFlowCountGlobal / 5.0 + " particles/s");

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
            } while (isOverlaped(x, y, diameter / 2));
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
            y = BELOW_SPACE + L - (2 * p.getRadius()) + ((tries / 10) * 2 * p.getRadius());
            tries++;
        } while (isOverlaped(x, y, p.getRadius()));
        p.setX(x);
        p.setY(y);
        p.setVelX(0);
        p.setVelY(0);
        p.setOut(false);
        GranularForce gf = new GranularForce(p, null);
        Beeman beeman = new Beeman(gf, dt);
        integrationMethod.put(p, beeman);
    }


    private boolean isOverlaped(double x, double y, double r) {
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
                if (Math.abs(p.getY() - BELOW_SPACE) < p.getRadius()) {
                    if (p.getX() <= (W / 2 - D / 2) || p.getX() >= (W / 2 + D / 2)) {
                        return new Particle(id++, p.getX(), BELOW_SPACE - p.getRadius(), 0, 0, p.getRadius(), p.getMass());
                    } else if (p.getX() - p.getRadius() <= (W / 2 - D / 2) && getDistance(p.getX(), p.getY(),
                            (W / 2 - D / 2), BELOW_SPACE) < p.getRadius()) {
                        return new Particle(id++, (W / 2 - D / 2), BELOW_SPACE, 0, 0, p.getRadius(), p.getMass());
                    } else if (p.getX() + p.getRadius() >= (W / 2 + D / 2) && getDistance(p.getX(), p.getY(),
                            (W / 2 + D / 2), BELOW_SPACE) < p.getRadius()) {
                        return new Particle(id++, (W / 2 + D / 2), BELOW_SPACE, 0, 0, p.getRadius(), p.getMass());
                    }
                }
            } else {
                return new Particle(id++, p.getX(), BELOW_SPACE - p.getRadius(), 0, 0, p.getRadius(), MASS);
            }
        }
        return null;
    }

    public void printState(int iteration, List<Particle> particles) {

        System.out.println(N);
        System.out.println(iteration);
        for (Particle p : particles) {
            System.out.println(p);
        }
    }

    public void printMeanParticlesEnergy(double t) {
        double ret = 0;
        for (Particle particle : particles) {
            ret += (1 / 2.0) * particle.getMass() * (Math.pow(particle.getVelX(), 2) + Math.pow(particle.getVelY(), 2));
        }
        systemEnergy.add(new Pair(t, ret));
    }

    private static void printFlowIntoFile(int n) {
        try{
            PrintWriter writer = new PrintWriter("flow-"+n+".tsv", "UTF-8");
            for (Integer i : flowint.keySet()) {
                writer.println(String.format(Locale.FRENCH, "%.3f",flowint.get(i)) + "\t" + i);
            }
            writer.println("------- System Energy --------");
            for (Pair p : systemEnergy) {
                writer.println(p.getKey() + "\t" + p.getValue());
            }

            writer.close();
        } catch (IOException e) {
            // do something
        }
    }

    public static void main(String[] args) {

        int L = 10;
        int W = 5;
        int D = 2;
        int PARTICLES = 50;
        int SIMULATION_TIME = 1;

        SiloGranularSimulation SGM = new SiloGranularSimulation(L, W, D, PARTICLES, SIMULATION_TIME);

        SGM.run();
        printFlowIntoFile(PARTICLES);
    }
}