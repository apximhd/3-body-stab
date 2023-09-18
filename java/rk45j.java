/*
 * rk45j.java
 * 
 * Copyright 2016, 2023 Alexey V. Pasechnik <Alexey.Pasechnik@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */

import java.lang.Math;
import java.io.*;
import java.util.*;

public class rk45j {

    //In the next string replace "gnuplot" to the full path of the gnuplot program in your system
    public static final String[] gnuplotCommand = {"gnuplot", "taei-j.plt"};

    // Common constants:
    public static final double pi = 4 * Math.atan(1); // 3.1415926...
    public static final double G = 6.67e-8;             // Newtonian gravitational constant, cm^3 g^-1 s^-2
    public static final double Msol = 1.99e33;          // Solar mass, g
    public static final double au = 1.496e13;           // Astronomy unit, cm
    public static final double kms = 1.0e5;             // cm/s per km/s
    public static final double cms = 2978000.0;         // Earth orbital velosity in cm/s
    public static final double deg = pi / 180.0;        // radians per degree 
    public static final double day = 3600 * 24;         // seconds per day
    public static final double year = 365.2422 * day;   // seconds per year

    // Initial conditions:
    public static final double Q = 7;                 // Q = a_3/a_12*(1-e_3)
    public static final double iAB = 0 * deg;           // A-B orbit inclination
    public static final double iAC = 140 * deg;         // C orbit inclination
    public static final double MassA = 1.5 * Msol;      // Mass of A-component
    public static final double MassB = 0.5 * Msol;      // Mass of B-component
    public static final double MassC = 0.5 * Msol;      // Mass of C-component
    public static final double eAB = 0.6;               // A-B orbit eccentricity
    public static final double eAC = 0.5;               // C orbit eccentricity
    public static final double a12 = 1 * au;            // A-B semimajor axis
    public static final double a3 = Q * a12 / (1 - eAC);// C orbit major semiaxe
    public static final double OmegaAB = 0 * deg;       // A-B orbit's node line longitude
    public static final double OmegaAC = 0 * deg;       // C orbit's node line longitude
    public static final double PeriArgAB = 60 * deg;    // AB orbit pericenter argument
    public static final double PeriArgAC = 155 * deg;   // C orbit pericenter argument
    public static final double AnomAB = 0 * deg;        // Mean anomaly at epoch t=0, component A
    public static final double AnomAC = 0 * deg;        // Mean anomaly at epoch t=0, component C

    public static final double n_out = 100000; // Number of periods outer binary;
    public static final double Nouts = 1000; // Number of string in output;
    public static final double criteria = 10.0 * a3;
    // public static final double precision = 1.0 * kms;
    public static final double precision = 1.0 * kms;

    // Osculating orbital elements 
    public static class Elements {
        double a12;
        double e12;
        double cosi12;
        double a3;
        double e3;
        double cosi3;
        double cosimut;
    }


    public static void main (String args[])
    {
        Locale.setDefault(new Locale.Builder().setLanguage("en").setRegion("US").build());

        double Tmax = n_out * Math.sqrt(4.0 * (pi*pi) * (a3*a3*a3) / (G * (MassA + MassB + MassC)));
        double T, eps, dt = 1.0; // integration step
        double [] x1 = new double [3];
        double [] v1 = new double [3];
        double [] x2 = new double [3];
        double [] v2 = new double [3];
        double [] x3 = new double [3];
        double [] v3 = new double [3];
        double [] y =  new double [18];
        Elements elements = new Elements();

        long clock = System.currentTimeMillis();
        int every = (int)(Tmax/Nouts);

        KeplerToCartesian ( iAB, iAC, eAB, eAC, PeriArgAB, PeriArgAC,
                            OmegaAB, OmegaAC, AnomAB, AnomAC,
                            a12, a3, MassA*G, MassB*G, MassC*G,
                            x1, v1, x2, v2, x3, v3);

        System.out.format ("Mass_A = %.1f Msol,   ", MassA/Msol);
        System.out.format ("Mass_B = %.1f Msol,   ", MassB/Msol);
        System.out.format ("Mass_C = %.1f Msol\n", MassC/Msol);
        System.out.format ("Inclination outer = %.2f deg\n", iAC/deg);
        System.out.format ("Q = %.1f,  a_12 = %.1f a.u.,  a_3 = %.1f a.u.\n", Q, a12/au, a3/au);
        System.out.format ("Eccentricity inner = %.2f,  Eccentricity outer = %.2f\n", eAB, eAC);
        System.out.format ("T inner = %.1f years\n", (Math.sqrt(4.0*(pi*pi)*(a12*a12*a12)/(G*(MassA+MassB)))/year));
        System.out.format ("T outer = %.1f years\n", (Math.sqrt(4.0*(pi*pi)*(a3*a3*a3)/(G*(MassA+MassB+MassC)))/year));
        System.out.format ("T max = %.1f years = %d periods of outer\n",
            Tmax/year,
            (int)(Tmax/(Math.sqrt(4.0*(pi*pi)*(a3*a3*a3)/(G*(MassA+MassB+MassC))))));
                            
        for (int i = 0; i < 3; i++)
        {
            y[i] = x1[i];
            y[i + 3] = x2[i];
            y[i + 6] = x3[i];
            y[i + 9] = v1[i];
            y[i + 12] = v2[i];
            y[i + 15] = v3[i];
        }

        PrintWriter XYZV, TAEI;
        try {    
            XYZV = new PrintWriter("XYZV_RK45_J");
            TAEI = new PrintWriter("TAEI_RK45_J");
            double step=0.0;
            double steps=0.0;
            int stable=1;
            for (T = 0; T<=Tmax; T+=dt) {
                eps = RK45(y, dt, MassA*G, MassB*G, MassC*G);
                dt = precision/eps;
                steps+=1.0;
                if (T > every*step) {
                    getElements (MassA/Msol, MassB/Msol, MassC/Msol, y, elements);
                    if(((int)(10000.0*T/Tmax))%10 == 0) {
                        System.out.format ("Complete %5.1f%%,  eps = %8.3f,  dt = %8.4f years,  step: %8.3f\u00D710\u2076\r", 100.0*T/Tmax, eps, dt/year, steps*0.000001);
                    }
                    if (elements.e12 > 1.0 ||
                        elements.e3 > 1.0 ||
                        elements.a12 < 0.0 ||
                        elements.a3 < 0.0 ||
                        elements.a12 > elements.a3 ||
                        elements.a3 > criteria) {
                        System.out.format ("\n\nSystem brakes after %6.0f revolutions\n", T/(Math.sqrt(4.0*(pi*pi)*(a3*a3*a3)/(G*(MassA+MassB+MassC)))));
                        stable=0;
                        break;
                    }
                    for (int i=0;i<9;i++)   XYZV.format ("%12.5e ", y[i]/au);
                    for (int i=9;i<18;i++)  XYZV.format ("%12.5e ", y[i]/kms);
                    XYZV.format ("%12.5e\n", eps);

                    TAEI.format ("%12.5e", T/(Math.sqrt(4.0*(pi*pi)*(a3*a3*a3)/G/(MassA+MassB+MassC))));
                    TAEI.format ("%12.5e", T/year);
                    TAEI.format ("%10.3f", elements.a12);
                    TAEI.format ("%10.3f", elements.e12);
                    TAEI.format ("%10.3f", elements.cosi12);
                    TAEI.format ("%10.3f", elements.a3);
                    TAEI.format ("%10.3f", elements.e3);
                    TAEI.format ("%10.3f", elements.cosi3);
                    TAEI.format ("%10.3f", elements.cosimut);
                    TAEI.format (" %d %d %d %d %d\n", 0,0,0,0,0); // Ending zeros for compatipility with G3R output

                    step+=1;
                }
            }
            clock = System.currentTimeMillis() - clock;
            if(stable!=0) System.out.format ("\n\nSystem is stable during %6.0f revolutions\n",
                                T/(Math.sqrt(4.0*(pi*pi)*(a3*a3*a3)/(G*(MassA+MassB+MassC)))));
            System.out.format ("t = %5.2f seconds.%n",clock/1000.);
            XYZV.close();
            TAEI.close();

            try {
                Runtime.getRuntime().exec(gnuplotCommand);
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }


    static double abs(double x) {
        if (x>=0) {
            return x;
        } else {
            return -x;
        }
    }


    static double dmax (double array[], int n) {
        double element;
        element = 0.0;
        for (int i = 0; i < n; i++)
            if (element < abs(array[i])) element = abs(array[i]);
        return element;
    }


    static double KeplerEquation(double A, double e) {
        double E = A;
        for (int i = 0; i < 20; i++) E = e * Math.sin(E) + A;
        return E;
    }

    static void KeplerToCartesian ( double iAB,
                                    double iAC,
                                    double eAB, double eAC,
                                    double PeriArgAB, double PeriArgAC,
                                    double OmegaAB, double OmegaAC,
                                    double AnomAB, double AnomAC,
                                    double a12, double a3,
                                    double M1, double M2, double M3,
                                    double[] x1, double[] v1,
                                    double[] x2, double[] v2,
                                    double[] x3, double[] v3) {

        double M12 = M1 + M2;
        double M123 = M1 + M2 + M3;
        double ecc12 = eAB;        // Eccentricity of inner orbit
        double ecc3 = eAC;         // Eccentricity of outer orbits

        double A1 = a12*(Math.cos(PeriArgAB)*Math.cos(OmegaAB) 
                        -Math.sin(PeriArgAB)*Math.sin(OmegaAB)*Math.cos(iAB));
        double A2 = a12*(Math.cos(PeriArgAB)*Math.sin(OmegaAB) 
                        +Math.sin(PeriArgAB)*Math.cos(OmegaAB)*Math.cos(iAB));
        double A3 = a12*Math.sin(PeriArgAB)*Math.sin(iAB);
        double B1 = a12*Math.sqrt(1-(ecc12*ecc12))*
                      (-Math.sin(PeriArgAB)*Math.cos(OmegaAB) 
                       -Math.cos(PeriArgAB)*Math.sin(OmegaAB)*Math.cos(iAB));
        double B2 = a12*Math.sqrt(1-(ecc12*ecc12))*
                      (-Math.sin(PeriArgAB)*Math.sin(OmegaAB) 
                       +Math.cos(PeriArgAB)*Math.cos(OmegaAB)*Math.cos(iAB));
        double B3 = a12*Math.sqrt(1-(ecc12*ecc12))*Math.cos(PeriArgAB)*Math.sin(iAB);
        double AC = Math.cos(KeplerEquation(AnomAB, ecc12)) - ecc12;
        double BC = Math.sin(KeplerEquation(AnomAB, ecc12));
        double AD = -Math.sin(KeplerEquation(AnomAB, ecc12))/(Math.pow(a12,1.5)*
                    (1 - ecc12*Math.cos(KeplerEquation(AnomAB, ecc12))));
        double BD = Math.cos(KeplerEquation(AnomAB, ecc12))/(Math.pow(a12,1.5)*
                    (1 - ecc12*Math.cos(KeplerEquation(AnomAB, ecc12))));
        double MW1 = -M2/M12;
        double MW2 = M1/M12;

        x1[0] = MW1*(A1*AC + B1*BC);
        x1[1] = MW1*(A2*AC + B2*BC);
        x1[2] = MW1*(A3*AC + B3*BC);
        x2[0] = MW2*(A1*AC + B1*BC);
        x2[1] = MW2*(A2*AC + B2*BC);
        x2[2] = MW2*(A3*AC + B3*BC);
        v1[0] = MW1*(A1*AD + B1*BD)*Math.sqrt(M12);
        v1[1] = MW1*(A2*AD + B2*BD)*Math.sqrt(M12);
        v1[2] = MW1*(A3*AD + B3*BD)*Math.sqrt(M12);
        v2[0] = MW2*(A1*AD + B1*BD)*Math.sqrt(M12);
        v2[1] = MW2*(A2*AD + B2*BD)*Math.sqrt(M12);
        v2[2] = MW2*(A3*AD + B3*BD)*Math.sqrt(M12);

        double A31 = a3*(Math.cos(PeriArgAC)*Math.cos(OmegaAC)
                        -Math.sin(PeriArgAC)*Math.sin(OmegaAC)*Math.cos(iAC));
        double A32 = a3*(Math.cos(PeriArgAC)*Math.sin(OmegaAC)
                        +Math.sin(PeriArgAC)*Math.cos(OmegaAC)*Math.cos(iAC));
        double A33 = a3*Math.sin(PeriArgAC)*Math.sin(iAC);
        double B31 = a3*Math.sqrt(1-(ecc3*ecc3))*
                      (-Math.sin(PeriArgAC)*Math.cos(OmegaAC)
                       -Math.cos(PeriArgAC)*Math.sin(OmegaAC)*Math.cos(iAC));
        double B32 = a3*Math.sqrt(1-(ecc3*ecc3))*
                      (-Math.sin(PeriArgAC)*Math.sin(OmegaAC)
                       +Math.cos(PeriArgAC)*Math.cos(OmegaAC)*Math.cos(iAC));
        double B33 = a3*Math.sqrt(1-(ecc3*ecc3))*Math.cos(PeriArgAC)*Math.sin(iAC);
        double A3C = Math.cos(KeplerEquation(AnomAC, ecc3)) - ecc3;
        double B3C = Math.sin(KeplerEquation(AnomAC, ecc3));
        double A3D = -Math.sin(KeplerEquation(AnomAC, ecc3))/(Math.pow(a3,1.5)*
                    (1 - ecc3*Math.cos(KeplerEquation(AnomAC, ecc3))));
        double B3D = Math.cos(KeplerEquation(AnomAC, ecc3))/(Math.pow(a3,1.5)*
                    (1 - ecc3*Math.cos(KeplerEquation(AnomAC, ecc3))));
 
        x3[0] = A31*A3C + B31*B3C;
        x3[1] = A32*A3C + B32*B3C;
        x3[2] = A33*A3C + B33*B3C;
        v3[0] = (A31*A3D + B31*B3D)*Math.sqrt(M123);
        v3[1] = (A32*A3D + B32*B3D)*Math.sqrt(M123);
        v3[2] = (A33*A3D + B33*B3D)*Math.sqrt(M123);

        for (int i=0;i<3; i++) {
            double vcm = (v1[i]*M1 + v2[i]*M2 + v3[i]*M3)/(M1 + M2 + M3);
            v1[i] -= vcm;
            v2[i] -= vcm;
            v3[i] -= vcm;
        }
    }


    static double gf(double[] y, double[] f, double M1, double M2, double M3, boolean flag) {
        double [] a= new double [9];
        double r12 = Math.sqrt( (y[3]-y[0])*(y[3]-y[0]) + 
                                (y[4]-y[1])*(y[4]-y[1]) + 
                                (y[5]-y[2])*(y[5]-y[2]) );
        r12=1.0/(r12*r12*r12);
        double r13 = Math.sqrt( (y[6]-y[0])*(y[6]-y[0]) +
                                (y[7]-y[1])*(y[7]-y[1]) +
                                (y[8]-y[2])*(y[8]-y[2]) );
        r13=1.0/(r13*r13*r13);
        double r23 = Math.sqrt( (y[3]-y[6])*(y[3]-y[6]) +
                                (y[4]-y[7])*(y[4]-y[7]) +
                                (y[5]-y[8])*(y[5]-y[8]) );
        r23=1.0/(r23*r23*r23);
        for (int i=0; i<3; i++) {
            a[i]   = (M2*(y[i+3]-y[i])*r12 + M3*(y[i+6]-y[i])*r13);
            a[i+3] = (M1*(y[i]-y[i+3])*r12 + M3*(y[i+6]-y[i+3])*r23);
            a[i+6] = (M1*(y[i]-y[i+6])*r13 + M2*(y[i+3]-y[i+6])*r23);
        }
        for (int i=0; i<9; i++) {
            f[i] = y[i + 9]; // f[0]..f[8] - velocities
            f[i+9] = a[i];   // f[9]..f[17] - accelerations
        }
        if(flag) {
            return dmax(a,9);
        } else {
            return 0;
        }
    }


    static double RK45 (double [] y, double dt, double M1, double M2, double M3) {
        double[] k1 = new double[18];
        double[] k2 = new double[18];
        double[] k3 = new double[18];
        double[] k4 = new double[18];
        double[] k5 = new double[18];
        double[] y1 = new double[18];
        double[] f =  new double[18];
        double eps;
        double c21=1.0/4.0;
        double c31=3.0/32.0, c32=9.0/32.0;
        double c41=1932.0/2197.0, c42=-7200.0/2197.0, c43=7296.0/2197.0;
        double c51=439.0/216.0, c52=-8.0, c53=3680.0/513.0, c54=-845.0/4104.0;
        double a1=25.0/216.0, a3=1408.0/2565.0, a4=2197.0/4104.0, a5=-1.0/5.0;
  
        eps = gf (y, f, M1, M2, M3, false);
        for (int i=0; i<18; i++)
        {
            k1[i] = dt*f[i];
            y1[i] = y[i] + c21*k1[i];
        }
        eps = gf (y1, f, M1, M2, M3, false);
        for (int i=0; i<18; i++)
        {
            k2[i] = dt*f[i];
            y1[i] = y[i] + c31*k1[i]+c32*k2[i];
        }
        eps = gf (y1, f, M1, M2, M3, false);
        for (int i=0; i<18; i++)
        {
            k3[i] = dt*f[i];
            y1[i] = y[i] + c41*k1[i]+c42*k2[i]+c43*k3[i];
        }
        eps = gf (y1,f, M1, M2, M3, false);  
        for (int i=0; i<18; i++)
        {
            k4[i] = dt*f[i];
            y1[i] = y[i] + c51*k1[i]+c52*k2[i]+c53*k3[i]+c54*k4[i];
        }
        eps = gf (y1,f, M1, M2, M3, true);
        for (int i=0; i<18; i++)
        {
            k5[i] = dt*f[i];
            y[i] += a1*k1[i]+a3*k3[i]+a4*k4[i]+a5*k5[i];
        }
        return eps;
    }


    static void getElements (double m1, double m2, double m3, double [] xyzv, Elements elements)
    {
        double x1 = xyzv[0] / au;
        double y1 = xyzv[1] / au;
        double z1 = xyzv[2] / au;
        double x2 = xyzv[3] / au;
        double y2 = xyzv[4] / au;
        double z2 = xyzv[5] / au;
        double x3 = xyzv[6] / au;
        double y3 = xyzv[7] / au;
        double z3 = xyzv[8] / au;
        double vx1 = xyzv[9] / cms;
        double vy1 = xyzv[10] / cms;
        double vz1 = xyzv[11] / cms;
        double vx2 = xyzv[12] / cms;
        double vy2 = xyzv[13] / cms;
        double vz2 = xyzv[14] / cms;
        double vx3 = xyzv[15] / cms;
        double vy3 = xyzv[16] / cms;
        double vz3 = xyzv[17] / cms;
        
        double x = (m1 + m2)/m2*(x1 - (x1*m1 + x2*m2)/(m1 + m2));
        double y = (m1 + m2)/m2*(y1 - (y1*m1 + y2*m2)/(m1 + m2));
        double z = (m1 + m2)/m2*(z1 - (z1*m1 + z2*m2)/(m1 + m2));
        double vx = (m1 + m2)/m2*(vx1 - (vx1*m1 + vx2*m2)/(m1 + m2));
        double vy = (m1 + m2)/m2*(vy1 - (vy1*m1 + vy2*m2)/(m1 + m2));
        double vz = (m1 + m2)/m2*(vz1 - (vz1*m1 + vz2*m2)/(m1 + m2));
    
        double mu = Math.sqrt(m1 + m2);
        double wx = vx/mu;
        double wy = vy/mu;
        double wz = vz/mu;
        double r = Math.sqrt(x*x + y*y + z*z);
        double w2 = wx*wx + wy*wy + wz*wz;
        double Alpha = 2/r - w2;
        double jx = y*wz - z*wy;
        double jy = z*wx - x*wz;
        double jz = x*wy - y*wx;
        double d = Math.sqrt(jx*jx + jy*jy + jz*jz);
        double ex = wy*jz - wz*jy - x/r;
        double ey = wz*jx - wx*jz - y/r;
        double ez = wx*jy - wy*jx - z/r;
        elements.e12 = Math.sqrt(ex*ex + ey*ey + ez*ez); // inner binaty eccentricity
        elements.cosi12 = jz/d;                          // cosine of inner inclination
        elements.a12 = 1/Alpha;                          // inner binary semimajor axis

        x = (m1 + m2 + m3)/(m1 + m2)*(x3 - (x1*m1 + x2*m2 + x3*m3)/(m1 + m2 + m3));
        y = (m1 + m2 + m3)/(m1 + m2)*(y3 - (y1*m1 + y2*m2 + y3*m3)/(m1 + m2 + m3));
        z = (m1 + m2 + m3)/(m1 + m2)*(z3 - (z1*m1 + z2*m2 + z3*m3)/(m1 + m2 + m3));
        vx = (m1 + m2 + m3)/(m1 + m2)*(vx3 - (vx1*m1 + vx2*m2 + vx3*m3)/(m1 + m2 + m3));
        vy = (m1 + m2 + m3)/(m1 + m2)*(vy3 - (vy1*m1 + vy2*m2 + vy3*m3)/(m1 + m2 + m3));
        vz = (m1 + m2 + m3)/(m1 + m2)*(vz3 - (vz1*m1 + vz2*m2 + vz3*m3)/(m1 + m2 + m3));

        mu = Math.sqrt(m1 + m2 + m3);
        wx = vx/mu;
        wy = vy/mu;
        wz = vz/mu;
        r = Math.sqrt(x*x + y*y + z*z);
        w2 = wx*wx + wy*wy + wz*wz;
        Alpha = 2/r - w2;
        double jx3 = y*wz - z*wy;
        double jy3 = z*wx - x*wz;
        double jz3 = x*wy - y*wx;
        double d3 = Math.sqrt(jx3*jx3 + jy3*jy3 + jz3*jz3);
        ex = wy*jz3 - wz*jy3 - x/r;
        ey = wz*jx3 - wx*jz3 - y/r;
        ez = wx*jy3 - wy*jx3 - z/r;
        elements.e3 = Math.sqrt(ex*ex + ey*ey + ez*ez);       // outer eccentricity
        elements.cosi3 = jz3/d3;                              // cosine of outer inclination
        elements.a3 = 1/Alpha;                                // outer semimajor axis
        elements.cosimut = (jx*jx3 + jy*jy3 + jz*jz3)/(d*d3); // cosine of mutual inclination
    }
}
