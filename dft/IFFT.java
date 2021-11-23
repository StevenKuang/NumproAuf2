package dft;

import java.util.Arrays;

/**
 * Schnelle inverse Fourier-Transformation
 *
 * @author Sebastian Rettenberger
 */
public class IFFT {
    /**
     * Schnelle inverse Fourier-Transformation (IFFT).
     *
     * Die Funktion nimmt an, dass die Laenge des Arrays c immer eine
     * Zweierpotenz ist. Es gilt also: c.length == 2^m fuer ein beliebiges m.
     */
    public static double power(double a, int n) {
        if (n == 0) 
            return 1;
        return a*power(a,n-1);
    }
    public static Complex[] ifft(Complex[] c) {
        // TODO: diese Methode ist zu implementieren
        int n = c.length;
        Complex[] v = new Complex[n];
        Complex omega;
        if (n == 1) {
            v[0] = c[0];
        } else {
            int m = n/2;
            Complex[] z1 = new Complex[m];
            Complex[] z2 = new Complex[m];
            for (int i = 0; i < m; i++) {
                z1[i] = c[2*i];
                z2[i] = c[2*i+1];  
            }
            z1 = ifft(z1);
            z2 = ifft(z2);
            omega = new Complex(Math.cos(2*Math.PI/n), Math.sin(2*Math.PI/n));
            for (int j = 0; j < m; j++) {
                // omega = Complex.fromPolar(1, j * 2 * Math.PI / n);   // no conj
                v[j] = z1[j].add((omega.power(j)).mul(z2[j]));
                v[m+j] = z1[j].sub((omega.power(j)).mul(z2[j]));
            }
        }
        return v;
    }
}
