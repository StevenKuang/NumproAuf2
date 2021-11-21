import java.util.Arrays;

/**
 * Die Klasse CubicSpline bietet eine Implementierung der kubischen Splines. Sie
 * dient uns zur effizienten Interpolation von aequidistanten Stuetzpunkten.
 *
 * @author braeckle
 */
public class CubicSpline implements InterpolationMethod {

    /**
     * linke und rechte Intervallgrenze x[0] bzw. x[n]
     */
    double a, b;

    /**
     * Anzahl an Intervallen
     */
    int n;

    /**
     * Intervallbreite
     */
    double h;

    /**
     * Stuetzwerte an den aequidistanten Stuetzstellen
     */
    double[] y;

    /**
     * zu berechnende Ableitunge an den Stuetzstellen
     */
    double yprime[];

    /**
     * {@inheritDoc} Zusaetzlich werden die Ableitungen der stueckweisen
     * Polynome an den Stuetzstellen berechnet. Als Randbedingungen setzten wir
     * die Ableitungen an den Stellen x[0] und x[n] = 0.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        this.a = a;
        this.b = b;
        this.n = n;
        h = ((double) b - a) / (n);

        this.y = Arrays.copyOf(y, n + 1);

        /* Randbedingungen setzten */
        yprime = new double[n + 1];
        yprime[0] = 0;
        yprime[n] = 0;

        /* Ableitungen berechnen. Nur noetig, wenn n > 1 */
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * getDerivatives gibt die Ableitungen yprime zurueck
     */
    public double[] getDerivatives() {
        return yprime;
    }

    /**
     * Setzt die Ableitungen an den Raendern x[0] und x[n] neu auf yprime0 bzw.
     * yprimen. Anschliessend werden alle Ableitungen aktualisiert.
     */
    public void setBoundaryConditions(double yprime0, double yprimen) {
        yprime[0] = yprime0;
        yprime[n] = yprimen;
        if (n > 1) {
            computeDerivatives();
        }
    }

    /**
     * Berechnet die Ableitungen der stueckweisen kubischen Polynome an den
     * einzelnen Stuetzstellen. Dazu wird ein lineares System Ax=c mit einer
     * Tridiagonalen Matrix A und der rechten Seite c aufgebaut und geloest.
     * Anschliessend sind die berechneten Ableitungen y1' bis yn-1' in der
     * Membervariable yprime gespeichert.
     * <p>
     * Zum Zeitpunkt des Aufrufs stehen die Randbedingungen in yprime[0] und yprime[n].
     * Speziell bei den "kleinen" Faellen mit Intervallzahlen n = 2
     * oder 3 muss auf die Struktur des Gleichungssystems geachtet werden. Der
     * Fall n = 1 wird hier nicht beachtet, da dann keine weiteren Ableitungen
     * berechnet werden muessen.
     */
    public void computeDerivatives() {
        /* TODO: diese Methode ist zu implementieren */

        if (n == 2) {
            // in this case we only have to look for the yprime[1]
            // according to the formula we don't have an LGS instead we just have one equality because just one unknown
            // yprime[0] + 4*yprime[1] + yprime[2] = ( 3/h * (y[2] - y[0]) - yprime[0]) <- I am not quite sure this works though

            yprime[1] = ((3 / h) * (y[2] - y[0]) -  yprime[0]) / 4;
            return;
        }
        if (n == 3) {
            // in this case we don't have enough dimension to build a Tridiagonalmatrix
            // we have to solve the lgs ourselves (the matrix is a 2x2 matrix)
            // wrote a helper method to solve it with Gauss (I did it first on paper manually)

            //our b vector
            double[] bVector = new double[2];
            bVector[0] = (3 / h) * (y[2] - y[0]) - yprime[0];
            bVector[1] = (3 / h) * (y[3] - y[1]) - yprime[3];
            double[] result = solveWithSpecificMatrix(bVector);
            yprime[1] = result[0];
            yprime[2] = result[1];
            return;
        }
        double[][] dVector = buildMatrix(n-1);
        TridiagonalMatrix triMatrix = new TridiagonalMatrix(dVector[0], dVector[1], dVector[2]);
        double[] bVector = new double[n-1];
        for (int i = 0; i < n-1; i++) {
            bVector[i] = (3/h) * ( y[i+2] - y[i] );
        }
        bVector[0] -= yprime[0];
        bVector[n-2] -= yprime[n];
        double[] result = triMatrix.solveLinearSystem(bVector);
        for (int i = 1; i < n ; i++) {
            yprime[i] = result[i-1];
        }
    }

    /**
     * {@inheritDoc} Liegt z ausserhalb der Stuetzgrenzen, werden die
     * aeussersten Werte y[0] bzw. y[n] zurueckgegeben. Liegt z zwischen den
     * Stuetzstellen x_i und x_i+1, wird z in das Intervall [0,1] transformiert
     * und das entsprechende kubische Hermite-Polynom ausgewertet.
     */
    @Override
    public double evaluate(double z) {
        /* TODO: diese Methode ist zu implementieren */
        if (z <= a) return y[0];
        if (z >= b) return y[n];
        /* Intervall finden, in dem z liegt */
        int interv = (int) ((z - a) / h);   // -> this will find our i ?
        double xi = a + interv * h;
        //double xip1 = a + (interv + 1) * h;
        // transform z to t
        double t = (z - xi) / h;
        double H0,H1,H2,H3;
        double tSquared = Math.sqrt(t);
        double tPower3 = t * Math.sqrt(t);
        H0 = 1 - 3 * tSquared + 2 * tPower3;
        H1 = 3 * tSquared  - 2 * tPower3;
        H2 = t - 2 * tSquared + tPower3;
        H3 = - tSquared + tPower3;

        return y[interv] * H0 + y[interv+1] * H1 + h * yprime[interv] * H2 + h * yprime[interv+1] * H3;

    }

    private double[] solveWithSpecificMatrix(double[] b) {
        double x, y;
        x = (4 * b[0] - b[1]) / 7;
        y = b[0] - 4 * x;
        return new double[]{x, y};
    }

    private double[][] buildMatrix(int dim) {
        // first we build the diagonal : it has the same dimension as our end-Matrix and all of its elements is the number 4
        double[] diag = new double[dim];
        double[] lower = new double[dim - 1];
        double[] upper = new double[dim - 1];
        Arrays.fill(diag, 4);
        Arrays.fill(upper, 1);
        Arrays.fill(lower, 1);
        return new double[][] {upper, diag, lower};
    }
}
