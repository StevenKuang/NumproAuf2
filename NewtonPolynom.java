import java.util.Arrays;

/**
 * Die Klasse Newton-Polynom beschreibt die Newton-Interpolation. Die Klasse
 * bietet Methoden zur Erstellung und Auswertung eines Newton-Polynoms, welches
 * uebergebene Stuetzpunkte interpoliert.
 *
 * @author braeckle
 *
 */
public class NewtonPolynom implements InterpolationMethod {

    /** Stuetzstellen xi */
    double[] x;

    /**
     * Koeffizienten/Gewichte des Newton Polynoms p(x) = a0 + a1*(x-x0) +
     * a2*(x-x0)*(x-x1)+...
     */
    double[] a;

    /**
     * die Diagonalen des Dreiecksschemas. Diese dividierten Differenzen werden
     * fuer die Erweiterung der Stuetzstellen benoetigt.
     */
    double[] f;

    /**
     * leerer Konstruktore
     */
    public NewtonPolynom() {
    };

    /**
     * Konstruktor
     *
     * @param x
     *            Stuetzstellen
     * @param y
     *            Stuetzwerte
     */
    public NewtonPolynom(double[] x, double[] y) {
        this.init(x, y);
    }

    /**
     * {@inheritDoc} Zusaetzlich werden die Koeffizienten fuer das
     * Newton-Polynom berechnet.
     */
    @Override
    public void init(double a, double b, int n, double[] y) {
        x = new double[n + 1];
        double h = (b - a) / n;

        for (int i = 0; i < n + 1; i++) {
            x[i] = a + i * h;
        }
        computeCoefficients(y);
    }

    /**
     * Initialisierung der Newtoninterpolation mit beliebigen Stuetzstellen. Die
     * Faelle "x und y sind unterschiedlich lang" oder "eines der beiden Arrays
     * ist leer" werden nicht beachtet.
     *
     * @param x
     *            Stuetzstellen
     * @param y
     *            Stuetzwerte
     */
    public void init(double[] x, double[] y) {
        this.x = Arrays.copyOf(x, x.length);
        computeCoefficients(y);
    }

    /**
     * computeCoefficients belegt die Membervariablen a und f. Sie berechnet zu
     * uebergebenen Stuetzwerten y, mit Hilfe des Dreiecksschemas der
     * Newtoninterpolation, die Koeffizienten a_i des Newton-Polynoms. Die
     * Berechnung des Dreiecksschemas soll dabei lokal in nur einem Array der
     * Laenge n erfolgen (z.B. spaltenweise Berechnung). Am Ende steht die
     * Diagonale des Dreiecksschemas in der Membervariable f, also f[0],f[1],
     * ...,f[n] = [x0...x_n]f,[x1...x_n]f,...,[x_n]f. Diese koennen spaeter bei
     * der Erweiterung der Stuetzstellen verwendet werden.
     *
     * Es gilt immer: x und y sind gleich lang.
     */
    private void computeCoefficients(double[] y) {
        int length = y.length;
        int count=0;
        this.a=Arrays.copyOf(y, y.length);
        f = new double[length];
        f[count++] = y[length-1];

        for(int j = 1 ; j < length ; j++ ) {
            for(int i = length-1; i > j-1; i--) {
                a[i] = (a[i] - a[i-1]) / (x[i]-x[i-j]);
            }
            f[count++] = a[length-1];
        }
        /*
        double dreieck[][] = new double [length][length];
        for(int i = 0 ; i < length ; i++) dreieck[0][i]=y[i];
        for(int i = 1 ; i < length ; i++) {
            for(int j = 0 ; j < length - i ; j++) {
                dreieck[i][j] = (dreieck[i-1][j+1] - dreieck[i-1][j])/(x[j+i] - x[j]);
            }
        }
        for(int i = 0; i < length ;i++) {
            a[i] = dreieck[i][0];
        }

         for(int i = 0,j = length - 1; i < length  ;i++,j--) {
            f[i] = dreieck[i][j];
        }
        //System.out.println(Arrays.deepToString(dreieck));
        //System.out.println(Arrays.toString(a));

        /* TODO: diese Methode ist zu implementieren */
        System.out.println(Arrays.toString(f));
    }

    /**
     * Gibt die Koeffizienten des Newton-Polynoms a zurueck
     */
    public double[] getCoefficients() {
        return a;
    }

    /**
     * Gibt die Dividierten Differenzen der Diagonalen des Dreiecksschemas f
     * zurueck
     */
    public double[] getDividedDifferences() {
        return f;
    }

    /**
     * addSamplintPoint fuegt einen weiteren Stuetzpunkt (x_new, y_new) zu x
     * hinzu. Daher werden die Membervariablen x, a und f vergoessert und
     * aktualisiert . Das gesamte Dreiecksschema muss dazu nicht neu aufgebaut
     * werden, da man den neuen Punkt unten anhaengen und das alte
     * Dreiecksschema erweitern kann. Fuer diese Erweiterungen ist nur die
     * Kenntnis der Stuetzstellen und der Diagonalen des Schemas, bzw. der
     * Koeffizienten noetig. Ist x_new schon als Stuetzstelle vorhanden, werden
     * die Stuetzstellen nicht erweitert.
     *
     * @param x_new
     *            neue Stuetzstelle
     * @param y_new
     *            neuer Stuetzwert
     */
    public void addSamplingPoint(double x_new, double y_new) {
        this.x=Arrays.copyOf(x, x.length+1);
        this.a=Arrays.copyOf(a, a.length+1);
        this.f=Arrays.copyOf(f, f.length+1);
        x[x.length-1] = x_new;
        double temp1 = y_new;
        double temp;
        int count=0;
        for(int i = 0; i < f.length ; i+=2) {
            temp = (temp1 - f[i]) / (x_new - x[x.length-2-i]);
            f[i]= temp1;
            count++;
            if(count==f.length-1) {
                f[count] = temp; break;
            }
            temp1 = (temp-f[i+1]) / (x_new - x[x.length-3-i]);
            f[i+1] = temp;
            count++;
            if(count==f.length-1) {
                f[count] = temp1; break;
            }
        }
        a[a.length-1]=f[f.length-1];
       // System.out.println(Arrays.toString(f));
        //System.out.println(Arrays.toString(a));
        /* TODO: diese Methode ist zu implementieren */
    }

    /**
     * {@inheritDoc} Das Newton-Polynom soll effizient mit einer Vorgehensweise
     * aehnlich dem Horner-Schema ausgewertet werden. Es wird davon ausgegangen,
     * dass die Stuetzstellen nicht leer sind.
     */
    @Override
    public double evaluate(double z) {
        /* TODO: diese Methode ist zu implementieren */
        double result=a[a.length-1];

        for(int i=a.length-2;i>=0;i--){
            result=a[i]+(z-x[i])*result;
        }
        return result;
    }
}

