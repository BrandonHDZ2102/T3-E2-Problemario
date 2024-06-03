
package problemariot3;

public class GaussJordan4 {
    public static void main(String[] args) {
        // Definir la matriz de coeficientes y el vector de términos independientes
        double[][] A = {
            {4, -2, 1},
            {1, 6, -1},
            {2, 1, 8}
        };

        double[] b = {7, 13, 17};

        // Resolver el sistema de ecuaciones
        double[] x = gaussJordan(A, b);

        // Imprimir la solución
        System.out.println("Solución:");
        for (int i = 0; i < x.length; i++) {
            System.out.printf("x%d = %.2f\n", i + 1, x[i]);
        }
    }

    public static double[] gaussJordan(double[][] A, double[] b) {
        int n = A.length;

        // Formar la matriz aumentada
        double[][] augmentedMatrix = new double[n][n + 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A[i][j];
            }
            augmentedMatrix[i][n] = b[i];
        }

        // Aplicar el método de Gauss-Jordan
        for (int i = 0; i < n; i++) {
            // Hacer que el pivote sea igual a 1
            double pivot = augmentedMatrix[i][i];
            for (int j = i; j < n + 1; j++) {
                augmentedMatrix[i][j] /= pivot;
            }

            // Hacer ceros debajo del pivote
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor = augmentedMatrix[k][i];
                    for (int j = i; j < n + 1; j++) {
                        augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                    }
                }
            }
        }

        // Extraer la solución del sistema
        double[] x = new double[n];
        for (int i = 0; i < n; i++) {
            x[i] = augmentedMatrix[i][n];
        }

        return x;
    }
}

/*
Entradas:
Matriz de coeficientes A:
4 -2  1
1  6 -1
2  1  8

Vector de términos independientes b:
7
13
17

Salidas:
Solución:
x1 = 2.00
x2 = 1.00
x3 = 1.00
*/
