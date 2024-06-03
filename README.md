# T6----E2----Problemario-Metodos

### Equipo

- Brandon Hernández Espinosa
- Italia Yoselin Lozada Olvera

  ## Eliminación Gaussiana.
### Descripción: 

La eliminación gaussiana es un algoritmo utilizado para resolver sistemas de ecuaciones lineales. El método consiste en aplicar operaciones elementales sobre las filas de una matriz aumentada para transformarla en una matriz triangular superior, de la cual se pueden resolver las variables mediante sustitución hacia atrás.

### Pseudocódigo:
```
Algoritmo EliminacionGaussiana(A, b)
    Entrada: 
        A: Matriz de coeficientes de tamaño n x n
        b: Vector de términos independientes de tamaño n
    Salida: 
        x: Vector solución de tamaño n

    AugmentedMatrix <- [A|b]  # Construir la matriz aumentada
    n <- número de filas de A

    Para k desde 1 hasta n hacer:
        # Encontrar el pivote
        Pivote <- AugmentedMatrix[k][k]
        Si Pivote = 0 entonces
            Intercambiar la fila k con una fila i donde AugmentedMatrix[i][k] != 0
            Pivote <- AugmentedMatrix[k][k]

        # Hacer cero los elementos debajo del pivote
        Para i desde k+1 hasta n hacer:
            Factor <- AugmentedMatrix[i][k] / Pivote
            Para j desde k hasta n+1 hacer:
                AugmentedMatrix[i][j] <- AugmentedMatrix[i][j] - Factor * AugmentedMatrix[k][j]

    # Sustitución hacia atrás
    x[n] <- AugmentedMatrix[n][n+1] / AugmentedMatrix[n][n]
    Para i desde n-1 hasta 1 hacer:
        Suma <- 0
        Para j desde i+1 hasta n hacer:
            Suma <- Suma + AugmentedMatrix[i][j] * x[j]
        x[i] <- (AugmentedMatrix[i][n+1] - Suma) / AugmentedMatrix[i][i]

    Retornar x
```

### Implementacion 
- Implementacion en C#

```
using System;

class Program
{
    static void Main(string[] args)
    {
        double[,] A = {
            { 2, 1, -1 },
            { -3, -1, 2 },
            { -2, 1, 2 }
        };

        double[] b = { 8, -11, -3 };

        double[] result = GaussianElimination(A, b);

        Console.WriteLine("Soluciones:");
        for (int i = 0; i < result.Length; i++)
        {
            Console.WriteLine($"x{i + 1} = {result[i]:0.0000}");
        }
    }

    static double[] GaussianElimination(double[,] A, double[] b)
    {
        int n = b.Length;

        // Construir la matriz aumentada
        double[,] augmentedMatrix = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmentedMatrix[i, j] = A[i, j];
            }
            augmentedMatrix[i, n] = b[i];
        }

        // Aplicar eliminación gaussiana
        for (int k = 0; k < n; k++)
        {
            // Encontrar el pivote
            double pivot = augmentedMatrix[k, k];
            if (pivot == 0)
            {
                // Intercambiar la fila k con una fila i donde augmentedMatrix[i, k] != 0
                for (int i = k + 1; i < n; i++)
                {
                    if (augmentedMatrix[i, k] != 0)
                    {
                        SwapRows(augmentedMatrix, k, i);
                        pivot = augmentedMatrix[k, k];
                        break;
                    }
                }
            }

            // Hacer cero los elementos debajo del pivote
            for (int i = k + 1; i < n; i++)
            {
                double factor = augmentedMatrix[i, k] / pivot;
                for (int j = k; j <= n; j++)
                {
                    augmentedMatrix[i, j] -= factor * augmentedMatrix[k, j];
                }
            }
        }

        // Sustitución hacia atrás
        double[] x = new double[n];
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < n; j++)
            {
                sum += augmentedMatrix[i, j] * x[j];
            }
            x[i] = (augmentedMatrix[i, n] - sum) / augmentedMatrix[i, i];
        }

        return x;
    }

    static void SwapRows(double[,] matrix, int row1, int row2)
    {
        int columns = matrix.GetLength(1);
        for (int i = 0; i < columns; i++)
        {
            double temp = matrix[row1, i];
            matrix[row1, i] = matrix[row2, i];
            matrix[row2, i] = temp;
        }
    }
}
```


### Ejmplos en java

[Ejemplo 1](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/EliminacionGaussiana1.java)

[Ejemplo 2](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/EliminacionGaussiana2.java)

[Ejemplo 3](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/EliminacionGaussiana3.java)

[Ejemplo 4](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/EliminacionGaussiana4.java)

[Ejemplo 5](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/EliminacionGaussiana5.java)


---------------------------------------------------------------------------


 ## Metodo de Gauss-Jordan.
### Descripción: 

La eliminación Gauss-Jordan es una extensión de la eliminación gaussiana, que no solo transforma la matriz aumentada en una forma triangular superior, sino que la reduce a la forma de matriz identidad. Este método se utiliza para resolver sistemas de ecuaciones lineales y para encontrar la inversa de una matriz.

### Pseudocódigo:
```
Algoritmo EliminacionGaussJordan(A, b)
    Entrada: 
        A: Matriz de coeficientes de tamaño n x n
        b: Vector de términos independientes de tamaño n
    Salida: 
        x: Vector solución de tamaño n

    AugmentedMatrix <- [A|b]  # Construir la matriz aumentada
    n <- número de filas de A

    Para k desde 1 hasta n hacer:
        # Encontrar el pivote
        Pivote <- AugmentedMatrix[k][k]
        Si Pivote = 0 entonces
            Intercambiar la fila k con una fila i donde AugmentedMatrix[i][k] != 0
            Pivote <- AugmentedMatrix[k][k]

        # Normalizar la fila del pivote
        Para j desde k hasta n+1 hacer:
            AugmentedMatrix[k][j] <- AugmentedMatrix[k][j] / Pivote

        # Hacer cero los elementos en las otras filas en la columna del pivote
        Para i desde 1 hasta n hacer:
            Si i != k entonces
                Factor <- AugmentedMatrix[i][k]
                Para j desde k hasta n+1 hacer:
                    AugmentedMatrix[i][j] <- AugmentedMatrix[i][j] - Factor * AugmentedMatrix[k][j]

    # Extraer la solución
    Para i desde 1 hasta n hacer:
        x[i] <- AugmentedMatrix[i][n+1]

    Retornar x

```

### Implementacion 
- Implementacion en C#

```
using System;

class Program
{
    static void Main(string[] args)
    {
        double[,] A = {
            { 2, 1, -1 },
            { -3, -1, 2 },
            { -2, 1, 2 }
        };

        double[] b = { 8, -11, -3 };

        double[] result = GaussJordanElimination(A, b);

        Console.WriteLine("Soluciones:");
        for (int i = 0; i < result.Length; i++)
        {
            Console.WriteLine($"x{i + 1} = {result[i]:0.0000}");
        }
    }

    static double[] GaussJordanElimination(double[,] A, double[] b)
    {
        int n = b.Length;

        // Construir la matriz aumentada
        double[,] augmentedMatrix = new double[n, n + 1];
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                augmentedMatrix[i, j] = A[i, j];
            }
            augmentedMatrix[i, n] = b[i];
        }

        // Aplicar eliminación Gauss-Jordan
        for (int k = 0; k < n; k++)
        {
            // Encontrar el pivote
            double pivot = augmentedMatrix[k, k];
            if (pivot == 0)
            {
                // Intercambiar la fila k con una fila i donde augmentedMatrix[i, k] != 0
                for (int i = k + 1; i < n; i++)
                {
                    if (augmentedMatrix[i, k] != 0)
                    {
                        SwapRows(augmentedMatrix, k, i);
                        pivot = augmentedMatrix[k, k];
                        break;
                    }
                }
            }

            // Normalizar la fila del pivote
            for (int j = k; j <= n; j++)
            {
                augmentedMatrix[k, j] /= pivot;
            }

            // Hacer cero los elementos en las otras filas en la columna del pivote
            for (int i = 0; i < n; i++)
            {
                if (i != k)
                {
                    double factor = augmentedMatrix[i, k];
                    for (int j = k; j <= n; j++)
                    {
                        augmentedMatrix[i, j] -= factor * augmentedMatrix[k, j];
                    }
                }
            }
        }

        // Extraer la solución
        double[] x = new double[n];
        for (int i = 0; i < n; i++)
        {
            x[i] = augmentedMatrix[i, n];
        }

        return x;
    }

    static void SwapRows(double[,] matrix, int row1, int row2)
    {
        int columns = matrix.GetLength(1);
        for (int i = 0; i < columns; i++)
        {
            double temp = matrix[row1, i];
            matrix[row1, i] = matrix[row2, i];
            matrix[row2, i] = temp;
        }
    }
}

```


### Ejmplos en java

[Ejemplo 1](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussJordan1.java)

[Ejemplo 2](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussJordan2.java)

[Ejemplo 3](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussJordan3.java)

[Ejemplo 4](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussJordan4.java)

[Ejemplo 5](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussJordan5.java)


---------------------------------------------------------------------------

 ## Método de Gauss-Seidel.
### Descripción: 

El método de Gauss-Seidel es un método iterativo utilizado para resolver sistemas de ecuaciones lineales. A diferencia de los métodos directos como la eliminación Gaussiana, los métodos iterativos construyen una sucesión de aproximaciones que convergen a la solución del sistema. El método de Gauss-Seidel mejora las aproximaciones sucesivas utilizando los valores más recientes de las variables.

### Pseudocódigo:
```
Algoritmo GaussSeidel(A, b, x, tol, maxIter)
    Entrada:
        A: Matriz de coeficientes de tamaño n x n
        b: Vector de términos independientes de tamaño n
        x: Vector inicial de tamaño n
        tol: Tolerancia para la convergencia
        maxIter: Número máximo de iteraciones
    Salida:
        x: Vector solución de tamaño n

    n <- tamaño de b
    iter <- 0
    error <- tol + 1  # Inicializar el error para entrar en el bucle

    Mientras error > tol Y iter < maxIter hacer
        x_old <- x  # Guardar la solución anterior

        Para i desde 1 hasta n hacer
            suma <- 0
            Para j desde 1 hasta n hacer
                Si j != i entonces
                    suma <- suma + A[i][j] * x[j]
                FinSi
            FinPara
            x[i] <- (b[i] - suma) / A[i][i]
        FinPara

        # Calcular el error como la norma del vector de diferencias
        error <- Norma(x - x_old)
        iter <- iter + 1
    FinMientras

    Retornar x

```

### Implementacion 
- Implementacion en C#

```
using System;

class Program
{
    static void Main(string[] args)
    {
        double[,] A = {
            { 4, 1, 2 },
            { 3, 5, 1 },
            { 1, 1, 3 }
        };

        double[] b = { 4, 7, 3 };
        double[] x = { 0, 0, 0 };  // Vector inicial
        double tol = 1e-10;
        int maxIter = 1000;

        double[] result = GaussSeidel(A, b, x, tol, maxIter);

        Console.WriteLine("Soluciones:");
        for (int i = 0; i < result.Length; i++)
        {
            Console.WriteLine($"x{i + 1} = {result[i]:0.0000000000}");
        }
    }

    static double[] GaussSeidel(double[,] A, double[] b, double[] x, double tol, int maxIter)
    {
        int n = b.Length;
        double[] x_old = new double[n];
        double error = tol + 1;
        int iter = 0;

        while (error > tol && iter < maxIter)
        {
            Array.Copy(x, x_old, n);

            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += A[i, j] * x[j];
                    }
                }
                x[i] = (b[i] - sum) / A[i, i];
            }

            error = 0;
            for (int i = 0; i < n; i++)
            {
                error += Math.Abs(x[i] - x_old[i]);
            }

            iter++;
        }

        return x;
    }
}


```


### Ejmplos en java

[Ejemplo 1](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussSeidel.java)

[Ejemplo 2](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussSeidel2.java)

[Ejemplo 3](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussSeidel3.java)

[Ejemplo 4](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussSeidel4.java)

[Ejemplo 5](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/GaussSeidel5.java)


---------------------------------------------------------------------------

 ## Método de Jacobi.
### Descripción: 

El método de Jacobi es un método iterativo para resolver sistemas de ecuaciones lineales. Es similar al método de Gauss-Seidel, pero a diferencia de este, el método de Jacobi utiliza solo los valores de la iteración anterior para actualizar las variables. Esto significa que todas las actualizaciones se realizan simultáneamente en cada paso de iteración.

### Pseudocódigo:
```
Algoritmo Jacobi(A, b, x, tol, maxIter)
    Entrada:
        A: Matriz de coeficientes de tamaño n x n
        b: Vector de términos independientes de tamaño n
        x: Vector inicial de tamaño n
        tol: Tolerancia para la convergencia
        maxIter: Número máximo de iteraciones
    Salida:
        x: Vector solución de tamaño n

    n <- tamaño de b
    x_new <- vector de tamaño n inicializado en 0
    iter <- 0
    error <- tol + 1  # Inicializar el error para entrar en el bucle

    Mientras error > tol Y iter < maxIter hacer
        Para i desde 1 hasta n hacer
            suma <- 0
            Para j desde 1 hasta n hacer
                Si j != i entonces
                    suma <- suma + A[i][j] * x[j]
                FinSi
            FinPara
            x_new[i] <- (b[i] - suma) / A[i][i]
        FinPara

        # Calcular el error como la norma del vector de diferencias
        error <- Norma(x_new - x)
        x <- x_new
        iter <- iter + 1
    FinMientras

    Retornar x


```

### Implementacion 
- Implementacion en C#

```
using System;

class Program
{
    static void Main(string[] args)
    {
        double[,] A = {
            { 4, 1, 2 },
            { 3, 5, 1 },
            { 1, 1, 3 }
        };

        double[] b = { 4, 7, 3 };
        double[] x = { 0, 0, 0 };  // Vector inicial
        double tol = 1e-10;
        int maxIter = 1000;

        double[] result = Jacobi(A, b, x, tol, maxIter);

        Console.WriteLine("Soluciones:");
        for (int i = 0; i < result.Length; i++)
        {
            Console.WriteLine($"x{i + 1} = {result[i]:0.0000000000}");
        }
    }

    static double[] Jacobi(double[,] A, double[] b, double[] x, double tol, int maxIter)
    {
        int n = b.Length;
        double[] x_new = new double[n];
        double error = tol + 1;
        int iter = 0;

        while (error > tol && iter < maxIter)
        {
            for (int i = 0; i < n; i++)
            {
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        sum += A[i, j] * x[j];
                    }
                }
                x_new[i] = (b[i] - sum) / A[i, i];
            }

            error = 0;
            for (int i = 0; i < n; i++)
            {
                error += Math.Abs(x_new[i] - x[i]);
            }

            Array.Copy(x_new, x, n);
            iter++;
        }

        return x;
    }
}


```


### Ejmplos en java

[Ejemplo 1](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/Jacobi1.java)

[Ejemplo 2](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/Jacobi2.java)

[Ejemplo 3](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/Jacobi3.java)

[Ejemplo 4](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/Jacobi4.java)

[Ejemplo 5](https://github.com/BrandonHDZ2102/T3-E2-Problemario/blob/main/ProblemarioT3/src/problemariot3/Jacobi5.java)


---------------------------------------------------------------------------
