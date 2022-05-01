#include <stdlib.h>
#include <stdio.h>

typedef struct {int num; int den;} Fraction;
typedef struct {int dim; Fraction** elements;} MatrixType;


/******* FUNCIONS PROTOTYPE *******/

void getMatrix(MatrixType* Matrix);
void printMatrix(MatrixType Matrix);
Fraction calcDet(MatrixType Matrix);
MatrixType findSubMatrix(MatrixType Matrix, int rowToDelete, int colToDelete);
MatrixType findCofactorsMatrix(MatrixType Matrix);
int invert(MatrixType* Matrix);

void simplify(Fraction* q);
int gcd(int max, int min);
Fraction sumFractions(Fraction n1, Fraction n2);
Fraction subtractFractions(Fraction n1, Fraction n2);


/******* MAIN *******/

int main(int argc, char *argv[])
{   
    MatrixType Matrix;

    getMatrix(&Matrix);
    printf("Matrix:\n");
    printMatrix(Matrix);

    // Invert the matrix if it is possible
    if (invert(&Matrix) == 1)
    {
        printf("Inverse matrix:\n");
        printMatrix(Matrix);
    }
    else printf("The given matrix is not invertible.\n");

    return 0;
}


/******* FUNCTIONS *******/

void getMatrix(MatrixType* Matrix)
{   
    int row, col;
    char answer;

    printf("Type the dimension of the matrix > ");
    scanf("%d", &(Matrix->dim));
    printf("\n");

    Matrix->elements = (Fraction**) malloc(Matrix->dim * sizeof(Fraction));  // Allocalte the N rows of the matrix
    for (row = 0; row < Matrix->dim; row++)
    {
        Matrix->elements[row] = (Fraction*) malloc(Matrix->dim * sizeof(Fraction));  // Allocate N elements for each row
    }

    printf("Only rational numbers are accepted (no irrational numbers or complex numbers)\n");
    printf("How to correctly type numebrs: \n");
    printf("1) Positive integer: 3\n");
    printf("2) Negative integer: -2\n");
    printf("3) Positive fraction: 1/3\n");
    printf("4) Negative fraction: -7/4\n");
    printf("Understood? (y) > ");
    scanf("%*c%c", &answer);
    printf("\n");

    // Fill the matrix with numbers (Fractions)
    printf("Type the element in:\n");
    for (row = 0; row < Matrix->dim; row++)
    {
        for (col = 0; col < Matrix->dim; col++)
        {   
            printf("row %d column %d > ", row+1, col+1);
            Matrix->elements[row][col].den = 1;
            scanf("%d/%d", &(Matrix->elements[row][col].num), &(Matrix->elements[row][col].den));
            // If only a number is typed, only the numerator is modified, and the denominator remains 1
        }
    }
    printf("\n");
}


void printMatrix(MatrixType Matrix)
{
    int row, col;
    int sign;  // To bring the denominator sign to the numerator

    for (row = 0; row < Matrix.dim; row++)
    {
        for (col = 0; col < Matrix.dim; col++)
        {   
            if (Matrix.elements[row][col].den < 0)  // to always keep the eventual negative sign at the numerator
            {
                sign = -1;
                Matrix.elements[row][col].den *= -1;
            }
            else sign = 1;
            Matrix.elements[row][col].num *= sign;  // adjust the numerator sign

            if (Matrix.elements[row][col].num >= 0) printf(" ");  // to nicely display all the numebrs
            if (Matrix.elements[row][col].den == 1 || Matrix.elements[row][col].num == 0)  // the number is an integer
            {
                printf("%d\t", Matrix.elements[row][col].num);
            }
            else  // the number is a fraction
            {   
                printf("%d/%d\t", Matrix.elements[row][col].num, Matrix.elements[row][col].den);
            }
        }
        printf("\n");
    }
    printf("\n");
}


/*Returns the determinant of the matrixed given as a parameter.*/
/*The determinant is calculated recursively using the Laplace method always on the first row.*/
Fraction calcDet(MatrixType Matrix)
{   
    MatrixType SubMatrix;
    int col;  // index of the column
    Fraction determinant;
    Fraction partialDeterminant;

    if (Matrix.dim == 1)
    {   
        determinant.num = Matrix.elements[0][0].num;
        determinant.den = Matrix.elements[0][0].den;
        return determinant;
    }

    determinant.num = 0;
    determinant.den = 1;

    for (col = 0; col < Matrix.dim; col++)
    {   
        SubMatrix = findSubMatrix(Matrix, 0, col);
        partialDeterminant = calcDet(SubMatrix); 

        // Moltiplico il partialDeterminant per l'elemento di cui avevo eliminato riga e colonna
        partialDeterminant.num *= Matrix.elements[0][col].num; 
        partialDeterminant.den *= Matrix.elements[0][col].den;

        // Sum (or subtract) the partialDeterminant to (or from) the total one
        if (col % 2 == 0) determinant = sumFractions(determinant, partialDeterminant);
        else determinant = subtractFractions(determinant, partialDeterminant);
    }

    return determinant;
}


/*Returns the matrix given as a parameter, without the specified row and colum.*/ 
MatrixType findSubMatrix(MatrixType Matrix, int rowToDelete, int colToDelete)
{
    MatrixType SubMatrix;
    int row, col;  // for the original matrix
    int new_row, new_col;  // for the new sub matrix

    SubMatrix.dim = Matrix.dim - 1;
    SubMatrix.elements = (Fraction**) malloc(SubMatrix.dim * sizeof(Fraction));  // Allocate the N-1 rows of the sub matrix
    for (new_row = 0; new_row < SubMatrix.dim; new_row++)
    {
        SubMatrix.elements[new_row] = (Fraction*) malloc(SubMatrix.dim * sizeof(Fraction));
        // Allocate N-1 elements for each row
    }

    /* Copy the correct elements from the original matrix in the sub matrix */
    new_row = 0;
    for (row = 0; row < Matrix.dim; row++)
    {   
        if (row != rowToDelete)
        {   
            new_col = 0;
            for (col = 0; col < Matrix.dim; col++)
            {   
                if (col != colToDelete)
                {
                    SubMatrix.elements[new_row][new_col].num = Matrix.elements[row][col].num;
                    SubMatrix.elements[new_row][new_col].den = Matrix.elements[row][col].den;
                    new_col++;
                }
            }
            new_row++;
        }
    }
    return SubMatrix;
}


/*Returns the matrix of the cofactors (already transposed) of the matrix given ad a parameter.*/ 
MatrixType findCofactorsMatrix(MatrixType Matrix)
{   
    MatrixType Cofactors;
    MatrixType SubMatrix;
    int row, col;  // for the cofactors matrix
    int sign;  // sign of the actual cofactor

    Cofactors.dim = Matrix.dim;
    Cofactors.elements = (Fraction**) malloc(Cofactors.dim * sizeof(Fraction));  // Allocalte the N rows of the cofacotrs matrix
    for (row = 0; row < Cofactors.dim; row++)
    {
        Cofactors.elements[row] = (Fraction*) malloc(Cofactors.dim * sizeof(Fraction));
        // Allocate N elements for each row
    }

    for (row = 0; row < Cofactors.dim; row++)
    {
        if (row % 2 == 0) sign = 1;
        else sign = -1;

        for (col = 0; col < Cofactors.dim; col++)
        {   
            SubMatrix = findSubMatrix(Matrix, row, col);
            Cofactors.elements[col][row] = calcDet(SubMatrix);
            Cofactors.elements[col][row].num *= sign;
            // row and col are inverted so that the matrix is already transposed
            sign = sign * (-1);
        }
    }
    return Cofactors;
}


/*Invert the matrix given as a parameter.*/
/*Returns 1 if the matrix has been inverted, 0 if the matrix was not invertible.*/
int invert(MatrixType* Matrix)
{
    Fraction determinant;
    MatrixType Cofactors;
    int row, col;

    determinant = calcDet(*Matrix);
    if (determinant.num == 0) return 0;

    if (Matrix->dim == 1)  // Matrix with a single element
    {   
        Matrix->elements[0][0].num = determinant.den;
        Matrix->elements[0][0].den = determinant.num;
        simplify(&(Matrix->elements[row][col]));
        return 1;
    }

    Cofactors = findCofactorsMatrix(*Matrix);
    for (row = 0; row < Matrix->dim; row++)
    {
        for (col = 0; col < Matrix->dim; col++)
        {
            // It is correct that the numerators and denominators are crossed
            Matrix->elements[row][col].num = determinant.den * Cofactors.elements[row][col].num;
            Matrix->elements[row][col].den = determinant.num * Cofactors.elements[row][col].den;
            simplify(&(Matrix->elements[row][col]));
        }
    }
    return 1;
}


/*Simplify the fraction given as a parameter.*/
void simplify(Fraction* q)
{
    int divisor;  // gcd of numerator and denominator
    int sign_num;
    int sign_den;

    if (q->num == 0) return;

    // Make the numbers positive in order to find the greatest common divisor
    if (q->num >= 0) sign_num = 1;
    else
    {
        sign_num = -1;
        q->num *= -1;
    }
    if (q->den >= 0) sign_den = 1;
    else
    {
        sign_den = -1;
        q->den *= -1;
    }

    divisor = gcd(q->num, q->den);
    q->num /= divisor;
    q->den /= divisor;

    // Adjust the signs
    q->num *= sign_num;
    q->den *= sign_den;
}


/*Find the greatest common divisor of the two numbers given as parameters.*/
int gcd(int max, int min)
{   
    int temp;
    if (min > max)
    {
        temp = max;
        max = min;
        min = temp;
    }
    if (max % min == 0) return min;
    return gcd(min, max % min);
}


/*Return n1 + n2*/
Fraction sumFractions(Fraction n1, Fraction n2)
{
    Fraction sum;
    int sign_den1;
    int sign_den2;

    // Make the numbers positive in order to find the greatest common divisor
    if (n1.den >= 0) sign_den1 = 1;
    else
    {
        sign_den1 = -1;
        n1.den *= -1;
    }
    if (n2.den >= 0) sign_den2 = 1;
    else
    {
        sign_den2 = -1;
        n2.den *= -1;
    }

    sum.den = (n1.den * n2.den) / gcd(n1.den, n2.den);  // least common multiple of the denominators
    sum.num = ((sum.den / n1.den) * n1.num * sign_den1) + ((sum.den / n2.den) * n2.num * sign_den2);
    simplify(&sum);
    return sum;
}


/*Return n1 - n2*/
Fraction subtractFractions(Fraction n1, Fraction n2)
{
    Fraction difference;
    int sign_den1;
    int sign_den2;

    // Make the numbers positive in order to find the greatest common divisor
    if (n1.den >= 0) sign_den1 = 1;
    else
    {
        sign_den1 = -1;
        n1.den *= -1;
    }
    if (n2.den >= 0) sign_den2 = 1;
    else
    {
        sign_den2 = -1;
        n2.den *= -1;
    }

    difference.den = (n1.den * n2.den) / gcd(n1.den, n2.den);  // least common multiple of the denominators
    difference.num = ((difference.den / n1.den) * n1.num * sign_den1) - ((difference.den / n2.den) * n2.num * sign_den2);
    simplify(&difference);
    return difference;
}