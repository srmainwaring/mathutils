//
// Created by frongere on 14/12/17.
//

#include "MathUtils.h"

using namespace mathutils;

int main(int argc, char* argv[]) {

    // Instatiating a matrix with its dimensions (it is set to zeros by default)
    auto myMatrix = MatrixMN<double>(4, 4);

//    // Displaying the matrix
    std::cout << myMatrix << "\n\n";

    // Getting dimensions of the matrix
    std::cout << "Nb Rows : " << myMatrix.GetNbRows() << std::endl << std::endl;
    std::cout << "Nb Cols : " << myMatrix.GetNbCols() << std::endl << std::endl;

    // Randomizing the coefficients
    myMatrix.Randomize();
    std::cout << myMatrix << "\n\n";

    // Getting individual coefficients (different means)
    std::cout << "Coeff 0, 0 : " << myMatrix.at(0, 0) << "\n\n";
    std::cout << "Coeff 0, 0 : " << myMatrix(0, 0) << "\n\n";

    myMatrix.at(0, 0) = 100;
    std::cout << "Coeff 0, 0 (100): " << myMatrix(0, 0) << "\n\n";

    myMatrix(1, 1) = 200;
    std::cout << "Coeff 1, 1 (200): " << myMatrix(1, 1) << "\n\n";

    std::cout << myMatrix << "\n\n";

    // Transpose the matrix inplace
    myMatrix.Transpose();
    std::cout << myMatrix << "\n\n";

    std::cout << "Nb Rows : " << myMatrix.GetNbRows() << std::endl << std::endl;
    std::cout << "Nb Cols : " << myMatrix.GetNbCols() << std::endl << std::endl;

    // TODO: montrer utilisation de la fonction transpose

    // Get the diagonal
    auto diag = myMatrix.GetDiag();
    std::cout << diag << "\n\n";

    // Get a column
    auto col2 = myMatrix.GetColumn(2);
    std::cout << col2 << "\n\n";

    // Get a row
    auto row1 = myMatrix.GetRow(1);
    std::cout << row1 << "\n\n";

    // Setting every coeffs to 0
    myMatrix.SetNull();
    std::cout << myMatrix << "\n\n";

    // Making it identity (even in square matrices)
    myMatrix.SetIdentity();
    std::cout << myMatrix << "\n\n";


    // Inverse of the matrix (inplace)
    auto initialMatrix = myMatrix;
    myMatrix.Inverse();
    std::cout << myMatrix << "\n\n";

    // Inverse of the Matrix (copy)
    auto invMatrix = myMatrix.GetInverse();
    std::cout << invMatrix << "\n\n";

    // Verification:
    std::cout << invMatrix - initialMatrix << "\n\n"; // Must be zero
    std::cout << invMatrix * myMatrix << "\n\n";   // Must be identity

    // Is the matrix positive semi definite
    auto isPSD = myMatrix.IsPositiveSemiDefinite();

    // QR Decomposition
    MatrixMN<double> Q, R;
    myMatrix.GetQRDecomposition(Q, R);

    std::cout << Q << "\n\n";
    std::cout << R << "\n\n";
    std::cout << myMatrix - Q*R << "\n\n";


    // LU decomposition
    MatrixMN<double> P, L, U;
    myMatrix.GetLUDecomposition(P, L, U);

    std::cout << P << "\n\n";
    std::cout << L << "\n\n";
    std::cout << U << "\n\n";
    std::cout << P * myMatrix - L*U << "\n\n"; // Must be null !


    // Cholesky decomposition
//    myMatrix.Randomize();
    myMatrix *= myMatrix;  // Making sure we have a matrix semi positive definite

    MatrixMN<double> Lchol;
    myMatrix.GetCholeskyDecomposition(Lchol);

    std::cout << Lchol << "\n\n";
    std::cout << Transpose(Lchol) << "\n\n";
    std::cout << myMatrix << "\n\n";
    std::cout << Lchol*Transpose(Lchol) << "\n\n";  // FIXME: ne donne pas le bon resultat !!!
//    std::cout << myMatrix - L * Transpose(L) << "\n\n";



    // SVD decomposition
    MatrixMN<double> Usvd, SingVal, Vsvd;
    myMatrix.GetSVDDecomposition(Usvd, SingVal, Vsvd);

    std::cout << Usvd << "\n\n";
    std::cout << SingVal << "\n\n";
    std::cout << Vsvd << "\n\n";

    std::cout << myMatrix - Usvd * SingVal.asDiagonal() * Vsvd.adjoint() << "\n\n";


    // Pseudo inverse
    auto matRect = MatrixMN<double>(10, 5);
    matRect.Randomize();
    auto matRect_pinv = matRect.GetPseudoInverse();

    std::cout << matRect_pinv * matRect << "\n\n";
    std::cout << matRect_pinv << "\n\n";


    // Fonction Pinv
    auto pinvMat = Pinv(matRect);
    std::cout << pinvMat << "\n\n";




    return 0;
}