//
// Created by frongere on 14/12/17.
//

#include "MathUtils.h"

using namespace mathutils;


void PrintHeader(std::string title) {
    std::cout << "\n=====================================================================" << std::endl;
    std::cout << "    " << title << std::endl;
    std::cout << "=====================================================================" << std::endl;
}

void PrintInfo(std::string info) {
    std::cout << info << ":" << std::endl;
}


int main(int argc, char* argv[]) {


    // ========================================================================
    //    Matrix creation, initialization and extractions (indexing)
    // ========================================================================
    PrintHeader("Instantiating a matrix with 6 rows and 3 columns");

    // Instatiating a matrix with its dimensions (it is set to zeros by default)
    auto myMatrix = MatrixMN<double>(6, 3);

//    // Displaying the matrix
    std::cout << myMatrix << "\n\n";

    // Getting dimensions of the matrix
    std::cout << "Nb Rows : " << myMatrix.GetNbRows() << std::endl << std::endl;
    std::cout << "Nb Cols : " << myMatrix.GetNbCols() << std::endl << std::endl;

    PrintInfo("Randomizing coefficients");
    // Randomizing the coefficients
    myMatrix.Randomize();
    std::cout << myMatrix << "\n\n";

    // Getting individual coefficients (different means)
    PrintInfo("Getting and setting individual coefficients by indices");
    std::cout << "Coeff 0, 0 : " << myMatrix.at(0, 0) << "\n\n";
    std::cout << "Coeff 0, 0 : " << myMatrix(0, 0) << "\n\n";

    myMatrix.at(0, 0) = 100;
    std::cout << "Coeff 0, 0 (100): " << myMatrix(0, 0) << "\n\n";

    myMatrix(1, 1) = 200;
    std::cout << "Coeff 1, 1 (200): " << myMatrix(1, 1) << "\n\n";

    std::cout << myMatrix << "\n\n";

    PrintInfo("Matrix transposition");
    // Transpose the matrix inplace
    myMatrix.Transpose();
    std::cout << myMatrix << "\n\n";

    std::cout << "Nb Rows : " << myMatrix.GetNbRows() << std::endl << std::endl;
    std::cout << "Nb Cols : " << myMatrix.GetNbCols() << std::endl << std::endl;

    // TODO: montrer utilisation de la fonction transpose

    PrintInfo("Diagonal extraction");
    // Get the diagonal
    auto diag = myMatrix.GetDiag();
    std::cout << diag << "\n\n";

    PrintInfo("Column extraction by index");
    // Get a column
    auto col2 = myMatrix.GetColumn(2);
    std::cout << col2 << "\n\n";

    PrintInfo("Row extraction by index");
    // Get a row
    auto row1 = myMatrix.GetRow(1);
    std::cout << row1 << "\n\n";

    PrintInfo("Setting coeffs to 0");
    // Setting every coeffs to 0
    myMatrix.SetNull();
    std::cout << myMatrix << "\n\n";

    PrintInfo("Set identity matrix");
    // Making it identity (even in square matrices)
    myMatrix.SetIdentity();
    std::cout << myMatrix << "\n\n";

    // ========================================================================
    //    Square matrix inversion
    // ========================================================================

    PrintHeader("Inversion of square matrices");
    MatrixMN<double> matSquare(6, 6);
    matSquare.Randomize();
    PrintInfo("The square Matrix to inverse");
    std::cout << matSquare << "\n\n";

    PrintInfo("Inverse of this matrix inplace");
    // Inverse of the matrix (inplace)
    auto initialMatrix = matSquare;
    matSquare.Inverse();
    std::cout << matSquare << "\n\n";

    PrintInfo("Inverse of the same matrix by copy (we come back to the original...)");
    // Inverse of the Matrix (copy)
    auto invMatrix = matSquare.GetInverse();
    std::cout << invMatrix << "\n\n"; // Must get back the initial matrix

    PrintInfo("Verifications ( (A^-1)^-1 - A)");
    std::cout << invMatrix - initialMatrix << "\n\n"; // Must be zero

    PrintInfo("Verifications ( A^-1 * A)");
    std::cout << invMatrix * matSquare << "\n\n";   // Must be identity

    // Is the matrix positive semi definite
    auto isPSD = matSquare.IsPositiveSemiDefinite();


    // ========================================================================
    //    QR Decomposition
    // ========================================================================

    PrintHeader("Decompositions");
    PrintHeader("QR Decomposition (A = QR)");
    myMatrix.Transpose();
    myMatrix.Randomize();

    PrintInfo("The matrix A to decompose");
    std::cout << myMatrix << "\n\n";

    // QR Decomposition
    MatrixMN<double> Q, R;
    myMatrix.GetQRDecomposition(Q, R);

    PrintInfo("Q");
    std::cout << Q << "\n\n";
    PrintInfo("R");
    std::cout << R << "\n\n";
    PrintInfo("Verification (A - Q*R)");
    std::cout << myMatrix - Q*R << "\n\n";

    // ========================================================================
    //    LU Decomposition
    // ========================================================================

    PrintHeader("LU Decomposition (PA = LU) for A square");
    PrintInfo("The matrix A to decompose");
    matSquare.Randomize();
    std::cout << matSquare << "\n\n";

    // LU decomposition
    MatrixMN<double> P, L, U;
    matSquare.GetLUDecomposition(P, L, U);

    PrintInfo("P");
    std::cout << P << "\n\n";
    PrintInfo("L");
    std::cout << L << "\n\n";
    PrintInfo("U");
    std::cout << U << "\n\n";
    PrintInfo("Verification (P*A - L*U)");
    std::cout << P * matSquare - L*U << "\n\n"; // Must be null !

    // ========================================================================
    //    Cholesky Decomposition
    // ========================================================================

//    PrintHeader("Cholesky Decomposition A = L*L^T");
//    // Cholesky decomposition
//    myMatrix.Randomize();
//    myMatrix *= myMatrix;  // Making sure we have a matrix semi positive definite
//    std::cout << myMatrix << "\n\n";
//
//    MatrixMN<double> Lchol;
//    myMatrix.GetCholeskyDecomposition(Lchol);
//
//    PrintInfo("L");
//    std::cout << Lchol << "\n\n";
//    std::cout << Transpose(Lchol) << "\n\n";
//    std::cout << myMatrix << "\n\n";
//    std::cout << Lchol*Transpose(Lchol) << "\n\n";  // FIXME: ne donne pas le bon resultat !!!
////    std::cout << myMatrix - L * Transpose(L) << "\n\n";

    // ========================================================================
    //    SVD Decomposition
    // ========================================================================

    PrintHeader("SVD Decomposition A = U * Sdiag * (V*)");
    // SVD decomposition
    PrintInfo("The matrix to decompose");
    myMatrix.Randomize();
    std::cout << myMatrix << "\n\n";

    MatrixMN<double> Usvd, SingVal, Vsvd;
    myMatrix.GetSVDDecomposition(Usvd, SingVal, Vsvd);

    PrintInfo("U");
    std::cout << Usvd << "\n\n";
    PrintInfo("Singular values");
    std::cout << SingVal << "\n\n";
    PrintInfo("V");
    std::cout << Vsvd << "\n\n";

    PrintInfo("Verification (A - U * Sdiag * V*)");
    std::cout << myMatrix - Usvd * SingVal.asDiagonal() * Vsvd.adjoint() << "\n\n";

    // ========================================================================
    //    Pseudo inverse computation
    // ========================================================================

    PrintHeader("Pseudo inverse of a rectangular matrix");
    // Pseudo inverse
    auto matRect = MatrixMN<double>(10, 5);
    matRect.Randomize();
    auto matRect_pinv = matRect.GetPseudoInverse();

    PrintInfo("A+");
    std::cout << matRect_pinv << "\n\n";

    PrintInfo("Using the function version");
    // Fonction Pinv
    auto pinvMat = Pinv(matRect);
    std::cout << pinvMat << "\n\n";


    PrintInfo("Verification (A+ * A)");
    std::cout << matRect_pinv * matRect << "\n\n";



    return 0;
}