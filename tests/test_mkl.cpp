//
// Created by frongere on 10/05/19.
//

#include "MathUtils/MathUtils.h"

#ifdef WITH_MKL
#include "mkl.h"
#include "mkl_df.h"
#endif


int main() {


    #ifdef WITH_MKL
    std::cout << "coucou" << std::endl;
    #endif








    return 0;
}
