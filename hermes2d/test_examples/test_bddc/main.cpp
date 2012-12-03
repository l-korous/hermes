#define HERMES_REPORT_ALL
#include "definitions.h"
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "mpi.h"
#include <boost/lexical_cast.hpp>

extern "C" {
#include "/home/lukas/dependencies/include/bddcml_interface_c.h"
}

int main(int argc, char* argv[])
{
    MPI_Init(&argc, &argv);

    int commFortran = MPI_Comm_c2f(MPI_COMM_WORLD);

    // Orient in the communicator
    int nProc;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nProc );

    MPI_Comm_rank(MPI_COMM_WORLD, &rank );

    // BDDCML data initialization.
    int la = 64;
    int nelem = 6, nnod = 12, ndof = 12, ndim = 2, meshdim = 2, nelems = 3, nnods = 8, ndofs = 8, linet = 12, lnnet = 3, lnndf = 8, lisngn = 8, lisvgvn = 8, lisegn = 3,
            lxyz1 = 8, lxyz2 = 2, lifix = 8, lfixv = 8, lrhs = 8, is_rhs_complete = 1, lsol = 8, matrixtype = 0, is_assembled_int = 1;
    int *nnet = new int[3], *nndf = new int[8], *inet = new int[12], *isngn = new int[8], *isvgvn = new int[8], *isegn = new int[3], *ifix = new int[8], *i_sparse = new int[64], *j_sparse = new int[64];
    double *fixv = new double[8], *rhs = new double[8], *sol = new double[8], *a_sparse = new double[64], *xyz = new double[16];
    nnet[0] = nnet[1] = nnet[2] = 4;
    nndf[0] = nndf[1] = nndf[2] = nndf[3] = nndf[4] = nndf[5] = nndf[6] = nndf[7] = 1;
    sol[0] = sol[1] = sol[2] = sol[3] = sol[4] = sol[5] = sol[6] = sol[7] = 0;

    // For constraints (not used).
    int number = 0;

    std::map<int, int> sub0globaltolocal;
    sub0globaltolocal.insert(std::pair<int, int>(3, 0));
    sub0globaltolocal.insert(std::pair<int, int>(2, 1));
    sub0globaltolocal.insert(std::pair<int, int>(4, 2));
    sub0globaltolocal.insert(std::pair<int, int>(10, 3));
    sub0globaltolocal.insert(std::pair<int, int>(8, 4));
    sub0globaltolocal.insert(std::pair<int, int>(7, 5));
    sub0globaltolocal.insert(std::pair<int, int>(6, 6));
    sub0globaltolocal.insert(std::pair<int, int>(11, 7));

    std::map<int, int> sub1globaltolocal;
    sub1globaltolocal.insert(std::pair<int, int>(0, 0));
    sub1globaltolocal.insert(std::pair<int, int>(1, 1));
    sub1globaltolocal.insert(std::pair<int, int>(4, 2));
    sub1globaltolocal.insert(std::pair<int, int>(9, 3));
    sub1globaltolocal.insert(std::pair<int, int>(3, 4));
    sub1globaltolocal.insert(std::pair<int, int>(2, 5));
    sub1globaltolocal.insert(std::pair<int, int>(5, 6));
    sub1globaltolocal.insert(std::pair<int, int>(10, 7));

    std::cout << " I am here. " << std::endl;

    ///// Tahle cast je jen assembling matice //
    Hermes::Hermes2D::Mesh mesh;
    // Load the mesh.
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("domain.xml", &mesh);

    Views::MeshView m;
    Hermes::Hermes2D::H1Space<double> space(&mesh, 1);

    // Initialize the weak formulation.
    CustomWeakFormPoisson wf(new Hermes::Hermes2DFunction<double>(-10));

    // Initialize the solution.
    Hermes::Hermes2D::Solution<double> sln;

    // Initialize linear solver.
    Hermes::Hermes2D::DiscreteProblemLinear<double> dp(&wf, &space);

    /// Jacobian.
    SparseMatrix<double>* jacobian = Algebra::create_matrix<double>();

    /// Residual.
    Vector<double>* residual = Algebra::create_vector<double>();

    // Solve the linear problem.
    dp.assemble(jacobian, residual);

    int verbose_level = 2;
    int numbase = 0;
    std::vector<int> subdomainsNumber(2);
    subdomainsNumber[0] = 2;
    subdomainsNumber[1] = 1;
    int subdomainsNumberTotal = 2;
    int levelsNumber = 2;
    int one = 1;

    bddcml_init(&levelsNumber, &(subdomainsNumber[0]), &subdomainsNumberTotal, &one, &commFortran, &verbose_level, &numbase);

    std::string numString = boost::lexical_cast<std::string>( rank );
    std::string rhsFile = "rhs" + numString;
    std::string iSpFile = "iSparse" + numString;
    std::string jSpFile = "jSparse" + numString;
    std::string aSpFile = "aSparse" + numString;
    std::ofstream outRhs(rhsFile.c_str(), std::ofstream::out );
    std::ofstream outISparse(iSpFile.c_str(), std::ofstream::out );
    std::ofstream outJSparse(jSpFile.c_str(), std::ofstream::out );
    std::ofstream outASparse(aSpFile.c_str(), std::ofstream::out );

    // subdomain 0
    if(rank == 0)
    {
        inet[0] = 0; inet[1] = 1; inet[2] = 5; inet[3] = 4; inet[4] = 1; inet[5] = 2; inet[6] = 6; inet[7] = 5; inet[8] = 2; inet[9] = 3; inet[10] = 7; inet[11] = 6;
        isngn[0] = 3; isngn[1] = 2; isngn[2] = 5; isngn[3] = 10; isngn[4] = 8; isngn[5] = 7; isngn[6] = 6; isngn[7] = 11;
        isvgvn[0] = 3; isvgvn[1] = 2; isvgvn[2] = 5; isvgvn[3] = 10; isvgvn[4] = 8; isvgvn[5] = 7; isvgvn[6] = 6; isvgvn[7] = 11;
        isegn[0] = 3; isegn[1] = 4; isegn[2] = 5;
        xyz[0] = 0; xyz[1] = 1; xyz[2] = 2; xyz[3] = 3; xyz[4] = 0; xyz[5] = 1; xyz[6] = 2; xyz[7] = 3; xyz[8] = 1; xyz[9] = 1; xyz[10] = 1; xyz[11] = 1; xyz[12] = 2; xyz[13] = 2; xyz[14] = 2; xyz[15] = 2;
        ifix[0] = 1; ifix[1] = 0; ifix[2] = 0; ifix[3] = 1; ifix[4] = 1; ifix[5] = 1; ifix[6] = 1; ifix[7] = 1;
        fixv[0] = 20; fixv[1] = 20; fixv[2] = 20; fixv[3] = 20; fixv[4] = 20; fixv[5] = 20; fixv[6] = 20; fixv[7] = 20;

        for(int i = 0; i < residual->length(); i++)
        {
            if(sub0globaltolocal.find(i) == sub0globaltolocal.end())
                continue;
            rhs[sub0globaltolocal.find(i)->second] = residual->get(i);
        }

        la = 0;
        for(int i = 0; i < jacobian->get_size(); i++)
        {
            if(sub0globaltolocal.find(i) == sub0globaltolocal.end())
                continue;
            for(int j = 0; j < jacobian->get_size(); j++)
            {
                if(sub0globaltolocal.find(j) == sub0globaltolocal.end())
                    continue;
                i_sparse[la] = sub0globaltolocal.find(i)->second;
                j_sparse[la] = sub0globaltolocal.find(j)->second;
                a_sparse[la++] = jacobian->get(i, j);
            }
        }
    }
    // subdomain 1
    else
    {
        inet[0] = 0; inet[1] = 1; inet[2] = 5; inet[3] = 4; inet[4] = 1; inet[5] = 2; inet[6] = 6; inet[7] = 5; inet[8] = 2; inet[9] = 3; inet[10] = 7; inet[11] = 6;
        isngn[0] = 0; isngn[1] = 1; isngn[2] = 4; isngn[3] = 9; isngn[4] = 3; isngn[5] = 2; isngn[6] = 5; isngn[7] = 10;
        isvgvn[0] = 0; isvgvn[1] = 1; isvgvn[2] = 4; isvgvn[3] = 9; isvgvn[4] = 3; isvgvn[5] = 2; isvgvn[6] = 5; isvgvn[7] = 10;
        isegn[0] = 0; isegn[1] = 1; isegn[2] = 2;
        xyz[0] = 0; xyz[1] = 1; xyz[2] = 2; xyz[3] = 3; xyz[4] = 0; xyz[5] = 1; xyz[6] = 2; xyz[7] = 3; xyz[8] = 0; xyz[9] = 0; xyz[10] = 0; xyz[11] = 0; xyz[12] = 1; xyz[13] = 1; xyz[14] = 1; xyz[15] = 1;
        ifix[0] = 1; ifix[1] = 1; ifix[2] = 1; ifix[3] = 1; ifix[4] = 1; ifix[5] = 0; ifix[6] = 0; ifix[7] = 1;
        fixv[0] = 20; fixv[1] = 20; fixv[2] = 20; fixv[3] = 20; fixv[4] = 20; fixv[5] = 20; fixv[6] = 20; fixv[7] = 20;

        for(int i = 0; i < residual->length(); i++)
        {
            if(sub1globaltolocal.find(i) == sub1globaltolocal.end())
                continue;
            rhs[sub1globaltolocal.find(i)->second] = residual->get(i);
        }

        la = 0;
        for(int i = 0; i < jacobian->get_size(); i++)
        {
            if(sub1globaltolocal.find(i) == sub1globaltolocal.end())
                continue;
            for(int j = 0; j < jacobian->get_size(); j++)
            {
                if(sub1globaltolocal.find(j) == sub1globaltolocal.end())
                    continue;
                i_sparse[la] = sub1globaltolocal.find(i)->second;
                j_sparse[la] = sub1globaltolocal.find(j)->second;
                a_sparse[la++] = jacobian->get(i, j);
            }
        }
    }

    la = 64;

    for(int i = 0; i < lrhs; i++)
        outRhs << rhs[i] << std::endl;
    outRhs.close();

    for(int i = 0; i < la; i++)
        outISparse << i_sparse[i] << std::endl;
    outISparse.close();

    for(int i = 0; i < la; i++)
        outJSparse << j_sparse[i] << std::endl;
    outJSparse.close();

    for(int i = 0; i < la; i++)
        outASparse << a_sparse[i] << std::endl;
    outASparse.close();

    std::cout << "I am " << rank << std::endl;
    std::cout << "rhsValues" << std::endl;
    for(int i = 0; i < 8; i++)
        std::cout << rhs[i] << std::endl;
    std::cout << std::endl;

    std::cout << "iSparse" << std::endl;
    for(int i = 0; i < la; i++)
        std::cout << i_sparse[i] << std::endl;
    std::cout << std::endl;

    std::cout << "jSparse" << std::endl;
    for(int i = 0; i < la; i++)
        std::cout << j_sparse[i] << std::endl;
    std::cout << std::endl;

    std::cout << "aSparse" << std::endl;
    for(int i = 0; i < la; i++)
        std::cout << a_sparse[i] << std::endl;

    std::cout << "I am rank: " << numString << std::endl;

    std::ifstream inRhs(rhsFile.c_str(), std::ifstream::in );
    for(int i = 0; i < lrhs; i++)
        inRhs >> rhs[i];
    inRhs.close();

    std::ifstream inISparse(iSpFile.c_str(), std::ifstream::in );
    for(int i = 0; i < la; i++)
        inISparse >> i_sparse[i];
    inISparse.close();

    std::ifstream inJSparse(jSpFile.c_str(), std::ifstream::in );
    for(int i = 0; i < la; i++)
        inJSparse >> j_sparse[i];
    inJSparse.close();

    std::ifstream inASparse(aSpFile.c_str(), std::ifstream::in );
    for(int i = 0; i < la; i++)
        inASparse >> a_sparse[i];
    inASparse.close();

    std::cout << "la = " << la << std::endl;
    for(int i = 0; i < la; i++)
        std::cout << i_sparse[i] << " " << j_sparse[i] << " " << a_sparse[i] << std::endl;

    double constraints = 1.0;
    bddcml_upload_subdomain_data(&nelem, &nnod, &ndof, &ndim, &meshdim, &rank, &nelems, &nnods, &ndofs, inet, &linet, nnet, &lnnet, nndf, &lnndf, isngn, &lisngn, isvgvn, &lisvgvn, isegn, &lisegn, xyz, &lxyz1, &lxyz2, ifix, &lifix, fixv, &lfixv, rhs, &lrhs, &is_rhs_complete, sol, &lsol, &matrixtype, i_sparse, j_sparse, a_sparse, &la, &is_assembled_int, &constraints, &number, &number);

    int use_defaults_int          = 1;
    int parallel_division_int     = 0;
    int use_arithmetic_int        = 0;
    int use_user_constraints_int  = 0;
    int use_adaptive_int = 0;

    std::cout << "About to setup preconditioner." << std::endl;

    // setup multilevel BDDC preconditioner
    bddcml_setup_preconditioner( & matrixtype, & use_defaults_int, & parallel_division_int, & use_arithmetic_int, & use_adaptive_int, & use_user_constraints_int );

    std::cout << "About to solve." << std::endl;

    int method = -1;
    int iterations;
    int covergenceReasons;
    double conditionNumber;
    double tol = 10e-6;
    bddcml_solve(&commFortran, &method, &tol, &method, &method, &iterations, &covergenceReasons, &conditionNumber);

    int length = 12;
    double * solution_vector;
    if ( rank == 0 )
        solution_vector = new double[length];

    bddcml_download_global_solution(solution_vector, &length);
    if ( rank == 0 )
    {
        for(int i_l = 0; i_l < length; i_l++)
            std::cout << solution_vector[i_l] << std::endl;

        Hermes2D::Solution<double>::vector_to_solution(solution_vector, &space, &sln);

        // Visualize the solution.
        Hermes::Hermes2D::Views::ScalarView viewS(rhsFile.c_str(), new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

        viewS.show(&sln, Hermes::Hermes2D::Views::HERMES_EPS_LOW);
        viewS.wait_for_close();
    }

    bddcml_finalize();

    MPI_Finalize();
    return 0;
}
