#include "fvCFD.H"
#include <xGEMS/ChemicalEngine.hpp>
#include <eigen3/Eigen/Dense>

#include <exception>

int main(int argc, char* argv[])
{
    if (argc != 2)
    {
        Foam::Info << "Usage: xgemsFoamProbe system-dat.lst" << Foam::endl;
        return 1;
    }

    try
    {
        xGEMS::ChemicalEngine engine(argv[1]);
        const Eigen::VectorXd bulk = engine.elementAmounts();
        const int status = engine.equilibrate
        (
            engine.temperature() + 1.0,
            engine.pressure(),
            bulk
        );

        Foam::Info << "xGEMS status: " << status
                   << ", converged: " << engine.converged()
                   << ", iterations: " << engine.numIterations()
                   << Foam::endl;

        return engine.converged() ? 0 : 2;
    }
    catch (const std::exception& error)
    {
        Foam::Info << "xGEMS error: " << error.what() << Foam::endl;
        return 3;
    }
}
