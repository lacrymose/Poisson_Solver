#include "Poisson.cpp"

int main ()
{
    deallog.depth_console(0);
    {
        PoissonProblem::Poisson<2> poisson_2d(1,4);
        poisson_2d.run();
    }

    return 0;
}
