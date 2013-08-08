#include "IPOptProblemInterface.h"
#include "IpIpoptApplication.hpp"

using namespace gnlp;

int main(int argc, char *argv[])
{
  // Set up the problem
  ConstrainedNonLinearProblem cnlp;

  // Create and add the four variables
  Variable *x1 = cnlp.AddVariable("x1", 1, 1, 5);
  Variable *x2 = cnlp.AddVariable("x2", 5, 1, 5);
  Variable *x3 = cnlp.AddVariable("x3", 5, 1, 5);
  Variable *x4 = cnlp.AddVariable("x4", 1, 1, 5);

  // Set the objective function
  Expression *objective =
      new BinarySum(&cnlp,
                    new BinaryProduct(&cnlp,
                                      new BinaryProduct(&cnlp, x1, x4),
                                      new BinarySum(&cnlp,
                                                    new BinarySum(&cnlp, x1, x2),
                                                    x3)),
                    x3);
  cnlp.SetObjective(objective);

  // Set the constraints
  Expression *c1 = new BinaryProduct(&cnlp,
                                     new BinaryProduct(&cnlp, x1, x2),
                                     new BinaryProduct(&cnlp, x3, x4));
  cnlp.AddConstraint(c1, "TEST1", 25, ConstrainedNonLinearProblem::UBINF);

  BigSum *c2 = new BigSum(&cnlp);
  c2->AddSummand(new Square(&cnlp, x1));
  c2->AddSummand(new Square(&cnlp, x2));
  c2->AddSummand(new Square(&cnlp, x3));
  c2->AddSummand(new Square(&cnlp, x4));
  cnlp.AddConstraint(c2, "TEST2", 40, 40);

  // Initialize the problem
  cnlp.SetupProblem(true);

  // Create the IPopt wrapper
  SmartPtr<IPOptProblemInterface> ip = new IPOptProblemInterface(&cnlp);

  // Set up the IPopt problem
  // Create a new instance of IpoptApplication
  //  (use a SmartPtr, not raw)
  // We are using the factory, since this allows us to compile this
  // example with an Ipopt Windows DLL
  SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

  // Change some options
  // Note: The following choices are only examples, they might not be
  //       suitable for your optimization problem.
  app->Options()->SetNumericValue("tol", 1e-9);
  app->Options()->SetStringValue("mu_strategy", "adaptive");
  app->Options()->SetStringValue("output_file", "ipopt.out");
  app->Options()->SetStringValue("derivative_test", "second-order");

  // Intialize the IpoptApplication and process the options
  ApplicationReturnStatus status;
  status = app->Initialize();
  if (status != Solve_Succeeded) {
    printf("\n\n*** Error during initialization!\n");
    return (int) status;
  }

  // Ask Ipopt to solve the problem
  status = app->OptimizeTNLP(GetRawPtr(ip));

  if (status == Solve_Succeeded) {
    printf("\n\n*** The problem solved!\n");
    printf("SOLUTION %12.6f %12.6f %12.6f %12.6f\n",
           x1->Evaluate(), x2->Evaluate(), x3->Evaluate(), x4->Evaluate());
  }
  else {
    printf("\n\n*** The problem FAILED!\n");
  }

  // As the SmartPtrs go out of scope, the reference count
  // will be decremented and the objects will automatically
  // be deleted.

  return (int) status;
}

