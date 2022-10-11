#ifndef ELASTIC_PROBLEM_H_
#define ELASTIC_PROBLEM_H_


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>

#include<fstream>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/lac/sparse_direct.h>
#include "material_class.h"
#include "perturbation_class.h"



// In this example, we need vector-valued finite elements. The support for
// these can be found in the following include file:
#include <deal.II/fe/fe_system.h>
// We will compose the vector-valued finite elements from regular Q1 elements
// which can be found here, as usual:
#include <deal.II/fe/fe_q.h>
#include <deal.II/base/tensor_function.h>

// This again is C++:
#include <fstream>
#include <iostream>

#define DIM 2

namespace Step8
{
using namespace dealii;




class ElasticProblem
{
public:
	ElasticProblem();

	void run_iterative();
	void run_criticalsurface();
	void run_forwardsolve();
	double get_Strain_Energy();

private:



	void setup_system();
	void assemble_nonlinear_system(const unsigned int);
	void solve();
	void output_results(const unsigned int cycle) const;
	void output_zero_results(const unsigned int cycle) const;
	void solve_linear_system();
	int solve_linear_system_critical(int);
	void output_perturbation(const unsigned int cycle) const;
	void setup_constraints();

	void feasibilitycalc();

	void As_function(const std::vector<Point<DIM>> &, std::vector<Tensor<2, DIM>> & );
	void Ab_function(const std::vector<Point<DIM>> &, std::vector<Tensor<2, DIM>> & );

	void Initialize_Internal_Variables();
	void Update_Internal_Variables(const int);

	void FlipConfigUD(Vector<double> &);



	Tensor<1,DIM> basis_vector(const int);
	Tensor<2,DIM> outer_product(const Tensor<1,DIM> &, const Tensor<1,DIM> &);
	double Tensor_Inner(const Tensor<2,DIM> &, const Tensor<2,DIM> &);
	double Det2D(const Tensor<2,DIM> & );
	Tensor<2,DIM> Tensor_Transpose(const Tensor<2,DIM> & );
	Tensor<2,DIM> cofactor2D(const Tensor<2,DIM> & F);
	double BilinearProduct(const Tensor<2,DIM> &, const Tensor<4,DIM> & , const Tensor<2,DIM> &);
	double Tensor_Trace(const Tensor<2,DIM> & );
	Tensor<2,DIM> Tensor_Sym(const Tensor<2,DIM> & );
	void solution_stabilization(const std::vector<Point<DIM>> &,
			std::vector<Tensor<1, DIM+1>> &  );
	void right_hand_side(const std::vector<Point<DIM>> &,
			std::vector<Tensor<1, DIM+5>> & );
	void Youngs_Modulus(const std::vector<Point<DIM>> &,
				std::vector<double> &  );

	std::vector<double> linspace(double , double , int );
	void OutputVector(std::vector<double> & , std::vector<double> &, std::string);
	Tensor<1,DIM+1> crossproduct(const Tensor<1,DIM+1> &,const Tensor<1,DIM+1> &);
	double TensorNorm1D(const Tensor<1,DIM+1> &);
	double TensorNorm2D(const Tensor<2,DIM> &);

	bool CharacteristicFunction(const Point<DIM> &);
	double stab_magnitude = 0.0;
	double perturbation_magnitude = 0.0;
	unsigned int number_dofs = 0;
	double alpha = 1.0;

	double StrainEnergy = 0.0;



	Tensor<1,DIM+1> Light_Direction;
	Tensor<2,DIM> DirectionTens;

	double T = 1.0;
	double dT = 0.01;

	double Light_Magnitude_As = 1.0;
	double Light_Magnitude_Ab = 1.0;

	Tensor<1,DIM> offset;

	const double nu_global = 0.49;
	const double hv = 0.01;

	double internal_drag = 1.0;
	Triangulation<DIM> triangulation;
	DoFHandler<DIM>    dof_handler;

	FESystem<DIM> fe;

	AffineConstraints<double> constraints;

	SparsityPattern      sparsity_pattern;
	SparseMatrix<double> system_matrix;

	Vector<double> solution;
	Vector<double> step_direction;
	Vector<double> system_rhs;
	Vector<double> residual;
	Vector<double> evaluation_point;
	Vector<double> solution_old;
	Vector<double> perturbation_direction;
	Vector<double> prev_solution;

	Point<DIM> corner1, corner2;
	std::vector<unsigned int>  grid_dimensions;

	std::vector<Material_Class> Material_Vector_InPlane;
	std::vector<Material_Class> Material_Vector_Bending;

	//double hsc = pow(hv,2.0)/(12*(1.0 - pow(nu_global,2.0)));
	double hsc = 1.458489626609542e-05;
	double hscinv = 1.0/hsc;
	double Eval = 1.e2;

	double kappa_Global = 0.0;

	double epsilon_Global = 0.0;




	Tensor<2,2> As_Global;
	Tensor<2,2> Ab_Global;

	std::vector<std::vector<Tensor<2,2>>> Ab_Internal;
	std::vector<std::vector<Tensor<2,2>>> Ab0_Internal;
	std::vector<std::vector<Tensor<2,2>>> As_Internal;
	std::vector<std::vector<Tensor<2,2>>> As0_Internal;



	double resnorm = 10.0;

	double rp = 1.;
	double rd = 1.;
	double kd = 1.;
	double tau = 10.;
	double tol = pow(10.0,-10.0);
	double drag_coeff = 1.0;
	double muval = 10000.*Eval;
	const double mumin = 10000.*Eval;
	double L = 1.0;

	void ReadCurvatureFiles(std::vector<double> &,std::vector<double> &,std::vector<double> &);
};


}

#endif /* ELASTIC_PROBLEM_H_ */
