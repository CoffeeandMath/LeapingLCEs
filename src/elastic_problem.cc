#ifndef ELASTIC_PROBLEM_CC_
#define ELASTIC_PROBLEM_CC_
#include "elastic_problem.h"




namespace Step8
{
using namespace dealii;





void ElasticProblem::run_forwardsolve(){

	// Setting the geometric parameters
	L = 1.0; // Length of side

	// Setting edge locations of the rectangle
	corner1(0) = -(L+.1)/2.0;
	corner1(1) = -L/2.0;
	corner2(0) = (L+.1)/2.0;
	corner2(1) = L/2.0;

	// Offsetting the defect center
	offset[0] = 0.0*L;
	offset[1] = 0.0*L;

	// Setting up the discretization to 32 elements per edge
	grid_dimensions.resize(2);
	int GridNum = 32;
	grid_dimensions[0] = GridNum;
	grid_dimensions[1] = GridNum;


	GridGenerator::subdivided_hyper_rectangle (triangulation, grid_dimensions, corner1, corner2, false);

	// Outputting basic system information
	std::cout << "Single Cycle Run with Grid Size: " << grid_dimensions[0] << " x " << grid_dimensions[1] << std::endl;


	std::cout << "   Number of active cells:       "
			<< triangulation.n_active_cells() << std::endl;
	setup_system();

	std::cout << "   Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;



	// Sets the spontaneous in plane strain
	double epsilonmax = -0.03;

	// Sets the initial spontaneous curvature
	double kappaneggoal = -0.1;

	// Number of steps in the initialization phase
	int Nsteps1 = 150;
	int NstepsInit = 250;
	// Setting internal drag (penalization of deformation from one time step to the other)
	internal_drag = 100.0;
	Tensor<1,DIM> LCN_Alignment;
	Tensor<2,DIM> Iden2D;
	Iden2D[0][0] = 1.0;
	Iden2D[1][1] = 1.0;

	// Sets the orientation of the liquid crystal network alignment in the r, theta coordinate system
	LCN_Alignment[0] = 1.0;


	DirectionTens = outer_product(LCN_Alignment,LCN_Alignment);
	DirectionTens = DirectionTens/TensorNorm2D(DirectionTens);




	Tensor<2,DIM> As_sc;

	Tensor<2,DIM> Ab_sc;


	stab_magnitude = 000.;


	std::vector<double> ratio1 = linspace(0.0,1.0,Nsteps1);
	std::vector<double> ratioInit = linspace(0.0,1.0,NstepsInit);



	int update_counter = 0;
	int stage;

	stage = 0;
	As_Global = 0.0*DirectionTens;
	Ab_Global = 0.0*DirectionTens;



	double internal_drag_init = internal_drag;
	evaluation_point = solution;
	assemble_nonlinear_system(2);

	double scaling_factor = 1.0;


	Update_Internal_Variables(stage);
	solve_linear_system();


	for (unsigned int i = 0; i < NstepsInit; i++){
		std::cout << "Stage: " << stage << ", STARTING STEP: " << i << "   --------------------------------------- " << std::endl;
		//As_Global = ratio1[i]*Light_Magnitude_As*DirectionTens;
		std::cout << "Epsilon Value: " << ratioInit[i]*epsilonmax << std::endl;
		std::cout << "Kappa Value: " << ratioInit[i]*kappaneggoal << std::endl;

		epsilon_Global = ratioInit[i]*epsilonmax;
		kappa_Global = (ratioInit[i]+0.05*(1.0 - ratioInit[i]))*kappaneggoal;

		Update_Internal_Variables(stage);
		solve_linear_system();

		output_results(update_counter);
		update_counter++;
		std::cout << "Updating Internal Variables" << std::endl;


	}
	double dkappa = 0.02*L;
	double current_kappa = kappaneggoal;
	Vector<double> snapped_solution;
	FlipConfigUD(snapped_solution);
	stage = 1;
	Vector<double> last_solution = solution;

	//First I need to solve for the first snapped through solution.
	solution = last_solution;


	epsilon_Global = epsilonmax;

	bool did_snapthrough_happen = false;
	Vector<double> last_iteration_solution = solution;
	double change_in_solution_norm = 1.0;
	bool firstsolve = true;
	while (!did_snapthrough_happen) {
		std::cout << "Stage: " << stage << ", STARTING STEP: " << 0 << "   --------------------------------------- " << std::endl;
		//As_Global = ratio1[i]*Light_Magnitude_As*DirectionTens;
		current_kappa += dkappa;
		std::cout << "Epsilon Value: " << epsilonmax << std::endl;
		std::cout << "Kappa Value: " << current_kappa << std::endl;
		std::cout << "Counter: " << update_counter << std::endl;


		kappa_Global = current_kappa;
		Update_Internal_Variables(stage);
		last_iteration_solution = solution;
		int cmax = 20;
		int solvenum = solve_linear_system_critical(cmax);
		double last_change_in_solution_norm = change_in_solution_norm;
		change_in_solution_norm = 0.0;
		for (unsigned int k = 0; k < solution.size(); k++){
			change_in_solution_norm += pow(solution(k) - last_iteration_solution(k),2.0);
		}
		change_in_solution_norm = sqrt(change_in_solution_norm)/solution.size();



		double relative_solution_change = change_in_solution_norm/last_change_in_solution_norm;
		std::cout << "Relative Change in Solution: " << relative_solution_change << std::endl;

		if ((!firstsolve && (relative_solution_change > 20.0)) || cmax == solvenum) {
			did_snapthrough_happen = true;
		}
		output_results(update_counter);
		update_counter++;
		std::cout << "Updating Internal Variables" << std::endl;
		firstsolve = false;

	}
	solution = snapped_solution;
	int solvenum2 = solve_linear_system_critical(20);
	snapped_solution = solution;


	for (unsigned int i = 1; i < 400; i++){
		std::cout << "Stage: " << stage << ", STARTING STEP: " << 0 << "   --------------------------------------- " << std::endl;
		//As_Global = ratio1[i]*Light_Magnitude_As*DirectionTens;
		current_kappa += dkappa;
		std::cout << "Epsilon Value: " << epsilonmax << std::endl;
		std::cout << "Kappa Value: " << current_kappa << std::endl;
		std::cout << "Counter: " << update_counter << std::endl;


		kappa_Global = current_kappa;
		Update_Internal_Variables(stage);
		last_iteration_solution = solution;
		int cmax = 20;
		int solvenum = solve_linear_system_critical(cmax);
		double last_change_in_solution_norm = change_in_solution_norm;
		change_in_solution_norm = 0.0;
		for (unsigned int k = 0; k < solution.size(); k++){
			change_in_solution_norm += pow(solution(k) - last_iteration_solution(k),2.0);
		}
		change_in_solution_norm = sqrt(change_in_solution_norm)/solution.size();



		double relative_solution_change = change_in_solution_norm/last_change_in_solution_norm;
		std::cout << "Relative Change in Solution: " << relative_solution_change << std::endl;

		if ((!firstsolve && (relative_solution_change > 20.0)) || cmax == solvenum) {
			did_snapthrough_happen = true;
		}
		output_results(update_counter);
		update_counter++;
		std::cout << "Updating Internal Variables" << std::endl;
		firstsolve = false;
	}

}


int ElasticProblem::solve_linear_system_critical(int cmax){

	resnorm = 2.0*tol;
	double stepsize = 2.0*tol;
	unsigned int cntr = 0;


	for (unsigned int i = 0; i < evaluation_point.size(); i++) {
		evaluation_point(i) = solution(i) + perturbation_magnitude*perturbation_direction(i);
	}
	solution_old = solution;

	bool drag_on = false;
	bool resolve = false;
	double internal_drag_init = internal_drag;
	SparseDirectUMFPACK a_direct;
	internal_drag = 0.0;
	int cntrlim = 20;
	while (resnorm > tol && stepsize > tol && cntr < cmax){
		cntr++;
		prev_solution = solution;

		if (cntr >= cntrlim ) {

			if (cntr == cntrlim) {
				//alpha = 0.9;
				internal_drag = 100.*((double) (cntr/10)) * internal_drag_init;
				solution = solution_old;
				evaluation_point = solution;

			}

			solution_old = solution;

			//resolve = true;
		}

		std::cout << "Assembling System" << std::endl;


		assemble_nonlinear_system(2); //Initializes the system matrix and right hand side vector
		std::cout << "Solving Linear System" << std::endl;

		a_direct.initialize(system_matrix);
		step_direction = 0.0;

		a_direct.vmult(step_direction,system_rhs);
		stepsize = sqrt(step_direction.norm_sqr())/step_direction.size();
		//residual = step_direction;
		constraints.distribute(step_direction);

		std::cout << "Iteration: " << cntr << std::endl;
		std::cout << "Step Size: " << stepsize<< std::endl;

		solution.add(-1.0*alpha,step_direction);

		evaluation_point = solution;
		assemble_nonlinear_system(1);
		resnorm = sqrt(residual.norm_sqr())/residual.size();
		std::cout << "Residual: " << resnorm << std::endl;


		feasibilitycalc();
		if (rp/rd > tau){
			muval = kd*muval;
		} else if (rd/rp > tau){
			muval = std::max(muval/kd,mumin);
		} else {

		}
		/*
		if (resnorm > 5.0) {
			drag_on = true;
			std::cout << "Adding Drag" << std::endl;
			solution = solution_old;
			evaluation_point = solution_old;

			drag_coeff = 100.0*internal_drag;

		} else if (resnorm < 0.0010 && drag_on) {
			std::cout << "Removing Drag" << std::endl;
			solution_old = solution;
			drag_coeff = 0.0;
			drag_on = false;
		}
		 */
	}
	internal_drag = internal_drag_init;
	alpha = 1.0;

	if (resolve) {
		solve_linear_system();
	}

	return cntr;
}


void ElasticProblem::As_function(const std::vector<Point<DIM>> &points,
		std::vector<Tensor<2, DIM>> &  values){
	Tensor<2,DIM> zerotensor;


	Tensor<1,DIM> n1;
	n1[0] = -1.0/sqrt(2);
	n1[1] = 1.0/sqrt(2);

	Tensor<1,DIM> n2;
	n2[0] = 1.0/sqrt(2);
	n2[1] = 1.0/sqrt(2);
	bool c1 = false;
	bool c2 = false;
	double theta = 3.14159/2.0;
	Tensor<2,DIM> R1;
	R1[0][0] = cos(theta);
	R1[0][1] = -1.0*sin(theta);
	R1[1][0] = sin(theta);
	R1[1][1] = cos(theta);
	Tensor<2,DIM> R2;
	R2[0][0] = cos(theta);
	R2[0][1] = sin(theta);
	R2[1][0] = -1.0*sin(theta);
	R2[1][1] = cos(theta);
	for (unsigned int i = 0; i < values.size(); i++){
		Tensor<1,DIM> x;
		x[0] = points[i][0];
		x[1] = points[i][1];


		if (x*n1 >= 0.0) {
			c1 = true;
		} else {
			c1 = false;
		}

		if (x*n2 >= 0.0) {
			c2 = true;
		} else {
			c2 = false;
		}

		if (c1) {
			if (c2) {
				// Upper section
				values[i] = As_Global;
			} else {
				// Left section
				values[i] = R1*As_Global*R2;
			}
		} else {
			if (c2) {
				// right section
				values[i] = R1*As_Global*R2;
			} else {
				// bottom section
				values[i] = As_Global;
			}
		}



	}



}

void ElasticProblem::Ab_function(const std::vector<Point<DIM>> &points,
		std::vector<Tensor<2, DIM>> &  values){




	for (unsigned int i = 0; i < values.size(); i++){
		values[i] = Ab_Global;
	}

}

void ElasticProblem::Initialize_Internal_Variables() {


	QGauss<DIM> quadrature_formula(fe.degree + 3);

	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);


	const unsigned int n_q_points    = quadrature_formula.size();


	Ab_Internal.resize(triangulation.n_active_cells());
	Ab0_Internal.resize(triangulation.n_active_cells());
	As_Internal.resize(triangulation.n_active_cells());
	As0_Internal.resize(triangulation.n_active_cells());
	std::vector<Point<DIM>> quadrature_points(n_q_points);
	quadrature_points = fe_values.get_quadrature_points();


	Tensor<2,2> zerotensor;
	zerotensor[0][0] = 0.0;
	zerotensor[1][0] = 0.0;
	zerotensor[0][1] = 0.0;
	zerotensor[1][1] = 0.0;
	double pi = 3.14159;
	Tensor<2,DIM> Abinit;
	// Now we can begin with the loop over all cells:
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		unsigned int cell_index = cell->active_cell_index();
		fe_values.reinit(cell);
		Ab_Internal[cell_index].resize(n_q_points);
		As_Internal[cell_index].resize(n_q_points);
		Ab0_Internal[cell_index].resize(n_q_points);
		As0_Internal[cell_index].resize(n_q_points);


		for (const unsigned int q_point : fe_values.quadrature_point_indices()){

			Ab_Internal[cell_index][q_point] = zerotensor;

			As_Internal[cell_index][q_point] = zerotensor;
			As0_Internal[cell_index][q_point] = zerotensor;

			double defmag = 0.1;
			Abinit[0][0] = -0.0*defmag*sin(pi*quadrature_points[q_point][0]/L);
			Abinit[0][1] = 0.0;
			Abinit[1][0] = 0.0;
			Abinit[1][1] = 0.0*defmag*sin(pi*quadrature_points[q_point][1]/L);



		}


	}


}

void ElasticProblem::Update_Internal_Variables(const int stage) {

	QGauss<DIM> quadrature_formula(fe.degree + 3);

	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();


	std::vector<Tensor<1,DIM>> u_q(n_q_points);
	std::vector<double> u3_q(n_q_points);
	std::vector<Tensor<2,DIM>> grad_u_q(n_q_points);
	std::vector<Tensor<1,DIM>> grad_u3_q(n_q_points);
	std::vector<Tensor<1,DIM>> w_q(n_q_points);
	std::vector<Tensor<2,DIM>> grad_w_q(n_q_points);
	std::vector<Tensor<1,DIM>> lagrange_multiplier_q(n_q_points);
	std::vector<Tensor<1,DIM>> u_old_q(n_q_points);
	std::vector<Point<DIM>> quadrature_points(n_q_points);
	std::vector<double> u3_old_q(n_q_points);
	std::vector<Tensor<2,DIM>> As_values(n_q_points);

	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);



	// Now we can begin with the loop over all cells:
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		unsigned int cell_index = cell->active_cell_index();


		const FEValuesExtractors::Vector u(0);
		const FEValuesExtractors::Scalar u3(DIM);
		const FEValuesExtractors::Vector w(DIM+1);
		const FEValuesExtractors::Vector lagrange_multiplier(DIM+3);


		fe_values.reinit(cell);



		fe_values[u].get_function_values(evaluation_point,u_q);
		fe_values[u3].get_function_values(evaluation_point,u3_q);
		fe_values[u].get_function_gradients(evaluation_point,grad_u_q);
		fe_values[u3].get_function_gradients(evaluation_point,grad_u3_q);
		fe_values[w].get_function_values(evaluation_point,w_q);
		fe_values[w].get_function_gradients(evaluation_point,grad_w_q);
		fe_values[lagrange_multiplier].get_function_values(evaluation_point,lagrange_multiplier_q);
		fe_values[u].get_function_values(solution_old,u_old_q);
		fe_values[u3].get_function_values(solution_old,u3_old_q);
		quadrature_points = fe_values.get_quadrature_points();

		As_function(fe_values.get_quadrature_points(),As_values);



		for (const unsigned int q_point : fe_values.quadrature_point_indices()){
			if (stage == 0 || stage == 1) {



				Tensor<1,DIM> x;
				x[0] = quadrature_points[q_point][0] - offset[0];
				x[1] = quadrature_points[q_point][1] - offset[1];

				double xnorm = sqrt(pow(x[0],2.0) + pow(x[1],2.0));

				if (xnorm > 10e-6) {
					Tensor<1,DIM> er = x/xnorm;
					Tensor<1,DIM> etheta;
					etheta[0] = -er[1];
					etheta[1] = er[0];

					Tensor<2,DIM> DirectionTensor = outer_product(etheta, etheta);

					As_Internal[cell_index][q_point] = epsilon_Global * DirectionTensor;
					Ab_Internal[cell_index][q_point] = kappa_Global * DirectionTensor;
				} else {
					Tensor<2,DIM> IdentityTens;
					IdentityTens[0][0] = 1.0;
					IdentityTens[1][1] = 1.0;
					As_Internal[cell_index][q_point] = epsilon_Global * IdentityTens;
					Ab_Internal[cell_index][q_point] = kappa_Global * IdentityTens;
				}





			}





			// Updating Ab

		}



	}

}

Tensor<1,DIM> ElasticProblem::basis_vector(const int i) {
	Tensor<1,DIM> vout;
	if (i < DIM) {
		vout[i] = 1.0;
	}
	return vout;

}

bool ElasticProblem::CharacteristicFunction(const Point<DIM> & point){
	return true;
}
Tensor<2,DIM> ElasticProblem::outer_product(const Tensor<1,DIM> & v1, const Tensor<1,DIM> & v2) {
	Tensor<2,DIM> Tout;
	for (unsigned int i = 0; i < DIM; i++) {
		for (unsigned int j = 0; j < DIM; j++) {
			Tout[i][j] = v1[i]*v2[j];
		}
	}
	return Tout;
}



double ElasticProblem::Tensor_Inner(const Tensor<2,DIM> & H1, const Tensor<2,DIM> & H2) {
	double sum = 0.0;
	for (unsigned int i = 0; i < DIM; i++) {
		for (unsigned int j = 0; j < DIM; j++) {
			sum += H1[i][j]*H2[i][j];
		}
	}
	return sum;
}

double ElasticProblem::Det2D(const Tensor<2,DIM> & H) {
	return H[0][0]*H[DIM-1][DIM-1] - H[0][DIM-1]*H[DIM-1][0];
}


Tensor<2,DIM> ElasticProblem::Tensor_Transpose(const Tensor<2,DIM> & H) {
	Tensor<2,DIM> Tout;
	for (unsigned int i = 0; i < DIM; i++) {
		for (unsigned int j = 0; j < DIM; j++) {
			Tout[j][i] = H[i][j];
		}
	}
	return Tout;
}


Tensor<2,DIM> ElasticProblem::cofactor2D(const Tensor<2,DIM> & F) {
	Tensor<2,DIM> CofacF;
	CofacF[0][0] = F[DIM-1][DIM-1];
	CofacF[0][DIM-1] = -1.0*F[DIM-1][0];
	CofacF[DIM-1][0] = -1.0*F[0][DIM-1];
	CofacF[DIM-1][DIM-1] = F[0][0];
	return CofacF;
}


double ElasticProblem::BilinearProduct(const Tensor<2,DIM> & F1, const Tensor<4,DIM> & C, const Tensor<2,DIM> & F2) {
	double sum = 0.0;

	for (unsigned int i = 0; i < DIM; i++) {
		for (unsigned int j = 0; j < DIM; j++) {
			for (unsigned int k = 0; k < DIM; k++) {
				for (unsigned int l = 0; l < DIM; l++) {
					sum += F1[i][j]*C[i][j][k][l]*F2[k][l];
				}
			}
		}
	}
	return sum;
}


double ElasticProblem::Tensor_Trace(const Tensor<2,DIM> & F) {
	double sum = 0.0;
	for (unsigned int i = 0; i < DIM; i++){
		sum += F[i][i];
	}
	return sum;
}


Tensor<2,DIM> ElasticProblem::Tensor_Sym(const Tensor<2,DIM> & F){
	return 0.5*(F + Tensor_Transpose(F));
}

Tensor<1,DIM+1> ElasticProblem::crossproduct(const Tensor<1,DIM+1> & v1,const Tensor<1,DIM+1> & v2) {
	Tensor<1,DIM+1> vout;

	vout[0] = v1[1]*v2[2] - v1[2]*v2[1];
	vout[1] = v1[2]*v2[0] - v1[0]*v2[2];
	vout[2] = v1[0]*v2[1] - v1[1]*v2[0];
	return vout;
}
double ElasticProblem::TensorNorm1D(const Tensor<1,DIM+1> & v){
	return sqrt(pow(v[0],2.0)+pow(v[1],2.0)+pow(v[2],2.0));
}

double ElasticProblem::TensorNorm2D(const Tensor<2,DIM> & H) {
	double sum = 0.0;

	for (unsigned int i = 0; i < DIM; i++) {
		for (unsigned int j = 0; j < DIM; j++) {
			sum += pow(H[i][j],2.0);
		}
	}
	return sqrt(sum);
}

void ElasticProblem::solution_stabilization(const std::vector<Point<DIM>> &points,
		std::vector<Tensor<1, DIM+1>> &  values)
{
	double max_radius = .1;




	for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
	{



		if ((points[point_n]).norm_square() < max_radius * max_radius ) {
			values[point_n][0] = stab_magnitude;;
			values[point_n][1] = stab_magnitude;;
			values[point_n][2] = stab_magnitude;;
		} else {
			values[point_n][0] = 0.0;
			values[point_n][1] = 0.0;
			values[point_n][2] = 0.0;
		}



	}


}


void ElasticProblem::right_hand_side(const std::vector<Point<DIM>> &points,
		std::vector<Tensor<1, DIM+5>> &  values)
{

	Point<DIM> point_1, point_2;
	point_1(0) = 0.5;
	point_2(0) = -0.5;

	for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
	{


		if (((points[point_n] - point_1).norm_square() < 0.2 * 0.2) ||
				((points[point_n] - point_2).norm_square() < 0.2 * 0.2))
			values[point_n][0] = 0.0*1.0;
		else
			values[point_n][0] = 0.0;


		if (points[point_n].norm_square() < 0.2 * 0.2)
			values[point_n][1] = 0.0*1.0;
		else
			values[point_n][1] = 0.0;

	}

	for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
	{


		for (unsigned int i = 0; i < 5; i++){
			values[point_n][i] = 0.0;
		}

	}
}

void ElasticProblem::Youngs_Modulus(const std::vector<Point<DIM>> &points,
		std::vector<double> &  values){

	for (unsigned int point_n = 0; point_n < points.size(); ++point_n)
	{
		if (CharacteristicFunction(points[point_n])) {
			values[point_n] = Eval;
		} else {
			values[point_n] = 10.0*Eval;
		}
	}

}

std::vector<double> ElasticProblem::linspace(double start_in, double end_in, int num_in)
{

	std::vector<double> linspaced;

	double start = static_cast<double>(start_in);
	double end = static_cast<double>(end_in);
	double num = static_cast<double>(num_in);

	if (num == 0) { return linspaced; }
	if (num == 1)
	{
		linspaced.push_back(start);
		return linspaced;
	}

	double delta = (end - start) / (num - 1);

	for(int i=0; i < num-1; ++i)
	{
		linspaced.push_back(start + delta * i);
	}
	linspaced.push_back(end); // I want to ensure that start and end
	// are exactly the same as the input
	return linspaced;
}


void ElasticProblem::OutputVector(std::vector<double> & v1, std::vector<double> & v2, std::string name) {
	std::ofstream myfile;
	myfile.open(name);
	myfile << "x,y\n";
	for (unsigned int i = 0; i < v1.size(); i++) {
		myfile << v1[i] << "," << v2[i] << "\n";
	}
	myfile.close();
}


ElasticProblem::ElasticProblem()
: dof_handler(triangulation)
, fe(FESystem<DIM>(FE_Q<DIM>(2),2), 1,FE_Q<DIM>(1), 1,FESystem<DIM>(FE_Q<DIM>(1),2), 1,FESystem<DIM>(FE_Q<DIM>(1),2), 1)
{}


void ElasticProblem::setup_system()
{
	dof_handler.distribute_dofs(fe);
	number_dofs = dof_handler.n_dofs();
	solution.reinit(number_dofs);
	solution_old.reinit(number_dofs);
	residual.reinit(number_dofs);
	evaluation_point.reinit(number_dofs);
	perturbation_direction.reinit(number_dofs);
	Material_Vector_InPlane.resize(triangulation.n_active_cells());
	Material_Vector_Bending.resize(triangulation.n_active_cells());


	system_rhs.reinit(number_dofs);

	setup_constraints();

	VectorTools::interpolate(dof_handler,PerturbationValues(), perturbation_direction);
	//VectorTools::interpolate(dof_handler,Functions::ZeroFunction<DIM>(DIM+5), perturbation_direction);



	DynamicSparsityPattern dsp(number_dofs, number_dofs);
	DoFTools::make_sparsity_pattern(dof_handler,
			dsp,
			constraints,
			/*keep_constrained_dofs = */ false);
	sparsity_pattern.copy_from(dsp);

	system_matrix.reinit(sparsity_pattern);

	std::ofstream out("sparsity_pattern2.svg");
	sparsity_pattern.print_svg(out);


	constraints.distribute(perturbation_direction);
	Initialize_Internal_Variables();


}





void ElasticProblem::assemble_nonlinear_system(const unsigned int solve_level)
{
	system_matrix = 0;
	system_rhs    = 0;
	QGauss<DIM> quadrature_formula(fe.degree + 3);

	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);

	std::vector<Tensor<1,DIM>> u_q(n_q_points);
	std::vector<double> u3_q(n_q_points);
	std::vector<Tensor<2,DIM>> grad_u_q(n_q_points);
	std::vector<Tensor<1,DIM>> grad_u3_q(n_q_points);
	std::vector<Tensor<1,DIM>> w_q(n_q_points);
	std::vector<Tensor<2,DIM>> grad_w_q(n_q_points);
	std::vector<Tensor<1,DIM>> lagrange_multiplier_q(n_q_points);
	std::vector<Tensor<1,DIM>> u_old_q(n_q_points);
	std::vector<double> u3_old_q(n_q_points);


	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


	std::vector<double> lambda_values(n_q_points);
	std::vector<double> mu_values(n_q_points);
	std::vector<double> nu_values(n_q_points);
	std::vector<double> Emod_values(n_q_points);


	Functions::ConstantFunction<DIM> lambda(1.), mu(1.), nu(nu_global);


	std::vector<Tensor<1, DIM+5>> rhs_values(n_q_points);

	std::vector<Tensor<1,DIM+1>> stab_values(n_q_points);



	StrainEnergy = 0.0;

	// Now we can begin with the loop over all cells:
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		unsigned int cell_index = cell->active_cell_index();


		const FEValuesExtractors::Vector u(0);
		const FEValuesExtractors::Scalar u3(DIM);
		const FEValuesExtractors::Vector w(DIM+1);
		const FEValuesExtractors::Vector lagrange_multiplier(DIM+3);

		cell_matrix = 0;
		cell_rhs    = 0;

		fe_values.reinit(cell);

		// Next we get the values of the coefficients at the quadrature
		// points. Likewise for the right hand side:
		lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
		mu.value_list(fe_values.get_quadrature_points(), mu_values);
		nu.value_list(fe_values.get_quadrature_points(), nu_values);
		right_hand_side(fe_values.get_quadrature_points(), rhs_values);
		solution_stabilization(fe_values.get_quadrature_points(), stab_values);
		Youngs_Modulus(fe_values.get_quadrature_points(), Emod_values);



		fe_values[u].get_function_values(evaluation_point,u_q);
		fe_values[u3].get_function_values(evaluation_point,u3_q);
		fe_values[u].get_function_gradients(evaluation_point,grad_u_q);
		fe_values[u3].get_function_gradients(evaluation_point,grad_u3_q);
		fe_values[w].get_function_values(evaluation_point,w_q);
		fe_values[w].get_function_gradients(evaluation_point,grad_w_q);
		fe_values[lagrange_multiplier].get_function_values(evaluation_point,lagrange_multiplier_q);
		fe_values[u].get_function_values(solution_old,u_old_q);
		fe_values[u3].get_function_values(solution_old,u3_old_q);



		for (const unsigned int q_point :
				fe_values.quadrature_point_indices())
		{
			const Tensor<2,DIM> In_Plane_Strain_q = Tensor_Sym(grad_u_q[q_point]) + 0.5*outer_product(grad_u3_q[q_point],grad_u3_q[q_point])- As_Internal[cell_index][q_point];
			const Tensor<2,DIM> Bending_Strain_q = grad_w_q[q_point] - Ab_Internal[cell_index][q_point];

			const double nu_loc = nu_values[q_point];
			//const Tensor<2,DIM> grad_u = fe_values[u].get_function_gradients(q_point);
			//const Tensor<2,DIM> sym_grad_u = 0.5*(grad_u + Tensor_Transpose<DIM>(grad_u));
			Material_Vector_InPlane[cell_index].set_Params(Emod_values[q_point], 0.49, In_Plane_Strain_q);
			Material_Vector_Bending[cell_index].set_Params(Emod_values[q_point], 0.49, Bending_Strain_q);
			//Material_Class Material_InPlane(Eval, nu_loc, In_Plane_Strain_q);
			//Material_Class Material_Bending(Eval, nu_loc, Bending_Strain_q);


			StrainEnergy += (Material_Vector_InPlane[cell_index].getQ2() + hsc*Material_Vector_Bending[cell_index].getQ2())*fe_values.JxW(q_point);

			for (const unsigned int i : fe_values.dof_indices())
			{

				const Tensor<1,DIM> U_i_q = fe_values[u].value(i,q_point);
				const double U3_i_q = fe_values[u3].value(i,q_point);
				//const double div_Phi_i_u = fe_values[u].divergence(i,q_point);
				const Tensor<2,DIM> grad_U_i_q = fe_values[u].gradient(i,q_point);
				const Tensor<2,DIM> symgrad_U_i_q = Tensor_Sym(grad_U_i_q);
				//const Tensor<2,DIM> cofac_grad_Phi_i_u = cofactor2D<DIM>(symgrad_Phi_i_u);

				const Tensor<1,DIM> grad_U3_i_q = fe_values[u3].gradient(i,q_point);

				const Tensor<1,DIM> W_i_q = fe_values[w].value(i,q_point);
				//const double div_w_i_u = fe_values[w].divergence(i,q_point);
				const Tensor<2,DIM> grad_W_i_q = fe_values[w].gradient(i,q_point);
				//const Tensor<2,DIM> cofac_grad_w_i_u = cofactor2D<DIM>(grad_w_i_u);
				const Tensor<1,DIM> Lagrange_multiplier_i_q = fe_values[lagrange_multiplier].value(i,q_point);

				const Tensor<2,DIM> dsymwotimesW_i_q = Tensor_Sym(outer_product(W_i_q,w_q[q_point]));
				const Tensor<2,DIM> dsym_grad_u3_otimes_grad_U3_i_q = Tensor_Sym(outer_product(grad_U3_i_q,grad_u3_q[q_point]));

				if (solve_level >= 2){
					for (const unsigned int j : fe_values.dof_indices())
					{
						const Tensor<1,DIM> U_j_q = fe_values[u].value(j,q_point);
						const double U3_j_q = fe_values[u3].value(j,q_point);


						//const double div_Phi_j_u = fe_values[u].divergence(j,q_point);
						const Tensor<2,DIM> grad_U_j_q = fe_values[u].gradient(j,q_point);
						const Tensor<2,DIM> symgrad_U_j_q = Tensor_Sym(grad_U_j_q);
						//const Tensor<2,DIM> cofac_grad_Phi_j_u = cofactor2D<DIM>(symgrad_Phi_j_u);

						const Tensor<1,DIM> grad_U3_j_q = fe_values[u3].gradient(j,q_point);
						const Tensor<1,DIM> W_j_q = fe_values[w].value(j,q_point);
						//const double div_w_j_u = fe_values[w].divergence(j,q_point);
						const Tensor<2,DIM> grad_W_j_q = fe_values[w].gradient(j,q_point);

						const Tensor<1,DIM> Lagrange_multiplier_j_q = fe_values[lagrange_multiplier].value(j,q_point);

						const Tensor<2,DIM> dsymwotimesW_j_q = Tensor_Sym(outer_product(W_j_q,w_q[q_point]));
						const Tensor<2,DIM> dsym_grad_u3_otimes_grad_U3_j_q = Tensor_Sym(outer_product(grad_U3_j_q,grad_u3_q[q_point]));

						const Tensor<2,DIM> ddsymWotimesW_ij_q = Tensor_Sym(outer_product(W_i_q,W_j_q));
						const Tensor<2,DIM> ddsym_grad_U3_otimes_grad_U3_ij_q = Tensor_Sym(outer_product(grad_U3_i_q,grad_U3_j_q));
						//double directsolve = (div_Phi_i_u*div_Phi_j_u - (1.0-nu_values[q_point]) * Tensor_Inner<DIM>(cofac_grad_Phi_i_u,symgrad_Phi_j_u) 
						//	)*fe_values.JxW(q_point); 


						//cell_matrix(i,j) += (w_i_u*w_j_u)*fe_values.JxW(q_point); 

						//cell_matrix(i,j) += (grad_phi_i_u3*grad_phi_j_u3)*fe_values.JxW(q_point); 

						//cell_matrix(i,j) += (Tensor_Inner<DIM>(Material_InPlane.getdQ2dF(),symgrad_Phi_j_u))*fe_values.JxW(q_point);

						// Internal strain energies
						cell_matrix(i,j) += hscinv*BilinearProduct(symgrad_U_i_q+dsym_grad_u3_otimes_grad_U3_i_q,Material_Vector_InPlane[cell_index].getddQ2ddF(),symgrad_U_j_q + dsym_grad_u3_otimes_grad_U3_j_q)*fe_values.JxW(q_point);
						cell_matrix(i,j) += hscinv*Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),ddsym_grad_U3_otimes_grad_U3_ij_q)*fe_values.JxW(q_point);
						cell_matrix(i,j) += BilinearProduct(grad_W_i_q,Material_Vector_Bending[cell_index].getddQ2ddF(),grad_W_j_q)*fe_values.JxW(q_point);

						//Constraint satisfaction

						cell_matrix(i,j) -= (Lagrange_multiplier_i_q*(W_j_q - grad_U3_j_q))*fe_values.JxW(q_point);

						cell_matrix(i,j) -= (Lagrange_multiplier_j_q)*(W_i_q - grad_U3_i_q)*fe_values.JxW(q_point);

						cell_matrix(i,j) += muval*((W_i_q - grad_U3_i_q)*(W_j_q - grad_U3_j_q))*fe_values.JxW(q_point);

						// Making the hessian matrix invertible
						cell_matrix(i,j) += internal_drag*(U_i_q*U_j_q + U3_i_q*U3_j_q)*fe_values.JxW(q_point);

						cell_matrix(i,j) += (U_i_q[0] * stab_values[q_point][0] * U_j_q[0] + U_i_q[1] * stab_values[q_point][1] * U_j_q[1] + U3_i_q * stab_values[q_point][2] * U3_j_q)*fe_values.JxW(q_point);

						cell_matrix(i,j) += drag_coeff*(U_i_q*U_j_q + U3_i_q*U3_j_q)*fe_values.JxW(q_point);

						//cell_matrix(i,j) += hsc*(div_w_i_u*div_w_j_u - (1.0-nu_values[q_point]) * Tensor_Inner<DIM>(cofac_grad_w_i_u,grad_w_j_u) 
						//	)*fe_values.JxW(q_point); 


					}
				}




				const unsigned int component_i = fe.system_to_component_index(i).first;

				Tensor<1,DIM> wforce;
				wforce[0] = 0.0;
				wforce[1] = 0.0;

				cell_rhs(i) += (W_i_q*wforce)*fe_values.JxW(q_point);
				cell_rhs(i) += fe_values.shape_value(i, q_point) *
						rhs_values[q_point][component_i] *
						fe_values.JxW(q_point);

				//



				cell_rhs(i) += hscinv*(Tensor_Inner(Material_Vector_InPlane[cell_index].getdQ2dF(),symgrad_U_i_q + dsym_grad_u3_otimes_grad_U3_i_q))*fe_values.JxW(q_point);
				cell_rhs(i) += (Tensor_Inner(Material_Vector_Bending[cell_index].getdQ2dF(),grad_W_i_q))*fe_values.JxW(q_point);
				cell_rhs(i) -= Lagrange_multiplier_i_q*(w_q[q_point] - grad_u3_q[q_point])*fe_values.JxW(q_point);
				cell_rhs(i) -= lagrange_multiplier_q[q_point]*(W_i_q - grad_U3_i_q)*fe_values.JxW(q_point);
				cell_rhs(i) += muval*((w_q[q_point] - grad_u3_q[q_point])*(W_i_q - grad_U3_i_q))*fe_values.JxW(q_point);
				cell_rhs(i) += internal_drag*((U_i_q*(u_q[q_point] - u_old_q[q_point])) + U3_i_q*(u3_q[q_point] - u3_old_q[q_point]))*fe_values.JxW(q_point);
				cell_rhs(i) += (U_i_q[0] * stab_values[q_point][0] *u_q[q_point][0] + U_i_q[1] * stab_values[q_point][1] *u_q[q_point][1] + U3_i_q* stab_values[q_point][2] *u3_q[q_point])*fe_values.JxW(q_point);

				cell_rhs(i) += drag_coeff*((u_q[q_point] - u_old_q[q_point])*U_i_q + (u3_q[q_point] - u3_old_q[q_point])*U3_i_q)*fe_values.JxW(q_point);
			}
		}



		// The transfer from local degrees of freedom into the global matrix
		// and right hand side vector does not depend on the equation under
		// consideration, and is thus the same as in all previous
		// examples.
		cell->get_dof_indices(local_dof_indices);
		for (unsigned int i = 0; i < local_dof_indices.size(); i++) {

			system_rhs[local_dof_indices[i]] += cell_rhs[i];
			for (unsigned int j = 0; j < local_dof_indices.size(); j++) {
				system_matrix.add(local_dof_indices[i],local_dof_indices[j],cell_matrix[i][j]);
			}
		}

		//constraints.distribute_local_to_global(
		//		cell_matrix, cell_rhs, local_dof_indices, system_matrix, system_rhs);

	}
	constraints.condense(system_matrix);
	constraints.condense(system_rhs);
	residual = system_rhs;
}


// @sect4{ElasticProblem::solve}



void ElasticProblem::solve()
{
	/*
	SolverControl            solver_control(1000, 1e-12);
	SolverCG<Vector<double>> cg(solver_control);

	PreconditionSSOR<SparseMatrix<double>> preconditioner;
	preconditioner.initialize(system_matrix, 1.2);

	cg.solve(system_matrix, solution, system_rhs, preconditioner);

	constraints.distribute(solution);
	 */
	SparseDirectUMFPACK a_direct;
	a_direct.initialize(system_matrix);
	a_direct.vmult(solution,system_rhs);
	constraints.distribute(solution);
	std::cout << solution.norm_sqr()<< std::endl;
}


double ElasticProblem::get_Strain_Energy()
{

	evaluation_point = solution;
	assemble_nonlinear_system(0);
	return StrainEnergy;

}





void ElasticProblem::solve_linear_system()
{
	resnorm = 2.0*tol;
	double stepsize = 2.0*tol;
	unsigned int cntr = 0;


	for (unsigned int i = 0; i < evaluation_point.size(); i++) {
		evaluation_point(i) = solution(i) + perturbation_magnitude*perturbation_direction(i);
	}
	solution_old = solution;

	bool drag_on = false;
	bool resolve = false;
	double internal_drag_init = internal_drag;
	SparseDirectUMFPACK a_direct;

	while (resnorm > tol && stepsize > tol){
		cntr++;
		prev_solution = solution;

		if (cntr >= 20 ) {

			if (cntr == 20) {
				alpha = 0.9;
				internal_drag = 20.*((double) (cntr/10)) * internal_drag_init;
				solution = solution_old;
				evaluation_point = solution;

			}

			solution_old = solution;

			resolve = true;
		}

		std::cout << "Assembling System" << std::endl;


		assemble_nonlinear_system(2); //Initializes the system matrix and right hand side vector
		std::cout << "Solving Linear System" << std::endl;

		a_direct.initialize(system_matrix);
		step_direction = 0.0;

		a_direct.vmult(step_direction,system_rhs);
		stepsize = sqrt(step_direction.norm_sqr())/step_direction.size();
		//residual = step_direction;
		constraints.distribute(step_direction);

		std::cout << "Iteration: " << cntr << std::endl;
		std::cout << "Step Size: " << stepsize<< std::endl;

		solution.add(-1.0*alpha,step_direction);

		evaluation_point = solution;
		assemble_nonlinear_system(1);
		resnorm = sqrt(residual.norm_sqr())/residual.size();
		std::cout << "Residual: " << resnorm << std::endl;
		feasibilitycalc();
		if (rp/rd > tau){

			muval = kd*muval;
		} else if (rd/rp > tau){
			muval = std::max(muval/kd,mumin);
		}
		/*
		if (resnorm > 5.0) {
			drag_on = true;
			std::cout << "Adding Drag" << std::endl;
			solution = solution_old;
			evaluation_point = solution_old;

			drag_coeff = 100.0*internal_drag;

		} else if (resnorm < 0.0010 && drag_on) {
			std::cout << "Removing Drag" << std::endl;
			solution_old = solution;
			drag_coeff = 0.0;
			drag_on = false;
		}
		 */
	}
	internal_drag = internal_drag_init;
	alpha = 1.0;

	if (resolve) {
		solve_linear_system();
	}
}
void ElasticProblem::feasibilitycalc(){

	QGauss<DIM> quadrature_formula(fe.degree + 3);

	FEValues<DIM> fe_values(fe,
			quadrature_formula,
			update_values | update_gradients |
			update_quadrature_points | update_JxW_values);

	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	const unsigned int n_q_points    = quadrature_formula.size();

	FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
	Vector<double>     cell_rhs(dofs_per_cell);

	std::vector<Tensor<1,DIM>> u_q(n_q_points);
	std::vector<double> u3_q(n_q_points);
	std::vector<Tensor<2,DIM>> grad_u_q(n_q_points);
	std::vector<Tensor<1,DIM>> grad_u3_q(n_q_points);
	std::vector<Tensor<1,DIM>> grad_u3_old_q(n_q_points);
	std::vector<Tensor<1,DIM>> w_q(n_q_points);
	std::vector<Tensor<2,DIM>> grad_w_q(n_q_points);
	std::vector<Tensor<1,DIM>> lagrange_multiplier_q(n_q_points);
	std::vector<Tensor<1,DIM>> u_old_q(n_q_points);
	std::vector<double> u3_old_q(n_q_points);


	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);


	std::vector<double> lambda_values(n_q_points);
	std::vector<double> mu_values(n_q_points);
	std::vector<double> nu_values(n_q_points);
	std::vector<double> Emod_values(n_q_points);


	Functions::ConstantFunction<DIM> lambda(1.), mu(1.), nu(nu_global);


	std::vector<Tensor<1, DIM+5>> rhs_values(n_q_points);

	std::vector<Tensor<1,DIM+1>> stab_values(n_q_points);



	StrainEnergy = 0.0;
	rp = 0.;
	rd = 0.;

	// Now we can begin with the loop over all cells:
	for (const auto &cell : dof_handler.active_cell_iterators())
	{
		unsigned int cell_index = cell->active_cell_index();


		const FEValuesExtractors::Vector u(0);
		const FEValuesExtractors::Scalar u3(DIM);
		const FEValuesExtractors::Vector w(DIM+1);
		const FEValuesExtractors::Vector lagrange_multiplier(DIM+3);

		cell_matrix = 0;
		cell_rhs    = 0;

		fe_values.reinit(cell);

		// Next we get the values of the coefficients at the quadrature
		// points. Likewise for the right hand side:
		lambda.value_list(fe_values.get_quadrature_points(), lambda_values);
		mu.value_list(fe_values.get_quadrature_points(), mu_values);
		nu.value_list(fe_values.get_quadrature_points(), nu_values);
		right_hand_side(fe_values.get_quadrature_points(), rhs_values);
		solution_stabilization(fe_values.get_quadrature_points(), stab_values);
		Youngs_Modulus(fe_values.get_quadrature_points(), Emod_values);



		fe_values[u].get_function_values(evaluation_point,u_q);
		fe_values[u3].get_function_values(evaluation_point,u3_q);
		fe_values[u].get_function_gradients(evaluation_point,grad_u_q);
		fe_values[u3].get_function_gradients(evaluation_point,grad_u3_q);
		fe_values[w].get_function_values(evaluation_point,w_q);
		fe_values[w].get_function_gradients(evaluation_point,grad_w_q);
		fe_values[lagrange_multiplier].get_function_values(evaluation_point,lagrange_multiplier_q);
		fe_values[u].get_function_values(prev_solution,u_old_q);
		fe_values[u3].get_function_values(prev_solution,u3_old_q);
		fe_values[u3].get_function_gradients(prev_solution,grad_u3_old_q);



		for (const unsigned int q_point :
				fe_values.quadrature_point_indices())
		{

			Tensor<1,DIM> du3F = grad_u3_q[q_point] - w_q[q_point];
			Tensor<1,DIM> du3old = grad_u3_q[q_point] - grad_u3_old_q[q_point];
			rp += scalar_product(du3F,du3F)*fe_values.JxW(q_point);
			rd += scalar_product(du3old,du3old)*fe_values.JxW(q_point);




		}



	}

	rp = sqrt(rp);
	rd = (muval/Eval)*sqrt(rd);
	std::cout << "Primal feasibility: " <<rp << std::endl;
	std::cout << "Dual feasibility: " <<rd << std::endl;




}

void ElasticProblem::setup_constraints(){
	constraints.clear();

	DoFTools::make_hanging_node_constraints(dof_handler, constraints);
	//VectorTools::interpolate_boundary_values(dof_handler,
	//		1,
	//		Functions::ZeroFunction<DIM>(DIM+5),
	//		constraints);

	// Constraint stuff
	std::vector<bool> u1_components = {true,false,false,false,false,false,false};
	ComponentMask u1_mask(u1_components);

	std::vector<bool> u2_components = {false,true,false,false,false,false,false};
	ComponentMask u2_mask(u2_components);

	std::vector<bool> u3_components = {false,false,true,false,false,false,false};
	ComponentMask u3_mask(u3_components);

	std::vector<bool> w1_components = {false,false,false,true,false,false,false};
	ComponentMask w1_mask(w1_components);

	std::vector<bool> w2_components = {false,false,false,false,true,false,false};
	ComponentMask w2_mask(w2_components);

	std::vector<bool> is_u1_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u1_mask, is_u1_comp);

	std::vector<bool> is_u2_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u2_mask, is_u2_comp);

	std::vector<bool> is_u3_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u3_mask, is_u3_comp);

	std::vector<bool> is_w1_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, w1_mask, is_w1_comp);

	std::vector<bool> is_w2_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, w2_mask, is_w2_comp);

	std::vector<Point<DIM>> support_points(number_dofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);

	unsigned int closestzeropointind = 0;
	double tempdist = sqrt(pow(support_points[0][0],2.0) + pow(support_points[0][1],2.0));

	unsigned int furthestrightmidpointind = 0;
	double rightdist = sqrt(pow(support_points[0][0] - L/2.0,2.0) + pow(support_points[0][1],2.0));

	for (unsigned int i = 0; i < support_points.size(); i++) {
		double idist = sqrt(pow(support_points[i][0],2.0) + pow(support_points[i][1],2.0));
		if (idist < tempdist){
			tempdist = idist;
			closestzeropointind = i;
		}

		double temprdist = sqrt(pow(support_points[i][0] - L/2.0,2.0) + pow(support_points[i][1],2.0));

		if (temprdist < rightdist) {
			rightdist = temprdist;
			furthestrightmidpointind = i;
		}
	}

	for (unsigned int i = 0; i < number_dofs; i++) {

		if (fabs(support_points[i][0] - support_points[closestzeropointind][0]) < 1.0e-6
				&& fabs(support_points[i][1]- support_points[closestzeropointind][1]) < 1.0e-6
				&& (is_u1_comp[i] || is_u2_comp[i] || is_u3_comp[i] || is_w1_comp[i] || is_w2_comp[i])) {
			constraints.add_line(i);
		}

		if (fabs(support_points[i][0] - support_points[furthestrightmidpointind][0]) < 1.0e-6
				&& fabs(support_points[i][1]- support_points[furthestrightmidpointind][1]) < 1.0e-6
				&& is_u2_comp[i]) {
			constraints.add_line(i);
		}

		/*
		if (fabs(support_points[i][0] + L/2) < 1.0e-6
				&& fabs(support_points[i][1] + L/2) < 1.0e-6
				&& (is_u1_comp[i] || is_u2_comp[i] || is_u3_comp[i])) {
			constraints.add_line(i);
		}

		if (fabs(support_points[i][0] + L/2) < 1.0e-6
				&& fabs(support_points[i][1] - L/2) < 1.0e-6
				&& (is_u1_comp[i] || is_u3_comp[i])) {
			constraints.add_line(i);

		}
		if (fabs(support_points[i][0] - L/2) < 1.0e-6
				&& fabs(support_points[i][1] - L/2) < 1.0e-6
				&& (is_u3_comp[i])) {
			constraints.add_line(i);

		}

		if (fabs(support_points[i][0] - L/2) < 1.0e-6
				&& fabs(support_points[i][1] + L/2) < 1.0e-6
				&& (is_u3_comp[i])) {
			constraints.add_line(i);

		}
		 */

	}


	constraints.close();
}



void ElasticProblem::output_results(const unsigned int cycle) const
{
	DataOut<DIM> data_out;
	data_out.attach_dof_handler(dof_handler);

	std::vector<std::string> solution_names;
	switch (DIM)
	{
	case 1:
		solution_names.emplace_back("displacement");
		break;
	case 2:
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		solution_names.emplace_back("z_displacement");
		solution_names.emplace_back("w_1_gradient");
		solution_names.emplace_back("w_2_gradient");
		solution_names.emplace_back("l_1");
		solution_names.emplace_back("l_2");
		break;
	case 3:
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		solution_names.emplace_back("z_displacement");
		break;
	default:
		Assert(false, ExcNotImplemented());
	}

	// After setting up the names for the different components of the
	// solution vector, we can add the solution vector to the list of
	// data vectors scheduled for output. Note that the following
	// function takes a vector of strings as second argument, whereas
	// the one which we have used in all previous examples accepted a
	// string there. (In fact, the function we had used before would
	// convert the single string into a vector with only one element
	// and forwards that to the other function.)
	data_out.add_data_vector(solution, solution_names);
	data_out.build_patches();

	std::ofstream output("solutions/solution-" + std::to_string(cycle) + ".vtk");
	data_out.write_vtk(output);
}

void ElasticProblem::output_zero_results(const unsigned int cycle) const
{
	DataOut<DIM> data_out;
	data_out.attach_dof_handler(dof_handler);

	std::vector<std::string> solution_names;
	switch (DIM)
	{
	case 1:
		solution_names.emplace_back("displacement");
		break;
	case 2:
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		solution_names.emplace_back("z_displacement");
		solution_names.emplace_back("w_1_gradient");
		solution_names.emplace_back("w_2_gradient");
		solution_names.emplace_back("l_1");
		solution_names.emplace_back("l_2");
		break;
	case 3:
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		solution_names.emplace_back("z_displacement");
		break;
	default:
		Assert(false, ExcNotImplemented());
	}

	// After setting up the names for the different components of the
	// solution vector, we can add the solution vector to the list of
	// data vectors scheduled for output. Note that the following
	// function takes a vector of strings as second argument, whereas
	// the one which we have used in all previous examples accepted a
	// string there. (In fact, the function we had used before would
	// convert the single string into a vector with only one element
	// and forwards that to the other function.)
	Vector<double> zero_solution(solution.size());
	for (unsigned int i = 0; i < zero_solution.size(); i++) {
		zero_solution(i) = 0.0;
	}
	data_out.add_data_vector(zero_solution, solution_names);
	data_out.build_patches();

	std::ofstream output("solutions/solution-" + std::to_string(cycle) + ".vtk");
	data_out.write_vtk(output);
}


void ElasticProblem::ReadCurvatureFiles(std::vector<double> & tvec,std::vector<double> &sigvec,std::vector<double> & kappavec){
	std::ifstream file1;
	std::ifstream file2;
	std::ifstream file3;
	/*
	 * This example assumes you are reading
	 * strings from the file, but you can
	 * change it to whatever data type that
	 * is stored in the file.
	 */
	tvec.clear();
	sigvec.clear();
	kappavec.clear();

	std::string inputString;

	file1.open("time.csv");
	while(file1>>inputString) //reads one string at a time
		tvec.push_back(std::stod(inputString)); //add it to data vector
	file1.close();


	file2.open("inplane.csv");
	while(file2>>inputString)
		sigvec.push_back(std::stod(inputString));
	file2.close();

	file3.open("bend.csv");
	while(file3>>inputString)
		kappavec.push_back(std::stod(inputString));
	file3.close();



}

void ElasticProblem::output_perturbation(const unsigned int cycle) const
{
	DataOut<DIM> data_out;
	data_out.attach_dof_handler(dof_handler);

	std::vector<std::string> solution_names;
	switch (DIM)
	{
	case 1:
		solution_names.emplace_back("displacement");
		break;
	case 2:
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		solution_names.emplace_back("z_displacement");
		solution_names.emplace_back("w_1_gradient");
		solution_names.emplace_back("w_2_gradient");
		solution_names.emplace_back("l_1");
		solution_names.emplace_back("l_2");
		break;
	case 3:
		solution_names.emplace_back("x_displacement");
		solution_names.emplace_back("y_displacement");
		solution_names.emplace_back("z_displacement");
		break;
	default:
		Assert(false, ExcNotImplemented());
	}

	// After setting up the names for the different components of the
	// solution vector, we can add the solution vector to the list of
	// data vectors scheduled for output. Note that the following
	// function takes a vector of strings as second argument, whereas
	// the one which we have used in all previous examples accepted a
	// string there. (In fact, the function we had used before would
	// convert the single string into a vector with only one element
	// and forwards that to the other function.)

	data_out.add_data_vector(perturbation_direction, solution_names);
	data_out.build_patches();

	std::ofstream output("solutions/perturbation.vtk");
	data_out.write_vtk(output);
}

void ElasticProblem::FlipConfigUD(Vector<double> & flippedvector){
	std::vector<bool> u1_components = {true,false,false,false,false,false,false};
	ComponentMask u1_mask(u1_components);

	std::vector<bool> u2_components = {false,true,false,false,false,false,false};
	ComponentMask u2_mask(u2_components);

	std::vector<bool> u3_components = {false,false,true,false,false,false,false};
	ComponentMask u3_mask(u3_components);

	std::vector<bool> w1_components = {false,false,false,true,false,false,false};
	ComponentMask w1_mask(w1_components);

	std::vector<bool> w2_components = {false,false,false,false,true,false,false};
	ComponentMask w2_mask(w2_components);

	std::vector<bool> is_u1_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u1_mask, is_u1_comp);

	std::vector<bool> is_u2_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u2_mask, is_u2_comp);

	std::vector<bool> is_u3_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, u3_mask, is_u3_comp);

	std::vector<bool> is_w1_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, w1_mask, is_w1_comp);

	std::vector<bool> is_w2_comp(number_dofs, false);
	DoFTools::extract_dofs(dof_handler, w2_mask, is_w2_comp);

	std::vector<Point<DIM>> support_points(number_dofs);
	MappingQ1<DIM> mapping;
	DoFTools::map_dofs_to_support_points(mapping, dof_handler, support_points);



	flippedvector.reinit(solution.size());
	for (unsigned int i = 0; i < number_dofs; i++) {


		if (!is_u1_comp[i] && !is_u2_comp[i]) {
			flippedvector[i] = -1.0*solution[i];
		} else {
			flippedvector[i] = solution[i];
		}


	}



}







}


#endif // ELASTIC_PROBLEM_CC_
