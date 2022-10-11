#ifndef PERTURBATION_CLASS_CC_
#define PERTURBATION_CLASS_CC_
#include "perturbation_class.h"

namespace Step8
{
using namespace dealii;




double PerturbationValues::value(const Point<DIM> &p,
		const unsigned int i) const
{

	Tensor<2,DIM> H;
	// /*
	// Diagonal Deformation
	H[0][0] = 0.0;
	H[0][1] = 0.0;
	H[1][0] = 0.0;
	H[1][1] = -1.0;
	// */

	/*
			// Aligned Deformation
			H[0][0] = 1.0;

	 */


	Tensor<1,DIM> x;
	x[0] = p(0);
	x[1] = p(1);
	//std::cout << i << std::endl;

	double u3out = 0.5*(x*(H*x));
	Tensor<1,DIM> wvecout = H*x;

	if (i==0) {
		// u1


	} else if (i == 1) {
		// u2


	} else if (i == 2) {
		// u3
		return u3out;

	} else if (i == 3) {
		// w1

		return wvecout[0];


	} else if (i == 4) {
		// w2
		return wvecout[1];



	} else if (i == 5) {
		// l1


	} else if (i == 6) {
		// l2


	}
	return 0.0;

}

void PerturbationValues::vector_value ( const Point<DIM> &p,
		Vector<double>   &values ) const
{
	for (unsigned int c = 0; c < this->n_components; ++c)
		values(c) = value (p, c);
}

}



#endif // PERTURBATION_CLASS_CC_
