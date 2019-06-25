#pragma once
#include "config.h"

// Declare

void findCross(leq leq_1, leq leq_2, db save);
void findCross(leq leq_1, ceq ceq_1, db save);
void findCross(ceq ceq_1, ceq ceq_2, db save);

void integral(leq leq_1, leq leq_2, db save);
void integral(leq leq_1, ceq ceq_1, db save);
void integral(ceq ceq_1, ceq ceq_2, db save);

// Define

void findCross(leq leq_1, leq leq_2, db save) {
	float diff_gradient = leq_1.gradient - leq_2.gradient;	// leq_1 - leq_2 (left term)
	float diff_basis = leq_2.basis - leq_1.basis;			// leq_2 - leq_1 (right term)

	// If those equations are parallel, it can be impossible or indeterminate
	if (diff_gradient == 0) {
		// impossible
		if (diff_basis == 0) {
			save.cross_type[save.ptr_bool] = IMPOSSIBLE;
			save.cross_corX[save.ptr_corX] = 0;
			save.ptr_bool++;
			save.ptr_corX++;
		}
		// indeterminate
		else {
			save.cross_type[save.ptr_bool] = INDETERMINATE;
			save.cross_corX[save.ptr_corX] = 0;
			save.ptr_bool++;
			save.ptr_corX++;
		}
	}
	// Having a root
	else {
		save.cross_type[save.ptr_bool] = POSSIBLE;
		save.cross_corX[save.ptr_corX] = diff_basis / diff_gradient;
		save.ptr_bool++;
		save.ptr_corX++;
	}
}

void findCross(leq leq_1, ceq ceq_1, db save) {
	float determinate = 
}