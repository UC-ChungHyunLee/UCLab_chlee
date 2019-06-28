#pragma once
#include "config.h"

// Declare

void checkRange(ceq ceq_t);

void findCross(leq leq_1, leq leq_2, db save);
void findCross(leq leq_1, ceq ceq_1, db save);
void findCross(ceq ceq_1, ceq ceq_2, db save);

void integral(leq leq_1, leq leq_2, db save);
void integral(leq leq_1, ceq ceq_1, db save);
void integral(ceq ceq_1, ceq ceq_2, db save);

// Define

void checkRange(ceq ceq_t) {
	float x_min = cos(getRadian(ceq_t.angle_min));
	float x_max = cos(getRadian(ceq_t.angle_max));
	float y_min = sin(getRadian(ceq_t.angle_min));
	float y_max = sin(getRadian(ceq_t.angle_max));

	for (int i = ceq_t.angle_min; i <= ceq_t.angle_max; i++) {
		// find corX_max
		if (cos(getRadian(i)) >= x_max) {
			x_max = cos(getRadian(i));
		}
		// find corX_min
		if (cos(getRadian(i) <= x_min)) {
			x_min = cos(getRadian(i));
		}
		// find corY_max
		if (sin(getRadian(i)) >= y_max) {
			y_max = sin(getRadian(i));
		}
		// find corY_min
		if (sin(getRadian(i)) <= y_min) {
			y_min = sin(getRadian(i));
		}	
	}

	ceq_t.corX_max = x_max;
	ceq_t.corX_min = x_min;
	ceq_t.corY_max = y_max;
	ceq_t.corY_min = y_min;
}


void findCross(leq leq_1, leq leq_2, db save) {
	float diff_gradient = leq_1.gradient - leq_2.gradient;	// leq_1 - leq_2 (left term)
	float diff_basis = leq_2.basis - leq_1.basis;			// leq_2 - leq_1 (right term)
	float temp_corX = 0.0;

	// If those equations are parallel, it can be impossible or indeterminate
	if (diff_gradient == 0) {
		// impossible
		if (diff_basis == 0) {
			save.cross_type[save.ptr] = IMPOSSIBLE;
			save.cross_corX[save.ptr] = 0;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = LINEAR;
			save.ptr++;
		}
		// indeterminate
		else {
			save.cross_type[save.ptr] = INDETERMINATE;
			save.cross_corX[save.ptr] = 0;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = LINEAR;
			save.ptr++;
		}
	}
	// Having a root
	else {
		temp_corX = diff_basis / diff_gradient;
		// Check that the cross coordinate is in the line segment
		if ((leq_1.corX_min <= temp_corX) && (leq_1.corX_max >= temp_corX) && (leq_2.corX_min <= temp_corX) && (leq_2.corX_max >= temp_corX)) {
			save.cross_type[save.ptr] = POSSIBLE;
			save.cross_corX[save.ptr] = temp_corX;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = LINEAR;
			save.ptr++;
		}
		// Else, it means that the cross coordinate is not in the line segment
		else {
			save.cross_type[save.ptr] = IMPOSSIBLE;
			save.cross_corX[save.ptr] = 0;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = LINEAR;
			save.ptr++;
		}
	}
}

void findCross(leq leq_1, ceq ceq_1, db save) {
	float temp_a = leq_1.gradient * leq_1.gradient + 1;
	float temp_b = ((-2.0) * ceq_1.corX) + (2 * leq_1.gradient * (leq_1.basis - ceq_1.corY));
	float temp_c = (ceq_1.corX * ceq_1.corX) + (leq_1.basis - ceq_1.corY) * (leq_1.basis - ceq_1.corY) - (ceq_1.radius * ceq_1.radius);
	float determinate = (temp_b * temp_b) - (4 * temp_a * temp_c);
	

	// impossible
	if (determinate < 0) {
		save.cross_type[save.ptr] = IMPOSSIBLE;
		save.cross_corX[save.ptr] = 0;
		save.cross_eq1[save.ptr] = LINEAR;
		save.cross_eq2[save.ptr] = CIRCULAR;
		save.ptr++;
	}
	// dual root
	else if (determinate = 0) {
		float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / 2 * temp_a * temp_c;
		float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / 2 * temp_a * temp_c;
		float tempA_corY = root_a * leq_1.gradient + leq_1.basis;
		float tempB_corY = root_b * leq_1.gradient + leq_1.basis;
				
		if (root_a == root_b) {
			// Check that the cross coordinate is in the condition range.
			if ((leq_1.corX_min <= root_a) && (leq_1.corX_max >= root_a) &&																		// Check that the cross coordinate is in the line segment.
				(ceq_1.corX_min <= root_a) && (ceq_1.corX_max >= root_a) && (ceq_1.corY_min <= tempA_corY) && (ceq_1.corY_max >= tempA_corY)) {	// Check that the cross coordinate is in the arc of the sector.
			
				save.cross_type[save.ptr] = DUALROOT;
				save.cross_corX[save.ptr] = root_a;
				save.cross_eq1[save.ptr] = LINEAR;
				save.cross_eq2[save.ptr] = CIRCULAR;
				save.ptr++;
			}
			// Else, it means that the cross coordinate is not in the condition range.
			else {
				save.cross_type[save.ptr] = IMPOSSIBLE;
				save.cross_corX[save.ptr] = root_a;
				save.cross_eq1[save.ptr] = LINEAR;
				save.cross_eq2[save.ptr] = CIRCULAR;
				save.ptr++;
			}
		}
		else {
			printf("Dual Root ERR (LEQ & CEQ)....\n");
			system("PUASE");
		}
	}
	
	// two different roots
	else {
		float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / 2 * temp_a * temp_c;
		float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / 2 * temp_a * temp_c;
		float tempA_corY = root_a * leq_1.gradient + leq_1.basis;
		float tempB_corY = root_b * leq_1.gradient + leq_1.basis;

		// Check that the cross coordinate is in the condition range.
		if ((leq_1.corX_min <= root_a) && (leq_1.corX_max >= root_a) &&																		// Check that the cross coordinate A is in the line segment.
			(ceq_1.corX_min <= root_a) && (ceq_1.corX_max >= root_a) && (ceq_1.corY_min <= tempA_corY) && (ceq_1.corY_max >= tempA_corY)) {	// Check that the cross coordinate A is in the arc of the sector.
			save.cross_type[save.ptr] = POSSIBLE;
			save.cross_corX[save.ptr] = root_a;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = CIRCULAR;
			save.ptr++;
		}
		else {
			save.cross_type[save.ptr] = IMPOSSIBLE;
			save.cross_corX[save.ptr] = 0;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = CIRCULAR;
			save.ptr++;
		}
		if ((leq_1.corX_min <= root_b) && (leq_1.corX_max >= root_b) && 																	// Check that the cross coordinate B is in the line segment.
			(ceq_1.corX_min <= root_b) && (ceq_1.corX_max >= root_b) && (ceq_1.corY_min <= tempB_corY) && (ceq_1.corY_max >= tempB_corY)) {	// Check that the cross coordinate B is in the arc of the sector.
			save.cross_type[save.ptr] = POSSIBLE;
			save.cross_corX[save.ptr] = root_b;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = CIRCULAR;
			save.ptr++;
		}
		else {
			save.cross_type[save.ptr] = IMPOSSIBLE;
			save.cross_corX[save.ptr] = 0;
			save.cross_eq1[save.ptr] = LINEAR;
			save.cross_eq2[save.ptr] = CIRCULAR;
			save.ptr++;
		}
		
	}
}

void findCross(ceq ceq_1, ceq ceq_2, db save) {
	// Since we set the 4 sectors already, we can define easily.
	// First, we find the parallel line between two circle.
	// Second, we determine whether the line has root with those circles.

	// Pre-condition : congruence
	if ((ceq_1.corX == ceq_2.corX) && (ceq_1.corY == ceq_2.corY) && (ceq_1.radius == ceq_2.radius)) {
		save.cross_type[save.ptr] = CONGRUENCE;
		save.cross_corX[save.ptr] = 0;
		save.cross_eq1[save.ptr] = CIRCULAR;
		save.cross_eq2[save.ptr] = CIRCULAR;
		save.ptr++;
	}
	// Pre-condition : Contain (is equal to impossible)
	else if ((ceq_1.corX == ceq_2.corX) && (ceq_1.corY == ceq_2.corY)) {
		save.cross_type[save.ptr] = IMPOSSIBLE;
		save.cross_corX[save.ptr] = 0;
		save.cross_eq1[save.ptr] = CIRCULAR;
		save.cross_eq2[save.ptr] = CIRCULAR;
		save.ptr++;
	}
	// Pre-condition : The locational relationship between the two circle is horizontal (it means that the coordinate X of the central of the two circle are same : y = k)
	else if (ceq_1.corX == ceq_1.corX) {
		float temp_gradient = (ceq_2.corX - ceq_1.corX) / (ceq_1.corY - ceq_2.corY); // The gradient of the line is 0;
		float temp_basis = (ceq_2.radius * ceq_2.radius - ceq_1.radius * ceq_1.radius) / 2 * (ceq_1.corY - ceq_2.corY);
		float temp_a = temp_gradient * temp_gradient + 1;
		float temp_b = ((-2.0) * ceq_1.corX) + (2 * temp_gradient * temp_basis) + ((-2.0) * temp_gradient * ceq_1.corY);
		float temp_c = (ceq_1.corX * ceq_1.corX) + ((-2.0) * ceq_1.corY * temp_basis) + (ceq_1.corY * ceq_1.corY) - (ceq_1.radius * ceq_1.radius);
		float determinate = (temp_b * temp_b) - (4 * temp_a * temp_c);

		// Impossible
		if (determinate < 0) {
			save.cross_type[save.ptr] = IMPOSSIBLE;
			save.cross_corX[save.ptr] = 0;
			save.cross_eq1[save.ptr] = CIRCULAR;
			save.cross_eq2[save.ptr] = CIRCULAR;
			save.ptr++;
		}
		// Dual root (In between the two circle, Dual root has no meaning for overlapping)
		else if (determinate = 0) {
			float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / 2 * temp_a * temp_c;
			float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / 2 * temp_a * temp_c;
			
			if (root_a == root_b) {
				save.cross_type[save.ptr] = DUALROOT;
				save.cross_corX[save.ptr] = root_a;
				save.cross_eq1[save.ptr] = CIRCULAR;
				save.cross_eq2[save.ptr] = CIRCULAR;
				save.ptr++;
			}
			else {
				printf("DUAL ROOT ERR... (CEQ and CEQ) \n");
				system("PAUSE");
			}
		}
		else if (determinate > 0) {
			
		}
	}
}