#pragma once
#include "config.h"
#define M_PI 3.14159265358979323846

// Declare

double getRadian(float angle);
float min(float fleq, float sleq, float tceq);
float min(float a, float b);
float max(float fleq, float sleq, float tceq);
float max(float a, float b);

void init_SU(SU* su);
void init_DB(db* save);

void checkRange(ceq* ceq_t);
void sortCoordinate(db* save);

void findCross(leq leq_1, leq leq_2, db* save);
void findCross(leq leq_1, ceq ceq_1, db* save);
void findCross(ceq ceq_1, ceq ceq_2, db* save);

void cal(SU su1, SU su2, db save[]);

void integral(leq leq_1, leq leq_2, db save);
void integral(leq leq_1, ceq ceq_1, db save);
void integral(ceq ceq_1, ceq ceq_2, db save);



// Define

double getRadian(float angle) {
	return angle * (M_PI / 180);
}

float min(float fleq, float sleq, float tceq) {
	if (fleq <= sleq) {
		if (fleq <= tceq) {
			return fleq;
		}
		else {
			return tceq;
		}
	}
	else {
		if (sleq <= tceq) {
			return sleq;
		}
		else {
			return tceq;
		}
	}
}

float min(float a, float b) {
	if (a <= b) {
		return a;
	}
	else {
		return b;
	}
}

float max(float fleq, float sleq, float tceq) {
	if (fleq >= sleq) {
		if (fleq >= tceq) {
			return fleq;
		}
		else {
			return tceq;
		}
	}
	else {
		if (sleq >= tceq) {
			return sleq;
		}
		else {
			return tceq;
		}
	}
}

float max(float a, float b) {
	if (a >= b) {
		return a;
	}
	else {
		return b;
	}
}

void init_SU(SU* su) {
	
	for (int i = 0; i < 4; i++) {
		su->beam[i].ch = 0;
		su->beam[i].root = su->index;
		su->beam[i].index = i;
		su->beam[i].fleq.root = su->index;
		su->beam[i].fleq.num = 1;
		su->beam[i].sleq.root = su->index;
		su->beam[i].sleq.num = 2;
		su->beam[i].tceq.root = su->index;
		su->beam[i].fleq.index = su->beam[i].index;
		su->beam[i].sleq.index = su->beam[i].index;
		su->beam[i].tceq.index = su->beam[i].index;

		su->beam[i].fleq.gradient = tan(getRadian(45 + 90 * i));
		su->beam[i].fleq.basis = su->corY - su->corX * su->beam[i].fleq.gradient;
		if (cos(getRadian(45 + 90 * i)) < 0) {
			su->beam[i].fleq.corX_min = su->corX + su->radius * cos(getRadian(45 + 90 * i));
			su->beam[i].fleq.corX_max = su->corX;
			su->beam[i].fleq.corY_min = su->corY;
			su->beam[i].fleq.corY_max = su->corY + su->radius * sin(getRadian(45 + 90 * i));
		}
		else {
			su->beam[i].fleq.corX_max = su->corX + su->radius * cos(getRadian(45 + 90 * i));
			su->beam[i].fleq.corX_min = su->corX;
			su->beam[i].fleq.corY_max = su->corY;
			su->beam[i].fleq.corY_min = su->corY + su->radius * sin(getRadian(45 + 90 * i));
		}

		su->beam[i].sleq.gradient = tan(getRadian(45 + 90 * (i + 1)));
		su->beam[i].sleq.basis = su->corY - su->corX * su->beam[i].sleq.gradient;
		if (cos(getRadian(45 + 90 * (i + 1))) < 0) {
			su->beam[i].sleq.corX_min = su->corX + su->radius * cos(getRadian(45 + 90 * (i + 1)));
			su->beam[i].sleq.corX_max = su->corX;
			su->beam[i].fleq.corY_min = su->corY;
			su->beam[i].fleq.corY_max = su->corY + su->radius * sin(getRadian(45 + 90 * (i + 1)));
		}
		else {
			su->beam[i].sleq.corX_max = su->corX + su->radius * cos(getRadian(45 + 90 * (i + 1)));
			su->beam[i].sleq.corX_min = su->corX;
			su->beam[i].fleq.corY_max = su->corY;
			su->beam[i].fleq.corY_min = su->corY + su->radius * sin(getRadian(45 + 90 * (i + 1)));
		}

		su->beam[i].tceq.angle_min = 45 + 90 * i;
		su->beam[i].tceq.angle_max = 45 + 90 * (i + 1);
		su->beam[i].tceq.corX = su->corX;
		su->beam[i].tceq.corY = su->corY;
		su->beam[i].tceq.radius = su->radius;
		checkRange(&su->beam[i].tceq);

		su->beam[i].x_min = min(su->beam[i].sleq.corX_min, su->beam[i].fleq.corX_min, su->beam[i].tceq.corX_min);
		su->beam[i].x_max = max(su->beam[i].sleq.corX_max, su->beam[i].fleq.corX_max, su->beam[i].tceq.corX_max);
		su->beam[i].y_min = min(su->beam[i].sleq.corY_min, su->beam[i].fleq.corY_min, su->beam[i].tceq.corY_min);
		su->beam[i].y_max = max(su->beam[i].sleq.corY_max, su->beam[i].fleq.corY_max, su->beam[i].tceq.corY_max);
	}
}

void init_DB(db* save) {
	save->ptr = 0;

	for (int i = 0; i < BUFFERSIZE; i++) {
		save->cross_corX[i] = 0.0;
		save->cross_corY[i] = 0.0;
		save->cross_type[i] = -1;
		save->cross_eq1[i] = -1;
		save->cross_eq2[i] = -1;
		save->cross_ceq1[i] = ceq_temp;
		save->cross_ceq2[i] = ceq_temp;
		save->cross_leq1[i] = leq_temp;
		save->cross_leq2[i] = leq_temp;		
	}
}

void checkRange(ceq* ceq_t) {
	float x_min = ceq_t->corX + ceq_t->radius * cos(getRadian(ceq_t->angle_min));
	float x_max = ceq_t->corX + ceq_t->radius * cos(getRadian(ceq_t->angle_max));
	float y_min = ceq_t->corY + ceq_t->radius * sin(getRadian(ceq_t->angle_min));
	float y_max = ceq_t->corY + ceq_t->radius * sin(getRadian(ceq_t->angle_max));

	for (int i = ceq_t->angle_min; i <= ceq_t->angle_max; i++) {
		// find corX_max
		if ((ceq_t->corX + ceq_t->radius * cos(getRadian(i))) >= x_max) {
			x_max = ceq_t->corX + ceq_t->radius * cos(getRadian(i));
		}
		// find corX_min
		if ((ceq_t->corX + ceq_t->radius * cos(getRadian(i))) <= x_min) {
			x_min = ceq_t->corX + ceq_t->radius * cos(getRadian(i));
		}
		// find corY_max
		if ((ceq_t->corY + ceq_t->radius * sin(getRadian(i))) >= y_max) {
			y_max = ceq_t->corX + ceq_t->radius * sin(getRadian(i));
		}
		// find corY_min
		if ((ceq_t->corY + ceq_t->radius * sin(getRadian(i))) <= y_min) {
			y_min = ceq_t->corX + ceq_t->radius * sin(getRadian(i));
		}	
	}

	ceq_t->corX_max = x_max;
	ceq_t->corX_min = x_min;
	ceq_t->corY_max = y_max;
	ceq_t->corY_min = y_min;
}

// Bubble Sorting
void sortCoordinate(db* save) {
	float temp_corX;
	float temp_corY;
	float temp_eq_1;
	float temp_eq_2;
	leq temp_leq_1;
	leq temp_leq_2;
	ceq temp_ceq_1;
	ceq temp_ceq_2;
	for (int i = 0; i < save->ptr; i++) {
		for (int k = i + 1; k < save->ptr; k++) {
			if (save->cross_corX[k] < save->cross_corX[i]) {
				temp_corX = save->cross_corX[i];
				temp_corY = save->cross_corY[i];
				temp_eq_1 = save->cross_eq1[i];
				temp_eq_2 = save->cross_eq2[i];
				temp_leq_1 = save->cross_leq1[i];
				temp_leq_2 = save->cross_leq2[i];
				temp_ceq_1 = save->cross_ceq1[i];
				temp_ceq_2 = save->cross_ceq2[i];
				save->cross_corX[i] = save->cross_corX[k];
				save->cross_corY[i] = save->cross_corY[k];
				save->cross_eq1[i] = save->cross_eq1[k];
				save->cross_eq2[i] = save->cross_eq2[k];
				save->cross_leq1[i] = save->cross_leq1[k];
				save->cross_leq2[i] = save->cross_leq2[k];
				save->cross_ceq1[i] = save->cross_ceq1[k];
				save->cross_ceq2[i] = save->cross_ceq2[k];
				save->cross_corX[k] = temp_corX;
				save->cross_corY[k] = temp_corY;
				save->cross_eq1[k] = temp_eq_1;
				save->cross_eq2[k] = temp_eq_2;
				save->cross_leq1[k] = temp_leq_1;
				save->cross_leq2[k] = temp_leq_2;
				save->cross_ceq1[k] = temp_ceq_1;				
				save->cross_ceq2[k] = temp_ceq_2;
			}
		}
	}
}

void findCross(leq leq_1, leq leq_2, db* save) {
	float diff_gradient = leq_1.gradient - leq_2.gradient;	// leq_1 - leq_2 (left term)
	float diff_basis = leq_2.basis - leq_1.basis;			// leq_2 - leq_1 (right term)
	float temp_corX = 0.0;

	// If those equations are parallel, it can be impossible or indeterminate
	if (diff_gradient == 0) {
		// impossible
		if (diff_basis == 0) {
			/*
			save->cross_type[save->ptr] = IMPOSSIBLE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = LINEAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_2;
			save->cross_ceq1[save->ptr] = ceq_temp;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
			*/
		}
		// indeterminate
		else {
			/*
			save->cross_type[save->ptr] = INDETERMINATE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = LINEAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_2;
			save->cross_ceq1[save->ptr] = ceq_temp;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
			*/
		}
	}
	// Having a root
	else {
		temp_corX = diff_basis / diff_gradient;
		// Check that the cross coordinate is in the line segment
		if ((leq_1.corX_min <= temp_corX) && (leq_1.corX_max >= temp_corX) && (leq_2.corX_min <= temp_corX) && (leq_2.corX_max >= temp_corX)) {
			save->cross_type[save->ptr] = POSSIBLE;
			save->cross_corX[save->ptr] = temp_corX;
			save->cross_corY[save->ptr] = leq_1.gradient * temp_corX + leq_1.basis;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = LINEAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_2;
			save->cross_ceq1[save->ptr] = ceq_temp;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
		}
		// Else, it means that the cross coordinate is not in the line segment
		else {
			/*
			save->cross_type[save->ptr] = IMPOSSIBLE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = LINEAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_2;
			save->cross_ceq1[save->ptr] = ceq_temp;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
			*/
		}
	}
}

void findCross(leq leq_1, ceq ceq_1, db* save) {
	float temp_a = leq_1.gradient * leq_1.gradient + 1;
	float temp_b = ((-2.0) * ceq_1.corX) + (2 * leq_1.gradient * (leq_1.basis - ceq_1.corY));
	float temp_c = (ceq_1.corX * ceq_1.corX) + (leq_1.basis - ceq_1.corY) * (leq_1.basis - ceq_1.corY) - (ceq_1.radius * ceq_1.radius);
	float determinate = (temp_b * temp_b) - (4 * temp_a * temp_c);
	

	// impossible
	if (determinate < 0) {
		/*
		save->cross_type[save->ptr] = IMPOSSIBLE;
		save->cross_corX[save->ptr] = 0;
		save->cross_eq1[save->ptr] = LINEAR;
		save->cross_eq2[save->ptr] = CIRCULAR;
		save->cross_leq1[save->ptr] = leq_1;
		save->cross_leq2[save->ptr] = leq_temp;
		save->cross_ceq1[save->ptr] = ceq_1;
		save->cross_ceq2[save->ptr] = ceq_temp;
		save->ptr++;
		*/
	}
	// dual root
	else if (determinate == 0) {
		float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / (2 * temp_a * temp_c);
		float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / (2 * temp_a * temp_c);
		float tempA_corY = root_a * leq_1.gradient + leq_1.basis;
		float tempB_corY = root_b * leq_1.gradient + leq_1.basis;
				
		if (root_a == root_b) {
			// Check that the cross coordinate is in the condition range.
			if ((leq_1.corX_min <= root_a) && (leq_1.corX_max >= root_a) &&																		// Check that the cross coordinate is in the line segment.
				(ceq_1.corX_min <= root_a) && (ceq_1.corX_max >= root_a) && (ceq_1.corY_min <= tempA_corY) && (ceq_1.corY_max >= tempA_corY)) {	// Check that the cross coordinate is in the arc of the sector.
			
				save->cross_type[save->ptr] = DUALROOT;
				save->cross_corX[save->ptr] = root_a;
				save->cross_corY[save->ptr] = root_a * leq_1.gradient + leq_1.basis;
				save->cross_eq1[save->ptr] = LINEAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_1;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_temp;
				save->ptr++;
			}
			// Else, it means that the cross coordinate is not in the condition range.
			else {
				/*
				save->cross_type[save->ptr] = IMPOSSIBLE;
				save->cross_corX[save->ptr] = root_a;
				save->cross_eq1[save->ptr] = LINEAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_1;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_temp;
				save->ptr++;
				*/
			}
		}
		else {
			printf("Dual Root ERR (LEQ & CEQ)....\n");
			system("PUASE");
		}
	}
	
	// two different roots
	else {
		float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / (2 * temp_a);
		float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / (2 * temp_a);
		float tempA_corY = root_a * leq_1.gradient + leq_1.basis;
		float tempB_corY = root_b * leq_1.gradient + leq_1.basis;

		// Check that the cross coordinate is in the condition range.
		if ((leq_1.corX_min <= root_a) && (leq_1.corX_max >= root_a) &&																		// Check that the cross coordinate A is in the line segment.
			(ceq_1.corX_min <= root_a) && (ceq_1.corX_max >= root_a) && (ceq_1.corY_min <= tempA_corY) && (ceq_1.corY_max >= tempA_corY)) {	// Check that the cross coordinate A is in the arc of the sector.
			save->cross_type[save->ptr] = POSSIBLE;
			save->cross_corX[save->ptr] = root_a;
			save->cross_corY[save->ptr] = root_a * leq_1.gradient + leq_1.basis;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
		}
		else {
			/*
			save->cross_type[save->ptr] = IMPOSSIBLE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
			*/
		}
		if ((leq_1.corX_min <= root_b) && (leq_1.corX_max >= root_b) && 																	// Check that the cross coordinate B is in the line segment.
			(ceq_1.corX_min <= root_b) && (ceq_1.corX_max >= root_b) && (ceq_1.corY_min <= tempB_corY) && (ceq_1.corY_max >= tempB_corY)) {	// Check that the cross coordinate B is in the arc of the sector.
			save->cross_type[save->ptr] = POSSIBLE;
			save->cross_corX[save->ptr] = root_b;
			save->cross_corY[save->ptr] = root_b * leq_1.gradient + leq_1.basis;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
		}
		else {
			/*
			save->cross_type[save->ptr] = IMPOSSIBLE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = LINEAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_1;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_temp;
			save->ptr++;
			*/
		}
		
	}
}

void findCross(ceq ceq_1, ceq ceq_2, db* save) {
	// Since we set the 4 sectors already, we can define easily.
	// First, we find the parallel line between two circle.
	// Second, we determine whether the line has root with those circles.

	// Pre-condition : 
	float det_ceq = sqrt(((ceq_1.corX - ceq_2.corX) * (ceq_1.corX - ceq_2.corX)) + ((ceq_1.corY - ceq_2.corY) * (ceq_1.corY - ceq_2.corY)));
	if ((det_ceq < (ceq_1.radius + ceq_2.radius)) && (det_ceq > 0)) {
		// Pre-condition : The locational relationship between the two circle is horizontal(it means that the coordinate Y of the central of the two circle are same : x = k)
		if ((ceq_1.corX != ceq_2.corX) && (ceq_1.corY == ceq_2.corY)) {
			float temp_corX = ((ceq_1.corX * ceq_1.corX - ceq_2.corX * ceq_2.corX) - (ceq_1.radius * ceq_1.radius - ceq_2.radius * ceq_2.radius)) / (2 * (ceq_1.corX - ceq_2.corX));
			if ((temp_corX >= ceq_1.corX_min) && (temp_corX <= ceq_1.corX_max) && (temp_corX >= ceq_1.corX_min) && (temp_corX <= ceq_1.corX_max)) {
				save->cross_type[save->ptr] = POSSIBLE;
				save->cross_corX[save->ptr] = temp_corX;
				save->cross_corY[save->ptr] = ceq_1.corY - sqrt(ceq_1.radius * ceq_1.radius - (temp_corX - ceq_1.corX) * (temp_corX - ceq_1.corX));
				save->cross_eq1[save->ptr] = CIRCULAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_temp;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_2;
				save->ptr++;
				save->cross_type[save->ptr] = POSSIBLE;
				save->cross_corX[save->ptr] = temp_corX;
				save->cross_corY[save->ptr] = ceq_1.corY + sqrt(ceq_1.radius * ceq_1.radius - (temp_corX - ceq_1.corX) * (temp_corX - ceq_1.corX));
				save->cross_eq1[save->ptr] = CIRCULAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_temp;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_2;
				save->ptr++;
			}
			else {
				/* void */
			}
		}
		// Pre-condition : The locational relationship between the two circle is vertical (it means that the coordinate X of the central of the two circle are same : y = k)
		else if ((ceq_1.corX == ceq_2.corX) && (ceq_1.corY != ceq_2.corY)) {
			float temp_corY = ((ceq_1.corY * ceq_1.corY - ceq_2.corY * ceq_2.corY) - (ceq_1.radius * ceq_1.radius - ceq_2.radius * ceq_2.radius)) / (2 * (ceq_1.corY - ceq_2.corY));
			float root_a = ceq_1.corX - sqrt((ceq_1.radius * ceq_1.radius) - (temp_corY - ceq_1.corY) * (temp_corY - ceq_1.corY));
			float root_b = ceq_1.corX + sqrt((ceq_1.radius * ceq_1.radius) - (temp_corY - ceq_1.corY) * (temp_corY - ceq_1.corY));
			if ((root_a >= ceq_1.corX_min) && (root_a <= ceq_1.corX_max) && (root_a >= ceq_1.corX_min) && (root_a <= ceq_1.corX_max)) {
				save->cross_type[save->ptr] = POSSIBLE;
				save->cross_corX[save->ptr] = root_a;
				save->cross_corY[save->ptr] = temp_corY;
				save->cross_eq1[save->ptr] = CIRCULAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_temp;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_2;
				save->ptr++;
			}
			else {
				/* void */
			}

			if ((root_b >= ceq_1.corX_min) && (root_b <= ceq_1.corX_max) && (root_b >= ceq_2.corX_min) && (root_b <= ceq_2.corX_max)) {
				save->cross_type[save->ptr] = POSSIBLE;
				save->cross_corX[save->ptr] = root_b;
				save->cross_corY[save->ptr] = temp_corY;
				save->cross_eq1[save->ptr] = CIRCULAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_temp;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_2;
				save->ptr++;
			}
			else {
				/* void */
			}

		}
		else if ((ceq_1.corX != ceq_2.corX) && (ceq_1.corY != ceq_2.corY)) {
			float temp_gradient = (ceq_2.corX - ceq_1.corX) / (ceq_1.corY - ceq_2.corY); // The gradient of the line is 0;
			float temp_basis = (((ceq_1.corX * ceq_1.corX) - (ceq_2.corX * ceq_2.corX)) + ((ceq_1.corY * ceq_1.corY) - (ceq_2.corY * ceq_2.corY)) + ((ceq_2.radius * ceq_2.radius) - (ceq_1.radius * ceq_1.radius))) / (2 * (ceq_1.corY - ceq_2.corY));
			float temp_a = temp_gradient * temp_gradient + 1;
			float temp_b = 2 * (temp_gradient * (temp_basis - ceq_1.corY) - ceq_1.corX);
			float temp_c = (ceq_1.corX * ceq_1.corX) + ((temp_basis - ceq_1.corY) * (temp_basis - ceq_1.corY)) - (ceq_1.radius * ceq_1.radius);
			float determinate = (temp_b * temp_b) - (4 * temp_a * temp_c);

			// Impossible
			if (determinate < 0) {
				/*
				save->cross_type[save->ptr] = IMPOSSIBLE;
				save->cross_corX[save->ptr] = 0;
				save->cross_eq1[save->ptr] = CIRCULAR;
				save->cross_eq2[save->ptr] = CIRCULAR;
				save->cross_leq1[save->ptr] = leq_temp;
				save->cross_leq2[save->ptr] = leq_temp;
				save->cross_ceq1[save->ptr] = ceq_1;
				save->cross_ceq2[save->ptr] = ceq_2;
				save->ptr++;
				*/
			}
			// Dual root (In between the two circle, Dual root has no meaning for overlapping)
			else if (determinate == 0) {
				float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / (2 * temp_a);
				float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / (2 * temp_a);

				if (root_a == root_b) {
					save->cross_type[save->ptr] = DUALROOT;
					save->cross_corX[save->ptr] = root_a;
					save->cross_corY[save->ptr] = root_a * temp_gradient + temp_basis;
					save->cross_eq1[save->ptr] = CIRCULAR;
					save->cross_eq2[save->ptr] = CIRCULAR;
					save->cross_leq1[save->ptr] = leq_temp;
					save->cross_leq2[save->ptr] = leq_temp;
					save->cross_ceq1[save->ptr] = ceq_1;
					save->cross_ceq2[save->ptr] = ceq_2;
					save->ptr++;
				}
				else {
					printf("DUAL ROOT ERR... (CEQ and CEQ) \n");
					system("PAUSE");
				}
			}
			else if (determinate > 0) {
				float root_a = (((-1.0) * temp_b) - sqrt(determinate)) / (2 * temp_a);
				float root_b = (((-1.0) * temp_b) + sqrt(determinate)) / (2 * temp_a);
				float tempA_corY = root_a * temp_gradient + temp_basis;
				float tempB_corY = root_b * temp_gradient + temp_basis;
				if ((root_a >= ceq_1.corX_min) && (root_a <= ceq_1.corX_max) && (root_a >= ceq_2.corX_min) && (root_a <= ceq_2.corX_max)
					&& (tempA_corY >= ceq_1.corY_min) && (tempA_corY <= ceq_1.corY_max) && (tempA_corY >= ceq_2.corY_min) && (tempA_corY <= ceq_2.corY_max)) {
					save->cross_type[save->ptr] = POSSIBLE;
					save->cross_corX[save->ptr] = root_a;
					save->cross_corY[save->ptr] = root_a * temp_gradient + temp_basis;
					save->cross_eq1[save->ptr] = CIRCULAR;
					save->cross_eq2[save->ptr] = CIRCULAR;
					save->cross_leq1[save->ptr] = leq_temp;
					save->cross_leq2[save->ptr] = leq_temp;
					save->cross_ceq1[save->ptr] = ceq_1;
					save->cross_ceq2[save->ptr] = ceq_2;
					save->ptr++;
				}
				else {
					/* void */
				}

				if ((root_b >= ceq_1.corX_min) && (root_b <= ceq_1.corX_max) && (root_b >= ceq_2.corX_min) && (root_b <= ceq_2.corX_max)
					&& (tempB_corY >= ceq_1.corY_min) && (tempB_corY <= ceq_1.corY_max) && (tempB_corY >= ceq_2.corY_min) && (tempB_corY <= ceq_2.corY_max)) {
					save->cross_type[save->ptr] = POSSIBLE;
					save->cross_corX[save->ptr] = root_b;
					save->cross_corY[save->ptr] = root_b * temp_gradient + temp_basis;
					save->cross_eq1[save->ptr] = CIRCULAR;
					save->cross_eq2[save->ptr] = CIRCULAR;
					save->cross_leq1[save->ptr] = leq_temp;
					save->cross_leq2[save->ptr] = leq_temp;
					save->cross_ceq1[save->ptr] = ceq_1;
					save->cross_ceq2[save->ptr] = ceq_2;
					save->ptr++;
				}
				else {
					/* void */
				}
			}
		}
	}
	else if (det_ceq == 0) {
		// Pre-condition : congruence
		if ((ceq_1.corX == ceq_2.corX) && (ceq_1.corY == ceq_2.corY) && (ceq_1.radius == ceq_2.radius)) {
			/*
			save->cross_type[save->ptr] = CONGRUENCE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = CIRCULAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_temp;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_2;
			save->ptr++;
			*/
		}
		// Pre-condition : Contain (is equal to impossible)
		else if ((ceq_1.corX == ceq_2.corX) && (ceq_1.corY == ceq_2.corY)) {
			/*
			save->cross_type[save->ptr] = IMPOSSIBLE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = CIRCULAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_temp;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_2;
			save->ptr++;
			*/
		}
		else {
			printf("CIRCULAR CONDITION ERR (DET_CEQ == 0, BUT COORDINATE IS NOT SAME...\n");
			// Pre-condition : Non-overlapping
			/*
			save->cross_type[save->ptr] = IMPOSSIBLE;
			save->cross_corX[save->ptr] = 0;
			save->cross_eq1[save->ptr] = CIRCULAR;
			save->cross_eq2[save->ptr] = CIRCULAR;
			save->cross_leq1[save->ptr] = leq_temp;
			save->cross_leq2[save->ptr] = leq_temp;
			save->cross_ceq1[save->ptr] = ceq_1;
			save->cross_ceq2[save->ptr] = ceq_2;
			save->ptr++;
			*/
		}
	}
	else if (det_ceq > (ceq_1.radius + ceq_2.radius)) {
		// Pre-condition : Non-overlapping
		/*
		save->cross_type[save->ptr] = IMPOSSIBLE;
		save->cross_corX[save->ptr] = 0;
		save->cross_eq1[save->ptr] = CIRCULAR;
		save->cross_eq2[save->ptr] = CIRCULAR;
		save->cross_leq1[save->ptr] = leq_temp;
		save->cross_leq2[save->ptr] = leq_temp;
		save->cross_ceq1[save->ptr] = ceq_1;
		save->cross_ceq2[save->ptr] = ceq_2;
		save->ptr++;
		*/
	}
	else {
		printf("CIRCULAR CONDITON ERR (DET_CEQ)...\n");
		/*
		save->cross_type[save->ptr] = IMPOSSIBLE;
		save->cross_corX[save->ptr] = 0;
		save->cross_eq1[save->ptr] = CIRCULAR;
		save->cross_eq2[save->ptr] = CIRCULAR;
		save->cross_leq1[save->ptr] = leq_temp;
		save->cross_leq2[save->ptr] = leq_temp;
		save->cross_ceq1[save->ptr] = ceq_1;
		save->cross_ceq2[save->ptr] = ceq_2;
		save->ptr++;
		*/
	}
}


void cal(SU su1, SU su2, db save[]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			// SU #1의 i번째 Beam 자가교점
			findCross(su1.beam[i].fleq, su1.beam[i].sleq, &save[(4 * i) + j]);
			findCross(su1.beam[i].fleq, su1.beam[i].tceq, &save[(4 * i) + j]);
			findCross(su1.beam[i].sleq, su1.beam[i].tceq, &save[(4 * i) + j]);
			// SU #2의 j번째 Beam 자가교점
			findCross(su2.beam[j].fleq, su2.beam[j].sleq, &save[(4 * i) + j]);
			findCross(su2.beam[j].fleq, su2.beam[j].tceq, &save[(4 * i) + j]);
			findCross(su2.beam[j].sleq, su2.beam[j].tceq, &save[(4 * i) + j]);
			// SU #1의 i번째 Beam 중 fleq 및 SU #2의 j번째 Beam 교점
			findCross(su1.beam[i].fleq, su2.beam[j].fleq, &save[(4 * i) + j]);
			findCross(su1.beam[i].fleq, su2.beam[j].sleq, &save[(4 * i) + j]);
			findCross(su1.beam[i].fleq, su2.beam[j].tceq, &save[(4 * i) + j]);
			// SU #1의 i번째 Beam 중 sleq 및 SU #2의 j번째 Beam 교점
			findCross(su1.beam[i].sleq, su2.beam[j].fleq, &save[(4 * i) + j]);
			findCross(su1.beam[i].sleq, su2.beam[j].sleq, &save[(4 * i) + j]);
			findCross(su1.beam[i].sleq, su2.beam[j].tceq, &save[(4 * i) + j]);
			// SU #1의 i번째 Beam 중 tceq 및 SU #2의 j번째 Beam 교점
			findCross(su2.beam[j].fleq, su1.beam[i].tceq, &save[(4 * i) + j]);
			findCross(su2.beam[j].sleq, su1.beam[i].tceq, &save[(4 * i) + j]);
			findCross(su1.beam[i].tceq, su2.beam[j].tceq, &save[(4 * i) + j]);

			// 교점을 구한 다음 공통 범위 내에 존재하는 것 외에는 전부 지움.
			for (int k = 0; k < save[(4 * i) + j].ptr; k++) {
				float x_min = max(su1.beam[i].x_min, su2.beam[j].x_min);
				float x_max = min(su1.beam[i].x_max, su2.beam[j].x_max);
				float y_min = max(su1.beam[i].y_min, su2.beam[j].y_min);
				float y_max = min(su1.beam[i].y_max, su2.beam[j].y_max);

				if ((save[(4 * i) + j].cross_corX[k] > x_max) || (save[(4 * i) + j].cross_corX[k] < x_min) || (save[(4 * i) + j].cross_corY[k] > y_max) || (save[(4 * i) + j].cross_corY[k] < y_min)) {
					for (int l = k; l < save[(4 * i) + j].ptr; l++) {
						if (l < BUFFERSIZE) {
							save[(4 * i) + j].cross_ceq1[l] = save[(4 * i) + j].cross_ceq1[l + 1];
							save[(4 * i) + j].cross_ceq2[l] = save[(4 * i) + j].cross_ceq2[l + 1];
							save[(4 * i) + j].cross_leq1[l] = save[(4 * i) + j].cross_leq1[l + 1];
							save[(4 * i) + j].cross_leq2[l] = save[(4 * i) + j].cross_leq2[l + 1];
							save[(4 * i) + j].cross_eq1[l] = save[(4 * i) + j].cross_eq1[l + 1];
							save[(4 * i) + j].cross_eq2[l] = save[(4 * i) + j].cross_eq2[l + 1];
							save[(4 * i) + j].cross_type[l] = save[(4 * i) + j].cross_type[l + 1];
							save[(4 * i) + j].cross_corX[l] = save[(4 * i) + j].cross_corX[l + 1];
							save[(4 * i) + j].cross_corY[l] = save[(4 * i) + j].cross_corY[l + 1];
						}
						else if (l == BUFFERSIZE) {
							save[(4 * i) + j].cross_ceq1[l] = ceq_temp;
							save[(4 * i) + j].cross_ceq2[l] = ceq_temp;
							save[(4 * i) + j].cross_leq1[l] = leq_temp;
							save[(4 * i) + j].cross_leq2[l] = leq_temp;
							save[(4 * i) + j].cross_eq1[l] = -1;
							save[(4 * i) + j].cross_eq2[l] = -1;
							save[(4 * i) + j].cross_type[l] = -1;
							save[(4 * i) + j].cross_corX[l] = 0.0;
							save[(4 * i) + j].cross_corY[l] = 0.0;
						}
						else {
							printf("INSERTING ERR....\n");
							system("PAUSE");
						}
						l--;
						save[(4 * i) + j].ptr--;
					}
				}
			}

			sortCoordinate(&save[(4 * i) + j]);
		}
	}
}