#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#define _USE_MATHE_DEFINES

#define BUFFERSIZE 10000

#define POSSIBLE 0
#define IMPOSSIBLE 1
#define INDETERMINATE 2

typedef struct dataSave {
	int ptr_corX; // Pointer of coordinate which crossed between equation_1 and equation_2
	int ptr_eq1; // Pointer of SU which is one of the crossed (equation_1)
	int ptr_eq2; // Pointer of SU which is one of the crossed (equation_2)
	int ptr_bool; // Pointer of relationship between equation_1 and equation_2 such as having a root, impossible and indeterminate
	float cross_corX[BUFFERSIZE];
	float cross_eq1[BUFFERSIZE][BUFFERSIZE];
	float cross_eq2[BUFFERSIZE][BUFFERSIZE];
	bool cross_type[BUFFERSIZE];
}db;

typedef struct LinearEquation {
	float gradient;
	float basis;
	float corX_min;
	float corX_max;
} leq;

typedef struct CircularEquation {
	float corX;
	float corY;
	float radius;
	float angle_min;
	float angle_max;
} ceq;

typedef struct Sector {
	float ch; // Channel State with Probability
	int index; // 0~3
	leq fleq; // First Linear Equation
	leq sleq; // Second Linear Equation
	ceq tceq; // Third Circular Equation
} sec;

typedef struct PrimaryUser {
	float corX;
	float corY;
	float radius;
	int flag;
} PU;

typedef struct SecondaryUser {
	int index;
	float corX;
	float corY;
	float radius;

	sec beam[4];
};
