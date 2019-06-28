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
#define DUALROOT 1
#define IMPOSSIBLE 2
#define INDETERMINATE 3
#define CONGRUENCE 4

#define LINEAR 0
#define CIRCULAR 1

typedef struct dataSave {
	int ptr;
	float cross_corX[BUFFERSIZE];
	float cross_eq1[BUFFERSIZE];
	float cross_eq2[BUFFERSIZE];
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
	float corX_min;
	float corX_max;
	float corY_min;
	float corY_max;
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
