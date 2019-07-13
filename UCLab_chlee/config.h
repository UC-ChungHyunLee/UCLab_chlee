#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#define _USE_MATHE_DEFINES

#define BUFFERSIZE 15

#define POSSIBLE 0 // 해를 가짐 (선분과 선분이면 1개, 선분과 호면 2개, 호와 호면 2개)
#define DUALROOT 1 // 이중근 (선분과 호 또는 호와 호에서 이중근을 가질 수 있음)
#define IMPOSSIBLE 2 // 불능 (해가 없음, 겹치지 않음)
#define INDETERMINATE 3 // 부정 (선분이 겹침)
#define CONGRUENCE 4 // 합동 (호가 겹침)

#define LINEAR 0
#define CIRCULAR 1

typedef struct LinearEquation {
	float gradient;
	float basis;
	float corX_min;
	float corX_max;
	float corY_min;
	float corY_max;
	int root;
	int index;
	int num;
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
	int root;
	int index;
} ceq;

typedef struct Sector {
	float ch; // Channel State with Probability
	int root; // SU index
	int index; // 0~3
	leq fleq; // First Linear Equation
	leq sleq; // Second Linear Equation
	ceq tceq; // Third Circular Equation
	float x_min;
	float x_max;
	float y_min;
	float y_max;
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
} SU;

typedef struct dataSave {
	int ptr;
	float cross_corX[BUFFERSIZE];
	float cross_corY[BUFFERSIZE];
	int cross_eq1[BUFFERSIZE];
	int cross_eq2[BUFFERSIZE];
	int cross_type[BUFFERSIZE];
	leq cross_leq1[BUFFERSIZE];
	leq cross_leq2[BUFFERSIZE];
	ceq cross_ceq1[BUFFERSIZE];
	ceq cross_ceq2[BUFFERSIZE];
}db;

leq leq_temp;
ceq ceq_temp;