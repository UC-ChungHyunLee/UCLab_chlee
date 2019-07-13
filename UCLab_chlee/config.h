#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#define _USE_MATHE_DEFINES

#define BUFFERSIZE 15

#define POSSIBLE 0 // �ظ� ���� (���а� �����̸� 1��, ���а� ȣ�� 2��, ȣ�� ȣ�� 2��)
#define DUALROOT 1 // ���߱� (���а� ȣ �Ǵ� ȣ�� ȣ���� ���߱��� ���� �� ����)
#define IMPOSSIBLE 2 // �Ҵ� (�ذ� ����, ��ġ�� ����)
#define INDETERMINATE 3 // ���� (������ ��ħ)
#define CONGRUENCE 4 // �յ� (ȣ�� ��ħ)

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