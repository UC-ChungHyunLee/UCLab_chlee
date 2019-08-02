#include "lib.h"


int main(void) {

	SU testSU1;
	testSU1.index = 1;
	testSU1.corX = 2.0;
	testSU1.corY = 2.0;
	testSU1.radius = 2.0;
	SU testSU2;
	testSU2.index = 2;
	testSU2.corX = 4.0;
	testSU2.corY = 4.0;
	testSU2.radius = 2.0;
	init_SU(&testSU1);
	init_SU(&testSU2);

	db save[16];
	for (int i = 0; i < 16; i++) {
		init_DB(&save[i]);
	}
	/*
	cal(testSU1, testSU2, save);
	for (int i = 0; i < 4; i++) {
		printf("SU #1�� %d��° Beam &\n ", i);
		for (int j = 0; j < 4; j++) {
			printf("\tSU #2�� %d��° Beam�� ���� : ", j);
			for (int k = 0; k < save[(4 * i) + j].ptr; k++) {
				printf(" %.4f, ", save[(4 * i) + j].cross_corX[k]);
			}
			printf("\n");
		}
	}
	*/
	
	save[0].cross_corX[0] = 2;
	save[0].cross_corY[0] = 4;
	save[0].cross_corX[1] = 2.5858;
	save[0].cross_corY[1] = 2.5858;

	float temp = integral(testSU1.beam->tceq, testSU2.beam->tceq, save[0], 0, 1);

	printf("%.4f\n", temp);
	
	/*
	printf("#%d SU�� #%d SU�� ���� �׽�Ʈ..\n", testSU1.index, testSU2.index);
	for (int i = 0; i < save.ptr; i++) {
		printf("���� %d�� :\n", i + 1);
		switch (save.cross_eq1[i]) {
		case LINEAR:
			switch (save.cross_eq2[i]) {
			case LINEAR:
				printf("\t#%d SU�� %d��° ��ä���� %d��° ���� - ���� %d, ���� %d\n", save.cross_leq1[i].root, save.cross_leq1[i].index, save.cross_leq1[i].num, (int)save.cross_leq1[i].gradient, (int)save.cross_leq1[i].basis);
				printf("\t#%d SU�� %d��° ��ä���� %d��° ���� - ���� %d, ���� %d\n", save.cross_leq2[i].root, save.cross_leq2[i].index, save.cross_leq2[i].num, (int)save.cross_leq2[i].gradient, (int)save.cross_leq2[i].basis);
				switch (save.cross_type[i]) {
				case POSSIBLE:
					printf("\t\t %.4f���� ������ ����\n\n", save.cross_corX[i]);
					break;
				case IMPOSSIBLE:
					printf("\t\t ������ ������ ����\n\n");
					break;
				case INDETERMINATE:
					printf("\t\t �յ���\n\n");
					break;
				default:
					printf("\t\t ������ �߻�\n\n");
					break;
				}
				break;
			case CIRCULAR:
				printf("\t#%d SU�� %d��° ��ä���� %d��° ���� - ���� %d, ���� %d\n", save.cross_leq1[i].root, save.cross_leq1[i].index, save.cross_leq1[i].num, (int)save.cross_leq1[i].gradient, (int)save.cross_leq1[i].basis);
				printf("\t#%d SU�� %d��° ��ä���� ȣ - �� %.4f���� �� %.4f\n", save.cross_ceq1[i].root, save.cross_ceq1[i].index, save.cross_ceq1[i].angle_min, save.cross_ceq1[i].angle_max);
				switch (save.cross_type[i]) {
				case POSSIBLE:
				printf("\t\t %.4f���� ������ ����\n\n", save.cross_corX[i]);
					break;
				case DUALROOT:
					printf("\t\t ���߱��� ����\n\n");
					break;
				case IMPOSSIBLE:
					printf("\t\t ������ ������ ����\n\n");
					break;
				default:
					printf("\t\t ������ �߻�\n\n");
					break;
				}
				break;
			default:
				printf("\t��� ���� �߻�\n\n");
				break;
			}
			break;
		case CIRCULAR:
			switch (save.cross_eq2[i]) {
			case CIRCULAR:
				printf("\t#%d SU�� %d��° ��ä���� ȣ - �� %.4f���� �� %.4f\n", save.cross_ceq1[i].root, save.cross_ceq1[i].index, save.cross_ceq1[i].angle_min, save.cross_ceq1[i].angle_max);
				printf("\t#%d SU�� %d��° ��ä���� ȣ - �� %.4f���� �� %.4f\n", save.cross_ceq2[i].root, save.cross_ceq2[i].index, save.cross_ceq2[i].angle_min, save.cross_ceq2[i].angle_max);
				switch (save.cross_type[i]) {
				case POSSIBLE:
					printf("\t\t %.4f���� ������ ����\n\n", save.cross_corX[i]);
					break;
				case CONGRUENCE:
					printf("\t\t ȣ�� ��Ȯ�ϰ� ��ħ, �� ������ SU 2���� �������� ������\n\n");
					break;
				case IMPOSSIBLE:
					printf("\t\t ������ ������ ����\n\n");
					break;
				default:
					printf("\t\t ������ �߻�\n\n");
					break;
				}
				break;
			default:
				printf("\tEquation 1�� ���̸� Equation 2 ���� ���� ���ۿ� �����Ƿ� ������\n");
				break;
			}
			break;
		default:
			printf("\t ���е� �ƴϰ� ���� �ƴ�, �� ������\n\n");
			break;
		}
	}
	*/
	system("PAUSE");

	return 0;
}