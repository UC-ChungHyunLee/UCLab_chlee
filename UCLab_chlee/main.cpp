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

	db save;
	init_DB(&save);
	
	cal(testSU1, testSU2, &save);
	
	printf("#%d SU�� #%d SU�� ���� �׽�Ʈ..\n", testSU1.index, testSU2.index);
	for (int i = 0; i < save.ptr; i++) {
		printf("���� %d�� :", i + 1);
		switch (save.cross_eq1[i]) {
		case LINEAR:
			printf("\t���� - ���� %d, ���� %d\n", (int)save.cross_leq1[i].gradient, (int)save.cross_leq1[i].basis);
			switch (save.cross_eq2[i]) {
			case LINEAR:
				printf("\t���� - ���� %d, ���� %d\n", (int)save.cross_leq2[i].gradient, (int)save.cross_leq2[i].basis);
				printf("#%d SU�� %d��° ��ä���� %d��° ���а� #%d SU�� %d��° ��ä���� %d��° ���а� ���� üũ\n", save.cross_leq1[i].root, save.cross_leq1[i].index, save.cross_leq1[i].num, save.cross_leq2[i].root, save.cross_leq2[i].index, save.cross_leq2[i].num);
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
				printf("\tȣ - �� %.4f���� �� %.4f\n", save.cross_ceq1[i].angle_min, save.cross_ceq1[i].angle_max);
				printf("#%d SU�� %d��° ��ä���� %d��° ���а� #%d SU�� %d��° ��ä���� ȣ �� ���� üũ\n", save.cross_leq1[i].root, save.cross_leq1[i].index, save.cross_leq1[i].num, save.cross_ceq1[i].root, save.cross_ceq1[i].index);
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
			printf("\tȣ - �� %.4f���� �� %.4f\n", save.cross_ceq1[i].angle_min, save.cross_ceq1[i].angle_max);
			switch (save.cross_eq2[i]) {
			case CIRCULAR:
				printf("#%d SU�� %d��° ��ä���� ȣ�� #%d SU�� %d��° ��ä���� ȣ�� ���� üũ\n", save.cross_ceq1[i].root, save.cross_ceq1[i].index, save.cross_ceq2[i].root, save.cross_ceq2[i].index);
				printf("\tȣ - �� %.4f���� �� %.4f\n", save.cross_ceq2[i].angle_min, save.cross_ceq2[i].angle_max);
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

	system("PAUSE");

	return 0;
}