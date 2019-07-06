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
	
	printf("#%d SU와 #%d SU간 교점 테스트..\n", testSU1.index, testSU2.index);
	for (int i = 0; i < save.ptr; i++) {
		printf("교점 %d번 :", i + 1);
		switch (save.cross_eq1[i]) {
		case LINEAR:
			printf("\t선분 - 기울기 %d, 절편 %d\n", (int)save.cross_leq1[i].gradient, (int)save.cross_leq1[i].basis);
			switch (save.cross_eq2[i]) {
			case LINEAR:
				printf("\t선분 - 기울기 %d, 절편 %d\n", (int)save.cross_leq2[i].gradient, (int)save.cross_leq2[i].basis);
				printf("#%d SU의 %d번째 부채꼴의 %d번째 선분과 #%d SU의 %d번째 부채꼴의 %d번째 선분간 교점 체크\n", save.cross_leq1[i].root, save.cross_leq1[i].index, save.cross_leq1[i].num, save.cross_leq2[i].root, save.cross_leq2[i].index, save.cross_leq2[i].num);
				switch (save.cross_type[i]) {
				case POSSIBLE:
					printf("\t\t %.4f에서 교점을 가짐\n\n", save.cross_corX[i]);
					break;
				case IMPOSSIBLE:
					printf("\t\t 교점을 가지지 않음\n\n");
					break;
				case INDETERMINATE:
					printf("\t\t 합동임\n\n");
					break;
				default:
					printf("\t\t 계산오류 발생\n\n");
					break;
				}
				break;
			case CIRCULAR:
				printf("\t호 - 각 %.4f에서 각 %.4f\n", save.cross_ceq1[i].angle_min, save.cross_ceq1[i].angle_max);
				printf("#%d SU의 %d번째 부채꼴의 %d번째 선분과 #%d SU의 %d번째 부채꼴의 호 간 교점 체크\n", save.cross_leq1[i].root, save.cross_leq1[i].index, save.cross_leq1[i].num, save.cross_ceq1[i].root, save.cross_ceq1[i].index);
				switch (save.cross_type[i]) {
				case POSSIBLE:
				printf("\t\t %.4f에서 교점을 가짐\n\n", save.cross_corX[i]);
					break;
				case DUALROOT:
					printf("\t\t 이중근을 가짐\n\n");
					break;
				case IMPOSSIBLE:
					printf("\t\t 교점을 가지지 않음\n\n");
					break;
				default:
					printf("\t\t 계산오류 발생\n\n");
					break;
				}
				break;
			default:
				printf("\t계산 오류 발생\n\n");
				break;
			}
			break;
		case CIRCULAR:
			printf("\t호 - 각 %.4f에서 각 %.4f\n", save.cross_ceq1[i].angle_min, save.cross_ceq1[i].angle_max);
			switch (save.cross_eq2[i]) {
			case CIRCULAR:
				printf("#%d SU의 %d번째 부채꼴의 호와 #%d SU의 %d번째 부채꼴의 호간 교점 체크\n", save.cross_ceq1[i].root, save.cross_ceq1[i].index, save.cross_ceq2[i].root, save.cross_ceq2[i].index);
				printf("\t호 - 각 %.4f에서 각 %.4f\n", save.cross_ceq2[i].angle_min, save.cross_ceq2[i].angle_max);
				switch (save.cross_type[i]) {
				case POSSIBLE:
					printf("\t\t %.4f에서 교점을 가짐\n\n", save.cross_corX[i]);
					break;
				case CONGRUENCE:
					printf("\t\t 호가 정확하게 겹침, 즉 동일한 SU 2개가 비교함으로 계산오류\n\n");
					break;
				case IMPOSSIBLE:
					printf("\t\t 교점을 가지지 않음\n\n");
					break;
				default:
					printf("\t\t 계산오류 발생\n\n");
					break;
				}
				break;
			default:
				printf("\tEquation 1이 원이면 Equation 2 또한 원인 경우밖에 없으므로 계산오류\n");
				break;
			}
			break;
		default:
			printf("\t 선분도 아니고 원도 아님, 즉 계산오류\n\n");
			break;
		}
	}

	system("PAUSE");

	return 0;
}