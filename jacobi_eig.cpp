#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
using namespace std;
#include "jacobi_eig.h"

int sgn(double val) 
{
		if (val > 0) { return 1; }else if (val < 0) { return -1;}else{return 0;} 
}


void get_jk_matrix(float **matrix /*float matrix*/, int P) {

	double Sin2, Cos2, Cos_, Sin_, Tan2, Cot2, tmp;
	double den = 0.0, num = 0.0;
	double sumsqr, prev_sumsqr = 0.0;
	const double eps = 1E-16;

	for (int count = 1; count <= 15; count++) {
		for (int j = 0; j < P - 1; j++) {
			for (int k = j + 1; k < P; k++) {
				num = 0.0;
				den = 0.0;
				for (int i = 0; i < P; i++) {
					num += 2 * matrix[i][j] * matrix[i][k];
					den += (matrix[i][j] + matrix[i][k]) * (matrix[i][j] - matrix[i][k]);
				}
				if (abs(num) < eps && den >= 0) {
					break; //exit for loop
				}
				// find sin(2*theta) and cos(2*theta)
				if (abs(num) <= abs(den)) {    // 2*sigma(xy) <= |sqrt(x2-y2)|
					Tan2 = abs(num) / abs(den);
					Cos2 = 1 / sqrt(1 + Tan2 * Tan2);
					Sin2 = Tan2 * Cos2;
				}
				else {
					Cot2 = abs(den) / abs(num);
					Sin2 = 1 / sqrt(1 + Cot2 * Cot2);
					Cos2 = Cot2 * Sin2;
				}
				// condition +/- for cos(theta) and sin(theta)
				Cos_ = sqrt((1 + Cos2) / 2);
				Sin_ = Sin2 / (2 * Cos_);
				if (den < 0) {
					tmp = Cos_;
					Cos_ = Sin_;
					Sin_ = tmp;
				}
				Sin_ = sgn(num) * Sin_;
				// perform rotation on matrix
				for (int i = 0; i < P; i++) {
					tmp = matrix[i][j];
					matrix[i][j] = tmp * Cos_ + matrix[i][k] * Sin_;
					matrix[i][k] = -tmp * Sin_ + matrix[i][k] * Cos_;
				}
			}// NEXT k
		}// NEXT j

		sumsqr = 0;
		for (int c = 0; c < P; c++) {
			sumsqr += matrix[0][c] * matrix[c][0];
		}
		// display some results from the iteration.
		if (count > 5) {
			cout << "\t\t" << "---------------- " << "count:" << count << " ---------------" << endl;
			cout << "\t\t" << "prev_sumsqr = " << prev_sumsqr << ";\n\t\t" << "sumsqr = " << sumsqr << ";\n";
			cout << "\t\t" << "Error:" << abs(sumsqr - prev_sumsqr) << endl;
		}
		// if error is small break loop else continue
		if (abs(sumsqr - prev_sumsqr) < eps) {
			break;
		}
		else { prev_sumsqr = sumsqr; }
	}// NEXT count
}

/*********************** Define Member Functions for Class myStack *****************/
void myStack::push(const StackElement & value) {
	if (myTop < STACK_CAPACITY - 1)
	{
		myArray[++myTop] = value;
	}
	else
	{
		cerr << "Stack is full! \n";
	}
}

void myStack::pop()
{
	if (myTop >= 0)
		myTop--;
	else
		cerr << "Stack is empty!  Cannot remove top! \n";
}


void myStack::display(ostream & out) const
{
	for (int i = myTop; i >= 0; i--)
	{
		out << myArray[i] << endl;
	}
}

StackElement myStack::top() const
{
	if (myTop >= 0)
	{
		return myArray[myTop];
	}
	else { cerr << " Stack is Empty ! \n"; return -1; }
}