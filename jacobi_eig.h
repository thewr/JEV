#ifndef jacobi_eig_h
#define jacobi_eig_h

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdlib>

using namespace std;



// get sign 
//int sgn(double val){ 
//	if (val > 0) { return 1; }else if (val < 0) { return -1;}else{return 0;} 
//}

//function to display 2D array
//void display(float **matrix,int rows, int cols){
//	for (int r = 0; r < rows; ++r){
//		for (int c = 0; c < cols; ++c){
//			cout << right << setw(8) << setfill(' ') << setprecision(3);
//			cout << matrix[r][c] << " ";
//
//		}cout << endl;
//	}cout << endl;
//}


////-------------------------------------------------------------
//// Functions
typedef vector <double> record_t;
typedef vector <record_t> data_t;
////-----------------------------------------------------------------------------
//// Let's overload the stream input operator to read a list of CSV fields (a CSV record).
//// Remember, a record is a list of doubles separated by commas ','.
istream& operator >> (istream & ins, record_t & record)
{
	// make sure that the returned record contains only the stuff we read now
	record.clear();

	//	// read the entire line into a string (a CSV record is terminated by a newline)
	string line;
	getline(ins, line);

	//	// now we'll use a stringstream to separate the fields out of the line
	stringstream ss(line);
	string field;
	while (getline(ss, field, ','))
	{
		// for each field we wish to convert it to a double
		// (since we require that the CSV contains nothing but floating-point values)
		stringstream fs(field);
		double f = 0.0;  // (default value is 0.0)
		fs >> f;

		// add the newly-converted field to the end of the record
		record.push_back(f);
	}

	// Now we have read a single line, converted into a list of fields, converted the fields
	// from strings to doubles, and stored the results in the argument record, so
	// we just return the argument stream as required for this kind of input overload function.
	return ins;
}

//-----------------------------------------------------------------------------
// Let's likewise overload the stream input operator to read a list of CSV records.
// This time it is a little easier, just because we only need to worry about reading
// records, and not fields.
istream& operator >> (istream& ins, data_t& data)
{
	// make sure that the returned data only contains the CSV data we read here
	data.clear();

	// For every record we can read from the file, append it to our resulting data
	record_t record;
	while (ins >> record)
	{
		data.push_back(record);
	}

	// Again, return the argument stream as required for this kind of input stream overload.
	return ins;
}


//------------------------------------------------------------
// Apply JK method to get eigen values/vectors from corr matrix
//void get_jk_matrix(float **matrix /*float matrix*/, int P) {
//
//	double Sin2, Cos2, Cos_, Sin_, Tan2, Cot2, tmp;
//	double den = 0.0, num = 0.0;
//	double sumsqr, prev_sumsqr = 0.0;
//	const double eps = 1E-16;
//
//	for (int count = 1; count <= 15; count++){
//		for (int j = 0; j < P - 1; j++)	{
//			for (int k = j + 1; k < P; k++){
//				num = 0.0;
//				den = 0.0;
//				for (int i = 0; i < P; i++){
//					num += 2 * matrix[i][j] * matrix[i][k];
//					den += (matrix[i][j] + matrix[i][k]) * (matrix[i][j] - matrix[i][k]);
//				}
//				if (abs(num) < eps && den >= 0)	{
//					break; //exit for loop
//				}
//				// find sin(2*theta) and cos(2*theta)
//				if (abs(num) <= abs(den)){    // 2*sigma(xy) <= |sqrt(x2-y2)|
//					Tan2 = abs(num) / abs(den);
//					Cos2 = 1 / sqrt(1 + Tan2 * Tan2);
//					Sin2 = Tan2 * Cos2;
//				}
//				else{
//					Cot2 = abs(den) / abs(num);
//					Sin2 = 1 / sqrt(1 + Cot2 * Cot2);
//					Cos2 = Cot2 * Sin2;
//				}
//				// condition +/- for cos(theta) and sin(theta)
//				Cos_ = sqrt((1 + Cos2) / 2);
//				Sin_ = Sin2 / (2 * Cos_);
//				if (den < 0){
//					tmp = Cos_;
//					Cos_ = Sin_;
//					Sin_ = tmp;
//				}
//				Sin_ = sgn(num) * Sin_;
//				// perform rotation on matrix
//				for (int i = 0; i < P; i++){
//					tmp = matrix[i][j];
//					matrix[i][j] = tmp * Cos_ + matrix[i][k] * Sin_;
//					matrix[i][k] = -tmp * Sin_ + matrix[i][k] * Cos_;
//				}
//			}// NEXT k
//		}// NEXT j
//
//		sumsqr = 0;
//		for (int c = 0; c < P; c++){
//			sumsqr += matrix[0][c] * matrix[c][0];
//		}
//		// display some results from the iteration.
//		if (count > 5) {
//			cout << "\t\t" << "---------------- " << "count:" << count << " ---------------" << endl;
//			cout << "\t\t" << "prev_sumsqr = " << prev_sumsqr << ";\n\t\t" << "sumsqr = " << sumsqr << ";\n";
//			cout << "\t\t" << "Error:" << abs(sumsqr - prev_sumsqr) << endl;
//		}
//		// if error is small break loop else continue
//		if (abs(sumsqr - prev_sumsqr) < eps){
//			break;}else{prev_sumsqr = sumsqr;}
//	}// NEXT count
//
//
//
//
//}


//h for myStack
typedef double StackElement;
const int STACK_CAPACITY = 30;


class myStack
{
public:
	/* --- Constructor ---
	// Precondition:  myStack has been declared
	// Postcondition:  myStack has been constructed as an empty stack
	******************************************************************/
	myStack() { myTop = -1; };
	bool empty() const { return{ myTop == -1 }; };
	void push(const StackElement & value);
	void display(ostream & out) const;
	StackElement top() const;
	void pop();

	/* Data Members
	******************************/
private:
	StackElement myArray[STACK_CAPACITY];
	int myTop;
};

//
//void myStack::push(const StackElement & value) {
//	if (myTop < STACK_CAPACITY - 1)
//	{
//		myArray[++myTop] = value;
//	}
//	else
//	{
//		cerr << "Stack is full! \n";
//	}
//}
//
//void myStack::pop()
//{
//	if (myTop >= 0)
//		myTop--;
//	else
//		cerr << "Stack is empty!  Cannot remove top! \n";
//}
//
//
//void myStack::display(ostream & out) const
//{
//	for (int i = myTop; i >= 0; i--)
//	{
//		out << myArray[i] << endl;
//	}
//}
//
//StackElement myStack::top() const
//{
//	if (myTop >= 0)
//	{
//		return myArray[myTop];
//	}
//	else { cerr << " Stack is Empty ! \n"; }
//}

#endif