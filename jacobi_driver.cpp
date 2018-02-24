#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cstdlib>
using namespace std;
#include "jacobi_eig.h"


typedef vector<double> vector_1d;
typedef vector < vector_1d > vector_2d;





//---------------------------------------------------------------------------
// We allocate, load from files, and process arrays.
// In the end, we get the jacobi matrix eigenvalues/vectors.
// Let's also print to files.
int main()
{
	// the data we want.
	data_t data_matrix;

	// here is the file containing the data. Read into data.
	ifstream infile("data.csv");
	infile >> data_matrix;

	// if error reading file.
	if (!infile.eof())
	{
		cout << "Error reading file!\n";
		return 1;
	}

	infile.close();

	// basic information about the file.
	cout << "Your CSV file contains " << data_matrix.size() << " records.\n";
	int ROWS = data_matrix.size();   //ROWS OF MATRIX
	unsigned max_record_size = 0;
	for (unsigned n = 0; n < data_matrix.size(); n++)
		if (max_record_size < data_matrix[n].size())
			max_record_size = data_matrix[n].size();
	cout << "The largest record has " << max_record_size << " fields.\n\n";
	int COLS = max_record_size;
	int N = COLS;

	for (unsigned i = 0; i < data_matrix.size(); i++)
	{
		for (unsigned j = 0; j < max_record_size; j++)
		{
			data_matrix[i][j] = data_matrix[i][j];
		}
	} // data_matrix loaded



	// allocate memory for 2D array
	// array: data_matrix_T[COLS][ROWS]
	float** data_matrix_T = new float*[COLS];
	for (int i = 0; i < COLS; ++i)
		data_matrix_T[i] = new float[ROWS];
	// array: cov_matrix[N][N]
	float** cov_matrix = new float*[N];
	for (int i = 0; i < N; ++i)
		cov_matrix[i] = new float[N];
	// array: corr_matrix[N][N]
	float **corr_matrix = new float*[N];
	for (int i = 0; i < N; i++) {
		corr_matrix[i] = new float[N];
	}


	// find transpose of data_matrix.
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < COLS; j++) {
			data_matrix_T[j][i] = data_matrix[i][j];
		}
	}

	vector_1d cov_diag(N);
	// find covariance matrix.
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			cov_matrix[j][i] = 0.0;
			for (int k = 0; k < ROWS; k++) {
				cov_matrix[j][i] += data_matrix_T[j][k] * data_matrix[k][i] / (ROWS - 1);
			}
		}
	}

	// display matrix
	//cout << "\n" << "Covariance Matrix" << endl;
	//display(cov_matrix, N, N);

	// find correlation matrix from covariance matrix.
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			corr_matrix[i][j] = cov_matrix[i][j] / (sqrt(cov_matrix[i][i])*sqrt(cov_matrix[j][j]));
		}
	}

	// display matrix
	//cout << "\n" << "Correlation Matrix" << endl;
	//display(corr_matrix, N, N);

	// ---------------------------------------------------
	// Free each sub-array
	// deallocate data_matrix_T
	for (int i = 0; i < COLS; ++i) {
		delete[] data_matrix_T[i];
	}	delete[] data_matrix_T;

	// deallocate corr_matrix
	for (int i = 0; i < N; ++i) {
		delete[] corr_matrix[i];
	}	delete[] corr_matrix;


	//----------------------------------------------------
	// do jakobi method to eig
	get_jk_matrix(cov_matrix, N);

	//-------------------------------------------------
	// get eigen values and vectors
	vector_1d Eig_Vals(N);
	vector_2d X;
	//now we have an empty 2D-matrix of size (0,0). Resizing it with one single command:
	X.resize(N /*COLUMNS*/, vector_1d(2 /*ROWS*/, 0.0 /*INIT VAL*/));
	// and we are good to go .



	// find eigen values
	for (int j = 0; j < N; j++)
	{
		for (int k = 0; k < N; k++)
		{
			Eig_Vals[j] += (cov_matrix[k][j])*(cov_matrix[k][j]);

		}
		Eig_Vals[j] = sqrt(Eig_Vals[j]);
	}

	// find eigen vectors
	for (int i = 0; i < 2; i++) //column
	{
		for (int j = 0; j < N; j++) //row
		{
			X[j][i] = cov_matrix[j][i] / Eig_Vals[i];
		}
	}



	ofstream out;
	out.open("eig_vectors.txt");
	string Name[4] = { "A","B","C","D" };
	for (int i = 0;i < 2;i++) {
		out << "Cells(k, feature" << Name[i] << ").Value = ";
		for (int j = 0; j < N; j++)
		{
			if (j < (N - 1)) {
				out << X[j][i] << "*" << Name[j] << " +";
			}
			else
			{
				out << X[j][i] << "*" << Name[j];
			}
		}  out << endl;
	}out.close();

	cout << "\nEig Values: ";
	for (int i = 0;i < N;i++)
		cout << Eig_Vals[i] << " ";
	cout << endl;


	cout << "\nE.Vector(#1): \n";
	for (int i = 0;i < N;i++) {
		for (int j = 0; j < 2; j++) {
			cout << right << setw(8) << setfill(' ') << setprecision(3);
			cout << X[i][j] << " ";
		}cout << endl;
	}

	//---------------------------------------------------------------------------------
	// calculate and store in 2D vector the components data_matrix



	float sum = 0.0;
	const int num_clusters = 7;
	vector_2d f_cmp_Vals; f_cmp_Vals.resize(ROWS, vector_1d(N, 0.0));
	vector_2d f_clust_Vals; f_clust_Vals.resize(num_clusters, vector_1d(2, 0.0));
	vector_2d f_clust_Dist; f_clust_Dist.resize(ROWS, vector_1d(num_clusters, 0.0));
	vector_2d w; w.resize(ROWS, vector_1d(num_clusters, 0.0)); //w[i][j]

	double min[2] = { 0.0, 0.0 };
	double max[2] = { 0.0, 0.0 };




	for (int rows = 0; rows < ROWS; rows++) {
		for (int k = 0; k < 2; k++) {
			sum = 0.0;
			for (int i = 0; i < N; i++) {
				sum += data_matrix[rows][i] * X[i][k];
			}
			f_cmp_Vals[rows][k] = sum;


			//find min & max of feature components of data_matrix
			if (rows == 0)
			{
				min[k] = f_cmp_Vals[rows][k];
				max[k] = f_cmp_Vals[rows][k];
			}
			else if (f_cmp_Vals[rows][k] < min[k])
			{
				min[k] = f_cmp_Vals[rows][k];
			}
			else if (f_cmp_Vals[rows][k] > max[k])
			{
				max[k] = f_cmp_Vals[rows][k];
			}
		}// Next k
	}// Next rows


	int i = 0;
	int k = 0;


	//--------------------------------------------------
	//  update centroids
	sum = 0.0;
	k = 0;
	int j = 0;
	int count;
	float den, num;
	const int M = 2.0;

	den = 0.0;
	num = 0.0;
	
	// generate initial clusters
	f_clust_Vals[0][0] = min[0]; f_clust_Vals[0][1] = min[1];
	for (int i = 1; i < num_clusters;i++) {
		count = 0;
		while (count < 2) {
			f_clust_Vals[i][count] = f_clust_Vals[i - 1][count] + (max[count] - min[count]) / (num_clusters);
			count++;
		}
	}

	//========= display clusters =============== 
	cout << "\n" << "\t\t x\t y" << endl;
	for (int rows = 0; rows < num_clusters; rows++) {
		cout << "cluster #" << rows << ":";
		for (int k = 0; k < 2; k++) {
			cout << right << setw(8) << setfill(' ') << setprecision(3);
			cout << f_clust_Vals[rows][k] << " ";
		}cout << endl;
	}
	//==========================================


	// get distance
	for (int rows = 0; rows < ROWS; rows++) {
		for (int k = 0; k < num_clusters; k++) {
			count = 0;
			while (count < 2) {
				f_clust_Dist[rows][k] += pow((f_cmp_Vals[rows][count] - f_clust_Vals[k][count]), 2.0);
				count++;
			}
		}
	}

	//cout << "\n" << "components:" << endl;
	//for (int rows = 15; rows < 25; rows++) {
	//	cout << "Player_#" << rows + 1 << ": ";
	//	for (int k = 0; k < num_clusters; k++) {
	//		cout << f_clust_Dist[rows][k] << " ";
	//	} cout << endl;
	//}


	cout << "\n\n";
	myStack s;

	int rank;
	double value;
	vector_1d rank_all(ROWS);


	// get ranks
	for (int rows = 0; rows < ROWS; rows++) {
		for (int k = num_clusters - 1; k >= 0; k--) s.push(f_clust_Dist[rows][k]);
		value = s.top();
		s.pop();
		rank = 0;
		while (!s.empty())
		{
			if (value > s.top())
			{
				value = s.top();
				s.pop();
				rank++;
			}
			else {
				s.pop();
			}
		}
		rank_all[rows] = rank;
	}

	//cout << "\n" << "\t\t x\t y" << endl;
	for (int rows = 0; rows < num_clusters; rows++) {
		cout << "cluster #" << rows << ":";
		for (int k = 0; k < 2; k++) {
			cout << right << setw(8) << setfill(' ') << setprecision(3);
			cout << f_clust_Vals[rows][k] << " ";
		}cout << endl;
	}

	// w[i][j] => num = 1; den = <sigm_c_k=1> [ ||xi - cj|| / ||xi-ck||)^(2/(m-1) ];
	// c[k] => num = <sigma_x> w[k] * (x)^m * x; den = <sigma_x> w[k] * (x)^m;

	
	//--------------- ITERATION 1 --------------------------------
	// get weights for data
	den = 0.0;
	num = 0.0;
	for (int i = 0; i < ROWS; i++) {
		for (int j = 0; j < num_clusters;j++) {
			sum = 0.0;
			k = 0;
			num = pow(f_clust_Dist[i][j], 1 / (M - 1));
			while (k < num_clusters)
			{
				
				den = pow(f_clust_Dist[i][k], 1 / (M - 1) ); //* f_clust_Dist[i][k]; //m = 2;
				sum += num / den;
				k++;
			}
				w[i][j] = 1 / sum;
		}
	}

	//// display weights
	cout << "\nweights:\n";
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < num_clusters; j++)
		{
			cout << right << setw(8) << setfill(' ') << setprecision(2);
			cout << w[i][j] << " ";
		}
		cout << endl;
	}


	//--------------------------------------------------
	//  update centroids
	k = 0;
	j = 0;

	den = 0.0;
	num = 0.0;
	while (k < num_clusters) // O O O O ...
	{
		j = 0;
		while (j < 2)  // |X|Y| ...  
		{
			num = 0.0;
			den = 0.0;
			for (int i = 0; i < ROWS; i++) //ROWS
			{
				num += pow(w[i][k], M) * (f_cmp_Vals[i][j]);
				den += pow(w[i][k], M);			
			}
			sum = num / den;
			f_clust_Vals[k][j] = sum;
			j++;
		}
		
		k++;	
	}

	//========= display clusters =============== 
	//cout << "\n" << "\t\t x\t y" << endl;
	//for (int rows = 0; rows < num_clusters; rows++) {
	//	cout << "cluster #" << rows << ":";
	//	j = 0;
	//	while (j < 2)
	//	{
	//		cout << right << setw(8) << setfill(' ') << setprecision(3);
	//		cout << f_clust_Vals[rows][j] << " ";
	//		j++;
	//	}
	//	cout << endl;
	//}
	//============== END OF ITERATION ===========


	float hold = 0.0;
	for (int index = 0; index < 50; index++) {
		// get weights for data

		// get distance
		for (int rows = 0; rows < ROWS; rows++) {
			for (int k = 0; k < num_clusters; k++) {
				count = 0;
				while (count < 2) {
					f_clust_Dist[rows][k] += pow((f_cmp_Vals[rows][count] - f_clust_Vals[k][count]), 2.0);
					count++;
				}
			}
		}



		for (int i = 0; i < ROWS; i++) {
			for (int j = 0; j < num_clusters;j++) {
				sum = 0.0;
				k = 0;
				num = pow(f_clust_Dist[i][j], 1 / ( M - 1 ));
				while (k < num_clusters)
				{

					den = pow(f_clust_Dist[i][k], 1 / (M - 1)); //* f_clust_Dist[i][k]; //m = 2;
					sum += num / den;
					k++;
				}

				w[i][j] = 1 / sum;

			}

		}

		//--------------------------------------------------
		//  update centroids
		sum = 0.0;
		k = 0;
		j = 0;

		while (k < num_clusters) // O O O O ...
		{
			j = 0;
			while (j < 2)  // |X|Y| ...  
			{
				num = 0.0;
				den = 0.0;
				for (int i = 0; i < ROWS; i++) //ROWS
				{
					num += pow(w[i][k], M) * (f_cmp_Vals[i][j]);
					den += pow(w[i][k], M);
				}
				sum = num / den;
				f_clust_Vals[k][j] = sum;
				j++;
			}

			k++;
		}

		sum = 0.0;
		for (int i = 0; i < ROWS; i++) {
				
			for (int j = 0; j < num_clusters; j++) {
				sum += pow(w[i][j], M)*f_clust_Dist[i][j];
			}
		}
			
		if (abs(sum - hold) < 1) {break;} else { hold = sum; }

	} // end index

	//========== display weights ===============
	cout << "\nweights:\n";
	for (int i = 0; i < 10; i++)
	{
		for (int j = 0; j < num_clusters; j++)
		{
			cout << right << setw(8) << setfill(' ') << setprecision(2);
			cout << w[i][j] << " ";
		}
		cout << endl;
	}
	//==========================================


	//========= display clusters =============== 
	cout << "\n" << "\t\t x\t y" << endl;
	for (int rows = 0; rows < num_clusters; rows++) {
		cout << "cluster #" << rows << ":";
		j = 0;
		while (j < 2)
		{
			cout << right << setw(8) << setfill(' ') << setprecision(3);
			cout << f_clust_Vals[rows][j] << " ";
			j++;
		}
		cout << endl;
	}


	for (int rows = 0; rows < ROWS; rows++) {
		for (int k = num_clusters - 1; k >= 0; k--) s.push(f_clust_Dist[rows][k]);
		value = s.top();
		s.pop();
		rank = 0;
		while (!s.empty())
		{
			if (value > s.top())
			{
				value = s.top();
				s.pop();
				rank++;
			}
			else {
				s.pop();
			}
		}
		rank_all[rows] = rank;
	}





	system("pause");
	return 0;

}