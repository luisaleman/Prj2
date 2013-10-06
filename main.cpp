#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstdio>
#include <sstream>
#include <fstream>
#include <string>
#include "time.h"


using namespace std;
/*I declare prototype to avoid the compiler believing that the functions declared
at the end don't exist*/
void jacobi_solver (double **A, double **R, int n);
void rotation(double **A, double **R, int k, int l, int n);
double maxoffdiagonal(double **A, int *k, int *l, int n);
void output(double *A, double *rho, double **eigenvector, int r, int n);
void quicksort(double *array, double **extra, int start, int end);
double divide(double *array, double **extra, int start, int end);

int main(int argc, char* argv[]){

    double rho_max,h, total_time;
    int number_steps,dimension, i,j;
    double omega;
    omega=5;
    cout << endl << "Please enter rho max: ";
    cin >> rho_max;
    cout << endl << "Please enter number of steps: ";
    cin >> number_steps;
    //Set the step length "h"
    h=(double)(rho_max)/(double)(number_steps);
    cout << endl<< "The step leght was: " << h << endl;
    dimension=number_steps-1;

    //Initialize rho
    double *rho;
    rho= new double [dimension];
    //Set up rho, we go from n=1 to n=nmax-1
    for(i=0;i<dimension;i++){
        rho[i]=h*(i+1);
    }

    //Initialize the system matrix
    double **system_matrix;
    system_matrix= new double*[dimension];
    for(i=0;i<dimension;i++){
        system_matrix[i]= new double[dimension];
        }
    //Set up the system matrix
    for(i=0;i<dimension;i++){
        system_matrix[i][i]=(2/pow(h,2))+pow(omega,2)*pow(rho[i],2)+(1/rho[i]);
        }
    for(i=0;i<dimension-1;i++){
        system_matrix[i][i+1]=(-1/pow(h,2));
        system_matrix[i+1][i]=(-1/pow(h,2));
        }

    //Initialize eigenvector matrix
    double **eigen_matrix;
    eigen_matrix=new double*[dimension];
    for(i=0;i<dimension;i++){
        eigen_matrix[i]=new double [dimension];
    }

    //To know the elapsed time for the jacobi solver
    clock_t start, finish;
    start = clock();
    //Solve the problem
    jacobi_solver(system_matrix, eigen_matrix, dimension);
    finish = clock();
    total_time=((finish-start)/(double)CLOCKS_PER_SEC);
    cout << "The elapsed time was " << total_time<< " seg"<< endl;

    /*Due to the fact that the eigenvalues are not ordered, we have to do so
    *just to make it easy to put the into a file or plot them
    *also the eigenvectors wouldn't be ordered.
    */
    double *eigenvalues;
    eigenvalues= new double[dimension];
    for(i=0;i<dimension;i++){
        eigenvalues[i]=system_matrix[i][i];
    }
    for(i=0;i<3;i++){
        printf("%i eigenvalue is %.4f \n", i+1, eigenvalues[i]);
    }
    quicksort(eigenvalues,eigen_matrix, 0,dimension-1);
    //Output into a file
    output(eigenvalues, rho, eigen_matrix ,rho_max,number_steps);
    for(i=0;i<3;i++){
        printf("%i eigenvalue is %.4f \n", i+1, eigenvalues[i]);
    }
    return 0;
}
#pragma endregion

#pragma region jacobi
/*Function of the jacobi solver, it uses both functions rotation and
 *maxoffdiagonal
*/
void jacobi_solver (double **A, double **R, int n)
{
    int k, l;
    double tolerance = 1.0e-10;                              //set the tolerance for the off-diagonal elements
    double max_iterations = (double)n*(double)n*(double)n;               //to set a limit on the number of iterations
    int iterations = 0;                                     //to count the iterations and set the other limit
    double max_off_diag;
    max_off_diag= maxoffdiagonal(A, &k, &l, n);

    // We need to create the matrix to store the eigenvectors (nxn)
    //at the begining of the process is the identity matrix
    for (int i=0;i<n;i++){
        for(int j= 0; j<n;j++){
            if(i==j){
                R[i][j]=1.0;
            }
            else{
                R[i][j]=0.0;
            }
        }
    }

    while (fabs(max_off_diag)>tolerance && (double)iterations < max_iterations){
        rotation (A, R, k , l , n);
        max_off_diag = maxoffdiagonal(A, &k, &l, n);
        iterations++;
}
    cout << "The number of iterations was: "<< iterations << endl;

    return ;
}
#pragma endregion

#pragma region rotation
/*This function computes tau, tan, sin and cos using the values given k and l
 *and also implements the similarity transformation to the matrix elements
 *and keeps track on the eigenvectors in the last part. So we can easily have both
 *eigenvectors and eigenvalues
*/
void rotation(double **A, double **R, int k, int l, int n){
    double sin,cos;
    if ( A[k][l]!= 0.0){
        double tan, tau;
        tau= (A[l][l]-A[k][k])/(2*A[k][l]);
        //we have to choose the smallest value for tan, the best way is:
        //if tau is positive we have to sum the sqrt to -tau
        if ( tau >0){
            tan=-tau + sqrt(1+pow(tau,2));
        }
        //if tau is negative we have to substract the sqrt to -tau
        else{
            tan=-tau - sqrt(1+pow(tau,2));
        }
        cos= 1/sqrt(1 + pow(tan,2));
        sin=cos*tan;
    }
    else{
        cos=1.0;
        sin=0.0;
    }
    //we proceed to implement the transformation to the matrix A
    double a_kk, a_ll, a_ik, a_il,a_kl, r_ik, r_il;
    a_kk = A[k][k];
    a_ll = A[l][l];
    a_kl = A[k][l];
    A[k][k]=pow(cos,2)*a_kk - 2*cos*sin*A[k][l] + pow(sin,2)*a_ll;
    A[l][l]=pow(sin,2)*a_kk + 2*cos*sin*A[k][l] + pow(cos,2)*a_ll;
    A[k][l]=0.0;//(a_kk-a_ll)*cos*sin+a_kl*(pow(cos,2)-pow(sin,2));
    A[l][k]=A[k][l];
    for(int i= 0; i <n; i++){
        if(i != k && i != l){
            a_ik=A[i][k];
            a_il=A[i][l];
            A[i][k]=cos*a_ik - sin*a_il;
            A[k][i]=A[i][k];
            A[i][l]=cos*a_il + sin*a_ik;
            A[l][i]=A[i][l];
        }
        //we compute the new aigenvectors
        r_ik=R[i][k];
        r_il=R[i][l];
        R[i][k]=cos*r_ik - sin*r_il;
        R[i][l]=cos*r_il + sin*r_ik;
    }
    return;
}
#pragma endregion

#pragma region maxoffdiagonal
/* The function to get the maximum element off the diagonal and get its coordinates
 *"l" "k"
*/
double maxoffdiagonal(double **A, int *k, int *l, int n){
    double max_value=0.00;
    for(int i=0; i<n; i++){
        for (int j=i+1;j<n;j++){
            if (max_value<fabs(A[i][j])){
                max_value=fabs(A[i][j]);
                *l=i;
                *k=j;
            }
        }
    }
    return max_value;
}
#pragma endregion


#pragma region output
/*The function output brings into a file the data of (in the next order):
 *rho[i], the first three eigenvectors, the square of the eigenvectors, the eigenvalues
 *
 *NOTICE:the eigenvalues were sorted in ascending order and so the eigenvectors
*/
void output(double *A, double *rho, double **eigenvector, int r, int n){
    //char *outputname;
    string outputname;
    stringstream str2, str1;
    str1 << r;
    str2 << n;
    outputname = str1.str();
    outputname += "rho_";
    outputname += str2.str();
    outputname += "steps.txt";

    ofstream myfile;
    myfile.open(outputname.c_str());
    for (int i=0;i<n-1;i++){
        myfile << setprecision(8) <<rho[i]<<setw(20);
        myfile << setprecision(8) <<eigenvector[i][0]<<setw(15);
        myfile << setprecision(8) <<eigenvector[i][1]<<setw(15);
        myfile << setprecision(8) <<eigenvector[i][2]<<setw(15);
        myfile << setprecision(8) <<eigenvector[i][0]*eigenvector[i][0]<<setw(15);
        myfile << setprecision(8) <<eigenvector[i][1]*eigenvector[i][1]<<setw(15);
        myfile << setprecision(8) <<eigenvector[i][2]*eigenvector[i][2]<<setw(15);
        myfile << setprecision(8)<<A[i]<< endl ;
    }
    myfile.close();
}
#pragma endregion

#pragma region quicksort
/* Recursive function to do the sorting in ascending order due that the eigenvalues
 *sometimes changed of place with others, and also, we sort the eigenvectors just to have
 *that the first eigenvalue corresponds to the first column, etc.
 */
void quicksort(double *array,double **extra, int start, int end)
{
    double pivot;

    if (start < end) {
        pivot = divide(array, extra, start, end);

        // Ordeno la lista de los menores
        quicksort(array,extra, start, pivot - 1);

        // Ordeno la lista de los mayores
        quicksort(array, extra, pivot + 1, end);
    }
}
// Function to divide the array and do the interchange
double divide(double *array, double **extra, int start, int end) {
    int left;
    int right;
    double pivot;
    double temp;
    double *temp2;
    temp2=new double[end+1];

    pivot = array[start];
    left = start;
    right = end;

    // while the indexes don't cross with each other
    while (left < right) {
        while (array[right] > pivot) {
            right--;
        }

        while ((left < right) && (array[left] <= pivot)) {
            left++;
        }

        // if they havent crossed, we continue
        if (left < right) {
            temp = array[left];
            array[left] = array[right];
            array[right] = temp;
            for(int i=0;i<end+1;i++){
                temp2[i]=extra[i][left];
                extra[i][left]=extra[i][right];
            extra[i][right]=temp2[i];
            }
        }
    }

    // if the indexes have already crossed, we stop and mark the pivot
    temp = array[right];
    array[right] = array[start];
    array[start] = temp;
    for(int i=0;i<end+1;i++){
    temp2[i]= extra[i][right];
    extra[i][right]=extra[i][start];
    extra[i][start]=temp2[i];
    }

    // the new position of the pivot
    return right;
}
#pragma endregion

