// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include "CImg.h"
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <Sift.h>
#include <tgmath.h>

//Use the cimg namespace to access functions easily
using namespace cimg_library;
using namespace std;
#define N 3

void draw_descriptor_image(CImg<double> image, const vector<SiftDescriptor> descriptors, const char *filename)
{
  for(unsigned int i=0; i < descriptors.size(); i++)
    {
      int tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
      double color_point[] = {255.0, 255.0, 0};
      for(int x=-2; x<3; x++)
	for(int y=-2; y<3; y++)
	  if(x==0 || y==0)
	    for(int c=0; c<3; c++){
	      //Find if coordinates are in workspace to draw crosshair
	      tx1 = (descriptors[i].col + y - 1);
	      ty1 = (descriptors[i].row + x - 1);
	      if (tx1 >= 0 && tx1 < image.width() && ty1 >= 0 && ty1 < image.height())
		image( tx1, ty1, 0, c) = color_point[c];				
	    }
    }
  image.get_normalize(0,255).save(filename);
}

void getCofactor(double A[N][N], double temp[N][N], int p, int q, int n)
{
    int i = 0, j = 0;
 
    // Looping for each element of the matrix
    for (int row = 0; row < n; row++)
    {
        for (int col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i][j++] = A[row][col];
 
                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

double determinant(double A[N][N], int n)
{
    double D = 0; // Initialize result
 
    //  Base case : if matrix contains single element
    if (n == 1)
        return A[0][0];
 
    double temp[N][N]; // To store cofactors
 
    int sign = 1;  // To store sign multiplier
 
     // Iterate for each element of first row
    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0][f] * determinant(temp, n - 1);
 
        // terms are to be added with alternate sign
        sign = -sign;
    }
 
    return D;
}

void adjoint(double A[N][N],double adj[N][N])
{
    if (N == 1)
    {
        adj[0][0] = 1;
        return;
    }
 
    // temp is used to store cofactors of A[][]
    double sign = 1, temp[N][N];
 
    for (int i=0; i<N; i++)
    {
        for (int j=0; j<N; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, N);
 
            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;
 
            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j][i] = (sign)*(determinant(temp, N-1));
        }
    }
}

void inverseMatrix(double A[N][N], double inverse[N][N])
{
    // Find determinant of A[][]
    double det = determinant(A, N);
    if (det == 0)
    {
        cout << "Singular matrix, can't find its inverse";
    }
 
    // Find adjoint
    double adj[N][N];
    adjoint(A, adj);
 
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
            inverse[i][j] = adj[i][j]/double(det);
}

void getWarpedImage(CImg<double> input_image,double projectiveTransform[][3])
{
	double inverse[3][3];
	inverseMatrix(projectiveTransform,inverse);
	CImg<double> output_image(1024,1024);
	for(int i=0;i<input_image.width();i++)
	{
		for(int j=0;j<input_image.height();j++)
		{
			float y = (inverse[0][0] * i + inverse[0][1] *j + inverse[0][2]*1);
			float x = (inverse[1][0] * i + inverse[1][1] *j + inverse[1][2]*1);
			float w = (inverse[2][0] * i + inverse[2][1] *j + inverse[2][2]*1);
//			cout<< y <<' '<<endl;
			if((x/w > 0) && (x/w<1024) && (y/w>0) && (y/w<1024))	
				output_image(i,j) = input_image(int(y/w),int(x/w));
		}
	}


	output_image.save("output.png");

}

CImg<double> getTransformationMatrix(int x1,int y1,int x2,int y2,int x3,int y3,int x4,int y4,int X1,int Y1,int X2,int Y2,int X3,int Y3,int X4,int Y4)
{
	CImg<double> A(8,8);
	A(0,0) = x1; A(0,1) = y1; A(0,2) = 1; A(0,3) = 0; A(0,4) = 0; A(0,5) = 0; A(0,6) = -x1*X1; A(0,7) = -y1*X1;
	
	A(1,0) = 0; A(1,1) = 0; A(1,2) = 0; A(1,3) = x1; A(1,4) = y1; A(1,5) = 1; A(1,6) = -x1*Y1; A(1,7) = -y1*Y1;

	A(2,0) = x2; A(2,1) = y2; A(2,2) = 1; A(2,3) = 0; A(2,4) = 0; A(2,5) = 0; A(2,6) = -x2*X2; A(2,7) = -y2*X2;

	A(3,0) = 0; A(3,1) = 0; A(3,2) = 0; A(3,3) = x2; A(3,4) = y2; A(3,5) = 1; A(3,6) = -x2*Y2; A(3,7) = -y2*Y2;
	
	A(4,0) = x3; A(4,1) = y3; A(4,2) = 1; A(4,3) = 0; A(4,4) = 0; A(4,5) = 0; A(4,6) = -x3*X3; A(4,7) = -y3*X3;
	
	A(5,0) = 0; A(5,1) = 0; A(5,2) = 0; A(5,3) = x3; A(5,4) = y3; A(5,5) = 1; A(5,6) = -x3*Y3; A(5,7) = -y3*Y3;

	A(6,0) = x4; A(6,1) = y4; A(6,2) = 1; A(6,3) = 0; A(6,4) = 0; A(6,5) = 0; A(6,6) = -x4*X4; A(6,7) = -y4*X4;
	
	A(7,0) = 0; A(7,1) = 0; A(7,2) = 0; A(7,3) = x4; A(7,4) = y4; A(7,5) = 1; A(7,6) = -x4*Y4; A(7,7) = -y4*Y4;

	CImg<double> B(8,1);

	B(0,0) = X1;
	B(1,0) = Y1;	
	B(2,0) = X2;
	B(3,0) = Y2;
	B(4,0) = X3;
	B(5,0) = Y3;
	B(6,0) = X4;
	B(7,0) = Y4;

//	CImg<double> transformMatrix(9,1) = B.solve(A);
	CImg<double>::B.solve(A);
	for(int i=0;i<B.width();i++)
		for(int j=0;j<B.height();j++)
			cout<<B(i,j)<<' '<<endl;

	return A;
}
int main(int argc, char **argv)
{
  try {
    
    /*
      TEST CODE - STARTS
    */
	string part = "";
	CImg<double> input_image("images/part1/lincoln.png");
	CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
//	vector<SiftDescriptor> input_descriptors = Sift::compute_sift(input_gray);
//	draw_descriptor_image(input_image, input_descriptors, "input_image.png");
	double projectiveTransform[3][3];
	double inverse[3][3];
	projectiveTransform[0][0] = 0.907;
	projectiveTransform[0][1] = 0.258;
	projectiveTransform[0][2] = -182;
	projectiveTransform[1][0] = -0.153;
	projectiveTransform[1][1] = 1.44;
	projectiveTransform[1][2] = 58;
	projectiveTransform[2][0] = 0.000306;
	projectiveTransform[2][1] = 0.000731;
	projectiveTransform[2][2] = 1;

	inverse[0][0] = 1.124;
	inverse[0][1] = -0.3146;
	inverse[0][2] = 222.95;
	inverse[1][0] = 0.1088;
	inverse[1][1] = 0.685;
	inverse[1][2] = -19.92;
	inverse[2][0] = 0.00026;
	inverse[2][1] = -0.000597;
	inverse[2][2] = 1.0827;

	inverseMatrix(projectiveTransform,inverse);
/*
	for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			cout<<inverse[i][j]<<' '<<endl;	
*/
	CImg<double> output_image(1024,1024);
	for(int i=0;i<input_image.width();i++)
	{
		for(int j=0;j<input_image.height();j++)
		{
			float y = (inverse[0][0] * i + inverse[0][1] *j + inverse[0][2]*1);
			float x = (inverse[1][0] * i + inverse[1][1] *j + inverse[1][2]*1);
			float w = (inverse[2][0] * i + inverse[2][1] *j + inverse[2][2]*1);
//			cout<< y <<' '<<endl;
			if((x/w > 0) && (x/w<1024) && (y/w>0) && (y/w<1024))	
				output_image(i,j) = input_image(int(y/w),int(x/w));
		}
	}


	output_image.save("output.png");
	getTransformationMatrix(318,256,534,372,316,670,73,473,141,131,480,159,493,630,64,601);	
			

    /*
      TEST CODE - ENDS
    */
    
    if(part == "part1"){
      // Billboard
    }	
    else if(part == "part2"){
      // Blending
    }
    else if(part == "part3"){
      // RANSAC
    }
    else if(part == "part4"){
      // Panorama
    }
    
    
    // feel free to add more conditions for other parts (e.g. more specific)
    //  parts, for debugging, etc.
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}
