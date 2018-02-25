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
// #include <tgmath.h>
#include <math.h>
#include <float.h>

//Use the cimg namespace to access functions easily
using namespace cimg_library;
using namespace std;
#define N 3

struct data_tuple
{
    int feature_point_1;
    int feature_point_2;
    double short_distance;
};

struct line{
    int x_1,y_1;
    int x_2,y_2;
};
CImg <double> draw_descriptor_image(CImg<double> image, const vector<SiftDescriptor> descriptors, const char *filename)
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
	      {	
		image( tx1, ty1, 0, c) = color_point[c];		
	      }				
	    }
    }
  image.get_normalize(0,255).save(filename);
  return image;
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
	
	cimg_forXY(input_image,i,j)
	{
		float x = (inverse[0][0] * i + inverse[0][1] *j + inverse[0][2]*1);
		float y = (inverse[1][0] * i + inverse[1][1] *j + inverse[1][2]*1);
		float w = (inverse[2][0] * i + inverse[2][1] *j + inverse[2][2]*1);
		if((x/w > 0) && (x/w<1024) && (y/w>0) && (y/w<1024))	
			output_image(i,j) = input_image(int(x/w),int(y/w));
	}

	output_image.save("output.png");

}

CImg<double> getTransformationMatrix(int x1,int y1,int x2,int y2,int x3,int y3,int x4,int y4,int X1,int Y1,int X2,int Y2,int X3,int Y3,int X4,int Y4)
{
	CImg<long double> A(8,8);
	A(0,0) = x1; A(0,1) = y1; A(0,2) = 1; A(0,3) = 0; A(0,4) = 0; A(0,5) = 0; A(0,6) = -x1*X1; A(0,7) = -y1*X1;
	
	A(1,0) = 0; A(1,1) = 0; A(1,2) = 0; A(1,3) = x1; A(1,4) = y1; A(1,5) = 1; A(1,6) = -x1*Y1; A(1,7) = -y1*Y1;

	A(2,0) = x2; A(2,1) = y2; A(2,2) = 1; A(2,3) = 0; A(2,4) = 0; A(2,5) = 0; A(2,6) = -x2*X2; A(2,7) = -y2*X2;

	A(3,0) = 0; A(3,1) = 0; A(3,2) = 0; A(3,3) = x2; A(3,4) = y2; A(3,5) = 1; A(3,6) = -x2*Y2; A(3,7) = -y2*Y2;
	
	A(4,0) = x3; A(4,1) = y3; A(4,2) = 1; A(4,3) = 0; A(4,4) = 0; A(4,5) = 0; A(4,6) = -x3*X3; A(4,7) = -y3*X3;
	
	A(5,0) = 0; A(5,1) = 0; A(5,2) = 0; A(5,3) = x3; A(5,4) = y3; A(5,5) = 1; A(5,6) = -x3*Y3; A(5,7) = -y3*Y3;

	A(6,0) = x4; A(6,1) = y4; A(6,2) = 1; A(6,3) = 0; A(6,4) = 0; A(6,5) = 0; A(6,6) = -x4*X4; A(6,7) = -y4*X4;
	
	A(7,0) = 0; A(7,1) = 0; A(7,2) = 0; A(7,3) = x4; A(7,4) = y4; A(7,5) = 1; A(7,6) = -x4*Y4; A(7,7) = -y4*Y4;

	CImg<long double> B(1,8);
//	cout<<B.width()<<endl;
	CImg<long double> X(8,1);
	B(0,0) = X1;
	B(0,1) = Y1;	
	B(0,2) = X2;
	B(0,3) = Y2;
	B(0,4) = X3;
	B(0,5) = Y3;
	B(0,6) = X4;
	B(0,7) = Y4;

//	CImg<double> transformMatrix(8,1) = A.solve(B);
	X = B.solve(A.transpose());
/*	cimg_forXY(X,i,j)
	{
		cout<<X(i,j)<<endl;
	}
*/
	return X;
}

CImg<double> drawLines(CImg<double> image,std::vector<line> lineVector){
    const unsigned char color[] = { 255,128,64 };
    // input_image.draw_line(40,40,80,70,color);
    for(int i=0;i<lineVector.size();i++){
        image.draw_line(lineVector[i].x_1,lineVector[i].y_1,lineVector[i].x_2,lineVector[i].y_2,color);
    }
    return image;

}

bool comparator_function(data_tuple obj_1,data_tuple obj_2){
    return (obj_1.short_distance<obj_2.short_distance);
}

double getEuclideanDistance(std::vector<float> datapoint_1,std::vector<float> datapoint_2){
double distance=0.0;
// cout<<datapoint_1.size()<<","<<datapoint_2.size()<<endl;
for(int i=0;i<int(datapoint_1.size());i++){
    distance+=(datapoint_1[i]-datapoint_2[i])*(datapoint_1[i]-datapoint_2[i]);
    // cout<<"("<<datapoint_1[i]<<","<<datapoint_2[i]<<")"<<endl;
}
// cout<<"D:"<<sqrtf(distance)<<endl;
return sqrt(distance);

}

CImg<double> stcichImages(CImg<double> image_1,CImg<double> image_2){
int max_width=image_1.width();
int max_height=image_1.height();

CImg<double> stiched_image(2*max_width,max_height,1,3);
for(int i=0;i<max_height;i++){
    for(int j=0;j<max_width;j++){
        for(int c=0;c<3;c++){
            stiched_image(j,i,0,c)=image_1(j,i,0,c);
        }
    }
}

for(int i=0;i<max_height;i++){
    for(int j=0;j<max_width;j++){
        for(int c=0;c<3;c++){
            stiched_image(j+max_width,i,0,c)=image_2(j,i,0,c);
        }
    }
}
return stiched_image;
}

CImg<double> getMatchingWithBruteForce(CImg<double> input_image_1,CImg<double> input_image_2,double threshold=0.65) {

// input_image_1=input_image_1.resize(input_image_2.width(),input_image_2.height());
// input_image_1.save("resized.png");
// temp.append_object3d(input_image_2).save("merged.png");
    cout<<"C:"<<input_image_1(1,1,0,1)<<endl;
    int max_width,max_height;
    if(input_image_1.height()<=input_image_2.height()){
        max_height=input_image_2.height();
    }
    else{
        max_height=input_image_1.height();
    }
     if(input_image_1.width()<=input_image_2.width()){
        max_width=input_image_2.width();
     }
     else{
        max_width=input_image_1.width();
     }
input_image_1=input_image_1.resize(max_width,max_height);
input_image_1.save("resize_1.png");

input_image_2=input_image_2.resize(max_width,max_height);
input_image_2.save("resize_2.png");

CImg <double> stiched_image=stcichImages(input_image_1,input_image_2);
stiched_image.save("stiched_image.png");

vector<SiftDescriptor> image_1_descriptors=Sift::compute_sift(input_image_1.get_RGBtoHSI().get_channel(2));
vector<SiftDescriptor> image_2_descriptors=Sift::compute_sift(input_image_2.get_RGBtoHSI().get_channel(2));

std::vector<data_tuple> matching_vectors;

for(int i=0;i<image_1_descriptors.size();i++){
    int point_1=-1,point_2=-1;
    double sd_1=DBL_MAX;
    double sd_2=DBL_MAX;
    double curr_distance;
    for(int j=0;j<image_2_descriptors.size();j++){
        curr_distance=getEuclideanDistance(image_1_descriptors[i].descriptor,image_2_descriptors[j].descriptor);
        if(curr_distance<sd_1){
            point_2=point_1;
            point_1=j;
            sd_2=sd_1;
            sd_1=curr_distance;
        }
        else{
            if(curr_distance<sd_2){
                point_2=j;
                sd_2=curr_distance;
            }
        }
    }
    double d_ratio=sd_1/sd_2;
    if((point_2!=-1)&&(d_ratio<=threshold)){
        matching_vectors.push_back({point_1,point_2,sd_1});
    }
}
for(int i=0;i<image_2_descriptors.size();i++){
    image_2_descriptors[i].col=image_2_descriptors[i].col+max_width;
}
// std::sort(matching_vectors.begin(),matching_vectors.end(),comparator_function);
CImg<double> sift_descriptor;
sift_descriptor=draw_descriptor_image(stiched_image,image_1_descriptors,"image_1_sift.png");
sift_descriptor=draw_descriptor_image(sift_descriptor,image_2_descriptors,"image_sift.png");
std::vector<line>lineVector;
for(int i=0;i<matching_vectors.size();i++){
    int row_1=image_1_descriptors[matching_vectors[i].feature_point_1].col;
    int col_1=image_1_descriptors[matching_vectors[i].feature_point_1].row;
    int row_2=image_2_descriptors[matching_vectors[i].feature_point_2].col;
    int col_2=image_2_descriptors[matching_vectors[i].feature_point_2].row;
    lineVector.push_back({row_1,col_1,row_2,col_2});
    }
CImg<double> sift_match_image=drawLines(sift_descriptor,lineVector);
sift_match_image.save("sift_match.png");
}

int main(int argc, char **argv)
{
  try {
        string part = argv[1];
    CImg<double> input_image("images/part1/lincoln.png");
    CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> input_descriptors = Sift::compute_sift(input_gray);
    draw_descriptor_image(input_image, input_descriptors, "input_image.png");
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

	// inverseMatrix(projectiveTransform,inverse);

    CImg<double> output_image(1024,1024);
    
    cimg_forXY(input_image,i,j)
    {
        float x = (inverse[0][0] * i + inverse[0][1] *j + inverse[0][2]*1);
        float y = (inverse[1][0] * i + inverse[1][1] *j + inverse[1][2]*1);
        float w = (inverse[2][0] * i + inverse[2][1] *j + inverse[2][2]*1);
        if((x/w > 0) && (x/w<1024) && (y/w>0) && (y/w<1024))    
            output_image(i,j) = input_image(int(x/w),int(y/w));
    }
    output_image.save("output.png");
    getTransformationMatrix(318,256,534,372,316,670,73,473,141,131,480,159,493,630,64,601);

    CImg<double> billboard_image("images/part1/billboard1.jpg");
    input_gray = billboard_image.get_RGBtoHSI().get_channel(2);

    input_descriptors = Sift::compute_sift(input_gray);
    draw_descriptor_image(billboard_image, input_descriptors, "sift.png");
/*  for(int i=0;i<input_gray.height();i++)
        for(int j=0;j<input_gray.width();j++)
        {
            if(input_gray(i,j) == 255)
            {
                cout<<i<<' '<<j<<endl;
                break;
            }   
        }
    
*/
    
    
    CImg<double> kernel(5,5);

    kernel(0,0) = 1/256 ; kernel(0,1) = 4/256 ; kernel(0,2) = 6/256 ;kernel(0,3) = 4/256 ;kernel(0,4) = 1/256 ;
    kernel(1,0) = 4/256 ; kernel(1,1) = 16/256 ; kernel(1,2) = 24/256 ;kernel(1,3) = 16/256 ;kernel(1,4) = 4/256 ;
    kernel(2,0) = 6/256 ; kernel(2,1) = 24/256 ; kernel(2,2) = 36/256 ;kernel(2,3) = 24/256 ;kernel(2,4) = 6/256 ;
    kernel(3,0) = 4/256 ; kernel(3,1) = 16/256 ; kernel(3,2) = 24/256 ;kernel(3,3) = 16/256 ;kernel(3,4) = 4/256 ;
    kernel(4,0) = 1/256 ; kernel(4,1) = 4/256 ; kernel(4,2) = 6/256 ;kernel(4,3) = 4/256 ;kernel(4,4) = 1/256 ;

    CImg<double> conv_image("images/part1/lincoln.png");
    CImg<double> conv_gray = conv_image.get_RGBtoHSI().get_channel(2);  
//  CImg<double> conv = input_gray.get_convolve(kernel,0);
        
    conv_gray.save("conv.png");
    
       
	if(part == "part1"){
	}
	else if(part == "part2"){
 
        }
	else if(part == "part3"){
	string file_1=argv[2];
	string file_2=argv[3];
	CImg<double> input_image_1(argv[2]);
	CImg<double> input_image_2(argv[3]);
	getMatchingWithBruteForce(input_image_1,input_image_2);
	
	}
	else if(part == "part4"){

	}
//    
//    
//feel free to add more conditions for other parts (e.g. more specific)
//parts, for debugging, etc.
  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}

