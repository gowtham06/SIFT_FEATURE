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
#include <cstdlib>
#include <ctime>

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

struct Model{
CImg <double> transform_matrix;
int maxInliners;
std::vector<line> match_lines;
};

CImg <double> draw_descriptor_image(CImg<double> image, const vector<SiftDescriptor> descriptors, const char *filename)
{
    // cout<<"In draw_descriptor_image()m_1 function"<<endl;
  for(unsigned int i=0; i < descriptors.size(); i++)
    {
      int tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
      double color_point[] = {255.0, 255.0, 0};
      for(int x=-2; x<3; x++)
	for(int y=-2; y<3; y++)
	  if(x==0 || y==0)
	    for(int c=0; c<3; c++){
	      //Find if coordinates are in workspace to draw crosshair
	      tx1 = (descriptors.at(i).col + y - 1);
	      ty1 = (descriptors.at(i).row + x - 1);
	      if (tx1 >= 0 && tx1 < image.width() && ty1 >= 0 && ty1 < image.height())
	      {	
		image( tx1, ty1, 0, c) = color_point[c];		
	      }				
	    }
    }
  image.get_normalize(0,255).save(filename);
  // cout<<"End of draw_descriptor_image()m_1 function"<<endl;
  return image;
}

CImg <double> draw_descriptor_image(CImg<double> image, const vector<line> match_lines, const char *filename)
{
    // cout<<"In draw_descriptor_image()m_2 function"<<endl;
    int i=0;
  for(i=0; i < int(match_lines.size()); i++)
    {
      int tx1 = 0, ty1 = 0, tx2 = 0, ty2 = 0;
      double color_point[] = {255.0, 255.0, 0};
      for(int x=-2; x<3; x++)
    for(int y=-2; y<3; y++)
      if(x==0 || y==0)
        for(int c=0; c<3; c++){
          //Find if coordinates are in workspace to draw crosshair
          tx1 = (match_lines.at(i).x_1 + y - 1);
          ty1 = (match_lines.at(i).y_1 + x - 1);

          tx2=(match_lines.at(i).x_2+y-1);
          ty2 = (match_lines.at(i).y_2 + x - 1);

          if (tx1 >= 0 && tx1 < image.width() && ty1 >= 0 && ty1 < image.height())
          { 
        image( tx1, ty1, 0, c) = color_point[c];        
          }             
         if (tx2 >= 0 && tx2 < image.width() && ty2 >= 0 && ty2 < image.height())
          { 
        image( tx2, ty2, 0, c) = color_point[c];        
          }             
           
        }
    }
  image.get_normalize(0,255).save(filename);
  // cout<<"End of draw_descriptor_image()m_2 function"<<endl;
  return image;
}



void getWarpedImage(CImg<double> input_image,double projectiveTransform[][3])
{
    double inverse[3][3];
    inverseMatrix(projectiveTransform,inverse);
    CImg<double> output_image(input_image.width(),input_image.height(),1,3,255);
    
    cimg_forXY(input_image,i,j)
    {
        float x = (inverse[0][0] * i + inverse[0][1] *j + inverse[0][2]*1);
        float y = (inverse[1][0] * i + inverse[1][1] *j + inverse[1][2]*1);
        float w = (inverse[2][0] * i + inverse[2][1] *j + inverse[2][2]*1);
        if((x/w > 0) && (x/w<input_image.width()) && (y/w>0) && (y/w<input_image.height())) 
        {
            output_image(i,j,0,0) = input_image(int(x/w),int(y/w),0,0);
            output_image(i,j,0,1) = input_image(int(x/w),int(y/w),0,1);
            output_image(i,j,0,2) = input_image(int(x/w),int(y/w),0,2);
        }
    }

    output_image.save("outputPart14_1.png");

}

CImg<double> getTransformationMatrix(int x1,int y1,int x2,int y2,int x3,int y3,int x4,int y4,int X1,int Y1,int X2,int Y2,int X3,int Y3,int X4,int Y4)
{
    CImg<long double> A(8,8);
    A(0,0) = double(x1); A(0,1) = double(y1); A(0,2) = double(1); A(0,3) = double(0); A(0,4) = double(0); A(0,5) = double(0); A(0,6) = double(-x1*X1); A(0,7) = double(-y1*X1);

    A(1,0) = double(0); A(1,1) = double(0); A(1,2) = double(0); A(1,3) = double(x1); A(1,4) = double(y1); A(1,5) = double(1); A(1,6) = double(-x1*Y1); A(1,7) = double(-y1*Y1);

    A(2,0) = double(x2); A(2,1) = double(y2); A(2,2) = double(1); A(2,3) = double(0); A(2,4) = double(0); A(2,5) = double(0); A(2,6) = double(-x2*X2); A(2,7) = double(-y2*X2);

    A(3,0) = double(0); A(3,1) = double(0); A(3,2) = double(0); A(3,3) = double(x2); A(3,4) = double(y2); A(3,5) = double(1); A(3,6) = double(-x2*Y2); A(3,7) = double(-y2*Y2);

    A(4,0) = double(x3); A(4,1) = double(y3); A(4,2) =double(1); A(4,3) = double(0); A(4,4) = double(0); A(4,5) = double(0); A(4,6) =double(-x3*X3); A(4,7) = double(-y3*X3);

    A(5,0) = double(0); A(5,1) = double(0); A(5,2) = double(0); A(5,3) = double(x3); A(5,4) = double(y3); A(5,5) = double(1); A(5,6) = double(-x3*Y3); A(5,7) = double(-y3*Y3);

    A(6,0) = double(x4); A(6,1) = double(y4); A(6,2) = double(1); A(6,3) = double(0); A(6,4) = double(0); A(6,5) = double(0); A(6,6) = double(-x4*X4); A(6,7) = double(-y4*X4);

    A(7,0) = double(0); A(7,1) = double(0); A(7,2) = double(0); A(7,3) = double(x4); A(7,4) = double(y4); A(7,5) = double(1); A(7,6) = double(-x4*Y4); A(7,7) = double(-y4*Y4);


    CImg<long double> B(1,8);
    //  cout<<B.width()<<endl;
    CImg<double> X(3,3);
    B(0,0) = X1;
    B(0,1) = Y1;    
    B(0,2) = X2;
    B(0,3) = Y2;
    B(0,4) = X3;
    B(0,5) = Y3;
    B(0,6) = X4;
    B(0,7) = Y4;

    X = B.solve(A.transpose());
    return X;
    /*CImg<double> result(3,3);
    cimg_forXY(result,i,j)
    {
        result(i,j) = X(i,j);
    }
    result(2,2) = 1;    
    return result;*/
}


CImg<double> drawLines(CImg<double> image,std::vector<line> lineVector){
    // cout<<"In drawLines () function"<<endl;
    const unsigned char color[] = { 255,128,64 };
    // input_image.draw_line(40,40,80,70,color);
    for(int i=0;i<int(lineVector.size());i++){
        image.draw_line(lineVector.at(i).x_1,lineVector.at(i).y_1,lineVector.at(i).x_2,lineVector.at(i).y_2,color);
    }
    // cout<<"End of drawLines() function"<<endl;
    return image;

}

bool comparator_function(data_tuple obj_1,data_tuple obj_2){
    // cout<<"Inside comparator_function()"<<endl;
    // cout<<"End of comparator_function() function"<<endl;
    
    return (obj_1.short_distance<obj_2.short_distance);
}


double getEuclideanDistance(std::vector<float> datapoint_1,std::vector<float> datapoint_2){
double distance=0.0;
// cout<<"getEuclideanDistanc() function"<<endl;
// cout<<datapoint_1.size()<<","<<datapoint_2.size()<<endl;
for(int i=0;i<int(datapoint_1.size());i++){
    distance+=(datapoint_1.at(i)-datapoint_2.at(i))*(datapoint_1.at(i)-datapoint_2.at(i));
    // cout<<"("<<datapoint_1[i]<<","<<datapoint_2[i]<<")"<<endl;
}
// cout<<"D:"<<sqrtf(distance)<<endl;
return sqrt(distance);

}

CImg<double> stcichImages(CImg<double> image_1,CImg<double> image_2){
int max_width=image_1.width();
int max_height=image_1.height();

CImg<double> stiched_image(2*max_width,max_height,1,3);
// cout<<"In stcichImages() function"<<endl;
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

void draw_sift_matches(CImg<double> input_image_1,CImg<double> input_image_2,std::vector<line> lineVector,int max_width,const char* output_file_name){
    // cout<<"In draw_sift_matches()_m2 function"<<endl;
    vector<SiftDescriptor> image_1_descriptors=Sift::compute_sift(input_image_1.get_RGBtoHSI().get_channel(2));
    vector<SiftDescriptor> image_2_descriptors=Sift::compute_sift(input_image_2.get_RGBtoHSI().get_channel(2));
    int i=0;
    std::vector<SiftDescriptor> v;
    for(i=0;i<int(image_2_descriptors.size());i++){
    image_2_descriptors[i].col=image_2_descriptors[i].col+max_width;
    }
    CImg <double> stiched_image=stcichImages(input_image_1,input_image_2);
    CImg<double> sift_descriptor;
    for(i=0;i<int(lineVector.size());i++){
        lineVector[i].x_2=lineVector[i].x_2+max_width;
    }
    sift_descriptor=draw_descriptor_image(stiched_image,lineVector,output_file_name);
    CImg<double> sift_match_image=drawLines(sift_descriptor,lineVector);
    sift_match_image.save(output_file_name);
    // cout<<"End of draw_sift_matches()_m2 function"<<endl;
}

std::vector<line> get_sift_matching(CImg<double> input_image_1,CImg<double> input_image_2,double threshold=0.85) {
    // cout<<"get_sift_matching() function"<<endl;
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
input_image_2=input_image_2.resize(max_width,max_height);
vector<SiftDescriptor> image_1_descriptors=Sift::compute_sift(input_image_1.get_RGBtoHSI().get_channel(2));
vector<SiftDescriptor> image_2_descriptors=Sift::compute_sift(input_image_2.get_RGBtoHSI().get_channel(2));
std::vector<data_tuple> matching_vectors;
int i=0,j=0;
for(i=0;i<int(image_1_descriptors.size());i++){
    int point_1=-1,point_2=-1;
    double sd_1=DBL_MAX;
    double sd_2=DBL_MAX;
    double curr_distance;
    for( j=0;j<int(image_2_descriptors.size());j++){
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
        matching_vectors.push_back({i,point_1,sd_1});
    }

}
    std::vector<line>lineVector;
    
for(i=0;i<int(matching_vectors.size());i++){
    int row_1=image_1_descriptors[matching_vectors[i].feature_point_1].col;
    int col_1=image_1_descriptors[matching_vectors[i].feature_point_1].row;
    int row_2=image_2_descriptors[matching_vectors[i].feature_point_2].col;
    int col_2=image_2_descriptors[matching_vectors[i].feature_point_2].row;
    lineVector.push_back({row_1,col_1,row_2,col_2});
    }
    draw_sift_matches(input_image_1,input_image_2,lineVector,max_width,"sift_match.png");
    // cout<<"End of get_sift_matching() function"<<endl;
   return lineVector;

}

Model generateRansacModel(CImg<double> image1, CImg<double> image2)
{
    
    // cout<<"In generateRansacModel() function"<<endl;
    struct Model ransac_model;
    vector<line> siftCorrespondence = get_sift_matching(image1, image2);
    int iterationCount = 30000, i, j, inliner, maxInline=-2147483648;
    double x_1, y_1, x_2, y_2, X1, Y1, h;
    line l1, l2, l3, l4;
    CImg<double> tranformMatrix, bestMatrix;
    std::vector<line> sift_match_lines;
    for(i=0;i<iterationCount;i++)
    {
        std::srand(unsigned(std::time(0)));
        std::random_shuffle(siftCorrespondence.begin(), siftCorrespondence.end());
        l1 = siftCorrespondence.at(0);
        l2 = siftCorrespondence.at(1);
        l3 = siftCorrespondence.at(2);
        l4 = siftCorrespondence.at(3);
        inliner = 0;
        std::vector<line> temp;
         tranformMatrix = getTransformationMatrix(l1.x_1,l1.y_1,l2.x_1,l2.y_1,l3.x_1,l3.y_1,l4.x_1,l4.y_1,l1.x_2,l1.y_2,l2.x_2,l2.y_2,l3.x_2,l3.y_2,l4.x_2,l4.y_2);
        for(j=4;j<int(siftCorrespondence.size());j++)
        {
            x_1 = siftCorrespondence.at(j).x_1;
            y_1 = siftCorrespondence.at(j).y_1;
            x_2 = siftCorrespondence.at(j).x_2;
            y_2 = siftCorrespondence.at(j).y_2;
            
            X1 = tranformMatrix(0,0)*x_1 + tranformMatrix(0,1)*y_1 + tranformMatrix(0,2)*1;
            Y1 = tranformMatrix(0,3)*x_1 + tranformMatrix(0,4)*y_1 + tranformMatrix(0,5)*1;
            h =  tranformMatrix(0,6)*x_1 + tranformMatrix(0,7)*y_1 + 1;
            
            X1/=h;
            Y1/=h;
            
            if(abs(X1-x_2) <= 10 and abs(Y1-y_2) <= 10){
                temp.push_back({x_1,y_1,x_2,y_2});
                inliner++;
            }
        }
        

        if(inliner > maxInline)
        {
            maxInline = inliner;
            bestMatrix = tranformMatrix;
            sift_match_lines=temp;
        }

    }
    ransac_model.transform_matrix=bestMatrix;
    ransac_model.maxInliners=maxInline;
    ransac_model.match_lines=sift_match_lines;
    // cout<<"End of generateRansacModel() function"<<endl;

    return ransac_model;


}

void  sift_match_pruning(CImg<double> input_image_1,CImg <double> input_image_2){
    // cout<<"sift_match_pruning() function"<<endl;
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
input_image_2=input_image_2.resize(max_width,max_height);
Model ransac_model=generateRansacModel(input_image_1,input_image_2);
// draw_sift_matches(input_image_1,input_image_2,ransac_model.match_lines,max_width,"sift_ransac_match.jpg");
 draw_sift_matches(input_image_1,input_image_2,ransac_model.match_lines,max_width,"sift_ransac_match.png");
// cout<<"End of sift_match_pruning() function"<<endl;


}


int main(int argc, char **argv)
{
  try {
        string part = argv[1];
    CImg<double> input_image("images/part1/lincoln.png");
    CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
    vector<SiftDescriptor> input_descriptors = Sift::compute_sift(input_gray);
    draw_descriptor_image(input_image, input_descriptors, "input_image.png");
    
    
       
	if(part == "part1"){
			CImg<double> input_image("images/part1/lincoln.png");
			CImg<double> input_gray = input_image.get_RGBtoHSI().get_channel(2);
			CImg<double> projectiveTransform1(3,3);

			projectiveTransform1(0,0) = 0.907;
			projectiveTransform1(1,0) = 0.258;
			projectiveTransform1(2,0) = -182;
			projectiveTransform1(0,1) = -0.153;
			projectiveTransform1(1,1) = 1.44;
			projectiveTransform1(2,1) = 58;
			projectiveTransform1(0,2) = -0.000306;
			projectiveTransform1(1,2) = 0.000731;
			projectiveTransform1(2,2) = 1;

			CImg<double> inverse1 = projectiveTransform1.invert();
			CImg<double> output_image(1024,1024,1,3,0);


			cimg_forXY(input_image,i,j)
			{
				float x = (inverse1(0,0) * i + inverse1(1,0) *j + inverse1(2,0)*1);
				float y = (inverse1(0,1) * i + inverse1(1,1) *j + inverse1(2,1)*1);
				float w = (inverse1(0,2) * i + inverse1(1,2) *j + inverse1(2,2)*1);
				if((x/w > 0) && (x/w<input_image.width()) && (y/w>0) && (y/w<input_image.height()))	
				{
					output_image(i,j,0,0) = input_image(int(x/w),int(y/w),0,0);
					output_image(i,j,0,1) = input_image(int(x/w),int(y/w),0,1);
					output_image(i,j,0,2) = input_image(int(x/w),int(y/w),0,2);	
				}	
			}
//			getWarpedImage(input_image,projectiveTransform1);
			output_image.save("lincoln_warped.png");
			CImg<double> transform(3,3);
			transform = getTransformationMatrix(141,131,480,159,493,630,64,601,318,256,534,372,316,670,73,473);

			
			cimg_forXY(transform,i,j)
				cout<<transform(i,j)<<endl;

			CImg<double> book("images/part1/book2.jpg");
			CImg<double> inverse2 = transform.invert();
			CImg<double> output_image1(book.width(),book.height(),1,3,0);


			cimg_forXY(book,i,j)
			{
				float x = (inverse2(0,0) * i + inverse2(1,0) *j + inverse2(2,0)*1);
				float y = (inverse2(0,1) * i + inverse2(1,1) *j + inverse2(2,1)*1);
				float w = (inverse2(0,2) * i + inverse2(1,2) *j + inverse2(2,2)*1);
				if((x/w > 0) && (x/w<book.width()) && (y/w>0) && (y/w<book.height()))	
				{
					output_image1(i,j) = book(int(x/w),int(y/w));
					output_image1(i,j,0,0) = book(int(x/w),int(y/w),0,0);
					output_image1(i,j,0,1) = book(int(x/w),int(y/w),0,1);
					output_image1(i,j,0,2) = book(int(x/w),int(y/w),0,2);	
				}	
			}
			output_image1.save("book_warped.png");

			CImg<double> transform1(3,3);
			transform1  = getTransformationMatrix(101,60,530,60,101,204,530,204,0,0,0,1024,1024,0,1024,1024);
		
			getWarpedImage(input_image,transform1);	
			CImg<double> transform2(3,3);
			transform2  = getTransformationMatrix(0,0,0,1024,1024,0,1024,1024,176,54,1107,260,147,625,1126,702);

			CImg<double> transform3(3,3);
			transform3  = getTransformationMatrix(141,131,480,159,493,630,64,601,318,256,534,372,316,670,73,473);



	}
	else if(part == "part2"){
		
			CImg<double> kernel(5,5);

			kernel(0,0) = 1.00/256.00 ; kernel(0,1) = 4.00/256.00; kernel(0,2) = 6.00/256.00 ;kernel(0,3) = 4.00/256.00 ;kernel(0,4) = 1.00/256.00 ;
			kernel(1,0) = 4.00/256.00 ; kernel(1,1) = 16.00/256.00 ; kernel(1,2) = 24.00/256.00 ;kernel(1,3) = 16.00/256.00 ;kernel(1,4) = 4.00/256.00;
			kernel(2,0) = 6.00/256.00 ; kernel(2,1) = 24.00/256.00 ; kernel(2,2) = 36.00/256.00 ;kernel(2,3) = 24.00/256.00 ;kernel(2,4) = 6.00/256.00 ;
			kernel(3,0) = 4.00/256.00 ; kernel(3,1) = 16.00/256.00 ; kernel(3,2) = 24.00/256.00 ;kernel(3,3) = 16.00/256.00 ;kernel(3,4) = 4.00/256.00 ;
			kernel(4,0) = 1.00/256.00 ; kernel(4,1) = 4.00/256.00 ; kernel(4,2) = 6.00/256.00 ;kernel(4,3) = 4.00/256.00 ;kernel(4,4) = 1.00/256.00 ;



			CImg<double> apple_image(argv[2]);
			CImg<double> orange_image(argv[3]);
			CImg<double> mask(argv[4]);	

			CImg<double> G0_apple = apple_image;
			CImg<double> G0_orange = orange_image;
			CImg<double> G0_mask(307,307,1,3);
			cimg_forXYC(G0_mask,i,j,k)
			{
				G0_mask(i,j,k) = mask(i,j);
			}
			G0_mask.normalize(0,1);	 

			CImg<double> G1_apple;
			CImg<double> G1_orange;
			CImg<double> G1_mask;
			CImg<double> G1_apple_temp = apple_image.get_convolve(kernel,0,false);
			CImg<double> G1_orange_temp = orange_image.get_convolve(kernel,0,false);
			CImg<double> G1_mask_temp = G0_mask.get_convolve(kernel,0,false);
			G1_apple = G1_apple_temp.get_resize_halfXY();
			G1_orange = G1_orange_temp.get_resize_halfXY();
			G1_mask = G1_mask_temp.get_resize_halfXY();

			CImg<double> G2_apple;
			CImg<double> G2_orange;
			CImg<double> G2_mask;
			CImg<double> G2_apple_temp = G1_apple.get_convolve(kernel,0,false);
			CImg<double> G2_orange_temp = G1_orange.get_convolve(kernel,0,false);
			CImg<double> G2_mask_temp = G1_mask.get_convolve(kernel,0,false);
			G2_apple = G2_apple_temp.get_resize_halfXY();
			G2_orange = G2_orange_temp.get_resize_halfXY();
			G2_mask = G2_mask_temp.get_resize_halfXY();
			G2_apple.save("G2_apple.jpg");

			CImg<double> G3_apple;
			CImg<double> G3_orange;
			CImg<double> G3_mask;
			CImg<double> G3_apple_temp = G2_apple.get_convolve(kernel,0,false);
			CImg<double> G3_orange_temp = G2_orange.get_convolve(kernel,0,false);
			CImg<double> G3_mask_temp = G2_mask.get_convolve(kernel,0,false);
			G3_apple = G3_apple_temp.get_resize_halfXY();
			G3_orange = G3_orange_temp.get_resize_halfXY();
			G3_mask = G3_mask_temp.get_resize_halfXY();
			G3_apple.save("G3_apple.jpg");	

			CImg<double> G4_apple;
			CImg<double> G4_orange;
			CImg<double> G4_mask;
			CImg<double> G4_apple_temp = G3_apple.get_convolve(kernel,0,false);
			CImg<double> G4_orange_temp = G3_orange.get_convolve(kernel,0,false);
			CImg<double> G4_mask_temp = G3_mask.get_convolve(kernel,0,false);
			G4_apple = G4_apple_temp.get_resize_halfXY();
			G4_orange = G4_orange_temp.get_resize_halfXY();
			G4_mask = G4_mask_temp.get_resize_halfXY();
			G4_apple.save("G4_apple.jpg");

			CImg<double> G5_apple;
			CImg<double> G5_orange;
			CImg<double> G5_mask;
			CImg<double> G5_apple_temp = G4_apple.get_convolve(kernel,0,false);
			CImg<double> G5_orange_temp = G4_orange.get_convolve(kernel,0,false);
			CImg<double> G5_mask_temp = G4_mask.get_convolve(kernel,0,false);
			G5_apple = G5_apple_temp.get_resize_halfXY();
			G5_orange = G5_orange_temp.get_resize_halfXY();
			G5_mask = G5_mask_temp.get_resize_halfXY();
			G5_apple.save("G5_apple.jpg");

			CImg<double> L0_apple = G5_apple;
			CImg<double> L0_orange = G5_orange;

			CImg<double> L1_apple;
			CImg<double> L1_orange;
			CImg<double> G5_apple_upscale = G5_apple.resize_doubleXY();
			CImg<double> G5_orange_upscale = G5_orange.resize_doubleXY();	
			L1_apple = G4_apple - G5_apple_upscale;
			L1_orange = G4_orange - G5_orange_upscale;
			cout<<G5_apple_upscale.width()<<' ' <<G5_apple_upscale.height()<<endl;
			L1_apple.normalize(0,255);
			L1_orange.normalize(0,255);
			L1_apple.save("L1_apple.jpg");


			CImg<double> L2_apple;
			CImg<double> L2_orange;
			CImg<double> G4_apple_upscale = G4_apple.resize_doubleXY();
			CImg<double> G4_orange_upscale = G4_orange.resize_doubleXY();
			L2_apple = G3_apple - G4_apple_upscale;
			L2_orange = G3_orange - G4_orange_upscale;
			L2_apple.normalize(0,255);
			L2_orange.normalize(0,255);
			L2_apple.save("L2_apple.jpg");
		
			
			CImg<double> L3_apple;
			CImg<double> L3_orange;
			CImg<double> G3_apple_upscale = G3_apple.resize_doubleXY();
			CImg<double> G3_orange_upscale = G3_apple.resize_doubleXY();
			L3_apple = G2_apple - G3_apple_upscale;
			L3_orange = G2_orange - G3_orange_upscale;
			L3_apple.normalize(0,255);
			L3_orange.normalize(0,255);
			L3_apple.save("L3_apple.jpg");

			CImg<double> L4_apple;
			CImg<double> L4_orange;
			CImg<double> G2_apple_upscale = G2_apple.resize_doubleXY();
			CImg<double> G2_orange_upscale = G2_orange.resize_doubleXY();
			L4_apple = G1_apple - G2_apple_upscale;
			L4_orange = G1_orange - G2_orange_upscale;
			L4_apple.normalize(0,255);
			L4_orange.normalize(0,255);
			L4_apple.save("L4_apple.jpg");

			CImg<double> L5_apple;
			CImg<double> L5_orange;
			CImg<double> G1_apple_upscale = G1_apple.resize_doubleXY();
			CImg<double> G1_orange_upscale = G1_orange.resize_doubleXY();
			L5_apple = G0_apple - G1_apple_upscale;
			L5_orange = G0_orange - G1_orange_upscale;
			L5_apple.normalize(0,255);
			L5_orange.normalize(0,255);
			L5_apple.save("L5_apple.jpg");

			CImg<double> L0_blended(10,10);
			CImg<double> L1_blended(20,20);
			CImg<double> L2_blended(40,40);
			CImg<double> L3_blended(80,80);
			CImg<double> L4_blended(160,160);
			CImg<double> L5_blended;
			CImg<double> G0_mask_ones(307,307,3);
			CImg<double> G1_mask_ones(G1_mask.width(),G1_mask.height(),3);
			G0_mask_ones.fill(1);




			L5_blended = G0_mask.mul(L5_apple) + (G0_mask_ones - G0_mask).mul(L5_orange);	
			L5_blended.normalize(0,255);
			L5_blended.save("L5_blended.jpg");

			L4_blended = G1_mask.mul(L4_apple) + (G0_mask_ones.resize(G1_mask.width(),G1_mask.height()) - G1_mask).mul(L4_orange);
			L4_blended.normalize(0,255);
			L4_blended.save("L4_blended.jpg");

			L3_blended = G2_mask.mul(L3_apple) + (G0_mask_ones.resize(G2_mask.width(),G2_mask.height()) - G2_mask).mul(L3_orange);
			L3_blended.normalize(0,255);
			L3_blended.save("L3_blended.jpg");

			L2_blended = G3_mask.mul(L2_apple) + (G0_mask_ones.resize(G3_mask.width(),G3_mask.height()) - G3_mask).mul(L2_orange);
			L2_blended.normalize(0,255);
			L2_blended.save("L2_blended.jpg");

			L1_blended = G4_mask.mul(L1_apple) + (G0_mask_ones.resize(G4_mask.width(),G4_mask.height()) - G4_mask).mul(L1_orange);
			L1_blended.normalize(0,255);
			L1_blended.save("L1_blended.jpg");

			L0_blended = G5_mask.mul(L0_apple) + (G0_mask_ones.resize(G5_mask.width(),G5_mask.height()) - G5_mask).mul(L0_orange);
			L0_blended.normalize(0,255);
			L0_blended.save("L0_blended.jpg");


			CImg<double> G0;
			CImg<double> G1;
			CImg<double> G2;
			CImg<double> G3;
			CImg<double> G4;
			CImg<double> G5;

			G0 = L0_blended;
			G1 = L1_blended + G0.resize(L1_blended.width(),L1_blended.height());
			G2 = L2_blended + G1.resize(L2_blended.width(),L2_blended.height());
			G3 = L3_blended + G2.resize(L3_blended.width(),L3_blended.height());
			G4 = L4_blended + G3.resize(L4_blended.width(),L4_blended.height());
			G5 = L5_blended + G4.resize(L5_blended.width(),L5_blended.height());
			G5.normalize(0,255);	
			G5.save("Final_blended.jpg");



 
        }
	else if(part == "part3"){
	string file_1=argv[2];
	string file_2=argv[3];
	CImg<double> input_image_1(argv[2]);
	CImg<double> input_image_2(argv[3]);
	get_sift_matching(input_image_1,input_image_2);
    sift_match_pruning(input_image_1,input_image_2);
    // cout<<"part3 end"<<endl;
 }
	else if(part == "part4"){
        string file_1=argv[2];
            string file_2=argv[3];
            CImg<double> input_image_1(argv[2]);
            CImg<double> input_image_2(argv[3]);
        Model ransac_model=generateRansacModel(input_image_1,input_image_2);
        CImg<double> transformationMatrix=ransac_model.transform_matrix;
        double transformedMatrix[3][3];

        transformedMatrix[0][0] = transformationMatrix(0,0);
                transformedMatrix[0][1] = transformationMatrix(0,1);
                transformedMatrix[0][2] = transformationMatrix(0,2);
                transformedMatrix[1][0] = transformationMatrix(0,3);
                transformedMatrix[1][1] = transformationMatrix(0,4);
                transformedMatrix[1][2] = transformationMatrix(0,5);
                transformedMatrix[2][0] = transformationMatrix(0,6);
                transformedMatrix[2][1] = transformationMatrix(0,7);
                transformedMatrix[2][2] = 1;
        
        getWarpedImage(input_image_2, transformedMatrix);
    }

  }
  catch(const string &err) {
    cerr << "Error: " << err << endl;
  }
}

