// kate: replace-tabs off; indent-width 4; indent-mode normal
// vim: ts=4:sw=4:noexpandtab
/*

Copyright (c) 2010--2012,
François Pomerleau and Stephane Magnenat, ASL, ETHZ, Switzerland
You can contact the authors at <f dot pomerleau at gmail dot com> and
<stephane at magnenat dot net>

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL ETH-ASL BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "pointmatcher/PointMatcher.h"
#include "pointmatcher/PointMatcherPrivate.h"
#include <cassert>
#include <iostream>
#include "boost/filesystem.hpp"

#include "pointmatcher/IO.h"
#include "pointmatcher/IOFunctions.h"
#include "pointmatcher/InspectorsImpl.h"
#include "Eigen/Core"

// For logging
#include "pointmatcher/PointMatcherPrivate.h"
#include <typeinfo>


using namespace std;
using namespace Eigen;

struct _rows {
    int cols_count;
    float *cols;
};

struct _unit {
    int rows_count;
    struct _rows *rows;
};
float ans[16];
float* test_icp(struct _unit *param1,int s_r, int t_r)
{
    int i,j;
	typedef PointMatcher<float> PM;
    typedef typename Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrix;
	typedef PM::DataPoints DP;

    DP loadedPoints1;
    DP loadedPoints2;
    Matrix sourcefeatures = Matrix(3, s_r);
    Matrix targetfeatures = Matrix(3, t_r);
    cout<<param1->rows_count<<endl;
    for(i=0;i<s_r+t_r;i++){
        for(j=0;j<3;j++){
            //printf("%d,%d = %f \n",i,j,param1->rows[i].cols[j]);
            if(i<s_r)
            	sourcefeatures(j,i) = param1->rows[i].cols[j];
            else
            	targetfeatures(j,i-s_r) = param1->rows[i].cols[j];
		}
    }
    loadedPoints1.addFeature("x", sourcefeatures.row(0));
    loadedPoints1.addFeature("y", sourcefeatures.row(1));
    loadedPoints1.addFeature("z", sourcefeatures.row(2));

    loadedPoints2.addFeature("x", targetfeatures.row(0));
    loadedPoints2.addFeature("y", targetfeatures.row(1));
    loadedPoints2.addFeature("z", targetfeatures.row(2));
    PM::ICP icp;
	icp.setDefault();
    PM::TransformationParameters T2 = icp(loadedPoints1, loadedPoints2);
    //cout << "Final transformation:" << endl << T2 << endl;
    float* T = &T2(0,0);
    std::vector<float> x2y2;
    x2y2.assign(T, T + T2.size());

    for(int k = 0;k<16;k++){
	    ans[k]= x2y2[k];
        //cout<<ans[k]<<endl;
	}
    return ans;
}

void validateArgs(int argc, char *argv[], bool& isCSV);

/**
  * Code example for ICP taking 2 points clouds (2D or 3D) relatively close 
  * and computing the transformation between them.
  */

int main(int argc, char *argv[])
{
	bool isCSV = true;
	validateArgs(argc, argv, isCSV);
	typedef PointMatcher<float> PM;
    typedef typename Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrix;
	typedef PM::DataPoints DP;
	
	// Load point clouds
	const DP ref(DP::load(argv[1]));
	const DP data(DP::load(argv[2]));
    DP loadedPoints1;
    DP loadedPoints2;
    Matrix features = Matrix(1, 10000);
    for(int i = 0 ; i<10000;i++){
		features(0,i) = 0.5f+i/10000.0;

	}
    cout<<"features :"<<endl;
    cout<<features<<endl;
    loadedPoints1.addFeature("x", features.row(0));
    loadedPoints1.addFeature("y", features.row(1));
    loadedPoints1.addFeature("z", features.row(2));
    //loadedPoints1.addFeature("pad",features.row(3));

    loadedPoints2.addFeature("y", -features.row(0));
    loadedPoints2.addFeature("z", -features.row(1));
    loadedPoints2.addFeature("x", -features.row(2));
    //loadedPoints1.addFeature("pad", features.row(3));



	// Create the default ICP algorithm
	PM::ICP icp;
	
	// See the implementation of setDefault() to create a custom ICP algorithm
	icp.setDefault();
    cout<<"loadedPoints1 :"<<endl;
    cout<<loadedPoints1.features.rows()<<endl;
    cout<<loadedPoints1.features.cols()<<endl;
    cout<<"loadedPoints2 :"<<endl;
    cout<<loadedPoints2.features.rows()<<endl;
    cout<<loadedPoints2.features.cols()<<endl;

    cout<<typeid(features.row(1)).name()<<endl;
    cout<<typeid(ref.features.row(1)).name()<<endl;

	// Compute the transformation to express data in ref
	PM::TransformationParameters T2 = icp(loadedPoints1, loadedPoints2);

	// Transform data to express it in ref
	DP data_out(data);
	icp.transformations.apply(data_out, T2);
	
	// Safe files to see the results
	ref.save("test_ref.vtk");
	data.save("test_data_in.vtk");
	data_out.save("test_data_out.vtk");
	cout << "Final transformation:" << endl << T2 << endl;

	return 0;
}

void validateArgs(int argc, char *argv[], bool& isCSV )
{
	if (argc != 3)
	{
		cerr << "Wrong number of arguments, usage " << argv[0] << " reference.csv reading.csv" << endl;
		cerr << "Will create 3 vtk files for inspection: ./test_ref.vtk, ./test_data_in.vtk and ./test_data_out.vtk" << endl;
		cerr << endl << "2D Example:" << endl;
		cerr << "  " << argv[0] << " ../../examples/data/2D_twoBoxes.csv ../../examples/data/2D_oneBox.csv" << endl;
		cerr << endl << "3D Example:" << endl;
		cerr << "  " << argv[0] << " ../../examples/data/car_cloud400.csv ../../examples/data/car_cloud401.csv" << endl;
		exit(1);
	}
}
extern "C" {
        float* Test_icp(struct _unit *param1,int r1, int r2){return test_icp(param1,r1,r2);}
}
