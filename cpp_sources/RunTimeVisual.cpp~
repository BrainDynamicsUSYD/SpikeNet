#include <vector>
#include <opencv2/contrib/contrib.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include "RunTimeVisual.h"

using namespace cv;
using namespace std;

void run_time_visual(const char* win_name, vector<double> &v, const int row, const int col, const double min, const double max){

	
	Mat image,image_255,image_jet;

	// scale data
	vector<double> v_01;
	int range = max - min;
	for (int i = 0;i < row*col; i++){
		v_01.push_back((v[i]-min)/range);
	}

	// read scaled data into image
	size_t step = Mat::AUTO_STEP;
	image = Mat(row, col, CV_64F, v_01.data(), step); // CV_64F for double

	// convert double unit8
    	image.convertTo(image_255, CV_8UC1, 255.0 , 0); 

        // apply colormap
        applyColorMap(image_255, image_jet, COLORMAP_JET);

	// show image
   	namedWindow( win_name, CV_WINDOW_NORMAL);// Create a window for display. CV_WINDOW_AUTOSIZE 
	imshow( win_name, image_jet ); 
	waitKey(2); // two threads, waitKey(1) gives 1 ms for the plotting thread to finish the picture
	// For a 10,000 neuron population and visualisation at each step, 2 ms gives all the frames. 
}


