// B657 assignment 2 skeleton code
//
// Compile with: "make"
//
// See assignment handout for command line and project specifications.


//Link to the header file
#include <ctime>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <math.h>
#include "CImg.h"
#include "Sift.h"
#include "Matrix.h"
#include "Transform.h"

//Use the cimg namespace to access the functions easily
using namespace cimg_library;
using namespace std;

// holds the name of the image file and the matching score of sift descriptors
class Image
{
    public:

    string name;
    int count; 

    void setParameters(string name, int count)
    {
        this->name = name;
        this->count = count;
    }

    int getCount()
    {
        return count;
    }

    string getName()
    {
        return name;
    }

    bool operator < (const Image& i1) const
    {
        return count > i1.count;
    }
};

int sift_matching(CImg<double> input1, CImg<double> input2, string part, int MINIMUM_SIFT_DISTANCE  = 100)
{
    int matches=0, maxHeight = 0;
    const unsigned char color[] = {0, 255, 0};
    
    CImg<double> greyScale1, greyScale2;
    
    // converting the images to grayscale
    if(input1.spectrum() == 1)  
        greyScale1 = input1;
    else    
        greyScale1 = input1.get_RGBtoHSI().get_channel(2);

    if(input2.spectrum() == 1)
    	greyScale2 = input2;
    else
	    greyScale2 = input2.get_RGBtoHSI().get_channel(2); 
    
    vector<SiftDescriptor> desc1 = Sift::compute_sift(greyScale1);
    vector<SiftDescriptor> desc2 = Sift::compute_sift(greyScale2);
    
    
    // creating a merged image using the two images
    maxHeight = max(input1.height(),input2.height());

    CImg<double> output(input1.width() + input2.width(), maxHeight, 1, input1.spectrum(), 0);

	// first image on the left of the merged image
	if(part == "part1.1")
    {
        for(int i = 0; i < input1.width(); i++)
    	{
    		for(int j = 0; j < input1.height(); j++)
    		{
    	        for(int k = 0; k < input1.spectrum(); k++)
    			{
    			    output(i, j, 0, k) = input1(i, j, 0, k);                    
    			}
    		}
    	}

    	// second image on the right of the merged image
    	for(int i = input1.width(); i < output.width(); i++)
    	{
    		for(int j = 0; j < input2.height(); j++)
    		{
    	        for(int k = 0; k < input2.spectrum(); k++)
    			{
    			    output(i, j, 0, k) = input2(i - input1.width(), j, 0, k);                    
    			}
    		}
    	}
    }
    
    // matching both the descriptors and obtaining the number of matches
    for(int i = 0; i < desc1.size(); i++)
    {
        for(int j = 0; j < desc2.size(); j++)
        {
            double sum = 0;

            for(int k = 0; k < 128; k++)
            {
                sum += pow((desc1[i].descriptor[k]-desc2[j].descriptor[k]),2);
            }
            
            sum = sqrt(sum);

            if(sum < MINIMUM_SIFT_DISTANCE)
            {
                if(part == "part1.1")
                {
                    output.draw_line(desc1[i].col, desc1[i].row, desc2[j].col + input1.width(), desc2[j].row, color);
                }
                matches++;
            }
        }
    }
    
    // saving the sift output to a image file
    if(part == "part1.1") 
    {
        output.save("sift.png");
    }

    return matches; 
}

void img_display(const CImg<double>& img) {
	cimg_library::CImgDisplay main_disp(img, "temp");
	while (!main_disp.is_closed() ) 
		main_disp.wait();

}

int temp() {
#if 0
	CImg<double> img("lincoln.png");
	SqMatrix mat = SqMatrix::identity(3);

	Transformation t(mat);
	t.translate(30,30);
	t.rotate(-M_PI/6.0);
	CImg<double> result(transform_image<double>(img, t));
	cimg_library::CImgDisplay d1(img, "temp1");
	cimg_library::CImgDisplay d2(result, "temp2");
	result.save("result.png");

	while (!d1.is_closed() || !d2.is_closed()) {
		d1.wait();
		d2.wait();
	}
#else
	SqMatrix mat(3);
	mat(0,0) = 0.9;
	mat(0,1) = 0.258;
	mat(0,2) = -182;
	mat(1,0) = -0.153;
	mat(1,1) = 1.44;
	mat(1,2) = 58;
	mat(2,0) = -0.000306;
	mat(2,1) = 0.000731;
	mat(2,2) = 1;

	//mat = mat.transpose();

	Transformation t(mat);
	CImg<double> img("lincoln.png");
	CImg<double> result = transform_image<double>(img, t);
	img_display(result);

	return 0;
#endif
}

//-----------------------------------------------------------------------------------------

bool nearest_neighbour(const vector< double> &input1, const vector< double> &input2)
{
	for (int i = 0; i < 10; i++)
	{
	if (input1[i] != input2[i])
	{
	return false;
	}
	}
	return true;
}

void projection_function(CImg<double> input1, CImg<double> input2)
{
	int matches=0, maxHeight = 0;
	    const unsigned char color[] = {0, 255, 0};

	    CImg<double> greyScale1, greyScale2;

	    // converting the images to grayscale
	    if(input1.spectrum() == 1)
	        greyScale1 = input1;
	    else
	        greyScale1 = input1.get_RGBtoHSI().get_channel(2);

	    if(input2.spectrum() == 1)
	    	greyScale2 = input2;
	    else
		    greyScale2 = input2.get_RGBtoHSI().get_channel(2);

	    vector<SiftDescriptor> desc1 = Sift::compute_sift(greyScale1);
	    vector<SiftDescriptor> desc2 = Sift::compute_sift(greyScale2);
	    int size_input1=desc1.size();
	    int size_input2=desc2.size();
	    std::vector< double> X(128);

	    /*for(int i = 0; i < desc1.size(); i++)
	       {
	    	vector<int> row;
	           for(int k = 0; k < 128; k++)
	               {
	        	   	   row.push_back(unifRand());
	               }
	           X.push_back(row);
	       }*/

	    for(int k = 0; k < desc1.size(); k++)
	    	{
	    		X.push_back(unifRand());
	    	}

	    int w=255;
	    std::vector< double> F1(size_input1);
	    //CImg<double> F1(size_input1,size_input2);
	    std::vector< double> F2(size_input2);
	    //CImg<double> F2(size_input1,size_input2);
	    vector< double> F12;
	    for(int i=0;i<size_input1;i++)
	    {
	    for(int j=0;j<128;j++)
	    {
	    	F1[i]=F1[i]+((X[j]*desc1[i])/w);
	    	//F1(j,i) += (X(j, i) * desc1[i]/w );

	    }
	    }

	    for(int i=0;i <size_input2;i++)
	    	    {
	    	for(int j=0;j<128;j++)
	    	{
	    			F2[i]=F2[i]+((X[j]*desc2[i])/w);
	    			//F2(j,i) += (X(j, i) * desc2[i]/w );
	    	}

	    	    }

	    if(F1.size()!=F2.size()){
	            cout << "Vectors are not of the same size and hence the scalar product cannot be calculated" << endl;
	            return -1;
	        }

	    //Summary vectors
	    double distance = 0.0;
	    int lowest_index = -1, second_lowest_index;
	    map< int, pair<int, double> > image_map;
	    for (int i = 0; i < size_input1; i++)
	    	{
	    	double minimum_distance = std::numeric_limits<double>::infinity(), second_dist;
	    	int j;
	    	for (j = 0; j < size_input2; j++)
	    	{
	    	if (nearest_neighbour(F1, F2))
	    	{
	    	double ed = euclidean_distance(input1[i], input2[j]);
	    	printf("%d<>%d:%f\n", i, j, diff);
	    	if (minimum_distance > ed)
	    	{
	    	second_dist = minimum_distance;
	    	second_index = lowest_index;
	    	minimum_distance = ed;
	    	lowest_index = j;
	    	}
	    	else if (second_dist > diff)
	    	{
	    	second_dist = ed;
	    	second_lowest_index = j;
	    	}
	    	}
	    	}
	    	if (minimum_distance == std::numeric_limits<double>::infinity())
	    	{
	    	for (j = 0; j < size_input2; j++)
	    	{
	    	double ed = euclidean_distance(input1[i], input2[j]);
	    	if (minimum_distance > ed)
	    	{
	    	second_dist = minimum_distance;
	    	second_index = lowest_index;
	    	minimum_distance = ed;
	    	lowest_index = j;
	    	}
	    	else if (second_dist > ed)
	    	{
	    	second_dist = ed;
	    	second_lowest_index = j;
	    	}
	    	}
	    	}

	    	printf("dist: %f, dist2: %f\n", minimum_distance, second_dist);
	    		distance += minimum_distance;
	    		if (second_dist - minimum_distance > 100)
	    		{
	    		map< int, pair<int, double> >::iterator pos = match_map.find(min_index);
	    		if (pos == image_map.end())
	    		{
	    			image_map[minimum_index] = pair<int, double>(i, minimum_distance);
	    		}
	    		else if (pos->second.second > minimum_distance)
	    		{
	    			image_map[minimum_index] = pair<int, double>(i, minimum_distance);
	    		}

	    		//pairs.push_back(pair<int, int>(i, min_index));
	    		}

	    		}

	    		for (map< int, pair<int, double> >::iterator it = match_map.begin();
	    		it != image_map.end(); it++)
	    		{
	    		pairs.push_back(pair<int, int>(it->second.first, it->first));
	    		}

	    		for (int i = 0; i < pairs.size(); i++)
	    		{
	    		printf("%d -->> %d\n", pairs[i].first, pairs[i].second);
	    		}

	    		return dist;


}


float euclidean_distance(const SiftDescriptor &input1, const SiftDescriptor &input2)
{
	double distance = 0;
		for (int i = 0; i < 128; i++)
		{
			distance=distance+(input1.descriptor[i] - input2.descriptor[i]) * (input1.descriptor[i] - input2.descriptor[i]);
		}
		distance=sqrt(distance);
		return distance;
}

//
// Generate a random number between 0 and 1
// return a uniform number in [0,1].
double unifRand()
{
    return rand() / double(RAND_MAX);
}


//----------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
	return temp();
	try {
		if(argc < 2)
		{
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    a2-p1 part_id ..." << endl;
			return -1;
		}

		string part = argv[1];

		// Sift detector
		if(part == "part1.1")
		{
			if(argc < 4)
			{
				cout << "Insufficent number of arguments." << endl;
				return -1;
			}

			CImg<double> input1(argv[2]);
			CImg<double> input2(argv[3]);

			double matches = sift_matching(input1, input2, part);
			cout<<"The Number of Sift Matches: "<<matches<<endl;            
		}
		else if(part == "part1.2")
		{
			if(argc < 4)
			{
				cout << "Insufficent number of arguments." << endl;
				return -1;
			}

			CImg<double> input(argv[2]);

			// creating a vector to store the sift matches for each of the image
			vector<Image> images;

			// take in all the images from the cammand line and calculate sift matches of the descriptors
			for(int i=3; i<argc; i++)
			{
				CImg<double> img(argv[i]);
				Image I;
				I.setParameters(argv[i], sift_matching(input, img, part));
				images.push_back(I);
			} 

			// sort according to the sift match count 
			std::sort(images.begin(), images.end());

			// printing the results
			cout<<"Top matches for Input Image: "<<argv[2]<<endl<<"Image Name: \t\t"<<"Count: "<<endl;
			for(int i = 0; i < images.size(); i++)
			{
				cout<<images[i].getName()<<"\t\t  "<<images[i].getCount()<<endl;
			}
		}
		else if(part == "part2")
		{
			// do something here!
		}
		else
			throw std::string("unknown part!");

		// feel free to add more conditions for other parts (e.g. more specific)
		//  parts, for debugging, etc.
	}
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}


