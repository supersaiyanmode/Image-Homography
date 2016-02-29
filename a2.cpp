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
#include <list>
#include "CImg.h"
#include "Sift.h"
#include "Matrix.h"
#include "Transform.h"
#include <map>
#include "MappedCoordinates.h"
#include "HomographyEstimator.h"

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



double euclidean_distance(const SiftDescriptor &input1, const SiftDescriptor &input2)
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
CImg<double> compute_X()
{
	CImg<double> rand_mat(5, 128);
	rand_mat.rand(0,1.0);
	return rand_mat;
}

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

bool compare(const CImg<double> &input1, const CImg<double> &input2,int size1, int size2)
{
	for (int i = 0; i < 5; i++)
	{
	if (input1(i, size1) != input2(i, size2))
	{
	return false;
	}
	}
	return true;
}

double projection_function(const vector<SiftDescriptor> &input1, const vector<SiftDescriptor> &input2,const CImg<double> &X,vector< pair<int, int> > &pairs)
{
	double w = 255;
	int sub_iter = 5;
	int input_size = input1.size();
	int target_size = input2.size();
	CImg<double> source_vectors(sub_iter, input_size);
	for (int y = 0; y < input_size; y++) {
		for (int j = 0; j < sub_iter; j++) {
			for (int i = 0; i < 128; i++) {
				source_vectors(j, y) += X(j, i) * input1[i].descriptor[i] / w;
			}
			source_vectors(j, y) = (int) source_vectors(j, y);
		}
	}

	CImg<double> target_vectors(sub_iter, target_size);
	for (int y = 0; y < target_size; y++) {
		for (int j = 0; j < sub_iter; j++) {
			for (int i = 0; i < 128; i++) {
				target_vectors(j, y) += X(j, i) * input2[i].descriptor[i] / w;
			}
			target_vectors(j, y) = (int) target_vectors(j, y);
		}
	}

	double dist = 0.0;
	int min_index = -1, second_min_index;

	std::map<int, pair<int, double> > match_map;

	for (int i = 0; i < input_size; i++) {
		double min_dist = std::numeric_limits<double>::infinity(), second_min_dist;
		int j;
		for (j = 0; j < target_size; j++) {
			if (compare(source_vectors, target_vectors, i, j)) {
				double diff = euclidean_distance(input1[i], input2[j]);
				if (min_dist > diff) {
					second_min_dist = min_dist;
					second_min_index = min_index;
					min_dist = diff;
					min_index = j;
				} else if (second_min_dist > diff) {
					second_min_dist = diff;
					second_min_index = j;
				}
			}
		}
		if (min_dist == std::numeric_limits<double>::infinity()) {
			for (j = 0; j < target_size; j++) {
				double diff = euclidean_distance(input1[i], input2[j]);
				if (min_dist > diff) {
					second_min_dist = min_dist;
					second_min_index = min_index;
					min_dist = diff;
					min_index = j;
				} else if (second_min_dist > diff) {
					second_min_dist = diff;
					second_min_index = j;
				}
			}
		}
		dist += min_dist;
		if (second_min_dist - min_dist > 100) {
			map<int, pair<int, double> >::iterator pos = match_map.find(
					min_index);
			if (pos == match_map.end()) {
				match_map[min_index] = pair<int, double>(i, min_dist);
			} else if (pos->second.second > min_dist) {
				match_map[min_index] = pair<int, double>(i, min_dist);
			}
		}

	}
	for (map<int, pair<int, double> >::iterator it = match_map.begin();
			it != match_map.end(); it++) {
		pairs.push_back(pair<int, int>(it->second.first, it->first));
	}

	return dist;
}

void duplicate(CImg<double> &target, const CImg<double> &source)
{
	for (int i = 0; i < source.height(); i++)
	{
	for (int j = 0; j < source.width(); j++)
	{
	target(j, i, 0, 0) = source(j, i, 0, 0);
	target(j, i, 0, 1) = source(j, i, 0, 1);
	target(j, i, 0, 2) = source(j, i, 0, 2);
	}
	}
}

CImg<double> image_join(const CImg<double> &input1, const CImg<double> &input2)
{
	int coloumns = input1.width() + input2.width();
	int rows = max(input1.height(), input2.height());

	CImg<double> new_image(coloumns, rows, 1, 3);
	new_image = 0.0;

	duplicate(new_image, input2);
	new_image.shift(input1.width());
	duplicate(new_image, input1);

	return new_image;
}

//----------------------------------------------------------------------------------------


int main(int argc, char **argv)
{
	//return temp();
	try {
		if(argc < 2)
		{
			cout << "Insufficent number of arguments; correct usage:" << endl;
			cout << "    a2-p1 part_id ..." << endl;
			return -1;
		}

		string part = argv[1];

		// Sift detector
		if(part == "part1")
		{
			string inputFile = argv[2];
			CImg<double> input_image(inputFile.c_str());
			CImg<double> gray = input_image.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> descriptors = Sift::compute_sift(gray);

			list< pair<double, int> > ordered;

			vector< CImg<double> > final_targets;
			vector< vector< SiftDescriptor > > records;
			vector< vector< pair< int, int > > > pair_of_pics;
			const int first_i_target = 2;
			for (int file = first_i_target; file < argc; file ++)
			{
			final_targets.push_back(CImg<double>(argv[file]));
			CImg<double> &image_target = final_targets.back();
			CImg<double> gray = image_target.get_RGBtoHSI().get_channel(2);
			vector<SiftDescriptor> target_descriptors = Sift::compute_sift(gray);
			records.push_back(target_descriptors);

			pair_of_pics.push_back(vector< pair<int, int> >());
			vector< pair<int, int> > &vector_pairs = pair_of_pics.back();
			CImg<double> X = compute_X();
			double distance = projection_function(descriptors, target_descriptors, X, vector_pairs);

			list< pair<double, int> >::iterator final_list;
			for (final_list = ordered.begin(); final_list != ordered.end(); final_list++)
			{
			if (final_list->first >= distance)
			{
			break;
			}
			}

			ordered.insert(final_list, pair<double, int>(distance, file - first_i_target));
			}

			for (list<pair<double, int> >::iterator it = ordered.begin();
				it != ordered.end(); it++)
				{
				int ind = it->second;
				CImg<double> &target_image = final_targets[ind];
				vector< pair<int, int> > &pairs = pair_of_pics[ind];
				gray = target_image.get_RGBtoHSI().get_channel(2);
				vector<SiftDescriptor> &target_descriptors = records[ind];
				CImg<double> pair_image = image_join(input_image, target_image);
				const unsigned char color[] = { 0,255,0 };
				for (int i = 0; i < pairs.size(); i++)
				{
				int x1 = descriptors[pairs[i].first].col;
				int y1 = descriptors[pairs[i].first].row;
				int x2 = target_descriptors[pairs[i].second].col;
				int y2 = target_descriptors[pairs[i].second].row;
				x2 += input_image.width();
				pair_image.draw_line(x1, y1, x2, y2, color);
				}

				CImgDisplay diaplay_pair(pair_image, "Joined");
				while (!diaplay_pair.is_closed())
				{
					diaplay_pair.wait();
				}

				pair_image.save("joined.png");
				}
		}
		else if(part == "part2")
		{


		}
		else
		{
			throw std::string("unknown part!");
		}
	}
	catch(const string &err) {
		cerr << "Error: " << err << endl;
	}
}


