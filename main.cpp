#include "bsplinecdata.h"
using namespace std;

int main()
{
	// 1. declare a Cdata object named cdat using the default constructor
	Cdata cdat;

	// 2. declare an ifstream object connected to the file containing the data points (points.dat)
	ifstream ifs("dataSet2.txt");

	// 3. declare an ofstream object connected to output file for writing the evaluated points (BsplineCurve.dat)
	ofstream ofs("ResdataSet5.dat");

	// 4. read the data points from the file into cdat
	ifs >> cdat;

	// 5. declare a BsplineCdata object named bsp and initialise with cdat from 1
	BsplineCdata bsp(cdat);

	// 6. Use the ofstream object to write to the output file BsplineCurve.dat the BsplineCurve object returned by calling the method getBspline() for bsp
	ofs << bsp.getBspline();

	return 0;
}




	