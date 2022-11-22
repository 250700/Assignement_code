#include "cdata.h"
#include <cmath>
#define M_PI

// constructors
Cdata::Cdata() : n(0), pts(), param(), tang() {}

Cdata::Cdata(int num, const Vector<Point3D>& data) : n(num), pts(data), param(n), tang(n)
{
	Foley();
	computeTangents();
}


// uniform parameterisation
void Cdata::uniform()
{
	// compute the parameter values in the vector data member param according to the uniform method (page 25)
	for (int i = 0; i<n; i++) {
		param[i] = (double)i;
		std::cout << param[i] << "\n";
	
	}
}

// chord length parameterisation
void Cdata::chord()
{
	// compute the parameter values in the vector data member param according to the chord length method (page 26)
	param[0] = 0.0;

	double sum = 0.0;

	for (int i = 0; i<n - 1; i++) {
		double dist = sqrt(pow(pts[i + 1].getX() - pts[i].getX(), 2) +
			pow(pts[i + 1].getY() - pts[i].getY(), 2) +
			pow(pts[i + 1].getZ() - pts[i].getZ(), 2));
		sum = sum + dist;
		param[i + 1] = sum;
		std::cout << param[i+1] << "\n";
	}
}

// centripetal parameterisation
void Cdata::centripetal()
{
	// compute the parameter values in the vector data member param according to the centripetal method (page 26)
	param[0] = 0.0;

	double sum = 0.0;

	for (int i = 0; i<n - 1; i++) {
		double dist = sqrt(pow(pts[i + 1].getX() - pts[i].getX(), 2) +
			pow(pts[i + 1].getY() - pts[i].getY(), 2) +
			pow(pts[i + 1].getZ() - pts[i].getZ(), 2));
		sum = sum + sqrt(dist);
		param[i + 1] = sum;
		std::cout << param[i + 1] << "\n";
	}
}
void Cdata::Foley()
{
	param[0] = 0.0;
	param[1] = 0.0;
	for (int i = 1; i < n-2; i++) {
		//std::cout << pts[i].getX()<< "\n";
		double current_dist= pow(pts[i + 1].getX() - pts[i].getX(), 2) +pow(pts[i + 1].getY() - pts[i].getY(), 2) +pow(pts[i + 1].getZ() - pts[i].getZ(), 2);
		//std::cout << current_dist << "\n";
		double previous_dist= pow(pts[i].getX() - pts[i-1].getX(), 2) +pow(pts[i].getY() - pts[i-1].getY(), 2) +pow(pts[i].getZ() - pts[i-1].getZ(), 2);
		//std::cout << previous_dist << "\n";
		double next_dist= pow(pts[i + 2].getX() - pts[i+1].getX(), 2) +pow(pts[i + 2].getY() - pts[i+1].getY(), 2) +pow(pts[i + 2].getZ() - pts[i+1].getZ(), 2);
		//std::cout << next_dist << "\n";
		double next_dot_product = pts[i + 1].getX() * pts[i].getX() + pts[i + 1].getY() * pts[i].getY() + pts[i + 1].getZ() * pts[i].getZ();
		//std::cout << next_dot_product << "\n";
		double current_dot_product = pts[i].getX() * pts[i-1].getX() + pts[i].getY() * pts[i-1].getY() + pts[i].getZ() * pts[i-1].getZ();
		//std::cout << current_dot_product << "\n";
		double previous_norm = sqrt(pow(pts[i - 1].getX(), 2) + pow(pts[i - 1].getY(), 2) + pow(pts[i - 1].getZ(), 2));
		//std::cout << previous_norm << "\n";
		double current_norm = sqrt(pow(pts[i].getX(), 2) + pow(pts[i].getY(), 2) + pow(pts[i].getZ(), 2));
		//std::cout << current_norm << "\n";
		double next_norm = sqrt(pow(pts[i + 1].getX(), 2) + pow(pts[i + 1].getY(), 2) + pow(pts[i + 1].getZ(), 2));
		//std::cout << next_norm << "\n";
		double current_teta = acos(current_dot_product / (current_norm * previous_norm));
		//std::cout << current_teta << "\n";
		double next_teta = acos(next_dot_product / (next_norm * current_norm));
		//std::cout << next_teta << "\n";
		double current_adjust_teta = (0.5 * M_PI + (M_PI - current_teta) - fabs(0.5 * M_PI - (M_PI - current_teta)) / 2);
		//std::cout << current_adjust_teta << "\n";
		double next_adjust_teta= (0.5 * M_PI + (M_PI - next_teta) - fabs(0.5 * M_PI - (M_PI - next_teta)) / 2);
		//std::cout << next_adjust_teta << "\n";
		double dif = current_dist * (1 + 1.5 * (current_adjust_teta * previous_dist) / (current_dist + previous_dist) + 1.5 * (next_adjust_teta * next_dist) / (current_dist + next_dist));
		std::cout << dif << "\n";
		param[i+2] = dif;
		//std::cout << dif << "\n";
		//param[i+2] = (double)i;
	}
	
}

// ACCESSOR METHODS
int Cdata::getNum() const
{
	return n;
}

Vector<Point3D> Cdata::getPoints() const
{
	return pts;
}

Vector<double> Cdata::getParam() const
{
	return param;
}


Vector<Point3D> Cdata::getTang() const
{
	return tang;
}


// READ/WRITE METHODS
void Cdata::readfile(std::ifstream& ifs)  // read from a file
{
	int num;
	ifs >> num;

	Vector<Point3D> data(num);
	ifs >> data;

	*this = Cdata(num, data);
}

void Cdata::writefile(std::ofstream& ofs) const // write to a file
{
	ofs << n << "\n";
	ofs << pts;
	ofs << param;
	ofs << tang;
}

void Cdata::write(std::ostream& os) const // write to the screen
{
	os << n << "\n";
	os << pts;
	os << param;
	os << tang;
}

void Cdata::read(std::istream &is)  // read from the keyboard
{
	std::cout << "number of points\n";
	int num;
	is >> num;

	std::cout << "input points\n";
	Vector<Point3D> data(num);
	is >> data;

	*this = Cdata(num, data);
}


// compute the tangent vectors according to Bessel's formula
void Cdata::computeTangents() 
{
	double h1, h2;
	Point3D temp;

	/* compute Bessel slope for the first data point */
	h1 = param[1] - param[0];
	h2 = param[2] - param[1];

	temp = (((pts[2] - pts[1]) / h2)*h1 + ((pts[1] - pts[0]) / h1)*h2) / (h1 + h2);
	tang[0] = (2 / h1)*(pts[1] - pts[0]) - temp;

	/* compute Bessel slope for the last data point */
	h1 = param[n - 2] - param[n - 3];
	h2 = param[n - 1] - param[n - 2];

	temp = (((pts[n - 1] - pts[n - 2]) / h2)*h1 + ((pts[n - 2] - pts[n - 3]) / h1)*h2) / (h1 + h2);
	tang[n-1] = 2 * (pts[n - 1] - pts[n - 2]) / h2 - temp;

	// compute remaining tangents
	for (int i = 1; i < n - 1; i++) {
		double hi1 = param[i + 1] - param[i];
		double hi2 = param[i] - param[i - 1];
		tang[i] = (((pts[i+1] - pts[i]) / hi1)*hi2 + ((pts[i] - pts[i-1]) / hi2)*hi1) / (hi1 + hi2);
	}
}
