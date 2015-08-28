//#include "TFormulaLite.h"
#include "TFormula.h"
#include "TStopwatch.h"
#include "TRandom.h"
#include "TSystem.h"
#include <complex>
#include <vector>
#include <typeinfo>
#include <TROOT.h>
#include "Vc/Vc"
#include <array>
#include <Vc/Allocator>
//#include <Vc/vector.h>
#define V_SIZE 2

int main()
{
gSystem->Load("libVcshared");
gROOT->ProcessLine("#include <Vc/Vc>");

const int n = 64;
const int n2 = n/Vc::double_v::Size;
const int n3 = n/Vc::float_v::Size;

std::vector<double> p={1.,10.,100.,1000.,1.,1.,1.,1.,1.,1.,100.,1.,1.,1.};
std::vector<float> p2={1.f,10.f,100.f,1000.f,1.f,1.f,1.f,1.f,1.f,1.f,100.f,1.f,1.f,1.f};

TStopwatch w;
Int_t vSize = 64;


std::vector<double> points(n);
for (int j = 0; j < n; ++j) 
{
	points[j] = gRandom->Rndm();
//	std::cout << "Random point: " << points[j] << std::endl;
}

std::vector<float> pointsF(n);
for (int j = 0; j < n; ++j) 
{
	pointsF[j] = gRandom->Rndm();
//	std::cout << "Random float point: " << pointsF[j] << std::endl;
}

std::cout << std::endl;
std::vector<Vc::double_v> points2(n2);
for (int j = 0; j < n2; ++j) 
{
	points2[j] = gRandom->Rndm();
//	std::cout << "Random point: " << points2[j] << std::endl;
}



Vc::double_v x;
for (int j = 0; j < 2; ++j) 
{
	x[j] = j+1;
//	std::cout << "x[" << j << "] = " << x[j] << std::endl;
}
TFormula<double> f("f", "x+2.");

double result = f.Eval(2.,p);
std::cout << "result = " << result << std::endl;


Vc::double_v y;
for (int j = 0; j < 2; ++j) 
{
	y[j] = j+1;
//	std::cout << "y[" << j << "] = " << y[j] << std::endl;
}
std::vector<float> z(4);
for (int j = 0; j < 4; ++j) 
{
	z[j] = 1;
//	std::cout << "z[" << j << "] = " << z[j] << std::endl;
}
std::vector<float> t(4);
for (int j = 0; j < 4; ++j) 
{
	t[j] = 2;
//	std::cout << "t[" << j << "] = " << t[j] << std::endl;
}
//std::cout<< "Size of x = " << x.size() << std::endl;
std::array<double,2> x2 = {1, 2};
std::vector<float> x3 = {1, 2, 3, 4, 5};


std::vector<double> var[2];
for (int i = 0; i < 4; ++i)
{
	var[i][0] = x[i];
}
for (int i = 0; i < 4; ++i)
{
	var[i][1] = y[i];
}

for (int i = 0; i < 2; ++i)
{
	for(int j = 0; j < 4; ++j)
	{
		std::cout << "var[" << i << "][" << j << "] = " << var[i][j] << std::endl;
	} 
}

Vc::double_v result2 = h.Eval(x);
for (int j = 0; j < 2; ++j) 
	std::cout << result2[j] << std::endl;
//std::cout << "result = " << result << std::endl;
std::cout << std::endl;



//sin(x)/x+x*x*x-2*x+cos(2*x)


TFormula<double> g("g", "x/2+10*x");
std::cout << "TEST TFormula" << std::endl;
w.Start();
for (int j = 0; j < n; ++j) 
	{
		volatile double y = g.Eval(points[j]);
//		std::cout << "g.Eval(points[j]) = " << g.Eval(points[j]) << std::endl;	
	}
w.Print();

/*
std::cout << std::endl;
std::cout << "TEST TFormula double" << std::endl;
TFormula<double> f1("f1", "x/2+10*x");
w.Start();
for (int j = 0; j < n; ++j) 
	{
		volatile double y = f1.Eval(points[j]);
//		std::cout << "f1.Eval(points[j]) = " << f1.Eval(points[j]) << std::endl;
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormula float" << std::endl;
TFormula<float> f2("f2", "x/2+10*x");
w.Start();
for (int j = 0; j < n; ++j) 
	{
		volatile float y = f2.Eval(pointsF[j]);
//		std::cout << "f2.Eval(pointsF[j]) = " << f2.Eval(pointsF[j]) << std::endl;
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormula Vc::double_v" << std::endl;
TFormula<Vc::double_v> f3("f3", "x/2+10*x");
std::vector<Vc::double_v> d(n2);
for (int j = 0; j < n2; j++) 
	{
		for (int k = 0; k < Vc::double_v::Size; k++) 
		{
			d[j][k] = points[j*Vc::double_v::Size+k];
		}
	}
w.Start();
for (int j = 0; j < n2; j++) 
	{
		volatile Vc::double_v y = f3.Eval(d[j]);
	}	
//		std::cout << "f3.Eval(d) = " << f3.Eval(d) << std::endl;	
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormula Vc::float_v" << std::endl;
TFormula<Vc::float_v> f4("f4", "x/2+10*x");
std::vector<Vc::float_v> d2(n2);
for (int j = 0; j < n3; j++) 
	{
		for (int k = 0; k < Vc::float_v::Size; k++) 
		{
			d2[j][k] = points[j*Vc::float_v::Size+k];
		}
	}
w.Start();
for (int j = 0; j < n3; j++) 
	{
		volatile Vc::float_v y = f4.Eval(d2[j]);
	}	
//		std::cout << "f3.Eval(d) = " << f3.Eval(d) << std::endl;	
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormula std::vector<double>" << std::endl;
TFormula<std::vector<double>> f5("f5", "x/2+10*x");
std::vector<std::vector<double> > d5(n/vSize);
Int_t n5 = d5.size();
for (int j = 0; j < n5; j++) 
	{
		d5[j].resize(vSize);
		for (int k = 0; k < vSize; k++) 
		{
			if( j*vSize+k >= n )
				break;
			d5[j][k] = points[j*vSize+k];
		}
	}
w.Start();
for (int j = 0; j < n5; j++) 
	{

		volatile std::vector<double> y = f5.Eval(d5[j]);		
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormula std::vector<float>" << std::endl;
TFormula<std::vector<float>> f6("f6", "x/2+10*x");
std::vector<std::vector<float> > d6(n/vSize);
Int_t n6 = d6.size();
for (int j = 0; j < n6; j++) 
	{
		d6[j].resize(vSize);
		for (int k = 0; k < vSize; k++) 
		{
			if( j*vSize+k >= n )
				break;
			d6[j][k] = points[j*vSize+k];
		}
	}
w.Start();
for (int j = 0; j < n6; j++) 
	{

		volatile std::vector<float> y = f6.Eval(d6[j]);		
	}
w.Print();
std::cout << std::endl;


/*
TFormula h("h", "[0]*sin(x)+[1]*x*x*x");
std::cout << "TEST TFormula parameteres" << std::endl;
h.SetParameter(0,1.);
h.SetParameter(1,10.);
w.Start();
for (int j = 0; j < n; ++j) 
	{
		volatile double y = h.Eval(points[j]);
//		std::cout << "h.Eval(points[j]) = " << h.Eval(points[j]) << std::endl;	
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormulaLite double parameteres" << std::endl;
TFormulaLite<double> h1("h1", "p[0]*sin(x)+p[1]*x*x*x");
w.Start();
for (int j = 0; j < n; ++j) 
	{
		volatile double y = h1.Eval(points[j],p);
//		std::cout << "f1.Eval(points[j]) = " << f1.Eval(points[j]) << std::endl;
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormulaLite float parameteres" << std::endl;
TFormulaLite<float> h2("h2", "p[0]*sin(x)+p[1]*x*x*x");
w.Start();
for (int j = 0; j < n; ++j) 
	{
		volatile float y = h2.Eval(pointsF[j],p);
//		std::cout << "f2.Eval(pointsF[j]) = " << f2.Eval(pointsF[j]) << std::endl;
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormulaLite Vc::double_v parameteres" << std::endl;
TFormulaLite<Vc::double_v> h3("h3", "p[0]*sin(x)+p[1]*x*x*x");
std::vector<Vc::double_v> l(n2);
for (int j = 0; j < n2; j++) 
	{
		for (int k = 0; k < Vc::double_v::Size; k++) 
		{
			l[j][k] = points[j*Vc::double_v::Size+k];
		}
	}
w.Start();
for (int j = 0; j < n2; j++) 
	{
		volatile Vc::double_v y = h3.Eval(l[j],p);
	}	
//		std::cout << "f3.Eval(d) = " << f3.Eval(d) << std::endl;	
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormulaLite Vc::float_v parameteres" << std::endl;
TFormulaLite<Vc::float_v> h4("h4", "p[0]*sin(x)+p[1]*x*x*x");
std::vector<Vc::float_v> l2(n2);
for (int j = 0; j < n3; j++) 
	{
		for (int k = 0; k < Vc::float_v::Size; k++) 
		{
			l2[j][k] = points[j*Vc::float_v::Size+k];
		}
	}
w.Start();
for (int j = 0; j < n3; j++) 
	{
		volatile Vc::float_v y = h4.Eval(l2[j],p2);
	}	
//		std::cout << "f3.Eval(d) = " << f3.Eval(d) << std::endl;	
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormulaLite std::vector<double> parameteres" << std::endl;
TFormulaLite<std::vector<double>> h5("h5", "p[0]*sin(x)+p[1]*x*x*x");
std::vector<std::vector<double> > l5(n/vSize);
Int_t n5 = l5.size();
for (int j = 0; j < n5; j++) 
	{
		l5[j].resize(vSize);
		for (int k = 0; k < vSize; k++) 
		{
			if( j*vSize+k >= n )
				break;
			l5[j][k] = points[j*vSize+k];
		}
	}
w.Start();
for (int j = 0; j < n5; j++) 
	{

		volatile std::vector<double> y = h5.Eval(l5[j],p);		
	}
w.Print();
std::cout << std::endl;
std::cout << "TEST TFormulaLite std::vector<float> parameteres" << std::endl;
TFormulaLite<std::vector<float>> h6("h6", "p[0]*sin(x)+p[1]*x*x*x");
std::vector<std::vector<float> > l6(n/vSize);
Int_t n6 = l6.size();
for (int j = 0; j < n6; j++) 
	{
		l6[j].resize(vSize);
		for (int k = 0; k < vSize; k++) 
		{
			if( j*vSize+k >= n )
				break;
			l6[j][k] = points[j*vSize+k];
		}
	}
w.Start();
for (int j = 0; j < n6; j++) 
	{

		volatile std::vector<float> y = h6.Eval(l6[j],p);		
	}
w.Print();
std::cout << std::endl;
*/

return 0;

}
