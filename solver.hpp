#pragma once

#include <iostream>
#include <complex>


using namespace std;

namespace solver{

	class RealVariable{
		public:
		double a;
		double b;	
		double c;
		RealVariable(){
			a=0;
			b=1;
			c=0;
		}
		RealVariable(double _a, double _b, double _c){
			a=_a;
			b=_b;
			c=_c;
		}
		
	};
	class ComplexVariable{
		public:
		complex<double> a;
		complex<double> b;
		complex<double> c;
		ComplexVariable(){
			a=0;
			b=1;
			c=0;
		}
		ComplexVariable(double _a, double _b, complex<double> _c){
			a=complex<double> (_a,0);
			b=complex<double> (_b,0);
			c=_c;
		}
		ComplexVariable(complex<double> _a, complex<double> _b, complex<double> _c){
			a=_a;
			b=_b;
			c=_c;
		}
	};

	double solve(RealVariable r);
	double solve(RealVariable r, int real);

	complex<double> solve(ComplexVariable c);

	//-------------------RealVariables------------

	RealVariable operator + (RealVariable r1, RealVariable r2);
	RealVariable operator + (double r1, RealVariable r2);
	RealVariable operator + (RealVariable r1, double r2);

	RealVariable operator - (RealVariable r1, RealVariable r2);
	RealVariable operator - (double r1, RealVariable r2);
	RealVariable operator - (RealVariable r1, double r2);
	RealVariable operator - (RealVariable r1);


	RealVariable operator * (RealVariable r1, RealVariable r2);
	RealVariable operator * (double r1, RealVariable r2);
	RealVariable operator * (RealVariable r1, double r2);


	RealVariable operator / (RealVariable r1, RealVariable r2);
	RealVariable operator / (double r1, RealVariable r2);
	RealVariable operator / (RealVariable r1, double r2);


	RealVariable operator ^ (RealVariable r1, RealVariable r2);
	RealVariable operator ^ (double r1, RealVariable r2);
	RealVariable operator ^ (RealVariable r1, double r2);

	RealVariable operator == (RealVariable r1, RealVariable r2);
	RealVariable operator == (double r1, RealVariable r2);
	RealVariable operator == (RealVariable r1, double r2);

	//-----------------------ComplexVariables-----------------

	ComplexVariable operator + (ComplexVariable c1, ComplexVariable c2);
	ComplexVariable operator + (double c1, ComplexVariable c2);
	ComplexVariable operator + (ComplexVariable c1, double c2);
	ComplexVariable operator + (complex<double> c1, ComplexVariable c2);
	ComplexVariable operator + (ComplexVariable c1, std::complex<double> c2);
	ComplexVariable operator + (complex<double> c1, int c2);
	ComplexVariable operator + (int c1, complex<double> c2);
	ComplexVariable operator + (complex<double> c1, double c2);
	ComplexVariable operator + (double c1, complex<double> c2);

	ComplexVariable operator - (ComplexVariable c1, ComplexVariable c2);
	ComplexVariable operator - (double c1, ComplexVariable c2);
	ComplexVariable operator - (ComplexVariable c1, double c2);
	ComplexVariable operator - (complex<double> c1, ComplexVariable c2);
	ComplexVariable operator - (ComplexVariable c1, complex<double> c2);
	ComplexVariable operator - (complex<double> c1, int c2);
	ComplexVariable operator - (int c1, complex<double> c2);
	ComplexVariable operator - (complex<double> c1, double c2);
	ComplexVariable operator - (double c1, complex<double> c2);
	ComplexVariable operator - (ComplexVariable c1);

	ComplexVariable operator * (ComplexVariable c1, ComplexVariable c2);
	ComplexVariable operator * (double c1, ComplexVariable c2);
	ComplexVariable operator * (ComplexVariable c1, double c2);
	ComplexVariable operator * (complex<double> c1, ComplexVariable c2);
	ComplexVariable operator * (ComplexVariable c1, complex<double> c2);
	ComplexVariable operator * (complex<double> c1, int c2);
	ComplexVariable operator * (int c1, complex<double> c2);
	ComplexVariable operator * (complex<double> c1, double c2);
	ComplexVariable operator * (double c1, complex<double> c2);

	ComplexVariable operator / (ComplexVariable c1, ComplexVariable c2);
	ComplexVariable operator / (double c1, ComplexVariable c2);
	ComplexVariable operator / (ComplexVariable c1, double c2);
	ComplexVariable operator / (complex<double> c1, ComplexVariable c2);
	ComplexVariable operator / (ComplexVariable c1, complex<double> c2);
	ComplexVariable operator / (complex<double> c1, int c2);
	ComplexVariable operator / (int c1, complex<double> c2);
	ComplexVariable operator / (complex<double> c1, double c2);
	ComplexVariable operator / (double c1, complex<double> c2);
	complex<double> operator / (complex<double> c, complex<double> b);

	ComplexVariable operator ^ (ComplexVariable c1, int c2);

	ComplexVariable operator == (ComplexVariable c1, ComplexVariable c2);
	ComplexVariable operator == (double c1, ComplexVariable c2);
	ComplexVariable operator == (ComplexVariable c1, double c2);
	ComplexVariable operator == (complex<double> c1, ComplexVariable c2);
	ComplexVariable operator == (ComplexVariable c1, complex<double> c2);

	bool operator == (complex<double> c1, complex<double> c2);
	bool operator == (double c1, complex<double> c2);
	bool operator == (int c1, complex<double> c2);
	bool operator == (complex<double> c1, double c2);
	bool operator == (complex<double> c1, int c2);

	bool isZero(complex<double> x);
	
	

}
