#include <iostream>
#include <complex>
#include <cmath>
#include "solver.hpp"


using namespace std;
namespace solver{

	double solve(RealVariable r){
		double a=r.a, b=r.b, c=r.c;
		double ans;
		
		if(a == 0 && b == 0){
			if(c==0)
				ans = 0;
			else
				throw runtime_error("there is no answer for this equation!");
		}
		else if(a!=0){
			double b2= pow(b,2);
			if((b2)-4*a*c < 0){
				throw runtime_error("answer not a real number!");
			}
			else if((b2)-4*a*c == 0)
				return (-b/(2*a));
			else{
				double delta =sqrt( (b2)-4*a*c );
				ans = (-b-delta)/(2*a);

			}
		}
		else
			ans = -c/b;
		r.a=0;	
		r.b=1;
		r.c=0;
		
		return ans;
	}
	

	complex<double> solve(ComplexVariable im){
		complex<double> a=im.a, b=im.b, c=im.c;
		complex<double> ans;
				
		if(isZero(a)){		
			if(isZero(b)){
				if(isZero(c))
					ans = 0;
				else
					throw runtime_error("there is no answer for this equation!");
			}
			else{
				if(imag(b)==0)
					ans = complex<double>(-real(c)/real(b), -imag(c)/real(b));
				else{
					ans = c/b;
				}
			}
		}
		else{	
			if(isZero(b)){
				complex<double> temp = c/a;
				if(imag(temp) == 0){
					if(real(temp)>=0)
						ans = complex<double>(sqrt(real(temp)), 0);
					else
						ans = complex<double>(0, sqrt(-real(temp)));
				}
			}
		}
		im.a=0;	
		im.b=1;
		im.c=0;

		return ans;
	}
	
	//--------------------RealVariables--------------------
	//--------------------Operator + ----------------------
	RealVariable operator + (RealVariable r1, RealVariable r2){
		RealVariable x(r1.a+r2.a, r1.b+r2.b, r1.c+r2.c);
		return x;
	}
	RealVariable operator + (double r1, RealVariable r2){
		RealVariable x(r2.a, r2.b, r2.c+r1);
		return x;
	}
	RealVariable operator + (RealVariable r1, double r2){
		return r2+r1;
	}


	//--------------------Operator - ----------------------
	RealVariable operator - (RealVariable r1, RealVariable r2){
		RealVariable x(r1.a-r2.a, r1.b-r2.b, r1.c-r2.c);
		
		return x;
	}
	RealVariable operator - (double r1, RealVariable r2){
		RealVariable x(-r2.a,-r2.b,r1-r2.c);
		return x;
	}
	RealVariable operator - (RealVariable r1, double r2){
		RealVariable x(r1.a,r1.b,r1.c-r2);
		return x;
	}
	RealVariable operator - (RealVariable r1){
		return 0-r1;
	}

	//--------------------Operator * ----------------------
	RealVariable operator * (RealVariable r1, RealVariable r2){
		if((r1.a>0 && (r2.a>0 || r2.b>0)) || (r2.a>0 && (r1.a>0 ||r1.b>0)))
			throw runtime_error("Cannot be more than x^2");
		RealVariable x;
		x.a=r1.a*r2.c + r1.b*r2.b + r1.c*r2.a;
		x.b=r1.b*r2.c + r2.b*r1.c;
		x.c=r1.c*r2.c;
		return x;
	}
	RealVariable operator * (double r1, RealVariable r2){
		RealVariable x(r1*r2.a, r1*r2.b, r1*r2.c);
		return x;
	}
	RealVariable operator * (RealVariable r1, double r2){
		return r2*r1;
	}


	//--------------------Operator / ----------------------
	RealVariable operator / (RealVariable r1, RealVariable r2){
		if(r2.a == 0 && r2.b == 0 && r2.c == 0)
			throw runtime_error("Cannot divide by 0!");
		if(r2.a!=0 || r2.b!=0)
			throw runtime_error("Illegal input");
		RealVariable x(r1.a/r2.c, r1.b/r2.c, r1.c/r2.c);
		return x;
	}
	RealVariable operator / (double r1, RealVariable r2){
		if(r2.a == 0 && r2.b == 0 && r2.c == 0)
			throw runtime_error("Cannot divide by 0!");
		RealVariable x(0,0,r1/r2.c);
		return x;
	}
	RealVariable operator / (RealVariable r1, double r2){
		if(r2 == 0)
			throw runtime_error("Cannot divide by 0!");
		RealVariable x(r1.a/r2, r1.b/r2, r1.c/r2);
		return x;
	}

	//--------------------Operator ^ ----------------------
	RealVariable operator ^ (RealVariable r1, RealVariable r2){

		if(r2.a > 0 || r2.b > 0 || r1.a > 0 || r2.c < 0 || r2.c > 2)
			throw runtime_error("not legal input!");
		if(r2.c == 0){
			RealVariable x(0,0,1);
			return x;	
		}
		else if(r2.c == 1){
			return r1;
		}
		else{
			RealVariable x;
			x.a = pow(r1.b,2);
			x.b = 2*r1.b*r1.c;
			x.c = pow(r1.c,2);
			return x;
		}
	}

	RealVariable operator ^ (double r1, RealVariable r2){
		if(r2.a>0 || r2.b>0)
			throw runtime_error("Not a legal input!");
		RealVariable x(0,0,r1);	
		return x^r2;
	}
	RealVariable operator ^ (RealVariable r1, double r2){
		RealVariable x(0,0,r2);	
		return r1^x;
	}


	//--------------------Operator == ----------------------
	RealVariable operator == (RealVariable r1, RealVariable r2){
		return r2 - r1;
	}
	RealVariable operator == (double r1, RealVariable r2){
		return r2 - r1;
	}
	RealVariable operator == (RealVariable r1, double r2){
		return r2 - r1;
	}

	//---------------------------------ComplexVariables--------------------
	//--------------------Operator + ----------------------	
	ComplexVariable operator + (ComplexVariable c1, ComplexVariable c2){
		ComplexVariable y(c1.a+c2.a, c1.b+c2.b, c1.c+c2.c);
		return y;
	}
	ComplexVariable operator + (double c1, ComplexVariable c2){
		complex<double> c3 = complex<double>(c1+real(c2.c), imag(c2.c));
		ComplexVariable y(c2.a, c2.b, c3);
		return y;
	}
	ComplexVariable operator + (ComplexVariable c1, double c2){
		return c2+c1;
	}
	ComplexVariable operator + (complex<double> c1, ComplexVariable c2){
		complex<double> c3 = complex<double>(real(c1)+real(c2.c), imag(c1)+imag(c2.c));
		ComplexVariable y(c2.a, c2.b, c3);
		return y;
	}
	ComplexVariable operator + (ComplexVariable c1, complex<double> c2){
		return c2+c1;
	}
	ComplexVariable operator + (complex<double> c1, int c2){
		return c2+c1;
	}
	ComplexVariable operator + (int c1, complex<double> c2){
		complex<double> c3 = complex<double>(c1+real(c2), imag(c2));
		ComplexVariable y(0.0, 0.0, c3);
		return y;
	}
	ComplexVariable operator + (complex<double> c1, double c2){
		complex<double> c3 = complex<double>(c2+real(c1), imag(c1));
		ComplexVariable y(0.0, 0.0, c3);
		return y;
	}
	ComplexVariable operator + (double c1, complex<double> c2){
		return c2+c1;
	}



	//--------------------Operator - ----------------------
	ComplexVariable operator - (ComplexVariable c1, ComplexVariable c2){
		ComplexVariable y(c1.a-c2.a, c1.b-c2.b, c1.c-c2.c);
		
		return y;
	}
	ComplexVariable operator - (double c1, ComplexVariable c2){
		complex<double> c3 = complex<double>(c1-real(c2.c), -imag(c2.c));
		ComplexVariable y(-c2.a,-c2.b,c3);
		return y;
	}
	ComplexVariable operator - (ComplexVariable c1, double c2){
		complex<double> c3 = complex<double>(real(c1.c)-c2, imag(c1.c));
		ComplexVariable y(c1.a, c1.b, c3);
		return y;
	}
	ComplexVariable operator - (complex<double> c1, ComplexVariable c2){
		ComplexVariable y(-c2.a,-c2.b,c1-c2.c);
		return y;
	}
	ComplexVariable operator - (ComplexVariable c1, complex<double> c2){
		ComplexVariable y(c1.a, c1.b, c1.c-c2);
		return y;
	}
	ComplexVariable operator - (complex<double> c1, int c2){
		complex<double> c3 = complex<double>(real(c1)-c2, imag(c1));
		ComplexVariable y(0.0,0.0,c3);
		return y;
	}
	ComplexVariable operator - (int c1, complex<double> c2){
		complex<double> c3 = complex<double>(c1-real(c2), -imag(c2));
		ComplexVariable y(0.0,0.0,c3);
		return y;
	}
	ComplexVariable operator - (complex<double> c1, double c2){
		complex<double> c3 = complex<double>(real(c1)-c2, imag(c1));
		ComplexVariable y(0.0,0.0,c3);
		return y;
	}
	ComplexVariable operator - (double c1, complex<double> c2){
		complex<double> c3 = complex<double>(c1-real(c2),-imag(c2));
		ComplexVariable y(0.0,0.0,c3);
		return y;
	}
	ComplexVariable operator - (ComplexVariable c1){
		ComplexVariable y(-c1.a , -c1.b, -c1.c);
		return y;
	}







	//--------------------Operator * ----------------------
	ComplexVariable operator * (ComplexVariable c1, ComplexVariable c2){
		if((!isZero(c1.a) && (!isZero(c2.a) || !isZero(c2.b))) ||   (!isZero(c2.a) && (!isZero(c1.a) || !isZero(c1.b)))) 
			throw runtime_error("Cannot be more than x^2");


		ComplexVariable y;
		y.a=c1.a*c2.c + c1.b*c2.b + c1.c*c2.a;
		y.b=c1.b*c2.c + c2.b*c1.c;
		y.c=c1.c*c2.c;
		return y;



	}
	ComplexVariable operator * (double c1, ComplexVariable c2){
		complex<double> _a = complex<double> (c1*real(c2.a), c1*imag(c2.a));
		complex<double> _b = complex<double> (c1*real(c2.b), c1*imag(c2.b));
		complex<double> _c = complex<double> (c1*real(c2.c), c1*imag(c2.c));
		ComplexVariable y(_a, _b, _c);
		return y;
	}
	ComplexVariable operator * (ComplexVariable c1, double c2){
		return c2*c1;
	}
	ComplexVariable operator * (complex<double> c1, ComplexVariable c2){
		complex<double> _a = c2.a*c1;
		complex<double> _b = c2.b*c1;
		complex<double> _c = c2.c*c1;
		ComplexVariable y(_a,_b,_c);
		return y;
	}
	ComplexVariable operator * (ComplexVariable c1, complex<double> c2){
		return c2*c1;
	}
	ComplexVariable operator * (complex<double> c1, int c2){
		complex<double> c3 = complex<double>(c2*real(c1) , c2*imag(c1));
		ComplexVariable y(0,0,c3);
		return y;
	}
	ComplexVariable operator * (int c1, complex<double> c2){
		return c2*c1;
	}
	ComplexVariable operator * (complex<double> c1, double c2){
		complex<double> c3 = complex<double>(c2*real(c1) , c2*imag(c1));
		ComplexVariable y(0,0,c3);
		return y;
	}
	ComplexVariable operator * (double c1, complex<double> c2){
		return c2*c1;
	}










	//--------------------Operator / ----------------------
	ComplexVariable operator / (ComplexVariable c1, ComplexVariable c2){
		if(isZero(c2.a) && isZero(c2.b) && isZero(c2.c))
			throw runtime_error("Cannot divide by 0!");
		if(!isZero(c2.a) || !isZero(c2.b))
			throw runtime_error("Illegal input");
					
		complex<double> _a = c1.a/c2.c;
		complex<double> _b = c1.b/c2.c;
		complex<double> _c = c1.c/c2.c;
		ComplexVariable y(_a,_b,_c);
		return y;
	}
	ComplexVariable operator / (double c1, ComplexVariable c2){	
		if(isZero(c2.a) && isZero(c2.b) && isZero(c2.c))	
			throw runtime_error("Cannot divide by 0!");	
		complex<double> c3 = complex<double>(c1/real(c2.c), 0);
		ComplexVariable y(0,0,c3);
		return y;	
		
	}
	ComplexVariable operator / (ComplexVariable c1, double c2){
		if(c2==0)	
			throw runtime_error("Cannot divide by 0!");	
		complex<double> _a = complex<double>(real(c1.a)/c2, imag(c1.a)/c2);
		complex<double> _b = complex<double>(real(c1.b)/c2, imag(c1.b)/c2);
		complex<double> _c = complex<double>(real(c1.c)/c2, imag(c1.c)/c2);
		
		ComplexVariable y(_a,_b,_c);
		return y;
	}
	ComplexVariable operator / (complex<double> c1, int c2){
		if(c2==0)	
			throw runtime_error("Cannot divide by 0!");	
		complex<double> c3 = complex<double>(real(c1)/c2, imag(c1)/c2);
		ComplexVariable y(0,0,c3);
		return y;
	}

	complex<double> operator / (complex<double> c, complex<double> b){
		complex<double> b1 = complex<double>(real(b), -imag(b));
		complex<double> m = -c*b1;
		
		double b2 = pow(real(b),2) + pow(imag(b),2);
		complex<double> ans =complex<double> (real(m)/b2, imag(m)/b2);
		return ans;
	}








	//--------------------Operator ^ ----------------------
	ComplexVariable operator ^ (ComplexVariable c1, int c2){
		if(c2<0 || c2>3)
			throw runtime_error("Not a legal input!");
		else if(c2 == 0){
			ComplexVariable y(0,0,1);
			return y;
		}
		else if(c2 == 1){
			return c1;
		}
		else{
			ComplexVariable y;
			y.a = complex<double>(pow(real(c1.b),2)-pow(imag(c1.b),2), 2*real(c1.b)*imag(c1.b));
			y.b = complex<double>(real(c1.b)*real(c1.c)-imag(c1.b)*imag(c1.c) , real(c1.b)*imag(c1.c)+real(c1.c)*imag(c1.c));
			y.c = complex<double>(pow(real(c1.c),2)-pow(imag(c1.c),2),2*real(c1.c)*imag(c1.c));
			return y;
		}
	}
	

	





	//--------------------Operator == ----------------------
	ComplexVariable operator == (ComplexVariable c1, ComplexVariable c2){
		return c1-c2;
	}
	ComplexVariable operator == (double c1, ComplexVariable c2){
		return c1-c2;
	}
	ComplexVariable operator == (ComplexVariable c1, double c2){
		return c1-c2;
	}
	ComplexVariable operator == (complex<double> c1, ComplexVariable c2){
		return c1-c2;
	}
	ComplexVariable operator == (ComplexVariable c1, complex<double> c2){
		return c1-c2;
	}



	bool operator == (complex<double> c1, complex<double> c2){
		return (real(c1) == real(c2) && imag(c1) == imag(c2));
	}
	bool operator == (double c1, complex<double> c2){
		return (c1 == real(c2));
	}
	bool operator == (complex<double> c1, double c2){
		return (c2 == real(c1));
	}
	bool operator == (int c1, complex<double> c2){
		return (c1 == real(c2));
	}
	bool operator == (complex<double> c1, int c2){
		return (c2 == real(c1));
	}

	bool isZero(complex<double> x){
		if(real(x) == 0 && imag(x) == 0)
			return true;	
		return false;		
	}




}
