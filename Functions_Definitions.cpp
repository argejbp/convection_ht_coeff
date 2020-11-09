#include <iostream>
#include <cmath>
#include "Functions_Declarations.h"
using namespace std;


/* 

	Functions for calculations of Heat Transfer Coefficients

*/

///////////////////////////////////////////////////////////////////////////

double Mean_Temperature(double T1, double T2){
	double Tm;
	Tm = (T1 + T2)/2;
	return Tm;
}

///////////////////////////////////////////////////////////////////////////

double Reynolds (double V, double rho, double v, double L_ref){
	double Re;
	Re = V*rho*L_ref/v;
	return Re;
}

///////////////////////////////////////////////////////////////////////////

double Inner_Friction_Factor(const double &Re) {
	double f;
	if (Re <= 2000) {
		f = 64 / Re;
	}
	else if (Re >= 3000) {
		f = pow((0.790*log(Re) - 1.64), -2); //Petukhov First Equation
	}
	return f;
}

///////////////////////////////////////////////////////////////////////////

double convection_heat_coefficient (const double &Nu, const double &k, const double &L_ref) {
	double h;
	h = Nu * k / L_ref;
	return h;
}

///////////////////////////////////////////////////////////////////////////

double CI_Pipes (const int &BC, const double &Re, const double &f, const double &L_ref, const double &Pr) {
	double Nu;
	if (Re <= 2000){
		switch(BC) {
			case 1: 			//Wall Heat Flow Boundary Condition
			Nu = 3.66;
			break;

			case 2:				//Wall Temperature Boundary Condition
			Nu = 4.36;
			break;
		}
	}
	else if (Re > 2000) {
		Nu = ((f / 8)*(Re - 1000)*Pr) / (1 + 12.7*pow((f / 8), 0.5)*(pow(Pr, 0.6666667) - 1));
	}
	return Nu;
}

///////////////////////////////////////////////////////////////////////////


double Reyleigh(const double &L_ref, const double &Pr, const double &v, const double &Tm, const double &T1, const double &T2, const double g = 9.81){
	double Ra;
	Ra = g * (1 / (Tm + 273.15))*(abs(T1 - T2))*pow(L_ref, 3)* Pr / pow(v, 2);
	return Ra;
}

///////////////////////////////////////////////////////////////////////////

double CN_Vertical_Plate(const double &L_ref, const double &Pr, const double &Ra){
	double Nu;
	if (Ra >= pow (10, 3) && Ra <= pow (10, 13)){
		Nu = pow((0.825 + 0.387*pow(Ra, 0.166667)) / (pow((1 + pow(0.492 / Pr, 0.5625)), 0.29629629629)), 2);
		return Nu;
	}

	else {
		cout << "There is no correlation for this Rayleigh Number. Ra = " << Ra;
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////

double CN_Horizontal_Plate_1(const int &shape, const double &Ra){
	double Nu;

	if (Ra >= pow(10, 4) && Ra <= pow(10, 7)) {
			Nu = 0.54*pow(Ra, 0.25);
			return Nu;
	}

	else if (Ra >= pow(10, 7) && Ra <= pow(10, 11)) {
		Nu = 0.15*pow(Ra, 0.333333);
		return Nu;
	}

	else {
		cout << "There is no correlation for this Rayleigh Number. Ra = " << Ra;
		return 0;
	}

}

///////////////////////////////////////////////////////////////////////////


double CN_Horizontal_Plate_2(const int &shape, const double &Ra){
	double Nu;

	if (Ra >= pow(10, 5) && Ra <= pow(10, 11)) {
		Nu = 0.27*pow(Ra, 0.25);
		return Nu;
	}

	else {
		cout << "There is no correlation for this Rayleigh Number. Ra = " << Ra;
		return 0;
	}


}

///////////////////////////////////////////////////////////////////////////

double CN_Sphere(const double &Pr, const double &Ra){
	double Nu;
	if (Ra <= 1E+11 || Pr >= 0.7) {
		Nu = 2 + (0.589*pow(Ra, 0.25)) / pow((1 + pow((0.469 / Pr), 0.5625)), 0.44444444);
	return Nu;
	}
	else {
		cout << "There is no correlation for this Rayleigh Number. Ra = " << Ra;
		return 0;
	}	
}

///////////////////////////////////////////////////////////////////////////

double Horizontal_Cylinder(const double &Ra, const double &Pr){
	double Nu;
	if (Ra <= 1E+12) {
		Nu = pow((0.6 + (0.387*pow(Ra, 0.166667)) / (pow(1 + pow(0.559 / Pr, 0.5625), 0.2962962963))), 2);
		return Nu;
	}

	else {
		cout << "There is no correlation for this Rayleigh Number. Ra = " << Ra;
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////

double CN_Vertical_Enclosure(const double &H, const double &L_Encl, const double &Pr, const double &Ra){
	double Nu;
	if (H / L_Encl >= 1 || H / L_Encl < 2) {
		if (Ra*Pr / (0.2 + Pr) > 1E+03) {
			Nu = 0.18*pow((Pr / (0.2 + Pr)*Ra), 0.29);
			return Nu;
		}
	}

	if (H / L_Encl >= 2 || H / L_Encl < 10) {
		if (Ra <= 1E+10) {
			Nu = 0.22*pow((Pr / (0.2 + Pr)*Ra), 0.28)*pow(H / L_Encl, -0.25);
			return Nu;
		}
	}

	else {
		cout << "There is no available correlation programmed for these conditions";
		cout << "\nRa = " << Ra;
		cout << "\nH/L = " << H / L_Encl;
		cout << "\nPr = " << Pr;
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////

double CN_Horizontal_Enclosure(const double &Pr, const double &Ra){
	double Nu;
	if (Ra <= 1E+12) {
		Nu = pow((0.6 + (0.387*pow(Ra, 0.166667)) / (pow(1 + pow(0.559 / Pr, 0.5625), 0.2962962963))), 2);
		return Nu;
	}

	else {
		cout << "There is no correlation for this Rayleigh Number. Ra = " << Ra;
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////

double CE_Flat_Plate(const double &Re, const double &Pr, const int &a){
	
	double Nu, Cf;
	
	if (Re < 5E+05 || Pr >0.6 || Pr <= 60) {
		Cf = 1.33 / pow(Re, 0.5);
		switch (a){
			case 0:
			Nu = 0.664*pow(Re, 0.5)*pow(Pr, 0.33333); // Constant temperature (Laminar)
			break;

			case 1:
			Nu = 0.453*pow(Re, 0.5)*pow(Pr, 0.33333); // Constant heat flux (Laminar)
			break;
		}

	}

	if (Re >= 5E+05) {
		Cf = 0.074 / pow(Re, 0.2);
		switch (a) {
		case 0:
			Nu = 0.037*pow(Re, 0.8)*pow(Pr, 0.33333); // Constant temperature (Turbulent)
			break;
		case 1:
			Nu = 0.0308*pow(Re, 0.8)*pow(Pr, 0.33333); // Constant heat flux (Turbulent)
			break;
		}

	}
	
	return Nu;
}

///////////////////////////////////////////////////////////////////////////

double CE_Cylinder(const double &Re, const double &Pr){
	double Nu;
	if (Re*Pr > 0.2) {
		Nu = 0.3 + 0.62*pow(Re, 0.5)*pow(Pr, 0.33333) / (pow((1 + pow((0.4 / Pr), 0.66667)), 0.25))*(pow(1 + pow((Re / 282000), 0.625), 0.8));
		return Nu;
	}

	else {
		cout << "There is no available correlation for those conditions";
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////

double CE_Sphere (const double &Re, const double &Pr, const double &miu_inf, const double &miu_s){
	double Nu;
	if (Re >= 3.5 || Re <= 80000 || Pr >= 0.7 || Pr <= 380 || miu_inf / miu_s >= 1 || miu_inf / miu_s <= 380) {
		Nu = 2 + (0.4*pow(Re, 0.5) + 0.06*pow(Re, 0.66667))*pow(Pr, 0.4)*(miu_inf / miu_s);
		return Nu;
	}

	else {
		cout << "There is no available correlation for those conditions";
		return 0;
	}
}

///////////////////////////////////////////////////////////////////////////