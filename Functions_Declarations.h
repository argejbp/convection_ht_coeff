using namespace std;

double Mean_Temperature(double T1, double T2);  //Listo
double convection_heat_coefficient (const double &Nu, const double &k, const double &L_ref);  //Listo

double Reynolds (double V, double rho, double v);  //Listo
double Inner_Friction_Factor (const double &Re);  //Listo
double CI_Pipes (const int &BC, const double &Re, const double &f, const double &L_ref, const double &Pr, double &Nu);  //Listo

double CE_Sphere (const double &Re); //Listo
double CE_Cylinder(const double &Re); //Listo
double CE_Flat_Plate(const double &Re, const double &Pr); //Listo

double Reyleigh(const double &L_ref, const double &Pr, const double &v, const double &Tm, const double &T1, const double &T2, const double g = 9.81); //Listo
double CN_Vertical_Plate(const double &L_ref, const double &Pr, const double &Ra); //Listo
double CN_Horizontal_Plate_1(const int &shape, const double &Ra); //Listo
double CN_Horizontal_Plate_2(const int &shape, const double &Ra); //Listo
double CN_Sphere(const double &Pr, const double &Ra); //Listo
double Horizontal_Cylinder(const double &Ra, const double &Pr); //Listo
double CN_Vertical_Enclosure(const double &H, const double &L_Encl, const double &Pr, const double &Ra); //Listo
double CN_Horizontal_Enclosure(const double &Pr, const double &Ra); //Listo