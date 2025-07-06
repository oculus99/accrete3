// version 0000.0004aaaaa

#include <iostream>
#include <fstream>   // For std::ofstream
#include <cstring>   // For strcmp (C-style string functions)
#include <cmath>     // For math functions like pow, log10, floor, sqrt, fabs, cbrt, M_PI
#include <cstdlib>   // For atoi, qsort, srand48, drand48, lrand48, time
#include <cstdio>    // For sprintf (needed if using sprintf, though std::to_string or string streams are preferred in C++)
#include <iomanip>   // For std::fixed and std::setprecision
#include <string>    // For std::string


// The original code included <rand48.h> which is non-standard.
// These functions are usually declared in <cstdlib> or <stdlib.h>
// for POSIX systems. We'll assume they are available.
extern "C" {
    double drand48();
    long    lrand48();
    void    srand48(long);
};

// --- Constants for conversions and calculations ---
const double SOLAR_MASS_TO_EARTH_MASS = 333000.0; // 1 Solar Mass = approx 333,000 Earth Masses
const double EARTH_RADIUS_METERS = 6.371e6;       // Earth Radius in meters
const double EARTH_MASS_KG = 5.972e24;            // Earth Mass in kilograms
const double JUPITER_DENSITY_KG_M3 = 1330.0;      // Jupiter density kg/m^3
const double EARTH_DENSITY_KG_M3 = 5510.0;        // Earth density kg/m^3
const double ICE_DENSITY_KG_M3 = 917.0;         // Density of ice (water ice) in kg/m^3
const double AU_TO_METERS = 1.496e11;             // 1 AU in meters
const double SUN_SURFACE_TEMP_K = 5778.0;         // Sun's effective surface temperature in Kelvin
const double SUN_RADIUS_METERS = 6.957e8;         // Sun's radius in meters

// Albedo estimates
const double ALBEDO_ROCKY_PLANET = 0.3;
const double ALBEDO_GAS_GIANT = 0.5;
const double ALBEDO_ICY_PLANET = 0.7; // Higher albedo for icy worlds

// New: Snow Line definition
const double SNOW_LINE_AU = 2.5; // Lumiraja AU:na (2.5 * ekosfäärin etäisyys, joka Auringolla 1 AU)

// Simulation parameters
struct DoleParams {
    double
    A,
    B,
    Alpha,
    Beta,
    Eccentricity,
    Gamma,
    K,
    MassSol,
    MassSun,
    MassNuclei,
    W,
    KMigration,           // Migration coefficient during accretion
    MigrationFinalFraction; // New: Final fraction for post-accretion migration (QMigration)

   DoleParams();
    DoleParams( double pA, double pAlpha,double pBeta,double pEccentricity,double pGamma,double pK,double pMassSol,double pMassSun,double pMassNuclei,double pB,double pW, double pKMigration, double pMigrationFinalFraction);
    
    double density(double au);
    double gasdensity(double au, double massratio);
    double masscritical(double au);
    double lowbound(double radius, double margin);
    double highbound(double radius, double margin);
};

// Initialize to defaults. See Sagan's article for insight into changing them.
DoleParams::DoleParams() {
    A = .00150;
    Alpha = 5;
    Beta = 0.5;
    Eccentricity = 0.15; //0.15
    Gamma = 1 / 3.0;
    K = 50;
    MassSol = 1;
    MassSun = MassSol * 1.0;
    MassNuclei = MassSun * 1e-15;
    B = MassSun * 1.2e-5;
    W = .2;
    KMigration = 0.0; // Default value: no migration
    MigrationFinalFraction = 1.0; // Default: no post-accretion migration (i.e., new_axis = old_axis * 1.0)
}

DoleParams::DoleParams( double pA, double pAlpha,double pBeta,double pEccentricity,double pGamma,double pK,double pMassSol,double pMassSun,double pMassNuclei,double pB,double pW, double pKMigration, double pMigrationFinalFraction)
    {
    A = pA ;
    Alpha =pAlpha;
    Beta = pBeta ;
    Eccentricity = pEccentricity;  
    Gamma = pGamma ;
    K = pK;
    MassSol = pMassSol;
    MassSun = pMassSun;
    MassNuclei = pMassNuclei;
    B = pB;
    W = pW;   
    KMigration = pKMigration; // Set migration coefficient
    MigrationFinalFraction = pMigrationFinalFraction; // Set QMigration
}


// Return disk density at orbital distance au
double DoleParams::density(double au) {
    return A * std::exp(-Alpha * std::pow(au, Gamma));
}

// Returns dust+gas density swept at orbital distance AU by a body of given
//  ratio (critical mass / mass).
double DoleParams::gasdensity(double au, double massratio) {
    return (K * density(au)) / (1 + std::sqrt(massratio) * (K-1));
}

// Return critical mass to form a gas giant at orbital distance au
double DoleParams::masscritical(double au) {
    return B * std::pow(au, -.75);
}

// I don't know what these functions really *represent*, but they're
//  repeatedly used.
double DoleParams::lowbound(double radius, double margin) {
    return radius - margin - W * (radius - margin) / (1 + W);
}

double DoleParams::highbound(double radius, double margin) {
    return radius + margin + W * (radius + margin) / (1 - W);
}

const double epsilon = .00001;
const int Maxint = 32767;

int imin(int a, int b = Maxint, int c = Maxint) {
    if (a < b)
    return (a < c) ? a : c;
    else
    return (b < c) ? b : c;
}

int imax(int a, int b = -Maxint, int c = -Maxint) {
    if (a > b)
    return (a > c) ? a : c;
    else
    return (b > c) ? b : c;
}

std::ofstream *Ps;
int PsPage = 1;
double
    PsBase = 50,           // base offset from lower left corner
    PsXscale = 1, PsYscale = 1,
    PsXoff = 0, PsYoff = 0;

void ps_end() {
    *Ps << "%%Trailer\nend" << std::flush;
    delete Ps;
    Ps = NULL;
}

void ps_beginpage(int page) {
    *Ps << "%%Page: " << page << ' ' << page << '\n'
    << PsXoff+PsBase << ' ' << PsYoff+PsBase << " translate\n"
    << PsXscale << ' ' << PsYscale << ' ' << " scale\n"
    << "/Helvetica findfont " << 9 / PsXscale << " scalefont setfont\n"
    << "0 setlinewidth\n";
}

void ps_begin(const char *file) {
    if (Ps != NULL) {
    ps_end();
    }

    Ps = new std::ofstream(file);
    *Ps << "%!PS-Adobe-2.1\n"
    << "%%Pages: 3\n"
    << "%%EndComments\n"
    << "/Helvetica findfont 12 scalefont setfont\n"
    << "0 setlinewidth\n"
    << "newpath\n"
    << "%%EndProlog\n";

    ps_beginpage(PsPage++);
}

void ps_showpage() {
    *Ps << "showpage\n";
    ps_beginpage(PsPage++);
}

void ps_window(double x1, double y1, double x2, double y2)
{
    const double width = 450;
    double
    xspan = x2 - x1, yspan = y2 - y1;

    PsXscale = width / xspan;
    PsYscale = PsXscale;          // could be width / yspan
    PsXoff   = -PsXscale * x1;
    PsYoff   = -PsYscale * y1;
}

void ps_circle(double x, double y, double radius, int fill)
{
    *Ps << x << ' ' << y << ' ' << radius << " 0 360 arc ";
    *Ps << (fill ? "fill" : "stroke") << std::endl;
}

void randomize(long seed)
{
    srand48(seed);
}

double rand01()
{
    return drand48();
}

double sqr(double x)
{
    return x * x;
}

// Return closest integer to arg
int cint(double arg)
{
    return static_cast<int>(std::floor(arg+0.5));
}

void ps_pset(double x, double y)
{
    ps_circle(x, y, 0.01, 0);
}

void ps_line(double x1, double y1, double x2, double y2)
{
    *Ps << x1 << ' ' << y1 << " moveto "
    << x2 << ' ' << y2 << " lineto stroke\n";
}

void ps_text(double x, double y, const char *s) {
    *Ps << x << ' ' << y << " moveto (" << s << ") show newpath\n";
}

// Draw scale on figure
void logscale(const char *xlabel = "", const char *ylabel = "") {
    ps_line(-1, -1,  3, -1);
    ps_line( 3, -1,  3,  1);
    ps_line( 3,  1,  3, -1);
    ps_line( 3, -1, -1, -1);

    ps_line(-1, 1, 3, 1);
    for (double au = 1; au <= 10 + epsilon; au += 1)  {
    ps_line(std::log10(au/10), 1, std::log10(au/10), .95);
    ps_line(std::log10(au), 1, std::log10(au), .95);
    ps_line(std::log10(au*10), 1, std::log10(au*10), .95);
    }

    ps_text(-1, 1, ".1");
    ps_text( 0, 1, "1");
    ps_text( 1, 1, "10");
    ps_text( 2, 1, "100");

    ps_text(2.3, 1, xlabel);
    ps_text(-1, .9, ylabel);
}

struct Nucleus {
    double
    axis,          // semimajor axis of the nuclei orbit
    eccen,         // eccentricity of the nuclei orbit
    mass,          // mass of the nuclei
    pRad,          // orbital radius at perigee
    aRad,          // orbital radius at apogee
    pAttr,         // grav attract dist at perigee
    aAttr;         // grav attract dist at apogee
    enum PlanetType { // Now explicitly named PlanetType
        Rock, Icy, GasGiant 
    }   type;          // type of planet

    Nucleus();
    void dump(int, std::ostream &, DoleParams *);
    friend int nucleusCompare(const void *, const void *);
    double lowbound(DoleParams *params);
    double highbound(DoleParams *params);

    // --- New helper functions for planet properties ---
    // Returns radius in kilometers
    double get_radius_km(); 
    double get_temperature_k();
};

Nucleus::Nucleus() {
    axis = eccen = mass = 0;
    pRad = aRad = 0;
    pAttr = aAttr = 0;
    type = Nucleus::Rock; // Default to Rock, will be determined later
}

double Nucleus::lowbound(DoleParams *params) {
    return params->lowbound(pRad, pAttr);
}

double Nucleus::highbound(DoleParams *params) {
    return params->highbound(aRad, aAttr);
}

// Comparison function used by qsort()
int nucleusCompare(const void *p1, const void *p2) {
    double   r1 = static_cast<const Nucleus *>(p1)->axis,
             r2 = static_cast<const Nucleus *>(p2)->axis;
    return (r1 < r2) ? -1 : (r1 > r2) ? 1 : 0;
}

// Dump nucleus stats to specified stream
void Nucleus::dump(int num, std::ostream &o, DoleParams *params) {
    double
    xplimit = pRad - pAttr,
    xalimit = aRad + aAttr,
    lowrange = lowbound(params),
    highrange = highbound(params),
    massCrit = params->masscritical(pRad);

    o << "Nucleus " << num << '\n';
    o << "\tRadius\t\t" << axis << '\n';
    o << "\tEccentricity\t" << eccen << '\n';
    o << "\tMass\t\t" << mass << '\n';
    o << "\tType\t\t";
    if (type == Nucleus::GasGiant) o << "Gas giant"; // Corrected: Nucleus::GasGiant
    else if (type == Nucleus::Icy) o << "Icy";      // Corrected: Nucleus::Icy
    else o << "Rocky";
    o << '\n';

    o << "\tRange          = [" << lowrange << "," << highrange << "]\n";
    o << "\tX[pa]limit     = [" << xplimit << "," << xalimit << "]\n";
    o << "\tR{peri,apo}gee = [" << pRad << "," << aRad << "]\n";

    o << "\tCritical Mass  = " << massCrit << ' '
      << '(' << mass / massCrit << "% of Critical)\n";
}

// --- Implementation of new helper functions for Nucleus ---
double Nucleus::get_radius_km() {
    double mass_kg = this->mass * SOLAR_MASS_TO_EARTH_MASS * EARTH_MASS_KG;
    double density_kg_m3;

    if (this->type == Nucleus::Rock) { // Corrected: Nucleus::Rock
        density_kg_m3 = EARTH_DENSITY_KG_M3;
    } else if (this->type == Nucleus::Icy) { // Corrected: Nucleus::Icy
        density_kg_m3 = ICE_DENSITY_KG_M3; 
    }
    else { // GasGiant
        density_kg_m3 = JUPITER_DENSITY_KG_M3; 
    }

    // Volume V = Mass / Density
    // Radius R = (3V / 4pi)^(1/3) in meters
    double radius_meters = std::cbrt((3 * mass_kg) / (4 * M_PI * density_kg_m3));
    return radius_meters / 1000.0; // Convert to kilometers
}

double Nucleus::get_temperature_k() {
    double albedo;
    if (this->type == Nucleus::Rock) { // Corrected: Nucleus::Rock
        albedo = ALBEDO_ROCKY_PLANET;
    } else if (this->type == Nucleus::Icy) { // Corrected: Nucleus::Icy
        albedo = ALBEDO_ICY_PLANET;
    }
    else { // GasGiant
        albedo = ALBEDO_GAS_GIANT;
    }

    // Standard black body equilibrium temperature (for fast rotating)
    // T_eq = T_sun * (R_sun / (2 * D))^(1/2) * (1 - A)^(1/4)
    // D is distance in meters, so axis_au * AU_TO_METERS
    double distance_meters = this->axis * AU_TO_METERS;
    if (distance_meters == 0) return 0; // Avoid division by zero for planets at axis 0

    double term_distance = SUN_RADIUS_METERS / (2 * distance_meters);
    if (term_distance < 0) term_distance = 0; // Ensure non-negative under sqrt
    
    double temp_eq = SUN_SURFACE_TEMP_K * std::sqrt(term_distance) * std::pow(1 - albedo, 0.25);

    return temp_eq;
}


int main(int ac, char *av[])
{
    // --- Help Option Check ---
    for (int i = 1; i < ac; i++) {
        if (strcmp(av[i], "-h") == 0 || strcmp(av[i], "--help") == 0) {
            std::cout << "Käyttö: " << av[0] << " [asetukset]\n";
            std::cout << "Simuloi planeettojen muodostumista ja migraatiota Dole-mallin pohjalta.\n\n";
            std::cout << "Asetukset:\n";
            std::cout << "  -h, --help           Näytä tämä ohjeviesti ja poistu.\n";
            std::cout << "  -dump                Tulosta yksityiskohtaisia tietoja konsoliin.\n";
            std::cout << "  -seed <luku>         Aseta satunnaislukugeneraattorin siemen (oletus: nykyinen aika).\n";
            std::cout << "  -maxplanets <luku>   Maksimi planeettojen määrä (oletus: 20).\n";
            std::cout << "  -aumin <AU>          Minimi etäisyys (AU) planeetan muodostumiselle (oletus: 0.3).\n";
            std::cout << "  -aumax <AU>          Maksimi etäisyys (AU) planeetan muodostumiselle (oletus: 50.0).\n";
            std::cout << "  -A <arvo>            Dole-mallin parametri A (oletus: 0.00150).\n";
            std::cout << "  -Alpha <arvo>        Dole-mallin parametri Alpha (oletus: 5.0).\n";
            std::cout << "  -Beta <arvo>         Dole-mallin parametri Beta (oletus: 0.5).\n";
            std::cout << "  -Eccentricity <arvo> Dole-mallin parametri Eccentricity (oletus: 0.15).\n";
            std::cout << "  -Gamma <arvo>        Dole-mallin parametri Gamma (oletus: 0.333333).\n";
            std::cout << "  -K <arvo>            Dole-mallin parametri K (oletus: 50.0).\n";
            std::cout << "  -MassSol <arvo>      Auringon massa (oletus: 1.0).\n";
            std::cout << "  -MassSun <arvo>      Auringon massa (oletus: 1.0).\n";
            std::cout << "  -MassNuclei <arvo>   Alkunukleuksen massa (oletus: 1e-15 * MassSun).\n";
            std::cout << "  -B <arvo>            Dole-mallin parametri B (oletus: 1.2e-5 * MassSun).\n";
            std::cout << "  -W <arvo>            Dole-mallin parametri W (oletus: 0.2).\n";
            std::cout << "  -KMigration <kerroin> Migraatiokerroin massan keräämisen aikana (oletus: 0.0, ei migraatiota).\n";
            std::cout << "                       Positiivinen arvo siirtää planeettaa sisäänpäin.\n";
            std::cout << "  -QMigration <kerroin> Migraatiokerroin massan keräämisen JÄLKEEN (oletus: 1.0, ei migraatiota).\n";
            std::cout << "                       Arvo 0.1 siirtää planeetan 10% alkuperäisestä etäisyydestään.\n";
            std::cout << "                       (Esim. 10 AU -> 1 AU, jos kerroin 0.1).\n";
            std::cout << "\nLumiraja asetettu arvoon: " << SNOW_LINE_AU << " AU.\n";
            return 0; // Poistu ohjelmasta ohjeen näyttämisen jälkeen
        }
    }
    // --- End Help Option Check ---

    int dumpflag = 0; // Default to no console dump unless -dump is specified

    long seed = time(0);
    int MaxPlanets = 20;    // Size of dynamic array containing generated planets

    double aumin = 0.3;  // Changed to double for precision
    double aumax = 50;   // Changed to double for precision

    // Default DoleParams values (as in the DoleParams constructor)
    double pA = .00150;
    double pAlpha = 5;
    double pBeta = 0.5;
    double pEccentricity = 0.15;
    double pGamma = 1 / 3.0;
    double pK = 50;
    double pMassSol = 1;
    double pMassSun = pMassSol * 1.0;
    double pMassNuclei = pMassSun * 1e-15;
    double pB = pMassSun * 1.2e-5;
    double pW = .2;
    double pKMigration = 0.0; // Default value: no migration
    double pMigrationFinalFraction = 1.0; // Default: no post-accretion migration

    // Command-line argument parsing
    for (int i = 1; i < ac; i++) { // Start from 1 to skip program name
        if (strcmp(av[i], "-dump") == 0) {
            dumpflag = 1;
        } else if (strcmp(av[i], "-seed") == 0) {
            if (i + 1 < ac) {
                seed = std::atol(av[++i]); // Use atol for long
            } else {
                std::cerr << "Error: -seed requires a value." << std::endl;
                return 1;
            }
        } else if (strcmp(av[i], "-maxplanets") == 0) {
            if (i + 1 < ac) {
                MaxPlanets = std::atoi(av[++i]);
            } else {
                std::cerr << "Error: -maxplanets requires a value." << std::endl;
                return 1;
            }
        } else if (strcmp(av[i], "-aumin") == 0) {
            if (i + 1 < ac) {
                aumin = std::atof(av[++i]);
            } else {
                std::cerr << "Error: -aumin requires a value." << std::endl;
                return 1;
            }
        } else if (strcmp(av[i], "-aumax") == 0) {
            if (i + 1 < ac) {
                aumax = std::atof(av[++i]);
            } else {
                std::cerr << "Error: -aumax requires a value." << std::endl;
                return 1;
            }
        }
        // DoleParams specific arguments
        else if (strcmp(av[i], "-A") == 0) {
            if (i + 1 < ac) pA = std::atof(av[++i]);
            else { std::cerr << "Error: -A requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-Alpha") == 0) {
            if (i + 1 < ac) pAlpha = std::atof(av[++i]);
            else { std::cerr << "Error: -Alpha requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-Beta") == 0) {
            if (i + 1 < ac) pBeta = std::atof(av[++i]);
            else { std::cerr << "Error: -Beta requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-Eccentricity") == 0) {
            if (i + 1 < ac) pEccentricity = std::atof(av[++i]);
            else { std::cerr << "Error: -Eccentricity requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-Gamma") == 0) {
            if (i + 1 < ac) pGamma = std::atof(av[++i]);
            else { std::cerr << "Error: -Gamma requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-K") == 0) {
            if (i + 1 < ac) pK = std::atof(av[++i]);
            else { std::cerr << "Error: -K requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-MassSol") == 0) {
            if (i + 1 < ac) pMassSol = std::atof(av[++i]);
            else { std::cerr << "Error: -MassSol requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-MassSun") == 0) {
            if (i + 1 < ac) pMassSun = std::atof(av[++i]);
            else { std::cerr << "Error: -MassSun requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-MassNuclei") == 0) {
            if (i + 1 < ac) pMassNuclei = std::atof(av[++i]);
            else { std::cerr << "Error: -MassNuclei requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-B") == 0) {
            if (i + 1 < ac) pB = std::atof(av[++i]);
            else { std::cerr << "Error: -B requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-W") == 0) {
            if (i + 1 < ac) pW = std::atof(av[++i]);
            else { std::cerr << "Error: -W requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-KMigration") == 0) { // Migration coefficient during accretion
            if (i + 1 < ac) pKMigration = std::atof(av[++i]);
            else { std::cerr << "Error: -KMigration requires a value." << std::endl; return 1; }
        }
        else if (strcmp(av[i], "-QMigration") == 0) { // Post-accretion migration coefficient
            if (i + 1 < ac) pMigrationFinalFraction = std::atof(av[++i]);
            else { std::cerr << "Error: -QMigration requires a value." << std::endl; return 1; }
        }
        else {
            std::cerr << "Tuntematon argumentti: " << av[i] << std::endl;
            std::cerr << "Käytä -h tai --help nähdäksesi ohjeet.\n";
            return 1; // Exit on unknown argument
        }
    }

    // Create DoleParams object with potentially overridden values
    DoleParams *params = new DoleParams(pA, pAlpha, pBeta, pEccentricity, pGamma, pK, pMassSol, pMassSun, pMassNuclei, pB, pW, pKMigration, pMigrationFinalFraction);

    // --- Dynamic aumin adjustment based on QMigration ---
    double original_aumin = aumin; // Store original aumin for comparison
    double effective_aumin_due_to_migration = original_aumin * params->MigrationFinalFraction;
    // If migration moves planets significantly inwards, adjust aumin lower for checks
    if (effective_aumin_due_to_migration < original_aumin) {
        aumin = effective_aumin_due_to_migration;
        // Ensure aumin doesn't become unrealistically small (e.g., minimum 0.01 AU)
        if (aumin < 0.01) aumin = 0.01;
        std::cout << "Huom: aumin säädetty migraation vuoksi arvoon: " << std::fixed << std::setprecision(4) << aumin << " AU\n";
    }
    // --- End Dynamic aumin adjustment ---

    double   nucleieccent, nucleiradius, masslast, eccent, mass, rperigee,
             rapogee, radius, mu, xperigee, xapogee, masscritical,
             density, bw, volume, da, dp, masspass;

    int      lowband, highband, dustcheck, iterate, hit;

    // Bands are measured in 1/5 au intervals from 0 to 50 AU
    const int
    MaxBand = 50 * 5,
    MaxNuclei = 1000;    // Max # of nuclei to inject


    enum { Mixed, Gas, Empty }
        band[MaxBand+1];    // Contents of band


    Nucleus *nuclei = new Nucleus[MaxPlanets];
    
    int nplanet = 0; 

    randomize(seed);

    ps_window(-1, -1, 2, 1);
    ps_begin("planets.ps");

    // Set up primordial cloud
    // Dust within .3 au is swept out by radiation pressure

    for (int n_band = 1; n_band <= MaxBand; n_band++) 
    band[n_band] = Mixed;

    // Draw scale on figure
    logscale("AU");

    // Plot density function; space samples even on log scale
    double max = 1.1 * params->density(0.3), lastau, lastrho,
             logaumin = std::log10(aumin),
             logaumax = std::log10(aumax),
             logaustep = (logaumax - logaumin) / 100;

    for (double logau = logaumin; logau <= logaumax; logau += logaustep) {
    double rho = params->density(std::pow(10,logau)) / max;

    if (logau > logaumin)
        ps_line(lastau, lastrho, logau, rho);
    else {
        char buf[100];
        std::sprintf(buf, "density [max = %g]", max);
        ps_text(logau, rho, buf);
    }

    lastau = logau;
    lastrho = rho;
    }

    double yBand = 0.9;
    // Inject nuclei into the cloud and acrete gas and dust
    int nuclei_count = 0; // Declare a new variable to keep track of nuclei processed
    for (nuclei_count = 0; nuclei_count < MaxNuclei; nuclei_count++) {
    // Check to see if all bands are taken
    int bandfull = 0;
    for (int i = 1; i <= MaxBand; i++) {
        if (band[i] == Mixed)  {
        bandfull = 1;
        break;
        }
    }

    if (bandfull == 0)
        break;

    nucleieccent = 1 - std::pow(1 - rand01(), .077);
    while ((nucleiradius = rand01() * 50) < .3)
        ;

    Nucleus body;
    body.axis = nucleiradius;
    body.eccen = nucleieccent;
    body.mass = params->MassNuclei;
    // Set initial type, will be determined accurately later
    body.type = Nucleus::Rock; 

    // This is only used on first pass of nuclei
    masslast = 0;
    iterate = 0;

    // --- Start of Mass accumulation loop ---
    int current_mass_iteration = 0;
    const int MAX_MASS_ACCUMULATION_ITERATIONS = 5000; 

    while (1) {
        current_mass_iteration++;
        if (current_mass_iteration > MAX_MASS_ACCUMULATION_ITERATIONS) {
            std::cerr << "Warning: Mass accumulation loop for nucleus " << nuclei_count
                      << " exceeded " << MAX_MASS_ACCUMULATION_ITERATIONS << " iterations. Breaking to prevent freeze." << std::endl;
            break; // Force break to prevent infinite loop
        }

        radius = body.axis;
        eccent = body.eccen;
        mass = body.mass;
        body.pRad = rperigee = nucleiradius * (1 - nucleieccent);
        body.aRad = rapogee = nucleiradius * (1 + nucleieccent);
        if (mass == 0)
        mass = 1e-15;
        mu = std::pow(mass / (1 + mass), .25);
        body.pAttr = xperigee = rperigee * mu;
        body.aAttr = xapogee = rapogee * mu;

        // Nuclei sweeps up band
        // Establish bounds on swept volume
        lowband = cint(body.lowbound(params) * 5);
        if (lowband < 1)
        lowband = 1;

        highband = cint(body.highbound(params) * 5);
        if (highband < 1)
        highband = 1;
        else if (highband > MaxBand)
        highband = MaxBand;

        if (lowband == highband)
        highband++;

        // calculate critical mass limits
        masscritical = params->masscritical(rperigee);

        /* check for bands with dust within range */
        if (iterate == 0) {
        dustcheck = 0;

        for (int bandno = lowband; bandno <= highband; bandno++) {
            if (masslast > 0 || band[bandno] == Mixed) {
            dustcheck = 1;
            break;
            }
        }

        // If no bands have dust in them (band[bandno] == Mixed), then dud
        // dustcheck == 1 means dust is in some band

        if (dustcheck == 0) {
            if (masslast == 0) {
            body.axis = 0;
            body.eccen = 0;
            body.mass = 0;
            }
            break; // exit mass accumulation loop
        }
        }

        // Calculate mass using Dole donut approximation
        if (mass == 0)
        mass = 1e-15;

        if (mass < masscritical)
        density = params->density(body.axis);
        else {
        // Sweeps dust and gas
        density = params->gasdensity(body.axis, masscritical / mass);
        }

        bw = 2 * body.axis * body.eccen +
             xapogee + xperigee +
             params->W * (rapogee + xapogee) / (1 - params->W) +
             params->W * (rperigee - xperigee) / (1 + params->W);

        volume = 2 * M_PI * body.axis * bw * (xapogee + xperigee);

        double sweepmass = volume * density;    // unused

        dp = body.lowbound(params);
        da = body.highbound(params);

        lowband = cint(dp * 5);
        highband = cint(da * 5);

        if (lowband == highband)
        highband++;
        if (highband > MaxBand)
        highband = MaxBand;
        masspass = 0;

        for (int bandno = lowband; bandno <= highband; bandno++) {
        double
            au = bandno / 5.0,
            xpnow, xanow,
            bandvol, bandvol2;

        // Calculate mass of wedge of doughnut
        xpnow = xperigee + (xapogee - xperigee) *
                  (bandno - lowband) / static_cast<double>(highband - lowband);
        xanow = xperigee + (xapogee - xperigee) *
                  (bandno + 1 - lowband) / static_cast<double>(highband - lowband);

        bandvol = 2 * M_PI * au * 0.2 * (xpnow + xanow);

        for (int i = 0; i < nplanet; i++) { 
            double
            dp2 = nuclei[i].lowbound(params),
            da2 = nuclei[i].highbound(params);

            if (da2 < au || dp2 > au + 0.2)
            continue;

            // Overlap exists, find bounds on overlap volume
            double
            bw2 = 2 * nuclei[i].axis * nuclei[i].eccen +
              nuclei[i].pAttr +
              nuclei[i].aAttr +
              params->W * (nuclei[i].aRad + nuclei[i].aAttr) / (1 - params->W) +
              params->W * (nuclei[i].pRad - nuclei[i].pAttr) / (1 + params->W);

            // At au, overlap has xp2now and xa2now as heights
            double xp2now = nuclei[i].pAttr +
                  (nuclei[i].aAttr - nuclei[i].pAttr) *
                    (au - dp2) / (da2 - dp2);
            double xa2now = nuclei[i].pAttr +
                  (nuclei[i].aAttr - nuclei[i].pAttr) *
                    (au + 0.2 - dp2) / (da2 - dp2);
            bandvol2 = 2 * M_PI * au * 0.2 * (xp2now + xa2now);

            // If previously swept band larger than this swept, no dust
            //    swept up.

            if (bandvol2 >= bandvol)
            bandvol = 0;
            else
            bandvol -= bandvol2;
        }

        masspass += bandvol * density;
        }

        body.mass = mass = masspass;

        // KMigration logic (during accretion)
        double delta_axis = 0.0;
        if (params->KMigration != 0.0) {
            delta_axis = -params->KMigration * body.axis; 
        }

        // Apply KMigration to the semi-major axis
        body.axis += delta_axis;

        // Ensure the planet stays within the defined AU bounds
        if (body.axis < aumin) body.axis = aumin;
        if (body.axis > aumax) body.axis = aumax;

        // Update other orbital parameters based on the new axis
        body.pRad = body.axis * (1 - body.eccen);
        body.aRad = body.axis * (1 + body.eccen);
        mu = std::pow(body.mass / (1 + body.mass), .25);
        body.pAttr = body.pRad * mu;
        body.aAttr = body.aRad * mu;

        // check for mass growth convergence
        if (mass == 0)
        mass = 1e-15;
        
        // **Revised convergence condition:**
        // Break if mass has converged. KMigration will occur for the next overall nucleus iteration.
        if (std::fabs(masslast / mass - 1) < 0.01) {
            break; 
        }
        
        masslast = mass;
        iterate = 1; 
    }    // --- End of mass accumulation loop ---

    // Set planet type based on final mass and position
    if (mass >= masscritical) {
        body.type = Nucleus::GasGiant;
    } else {
        // If not a gas giant, check against snow line for Rocky or Icy
        if (body.axis >= SNOW_LINE_AU) {
            body.type = Nucleus::Icy;
        } else {
            body.type = Nucleus::Rock;
        }
    }


    // Clear out bands emptied of dust
    if (dustcheck == 1) {
        std::cout << "event " << nuclei_count << " completed mass growth; swept from "
                  << dp << " to " << da;
        if (body.type == Nucleus::GasGiant)
        std::cout << "(gas giant)";
        else if (body.type == Nucleus::Icy)
        std::cout << "(icy)";
        else
        std::cout << "(rocky)";
        std::cout << std::endl;
    }
    if (lowband == highband)
        highband++;

    for (int bandno = lowband; bandno <= highband; bandno++) {
        if (mass < masscritical) { // For small planets, clear as Gas or Mixed to Empty based on type?
            // This logic is a bit ambiguous in the original code. 
            // For simplicity, let's assume if it's a small planet,
            // it still clears the band as "Gas" (if mixed) or "Empty" if it consumes everything.
            // The original logic here seemed to set band to Gas if Mixed for non-giants.
            if (band[bandno] == Mixed)
                band[bandno] = Gas; 
        }

        if (mass >= masscritical) {
        // Clear out bands emptied by gas giant formation
        if (band[bandno] == Gas || band[bandno] == Mixed)
            band[bandno] = Empty;
        }
    }

    // check for orbital and  gravitational overlap and coalescence
    if (body.axis == 0)
        continue;

    hit = 0;
    for (int i = 0; i < nplanet; i++) { 
        double  newradius, newmass, neweccent, term1, term2, term3, munew;

        if ((nuclei[i].aRad < (rperigee - xperigee) &&
           (nuclei[i].aRad + nuclei[i].aAttr) < rperigee) ||
           ((rapogee + xapogee) < nuclei[i].pRad &&
           rapogee < (nuclei[i].pRad - nuclei[i].pAttr)))
        continue;

        std::cout << "coalescence of nuclei  " << nuclei_count << ' ' << i << std::endl;

        hit = 1;
        newradius = (body.mass + nuclei[i].mass) /
            (body.mass / body.axis + nuclei[i].mass / nuclei[i].axis);

        term1 = body.mass * std::sqrt(body.axis) * std::sqrt(1 - sqr(body.eccen));
        term2 = nuclei[i].mass * std::sqrt(nuclei[i].axis) *
                      std::sqrt(1 - sqr(nuclei[i].eccen));
        term3 = (body.mass + nuclei[i].mass) * std::sqrt(newradius);

        neweccent = sqr(std::abs(1 - std::pow((term1 + term2) / term3, 2)));
        newmass = body.mass + nuclei[i].mass;

        // When merging, the new planet's type should ideally reflect its final state.
        // For simplicity, if either was a gas giant, the result is a gas giant.
        // Otherwise, re-evaluate based on new mass and position.
        Nucleus::PlanetType merged_type; 
        if (body.type == Nucleus::GasGiant || nuclei[i].type == Nucleus::GasGiant) {
            merged_type = Nucleus::GasGiant;
        } else {
            // Re-evaluate based on new_axis
            if (newradius >= SNOW_LINE_AU) {
                merged_type = Nucleus::Icy;
            } else {
                merged_type = Nucleus::Rock;
            }
        }


        nuclei[i].axis = 0;
        nuclei[i].eccen = 0;
        nuclei[i].mass = 0;
        body.axis = newradius;
        body.eccen = neweccent;
        body.mass = newmass;
        body.type = merged_type; // Set merged type
        body.pRad = newradius * (1 - neweccent);
        body.aRad = newradius * (1 + neweccent);
        munew = std::pow(newmass / (1 + newmass), .25);
        body.pAttr = body.pRad * munew;
        body.aAttr = body.aRad * munew;
    }

    mass = body.mass;

    // Show depletion of bands
    int
        lowband2 = cint(body.lowbound(params) * 5),
        highband2 = cint(body.highbound(params) * 5);

    lowband2 = imax(lowband2, 6, lowband);
    highband2 = imax(imin(highband2, MaxBand, highband), 6);
    if (lowband2 == highband2)
        highband2++;

    ps_line(std::log10(lowband2 * 0.2), yBand, std::log10(highband2 * 0.2), yBand);
    yBand -= 0.01;

    // Add new planet
    if (nplanet >= MaxPlanets) { 
        Nucleus *newplanet = new Nucleus[MaxPlanets*2];

        // Copy old planets to new
        for (int i = 0; i < nplanet; i++) 
        newplanet[i] = nuclei[i];

        delete [] nuclei;
        nuclei = newplanet;
        MaxPlanets *= 2;
    }
    nuclei[nplanet++] = body; 

    // Sweep planet array, removing merged bodies
    for (int i = 0; i < nplanet; i++) { 
        if (nuclei[i].axis == 0) {
        std::cerr << "Nuking body " << i << "; " << nplanet << " remaining" << std::endl;
          nuclei[i] = nuclei[nplanet-1];
          nuclei[nplanet-1].axis = 0;
          nplanet--;
        }
    }
    } // End of the main loop that iterates MaxNuclei times

    std::cout << " all bands taken....." << nuclei_count << "     nuclei used       " << std::endl;
    ps_text(2.2, 0.9, "emptied bands");

    // --- Apply QMigration AFTER all planets have been formed and merged ---
    if (params->MigrationFinalFraction != 1.0) {
        std::cout << "Applying post-accretion migration (QMigration)...\n";
        for (int i = 0; i < nplanet; ++i) {
            if (nuclei[i].axis > 0) { // Only migrate valid planets
                double old_axis = nuclei[i].axis;
                nuclei[i].axis *= params->MigrationFinalFraction;
                
                // Ensure it stays within bounds after migration
                if (nuclei[i].axis < aumin) nuclei[i].axis = aumin;
                if (nuclei[i].axis > aumax) nuclei[i].axis = aumax;

                // Update orbital parameters as eccentricity remains the same but axis changes
                nuclei[i].pRad = nuclei[i].axis * (1 - nuclei[i].eccen);
                nuclei[i].aRad = nuclei[i].axis * (1 + nuclei[i].eccen);
                double munew = std::pow(nuclei[i].mass / (1 + nuclei[i].mass), .25);
                nuclei[i].pAttr = nuclei[i].pRad * munew;
                nuclei[i].aAttr = nuclei[i].aRad * munew;

                // Re-evaluate type after QMigration if it's not a gas giant
                if (nuclei[i].type != Nucleus::GasGiant) { 
                     if (nuclei[i].axis >= SNOW_LINE_AU) {
                         nuclei[i].type = Nucleus::Icy; 
                     } else {
                         nuclei[i].type = Nucleus::Rock; 
                     }
                }
                // For debugging:
                // std::cout << "Planet " << i << ": Old Axis " << old_axis << " AU -> New Axis " << nuclei[i].axis << " AU\n";
            }
        }
    }

    // Sort nuclei by radius (again, as QMigration might change order)
    std::qsort(static_cast<void *>(nuclei), nplanet, sizeof(Nucleus), nucleusCompare); 

    // --- Hill-kriteerin mukainen suodatus ---
    std::cout << "Suoritetaan Hill-kriteerin mukainen suodatus...\n";
    bool unstable_orbits_found = true;
    while (unstable_orbits_found) {
        unstable_orbits_found = false;
        // Re-sort after any potential previous removal to ensure adjacent planets are checked correctly
        std::qsort(static_cast<void *>(nuclei), nplanet, sizeof(Nucleus), nucleusCompare);

        for (int i = 0; i < nplanet; ++i) {
            if (nuclei[i].axis == 0 || nuclei[i].mass == 0) continue; // Skip removed planets

            for (int j = i + 1; j < nplanet; ++j) {
                if (nuclei[j].axis == 0 || nuclei[j].mass == 0) continue; // Skip removed planets

                // Calculate Hill radii for both planets
                // Mass needs to be in terms of Solar Masses for consistency with Dole's units
                // (or ensure M_sun is in same units as planet mass)
                double M_sun_for_hill = params->MassSun; // Typically 1.0 Solar Mass (Dole's units)

                // Convert planet masses to Solar Mass for calculation
                double m_i_sol = nuclei[i].mass; 
                double m_j_sol = nuclei[j].mass;

                // Hill radius calculation: a * cbrt(m / (3 * M_sun))
                double rh_i = nuclei[i].axis * std::cbrt(m_i_sol / (3.0 * M_sun_for_hill));
                double rh_j = nuclei[j].axis * std::cbrt(m_j_sol / (3.0 * M_sun_for_hill));

                // Distance between planets
                double distance_between_planets = std::fabs(nuclei[i].axis - nuclei[j].axis);

                // Check for overlap. A typical stability criterion is 2*R_H.
                // For removal, we can use a stricter criterion, e.g., if distance < sum of Hill radii.
                // Or for robust stability, distance > some multiple of (RH_i + RH_j)
                // Using (rh_i + rh_j) as the threshold for removal.
                if (distance_between_planets < (rh_i + rh_j)) {
                    std::cout << "Instability found between Planet " << i+1 << " (" << std::fixed << std::setprecision(2) << nuclei[i].axis << " AU) and Planet " << j+1 << " (" << nuclei[j].axis << " AU).\n";
                    std::cout << "  Hill radii: " << rh_i << " AU (P" << i+1 << "), " << rh_j << " AU (P" << j+1 << "). Sum: " << (rh_i + rh_j) << " AU.\n";
                    std::cout << "  Distance: " << distance_between_planets << " AU.\n";

                    // Remove the less massive planet
                    if (nuclei[i].mass <= nuclei[j].mass) {
                        std::cout << "  Removing Planet " << i+1 << " (mass: " << nuclei[i].mass * SOLAR_MASS_TO_EARTH_MASS << " Earth Masses).\n";
                        nuclei[i].axis = 0; // Mark for removal
                        nuclei[i].mass = 0;
                    } else {
                        std::cout << "  Removing Planet " << j+1 << " (mass: " << nuclei[j].mass * SOLAR_MASS_TO_EARTH_MASS << " Earth Masses).\n";
                        nuclei[j].axis = 0; // Mark for removal
                        nuclei[j].mass = 0;
                    }
                    unstable_orbits_found = true; // Flag that changes occurred, need another pass
                    break; // Break inner loop, as the array structure changed
                }
            }
            if (unstable_orbits_found) break; // If a planet was removed, break outer loop to re-clean and restart check
        }

        // After a pass, re-clean the nuclei array if any were removed
        if (unstable_orbits_found) {
            int current_nplanet = 0;
            for (int k = 0; k < nplanet; ++k) {
                if (nuclei[k].axis != 0 && nuclei[k].mass != 0) {
                    nuclei[current_nplanet++] = nuclei[k];
                }
            }
            nplanet = current_nplanet; // Update the count of valid planets
            // Note: Array is re-sorted at the beginning of the while loop before checks.
        }
    }
    std::cout << "Hill-kriteerin suodatus valmis. Planeettoja jäljellä: " << nplanet << "\n";
    // --- Hill-kriteerin mukainen suodatus loppu ---


    // Draw planets, gas giants as filled circles
    double massSum = 0;
    for (int i = 0; i < nplanet; i++) { 
    massSum += nuclei[i].mass;

    double au = std::log10(nuclei[i].axis);
    // Use a dummy radius for plotting on log-scale for visibility
    // This radius is not the physical radius, just for visualization on the PS plot.
    double r_plot = std::cbrt(nuclei[i].mass * 1e5); // Scale for better visibility if needed

    // Fill circle for gas giants, outline for rocky, different fill for icy?
    // For PS, all non-GasGiant (Rocky, Icy) will be outlined for simplicity,
    // as PostScript doesn't easily support multiple fill patterns/colors without more complex code.
    bool fill_plot = (nuclei[i].type == Nucleus::GasGiant); 
    ps_circle(au, 0, r_plot, fill_plot);
    }
    ps_showpage();

    std::cout << "Total mass = " << massSum << std::endl;

    if (dumpflag) {
    std::cout << "Random number seed =" << seed << std::endl;

    for (int i = 0; i < nplanet; i++) { 
        if (nuclei[i].axis > 0)
        nuclei[i].dump(i, std::cout, params);
    }

    std::cout << "Bands with dust still in them:\n";
    for (int i = 1; i <= MaxBand; i++)
        if (band[i] == Mixed)
        std::cout << i << ' ';
    std::cout << std::endl;
    }

    // --- CSV Output Section ---
    std::ofstream csv_file("planets_output.csv");
    if (csv_file.is_open()) {
        // Updated header for planet_radius_km
        csv_file << "id,distance_au,ecc,mass_mearths,type,planet_radius_km,planet_temperature_k\n";
        csv_file << std::fixed << std::setprecision(8); // Set precision for output

        for (int i = 0; i < nplanet; i++) {
            if (nuclei[i].axis > 0 && nuclei[i].mass > 0) { // Only valid planets
                csv_file << i + 1 << ","; // Planet ID (1-based)
                csv_file << nuclei[i].axis << ","; // distance_au
                csv_file << nuclei[i].eccen << ","; // ecc
                csv_file << nuclei[i].mass * SOLAR_MASS_TO_EARTH_MASS << ","; // mass_mearths

                std::string planet_type_str;
                if (nuclei[i].type == Nucleus::GasGiant) { 
                    planet_type_str = "GasGiant";
                } else if (nuclei[i].type == Nucleus::Icy) { 
                    planet_type_str = "Icy";
                } else {
                    planet_type_str = "Rocky";
                }
                csv_file << planet_type_str << ","; // type

                // Use the new get_radius_km function
                csv_file << nuclei[i].get_radius_km() << ","; // planet_radius in KM

                double temperature_k = nuclei[i].get_temperature_k();
                csv_file << temperature_k << "\n"; // planet_temperature_k
            }
        }
        csv_file.close();
        std::cout << "Planet data saved to planets_output.csv\n";
    } else {
        std::cerr << "Error: Unable to open planets_output.csv for writing.\n";
    }
    // --- End CSV Output Section ---


    ps_end();
    delete params; // Clean up dynamically allocated DoleParams object
    delete[] nuclei; // Clean up dynamically allocated nuclei array

    return (0);
}
