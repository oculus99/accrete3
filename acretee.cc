
// v 0002 own density distribution

#include <iostream>
#include <fstream>    // For std::ofstream
#include <cstring>    // For strcmp (C-style string functions)
#include <cmath>      // For math functions like pow, log10, floor, sqrt, fabs, cbrt, M_PI
#include <cstdlib>    // For atoi, qsort, srand48, drand48, lrand48, time
#include <cstdio>     // For sprintf (needed if using sprintf, though std::to_string or string streams are preferred in C++)
#include <iomanip>    // For std::fixed and std::setprecision
#include <string>     // For std::string
#include <vector>     // For std::vector to store planets
#include <algorithm> // For std::sort

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
const double EARTH_RADIUS_METERS = 6.371e6;         // Earth Radius in meters
const double EARTH_MASS_KG = 5.972e24;              // Earth Mass in kilograms
const double JUPITER_DENSITY_KG_M3 = 1330.0;       // Jupiter density kg/m^3
const double EARTH_DENSITY_KG_M3 = 5510.0;          // Earth density kg/m^3 (approximate for rocky planets)
const double ICE_DENSITY_KG_M3 = 917.0;          // Density of ice (water ice) in kg/m^3
const double IRON_DENSITY_KG_M3 = 7874.0;          // Density of iron in kg/m^3 (approximate for iron-rich core)

const double AU_TO_METERS = 1.496e11;              // 1 AU in meters
const double SUN_SURFACE_TEMP_K = 5778.0;          // Sun's effective surface temperature in Kelvin
const double SUN_RADIUS_METERS = 6.957e8;          // Sun's radius in meters

// Albedo estimates
const double ALBEDO_ROCKY_PLANET = 0.3;
const double ALBEDO_GAS_GIANT = 0.5;
const double ALBEDO_ICY_PLANET = 0.7; // Higher albedo for icy worlds

// New: Snow Line definition
const double SNOW_LINE_AU = 5.0; // Lumiraja AU:na (2.5 * ekosfäärin etäisyys, joka Auringolla 1 AU)

// --- Constants for material composition ---
// Proportions for rocky/iron dust inside snow line
const double P_IRON_ROCKY = 0.3;
const double P_ROCK_ROCKY = 0.7;

// Proportions for icy/rocky dust outside snow line
const double P_ICE_ICY = 0.75;
const double P_ROCK_ICY = 0.25;

// Snow line factor: (p_ir/ Assuming SOLAR_MASS_TO_EARTH_MASS is defined elsewhere, e.g.:
// const double SOLAR_MASS_TO_EARTH_MASS = 333000.0;



// Based on your default proportions: (0.3+0.7) / 0.75 = 1 / 0.75 = 1.333...
// If you intended this to be exactly 4.0, please clarify.
const double SNOW_LINE_DUST_FACTOR = (P_IRON_ROCKY + P_ROCK_ROCKY) / P_ICE_ICY;

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
    KMigration,         // Migration coefficient during accretion
    MigrationFinalFraction; // New: Final fraction for post-accretion migration (QMigration)

    DoleParams();
    DoleParams( double pA, double pAlpha,double pBeta,double pEccentricity,double pGamma,double pK,double pMassSol,double pMassSun,double pMassNuclei,double pB,double pW, double pKMigration, double pMigrationFinalFraction);
    
    // Original density function, now uses total_dust_density
    double density(double au);
    double gasdensity(double au, double massratio);
    double masscritical(double au);
    double lowbound(double radius, double margin);
    double highbound(double radius, double margin);

    // New: Functions to get dust density for specific components
    double base_dust_density(double au); // Private helper
    double dust_density_iron(double au);
    double dust_density_rock(double au);
    double dust_density_ice(double au);
    double total_dust_density(double au); // Total dust density including snow line factor
};

// Initialize to defaults. See Sagan's article for insight into changing them.
DoleParams::DoleParams() {
    A = .00150;  // orig  .00150
    Alpha = 5;
    Gamma = 1 / 3.0;
    Beta = 0.5;
    Eccentricity = 0.15; //0.15
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

// Private helper to get the base density before composition or snow line factor
double DoleParams::base_dust_density(double au) {
  //  return A * std::exp(-Alpha * std::pow(au, Gamma));
    
    return A* std::pow(au/5,-0.333)*0.0025 ;
}

// Return dust density for iron component
double DoleParams::dust_density_iron(double au) {
    if (au < SNOW_LINE_AU) {
        return base_dust_density(au) * P_IRON_ROCKY;
    }
    return 0.0;
}

// Return dust density for rock component
double DoleParams::dust_density_rock(double au) {
    if (au < SNOW_LINE_AU) {
        return base_dust_density(au) * P_ROCK_ROCKY;
    } else {
        return base_dust_density(au) * P_ROCK_ICY;
    }
}

// Return dust density for ice component
double DoleParams::dust_density_ice(double au) {
    if (au >= SNOW_LINE_AU) {
        return base_dust_density(au) * P_ICE_ICY;
    }
    return 0.0;
}

// Return total dust density at orbital distance au, applying snow line factor
double DoleParams::total_dust_density(double au) {
    double base_rho = base_dust_density(au);
    if (au >= SNOW_LINE_AU) {
        return base_rho * SNOW_LINE_DUST_FACTOR;
    } else {
        return base_rho;
    }
}

// Return disk density at orbital distance au (now represents total accretable solid material)
double DoleParams::density(double au) {
    return total_dust_density(au);
}

// Returns dust+gas density swept at orbital distance AU by a body of given
// ratio (critical mass / mass). Now uses total_dust_density.
double DoleParams::gasdensity(double au, double massratio) {
    return (K * total_dust_density(au)) / (1 + std::sqrt(massratio) * (K-1));
}

// Return critical mass to form a gas giant at orbital distance au
double DoleParams::masscritical(double au) {
    return B * std::pow(au, -.75);
}

// I don't know what these functions really *represent*, but they're
// repeatedly used.
double DoleParams::lowbound(double radius, double margin) {
    return radius - margin - W * (radius - margin) / (1 + W);
}

double DoleParams::highbound(double radius, double margin) {
    return radius + margin + W * (radius + margin) / (1 - W);
}



// Assuming SOLAR_MASS_TO_EARTH_MASS is defined elsewhere, e.g.:
// const double SOLAR_MASS_TO_EARTH_MASS = 333000.0;



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
    PsBase = 50,          // base offset from lower left corner
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
    PsYscale = PsXscale;           // could be width / yspan
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
    axis,            // semimajor axis of the nuclei orbit
    eccen,            // eccentricity of the nuclei orbit
    mass,            // mass of the nuclei
    pRad,            // orbital radius at perigee
    aRad,            // orbital radius at apogee
    pAttr,            // grav attract dist at perigee
    aAttr;            // grav attract dist at apogee
    
    // Updated enum with new planet types
    enum PlanetType {
        iron_p,      // Iron-rich planet
        rock_p,      // Rocky planet
        ice_p,        // Icy planet
        gas_p         // Gas giant
    }    type;              // type of planet

    // New: Track accumulated material for more accurate typing
    double accumulated_iron_mass; // in solar masses
    double accumulated_rock_mass; // in solar masses
    double accumulated_ice_mass;  // in solar masses


    Nucleus();
    void dump(int, std::ostream &, DoleParams *);
    friend int nucleusCompare(const void *, const void *);
    double lowbound(DoleParams *params);
    double highbound(DoleParams *params);

    // --- New helper functions for planet properties ---
    // Returns radius in kilometers
    double get_radius_km();
    double get_temperature_k();
    void determine_planet_type(DoleParams* params); // New function to set planet type
        bool operator<(const Nucleus& other) const {
        return axis < other.axis;
    }
    
};

void export_planets_to_csv(const std::vector<Nucleus>& planets, const std::string& filename);

Nucleus::Nucleus() {
    axis = eccen = mass = 0;
    pRad = aRad = 0;
    pAttr = aAttr = 0;
    type = Nucleus::rock_p; // Default to rock_p
    accumulated_iron_mass = 0.0;
    accumulated_rock_mass = 0.0;
    accumulated_ice_mass = 0.0;
}

double Nucleus::lowbound(DoleParams *params) {
    return params->lowbound(pRad, pAttr);
}

double Nucleus::highbound(DoleParams *params) {
    return params->highbound(aRad, aAttr);
}

// Comparison function used by qsort()
int nucleusCompare(const void *p1, const void *p2) {
    double    r1 = static_cast<const Nucleus *>(p1)->axis,
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
    o << "\tRadius\t\t" << std::fixed << std::setprecision(4) << axis << " AU\n";
    o << "\tEccentricity\t" << std::fixed << std::setprecision(4) << eccen << '\n';
    
    // Convert mass from Solar Masses to Earth Masses
    o << "\tMass\t\t" << std::fixed << std::setprecision(4) << (mass * SOLAR_MASS_TO_EARTH_MASS) << " Earth Masses\n";
    
    o << "\tType\t\t";
    if (type == Nucleus::gas_p) o << "Gas giant";
    else if (type == Nucleus::ice_p) o << "Icy";
    else if (type == Nucleus::iron_p) o << "Iron";
    else o << "Rocky";
    o << '\n';

    o << "\tRange         = [" << std::fixed << std::setprecision(4) << lowrange << "," << std::fixed << std::setprecision(4) << highrange << "]\n";
    o << "\tX[pa]limit    = [" << std::fixed << std::setprecision(4) << xplimit << "," << std::fixed << std::setprecision(4) << xalimit << "]\n";
    o << "\tR{peri,apo}gee = [" << std::fixed << std::setprecision(4) << pRad << "," << std::fixed << std::setprecision(4) << aRad << "]\n";

    o << "\tCritical Mass  = " << std::fixed << std::setprecision(4) << (massCrit * SOLAR_MASS_TO_EARTH_MASS) << " Earth Masses "
      << '(' << std::fixed << std::setprecision(2) << (mass / massCrit * 100) << "% of Critical)\n";

    // Calculate and print material mass fractions
    double total_solid_mass = accumulated_iron_mass + accumulated_rock_mass + accumulated_ice_mass;
    double gas_mass = mass - total_solid_mass; // Assuming any mass beyond solids is gas
    if (gas_mass < 0) gas_mass = 0; // Ensure gas mass is not negative due to precision

    double total_mass_for_fractions = total_solid_mass + gas_mass; // This should be equal to 'mass'

    o << "\tMaterial Mass Fractions:\n";
    if (total_mass_for_fractions > 0) {
        o << "\t\tIron:\t" << std::fixed << std::setprecision(4) << (accumulated_iron_mass / total_mass_for_fractions) << '\n';
        o << "\t\tRock:\t" << std::fixed << std::setprecision(4) << (accumulated_rock_mass / total_mass_for_fractions) << '\n';
        o << "\t\tIce:\t" << std::fixed << std::setprecision(4) << (accumulated_ice_mass / total_mass_for_fractions) << '\n';
        o << "\t\tGas:\t" << std::fixed << std::setprecision(4) << (gas_mass / total_mass_for_fractions) << '\n';
    } else {
        o << "\t\tNo accumulated mass to calculate fractions.\n";
    }
}


// --- Implementation of new helper functions for Nucleus ---
double Nucleus::get_radius_km() {
    double mass_kg = this->mass * SOLAR_MASS_TO_EARTH_MASS * EARTH_MASS_KG;
    double density_kg_m3;

    if (this->type == Nucleus::rock_p) {
        density_kg_m3 = EARTH_DENSITY_KG_M3;
    } else if (this->type == Nucleus::ice_p) {
        density_kg_m3 = ICE_DENSITY_KG_M3;
    } else if (this->type == Nucleus::iron_p) {
        density_kg_m3 = IRON_DENSITY_KG_M3;
    }
    else { // gas_p
        density_kg_m3 = JUPITER_DENSITY_KG_M3;
    }

    // Volume V = Mass / Density
    // Radius R = (3V / 4pi)^(1/3) in meters
    // Handle potential division by zero or negative density if mass is 0 or density is not correctly set.
    if (density_kg_m3 <= 0 || mass_kg <= 0) {
        return 0.0;
    }
    double radius_meters = std::cbrt((3 * mass_kg) / (4 * M_PI * density_kg_m3));
    return radius_meters / 1000.0; // Convert to kilometers
}

double Nucleus::get_temperature_k() {
    double albedo;
    if (this->type == Nucleus::rock_p || this->type == Nucleus::iron_p) {
        albedo = ALBEDO_ROCKY_PLANET;
    } else if (this->type == Nucleus::ice_p) {
        albedo = ALBEDO_ICY_PLANET;
    }
    else { // gas_p
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

// New: Function to determine planet type based on accumulated mass and critical mass
void Nucleus::determine_planet_type(DoleParams* params) {
    double mass_critical_for_gas = params->masscritical(this->axis); // Critical mass at current orbit

    // Check for Gas Giant first
    if (this->mass >= mass_critical_for_gas) {
        this->type = Nucleus::gas_p;
        return;
    }

    // If not a gas giant, determine based on solid composition
    double total_solid_mass = accumulated_iron_mass + accumulated_rock_mass + accumulated_ice_mass;

    if (total_solid_mass > 0) {
        double iron_proportion = accumulated_iron_mass / total_solid_mass;
        double ice_proportion = accumulated_ice_mass / total_solid_mass;

        // Arbitrary thresholds for now, these can be tuned.
        // For instance, if iron > 50% of total solid, it's iron_p
        if (iron_proportion > 0.5) {
            this->type = Nucleus::iron_p;
        }
        // If ice is significant (> 50% for example), it's ice_p
        else if (ice_proportion > 0.5) {
            this->type = Nucleus::ice_p;
        }
        // Otherwise, it's a rocky planet
        else {
            this->type = Nucleus::rock_p;
        }
    } else {
        // If no solid mass accumulated (shouldn't happen for a formed nucleus but for safety)
        this->type = Nucleus::rock_p; // Default to rock if no material accreted
    }
}


// Define a struct to hold material availability for each band
struct BandMaterial {
    double dust_iron;
    double dust_rock;
    double dust_ice;
    // double gas; // Gas is handled implicitly via masscritical / gasdensity
    //           // If explicit gas tracking is needed, this would be uncommented
    //           // and its usage defined.

    BandMaterial() : dust_iron(0.0), dust_rock(0.0), dust_ice(0.0) {} // Initialize all to 0
};


//void export_planets_to_csv(const std::vector<Nucleus>& planets, const std::string& filename); 

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
            std::cout << "Pölyn koostumus lumirajan sisäpuolella: rauta " << P_IRON_ROCKY << ", kivi " << P_ROCK_ROCKY << ".\n";
            std::cout << "Pölyn koostumus lumirajan ulkopuolella: jää " << P_ICE_ICY << ", kivi " << P_ROCK_ICY << ".\n";
            std::cout << "Lumirajan läheisyydessä pölyn määrä kasvaa kertoimella: " << SNOW_LINE_DUST_FACTOR << ".\n";
            return 0; // Poistu ohjelmasta ohjeen näyttämisen jälkeen
        }
    }
    // --- End Help Option Check ---

    int dumpflag = 0; // Default to no console dump unless -dump is specified

    long seed = time(0);
    int MaxPlanets = 1000;     // orig 20 Size of dynamic array containing generated planets

    double aumin = 0.3;  // Changed to double for precision
    double aumax = 50;   // Changed to double for precision

    // Default DoleParams values (as in the DoleParams constructor)
    double pA = .00150; // 0.00150;
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

    double    nucleieccent, nucleiradius, masslast, eccent, mass, rperigee,
              rapogee, radius, mu, xperigee, xapogee, masscritical,
              density, bw, volume, da, dp; // masspass removed, now mass_acc_iron/rock/ice

    int      lowband, highband, dustcheck, iterate, hit;

    // Bands are measured in 1/5 au intervals from 0 to 50 AU
    const int
    MaxBand = 50 * 5,
    MaxNuclei = 1000;    // Max # of nuclei to inject

    // Replaced enum with struct for detailed material tracking
    BandMaterial band[MaxBand+1];

    // Changed to std::vector for dynamic sizing and easier management
    std::vector<Nucleus> planets;
    planets.reserve(MaxPlanets); // Pre-allocate memory

    int nplanet = 0;

    randomize(seed);

    ps_window(-1, -1, 2, 1);
    ps_begin("planets.ps");

    // Set up primordial cloud
    // Dust within .3 au is swept out by radiation pressure

    // Initialize bands with specific material compositions
    for (int n_band = 1; n_band <= MaxBand; n_band++) {
        double au = n_band / 5.0; // Convert band number back to AU

        // Initial total dust amount for this band segment (0.2 AU width)
        double initial_base_dust_in_band_segment = params->base_dust_density(au) * 0.2;

        if (au < SNOW_LINE_AU) {
            band[n_band].dust_iron = initial_base_dust_in_band_segment * P_IRON_ROCKY;
            band[n_band].dust_rock = initial_base_dust_in_band_segment * P_ROCK_ROCKY;
            band[n_band].dust_ice = 0.0;
        } else {
            // Apply the snow line factor to the total dust available in this region initially
            // This makes the *total* dust amount in these bands higher from the start
            double scaled_total_dust_for_ice_region = initial_base_dust_in_band_segment * SNOW_LINE_DUST_FACTOR;
            band[n_band].dust_iron = 0.0; // No distinct iron dust outside snow line
            band[n_band].dust_rock = scaled_total_dust_for_ice_region * P_ROCK_ICY;
            band[n_band].dust_ice = scaled_total_dust_for_ice_region * P_ICE_ICY;
        }
        // The `gas` field in BandMaterial is not explicitly used for initial filling here,
        // as gas accretion is currently handled via the `masscritical` and `gasdensity` function.
        // It implies gas is always "available" when a certain mass is reached.
    }


    // Draw scale on figure
    logscale("AU");

    // Plot density function; space samples even on log scale
    // This plot now represents the *total* initial dust density, including the snow line factor.
    double max_density_value = 0.0;
    for (double test_au = aumin; test_au <= aumax; test_au += 0.1) {
        double current_rho = params->total_dust_density(test_au);
        if (current_rho > max_density_value) {
            max_density_value = current_rho;
        }
    }

    double lastau, lastrho;
    double logaumin = std::log10(aumin);
    double logaumax = std::log10(aumax);
    double logaustep = (logaumax - logaumin) / 100;

    for (double logau = logaumin; logau <= logaumax + epsilon; logau += logaustep) {
        double au_val = std::pow(10,logau);
        double rho = params->total_dust_density(au_val) / max_density_value;

        if (logau > logaumin)
            ps_line(lastau, lastrho, logau, rho);
        else {
            char buf[100];
            std::sprintf(buf, "density [max = %g]", max_density_value);
            ps_text(logau, rho, buf);
        }

        lastau = logau;
        lastrho = rho;
    }

    double yBand = 0.9;
    // Inject nuclei into the cloud and acrete gas and dust
    int nuclei_count = 0; // Declare a new variable to keep track of nuclei processed
    for (nuclei_count = 0; nuclei_count < MaxNuclei; nuclei_count++) {
        // Check to see if all bands are taken (i.e., no dust left)
        int band_has_dust = 0;
        for (int i = 1; i <= MaxBand; i++) {
            if (band[i].dust_iron > 0 || band[i].dust_rock > 0 || band[i].dust_ice > 0)  {
                band_has_dust = 1;
                break;
            }
        }

        if (band_has_dust == 0)
            break; // All dust collected, exit

        nucleieccent = 1 - std::pow(1 - rand01(), .077);
        while ((nucleiradius = rand01() * 50) < .3)
            ;

        Nucleus body;
        body.axis = nucleiradius;
        body.eccen = nucleieccent;
        body.mass = params->MassNuclei;
        body.type = Nucleus::rock_p; // Initial default type, will be refined later
        body.accumulated_iron_mass = 0.0;
        body.accumulated_rock_mass = 0.0;
        body.accumulated_ice_mass = 0.0;


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

            mu = mass / params->MassSun;
            xperigee = rperigee;
            xapogee = rapogee;
            body.pAttr = rperigee * std::pow(mu, 1 / 3.0);
            body.aAttr = rapogee * std::pow(mu, 1 / 3.0);

            lowband = cint(body.lowbound(params) * 5.0);
            highband = cint(body.highbound(params) * 5.0);

            lowband = imax(1, lowband);
            highband = imin(MaxBand, highband);

            dustcheck = 0;
            // Iterate over bands and accrete material
            for (int bandnum = lowband; bandnum <= highband; bandnum++) {
                double band_au = bandnum / 5.0; // AU for the center of this band

                double current_iron_dust = band[bandnum].dust_iron;
                double current_rock_dust = band[bandnum].dust_rock;
                double current_ice_dust = band[bandnum].dust_ice;

                double total_dust_in_band = current_iron_dust + current_rock_dust + current_ice_dust;

                if (total_dust_in_band > 0) {
                    dustcheck = 1;

                    // Accrete a fraction of the dust from this band
                    double mass_accreted_from_band_factor = rand01() * params->W;
                    
                    double accreted_iron = current_iron_dust * mass_accreted_from_band_factor;
                    double accreted_rock = current_rock_dust * mass_accreted_from_band_factor;
                    double accreted_ice = current_ice_dust * mass_accreted_from_band_factor;

                    // Update planet's accumulated material
                    body.accumulated_iron_mass += accreted_iron;
                    body.accumulated_rock_mass += accreted_rock;
                    body.accumulated_ice_mass += accreted_ice;

                    // Update band's remaining dust
                    band[bandnum].dust_iron -= accreted_iron;
                    band[bandnum].dust_rock -= accreted_rock;
                    band[bandnum].dust_ice -= accreted_ice;

                    body.mass += (accreted_iron + accreted_rock + accreted_ice);
                }
            }

            if (dustcheck == 0) { // If no dust was found in any band
                break;
            }

            if (body.mass > params->masscritical(body.axis) && masslast < params->masscritical(body.axis)) {
                body.mass = params->masscritical(body.axis) * params->K; // Accrete gas!
            }

            if (body.mass == masslast) {
                if (iterate > 100)
                    break;
                else
                    iterate++;
            } else {
                masslast = body.mass;
                iterate = 0;
            }

            // Apply migration during accretion if KMigration is set
            if (params->KMigration > 0 && body.mass > params->MassNuclei) { // Only migrate if it's grown
                double migration_factor = 1.0 - (params->KMigration * std::log10(body.mass / params->MassNuclei));
                if (migration_factor < 0.01) migration_factor = 0.01; // Avoid migration to near zero AU
                body.axis = nucleiradius * migration_factor;
                // If migration changes axis, adjust perigee/apogee/eccentricity as well to maintain orbit shape relative to new axis
                body.pRad = body.axis * (1 - body.eccen);
                body.aRad = body.axis * (1 + body.eccen);
                // Re-calculate attractive distances based on new axis/mass
                body.pAttr = body.pRad * std::pow(body.mass / params->MassSun, 1 / 3.0);
                body.aAttr = body.aRad * std::pow(body.mass / params->MassSun, 1 / 3.0);
            }
        } // End of mass accretion loop

        // After accretion, apply final migration if QMigration (MigrationFinalFraction) is set
        if (params->MigrationFinalFraction != 1.0) {
            body.axis = body.axis * params->MigrationFinalFraction;
            body.pRad = body.axis * (1 - body.eccen);
            body.aRad = body.axis * (1 + body.eccen);
            body.pAttr = body.pRad * std::pow(body.mass / params->MassSun, 1 / 3.0);
            body.aAttr = body.aRad * std::pow(body.mass / params->MassSun, 1 / 3.0);
        }

        // Only add to planets if it's within the valid AU range after potential migration
        if (body.axis >= aumin && body.axis <= aumax) {
            body.determine_planet_type(params); // Determine the final type after all accretion
            planets.push_back(body); // Add the formed nucleus to the vector of planets
            nplanet++;

            if (nplanet >= MaxPlanets) break; // Stop if max planets reached
        }
    } // End of nuclei injection loop

    // Sort planets by their semi-major axis
 
    qsort(planets.data(), nplanet, sizeof(Nucleus), nucleusCompare);

    std::cout << "\nSimulaatio valmis. Löydetty " << nplanet << " planeettaa.\n\n";


	export_planets_to_csv(planets, "planet_data.csv");
	
	exit(-1);

    // Output all formed planets
    for (int i = 0; i < nplanet; i++) {
        planets[i].dump(i + 1, std::cout, params);
        std::cout << "\n----------------------------------------\n\n";
    }

    // Draw the generated planets
    ps_text(std::log10(aumin), 0.8, "Planets");
    for (int i = 0; i < nplanet; i++) {
        ps_pset(std::log10(planets[i].axis), yBand);
        char buf[100];
        std::sprintf(buf, "%d", i+1);
        ps_text(std::log10(planets[i].axis), yBand + 0.02, buf);
    }
    ps_showpage();

    // Plot Planet Radii vs. AU
    ps_beginpage(PsPage++);
    ps_window(0, 0, 5, 20); // Adjust window for radii (e.g., 0-20 Earth radii or more if larger planets possible)
    ps_text(0.1, 19, "Planet Radius (1000s of km) vs. Orbital Distance (AU)");
    ps_line(0.1, 0.1, 4.9, 0.1); // X-axis
    ps_line(0.1, 0.1, 0.1, 19.9); // Y-axis

    // X-axis labels (AU)
    for (int i = 1; i <= 5; ++i) {
        ps_line(0.1 + (i-1)*0.96, 0.1, 0.1 + (i-1)*0.96, 0.08); // Tick marks
        char label[10];
        std::sprintf(label, "%d AU", i*10);
        ps_text(0.1 + (i-1)*0.96, 0.05, label);
    }

    // Y-axis labels (Radius in 1000s of km)
    for (int i = 0; i <= 19; ++i) {
        ps_line(0.1, 0.1 + i, 0.08, 0.1 + i); // Tick marks
        char label[10];
        std::sprintf(label, "%d", i);
        ps_text(0.05, 0.1 + i, label);
    }


    for (int i = 0; i < nplanet; i++) {
        double radius_km = planets[i].get_radius_km();
        // Scale radius for plotting (e.g., if max radius is 70,000 km, divide by 1000 for y-axis range 0-70)
        double plot_radius_y = radius_km / 1000.0; // Plotting in 1000s of km

        // Scale AU for plotting (e.g., to fit 0-50 AU into x-axis 0.1-4.9)
        double plot_au_x = 0.1 + (planets[i].axis / aumax) * (4.9 - 0.1);

        ps_circle(plot_au_x, plot_radius_y, 0.05, 0); // Plot planet as circle
        char buf[100];
        std::sprintf(buf, "%d", i+1);
        ps_text(plot_au_x + 0.05, plot_radius_y + 0.05, buf);
    }
    ps_showpage();


    // Plot Planet Temperature vs. AU
    ps_beginpage(PsPage++);
    ps_window(0, 0, 5, 400); // Adjust window for temperature (e.g., 0-400K)
    ps_text(0.1, 390, "Planet Temperature (K) vs. Orbital Distance (AU)");
    ps_line(0.1, 0.1, 4.9, 0.1); // X-axis
    ps_line(0.1, 0.1, 0.1, 399.9); // Y-axis

    // X-axis labels (AU) - reuse previous AU scaling
    for (int i = 1; i <= 5; ++i) {
        ps_line(0.1 + (i-1)*0.96, 0.1, 0.1 + (i-1)*0.96, 0.08); // Tick marks
        char label[10];
        std::sprintf(label, "%d AU", i*10);
        ps_text(0.1 + (i-1)*0.96, 0.05, label);
    }

    // Y-axis labels (Temperature in K)
    for (int i = 0; i <= 400; i += 50) {
        ps_line(0.1, (double)i, 0.08, (double)i); // Tick marks
        char label[10];
        std::sprintf(label, "%d K", i);
        ps_text(0.05, (double)i, label);
    }


    for (int i = 0; i < nplanet; i++) {
        double temp_k = planets[i].get_temperature_k();
        double plot_au_x = 0.1 + (planets[i].axis / aumax) * (4.9 - 0.1);

        ps_circle(plot_au_x, temp_k, 0.05, 0); // Plot planet as circle
        char buf[100];
        std::sprintf(buf, "%d", i+1);
        ps_text(plot_au_x + 0.05, temp_k + 5, buf);
    }
    ps_showpage();


    ps_end();
    delete params; // Clean up DoleParams object

    return 0;
}



void export_planets_to_csv(const std::vector<Nucleus>& planets, const std::string& filename) {
    std::ofstream csv_file(filename);

    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open CSV file " << filename << " for writing.\n";
        return;
    }

    // Write CSV header
    csv_file << "radius_au,eccentricity,mass_mearths,planet_type,p_iron,p_rock,p_ice,p_gas\n";

    // Create a mutable copy of the planets vector to sort
    std::vector<Nucleus> sorted_planets = planets;

    // Sort planets by radius (axis) in ascending order
    std::sort(sorted_planets.begin(), sorted_planets.end());

    // Write planet data
    for (const auto& planet : sorted_planets) {
        double mass_earth_masses = planet.mass * SOLAR_MASS_TO_EARTH_MASS;

        double total_solid_mass = planet.accumulated_iron_mass + planet.accumulated_rock_mass + planet.accumulated_ice_mass;
        double gas_mass = planet.mass - total_solid_mass;
        if (gas_mass < 0) gas_mass = 0; // Ensure non-negative gas mass

        double total_mass_for_fractions = total_solid_mass + gas_mass;

        double p_iron = 0.0;
        double p_rock = 0.0;
        double p_ice = 0.0;
        double p_gas = 0.0;

        if (total_mass_for_fractions > 0) {
            p_iron = planet.accumulated_iron_mass / total_mass_for_fractions;
            p_rock = planet.accumulated_rock_mass / total_mass_for_fractions;
            p_ice = planet.accumulated_ice_mass / total_mass_for_fractions;
            p_gas = gas_mass / total_mass_for_fractions;
        }

        std::string planet_type_str;
        switch (planet.type) {
            case Nucleus::iron_p: planet_type_str = "Iron"; break;
            case Nucleus::rock_p: planet_type_str = "Rocky"; break;
            case Nucleus::ice_p:  planet_type_str = "Icy"; break;
            case Nucleus::gas_p:  planet_type_str = "Gas Giant"; break;
        }

        csv_file << std::fixed << std::setprecision(6)
                 << planet.axis << ","
                 << planet.eccen << ","
                 << mass_earth_masses << ","
                 << planet_type_str << ","
                 << p_iron << ","
                 << p_rock << ","
                 << p_ice << ","
                 << p_gas << "\n";
    
    
    if( mass_earth_masses>1e-4) {
    std::cout<<planet.axis<< ","    << mass_earth_masses << ","<< planet_type_str << "\n";
	
		}
    
    }

    csv_file.close();
    std::cout << "Planet data exported to " << filename << std::endl;
}


