

//
// acretee  v 0000.0000.0008
// based on acrete
//
// use own surface density profile, pov ray output
//

#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>



extern "C" {
    double drand48();
    long     lrand48();
    void     srand48(long);
};

const double SOLAR_MASS_TO_EARTH_MASS = 333000.0;
const double EARTH_RADIUS_METERS = 6.371e6;
const double EARTH_MASS_KG = 5.972e24;
const double JUPITER_DENSITY_KG_M3 = 1330.0;
const double EARTH_DENSITY_KG_M3 = 5510.0;
const double ICE_DENSITY_KG_M3 = 917.0;
const double IRON_DENSITY_KG_M3 = 7874.0;

const double AU_TO_METERS = 1.496e11;
const double SUN_SURFACE_TEMP_K = 5778.0;
const double SUN_RADIUS_METERS = 6.957e8;

const double ALBEDO_ROCKY_PLANET = 0.3;
const double ALBEDO_GAS_GIANT = 0.5;
const double ALBEDO_ICY_PLANET = 0.7;

const double SNOW_LINE_AU = 5.0;

const double P_IRON_ROCKY = 0.3;
const double P_ROCK_ROCKY = 0.7;

const double P_ICE_ICY = 0.75;
const double P_ROCK_ICY = 0.25;

const double SNOW_LINE_DUST_FACTOR = (P_IRON_ROCKY + P_ROCK_ROCKY) / P_ICE_ICY; // (0.3+0.7) / 0.75 = 1.333...

struct DoleParams {
    double A, B, Alpha, Beta, Eccentricity, Gamma, K, MassSol, MassSun, MassNuclei, W, KMigration, MigrationFinalFraction;

    DoleParams();
    DoleParams(double pA, double pAlpha, double pBeta, double pEccentricity, double pGamma, double pK, double pMassSol, double pMassSun, double pMassNuclei, double pB, double pW, double pKMigration, double pMigrationFinalFraction);

    double density(double au);
    double gasdensity(double au, double massratio);
    double masscritical(double au);
    double lowbound(double radius, double margin);
    double highbound(double radius, double margin);

    double base_dust_density(double au);
    double dust_density_iron(double au);
    double dust_density_rock(double au);
    double dust_density_ice(double au);
    double total_dust_density(double au);
};

DoleParams::DoleParams() {
    A = .00150 * 100;
    Alpha = 5;
    Gamma = 1 / 3.0;
    Beta = 0.5;
    Eccentricity = 0.15;
    K = 50;
    MassSol = 1;
    MassSun = MassSol * 1.0;
    MassNuclei = MassSun * 1e-15;
    B = MassSun * 1.2e-5 * 3;
    W = .2;
    KMigration = 0.0;
    MigrationFinalFraction = 1.0;
}

DoleParams::DoleParams(double pA, double pAlpha, double pBeta, double pEccentricity, double pGamma, double pK, double pMassSol, double pMassSun, double pMassNuclei, double pB, double pW, double pKMigration, double pMigrationFinalFraction)
{
    A = pA;
    Alpha = pAlpha;
    Beta = pBeta;
    Eccentricity = pEccentricity;
    Gamma = pGamma;
    K = pK;
    MassSol = pMassSol;
    MassSun = pMassSun;
    MassNuclei = pMassNuclei;
    B = pB;
    W = pW;
    KMigration = pKMigration;
    MigrationFinalFraction = pMigrationFinalFraction;
}

double DoleParams::base_dust_density(double au) {
	// we use own dist densioty not dole's
    double den = 0.0;
    if (au < 0.001) au = 0.001;
    den = A * std::pow(au / 2.5, -1) * 0.0020;
    if (au < 4) { den = A * 1e-6; }
    if (au < 2) { den = A * 0.1; }
    return den;
}

double DoleParams::dust_density_iron(double au) {
    if (au < SNOW_LINE_AU) {
        return base_dust_density(au) * P_IRON_ROCKY;
    }
    return 0.0;
}

double DoleParams::dust_density_rock(double au) {
    if (au < SNOW_LINE_AU) {
        return base_dust_density(au) * P_ROCK_ROCKY;
    } else {
        return base_dust_density(au) * P_ROCK_ICY;
    }
}

double DoleParams::dust_density_ice(double au) {
    if (au >= SNOW_LINE_AU) {
        return base_dust_density(au) * P_ICE_ICY;
    }
    return 0.0;
}

double DoleParams::total_dust_density(double au) {
    double base_rho = base_dust_density(au);
    if (au >= SNOW_LINE_AU) {
        return base_rho * SNOW_LINE_DUST_FACTOR;
    } else {
        return base_rho;
    }
}

double DoleParams::density(double au) {
    double val = total_dust_density(au);
    if (val < 1e-30) return 1e-30;
    return val;
}

double DoleParams::gasdensity(double au, double massratio) {
    double total_dust = total_dust_density(au);
    if (total_dust < 1e-30) total_dust = 1e-30;

    double denominator = (1 + std::sqrt(massratio) * (K - 1));
    if (denominator == 0) denominator = 1e-30;

    return (K * total_dust) / denominator;
}

double DoleParams::masscritical(double au) {
    if (au <= 0) au = 0.001;
    return B * std::pow(au, -.75);
}

double DoleParams::lowbound(double radius, double margin) {
    double denominator = (1 + W);
    if (denominator == 0) denominator = 1e-30;
    return radius - margin - W * (radius - margin) / denominator;
}

double DoleParams::highbound(double radius, double margin) {
    double denominator = (1 - W);
    if (denominator == 0) denominator = 1e-30;
    return radius + margin + W * (radius + margin) / denominator;
}

const double epsilon = .00001;
const int Maxint = 32767;

int imin(int a, int b = Maxint, int c = Maxint) {
    if (a < b)
    return (a < c) ? a : c;
    else
    return (b < c) ? b : c;
}

// Corrected imax function logic
int imax(int a, int b = -Maxint, int c = -Maxint) {
    if (a > b) {
        return (a > c) ? a : c;
    } else {
        return (b > c) ? b : c; // Corrected from (b < c)
    }
}

std::ofstream *Ps;
int PsPage = 1;
double PsBase = 50, PsXscale = 1, PsYscale = 1, PsXoff = 0, PsYoff = 0;

void ps_end() {
    if (Ps != NULL) {
        *Ps << "%%Trailer\nend" << std::flush;
        delete Ps;
        Ps = NULL;
    }
}

void ps_beginpage(int page) {
    if (Ps == NULL) return;
    *Ps << "%%Page: " << page << ' ' << page << '\n'
    << PsXoff + PsBase << ' ' << PsYoff + PsBase << " translate\n"
    << PsXscale << ' ' << PsYscale << ' ' << " scale\n"
    << "/Helvetica findfont " << 9 / PsXscale << " scalefont setfont\n"
    << "0 setlinewidth\n";
}

void ps_begin(const char *file) {
    if (Ps != NULL) {
        ps_end();
    }

    Ps = new std::ofstream(file);
    if (!Ps->is_open()) {
        std::cerr << "VIRHE: PostScript-tiedostoa '" << file << "' ei voitu avata." << std::endl;
        delete Ps;
        Ps = NULL;
        return;
    }

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
    if (Ps == NULL) return;
    *Ps << "showpage\n";
    ps_beginpage(PsPage++);
}

void ps_window(double x1, double y1, double x2, double y2)
{
    const double width = 450;
    double xspan = x2 - x1, yspan = y2 - y1;

    if (xspan == 0) xspan = 1.0;
    if (yspan == 0) yspan = 1.0;

    PsXscale = width / xspan;
    PsYscale = PsXscale;
    PsXoff = -PsXscale * x1;
    PsYoff = -PsYscale * y1;
}

void ps_circle(double x, double y, double radius, int fill)
{
    if (Ps == NULL) return;
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

int cint(double arg)
{
    return static_cast<int>(std::floor(arg + 0.5));
}

void ps_pset(double x, double y)
{
    ps_circle(x, y, 0.01, 0);
}

void ps_line(double x1, double y1, double x2, double y2)
{
    if (Ps == NULL) return;
    *Ps << x1 << ' ' << y1 << " moveto "
    << x2 << ' ' << y2 << " lineto stroke\n";
}

void ps_text(double x, double y, const char *s) {
    if (Ps == NULL) return;
    *Ps << x << ' ' << y << " moveto (" << s << ") show newpath\n";
}

void logscale(const char *xlabel = "", const char *ylabel = "") {
    ps_line(-1, -1, 3, -1);
    ps_line(3, -1, 3, 1);
    ps_line(3, 1, 3, -1);
    ps_line(3, -1, -1, -1);

    ps_line(-1, 1, 3, 1);
    for (double au = 1; au <= 10 + epsilon; au += 1) {
        if (au / 10.0 > 0) ps_line(std::log10(au / 10.0), 1, std::log10(au / 10.0), .95);
        if (au > 0) ps_line(std::log10(au), 1, std::log10(au), .95);
        if (au * 10.0 > 0) ps_line(std::log10(au * 10.0), 1, std::log10(au * 10.0), .95);
    }

    ps_text(-1, 1, ".1");
    ps_text(0, 1, "1");
    ps_text(1, 1, "10");
    ps_text(2, 1, "100");

    ps_text(2.3, 1, xlabel);
    ps_text(-1, .9, ylabel);
}

struct Nucleus {
    double axis, eccen, mass, pRad, aRad, pAttr, aAttr;

    enum PlanetType {
        iron_p, rock_p, ice_p, gas_p
    } type;

    double accumulated_iron_mass;
    double accumulated_rock_mass;
    double accumulated_ice_mass;

    Nucleus();
    void dump(int, std::ostream &, DoleParams *) const;
    friend int nucleusCompare(const void *, const void *);
    double lowbound(DoleParams *params) const;
    double highbound(DoleParams *params) const;

    double get_radius_km() const;
    double get_temperature_k() const;
    void determine_planet_type(DoleParams *params);
    bool operator<(const Nucleus &other) const {
        return axis < other.axis;
    }
};



Nucleus::Nucleus() {
    axis = eccen = mass = 0;
    pRad = aRad = 0;
    pAttr = aAttr = 0;
    type = Nucleus::rock_p;
    accumulated_iron_mass = 0.0;
    accumulated_rock_mass = 0.0;
    accumulated_ice_mass = 0.0;
}

double Nucleus::lowbound(DoleParams *params) const {
    return params->lowbound(pRad, pAttr);
}

double Nucleus::highbound(DoleParams *params) const {
    return params->highbound(aRad, aAttr);
}

int nucleusCompare(const void *p1, const void *p2) {
    double r1 = static_cast<const Nucleus *>(p1)->axis,
           r2 = static_cast<const Nucleus *>(p2)->axis;
    return (r1 < r2) ? -1 : (r1 > r2) ? 1 : 0;
}

void Nucleus::dump(int num, std::ostream &o, DoleParams *params) const {
    double xplimit = pRad - pAttr,
           xalimit = aRad + aAttr,
           lowrange = lowbound(params),
           highrange = highbound(params),
           massCrit = params->masscritical(pRad);

    o << "Nucleus " << num << '\n';
    o << "\tRadius\t\t" << std::fixed << std::setprecision(4) << axis << " AU\n";
    o << "\tEccentricity\t" << std::fixed << std::setprecision(4) << eccen << '\n';

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

    double total_solid_mass = accumulated_iron_mass + accumulated_rock_mass + accumulated_ice_mass;
    double gas_mass = mass - total_solid_mass;
    if (gas_mass < 0) gas_mass = 0;

    double total_mass_for_fractions = total_solid_mass + gas_mass;

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

double Nucleus::get_radius_km() const {
    double mass_kg = this->mass * SOLAR_MASS_TO_EARTH_MASS * EARTH_MASS_KG;
    double density_kg_m3;

    if (this->type == Nucleus::rock_p) {
        density_kg_m3 = EARTH_DENSITY_KG_M3;
    } else if (this->type == Nucleus::ice_p) {
        density_kg_m3 = ICE_DENSITY_KG_M3;
    } else if (this->type == Nucleus::iron_p) {
        density_kg_m3 = IRON_DENSITY_KG_M3;
    } else {
        density_kg_m3 = JUPITER_DENSITY_KG_M3;
    }

    if (density_kg_m3 <= 0 || mass_kg <= 0) {
        return 0.0;
    }
    double radius_meters = std::cbrt((3 * mass_kg) / (4 * M_PI * density_kg_m3));
    return radius_meters / 1000.0;
}

double Nucleus::get_temperature_k() const {
    double albedo;
    if (this->type == Nucleus::rock_p || this->type == Nucleus::iron_p) {
        albedo = ALBEDO_ROCKY_PLANET;
    } else if (this->type == Nucleus::ice_p) {
        albedo = ALBEDO_ICY_PLANET;
    } else {
        albedo = ALBEDO_GAS_GIANT;
    }

    double distance_meters = this->axis * AU_TO_METERS;
    if (distance_meters <= 0) return 0;

    double term_distance = SUN_RADIUS_METERS / (2 * distance_meters);
    if (term_distance < 0) term_distance = 0;

    double temp_eq = SUN_SURFACE_TEMP_K * std::sqrt(term_distance) * std::pow(1 - albedo, 0.25);

    return temp_eq;
}

void Nucleus::determine_planet_type(DoleParams *params) {
    double mass_critical_for_gas = params->masscritical(this->axis);

    if (this->mass >= mass_critical_for_gas * (1.0 - epsilon)) {
        this->type = Nucleus::gas_p;
        return;
    }

    double total_solid_mass = accumulated_iron_mass + accumulated_rock_mass + accumulated_ice_mass;

    if (total_solid_mass > 0) {
        double iron_proportion = accumulated_iron_mass / total_solid_mass;
        double ice_proportion = accumulated_ice_mass / total_solid_mass;

        if (iron_proportion > 0.5) {
            this->type = Nucleus::iron_p;
        } else if (ice_proportion > 0.5) {
            this->type = Nucleus::ice_p;
        } else {
            this->type = Nucleus::rock_p;
        }
    } else {
        this->type = Nucleus::rock_p;
    }
}

struct BandMaterial {
    double dust_iron;
    double dust_rock;
    double dust_ice;

    BandMaterial() : dust_iron(0.0), dust_rock(0.0), dust_ice(0.0) {}
};

int generate_and_render_povray(std::vector<Nucleus> planets, int planetnum) ;

void export_planets_to_csv(std::vector<Nucleus>  planets, std::string filename);

int main(int ac, char *av[])
{
    for (int i = 1; i < ac; i++) {
        if (strcmp(av[i], "-h") == 0 || strcmp(av[i], "--help") == 0) {
            std::cout << "Käyttö: " << av[0] << " [asetukset]\n";
            std::cout << "Simuloi planeettojen muodostumista ja migraatiota Dole-mallin pohjalta.\n\n";
            std::cout << "Asetukset:\n";
            std::cout << "  -h, --help           Näytä tämä ohjeviesti ja poistu.\n";
            std::cout << "  -dump                Tulosta yksityiskohtaisia tietoja konsoliin.\n";
            std::cout << "  -seed <luku>         Aseta satunnaislukugeneraattorin siemen (oletus: nykyinen aika).\n";
            std::cout << "  -maxplanets <luku>   Maksimi planeettojen määrä (oletus: 1000).\n";
            std::cout << "  -aumin <AU>          Minimi etäisyys (AU) planeetan muodostumiselle (oletus: 0.3).\n";
            std::cout << "  -aumax <AU>          Maksimi etäisyys (AU) planeetan muodostumiselle (oletus: 50.0).\n";
            std::cout << "  -A <arvo>            Dole-mallin parametri A (oletus: 0.00150).\n";
            std::cout << "  -Alpha <arvo>        Dole-mallin parametri Alpha (oletus: 5.0).\n";
            std::cout << "  -Beta <arvo>         Dole-mallin parametri Beta (oletus: 0.5).\n";
            std::cout << "  -Eccentricity <arvo> Dole-mallin parametri Eccentricity (oletus: 0.15).\n";
            std::cout << "  -Gamma <arvo>        main Dole-mallin parametri Gamma (oletus: 0.333333).\n";
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
            std::cout << "  -AccretionFrac <kerroin> Fraktio pölystä, joka kerätään per iteraatio (oletus: 0.001).\n";
            std::cout << "  -MaxAccrIter <luku> Maksimi akkretion iteraatiot per planeetta (oletus: 1000000).\n";

            std::cout << "\nLumiraja asetettu arvoon: " << SNOW_LINE_AU << " AU.\n";
            std::cout << "Pölyn koostumus lumirajan sisäpuolella: rauta " << P_IRON_ROCKY << ", kivi " << P_ROCK_ROCKY << ".\n";
            std::cout << "Pölyn koostumus lumirajan ulkopuolella: jää " << P_ICE_ICY << ", kivi " << P_ROCK_ICY << ".\n";
            std::cout << "Lumirajan läheisyydessä pölyn määrä kasvaa kertoimella: " << SNOW_LINE_DUST_FACTOR << ".\n";
            return 0;
        }
    }

    int dumpflag = 0;
    long seed = time(0);
    int MaxPlanets = 1000;
    double aumin = 0.3;
    double aumax = 50;

    double pA = .00150 * 20;
    double pAlpha = 5;
    double pBeta = 0.5;
    double pEccentricity = 0.15;
    double pGamma = 1 / 3.0;
    double pK = 50;
    double pMassSol = 1;
    double pMassSun = pMassSol * 1.0;
    double pMassNuclei = pMassSun * 1e-15;
    double pB = pMassSun * 1.2e-5 * 3;
    double pW = .2;
    double pKMigration = 0.0;
    double pMigrationFinalFraction = 1.0;

    double ACCRETION_FRACTION_PER_ITERATION = 0.001;
    int MAX_ACCRETION_ITERATIONS = 1000000;

    for (int i = 1; i < ac; i++) {
        if (strcmp(av[i], "-dump") == 0) {
            dumpflag = 1;
        } else if (strcmp(av[i], "-seed") == 0) {
            if (i + 1 < ac) { seed = std::atol(av[++i]); } else { std::cerr << "Error: -seed requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-maxplanets") == 0) {
            if (i + 1 < ac) { MaxPlanets = std::atoi(av[++i]); } else { std::cerr << "Error: -maxplanets requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-aumin") == 0) {
            if (i + 1 < ac) { aumin = std::atof(av[++i]); } else { std::cerr << "Error: -aumin requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-aumax") == 0) {
            if (i + 1 < ac) { aumax = std::atof(av[++i]); } else { std::cerr << "Error: -aumax requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-A") == 0) {
            if (i + 1 < ac) pA = std::atof(av[++i]); else { std::cerr << "Error: -A requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-Alpha") == 0) {
            if (i + 1 < ac) pAlpha = std::atof(av[++i]); else { std::cerr << "Error: -Alpha requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-Beta") == 0) {
            if (i + 1 < ac) pBeta = std::atof(av[++i]); else { std::cerr << "Error: -Beta requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-Eccentricity") == 0) {
            if (i + 1 < ac) pEccentricity = std::atof(av[++i]); else { std::cerr << "Error: -Eccentricity requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-Gamma") == 0) {
            if (i + 1 < ac) pGamma = std::atof(av[++i]); else { std::cerr << "Error: -Gamma requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-K") == 0) {
            if (i + 1 < ac) pK = std::atof(av[++i]); else { std::cerr << "Error: -K requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-MassSol") == 0) {
            if (i + 1 < ac) pMassSol = std::atof(av[++i]); else { std::cerr << "Error: -MassSol requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-MassSun") == 0) {
            if (i + 1 < ac) pMassSun = std::atof(av[++i]); else { std::cerr << "Error: -MassSun requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-MassNuclei") == 0) {
            if (i + 1 < ac) pMassNuclei = std::atof(av[++i]); else { std::cerr << "Error: -MassNuclei requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-B") == 0) {
            if (i + 1 < ac) pB = std::atof(av[++i]); else { std::cerr << "Error: -B requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-W") == 0) {
            if (i + 1 < ac) pW = std::atof(av[++i]); else { std::cerr << "Error: -W requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-KMigration") == 0) {
            if (i + 1 < ac) pKMigration = std::atof(av[++i]); else { std::cerr << "Error: -KMigration requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-QMigration") == 0) {
            if (i + 1 < ac) pMigrationFinalFraction = std::atof(av[++i]); else { std::cerr << "Error: -QMigration requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-AccretionFrac") == 0) {
            if (i + 1 < ac) ACCRETION_FRACTION_PER_ITERATION = std::atof(av[++i]); else { std::cerr << "Error: -AccretionFrac requires a value." << std::endl; return 1; }
        } else if (strcmp(av[i], "-MaxAccrIter") == 0) {
            if (i + 1 < ac) MAX_ACCRETION_ITERATIONS = std::atoi(av[++i]); else { std::cerr << "Error: -MaxAccrIter requires a value." << std::endl; return 1; }
        } else {
            std::cerr << "Tuntematon argumentti: " << av[i] << std::endl;
            std::cerr << "Käytä -h tai --help nähdäksesi ohjeet.\n";
            return 1;
        }
    }

    DoleParams *params = new DoleParams(pA, pAlpha, pBeta, pEccentricity, pGamma, pK, pMassSol, pMassSun, pMassNuclei, pB, pW, pKMigration, pMigrationFinalFraction);

    double original_aumin = aumin;
    double effective_aumin_due_to_migration = original_aumin * params->MigrationFinalFraction;
    if (effective_aumin_due_to_migration < original_aumin) {
        aumin = effective_aumin_due_to_migration;
        if (aumin < 0.01) aumin = 0.01;
        std::cerr << "Huom: aumin säädetty migraation vuoksi arvoon: " << std::fixed << std::setprecision(4) << aumin << " AU\n";
    }

    double nucleieccent, nucleiradius, masslast, eccent, mass, rperigee,
              rapogee, radius, mu, xperigee, xapogee, masscritical,
              density, bw, volume, da, dp;

    int lowband, highband, dustcheck, iterate, hit;

    const int MaxBand = 50 * 5;
    BandMaterial band[MaxBand + 1];
    std::vector<Nucleus> planets;
    planets.reserve(MaxPlanets);
    int nplanet = 0;

    randomize(seed);

    ps_window(-1, -1, 2, 1);
    ps_begin("planets.ps");

    if (Ps == NULL) {
        std::cerr << "VIRHE: PostScript-tiedoston alustus epäonnistui. Lopetetaan." << std::endl;
        delete params;
        return 1;
    }

    for (int n_band = 1; n_band <= MaxBand; n_band++) {
        double au = n_band / 5.0;
        double initial_base_dust_in_band_segment = params->base_dust_density(au) * 0.2;

        if (au < SNOW_LINE_AU) {
            band[n_band].dust_iron = initial_base_dust_in_band_segment * P_IRON_ROCKY;
            band[n_band].dust_rock = initial_base_dust_in_band_segment * P_ROCK_ROCKY;
            band[n_band].dust_ice = 0.0;
        } else {
            double scaled_total_dust_for_ice_region = initial_base_dust_in_band_segment * SNOW_LINE_DUST_FACTOR;
            band[n_band].dust_iron = 0.0;
            band[n_band].dust_rock = scaled_total_dust_for_ice_region * P_ROCK_ICY;
            band[n_band].dust_ice = scaled_total_dust_for_ice_region * P_ICE_ICY;
        }
    }

    std::cerr << "\n--- Alkuperäinen pölytilanne (BandMaterial) ---\n";
    for (int n_band = 1; n_band <= MaxBand; n_band++) {
        double au = n_band / 5.0;
        if (band[n_band].dust_iron > 1e-15 || band[n_band].dust_rock > 1e-15 || band[n_band].dust_ice > 1e-15) {
            std::cerr << "Band " << n_band << " (" << std::fixed << std::setprecision(2) << au << " AU): "
                      << "Rauta: " << std::scientific << band[n_band].dust_iron
                      << ", Kivi: " << std::scientific << band[n_band].dust_rock
                      << ", Jää: " << std::scientific << band[n_band].dust_ice << std::endl;
        }
    }
    std::cerr << "--------------------------------------------------\n";

    logscale("AU");

    double max_density_value = 0.0;
    for (double test_au = aumin; test_au <= aumax; test_au += 0.1) {
        double current_rho = params->total_dust_density(test_au);
        if (current_rho > max_density_value) {
            max_density_value = current_rho;
        }
    }
    if (max_density_value == 0) max_density_value = 1.0;

    double lastau, lastrho;
    double logaumin;

    lastau = 0.1;
    lastrho = params->total_dust_density(lastau);

    double log10_lastau = (lastau > 0) ? std::log10(lastau) : std::log10(1e-10);
    double log10_lastrho_norm = (lastrho / max_density_value > 0) ? std::log10(lastrho / max_density_value) : std::log10(1e-10);

    ps_line(log10_lastau, log10_lastrho_norm,
            log10_lastau, log10_lastrho_norm);

    for (double au_plot = 0.1 + 0.01; au_plot <= aumax + epsilon; au_plot += 0.01) {
        double current_rho_plot = params->total_dust_density(au_plot);
        double log10_au_plot = (au_plot > 0) ? std::log10(au_plot) : std::log10(1e-10);
        double log10_current_rho_norm = (current_rho_plot / max_density_value > 0) ? std::log10(current_rho_plot / max_density_value) : std::log10(1e-10);

        ps_line(log10_lastau, log10_lastrho_norm,
                log10_au_plot, log10_current_rho_norm);
        lastau = au_plot;
        lastrho = current_rho_plot;
        log10_lastau = log10_au_plot;
        log10_lastrho_norm = log10_current_rho_norm;
    }

    ps_text(-.5, .5, "Dust Density");

    for (int n = 0; n < MaxPlanets; n++) {
        std::cerr << "Aloitetaan uuden nukleuksen muodostus (nro " << n + 1 << "/" << MaxPlanets << ")"<< std::endl;
        Nucleus new_planet;

        bool found_dust_band = false;
        int do_while_attempts = 0;
        const int MAX_DO_WHILE_ATTEMPTS = 1000;

        do {
            do_while_attempts++;
            if (do_while_attempts > MAX_DO_WHILE_ATTEMPTS) {
                std::cerr << "VAROITUS: Ylitetty yritysmäärä löytää pölyä. Hypätään nukleus " << n + 1 << " yli." << std::endl;
                found_dust_band = false;
                break;
            }

            new_planet.axis = rand01() * (aumax - aumin) + aumin;
            new_planet.eccen = params->Eccentricity * (1 - rand01() * rand01());
            new_planet.mass = params->MassNuclei;

            new_planet.pRad = new_planet.axis * (1 - new_planet.eccen);
            new_planet.aRad = new_planet.axis * (1 + new_planet.eccen);

            if (new_planet.pRad <= 0) new_planet.pRad = 0.001;
            if (new_planet.aRad <= 0) new_planet.aRad = 0.001;

            mu = new_planet.mass * params->MassSol;
            if (mu <= 0) mu = 1e-30; // Ensure mu is positive for cbrt

            double density_pRad = params->density(new_planet.pRad);
            double density_aRad = params->density(new_planet.aRad);

            if (density_pRad < 1e-30) density_pRad = 1e-30;
            if (density_aRad < 1e-30) density_aRad = 1e-30;

            double val_pAttr = mu / (2.0 * density_pRad);
            double val_aAttr = mu / (2.0 * density_aRad);

            // Clamp values before cbrt to prevent extreme results if they are out of reasonable range
            if (val_pAttr < 0) val_pAttr = 1e-30; // Should not be negative if mu and density are positive
            if (val_aAttr < 0) val_aAttr = 1e-30;

            new_planet.pAttr = new_planet.pRad * std::cbrt(val_pAttr);
            new_planet.aAttr = new_planet.aRad * std::cbrt(val_aAttr);

            // Additional clamping for pAttr and aAttr to prevent extremely large values
            const double MAX_ATTR_AU = 1000.0; // Define a reasonable max attraction distance in AU
            if (new_planet.pAttr > MAX_ATTR_AU) new_planet.pAttr = MAX_ATTR_AU;
            if (new_planet.aAttr > MAX_ATTR_AU) new_planet.aAttr = MAX_ATTR_AU;
            if (new_planet.pAttr < 0) new_planet.pAttr = 0.0; // Should not be negative
            if (new_planet.aAttr < 0) new_planet.aAttr = 0.0; // Should not be negative

            lowband = cint(new_planet.lowbound(params) * 5.0);
            highband = cint(new_planet.highbound(params) * 5.0);

            lowband = imax(1, lowband);
            highband = imin(MaxBand, highband);

            std::cerr << "  Nukleus " << n + 1 << ": Etsii pölyä bandeilta " << lowband << " - " << highband << std::endl;

            dustcheck = 0;
            for (int j = lowband; j <= highband; j++) {
                if (band[j].dust_iron > 0 || band[j].dust_rock > 0 || band[j].dust_ice > 0) {
                    dustcheck = 1;
                    break;
                }
            }
            if (dustcheck == 1) found_dust_band = true;
            else {
                std::cerr << "  Nukleus " << n + 1 << ": Ei pölyä etäisyydellä " << std::fixed << std::setprecision(4) << new_planet.axis << " AU. Uusi yritys." << std::endl;
            }

        } while (!found_dust_band);

        if (!found_dust_band) {
            continue;
        }

        hit = 1;
        iterate = 0;

        new_planet.accumulated_iron_mass = 0.0;
        new_planet.accumulated_rock_mass = 0.0;
        new_planet.accumulated_ice_mass = 0.0;

        while (hit && iterate < MAX_ACCRETION_ITERATIONS && (nplanet < MaxPlanets)) {
            iterate++;
            masslast = new_planet.mass;

            if (iterate % 10000 == 0 || iterate == 1 || iterate == MAX_ACCRETION_ITERATIONS - 1) {
                std::cerr << "    Akkretion iteraatio " << iterate
                          << ", nukleus " << n + 1
                          << ", massa: " << std::fixed << std::setprecision(10) << new_planet.mass * SOLAR_MASS_TO_EARTH_MASS << " M_Earth"
                          << ", akseli: " << std::fixed << std::setprecision(10) << new_planet.axis << " AU" << std::endl;
            }

            mu = new_planet.mass * params->MassSol;
            if (mu <= 0) mu = 1e-30; // Ensure mu is positive for cbrt

            if (new_planet.pRad <= 0) new_planet.pRad = 0.001;
            if (new_planet.aRad <= 0) new_planet.aRad = 0.001;

            double current_density_pRad = params->density(new_planet.pRad);
            double current_density_aRad = params->density(new_planet.aRad);
            if (current_density_pRad < 1e-30) current_density_pRad = 1e-30;
            if (current_density_aRad < 1e-30) current_density_aRad = 1e-30;

            double val_pAttr_accr = mu / (2.0 * current_density_pRad);
            double val_aAttr_accr = mu / (2.0 * current_density_aRad);

            if (val_pAttr_accr < 0) val_pAttr_accr = 1e-30;
            if (val_aAttr_accr < 0) val_aAttr_accr = 1e-30;

            new_planet.pAttr = new_planet.pRad * std::cbrt(val_pAttr_accr);
            new_planet.aAttr = new_planet.aRad * std::cbrt(val_aAttr_accr);

            // Additional clamping for pAttr and aAttr during accretion loop
            const double MAX_ATTR_AU_ACCR = 1000.0;
            if (new_planet.pAttr > MAX_ATTR_AU_ACCR) new_planet.pAttr = MAX_ATTR_AU_ACCR;
            if (new_planet.aAttr > MAX_ATTR_AU_ACCR) new_planet.aAttr = MAX_ATTR_AU_ACCR;
            if (new_planet.pAttr < 0) new_planet.pAttr = 0.0;
            if (new_planet.aAttr < 0) new_planet.aAttr = 0.0;

            lowband = cint(new_planet.lowbound(params) * 5.0);
            highband = cint(new_planet.highbound(params) * 5.0);

            lowband = imax(1, lowband);
            highband = imin(MaxBand, highband);

            for (int j = lowband; j <= highband; j++) {
                if (band[j].dust_iron > 0) {
                    double actual_accreted_iron = band[j].dust_iron * ACCRETION_FRACTION_PER_ITERATION;
                    new_planet.accumulated_iron_mass += actual_accreted_iron;
                    band[j].dust_iron -= actual_accreted_iron;
                    if (band[j].dust_iron < 0) band[j].dust_iron = 0;
                }

                if (band[j].dust_rock > 0) {
                    double actual_accreted_rock = band[j].dust_rock * ACCRETION_FRACTION_PER_ITERATION;
                    new_planet.accumulated_rock_mass += actual_accreted_rock;
                    band[j].dust_rock -= actual_accreted_rock;
                    if (band[j].dust_rock < 0) band[j].dust_rock = 0;
                }

                if (band[j].dust_ice > 0) {
                    double actual_accreted_ice = band[j].dust_ice * ACCRETION_FRACTION_PER_ITERATION;
                    new_planet.accumulated_ice_mass += actual_accreted_ice;
                    band[j].dust_ice -= actual_accreted_ice;
                    if (band[j].dust_ice < 0) band[j].dust_ice = 0;
                }
            }

            double total_solids_this_iteration = new_planet.accumulated_iron_mass + new_planet.accumulated_rock_mass + new_planet.accumulated_ice_mass;
            new_planet.mass = params->MassNuclei + total_solids_this_iteration;

            masscritical = params->masscritical(new_planet.axis);
            if (new_planet.mass >= masscritical) {
                double current_mass_ratio = new_planet.mass / masscritical;
                if (current_mass_ratio <= 0) current_mass_ratio = 1.0;

                new_planet.mass = new_planet.mass + params->gasdensity(new_planet.axis, masscritical / new_planet.mass) * (new_planet.aRad - new_planet.pRad);
                if (new_planet.mass > params->MassSun * 2) {
                    new_planet.mass = params->MassSun * 2;
                }
            }

            if (new_planet.mass > 1e-20 && (std::fabs(new_planet.mass - masslast) < (0.000001 * new_planet.mass)) && iterate > 100) {
                hit = 0;
                std::cerr << "    Akkretio pysäytetty: massa stabiloitui. \n";
            } else if (new_planet.mass <= params->MassNuclei * 10 && iterate > MAX_ACCRETION_ITERATIONS / 10) {
                hit = 0;
                std::cerr << "    Akkretio pysäytetty: massa pysyi erittäin pienenä pitkään. \n";
            } else if (iterate >= MAX_ACCRETION_ITERATIONS) {
                hit = 0;
                std::cerr << "    Akkretio pysäytetty: maksimi iteraatiot saavutettu. \n";
            }

            if (params->KMigration != 0.0) {
                new_planet.axis *= (1.0 - params->KMigration * (new_planet.mass / params->MassSun));
                if (new_planet.axis <= 0) {
                    new_planet.axis = 0.01;
                    std::cerr << "    VAROITUS: Akseli muuttui nollaksi tai negatiiviseksi migraation takia. Asetettu: " << new_planet.axis << std::endl;
                    hit = 0;
                }
                new_planet.pRad = new_planet.axis * (1 - new_planet.eccen);
                new_planet.aRad = new_planet.axis * (1 + new_planet.eccen);
            }
        }
        std::cerr << "Nukleuksen " << n + 1 << " akkretiovaihe päättynyt." << std::endl;

        if (params->MigrationFinalFraction != 1.0) {
            new_planet.axis *= params->MigrationFinalFraction;
            if (new_planet.axis <= 0) {
                new_planet.axis = 0.01;
                std::cerr << "    VAROITUS: Akseli muuttui nollaksi tai negatiiviseksi loppumigraation takia. Asetettu: " << new_planet.axis << std::endl;
            }
            new_planet.pRad = new_planet.axis * (1 - new_planet.eccen);
            new_planet.aRad = new_planet.axis * (1 + new_planet.eccen);
        }

        if (new_planet.axis >= aumin && new_planet.axis <= aumax && (new_planet.mass * SOLAR_MASS_TO_EARTH_MASS) > 0.0001) {
            new_planet.determine_planet_type(params);
            planets.push_back(new_planet);
            nplanet++;
            std::cerr << "  Nukleus " << n + 1 << " hyväksytty planeetaksi (Massa: "
                      << std::fixed << std::setprecision(4) << new_planet.mass * SOLAR_MASS_TO_EARTH_MASS << " M_Earth, Akseli: "
                      << std::fixed << std::setprecision(4) << new_planet.axis << " AU). Planeettoja nyt: " << nplanet << std::endl;
        } else {
            std::cerr << "  Nukleus " << n + 1 << " hylätty (Massa: "
                      << std::fixed << std::setprecision(4) << new_planet.mass * SOLAR_MASS_TO_EARTH_MASS << " M_Earth, Akseli: "
                      << std::fixed << std::setprecision(4) << new_planet.axis << " AU). Syy: Ei täytä vaatimuksia (etäisyys/massa)." << std::endl;
        }

        if (nplanet >= MaxPlanets) {
            std::cerr << "Maksimimäärä planeettoja (" << MaxPlanets << ") saavutettu. Lopetetaan planeettojen muodostus." << std::endl;
            break;
        }
    }

    std::sort(planets.begin(), planets.end());



	int planetnum=planets.size();

//  printf(" %i ", planetnum);
  
 //   printf("%lf %lf ", planets[4].axis, planets[4].mass*333000);
//exit(-1);
//	char romans1[50][128]={"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XI", "XII", "XIII", "XIV", "XV", "XVI"};
//	char alphabets1[50][128]={"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o"};
//	printf(" %s", romans1[2]);
	
//	exit(-1);
	
   generate_and_render_povray(planets, planetnum);

	export_planets_to_csv(planets,(std::string) "planets.csv");

    std::cout << "\n--- Planeet system ---\n";
    std::cout << "Seed of simulation: " << seed << "\n";
    std::cout << "Generated planets: " << planets.size() << "\n";
    std::cout << "------------------------------------------\n\n";

    std::cout << std::left << std::setw(5) << "Planet # "
              << std::left << std::setw(12) << "Distance (AU) "
              << std::left << std::setw(12) << "Mass (M_Earth) "
              << std::left << std::setw(12) << "Type "
              << std::left << std::setw(12) << "Radius (km) "
              << std::left << std::setw(12) << "Temperature (K) "
              << std::left << std::setw(10) << "Fe (%) "
              << std::left << std::setw(10) << "Rock (%) "
              << std::left << std::setw(10) << "Ice (%) "
              << std::left << std::setw(10) << "Gas (%) "
              << std::endl;
    std::cout << std::string(110, '-') << std::endl;
    for (size_t i = 0; i < planets.size(); ++i) {
        const Nucleus &p = planets[i];

        if ((p.mass * SOLAR_MASS_TO_EARTH_MASS) < 0.03) {
            continue;
        }

        double total_solid_mass_planet = p.accumulated_iron_mass + p.accumulated_rock_mass + p.accumulated_ice_mass;
        double gas_mass_planet = p.mass - total_solid_mass_planet;
        if (gas_mass_planet < 0) gas_mass_planet = 0;

        double total_mass_for_fractions_planet = total_solid_mass_planet + gas_mass_planet;
        if (total_mass_for_fractions_planet == 0) total_mass_for_fractions_planet = p.mass;

        std::cout << std::left << std::setw(5) << (i + 1)
                  << std::left << std::setw(12) << std::fixed << std::setprecision(3) << p.axis
                  << std::left << std::setw(15) << std::fixed << std::setprecision(3) << (p.mass * SOLAR_MASS_TO_EARTH_MASS);

        std::string planet_type_str;
        if (p.type == Nucleus::gas_p) planet_type_str = "gas";
        else if (p.type == Nucleus::ice_p) planet_type_str = "ice";
        else if (p.type == Nucleus::iron_p) planet_type_str = "fe";
        else planet_type_str = "rock";
        std::cout << std::left << std::setw(12) << planet_type_str
                  << std::left << std::setw(12) << std::fixed << std::setprecision(1) << p.get_radius_km()
                  << std::left << std::setw(12) << std::fixed << std::setprecision(1) << p.get_temperature_k();

        if (total_mass_for_fractions_planet > 0) {
            std::cout << std::left << std::setw(10) << std::fixed << std::setprecision(1) << (p.accumulated_iron_mass / total_mass_for_fractions_planet * 100.0)
                      << std::left << std::setw(10) << std::fixed << std::setprecision(1) << (p.accumulated_rock_mass / total_mass_for_fractions_planet * 100.0)
                      << std::left << std::setw(10) << std::fixed << std::setprecision(1) << (p.accumulated_ice_mass / total_mass_for_fractions_planet * 100.0)
                      << std::left << std::setw(10) << std::fixed << std::setprecision(1) << (gas_mass_planet / total_mass_for_fractions_planet * 100.0);
        } else {
            std::cout << std::left << std::setw(10) << "N/A"
                      << std::left << std::setw(10) << "N/A"
                      << std::left << std::setw(10) << "N/A"
                      << std::left << std::setw(10) << "N/A";
        }
        std::cout << std::endl;

        double planet_mass_earth = p.mass * SOLAR_MASS_TO_EARTH_MASS;
        if (planet_mass_earth > 0.03) {
            double rendering_radius = 0.05 * std::pow(planet_mass_earth, 1.0 / 6.0);
            if (rendering_radius <= 0) rendering_radius = 0.01;

            double log10_axis = (p.axis > 0) ? std::log10(p.axis) : std::log10(1e-10);

            if (Ps != NULL) {
                *Ps << "newpath " << log10_axis << ' ' << 0 << ' ' << rendering_radius << " 0 360 arc ";
                if (p.type == Nucleus::gas_p) {
                    *Ps << "0 0.5 1 setrgbcolor fill\n";
                } else if (p.type == Nucleus::ice_p) {
                    *Ps << "0.7 0.9 1 setrgbcolor fill\n";
                } else if (p.type == Nucleus::iron_p) {
                    *Ps << "0.6 0.3 0 setrgbcolor fill\n";
                } else {
                    *Ps << "0.5 0.5 0.5 setrgbcolor fill\n";
                }
            }
        }

        if (dumpflag) {
            p.dump(static_cast<int>(i + 1), std::cout, params);
            std::cout << std::endl;
        }
    }
    std::cout << std::string(110, '-') << std::endl;

    ps_showpage();
    ps_end();



    delete params;

    return 0;
}



int generate_and_render_povray(std::vector<Nucleus> planets, int planetnum) {

    int counter=0;
    int n=0;

	double massaa=0.0;
    double star_radius_pov = 1.5; // POV-Ray units
    double planet_radius_scale =0.5; // 1 Km in simulation = 0.0005 POV-Ray units
    double distance_scale = 3.0;        // 1 AU in simulation = 5.0 POV-Ray units

    double planet_x =0;  
	double planet_y = 0.0; // Assume they are on the x-z plane (y=0)
    double planet_z = 0.0;    
    double planet_pov_radius = 0;

    FILE *fp_pov=NULL;
    char pov_filename[] = "system_scene.pov";
    char output_image_filename[] = "system_render.png";
	char buffu [256]="";
    char command[512]="";
     char charnames[256]="";
	char romans1[50][128]={"I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XI", "XII", "XIII", "XIV", "XV", "XVI"};

    
    
    memset(buffu,0,256);
    memset(command,0,256);
 

    memset(charnames, 0, 256);
 

    // Luodaan char-taulukko löytöjärjestyksessä oleville nimille
    // Oletetaan, että charnames on riittävän suuri (esim. 256)
    // --- Ensimmäinen lajittelu: Massan mukaan laskevasti (massiivisin ensin) ---
    // Tämä järjestää planets-vektorin niin, että massiivisin planeetta on ensimmäisenä.
  //  std::sort(planets.begin(), planets.end(), [](const Nucleus& a, const Nucleus& b) {
  //      return a.mass > b.mass; // Lajittelee massan mukaan, suurin ensin
  //  });
        // --- Massalajittelu ---
   /*
   
    std::sort(planets.begin(), planets.end(), [](Nucleus& a, Nucleus& b) {
        return a.mass > b.mass; // Lajittelee massan mukaan, suurin ensin
    });

    */
    
    
    int i = 0;
    for (const auto& p : planets) {
        if (i < 256) {
            charnames[i] = static_cast<char>('a' + i);
            std::cout << "Planeetta: " << p.axis << ", Massa: " << p.mass
                      << ", Löytöjärjestys: " << charnames[i] << std::endl;
            i++;
        } else {
            std::cerr << "Varoitus: charnames-taulukko täynnä, kaikkia planeettoja ei voitu käsitellä." << std::endl;
            break;
        }
    }
 
 
 printf("\n>%s<", charnames);
/*
  std::sort(planets.begin(), planets.end(), [](const Nucleus& a, const Nucleus& b) {
        return a.axis > b.axis; // Lajittelee axis-arvon mukaan, suurin ensin
    });
*/

// exit(-1);

    //planet_pointer node1;


 //printf("%lf %lf ", planets[4].axis, planets[4].mass*333000);
 
 //exit(-1);

    // --- POV-Ray File Generation ---
    fp_pov = fopen(pov_filename, "w");
    if (fp_pov == NULL) {
        perror("Error: Unable to open POV-Ray file");
        return -1;
    }

    // POV-Ray header and basic scene setup
    fprintf(fp_pov, "// POV-Ray Scene file generated by System Simulator\n");
    fprintf(fp_pov, "// Generated on: %s\n", __DATE__);
    fprintf(fp_pov, "\n");
    fprintf(fp_pov, "#include \"colors.inc\"\n"); // Standard POV-Ray color definitions
    fprintf(fp_pov, "#include \"functions.inc\"\n");
    // Camera setup (adjust position and look_at for best view)
    // Positioned to look down slightly on the orbital plane
    fprintf(fp_pov, "camera {\n");
    fprintf(fp_pov, "  location <14, 0, -100> // x, y (up), z (depth)\n");
    fprintf(fp_pov, "  look_at <14, 0, 0>\n");
    fprintf(fp_pov, "  right x * image_width / image_height\n");
    fprintf(fp_pov, "  angle 20 // Field of view\n");
    fprintf(fp_pov, "}\n\n");

    // Light source - a strong white light to illuminate the scene
    // This could also be your central star itself
    fprintf(fp_pov, "light_source { <0, 0, -10> color White * 1.5 }\n\n");


/*
    fprintf(fp_pov, "sphere {\n");
    fprintf(fp_pov, "  <0, 0, 0>, %f // Center and radius\n", star_radius_pov);
    fprintf(fp_pov, "  pigment {  ");
  	fprintf(fp_pov, " crackle scale 0.3 \n");	
	fprintf(fp_pov, " color_map { \n");			
			
	fprintf(fp_pov, "  [0.0 rgb <1,1,0.5>*0.8 ] \n");
	fprintf(fp_pov, "  [1.0 rgb <1,1,0.75> ] \n");	
				
	fprintf(fp_pov, "  }// ... color map \n");	  
     fprintf(fp_pov, "   } // Yellow, semi-transparent for glow effect\n");
    fprintf(fp_pov, "  // Add an emissive finish for a glowing effect\n");
    fprintf(fp_pov, "  finish { ambient 1 diffuse 0 emission 1 }\n");
    fprintf(fp_pov, "  // Add a halo for a more realistic star glow (requires photons in render settings)\n");
    // fprintf(fp_pov, "  photons { emission 1 }\n"); // Uncomment if you want to use photons for glow
    fprintf(fp_pov, "}\n\n");
*/



fprintf(fp_pov, "background { color rgb <0, 0, 0.0> } // Tumma avaruus \n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "#declare bright_star= union {\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "light_source {\n");
fprintf(fp_pov, "    <0, 0, 0>\n");
fprintf(fp_pov, "    color rgb <1, 1, 1>\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "looks_like \n");
fprintf(fp_pov, "{\n");
fprintf(fp_pov, "// Tähti (kirkas pallo)\n");
fprintf(fp_pov, "sphere {\n");
fprintf(fp_pov, "    <0, 0, 0>, 1\n ");
fprintf(fp_pov, "    texture {\n");
fprintf(fp_pov, "        pigment { color rgb <1, 1, 0.8> } // Keltainen/valkoinen sävy \n");
fprintf(fp_pov, "        finish {\n");
fprintf(fp_pov, "            emission 1\n");
fprintf(fp_pov, "            diffuse 0.2\n");
fprintf(fp_pov, "            specular 0.5\n");
fprintf(fp_pov, "            roughness 0.01\n");
fprintf(fp_pov, "        }\n");
fprintf(fp_pov, "    }\n");
fprintf(fp_pov, "}\n");
fprintf(fp_pov, "}\n");
fprintf(fp_pov, "}\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "// Mediheheku (säteilevä vaikutus tähden ympärillä)\n ");
fprintf(fp_pov, "sphere {\n");
fprintf(fp_pov, "    <0, 0, 0>, 1 // Hieman suurempi pallo hehkua varten\n");
fprintf(fp_pov, "    hollow\n");
fprintf(fp_pov, "    material {\n");
fprintf(fp_pov, "        texture {\n");
fprintf(fp_pov, "                  pigment { color rgbt <1, 1, 1, 1> } // Lähes läpinäkyvä  \n ");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "        }\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "        interior {\n");
fprintf(fp_pov, "            media {\n");
fprintf(fp_pov, "                scattering { 3, rgb <1, 1, 0.8> * 1/20 }\n");
fprintf(fp_pov, "                //emission 1/10000\n");
fprintf(fp_pov, "                density { spherical poly_wave 2 density_map {\n");
fprintf(fp_pov, "                    [0 color rgbt <0, 0, 0,1>]\n");
fprintf(fp_pov, "                    [1 color rgbt <1, 1, 1,0>]\n");
fprintf(fp_pov, "                }}\n");
fprintf(fp_pov, "                samples 10,20\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "            }\n");
fprintf(fp_pov, "        }\n");
fprintf(fp_pov, "    }\n");
fprintf(fp_pov, "scale 3\n");
fprintf(fp_pov, "}\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "}\n");
fprintf(fp_pov, "\n");
fprintf(fp_pov, "object {bright_star translate x*-2}\n");
fprintf(fp_pov, "\n");



    // --- Planets ---
    // Scale factor for planet radii and distances in POV-Ray units
    // Adjust these to make planets visible and to fit within the view


  //  node1 = planet_head;
    counter = 0;
    for (n=0; n<planetnum; n++) {
		
	//	counter=n;
		
        // Calculate planet position in POV-Ray coordinates
        // Simple circular orbit for visualization. For eccentricity, you'd need more complex math.
        // We'll place them along the X-axis for simplicity in this example.
        //double planet_x =2+ log(node1->a) * distance_scale;
		planet_x =3+counter*distance_scale;  
		planet_y = 0.0; // Assume they are on the x-z plane (y=0)
        planet_z = 0.0;

        //double planet_pov_radius = sqrt(node1->radius) * planet_radius_scale;
        
        massaa=planets[n].mass*333000;
        
		if(massaa>0.02)
		{
        
        planet_pov_radius = pow((double)massaa, 0.15) * planet_radius_scale;
		printf("\n %i %lf %lf ", counter, massaa,  planet_pov_radius);

 
      
  //    exit(-1);
    


		memset(buffu,0,256);

        fprintf(fp_pov, "// Planet #%d\n", counter);
        fprintf(fp_pov, "sphere {\n");
        fprintf(fp_pov, "  <%f, %f, %f>, %f // Position and radius\n",
                planet_x, planet_y, planet_z, planet_pov_radius);


     if(planets[n].type == Nucleus::PlanetType::gas_p)
		  {
			fprintf(fp_pov, "  pigment {  \n");
			
					fprintf(fp_pov, " gradient y  sine_wave frequency 1.5 scale 5 warp {turbulence 0.5 } scale 1/5 turbulence 0.1   \n");	
			fprintf(fp_pov, " color_map { \n");			
			
			fprintf(fp_pov, "  [0.0 rgb <0.533333, 0.427451, 0.352941> ] \n");
			fprintf(fp_pov, "  [0.5 rgb <0.917647, 0.631373, 0.454902> ] \n");	
			fprintf(fp_pov, "  [1.0 rgb <0.992157, 0.952941, 0.847059> ] \n");	
				
			fprintf(fp_pov, "  }// ... color map \n");				
			
			fprintf(fp_pov, " } //...pigment \n  ");
	        fprintf(fp_pov, "  finish {  diffuse 0.65 ambient 0 } // Shiny finish\n");
		  //  fprintf(fp_pov, "  normal { wrinkles scale y/10 scale 3 warp {turbulence 0.1} scale 0.1 bump_size 0.1 } // \n");
		  }
      if(planets[n].type == Nucleus::PlanetType::ice_p)
		  {
			fprintf(fp_pov, "  pigment { color rgb <1, 1, 1> } \n");
	        fprintf(fp_pov, "  finish { phong 0.8 ambient 0} // Shiny finish\n");
		        fprintf(fp_pov, "  normal { agate scale 0.1  turbulence 0.2 bump_size -0.4  } // \n");
		  }
           
      if(planets[n].type == Nucleus::PlanetType::rock_p)
		  {
			int earthlike=0;
			double tempera=planets[n].get_temperature_k();
			if( (planets[n].mass*333000)>0.8) {
				if( (planets[n].mass*333000)<2.0) {
		 
		 
					if(tempera>273.0) {
							   if(tempera<320.0){
											earthlike=1;
										}
								}
		 
			
							} 
		 
					}
		 
		//earthlike=1;
		if(earthlike==1)
		{
				fprintf(fp_pov, "  texture {  \n");
				fprintf(fp_pov, "  pigment {  \n");
			
				fprintf(fp_pov, " wrinkles  scale 5 warp {turbulence 0.01 } scale 0.2  scale 0.5  \n");	
				fprintf(fp_pov, " color_map { \n");			
			
				fprintf(fp_pov, "  [0.0 rgb <0,0,1> ] \n");
				fprintf(fp_pov, "  [0.5 rgb <0,0,1> ] \n");
				fprintf(fp_pov, "  [0.5 rgb <0,1,0> ] \n");	
			    fprintf(fp_pov, "  [1.0 rgb <0.796078, 0.545098, 0.345098> ] \n");	
								
				fprintf(fp_pov, "  }// ... color map \n");				
			
				fprintf(fp_pov, " } //...pigment \n  ");
				fprintf(fp_pov, "  finish { phong 0.8 } // Shiny finish\n");
				fprintf(fp_pov, "  normal { wrinkles scale y/10 scale 3 warp {turbulence 0.1} scale 0.1 bump_size 0.1 } // \n");
					fprintf(fp_pov, " } \n");
		
					fprintf(fp_pov, "  texture { // clouds  \n");
				fprintf(fp_pov, "  pigment {  \n");
			
				fprintf(fp_pov, " granite  turbulence 1 \n");	
				fprintf(fp_pov, " color_map { \n");			
			
				fprintf(fp_pov, "  [0.0 rgbt <0,0,0,1> ] \n");
				fprintf(fp_pov, "  [0.3 rgbt <0,0,1,1> ] \n");
				fprintf(fp_pov, "  [0.5 rgbt <1,1,10> ] \n");	
			    fprintf(fp_pov, "  [1.0 rgbt <1,1,1,0> ] \n");	
								
				fprintf(fp_pov, "  }// ... color map \n");				
			
				fprintf(fp_pov, " } //...pigment \n  ");
				fprintf(fp_pov, "  finish { phong 0.8 ambient 0} // Shiny finish\n");
			//	fprintf(fp_pov, "  normal { wrinkles scale y/10 scale 3 warp {turbulence 0.1} scale 0.1 bump_size 0.1 } // \n");
					fprintf(fp_pov, " } \n");	
			}
			else
				{
				fprintf(fp_pov, "  pigment {  \n");
			
				fprintf(fp_pov, " wrinkles  scale 5 warp {turbulence 0.01 } scale 0.2  scale 0.5  \n");	
				fprintf(fp_pov, " color_map { \n");			
			
				fprintf(fp_pov, "  [0.0 rgb <0.266667, 0.25098, 0.203922> ] \n");	
			    fprintf(fp_pov, "  [1.0 rgb <0.87451, 0.729412, 0.541176> ] \n");	
								
				fprintf(fp_pov, "  }// ... color map \n");				
			
				fprintf(fp_pov, " } //...pigment \n  ");				
				
					fprintf(fp_pov, "  finish { diffuse 0.5 roughness 0.05 ambient 0}\n"); // Less shiny, more diffuse
					fprintf(fp_pov, "  normal { granite scale 0.5 turbulence 0.2 bump_size -0.3 } // \n");	
				}
		
		
	//	 }
		 
		 
		 
		  }
           
/*
        if (node1->gas_giant) {
            fprintf(fp_pov, "  pigment { color rgb <0.4, 0.8, 0.6> } // Geenish for gas giants\n");
            fprintf(fp_pov, "  finish { phong 0.8 } // Shiny finish\n");
        } else {
            // Terrestrial planet coloring based on features (simplified)
            if (node1->hydrosphere > 0.5 && node1->cloud_cover < 0.5) {
                fprintf(fp_pov, "  pigment { color rgb <0.2, 0.4, 0.8> } // Water world (blue)\n");
            } else if (node1->ice_cover > 0.3) {
                fprintf(fp_pov, "  pigment { color White } // Ice world\n");
            } else if (node1->surface_pressure > 0.001) { // Basic atmospheric planet
                fprintf(fp_pov, "  pigment { color rgb <0.2, 0.6, 0.2> } // Green/Earth-like\n");
            } else {
                fprintf(fp_pov, "  pigment { color rgb <0.7, 0.4, 0.2> } // Rocky/desert planet\n");
            }
            fprintf(fp_pov, "  finish { phong 0.6 roughness 0.05 ambient 0}\n"); // Less shiny, more diffuse
        }
        */
        fprintf(fp_pov, "}\n\n");


       fprintf(fp_pov, "text { \n\n");
      fprintf(fp_pov, "ttf \"timrom.ttf\" ");

    memset(buffu,0,256);
  sprintf( buffu ,"\"%d\" ", counter+1 );
 // fprintf(fp_pov,"%s", buffu);
 fprintf(fp_pov," \"%s\" ", romans1[counter]);
 
  fprintf(fp_pov, " 0.15,0 \n pigment {color rgb <1,1,1> }\n");
      fprintf(fp_pov, " translate y*-2.25 \n\n");
    fprintf(fp_pov, " translate x*%f \n\n", planet_x-0.25); // planet_x
      fprintf(fp_pov, " }\n\n");

counter++;
} // ... mass over limit
else {
//	n--;
	}

    //    counter++;
    //    node1 = node1->next_planet;
    }



    fclose(fp_pov);

//	 exit(-1);
    printf("POV-Ray scene file '%s' generated successfully.\n", pov_filename);

    // --- Call POV-Ray to Render ---

    memset(command,0,512);
    // Adjust the path to povray.exe if it's not in your system's PATH
    // On Linux/macOS: "povray -W1200 -H400 +A +R -O%s %s"
    // On Windows: "pvengine /RENDER /EXIT /W1200 /H400 /O%s %s"
    // I'll use a generic command that works on many systems if 'povray' is in PATH.
    // Make sure to replace `povray` with `pvengine` if you're on Windows and it's not in your PATH.
    sprintf(command, "/usr/bin/povray -W1800 -H300 -O%s %s", output_image_filename, pov_filename);

    printf("Executing POV-Ray command: %s\n", command);
    int result = system(command);

    if (result == 0) {
        printf("POV-Ray rendering successful. Image saved to '%s'\n", output_image_filename);
    } else {
        fprintf(stderr, "Error: POV-Ray rendering failed with code %d. Make sure POV-Ray is installed and in your system's PATH.\n", result);
        fprintf(stderr, "You might need to adjust the 'povray' command to 'pvengine' or specify its full path.\n");
    }

 return(0);
}



void export_planets_to_csv(std::vector<Nucleus> planets,  std::string filename) {

	Nucleus planet;
	
    int n=0;
	int planetnum=planets.size();

    std::ofstream csv_file(filename);


    std::vector<Nucleus> sorted_planets = planets;

    std::sort( sorted_planets.begin(), sorted_planets.end() );

    if (!csv_file.is_open()) {
        std::cerr << "Error: Could not open CSV file " << filename << " for writing.\n";
        return;
    }

    // Write CSV header
    csv_file << "radius_au,eccentricity,mass_mearths,planet_type,p_iron,p_rock,p_ice,p_gas, radius_km,teq_k\n";
	
	printf("\n %i ", planetnum);
	
	for (n=0;n<planetnum;n++)
	{
	 planet=planets[n];
	
		  //      csv_file << std::fixed << std::setprecision(6)
           //      << planet.axis << ","
            //     << planet.eccen << ",";
     
             double mass_earth_masses = planet.mass * SOLAR_MASS_TO_EARTH_MASS;

        double total_solid_mass = planet.accumulated_iron_mass + planet.accumulated_rock_mass + planet.accumulated_ice_mass;
        double gas_mass = planet.mass - total_solid_mass;
        if (gas_mass < 0) gas_mass = 0; // Ensure non-negative gas mass

        double total_mass_for_fractions = total_solid_mass + gas_mass;

        double p_iron = 0.0;
        double p_rock = 0.0;
        double p_ice = 0.0;
        double p_gas = 0.0;
		double planet_radius_km=planets[n].get_radius_km();
		double planet_equilibrium_temperature=planets[n].get_temperature_k();

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
            case Nucleus::gas_p:  planet_type_str = "Giant"; break;
        }

        csv_file << std::fixed << std::setprecision(6)
                 << planet.axis << ","
                 << planet.eccen << ","
                 << mass_earth_masses << ","
                 << planet_type_str << ","
                 << p_iron << ","
                 << p_rock << ","
                 << p_ice << ","
                 << p_gas << ","
				 << std::setprecision(2)
				 <<	round(planet_radius_km)<<","
				 <<	round(planet_equilibrium_temperature)<< "\n";
     
     
		}
	
	

    csv_file.close();
    std::cout << "Planet data exported to " << filename << std::endl;
	

	




}
