
/*%CC acrete.c -o acrete -lm

* This will also work with g++ 1.40.3: g++ acrete.c -o acrete -lm
*
* Dole accretion model for planetary formation.
* Adapted from Wales Larrison's BASIC code.
* References:
*
* Dole, S.H. "Formation of Planetary Systems by Aggregation:
* a Computer Simulation" Icarus, vol 13, p 494, 1970
* Isaacman, R. & Sagan, C. "Computer Simulation of Planetary
* Accretion Dynamics: Sensitivity to Initial Conditions"
* Icarus, vol 31, p 510, 1977
*
* Usage:
* acrete [-seed #] [-dump]
*
* -seed specifies initial value to random number generator (uses time otherwise)
* -dump causes a dump of the generated system on stdout
* Produces a PostScript file "planets.ps"
*
* Jon Leech (leech@cs.unc.edu)
*/

#include <iostream>
#include <fstream>   // For std::ofstream
#include <cstring>   // For strcmp (C-style string functions)
#include <cmath>     // For math functions like pow, log10, floor, sqrt, fabs, cbrt
#include <cstdlib>   // For atoi, qsort, srand48, drand48, lrand48, time
#include <cstdio>    // For sprintf (needed if using sprintf, though std::to_string or string streams are preferred in C++)


// The original code included <rand48.h> which is non-standard.
// These functions are usually declared in <cstdlib> or <stdlib.h>
// for POSIX systems. We'll assume they are available.
extern "C" {
    double drand48();
    long   lrand48();
    void   srand48(long);
};

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
    W;

    DoleParams();
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
    return (b > c) ? b : c; // Corrected this line
}

std::ofstream *Ps;
int PsPage = 1;
double
    PsBase = 50,             // base offset from lower left corner
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
    PsYscale = PsXscale;         // could be width / yspan
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
    enum {
    Rock, GasGiant
    }   type;          // type of planet

    Nucleus();
    void dump(int, std::ostream &, DoleParams *);
    friend int nucleusCompare(const void *, const void *);
    double lowbound(DoleParams *params);
    double highbound(DoleParams *params);
};

Nucleus::Nucleus() {
    axis = eccen = mass = 0;
    pRad = aRad = 0;
    pAttr = aAttr = 0;
    type = Rock;
}

double Nucleus::lowbound(DoleParams *params) {
    return params->lowbound(pRad, pAttr);
}

double Nucleus::highbound(DoleParams *params) {
    return params->highbound(aRad, aAttr);
}

// Comparison function used by qsort()
int nucleusCompare(const void *p1, const void *p2) {
    double  r1 = static_cast<const Nucleus *>(p1)->axis,
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
    o << "\tType\t\t" << ((type == GasGiant) ? "Gas giant" : "Rocky") << '\n';

    o << "\tRange          = [" << lowrange << "," << highrange << "]\n";
    o << "\tX[pa]limit     = [" << xplimit << "," << xalimit << "]\n";
    o << "\tX{peri,apo}gee = [" << pAttr << "," << aAttr << "]\n";
    o << "\tR{peri,apo}gee = [" << pRad << "," << aRad << "]\n"; // Corrected this line, was 'pRad' twice

    o << "\tCritical Mass  = " << massCrit << ' '
      << '(' << mass / massCrit << "% of Critical)\n";
}

int main(int ac, char *av[])
{
    int dumpflag = 0;
    long seed = time(0);

    for (int i = 0; i < ac; i++) {
    if (strcmp(av[i], "-dump") == 0) // Use strcmp directly, not std::strcmp
        dumpflag = 1;
    else if (strcmp(av[i], "-seed") == 0)
        seed = std::atoi(av[++i]);
    }

    double  nucleieccent, nucleiradius, masslast, eccent, mass, rperigee,
            rapogee, radius, mu, xperigee, xapogee, masscritical,
            density, bw, volume, da, dp, masspass;

    int      lowband, highband, dustcheck, iterate, hit;

    // Bands are measured in 1/5 au intervals from 0 to 50 AU
    const int
    MaxBand = 50 * 5,
    MaxNuclei = 1000;   // Max # of nuclei to inject

    enum { Mixed, Gas, Empty }
        band[MaxBand+1];    // Contents of band

    int MaxPlanets = 20;    // Size of dynamic array containing generated planets
    Nucleus *nuclei = new Nucleus[MaxPlanets];
    DoleParams *params = new DoleParams;

    // Declare nplanet here and initialize it to 0
    int nplanet = 0; // <--- ADDED THIS LINE

    randomize(seed);

    ps_window(-1, -1, 2, 1);
    ps_begin("planets.ps");

    // Set up primordial cloud
    // Dust within .3 au is swept out by radiation pressure

    for (int n_band = 1; n_band <= MaxBand; n_band++) // Renamed n to n_band to avoid conflict
    band[n_band] = Mixed;

    // Draw scale on figure
    logscale("AU");

    // Plot density function; space samples even on log scale
    double max = 1.1 * params->density(0.3), lastau, lastrho,
            logaumin = std::log10(0.3),
            logaumax = std::log10(50),
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

    // This is only used on first pass of nuclei
    masslast = 0;
    iterate = 0;

    // Mass accumulation loop - continue until minimal accumulation
    while (1) {
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
            // std::cout << "no dust; " << nuclei_count << "is a dud" << std::endl;
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

        for (int i = 0; i < nplanet; i++) { // nplanet used here
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
            // break; // Original code had a break here, but it would only check against the first planet.
                     // Removing it allows checking against all existing planets. If this causes
                     // different simulation results, you might need to re-evaluate the original intent.
        }

        masspass += bandvol * density;
        }

        body.mass = mass = masspass;

        // std::cout << mass << ' ' << nuclei_count << std::endl;

        // check for mass growth convergence
        if (mass == 0)
        mass = 1e-15;
        if (mass >= masscritical) {
        body.type = Nucleus::GasGiant;
        }
        if (std::fabs(masslast / mass - 1) < 0.01)
        break;
        masslast = mass;

        iterate = 1;
    }    // end mass accumulation loop

    // Clear out bands emptied of dust
    if (dustcheck == 1) {
        std::cout << "event " << nuclei_count << " completed mass growth; swept from "
                  << dp << " to " << da;
        if (mass > masscritical)
        std::cout << "(gas giant)";
        std::cout << std::endl;
    }
    if (lowband == highband)
        highband++;

    for (int bandno = lowband; bandno <= highband; bandno++) {
        if (mass < masscritical) {
        if (band[bandno] == Mixed)
            band[bandno] = Gas;
        }

        if (mass >= masscritical) {
        // Clear out bands emptied by gas
        body.type = Nucleus::GasGiant;
        if (band[bandno] == Gas || band[bandno] == Mixed)
            band[bandno] = Empty;
        }
    }

    // std::cout << "check nucleus " << nuclei_count << " for overlap and coalescence\n";

    // check for orbital and  gravitational overlap and coalescence
    // std::cout << body.pRad - body.pAttr
    //     << body.aRad+body.aAttr << std::endl;

    if (body.axis == 0)
        continue;

    hit = 0;
    for (int i = 0; i < nplanet; i++) { // nplanet used here
        double  newradius, newmass, neweccent, term1, term2, term3, munew;

        // std::cout << nuclei_count << '\t' << i << '\t'
        //          << nuclei[i].aRad << '\t' << rperigee-xperigee << '\t'
        //          << nuclei[i].pRad << '\t' << rapogee+xapogee << std::endl;

        if ((nuclei[i].aRad < (rperigee - xperigee) &&
           (nuclei[i].aRad + nuclei[i].aAttr) < rperigee) ||
           ((rapogee + xapogee) < nuclei[i].pRad &&
           rapogee < (nuclei[i].pRad - nuclei[i].pAttr)))
        continue;

        std::cout << "coalesence of nuclei  " << nuclei_count << ' ' << i << std::endl;

        hit = 1;
        newradius = (body.mass + nuclei[i].mass) /
            (body.mass / body.axis + nuclei[i].mass / nuclei[i].axis);

        term1 = body.mass * std::sqrt(body.axis) * std::sqrt(1 - sqr(body.eccen));
        term2 = nuclei[i].mass * std::sqrt(nuclei[i].axis) *
                      std::sqrt(1 - sqr(nuclei[i].eccen));
        term3 = (body.mass + nuclei[i].mass) * std::sqrt(newradius);

        neweccent = sqr(std::abs(1 - std::pow((term1 + term2) / term3, 2)));
        newmass = body.mass + nuclei[i].mass;

        // std::cerr << "Nuking body " << i << std::endl;
        nuclei[i].axis = 0;
        nuclei[i].eccen = 0;
        nuclei[i].mass = 0;
        body.axis = newradius;
        body.eccen = neweccent;
        body.mass = newmass;
        body.pRad = newradius * (1 - neweccent);
        body.aRad = newradius * (1 + neweccent);
        munew = std::pow(newmass / (1 + newmass), .25);
        body.pAttr = body.pRad * munew;
        body.aAttr = body.aRad * munew;

        //    std::cout << "new mass = " << newmass << " new radius = " << newradius
        //              << " new eccentricity = " << neweccent << std::endl;
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

    // iterate for mass captured
    //  std::cout << "Start of iteration for mass" << N
    //        << " mass = " << mass
    //        << " masslast = " << masslast << std::endl;

    if (body.mass >= masscritical)
        body.type = Nucleus::GasGiant;

    // Add new planet
    if (nplanet >= MaxPlanets) { // nplanet used here
        Nucleus *newplanet = new Nucleus[MaxPlanets*2];

        // Copy old planets to new
        for (int i = 0; i < nplanet; i++) // nplanet used here
        newplanet[i] = nuclei[i];

        delete [] nuclei;
        nuclei = newplanet;
        MaxPlanets *= 2;
    }
    nuclei[nplanet++] = body; // nplanet used here

    // Sweep planet array, removing merged bodies
    for (int i = 0; i < nplanet; i++) { // nplanet used here
        if (nuclei[i].axis == 0) {
        std::cerr << "Nuking body " << i << "; " << nplanet << " remaining" << std::endl;
          nuclei[i] = nuclei[nplanet-1];
          nuclei[nplanet-1].axis = 0;
          nplanet--;
        }
    }
    } // End of the main loop that iterates MaxNuclei times

    std::cout << " all bands taken....." << nuclei_count << "    nuclei used      " << std::endl; // Use nuclei_count
    ps_text(2.2, 0.9, "emptied bands");

    // Sort nuclei by radius
    std::qsort(static_cast<void *>(nuclei), nplanet, sizeof(Nucleus), nucleusCompare); // nplanet used here

    // Draw planets, gas giants as filled circles
    double massSum = 0;
    for (int i = 0; i < nplanet; i++) { // nplanet used here
    massSum += nuclei[i].mass;

    double au = std::log10(nuclei[i].axis);
    double r = std::cbrt(nuclei[i].mass);

    ps_circle(au, 0, r, nuclei[i].type == Nucleus::GasGiant);
    }
    ps_showpage();

    std::cout << "Total mass = " << massSum << std::endl;

    if (dumpflag) {
    std::cout << "Random number seed =" << seed << std::endl;

    for (int i = 0; i < nplanet; i++) { // nplanet used here
        if (nuclei[i].axis > 0)
        nuclei[i].dump(i, std::cout, params);
    }

    std::cout << "Bands with dust still in them:\n";
    for (int i = 1; i <= MaxBand; i++)
        if (band[i] == Mixed)
        std::cout << i << ' ';
    std::cout << std::endl;
    }

    ps_end();

    return (0);
}
