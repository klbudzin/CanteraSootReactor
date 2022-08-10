//! @file IdealGasReactor.h

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#ifndef CT_IDEALGASREACTORSOOT_H
#define CT_IDEALGASREACTORSOOT_H

#include "Reactor.h"

namespace Cantera
{

/**
 * Class IdealGasReactorSoot is a class for stirred reactors that is specifically
 * optimized for ideal gases. In this formulation, temperature replaces the
 * total internal energy as a state variable.
 */
class IdealGasReactorSoot : public Reactor
{
public:
    IdealGasReactorSoot(double rhoTotal,double NDDInitial = 0., double m_MassCarbonInitial = 0.,bool output = false);

    virtual std::string type() const {
        return "IdealGasReactorSoot";
    }

    virtual void setThermoMgr(ThermoPhase& thermo);

    virtual void getState(doublereal* y);

    virtual void initialize(doublereal t0 = 0.0);

    virtual void eval(double t, double* LHS, double* RHS);

    virtual void updateState(doublereal* y);

    //! Return the index in the solution vector for this reactor of the
    //! component named *nm*. Possible values for *nm* are "mass",
    //! "volume", "temperature", the name of a homogeneous phase species, or the
    //! name of a surface species.
    virtual size_t componentIndex(const std::string& nm) const;
    std::string componentName(size_t k);
    doublereal GetCarbonMassFrac();
    doublereal GetNd();
protected:
    //!Calculate the Solid carbons specific heat at constant volume 
    doublereal CalculateSootCv(doublereal T);
    //!Calculate the Solid carbon's specific heat at constant pressure
    doublereal CalculateSootCp_Ru(doublereal T);
    //!Calculate The Solid Carbon's specific enthalpy
    doublereal CalculateSootH_k_RuT(doublereal T);
    //!Calculate The Solid Carbon's internal energy
    doublereal CalculateSootU_k(doublereal T);
    //!Calculate The total energy of the system
    doublereal CalculateTotalEnergy(doublereal Temperature);
    //!Calculate the Temperature From Total Energy
    void SolveTemperatureAndSetState(doublereal gasDensity);

    vector_fp m_uk; //!< Species molar internal energies
    vector_fp m_Concentration_k;
    double *Yk_gas; //!< Species mass fractions w/ respect to the total mass = mass of solid + mass of gas
    doublereal TotalDensity;
    doublereal TotalEnergy;
    doublereal Y_Carbon;//!< Mass fraction of solid carbon
    doublereal Nd; //!< Soot Number Density
    doublereal const ConstantPressure = OneAtm; //HardCoded Constant Pressure for now
   // doublereal const StartingRho; //Possibility if we want to leave constant pressure term for the simit side of things

    //Soot reaction rates Calculator Methods
    void CalculateNucleationRate(doublereal T, doublereal C2H2Conc);
    void CalculateSurfaceGrowthRate(doublereal T, doublereal C2H2Conc);
    void CalculateAgglomerationRate(doublereal T);
    void CalculateO2OxidationRate(doublereal T, doublereal O2Conc);
    void CalculateOOxidationRate(doublereal T, doublereal OConc);
    void CalculateOHOxidationRate(doublereal T, doublereal OHConc);
    void CalculateSootRates(doublereal T);
    void CalculateO2OxidationRateLeungEtAl(doublereal T, doublereal O2Conc);
    //Soot Reaction Rate variable members
    doublereal NucleationRate = 0;
    doublereal SurfaceGrowthRate = 0;
    doublereal AgglomerationRate = 0;
    doublereal O2OxidationRate = 0;
    doublereal OOxidationRate = 0;
    doublereal OHOxidationRate = 0;
    //Species Indices we would like to know to ease calculation of Soot Reaction Rates set in the initialize step
    int O2Ind;
    int OHInd;
    int OInd;
    int C2H2Ind;
    int COInd;
    int H2Ind;
    int HInd;

    //Soot reaction Rate Constants and parameters
    static doublereal Cmin; // the number of carbon atoms per carbon particles
    static doublereal Ca;
    static doublereal OxidationCollisionEfficiency; // The collision efficiency for the OH and O oxidiation rates
    doublereal SurfGrowthRateConstant;
    doublereal SurfaceAreaConstant; // Parts of the Surface Area term that don't depend on Yc and Nd
    doublereal AgglomerationConstant;
    static doublereal NdNucleationConversionTerm; //2/Cmin*Avagadro's Number -> relates formation of solid carbon moles -> formation of solid carbon atoms -> formation of solid carbon particles (1 soot)
 
    //Soot diameter at current timeStep
    doublereal CurrentSootAvgDiameter;


    //NASA7Coefficients or C(S)
    static vector_fp CSNasa7TLow;
    static vector_fp CSNasa7THigh;
    static double MW_Carbon;
    static double CarbonDensity;
    bool printOutInitialState;
};

}

#endif
