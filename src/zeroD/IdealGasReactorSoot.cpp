//! @file IdealGasReactor.cpp A zero-dimensional reactor

// This file is part of Cantera. See License.txt in the top-level directory or
// at https://cantera.org/license.txt for license and copyright information.

#include "cantera/zeroD/IdealGasReactorSoot.h"
#include "cantera/zeroD/FlowDevice.h"
#include "cantera/zeroD/Wall.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/base/utilities.h"

#include <boost/math/tools/roots.hpp>
namespace bmt = boost::math::tools;
using namespace std;


namespace Cantera
{
    //Thermo-Properties of Solid Carbon
    vector_fp IdealGasReactorSoot::CSNasa7TLow = { -3.108720720e-01, 4.403536860e-03, 1.903941180e-06,-6.385469660e-09, 2.989642480e-12,
       -1.086507940e+02, 1.113829530e+00 };
    vector_fp IdealGasReactorSoot::CSNasa7THigh = { 1.455718290e+00, 1.717022160e-03,-6.975627860e-07, 1.352770320e-10,-9.675906520e-15,
       -6.951388140e+02,-8.525830330e+00 };
    doublereal IdealGasReactorSoot::MW_Carbon = 12.0107;
    doublereal IdealGasReactorSoot::CarbonDensity = 2000;
    doublereal IdealGasReactorSoot::Cmin = 700.; //Average Number of carbon atoms per solid carbon soot particle
    doublereal IdealGasReactorSoot::Ca = 3.; //This is the Agglomeration rate constant currently defined from What other investigators used as said by (Leung and Linsdstedt in 1991 (They used 9) )
    doublereal IdealGasReactorSoot::OxidationCollisionEfficiency = .2;
    doublereal IdealGasReactorSoot::NdNucleationConversionTerm = 2. / Cmin * Avogadro;

    


    //Constructor, Set the Starting Values of NDD and the Carbon Mass Fractions, and the constant total density
    //The gas phase is introduced by the addition of a solution using Reactor.add(Solution)
    IdealGasReactorSoot::IdealGasReactorSoot(double rhoT,double NDDInitial, double Y_MassCarbonInitial,bool Output) : Reactor() {
        this->Nd = NDDInitial;
        this->Y_Carbon = Y_MassCarbonInitial;
        this->TotalDensity = rhoT;
        this->printOutInitialState = Output;
        this->TotalEnergy = 0.;

        //Set The Soot Reaction Rate Constants
        this->SurfGrowthRateConstant = (std::pow(Pi, 1. / 6.) * std::pow(rhoT, .5) *
            std::pow(6. / CarbonDensity, 1. / 3.)) * 6000.;
        this->AgglomerationConstant = 2. * Ca * std::pow(6./ Pi / CarbonDensity, 1. / 6.) * std::pow(6. * Boltzmann / CarbonDensity, 1. / 2.) * rhoT * rhoT;
        this->SurfaceAreaConstant = Pi * std::pow(6. / Pi / CarbonDensity, 2. / 3.) * rhoT;
    }

    //This Method just ensures that the thermoPhase or solution is an idealGas
    void IdealGasReactorSoot::setThermoMgr(ThermoPhase& thermo)
    {
        if (thermo.type() != "IdealGas") {
            throw CanteraError("IdealGasReactorSoot::setThermoMgr",
                               "Incompatible phase type provided");
        }
        Reactor::setThermoMgr(thermo);
    }
    
    //Below Sets the state in Y given the current m_thermo variable and Y_Carbon/NDD states
    void IdealGasReactorSoot::getState(double* y)
    {
        if (m_thermo == 0) {
            throw CanteraError("IdealGasReactor::getState",
                "Error: reactor is empty.");
        }
        m_thermo->restoreState(m_state);
        //If TotalEnergyNot Set, set it
        if (TotalEnergy == 0)
            TotalEnergy = CalculateTotalEnergy(m_thermo->temperature());

        // To Be Consitent with the rest of the reactors, The total mass and volume are the first and second
        // Solution Variable. In the current case they are assumed not to change
        // set the first component to the gaseous Mass (Density defined as (m_gas /Vol))
        //double mass_gas = m_thermo->density() * m_vol;
        m_mass = TotalDensity * m_vol;
        y[0] = m_mass; //Represents the Total Mass
        //if (m_mass < mass_gas)
        //    std::cout << "Something whent wrong and the gaseous mass is larger than the total mass....";
        y[1] = m_vol;  //set the second component to the total volume
        // Set the third component to the temperature
        y[2] = m_thermo->temperature();
        // set components y+3 ... y+K+2 to the mass fractions of each species (this is Y_kgas)
        m_thermo->getMassFractions(Yk_gas);
        //Now we need to scale them by the Carbon MassFraction
        double oneMYc = 1 - Y_Carbon;
        for (size_t i = 0; i < m_nsp; i++) 
            y[i + 3] = Yk_gas[i] * oneMYc;

        //Remaining state space is mass carbon, Ndd
        y[m_nv-2] = Y_Carbon;
        y[m_nv - 1] = Nd;
    }

    void IdealGasReactorSoot::initialize(doublereal t0)
    {
        Reactor::initialize(t0);
        //Add in two more equations for Carbon Mass and Soot Number Density
        m_nv += 2;
        m_uk.resize(m_nsp, 0.0);
        m_Concentration_k.resize(m_nsp, 0.0);
        Yk_gas = new double[m_nsp];
        //Grab the indices location of species of interest
        O2Ind = m_thermo->speciesIndex("O2");
        OHInd = m_thermo->speciesIndex("OH");
        OInd = m_thermo->speciesIndex("O");
        C2H2Ind = m_thermo->speciesIndex("C2H2");
        COInd = m_thermo->speciesIndex("CO");
        H2Ind = m_thermo->speciesIndex("H2");
        HInd = m_thermo->speciesIndex("H");
    }

    void IdealGasReactorSoot::updateState(doublereal* y)
    {
        // The components of y are [0] the total mass, [1] the total volume,
        // [2] the temperature, [3...K+3] are the mass fractions of each species,
        // and [K+3...] are the coverages of surface species on each wall.
        //Mass and volume should not be chaning for me -Kenny
        m_mass = y[0];//Total mass
        m_vol = y[1];//Total Volume

        //Grab the Carbon mass fraction and soot number density
        Y_Carbon = std::max(0.0,y[m_nv - 2]); //Carbon total mass fraction
        Nd = std::max(0.0,y[m_nv - 1]); //Soot number density

        //Need to scale back the mass fractions for the thermoState, Update our Gaseous Mass Fraction Vector
        double one_Yc = 1. / (1. - Y_Carbon);
        for (size_t i = 0; i < m_nsp; i++)
            Yk_gas[i] = y[i + 3] * (one_Yc);
        //Set the State of the thermophase with YTP
        m_thermo->setMassFractions_NoNorm(Yk_gas);
       
        //Now Set the thermoState
        //Option A) Set the State using the constant Pressure value
        //m_thermo->setState_TP(y[2], ConstantPressure);

        //Option B) Calculate New gaseous density from total density and set it
        doublereal gaseousdensity = TotalDensity * (1. - Y_Carbon) / (1 - Y_Carbon * TotalDensity / CarbonDensity);
        SolveTemperatureAndSetState(gaseousdensity);
        //m_thermo->setState_TR(y[2], gaseousdensity);

        //Update to connected reactors, currently shouldn't do anything... If anyone tries to update this hopefully
        //They worked on the through reactor equations
        //What it does do though is save the thermo state m_state to the thermo object!
        updateConnected(true);
    }

    void IdealGasReactorSoot::eval(double time, double* LHS, double* RHS)
    {
        //Grab Parts of the ODE Solver we will be changing as a reference to make things easier to read
        double& mcvdTdt = RHS[2]; // m * c_v * dT/dt
        double* mdYdt = RHS + 3; // mass * dY/dt

        //Make Sure the thermo object actually represents the saved state
        m_thermo->restoreState(m_state);
        //Need Temperature and Concentrations for Soot Reactions
        doublereal T = m_thermo->temperature();
        m_thermo->getConcentrations(&m_Concentration_k[0]);
        //Molecular Weights needed for Species Production Source
        const vector_fp& mw = m_thermo->molecularWeights();

        //Energy Equation Requires knowing the partial molar internal energies
       // m_thermo->getPartialMolarIntEnergies(&m_uk[0]);
        //CvTotal = Cv_g*(1-Yc) + Cvc*Yc
       // double CvTot = m_thermo->cv_mass() * (1 - Y_Carbon) + CalculateSootCv(T) * Y_Carbon;
        //double Uc = CalculateSootU_k(T);

        //Print Debug Statements for log
        /*
        if(!printOutInitialState) {
            doublereal e_s = GasConstant / MW_Carbon * (T * CalculateSootH_k_RuT(T)) - OneAtm / CarbonDensity;//sensible + chemical
            doublereal e_g = m_thermo->intEnergy_mass();// sensible + chemical
            std::cout << "Temperature = " << T << std::endl;
            std::cout << "Solid soot total energy: " << e_s << std::endl;
            std::cout << "Gaseous total energy: " << e_g << std::endl;
            std::cout << " ------------------------------------> TOTAL ENERGY : __ : = " << CalculateTotalEnergy(T) << std::endl;
            //std::cout << "Solid Simit sensible energy = H_C,sens(J/Kg) -P/rho_s = " << e_s - GasConstant / MW_Carbon * 298.15 * CalculateSootH_k_RuT(298.15) << std::endl;
            for (size_t n = 0; n < m_nsp; n++) {
                e_g -= Yk_gas[n] * m_thermo->Hf298SS(n) / mw[n];
            }
            //std::cout << "Gaseous Simit sensible energy = ct->intEnergy_mass() = " << e_g << std::endl;
        }
        if (printOutInitialState) {
            doublereal e_s = GasConstant/MW_Carbon * ( T * CalculateSootH_k_RuT(T) ) - OneAtm / CarbonDensity;//sensible + chemical
            doublereal e_g = m_thermo->intEnergy_mass();// sensible + chemical
            std::cout << "_____________________________________________________________" << std::endl;
            std::cout << " new Time Step" << std::endl;
            std::cout << "_____________________________________________________________" << std::endl;
            std::cout << "INITIALLY:::::" << std::endl;
            std::cout << "Temperature Initially = " << T << std::endl;
            std::cout << "Total Density = " << TotalDensity << std::endl;
            std::cout << "Solid soot total energy: " << e_s << std::endl;
            std::cout << "Gaseous total energy: " << e_g << std::endl;
            //std::cout << "Soot Cv = "<< CalculateSootCv(T) << std::endl;
            std::cout << "Total Cv = " << CvTot << std::endl;
            //std::cout << "Soot U_k = " << Uc << std::endl;
            std::cout << "Volume Total = " << m_vol << std::endl;
            //std::cout << "Volume gas = " << m_vol * (1 - TotalDensity / CarbonDensity * Y_Carbon) << std::endl;
          //  std::cout << "Solid Simit sensible energy = H_C,sens(J/Kg) -P/rho_s = " << e_s - GasConstant/MW_Carbon*298.15*CalculateSootH_k_RuT(298.15)<< std::endl;
          //  for (size_t n = 0; n < m_nsp; n++) {
           //     e_g -= Yk_gas[n] * m_thermo->Hf298SS(n) / mw[n];
           // }
           // std::cout << "Gaseous Simit sensible energy = ct->intEnergy_mass() = " << e_g << std::endl;
            printOutInitialState = false;
        }
        */

        if (m_chem) { //If chemistry is enabled, grab/solve for the species production terms (kmol/m^3/s)
            m_kin->getNetProductionRates(&m_wdot[0]); // "omega dot"
        }
        
        //Now add in all the Normal Reaction Species Terms
        for (size_t n = 0; n < m_nsp; n++) {
            // production in gas phase and from surfaces
            mdYdt[n] = (m_wdot[n] * m_vol) * mw[n] * (1 - TotalDensity / CarbonDensity * Y_Carbon);
            //Assign left-hand side of dYdt ODE as total mass
            LHS[n + 3] = m_mass;

            // heat release from gas phase and surface reactions
           // mcvdTdt -= m_wdot[n] * m_uk[n] * m_vol * (1 - TotalDensity / CarbonDensity * Y_Carbon) * factor;
        }
        
        //Now we Compute and add in the Soot Terms
        CalculateSootRates(T);
        if (printOutInitialState) {
            std::cout << "Current Temperature is            : " << T << std::endl;
            std::cout << "The carbon mass fraction is       : " << Y_Carbon << std::endl;
            std::cout << "The Actual Soot Number density is : " << Nd << std::endl;
            std::cout << "The Nucleation Rate is            : " << NucleationRate << std::endl;
            std::cout << "The Surface Growth Rate is        : " << SurfaceGrowthRate << std::endl;
            std::cout << "The Agglommeration Rate is        : " << AgglomerationRate << std::endl;
            std::cout << "The O2 Oxidation Rate is          : " << O2OxidationRate << std::endl;
            std::cout << "The O Oxidation Rate is           : " << OOxidationRate << std::endl;
            std::cout << "The OH Oxidation Rate is          : " << OHOxidationRate << std::endl;
            std::cout << "The C2H2 Concentration is         : " << m_Concentration_k[C2H2Ind] << std::endl;
            
        }
        mdYdt[C2H2Ind] += m_vol * mw[C2H2Ind] * (-NucleationRate - SurfaceGrowthRate);
        mdYdt[OInd] += m_vol * mw[OInd] * (-OOxidationRate);
        mdYdt[O2Ind] += m_vol * mw[O2Ind] * (-.5 * O2OxidationRate);
        mdYdt[OHInd] += m_vol * mw[OHInd] * (-OHOxidationRate);
        mdYdt[COInd] += m_vol * mw[COInd] * (O2OxidationRate + OHOxidationRate + OOxidationRate);
        mdYdt[H2Ind] += m_vol * mw[H2Ind] * (NucleationRate + SurfaceGrowthRate);
        mdYdt[HInd] += m_vol * mw[HInd] * (OHOxidationRate);
        //Solid Caron Source Term
        LHS[m_nv - 2] = TotalDensity;
        RHS[m_nv - 2] = MW_Carbon*(2*(NucleationRate+SurfaceGrowthRate) - O2OxidationRate - OOxidationRate - OHOxidationRate);
        //Nd Source Term
        LHS[m_nv - 1] = TotalDensity;
        RHS[m_nv - 1] = NdNucleationConversionTerm * NucleationRate - AgglomerationRate;

        //Now Add Terms To Energy Equation (TODO:: It wasn't working,switched to constant energy solver)
        // compression work and external heat transfer
        //mcvdTdt += -m_pressure * m_vdot + m_Qdot;
        //Add Solid Carbon Contribution to energy
        //mcvdTdt -= 0000000 * Uc * m_vol;
        mcvdTdt = 1;
  


        RHS[1] = 0;
        if (m_energy) {
            //TotalMass*Total Cv T\dot = [RHS](-sum_allspecies (Total Mass creation rate) *(internal species energy) )
           // LHS[2] = m_mass * CvTot; 

            LHS[2] = 1;
            //This below is to check that things make sense with the addition of carbon and will pull
            //out the carbon part out of the temperature equation for now
            //LHS[2] = m_mass * (1 - Y_Carbon) * m_thermo->cv_mass();// m_g C_v_g T\dot = RHS (this made sense
        }
        else {
            RHS[2] = 0;
        }
    }

    void IdealGasReactorSoot::SolveTemperatureAndSetState(doublereal gasDensity) {
        auto Etot_err = [this, gasDensity](double T) {
            m_thermo->setState_TR(T, gasDensity);
            return CalculateTotalEnergy(T) - TotalEnergy;
        };
        double TGuess = m_thermo->temperature();
        boost::uintmax_t maxiter = 100;
        std::pair<double, double> TT;
        try {
            TT = bmt::bracket_and_solve_root(
                Etot_err, TGuess, 1.2, true, bmt::eps_tolerance<double>(48), maxiter);
        }
        catch (std::exception&) {
            // Try full-range bisection if bracketing fails (for example, near
            // temperature limits for the phase's equation of state)
            try {
                TT = bmt::bisect(Etot_err, m_thermo->minTemp(), m_thermo->maxTemp(),
                    bmt::eps_tolerance<double>(48), maxiter);
            }
            catch (std::exception& err2) {
                // Set m_thermo back to a reasonable state if root finding fails
                m_thermo->setState_TR(TGuess, m_mass / m_vol);
                throw CanteraError("Reactor::updateState",
                    "{}\nat U = {}, rho = {}", err2.what(), TotalEnergy, m_mass / m_vol);
            }
        }
        if (fabs(TT.first - TT.second) > 1e-7 * TT.first) {
            throw CanteraError("Reactor::updateState", "root finding failed");
        }
        m_thermo->setState_TR(TT.second, gasDensity);
    }

    doublereal IdealGasReactorSoot::CalculateTotalEnergy(doublereal T) {
        double e_s = GasConstant / MW_Carbon * (T * CalculateSootH_k_RuT(T)) - OneAtm / CarbonDensity;
        double e_g = m_thermo->intEnergy_mass();
        return Y_Carbon * e_s + (1 - Y_Carbon) * e_g;
    }

    //Calculating The reactions having to do with soot.
    //!This method manages solving all of the soot reaction rates
    void IdealGasReactorSoot::CalculateSootRates(doublereal T) {
        CalculateNucleationRate(T, m_Concentration_k[C2H2Ind]);
        CalculateSurfaceGrowthRate(T, m_Concentration_k[C2H2Ind]);
        CalculateAgglomerationRate(T);
        CalculateO2OxidationRate(T, m_Concentration_k[O2Ind]);
        CalculateOOxidationRate(T, m_Concentration_k[OInd]);
        CalculateOHOxidationRate(T, m_Concentration_k[OHInd]);
    }
    //! The first step in the utilized pathway is the nucleation of Soot for: C2H2 -> 2C(s) + H2
    //! Here 2 C(s) atoms are formed from one C2H2 molecule, Along with a 2*Avagodro's Number/ Number of carbon atoms in the resulting carbon particale worth of carbon particles
    void IdealGasReactorSoot::CalculateNucleationRate(doublereal T, doublereal C2H2Conc)
    {   //Assumes C2H2Conc given in kmol/m^3 (Possibly need to add a volume correction) //Below is the nucleation rate given by Leung et al
        //NucleationRate = std::max(0.0, C2H2Conc) * 1350000. * std::exp(-21100. / T) * (1 - TotalDensity / CarbonDensity * Y_Carbon);
        //Below is the nucleation rate given by liu et al, and also what is used by zimmer et al
        NucleationRate = std::max(0.0, C2H2Conc) * 1000 * std::exp(-16103 / T) * (1 - TotalDensity / CarbonDensity * Y_Carbon);
    }
    //! The second step is the surface growth of an existing carbon particle: C2H2 + (n) C(s) -> (n+2) C(s) + H2, 
    //! where more C2H2 deposites onto an existing carbon particle. Here no new carbon particles are formed
    //! So The Nd of carbon particles doesn't change
    void IdealGasReactorSoot::CalculateSurfaceGrowthRate(doublereal T, doublereal C2H2Conc)
    {       //From Linstedt and Leung, Still assumes concentrations are given in kmol/m^3
            //R2 = k_2 * sqrt( Pi * (6MW_c/Pi/carbonDenstiy)^(2/3) )* (totalDensity^1/2) * (N)^1/6 *(Y_Carbon)^1/3
        //SurfaceGrowthRate = SurfGrowthRateConstant * std::max(0.0, C2H2Conc) * std::pow(Y_Carbon, 1. / 3.) * std::pow(Nd, 1. / 6.) * std::exp(-12100. / T) * (1 - TotalDensity / CarbonDensity * Y_Carbon); //Inherently always positive.... so it damn better be
            //Below is the nucleation rate given by liu et al, and also what is used by zimmer et al
        SurfaceGrowthRate = 700./6000. * SurfGrowthRateConstant * std::max(0.0, C2H2Conc) * std::pow(Y_Carbon, 1. / 3.) * std::pow(Nd, 1. / 6.) * std::exp(-10064. / T) * (1 - TotalDensity / CarbonDensity * Y_Carbon);
    }
    //! The Next step considers the agglomeration (Combination) of existing soot particles 2Particles of carbon -> 1 particle of carbon
    //! In this step the soot number density decreases as there are now less particles overall, even though the size of the avg carbon particle has increased
    void IdealGasReactorSoot::CalculateAgglomerationRate(doublereal T) {
        AgglomerationRate = AgglomerationConstant * std::pow(T, 1. / 2.) * std::pow(Y_Carbon, 1. / 6.) * std::pow(Nd, 11. / 6.);
        //AgglomerationRate = 0;
    }
    //! The Next 3 reactions have to do with the oxidation pathways of the carbon particles, i.e. the Solid carbon oxidizes through these reactions to from CO and possible Hydrogen
    //! The first is the O2 Pathway C(s) + .5 O2 -> CO
    void IdealGasReactorSoot::CalculateO2OxidationRateLeungEtAl(doublereal T, doublereal O2Conc)
    {
      //  O2OxidationRate = 17800. * std::pow(T, .5) * std::exp(39000. * 4184 / GasConstant / T) * O2Conc * ( Pi * std::pow(6 * MW_Carbon / Pi / CarbonDensity, 2. / .3) ) / MW_Carbon;
        //return O2OxidationConstant * O2Conc * std::pow(T, .5) * std::exp(-19680. / T) *
          //  std::pow(Y_Carbon, 2. / 3.) * std::pow(Nd, 1. / 3.);
    }

    void IdealGasReactorSoot::CalculateO2OxidationRate(doublereal T, doublereal O2Conc) {
        double ka = 200. * std::exp(-15098. / T);
        double kz = 21.3 * std::exp(2063. / T);
        double kb = 4.46e-2 * std::exp(-7650. / T);
        double kT = 1.51e6 * std::exp(-48817. / T);
        double PO2 = O2Conc * GasConstant * T/OneAtm;
        double xA = 1. / (1 + kT / (kb * PO2));
        O2OxidationRate = (ka * PO2 * xA / (1 + kz * PO2) + kb * PO2 * (1 - xA) ) * SurfaceAreaConstant * std::pow(Nd,1./3.) * std::pow(Y_Carbon,2./3.);
        //O2OxidationRate = 0;
    }

    //! The second is the O pathway C(s) + O -> CO
    void IdealGasReactorSoot::CalculateOOxidationRate(doublereal T, doublereal OConc) {
        //O Concentration assumed to be in kmol/m^3
        OOxidationRate =  0.001094 * OxidationCollisionEfficiency * std::max(0.0, OConc) * std::pow(T, 1. / 2.) * GasConstant * SurfaceAreaConstant * std::pow(Nd, 1. / 3.) * std::pow(Y_Carbon, 2. / 3.) * (1 - TotalDensity / CarbonDensity * Y_Carbon);
        //OOxidationRate = 0.;
    }
    //! The third is the OH pathway C(s) + OH -> CO + H
    void IdealGasReactorSoot::CalculateOHOxidationRate(doublereal T, doublereal OHConc) {
        //OH Concentration assumed to be in kmol/m^3
        OHOxidationRate = 0.001044 * OxidationCollisionEfficiency * std::max(0.0, OHConc) * std::pow(T, 1. / 2.) * GasConstant * SurfaceAreaConstant * std::pow(Nd, 1. / 3.) * std::pow(Y_Carbon, 2. / 3.) * (1 - TotalDensity / CarbonDensity * Y_Carbon);
        //OHOxidationRate = 0;
    }
    



    size_t IdealGasReactorSoot::componentIndex(const string& nm) const
    {
        size_t k = speciesIndex(nm);
        if (k != npos) {
            return k + 3;
        } else if (nm == "mass") {
            return 0;
        } else if (nm == "volume") {
            return 1;
        } else if (nm == "temperature") {
            return 2;
        } else if (nm == "massCarbon") {
            return m_nv - 2;
        }
        else if (nm == "sootNumberDensity") {
            return m_nv - 1;
        } else {
            return npos;
        }
    }

    std::string IdealGasReactorSoot::componentName(size_t k) {
        if (k == 2) {
            return "temperature";
        } else if (k == (m_nv - 2)) {
            return "massCarbon";
        } else if (k == (m_nv - 1)) {
            return "sootNumberDensity";
        } 
        else {
            return Reactor::componentName(k);
        }
    }

    doublereal IdealGasReactorSoot::CalculateSootCv(doublereal Temp) {
        //Cv_mass = Cv_mole/MW_Carbon = 1/MWCarbon Ru*(Cp/Ru - 1)
       // double Cp_Ru = 1000. / GasConstant;
        double Cp_Ru = CalculateSootCp_Ru(Temp);
        return GasConstant * (Cp_Ru- 1.) / MW_Carbon;
    }
    doublereal IdealGasReactorSoot::CalculateSootCp_Ru(doublereal Temp) {
        if ((Temp < 200) )
        {
            double t200 = CSNasa7TLow.at(0) + 200. * (CSNasa7TLow.at(1) + 200. * (CSNasa7TLow.at(2) + 200. * (CSNasa7TLow.at(3) + 200. * CSNasa7TLow.at(4))));
            double t300 = CSNasa7TLow.at(0) + 300. * (CSNasa7TLow.at(1) + 300. * (CSNasa7TLow.at(2) + 300. * (CSNasa7TLow.at(3) + 300. * CSNasa7TLow.at(4))));
            return t200 - (t300 - t200) / 100. * (200 - Temp);
        }
        else if (Temp <= 1000)
            return CSNasa7TLow.at(0) + Temp * (CSNasa7TLow.at(1) + Temp * (CSNasa7TLow.at(2) + Temp * (CSNasa7TLow.at(3) + Temp * CSNasa7TLow.at(4))));
        else if (Temp <= 5000)
            return CSNasa7THigh.at(0) + Temp * (CSNasa7THigh.at(1) + Temp * (CSNasa7THigh.at(2) + Temp * (CSNasa7THigh.at(3) + Temp * CSNasa7THigh.at(4))));
        else {
            double t5000 = CSNasa7THigh.at(0) + 5000. * (CSNasa7THigh.at(1) + 5000. * (CSNasa7THigh.at(2) + 5000. * (CSNasa7THigh.at(3) + 5000. * CSNasa7THigh.at(4))));
            double t4900 = CSNasa7THigh.at(0) + 4900. * (CSNasa7THigh.at(1) + 4900. * (CSNasa7THigh.at(2) + 4900. * (CSNasa7THigh.at(3) + 4900. * CSNasa7THigh.at(4))));
            return t5000 + (t5000 - t4900) / 100. * (Temp - 5000.);
        }
    }
    
    doublereal IdealGasReactorSoot::CalculateSootU_k(doublereal Temp) {
        return GasConstant * Temp * (CalculateSootH_k_RuT(Temp) - 1.);
    }
    doublereal IdealGasReactorSoot::CalculateSootH_k_RuT(doublereal Temp) {
      //  return 1000. / GasConstant;
        if ( (Temp < 200.) )
        {
            double t200 = CSNasa7TLow.at(5) / 200. + CSNasa7TLow.at(0) + 200. * (CSNasa7TLow.at(1) / 2. +
                200. * (CSNasa7TLow.at(2) / 3. + 200. * (CSNasa7TLow.at(3) / 4. +
                    200. * CSNasa7TLow.at(4) / 5.)));
            double t300 = CSNasa7TLow.at(5) / 300. + CSNasa7TLow.at(0) + 300. * (CSNasa7TLow.at(1) / 2. +
                300. * (CSNasa7TLow.at(2) / 3. + 300. * (CSNasa7TLow.at(3) / 4. +
                    300. * CSNasa7TLow.at(4) / 5.)));
            return t200 - (t300 - t200) / 100. * (200. - Temp);
        }
        else if (Temp <= 1000.)
            return CSNasa7TLow.at(5)/Temp + CSNasa7TLow.at(0) + Temp * (CSNasa7TLow.at(1)/2. +
                                  Temp * (CSNasa7TLow.at(2)/3. + Temp * (CSNasa7TLow.at(3)/4. +
                                  Temp * CSNasa7TLow.at(4)/5. ) ) );
        else if( Temp <= 5000.)
            return CSNasa7THigh.at(5)/Temp + CSNasa7THigh.at(0) + Temp * (CSNasa7THigh.at(1)/2. + 
                                  Temp * (CSNasa7THigh.at(2)/3. + Temp * (CSNasa7THigh.at(3)/4. + 
                                  Temp * CSNasa7THigh.at(4)/5. ) ) );
        else {
            double t5000 = CSNasa7THigh.at(5) / 5000. + CSNasa7THigh.at(0) + 5000. * (CSNasa7THigh.at(1) / 2. +
                5000. * (CSNasa7THigh.at(2) / 3. + 5000. * (CSNasa7THigh.at(3) / 4. +
                    5000. * CSNasa7THigh.at(4) / 5.)));
            double t4900 = CSNasa7THigh.at(5) / 4900. + CSNasa7THigh.at(0) + 4900. * (CSNasa7THigh.at(1) / 2. +
                4900. * (CSNasa7THigh.at(2) / 3. + 4900. * (CSNasa7THigh.at(3) / 4. +
                    4900. * CSNasa7THigh.at(4) / 5.)));
            return t5000 + (t5000 - t4900) / 100. * (Temp - 5000.);
        }

    }
    

    doublereal IdealGasReactorSoot::GetCarbonMassFrac() {
        return Y_Carbon;
    }
    doublereal IdealGasReactorSoot::GetNd() {
        return Nd;
    }



}
