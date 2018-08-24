// Copyright (c) 2018 Stig Rune Sellevag
//
// This file is distributed under the MIT License. See the accompanying file
// LICENSE.txt or http://www.opensource.org/licenses/mit-license.php for terms
// and conditions.

#ifndef NUMLIB_CONSTANTS_H
#define NUMLIB_CONSTANTS_H

namespace Numlib {

//------------------------------------------------------------------------------
//
// Provides definitions of mathematical constants, metric prefixes, physical
// constants, metric prefixes, and conversion factors.

// Mathematical constants:

constexpr double e = 2.7182818284590452354;
constexpr double pi = 3.14159265358979323846;
constexpr double sqrt2 = 1.41421356237309504880;

// Metric prefixes:

constexpr double yotta = 1.0e+24;
constexpr double zetta = 1.0e+21;
constexpr double exa = 1.0e+18;
constexpr double peta = 1.0e+15;
constexpr double tera = 1.0e+12;
constexpr double giga = 1.0e+9;
constexpr double mega = 1.0e+6;
constexpr double kilo = 1.0e+3;
constexpr double hecto = 1.0e+2;
constexpr double deca = 10.0;
constexpr double one = 1.0;
constexpr double deci = 1.0e-1;
constexpr double centi = 1.0e-2;
constexpr double milli = 1.0e-3;
constexpr double micro = 1.0e-6;
constexpr double nano = 1.0e-9;
constexpr double pico = 1.0e-12;
constexpr double femto = 1.0e-15;
constexpr double atto = 1.0e-18;
constexpr double zepto = 1.0e-21;
constexpr double yocto = 1.0e-24;

// CODATA recommended 2014 values for physical constants:
// Source: http://physics.nist.gov/constants
constexpr double m_u = 1.66053904000e-27;    // Atomic mass constant (kg)
constexpr double N_A = 6.02214085700e+23;    // Avogadro constant (mol^-1)
constexpr double a_0 = 0.529177210670;       // Bohr radius (1.0e-10 m)
constexpr double k = 1.38064852000e-23;      // Boltzmann constant (J K^-1)
constexpr double G_0 = 7.74809173100e-5;     // Conductance quantum (S)
constexpr double eps_0 = 8.85418781700e-12;  // Electric constant (F m^-1)
constexpr double m_e = 9.10938356000e-31;    // Electron mass (kg)
constexpr double eV = 1.60217662080e-19;     // Electron volt (J)
constexpr double ec = 1.60217662080e-19;     // Elementary charge (C)
constexpr double F = 9.64853328900e+4;       // Faraday constant (C mol^-1)
constexpr double alpha = 7.29735256640e-3;   // Fine-structure constant
constexpr double R = 8.31445980000;          // Gas constant (J mol^-1 K^-1)
constexpr double E_h = 4.35974465000e-18;    // Hartree energy (J)
constexpr double c_0 = 2.99792458e+8;        // Speed of light (m s^-1)
constexpr double mu_0 = 1.25663706140e-6;    // Magnetic constant (N A^-2)
constexpr double phi_0 = 2.06783383100e-15;  // Magnetic flux quantum (Wb)
constexpr double m_p_m_e = 1.83615267389e+3; // Proton-electron mass ratio
constexpr double G = 6.67408000000e-11;      // Newtonian grav. (m^3/(kg s^2))
constexpr double h = 6.62607004000e-34;      // Planck constant (J s)
constexpr double h_bar = 1.05457180000e-34;  // Planck-bar constant (J s)
constexpr double m_p = 1.67262189800e-27;    // Proton mass (kg)
constexpr double R_inf = 1.09737315685e+7;   // Rydberg constant (m^-1)
constexpr double std_atm = 1.01325000000e+5; // Std. atmospheric pressure (Pa)
constexpr double sigma = 5.67036700000e-8;   // Stefan-Boltzmann (W m^-2 K^-4)

// Conversion factors:

constexpr double cal_to_J = 4.184;                       // cal to J
constexpr double icm_to_kJ = 1.19626564e-02;             // cm^-1 to kJ/mol
constexpr double icm_to_K = 100.0 * h * c_0 / k;         // cm^-1 to K
constexpr double J_to_icm = 100.0 * h * c_0;             // J to cm^-1
constexpr double au_to_cm = a_0 * 1.0e-8;                // au to cm
constexpr double au_to_icm = E_h / (h * c_0 * 100.0);    // au to cm^-1
constexpr double au_to_s = h_bar / E_h;                  // au to s
constexpr double au_to_K = E_h / k;                      // au to K
constexpr double au_to_kg = m_e;                         // au to kg
constexpr double au_to_kgm2 = m_u * a_0 * a_0 * 1.0e-20; // amu bohr^2 to kg m^2
constexpr double GHz_to_K = giga * 4.79924470000e-11;    // GHz to K

} // namespace Numlib

#endif // NUMLIB_CONSTANTS_H
