### Python prototype for THMPH procedure : goal : make it more robust.

import numpy as np
from .material_mixture import AVOGADRO, CM2_TO_BARN

try:
    from iapws import IAPWS97
    IAPWS_AVAILABLE = True
except ImportError:
    IAPWS_AVAILABLE = False
    print("Warning: iapws library not found. Water properties functions will not work. Please install it using 'pip install iapws'.")


def THMPH(Pressure, Enthalpy):
    """
    Backwards inversion of steam tables to find water density and temperature.
    based on A. Hebert and P. Gallet's THMPH function. 

    Parmeters :
    ----------
        - Pressure: Pressure in Pa
        - Enthalpy: Enthalpy in J/Kg
    Returns:
    ----------
        - Temperature: Temperature in Kelvin
        - Density: Density in Kg/m^3
    """

    # Create an instance of the IAPWS97 class for water
    print(f"THMPH : Pressure: {Pressure} Pa, Enthalpy: {Enthalpy} J/Kg")
    water = IAPWS97(P=Pressure/1e6, h=Enthalpy/1e3)

    # Extract temperature and density from the water object
    Temperature = water.T
    Density = water.rho

    return Temperature, Density

def THMPT(Pressure, Temperature):
    """
    Parameters :
    ----------
        - Pressure: Pressure in Pascal
        - Temperature: Temperature in Kelvin

    Returns:
    ----------
        - Density: Density in Kg/m^3
        - Enthalpy: Enthalpy in J/Kg
        - Conductivity (k): Thermal conductivity in W/m/K
        - Viscosity (mu): Dynamic viscosity in Pa.s
        - Cp: Specific heat capacity in J/Kg/K
    """

    # Create an instance of the IAPWS97 class for water
    water = IAPWS97(P=Pressure/1e6, T=Temperature)

    # Extract properties from the water object
    Density = water.rho
    Enthalpy = water.h
    Conductivity = water.k
    Viscosity = water.mu
    Cp = water.cp

    return Density, Enthalpy*1e3, Conductivity, Viscosity, Cp

def THMSAT(Pressure):
    """
    Get saturation temperature at Pressure
    Parameters :
    ----------
        - Pressure: Pressure in Pascal 
    Returns:
    ----------
        - Temperature: Saturation temperature in Kelvin
    """
    # Create an instance of the IAPWS97 class for water
    water = IAPWS97(P=Pressure/1e6, x=0.0)

    # Extract saturation temperature from the water object
    Temperature = water.T

    return Temperature

def THMTX(Temperature, Quality):
    """
    Get Density, enthalpy, conductivity, dynamic viscosity and specific heat at Temperature and Quality
    Parameters :
    ----------
        - Temperature: Temperature in Kelvin
        - Quality: Quality of the mixture (0 for liquid, 1 for vapor)
    Returns:
    ----------
        - Density (rho): Density in Kg/m^3 at 
        - Enthalpy (h): Enthalpy in J/Kg
        - Conductivity (k): Thermal conductivity in W/m/K
        - Viscosity (mu): Dynamic viscosity in Pa.s
        - Specific heat capacity (Cp) in J/Kg/K

    """
    # Create an instance of the IAPWS97 class for water
    print(f"THMTX : Temperature: {Temperature} K, Quality: {Quality}")
    water = IAPWS97(T=Temperature, x=Quality)

    # Extract properties from the water object
    Density = water.rho
    Enthalpy = water.h
    Conductivity = water.k
    Viscosity = water.mu
    Cp = water.cp

    return Density, Enthalpy*1e3, Conductivity, Viscosity, Cp

def THMSAP(Temperature):
    """
    return the saturation pressure (Pa) as a function of the temperature (K)
    Parameters :
    ----------
        - Temperature: Temperature in Kelvin
    Returns:
    ----------
        - Pressure: Saturation pressure in Pascal

    """

    # Create an instance of the IAPWS97 class for water
    water = IAPWS97(T=Temperature, x=0)

    # Extract saturation pressure from the water object
    Pressure = water.P

    return Pressure

def compute_water_iso_densities_at_densities(densities):
    """Compute H and O isotopic number densities for light water at
    multiple mass densities.

    Parameters
    ----------
    densities : list[float]
        Mass densities in g/cm³.

    Returns
    -------
    list[dict]
        One ``{"H1": N_H, "O16": N_O}`` dict per density, in
        atoms/barn·cm.
    """
    M_H2O = 15.9994 + 2.0 * 1.00794
    results = []
    for rho in densities:
        N_MAT = rho * AVOGADRO / (M_H2O * CM2_TO_BARN)
        results.append({"H1": 2.0 * N_MAT, "O16": N_MAT})
    return results


if __name__ == "__main__":
    # Example usage
    Pressure = 7e6  # Pa
    
    # saturation temperature at 7 MPa
    T_sat = THMSAT(Pressure)
    print(f"Saturation temperature at {Pressure} Pa: {T_sat:.2f}K")

    # density of saturated liquid at 7 MPa
    rho_sat_liquid, _, _, _, _ = THMTX(T_sat, Quality=0)
    print(f"Density of saturated liquid at {Pressure} Pa: {rho_sat_liquid:.2f} kg/m^3")
    # density of saturated vapor at 7 MPa
    rho_sat_vapor, _, _, _, _ = THMTX(T_sat, Quality=1)
    print(f"Density of saturated vapor at {Pressure} Pa: {rho_sat_vapor:.2f} kg/m^3")

    rho_list = []
    for alpha in [0.0, 0.40, 0.80]:
        rho_m = alpha*rho_sat_vapor + (1-alpha)*rho_sat_liquid
        rho_list.append(rho_m*1e-3)  # convert to g/cm^3
        print(f"Density of mixture at {Pressure} Pa and quality {alpha}: {rho_m:.2f} kg/m^3")

    results = compute_water_iso_densities_at_densities(rho_list)

    for rho, res in zip(rho_list, results):
        print(f"At density {rho:.4f} g/cm³, alpha = {0.0 if rho == rho_sat_liquid*1e-3 else 1.0 if rho == rho_sat_vapor*1e-3 else (rho - rho_sat_liquid*1e-3)/(rho_sat_vapor*1e-3 - rho_sat_liquid*1e-3):.2f}:")
        print(f"H1 = {res['H1']:.6E} atoms/barn·cm, O16 = {res['O16']:.6E} atoms/barn·cm")