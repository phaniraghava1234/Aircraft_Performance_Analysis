import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math

# ---  Implement ISA Model Function (isa_model) ---
# Reference for ISA equations:
# Anderson, John D. Jr. "Aircraft Performance and Design." Chapter 2, "The Standard Atmosphere." McGraw-Hill, 1999.
# The implementation covers the troposphere (0-11km) and the lower stratosphere (11-20km, isothermal layer).
# Beyond 20km, it continues the isothermal trend as a simplification for this project's scope.
def isa_model(altitude_m):
    """
    Calculates International Standard Atmosphere (ISA) properties at a given altitude.

    Args:
        altitude_m (float): Altitude in meters.

    Returns:
        dict: A dictionary containing:
            'rho_kgm3': Air density (kg/m^3)
            'p_pa': Static pressure (Pa)
            'T_k': Static temperature (K)
            'a_ms': Speed of sound (m/s)
    """
    # ISA Fundamental Constants
    g0 = 9.80665    # Standard gravity (m/s^2)
    R = 287.05287   # Specific gas constant for dry air (J/(kg·K))
    gamma = 1.4     # Ratio of specific heats for air (dimensionless)

    # Sea Level ISA (Troposphere base, H = 0 km)
    T_sl_K = 288.15        # Temperature at Sea Level (15 degrees Celsius)
    P_sl_Pa = 101325       # Pressure at Sea Level (1 atmosphere)
    rho_sl_kgm3 = 1.225    # Density at Sea Level

    # Troposphere (0 to 11,000 m) Constants
    L_troposphere_K_per_m = -0.0065  # Temperature lapse rate in troposphere (K/m)
    H_tropopause_m = 11000   # Altitude of tropopause (m)

    # Lower Stratosphere (11,000 m to 20,000 m) Constants (Isothermal layer)
    H_stratosphere_top_m = 20000 # Top of the isothermal layer considered for this model

    if altitude_m <= H_tropopause_m:
        # Troposphere equations (temperature decreases linearly, pressure and density follow power laws)
        T_k = T_sl_K + L_troposphere_K_per_m * altitude_m
        P_pa = P_sl_Pa * (T_k / T_sl_K)**(-g0 / (L_troposphere_K_per_m * R))
        rho_kgm3 = rho_sl_kgm3 * (T_k / T_sl_K)**(-(g0 / (L_troposphere_K_per_m * R)) - 1)
    elif altitude_m <= H_stratosphere_top_m:
        # Isothermal Stratosphere (11km to 20km)
        # First, calculate properties at the tropopause (base of this layer)
        T_tropopause_K = T_sl_K + L_troposphere_K_per_m * H_tropopause_m
        P_tropopause_Pa = P_sl_Pa * (T_tropopause_K / T_sl_K)**(-g0 / (L_troposphere_K_per_m * R))

        T_k = T_tropopause_K # Temperature is constant in this layer
        P_pa = P_tropopause_Pa * np.exp(-g0 * (altitude_m - H_tropopause_m) / (R * T_k))
        rho_kgm3 = P_pa / (R * T_k) # Density from Ideal Gas Law (P = rho * R * T)
    else:
        # Beyond 20km (simplified: continue isothermal trend from 20km values for this project's scope)
        # A more complete ISA model defines further layers (e.g., 20-32km with positive lapse, etc.).
        # For our purposes, this simplified extension is sufficient.
        # Calculate properties at the top of the considered isothermal layer (20km)
        temp_at_20km = isa_model(H_stratosphere_top_m)['T_k']
        pressure_at_20km = isa_model(H_stratosphere_top_m)['p_pa']

        T_k = temp_at_20km # Assume temperature remains constant
        P_pa = pressure_at_20km * np.exp(-g0 * (altitude_m - H_stratosphere_top_m) / (R * T_k))
        rho_kgm3 = P_pa / (R * T_k) # Density from Ideal Gas Law

    # Speed of sound (a = sqrt(gamma * R * T))
    a_ms = np.sqrt(gamma * R * T_k)

    return {
        'rho_kgm3': rho_kgm3,
        'p_pa': P_pa,
        'T_k': T_k,
        'a_ms': a_ms
    }