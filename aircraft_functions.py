import numpy as np
# import math # Not strictly needed here if np handles everything, but can keep if desired for other functions.

# --- ISA Constants (moved from isa_model function for broader use) ---
# Reference: Anderson, John D. Jr. "Aircraft Performance and Design." Chapter 2, "The Standard Atmosphere." McGraw-Hill, 1999.
G0 = 9.80665    # Standard gravity (m/s^2)
R_AIR = 287.05287 # Specific gas constant for dry air (J/(kg·K))
GAMMA = 1.4     # Ratio of specific heats for air (dimensionless)

# Sea Level ISA (Troposphere base, H = 0 km)
T_SL_K = 288.15        # Temperature at Sea Level (15 degrees Celsius)
P_SL_PA = 101325       # Pressure at Sea Level (1 atmosphere)
RHO_SL_KGM3 = 1.225    # Density at Sea Level

# Troposphere (0 to 11,000 m) Constants
L_TROPOSPHERE_K_PER_M = -0.0065  # Temperature lapse rate in troposphere (K/m)
H_TROPOPAUSE_M = 11000   # Altitude of tropopause (m)

# Lower Stratosphere (11,000 m to 20,000 m) Constants (Isothermal layer)
H_STRATOSPHERE_TOP_M = 20000 # Top of the isothermal layer considered for this model


# --- Task 1.2: Implement ISA Model Function (isa_model) ---
def isa_model(altitude_m):
    """
    Calculates International Standard Atmosphere (ISA) properties at a given altitude.
    Considers troposphere (0-11km) and lower stratosphere (11-20km, isothermal).
    Beyond 20km, it continues the isothermal trend as a simplification for this project's scope.

    Args:
        altitude_m (float): Altitude in meters.

    Returns:
        dict: A dictionary containing:
            'rho_kgm3': Air density (kg/m^3)
            'p_pa': Static pressure (Pa)
            'T_k': Static temperature (K)
            'a_ms': Speed of sound (m/s)
    """
    if altitude_m <= H_TROPOPAUSE_M:
        # Troposphere equations
        T_k = T_SL_K + L_TROPOSPHERE_K_PER_M * altitude_m
        P_pa = P_SL_PA * (T_k / T_SL_K)**(-G0 / (L_TROPOSPHERE_K_PER_M * R_AIR))
        rho_kgm3 = RHO_SL_KGM3 * (T_k / T_SL_K)**(-(G0 / (L_TROPOSPHERE_K_PER_M * R_AIR)) - 1)
    elif altitude_m <= H_STRATOSPHERE_TOP_M:
        # Isothermal Stratosphere (11km to 20km)
        # Calculate properties at the tropopause (base of this layer)
        T_tropopause_K = T_SL_K + L_TROPOSPHERE_K_PER_M * H_TROPOPAUSE_M
        P_tropopause_Pa = P_SL_PA * (T_tropopause_K / T_SL_K)**(-G0 / (L_TROPOSPHERE_K_PER_M * R_AIR))

        T_k = T_tropopause_K # Temperature is constant in this layer
        P_pa = P_tropopause_Pa * np.exp(-G0 * (altitude_m - H_TROPOPAUSE_M) / (R_AIR * T_k))
        rho_kgm3 = P_pa / (R_AIR * T_k) # Density from Ideal Gas Law
    else:
        # Beyond 20km (simplified: continue isothermal trend from 20km values)
        # Recursively call isa_model to get values at 20km boundary without recomputing
        # (Note: This is a slight simplification; normally, you'd calculate these boundary values once)
        temp_at_20km = isa_model(H_STRATOSPHERE_TOP_M)['T_k']
        pressure_at_20km = isa_model(H_STRATOSPHERE_TOP_M)['p_pa']

        T_k = temp_at_20km # Assume temperature remains constant
        P_pa = pressure_at_20km * np.exp(-G0 * (altitude_m - H_STRATOSPHERE_TOP_M) / (R_AIR * T_k))
        rho_kgm3 = P_pa / (R_AIR * T_k) # Density from Ideal Gas Law

    # Speed of sound (a = sqrt(gamma * R * T))
    a_ms = np.sqrt(GAMMA * R_AIR * T_k)

    return {
        'rho_kgm3': rho_kgm3,
        'p_pa': P_pa,
        'T_k': T_k,
        'a_ms': a_ms
    }

# Placeholder for future functions (Task 1.3 - 1.9, etc.)
# aircraft_params = {} # This will be defined in the main script now
# def thrust_available(...): pass
# def fuel_burn_rate(...): pass
# def calculate_lift(...): pass
# def calculate_drag(...): pass
# def thrust_required(...): pass
# def calculate_v_stall(...): pass


# --- TEMPORARY DEBUGGING BLOCK (DELETE AFTER USE) ---
# Used to get precise expected values for unit tests.
# Run this file directly: python aircraft_functions.py
# Copy the exact output values into your test_aircraft_functions.py file.
print("\n--- DEBUGGING ISA MODEL OUTPUTS FOR UNIT TESTS ---")
debug_altitudes = [5000, 11000, 15000] # Altitudes that failed

for alt_m in debug_altitudes:
    result = isa_model(alt_m)
    print(f"\n--- Altitude: {alt_m}m ---")
    print(f"Expected T_k = {result['T_k']:.15f}")
    print(f"Expected P_pa = {result['p_pa']:.15f}")
    print(f"Expected rho_kgm3 = {result['rho_kgm3']:.15f}")
    print(f"Expected a_ms = {result['a_ms']:.15f}")
print("--------------------------------------------------\n")