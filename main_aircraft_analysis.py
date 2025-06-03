import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math # For mathematical functions like log, pi, etc.

# --- Import functions and constants from our custom module ---
# Reference: Standard Python module import practice.
# This allows using functions like isa_model and constants like RHO_SL_KGM3 defined in aircraft_functions.py
from aircraft_functions import isa_model, RHO_SL_KGM3, G0, R_AIR, GAMMA

#Define Comprehensive Aircraft Parameters ---
# Reference:
# - Values are illustrative, inspired by a light general aviation aircraft (e.g., Cessna 172) for the 'prop' type.
# - Parameters chosen based on typical inputs for aircraft performance analysis as described in:
#   - Anderson, John D. Jr. "Aircraft Performance and Design." McGraw-Hill, 1999.
#   - Raymer, Daniel P. "Aircraft Design: A Conceptual Approach." AIAA Education Series.
# All parameters are defined in SI units (meters, kilograms, seconds, Newtons, Watts, Pascals, Kelvin)
aircraft_params = {
    # 1. Weights & Capacities (Newtons, kg)
    # Note: Masses are converted to Newtons (mass * g0) for consistency in force calculations.
    'weight_ref_N': 1100 * G0, # Reference weight for calculations (e.g., typical operating weight in Newtons)
                               # Using a value slightly less than MTOW for typical cruise/climb analysis.
    'weight_MTOW_N': 1111 * G0, # Maximum Takeoff Weight (N) - Cessna 172 example: 2450 lbs = 1111 kg
    'weight_MLW_N': 1111 * G0,  # Maximum Landing Weight (N) - For C172, MLW = MTOW
    'weight_OEW_N': 767 * G0,   # Operating Empty Weight (N) - Cessna 172 example: 1697 lbs = 767 kg
    'fuel_capacity_kg': 159,    # Maximum usable fuel load capacity (kg) - Cessna 172 example: 53 gal = 159 kg (at 0.72 kg/L)
    'payload_capacity_kg': 344, # Max Payload (MTOW - OEW - Max Fuel) for range diagram purposes.

    # 2. Dimensions & Geometry (m, m^2, dimensionless)
    'wing_area_m2': 16.2,         # Wing reference area (m^2) - Cessna 172 example: 174 sq ft = 16.2 m^2
    'wing_span_m': 11.0,          # Wing span (m) - Cessna 172 example: 36 ft 1 in = 11.0 m
    'aspect_ratio': (11.0**2) / 16.2, # Aspect Ratio (dimensionless) - calculated from b and S

    # 3. Aerodynamic Coefficients (Dimensionless, for different configurations)
    # Oswald efficiency factor (e) - affects induced drag (K = 1 / (pi * AR * e))
    'oswald_efficiency_clean': 0.75, # Typical value for general aviation aircraft
    'oswald_efficiency_TO': 0.70,    # Slightly lower efficiency with takeoff flaps
    'oswald_efficiency_Land': 0.65,  # Even lower efficiency with full landing flaps

    # Zero-Lift Drag Coefficient (CD0)
    'cd0_clean': 0.025,           # Typical value for clean configuration (no gear/flaps)
    'cd0_TO': 0.040,              # Higher with landing gear and takeoff flaps deployed
    'cd0_Land': 0.080,            # Highest with landing gear and full landing flaps

    # Induced Drag Factor (K) - calculated based on aspect_ratio and oswald_efficiency
    # These will be used to define the true drag polar for simulating data, and for performance calcs
    'k_clean': 1 / (np.pi * ((11.0**2) / 16.2) * 0.75),
    'k_TO': 1 / (np.pi * ((11.0**2) / 16.2) * 0.70),
    'k_Land': 1 / (np.pi * ((11.0**2) / 16.2) * 0.65),

    # Maximum Lift Coefficient (CL_max) for positive stall
    'cl_max_clean': 1.4,          # Max CL clean (no flaps)
    'cl_max_TO': 1.8,             # Max CL with takeoff flaps
    'cl_max_Land': 2.2,           # Max CL with full landing flaps

    # Minimum (most negative) Lift Coefficient for negative stall (optional, for V-n)
    'cl_min_clean': -0.8,         # Typical value for inverted flight/negative G stall

    # Drag coefficient increments (added to CD0_clean for specific states, or already part of TO/Land CD0s)
    # If CD0_TO and CD0_Land already include these, these increments are redundant for total CD0.
    # They are useful for understanding individual contributions, or if calculating CD dynamically.
    'cd_gear_increment': 0.015,   # Typical increment for landing gear extended
    'cd_flaps_TO_increment': 0.010, # Typical increment for takeoff flap setting
    'cd_flaps_Land_increment': 0.040, # Typical increment for landing flap setting (beyond TO)
    'cd_speedbrake_increment': 0.050, # Typical increment for speed brakes deployed

    # 4. Structural Limits (dimensionless load factor 'g')
    'n_max_limit': 3.8,           # Max positive load factor (e.g., for normal category aircraft)
    'n_min_limit': -1.5,          # Min negative load factor

    # 5. Speed Limits (Design Speeds - True Airspeed in m/s)
    # Conversions: 1 knot = 0.514444 m/s, 1 mph = 0.44704 m/s
    'v_mo_ms': 150 * 0.514444,    # Maximum Operating Speed (Vmo) in m/s (approx. 150 knots for C172)
    'v_c_ms': 120 * 0.514444,     # Design Cruise Speed (Vc) in m/s (approx. 120 knots for C172)
    'v_d_ms': 180 * 0.514444,     # Design Dive Speed (Vd) in m/s (approx. 180 knots for C172)
    'v_max_limit_ms': 180 * 0.514444, # The highest speed shown on the V-n diagram (using Vd)

    # Reference speed factors (dimensionless) - applied to stall speed (Vs) at relevant configurations
    # These factors are often derived from regulations (e.g., FAR Part 23)
    'v_r_factor': 1.1,            # Rotation speed factor (Vr = 1.1 * Vs_TO)
    'v_lof_factor': 1.2,          # Liftoff speed factor (Vlof = 1.2 * Vs_TO, minimum)
    'v_td_factor': 1.15,          # Touchdown speed factor (Vtd = 1.15 * Vs_Land)

    # 6. Engine/Propulsion Model Parameters
    'engine_type': 'prop',        # String: 'jet' or 'prop' - defines which engine model to use

    # --- Propeller Engine Example Parameters (Cessna 172-like) ---
    'max_power_sl_W': 119310, # Max Power at Sea Level, ISA (Watts) - approx 160 hp (converted)
                              # Reference: Cessna 172 engines are typically 160-180 hp
    'power_lapse_exponent': 0.7,   # Exponent for altitude effect on power (Power ~ (rho/rho_sl)^exp)
                                   # Typical range 0.6-0.8 for normally aspirated piston engines
    'power_speed_factor': 0.0005,  # Parameter for linear speed effect on power (P = P_max * (1 - k*V))
                                   # Simplified: Represents slight power decrease at higher speeds for fixed-pitch props
    'sfc_sl_per_sec': 8.35e-8, # Specific Fuel Consumption at SL, ISA (kg / (W*s))
                               # Reference: 0.5 lbs/hp/hr for piston, converted to SI: (0.5 * 0.453592) / (745.7 * 3600) = 8.35e-8 kg/W/s

    # --- Jet Engine Example Parameters (Uncomment and set 'engine_type' to 'jet' to use) ---
    # 'max_thrust_sl_N': 20000,      # Max Static Thrust at Sea Level, ISA (Newtons) - e.g., for a small business jet (20 kN)
    # 'thrust_lapse_exponent': 0.8,  # Exponent for altitude effect on thrust (Thrust ~ (rho/rho_sl)^exp)
    # 'thrust_speed_factor': 0.001,  # Parameter for linear speed effect on thrust (T = T_max * (1 - k*V))
    # 'tsfc_sl_per_sec': 1.4e-5, # Thrust Specific Fuel Consumption at SL, ISA (kg / (N*s))
    #                            # Reference: 0.5 lbs/(lbf*hr) for jet, converted to SI: (0.5 * 0.453592) / (4.44822 * 3600) = 1.4e-5 kg/(N*s)
    # Note: For simplicity, fuel consumption lapse with altitude/speed is often ignored or very simplified in these models.

    # 7. Ground Operations Parameters (dimensionless)
    # Reference: Typical values for runway conditions.
    'mu_rolling': 0.02,           # Coefficient of rolling friction (e.g., asphalt runway)
    'mu_braking': 0.40,           # Coefficient of braking friction (e.g., dry asphalt, full braking)
    'cl_ground_roll_TO': 0.3,     # Representative Lift Coefficient during takeoff ground roll (for friction calculation)
    'cl_ground_roll_Land': 0.8,   # Representative Lift Coefficient during landing ground roll (for friction calculation)

    # 8. Atmospheric Parameters (ISA Standard - Fixed, used for consistency if not imported from constants directly)
    # These are already defined globally in aircraft_functions.py, but included here for completeness of parameter list
    # 'isa_sl_rho_kgm3': 1.225,
    # 'isa_sl_temp_K': 288.15,
    # 'isa_lapse_rate_K_per_m': -0.0065,
    # 'isa_tropopause_alt_m': 11000,
    # 'isa_tropopause_temp_K': 216.65,
}
