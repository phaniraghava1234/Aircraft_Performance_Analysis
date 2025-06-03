# Aircraft Performance & Envelope Analysis Tool

**Author:** [Your Name]  
**Date:** [Insert Current Date]  

## 🛫 Overview

The **Aircraft Performance & Envelope Analysis Tool** is a Python-based application designed to analyze and visualize key aspects of aircraft performance using simplified yet educational aerodynamic and propulsion models. It emphasizes conceptual clarity and demonstrates fundamental aircraft performance engineering principles with practical implementation in code.

This project serves as a compact yet comprehensive demonstration of how aircraft design parameters influence performance, while acknowledging simplifications necessary for rapid prototyping.

---

## 🎯 Project Goals

- Implement the **International Standard Atmosphere (ISA)** model for atmospheric properties across altitudes.
- Define a comprehensive aircraft parameter model including aerodynamics, structure, and propulsion.
- Simulate and fit drag polar data using a parabolic model.
- Calculate and plot the **V-n Diagram** to visualize structural and speed limits.
- Determine and plot **Maximum Rate of Climb** and **Maximum Climb Angle** versus altitude.
- Estimate **Service Ceiling** from ROC profile.
- Plot **Thrust Available vs. Speed** and **Thrust Required vs. Speed** across representative altitudes.
- Implement **Breguet Range and Endurance equations** for jet and propeller aircraft.
- Calculate simplified **Takeoff and Landing Ground Roll distances**.
- Simulate the **Effect of Wind** on ground speed and fuel burn.
- *(Optional)* Outline logic for a **Payload-Range Diagram**.

---

## 🧠 Skills Demonstrated

- **Aerodynamics**: Lift, Drag, Drag Polar, Stall, Induced Drag modeling.
- **Performance Estimation**: Climb, Range, Endurance, Takeoff/Landing.
- **Propulsion Modeling**: Jet and propeller thrust/power variation.
- **Flight Envelope Analysis**: V-n Diagram with structural and aerodynamic limits.
- **Atmospheric Modeling**: ISA model for standard atmospheric conditions.
- **Data Fitting**: Synthetic flight test data and `curve_fit` for drag estimation.
- **Numerical Analysis**: Iterative computation across altitudes and speeds.
- **Visualization**: Clear plots using Matplotlib for all major performance metrics.
- **Python Proficiency**: Modular, well-documented code with NumPy/SciPy.

---

## ⚙️ Features

### ✅ Core Implementations

- ISA Model (`isa_model`)
- Aircraft parameter structure (dictionary or class)
- Thrust/Power models with lapse rate
- Fuel consumption estimation
- Synthetic drag polar data simulation
- Drag polar fitting with `scipy.optimize.curve_fit`
- Lift & drag helper functions
- Thrust required for level/climb flight
- V-n Diagram calculation and plotting
- Rate of Climb & Climb Angle vs. Altitude
- Service Ceiling estimation
- Thrust/Drag vs. Speed curves
- Breguet Range & Endurance (Jet/Prop)
- Ground roll distances (Takeoff & Landing)
- Wind effect simulation

---

## 📊 Example Plots

- **Drag Polar Curve**
- **V-n Diagram**
- **Rate of Climb vs Altitude**
- **Thrust Available vs Thrust Required**
- **Breguet Range/Endurance Outputs**
- **Ground Roll Estimates**

---

## 📁 File Structure

