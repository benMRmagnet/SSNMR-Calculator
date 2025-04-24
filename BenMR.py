# Add note to Reference calculator for IUPAC

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from PIL import Image, ImageTk
import os
import sys
import numpy as np

# Gyromagnetic ratios (γ) in 10^7*rad/s/T
GYRO_MAGNETIC_RATIOS = {
"1H": 26.7522128,
"2H": 4.10662791,
"3H": 28.5349779,
"3He": -20.3801587,
"6Li": 3.9371709,
"7Li": 10.3977013,
"8Li": 3.9593,
"9Be": -3.759666,
"10B": 2.8746786,
"11B": 8.5847044,
"13C": 6.728284,
"14N": 1.9337792,
"15N": -2.71261804,
"17O": -3.62808,
"19F": 25.18148,
"21Ne": -2.11308,
"23Na": 7.0808493,
"25Mg": -1.63887,
"27Al": 6.9762715,
"29Si": -5.319,
"31P": 10.8394,
"33S": 2.055685,
"35Cl": 2.624198,
"37Cl": 2.184368,
"37Ar": 3.656,
"39Ar": 2.17,
"39K": 1.2500608,
"40K": -1.5542854,
"41K": 0.68606808,
"43Ca": -1.803069,
"45Sc": 6.5087973,
"47Ti": -1.5105,
"49Ti": -1.51095,
"50V": 2.670649,
"51V": 7.0455117,
"53Cr": -1.5152,
"55Mn": 6.6452546,
"57Fe": 0.8680624,
"59Co": 6.332,
"61Ni": -2.3948,
"63Cu": 7.111789,
"65Cu": 7.60435,
"67Zn": 1.676688,
"69Ga": 6.438855,
"71Ga": 8.181171,
"73Ge": -0.9360303,
"75As": 4.596163,
"77Se": 5.1253857,
"79Br": 6.725616,
"81Br": 7.249776,
"83Kr": -1.0331,
"85Rb": 2.592705,
"87Rb": 8.7864,
"87Sr": -1.1639376,
"89Y": -1.3162791,
"91Zr": -2.49743,
"93Nb": 6.5674,
"95Mo": -1.751,
"97Mo": -1.788,
"99Tc": 6.046,
"99Ru": -1.229,
"101Ru": -1.377,
"103Rh": -0.8468,
"105Pd": -1.23,
"107Ag": -1.0889181,
"109Ag": -1.2518634,
"111Cd": -5.6983131,
"113Cd": -5.9609155,
"113In": 5.8845,
"115In": 5.8972,
"115Sn": -8.8013,
"117Sn": -9.58879,
"119Sn": -10.0317,
"121Sb": 6.4435,
"123Sb": 3.4892,
"123Te": -7.059098,
"125Te": -8.5108404,
"127I": 5.389573,
"129Xe": -7.452103,
"131Xe": 2.209076,
"133Cs": 3.5332539,
"135Ba": 2.6755,
"137Ba": 2.99295,
"138La": 3.557239,
"139La": 3.8083318,
"137Ce": 3.07,
"139Ce": 3.39,
"141Ce": 1.49,
"141Pr": 8.1907,
"143Nd": -1.457,
"145Nd": -0.898,
"143Pm": 7.282,
"147Pm": 3.53,
"147Sm": -1.115,
"149Sm": -0.9192,
"151Eu": 6.651,
"153Eu": 2.9369,
"155Gd": -0.82132,
"157Gd": -1.0769,
"159Tb": 6.431,
"161Dy": -0.9201,
"163Dy": 1.289,
"165Ho": 5.71,
"167Er": -0.77157,
"169Tm": -2.218,
"171Yb": 4.7288,
"173Yb": -1.3025,
"175Lu": 3.0552,
"176Lu": 2.1684,
"177Hf": 1.086,
"179Hf": -0.6821,
"181Ta": 3.2438,
"183W": 1.1282403,
"185Re": 6.1057,
"187Re": 6.1682,
"187Os": 0.6192895,
"189Os": 2.10713,
"191Ir": 0.4812,
"193Ir": 0.5227,
"195Pt": 5.8385,
"197Au": 0.47306,
"199Hg": 4.8457916,
"201Hg": -1.788769,
"203Tl": 15.5393338,
"205Tl": 15.6921808,
"207Pb": 5.58046,
"209Bi": 4.375,
"209Po": 7.35,
"211Rn": 5.76,
"223Fr": 3.85,
"223Ra": 0.86369,
"225Ra": 7.029,
"227Ac": 3.5,
"229Th": 0.4,
"231Pa": 3.21,
"235U": -0.52,
"237Np": 6.01,
"239Pu": 0.974,
"243Am": 1.54,
"247Cm": 0.2
}


# --- Dipolar Tool Window ---
def open_dipolar_tool():
    dipolar_window = tk.Toplevel(main_window)
    dipolar_window.title("Dipolar Distance Tool")
    dipolar_window.geometry("400x250")

    def calculate_dipolar():
        try:
            nucleus1 = nucleus1_var.get().strip()
            nucleus2 = nucleus2_var.get().strip()
            if nucleus1 not in GYRO_MAGNETIC_RATIOS or nucleus2 not in GYRO_MAGNETIC_RATIOS:
                raise KeyError(f"Invalid nucleus: {nucleus1} or {nucleus2}")

            gamma1 = GYRO_MAGNETIC_RATIOS[nucleus1] * 1e7  
            gamma2 = GYRO_MAGNETIC_RATIOS[nucleus2] * 1e7  
            mu_0 = 4 * np.pi * 1e-7
            h_bar = 1.0545718e-34

            mode = dipolar_mode.get()
            if mode == "Distance to Coupling":
                r_angstrom = float(distance_entry.get())
                if r_angstrom <= 0:
                    raise ValueError("Distance must be positive")
                r_meters = r_angstrom * 1e-10  
                D_rad_s = (mu_0 * gamma1 * gamma2 * h_bar) / (4 * np.pi * (r_meters ** 3))
                D_hz = abs(D_rad_s) / (2 * np.pi) 
                coupling_entry.delete(0, tk.END)
                coupling_entry.insert(0, f"{D_hz:.3f}")
            else:  # Coupling to Distance
                D_hz = float(coupling_entry.get())
                if D_hz <= 0:
                    raise ValueError("Coupling must be positive")
                D_rad_s = D_hz * 2 * np.pi  
                r_meters = (abs(mu_0 * gamma1 * gamma2 * h_bar) / (4 * np.pi * D_rad_s)) ** (1/3)
                r_angstrom = r_meters * 1e10  
                distance_entry.delete(0, tk.END)
                distance_entry.insert(0, f"{r_angstrom:.3f}")
        except KeyError as e:
            messagebox.showerror("Error", str(e))
            print(f"KeyError: {e}")
        except ValueError as e:
            messagebox.showerror("Error", str(e))
            print(f"ValueError: {e}")
        except ZeroDivisionError as e:
            messagebox.showerror("Error", "Division by zero occurred")
            print(f"ZeroDivisionError: {e}")
        except Exception as e:
            messagebox.showerror("Error", f"Unexpected error: {str(e)}")
            print(f"Unexpected error: {e}")

    def reset_all():
        distance_entry.delete(0, tk.END)
        coupling_entry.delete(0, tk.END)
        nucleus1_var.set("1H")
        nucleus2_var.set("1H")
        dipolar_mode.set("Distance to Coupling")

    dipolar_mode = tk.StringVar(value="Distance to Coupling")
    tk.Radiobutton(dipolar_window, text="Distance to Coupling", variable=dipolar_mode, value="Distance to Coupling").grid(row=0, column=0, padx=5, pady=5)
    tk.Radiobutton(dipolar_window, text="Coupling to Distance", variable=dipolar_mode, value="Coupling to Distance").grid(row=0, column=1, padx=5, pady=5)

    tk.Label(dipolar_window, text="Bond Distance (Å):").grid(row=1, column=0, padx=5, pady=5)
    distance_entry = tk.Entry(dipolar_window)
    distance_entry.grid(row=1, column=1, padx=5, pady=5)

    tk.Label(dipolar_window, text="Dipolar Coupling (Hz):").grid(row=2, column=0, padx=5, pady=5)
    coupling_entry = tk.Entry(dipolar_window)
    coupling_entry.grid(row=2, column=1, padx=5, pady=5)

    tk.Label(dipolar_window, text="Nucleus 1:").grid(row=3, column=0, padx=5, pady=5)
    nucleus1_var = tk.StringVar(value="1H")
    nucleus1_dropdown = ttk.Combobox(dipolar_window, textvariable=nucleus1_var, values=list(GYRO_MAGNETIC_RATIOS.keys()), state="readonly")
    nucleus1_dropdown.grid(row=3, column=1, padx=5, pady=5)

    tk.Label(dipolar_window, text="Nucleus 2:").grid(row=4, column=0, padx=5, pady=5)
    nucleus2_var = tk.StringVar(value="1H")
    nucleus2_dropdown = ttk.Combobox(dipolar_window, textvariable=nucleus2_var, values=list(GYRO_MAGNETIC_RATIOS.keys()), state="readonly")
    nucleus2_dropdown.grid(row=4, column=1, padx=5, pady=5)

    tk.Button(dipolar_window, text="Calculate", command=calculate_dipolar).grid(row=5, column=1, padx=5, pady=10)
    tk.Button(dipolar_window, text="Reset All", command=reset_all).grid(row=6, column=1, padx=5, pady=10)

# --- Chemical Shift Conversion Tool ---
def open_shift_tool():
    shift_window = tk.Toplevel(main_window)
    shift_window.title("Chemical Shift Converter")
    shift_window.geometry("350x500")

    def convert_to_haeberlen():
        try:
            d11 = float(d11_entry.get())
            d22 = float(d22_entry.get())
            d33 = float(d33_entry.get())
            deltas = [d11, d22, d33]
            delta_iso = (d11 + d22 + d33) / 3
            deviations = sorted([(abs(d - delta_iso), d) for d in deltas], reverse=True)
            delta_zz = deviations[0][1]
            delta_xx = deviations[2][1]
            delta_yy = deviations[1][1]
            aniso = delta_zz - delta_iso
            eta = abs((delta_yy - delta_xx) / aniso) if aniso != 0 else 0  # Take absolute value
            eta = max(0, min(1, eta))
            iso_entry.delete(0, tk.END); iso_entry.insert(0, f"{delta_iso:.2f}")
            aniso_entry.delete(0, tk.END); aniso_entry.insert(0, f"{aniso:.2f}")
            eta_entry.delete(0, tk.END); eta_entry.insert(0, f"{eta:.3f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input for δ₁₁, δ₂₂, or δ₃₃")

    def convert_to_herzfeld():
        try:
            d11 = float(d11_entry.get())
            d22 = float(d22_entry.get())
            d33 = float(d33_entry.get())
            deltas = sorted([d11, d22, d33], reverse=True)
            delta_iso = (d11 + d22 + d33) / 3
            omega = deltas[0] - deltas[2]
            kappa = 3 * (d22 - delta_iso) / omega if omega != 0 else 0
            kappa = max(-1, min(1, kappa))
            delta_iso_entry.delete(0, tk.END); delta_iso_entry.insert(0, f"{delta_iso:.2f}")
            omega_entry.delete(0, tk.END); omega_entry.insert(0, f"{omega:.2f}")
            kappa_entry.delete(0, tk.END); kappa_entry.insert(0, f"{kappa:.3f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input for δ₁₁, δ₂₂, or δ₃₃")

    def haeberlen_to_standard():
        try:
            delta_iso = float(iso_entry.get())
            aniso = float(aniso_entry.get())
            eta = float(eta_entry.get())
            delta_zz = delta_iso + aniso
            delta_xx_plus_delta_yy = 3 * delta_iso - delta_zz
            delta_yy_minus_delta_xx = eta * aniso
            delta_xx = (delta_xx_plus_delta_yy - delta_yy_minus_delta_xx) / 2
            delta_yy = (delta_xx_plus_delta_yy + delta_yy_minus_delta_xx) / 2
            deltas = sorted([delta_zz, delta_xx, delta_yy], reverse=True)
            d11 = deltas[0]
            d22 = deltas[1]
            d33 = deltas[2]
            d11_entry.delete(0, tk.END); d11_entry.insert(0, f"{d11:.2f}")
            d22_entry.delete(0, tk.END); d22_entry.insert(0, f"{d22:.2f}")
            d33_entry.delete(0, tk.END); d33_entry.insert(0, f"{d33:.2f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input for δ_iso, aniso, or η")

    def haeberlen_to_herzfeld():
        try:
            delta_iso = float(iso_entry.get())
            aniso = float(aniso_entry.get())
            eta = float(eta_entry.get())
            delta_zz = delta_iso + aniso
            delta_yy_minus_delta_xx = eta * aniso
            delta_xx_plus_delta_yy = 3 * delta_iso - delta_zz
            delta_xx = (delta_xx_plus_delta_yy - delta_yy_minus_delta_xx) / 2
            delta_yy = (delta_xx_plus_delta_yy + delta_yy_minus_delta_xx) / 2
            deltas = sorted([delta_zz, delta_xx, delta_yy], reverse=True)
            delta_iso_new = (delta_zz + delta_xx + delta_yy) / 3
            omega = deltas[0] - deltas[2]
            kappa = 3 * (delta_yy - delta_iso_new) / omega if omega != 0 else 0
            kappa = max(-1, min(1, kappa))
            delta_iso_entry.delete(0, tk.END); delta_iso_entry.insert(0, f"{delta_iso_new:.2f}")
            omega_entry.delete(0, tk.END); omega_entry.insert(0, f"{omega:.2f}")
            kappa_entry.delete(0, tk.END); kappa_entry.insert(0, f"{kappa:.3f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input for δ_iso, aniso, or η")

    def herzfeld_to_standard():
        try:
            delta_iso = float(delta_iso_entry.get())
            omega = float(omega_entry.get())
            kappa = float(kappa_entry.get())
            d11 = delta_iso + (omega / 2) * (1 - kappa / 3)
            d22 = delta_iso + (omega * kappa / 3)
            d33 = delta_iso - (omega / 2) * (1 + kappa / 3)
            d11_entry.delete(0, tk.END); d11_entry.insert(0, f"{d11:.2f}")
            d22_entry.delete(0, tk.END); d22_entry.insert(0, f"{d22:.2f}")
            d33_entry.delete(0, tk.END); d33_entry.insert(0, f"{d33:.2f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input for δ_iso, Ω, or κ")

    def herzfeld_to_haeberlen():
        try:
            delta_iso = float(delta_iso_entry.get())
            omega = float(omega_entry.get())
            kappa = float(kappa_entry.get())
            # Calculate principal components
            delta_11 = delta_iso + (omega / 2) * (1 - kappa / 3)
            delta_22 = delta_iso + (omega * kappa / 3)
            delta_33 = delta_iso - (omega / 2) * (1 + kappa / 3)
            deltas = [delta_11, delta_22, delta_33]
            # Assign Haeberlen components based on deviation from delta_iso
            deviations = sorted([(abs(d - delta_iso), d) for d in deltas], reverse=True)
            delta_zz = deviations[0][1]
            delta_xx = deviations[2][1]
            delta_yy = deviations[1][1]
            aniso = delta_zz - delta_iso
            eta = abs((delta_yy - delta_xx) / aniso) if aniso != 0 else 0
            eta = max(0, min(1, eta))
            iso_entry.delete(0, tk.END); iso_entry.insert(0, f"{delta_iso:.2f}")
            aniso_entry.delete(0, tk.END); aniso_entry.insert(0, f"{aniso:.2f}")
            eta_entry.delete(0, tk.END); eta_entry.insert(0, f"{eta:.3f}")
        except ValueError:
            messagebox.showerror("Error", "Invalid input for δ_iso, Ω, or κ")

    def reset_all():
        d11_entry.delete(0, tk.END)
        d22_entry.delete(0, tk.END)
        d33_entry.delete(0, tk.END)
        iso_entry.delete(0, tk.END)
        aniso_entry.delete(0, tk.END)
        eta_entry.delete(0, tk.END)
        delta_iso_entry.delete(0, tk.END)
        omega_entry.delete(0, tk.END)
        kappa_entry.delete(0, tk.END)

    tk.Label(shift_window, text="Standard Convention", font=("Arial", 10, "bold")).grid(row=0, column=0, columnspan=2, pady=5)
    tk.Label(shift_window, text="δ₁₁ (ppm):").grid(row=1, column=0, padx=5, pady=2)
    d11_entry = tk.Entry(shift_window)
    d11_entry.grid(row=1, column=1, padx=5, pady=2)
    tk.Label(shift_window, text="δ₂₂ (ppm):").grid(row=2, column=0, padx=5, pady=2)
    d22_entry = tk.Entry(shift_window)
    d22_entry.grid(row=2, column=1, padx=5, pady=2)
    tk.Label(shift_window, text="δ₃₃ (ppm):").grid(row=3, column=0, padx=5, pady=2)
    d33_entry = tk.Entry(shift_window)
    d33_entry.grid(row=3, column=1, padx=5, pady=2)
    tk.Button(shift_window, text="To Haeberlen", command=convert_to_haeberlen).grid(row=4, column=0, padx=5, pady=5)
    tk.Button(shift_window, text="To Herzfeld-Berger", command=convert_to_herzfeld).grid(row=4, column=1, padx=5, pady=5)

    tk.Label(shift_window, text="Haeberlen Convention", font=("Arial", 10, "bold")).grid(row=5, column=0, columnspan=2, pady=5)
    tk.Label(shift_window, text="δ_iso (ppm):").grid(row=6, column=0, padx=5, pady=2)
    iso_entry = tk.Entry(shift_window)
    iso_entry.grid(row=6, column=1, padx=5, pady=2)
    tk.Label(shift_window, text="δ_aniso (ppm):").grid(row=7, column=0, padx=5, pady=2)
    aniso_entry = tk.Entry(shift_window)
    aniso_entry.grid(row=7, column=1, padx=5, pady=2)
    tk.Label(shift_window, text="η:").grid(row=8, column=0, padx=5, pady=2)
    eta_entry = tk.Entry(shift_window)
    eta_entry.grid(row=8, column=1, padx=5, pady=2)
    tk.Button(shift_window, text="To Standard", command=haeberlen_to_standard).grid(row=9, column=0, padx=5, pady=5)
    tk.Button(shift_window, text="To Herzfeld-Berger", command=haeberlen_to_herzfeld).grid(row=9, column=1, padx=5, pady=5)

    tk.Label(shift_window, text="Herzfeld-Berger Convention", font=("Arial", 10, "bold")).grid(row=10, column=0, columnspan=2, pady=5)
    tk.Label(shift_window, text="δ_iso (ppm):").grid(row=11, column=0, padx=5, pady=2)
    delta_iso_entry = tk.Entry(shift_window)
    delta_iso_entry.grid(row=11, column=1, padx=5, pady=2)
    tk.Label(shift_window, text="Ω (ppm):").grid(row=12, column=0, padx=5, pady=2)
    omega_entry = tk.Entry(shift_window)
    omega_entry.grid(row=12, column=1, padx=5, pady=2)
    tk.Label(shift_window, text="κ:").grid(row=13, column=0, padx=5, pady=2)
    kappa_entry = tk.Entry(shift_window)
    kappa_entry.grid(row=13, column=1, padx=5, pady=2)
    tk.Button(shift_window, text="To Standard", command=herzfeld_to_standard).grid(row=14, column=0, padx=5, pady=5)
    tk.Button(shift_window, text="To Haeberlen", command=herzfeld_to_haeberlen).grid(row=14, column=1, padx=5, pady=5)

    tk.Button(shift_window, text="Reset All", command=reset_all).grid(row=15, column=0, columnspan=2, pady=10)
    
# --- Bloch-Siegert Calculator Window ---
def open_bloch_siegert_tool():
    bloch_window = tk.Toplevel(main_window)
    bloch_window.title("Bloch-Siegert Calculator")
    bloch_window.geometry("400x350")

    # Mapping of display fractions to decimal values
    spin_state_map = {
        "1/2": 0.5,
        "1": 1.0,
        "3/2": 1.5,
        "2": 2.0,
        "5/2": 2.5,
        "3": 3.0,
        "7/2": 3.5,
        "4": 4.0,
        "9/2": 4.5
    }

    def calculate_bloch_siegert():
        try:
            # Get Larmor frequencies from entry boxes (in MHz, then convert to Hz)
            wH = float(larmor_H_entry.get())  # MHz
            wX = float(larmor_X_entry.get())  # MHz
            if wH <= 0 or wX <= 0:
                raise ValueError("Larmor frequencies must be positive")
            wH *= 1e6  # Convert to Hz
            wX *= 1e6  # Convert to Hz

            # Get pulse length and spin state (convert display value to decimal)
            pn = float(pulse_length_entry.get())  # in microseconds
            if pn <= 0:
                raise ValueError("Pulse length must be positive")
            i = spin_state_map[spin_state_var.get()]  # Use decimal from mapping

            mode = bloch_mode.get()
            if mode == "null":
                vbs = 1 / (4 * pn * 1e-6)  # Bloch-Siegert shift in Hz (pn to seconds)
            elif mode == "inverted":
                vbs = 1 / (2 * pn * 1e-6)  # Bloch-Siegert shift in Hz (pn to seconds)
            else:
                raise ValueError("Invalid mode selected")

            # Calculate RF field (vx) in kHz
            vx = (wX / wH) * np.sqrt((wH * vbs) * (1 - (wX / wH) ** 2)) / 1e3  # kHz
            if np.isnan(vx) or np.iscomplex(vx):
                raise ValueError("Invalid calculation: check frequency ratio")

            # Calculate 90-degree pulse length (v) in microseconds
            v = 1e6 / ((0.5 + i) * 4 * (vx * 1e3))  # microseconds

            # Update output fields (temporarily enable them)
            rf_field_entry.config(state="normal")
            rf_field_entry.delete(0, tk.END)
            rf_field_entry.insert(0, f"{vx:.3f}")
            rf_field_entry.config(state="disabled")

            pulse_90_entry.config(state="normal")
            pulse_90_entry.delete(0, tk.END)
            pulse_90_entry.insert(0, f"{v:.3f}")
            pulse_90_entry.config(state="disabled")

        except ValueError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Unexpected error: {str(e)}")

    def reset_all():
        larmor_H_entry.delete(0, tk.END)
        larmor_X_entry.delete(0, tk.END)
        pulse_length_entry.delete(0, tk.END)
        spin_state_var.set("1/2")  # Reset to 1/2
        rf_field_entry.config(state="normal")
        rf_field_entry.delete(0, tk.END)
        rf_field_entry.config(state="disabled")
        pulse_90_entry.config(state="normal")
        pulse_90_entry.delete(0, tk.END)
        pulse_90_entry.config(state="disabled")
        bloch_mode.set("null")

    # Mode selection
    bloch_mode = tk.StringVar(value="null")
    tk.Radiobutton(bloch_window, text="Null Signal", variable=bloch_mode, value="null").grid(row=0, column=0, padx=5, pady=5)
    tk.Radiobutton(bloch_window, text="Inverted Signal", variable=bloch_mode, value="inverted").grid(row=0, column=1, padx=5, pady=5)

    # Input fields
    tk.Label(bloch_window, text="Observed Larmor Freq (MHz):").grid(row=1, column=0, padx=5, pady=5)
    larmor_H_entry = tk.Entry(bloch_window)
    larmor_H_entry.grid(row=1, column=1, padx=5, pady=5)

    tk.Label(bloch_window, text="X Larmor Freq (MHz):").grid(row=2, column=0, padx=5, pady=5)
    larmor_X_entry = tk.Entry(bloch_window)
    larmor_X_entry.grid(row=2, column=1, padx=5, pady=5)

    tk.Label(bloch_window, text="Pulse Length (µs):").grid(row=3, column=0, padx=5, pady=5)
    pulse_length_entry = tk.Entry(bloch_window)
    pulse_length_entry.grid(row=3, column=1, padx=5, pady=5)

    tk.Label(bloch_window, text="X Spin State:").grid(row=4, column=0, padx=5, pady=5)
    spin_state_var = tk.StringVar(value="1/2")
    spin_state_dropdown = ttk.Combobox(bloch_window, textvariable=spin_state_var, 
                                       values=["1/2", "1", "3/2", "2", "5/2", "3", "7/2", "4", "9/2"], 
                                       state="readonly")
    spin_state_dropdown.grid(row=4, column=1, padx=5, pady=5)

    # Output fields (disabled by default)
    tk.Label(bloch_window, text="RF Field (kHz):").grid(row=5, column=0, padx=5, pady=5)
    rf_field_entry = tk.Entry(bloch_window, state="disabled")
    rf_field_entry.grid(row=5, column=1, padx=5, pady=5)

    tk.Label(bloch_window, text="90° Pulse Length (µs):").grid(row=6, column=0, padx=5, pady=5)
    pulse_90_entry = tk.Entry(bloch_window, state="disabled")
    pulse_90_entry.grid(row=6, column=1, padx=5, pady=5)

    # Buttons
    tk.Button(bloch_window, text="Calculate", command=calculate_bloch_siegert).grid(row=7, column=1, padx=5, pady=10)
    tk.Button(bloch_window, text="Reset All", command=reset_all).grid(row=8, column=1, padx=5, pady=10)
    
# --- RF Field Calculator Window ---
def open_rf_field_tool():
    rf_window = tk.Toplevel(main_window)
    rf_window.title("RF Field Calculator")
    rf_window.geometry("350x250")

    def calculate_rf_field():
        try:
            # Get input values from boxes
            rf1_input = rf1_entry.get().strip()
            rf2_input = rf2_entry.get().strip()
            w1_input = w1_entry.get().strip()
            w2_input = w2_entry.get().strip()

            # Count how many fields are provided
            provided = sum(1 for x in [rf1_input, rf2_input, w1_input, w2_input] if x)

            if provided != 3:
                raise ValueError("Provide exactly three of: First RF, Second RF, First Power, Second Power")

            # Convert provided inputs to floats
            if rf1_input:
                rf1 = float(rf1_input)
                if rf1 <= 0:
                    raise ValueError("First RF must be positive")
            if rf2_input:
                rf2 = float(rf2_input)
                if rf2 <= 0:
                    raise ValueError("Second RF must be positive")
            if w1_input:
                w1 = float(w1_input)
                if w1 < 0:
                    raise ValueError("First Power must be non-negative")
            if w2_input:
                w2 = float(w2_input)
                if w2 < 0:
                    raise ValueError("Second Power must be non-negative")

            # Calculate the missing value based on w1/w2 = (rf1/rf2)^2
            if not rf1_input:  # Calculate rf1
                rf2 = float(rf2_input)
                w1 = float(w1_input)
                w2 = float(w2_input)
                s2 = w1 / w2
                s = np.sqrt(s2)
                rf1 = rf2 * s
                rf1_entry.delete(0, tk.END)
                rf1_entry.insert(0, f"{rf1:.3f}")
                print(f"Calculated: rf1 = {rf1:.3f} Hz")
            elif not rf2_input:  # Calculate rf2
                rf1 = float(rf1_input)
                w1 = float(w1_input)
                w2 = float(w2_input)
                s2 = w1 / w2
                s = np.sqrt(s2)
                rf2 = rf1 / s
                rf2_entry.delete(0, tk.END)
                rf2_entry.insert(0, f"{rf2:.3f}")
                print(f"Calculated: rf2 = {rf2:.3f} Hz")
            elif not w1_input:  # Calculate w1
                rf1 = float(rf1_input)
                rf2 = float(rf2_input)
                w2 = float(w2_input)
                s = rf1 / rf2
                s2 = s * s
                w1 = w2 * s2
                w1_entry.delete(0, tk.END)
                w1_entry.insert(0, f"{w1:.3f}")
                print(f"Calculated: w1 = {w1:.3f} W")
            elif not w2_input:  # Calculate w2
                rf1 = float(rf1_input)
                rf2 = float(rf2_input)
                w1 = float(w1_input)
                s = rf1 / rf2
                s2 = s * s
                w2 = w1 / s2
                w2_entry.delete(0, tk.END)
                w2_entry.insert(0, f"{w2:.3f}")
                print(f"Calculated: w2 = {w2:.3f} W")

        except ValueError as e:
            messagebox.showerror("Error", str(e))
        except Exception as e:
            messagebox.showerror("Error", f"Unexpected error: {str(e)}")

    def reset_all():
        rf1_entry.delete(0, tk.END)
        rf2_entry.delete(0, tk.END)
        w1_entry.delete(0, tk.END)
        w2_entry.delete(0, tk.END)

    # Input/Output fields
    tk.Label(rf_window, text="First Radiofrequency (Hz):").grid(row=0, column=0, padx=5, pady=5)
    rf1_entry = tk.Entry(rf_window)
    rf1_entry.grid(row=0, column=1, padx=5, pady=5)

    tk.Label(rf_window, text="Second Radiofrequency (Hz):").grid(row=1, column=0, padx=5, pady=5)
    rf2_entry = tk.Entry(rf_window)
    rf2_entry.grid(row=1, column=1, padx=5, pady=5)

    tk.Label(rf_window, text="First Power (W):").grid(row=2, column=0, padx=5, pady=5)
    w1_entry = tk.Entry(rf_window)
    w1_entry.grid(row=2, column=1, padx=5, pady=5)

    tk.Label(rf_window, text="Second Power (W):").grid(row=3, column=0, padx=5, pady=5)
    w2_entry = tk.Entry(rf_window)
    w2_entry.grid(row=3, column=1, padx=5, pady=5)

    # Buttons
    tk.Button(rf_window, text="Calculate", command=calculate_rf_field).grid(row=4, column=1, padx=5, pady=10)
    tk.Button(rf_window, text="Reset All", command=reset_all).grid(row=5, column=1, padx=5, pady=10)
    
# --- Reference Frequencies Tool (IUPAC Recommendations 2001) ---
def open_reference_frequencies_tool():
    ref_window = tk.Toplevel(main_window)
    ref_window.title("IUPAC Reference Frequencies")
    ref_window.geometry("350x250")

    # Frequency ratios in MHz
    FREQUENCY_RATIOS = {
  "1H": 100,
"2H": 15.350609,
"3He": 76.179424,
"6Li": 14.716147,
"7Li": 38.863797,
"9Be": 14.051825,
"10B": 10.743658,
"11B": 32.083974,
"13C": 25.14502,
"14N": 7.226317,
"15N": 10.136767,
"17O": 13.556429,
"19F": 94.094011,
"21Ne": 7.894126,
"23Na": 26.4519,
"25Mg": 6.118849,
"27Al": 26.056859,
"29Si": 8.465314,
"31P": 40.480747,
"33S": 7.676033,
"35Cl": 9.797909,
"37Cl": 8.155391,
"39K": 4.666345,
"40K": 5.806311,
"41K": 2.561281,
"43Ca": 6.730093,
"45Sc": 24.291745,
"47Ti": 5.637277,
"49Ti": 5.639134,
"51V": 26.302873,
"53Cr": 5.652001,
"55Mn": 24.789011,
"57Fe": 3.237778,
"59Co": 23.727074,
"61Ni": 8.936,
"63Cu": 26.515555,
"65Cu": 28.403567,
"67Zn": 6.256663,
"69Ga": 24.001351,
"71Ga": 30.496008,
"73Ge": 3.488315,
"75As": 17.122614,
"77Se": 19.071513,
"79Br": 25.05398,
"81Br": 27.006572,
"83Kr": 3.847614,
"85Rb": 9.654591,
"87Rb": 32.720133,
"87Sr": 4.333379,
"89Y": 4.90005,
"91Zr": 9.386566,
"93Nb": 24.476692,
"95Mo": 6.516655,
"97Mo": 6.653054,
"99Ru": 4.604071,
"101Ru": 5.162029,
"103Rh": 3.186447,
"105Pd": 4.58614,
"107Ag": 4.047614,
"109Ag": 4.653783,
"111Cd": 21.21457,
"113Cd": 22.193175,
"113In": 21.910582,
"115In": 21.912653,
"117Sn": 35.632295,
"119Sn": 37.290665,
"121Sb": 23.930265,
"123Sb": 12.991345,
"123Te": 26.168297,
"125Te": 31.548406,
"127I": 20.007467,
"129Xe": 27.810468,
"131Xe": 8.243394,
"133Cs": 13.116255,
"135Ba": 9.949615,
"137Ba": 11.232013,
"139La": 14.125618,
"141Pr": 30.62,
"143Nd": 5.45,
"145Nd": 3.36,
"147Sm": 4.17,
"149Sm": 3.44,
"151Eu": 24.86,
"153Eu": 10.98,
"155Gd": 3.07,
"157Gd": 4.03,
"159Tb": 24.04,
"161Dy": 3.44,
"163Dy": 4.82,
"165Ho": 21.34,
"167Er": 2.88,
"169Tm": 8.29,
"171Yb": 17.499306,
"173Yb": 4.821,
"175Lu": 11.404,
"176Lu": 8.131,
"177Hf": 4.015,
"179Hf": 2.51,
"181Ta": 12.003087,
"183W": 4.166159,
"185Re": 22.687466,
"187Re": 22.867766,
"187Os": 2.304362,
"189Os": 7.803362,
"191Ir": 1.793054,
"193Ir": 1.938614,
"195Pt": 21.496786,
"197Au": 1.729391,
"199Hg": 17.910879,
"201Hg": 6.611135,
"203Tl": 57.633818,
"205Tl": 57.983785,
"207Pb": 20.920597,
"235U": 1.8414

    }

    def update_frequency_ratio():
        """Update the frequency ratio display when nucleus is selected."""
        nucleus = nucleus_var.get()
        freq_ratio = FREQUENCY_RATIOS.get(nucleus, 0.0)
        freq_ratio_entry.config(state="normal")
        freq_ratio_entry.delete(0, tk.END)
        freq_ratio_entry.insert(0, f"{freq_ratio:.6f}")
        freq_ratio_entry.config(state="disabled")

    def calculate_reference_frequency():
        """Calculate the reference frequency based on SF and frequency ratio."""
        try:
            sf = float(sf_entry.get())
            nucleus = nucleus_var.get()
            freq_ratio = FREQUENCY_RATIOS.get(nucleus, 0.0)
            ref_freq = sf * (freq_ratio / 100.0)
            ref_freq_entry.delete(0, tk.END)
            ref_freq_entry.insert(0, f"{ref_freq:.6f}")
        except ValueError:
            messagebox.showerror("Error", "Please enter a valid SF value")
            ref_freq_entry.delete(0, tk.END)
            ref_freq_entry.insert(0, "0.000000")

    def reset_all():
        """Reset all input and output fields."""
        sf_entry.delete(0, tk.END)
        nucleus_var.set("1H")
        freq_ratio_entry.config(state="normal")
        freq_ratio_entry.delete(0, tk.END)
        freq_ratio_entry.insert(0, f"{FREQUENCY_RATIOS['1H']:.6f}")
        freq_ratio_entry.config(state="disabled")
        ref_freq_entry.delete(0, tk.END)
        ref_freq_entry.insert(0, "0.000000")

    # SF input
    tk.Label(ref_window, text="Enter SF for 1H here (MHz):").grid(row=0, column=0, padx=5, pady=5, sticky="w")
    sf_entry = tk.Entry(ref_window)
    sf_entry.grid(row=0, column=1, padx=5, pady=5)

    # Nucleus selection
    tk.Label(ref_window, text="Select Nucleus:").grid(row=1, column=0, padx=5, pady=5, sticky="w")
    nucleus_var = tk.StringVar(value="1H")
    nucleus_dropdown = ttk.Combobox(ref_window, textvariable=nucleus_var, 
                                   values=["1H",
"2H",
"3He",
"6Li",
"7Li",
"9Be",
"10B",
"11B",
"13C",
"14N",
"15N",
"17O",
"19F",
"21Ne",
"23Na",
"25Mg",
"27Al",
"29Si",
"31P",
"33S",
"35Cl",
"37Cl",
"39K",
"40K",
"41K",
"43Ca",
"45Sc",
"47Ti",
"49Ti",
"51V",
"53Cr",
"55Mn",
"57Fe",
"59Co",
"61Ni",
"63Cu",
"65Cu",
"67Zn",
"69Ga",
"71Ga",
"73Ge",
"75As",
"77Se",
"79Br",
"81Br",
"83Kr",
"85Rb",
"87Rb",
"87Sr",
"89Y",
"91Zr",
"93Nb",
"95Mo",
"97Mo",
"99Ru",
"101Ru",
"103Rh",
"105Pd",
"107Ag",
"109Ag",
"111Cd",
"113Cd",
"113In",
"115In",
"117Sn",
"119Sn",
"121Sb",
"123Sb",
"123Te",
"125Te",
"127I",
"129Xe",
"131Xe",
"133Cs",
"135Ba",
"137Ba",
"139La",
"141Pr",
"143Nd",
"145Nd",
"147Sm",
"149Sm",
"151Eu",
"153Eu",
"155Gd",
"157Gd",
"159Tb",
"161Dy",
"163Dy",
"165Ho",
"167Er",
"169Tm",
"171Yb",
"173Yb",
"175Lu",
"176Lu",
"177Hf",
"179Hf",
"181Ta",
"183W",
"185Re",
"187Re",
"187Os",
"189Os",
"191Ir",
"193Ir",
"195Pt",
"197Au",
"199Hg",
"201Hg",
"203Tl",
"205Tl",
"207Pb",
"235U"
], state="readonly")
    nucleus_dropdown.grid(row=1, column=1, padx=5, pady=5)
    nucleus_dropdown.bind("<<ComboboxSelected>>", lambda event: update_frequency_ratio())

    # Frequency ratio display (disabled)
    tk.Label(ref_window, text="Frequency Ratio:").grid(row=2, column=0, padx=5, pady=5, sticky="w")
    freq_ratio_entry = tk.Entry(ref_window, state="disabled")
    freq_ratio_entry.grid(row=2, column=1, padx=5, pady=5)
    # Initialize with 1H frequency ratio
    freq_ratio_entry.config(state="normal")
    freq_ratio_entry.insert(0, f"{FREQUENCY_RATIOS['1H']:.6f}")
    freq_ratio_entry.config(state="disabled")

    # Reference frequency output
    tk.Label(ref_window, text="Reference Frequency (MHz):").grid(row=3, column=0, padx=5, pady=5, sticky="w")
    ref_freq_entry = tk.Entry(ref_window)
    ref_freq_entry.grid(row=3, column=1, padx=5, pady=5)
    ref_freq_entry.insert(0, "0.000000")

    # Buttons
    tk.Button(ref_window, text="Calculate", command=calculate_reference_frequency).grid(row=4, column=1, padx=5, pady=10)
    tk.Button(ref_window, text="Reset All", command=reset_all).grid(row=5, column=1, padx=5, pady=10)
    
    # IUPAC Recommendations Label 
    tk.Label(ref_window, text="IUPAC Recommendations 2001").grid(row=6, column=1, padx=5, pady=5)
    
# --- Main Window ---
main_window = tk.Tk()
main_window.title("BenMR")
main_window.geometry("300x350")

tk.Label(main_window, text="Select Tool:").pack(pady=5)
tool_var = tk.StringVar()
tool_dropdown = ttk.Combobox(main_window, textvariable=tool_var, 
                            values=["Dipolar Distance Tool", "Chemical Shift Converter", 
                                    "Bloch-Siegert Calculator", "RF Field Calculator", 
                                    "Reference Frequencies"], 
                            state="readonly")
tool_dropdown.pack(pady=5)
tool_dropdown.set("Dipolar Distance Tool")

def open_tool():
    tool = tool_var.get()
    if tool == "Dipolar Distance Tool":
        open_dipolar_tool()
    elif tool == "Chemical Shift Converter":
        open_shift_tool()
    elif tool == "Bloch-Siegert Calculator":
        open_bloch_siegert_tool()
    elif tool == "RF Field Calculator":
        open_rf_field_tool()
    elif tool == "Reference Frequencies":
        open_reference_frequencies_tool()

tk.Button(main_window, text="Open Tool", command=open_tool).pack(pady=10)

logo_image = None
try:
    if getattr(sys, 'frozen', False):
        script_dir = sys._MEIPASS
    else:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    image_path = os.path.join(script_dir, "NMR_image_modified.png")
    if os.path.exists(image_path):
        image = Image.open(image_path)
        image = image.resize((200, 200), Image.Resampling.LANCZOS)
        logo_image = ImageTk.PhotoImage(image)
        logo_label = tk.Label(main_window, image=logo_image)
        logo_label.pack(pady=5)
    else:
        tk.Label(main_window, text="NMR image file not found").pack(pady=5)
except Exception as e:
    tk.Label(main_window, text=f"Error loading logo: {str(e)}").pack(pady=5)

main_window.logo_image = logo_image

main_window.mainloop()