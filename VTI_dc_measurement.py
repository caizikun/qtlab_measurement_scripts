# keithley measurement script for the VTI setup. 
# Ludo Cornelissen
# Last updated: 12/08/2014

# Use up to four keithleys simultaneously. 

import numpy as np
import qt
import measurement_module_dc as mmdc

################################################################# 
#Get the instruments from qtlab. Naming conventions are:
#Magnet power supply                                 : magnet_supply
#Lockin 1 (SR830)                                    : lockin1
#Lockin 2 (SR830)                                    : lockin2
#Lockin 3 (SR830)                                    : lockin3
#Lockin 4 (SR830)                                    : lockin4
#temperature controller                              : temp_control
#rotating sample holder                              : sample_holder
#Keithley 1                                          : keith1
#Keithley 2                                          : keith2
#Keithley 3                                          : keith3
#################################################################

magnet_supply = qt.instruments.get('magnet_supply')
keith1 = qt.instruments.get('keith1')
# keith2 = qt.instruments.get('keith2')
# lockin3 = qt.instruments.get('lockin3')
# lockin4 = qt.instruments.get('lockin4')
temp_control = qt.instruments.get('temp_control')
# sample_holder = qt.instruments.get('sample_holder')

# --------------------------- x --------------------------------    
#             Data directory and number of keithleys
# --------------------------- x --------------------------------
# Set the directory in which the data will be stored:
qt.config['datadir'] = 'C:/Data/Sander K/04_3T_Tesla/04_A4_Tesla'

# Number of keithleys to use in the measurement:
Nkeithleys = 1

# --------------------------- x ----------------------------
#                   IV meetkast parameters:
# --------------------------- x ----------------------------
# Pre-amp gains of the I/V meetkast (in V/A)
Gain1 = 1.0e0       # Gain for keith1
Gain2 = 1.0e0       # Gain for keith2
Gain3 = 1.0e0       # Gain for keith3
Gain4 = 1.0e0       # Gain for keith4

# Put the gains in a dict to pass onto the measurement class
Gains = {'keith1': Gain1, 'keith2': Gain2, 'keith3': Gain3, 'keith4': Gain4}

# --------------------------- x ----------------------------
#               Measurement time control: delays
# --------------------------- x ----------------------------
# Delay before starting the measurement, give instruments time to settle. Units are seconds.
delay1 = 20             
# Delay after setting a new value during sweep. This delay is waited before taking a datapoint
delay2 = 0.1
                          
# --------------------------- x -------------------------------
#               Measurement class initialization
# --------------------------- x -------------------------------
# The actual measurements are performed here, as the functions from the measure() class are called. 
# This is the class initialization:
m = mmdc.measure( Nkeithleys=Nkeithleys, 
             gains=Gains, 
             delay1=delay1, 
             delay2=delay2,
             plot_mode = 'Raw')
# Valid plot modes are: 'Processed', 'Raw', or 'Off'

# --------------------------- x -------------------------------
#                 Performing the measurements
# --------------------------- x -------------------------------
'''
# Some examples of how to use the measurement functions:
Vvector = np.linspace(0.004, 1.0, 50) - first make a sweep vector, in this case 50 values. 
m.lockin_output_sweep(lockin1, Vvector, trace_retrace=True) - put it in the measurement function you want to use. Sweep lockin output voltage

# Same for a 2D scan:
Avector = np.linspace(-180, 180, 100) - sweep vector for angle
Tvector = np.linspace(4, 200, 15) - sweep vector for temperature
m.angle_vs_temperature_scan(sample_holder, Avector, temp_control, Tvector, trace_retrace=False) - Sweep angle vs temperature

# For an overview of all implemented functions, input: 'run measurement_module_dc.py' in the qtlab command line.

# Happy measuring!
'''

# --------------------------- x -------------------------------
#            I-V Section
# --------------------------- x -------------------------------
# keith1.set_parameter_rate('source_value',0.000001,1000)
# keith2.set_parameter_rate('source_value',0.03,1000)

# Source V measure I
# keith1.set_source_value(0)
# keith1.set_source_mode('VOLT') 
# keith1.set_source_range(2)
# keith1.set_source_rate(0.25)
# keith1.set_output(1)
# keith1.set_source_compliance(1.5e-2)
# keith1.set_sense_mode('CURR')
# keith1.set_autorange(1)
# keith1.set_averaging_count(1)
# min, max = 0, 1.1
# V_sweep = np.linspace(0, max, 201)
# V_sweep = np.hstack((V_sweep, np.linspace(max, min, 201)))
# V_sweep = np.hstack((V_sweep, np.linspace(min, 0, 101)))
# m.keithley_sweep(keith1, 'Voltage', V_sweep, trace_retrace=False)

# Source I measure V
# keith1.set_source_value(0)
# keith1.set_source_mode('CURR')
# keith1.set_source_range(20e-3)
# keith1.set_source_rate(200e-6)
# keith1.set_source_compliance(10)
# keith1.set_output(1)
# keith1.set_sense_mode('VOLT')
# keith1.set_autorange(1)
# keith1.set_averaging_count(1)
# I_sweep = np.linspace(0, 5e-3, 101)
# m.keithley_sweep(keith1, 'Current', I_sweep, trace_retrace=True)

# Vc_sweep1=np.arange(0, 0.26, 0.01)
# Vc_sweep2=np.arange(-0.25,0.26, 0.01)
# Vc_sweep3=np.arange(-0.25,0.01, 0.01)
# Vc = np.append(Vc_sweep1,-Vc_sweep2)
# Vc = np.append(Vc,Vc_sweep3)

# Bsweep1=np.arange(0,1.5,0.05)
# Bsweep2=np.arange(2,4.5,0.5)
# Btotal=np.append(Bsweep1,Bsweep2)

# keith1.set_source_value(60)
# keith2.set_source_value(-2)
# magnet_supply.set_mode('Resistive')
# m.magnet_field_sweep(magnet_supply, Btotal, trace_retrace=True)

# Vg=np.arange(0, 60+0.3, 0.3)

# Ivector = np.arange(-1e-3, 1e-3+5e-6, 5e-6)
# I3 = np.arange(1e-3, 0.0-5e-6, -5e-6)
# Itotal = np.append(I1, Ivector)
# Itotal = np.append(Itotal, I3)


# keith2.set_source_value(-2)
# m.keithley_sweep(keith1, 'Voltage', Vg, trace_retrace=True)
# keith2.set_source_value(-1)
# m.keithley_sweep(keith1, 'Voltage', Vg, trace_retrace=True)
# keith2.set_source_value(1)
# m.keithley_sweep(keith1, 'Voltage', Vg, trace_retrace=True)
# keith2.set_source_value(2)
# m.keithley_sweep(keith1, 'Voltage', Vg, trace_retrace=True)


# magnet_supply.set_field(1.0)
# magnet_supply.set_mode('Persistent')
# m.keithley_sweep(keith1, 'Voltage', V1, trace_retrace=True)
# --------------------------- x -------------------------------
#            Magnetic Field sweep section
# --------------------------- x -------------------------------

# --------------- DC Keithley SMU settings for B sweep ---------------

# keith1.set_averaging_count(60)
# keith1.set_averaging_mode('Repeat')
# keith1.set_parameter_rate('source_value',1e-5,100)

# --------------- Field Sweep Vector definition ---------------
# Large field scan
Bsweep1 = np.arange(-2.5,-2.00+0.5,+0.5)
Bsweep2 = np.arange(-1.75,-1.00+0.25,+0.25)
Bsweep3 = np.arange(-0.8,0.8+0.1,+0.1)
Bsweep4 = np.arange(1,1.75+0.25,+0.25)
Bsweep5 = np.arange(2,2.5+0.5,+0.5)
Btotal = np.append(Bsweep1,Bsweep2)
Btotal = np.append(Btotal,Bsweep3)
Btotal = np.append(Btotal,Bsweep4)
Btotal = np.append(Btotal,Bsweep5)

# # Detailed small field scan
# Bsweep1 = np.linspace(-1, 1, 21)
# Bsweep2 = np.linspace(1.8, -2, 20)
# Btotal3 = np.append(Bsweep1,Bsweep2)

#Detailed small field scan
Bsweep1 = np.arange(-1.4,-0.40+0.1,+0.1)
Bsweep2 = np.arange(-0.35,0.35+0.025,+0.025)
Bsweep3 = np.arange(0.4,1.4+0.1,+0.1)
Btotal2 = np.append(Bsweep1,Bsweep2)
Btotal2 = np.append(Btotal2,Bsweep3)

#---- Ultra fine Field sweep ----
# Bsweep1 = np.arange(-0.2,-0.10,+0.05)
# Bsweep2 = np.arange(-0.1,0.1,+0.01)
# Bsweep3 = np.arange(0.1,0.2+0.05,+0.05)
# Btotal = np.append(Bsweep1,Bsweep2)
# Btotal = np.append(Btotal,Bsweep3)
# --------------- Magnetic field measurements ---------------

# Initialisation
# keith1.set_source_value(0)
# keith1.set_source_mode('CURR') #Source current
# keith1.set_source_rate(150e-6)
# keith1.set_source_compliance(1.3)
# keith1.set_output(1)
# keith1.set_sense_mode('VOLT') # Sense voltage
keith1.set_autorange(1)
# keith1.set_source_range(20e-3)
keith1.set_averaging_mode(2)
keith1.set_averaging_count(60)

# Measurement
# Hanle -1deg (But this is often wrongly mentioned., check logbook for Bper or Bpar)

keith_current = [1.5e-2]
keith_range = [2e-2]
# SET DELAY TO 1.5s
for n, i in enumerate(keith_current):
    keith1.set_source_range(keith_range[n])
    keith1.set_source_value(i)
    qt.msleep(1)
    # sample_holder.set_position(0)
    magnet_supply.set_heater('on')
    qt.msleep(10)
    m.repeated_magnet_sweep(magnet_supply, Btotal2,6, trace_retrace=True)
    # sample_holder.set_position(90)
    # m.repeated_magnet_sweep(magnet_supply, Btotal2,4, trace_retrace=True)
    magnet_supply.set_field(0.0)
    magnet_supply.set_heater('off')

    
# keith1.set_source_range(2e-2)
# keith1.set_source_value(5.4e-3)  
# qt.msleep(1)
# magnet_supply.set_heater('on')
# qt.msleep(10)
# m.repeated_magnet_sweep(magnet_supply, Btotal2,4, trace_retrace=True)
# magnet_supply.set_field(0.0)
# magnet_supply.set_heater('off') 

# sample_holder.set_position(90)    
# magnet_supply.set_heater('on')
# qt.msleep(10)
# m.repeated_magnet_sweep(magnet_supply, Btotal2,4, trace_retrace=True)
# magnet_supply.set_field(0.0)
# magnet_supply.set_heater('off') 

# keith1.set_source_value(0.0045)
# qt.msleep(60)
# magnet_supply.set_heater('on')
# qt.msleep(20)
# m.repeated_magnet_sweep(magnet_supply, Btotal,8, trace_retrace=True)
# magnet_supply.set_field(0.0)
# magnet_supply.set_heater('off')

# Angle set to 92.5 deg for ip-plane measurement
# sample_holder.set_position(92.5)
# keith1.set_source_value(-1e-3)
# qt.msleep(300)
# magnet_supply.set_heater('on')
# qt.msleep(20)
# m.repeated_magnet_sweep(magnet_supply, Btotal2,4, trace_retrace=True)


# m.repeated_magnet_sweep(magnet_supply, Btotal,2, trace_retrace=True)

# --------------------------- x -------------------------------
#            Return instruments to zero if needed
# --------------------------- x -------------------------------
# keith1.set_source_value(0)
# keith2.set_source_value(0)
# Return the field to zero
# magnet_supply.set_mode('Resistive')
# magnet_supply.set_field(0.0)
# magnet_supply.set_mode('persistent')
# magnet_supply.set_heater('off')