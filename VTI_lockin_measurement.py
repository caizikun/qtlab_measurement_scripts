# Lockin measurement script for the VTI setup. 
# Ludo Cornelissen
# Last updated: 12/08/2014

# Use up to four lockins simultaneously. 

import numpy as np
import qt
import measurement_module_ac as mma
reload(mma)

################################################################# 
#Get the instruments from qtlab. Naming conventions are:
#Magnet power supply                                 : magnet_supply
#Lockin 1 (SR830)                                    : lockin1
#Lockin 2 (SR830)                                    : lockin2
#Lockin 3 (SR830)                                    : lockin3
#Lockin 4 (SR830)                                    : lockin4
#temperature controller                              : temp_control
#rotating sample holder                              : sample_holder
#keithley                                            : keith1
#################################################################

magnet_supply = qt.instruments.get('magnet_supply')
lockin1 = qt.instruments.get('lockin1')
lockin2 = qt.instruments.get('lockin2')
# lockin3 = qt.instruments.get('lockin3')
# lockin4 = qt.instruments.get('lockin4')
temp_control = qt.instruments.get('temp_control')
sample_holder = qt.instruments.get('sample_holder')

# --------------------------- x --------------------------------    
#             Data directory and number of lockins
# --------------------------- x --------------------------------
# Set the directory in which the data will be stored:
qt.config['datadir'] = 'C:/Data/Juan/NFO_1'

# Number of lockins to use in the measurement:
Nlockins = 2

# --------------------------- x ----------------------------
#                   IV meetkast parameters:
# --------------------------- x ----------------------------
# Amplification factor of the lockin AC output (in A/V)
VIconv = 100e-6

# Pre-amp gains of the I/V meetkast (in V/A)
Gain1 = 1.0e0       # Gain for lockin1
Gain2 = 1.0e0       # Gain for lockin2
Gain3 = 1.0e0       # Gain for lockin3
Gain4 = 1.0e0       # Gain for lockin4

# Put the gains in a dict to pass onto the measurement class
Gains = {'lockin1': Gain1, 'lockin2': Gain2, 'lockin3': Gain3, 'lockin4': Gain4}
# Set the harmonics of the lockins. Dict to pass on to measurement class
Harmonics = {'lockin1':1, 'lockin2':2, 'lockin3':1, 'lockin4':2}

# --------------------------- x ----------------------------
#               Measurement time control: delays
# --------------------------- x ----------------------------
# Delay before starting the measurement, give instruments time to settle. Units are seconds.
delay1 = 15                                
# Delay after setting the a new value during sweep. This delay is waited before taking a datapoint
# Should be around 3-5*tau for lockin measurements.
delay2 = 10
                          
# --------------------------- x -------------------------------
#               Measurement class initialization
# --------------------------- x -------------------------------
# The actual measurements are performed here, as the functions from the measure() class are called. 
# This is the class initialization:
m = mma.measure( Nlockins=Nlockins, 
             harmonics=Harmonics, 
             VIconv=VIconv, 
             gains=Gains, 
             delay1=delay1, 
             delay2=delay2,
             plot_mode = 'Raw')
# Valid plot modes are: 'Processed', 'Raw', or 'Off'

# set the current which is applied to the sample.
Current_on_sample = 100e-6
V_out = Current_on_sample / VIconv
lockin1.set_amplitude(V_out)

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

# For an overview of all implemented functions, input: 'run measurement_module_ac.py' in the qtlab command line.

# Happy measuring!
'''

Bvector = np.arange(-2, 2, 0.05)
m.magnet_field_sweep(magnet_supply, Bvector, trace_retrace=True)

Bvector = np.arange(-1, 1, 0.02)
m.magnet_field_sweep(magnet_supply, Bvector, trace_retrace=True)

sample_holder.set_position(10)

Bvector = np.arange(-2, 2, 0.05)
m.magnet_field_sweep(magnet_supply, Bvector, trace_retrace=True)

Bvector = np.arange(-1, 1, 0.02)
m.magnet_field_sweep(magnet_supply, Bvector, trace_retrace=True)

# --------------------------- x -------------------------------
#            Return instruments to zero if needed
# --------------------------- x -------------------------------
''' Return lockin excitation to zero (0.004 is lowest possible excitation) '''
lockin1.set_amplitude(0.004)

''' Return the field to zero '''
#magnet_supply.set_mode('Resistive')
magnet_supply.set_field(0.0)
magnet_supply.set_heater('off')