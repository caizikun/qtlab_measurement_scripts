# Module containing measurement functions for use with up to 4 keithleys and 
# a temperature controller. 
# Ludo Cornelissen
# Last updated: 12/08/2014

# Use up to four keithleys for measurements. Implementing more keithleys
# can be done by changing the _take_data and _create_data functions.

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

from numpy import pi, random, arange, size, array, sin, cos, linspace, sinc, sqrt
from shutil import copyfile
from os import mkdir
from os.path import exists
from lib.file_support.spyview import SpyView
import inspect

import qt
import sys
import numpy as np
import data as d
import shutil
import os
import time
import datetime

################################################################# 
#Get the instruments from qtlab. It is very important that intrument
#names are not changed! Naming conventions are:
#Keithley 1                                   : keith1
#Keithley 2                                   : keith2
#Keithley 3                                   : keith3
#Keithley 4                                   : keith4
#temperature controller                       : temp_control
#
#Deviation from these conventions will cause the script to fail!
#################################################################
# INSTUMENTS are defined here!

try:    
    keith1 = qt.instruments.get('keith1')
except:
    keith1 = None
try:
    keith2 = qt.instruments.get('keith2')
except:
    keith2 = None
try:
    keith3 = qt.instruments.get('keith3')
except:
    keith3 = None
try:
    keith4 = qt.instruments.get('keith4')
except:
    keith4 = None
try:
    temp_control = qt.instruments.get('temp_control')
except:
    temp_control = None


# ----------------------------------------------------------------
##################### measurement class ############################
# ----------------------------------------------------------------
class measure():
    '''
    The main class which comprehends all the measurement data handling,
    as well as measurement functions that can be used in the script.
    '''
    
    def __init__(self, Nkeithleys, gains, delay1, delay2, plot_mode):
        '''
        Initialize the measure() class.
        Creates an incremental generator for keeping track of data file number.
        Input:
            Nkeithleys (int)          : Number of lockins to use. 1, 2 or 3
            gains (dict)            : dictionary containing the gains for keith1,
                                      keith2, keith3, keith4. With this, the pre-amp gains
                                      are meant.
            delay1 (float)          : delay to wait at the start of a sweep, in seconds
            delay2 (float)          : delay to wait before taking a datapoint during
                                      a sweep. This should be 5-15*tau for a lockin 
                                      measurement, depending on filter settings. Value
                                      here is in seconds.
            plot_mode (string)      : Valid options here are:
                                      "Raw"         : plots the raw voltage measured by the lockins
                                      "Processed"   : plots the resistance measured by the lockins,
                                                      so we already take the gains and IV conversion
                                                      into account.
                                      "Off"         : Use this when no plotting is desired. Can be
                                                      handy when precise timing is important, since
                                                      plotting introduces some delay.
        Output:
            None
        '''
        self.filename='data'
        self.generator=d.IncrementalGenerator(qt.config['datadir']+'\\'+self.filename,1);
        print 'Initializing dc measurements.'
        
        self._Nkeithleys = Nkeithleys
        self._gains = gains
        self._delay1 = delay1
        self._delay2 = delay2
        self._plot_mode = plot_mode
        self._plot_update_count = 0
        
        # For the VTI setup, we should have the temperature controller hooked up
        # to the PC and initialized in the script. check for this using the try-
        # except case below:
        try:
            temp_control
        except NameError:
            temp_control = None
        
        # Get all parameter values for all connected instruments, to write the most up-to-date 
        # settings file as possible
        for ins in qt.instruments.get_instruments().values():
            ins.get_all()
            
        
    def _create_data(self,x_vector,x_instrument,x_parameter,y_vector,y_instrument,y_parameter):
        '''
        Creates the qtlab data object. Also copies the python script,
        generates a spyview meta file and generates a settings file
        which shows the instrument parameters at the start of the 
        measurement.
        '''
        qt.Data.set_filename_generator(self.generator)
        data = qt.Data(name=self.filename)
        
        # We allow for 2D scans, thus adding two parameters to the data.
        data.add_coordinate(x_instrument + ' (' + x_parameter + ')',
                            size=len(x_vector),
                            start=x_vector[0],
                            end=x_vector[-1]) 
        data.add_coordinate(y_instrument +' (' + y_parameter +')',
                            size=len(y_vector),
                            start=y_vector[0],
                            end=y_vector[-1])
        
        # add a fancy trick here so that you can automatically
        # determine to do an IV or VI measurement with the keithley
        self._keithley_sense={}
        for ins in qt.instruments.get_instrument_names():
            if 'keith' in ins:
                sensemod = eval('%s.get_sense_mode()\n' % ins)
                self._keithley_sense.update({ins: sensemod})
            
        
        proc_unit = 'Ohm'
        
        if self._keithley_sense['keith1'] == 0:
            unit = 'V'
        elif self._keithley_sense['keith1'] == 1:
            unit = 'A'
        elif self._keithley_sense['keith1'] == 2:
            unit = 'Ohm'
        
        data.add_value('keith1',
                        size=len(x_vector),
                        units=unit)     # Read out keith1
        data.add_value('keith1_processed',
                        size=len(x_vector),
                        units=proc_unit)   # processed data column    
                        
        if self._Nkeithleys >= 2:
        
            if self._keithley_sense['keith2'] == 0:
                unit = 'V'
            elif self._keithley_sense['keith2'] == 1:
                unit = 'A'
            elif self._keithley_sense['keith2'] == 2:
                unit = 'Ohm'
                
            data.add_value('keith2',
                        size=len(x_vector),
                        units=unit)     # Read out keith2
            data.add_value('keith2_processed',
                        size=len(x_vector),
                        units=proc_unit)   # processed data column
        
        if self._Nkeithleys >= 3:
            if self._keithley_sense['keith3'] == 0:
                unit = 'V'
            elif self._keithley_sense['keith3'] == 1:
                unit = 'A'
            elif self._keithley_sense['keith3'] == 2:
                unit = 'Ohm'        
        
            data.add_value('keith3',
                        size=len(x_vector),
                        units=unit)     # Read out lockin 3 X
            data.add_value('keith3_processed',
                        size=len(x_vector),
                        units=proc_unit)   # processed data column                         

        if self._Nkeithleys == 4:
            if self._keithley_sense['keith4'] == 0:
                unit = 'V'
            elif self._keithley_sense['keith4'] == 1:
                unit = 'A'
            elif self._keithley_sense['keith4'] == 2:
                unit = 'Ohm'
                
            data.add_value('keith4',
                        size=len(x_vector),
                        units=unit)     # Read out lockin 4 X
            data.add_value('keith4_processed',
                        size=len(x_vector),
                        units=proc_unit)   # processed data column          
            
        # If the temperature controller is found, measure the temperature for every
        # datapoint.
        if temp_control is not None:
            data.add_value('measured_temperature',
                            size=len(x_vector),
                            units='K') # Read the temperature from the controller.                 
        
        data.add_value('Timestamp', size=len(x_vector), units='s')
        
        data.create_file()                                  # Create data file
        
        script_path = os.path.abspath(inspect.stack()[2][1])
        copy_path = os.path.abspath(qt.config['datadir'])+'\\'+self.filename + '_'+str(self.generator._counter-1) + '_script.py'
        shutil.copyfile(script_path, copy_path)           # Copy the python script into the data folder
        # SpyView(data).write_meta_file()                     # Create the spyview meta.txt file
        print 'Data file ' + self.filename + '_' + str(self.generator._counter-1) + '.dat created in directory ' + os.path.abspath(qt.config['datadir'])
        return data
    
    def _take_data(self):
        ''' 
        Reads data from the keithleys. Use up to four keithleys.
        '''
        
        qt.msleep(self._delay2)
        
        if self._Nkeithleys >= 1:
            X1 = keith1.get_sense_value()
            
            sense = self._keithley_sense['keith1']
            sourceval = keith1.get_source_value()
            gain = self._gains['keith1']
            if sense == 0:
                if sourceval != 0.0:
                    X1p = X1 / (gain * sourceval)
                else:
                    X1p = 0.0
            elif sense == 1:
                if X1 != 0.0:
                    X1p = sourceval/(gain*X1)
                else:
                    X1p = 0.0
            elif sense == 2:
                X1p = X1
            datavalues = [X1, X1p]
            
        if self._Nkeithleys >= 2:
            X2 = keith2.get_sense_value()
            
            sense = self._keithley_sense['keith2']
            sourceval = keith2.get_source_value()
            gain = self._gains['keith2']
            if sense == 0:
                if sourceval != 0.0:
                    X2p = X2 / (gain * sourceval)
                else:
                    X2p = 0.0
            elif sense == 1:
                if X2 != 0.0:
                    X2p = sourceval/(gain*X1)
                else:
                    X2p = 0.0
            elif sense == 2:
                X2p = X2
                
            for item in [X2, X2p]:
                datavalues.append(item)

        if self._Nkeithleys >= 3:
        
            sense = self._keithley_sense['keith3']
            sourceval = keith3.get_source_value()
            gain = self._gains['keith3']
            if sense == 0:
                if sourceval != 0.0:
                    X1p = X1 / (gain * sourceval)
                else:
                    X1p = 0.0
            elif sense == 1:
                if X1 != 0.0:
                    X1p = sourceval/(gain*X1)
                else:
                    X1p = 0.0
            elif sense == 2:
                X1p = X1
            for item in [X3, X3p]:
                datavalues.append(item)
        
        if self._Nkeithleys == 4:
        
            sense = self._keithley_sense['keith4']
            sourceval = keith4.get_source_value()
            gain = self._gains['keith4']
            if sense == 0:
                if sourceval != 0.0:
                    X1p = X1 / (gain * sourceval)
                else:
                    X1p = 0.0
            elif sense == 1:
                if X1 != 0.0:
                    X1p = sourceval/(gain*X1)
                else:
                    X1p = 0.0
            elif sense == 2:
                X1p = X1
            for item in [X4, X4p]:
                datavalues.append(item)
                
        else:
            pass
        
        if temp_control is not None:
            #T = temp_control.get_temperature()
            T = temp_control.get_temperatureA() # for the VTI we have to get temperature A
            datavalues.append(T)
        else:
            pass
        # add a timestamp to the datapoint
        datavalues.append(time.time())
        return datavalues
    
    def _time_output(self, x, x_vector, delay, Nkeithleys):
        '''
        Function to estimate how long one sweep takes, 
        and indicates progress during the measurement.
      
        Input:
            x (float)                     : current sweep value
            xvector (list of floats)      : the sweep to perform
            delay   (float)               : delay time for each datapoint,
                                            in ms
            Nkeithleys (int)              : Number of keithleys used in the
                                            measurement.
        '''
        global t_start
        xvector = x_vector.tolist() 
        Ntot = len(xvector) # Total number of datapoints to take
        if temp_control is not None:
            t_temp = 30.00
        else:
            t_temp = 0.0
        overhead = 100.0+Nkeithleys*80.0+t_temp  #Additional time for each datapoint, in ms
        delay_ms = delay*1000.0
      
        if x==xvector[0]:
            estimate = Ntot*(delay_ms + overhead)
            t_start = time.time()
        elif x==xvector[-1]:
            print "Estimated time remaining: 00:00:00"
            t_elapsed = time.time()-t_start
            m, s = divmod(t_elapsed, 60)
            h, m = divmod(m, 60)
            print "Total measurement time: %d:%02d:%02d" % (h, m, s)
            return
        else:
            estimate = (Ntot - xvector.index(x))*(delay_ms + overhead) 
        seconds = int(estimate/1000.0)
        m, s = divmod(seconds, 60)
        h, m = divmod(m, 60)
        print "Estimated time remaining: %d:%02d:%02d\r" % (h, m, s),
        sys.stdout.flush()
    
    def _create_plots_1D(self, data):
        plotlist = []
        self._plot_update_count = 0
        if self._plot_mode == 'Off':
            return plotlist
        elif self._plot_mode == 'Processed':
            startx = 3
        elif self._plot_mode == 'Raw':
            startx = 2
        else:
            print 'Invalid plot_mode selected. Choose from "Raw", "Processed" or "Off".'
            print 'Defaulting to "Off".'
            return plotlist
            
        inslist = qt.instruments.get_instrument_names()
        counter = 0
        ploty = False
        
        for ins in inslist:
            if 'keith' in ins and (counter+1) <= self._Nkeithleys:
                dim = startx+2*(counter)
                plotr = qt.Plot2D(data, name=ins, coordim=1, valdim=dim)
                plotlist.append(plotr)
                plotr.set_style('linespoints')
                plotr.set_ylabel(data.format_label(dim))
                counter += 1  
        return plotlist
    
    def _update_plots(self, plotlist):
        if self._plot_update_count % 2 == 0:
            for plot in plotlist:
                plot.update()
        self._plot_update_count += 1
	
    # ----------------------------------------------------------
    ################ Measurement functions #####################    
    # ----------------------------------------------------------
    #All measurement functions that can be called in this script
    #are defined here. We start with 1D scans, followed by 2D scans.
    
    # ------------- x --------------
    # 1D Scans
    # ------------- x --------------

    def time_trace(self, Npoints):
        '''
        Takes a trace as a function of time. Amount of datapoints is
        determined by Npoints. Total time is determined by Npoints*delay2
     
        Input:
        Npoints        :       Number of datapoints to record
     
        '''
        qt.mstart()
        tstart = 0
         
        x_vector = np.arange(Npoints)
        y_vector = [0]
        data = self._create_data(x_vector,'time','seconds',
                                y_vector,'none','y_parameter')
        print 'Recording %i datapoints as a function of time.' % Npoints
        
        plotlist = self._create_plots_1D(data)
        qt.msleep(self._delay1)
        
        st_time = time.time()
        
        for x in x_vector:
            t1 = time.time()
            datavalues = self._take_data()
            t2 = time.time()
            tavg = (t1+t2)/2 - st_time
            y = 0.0           
            data.add_data_point(tavg,y, *datavalues )
            self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
            self._update_plots(plotlist)
            
        end_time = time.time()
        data.new_block()
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'time_trace finished.'
        elapsed = end_time-st_time
        print 'Elapsed measurement time = %2.2f seconds.' % elapsed
    
    def lockin_output_sweep(self, lockin, Vvector, trace_retrace=False):
        '''
        Sweeps the output voltage of the lockin.
        
        Input:
            lockin (qtlab.instrument)       : lockin amplifier to sweep
            Vvector (np.array of floats)    : Vector of output voltages to sweep over
            trace_retrace (bool)            : True for trace/retrace, false trace only
        
        Output:
            None
        '''
        
        qt.mstart()
        # Create sweep vector
        x_vector = Vvector
        xstart = Vvector[0]
        xend = Vvector[-1]
        
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        # create data file, spyview metafile, copy script
        data = self._create_data(x_vector,lockin.get_name(),'Output [V]',
                                y_vector,'sweep direction','tr_1/retr_-1')
        print 'Sweeping lockin output voltage from %2.4f to %2.4f V.' % (xstart, xend)
        
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        plotlist = self._create_plots_1D(data)
        lockin.set_amplitude(xstart)
        
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Performing ' + tr_retr_map[y]
        
            for x in x_vector:
                lockin.set_amplitude(x)
                datavalues = self._take_data()
                data.add_data_point(x,y, *datavalues )
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            x_vector = x_vector[::-1]
            data.new_block()
            
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'lockin_output_sweep finished.'
    
    def frequency_sweep(self, lockin, Fvector, trace_retrace=False):
        '''
        Sweeps the frequency of the internal oscillator of lockin1.
        
        Input:
            lockin (qtlab.instrument)       : lockin amplifier to sweep
            Fvector (np.array of floats)    : Vector of frequencies to sweep over
            trace_retrace (bool)            : True for trace/retrace, false trace only
        
        Output:
            None
        '''
        qt.mstart()
        # Create sweep vector
        x_vector = Fvector
        fstart = x_vector[0]
        fend = x_vector[-1]
        
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        # create data file, spyview metafile, copy script
        data = self._create_data(x_vector,lockin.get_name(),'Frequency [Hz]',
                                y_vector,'sweep direction','tr_1/retr_-1')
        print 'Sweeping lockin frequency from %2.4f to %2.4f Hz.' % (fstart, fend)
        
        lockin.set_frequency(fstart)
        plotlist = self._create_plots_1D(data)
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Performing ' + tr_retr_map[y]
        
            for x in x_vector:
                lockin.set_frequency(x)
                datavalues = self._take_data()        
                data.add_data_point(x,y, *datavalues ) 
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            x_vector = x_vector[::-1]
            data.new_block()
        
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'frequency_sweep finished.'
    
    def lockin_current_sweep(self, lockin, Ivector, trace_retrace=False):
        '''
        Sweeps the current applied to the sample, by sweeping the output
        voltage of the lockin.
         
        Input:
            lockin (qtlab.instrument)       : lockin amplifier to sweep
            Ivector (np.array of floats)    : Vector of currents to sweep over
            trace_retrace (bool)            : True for trace/retrace, false trace only
         
        Output:
            None
        '''
        qt.mstart()
 
        # Create sweep vector
        x_vector = Ivector
        xstart = x_vector[0]
        xend = x_vector[-1]
        
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
         
        # create data file, spyview metafile, copy script
        data = self._create_data(x_vector,lockin.get_name(),'Current [A]',
                                y_vector,'sweep direction','tr_1/retr_-1')
        print 'Sweeping on sample current from %2.4e to %2.4e A.' % (xstart, xend)
        
        lockin1.set_amplitude(x_vector[0]/self._VIconv)
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        plotlist = self._create_plots_1D(data)
        
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Performing ' + tr_retr_map[y]            
         
            for x in x_vector:
                # Convert the current values given by the user to voltage values to
                # send to the lockin.
                x_V = x / self._VIconv
                lockin.set_amplitude(x_V)
                datavalues = self._take_data()          
                data.add_data_point(x,y, *datavalues )
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            x_vector = x_vector[::-1]
        
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'lockin_current_sweep finished.'
    
    def magnet_current_sweep(self, supply, Ivector, trace_retrace=True):
        '''
        Sweeps the current of the magnet power supply, 
        which should be named as magnet_supply at the top
        of this script. Units are Amps
        
        This function takes current values as its input. These
        values have to be specified in amps!
        
        Input:
            supply (qtlab.instrument)       : magnet supply to sweep
            Ivector (np.array of floats)    : vector of currents to sweep over
            trace_retrace (bool)            : Take a trace and retrace?
        
        Output:
            None
        '''
        qt.mstart()
        # Create sweep vector
        x_vector = Ivector
        Istart = x_vector[0]
        Iend = x_vector[-1]
        
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        data = self._create_data(x_vector,supply.get_name(),'Current [A]',
                                y_vector,'sweep direction','tr_1/retr_-1') # create data file, spyview metafile, copy script
                                
        print 'Sweeping magnet current from %2.3fA to %2.3fA.' % (Istart, Iend)
        
        supply.set_current(x_vector[0])
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        plotlist = self._create_plots_1D(data)
        
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Performing ' + tr_retr_map[y]

            for x in x_vector:
                supply.set_current(x)
                datavalues = self._take_data()            
                data.add_data_point(x,y, *datavalues)
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            x_vector = x_vector[::-1]
                
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'magnet_current_sweep finished.'
    
    def magnet_field_sweep(self, supply, Bvector, trace_retrace=True):
        '''
        Sweeps the magnet power supply , 
        which should be named as magnet_supply at the top
        of this script.Units are Tesla!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        Input:
            supply (qtlab.instrument)     : magnet supply to sweep
            Bvector (np.array of floats)  : Vector of B values to use
            trace_retrace (bool)          : Take a trace and retrace?
            
        Output:
            None
        '''
        qt.mstart()
        # Create sweep vector
        x_vector = Bvector
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        # initialize data structure
        data = self._create_data(x_vector,supply.get_name(),'Field [T]',
                                y_vector,'sweep direction','tr_1/retr_-1')  # create data file, spyview metafile, copy script
        print 'Sweeping field from %2.3fmT to %2.3fmT.' % (x_vector[0]*1000.0, x_vector[-1]*1000.0)
        
        supply.set_field(x_vector[0])
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        
        plotlist = self._create_plots_1D(data)
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Performing ' + tr_retr_map[y]
            for x in x_vector:
                supply.set_field(x)
                datavalues = self._take_data()
                data.add_data_point(x, y, *datavalues)
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            x_vector = x_vector[::-1]
            
        data.new_block()            
        data._write_settings_file()  # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'magnet_field_sweep finished.'    
    
    def angle_sweep(self, sample_hldr, Avector, trace_retrace=False):
        '''
        Sweeps the angle of the sample with respect to the magnetic field.
                
        Input:
            sample_hldr (qtlab.instrument)   : Sample holder to rotate
            Avector (np.array of floats)     : Vector of angles to sweep over
            trace_retrace (bool)             : if true, takes also a retrace
        
        Output:
            None
        '''
        qt.mstart()
        # Create sweep vector
        x_vector = Avector
        astart = x_vector[0]
        aend = x_vector[-1]
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        data = self._create_data(x_vector,'angle','degrees',
                                y_vector,'sweep direction','tr_1/retr_-1') # create data file, spyview metafile, copy script  
        print 'Sweeping angle from %2.3fdegrees to %2.3fdegrees.' % (astart, aend)
        
        sample_hldr.set_position(astart)
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        plotlist = self._create_plots_1D(data)
        
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Performing ' + tr_retr_map[y]
            for x in x_vector:
                sample_hldr.set_position(x)
                datavalues = self._take_data()            
                data.add_data_point(x,y, *datavalues )
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            x_vector = x_vector[::-1]

        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'angle_sweep finished.'
    
    def temperature_sweep(self, t_control, Tvector, trace_retrace=False):
        '''
        Sweeps the sample temperature.
                
        Input:
            t_control (qtlab.instrument)     : Temperature controller used for sweep
            Tvector (np.array of floats)     : Vector of temperatures to sweep over
            trace_retrace (bool)             : if true, takes also a retrace
            
        
        Output:
            None
        '''        
        qt.mstart()
        # Create sweep vector
        x_vector = Tvector
        Tstart = x_vector[0]
        Tend = x_vector[-1]
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        data = self._create_data(x_vector,'temperature','K',
                                y_vector,'sweep direction','tr_1/retr_-1') # create data file, spyview metafile, copy script
        print 'Sweeping temperature from %2.3fK to %2.3fK.' % (Tstart, Tend)
        
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        plotlist = self._create_plots_1D(data)
        t_control.ramp_temperature(Tstart)
        
        for y in y_vector:
            print 'Performing ' + tr_retr_map[y]
            qt.msleep(self._delay1)
            
            for x in x_vector:
                t_control.ramp_temperature(x)
                datavalues = self._take_data()            
                data.add_data_point(x,y, *datavalues )
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            x_vector=x_vector[::-1]
            
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'temperature_sweep finished.'               
        
    def keithley_sweep(self, keithley, source_parameter, Svector, trace_retrace=False):
        '''
        Sweeps the output of a keithley source meter unit.
                
        Input:
            keithley (qtlab.instrument)      : keithley used for sweep
            Svector (np.array of floats)     : Output values to sweep over
            source_parameter (string)        : 'Current' or 'Voltage', depending on
                                               what the keithley output should be.
            trace_retrace (bool)             : if true, takes also a retrace
            
        Output:
            None
        '''        
        qt.mstart()
        # Create sweep vector
        x_vector = Svector
        Kstart = x_vector[0]
        Kend = x_vector[-1]
        if trace_retrace:
            y_vector = [1.0, -1.0]
        else:
            y_vector = [1.0]
        
        if source_parameter == 'Current':
            keithley.set_source_mode('CURR')
            unit = 'A'
        elif source_parameter == 'Voltage':
            keithley.set_source_mode('VOLT')
            unit = 'V'
        else:
            print 'Invalid source parameter... choose "Current" or "Voltage".'
            return
        
        data = self._create_data(x_vector,source_parameter,unit,
                                y_vector,'sweep direction','tr_1/retr_-1') # create data file, spyview metafile, copy script
        print 'Sweeping keithley %s from %1.2e%s to %1.2e%s.' % (source_parameter, Kstart, unit, Kend, unit)
        
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        plotlist = self._create_plots_1D(data)
        keithley.set_source_value(Kstart)
        
        for y in y_vector:
            print 'Performing ' + tr_retr_map[y]
            qt.msleep(self._delay1)
            
            for x in x_vector:
                keithley.set_source_value(x)
                datavalues = self._take_data()            
                data.add_data_point(x,y, *datavalues )
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            x_vector=x_vector[::-1]
            
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file()
        qt.mend()
        print 'keithley_sweep finished.'         
        
    def manual_angle_sweep(self, trace_vector, retrace_vector=[]):
        '''
        Manually rotate the sample. Function waits for user confirmation
        whether the sample has been rotated to the right position.
        
        Input:
            trace_vector (np.array of floats)   :   vector of angles to rotate to
                                                    for the trace. 
            retrace_vector (np.array of floats) :   vector of angles to rotate to for
                                                    the retrace. If left empty, only
                                                    trace will be performed.
        '''
        qt.mstart()
        # Create sweep vector
        if not retrace_vector:
            y_vector = [1.0]
        else:
            y_vector = [1.0, -1.0]
            
        data = self._create_data(trace_vector,'angle','degrees',
                                y_vector,'sweep direction','tr_1/retr_-1') # create data file, spyview metafile, copy script    
        plotlist = self._create_plots_1D(data)
        tr_retr_map = {1.0 : 'trace', -1.0 : 'retrace'}
        
        for y in y_vector:
            print 'Performing ' + tr_retr_map[y]
            
            for x in trace_vector:
                print 'Rotate sample to %.1f degrees.' % x
                ans = raw_input('Press <enter> when sample is at correct position:')
                print 'Measuring...\r',
                sys.stdout.flush()
                datavalues = self._take_data()            
                data.add_data_point(x,y, *datavalues)
                self._update_plots(plotlist)
        
            if y == -1.0:
            
                print 'Performing ' + tr_retr_map[y]
                for x in retrace_vector:
                    print 'Rotate sample to %.1f degrees.' % x
                    ans = raw_input('Press <enter> when sample is at correct position:')
                    print 'Measuring...\r',
                    sys.stdout.flush()
                    datavalues = self._take_data()            
                    data.add_data_point(x,y, *datavalues)
                    self._update_plots(plotlist)
        
        data._write_settings_file()
        data.close_file()
        qt.mend()
        print 'manual_angle_sweep finished.'
        
    def repeated_magnet_sweep(self, supply, Bvector, repetitions, trace_retrace=False):
        '''
        Sweeps the field of the magnet power supply. Units are Tesla!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        Performs 'repetitions' number of magnet supply sweeps.
        
        If the argument trace_retrace is set to True, (default is False) 
        the sweep is performed back and forth, so using
        repetitions = 10 will result in five forward and five
        backward traces, all stored in the same data file.
        
        Input:
            supply (qtlab.instrument)        : magnet supply to sweep
            Bvector (np.array of floats)     : Vector of field values to sweep over
            repetitions (int)                : Number of sweeps to perform
            trace_retrace (bool)             : Sweep the magnet back and forth. If
                                               False, direction for all sweeps will
                                               be the same.
        Output:
            None
        '''    
        qt.mstart()
        
        # Create sweep vectors
        x_vector = Bvector
        Bstart = x_vector[0]
        Bend = x_vector[-1]
        
        y_vector = np.arange(1, repetitions + 1)
        data = self._create_data(x_vector,supply.get_name(),'Field [T]',
                                y_vector,'repetitons','n') # create data file, spyview metafile, copy script
        
        print 'Sweeping field from %2.3fmT to %2.3fmT. %i sweeps.' % (Bstart*1000, Bend*1000, repetitions)
        supply.set_field(x_vector[0])
        plotlist = self._create_plots_1D(data)
        
        for y in y_vector:
            qt.msleep(self._delay1)
            print 'Field sweep number %i.' % y
            for x in x_vector:
                supply.set_field(x)
                datavalues = self._take_data()
                data.add_data_point(x, y, *datavalues)
                self._time_output(x, x_vector, self._delay2, self._Nkeithleys)
                self._update_plots(plotlist)
            data.new_block()
            if trace_retrace:
                x_vector = x_vector[::-1]
            else:
                supply.set_field(x_vector[0])
                
        data._write_settings_file() # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'repeated_magnet_sweep finished.'

    # ------------- x --------------
    # 2D Scans
    # ------------- x --------------
    
    def magnet_vs_angle_scan(self, supply, Bvector, sample_hldr, Avector, trace_retrace=False):
        '''
        Sweeps the field of the magnet power supply, vs angle of the sample holder. Units are Tesla!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        This function performs a magnetic field sweep for every angle in Avector.
        
        Generates a 2d data file.
        
        Input:
            supply (qtlab.instrument)        : magnet supply to sweep
            Bvector (np.array of floats)     : Vector of field values to sweep over
            sample_hldr (qtlab.instrument)   : sample_holder to sweep angle
            Avector (np.array of floats)     : Vector of angles to sweep over
            trace_retrace (bool)             : Sweep the lockin back and forth. If
                                               False, direction for all sweeps will
                                               be the same.
                                      
        Output:
            None
        '''    
        qt.mstart()
        # Create sweep vectors
        B_vector = Bvector
        a_vector = Avector
        
        data = self._create_data(B_vector, supply.get_name(), 'Field [T]',
                                a_vector,sample_hldr.get_name(),'angle [degrees]')                                # create data file, spyview metafile, copy script
        print 'Scanning field vs angle. Recording %i datapoints in total.' % (len(B_vector)*len(a_vector))

        plotlist = self._create_plots_1D(data)
        
        for a in a_vector:
            supply.set_field(B_vector[0])
            sample_hldr.set_position(a)
            print 'Angle a = %2.2f degrees. Sweeping field.' %  (a)
            qt.msleep(self._delay1)
            
            for B in B_vector:
                supply.set_field(B)
                datavalues = self._take_data()
                data.add_data_point(a, B, *datavalues)
                self._update_plots(plotlist)
            data.new_block()
            
            if trace_retrace:
                B_vector_retr = B_vector[::-1]
                for B in B_vector_retr:
                    supply.set_field(B)
                    datavalues = self._take_data()
                    data.add_data_point(a, B, *datavalues)
                    self._update_plots(plotlist)
                data.new_block()
            else:
                pass
                
        data._write_settings_file()   # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'magnet_vs_angle_scan finished.'

    def angle_vs_magnet_scan(self, sample_hldr, Avector, supply, Bvector, trace_retrace=False, switch_persistent=False):
        '''
        Sweeps the angle of the sample holder vs the field of the magnet power supply. Units are Tesla!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        This function performs an angle sweep for every field value in Bvector. The angle sweeps
        will be performed with the magnet in persistent mode: switching to persistent mode takes
        approx 2 minutes to stabilize.
        
        Generates a 2d data file. If trace_retrace = True, will do a trace/retrace at all field values.
        
        Input:
            sample_hldr (qtlab.instrument)   : sample_holder to sweep angle
            Avector (np.array of floats)     : Vector of angles to sweep over
            supply (qtlab.instrument)        : magnet supply to sweep
            Bvector (np.array of floats)     : Vector of field values to sweep over
            trace_retrace (bool)             : Sweep the lockin back and forth. If
                                               False, direction for all sweeps will
                                               be the same.
            switch_persistent (bool)         : If true, switches the magnet to pers-
                                               istent mode before sweeping the angle.
                                               If false, magnet switch heater is left
                                               on during the angle sweep. Use true to
                                               conserve Helium.
                                      
        Output:
            None
        '''    
        qt.mstart()
        # Create sweep vectors
        B_vector = Bvector
        a_vector = Avector
        
        data = self._create_data(a_vector,sample_hldr.get_name(),'angle [degrees]',
                                 B_vector, supply.get_name(), 'Field [T]') # create data file, spyview metafile, copy script
        print 'Scanning angle vs field. Recording %i datapoints in total.' % (len(B_vector)*len(a_vector))
        
        plotlist = self._create_plots_1D(data)
        
        for B in B_vector:
            if switch_persistent:
                supply.set_mode('Resistive')
                supply.set_field(B)
                supply.set_mode('Persistent')
            else:
                supply.set_field(B)
            print 'Field B = %2.3fT. Sweeping angle.' %  (B)
            sample_hldr.set_position(a_vector[0])
            qt.msleep(self._delay1)
            
            for a in a_vector:
                sample_hldr.set_position(a)
                datavalues = self._take_data()
                data.add_data_point(a, B, *datavalues)
                self._update_plots(plotlist)
            data.new_block()
            
            if trace_retrace:
                a_vector_retr = a_vector[::-1]
                for a in a_vector_retr:
                    sample_hldr.set_position(a)
                    datavalues = self._take_data()
                    data.add_data_point(a, B, *datavalues)
                    self._update_plots(plotlist)
                data.new_block()
            else:
                sample_hldr.set_position(a_vector[0])
                
        data._write_settings_file()   # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'angle_vs_magnet_scan finished.'
        
    def magnet_vs_temperature_scan(self, supply, Bvector, t_control, Tvector, trace_retrace=False, Tramp_parameters=[0.01, 20, 1800]):
        '''
        Sweeps the magnet power supply, 
        which should be named as magnet_supply at the top
        of this script, vs temperature. Units are Tesla!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        Generates a 2d data file.
        
        Input:
            supply (qtlab.instrument)         : magnet supply to sweep
            Bvector (np.array of floats)      : Vector of field values to sweep over
            t_control (qtlab.instrument)      : temperature controller to sweep
            Tvector (np.array of floats)      : Vector of temperatures to sweep over
            trace_retrace (bool)              : Sweep the field back and forth. If
                                                False, direction for all sweeps will
                                                be the same.
            Tramp_parameters (list of floats) : List containing parameters for the tem-
                                                pareture ramping. Parameters are as follows:
                                                [precision, timestep, timeout]. Precision
                                                is the allowed temperature error (default
                                                10 mK). Timestep is the # of seconds of one
                                                step (default 20), Timeout is the # of 
                                                seconds after which the script will 
                                                automatically resume if no stable temperature 
                                                is reached (default is 1800).
                                      
        Output:
            None
        '''    
        qt.mstart()
        # Create sweep vectors
        B_vector = Bvector
        T_vector = Tvector
        
        data = self._create_data(B_vector, supply.get_name(), 'Field [T]',
                                T_vector, t_control.get_name(),'Temperature [K]') # create data file, spyview metafile, copy script
        print 'Scanning field vs temperature. Recording %i datapoints in total.' % (len(B_vector)*len(T_vector))
        plotlist = self._create_plots_1D(data)
        
        for T in T_vector:
            supply.set_field(B_vector[0])
            t_control.ramp_temperature(T, *Tramp_parameters)
            print 'Temperature T = %2.2fK. Sweeping field.' %  (T)
            qt.msleep(self._delay1)
            
            for B in B_vector:
                supply.set_field(B)
                datavalues = self._take_data()
                data.add_data_point(B, T, *datavalues)  
                self._update_plots(plotlist)
            data.new_block()
            
            if trace_retrace:
                B_vector_retr = B_vector[::-1]
                for B in B_vector_retr:
                    supply.set_field(B)
                    datavalues = self._take_data()
                    data.add_data_point(B, T, *datavalues)
                    self._update_plots(plotlist)
                data.new_block()
            else:
                pass
                
        data._write_settings_file()   # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'magnet_vs_temperature_scan finished.'        

    def angle_vs_temperature_scan(self, sample_hldr, Avector, t_control, Tvector, trace_retrace=False, Tramp_parameters=[0.01, 20, 1800]):
        '''
        Sweeps the sample angle of sample_hldr vs temperature.
        
        Generates a 2d data file. If trace_retrace = True, will do a trace/retrace
        at all temperatures.
        
        Input:
            sample_hldr (qtlab.instrument)    : sample holder to sweep
            Avector (np.array of floats)      : Vector of angles to sweep over
            t_control (qtlab.instrument)      : temperature controller to sweep
            Tvector (np.array of floats)      : Vector of temperatures to sweep over
            trace_retrace (bool)              : Sweep the field back and forth. If
                                                False, direction for all sweeps will
                                                be the same.
            Tramp_parameters (list of floats) : List containing parameters for the tem-
                                                pareture ramping. Parameters are as follows:
                                                [precision, timestep, timeout]. Precision
                                                is the allowed temperature error (default
                                                10 mK). Timestep is the # of seconds of one
                                                step (default 20), Timeout is the # of 
                                                seconds after which the script will 
                                                automatically resume if no stable temperature 
                                                is reached (default is 1800).                                               
        Output:
            None
        '''    
        qt.mstart()
        # Create sweep vectors
        A_vector = Avector
        T_vector = Tvector
        
        data = self._create_data(A_vector, sample_hldr.get_name(), 'Angle [degrees]',
                                T_vector, t_control.get_name(),'Temperature [K]') # create data file, spyview metafile, copy script
        print 'Scanning angle vs temperature. Recording %i datapoints in total.' % (len(A_vector)*len(T_vector))
        plotlist = self._create_plots_1D(data)
        
        for T in T_vector:
            sample_hldr.set_position(A_vector[0])
            t_control.ramp_temperature(T, *Tramp_parameters)
            print 'Temperature T = %2.2fK. Sweeping angle.' %  (T)
            qt.msleep(self._delay1)
            
            for A in A_vector:
                sample_hldr.set_position(A)
                datavalues = self._take_data()
                data.add_data_point(A, T, *datavalues)  
                self._update_plots(plotlist)
            data.new_block()
            
            if trace_retrace:
                A_vector_retr = A_vector[::-1]
                for A in A_vector_retr:
                    sample_hldr.set_position(A)
                    datavalues = self._take_data()
                    data.add_data_point(A, T, *datavalues)
                    self._update_plots(plotlist)
                data.new_block()
            else:
                pass
                
        data._write_settings_file()   # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'angle_vs_temperature_scan finished.'         
    
    def magnet_vs_lockin_voltage_scan(self, supply, Bvector, lockin, Vvector, trace_retrace=False):
        '''
        Sweeps the magnet power supply, vs the lockin output voltage of the
        lockin1. Units are Tesla!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        Generates a 2d data file.
        
        Input:
            supply (qtlab.instrument)        : magnet supply to sweep
            Bvector (np.array of floats)     : Vector of field values to sweep over
            lockin (qtlab.instrument)        : lockin amplifier to sweep
            Vvector (np.array of floats)     : Vector of output values to sweep over
            trace_retrace (bool)             : Sweep the lockin back and forth. If
                                               False, direction for all sweeps will
                                               be the same.
                                      
        Output:
            None
        '''    
        qt.mstart()
        # Create sweep vectors
        B_vector = Bvector
        x_vector = Vvector
        
        data = self._create_data(B_vector, supply.get_name(), 'Field [T]',
                                x_vector,lockin.get_name(),'Output voltage [V]')                                # create data file, spyview metafile, copy script
        print 'Scanning field vs lockin amplitude. Recording %i datapoints in total.' % (len(B_vector)*len(x_vector))
        plotlist = self._create_plots_1D(data)
        lockin.set_amplitude(x_vector[0])
        
        for B in B_vector:
            supply.set_field(B)
            qt.msleep(self._delay1)
            print 'Field B= %2.3fmT. Sweeping lockin amplitude.' %  (B*1000.0)
            for x in x_vector:
                lockin.set_amplitude(x)
                datavalues = self._take_data()
                data.add_data_point(B, x, *datavalues)                               
                self._update_plots(plotlist)
            data.new_block()
            if trace_retrace:
                x_vector = x_vector[::-1]
            else:
                pass
                
        data._write_settings_file()                                                                                                             # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'magnet_vs_lockin_voltage_scan finished.'
        
    def magnet_vs_frequency_scan(self, supply, Bvector, lockin, Fvector, trace_retrace=False):
        '''
        Sweeps the magnet power supply vs the frequency of the
        lockin. Units are Tesla and Hertz!!!
        
        This function takes field values as its input. These
        values have to be specified in Tesla!
        
        Generates a 2d data file.
        
        Input:
            supply (qtlab.instrument)        : magnet supply to sweep
            Bvector (np.array of floats)     : Vector of field values to sweep over
            lockin (qtlab.instrument)        : lockin amplifier to sweep
            Fvector (np.array of floats)     : Vector of frequencies to sweep over
            trace_retrace (bool)             : Sweep the lockin back and forth. If
                                               False, direction for all sweeps will
                                               be the same.
                                      
        Output:
            None
        '''    
        qt.mstart()
        # Create sweep vectors
        B_vector = Bvector
        x_vector = Fvector
        
        data = self._create_data(B_vector, supply.get_name(), 'Field [T]',
                                x_vector, lockin.get_name(), 'Frequency [Hz]')                                # create data file, spyview metafile, copy script
        print 'Scanning field vs lockin frequency. Recording %i datapoints in total.' % (len(B_vector)*len(x_vector))
        
        supply.set_field(B_vector[0])
        plotlist = self._create_plots_1D(data)
        
        for x in x_vector:
            lockin.set_frequency(x)
            qt.msleep(self._delay1)
            print 'Frequency F= %2.2fHz. Sweeping field.' % x
            for B in B_vector:
                supply.set_field(B)
                datavalues = self._take_data()
                data.add_data_point(B, x, *datavalues)                               
                self._update_plots(plotlist)
            data.new_block()
            if trace_retrace:
                x_vector = x_vector[::-1]
            else:
                pass
                
        data._write_settings_file()                                                                                                             # Overwrite the settings file created at the beginning, this ensures updating the sweep variable with the latest value
        data.close_file() 
        qt.mend()
        print 'magnet_vs_frequency_scan finished.'
                
if __name__ == '__main__':
    m = measure
    class_objects = dir(m)
    print 'Module containing the measure class.'
    print 'Associated 1D measurement functions are:\n'
    for function in class_objects:
        if function[0] != '_' and not '_vs_' in function:
            print function
    print '\nAssociated 2D measurement functions are:\n'
    for function in class_objects:
        if function[0] != '_' and '_vs_' in function:
            print function
    print '\nInitialize with measure(Nkeithleys, gains, delay1, delay2, plot_mode)'
    print 'Function help is available by calling help(measure.<function name>)'
    
