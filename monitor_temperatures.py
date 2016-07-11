# Temperature monitor script for VTI. Continuously monitors insert and bath temperatures.
# Ludo Cornelissen

import qt
import time
import data as d

temp_control = qt.instruments.get('temp_control')
try: 
    magnet_supply = qt.instruments.get('magnet_supply')
except:
    magnet_supply = None

generator=d.IncrementalGenerator(qt.config['datadir']+'\\'+'temperaturedata',1);
qt.Data.set_filename_generator(generator)

data = qt.Data('')
data.add_coordinate('time (s)')
data.add_value('Tinsert (K)')
data.add_value('Tbath (K)')

plotlist = []
plotlist.append(qt.Plot2D(data, name='Insert temperature', coordim=1, valdim=1))
plotlist.append(qt.Plot2D(data, name='Bath temperature', coordim=1, valdim=2))

if magnet_supply is not None:
    data.add_value('He level (mm)')
    plotlist.append(qt.Plot2D(data, name='He level', coordim=1, valdim=3))

timestep = 15
for plot in plotlist:
    plot.set_style('linespoints')

print 'Monitoring temperature, press <CTRL + c> to stop.'
qt.mstart()
tstart = time.time()

while True:
    try:
        time1 = time.time()
        tavg = time1 - tstart
        Ti = temp_control.get_temperatureA()
        Tb = temp_control.get_temperatureB()
        temp_control.get_heater_output1()
        temp_control.get_heater_output2()
        
        if magnet_supply is not None:
            He = magnet_supply.get_He_level()
            data.add_data_point(tavg, Ti, Tb, He)
        else:
            data.add_data_point(tavg, Ti, Tb)
        
        for plot in plotlist:
            plot.update()

        qt.msleep(timestep)
    except KeyboardInterrupt:
        break
        
tend = time.time()
print 'Monitoring stopped after %.1fs.' % (tend-tstart)
qt.mend()
data.close_file()

for plot in plotlist:
    plot.save_png()