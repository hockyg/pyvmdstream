###############################################################################
#    Copyright (C) 2013  Glen M. Hocky and Aaron S. Keys
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################

"""
        pyvmdstream
        A python library to open a TCP connection to VMD 
        
        based on the C++ vmdstream library by Aaron Keys

"""

import os
import sys
import socket
import subprocess
import time

import numpy as np

# These are hard coded values from vmd version 1.9
VMDSTARTCOLOR=33
VMDNCOLORS=1024

# note the following function comes from matplotlib documentation, Jan 17, 2012
# http://www.scipy.org/Cookbook/Matplotlib/ColormapTransformations
# cmaps shown on http://www.scipy.org/Cookbook/Matplotlib/Show_colormaps
def cmap_discretize(N,cmap="jet"):
    import matplotlib
    from matplotlib import cm
    """Return a discrete colormap from the continuous colormap cmap.
    
        cmap: colormap instance, eg. cm.jet. 
        N: number of colors.
    
    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    if type(cmap) == str:
        cmap = matplotlib.cm.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)

class VMDStream():
    """ 

        This is a class which facilitates controlling VMD via pyvmdstream.

        Though it is possible, and probably faster, to interact with VMD
          directly, a class such as this can make plugging this library into
          existing analysis code very clean.

        Example:
            vmdstream = pyvmdstream.VMDStream()
            vmdstream.draw_atomic(configuration)
            vmdstream.close()

            See examples/draw_test for complete and working usage examples.

    """
    def __init__(self,command_script="remote_ctl.tcl",port=5555):
        """

            Called when a new VMDStream object is created.
            Will create the file 'command_script' and open a connection to 
              VMD on port 'port'

        """
        self.s = vmdstart(command_script,port)

    def close(self):
        """ Close the connection to VMD """
        vmdstop(self.s)

    def send(self,sendstring):
        """ Send a command directly to VMD """
        self.s.send(sendstring)

    def draw_atomic( self, configuration, atomtypes=None, default_radius=0.5,
                     radii=None, radius_list=None, color_list=None,
                     color_value_list=None,sphere_resolution=30,
                     connecting_segments_types=None, bond_list=None, cylinder_radius_fraction=0.5, reset_view=True):
        """
        
            This is an example function for drawing an atomic configuration.
   
            Arguments:
                Configuration     - An (N,3) numpy array giving the 
                                      atomic coordinates
                Atomtypes         - A length N numpy array specifying atomtypes 
                                      in the range 0-natomtypes
                Default_radius    - default size radius
                Radii             - Specify the radius of each particle type
                                      in a numpy array of length (natomtypes)
                Radius_list       - Specify the radius of each particle in a
                                      numpy array of length N
                Color list        - Specify the color of each atom in a numpy
                                      array of length N integers in the range 
                                      VMDSTARTCOLOR to VMDNCOLORS,
                Color value       - Specify the color of each atom in a numpy
                                      array of length N floats in the range 0-1
                Sphere resolution - VMD sphere resolution, integer
                Connecting
                    Segments Types - Connect atomic coordinates of types in this
                                       list with cylinders if the same type
                Bond List          - List of bead pairs to connect with cylinders
                                     (0 indexed)
                Cylinder radius 
                    fraction      - Fraction of bead size for cylinder radius
                                     (default: 0.5)
                Reset view        - recenter view every time step (True/False)

        """
        natoms = len(configuration)

        #example initial setup
        self.s.send('axes location off\n')
        self.s.send('display projection orthographic\n')
        self.s.send('display resize 800 800\n')
        self.s.send("draw delete all\n")
        self.s.send('draw materials on\n')
        self.s.send('draw material "HardPlastic"\n')

        if color_value_list is not None:
            # readjust color value list
            color_list = np.array( np.floor( color_value_list*VMDNCOLORS ),dtype=int)

        if bond_list is not None:
            for link_idx in range(len(bond_list)):
                i,j = bond_list[link_idx]
                if color_list is not None:
                    self.s.send("draw color %i\n"%(color_list[i]+VMDSTARTCOLOR))
                elif atomtypes is not None:
                    self.s.send("draw color %i\n"%(atomtypes[i]))
                if radii is not None and atomtypes is not None:
                    this_radius = radii[atomtypes[i]]        
                elif radius_list is not None:
                    this_radius = radius_list[i]
                else:
                    this_radius = default_radius
                self.s.send("draw cylinder {%f %f %f} {%f %f %f} radius %f resolution %i filled yes\n"%( configuration[i,0],configuration[i,1],configuration[i,2],configuration[j,0],configuration[j,1],configuration[j,2],this_radius*cylinder_radius_fraction,sphere_resolution) )

        for i in range(len(configuration)):
            if color_list is not None:
                self.s.send("draw color %i\n"%(color_list[i]+VMDSTARTCOLOR))
            elif atomtypes is not None:
                self.s.send("draw color %i\n"%(atomtypes[i]))

            if radii is not None and atomtypes is not None:
                this_radius = radii[atomtypes[i]]        
            elif radius_list is not None:
                this_radius = radius_list[i]
            else:
                this_radius = default_radius
    
            self.s.send("draw sphere {%f %f %f} radius %f resolution %i\n"%( configuration[i,0],configuration[i,1],configuration[i,2],this_radius,sphere_resolution) )
            if i+1<len(configuration) and connecting_segments_types is not None:
                for connecting_segments_type in connecting_segments_types:
                    if atomtypes[i] == connecting_segments_type and atomtypes[i+1] == connecting_segments_type:
                        self.s.send("draw cylinder {%f %f %f} {%f %f %f} radius %f resolution %i filled yes\n"%( configuration[i,0],configuration[i,1],configuration[i,2],configuration[i+1,0],configuration[i+1,1],configuration[i+1,2],this_radius*cylinder_radius_fraction,sphere_resolution) )
 

        if reset_view is True:
            self.s.send('display resetview\n')
            self.s.send('scale by 1.3\n')

    def set_colorscale(self,colormap="jet",ncolors=VMDNCOLORS,startcolorid=VMDSTARTCOLOR):
        """

            This allows one to set the 'color scale' of VMD ( the colors from 
              %i to %i ) to colors from a matplotlib color map.

            Giving no colormap should reset the colors to VMD defaults.

        """
        # reset vmd color scale
        # an alternative would be to change the colorscale in vmd, e.g.
        #    self.s.send("color scale method GWR\n")
        if colormap == "":
            self.s.send("color scale method RGB\n")
        else:
            # this is an obj that takes a number in zero to one and returns
            #   r g b intensity values
            colors = cmap_discretize(ncolors,cmap=colormap)
            for i in range(ncolors):
                color_frac = float(i)/float(ncolors)
                color_vals = colors(color_frac)[:3]
                # set color startcolor+i to r g b value
                self.s.send("color change rgb %i %0.3f %0.3f %0.3f\n"%( startcolorid+i,color_vals[0],color_vals[1],color_vals[2]) )

    def render_tachyon(self,file_prefix="test_frame"):
        """ Render file via tachyon with current VMD defaults """
        #self.s.send('render Tachyon %(FILE_PREFIX)s "/usr/local/lib/vmd/tachyon_LINUXAMD64" -aasamples 12 %%s -format TARGA -o %%s.tga\n'%{'FILE_PREFIX':file_prefix+'.dat'})
        #self.s.send('render Tachyon %(FILE_PREFIX)s "tachyon" -aasamples 12 %%s -format TARGA -o %%s.tga\n'%{'FILE_PREFIX':file_prefix+'.dat'})
        self.s.send('render Tachyon %(FILE_PREFIX)s "/software/vmd-1.9.2-x86_64/lib/tachyon_LINUXAMD64" -aasamples 12 %%s -format TARGA -o %%s.tga\n'%{'FILE_PREFIX':file_prefix+'.dat'})

def ctl_script(port):
    """
        Return a vmd startup script
    
        This function returns a tcl startup script for vmd which opens TCP 
          connection on port PORT
    
        Arguments:
            port=PORT - this specifies TCP port VMD will open
    
    """

    return """\
#------------------------------------------------------------------
# $Id: remote_ctl.tcl,v 1.6 2003/02/12 21:33:11 oliver Exp $
# based on bounce.tcl and vmdcollab.tcl
# from http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/vmdcollab/
#
# start this in VMD and send commands to the listening port to have VMD
# execute them remotely
# (also see http://www.tcl.tk/scripting/netserver.html)
#
# Usage: vmd -e remote_ctl.tcl
# or vmd> source remote_ctl.tcl
#
# Security: we only allow connections from localhost (see acpt)
#
# Bugs:
# * once a wrong command was sent, the connection appears
# to 'block' and does not accept correct commands later
# * does not write result back to socket (one way connection...) so
# there is no way to inquire objects in vmd

namespace eval remote_ctl {
variable main
variable clients
variable default_vmd_port
set default_vmd_port 5555
# I am too dumb to set the default value for port from
# $default_vmd_port so I put 5555 in there literally
proc start { {port %(PORT)s} } {
variable main
set main [socket -server remote_ctl::acpt $port]
putlog "Listening on port $port"
}
proc acpt { sock addr port } {
variable clients
if {[string compare $addr "127.0.0.1"] != 0} {
putlog "Unauthorized connection attempt from $addr port $port"
close $sock
return
}
putlog "Accept $sock from $addr port $port"
set clients($sock) 1
fconfigure $sock -buffering line
fileevent $sock readable [list remote_ctl::recv $sock]
}
proc recv { sock } {
variable main
variable clients
if { [eof $sock] || [catch {gets $sock line}]} {
# end of file or abnormal connection drop:
# shut down this connection
close $sock
putlog "Closing $sock"
unset clients($sock)
} else {
if {[string compare $line "quit"] == 0} {
# prevent new connections
# existing connections stay open
# No -- Bug(?): 'quit' closes VMD...
putlog "Disallowing incoming connections by request of $sock"
close $main
}
# execute the received commands
# should check for runtime errors which otherwise leave the connection
# in an unusable state
# eval $line
set rc [catch $line result]
if { $rc } {
#puts $sock "Error executing comand '$line': n$result"
puts "Error executing comand '$line': n$result"
} else {
#puts $sock $result
#puts $result
}
}
}
###### would like the last line from stdout in line ###########
# (or any working solution....)
proc send { sock line} {
variable clients
# send reply to connecting client
putlog "send '$line' to $sock"
puts $sock $line
}

proc putlog { text } {
puts $text
return
}
}

remote_ctl::putlog "Starting remote_ctl server in vmd: connect with something like"
remote_ctl::putlog "telnet localhost %(PORT)s"
remote_ctl::start
"""%{'PORT':port}

def vmdstart(command_script="remote_ctl.tcl",port=5555):
    """
        Start a VMD session

        Open VMD and generate/run the file COMMAND_SCRIPT, which will open a 
          TCP connection on port PORT

        Arguments:
            command_script=COMMAND_SCRIPT - this specifies a file name where 
              VMD startup commands will be written out. Note: this file will be
              overwritten. The default is 'remote_ctl.tcl'
            port=PORT -- this specifies TCP port VMD will open

    """
    open(command_script,'w').write( ctl_script(port) )
    cmd = subprocess.Popen("vmd -e "+command_script,shell=True)
    time.sleep(5)
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    for i in range(5):
        try:
	    s.connect(("127.0.0.1",port))
        except Exception, e:
            if e.errno == 106 or e.errno == 56:
                #already connected
                break
            else:
                print e
        time.sleep(2)
    return s

def vmdstop(s):
    """
        Close a VMD session

        Sends the exit command to VMD and closes the python TCP socket which 
          was used for communication with VMD

        Arguments:
            s=SOCKET - an open python socket object

    """
    s.send("exit\n")
    s.close()

