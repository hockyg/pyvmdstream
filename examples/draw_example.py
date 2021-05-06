import os
import sys
import numpy as np
import time

import pyvmdstream

def load_xyzfile(xyzfile):
    """
        
        load an atomic configuration in a very simple xyz format

        Note: this function is for exapmle purposes and does not have any error
          checking and so will crash if syntax is not specified correctly

        Arguments:
            xyzfile - a plain text file in the format specified below. This
              may contain more than one configuration concatenated together

        Return values:
            trajectory - a list of numpy arrays, each of shape (N,3), and each
              containing a single configuration
            atomtype_list - a list of numpy arrays, each of size N, specifying
              an atomtype for each atom in the file
            boxsize_list - a list of numpy arrays, each specifying a set
              (Lx,Ly,Lz) prepresenting simulation box size dimensions 

        file format:
            N ANYTHING
            atomid1 atomtype1 x1 y1 z1
            atomid2 atomtype2 x2 y2 z2
            ...
            atomidN atomtypeN xN yN zN
            Lx Ly Lz

    """
    trajectory = []
    atomtype_list = []
    boxsize_list = []

    fh = open(xyzfile,'r')
    line = fh.readline()

    atomtype_label_list = []
    N = int( line.split()[0] )
    while line:
        configuration = np.zeros((N,3))
        atomtypes = np.zeros(N,dtype=int)
        for i in range(N):
            line = fh.readline()
            atomid,atomtype,x,y,z = line.split()[:5]

            # this converts atomtype to a numeric atomtype, starting with zero,
            #   regardless of format

            if atomtype in atomtype_label_list:
                atomtype_numeric = atomtype_label_list.index(atomtype) 
            else:
                atomtype_numeric = len(atomtype_label_list)
                atomtype_label_list.append(atomtype)

            atomtypes[i] = atomtype_numeric

            pos = np.array((x,y,z),dtype=float)
            # set the i'th row of configuration to be the xyz position of atomi
            configuration[i,:] = pos
        line = fh.readline()
        boxsize = np.array( line.split(),dtype=float)

        trajectory.append(configuration)
        atomtype_list.append(atomtypes)
        boxsize_list.append(boxsize)

        line = fh.readline()

    return trajectory, atomtype_list, boxsize_list

def divider():
    print("-"*80)

def main():
    """

        This function has a series of examples of how to use the pyvmdstream library

    """
    print("Warning, this program should write out files named test_frame_%%i.* in the current directory")
    kbd_input = input("Press Enter to continue, To cancel, press Ctrl-C or type anything else before pressing enter\n")
    if not kbd_input == "":
        print("Execution canceled")
        sys.exit()

    ##################
    # simple example #
    ##################
    divider()
    print("First, a very simple example of how to use the interface")
    divider()

    # create a new connection to vmd 
    s = pyvmdstream.vmdstart(port=5556)
    # send a command to vmd
    s.send("draw sphere {0 0 0} radius 1 resolution 30\n".encode())
    # wait so you can see the result
    time.sleep(3)
    # stop the connection to vmd
    pyvmdstream.vmdstop(s)

    ######################
    # VMDStream examples #
    ######################

    print("Now, examples using a sample VMDStream class")
    divider()

    # load configuration
    examplefile = os.path.join( os.path.dirname(sys.argv[0]), "example_2d_configurations.xyz")
    trajectory, atomtype_list, boxsize_list = load_xyzfile(examplefile)

    # create a new VMDStream object
    vmdstream = pyvmdstream.VMDStream()

    print("Drawing plain configuration")
    divider()

    # draw a very plain rendering of configuration 0
    vmdstream.draw_atomic( trajectory[0] )
    time.sleep(2)

    # set a pair of radii for the two atomtypes in the sample input
    radii = np.array((0.45,0.35))

    print("Drawing configuration colored by atomtype, and with correct sizes")
    divider()
    # render the scene specifying atomtypes for each atom, and radii for the two atomtypes 
    vmdstream.draw_atomic( trajectory[0], atomtype_list[0],radii=radii )
    time.sleep(2)

    print("Coloring randomly by color id")
    divider()
    # Pick random colors in the range 0-1023 for the N atoms
    random_color_list = np.random.randint(1024,size=len(trajectory[0]))
    vmdstream.draw_atomic( trajectory[0], atomtype_list[0],color_list=random_color_list,radii=radii )
    time.sleep(2)

    print("Coloring randomly by color values in 0-1")
    divider()
    # Pick random values in the range 0-1 ( which will be converted to color id's in the range 0-1023 ) for the N atoms
    random_values_list = np.random.random(size=len(trajectory[0]))
    vmdstream.draw_atomic( trajectory[0], atomtype_list[0],color_value_list=random_values_list,radii=radii )
    time.sleep(2)

    print("Coloring by atom type but random radius sizes")
    divider()
    # Choose random sizes for all N particles, color by atomtype "
    random_radii = 0.2+0.4*np.random.random(size=len(trajectory[0]))
    vmdstream.draw_atomic( trajectory[0], atomtype_list[0],radius_list=random_radii )
    time.sleep(2)

    print("Now make a small movie with frame 2, changing colors and radii with postition and time")
    divider()

    # set the color scheme to one from matplotlib, "jet" by default
    vmdstream.set_colorscale()
    vmdstream.set_colorscale(colormap="")
    nframes = 16

    # we will now use the second configuration from the sample file
    configuration = trajectory[1]
    atomtypes = atomtype_list[1]

    # generate a list of radii of length N, specifying the sizes of each atom individually
    radius_list = radii[atomtypes]
   
    # calculate the "x" distance from the side of the box
    dx = ( configuration[:,0] - boxsize_list[1][0] )
    # scale that distance to the range 0-1
    dx_scaled = (dx+dx.min())/(dx.max()-dx.min())

    for i in range(nframes):
        # sinusoidally color the particles based on distance from the side
        color_value_list = 0.5+0.5*np.sin(dx_scaled*np.pi + 2*i*np.pi/nframes-np.pi/2)
        # scale down the radius of eaach atom, then sinusoidally modify it based on the distance from the bottom
        radius_list_scaled = radius_list*0.8 + 0.2*np.sin(dx_scaled*np.pi + 2*i*np.pi/nframes-np.pi/2)
        vmdstream.draw_atomic( configuration, atomtypes,color_value_list=color_value_list,radius_list=radius_list_scaled )
        # save the configuration and render it using the "tachyon" plugin in vmd
        vmdstream.render_tachyon(file_prefix="test_frame_%02i"%i)

    # close the vmd connection
    vmdstream.close()

if __name__ == "__main__":
    main()
