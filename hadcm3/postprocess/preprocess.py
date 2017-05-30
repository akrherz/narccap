"""
 We need to convert the output files and save as much space in the 
 process...
"""
import sys
import os
import glob
import mm5_class
import mx.DateTime
import numpy
import lib

def process_MMOUT(datadir):
    os.chdir(datadir)
    # Figure out a list of MMOUT files
    files = glob.glob("MMOUT_DOMAIN1_????")
    files.sort()
    # Move us to interpb
    for file in files:
        # We don't wish to convert this file, it is just along for the ride
        if file == "MMOUT_DOMAIN1_000":
            continue
        os.chdir("/tera10/akrherz/narccap/INTERPB")
        # Figure out time axis
        taxis = lib.extract_times(datadir+file)

        # Setup variable substitution values
        vars = {}
        vars['mm5file'] = datadir+file
        vars['syear'] = taxis[0].year
        vars['smonth'] = taxis[0].month
        vars['sday'] = taxis[0].day
        vars['shour'] = taxis[0].hour
        if taxis[0].day == 21:
            taxis[-1] = taxis[0] + mx.DateTime.RelativeDateTime(months=1,day=1)
        vars['eyear'] = taxis[-1].year
        vars['emonth'] = taxis[-1].month
        vars['eday'] = taxis[-1].day
        vars['ehour'] = taxis[-1].hour
        

        # Edit the namelist.input for interb
        data = open('namelist.tpl', 'r').read()
        out = open('namelist.input', 'w')
        out.write( data % vars )
        out.close()
        
        # Run interb for each file
        print "Running INTERPB for %s [%s - %s]" % (file,
              taxis[0].strftime("%Y-%m-%d %H"), 
              taxis[-1].strftime("%Y-%m-%d %H"))
        os.system("./interpb >& interpb.log")

        # Move output file to right location
        os.rename("MMOUTP_DOMAIN1", datadir + file.replace("UT", "UTP"))
        # Cleanup
        os.system("rm -rf FILE_*")

        # Gzip
        os.chdir( datadir )
        os.system("gzip %s" % (file,))

        # Convert to NetCDF
        file = file.replace("UT", "UTP")
        mm5 = mm5_class.mm5(file)
        cmd = "archiver %s 0 %s" % (file, mm5.tsteps)
        print "Converting %s to NetCDF %s tsteps" % (file, mm5.tsteps)
        si,so = os.popen4( cmd )
        a = so.read() # Necessary to keep things blocking?

        # Remove MMOUTP files
        if file[:6] == "MMOUTP":
            os.remove(file)

def process_NCOUT(datadir):

    # Change directory
    os.chdir( datadir )
    # Look for any MMOUT and NCOUT files
    files = glob.glob('NCOUT_DOMAIN1_[0-9][0-9][0-9][0-9]')
    files.sort()
    # Loop over the files
    for file in files:
        # Skip the initial output files, we don't care about these
        if file == "NCOUT_DOMAIN1_000":
            continue
        mm5 = mm5_class.mm5(file)
        cmd = "archiver %s 0 %s" % (file, mm5.tsteps)
        print "Converting %s to NetCDF %s tsteps" % (file, mm5.tsteps)
        si,so = os.popen4( cmd )
        a = so.read() # Necessary to keep things blocking?
        
        # Gzip
        os.system("gzip %s" % (file,))

if __name__ == '__main__':
    datadir = sys.argv[1]
    #process_NCOUT(datadir)
    process_MMOUT(datadir)
