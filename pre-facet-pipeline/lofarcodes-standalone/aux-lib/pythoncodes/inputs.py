import getpass
from pyslalib import slalib, sladoc

def find_username():
    username = getpass.getuser()
    return username

#-----------------------------------------------------------

def yes_no(question):
    """ 
    Question requiring y/n input
    """

    input = str(raw_input(question))
    while input != 'y' and input != 'n':
        print 'Enter y or n '
        input = str(raw_input(question))
    return input
#-----------------------------------------------------------

def split_remove_blanks(line):
    """
    Remove the blanks from a line and return a string containing all elements of
 the line
    """
    line = line[:-1]
    line = line.split(' ')
    while '' in line:
        line.remove('')
    return line


#-----------------------------------------------------------

def enter_coords(input):

    """
    Allows user to enter coordinates in degrees or J2000 and returns
    radians

    """


    if input == 'J2000':

        RA = str(raw_input("Enter RA (hh mm ss)   "))
        RA = RA.split(' ')
        (rarad,status)  = slalib.sla_dtf2r(RA[0],RA[1],RA[2]) # hms2rad


        DEC= str(raw_input("Enter Dec (deg mm ss)    "))
        DEC = DEC.split()
        (decrad,status) = slalib.sla_daf2r(DEC[0],DEC[1],DEC[2])   # dms2rad
        

    elif input == 'DEG':
        RA = float(raw_input("Enter RA (degrees)   "))
        rarad = RA*deg2rad
        
        DEC= float(raw_input("Enter Dec (degrees) "))
        decrad = DEC*deg2rad
        


    return rarad,decrad


#-----------------------------------------------------------

