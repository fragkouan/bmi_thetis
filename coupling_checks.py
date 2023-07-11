"""
Check various stuff before coupling the models to ensure that no stupid mistake
will take place
"""
from colours import clrtxt


def check_coupling_status(thetis, swan_int):
    """
    """

    if (thetis=='2-way' and swan_int==2):
        stat = 'Fully coupled'
    elif (thetis=='Thetis-to-SWAN' and swan_int==1):
        stat = 'Thetis-to-SWAN'
    elif (thetis=='SWAN-to-Thetis' and swan_int==0):
        stat = 'SWAN-to-Thetis'
    elif (thetis=='No coupling' and swan_int==-1):
        stat = 'No coupling'
    else:
        swan = convert_int_to_swan_coupling_stat(swan_int)
        message = clrtxt.bold + clrtxt.red + "Terminating program: " +\
            clrtxt.end + "Thetis and SWAN do not have the same coupling " +\
            "status. Thetis has '" + clrtxt.bold + clrtxt.green + thetis + clrtxt.end +\
            "' as coupling status and SWAN '" + clrtxt.bold + clrtxt.green + swan +\
            clrtxt.end +"'. The two models need to have the same coupled status"
        raise SystemExit(message)

    return stat

def convert_int_to_swan_coupling_stat(swan_int):
    """
    """
    if swan_int==2:
        swan = 'Fully coupled'
    elif swan_int==1:
        swan = 'Thetis-to-SWAN'
    elif swan_int==0:
        swan = 'SWAN-to-Thetis'
    else:
        swan = 'No coupling'

    return swan

def check_coupling_dt(thetis, swan):
    if abs(thetis-swan)>0.0001:
        message = clrtxt.bold + clrtxt.red + "Terminating program: " +\
            clrtxt.end + "Thetis coupling timestep does not have the same " +\
            "value as the coupling timestep of SWAN. The former is " +\
            str(thetis) + " sec, while the latter is " + str(swan) + " sec."
        raise SystemExit(message)

def check_simulation_duration(thetis, swan):
    if abs(thetis-swan)>0.0001:
        message = clrtxt.bold + clrtxt.red + "Terminating program: " +\
            clrtxt.end + "Thetis simulation duration does not have the same " +\
            "value as the simulation duration of SWAN. The former is " +\
            str(thetis) + " sec, while the latter is " + str(swan) + " sec."
        raise SystemExit(message)
