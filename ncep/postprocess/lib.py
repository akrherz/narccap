
import mx.DateTime
import mm5_class

def extract_times(mm5file):
    """
    Function to return an array of mx.DateTime instances for the timestamps
    in this MM5 file
    """
    mm5 = mm5_class.mm5(mm5file)
    # Sample for first timestamp with variable u
    # Requires a modification to mm5_class.py to read header
    tstr = mm5.get_field_list()['u']['time'][:13]
    ts0 = mx.DateTime.strptime(tstr, "%Y-%m-%d_%H")
    interval = mm5.timeincrement
    ts1 = ts0 + mx.DateTime.RelativeDateTime(
                seconds=(interval * (mm5.tsteps-1) ))
    taxis = []
    now = ts0
    while (now <= ts1):
        taxis.append( now )
        now += mx.DateTime.RelativeDateTime(seconds=interval)
    del(mm5)
    return taxis
