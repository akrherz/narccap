"""
 NCOUT files did not contain a bigheader, nor sigmah information necessary
 for archiver to work.  This attempts to fix it!
"""

import struct
import os

sigmah = open("sigmah", "rb").read()
os.chdir("Run.NCEP")

for i in range(1, 1226):

    bigheader = open("MMOUT_DOMAIN1_%04i" % (i,), "rb").read()[
        :117620
    ]  # 117600 + 5x4

    a = open("NCOUT_DOMAIN1_%04i" % (i,), "rb")
    out = open("NCOUT_DOMAIN1_%04i.new" % (i,), "wb")
    out.write(bigheader)
    newtime = True
    while True:
        if newtime:
            print "WRITE SIGMAH!"
            out.write(sigmah)
            newtime = False
        data = a.read(4)
        if len(data) == 0:
            print "Done!"
            break
        out.write(data)
        flagbytes = struct.unpack(">I", data)[0]  # always 4?

        data = a.read(4)
        flag = struct.unpack(">I", data)[0]
        out.write(data)

        data = a.read(4)
        flagbytes2 = struct.unpack(">I", data)[0]
        out.write(data)

        print "Flag: %s flagbytes: %s flagbytes2: %s" % (
            flag,
            flagbytes,
            flagbytes2,
        )
        if flag == 2:
            newtime = True
            continue

        data = a.read(4)
        readbytes = struct.unpack(">I", data)[0]
        out.write(data)

        header = a.read(readbytes)
        out.write(header)
        if readbytes == 152:
            h = struct.unpack(">i4i4if4s4s24s9s25s46s", header)
            print h[12], h[13]

        data = a.read(4)
        readbytes2 = struct.unpack(">I", data)[0]
        out.write(data)
        print "readbytes: %s readbytes2: %s" % (readbytes, readbytes2)
        if flag == 1:
            data = a.read(4)
            out.write(data)
            datasz = struct.unpack(">I", data)[0]

            data = a.read(datasz)
            out.write(data)

            data = a.read(4)
            datasz2 = struct.unpack(">I", data)[0]
            out.write(data)
            print "Read datasize: %s datasize2: %s" % (datasz, datasz2)

    out.close()
