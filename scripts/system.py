from __future__ import division, print_function; __metaclass__ = type
import os, sys, math, itertools
import numpy
import west
from west import WESTSystem
from westpa.binning import RectilinearBinMapper, RecursiveBinMapper

import logging
log = logging.getLogger(__name__)
log.debug('loading module %r' % __name__)

class System(WESTSystem):
    def initialize(self):
        self.pcoord_ndim            = 2
        self.pcoord_len             = 2
        self.pcoord_dtype           = numpy.float32
        minfrag                     = 90
        outer_mapper                = RectilinearBinMapper([[0,minfrag-0.5,float('inf')],[0,float('inf')]])
        binbounds1                  = [1.0*i-0.5 for i in xrange(minfrag,121)] + [float('inf')]
        binbounds2                  = [15.0*i*(i-1)/2 for i in xrange(1,20)] + [float('inf')]
        self.bin_mapper             = RecursiveBinMapper(outer_mapper)
        self.bin_mapper.add_mapper(RectilinearBinMapper([binbounds1,binbounds2]),[minfrag,1])
        self.bin_target_counts      = numpy.empty((self.bin_mapper.nbins,), numpy.int)
        self.bin_target_counts[...] = 4
        #this is intended to change the bins where distance < 200 A to 16 simulations per bin
        for q1 in xrange(minfrag,121):
             for q2 in [0.5+15.0*i*(i-1)/2 for i in xrange(1,8)]:
                   print(q1,q2)
                   bin=self.bin_mapper.assign(numpy.array([(q1,q2)])) #tuple?
                   self.bin_target_counts[bin] = 16

def coord_loader(fieldname, coord_file, segment, single_point=False):
    coord_raw = numpy.loadtxt(coord_file, dtype=numpy.float32) 
    npts = len(coord_raw)
    assert coord_raw.shape[1] % 3 == 0
    ngrps = coord_raw.shape[1] // 3

    coords = numpy.empty((ngrps, npts, 3), numpy.float32)
    for igroup in xrange(ngrps):
        for idim in xrange(3):
            coords[igroup,:,idim] = coord_raw[:,igroup*3+idim]

    segment.data[fieldname] = coords
