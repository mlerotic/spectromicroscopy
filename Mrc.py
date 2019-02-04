"""Provide methods for reading and writing files in the MRC
format.

Requires NumPy.  This module has been imported successfully
when used with the following combinations of Python and
NumPy:
    Python 2.5.6 ; NumPy 1.3.0.dev6083
    Python 2.6.6 ; NumPy 1.4.1
    Python 2.7.9 ; NumPy 1.8.2
    Python 2.7.10 ; NumPy 1.8.0rc1
    Python 3.4.2 ; NumPy 1.8.2
.  Other combinations for Python versions greater than or
equal to 2.5 and NumPy versions greater than or equal to
1.3.0 likely work.

Implements the MRC file format as described at
http://msg.ucsf.edu/IVE/IVE4_HTML/IM_ref2.html .

The Mrc class is likely easiest to use if you want read-only
access to an Mrc file or want read-write access and the
modifications that you'll make do not affect the size or
format of the extended header or image data.  The Mrc class
does use memory mapping of the file.  Depending on the
system, that may make it unusable for large files.

An example of how to use the Mrc class is this:

import numpy
import Mrc
a = Mrc.bindFile('somefile.mrc')
# a is a NumPy array with the image data memory mapped from
# somefile.mrc.  You can use it directly with any function
# that will take a NumPy array.
hist = numpy.histogram(a, bins=200)
# Print out key information from the header.
a.Mrc.info()
# Use a.Mrc.hdr to access the MRC header fields.
wavelength0_nm = a.Mrc.hdr.wave[0]

If you only want a copy of all the highest resolution data
set from a MRC file as a NumPy array, you can use load():

import numpy
import Mrc
a = Mrc.load('somefile.mrc')

The Mrc2 class is what you would use if you wanted to create
a MRC file from scratch.  An example of that is:

import numpy
import Mrc
a = numpy.reshape(numpy.asarray(
                  numpy.linspace(0, 5999, 6000), dtype='i2'),
                  (10,20,30))
m = Mrc.Mrc2('newfile.mrc', mode='w')
# If you want the header to use the MRC 2014 format rather than
# the Priism format, insert
# m.hdr = Mrc.hdrChangeToMrc2014Format(m.hdr)
# Set the header size and pixel type fields.  You could also use
# m.initHdrForArr(a) instead.
m.setHdrForShapeType(a.shape, a.dtype)
# Set other fields in the header.
m.hdr.setSpacing(0.1, 0.1, 0.3)
m.hdr.wave[0] = 540
m.hdr.mmm1 = (0.0, 5999.0, 2999.5)
m.hdr.setTitle('Written by Mrc2 of Mrc.py')
m.writeHeader()
m.writeStack(a)
m.close()

A short cut for writing a NumPy array as an Mrc file is to use
save().  If the only fields that you want to set in the header
are the basic size and pixel type information, it is very simple
to use:

import numpy
import Mrc
a = numpy.reshape(numpy.asarray(
                  numpy.linspace(0, 5999, 6000), dtype='i2'),
                  (10,20,30))
Mrc.save(a, 'newfile.mrc', ifExists='overwrite')

Known limitations:
1) Does not support the 4-bit unsigned integer sample format (mode 101).
2) Does not provide any mechanism for accessing the lower resolution
versions of the image data allowed by Priism's version of the MRC format.
3) Uses a NumPy structured array type to represent mode 3 (complex values
represented by two signed 16-bit integers) image data.  Nothing is provided
to make that data act more like standard complex arrays in NumPy.

Release notes:
Version 0.1.1 corrected Mrc2.makeSymmetryInfo() since it didn't use self
when accessing the header information.  Also changed the logic in
getExtHeaderFormat() to match Priism's when working with a file without
the extended header type field but with the map field.  In that case assume
a non-crystallographic extended header if the space group is 0, 1, or 401.
Version 0.1.0 as included with Priism made the following changes which
could affect compatibility with client code:
1) Added FROM_PRIISM to the module to distinguish this version from other
versions of Mrc.py derived from Priithon.
2) Added HdrBase, HdrPriism, Hdr2014, Hdr2014Priism, ManagedTitleArray, and
ReorderedArray classes.  The combination of HdrBase and HdrPriism replaced
a class defined within implement_hdr().  They do not have the type data
attribute (2 byte field immediately after mmm1) that class had and change
the interpretation of nspg (absorbing what had been in the type data
attribute) and blank (split by adding ntst).  The title attribute changed
from being backed by a a 10a80 element in the structured array to a
(10,80)i1 element.  Since interactions with the title attribute are mediated
through a ManagedTitleArray, assigning strings to the titles should work
as before.
3) Mode 0 image data is interpreted as signed 8-bit integers for
compatibility with the MRC 2014 standard.
4) Mode 3 image data (complex values stored as two signed 16-bit integers)
was treated as one 4-byte floating-point value.  Now it's represented by a
NumPy structured array with two 16-bit signed integer fields.  The first
field is called 'real' and the second is called 'imag'.
5) The extInts and extFloats data attributes of the Mrc and Mrc2 classes are
now always set.  They may be None, depending on the size and format of the
extended header.
6) Removed insertExtHdr() from the Mrc class.
7) Removed extHdrSize, extHdrnint, and extHdrnfloat keywords from __init__()
for the Mrc class.
8) Changed the initial values in the header when an instance of Mrc2 is
created for a new file.
9) Changed the return value of MrcMode2dtype() from the generic Python type
equivalent to the sample representation to a NumPy dtype.
10) As a result of the changes for (2), changed the type returned by
implement_hdr() and makeHdrArray().
11) Changed initHdrArrayFrom() to have a return value.
"""

__author__ = 'Sebastian Haase <haase@msg.ucsf.edu>'
__license__ = 'BSD license - see PRIITHON_LICENSE.txt'
__version__ = '0.1.1'

import numpy as N
import sys
import os
# Python 3 renamed __builtin__ to builtins.  Work around that without using
# from __future__ so that this will work with older versions of Python 2.
try:
    import __builtin__ as builtins
except ImportError:
    import builtins
import string
import tempfile
import weakref

# Mark this as a Mrc.py originating from Priism in case a client wants to
# try to distinguish between this version and versions from elsewhere.
FROM_PRIISM = True


def bindFile(fn, writable=0):
    """Return a NumPy array memory mapped from an existing MRC file.

    The returned NumPy array will have an attribute, Mrc, that is
    an instance of the Mrc class.  You can use that to access the
    header or extended header of the file.  For instance, if x
    was returned by bindFile(), x.Mrc.hdr.Num is the number of x
    samples, number of y samples, and number of sections from the
    file's header.

    Positional parameters:
    fn -- Is the name of the MRC file to bind.

    Keyword parameters:
    writable -- If True, the returned array will allow the elements of
    the array to be modified.  The Mrc instance packaged with the
    array will also allow modification of the header and extended
    header entries.  Changes made to the array elements or the header
    will affect the file to which the array is bound.
    """

    mode = 'r'
    if writable:
        mode = 'r+'
    a = Mrc(fn, mode)

    return a.data_withMrc(fn)


class Mrc:
    """Provide memory mapped access to existing MRC files.

    If x is a Mrc instance, x.data has the highest resolution image
    data from the file, x.hdr has the header fields, and x.e has the
    extended header as an array of unsigned 8-bit integers.  If the
    extended header size is greater than zero, the extended header
    format is Priism's format, and the number of integers or
    floating-point values per section in the extended header is
    greater than zero, x.extInts has the integer values from
    the extended header, x.extFloats has the floating-point values
    from the extended header, and x.extSym is None.  If the
    extended header has symmetry information, x.extInts is None,
    x.extFloats is None, and x.extSym is an array of 80 character
    records for the symmetry information.  Any other cases will
    have x.extInts equal to None, x.extFloats equal to None, and
    x.extSym equal to None.

    If used to modify a file, the Mrc class does not provide
    public methods to change the format or size of the image data
    or change the size of the extended header data and then remap
    the image data and extended header data into memory.  Because of
    that, it is best to use the Mrc class for read-only access or
    for read-write access where the modifications do not change
    the layout or format of the MRC file.  The Mrc2 class can
    handle more general modifications to an existing MRC file.

    Depending on the system and the handling of memory mapping
    in Python, it may not be possible to memory map files which
    are too large (around 1 gigabyte was a frequent barrier with
    Python 2.5 and earlier).  The Mrc2 class which does not use
    memory mapping could be useful for those files.

    If the size of the extended header is not a multiple of the
    size of the data type used to represent one image sample,
    the memory mapped image data will be misaligned.  Use the
    Mrc2 class for files like that.

    Version 0.1.0 of Mrc.py as included with Priism removed the
    insertExtHdr() function from this class.  It also changed the
    conditions for when the extInts and extFloats attributes are
    set.
    """

    def __init__(self, path, mode='r'):
        """Initialize the Mrc object.

        Maps the entire file into memory.

        Verion 0.1.0 of Mrc.py as included with Priism removed
        the extHdrSize, extHdrnint, and extHdrnfloat keywords.

        Positional parameters:
        path -- Is the file to map into memory.

        Keyword parameters:
        mode -- If mode is 'r', requests read-only access
        to the file.  If mode is 'r+', requests read and
        write access to the file.
        """

        self.path = os.path.abspath(path)
        self.filename = os.path.basename(path)

        self.m = N.memmap(path, mode=mode)
        self.h = self.m[:1024]

        self.hdr = makeHdrArray(self.h)

        if hdrIsByteSwapped(self.hdr):
            self.hdr._array.dtype = self.hdr._array.dtype.newbyteorder()
            self.isByteSwapped = True
        else:
            self.isByteSwapped = False

        self.data_offset = 1024 + max(0, self.hdr.next)
        self.d = self.m[self.data_offset:]

        self.e = self.m[1024:self.data_offset]

        self.doDataMap()

        self.numInts = max(0, self.hdr.NumIntegers)
        self.numFloats = max(0, self.hdr.NumFloats)

        if self.hdr.next > 0:
            fmt = getExtHeaderFormat(self.hdr)
        else:
            fmt = -1
        if fmt == 0:
            self.doSymMap()
        elif fmt == 1 and (self.numInts > 0 or self.numFloats > 0):
            self.doExtHdrMap()
        else:
            self.extHdrArray = None
            self.extInts = None
            self.extFloats = None
            self.extSym = None


    def doExtHdrMap(self, nz=0):
        """Map a NumPy structured array to the Priism-style extended header.

        Creates self.extHdrArray, the structured array to represent the
        extended header.  Also generates self.extInts, a view of
        self.extHdrArray with the integer entries, and self.extFloats,
        a view of self.extHdrArray with the floating-point entries.
        Sets self.extSym to None.

        Keyword parameters:
        nz -- Is the number of sections of data to include in
        self.extHdrArray.  If nz is zero, the number of sections
        included will be the maximum of zero and the number of
        sections from the header (self.hdr.Num[2]).  If nz is
        less than zero, the number of sections will be the maximum
        possible given the size of the extended header and the
        number of integer and floating-point values per section.
        """

        if nz == 0:
            nz = max(0, self.hdr.Num[-1])

        maxnz = len(self.e) // ((self.numInts + self.numFloats) * 4)
        if nz < 0 or nz > maxnz:
            nz=maxnz

        byteorder = '='
        type_descr = [('int', '%s%di4' % (byteorder, self.numInts)),
                      ('float', '%s%df4' % (byteorder, self.numFloats))]

        self.extHdrArray = N.recarray(shape=nz, dtype=type_descr, buf=self.e)
        if self.isByteSwapped:
            self.extHdrArray = self.extHdrArray.newbyteorder()
        self.extInts = self.extHdrArray.field('int')
        self.extFloats = self.extHdrArray.field('float')
        self.extSym = None


    def doSymMap(self):
        """Map a NumPy structured array to the symmetry information.

        Creates self.extHdrArray, a structured array to represent the
        extended header.  Also generates self.extSym, an array of 80
        character strings mapped to as much of the extended header
        as possible.  Sets self.extInts and self.extFloats to None.
        """

        nrec = self.hdr.next // 80
        nrem = self.hdr.next - 80 * nrec
        type_descr = [('records', '(%d,80)i1' % nrec),
                      ('extra', '%di1' % nrem)]
        self.extHdrArray = N.recarray(shape=1, dtype=type_descr,
                                      buf=self.e)
        self.extSym = ManagedTitleArray(self.extHdrArray.field('records')[0])
        self.extInts = None
        self.extFloats = None


    def doDataMap(self):
        """Map a NumPy array to the highest resolution data set in the file.

        Creates self.data as the mapped array.
        """

        dtype = MrcMode2dtype(self.hdr.PixelType)
        shape = shapeFromHdr(self.hdr)

        self.data = self.d.view()
        self.data.dtype = dtype
        n0 = self.data.shape[0]
        n1 = N.prod(shape)
        if n0 < n1:
            # The file is smaller then the space needed for the lowest
            # resolution data set.  Truncate the slowest varying dimension.
            print('** WARNING **: file truncated - shape from header: %s '
                  'expected to get %s but got %s' %
                  (str(shape), str(N.prod(shape)), str(n0)))
            n1 = n1 // shape[0]
            s0 =  n0 // n1
            shape = (s0,) + shape[1:]
            self.data = self.data[:(s0*n1)]
        elif n0 > n1:
            # The file is larger than the space needed for the highest
            # resolution data set.  Ignore the excess.
            self.data = self.data[:n1]
        self.data.shape = shape

        if self.isByteSwapped:
            self.data = self.data.newbyteorder()


    def setTitle(self, s, i=-1, push=False, truncate=False):
        """Set a title in the MRC header.

        Provided for compatibility with previous versions of Mrc.py.
        In the version, you can use self.hdr.setTitle().

        Positional parameters:
        s -- Is the character string for the title.  If s is longer
        than 80 characters and truncate is False, a ValueError
        exception will be raised.  Since no byte swapping is done
        for the titles in the header, s should be encoded in ASCII
        or another format that does not use multibyte characters.

        Keyword parameters:
        i -- Is the index of the title to set.  If i is less than
        zero, the last title not in use will be set.  If i is less
        than zero and all the titles are in use and push is False
        or i is greater than 9, a ValueError exception will be
        raised.

        push -- If True, i is less than zero, and all titles are
        in use, titles will be pushed down before assigning s
        to the last title.  That will discard the first title and
        title[k] (for k greater than and equal to 0 and less than
        9) will be title[k+1] from before the change.  The push
        keyword was added in version 0.1.0 of Mrc.py as included
        with Priism

        truncate -- If True, only use the first 80 characters from
        s.  The truncate keyword was added in version 0.1.0 of
        Mrc.py as included with Priism.
        """

        self.hdr.setTitle(s, i=i, push=push, truncate=truncate)


    def axisOrderStr(self, onlyLetters=True):
        """Return a string indicating the ordering of dimensions.

        x, y, z, w, and t will appear at most once in the string, and
        at least three of them will be present.  The letters that do
        appear will be present in order from slowest varying to
        fastest varying.  The values for the axis field in the header
        do not affect the result.

        Keyword parameters:
        onlyLetters -- If True, only the letters for the dimensions
        will appear in the string.  If False, the first character
        of the string will be '[', the last character of the string
        will be ']', and a letter for a dimension will be preceded
        by a comma if it is not the first, slowest-varying, dimension.
        """

        return axisOrderStr(self.hdr, onlyLetters)


    def looksOK(self, verbose=1):
        """Perform basic tests on file.

        Currently tests the file size against what is expected from
        the header information.

        Keyword parameters:
        verbose -- If greater than or equal to one, some diagnostic
        information will be printed.  Larger values generate more
        diagnostic information.  Currently, values larger than three
        have the same effect as a value of three does.

        Return value:
        Returns True if all the tests pass.  Returns False if one or
        more fail.
        """

        shape = self.data.shape
        b = self.data.dtype.itemsize
        eb = N.prod(shape) * b
        ab = len(self.d)
        secb = N.prod(shape[-2:]) * b
        if self.hdr.sub > 1:
            nres = self.hdr.sub
            if self.hdr.zfac > 1:
                zfac = self.hdr.zfac
                if self.hdr.ImgSequence != 1:
                    # Treat an invalid image sequence value as if it was
                    # zero.  An image sequence value or zero or two has
                    # z as the fastest varying dimension after y and x.
                    if len(shape) >= 3:
                        resChangesNz = True
                        zind = -3
                    else:
                        resChangesNz = False
                else:
                    # WZT order has z changing next fastest after
                    # wavelength, y, and x.
                    if len(shape) >= 4:
                        resChangesNz = True
                        zind = -4
                    else:
                        resChangesNz = False
            else:
                resChangesNz = False
        else:
            nres = 1
            resChangesNz = False
        elb = 0
        resShape = list(shape)
        for i in range(1, nres):
            if resChangesNz:
                resShape[zind] = (resShape[zind] + zfac - 1) // zfac
            resShape[-1] = resShape[-1] // 2
            resShape[-2] = resShape[-2] // 2
            elb += N.prod(resShape)
        elb *= b

        # Computing the number of sections of data present only makes
        # sense if there is no additional resolutions present.
        if ab < eb or elb == 0:
            displaySectionInfo = True
            anSecs = ab / float(secb)
            enSecs = eb / float(secb)
        else:
            displaySectionInfo = False

        etb = eb + elb
        if verbose >= 3:
            print('expected total data bytes:            %s' % str(etb))
            print('expected number of resolutions:       %s' % str(nres))
            print('expected bytes in highest resolution: %s' % str(eb))
            print('expected bytes in other resolutions:  %s' % str(elb))
            print('data bytes in file:                   %s' % str(ab))
            if displaySectionInfo:
                print('expected total secs:                  %s' % str(enSecs))
                print('file has total secs:                  %s' % str(anSecs))

        if etb == ab:
            if verbose >= 2:
                print('OK')
            return 1
        elif etb < ab:
            if verbose >= 1:
                print('* have %s extra bytes in file' % str(ab - etb))
                if displaySectionInfo:
                    print('* have %.2f extra hidden sections in file' %
                          (anSecs - enSecs))
            return 0
        elif eb <= ab:
            if verbose >= 1:
                print('* lower resolution data truncated by %s bytes' %
                      str(etb - ab))
            return 0
        else:
            if verbose >= 1:
                print('* highest resolution data truncated by %s bytes' %
                      str(eb - ab))
                print('* file missing %.2f sections of highest resolution' %
                      (enSecs - anSecs))
                print('PLEASE SET shape to  %s sections !!! ' %
                      str(int(anSecs)))
            return 0


    def info(self):
        """Print useful information from header."""

        hdrInfo(self.hdr)


    def data_withMrc(self, fn):
        """Return the image data as a NumPy array.  The Mrc attribute of
        the returned array is this Mrc object.

        Positional parameters:
        fn -- Is for compatibility with bindFile() and previous versions
        of Mrc.py; its value is not used.
        """

        class ndarray_inMrcFile(N.ndarray):
            def __array_finalize__(self, obj):
                self.Mrc = getattr(obj, 'Mrc', None)

        data = self.data
        data.__class__ = ndarray_inMrcFile
        ddd = weakref.proxy(data)
        self.data = ddd
        data.Mrc = self

        return data


    def close(self):
        """Deletes the memory map of the MRC file.

        Implicitly commits to disk any changes made.
        """

        # As of NumPy 1.9, memmap no longer has a close method.  Instead
        # use del for all versions.
        if hasattr(self, 'm'):
            del self.m


###########################################################################
###########################################################################
###########################################################################
###########################################################################

def open(path, mode='r'):
    """Return a Mrc2 object for a MRC file with the given file name.

    Positional parameters:
    path -- Is the name of the file to open.  May be None to create
    a temporary file.  For more information, look at the documentation
    for Mrc2.__init__ .

    Keyword parameters:
    mode -- Controls how the file is opened.  Use 'r' for read-only
    access, 'r+' for read and write access to an existing file, 'w'
    for write access which will ignore the contents of the file if it
    already exists, and 'w+' for read and write access which will
    ignore the contents of the file if it already exists.  For more
    information, look at the documentation for Mrc2.__init__ .

    Return value:
    Returns a Mrc2 object for the file.  When mode is 'r' or 'r+',
    the header and extended header, if it is present, have been read,
    and the file is positioned at the start of the image data for
    the first image.  When mode is 'w' or 'w+', the header fields
    have been set to default values, the file is positioned at the
    start of the header, and a call to setHdrForShapeType() or
    initHdrForArr() for the object will be necessary before reading
    or writing image data.
    """

    return Mrc2(path, mode)


def load(fn):
    """Copy the highest resolution data set from a MRC file to memory.

    Positional parameters:
    fn -- Is the name of the MRC file from which to read.

    Return value:
    Returns a three-dimensional NumPy array with the data from the file.
    The array is dimensioned (nsections, ny, nx).  The ordering of elements
    in the first dimension is the same as in the file.  The array is not
    mapped to the file so changes to the array will not affect the file
    and changes to the file will not affect the array.
    """

    m = open(fn)
    a = m.readStack(m.hdr.Num[2])
    return a


def save(a, fn, ifExists='ask', zAxisOrder=None,
         hdr=None, hdrEval='',
         calcMMM=True,
         extInts=None, extFloats=None):
    """Save the contents of a NumPy array as a MRC-formatted file.

    Positional parameters:
    a -- Is the NumPy array to save.  a must have one to five dimensions.
    The type of a must be natively supported by the MRC format.
    Acceptable types include 32-bit floating point, unsigned and signed
    16-bit integer, 64-bit complex, and 8-bit signed integer.
    fn -- Is the name of the file to be written.

    Keyword parameters:
    ifExists -- Controls what happens if there is already a file with
    the given name.  ifExists should be 'ask', 'overwrite', or 'raise'.
    If the first letter of ifExists is 'o', the existing file will be
    overwritten.  If the first ifExists is 'a', execution will be
    suspended until there is a response to an interactive prompt about
    overwriting the file.  If the prompt is anwered with 'y' or 'Y'
    followed by a return, the file will be overwritten.  Any other
    answer followed by a return will generate an exception.  If the
    first letter of ifExists is neither 'o' nor 'a', an exception will
    be generated if the file already exists.

    zAxisOrder -- Controls how the dimensions besides the last two
    of the array are translated to the z, wavelength, and time axes
    of the file.  The ordering of the dimensions in zAxisOrder is
    from slowest varying (the first dimension of a), to fastest
    varying.  When zAxisOrder is None, it is equivalent to 'z' if
    the array has three dimensions, 'tz' if the array has four
    dimensions and to 'tzw' in all other cases.  Any ' ', '-', '.',
    or ',' characters in zAxisOrder are treated as delimiters and
    are stripped out.  The remaining characters are converted to
    lower case.  In the case where the array has three dimensions,
    only the last character in zAxisOrder after the delimiters have
    been stripped has an effect.  It affects the number of z
    samples (nz), number of wavelengths (nw), number of time points
    (nt), and image sequence as follows:
           nz          nw          nt          Priism sequence
    'z'    a.shape[0]  1           1           ZTW
    'w'    1           a.shape[0]  1           WZT
    't'    1           1           a.shape[0]  ZTW

    If the array has four dimensions, only the last two characters
    of zAxisOrder have an effect:
           nz          nw          nt          Priism sequence
    'zw'   a.shape[0]  a.shape[1]  1           WZT
    'zt'   not supported by Priism sequence types
    'wz'   a.shape[1]  a.shape[0]  1           ZTW
    'wt'   1           a.shape[0]  a.shape[1]  ZTW
    'tz'   a.shape[1]  1           a.shape[0]  ZTW
    'tw'   1           a.shape[1]  a.shape[0]  WZT

    If the array has five dimensions, only the last three characters
    of zAxisOrder have an effect:
           nz          nw          nt          Priism sequence
    'zwt'  not supported by Priism sequence types
    'ztw'  not supported by Priism sequence types
    'wtz'  a.shape[2]  a.shape[0]  a.shape[1]  ZTW
    'wzt'  not supported by Priism sequence types
    'tzw'  a.shape[1]  a.shape[2]  a.shape[0]  WZT
    'twz'  a.shape[2]  a.shape[1]  a.shape[0]  ZWT

    hdr -- If not None, should act like an instance of the HdrPriism
    or Hdr2014Priism classes, i.e. a value returned by makeHdrArray()
    or implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.  The values of
    all fields in hdr will be copied to the output file except for
    the number of samples ('Num' field; bytes 1 - 12), the pixel
    type ('PixelType' field; bytes 13 - 16), and the number of
    bytes in the extended header ('next' field; bytes 93 - 96).
    If both hdr and zAxisOrder are not None, the values for the
    number of wavelengths ('NumWaves' field; bytes 197 - 198),
    number of time points ('NumTimes' field; bytes 181 - 182),
    and image sequence ('ImgSequence' field; bytes 183 - 184)
    will be overwritten by the values in hdr.

    calcMMM -- If True, the minimum value for each wavelength,
    the maximum value for each wavelength, and the median value for
    first wavelength will be calculated and stored in the header.
    Those calculated values will overwrite the values set
    by hdr.

    hdrEval -- If not an empty string or None, hdrEval will be
    executed with Python's exec in a context where all global
    variables are accessible and the local variables are those of
    the calling function augmented with a variable named 'hdr'.
    That local variable represents the header to be saved to
    the file.  It is a header as returned by makeHdrArray() or
    implement_hdr().  hdrEval is executed after any changes to
    the header made due to the zAxisOrder, hdr, calcMMM,
    extInts, and extFloats parameters.

    extInts -- If not None, will be used to initialize the
    integer fields in the extended header.  When not None, extInts
    must be a NumPy array.  The size of the last dimension will
    be used as the number of integers per section.  The product of
    the sizes for the remaining dimensions will be used as the
    number of sections of data available.  When extInts is None
    and extFloats is not None, there will be no integer entries
    in the extended header.

    extFloats -- If not None, will be used to initialize the
    floating-point fields in the extended header.  When not None,
    extFloats must be a NumPy array.  The size of the last dimension
    will be used as the number of floating-point values per section.
    The product of the sizes for the remaining dimensions will be
    used as the number of sections of data available.  When extFloats
    is None and extInts is not None, there will be no floating-point
    entries in the extended header.  If there are a different number
    of sections of data available from extInts and extFloats, the
    size of the extended header will be set based on the larger
    number of sections and the unspecified values will be filled with
    zeros.
    """

    if os.path.exists(fn):
        if ifExists[0] == 'o':
            pass
        elif ifExists[0] == 'a':
            try:
                # First try the Python 2 way for this.
                answer = raw_input('overwrite?')
            except NameError:
                # Now try the Python 3 one; note that input() has a different
                # meaning in Python 2.
                answer = input('overwrite?')
            yes = answer.lower() == 'y'
            if not yes:
                raise RuntimeError('not overwriting existing file "%s"' % fn)
        else:
            raise RuntimeError('not overwriting existing file "%s"' % fn)

    m = Mrc2(fn, mode='w')

    m.initHdrForArr(a, zAxisOrder)
    if hdr is not None:
        m.hdr = initHdrArrayFrom(m.hdr, hdr)

    if calcMMM:
        wAxis = axisOrderStr(m.hdr).find('w')
        if wAxis < 0:
            m.hdr.mmm1 = computeMinMaxMedian(N.real(a))
        else:
            nw = m.hdr.NumWaves
            m.hdr.mmm1 = computeMinMaxMedian(N.real(a.take((0,), wAxis)))
            if nw >=2:
                m.hdr.mm2 = computeMinMax(N.real(a.take((1,), wAxis)))
            if nw >=3:
                m.hdr.mm3 = computeMinMax(N.real(a.take((2,), wAxis)))
            if nw >=4:
                m.hdr.mm4 = computeMinMax(N.real(a.take((3,), wAxis)))
            if nw >=5:
                m.hdr.mm5 = computeMinMax(N.real(a.take((4,), wAxis)))

    if extInts is not None:
        numints = extInts.shape[-1]
        numextsec_int = extInts.size // numints
    else:
        numints = 0
        numextsec_int = 0
    if extFloats is not None:
        numfloats = extFloats.shape[-1]
        numextsec_float = extFloats.size // numfloats
    else:
        numfloats = 0
        numextsec_float = 0
    if ((numints > 0 and numextsec_int > 0) or
        (numfloats > 0 and numextsec_float > 0)):
        m.makeExtendedHdr(numints, numfloats,
                          nSecs=max(numextsec_int, numextsec_float))
        if numints == 1:
            m.extInts[0:numextsec_int] = N.ravel(extInts)
            m.extInts[numextsec_int:] = 0
        elif numints > 1:
            m.extInts[0:numextsec_int, 0:numints] = (
                N.reshape(extInts, (numextsec_int, numints)))
            m.extInts[numextsec_int:, 0:numints] = 0
        if numfloats == 1:
            m.extFloats[0:numextsec_float] = N.ravel(extFloats)
            m.extFloats[numextsec_float:] = 0.0
        elif numfloats > 1:
            m.extFloats[0:numextsec_float, 0:numfloats] = (
                N.reshape(extFloats, (numextsec_float, numfloats)))
            m.extFloats[numextsec_float:, 0:numfloats] = 0.0

    if hdrEval:
        fr = sys._getframe(1)
        loc = {'hdr':m.hdr}
        loc.update(fr.f_locals)
        glo = fr.f_globals
        exec(hdrEval, glo, loc)
    m.writeHeader()
    m.writeExtHeader(seekTo0=True)
    m.writeStack(a)
    m.close()


def computeMinMaxMedian(array):
    """Compute statistics for save().

    This function was added in version 0.1.0 of Mrc.py as included
    with Priism.

    Positional parameters:
    array -- Is a NumPy array.

    Return value:
    Returns a tuple with the minimum, maximum, and median of the array.
    If array is structured, the values returned are for the first
    component.
    """

    if array.dtype.fields is None:
        return (N.min(array), N.max(array), N.median(array))
    else:
        return (N.min(array[array.dtype.names[0]]),
                N.max(array[array.dtype.names[0]]),
                N.median(array[array.dtype.names[0]]))


def computeMinMax(array):
    """Compute statistics for save().

    This function was added in version 0.1.0 of Mrc.py as included
    with Priism.

    Positional parameters:
    array -- Is a NumPy array.

    Return value:
    Returns a tuple with the minimum and maximum of the array.  If
    array is structured, the values returned are for the first
    component.
    """

    if array.dtype.fields is None:
        return (N.min(array), N.max(array))
    else:
        return (N.min(array[array.dtype.names[0]]),
                N.max(array[array.dtype.names[0]]))


###########################################################################
###########################################################################
###########################################################################
###########################################################################

class Mrc2:
    """Provide access to MRC files.

    Unlike the Mrc class, does not use a memory map to access the file.
    Actively manages the header and extended header data.

    If x is a Mrc2 instance. x.hdr has the header fields, and x.e has
    the extended header as an array of unsigned 8-bit integers.  If the
    extended header size is greater than zero, the extended header is
    in Priism's format, and the number of integers or floating-point
    values per section in the extended header is greater than zero,
    x.extInts has the integer values from the extended header,
    x.extFloats has the floating-point values from the extended header,
    and x.extSym is None.  If the extended header has symmetry
    information, x.extInts is None, x.extFloats is None, and x.extSym
    is an array of 80 character records for the symmetry information.
    Any other cases will have x.extInts equal to None, x.extFloats
    equal to None, and x.extSym equal to None.

    Some internal state depends on the values in hdr.Num and
    hdr.PixelType.  To modify those fields, the recommended procedure
    is to do so indirectly with initHdrForArr() or setHdrForShapeType().
    If you do modify those fields directly, call _initWhenHdrArraySet()
    so that the internal state is consistent with your changes.  For
    the fields related to the extended header, hdr.NumIntegers,
    hdr.NumFloats, and hdr.next, there's no public way to modify
    those directly while maintaining a consistent internal state
    for the extended header.  Use makeExtendedHdr() or
    makeSymmetryInfo() to modify those fields.

    For the image data, provides the functions seekSec(), readSec(),
    readStack(), writeSec(), and writeStack() to position the file
    at a given section, read image data, or write image data.

    Version 0.1.0 of Mrc.py as included with Priism changed the
    conditions for when the extInts and extFloats attributes are
    set.
    """

    def __init__(self, path, mode='r'):
        """Initialize the Mrc2 object.

        If mode is 'r' or 'r+', reads the header and, if present, the
        extended header, and positions the file at the start of the image
        data for the first section.  When mode is 'w' or 'w+', the
        header fields are set to the default values, the file is
        positioned at the start of the header, and a call to
        setHdrForShapeType() or initHdrForArr() will be necessary
        before reading or writing image data.

        Positional parameters:
        path -- Is the name of the file to use.  If path is None, a
        temporary file will be generated, and that file will be deleted
        when the close() method is called.  A value of None for path
        will only work if the mode parameter is set to 'w' or 'w+'.
        The _path attribute of the created object will be set to path
        if path is not None; otherwise, it will be set to the name of
        the temporary file.  Allowing None for path was added in
        version 0.1.0 of Mrc.py as included with Priism.

        Keyword parameters:
        mode -- Specifies how the file should be opened.  It has
        similar semantics to the mode parameter for Python's open()
        except that all modes implicitly include 'b' for working
        with binary files.  Allowed values for mode are:

              reading    writing     notes
        'r'   allowed    forbidden   path must exist when __init__ called
        'r+'  allowed    allowed     path must exist when __init__ called
        'w'   forbidden  allowed     path overwritten if it exists
        'w+'  allowed    allowed     path overwritten if it exists
        """

        if path is None:
            self._f, self._path = tempfile.mkstemp()
            self._f = os.fdopen(self._f, mode + 'b')
            self._delete_on_close = True
        else:
            self._f = builtins.open(path, mode + 'b')
            self._path = path
            self._delete_on_close = False
        self._name = os.path.basename(self._path)
        self._mode = mode

        self._hdrSize = 1024
        self._dataOffset = self._hdrSize

        self._fileIsByteSwapped = False

        if mode in ('r', 'r+'):
            self._initFromExistingFile()

            self.seekSec(0)
        else:
            self.hdr = makeHdrArray()
            self.hdr.Num = (0, 0, 0)
            self.hdr.PixelType = 1
            self.hdr.mst = (0, 0, 0)
            self.hdr.m = (1, 1, 1)
            self.hdr.d = (1.0, 1.0, 1.0)
            self.hdr.angle = (90.0, 90.0, 90.0)
            self.hdr.axis = (1, 2, 3)
            self.hdr.mmm1 = (0.0, 0.0, 0.0)
            self.hdr.nspg = 0
            self.hdr.next = 0
            self.hdr.dvid = 0xc0a0
            self.hdr.nblank = 0
            self.hdr.ntst = 0
            self.hdr.blank = 0
            self.hdr.NumIntegers = 0
            self.hdr.NumFloats = 0
            self.hdr.sub = 1
            self.hdr.zfac = 1
            self.hdr.mm2 = (0.0, 0.0)
            self.hdr.mm3 = (0.0, 0.0)
            self.hdr.mm4 = (0.0, 0.0)
            self.hdr.ImageType = 0
            self.hdr.LensNum = 0
            self.hdr.n1 = 0
            self.hdr.n2 = 0
            self.hdr.v1 = 0
            self.hdr.v2 = 0
            self.hdr.mm5 = (0.0, 0.0)
            self.hdr.NumTimes = 1
            self.hdr.ImgSequence = 0
            self.hdr.tilt = (0.0, 0.0, 0.0)
            self.hdr.NumWaves = 1
            self.hdr.wave = (0, 0, 0, 0, 0)
            self.hdr.zxy0 = (0.0, 0.0, 0.0)
            self.hdr.NumTitles = 0
            self.hdr.title = ' ' * 80

            self._shape = None
            self._shape2d = None
            self._dtype = None    # scalar data type of pixels
            self._secByteSize = 0
            self.e = N.zeros(0, dtype='u1')
            self._extHdrArray = None
            self.extInts = None
            self.extFloats = None
            self.extSym = None


    def initHdrForArr(self, arr, zAxisOrder=None):
        """Initialize the MRC header from the shape and type of a NumPy
        array.

        Positional parameters:
        arr -- Is the NumPy array whose shape and type are to be used.

        zAxisOrder -- Controls how the dimensions besides the last two
        of the array are translated to the z, wavelength, and time axes
        of the file.  The ordering of the dimensions in zAxisOrder is
        from slowest varying (the first dimension of a), to fastest
        varying.  When zAxisOrder is None, it is equivalent to 'z' if
        the array has three dimensions, 'tz' if the array has four
        dimensions and to 'tzw' in all other cases.  Any ' ', '-',
        '.', or ',' characters in zAxisOrder are treated as delimiters
        and are stripped out.  The remaining characters are converted
        to lower case.  The documentation for save() in this module
        has the details for how the zAxisOrder will set the header
        values for the image sequence type and the number of samples
        in z, wavelength, and time.
        """

        if zAxisOrder is None:
            if   arr.ndim ==3:
                zAxisOrder = 'z'
            elif arr.ndim ==4:
                zAxisOrder = 'tz'
            else:
                zAxisOrder = 'tzw'
        else:
            # remove delimiter characters '-., '
            zAxisOrder = zAxisOrder.translate(
                string.join([chr(i) for i in range(256)], ''), '-., ').lower()

        mrcmode = dtype2MrcMode(arr.dtype)
        init_simple(self.hdr, mrcmode, arr.shape,
                    isByteSwapped=self._fileIsByteSwapped)
        if arr.ndim == 1 or arr.ndim == 2:
            pass
        elif arr.ndim == 3:
            if   zAxisOrder[-1] == 'z':
                self.hdr.ImgSequence = 0
            elif zAxisOrder[-1] == 'w':
                self.hdr.ImgSequence = 1
                self.hdr.NumWaves = arr.shape[-3]
            elif zAxisOrder[-1] == 't':
                self.hdr.ImgSequence = 0
                self.hdr.NumTimes = arr.shape[-3]
            else:
                raise ValueError('unsupported axis order')
        elif arr.ndim == 4:
            if   zAxisOrder[-2:] == 'zt':
                raise ValueError('unsupported axis order; time varies '
                                 'faster than z')
            elif zAxisOrder[-2:] == 'tz':
                self.hdr.ImgSequence = 0
                self.hdr.NumTimes = arr.shape[-4]
            elif zAxisOrder[-2:] == 'wz':
                self.hdr.ImgSequence = 0
                self.hdr.NumWaves = arr.shape[-4]
            elif zAxisOrder[-2:] == 'zw':
                self.hdr.ImgSequence = 1
                self.hdr.NumWaves = arr.shape[-3]
            elif zAxisOrder[-2:] == 'tw':
                self.hdr.ImgSequence = 1
                self.hdr.NumWaves = arr.shape[-3]
                self.hdr.NumTimes = arr.shape[-4]
            elif zAxisOrder[-2:] == 'wt':
                self.hdr.ImgSequence = 0
                self.hdr.NumWaves = arr.shape[-4]
                self.hdr.NumTimes = arr.shape[-3]
            else:
                raise ValueError('unsupported axis order')
        elif arr.ndim == 5:
            if zAxisOrder[-3:] == 'wtz':
                self.hdr.ImgSequence = 0
                self.hdr.NumWaves = arr.shape[-5]
                self.hdr.NumTimes = arr.shape[-4]
            elif zAxisOrder[-3:] == 'tzw':
                self.hdr.ImgSequence = 1
                self.hdr.NumWaves = arr.shape[-3]
                self.hdr.NumTimes = arr.shape[-5]
            elif zAxisOrder[-3:] == 'twz':
                self.hdr.ImgSequence = 2
                self.hdr.NumWaves = arr.shape[-4]
                self.hdr.NumTimes = arr.shape[-5]
            else:
                raise ValueError('unsupported axis order')
        else:
             raise ValueError('unsupported array ndim')
        if self.hdr.NumWaves > 5:
            print('WARNING: more than 5 wavelengths for MRC file')

        self._initWhenHdrArraySet()


    def _initFromExistingFile(self):
        """Initialize the header for __init__ from the contents of the file."""

        self.seekHeader()
        buffer = N.fromfile(self._f, dtype='u1', count=1024)
        self.hdr = makeHdrArray(buffer, makeWeak=False)

        if hdrIsByteSwapped(self.hdr):
            self.hdr._array.dtype = self.hdr._array.dtype.newbyteorder()
            self._fileIsByteSwapped = True

        self._extHdrSize = self.hdr.next
        self._extHdrNumInts = max(0, self.hdr.NumIntegers)
        self._extHdrNumFloats = max(0, self.hdr.NumFloats)
        self._extHdrBytesPerSec = (
            (self._extHdrNumInts + self._extHdrNumFloats) * 4)
        self._dataOffset = self._hdrSize + self._extHdrSize

        if self._extHdrSize > 0:
            self.e = N.fromfile(self._f, dtype='u1', count=self._extHdrSize)
            fmt = getExtHeaderFormat(self.hdr)
        else:
            self.e = N.zeros(0, dtype='u1')
            fmt = -1
        if fmt == 0:
            nrec = self._extHdrSize // 80
            nrem = self._extHdrSize - 80 * nrec
            type_descr = [('records', '(%d,80)i1' % nrec),
                          ('extra', '%di1' % nrem)]
            self.extHdrArray = N.recarray(shape=1, dtype=type_descr,
                                          buf=self.e)
            self.extSym = ManagedTitleArray(
                self.extHdrArray.field('records')[0])
            self.extInts = None
            self.extFloats = None
        elif (fmt == 1 and
              (self._extHdrNumInts > 0 or self._extHdrNumFloats > 0)):
            nSecs = self._extHdrSize // self._extHdrBytesPerSec
            byteorder = '='
            type_descr = [
                ('int', '%s%di4' % (byteorder, self._extHdrNumInts)),
                ('float', '%s%df4' % (byteorder, self._extHdrNumFloats))]
            self._extHdrArray = N.recarray(shape=nSecs, dtype=type_descr,
                                           buf=self.e)
            if self._fileIsByteSwapped:
                self._extHdrArray = self._extHdrArray.newbyteorder()
            self.extInts = self._extHdrArray.field('int')
            self.extFloats = self._extHdrArray.field('float')
            self.extSym = None
        else:
            self._extHdrArray = None
            self.extInts = None
            self.extFloats = None
            self.extSym = None

        self._initWhenHdrArraySet()


    def _initWhenHdrArraySet(self):
        """Reset internal attributes based on size and pixel type in header."""

        nx, ny, nsecs =  self.hdr.Num
        if nx < 0:
            nx = 0
        if ny < 0:
            ny = 0
        if nsecs < 0:
            nsecs = 0
        self._shape = (nsecs, ny, nx) # todo: wavelenths , times
        self._shape2d = self._shape[-2:]
        self._dtype = MrcMode2dtype(self.hdr.PixelType)
        if self._fileIsByteSwapped:
            self._dtype = self._dtype.newbyteorder()
        self._secByteSize = self._dtype.itemsize * N.prod(self._shape2d)


    def setHdrForShapeType(self, shape, type):
        """Set the size and pixel type fields in the header.

        For a file opened in 'w' or 'w+' mode, this and initHdrForArr()
        are the two ways to make a Mrc2 object ready to read or write
        image data.

        As currently implemented, only uses the last two elements of
        shape and the product of all the remaining elements of shape
        to set the size fields in the header.  It does not modify the
        fields for the number of time points, number of wavelengths, or
        image sequence.

        Positional parameters:
        shape -- Is a tuple to specify the shape of the data to be stored
        in the MRC file.  The ith element of the tuple is the number of
        samples for the ith dimension.  The 0th dimension is the slowest
        varying. The fastest varying dimension, usually called x, is the
        last element in the tuple.  Shape should have at least two
        elements.

        type -- Is the NumPy dtype or Python type that will be used to
        represent each pixel value in the file.  If the type is not
        equivalent to one of the pixel formats supported by MRC, an
        exception will be raised.
        """

        mrcmode = dtype2MrcMode(type)
        self.hdr.PixelType =  mrcmode
        self.hdr.Num = shape[-1], shape[-2], N.prod(shape[:-2])
        self._initWhenHdrArraySet()


    def makeExtendedHdr(self, numInts, numFloats, nSecs=None):
        """Create a Priism extended header or remove the extended header.

        Will remove the extended header if nSecs is zero or both numInts
        and numFloats are zero.

        If header is in Priism's format, sets the space group to zero.
        Also sets the space group to zero if the header does not claim
        to support the exttyp field and the space group is different
        than 0, 1, or 401.

        The entries for a new header will all be zero.  If there
        already was an extended header, the resources for the previous
        extended header are released, and no attempt is made to copy
        the previous values to the new header.

        When a new extended header is created, the integer values can
        be accessed with self.extInts.  The floating point values can
        be accessed with self.extFloats.  Both are NumPy array views.
        If numInts is greater than one or is zero, the shape for
        self.extInts will be (nSecs, numInts).  If numInts is one,
        the shape for self.extFloats will be (nSecs,).  If numFloats
        is greater than one or is zero, the shape for self.extFloats
        will be (nSecs, numFloats).  If numFloats is one, the shape
        for self.extFloats will be (nSecs,).

        makeExtendedHdr() does not change the contents of the file.
        To commit changes made to the shape of the extended header,
        call writeHeader().  To commit changes made to the values
        in the extended header, call writeExtHeader().

        Positional parameters:
        numInts -- Is the number of integer values to store per section
        in the extended header.  Must be non-negative.

        numFloats -- Is the number of floating-point values to store
        per section in the extended header.  Must be non-negative.

        Keyword parameters:
        nSecs -- If not None, nSecs is the number of sections of
        storage to allocate for the extended header.  If nSecs is
        None, the number of sections allocated will be the number
        of sections from the header.
        """

        if numInts < 0 or numFloats < 0:
            raise ValueError('Number of integers or floating point '
                             'values is negative')
        if numInts > 32767 or numFloats > 32767:
            raise ValueError('Number of integers or floating point '
                             'values is too large to store in header fields')
        if nSecs is not None and nSecs < 0:
            raise ValueError('Number of sections is negative')
        if nSecs is None:
            if self._shape is None:
                nSecs = 0
            else:
                nSecs = self._shape[0]
        bytesPerSec = (numInts + numFloats) * 4
        ntot = minExtHdrSize(nSecs, bytesPerSec)
        if ntot > 2147483647:
            raise ValueError('Requested extended header size is too '
                             'large for the extended header size field')

        self._extHdrNumInts = self.hdr.NumIntegers = numInts
        self._extHdrNumFloats = self.hdr.NumFloats = numFloats
        if hdrHasExtType(self.hdr):
            self.hdr.exttyp = N.fromstring('AGAR', dtype='i1')
        else:
            if (hdrIsInPriismFormat(self.hdr) or
                (self.hdr.nspg != 0 and self.hdr.nspg != 1 and
                 self.hdr.nspg != 401)):
                self.hdr.nspg = 0
        self._extHdrBytesPerSec = bytesPerSec
        self._extHdrSize = self.hdr.next = ntot
        self._dataOffset = self._hdrSize + self._extHdrSize

        self.e = N.zeros(self._extHdrSize, dtype='u1')
        if self._extHdrSize > 0 and self._extHdrBytesPerSec > 0:
            nSecs = self._extHdrSize // self._extHdrBytesPerSec
            byteorder = '='
            type_descr = [
                ('int', '%s%di4' % (byteorder, self._extHdrNumInts)),
                ('float', '%s%df4' % (byteorder, self._extHdrNumFloats))]
            self._extHdrArray = N.recarray(nSecs, dtype=type_descr,
                                           buf=self.e)
            self.extInts = self._extHdrArray.field('int')
            self.extFloats = self._extHdrArray.field('float')
        else:
            self._extHdrArray = None
            self.extInts = None
            self.extFloats = None
        self.extSym = None


    def makeSymmetryInfo(self, nbytes, nspg=None):
        """Create the extended header for symmetry information.

        If the header is in Priism's format and the space group,
        after applying the nspg keyword, is zero, will raise a
        RuntimeError exception since files with zero for the
        space group are assumed to use Priism-style extended
        headers.

        Sets the NumIntegers and NumFloats fields to zero.

        The new extended header is filled with spaces.  If there
        already was an extended header, the resources for the
        previous extended header are released, and no attempt
        is made to copy the previous values to the new header.

        When a new extended header is created, the symmetry
        information, as an array of 80 character records, can
        be accessed with self.extSym.  self.extInts and
        self.extFloats are set to None.

        makeSymmetryInfo() does not change the contents of the file.
        To commit changes made to the shape of the extended header,
        call writeHeader().  To commit changes made to the values
        in the extended header, call writeExtHeader().

        This function was added in version 0.1.0 of Mrc.py as
        included with Priism.

        Positional parameters:
        nbytes -- Is the number of bytes to allocate.  A ValueError
        exception will be raised if nbytes is less than zero.  A
        value of zero will remove the extended header.  Note that
        values of nbytes that are not multiples of eight could lead
        to misalignment of image data if the file is memory mapped
        with the Mrc class.

        Keyword parameters:
        nspg -- If not None, the space group in the header will be
        set to the specified value.
        """

        if nbytes < 0:
            raise ValueError('Negative number of bytes requested for '
                             'extended header')
        if nbytes > 2147483647:
            raise ValueError('Requested number of bytes is too large '
                             'for the extended header size field')
        if nspg is not None:
            self.hdr.nspg = nspg
        if hdrIsInPriismFormat(self.hdr) and self.hdr.nspg == 0:
            raise RuntimeError('Used makeSymmetryInfo() when the space '
                               'group is zero')
        self.hdr.next = nbytes
        self.hdr.NumIntegers = 0
        self.hdr.NumFloats = 0
        if hdrHasExtType(self.hdr):
            self.hdr.exttyp = N.fromstring('MRC0', dtype='i1')
        self._extHdrSize = nbytes
        self._extHdrNumInts = 0
        self._extHdrNumFloats = 0
        self._extHdrBytesPerSec = 0
        self._dataOffset = self._hdrSize + self._extHdrSize
        if self._extHdrSize > 0:
            self.e = N.empty(self._extHdrSize, dtype='u1')
            # ASCII for space.
            self.e[:] = 32
            nrec = self._extHdrSize // 80
            nrem = self._extHdrSize - 80 * nrec
            type_descr = [('records', '(%d,80)i1' % nrec),
                          ('extra', '%di1' % nrem)]
            self._extHdrArray = N.recarray(shape=1, dtype=type_descr,
                                           buf=self.e)
            self.extSym = ManagedTitleArray(
                self._extHdrArray.field('records')[0])
        else:
            self.e = N.zeros(0, dtype='u1')
            self._extHdrArray = None
            self.extSym = None
        self.extInts = None
        self.extFloats = None


    def makeGenericExtendedHdr(self, nbytes, fmt):
        """Allocate space for an extended header that is not
        for symmetry information and does not use Priism's
        format.

        The bytes in the new extended header are set to zero.
        If there already was an extended header, the resources
        for the previous extended header are released, and no
        attempt is made to copy the previous values to the new
        header.

        The new extended header, as an array of unsigned bytes,
        can be accessed through self.e.  Sets self.extInts,
        self.extFloats, and self.extSym to None.

        makeGenericExtendedHdr() does not change the contents
        of the file.  To commit the changes made to the shape
        or format of the extended header, call writeHeader().
        To commit changes made to the values in the extended
        header, call writeExtHeader().

        This function was added in version 0.1.0 of Mrc.py as
        included with Priism.

        Positional parameters:
        nbytes -- Is the number of bytes to allocate.  If
        nbytes is zero, the extended header will be removed.
        A ValueError exception will be raised if nbytes is
        less than zero.  Note that values of nbytes that are
        not multiples of eight could lead to misalignment of
        image data if the file is memory mapped with the Mrc
        class.
        fmt -- Is a four character ASCII string describing
        the format of the extended header.  A ValueError
        exception will be raised if fmt is 'AGAR', 'MRC0',
        of 'CCP4'.  Use makeExtendedHdr() or makeSymmetryInfo()
        to create extended headers for those formats.  A
        RuntimeError exception will be raised if the header
        format does not store a string describing the
        extended header format.
        """

        if nbytes < 0:
            raise ValueError('Negative number of bytes requested for '
                             'extended header')
        if nbytes > 2147483647:
            raise ValueError('Requested number of bytes is too large '
                             'for the extended header size field')
        if fmt == 'AGAR':
            raise ValueError('Use makeExtendedHdr() to create an extended '
                             'with the Priism format')
        if fmt == 'MRC0' or fmt == 'CCP4':
            raise ValueError('Use makeSymmetryInfo() to create an extended '
                             'header with symmetry information')
        if len(fmt) != 4:
            raise ValueError('fmt is not a four character string')
        if not hdrHasExtType(self.hdr):
            raise RuntimeError('Header does not store a string code for the '
                               'extended header format')

        self.hdr.exttyp = N.fromstring(fmt, dtype='i1')
        self.hdr.next = nbytes
        self._extHdrSize = nbytes
        self._extHdrNumInts = 0
        self._extHdrNumFloats = 0
        self._extHdrBytesPerSec = 0
        self._dataOffset = self._hdrSize + self._extHdrSize
        self.e = N.zeros(nbytes, dtype='u1')
        self.extInts = None
        self.extFloats = None
        self.extSym = None


    def setTitle(self, s, i=-1, push=False, truncate=False):
        """Set a title in the MRC header.

        This function was added in version 0.1.0 of Mrc.py as
        included with Priism.  That version also allows
        calling setTitle() directly on the header:
        self.hdr.setTitle().

        Positional parameters:
        s -- Is the character string for the title.  If s is longer
        than 80 characters and truncate is False, a ValueError
        exception will be raised.  Since no byte swapping is done
        for the titles in the header, s should be encoded in ASCII
        or another format that does not use multibyte characters.

        Keyword parameters:
        i -- Is the index of the title to set.  If i is less than
        zero, the last title not in use will be set.  If i is less
        than zero and all the titles are in use and push is False
        or i is greater than 9, a ValueError exception will be
        raised.

        push -- If True, i is less than zero, and all titles are
        in use, titles will be pushed down before assigning s
        to the last title.  That will discard the first title and
        title[k] (for k greater than and equal to 0 and less than
        9) will be title[k+1] from before the change.

        truncate -- If True, only use the first 80 characters from
        s.
        """

        self.hdr.setTitle(s, i=i, push=push, truncate=truncate)


    def axisOrderStr(self, onlyLetters=True):
        """Return a string indicating the ordering of dimensions.

        x, y, z, w, and t will appear at most once in the string, and
        at least three of them will be present.  The letters that do
        appear will be present in order from slowest varying to
        fastest varying.  The values for the axis field in the header
        do not affect the result.

        This function was added in version 0.1.0 of Mrc.py as
        included with Priism.

        Keyword parameters:
        onlyLetters -- If True, only the letters for the dimensions
        will appear in the string.  If False, the first character
        of the string will be '[', the last character of the string
        will be ']', and a letter for a dimension will be preceded
        by a comma if it is not the first, slowest-varying dimension.
        """

        return axisOrderStr(self.hdr, onlyLetters)


    def info(self):
        """Print useful information from header."""

        hdrInfo(self.hdr)


    def close(self):
        """Close the file associated with the Mrc2 object.  Delete that
        file if it was created as a temporary file.
        """

        self._f.close()
        if self._delete_on_close:
            os.remove(self._path)


    def flush(self):
        """Flush any changes to the Mrc2 object's file to disk."""

        self._f.flush()


    def seekSec(self, i):
        """Seek to the start of the image data for a given section.

        Positional parameters:
        i -- Is the 0-based section index.
        """

        if self._secByteSize == 0:
            raise ValueError('not inited yet - unknown shape, type')
        self._f.seek(self._dataOffset + i * self._secByteSize)


    def seekHeader(self):
        """Seek to the start of the MRC header."""

        self._f.seek(0)


    def seekExtHeader(self):
        """Seek to the start of the MRC extended header."""

        self._f.seek(self._hdrSize)


    def readSec(self, i=None):
        """Read one image from the MRC file.

        Keyword parameters:
        i -- If i is None, starts the read at the current position for the
        file.  If i is not None, seeks to the start of section i and reads
        from that location.

        Return value:
        Returns a two-dimensional NumPy array with the image data.  The
        shape of the array is (self.hdr.Num[1], self.hdr.Num[0]).  The
        format for each image element is the same as in the file.
        """

        if i is not None:
            self.seekSec(i)

        a = N.fromfile(self._f, self._dtype, N.prod(self._shape2d))
        a.shape = self._shape2d
        return a


    def writeSec(self, a, i=None):
        """Write image data to the MRC file.

        Positional parameters:
        a -- Is a NumPy array with the data to write.  The format for each
        sample should match, up to byte order, the format for data values
        in the file.  No checks are done on the dimensions of a, so any
        amount of data can be written.  To best ensure compatibility with
        future versions of Mrc2, a should be two-dimensional with a shape
        of (self.Hdr.Num[1], self.hdr.Num[0]).

        Keyword parameters:
        i -- If i is None, the write starts at the current position for
        the file.  If i is not None, seeks to the start of section i and
        starts the write at that location.
        """

        if a.dtype.type != self._dtype.type:
            raise TypeError('type of data, %s, to write does not match '
                            'type, %s, in header' %
                            (a.dtype.name, self._dtype.name))
        if self._fileIsByteSwapped != isSwappedDtype(a.dtype):
            v = a.byteswap()
        else:
            v = a
        if i is not None:
            self.seekSec(i)

        return v.tofile(self._f)


    def readStack(self, nz, i=None):
        """Read nz images from the MRC file.

        Positional parameters:
        nz -- Is the number of images to read.

        Keyword parameters:
        i -- If i is None, the read starts at the current position for
        the file.  If i is not None, seeks to the start of section i and
        starts the read there.

        Return value:
        Returns a three-dimensional NumPy array with the image data.
        The shape of the array is (nz, self.hdr.Num[1], self.hdr.Num[0]).
        The format for each image element is the same as in the file.
        """

        if i is not None:
            self.seekSec(i)

        a = N.fromfile(self._f, self._dtype, nz * N.prod(self._shape2d))
        a.shape = (nz,) + self._shape2d
        return a


    def writeStack(self, a, i=None):
        """Write image data to the MRC file.

        Positional parameters:
        a -- Is a NumPy array with the data to write.  The format for each
        sample should match, up to byte order, the format for data values
        in the file.  No checks are done on the dimensions of a, so any
        amount of data can be written.  To best ensure compatibility with
        future versions of Mrc2, a should have at least two dimensions and
        the shape of the last two dimensions should be (self.Hdr.Num[1],
        self.Hdr.Num[0]).

        Keyword parameters:
        i -- If i is None, the write starts at the current position for
        the file.  If i is not None, seeks to the start of section i and
        starts the write at that location.
        """

        if a.dtype.type != self._dtype.type:
            raise TypeError('type of data, %s, to write does not match '
                            'type, %s, in header' %
                            (a.dtype.name, self._dtype.name))
        if self._fileIsByteSwapped != isSwappedDtype(a.dtype):
            v = a.byteswap()
        else:
            v = a
        if i is not None:
            self.seekSec(i)

        return v.tofile(self._f)


    def writeHeader(self, seekTo0=False):
        """Write the 1024 byte MRC header to the file.

        Keyword parameters:
        seekTo0 -- If True, the file position will be set after writing the
        header to the start of the first section's image data.
        """

        self.seekHeader()
        self.hdr._array.tofile(self._f)
        if seekTo0:
            self.seekSec(0)


    def writeExtHeader(self, seekTo0=False):
        """Write the extended header to the file.

        Keyword parameters:
        seekTo0 -- If True, the file position will be set after writing the
        extended header to the start of the first section's image data.
        """

        self.seekExtHeader()
        self.e.tofile(self._f)
        if seekTo0:
            self.seekSec(0)


###########################################################################
###########################################################################
###########################################################################
###########################################################################

class HdrBase(object):
    """Represents a MRC header without extensions.

    Only provides access to the parts of the header that are in the
    original formulation of the MRC format.  Those fields are listed
    below and can be accessed as x.name_of_field where x is an instance
    of this class.  For example, x.Num would access the number of
    columns, rows, and sections.

    This class was added in version 0.1.0 of Mrc.py as included with
    Priism.
    """

    # Use __slots__ so any assignments to unspecified attributes or
    # properties will give an AttributeError exception.
    __slots__ = ('_array',)


    def __init__(self, hdrArray):
        """Initialize the HdrBase object.

        Positional parameters:
        hdrArray -- Is a MRC header as a NumPy structured array with
        one element.  The dtype of the array will have to compatible
        with whatever subclass of HdrBase is being instantiated.
        A NumPy structured array with a dtype of numpy.dtype(mrcHdr_dtype)
        is compatible with HdrBase and HdrPriism.  A NumPy structured
        array with a dtype of numpy.dtype(mrc2014Hdr_dtype) is compatible
        with HdrBase, Hdr2014 and Hdr2014Priism.
        """

        self._array = hdrArray


    def _getNum(self):
        """Is three integers as a NumPy array.  The first is the number
        of columns; i.e. the number of samples in the fastest varying
        dimension as stored in the file.  The second is the number of
        rows; i.e. the number of samples in the second-fastest varying
        dimension as stored in the file.

        Occupies bytes 1 - 12 (numbered from one) in the header.
        """
        return self._array['Num'][0]
    def _setNum(self, value):
        self._array['Num'][0] = value
    Num = property(_getNum, _setNum)


    def _getPixelType(self):
        """Is an integer code for the format of the image data.
        For supported Python types or NumPy dtypes, Mrc.dtype2MrcMode()
        can compute the code.  Codes recognized by this software are (zero
        through 4 are part of the original formulation of the format):
        0:    signed (two's complement) 8-bit integer
        1:    signed (two's complement) 16-bit integer
        2:    32-bit IEEE floating-point value
        3:    real and imaginary parts of a complex value; both as signed
              (two's complement) 16-bit integers
        4:    real and imaginary parts of a complex value; both as 32-bit
              IEEE floating-point values
        5:    nothing currently standardized; treated as signed (two's
              complement) 16-bit integer
        6:    unsigned 16-bit integer
        7:    nothing currently standardized; treated as signed (two's
              complement) 32-bit integer
        101:  nothing currently standardized; treated as unsigned 4-bit
              integer; if the number of columns is odd, each row has an
              extra unused 4 bits at the end so rows start on 8-bit
              boundaries

        Occupies bytes 13 - 16 (numbered from one) in the header.
        """
        return self._array['PixelType'][0]
    def _setPixelType(self, value):
        self._array['PixelType'][0] = value
    PixelType = property(_getPixelType, _setPixelType)


    def _getmst(self):
        """Is three integers as a NumPy array.  For crystallographic
        data, these are the location of the first column, first row,
        and first section in the unit cell.  For non-crystallographic
        data, these are usually used to describe the dataset's
        relationship to another, larger, dataset from which it was
        drawn; these values might then be the indices in that larger
        dataset for the first sample in this dataset.

        Occupies bytes 17 - 28 (numbered from one) in the header.
        """
        return self._array['mst'][0]
    def _setmst(self, value):
        self._array['mst'][0] = value
    mst = property(_getmst, _setmst)


    def _getm(self):
        """Is three integers as a NumPy array.  For crystallographic
        data, holds the number of samples along the x, y, and z axes,
        respectively, of the unit cell.  For non-crystallographic data,
        it is common to use either ones for the values in m or to set
        the values equal to the values in Num.  An exception to that
        is for MRC 2014 files.  Then, m[2] is used to store the number
        of z slices per volume when the space group is 1 or 401.

        Occupies bytes 29 - 40 (numbered from one) in the header.
        """
        return self._array['m'][0]
    def _setm(self, value):
        self._array['m'][0] = value
    m = property(_getm, _setm)


    def _getd(self):
        """Is three floating-point values as a NumPy array.  For
        crystallographic data, they are the dimensions of the unit
        cell in Angstroms.  For non-crystallographic data, the
        elements of d are the elements of m times the spacing
        between samples in that dimension so that d[axis[0]] /
        m[axis[0]] is the spacing between samples in a column,
        d[axis[1]] / m[axis[1]] is the spacing between samples
        in a row, and d[axis[2]] / m[axis[2]] is the spacing
        between sections.  By convention, EM data uses
        Angstroms for the spacing and optical data uses microns.

        Occupies bytes 41 - 52 (numbered from one) in the header.
        """
        return self._array['d'][0]
    def _setd(self, value):
        self._array['d'][0] = value
    d = property(_getd, _setd)


    def _getangle(self):
        """Is three floating-point values as a NumPy array.  They are
        the angles, in degrees, between the axes of the unit cell.  For
        non-crystallographic data, these are all set to 90.

        Occupies bytes 53 - 64 (numbered from one) in the header.
        """
        return self._array['angle'][0]
    def _setangle(self, value):
        self._array['angle'][0] = value
    angle = property(_getangle, _setangle)


    def _getaxis(self):
        """Is three integers as a NumPy array.  The first is the axis
        (1 for x; 2 for y; 3 for z) that corresponds to columns in the
        file.  The second is the axis that corresponds to rows in the
        file.  The third is the axis that corresponds to sections in
        the file.

        Occupies bytes 65 - 76 (numbered from one) in the header.
        """
        return self._array['axis'][0]
    def _setaxis(self, value):
        self._array['axis'][0] = value
    axis = property(_getaxis, _setaxis)


    def _getmmm1(self):
        """Is three floating-point values as a NumPy array.  The first
        is the minimum density value.  The second is the maximum density
        value.  The third is the mean density value.  For files using
        Priism's format, these statistics are for the data from the
        first wavelength.

        Occupies bytes 77 - 88 (numbered from one) in the header.
        """
        return self._array['mmm1'][0]
    def _setmmm1(self, value):
        self._array['mmm1'][0] = value
    mmm1 = property(_getmmm1, _setmmm1)


    def _getnspg(self):
        """Is an integer to hold the space group for crystallographic
        data.  For non-crystallographic data, the convention is to
        set the space group to zero, one, or 401.  Priism's version
        of MRC uses zero.  Image2000 and MRC 2014 use zero for single
        images or image sequences, one for cases where the images
        represent a volume, and 401 for volume stacks.

        Occupies bytes 89 - 92 (numbered from one) in the header.
        """
        return self._array['nspg'][0]
    def _setnspg(self, value):
        self._array['nspg'][0] = value
    nspg = property(_getnspg, _setnspg)


    def _getnext(self):
        """Is the number of bytes in the extended header.

        Occupies bytes 93 - 96 (numbered from one) in the header.
        """
        return self._array['next'][0]
    def _setnext(self, value):
        self._array['next'][0] = value
    next = property(_getnext, _setnext)


    def _getNumTitles(self):
        """Is an integer for the number of titles used.  Up to 10
        titles can be stored.

        Occupies bytes 221 - 224 (numbered from one) in the header.
        """
        return self._array['NumTitles'][0]
    def _setNumTitles(self, value):
        self._array['NumTitles'][0] = value
    NumTitles = property(_getNumTitles, _setNumTitles)


    def _gettitle(self):
        """Are the titles as a ten element array where each element is
        an eighty character string.

        Occupies bytes 225 - 1024 (numbered from one) in the header.
        """
        return ManagedTitleArray(self._array['title'][0])
    def _settitle(self, value):
        a = ManagedTitleArray(self._array['title'][0])
        if isStringLike(value):
            nper = self._array['title'][0].shape[1]
            if len(value) >= len(a) * nper:
                for i in range(0, len(a)):
                    a[i] = value[i*80:(i+1)*80]
            else:
                for i in range(0, len(a)):
                    a[i] = value
        elif hasattr(value, '__len__'):
            if len(value) == len(a):
                for i in range(0, len(a)):
                    a[i] = value[i]
            else:
                raise ValueError('collection assigned to title does not have '
                                 ' %d elements' % len(a))
        else:
            raise TypeError('invalid type for assignment to title')
    title = property(_gettitle, _settitle)


    def getSpacing(self):
        """Return spacing between samples.

        By convention, the units for the spacing are microns for
        optical microscope data and Angstroms for electron microscope
        data.

        If the axis values in the header are invalid or a value in
        m in the header is zero, the spacing values are not well
        defined.

        Returns a three element NumPy array.  The first value is
        the spacing in the direction set by axis[0], the second value
        is the spacing in the direction set by axis[1], and the
        third value is the spacing in the direction set by axis[2].
        """

        r = N.empty(3, dtype='f4')
        m = self.m
        d = self.d
        ax = self.axis
        for i in range(0, 3):
            j = ax[i] - 1
            if (j < 0 or j > 2 or m[j] == 0):
                r[i] = 0.0
            else:
                r[i] = d[j] / m[j]
        return r


    def setSpacing(self, d0, d1, d2):
        """Set spacing between samples.

        By convention, the units for the spacing are microns for
        optical microscope data and Angstroms for electron microscope
        data.

        If the axis values in the header are invalid or a value in
        m in the header is zero, the spacing values are not well
        defined.

        Positional parameters:
        d0 -- Is the spacing between samples in the direction set by
        axis[0].

        d1 -- Is the spacing between samples in the direction set by
        axis[1].

        d2 -- Is the spacing between samples in the direction set by
        axis[2].
        """

        m = self.m
        d = self.d
        ax = self.axis
        j = ax[0] - 1
        if j >= 0 and j < 3:
            d[j] = d0 * m[j]
        j = ax[1] - 1
        if j >= 0 and j < 3:
            d[j] = d1 * m[j]
        j = ax[2] - 1
        if j >= 0 and j < 3:
            d[j] = d2 * m[j]


    def clearTitles(self):
        """Set the number of titles to zero and fill the titles with spaces."""

        self.NumTitles = 0
        self.title = ' ' * 80


    def setTitle(self, s, i=-1, push=False, truncate=False):
        """Set a title in the MRC header.

        Positional parameters:
        s -- Is the character string for the title.  If s is longer
        than 80 characters and truncate is False, a ValueError
        exception will be raised.  Since no byte swapping is done
        for the titles in the header, s should be encoded in ASCII
        or another format that does not use multibyte characters.

        Keyword parameters:
        i -- Is the index of the title to set.  If i is less than
        zero, the last title not in use will be set.  If i is less
        than zero and all the titles are in use and push is False
        or i is greater than 9, a ValueError exception will be
        raised.

        push -- If True, i is less than zero, and all titles are
        in use, titles will be pushed down before assigning s
        to the last title.  That will discard the first title and
        title[k] (for k greater than and equal to 0 and less than
        9) will be title[k+1] from before the change.

        truncate -- If True, only use the first 80 characters from
        s.
        """

        n = max(0, min(self.NumTitles, 10))
        b = N.fromstring(s, dtype='i1')
        if b.shape[0] > 80:
            if not truncate:
                raise ValueError('Mrc only support title up to 80 characters')
            b = b[0:80]
        elif b.shape[0] < 80:
            # Pad with spaces to match what is done by Priism libraries.
            b = N.concatenate((b, 32 * N.ones(80 - b.shape[0], dtype='i1')))

        if i < 0:
            i = n
            if n == 10 and push:
                for i in range(0, 9):
                    self.title[i] = self.title[i + 1]
                i = 9
                n = 9
        if i > 9:
            raise ValueError('Mrc only support up to 10 titles (0<=i<10)')
        if i >= n:
            self.NumTitles = i + 1

        self.title[i] = b


class HdrPriism(HdrBase):
    """Represents a MRC header with the Priism extensions.

    Provides access to all parts of the header.  The fields that
    are in common between the basic MRC format and Priism are
    described in HdrBase.  The fields that are Priism extensions
    listed below.

    Some properties are provided to make it easier to use instances
    of this class interchangeably with instances of Hdr2014Priism,
    especially for retrieving values without modifying them.

    This class was added in version 0.1.0 of Mrc.py as included
    with Priism.  The dynamically declared class it replaced
    had data attributes with the same names except that it did
    not have nblank, ntst, rms, and origin and it had a data
    attribute called type that was the the first two bytes of the
    nspg attribute in this class (the nspg attribute in that class
    class was only two bytes and was immediately after the type
    attribute).  The blank attribute in that class was larger,
    spanning the nblank, ntst, and blank attributes in this class.
    That class did not have getSpacing(), setSpacing(),
    clearTitles(), setTitle(), index2zwt() or zwt2index functions.
    """

    # Use __slots__ so any assignments to unspecified attributes or
    # properties will give an AttributeError exception.
    __slots__ = ()


    def index2zwt(self, i, over=False):
        """Convert section index to 3D index in z, wavelength, and time.

        Positional parameters:
        i -- Is the section index.  If less than zero, it is treated as a
        displacement from one past the end (i.e. -1 is the last valid
        section index, and -2 is the next-to-last valid section index).

        Keyword parameters:
        over -- If True, a value of i past the end will cause a ValueError
        exception to be raised.  If False, a value of i past the end will
        lead to a value for the index in the slowest-varying dimension to
        be past the end.

        Return value:
        Returns a three element tuple of the zero-based indices for z,
        wavelength, and time.
        """

        nw = max(1, self.NumWaves)
        nt = max(1, self.NumTimes)
        nz = max(1, self.Num[2] // (nw * nt))
        seq = self.ImgSequence
        if seq < 0 or seq > 2:
            seq = 0
        return index2zwt(i, nz, nw, nt, seq, over=over)


    def zwt2index(self, z, w, t, over=False):
        """Convert a 3D index in z, wavelength, and time to a section index.

        Positional parameters:
        z -- Is the index in z.  If less than zero, it is treated as a
        displacement from one past the end i.e. -1 is the last valid
        section index, and -2 is the next-to-last valid section index).

        w -- Is the wavelength index.  If less than zero, it is treated as
        a displacement from one past the end.

        t -- Is the time index.  If less than zero, it is treated as a
        displacement from one past the end.

        Keyword parameters:
        over -- If True, a ValueError exception will be raised if any of
        z, w, or t are past the end in their respective dimesions.  If
        False, a value for z, w, or t that is past the end will only
        result in a ValueError exception if that dimension is not the
        slowest-varying non-singleton dimension.

        Return value:
        Returns a zero-based index for the section.
        """

        nw = max(1, self.NumWaves)
        nt = max(1, self.NumTimes)
        nz = max(1, self.Num[2] // (nw * nt))
        seq = self.ImgSequence
        if seq < 0 or seq > 2:
            seq = 0
        return zwt2index(z, w, t, nz, nw, nt, seq, over=over)


    def _getdvid(self):
        """Is a 16-bit integer that is a code for the originator of the
        data.  Files from Priism and John Sedat's microscopes put
        0xc0a0 (-16224 as a signed integer) in this value.

        Occupies bytes 97 - 98 (numbered from one) in the header.
        """
        return self._array['dvid'][0]
    def _setdvid(self, value):
        self._array['dvid'][0] = value
    dvid = property(_getdvid, _setdvid)


    def _getnblank(self):
        """Is a 16-bit integer that is not currently reserved for any
        specific purpose.

        Occupies bytes 99 - 100 (numbered from one) in the header.
        """
        return self._array['nblank'][0]
    def _setnblank(self, value):
        self._array['nblank'][0] = value
    nblank = property(_getnblank, _setnblank)


    def _getntst(self):
        """Is an integer used to store the starting time index, i.e.
        the index in a longer time series for the first time point in
        this dataset.

        Occupies bytes 101 - 104 (numbered from one) in the header.
        """
        return self._array['ntst'][0]
    def _setntst(self, value):
        self._array['ntst'][0] = value
    ntst = property(_getntst, _setntst)


    def _getblank(self):
        """Is 24 bytes that is not currently reserved for any specific
        purpose.  Priism and this implementation do not perform any
        byte-swapping on this data.  If you store multi-byte quantities
        here, you will need to handle the byte-swapping.

        Occupies bytes 105 - 128 (numbered from one) in the header.
        """
        return self._array['blank'][0]
    def _setblank(self, value):
        self._array['blank'][0] = value
    blank = property(_getblank, _setblank)


    def _getNumIntegers(self):
        """Is the number of 4-byte integers to store per section in
        the extended header.

        Occupies bytes 129 - 130 (numbered from one) in the header.
        """
        return self._array['NumIntegers'][0]
    def _setNumIntegers(self, value):
        self._array['NumIntegers'][0] = value
    NumIntegers = property(_getNumIntegers, _setNumIntegers)


    def _getNumFloats(self):
        """Is the number of 4-byte floating-point values to store
        per section in the extended header.

        Occupies bytes 131 - 132 (numbered from one) in the header.
        """
        return self._array['NumFloats'][0]
    def _setNumFloats(self, value):
        self._array['NumFloats'][0] = value
    NumFloats = property(_getNumFloats, _setNumFloats)


    def _getsub(self):
        """Is the number of different resolutions stored in the
        dataset.

        Occupies bytes 133 - 134 (numbered from one) in the header.
        """
        return self._array['sub'][0]
    def _setsub(self, value):
        self._array['sub'][0] = value
    sub = property(_getsub, _setsub)


    def _getzfac(self):
        """For multiple resolutions, the number of z samples in a
        resolution will be the number of z samples in the next higher
        resolution divided by this integer and rounding up any remainder.

        Occupies bytes 135 - 136 (numbered from one) in the header.
        """
        return self._array['zfac'][0]
    def _setzfac(self, value):
        self._array['zfac'][0] = value
    zfac = property(_getzfac, _setzfac)


    def _getmm2(self):
        """Is two floating-point values as a NumPy array.  The first
        is the minimum density for the second wavelength, and the
        second is the maximum density for the second wavelength.

        Occupies bytes 137 - 144 (numbered from one) in the header.
        """
        return self._array['mm2'][0]
    def _setmm2(self, value):
        self._array['mm2'][0] = value
    mm2 = property(_getmm2, _setmm2)


    def _getmm3(self):
        """Is two floating-point values as a NumPy array.  The first
        is the minimum density for the third wavelength, and the
        second is the maximum density for the third wavelength.

        Occupies bytes 145 - 152 (numbered from one) in the header.
        """
        return self._array['mm3'][0]
    def _setmm3(self, value):
        self._array['mm3'][0] = value
    mm3 = property(_getmm3, _setmm3)


    def _getmm4(self):
        """Is two floating-point values as a NumPy array.  The first
        is the minimum density for the fourth wavelength, and the
        second is the maximum density for the fourth wavelength.

        Occupies bytes 153 - 160 (numbered from one) in the header.
        """
        return self._array['mm4'][0]
    def _setmm4(self, value):
        self._array['mm4'][0] = value
    mm4 = property(_getmm4, _setmm4)


    def _getImageType(self):
        """Is a 16-bit integer which is a code for the type of data in
        the file.  http://msg.ucsf.edu/IVE/IVE4_HTML/IM_ref2.html#ImageTypes
        describes the types that Priism defines and how the n1, n2, v1,
        and v2 fields are used for each type.

        Occupies bytes 161 - 162 (numbered from one) in the header.
        """
        return self._array['ImageType'][0]
    def _setImageType(self, value):
        self._array['ImageType'][0] = value
    ImageType = property(_getImageType, _setImageType)


    def _getLensNum(self):
        """For optical data, this 16-bit integer is used to indicate
        which microscope configuration was used when the data was
        collected.  Each microscope system would have its own table
        of configurations and corresponding integer codes.

        Occupies bytes 163 - 164 (numbered from one) in the header.
        """
        return self._array['LensNum'][0]
    def _setLensNum(self, value):
        self._array['LensNum'][0] = value
    LensNum = property(_getLensNum, _setLensNum)


    def _getn1(self):
        """Is a 16-bit integer whose interpretation depends on the value
        for ImageType.

        Occupies bytes 165 - 166 (numbered from one) in the header.
        """
        return self._array['n1'][0]
    def _setn1(self, value):
        self._array['n1'][0] = value
    n1 = property(_getn1, _setn1)


    def _getn2(self):
        """Is a 16-bit integer whose interpretation depends on the value
        for ImageType.

        Occupies bytes 167 - 168 (numbered from one) in the header.
        """
        return self._array['n2'][0]
    def _setn2(self, value):
        self._array['n2'][0] = value
    n2 = property(_getn2, _setn2)


    def _getv1(self):
        """Is a 16-bit integer which is used to store a floating point value,
        f, as round_to_nearest(f * 100.0).  What f is depends on the value
        for ImageType.  This Python implementation leaves the conversion
        between f and v1 to the caller.

        Occupies bytes 169 - 170 (numbered from one) in the header.
        """
        return self._array['v1'][0]
    def _setv1(self, value):
        self._array['v1'][0] = value
    v1 = property(_getv1, _setv1)


    def _getv2(self):
        """Is a 16-bit integer which is used to store a floating point value,
        f, as round_to_nearest(f * 100.0).  What f is depends on the value
        for ImageType.  This Python implementation leaves the conversion
        between f and v2 to the caller.

        Occupies bytes 171 - 172 (numbered from one) in the header.
        """
        return self._array['v2'][0]
    def _setv2(self, value):
        self._array['v2'][0] = value
    v2 = property(_getv2, _setv2)


    def _getmm5(self):
        """Is two floating-point values as a NumPy array.  The first
        is the minimum density for the fifth wavelength, and the
        second is the maximum density for the fifth wavelength.

        Occupies bytes 173 - 180 (numbered from one) in the header.
        """
        return self._array['mm5'][0]
    def _setmm5(self, value):
        self._array['mm5'][0] = value
    mm5 = property(_getmm5, _setmm5)


    def _getNumTimes(self):
        """Is an integer for the number of time points stored in the file.

        Occupies bytes 181 - 182 (numbered from one) in the header.
        """
        return self._array['NumTimes'][0]
    def _setNumTimes(self, value):
        self._array['NumTimes'][0] = value
    NumTimes = property(_getNumTimes, _setNumTimes)


    def _getImgSequence(self):
        """Is an integer code for how the sections are arranged into z,
        wavelength, and time points.  Three values are undersood
        by Priism:
        0:  Z varies fastest, followed by time, followed by wavelength
        1:  wavelength varies fastest, followed by z, followed by time
        2:  z varies fastest, followed by wavelength, followed by time

        Occupies bytes 183 - 184 (numbered from one) in the header.
        """
        return self._array['ImgSequence'][0]
    def _setImgSequence(self, value):
        self._array['ImgSequence'][0] = value
    ImgSequence = property(_getImgSequence, _setImgSequence)


    def _gettilt(self):
        """Is three floating-point values as a NumPy array.  The three
        values are a trio of rotation angles in degrees.  For image
        windows, Priism uses those angles to rotate the coordinates of
        overlayed objects to the coordinates aligned to the image axes.
        The transformation uses the first angle as rotation about the
        original +x axis of the objects with a positive angle rotating
        the +y axis towards the +z axis.  The second angle rotates about
        the +y axis from the first rotation with a positive angle
        rotating the +z axis towards the +x axis.  The third angle
        rotates about +z axis from the second rotation with a positive
        angle rotating the +x axis towards the +y axis.

        Occupies bytes 185 - 196 (numbered from one) in the header.
        """
        return self._array['tilt'][0]
    def _settilt(self, value):
        self._array['tilt'][0] = value
    tilt = property(_gettilt, _settilt)


    def _getNumWaves(self):
        """Is an integer for the number of wavelengths stored in the file.
        The wavelength dimension is handled differently than z or time
        in that other header values store per-wavelength metadata.  That
        storage allows for information about five or less wavelengths.
        Because of that limitation and because many Priism applications
        assume a maximum of five wavelengths, you would normally restrict
        the number of wavelengths stored to be five or less.

        Occupies bytes 197 - 198 (numbered from one) in the header.
        """
        return self._array['NumWaves'][0]
    def _setNumWaves(self, value):
        self._array['NumWaves'][0] = value
    NumWaves = property(_getNumWaves, _setNumWaves)


    def _getwave(self):
        """Is five 16-bit integers as a NumPy array.  They are the
        wavelengths for the emitted or transmitted light in
        nanometers rounded to the nearest integer.  For a broad
        passband, you would normally use some measure of the center
        of the passband as the wavelength to store in the header.

        Occupies bytes 199 - 208 (numbered from one) in the header.
        """
        return self._array['wave'][0]
    def _setwave(self, value):
        self._array['wave'][0] = value
    wave = property(_getwave, _setwave)


    def _getzxy0(self):
        """Is three floating-point values as a NumPy array.  The three
        values specify an origin to use in coordinate transformations.
        The first value is the z coordinate, the second value is the
        x coordinate, and the third value is the y coordinate.  The
        units for each are, by convention, the same as used for the
        spacing between samples:  Angstroms for EM data and microns
        for optical data.

        Occupies bytes 209 - 220 (numbered from one) in the header.
        """
        return self._array['zxy0'][0]
    def _setzxy0(self, value):
        self._array['zxy0'][0] = value
    zxy0 = property(_getzxy0, _setzxy0)


    def _getrms(self):
        """Is part of the MRC 2014 format but is not in the Priism format.
        Accesses will return -1.  Attempts to set this to something other
        than -1 will generate an AttributeError exception.
        """
        return -1.0
    def _setrms(self, value):
        if value != -1.0:
            raise AttributeError('Priism MRC header does not store RMS')
    rms = property(_getrms, _setrms)


    def _getorigin(self):
        """Is three floating-point values that are part of the MRC 2014
        format.  The Priism format zxy0 field has the same purpose but
        is in a different part of the header and has a different ordering
        for the coordinates.  Reorder and return those values when accessed.
        When set, reorder the values and set the fields used by the Priism
        format.
        """
        return ReorderedArray(self._array['zxy0'][0], (1, 2, 0))
    def _setorigin(self, value):
        self._array['zxy0'][0] = (value[2], value[0], value[1])
    origin = property(_getorigin, _setorigin,
)


class Hdr2014(HdrBase):
    """Represents a MRC 2014 header without extensions.

    Only provides access to the parts of the header that are in the
    MRC 2014 specification.  The fields that are part of the original
    MRC formulation are described in HdrBase.  The additions for
    MRC 2014 are listed below.

    This class was added in version 0.1.0 of Mrc.py as included with
    Priism.
    """

    # Use __slots__ so any assignments to unspecified attributes or
    # properties will give an AttributeError exception.
    __slots__ = ()


    def index2zwt(self, i, over=False):
        """Convert section index to 3D index in z, wavelength, and time.

        Positional parameters:
        i -- Is the section index.  If less than zero, it is treated as a
        displacement from one past the end (i.e. -1 is the last valid
        section index, and -2 is the next-to-last valid section index).

        Keyword parameters:
        over -- If True, a value of i past the end will cause a ValueError
        exception to be raised.  If False, a value of i past the end will
        lead to a value for the index in the slowest-varying dimension to
        be past the end.

        Return value:
        Returns a three element tuple of the zero-based indices for z,
        wavelength, and time.
        """

        if self.nspg == 401:
            if self.m[2] > 0:
                nz = self.m[2]
                nt = max(1, self.Num[2] // nz)
            else:
                nz = 1
                nt = max(self.Num[2], 1)
        else:
            nz = max(self.Num[2], 1)
            nt = 1
        nw = 1
        # All the available image sequence types have z varying faster
        # than time; use the one with wavelength varying fastest since
        # any overflow will not go into the wavelength dimension.
        return index2zwt(i, nz, nw, nt, 1, over=over)


    def zwt2index(self, z, w, t, over=False):
        """Convert a 3D index in z, wavelength, and time to a section index.

        Positional parameters:
        z -- Is the index in z.  If less than zero, it is treated as a
        displacement from one past the end i.e. -1 is the last valid
        section index, and -2 is the next-to-last valid section index).

        w -- Is the wavelength index.  If less than zero, it is treated as
        a displacement from one past the end.

        t -- Is the time index.  If less than zero, it is treated as a
        displacement from one past the end.

        Keyword parameters:
        over -- If True, a ValueError exception will be raised if any of
        z, w, or t are past the end in their respective dimesions.  If
        False, a value for z, w, or t that is past the end will only
        result in a ValueError exception if that dimension is not the
        slowest-varying non-singleton dimension.

        Return value:
        Returns a zero-based index for the section.
        """

        if self.nspg == 401:
            if self.m[2] > 0:
                nz = self.m[2]
                nt = max(1, self.Num[2] // nz)
            else:
                nz = 1
                nt = max(self.Num[2], 1)
        else:
            nz = max(self.Num[2], 1)
            nt = 1
        nw = 1
        return zwt2index(z, w, t, nz, nw, nt, 1, over=over)


    def _getexttyp(self):
        """Is four bytes treated a four character string that is a code
        for the layout of the extended header.  This implementation
        understands three different values, all encoded in ASCII, for
        this field:  'MRC0', 'CCP4' (treated as a synonym for 'MRC0'),
        and 'AGAR'.

        This field was introduced by the MRC 2014 standard and was not
        part of the earlier Image 2000 standard.

        Occupies bytes 105 - 108 (numbered from one) in the header.
        """
        return self._array['exttyp'][0]
    def _setexttyp(self, value):
        self._array['exttyp'][0] = value
    exttyp = property(_getexttyp, _setexttyp)


    def _getnversion(self):
        """Is a 4-byte integer which stores the version number of the
        MRC format used by this file.  The version number is 10 times
        the Gregorian year when the specification was issued plus a
        zero-based (up to 9) version number within a year.

        Files which use the MRC 2014 format would have 20140 in this
        field.

        This field was introduced by the MRC 2014 standard and was not
        part of the earlier Image 2000 standard.

        Occupies bytes 109 - 112 (numbered from one) in the header.
        """
        return self._array['nversion'][0]
    def _setnversion(self, value):
        self._array['nversion'][0] = value
    nversion = property(_getnversion, _setnversion)


    def _getorigin(self):
        """Is three floating-point values as a NumPy array.  The three
        values specify an origin to use in coordinate transformations.
        The first value is the x coordinate, the second value is the
        y coordinate, and the third value is the z coordinate.  The
        units for each are, by convention, the same as used for the
        spacing between samples:  Angstroms for EM data and microns
        for optical data.

        Occupies bytes 197 - 208 (numbered from one) in the header.
        """
        return self._array['origin'][0]
    def _setorigin(self, value):
        self._array['origin'][0] = value
    origin = property(_getorigin, _setorigin)


    def _getmap(self):
        """Is four bytes treated as a four character string that
        identifies this file as an MRC file.  The MRC 2014 and
        Image 2000 standards specify that this field should be
        set to 'MAP ' encoded in ASCII.

        Occupies bytes 209 - 212 (numbered from one) in the header.
        """
        return self._array['map'][0]
    def _setmap(self, value):
        self._array['map'][0] = value
    map = property(_getmap, _setmap)


    def _getmachst(self):
        """Is four bytes that are a code for how floating-point, complex,
        integer, and character values are stored.  In practice, two
        combinations of values are used:  0x44 and 0x41 in the first two
        bytes (68 and 65 in decimal; the MRC 2014 documentation has 0x44
        in the second byte as well) and unspecified values (typically zero)
        in the second two bytes for little-endian and 0x11 (17 in decimal)
        in the first two bytes and unspecified values (typically zero) in
        the second two bytes for big-endian.

        Occupies bytes 213 - 216 (numbered from one) in the header.
        """
        return self._array['machst'][0]
    def _setmachst(self, value):
        self._array['machst'][0] = value
    machst = property(_getmachst, _setmachst)


    def _getrms(self):
        """Is a floating-point value for the RMS deviation of the densities
        from the mean density.  If the RMS deviation has not been computed,
        the convention is to put a value less than zero in this field.

        Occupies bytes 217 - 220 (numbered from one) in the header.
        """
        return self._array['rms'][0]
    def _setrms(self, value):
        self._array['rms'][0] = value
    rms = property(_getrms, _setrms)


class Hdr2014Priism(Hdr2014):
    """Represents a MRC 2014 header with fields from Priism where they
    do not conflict with the MRC 2014 standard.

    Provides access to all parts of the header.  The fields that are
    part of the original MRC formulation are described in HdrBase.  The
    fields specific to MRC 2014 are described in Hdr2014.  The extensions
    from Priism are listed below.

    Some properties are provided to make it easier to use instances of
    this class interchangeably with instances of HdrPriism, especially
    for retrieving values without modifying them.

    This class was added in version 0.1.0 of Mrc.py as included with
    Priism.
    """

    # Use __slots__ so any assignments to unspecified attributes or
    # properties will give an AttributeError exception.
    __slots__ = ()


    def _getextra0(self):
        """Is 4 bytes that is not currently reserved for any specific
        purpose.  Priism and this implementation do not perform any
        byte-swapping on this data.  If you store multi-byte quantities
        here, you will need to handle the byte-swapping.

        Occupies bytes 97 - 100 (numbered from one) in the header.
        """
        return self._array['extra0'][0]
    def _setextra0(self, value):
        self._array['extra0'][0] = value
    extra0 = property(_getextra0, _setextra0)


    def _getntst(self):
        """Is an integer used to store the starting time index, i.e.
        the index in a longer time series for the first time point in
        this dataset.

        Occupies bytes 101 - 104 (numbered from one) in the header.
        """
        return self._array['ntst'][0]
    def _setntst(self, value):
        self._array['ntst'][0] = value
    ntst = property(_getntst, _setntst)


    def _getextra1(self):
        """Is 16 bytes that is not currently reserved for any specific
        purpose.  Priism and this implementation do not perform any
        byte-swapping on this data.  If you store multi-byte quantities
        here, you will need to handle the byte-swapping.

        Occupies bytes 113 - 128 (numbered from one) in the header.
        """
        return self._array['extra1'][0]
    def _setextra1(self, value):
        self._array['extra1'][0] = value
    extra1 = property(_getextra1, _setextra1)


    def _getNumIntegers(self):
        """Is the number of 4-byte integers to store per section in
        the extended header.

        Occupies bytes 129 - 130 (numbered from one) in the header.
        """
        return self._array['NumIntegers'][0]
    def _setNumIntegers(self, value):
        self._array['NumIntegers'][0] = value
    NumIntegers = property(_getNumIntegers, _setNumIntegers)


    def _getNumFloats(self):
        """Is the number of 4-byte floating-point values to store
        per section in the extended header.

        Occupies bytes 131 - 132 (numbered from one) in the header.
        """
        return self._array['NumFloats'][0]
    def _setNumFloats(self, value):
        self._array['NumFloats'][0] = value
    NumFloats = property(_getNumFloats, _setNumFloats)


    def _getsub(self):
        """Is the number of different resolutions stored in the
        dataset.

        Occupies bytes 133 - 134 (numbered from one) in the header.
        """
        return self._array['sub'][0]
    def _setsub(self, value):
        self._array['sub'][0] = value
    sub = property(_getsub, _setsub)


    def _getzfac(self):
        """For multiple resolutions, the number of z samples in a
        resolution will be the number of z samples in the next higher
        resolution divided by this integer and rounding up any remainder.

        Occupies bytes 135 - 136 (numbered from one) in the header.
        """
        return self._array['zfac'][0]
    def _setzfac(self, value):
        self._array['zfac'][0] = value
    zfac = property(_getzfac, _setzfac)


    def _getextra2(self):
        """Is 24 bytes that is not currently reserved for any specific
        purpose.  Priism and this implementation do not perform any
        byte-swapping on this data.  If you store multi-byte quantities
        here, you will need to handle the byte-swapping.

        Occupies bytes 137 - 160 (numbered from one) in the header.
        """
        return self._array['extra2'][0]
    def _setextra2(self, value):
        self._array['extra2'][0] = value
    extra2 = property(_getextra2, _setextra2)


    def _getImageType(self):
        """Is a 16-bit integer which is a code for the type of data in
        the file.  http://msg.ucsf.edu/IVE/IVE4_HTML/IM_ref2.html#ImageTypes
        describes the types that Priism defines and how the n1, n2, v1,
        and v2 fields are used for each type.

        Occupies bytes 161 - 162 (numbered from one) in the header.
        """
        return self._array['ImageType'][0]
    def _setImageType(self, value):
        self._array['ImageType'][0] = value
    ImageType = property(_getImageType, _setImageType)


    def _getLensNum(self):
        """For optical data, this 16-bit integer is used to indicate
        which microscope configuration was used when the data was
        collected.  Each microscope system would have its own table
        of configurations and corresponding integer codes.

        Occupies bytes 163 - 164 (numbered from one) in the header.
        """
        return self._array['LensNum'][0]
    def _setLensNum(self, value):
        self._array['LensNum'][0] = value
    LensNum = property(_getLensNum, _setLensNum)


    def _getn1(self):
        """Is a 16-bit integer whose interpretation depends on the value
        for ImageType.

        Occupies bytes 165 - 166 (numbered from one) in the header.
        """
        return self._array['n1'][0]
    def _setn1(self, value):
        self._array['n1'][0] = value
    n1 = property(_getn1, _setn1)


    def _getn2(self):
        """Is a 16-bit integer whose interpretation depends on the value
        for ImageType.

        Occupies bytes 167 - 168 (numbered from one) in the header.
        """
        return self._array['n2'][0]
    def _setn2(self, value):
        self._array['n2'][0] = value
    n2 = property(_getn2, _setn2)


    def _getv1(self):
        """Is a 16-bit integer which is used to store a floating point value,
        f, as round_to_nearest(f * 100.0).  What f is depends on the value
        for ImageType.  This Python implementation leaves the conversion
        between f and v1 to the caller.

        Occupies bytes 169 - 170 (numbered from one) in the header.
        """
        return self._array['v1'][0]
    def _setv1(self, value):
        self._array['v1'][0] = value
    v1 = property(_getv1, _setv1)


    def _getv2(self):
        """Is a 16-bit integer which is used to store a floating point value,
        f, as round_to_nearest(f * 100.0).  What f is depends on the value
        for ImageType.  This Python implementation leaves the conversion
        between f and v2 to the caller.

        Occupies bytes 171 - 172 (numbered from one) in the header.
        """
        return self._array['v2'][0]
    def _setv2(self, value):
        self._array['v2'][0] = value
    v2 = property(_getv2, _setv2)


    def _getextra3(self):
        """Is 12 bytes that is not currently reserved for any specific
        purpose.  Priism and this implementation do not perform any
        byte-swapping on this data.  If you store multi-byte quantities
        here, you will need to handle the byte-swapping.

        Occupies bytes 173 - 184 (numbered from one) in the header.
        """
        return self._array['extra3'][0]
    def _setextra3(self, value):
        self._array['extra3'][0] = value
    extra3 = property(_getextra3, _setextra3)


    def _gettilt(self):
        """Is three floating-point values as a NumPy array.  The three
        values are a trio of rotation angles in degrees.  For image
        windows, Priism uses those angles to rotate the coordinates of
        overlayed objects to the coordinates aligned to the image axes.
        The transformation uses the first angle as rotation about the
        original +x axis of the objects with a positive angle rotating
        the +y axis towards the +z axis.  The second angle rotates about
        the +y axis from the first rotation with a positive angle
        rotating the +z axis towards the +x axis.  The third angle
        rotates about +z axis from the second rotation with a positive
        angle rotating the +x axis towards the +y axis.

        Occupies bytes 185 - 196 (numbered from one) in the header.
        """
        return self._array['tilt'][0]
    def _settilt(self, value):
        self._array['tilt'][0] = value
    tilt = property(_gettilt, _settilt)


    def _getmm2(self):
        """Is not part of the MRC 2014 format with Priism extensions.
        Is part of the Priism format.  On access, will return a
        NumPy array with two zero values.  Any attempt to set to
        something which is not two zero values will generate an
        AttributeError exception.
        """
        return N.zeros(2, dtype='f4')
    def _setmm2(self, value):
        if value[0] != 0.0 or value[1] != 0.0:
            raise AttributeError('MRC 2014 header does not store minimum '
                                 'and maximum values for wavelengths past '
                                 'the first')
    mm2 = property(_getmm2, _setmm2)


    def _getmm3(self):
        """Is not part of the MRC 2014 format with Priism extensions.
        Is part of the Priism format.  On access, will return a
        NumPy array with two zero values.  Any attempt to set to
        something which is not two zero values will generate an
        AttributeError exception.
        """
        return N.zeros(2, dtype='f4')
    def _setmm3(self, value):
        if value[0] != 0.0 or value[1] != 0.0:
            raise AttributeError('MRC 2014 header does not store minimum '
                                 'and maximum values for wavelengths past '
                                 'the first')
    mm3 = property(_getmm3, _setmm3)


    def _getmm4(self):
        """Is not part of the MRC 2014 format with Priism extensions.
        Is part of the Priism format.  On access, will return a
        NumPy array with two zero values.  Any attempt to set to
        something which is not two zero values will generate an
        AttributeError exception.
        """
        return N.zeros(2, dtype='f4')
    def _setmm4(self, value):
        if value[0] != 0.0 or value[1] != 0.0:
            raise AttributeError('MRC 2014 header does not store minimum '
                                 'and maximum values for wavelengths past '
                                 'the first')
    mm4 = property(_getmm4, _setmm4)


    def _getmm5(self):
        """Is not part of the MRC 2014 format with Priism extensions.
        Is part of the Priism format.  On access, will return a
        NumPy array with two zero values.  Any attempt to set to
        something which is not two zero values will generate an
        AttributeError exception.
        """
        return N.zeros(2, dtype='f4')
    def _setmm5(self, value):
        if value[0] != 0.0 or value[1] != 0.0:
            raise AttributeError('MRC 2014 header does not store minimum '
                                 'and maximum values for wavelengths past '
                                 'the first')
    mm5 = property(_getmm5, _setmm5)


    def _getNumTimes(self):
        """The MRC 2014 format store the number of volumes implicitly
        and uses a space group of 401 to indicate that multiple volumes
        may be present.  Handle Priism's notion of number of time points
        with the MRC 2014 mechanism.
        """
        if self._array['nspg'][0] == 401:
            # The file represents a volume stack.  The third component
            # of m is the number of z slices per volume.
            if self._array['m'][0][2] > 0:
                return self._array['Num'][0][2] // self._array['m'][0][2]
            # The number of z slices per volume is invalid, assume one
            # slice per volume.
            return self._array['Num'][0][2]
        return 1
    def _setNumTimes(self, value):
        if value > 0:
            if self._array['nspg'][0] == 401:
                oldnz = self._array['m'][0][2]
            else:
                if value > 1:
                    self._array['nspg'][0] = 401
                oldnz = self._array['Num'][0][2]
            newnz = self._array['Num'][0][2] // value
            # Preserve pixel spacing.
            if oldnz != 0:
                self._array['d'][0][2] *= newnz / float(oldnz)
            else:
                self._array['d'][0][2] = 0.0
            self._array['m'][0][2] = newnz
        else:
            raise ValueError('Attempts to set the number of time points to '
                             'a nonpositive value')
    NumTimes = property(_getNumTimes, _setNumTimes)


    def _getImgSequence(self):
        """Is not part of the MRC 2014 format with Priism extensions.  All
        MRC 2014 files implicitly have z varying faster than time.  Return
        zero on accesses.  Any attempts to set the value to something
        other than zero will generate an AttributeError exception.
        """
        # This is Priism's code for the ZTW ordering.  All of Priism's
        # defined orderings are the equivalent for this purpose, though,
        # since the number of wavelengths is one and z varies faster than
        # time.
        return 0
    def _setImgSequence(self, value):
        if value != 0:
            raise AttributeError('The MRC 2014 header does not store the '
                                 'image sequence code')
    ImgSequence = property(_getImgSequence, _setImgSequence)


    def _getNumWaves(self):
        """Is not part of the MRC 2014 format with Priism extensions.
        All such files implicitly have a single wavelength.  Return
        one for all accesses.  Any attempts to set the value to something
        other than one will generate an AttributeError exception.
        """
        return 1
    def _setNumWaves(self, value):
        if value != 1:
            raise AttributeError('MRC 2014 format does not store the '
                                 'number of wavelengths')
    NumWaves = property(_getNumWaves, _setNumWaves)


    def _getwave(self):
        """The MRC 2014 format with Priism extensions does not store
        wavelength values.  Return a NumPy array with five zeros for
        any access.  Any attempt to set the wavelength values to
        something other than zeros will generate an AttributeError
        exception.
        """
        return N.zeros(5, dtype='i2')
    def _setwave(self, value):
        if (value[0] != 0 or value[1] != 0 or value[2] != 0 or
            value[3] != 0 or value[4] != 0):
            raise AttributeError('MRC 2014 header does not store '
                                 'wavelength values')
    wave = property(_getwave, _setwave)


    def _getzxy0(self):
        """Is three floating-point values that are part of the Priism
        format.  The MRC 2014 format origin field has the same purpose but
        is in a different part of the header and has a different ordering
        for the coordinates.  Reorder and return those values when accessed.
        When set, reorder the values and set the fields used by the MRC
        2014 format.
        """
        return ReorderedArray(self._array['origin'][0], (2, 0, 1))
    def _setzxy0(self, value):
        self._array['origin'][0] = (value[1], value[2], value[0])
    zxy0 = property(_getzxy0, _setzxy0)


###########################################################################
###########################################################################
###########################################################################
###########################################################################

class ManagedTitleArray(object):
    """Since strings in the header are represented by NumPy arrays of
    signed one byte integers, provide a more convenient interface so that
    a caller can work with Python strings.  This class represents an array
    of strings from the header or extended header.

    This class was added in version 0.1.0 of Mrc.py as included with
    Priism.
    """

    # Use __slots__ so any assignments to unspecified attributes or
    # properties will give an AttributeError exception.
    __slots__ = ('_array',)


    def __init__(self, array):
        """Initializes a ManagedTitleArray.

        Positional parameters:
        array -- Is a NumPy two-dimensional array of one-byte integers.
        The first dimension is the number of strings; the second
        dimension represents the characters in the string.
        """

        self._array = array


    def __repr__(self):
        s = 'ManagedTitleArray(' + repr(self._array) + ')'
        return s


    def __str__(self):
        s = '['
        for i in range(0, self._array.shape[0] - 1):
            s = s + "'" + self._array[i].tostring() + "' ,"
        if self._array.shape[0] > 0:
            s = s + "'" + self._array[-1].tostring() + "']"
        else:
            s = s + ']'
        return s


    def __len__(self):
        return self._array.shape[0]


    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self._array[i].tostring() for i in
                    range(*key.indices(self._array.shape[0]))]
        else:
            return self._array[key].tostring()


    def __setitem__(self, key, value):
        if isinstance(key, slice):
            if isStringLike(value):
                for i in range(*key.indices(self._array.shape[0])):
                    self[i] = value
            else:
                r = range(*key.indices(self._array.shape[0]))
                if len(r) != len(value):
                    raise ValueError('Number of elements on right hand side does '
                                     'not match number on left')
                j = 0
                for i in r:
                    self[i] = value[j]
                    j += 1
        else:
            b = N.fromstring(value, dtype='i1')
            if b.shape[0] >= self._array.shape[1]:
                self._array[key][:] = b[0:self._array.shape[1]]
            else:
                # Pad with spaces.
                self._array[key][:] = N.concatenate(
                    (b, 32 * N.ones(self._array.shape[1] - b.shape[0],
                                    dtype='i1')))


###########################################################################
###########################################################################
###########################################################################
###########################################################################

class ReorderedArray(object):
    """Provide an interface that acts like a reordered view for the first
    dimension of a NumPy array.

    This class was added in version 0.1.0 of Mrc.py as included with
    Priism.
    """

    # Use __slots__ so any assignments to unspecified attributes or
    # properties will give an AttributeError exception.
    __slots__ = ('_array', '_indices')


    def __init__(self, array, indices):
        """Initializes a ReorderedArray.

        Positional parameters:
        array -- Is a NumPy array with at least one dimension.
        indices -- Is an iterable with the new indices (from first to
        last) to use for first dimension of the array.  The length
        of indices should not exceed the size of the first dimension
        of array.
        """

        self._array = array
        self._indices = indices


    def __repr__(self):
        s = ('ReorderedArray(' + repr(self._array) + ',' +
             repr(self._indices) + ')')
        return s


    def __str__(self):
        s = '['
        for i in range(0, len(self._indices) - 1):
            s = s + str(self._array[self._indices[i]]) + ', '
        if len(self._indices) > 0:
            s = s + str(self._array[self._indices[-1]]) + ']'
        else:
            s = s + ']'
        return s


    def __len__(self):
        return len(self._indices)


    def __getitem__(self, key):
        if isinstance(key, slice):
            return [self._array[self._indices[i]] for i in
                    range(*key.indices(self._array.shape[0]))]
        else:
            return self._array[self._indices[key]]


    def __setitem__(self, key, value):
        if isinstance(key, slice):
            r = range(*key.indices(self._array.shape[0]))
            if hasattr(value, '__len__'):
                if len(r) != len(value):
                    raise ValueError('Number of elements on right hand side '
                                     'does not match number on left')
                j = 0
                for i in r:
                    self[self._indices[i]] = value[j]
                    j += 1
            else:
                for i in r:
                    self[self._indices[i]] = value
        else:
            self._array[self._indices[key]] = value


###########################################################################
###########################################################################
###########################################################################
###########################################################################

def minExtHdrSize(nSecs, bytesPerSec):
    """Return the smallest multiple of 1024 capable of holding the
    extended header data.

    Positional parameters:
    nSecs -- Is the number of sections of extended header data.  It is
    assumed to be a non-negative value.

    bytesPerSec -- Is the number of bytes per section to store.  It is
    assumed to be a non-negative value.
    """

    t = nSecs * bytesPerSec
    r = t % 1024
    if r != 0:
        t += 1024 - r
    return t


def MrcMode2dtype(mode):
    """Return a NumPy dtype equivalent to a MRC pixel type code.

    Raises a RuntimeError exception if the given mode is not known to this
    implementation.  Raises a NotImplementedError is mode is part of the
    Priism or MRC 2014 specifications but is not handled by this
    implementation.

    Version 0.1.0 of Mrc.py as included with Priism changed the return
    value from a Python type to a NumPy dtype.

    Positional parameters:
    mode -- Is the integer code for how sample values are represented, i.e.
    the value of the 'PixelType' field in the header.
    """

    if mode == 0:
        dt = N.dtype('i1')
    elif mode == 1:
        dt = N.dtype('i2')
    elif mode == 2:
        dt = N.dtype('f4')
    elif mode == 3:
        dt = N.dtype({'names':['real', 'imag'], 'formats':['i2', 'i2']})
    elif mode == 4:
        dt = N.dtype('c8')
    elif mode == 5:
        dt = N.dtype('i2')
    elif mode == 6:
        dt = N.dtype('u2')
    elif mode == 7:
        dt = N.dtype('i4')
    elif mode == 101:
        raise NotImplementedError('Mrc.py does not handle the unsigned '
                                  '4-bit data type')
    else:
        raise RuntimeError('Unknown pixel type code, %d' % mode)
    return dt


def dtype2MrcMode(dtype):
    """Return the MRC sample format code number equivalent to a Python
    type or NumPy dtype.

    Positional parameters:
    dtype -- Is the Python type or NumPy dtype used for each sample.

    Return value:
    Returns an integer for the field labeled 'PixelType' (bytes 13 - 16)
    in the MRC header.  If dtype is not equivalent to one of the types
    supported by the MRC format, a ValueError exception will be raised.
    """

    hastype = hasattr(dtype, 'type')
    if dtype == N.int8 or (hastype and dtype.type == N.int8):
        return 0
    if dtype == N.int16 or (hastype and dtype.type == N.int16):
        return 1
    if dtype == N.float32 or (hastype and dtype.type == N.float32):
        return 2
    if dtype == N.complex64 or (hastype and dtype.type == N.complex64):
        return 4
    if dtype == N.uint16 or (hastype and dtype.type == N.uint16):
        return 6
    if dtype == N.int32 or (hastype and dtype.type == N.int32):
        return 7
    if (hasattr(dtype, 'fields') and dtype.fields is not None and
        len(dtype.fields) == 2):
        fields_ok = 0
        for k in dtype.fields:
            if dtype.fields[k][0] == N.int16:
                if dtype.fields[k][1] == 0:
                    fields_ok |= 1
                elif dtype.fields[k][1] == 2:
                    fields_ok |= 2
        if fields_ok == 3:
            return 3
    if hasattr(dtype, 'name'):
        name = dtype.name
    else:
        name = str(dtype)
    raise ValueError('MRC does not support %s' % name)


def shapeFromHdr(hdr, verbose=0):
    """Return a tuple of array dimensions equivalent to the sizes set in a
    MRC header.

    Positional parameters:
    hdr -- Is the MRC header to use.  It is expected to be a header
    as returned by makeHdrArray() or implement_hdr().  If x is an
    instance of the Mrc or Mrc2 classes, you can use x.hdr for this
    parameter.  Non-positive values for the number of wavelengths
    or time points in the header will be treated as if they were
    equal to one.  Negative values in the header for the number
    of x samples, number of y samples, or number of sections will
    be treated as if they were zero.

    Keyword parameters:
    verbose -- If true, a string giving the ordering of the dimensions
    from slowest to fastest will be printed.  Each dimension will
    be represented by a single letter, 'x', 'y', 'z', 'w', or 't',
    and those letters will be separated by commas.

    Return value:
    Returns a tuple of integers which are the size for each dimension,
    from slowest to fastest varying, specified by the MRC header.
    """

    zOrder = hdr.ImgSequence # , 'Image sequence. 0=ZTW, 1=WZT, 2=ZWT. '),
    if zOrder < 0 or zOrder > 2:
        # The value is invalid; use the default, ZTW, instead.
        zOrder = 0
    nt, nw = hdr.NumTimes, hdr.NumWaves
    nx, ny, nsecs = hdr.Num
    if nt <= 0:
        nt=1
    if nw <= 0:
        nw=1
    if nx < 0:
        nx = 0
    if ny < 0:
        ny = 0
    if nsecs < 0:
        nsecs = 0
    nz = nsecs // (nt * nw)

    if nt == nw == 1:
        shape = (nz, ny, nx)
        orderLetters = 'zyx'
    elif nz == 1 == nw:
        shape = (nt, ny, nx)
        orderLetters = 'tyx'
    elif nt == 1 or nw == 1:
        if zOrder == 0 or zOrder == 2:
            nn = nt
            if nt == 1:
                nn = nw
                orderLetters = 'wyx'
            else:
                orderLetters = 'tyx'
            shape = (nn, nz, ny, nx)
        else: # if zOrder == 1:
            if nt == 1:
                shape = (nz, nw, ny, nx)
                orderLetters = 'zwyx'
            else:
                shape = (nt, nz, ny, nx)
                orderLetters = 'tzyx'

    else: # both nt and nw > 1
        if zOrder == 0:
            shape = (nw, nt, nz, ny, nx)
            orderLetters = 'wtzyx'
        elif zOrder == 1:
            shape = (nt, nz, nw, ny, nx)
            orderLetters = 'tzwyx'
        else: # zOrder == 2:
            shape = (nt, nw, nz, ny, nx)
            orderLetters = 'twzyx'

    if verbose:
        print(','.join(orderLetters))
    return shape


def implement_hdr(hdrArray, hasMap=False):
    """Return a HdrPriism or Hdr2014Priism instance to wrap the given
    NumPy structured array.

    If h is an object returned by this function, it can be used in
    statements like

    h.d = (1,2,3)

    or

    h.LensNum = 13

    to modify values in the header or in statements like

    shape = (h.Num[2], h.Num[1], h.Num[0])

    or

    mean = h.mmm1[2]

    to retrieve values from the header.  The documentation for HdrBase
    and HdrPriism describes the fields available in a HdrPriism
    instance.  The documentation for HdrBase, Hdr2014, and Hdr2014Priism
    describes the fields available in a Hdr2014Priism instance.

    The return value also has convenience methods for the spacing
    between samples and for the titles in the header.  Those are
    described in the documentation for HdrBase.

    To get the original NumPy structured array from the return value,
    use

    h._array

    Positional parameters:
    hdrArray -- Is a MRC header as a NumPy structured array with
    one element.  The dtype of the array should be
    numpy.dtype(mrcHdr_dtype) or numpy.dtype(mrc2014Hdr_dtype).

    Keyword parameters:
    hasMap --- If True, return an instance of the Hdr2014Priism class.
    If False, return an instance of HdrPriism class.  The hasMap
    keyword was added in version 0.1.0 of Mrc.py as included with
    Priism.

    Return value:
    If hasMap is False, the return value is an instance of the HdrPriism
    class.  If hasMap is True, the return value is an instance of
    the Hdr2014Priism class.  Version 0.1.0 of Mrc.py as included with
    Priism changed the return value from an instance of a class
    defined within implement_hdr() to an instance of HdrPriism or
    Hdr2014Priism.
    """

    if hasMap:
        return Hdr2014Priism(hdrArray)
    return HdrPriism(hdrArray)


def makeHdrArray(buffer=None, makeWeak=True):
    """Create a NumPy structured array for the header and wrap it for easy
    access to the header fields.

    Keyword parameters:
    buffer -- If not None, the structured array will be overlayed on top
    of the first 1024 bytes of buffer.  Typically, buffer would be a NumPy
    array of unsigned 8-bit integers read from the start of a MRC file.

    makeWeak -- Only has an effect if buffer is not None.  In that case,
    the returned object will be marked as a weak reference when makeWeak
    is True.  The makeWeak keyword was added in version 0.1.0 of
    Mrc.py as included with Priism.

    Return value:
    Returns an object representing the header.  If the header was not
    generated from an existing buffer, the header fields have not been
    initialized.  The documentation for the return value of implement_hdr()
    describes how the returned object may be used.  The NumPy structured
    array embedded as the _array attribute of the returned object either
    has a NumPy dtype of numpy.dtype(mrcHdr_dtype) or
    numpy.dtype(mrc2014Hdr_dtype).
    """

    if buffer is not None:
        h=buffer
        if hdrHasMapField(buffer):
            h.dtype = mrc2014Hdr_dtype
            hasmap = True
        else:
            h.dtype = mrcHdr_dtype
            hasmap = False
        if makeWeak:
            h = weakref.proxy(h)
    else:
        h = N.recarray(1, mrcHdr_dtype)
        hasmap = False
    return implement_hdr(h, hasMap=hasmap)


def hdrInfo(hdr):
    """Print a subset of information from a MRC header.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.
    """

    shape = hdr.Num[::-1]
    nz = shape[0]
    numInts = max(0, hdr.NumIntegers)
    numFloats = max(0, hdr.NumFloats)

    print('width:                       %s' % str(shape[2]))
    print('height:                      %s' % str(shape[1]))
    print('# total slices:              %s' % str(shape[0]))

    nt, nw = hdr.NumTimes, hdr.NumWaves

    if nt <= 0  or nw <= 0:
        print(' ** ERROR ** : NumTimes or NumWaves is not positive')
        print('NumTimes: %s' % str(nt))
        print('NumWaves: %s' % str(nw))
    else:
        if nt == 1  and  nw == 1:
            print()
        elif nw == 1:
            print('  (%d times for %d zsecs)' % (nt, nz//nt))
        elif nt == 1:
            print('  (%d waves in %d zsecs)' % (nw, nz//nw))
        else:
            print('  (%d times for %d waves in %d zsecs)'% (nt,
                                                           nw,
                                                           nz // (nw * nt)))

    if nt != 1  or  nw != 1:
        print('# slice order:        %d (0,1,2 = (ZTW or WZT or ZWT)' %
              hdr.ImgSequence)

    d = hdr.getSpacing()
    print('pixel width x    (um):       %s' % str(d[0]))
    print('pixel width y    (um):       %s' % str(d[1]))
    print('pixel height     (um):       %s' % str(d[2]))

    print('# wavelengths:               %s' % str(nw))
    print('   wavelength 1  (nm):       %s' % str(hdr.wave[0]))
    print('    intensity min/max/mean:  %s %s %s' %
          (str(hdr.mmm1[0]), str(hdr.mmm1[1]), str(hdr.mmm1[2])))
    if nw > 1:
        print('   wavelength 2  (nm):       %s' % str(hdr.wave[1]))
        print('    intensity min/max:       %s %s' %
              (str(hdr.mm2[0]), str(hdr.mm2[1])))
    if nw > 2:
        print('   wavelength 3  (nm):       %s' % str(hdr.wave[2]))
        print('    intensity min/max:       %s %s' %
              (str(hdr.mm3[0]), str(hdr.mm3[1])))
    if nw > 3:
        print('   wavelength 4  (nm):       %s' % str(hdr.wave[3]))
        print('    intensity min/max:       %s %s' %
              (str(hdr.mm4[0]), str(hdr.mm4[1])))
    if nw > 4:
        print('   wavelength 5  (nm):       %s' % str(hdr.wave[4]))
        print('    intensity min/max:       %s %s' %
              (str(hdr.mm5[0]), str(hdr.mm5[1])))

    if hdr.LensNum == 12:
        name = ' (60x)'
    elif hdr.LensNum == 13:
        name = ' (100x)'
    else:
        name = '(??)'
    print('lens type:                   %s %s' % (str(hdr.LensNum), name))

    print('origin   (um) x/y/z:         %s %s %s' %
          (str(hdr.zxy0[1]), str(hdr.zxy0[2]), str(hdr.zxy0[0])))

    if hdr.PixelType == 0:
        name = '8 bit (signed)'
    elif hdr.PixelType == 1:
        name = '16 bit (signed)'
    elif hdr.PixelType == 2:
        name = '32 bit (signed real)'
    elif hdr.PixelType == 3:
        name = '16 bit (signed complex integer)'
    elif hdr.PixelType == 4:
        name = '32 bit (signed complex real)'
    elif hdr.PixelType == 5:
        name = '16 bit (signed) IW_EMTOM'
    elif hdr.PixelType == 6:
        name = '16 bit (unsigned short)'
    elif hdr.PixelType == 7:
        name = '32 bit (signed long)'
    elif hdr.PixelType == 101:
        name = 'unsigned 4-bit'
    else:
        name = ' ** undefined ** '
    print('# pixel data type:             %s' % name)

    if hdr.next > 0:
        n = numInts + numFloats
        if n > 0:
            name = ' (%d secs)' % (hdr.next // (4. * n))
        else:
            name = ' (??? secs)'
        name2 = '  (%d ints + %d reals per section)' % (numInts, numFloats)
    else:
        name = ''
        name2 = None
    print('# extended header size:        %s %s' % (str(hdr.next), name))
    if name2 is not None:
        print(name2)

    if hdr.NumTitles < 0:
        print(' ** ERROR ** : NumTitles less than zero (NumTitles = %s )' %
              str(hdr.NumTitles))
    elif hdr.NumTitles > 0:
        n = hdr.NumTitles
        if n > 10:
            print(' ** ERROR ** : NumTitles larger than 10 (NumTitles = %s )' %
                  hdr.NumTitles)
            n=10
        for i in range(n):
            print('title %d: %s'%(i, hdr.title[i]))


def axisOrderStr(hdr, onlyLetters=True):
    """Return a string indicating the ordering of dimensions.

    x, y, z, w, and t will appear at most once in the string, and at least
    three of them will be present.  The letters that do appear will be
    present in order from slowest varying to fastest varying.  The values
    for the axis field in the header do not affect the result.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.

    Keyword parameters:
    onlyLetters -- If True, only the letters for the dimensions will appear
    in the string.  If False, the first character of the string will be '[',
    the last character of the string will be ']', and a letter for a
    dimension will be preceded by a comma if it is not the first, slowest-
    varying dimension.
    """

    # 'Image sequence. 0=ZTW, 1=WZT, 2=ZWT.'  Given from fastest to slowest.
    zOrder = hdr.ImgSequence
    if zOrder < 0 or zOrder > 2:
        # The value is invalid; use the ZTW instead.
        zOrder = 0
    nt, nw = max(1, hdr.NumTimes), max(1, hdr.NumWaves)
    if nt == nw == 1:
        orderLetters= 'zyx'
    elif nt == 1:
        orderLetters= ('wzyx', 'zwyx', 'wzyx')[zOrder]
    elif nw == 1:
        orderLetters= ('tzyx', 'tzyx', 'tzyx')[zOrder]
    else:
        orderLetters= ('wtzyx', 'tzwyx', 'twzyx')[zOrder]

    if onlyLetters:
        return orderLetters
    else:
        return '[' + ','.join(orderLetters) + ']'


def index2zwt(i, nz, nw, nt, seq, over=False):
    """Convert section index to 3D index in z, wavelength, and time.

    This function was added in version 0.1.0 of Mrc.py as included
    with Priism.

    Positional parameters:
    i -- Is the section index.  If less than zero, it is treated as a
    displacement from one past the end (i.e. -1 is the last valid
    section index, and -2 is the next-to-last valid section index).

    nz -- Is the number of samples in z.  Assumed to be positive.

    nw -- Is the number of wavelengths.  Assumed to be positive.

    nt -- Is the number of time points.  Assumed to be positive.

    seq -- Is the code for the interleaving of z, wavelength, and time.
    If seq is 0, z varies fastest, followed by time, then followed by
    wavelength.  If seq is 1, wavelength variest fastest, followed by z,
    then followed by time.  If seq is 2, z varies fastest, followed by
    wavelength, then followed by time.  For any other value, a
    ValueError exception will be raised.

    Keyword parameters:
    over -- If True, a value of i past the end will cause a ValueError
    exception to be raised.  If False, a value of i past the end will
    lead to a value for the index in the slowest-varying dimension to
    be past the end.

    Return value:
    Returns a three element tuple of the zero-based indices for z,
    wavelength, and time.
    """

    if i < 0:
        ri = i + nz * nw * nt
        if ri < 0:
            raise ValueError('section index, %d, is before first' % i)
    else:
        ri = i
    if over and ri >= nz * nw * nt:
        raise ValueError('index is greater than or equal to nz * nw * nt')
    if seq == 0:
        result = (ri % nz, ri // (nz * nt), (ri // nz) % nt)
    elif seq == 1:
        result = ((ri // nw) % nz, ri % nw, ri // (nw * nz))
    elif seq == 2:
        result = (ri % nz, (ri // nz) % nw, ri // (nz * nw))
    else:
        raise ValueError('invalid code, %d, for sequence arrangement' %
                         seq)
    return result


def zwt2index(z, w, t, nz, nw, nt, seq, over=False):
    """Convert a 3D index in z, wavelength, and time to a section index.

    This function was added in version 0.1.0 of Mrc.py as included
    with Priism.

    Positional parameters:
    z -- Is the index in z.  If less than zero, it is treated as a
    displacement from one past the end i.e. -1 is the last valid
    section index, and -2 is the next-to-last valid section index).

    w -- Is the wavelength index.  If less than zero, it is treated as
    a displacement from one past the end.

    t -- Is the time index.  If less than zero, it is treated as a
    displacement from one past the end.

    nz -- Is the number of samples in z.  Assumed to be positive.

    nw -- Is the number of wavelengths.  Assumed to be positive.

    nt -- Is the number of time points.  Assumed to be positive.

    Keyword parameters:
    over -- If True, a ValueError exception will be raised if any of
    z, w, or t are past the end in their respective dimesions.  If
    False, a value for z, w, or t that is past the end will only
    result in a ValueError exception if that dimension is not the
    slowest-varying non-singleton dimension.

    Return value:
    Returns a zero-based index for the section.
    """

    if z < 0:
        rz = z + nz
        if rz < 0:
            raise ValueError('%d is before first z index' % z)
    else:
        rz = z
    if w < 0:
        rw = w + nw
        if rw < 0:
            raise ValueError('%d is before first wavelength index' % w)
    else:
        rw = w
    if t < 0:
        rt = t + nt
        if rt < 0:
            raise ValueError('%d is before first time point index' % t)
    else:
        rt = t
    if over and (rz >= nz or rw >= nw or rt >= nt):
        raise ValueError('z, wavelength, or time indices out of bounds')
    if seq == 0:
        if (rz >= nz and (nt > 1 or nw > 1 or t > 0 or w > 0)):
            raise ValueError('z index past end and z is not slowest-varying '
                             'non-singleton dimension')
        if (rt >= nt and (nw > 1 or w > 0)):
            raise ValueError('time point index past end and time is not '
                             'slowest-varying non-singleton dimension')
        result = rz + nz * (rt + nt * rw)
    elif seq == 1:
        if (rw >= nw and (nz > 1 or nt > 1 or z > 0 or t > 0)):
            raise ValueError('wavelength index past end and wavelength is not '
                             'slowest-varying non-singleton dimension')
        if (rz >= nz and (nt > 1 or t > 0)):
            raise ValueError('z index past end and z is not slowest-varying '
                             'non-singleton dimension')
        result = rw + nw * (rz + nz * rt)
    elif seq == 2:
        if (rz >= nz and (nw > 1 or nt > 1 or w > 0 or t > 0)):
            raise ValueError('z index past end and z is not slowest-varying '
                             'non-singleton dimension')
        if (rw >= nw and (nt > 1 or t > 0)):
            raise ValueError('wavelength index past end and wavelength is not '
                             'slowest-varying non-singleton dimension')
        result = rz + nz * (rw + nw * rt)
    else:
        raise ValueError('invalid code, %d, for sequence arrangement' %
                         seq)
    return result


def init_simple(hdr, mode, nxOrShape, ny=None, nz=None, isByteSwapped=None):
    """Initialize a MRC header for a given size and pixel type code.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.

    mode -- Is the MRC pixel type code to use.

    nxOrShape -- If ny or nz is not None, this should be a scalar which
    is the number of samples in x.  If ny and nz are None, it must be an
    iterable.  If it has one element, that will be the number of samples
    in x; ny and nz will be one.  If it has two elements, the last will
    be the number of samples in x, the second to last will be the number
    of samples in y, and nz will be one.  If has more than two elements,
    the last will be the number of samples in x, the second to last will
    be the number of samples in y, and the product of the remaining
    elements will be the number of samples in z.

    Keyword parameters:
    ny -- If not None, this will be the number of samples in y.  If not
    None, nxOrShape should be a scalar.

    nz -- If not None, this will be the number of samples in z.  If not
    None, nxOrShape should be a scalar.

    isByteSwapped -- If None, the stamps for byte order in the header
    will be set to match the byte ordering of hdr.  If True, the stamps
    for byte order in the header will be set to be opposite of the
    native byte ordering.  If False, the stamps for byte order in the
    header will match the native byte ordering.  The isByteSwapped keyword
    was added in version 0.1.0 of Mrc.py as included with Priism.
    """

    if ny is None and nz is None:
        if len(nxOrShape) == 2:
            nz, (ny, nx) = 1, nxOrShape
        elif len(nxOrShape) == 1:
            nz, ny, (nx,) = 1, 1, nxOrShape
        elif len(nxOrShape) == 3:
            nz, ny, nx = nxOrShape
        else:
            ny, nx = nxOrShape[-2:]
            nz = N.prod(nxOrShape[:-2])
    else:
        if ny is None:
            ny = 1
        if nz is None:
            nz = 1
        nx = nxOrShape

    hdr.Num = (nx, ny, nz)
    hdr.PixelType = mode
    hdr.mst = (0, 0, 0)
    hdr.m = (1, 1, 1)
    hdr.d = (1, 1, 1)
    hdr.angle = (90, 90, 90)
    hdr.axis = (1, 2, 3)
    hdr.mmm1 = (0, 100000, 5000)
    hdr.nspg = 0
    hdr.next = 0
    hdr.ntst = 0
    hdr.NumIntegers = 0
    hdr.NumFloats = 0
    hdr.sub = 1
    hdr.zfac = 1
    hdr.ImageType = 0
    hdr.LensNum = 0
    hdr.n1 = 0
    hdr.n2 = 0
    hdr.v1 = 0
    hdr.v2 = 0
    hdr.tilt = (0, 0, 0)
    hdr.zxy0 = (0, 0, 0)
    hdr.NumTitles = 0
    hdr.title = ' ' * 80
    if hdrIsInPriismFormat(hdr):
        hdr.dvid = 0xc0a0
        hdr.nblank = 0
        hdr.blank = 0
        hdr.mm2 = (0, 10000)
        hdr.mm3 = (0, 10000)
        hdr.mm4 = (0, 10000)
        hdr.mm5 = (0, 10000)
        hdr.NumTimes = 1
        # Zero means that z changes fastest, then time, and then wavelength.
        # One means that wavelength changes fastest, then z, and then time.
        # Two means that z changes fastest, then wavelength, and then time.
        hdr.ImgSequence = 0
        hdr.NumWaves = 1
        hdr.wave = (999, 0, 0, 0, 0)
    else:
        hdr.extra0 = 0
        hdr.exttyp = N.fromstring('AGAR', dtype='i1')
        hdr.nversion = 20140
        hdr.extra1 = 0
        hdr.extra2 = 0
        hdr.extra3 = 0
        hdr.map = N.fromstring('MAP ', dtype='i1')
        if isByteSwapped is None:
            isByteSwapped = isSwappedDtype(hdr._array.dtype)
        if isByteSwapped:
            hdr.machst = (_SYSSWAPST0, _SYSSWAPST1, 0, 0)
        else:
            hdr.machst = (_SYSMACHST0, _SYSMACHST1, 0, 0)
        hdr.rms = -1.0


def initHdrArrayFrom(hdrDest, hdrSrc):
    """Copy all fields from hdrSrc to hdrDest except the number of
    samples ('Num' field; bytes 1 - 12), pixel type ('PixelType' field;
    bytes 13 - 16), and the number of bytes in the extended header
    ('next' field; bytes 93 - 96.  Do not change the byte order for
    hdrDest.

    Positional parameters:
    hdrDest -- Is the header to change.  Assumed to be a MRC header
    as returned by makeHdrArray() or implement_hdr().  If x is an
    instance of the Mrc or Mrc2 classes, you can use x.hdr for this
    parameter.
    hdrSrc -- Is the header from which to copy.  Assumed to be a
    MRC header as returned by makeHdrArray() or implement_hdr().
    If hdrIsInPriismFormat(hdrSrc) is different than
    hdrIsInPriismFormat(hdrDest), the format of the destination
    will be changed to match the source.

    Return value:
    Returns the header, wrapped with a HdrPriism or Hdr2014Priism
    instance, as appropriate.  The return value was added in version
    0.1.0 of Mrc.py as included with Priism.

    """

    srcIsPriism = hdrIsInPriismFormat(hdrSrc)
    destIsPriism = hdrIsInPriismFormat(hdrDest)
    if srcIsPriism:
        if not destIsPriism:
            if isSwappedDtype(hdrDest._array.dtype):
                hdrDest._array.dtype = N.dtype(mrcHdr_dtype).newbyteorder()
            else:
                hdrDest._array.dtype = N.dtype(mrcHdr_dtype)
            hdrDest = HdrPriism(hdrDest._array)
        hdrDest.mst = hdrSrc.mst
        hdrDest.m = hdrSrc.m
        hdrDest.d = hdrSrc.d
        hdrDest.angle = hdrSrc.angle
        hdrDest.axis = hdrSrc.axis
        hdrDest.mmm1 = hdrSrc.mmm1
        hdrDest.nspg = hdrSrc.nspg
        hdrDest.next = 0
        hdrDest.dvid = hdrSrc.dvid
        hdrDest.nblank = hdrSrc.nblank
        hdrDest.ntst = hdrSrc.ntst
        hdrDest.blank = hdrSrc.blank
        hdrDest.NumIntegers = 0
        hdrDest.NumFloats = 0
        hdrDest.sub = hdrSrc.sub
        hdrDest.zfac = hdrSrc.zfac
        hdrDest.mm2 = hdrSrc.mm2
        hdrDest.mm3 = hdrSrc.mm3
        hdrDest.mm4 = hdrSrc.mm4
        hdrDest.ImageType = hdrSrc.ImageType
        hdrDest.LensNum = hdrSrc.LensNum
        hdrDest.n1 = hdrSrc.n1
        hdrDest.n2 = hdrSrc.n2
        hdrDest.v1 = hdrSrc.v1
        hdrDest.v2 = hdrSrc.v2
        hdrDest.mm5 = hdrSrc.mm5
        hdrDest.NumTimes = hdrSrc.NumTimes
        hdrDest.ImgSequence = hdrSrc.ImgSequence
        hdrDest.tilt = hdrSrc.tilt
        hdrDest.NumWaves = hdrSrc.NumWaves
        hdrDest.wave = hdrSrc.wave
        hdrDest.zxy0 = hdrSrc.zxy0
        hdrDest.NumTitles = hdrSrc.NumTitles
        hdrDest.title = hdrSrc.title
    else:
        if destIsPriism:
            if isSwappedDtype(hdrDest._array.dtype):
                hdrDest._array.dtype = N.dtype(mrc2014Hdr_dtype).newbyteorder()
            else:
                hdrDest._array.dtype = N.dtype(mrc2014Hdr_dtype)
            hdrDest = Hdr2014Priism(hdrDest._array)
        hdrDest.mst = hdrSrc.mst
        hdrDest.m = hdrSrc.m
        hdrDest.d = hdrSrc.d
        hdrDest.angle = hdrSrc.angle
        hdrDest.axis = hdrSrc.axis
        hdrDest.mmm1 = hdrSrc.mmm1
        hdrDest.nspg = hdrSrc.nspg
        hdrDest.next = 0
        hdrDest.extra0 = hdrSrc.extra0
        hdrDest.ntst = hdrSrc.ntst
        hdrDest.exttyp = hdrSrc.exttyp
        hdrDest.nversion = hdrSrc.nversion
        hdrDest.extra1 = hdrSrc.extra1
        hdrDest.NumIntegers = 0
        hdrDest.NumFloats = 0
        hdrDest.sub = hdrSrc.sub
        hdrDest.zfac = hdrSrc.zfac
        hdrDest.extra2 = hdrSrc.extra2
        hdrDest.ImageType = hdrSrc.ImageType
        hdrDest.LensNum = hdrSrc.LensNum
        hdrDest.n1 = hdrSrc.n1
        hdrDest.n2 = hdrSrc.n2
        hdrDest.v1 = hdrSrc.v1
        hdrDest.v2 = hdrSrc.v2
        hdrDest.extra3 = hdrSrc.extra3
        hdrDest.tilt = hdrSrc.tilt
        hdrDest.origin = hdrSrc.origin
        hdrDest.map = hdrSrc.map
        hdrDest.machst = hdrSrc.machst
        hdrDest.rms = hdrSrc.rms
        hdrDest.NumTitles = hdrSrc.NumTitles
        hdrDest.title = hdrSrc.title
    return hdrDest


def setTitle(hdr, s, i=-1, push=False, truncate=False):
    """Set a title in the MRC header.

    Provided for compatibility with previous versions of Mrc.py.
    In this version, you can use hdr.setTitle().

    Positional parameters:
    hdr -- Is the MRC header to modify.  It is expected to be a
    header returned by makeHdrArray() or implement_hdr().  If x
    is an instance of the Mrc or Mrc2 classes, you can use x.hdr
    for this parameter.

    s -- Is the character string for the title.  If s is longer
    than 80 characters and truncate is False, a ValueError
    exception will be raised.  Since no byte swapping is done
    for the titles in the header, s should be encoded in ASCII
    or another format that does not use multibyte characters.

    Keyword parameters:
    i -- Is the index of the title to set.  If i is less than
    zero, the last title not in use will be set.  If i is less
    than zero and all the titles are in use and push is False
    or i is greater than 9, a ValueError exception will be
    raised.

    push -- If True, i is less than zero, and all titles are
    in use, titles will be pushed down before assigning s
    to the last title.  That will discard the first title and
    title[k] (for k greater than and equal to 0 and less than
    9) will be title[k+1] from before the change.  The push
    keyword was added in version 0.1.0 of Mrc.py as included
    with Priism.

    truncate -- If True, only use the first 80 characters from
    s.  The truncate keyword was added in version 0.1.0 of
    Mrc.py as included with Priism.
    """

    hdr.setTitle(s, i=i, push=push, truncate=truncate)


def hdrChangeToMrc2014Format(hdr):
    """Change the header to have the format from MRC 2014.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.

    Return value:
    Returns the header wrapped with a Hdr2014Priism instance.
    """

    if not hdrIsInPriismFormat(hdr):
        return hdr
    nw = max(hdr.NumWaves, 1)
    nt = max(hdr.NumTimes, 1)
    if nw > 1:
        if hdr.ImgSequence == 1 or (hdr.ImgSequence == 2 and nt > 1):
            print('WARNING: MRC 2014 header does not support multiple '
                  'wavelengths; any image data present will be out of sequence')
        else:
            print('WARNING: MRC 2014 header does not support multiple wavelengths')
    dx, dy, dz = hdr.getSpacing()
    origz, origx, origy = hdr.zxy0
    if isSwappedDtype(hdr._array.dtype):
        hdr._array.dtype = N.dtype(mrc2014Hdr_dtype).newbyteorder()
        isswapped = True
    else:
        hdr._array.dtype = N.dtype(mrc2014Hdr_dtype)
        isswapped = False
    hdr = Hdr2014Priism(hdr._array)
    hdr.map = N.fromstring('MAP ', dtype='i1')
    hdr.nversion = 20140
    if nt > 1:
        hdr.nspg = 401
        hdr.m[2] = max(hdr.Num[2], 0) // nt
        hdr.setSpacing(dx, dy, dz)
    if hdr.nspg == 0 or hdr.nspg == 401:
        hdr.exttyp = N.fromstring('AGAR', dtype='i1')
    else:
        hdr.exttyp = N.fromstring('MRC0', dtype='i1')
    hdr.origin = (origx, origy, origz)
    if isswapped:
        hdr.machst = (_SYSSWAPST0, _SYSSWAPST1, 0, 0)
    else:
        hdr.machst = (_SYSMACHST0, _SYSMACHST1, 0, 0)
    hdr.rms = -1.0
    return hdr


def hdrChangeToPriismFormat(hdr):
    """Change the header to have the format from Priism.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.

    Return value:
    Returns the header wrapped with a HdrPriism instance.

    """

    if hdrIsInPriismFormat(hdr):
        return hdr
    origx, origy, origz = hdr.origin
    spg = hdr.npsg
    if isSwappedDtype(hdr._array.dtype):
        hdr._array.dtype = N.dtype(mrcHdr_dtype).newbyteorder()
    else:
        hdr._array.dtype = N.dtype(mrcHdr_dtype)
    hdr = HdrPriism(hdr._array)
    if spg == 1 or spg == 401:
        hdr.nspg = 0
        if spg == 401:
            hdr.NumTimes = hdr.Num[2] // max(hdr.m[2], 1)
    hdr.dvid = 0xc0a0
    hdr.nblank = 0
    hdr.blank = 0
    hdr.mm2 = (0.0, 0.0)
    hdr.mm3 = (0.0, 0.0)
    hdr.mm4 = (0.0, 0.0)
    hdr.ImgSequence = 0
    hdr.NumWaves = 1
    hdr.mm5 = (0.0, 0.0)
    hdr.wave = (0, 0, 0, 0, 0)
    hdr.zxy0 = (origz, origx, origy)
    return hdr


def hdrHasMapField(hdr):
    """Return True if the header has a properly formatted map field.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is the header as an array of 1024 unsigned 8-bit integers.
    """

    # The map field is bytes 208 through 211 (numbered from zero).  It
    # has 'MAP ' in ASCII when properly formatted.
    return (hdr[208] == 77 and hdr[209] == 65 and hdr[210] == 80 and
            hdr[211] == 32)


def hdrIsByteSwapped(hdr):
    """Return True if the header is swapped from the native byte ordering.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.

    Return value:
    Returns True if the machine stamp (which is insensitive to whether
    hdr.dtype has been byte swapped) or the Num field (which is sensitive
    to whether hdr.dtype has been byte swapped) indicates that the header
    is not in the machine byte order.  Otherwise, returns False.  For
    ease of interpreting the return value, isSwappedDtype(hdr.dtype)
    should be False.
    """

    if hasattr(hdr, 'machst'):
        if (hdr.machst[0] == _SYSMACHST0 and hdr.machst[1] >= _SYSMACHST1LO and
            hdr.machst[1] <= _SYSMACHST1HI):
            return False
        elif (hdr.machst[0] == _SYSSWAPST0 and
              hdr.machst[1] >= _SYSSWAPST1LO and
              hdr.machst[1] <= _SYSSWAPST1HI):
            return True
    # Use the test employed by Priism.  Assumes that the actual number of
    # samples in x and y are both positive and less than 65536.  Under those
    # conditions a byte-swapped header will have both either less than zero or
    # greater than 65535.
    nx = hdr.Num[0]
    ny = hdr.Num[1]
    return (nx < 0 or nx > 65535) and (ny < 0 or ny > 65535)


def hdrIsInPriismFormat(hdr):
    """Return True if the header uses the format from Priism.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.
    """

    return not hasattr(hdr, 'map')


def hdrHasExtType(hdr):
    """Return True if the header claims to handle the extended header
    type field.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.
    """

    has = False
    if hasattr(hdr, 'nversion'):
        # Check for a version number that looks reasonable.
        if (hdr.nversion >= 20140 and hdr.nversion <= 20509 and
            hasattr(hdr, 'exttyp')):
            has = True
    return has


def getExtHeaderFormat(hdr):
    """Return a code for the format of the extended header.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    hdr -- Is a MRC header as returned by makeHdrArray() or
    implement_hdr().  If x is an instance of the Mrc or Mrc2
    classes, you can use x.hdr for this parameter.

    Return value:
    Returns zero if the extended header is to be interpreted
    as symmetry information.  Returns one if the extended
    header is to be interpreted as the Priism format.  Returns
    -1 if the extended header format is not understood by this
    implementation.
    """

    if hdrHasExtType(hdr):
        if (N.array_equal(hdr.exttyp, N.fromstring('MRC0', dtype='i1')) or
            N.array_equal(hdr.exttyp, N.fromstring('CCP4', dtype='i1'))):
            fmt = 0
        elif N.array_equal(hdr.exttyp, N.fromstring('AGAR', dtype='i1')):
            fmt = 1
        else:
            fmt = -1
    elif hasattr(hdr, 'map'):
        if hdr.nspg == 0 or hdr.nspg == 1 or hdr.nspg == 401:
            fmt = 1
        else:
            fmt = 0
    else:
        if hdr.nspg == 0:
            fmt = 1
        else:
            fmt = 0
    return fmt


def isStringLike(value):
    """Return True if value is a bytearray, bytes, str, or unicode.  Otherwise
    return False.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.
    """

    if isinstance(value, str):
        return True
    # bytes is only defined in Python 3 and newer versions of Python 2.
    try:
        if isinstance(value, bytes):
            return True
    except:
        pass
    # unicode is only defined in Python 2.
    try:
        if isinstance(value, unicode):
            return True
    except:
        pass
    # bytearray was introduced in Python 2.6.
    try:
        if isinstance(value, bytearray):
            return True
    except:
        pass
    return False


def isSwappedDtype(dtype):
    """Return True if the given NumPy dtype is not in machine byte order.

    Will raise a ValueError exception if dtype is a structured type and has
    components which have different byte orders.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    dtype -- Is the NumPy dtype to test.
    """

    swappingtype = checkDtypeSwapping(dtype)
    if swappingtype == -2:
        raise ValueError('Structured dtype has commponents with different '
                         'endianness')
    return swappingtype < 0


def checkDtypeSwapping(dtype):
    """Check the byte swapping of the given NumPy dtype.

    This function was added in version 0.1.0 of Mrc.py as
    included with Priism.

    Positional parameters:
    dtype -- Is the NumPy dtype to check.

    Return value:
    The return value will be one of the following:

    1 -- The dtype is compatible with the native byte order.  It has
    one or more components that are in the native byte order and all
    remaining components are not sensitive to the byte order.

    0 -- The dtype is compatible with the native byte order.  It only
    has components that are not sensitive to the byte order.

    -1 -- The dtype requires byte swapping to be compatible with the
    native byte order.  It has one or more components that are not in
    the native byte order and all remaining components are not sensitive
    to the byte order.

    -2 -- The dtype does not have a consistent byte ordering.  At least
    one component is in the native byte order, and at least one component
    is byte swapped relative to the native byte order.
    """

    if dtype.fields is None:
        if dtype.byteorder == '=' or dtype.byteorder == _SYSBYTEORDER:
            rval = 1
        elif dtype.byteorder == '|':
            rval = 0
        else:
            rval = -1
    else:
        rval = 0
        for f in dtype.fields:
            rval_child = checkDtypeSwapping(dtype.fields[f][0])
            if rval_child != 0 and rval != rval_child:
                if rval == 0:
                    rval = rval_child
                else:
                    rval = -2
                    break
    return rval


# The second element in each tuple is the name by which that field will be
# accessed.  If changed, similar changes will be necessary in the
# property definitions for the HdrBase or HdrPriism classes.
mrcHdrFields = [
    ('3i4', 'Num'),
    ('1i4', 'PixelType'),
    ('3i4', 'mst'),
    ('3i4', 'm'),
    ('3f4', 'd'),
    ('3f4', 'angle'),
    ('3i4', 'axis'),
    ('3f4', 'mmm1'),
    ('1i4', 'nspg'),
    ('1i4', 'next'),
    ('1i2', 'dvid'),
    ('1i2', 'nblank'),
    ('1i4', 'ntst'),
    ('24i1', 'blank'),
    ('1i2', 'NumIntegers'),
    ('1i2', 'NumFloats'),
    ('1i2', 'sub'),
    ('1i2', 'zfac'),
    ('2f4', 'mm2'),
    ('2f4', 'mm3'),
    ('2f4', 'mm4'),
    ('1i2', 'ImageType'),
    ('1i2', 'LensNum'),
    ('1i2', 'n1'),
    ('1i2', 'n2'),
    ('1i2', 'v1'),
    ('1i2', 'v2'),
    ('2f4', 'mm5'),
    ('1i2', 'NumTimes'),
    ('1i2', 'ImgSequence'),
    ('3f4', 'tilt'),
    ('1i2', 'NumWaves'),
    ('5i2', 'wave'),
    ('3f4', 'zxy0'),
    ('1i4', 'NumTitles'),
    ('(10,80)i1', 'title'),
]

mrcHdrNames = []
mrcHdrFormats = []
for ff in mrcHdrFields:
    mrcHdrFormats.append(ff[0])
    mrcHdrNames.append(ff[1])
del ff
del mrcHdrFields
mrcHdr_dtype = list(zip(mrcHdrNames, mrcHdrFormats))

# This describes the MRC 2014 , http://www.ccpem.ac.uk/mrc_format/mrc2014.php ,
# format with the addition of fields from Priism that do not conflict with the
# MRC 2014 format.  The fields that are Priism extensions are ntst,
# NumIntegers, NumFloats, sub, zfac, ImageType, LensNum, n1, n2, v1, v2, and
# tilt.
mrc2014HdrFields = [
    ('3i4', 'Num'),
    ('1i4', 'PixelType'),
    ('3i4', 'mst'),
    ('3i4', 'm'),
    ('3f4', 'd'),
    ('3f4', 'angle'),
    ('3i4', 'axis'),
    ('3f4', 'mmm1'),
    ('1i4', 'nspg'),
    ('1i4', 'next'),
    ('4i1', 'extra0'),
    ('1i4', 'ntst'),
    ('4i1', 'exttyp'),
    ('1i4', 'nversion'),
    ('16i1', 'extra1'),
    ('1i2', 'NumIntegers'),
    ('1i2', 'NumFloats'),
    ('1i2', 'sub'),
    ('1i2', 'zfac'),
    ('24i1', 'extra2'),
    ('1i2', 'ImageType'),
    ('1i2', 'LensNum'),
    ('1i2', 'n1'),
    ('1i2', 'n2'),
    ('1i2', 'v1'),
    ('1i2', 'v2'),
    ('12i1', 'extra3'),
    ('3f4', 'tilt'),
    ('3f4', 'origin'),
    ('4i1', 'map'),
    ('4u1', 'machst'),
    ('1f4', 'rms'),
    ('1i4', 'NumTitles'),
    ('(10,80)i1', 'title'),
]

mrc2014HdrNames = []
mrc2014HdrFormats = []
for ff in mrc2014HdrFields:
    mrc2014HdrFormats.append(ff[0])
    mrc2014HdrNames.append(ff[1])
del ff
del mrc2014HdrFields
mrc2014Hdr_dtype = list(zip(mrc2014HdrNames, mrc2014HdrFormats))

# Set character used for NumPy dtypes that corresponds to the native byte
# ordering.  Also record the values to use in the MRC 2014 machine stamp
# field if a file is in the native ordering or if it is byte-swapped relative
# to the native ordering.  Since the MRC 2014 documentation says that
# a machine stamp of 68 and 68 is safe for specifying little-endian, ignore the
# least significant 4 bits of the second byte (they specify how characters are
# encoded) when testing for a little-endian stamp..  For writing, follow the
# CCP4 documentation and use 68 and 65 for a little-endian machine stamp.
if sys.byteorder == 'little':
    _SYSBYTEORDER = '<'
    _SYSMACHST0 = 68
    _SYSMACHST1 = 65
    _SYSMACHST1LO = 64
    _SYSMACHST1HI = 79
    _SYSSWAPST0 = 17
    _SYSSWAPST1 = 17
    _SYSSWAPST1LO = 17
    _SYSSWAPST1HI = 17
else:
    _SYSBYTEORDER = '>'
    _SYSMACHST0 = 17
    _SYSMACHST1 = 17
    _SYSMACHST1LO = 17
    _SYSMACHST1HI = 17
    _SYSSWAPST0 = 68
    _SYSSWAPST1 = 65
    _SYSSWAPST1LO = 64
    _SYSSWAPST1HI = 79
