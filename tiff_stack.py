'''
Created on Nov 19, 2012

@author: Mirna Lerotic
'''

import numpy as np
from PIL import Image        
        
class TiffStackWrapper:
    def __init__(self,fname):
        '''fname is the full path '''
        self.im  = Image.open(fname)
        
        self.im.seek(0)
        # get image dimensions from the meta data the order is flipped
        # due to row major v col major ordering in tiffs and numpy
        self.im_sz = [self.im.tag[0x101][0],
                      self.im.tag[0x100][0]]
        self.cur = self.im.tell()

    def get_mode(self):
        return self.im.mode
        
    def get_frame(self,j):
        '''Extracts the jth frame from the image sequence.
        if the frame does not exist return None'''
        try:
            self.im.seek(j)
        except EOFError:
            return None

        self.cur = self.im.tell()
        return np.reshape(self.im.getdata(),self.im_sz)
    
    def __iter__(self):
        self.im.seek(0)
        self.old = self.cur
        self.cur = self.im.tell()
        return self

    def next(self):
        try:
            self.im.seek(self.cur)
            self.cur = self.im.tell()+1
        except EOFError:
            self.im.seek(self.old)
            self.cur = self.im.tell()
            raise StopIteration
        return np.reshape(self.im.getdata(),self.im_sz)        
        
