
from PyQt5 import QtCore
from PyQt5.QtCore import pyqtSignal, pyqtSlot
from queue import SimpleQueue, Empty
import os
import time
import numpy as np
from scipy import ndimage
from skimage import filters
from skimage.registration import phase_cross_correlation

class GeneralPurposeSignals(QtCore.QObject):
    finished = pyqtSignal()
    ithetaprogress = pyqtSignal(int)

class GeneralPurposeProcessor(QtCore.QRunnable):
    def __init__(self, parent, queue):
        super(GeneralPurposeProcessor, self).__init__()
        self.signals = GeneralPurposeSignals()
        self.funcdict = {"ShiftImg": self.ShiftImg, "AlignReferenced": self.AlignReferenced}
        self.parent = parent
        self.queue = queue
        self.current_itheta = 0

    @pyqtSlot()
    def run(self):
        #print('worker', threading.get_ident())
        while True:
            #print(self.parent.pool.pool.activeThreadCount())
            #print(QtCore.QThreadPool.activeThreadCount())
            try:
                workerfunc, *args = self.queue.get(False)
                self.funcdict[workerfunc](*args[0])
                #print("busy with task {}".format(workerfunc))
            except Empty: # if queue empty
                break

    def AlignReferenced(self, data, itheta):
        if self.current_itheta != itheta:
            self.current_itheta = itheta
            self.signals.ithetaprogress.emit(itheta)
        ref_img = self.EdgeDetect(self.Gauss(self.parent.stack.absdata_cropped[:, :, data[0],itheta]))
        mov_img = self.EdgeDetect(self.Gauss(self.parent.stack.absdata_cropped[:, :, data[1],itheta]))
        upsampling = self.UpsamplingFactor()
        drift, error, _ = phase_cross_correlation(ref_img, mov_img,upsample_factor=upsampling,normalization=None)
        self.parent.stack.shiftsdict[itheta]["errors"][data[0]] = round(error,4)
        if self.parent.cb_upsampling.isChecked():
            self.parent.stack.shiftsdict[itheta]["xdots"][data[0]] = round(drift[0],2)
            self.parent.stack.shiftsdict[itheta]["ydots"][data[0]] = round(drift[1],2)
        else:
            self.parent.stack.shiftsdict[itheta]["xdots"][data[0]] = round(drift[0],0)
            self.parent.stack.shiftsdict[itheta]["ydots"][data[0]] = round(drift[1],0)
    def ShiftImg(self, row,x,y,itheta):
        borders, padded = self.PadImg(self.parent.stack.absdata4d[:, :, row,itheta],-x,-y)
        if self.parent.cb_upsampling.isChecked():
            shifted = ndimage.fourier_shift(np.fft.fft2(padded), [float(-x),float(-y)])
        else:
            shifted = ndimage.fourier_shift(np.fft.fft2(padded), [int(round(-x,0)),int(round(-y,0))])
        shifted = np.fft.ifft2(shifted)
        self.parent.stack.absdata4d_shifted[:, :, row, itheta] = shifted.real[borders[0]:padded.shape[0]-borders[1],borders[2]:padded.shape[1]-borders[3]]
        return

    def PadImg(self, img, x, y):
        default = 10 # minimum expansion of image
        borders = 4*[default]
        if x<0:
            borders[1] = abs(int(np.floor(x)))+default
        elif x>0:
            borders[0] = abs(int(np.ceil(x)))+default
        if y<0:
            borders[3] = abs(int(np.floor(y)))+default
        elif y>0:
            borders[2] = abs(int(np.ceil(y)))+default
        padded = np.pad(img,((borders[0],borders[1]),(borders[2],borders[3])),mode = "edge")
        return (borders, padded)
    def Gauss(self, im):
        gaussed = ndimage.gaussian_filter(im, self.parent.spinBoxGauss.value())
        return gaussed
    def EdgeDetect(self, im):
        if self.parent.cb_edgedetect.isChecked():
            im = filters.farid(im)
        return im
    def UpsamplingFactor(self):
        if self.parent.cb_upsampling.isChecked():
            fac = 20 ## equal to 0.05 px precision
        else:
            fac = 5 ## equal to 0.2 px precision
        return fac

class TaskDispatcher(QtCore.QObject):
    def __init__(self,parent):
        print("0 - Task dispatcher called, pool initiated.")
        super(TaskDispatcher, self).__init__()
        try:
            self.worker.signals.ithetaprogress.disconnect()
            self.worker.signals.finished.disconnect()
        except:
            pass
        self.pool = QtCore.QThreadPool.globalInstance()
        try:
            cpus = len(os.sched_getaffinity(0)) # number of cpu threads. not supported on some platforms.
        except:
            cpus = os.cpu_count()
        self.pool.setMaxThreadCount(cpus)
        self.queue = SimpleQueue()
        self.parent = parent

    @pyqtSlot()
    def run(self):
        #print(self.queue.qsize())
        #print('pool', threading.get_ident())
        # if self.parent.com.stack_4d == 0:
        qsize = int(self.queue.qsize())
        maxthreads = self.pool.maxThreadCount()
        preferred_thread_number = min(maxthreads, qsize)
        print("2 - Starting threads: ",preferred_thread_number, "threads needed, number of tasks:",qsize)
        while int(self.queue.qsize()) and self.pool.activeThreadCount() < preferred_thread_number: #start as many threads as needed.
            #print("active threads"+str(self.pool.activeThreadCount())+" qsize "+str(int(self.queue.qsize())))
            worker = GeneralPurposeProcessor(self.parent,self.queue)
            worker.signals.ithetaprogress.connect(self.parent.IThetaProgress)
            self.pool.start(worker)
            time.sleep(0.02) # artifical delay to start threads with a little time separation. Otherwise ShiftImgs freezes in Win10 and MacOS for small stacks!
        self.pool.waitForDone()
