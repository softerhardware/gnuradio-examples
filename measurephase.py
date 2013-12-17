#!/usr/bin/python2
#

import os
from gnuradio import gr
from gnuradio import audio
from gnuradio import blocks
from gnuradio import digital
from gnuradio import analog
from gnuradio import filter
from gnuradio import fft
import scikits.audiolab as al

from gnuradio.eng_option import eng_option
from gnuradio import analog
import sys
import numpy as np

try:
    from gnuradio import qtgui
    from PyQt4 import QtGui, QtCore
    import sip
except ImportError:
    sys.stderr.write("Error: Program requires PyQt4 and gr-qtgui.\n")
    sys.exit(1)


rate = 96000




class readwave(gr.sync_block):
  "Reads 24bit 4 channel wave file"
  def __init__(self,fn):
    gr.sync_block.__init__(
      self,
      name = "readwave",
      in_sig = None,
      out_sig =[np.float32,np.float32,np.float32,np.float32],
    )
    self.fn = fn
    self.f = al.Sndfile(self.fn,'r')
    #self.fl = self.f.nframes
    #self.data = self.f.read_frames(self.fl)
    #self.ci = 0
    #f.close()


  def work(self, input_items, output_items):

    try:
      data = self.f.read_frames(len(output_items[0]))
  
      output_items[0][:] = data[:,0]
      output_items[1][:] = data[:,1]
      output_items[2][:] = data[:,1]
      output_items[3][:] = data[:,3]

      
      return len(output_items[0])
    except:
      print("Done!")
      self.f.close()
      return -1





class findphase_c(gr.sync_block):
  "Finds phase in a vector"
  def __init__(self,inpvecsz):
    gr.sync_block.__init__(
      self,
      name = "findphase_c",
      in_sig = [(np.complex64,inpvecsz),(np.complex64,inpvecsz)],
      out_sig =[(np.float32,20)],
    )

  def work(self, input_items, output_items):
    v0 = input_items[0][0]
    v1 = input_items[1][0]

    v0inds = np.argsort(np.abs(v0))
    #v1inds = np.argsort(np.absolute(v1))

    ## Use 20 largest values from v0
    for i in range(0,20):
      j = v0inds[-(i+1)]
      #print v0[j],v1[j]
      ang = np.angle(v0[j] * np.conj(v1[j]),deg=True)
      output_items[0][0][i] = ang
      #print ang

    return len(output_items)




class RX(gr.hier_block2):
  def __init__(self,i):
    gr.hier_block2.__init__(self, "RX",
      gr.io_signature(2,2,gr.sizeof_float),
      gr.io_signature(1,1,gr.sizeof_gr_complex))


    self.i = i

    self.BuildXlate(True,True)
    

  def BuildXlate(self,scope=True,tscope=True):

    self.f2c = blocks.float_to_complex()
    
    taps = filter.firdes.low_pass ( 1.0, rate, 2000, 750, filter.firdes.WIN_HAMMING )
    self.xlate = filter.freq_xlating_fir_filter_ccf ( 12, taps, 10000, rate )
    self.connect( self,self.f2c, self.xlate, self )
    self.connect( (self,1), (self.f2c,1) )

    if scope:
      self.scope = qtgui.sink_c(1024, filter.firdes.WIN_BLACKMAN_hARRIS, fc=0, bw=rate/12, name="FFT%d" % self.i, plotfreq=True, 
        plotwaterfall=True, plottime=True, plotconst=True)

      self.connect( self.xlate, self.scope )

      # Get the reference pointer to the SpectrumDisplayForm QWidget
      self.pyobj = sip.wrapinstance(self.scope.pyqwidget(), QtGui.QWidget)
      self.pyobj.show()


    if tscope:
      self.tscope = qtgui.sink_c(1024, filter.firdes.WIN_BLACKMAN_hARRIS, fc=0, bw=rate, name="Top FFT%d" % self.i, plotfreq=True, 
        plotwaterfall=True, plottime=False, plotconst=False)
   
      self.connect(self.f2c,self.tscope)

      self.tpyobj = sip.wrapinstance(self.tscope.pyqwidget(), QtGui.QWidget)
      self.tpyobj.show() 




class MeasurePhase(gr.top_block):
  def __init__(self):
    gr.top_block.__init__(self)

    self.qapp = QtGui.QApplication(sys.argv)

    self.sdrsource = audio.source(rate,"MeasurePhase")
    #11556698
    #self.sdrsource = readwave("/home/shaynal/beamforming/wspr/wspr11556698.wav")
    #self.sdrsource = blocks.throttle(4,96000)
    #self.connect(self.sdrsourcea,self.sdrsource) 

    self.rx0 = RX(0)
    self.rx1 = RX(1)

    self.connect( (self.sdrsource,0), (self.rx0,0) )
    self.connect( (self.sdrsource,1), (self.rx0,1) )
    self.connect( (self.sdrsource,2), (self.rx1,0) )
    self.connect( (self.sdrsource,3), (self.rx1,1) )

    self.BuildConjMult()
    self.BuildFFT()



  def BuildConjMult(self,scope=True):
    ## Compute angle difference
    #self.conj = blocks.conjugate_cc()
    self.mult = blocks.multiply_conjugate_cc()

    self.connect(self.rx0, blocks.multiply_const_cc(10000), (self.mult,0))
    self.connect(self.rx1, blocks.multiply_const_cc(10000), (self.mult,1))

    self.histo = qtgui.histogram_sink_f(1000,360,-179,180,"Histogram")
    #self.histo.enable_autoscale(False)
    self.histo.enable_accumulate(True)
    self.histo.enable_grid(True)
    #self.histo.enable_menu(True)
    self.connect(self.mult,blocks.complex_to_arg(), blocks.multiply_const_ff(180.0/np.pi),self.histo)

    self.pyobj = sip.wrapinstance(self.histo.pyqwidget(), QtGui.QWidget)
    self.pyobj.show()


  def BuildFFT(self,scope=True):
    fft0w = filter.window.blackmanharris(2048)
    fft0 = fft.fft_vcc(2048, True, fft0w, True)

    fft1w = filter.window.blackmanharris(2048)
    fft1 = fft.fft_vcc(2048, True, fft1w, True)

    self.connect(self.rx0,blocks.stream_to_vector(gr.sizeof_gr_complex, 2048),fft0)
    self.connect(self.rx1,blocks.stream_to_vector(gr.sizeof_gr_complex, 2048),fft1)

    v2s = blocks.vector_to_stream(gr.sizeof_float,20)

    fp2 = findphase_c(2048)
    self.connect(fft0,(fp2,0))
    self.connect(fft1,(fp2,1))
    self.connect(fp2,v2s)

    if scope:
      self.ffth = qtgui.histogram_sink_f(100,360,-179,180,"FFT Phase Histogram")
      #self.ffth.enable_autoscale(False)
      self.ffth.enable_accumulate(True)
      self.ffth.enable_grid(True)
      #self.histo.enable_menu(True)
      self.connect(v2s,self.ffth)

      self.ffthqt = sip.wrapinstance(self.ffth.pyqwidget(), QtGui.QWidget)
      self.ffthqt.show()




  
  def JackConnect(self):

    for i in range(0,4):
      os.system("jack_connect system:capture_%d MeasurePhase:in%d" % (i+1,i))


    
 

if __name__ == '__main__':
  fg = MeasurePhase();
  fg.start()
  fg.JackConnect()  
  fg.qapp.exec_()
  fg.stop()
