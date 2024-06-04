# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 15:19:35 2021

Class for handling hdf5 file I/O for SA based RFI monitor
This was derived from the think RF based monitor
@author: G. Hovey, 6 Nov 2022;Onsala Space Observatory, Sweden
"""
import datetime
import h5py
import numpy as np


class RfiMonFileIO(object):
        
    def __init__(self, parmD, dirspec='./'):
        self.parmD = parmD
        self.fd = None
        self.filespec = ''
        self.dirspec = dirspec
        #locTime = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H_%M_%S")

        
    def openFile(self, filespec='', modeStr='r'):
        """ open a filespec according to mode string
            if filespec is '' then a standard file spec is created from
            dirspec and parmD in __init__
            
            Note if a file is created modeStr = 'w' then the file is empty and
            headers needs to be created
        """
        if filespec =='': #if filespec not given
            filespec = self.makeFilespec() #no so create it
        if self.fd == None:
            try:
                self.fd = h5py.File(filespec, modeStr)
                self.filespec = filespec #If success update internal filespec
            except:
                print('error trying to open ', filespec)
                self.closeFile()
                return ''
        else:
            print('File already open')
        
        return filespec


    def makeFilespec(self):
        """ create a standardized timestamped Filespec"""
        locTime = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H_%M_%S")
        filespec = self.dirspec
        filespec += str(self.parmD['antenna']).replace(' ','')
        filespec += str(self.parmD['band']).replace(' ','') 
        filespec += '-' + locTime + '.hdf5'
        return filespec
                                    
    def reopenFile(self, modeStr):
        if self.fd != None:
            self.closeFile()
        self.openFile(self.filespec, modeStr)
    

    def rdHdr(self):
        if self.filespec == '':
            print('Error filespec not specified')
            return None
        if self.fd == None:        
            self.openFile(self.filespec, modeStr='r')
        
        hdrD = {}
        hdrD.update(self.fd.attrs)
        return hdrD
    

    def closeFile(self):
        if self.fd != None:
            self.fd.close()
            self.fd = None

    def rdField(self, fieldName):
        if self.filespec == '':
            print('Error filespec not specified')
            return None
        if self.fd == None:        
            self.openFile(self.filespec, modeStr='r')
        return self.fd[fieldName]

              
