#!/usr/bin/env python
"""
Author: M. Fahad Syed (fahad.syed@vtt.fi)
"""



import os
import sys
import subprocess
import operator
import traceback
import time
import shutil

import threading


traceLog     = "../TraceLog.txt"
executionLog = "../ExecutionLog.txt"
errorLog     = "../ErrorLog.txt"


def createFilePath(dirPath,fileName):

    try:
        path = ""
        
        length = len(dirPath)
        
        if dirPath.find("/", length-1) == -1:
            path = dirPath + "/" + fileName
        else:
            path = dirPath + fileName
            
        return path
        
    except Exception:
        print traceback.print_exc()
        return None

def createDirectoryPath(path,dirName):

    try:
        dirPath = ""
        
        length = len(path)
        
        if path.find("/", length-1) == -1:
            dirPath = path + "/" + dirName + "/"
        else:
            dirPath = path + dirName + "/"
            
        return dirPath
        
    except Exception:
        print traceback.print_exc()
        return None


def createDirectory(dirPath):

    try:
    
        if (dirPath != None):
            
            if not os.path.exists(dirPath):
                            
                os.makedirs(dirPath)
                
            return dirPath
        
        return None

    except Exception:
        print traceback.print_exc()
        return None


def logExecutionCall(call):

    try:
        
        print call
    
        fw = open(executionLog,"a")
        
        fw.write(call + "\n\n")

        fw.close()
        
    except Exception:
        print traceback.print_exc()


def logCall(call):

    try:
   
        fw = open(traceLog,"a")
        
        fw.write(call + "\n\n")

        fw.close()
        
    except Exception:
        print traceback.print_exc()


def executeCall(call):
    
    try:
                
        print call + "\n\n"

        executionLogFile = open(executionLog,"a")

        errorLogFile = open(errorLog,"a")
        
        executionLogFile.write(call + "\n\n")

        errorLogFile.write(call + "\n\n")
        
        subprocess.call(call, stdin=executionLogFile, stderr=errorLogFile, shell=True)
        
        errorLogFile.close()
        
        executionLogFile.close()
    
    except Exception:

        print traceback.print_exc()


def moveDirectory(src,dest):

    try:

        if ( os.path.exists(src) and os.path.exists(dest) ):
        
            shutil.move(src, dest)
                    
    except Exception:
        print traceback.print_exc()



def zipDirectory(zdir):

    try:
        
        if os.path.exists( zdir ):
        
            call = "gzip -r " + zdir
        
            executeCall(call)
        
    except Exception:
        print traceback.print_exc()


def unzipFileToDifferentLocation(zfile, uzFile):

    try:
        if ( (os.path.exists(zfile)) and (not os.path.exists(uzFile)) ):
            
            call = "gzip -d -c " + zfile + " > " +  uzFile
        
            executeCall(call)
        
    except Exception:
        print traceback.print_exc()


def unzipFile(zfile):

    try:
        print zfile + "\n"


        if os.path.exists(zfile):

            call = "gzip -d " + zfile
        
            executeCall(call)
        
    except Exception:
        print traceback.print_exc()


####################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################

   
    
def unzipDirectory(zdir):

    try:
        
        if os.path.exists( zdir ):
        
            call = "gzip -d " + zdir
        
            executeCall(call)
        
    except Exception:
        print traceback.print_exc()

    
def findMissingSampleLines(inFile, fromFile):
    
    try:
        
        fromFile = open(fromFile, "r")        
        inFile   = open(inFile,"r")
        
        fromFileList = fromFile.readlines()
        inFileList   = inFile.readlines()
        
        localSampleName = ""
        
        fromFileDataList = []
        toAppend    = []
        
                
        for fromFileLine in fromFileList:
            
            found = False
            
            fromFileDataList = fromFileLine.split("\t")
            
            fromFileSampleName = fromFileDataList[0]
            
                        
            for inFileLine in inFileList:
                
                inFileLineSampleName = inFileLine.split("\t")[0]
                
                if fromFileSampleName == inFileLineSampleName:
                    found = True
                    break
            
            if not found:
                                    
                toAppend.append(fromFileLine)

        
        fromFile.close()
        inFile.close()
        
        return toAppend
    
    except Exception:

        print traceback.print_exc()
        
        return []

    
def removeFile(file):

    try:
    
        if (file != None):
            
            if os.path.exists(file):
                            
                os.remove(file)
                
            return True
        
        return False

    except Exception:
        print traceback.print_exc()
        return False
    
    
def fileWriteLine(file, text):

    try:
    
        if (file != None):
            
            fileHandle = open(file,"a")
            fileHandle.write(text + "\n")
            fileHandle.close()
                
            return True
    
        return False
        

    except Exception:
        print traceback.print_exc()
        return False


def zipFile(file):

    try:
        if os.path.exists(file):

            call = "gzip " + file
        
            executeCall(call)
        
    except Exception:
        print traceback.print_exc()



def copyFile(src,dest):

    try:

        if ( os.path.exists(src) ):
        
            shutil.copy(src, dest)
                    
    except Exception:
        print traceback.print_exc()



def executeThreadedCall(call):
    
    try:
                
        print call + "\n\n"
        
        lock = threading.Lock()
 
        lock.acquire()
        
        executionLogFile = open(executionLog,"a")
       
        executionLogFile.write(call + "\n\n")

        executionLogFile.close()

        lock.release()
        

        subprocess.call(call, shell=True)
        
    
    except Exception:

        print traceback.print_exc()



    
