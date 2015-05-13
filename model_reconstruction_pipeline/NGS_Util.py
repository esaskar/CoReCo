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
    path = ""

    length = len(dirPath)

    if dirPath.find("/", length-1) == -1:
        path = dirPath + "/" + fileName
    else:
        path = dirPath + fileName

    return path

def createDirectoryPath(path,dirName):
    dirPath = ""

    length = len(path)

    if path.find("/", length-1) == -1:
        dirPath = path + "/" + dirName + "/"
    else:
        dirPath = path + dirName + "/"

    return dirPath

def createDirectory(dirPath):
    if (dirPath != None):

        if not os.path.exists(dirPath):

            os.makedirs(dirPath)

        return dirPath

    return None

def logExecutionCall(call):
    print call

    fw = open(executionLog,"a")

    fw.write(call + "\n\n")

    fw.close()

def logCall(call):
    fw = open(traceLog,"a")

    fw.write(call + "\n\n")

    fw.close()

def executeCall(call):
    executionLogFile = open(executionLog,"a")
    errorLogFile = open(errorLog,"a")
    executionLogFile.write(call + "\n\n")
    errorLogFile.write(call + "\n\n")
    rc = subprocess.call(call, stdin=executionLogFile, stderr=errorLogFile, shell=True)
    if rc > 0:
        raise Exception("System call failed (%d): %s" % (rc, call))

def moveDirectory(src,dest):
    if ( os.path.exists(src) and os.path.exists(dest) ):

        shutil.move(src, dest)


def zipDirectory(zdir):
    if os.path.exists( zdir ):

        call = "gzip -r " + zdir

        executeCall(call)

def unzipFileToDifferentLocation(zfile, uzFile):
    if ( (os.path.exists(zfile)) and (not os.path.exists(uzFile)) ):

        call = "gzip -d -c " + zfile + " > " +  uzFile

        executeCall(call)

def unzipFile(zfile):
    print zfile + "\n"


    if os.path.exists(zfile):

        call = "gzip -d " + zfile

        executeCall(call)


####################################################################################################################################################################################################################################################################################################################################
####################################################################################################################################################################################################################################################################################################################################

   
    
def unzipDirectory(zdir):
    if os.path.exists( zdir ):

        call = "gzip -d " + zdir

        executeCall(call)
    
def findMissingSampleLines(inFile, fromFile):
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
        
def removeFile(file):
    if (file != None):

        if os.path.exists(file):

            os.remove(file)

        return True

    return False
    
def fileWriteLine(file, text):
    if (file != None):

        fileHandle = open(file,"a")
        fileHandle.write(text + "\n")
        fileHandle.close()

        return True

    return False
        

def zipFile(file):
    if os.path.exists(file):

        call = "gzip " + file

        executeCall(call)


def copyFile(src,dest):
    if ( os.path.exists(src) ):

        shutil.copy(src, dest)

def executeThreadedCall(call):
    print call + "\n\n"

    lock = threading.Lock()

    lock.acquire()

    executionLogFile = open(executionLog,"a")

    executionLogFile.write(call + "\n\n")

    executionLogFile.close()

    lock.release()
    subprocess.call(call, shell=True)



    
