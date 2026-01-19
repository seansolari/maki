#!/usr/bin/env python3
import gzip
import os
import queue
import re
import shutil
import subprocess
from argparse import ArgumentParser
from multiprocessing import Process, Queue
from shutil import which
from typing import List

COMPRESSION_EXTNS = ['.tar.gz', '.tar', '.gz']
SEQ_EXTNS = [".fna", ".fasta", ".fa"]
EXTS_TO_DELETE = [".embl", ".fna", ".faa", ".ffn", ".gbff", ".hypotheticals.faa", ".hypotheticals.tsv"]
EXTS_TO_ZIP = [".json", ".inference.tsv", ".log", ".tsv"]
EXTS_TO_COPY = [".txt"]
EXTS_TO_FILTER = [".gff3"]

def checkBakta() -> bool:
    return which("bakta") is not None

def parseArgs():
    parser = ArgumentParser()
    parser.add_argument("-f", "--file", dest="inputFile", type=str, required=True, help="File containing paths to genome sequences")
    parser.add_argument("-o", "--out", dest="outputDir", type=str, required=True, help="Output folder")
    parser.add_argument("--db", dest="db", type=str, required=True, help="Bakta database")
    parser.add_argument("--tasks", dest="tasks", type=int, help="Number of parallel Bakta instances to run")
    parser.add_argument("--task-threads", dest="threadsPerTask", type=int, help="Number of threads per Bakta instance")
    parser.add_argument('--keep', dest='keepTemp', default=False, action='store_true')
    parser.add_argument("--filter", dest="filterStrings", action="append", help="Regex filters")
    parser.add_argument("--tmp-dir", dest="tempDir", type=str, help="Temporary directory path")
    return parser.parse_args()

def readLines(file: str) -> List[str]:
    with open(file, 'r') as f:
        return [line.strip() for line in f]

def stripExtensions(s: str) -> str:
    for suffix in COMPRESSION_EXTNS:
        if s.endswith(suffix):
            s = s[:-len(suffix)]
            break

    for suffix in SEQ_EXTNS:
        if s.endswith(suffix):
            s = s[:-len(suffix)]
            break

    return s

class BaktaProcess(Process):
    def __init__(self, workQueue: Queue, outputDir: str, dbFile: str, threads: int, keepTemp: bool, filterStrings: List[str], tempDir: str, *args, **kwargs):
        super(BaktaProcess, self).__init__(*args, **kwargs)
        self.queue = workQueue
        self.outputDir = outputDir
        self.dbFile = dbFile
        self.threads = threads
        self.keepTemp = keepTemp
        self.filterStrings = filterStrings
        self.tempDir = tempDir

    def run(self):
        while 1:
            try:
                outputName, fastaFile = self.queue.get(block=True, timeout=0.1)
            except queue.Empty:
                break

            # run Bakta and place in temp folder

            outDir = os.path.join(self.outputDir, outputName)
            tempOutDir = outDir + ".temp"

            try:
                self.runBakta(tempOutDir, outputName, fastaFile)
            except subprocess.CalledProcessError:
                print("[error] failed %s, skipping..." % outputName)
                continue

            print("[main] compressing and filtering bakta output")

            # create output folder

            os.mkdir(outDir)
            
            src = os.path.join(tempOutDir, outputName)
            dest = os.path.join(outDir, outputName)
            
            # (re)move auxilliary files
            
            if not self.keepTemp:
                self.removeFiles(src, *EXTS_TO_DELETE)
            else:
                self.moveAndGZip(src, dest, *EXTS_TO_DELETE)
            
            # copy and zip files
            
            self.moveAndGZip(src, dest, *EXTS_TO_ZIP)
            
            # copy files
            
            self.moveFiles(src, dest, *EXTS_TO_COPY)
            
            # filter file
            
            self.moveAndFilter(src, dest, self.filterStrings, *EXTS_TO_FILTER)
            
            # remove temp dir
            
            shutil.rmtree(tempOutDir)
    
    def runBakta(self, outDir: str, outPrefix: str, inputFile: str) -> None:
        baktaCommand = [
            "bakta",
            "--db", self.dbFile,
            "--verbose",
            "--output", outDir,
            "--prefix", outPrefix,
            "--threads", str(self.threads),
            "--keep-contig-headers",
            "--skip-plot",
            "--tmp-dir", self.tempDir,
            inputFile
        ]
        
        print("[main] running command: %s" % " ".join(baktaCommand))
        subprocess.check_call(baktaCommand)

    @staticmethod
    def removeFiles(srcPrefix: str, *exts: str) -> None:
        for sfx in exts:
            os.remove(srcPrefix + sfx)

    @staticmethod
    def moveAndGZip(srcPrefix: str, destPrefix: str, *exts: str) -> None:
        for sfx in exts:
            srcFile = srcPrefix + sfx
            destFile = destPrefix + sfx + ".gz"
            with open(srcFile, 'rb') as src, gzip.open(destFile, 'wb') as dst:
                dst.writelines(src)
            os.remove(srcFile)

    @staticmethod
    def moveFiles(srcPrefix: str, destPrefix: str, *exts: str) -> None:
        for sfx in exts:
            shutil.move(srcPrefix + sfx, destPrefix + sfx)

    @staticmethod
    def moveAndFilter(srcPrefix: str, destPrefix: str, filters: List[str], *exts: str):
        for sfx in exts:
            srcFile = srcPrefix + sfx
            destFile = destPrefix + sfx
            with open(srcFile, 'r') as src, open(destFile, 'w') as dst:
                for line in src:
                    if line.startswith("#"):
                        dst.write(line)
                    else:
                        for filterStr in filters:
                            if re.search(filterStr, line) is not None:
                                print("[main] line caught by filter %s\n\t%s" % (filterStr, line))
                                dst.write("#" + line)
                                break
                        else:
                            dst.write(line)

def main():
    if not checkBakta():
        print("[error] could not find bakta")
        exit(1)

    args = parseArgs()

    # find queries

    inputFiles = {
        (stripExtensions(os.path.basename(f)), f)
        for f in readLines(args.inputFile)
    }

    outputDir = os.path.abspath(args.outputDir)

    if not os.path.exists(outputDir):
        os.mkdir(outputDir)

    if not os.path.exists(args.tempDir):
        os.mkdir(args.tempDir)

    # remove any that have already been processed

    existingGenomes = {
        d
        for d in os.listdir(outputDir)
        if os.path.isdir(os.path.join(outputDir, d))
    }

    if existingGenomes:
        oldSize = len(inputFiles)
        inputFiles = {f for f in inputFiles if f[0] not in existingGenomes}
        newSize = len(inputFiles)

        if newSize < oldSize:
            print("[main] %d inputs already processed, %d remaining" % (oldSize - newSize, newSize))
    
    # processing

    inputQueue = Queue()
    
    for data in inputFiles:
        inputQueue.put(data)

    processes = [
        BaktaProcess(inputQueue, outputDir, args.db, args.threadsPerTask, args.keepTemp, args.filterStrings, args.tempDir)
        for _ in range(args.tasks)
    ]

    for p in processes:
        p.start()

    for p in processes:
        p.join()

if __name__ == "__main__":
    main()
