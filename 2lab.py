from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import numpy as np
import matplotlib.pyplot as plt
from typing import List
from scipy.signal import find_peaks

class Encoding:
    def __init__(self, name, minVal, maxVal):
        self.name = name
        self.maxVal = maxVal
        self.minVal = minVal

class Sequence:
    def __init__(self, name, sequence, quality):
        self.name = name
        self.sequence = sequence
        self.quality = quality

class Group:
    def __init__(self, start, stop, sum):
        self.start = start
        self.stop = stop

def readFastaq(fileName):
    result = []
    with open(fileName) as handle:
        for (name, sequence , quality) in FastqGeneralIterator(handle):
            result.append(Sequence(name, sequence, quality))
    return result

def getMinMax(sequences: List[Sequence]):
    minVal = 256
    maxVal = -1
    for sequence in sequences:
        for char in sequence.quality:
            minVal = min(minVal, ord(char))
            maxVal = max(maxVal, ord(char))

    return (minVal, maxVal)

def getEncoding(sequences: List[Sequence]):
    (minVal, maxVal) = getMinMax(sequences)
    encodings = [
        Encoding("Sanger Phred+33", 33, 73),
        Encoding("Solexa Solexa+64", 59, 104),
        Encoding("Illumina 1.3+ Phred+64", 64, 104),
        Encoding("Illumina 1.5+ Phred+64", 66, 105),
        Encoding("Illumina 1.8+ Phred+33", 33, 74)
    ]
    
    for encoding in encodings:
        if minVal >= encoding.minVal and maxVal <= encoding.maxVal:
            return encoding.name
    return ""

def getFrequencies(sequences: List[Sequence]):
    frequences = [0] * 100
    for sequence in sequences:
        percentage = (sequence.sequence.count("G") + sequence.sequence.count("C")) / len(sequence.sequence) * 100
        frequences[round(percentage)] += 1
    return frequences

def drawGraph(frequencies):
    plt.bar(np.arange(0, 100), frequencies)
    plt.ylabel("Number of Reads")
    plt.xlabel("Percentage of GC")
    plt.show()

def countGC(sequence):
    count = sequence.count("G") + sequence.count("C")
    return count / len(sequence) * 100.0 

def drawPlot():
    with open(filename) as handle:
        for (_, sequence, _) in FastqGeneralIterator(handle):
            gcPercentages.append((sequence.count("G") + sequence.count("C")) / len(sequence) * 100)
            sequenceArr.append(sequence)

    for i in range(0, 100):
        cnt = 0
        for gcPercentage in gcPercentages:
            if round(gcPercentage) == i:
                cnt += 1
        readsNumberArr.append(cnt)

    plt.bar(np.arange(0, 100), readsNumberArr)
    plt.ylabel("Readu skaicius")
    plt.xlabel("Kiek GC yra reade procentais")
    plt.show()

def getFiveSequences(gcPercentage, sequences: List[Sequence]):
    result = []
    for sequence in sequences:
        percentage = (sequence.sequence.count("G") + sequence.sequence.count("C")) / len(sequence.sequence) * 100
        if percentage == gcPercentage:
            result.append(sequence)
        if len(result) == 5:
            break
    return result

result = ""

def doBlastSearch(frequencies, sequences: List[Sequence]):
    peaks, _ = find_peaks(frequencies, height=100, distance=15)
    peakSequences = []
    for peak in peaks:
        peakSequences.extend(getFiveSequences(peak, sequences))

    for i in range(len(peakSequences)):
        print("Performing blast search: " + str(i + 1))
        result = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            hitlist_size=1,
            sequence=peakSequences[i].sequence,
            entrez_query='Bacteria [Organism]')
        
        parsedResults = NCBIXML.parse(result)
        for parsedResult in parsedResults:
            for alignment in parsedResult.alignments:
                print("Found " + str(len(alignment.hsps)) + " matches.")
                for hsp in alignment.hsps:
                    print("Title:", end=' ')
                    print(alignment.title)
                    print("Subject:", end=' ')
                    print(hsp.sbjct)
                    print("Read ID:", end=' ')
                    print(peakSequences[i].name)

fastaq = readFastaq("reads_for_analysis.fastq")
print(getEncoding(fastaq))

frequencies = getFrequencies(fastaq)
drawGraph(frequencies)

doBlastSearch(frequencies, fastaq)