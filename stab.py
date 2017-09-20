#!/usr/bin/python

import sys
import time
import hashlib
import argparse

def listOfSamplesForSNP(snpList, snpPos):
    resList = []
    for itm in snpList:
        if itm[0] == snpPos:
            resList.append(itm[1])
    return resList

def listOflowCovSNPSamples(calledSNPlist, lowCovListAsTuple, snpPos):
    lowCovRes = []
    for s in lowCovListAsTuple:
        if s[0] in calledSNPlist:
            continue
        else:
            if s[1] == snpPos:
                lowCovRes.append(s[0])
    # print(lowCovRes)
    # time.sleep(4)
    return lowCovRes

def listOfCoverageFromMpileup(cSNP, lcSNP, mpileupDict, snpPos):
    mpileupRes = []
    mpileupRet = []
    for s_key, s_value in mpileupDict.items():
        for m_tuple in s_value: # scan tuple for match to pos, then append if cov > 4
            if m_tuple[0] == snpPos:
                m_tupleInt = int(m_tuple[1])
                if m_tupleInt > 4:
                    mpileupRes.append(s_key)
    if mpileupRes: # if list has been populated, scan entries for matches in prev two lists
        for i in mpileupRes:
            if i in cSNP:
                continue
            else:
                if i in lcSNP:
                    continue
                else:
                    mpileupRet.append(i) # only append if a match has not been made
    else:
        return mpileupRes

    if not mpileupRet: # if there are no matches, return nil
        return mpileupRes
    else:
        return mpileupRet

def parseSNPline(sList, sDict):
    pos = 0
    posStart = True
    retString = ""
    sData = ""
    for stringItem in stringList:
        if stringItem in jDict:
            kStringAsDict = jDict[stringItem]
            if posStart:
                pos = int(kStringAsDict['snpPos'])
                posStart = False
            else:
                if pos != int(kStringAsDict['snpPos']):
                    print('Error, check SNP pos!')
                    print(pos)
                    print(kStringAsDict)
                    sys.exit()
                else:
                    retString += str(kStringAsDict['snpType']) + '\t'
                    if not sData:
                        sData = str(kStringAsDict['snpData']) # will automatically copy over noinfo entry
        else:
            continue
    retString += sData
    return retString

def main():
    parser = argparse.ArgumentParser(description='Snp TAble Builder using python')
    parser.add_argument('-c', '--concatenatedList', help='concatenated list of SNP positions from all samples using samtools, must be concatenated into a single file')
    parser.add_argument('-l', '--lowcovlist', help='low coverage output file from SNP analysis, preferably from DNASTAR.  See README for format instructions')
    parser.add_argument('-u', '--highOrNormalCovList', help='normal or high coverage output file from SNP analysis, preferably from DNASTAR.  See README for format instructions')

    args = parser.parse_args()

    if not args.concatenatedList:
        print('You must provide a concatenatedList file')
        sys.exit()
    if not args.lowcovlist:
        print('You must provide a lowcovlist file')
        sys.exit()
    if not args.highOrNormalCovList:
        print('You must provide a highOrNormalCovList file')
        sys.exit()


    # open and parse lowCov SNP list as [(sampleName, snpPos)]
    # also create list of all sample names
    snpDepthOf1List = []
    sampleNamesAsDict = dict()
    sampleNamesAsList = []
    snpDepthOf1_IDandPosition = []
    lowDItemAsInt = 0
    lowDItemAsString = ""
    with open(args.lowcovlist, 'r') as snpLowC:
        for lowCovItm in snpLowC:
            lowCovItm = lowCovItm.rstrip('\r\n')
            snpDepthOf1List.append(lowCovItm)
    del snpDepthOf1List[0]
    for lowDItem in snpDepthOf1List:
        splitlowDItem = []
        try:
            splitlowDItem = lowDItem.split('\t')
        except AttributeError:
            print(lowDItem)
            sys.exit()
        if not splitlowDItem[1]:
            continue # eliminates empty lines
        else:
            lowDItemAsInt = int(splitlowDItem[1])
            lowDItemAsString = splitlowDItem[0]
            snpDepthOf1_IDandPosition.append((lowDItemAsString, lowDItemAsInt))
            if lowDItemAsString not in sampleNamesAsDict:
                sampleNamesAsDict[lowDItemAsString] = 1
            else:
                sampleNamesAsDict[lowDItemAsString] += 1
    sampleNamesAsList = sampleNamesAsDict.keys()

    # open and parse concatenated SNP file as {'sampleName': [(snpPos, cov), (snpPos, cov),...]}
    # also create list of all sampleNames
    snpConcatenatedList = []
    with open(args.concatenatedList, 'r') as snpCList:
        for i in snpCList:
            i = i.rstrip('\r\n')
            snpConcatenatedList.append(i)


    snpMpileupCoverageDict = dict()
    d_key = ""
    d_value = []
    start = 1
    for lItm in snpConcatenatedList:
        lItm4Dict = lItm.split('#')
        if (len(lItm4Dict) > 1):
            if not start:
                snpMpileupCoverageDict[d_key] = d_value
                d_key = ""
                d_value = []
            d_key = lItm4Dict[4]
            continue
        else:
            start = 0
            lItm = lItm.split('\t')
            d_value.append((int(lItm[1]),int(lItm[3])))
    snpMpileupCoverageDict[d_key] = d_value
    # sampleNamesAsList = snpMpileupCoverageDict.keys()

    # DNASTAR data
    listOfResults = []
    snpString = ""
    with open(args.highOrNormalCovList, 'r') as DNASTARsnp:
        for p in DNASTARsnp:
            p = p.rstrip('\n')
            snpString += p
    listOfResults = snpString.split('\r')

    headerLine = listOfResults[0]
    del listOfResults[0]
    listOfResults.pop()
    # listOfUniqueSnpPos = [(uniqueSnpPos1, idx1),(uniqueSnpPos2, idx2),...]
    snpTracker = 73
    listOfUniqueSnpPos = []
    listOfUniqueSnpPos.append((snpTracker, 0))
    ct = 0
    for itmLine in listOfResults:
        itmLine = itmLine.split('\t')
        itmLineValueAsInt = int(itmLine[3])
        if itmLineValueAsInt != snpTracker:
            listOfUniqueSnpPos.append((itmLineValueAsInt, ct))
        snpTracker = itmLineValueAsInt
        ct += 1

    # snpFromDNASTARasTuple = [(snpPos, sampleName)]
    snpFromDNASTARasTuple = []
    snpCt = 0
    for otherItemLine in listOfResults:
        otherItemLine = otherItemLine.split('\t')
        otherItemLineAsInt = int(otherItemLine[3])
        snpFromDNASTARasTuple.append((otherItemLineAsInt, otherItemLine[1]))

    # snpEntryAsTuple = [(snpPos, snpEntryAsString)
    snpEntryAsTuple = []
    for idxVal in listOfUniqueSnpPos:
        indexValue = idxVal[1]
        snpEntryLine = listOfResults[indexValue]
        splitSnpEntryLine = snpEntryLine.split('\t')
        splitSnpEntryAsList = splitSnpEntryLine[3:]
        splitSnpEntryAsString = '\t'.join(splitSnpEntryAsList)
        snpEntryAsTuple.append((idxVal[0], splitSnpEntryAsString))

    namesHeader = ""
    for n in sampleNamesAsList:
        namesHeader += n + '\t'
    print(namesHeader + headerLine)

    # orderedAsEntries
    # [{SampleName : {snpType: snp_lowCov_none, snpData: snpEntryLine, useThisLine: BoolTrueOrFalse}}]

    sampleOrder = []
    ct = 0
    for c in sampleNamesAsList:
        sampleOrder.append((ct, c))
        ct += 1

    snpLL = []
    x = 0
    for snp in listOfUniqueSnpPos:
        llDict = dict()
        llString = hashlib.md5(str(snp[1])).hexdigest()
        snpPos = snp[1]
        sampleDict = dict()
        calledSnpSamples = listOfSamplesForSNP(snpFromDNASTARasTuple, snp[0])
        lowCovSamples = listOflowCovSNPSamples(calledSnpSamples, snpDepthOf1_IDandPosition, snp[0])
        mpileupSamples = listOfCoverageFromMpileup(calledSnpSamples, lowCovSamples, snpMpileupCoverageDict, snp[0])
        for each_snpSample in sampleNamesAsList:
            if each_snpSample in calledSnpSamples:
                sampleDict[each_snpSample] = 'SNP'
        for each_snpSample in sampleNamesAsList:
            if each_snpSample not in sampleDict:
                if each_snpSample in lowCovSamples:
                    sampleDict[each_snpSample] = 'NoCov-SNP'
        for each_snpSample in sampleNamesAsList:
            if each_snpSample not in sampleDict:
                if each_snpSample in lowCovSamples:
                    sampleDict[each_snpSample] = 'NoCov-WT'
        for each_snpSample in sampleNamesAsList:
            if each_snpSample not in sampleDict:
                sampleDict[each_snpSample] = 'NoInfo'
        llDict = {llString : {'snpPos': snpPos, 'sampleData': sampleDict, 'snpData': snpEntryAsTuple[x][1], 'sampleOrder': sampleOrder}}
        snpLL.append(llDict)
        x += 1

    y = 0 # iterator
    for h in snpLL:
        snpPos_unhashed = listOfUniqueSnpPos[y]
        snpPos_hashed = hashlib.md5(str(snpPos_unhashed[1])).hexdigest()
        snpPosData = h[snpPos_hashed]
        orderOfSamples = snpPosData['sampleOrder']
        orderOfSamples.sort()
        sample_dataAsDict = snpPosData['sampleData']
        snpLine_data = snpPosData['snpData']
        outputString = ""
        for i in orderOfSamples:
            sampleNameFromSampleOrder = i[1]
            sampleNameOutput = sample_dataAsDict[sampleNameFromSampleOrder]
            outputString += sampleNameOutput + '\t'
        outputString += snpLine_data
        print(outputString)
        y += 1
if __name__ == "__main__":
    main()
