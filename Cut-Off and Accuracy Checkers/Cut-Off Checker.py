import Bio.motifs as motifs
import csv
from operator import itemgetter
import time

class Promoter_Sequence(object):
    def __init__(self):
        self.reset()

    def reset(self):
        self.ID = ""
        self.Cell = ""
        self.Pathway = ""
        self.Reaction = ""
        self.Enzyme = ""
        self.Sequence = ""
        self.Putative_Site = []
        self.Hash = 0
    
    def count_site(self, DatabaseDir):
        Site_Counts = {}
        with open(DatabaseDir) as Database:
            for Motif in motifs.parse(Database, "jaspar"):
                Count = 0
                for Sites in self.Putative_Site:
                    if Motif.name == Sites[0]:
                        Count = Count + 1
                if Count == 0:
                    continue;
                Site_Counts[Motif.name] = Count
        return Site_Counts
    
    def write(self, DatabaseDir):
        Site_Counts = self.count_site(DatabaseDir)
        Total_Count = sum(Site_Counts.values())
        if len(self.Sequence) == 0:
            Site_Index = 0
        else:
            Site_Index = Total_Count/len(self.Sequence)
        Record = {"Cell ID": self.Cell, "Pathway": self.Pathway, "Enzyme": self.Enzyme, 
                  "Reaction": self.Reaction, "SubUnit": self.ID, "Promoter": self.Sequence, 
                  "Sites": self.Putative_Site, "Site Counts": Site_Counts, "Site Index": Site_Index, "Hash" : self.Hash}
        return Record

def Pseudocount_Adder(Motif):
    Experiments = 0
    for Index in range(0,4):
        Experiments += Motif[Index][0]
    Pseudocount = 2.5/100 if Experiments == 1 else 0.25 * (Experiments**0.5)
    return Pseudocount

def core_Pssm(Motif, background={"A":0.25, "C":0.25, "G":0.25, "T":0.25}):
    fqCore = {"A":[], "C":[], "G":[], "T":[]}
    for position in range(Motif.length):  
        pwmColumn = []  
        for index in range(0,4):
            pwmColumn.append(Motif[index][position])
        if max(pwmColumn)/sum(pwmColumn) > 0.95:
            for row, item in zip("ACGT", pwmColumn):
                fqCore[row].append(float(item))
        else:
            for row, item in zip("ACGT", background.values()):
                fqCore[row].append(float(item))  

    coreMatrix = motifs.Motif("ACGT", counts=fqCore)
    corePWM = coreMatrix.counts.normalize()  
    corePSSM = corePWM.log_odds(background={"A":0.31,"C":0.19, "G":0.19, "T":0.31})
    return corePSSM 

def Relative_Score(Score, Max_Score, Min_Score):
    if Max_Score - Min_Score == 0:
        return 0
    Ratio = (Score - Min_Score)/(Max_Score - Min_Score)
    return Ratio

def reverseComplement(string):
    newString = ""
    complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for letter in string:
        newString += complement[letter.upper()]
    
    return newString[::-1]

start = time.time()
databaseDict = {}
pssmDict = {}
coreDict = {}
with open("fullDatabase.txt") as database:
    for Motif in motifs.parse(database, "Jaspar"):
            PWM = Motif.counts.normalize(pseudocounts=Pseudocount_Adder(Motif.counts))
            PSSM = PWM.log_odds(background={"A":0.31,"C":0.19, "G":0.19, "T":0.31})
            CORE = core_Pssm(Motif.counts, background={"A":0.31,"C":0.19, "G":0.19, "T":0.31})
            
            try:
                databaseDict[Motif.name].append(Motif.base_id)
                pssmDict[Motif.base_id] = PSSM
                coreDict[Motif.base_id] = CORE
            except KeyError:
                databaseDict[Motif.name] = [Motif.base_id]
                pssmDict[Motif.base_id] = PSSM
                coreDict[Motif.base_id] = CORE

exonDict = {}
with open("exonSequences.txt") as exonDatabase:
    for line in exonDatabase:
        key, value = line.split(";")
        exonDict[key] = value[:-1] 

Promoter_List = []
with open("Inputs.csv") as Promoter_Data:
    Promoter_Reader = csv.DictReader(Promoter_Data, delimiter=';')
    Temporary_List = ["","","","","","",""]
    Hash_Tracker = [] 
    for Row in Promoter_Reader: 
        Temporary_List[0] = Row["Reference Cell Line"] if Temporary_List[0] == "" else Temporary_List[0]     
        if Row["Hash"] in Hash_Tracker:
            continue;
        Temporary_List[4] = Row["Reference Gene Name"]
        Temporary_List[5] = Row["Reference Sequence"]
        Temporary_List[6] = Row["Hash"]
        if Row["Pathway"] != "":
            Temporary_List[1] = Row["Pathway"]
        if Row["Reaction"] != "":
            Temporary_List[2] = Row["Reaction"]
        if Row["Reference Enzyme"] != "":
            Temporary_List[3] = Row["Reference Enzyme"]
        Promoter_List.append([Temporary_List[0],Temporary_List[1],Temporary_List[2],
                              Temporary_List[3],Temporary_List[4],Temporary_List[5],Temporary_List[6]])
        Hash_Tracker.append(Row["Hash"])

logFile = open("logs.txt", "w")
results = open("rsts.txt", "w")

tfMotifDict = {}
cutoffDict = {}
tfInteractsWith = {}
reversedDatabase = {}
for factor in databaseDict.keys():
    results.write("{} Results:\n".format(factor))
    for factorID in databaseDict[factor]:
        reversedDatabase[factorID] = factor
        cutoffStart = time.time()
        logFile.write("Cut-off check for {} has now started.\n".format(factorID))

        PSSM = pssmDict[factorID]
        CORE = coreDict[factorID]
        tfMotifDict[factorID] = []
        for exon in exonDict.keys():
            if len(exonDict[exon]) < PSSM.length:
                continue
            SCAN = PSSM.search(sequence=exonDict[exon], both=True) 
            for HIT in SCAN:
                sequence = exonDict[exon][HIT[0]:HIT[0]+PSSM.length]
                hitSequence = reverseComplement(sequence) if HIT[0] < 0 else sequence
                coreSCAN = CORE.search(sequence=hitSequence, threshold=0.0, both=False)
                coreHITS = list(coreSCAN)
                if coreHITS == []:
                    continue;
                elif Relative_Score(coreHITS[0][1], CORE.max, 0) > 0.95:
                    tfMotifDict[factorID].append([factorID , HIT[0] , HIT[1] , Relative_Score(HIT[1], PSSM.max, PSSM.min) , Relative_Score(coreHITS[0][1], CORE.max, 0) , exonDict[exon][HIT[0]:HIT[0]+PSSM.length]])
        sortedMotif = sorted(tfMotifDict[factorID], key=itemgetter(3), reverse=True)
        if len(sortedMotif) > 236:
            tfCutOff = sortedMotif[236][3]
        elif len(sortedMotif) > 0:
            tfCutOff = sortedMotif[-1][3]
        else: 
            tfCutOff = 0.75

        cutoffDict[factorID] = tfCutOff
        cutoffEnd = time.time()
        results.write("{} Cut-off Value = {}\n".format(factorID, tfCutOff))
        logFile.write("{} Cut-off check completed in {} secs.\n".format(factorID, cutoffEnd-cutoffStart))

        scanStart = time.time()

        for Item in Promoter_List:
            Promoter = Promoter_Sequence()
            Promoter.Cell = Item[0]
            Promoter.Pathway = Item[1]
            Promoter.Reaction = Item[2]
            Promoter.Enzyme = Item[3]
            Promoter.ID = Item[4]
            Promoter.Sequence = Item[5].lower()
            Promoter.Hash = Item[6]

            SCAN = PSSM.search(sequence=Promoter.Sequence, threshold=0.0, both=True) 
            for HIT in SCAN:
                threshold = tfCutOff
                if Relative_Score(HIT[1], PSSM.max, PSSM.min) > threshold:
                    sequence = Promoter.Sequence[HIT[0]:HIT[0]+PSSM.length]
                    hitSequence = reverseComplement(sequence) if HIT[0] < 0 else sequence
                    coreSCAN = CORE.search(sequence=hitSequence, threshold=0.0, both=False)
                    coreHITS = list(coreSCAN)
                    if coreHITS == []:
                        continue;
                    elif Relative_Score(coreHITS[0][1], CORE.max, 0) > 0.9:
                        try:
                            tfInteractsWith[factorID].append(Promoter.ID)
                            break;
                        except KeyError:
                            tfInteractsWith[factorID] = [Promoter.ID]
                            break;

        scanEnd = time.time()
        logFile.write("{} scan check completed in {} secs.\n".format(factorID, scanEnd-scanStart))

    results.flush()
    logFile.flush()
