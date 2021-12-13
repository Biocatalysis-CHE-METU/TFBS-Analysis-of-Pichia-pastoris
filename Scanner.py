import Bio.motifs as motifs
import csv
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

def Relative_Score(Score, Max_Score, Min_Score):
    if Max_Score - Min_Score == 0:
        return 0
    Ratio = (Score - Min_Score)/(Max_Score - Min_Score)
    return Ratio

def read_pathway(Pathway_Dir):
    with open(Pathway_Dir) as csv_file:
        Pathway_Dict = {}
        Reader = csv.DictReader(csv_file, delimiter = ";")
        for Row in Reader:
            if Row["Pathway"] != "":
                Pathway_Dict[Row["Pathway"]] = {}
            else:
                Row["Pathway"] = list(Pathway_Dict.keys())[-1] 
            if Row["Reference Enzyme"] != "":       
                Pathway_Dict[Row["Pathway"]][Row["Reference Enzyme"]] = {}
            else:
                Row["Reference Enzyme"] = list(Pathway_Dict[Row["Pathway"]].keys())[-1]
            Pathway_Dict[Row["Pathway"]][Row["Reference Enzyme"]][Row["Reference Gene Name"]] = {}
        return Pathway_Dict 

def measure_key_frequency(dict_of_dicts):
    frequency = {}
    for d in dict_of_dicts.values():
        for key in d.keys():
            frequency[key] = 1 if key not in frequency else frequency[key] + 1
    return frequency

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

def reverseComplement(string):
    newString = ""
    complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
    for letter in string:
        newString += complement[letter.upper()]

    return newString[::-1]  

Comparison_Directory = input("Please enter Comparison List directory\n")
Database_Directory = input("Please enter Putative Site database directory\n")
Thresholds_Directory = input("Please enter thresholds directory\n")
start = time.time()

Promoter_List = []
with open(Comparison_Directory) as Promoter_Data:
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

thresholdDict = {}
with open(Thresholds_Directory) as Threshold_Data:
    for line in Threshold_Data:
        lineListed = line.split()
        if lineListed[-1].isnumeric():
            thresholdDict[lineListed[0]] = float(lineListed[-1])
        else:
            continue

Results = open("Results.csv", "w")
write_sites = csv.writer(Results, delimiter=";")
Result_Dictionary = {}

EnzymeResults = open("ReferenceEnzymeWideAnalysis.csv", "w")
EnzymeResult_Writer = csv.writer(EnzymeResults, delimiter=";")

PathwayResults = open("ReferencePathwayWideAnalysis.csv", "w")
PathwayResult_Writer = csv.writer(PathwayResults, delimiter=";")

for Item in Promoter_List:
    Promoter = Promoter_Sequence()
    Promoter.Cell = Item[0]
    Promoter.Pathway = Item[1]
    Promoter.Reaction = Item[2]
    Promoter.Enzyme = Item[3]
    Promoter.ID = Item[4]
    Promoter.Sequence = Item[5].lower()
    Promoter.Hash = Item[6]
    with open(Database_Directory) as Database:
        for Motif in motifs.parse(Database, "jaspar"):
            PWM = Motif.counts.normalize(pseudocounts=Pseudocount_Adder(Motif.counts))
            PSSM = PWM.log_odds(background={"A":0.31,"C":0.19, "G":0.19, "T":0.31})
            CORE = core_Pssm(Motif.counts, background={"A":0.31,"C":0.19, "G":0.19, "T":0.31})
            SCAN = PSSM.search(sequence=Promoter.Sequence, threshold=0.0, both=True) 
            for HIT in SCAN:
                try:
                    threshold = thresholdDict[Motif.base_id]
                except KeyError:
                    threshold = 0.75
                if Relative_Score(HIT[1], PSSM.max, PSSM.min) > threshold:
                    sequence = Promoter.Sequence[HIT[0]:HIT[0]+Motif.length]
                    hitSequence = reverseComplement(sequence) if HIT[0] < 0 else sequence
                    coreSCAN = CORE.search(sequence=hitSequence, threshold=0.0, both=False)
                    coreHITS = list(coreSCAN)
                    if coreHITS == []:
                        continue;
                    elif Relative_Score(coreHITS[0][1], CORE.max, 0) > 0.9:
                        Promoter.Sequence = Promoter.Sequence[:HIT[0]] + Promoter.Sequence[HIT[0]:HIT[0]+Motif.length].upper() + Promoter.Sequence[HIT[0]+Motif.length:]
                        Promoter.Putative_Site.append([Motif.name , HIT[0] , HIT[1] , Relative_Score(HIT[1], PSSM.max, PSSM.min) , Relative_Score(coreHITS[0][1], CORE.max, 0) , Promoter.Sequence[HIT[0]:HIT[0]+Motif.length]])
    
    write_sites.writerow(["Cell ID" , "Pathway", "Reaction", "Enzyme", "Sequence ID" , "Updated Sequence" , "Number of Putative Sites", "Hash"])
    write_sites.writerow([Promoter.Cell, Promoter.Pathway, Promoter.Reaction, Promoter.Enzyme, Promoter.ID , Promoter.Sequence , len(Promoter.Putative_Site), Promoter.Hash])
    write_sites.writerow(["HIT ID" , "HIT Location" , "HIT Score" , "HIT Relative Score" , "HIT Core Similarity" , "HIT Identity"])
    for Item in Promoter.Putative_Site:
        write_sites.writerow(Item)
    write_sites.writerow(["NULL"])  

    Resulting_Object = Promoter.write(Database_Directory)
    if Promoter.ID not in Result_Dictionary.keys():
        Result_Dictionary[Promoter.ID] = Resulting_Object
    else:
        Result_Dictionary[Promoter.ID]["Sites"].append(Resulting_Object["Sites"])
        Result_Dictionary[Promoter.ID]["Hash"] = list(Result_Dictionary[Promoter.ID]["Hash"]).append(Resulting_Object["Hash"])
        Result_Dictionary[Promoter.ID]["Site Counts"] = {k: Result_Dictionary[Promoter.ID]["Site Counts"].get(k, 0) + Resulting_Object["Site Counts"].get(k, 0) for k in set(Result_Dictionary[Promoter.ID]["Site Counts"]) | set(Resulting_Object["Site Counts"])}
        Result_Dictionary[Promoter.ID]["Site Index"] = sum(Result_Dictionary[Promoter.ID]["Site Counts"].values())/len(Result_Dictionary[Promoter.ID]["Promoter"])
  
Pathway_Dictionary = read_pathway(Comparison_Directory)

for Pathway_Key in Pathway_Dictionary.keys():
    for Enzyme_Key in Pathway_Dictionary[Pathway_Key].keys():
        for SubUnit_Key in Pathway_Dictionary[Pathway_Key][Enzyme_Key].keys():
                Pathway_Dictionary[Pathway_Key][Enzyme_Key][SubUnit_Key] = Result_Dictionary[SubUnit_Key]
                
for Pathway_Key in Pathway_Dictionary.keys():
    Pathway_Sites_Dict = {}
    Pathway_Common_Keys = []
    for Enzyme_Key in Pathway_Dictionary[Pathway_Key].keys():
        Enzyme_Sites_Dict = {}
        Enzyme_Common_Keys = []
        for SubUnit_Key in Pathway_Dictionary[Pathway_Key][Enzyme_Key].keys():
            Pathway_Sites_Dict[SubUnit_Key] = (Pathway_Dictionary[Pathway_Key][Enzyme_Key][SubUnit_Key]["Site Counts"])
            Enzyme_Sites_Dict[SubUnit_Key] = (Pathway_Dictionary[Pathway_Key][Enzyme_Key][SubUnit_Key]["Site Counts"])
        Enzyme_Key_Frequency = measure_key_frequency(Enzyme_Sites_Dict)
        for Key in Enzyme_Key_Frequency:
            if Enzyme_Key_Frequency[Key] >= 0.15*len(Enzyme_Sites_Dict):
                Enzyme_Common_Keys.append([Key , Enzyme_Key_Frequency[Key]/len(Enzyme_Sites_Dict)])
        EnzymeResult_Writer.writerow([Pathway_Key, Enzyme_Key, Enzyme_Common_Keys])

    Pathway_Key_Frequency = measure_key_frequency(Pathway_Sites_Dict)
    for Key in Pathway_Key_Frequency:
        if Pathway_Key_Frequency[Key] >= 0.15*len(Pathway_Sites_Dict):
            Pathway_Common_Keys.append([Key, Pathway_Key_Frequency[Key]/len(Pathway_Sites_Dict)])
    PathwayResult_Writer.writerow([Pathway_Key, Pathway_Common_Keys])

end = time.time()
print("Run Time: %.3f" %(end-start))
