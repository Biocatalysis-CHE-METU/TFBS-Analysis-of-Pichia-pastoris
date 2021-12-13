import Bio.pairwise2 as pairwise2
from Bio.pairwise2 import format_alignment
import Bio.motifs as motifs
import csv
import time
import re

class SubUnit_Object(object): 
    def __init__(self, Cell_ID, Pathway, Reaction, Enzyme, SubUnit, Promoter_Sequence, Hash):
        self.Cell_Line = Cell_ID #The cell line this enzyme is found at. 
        self.Pathway_ID = Pathway #Identity of the pathway this subunit serves.
        self.Reaction_ID = Reaction #Identity of the reaction this subunit serves.
        self.Enzyme_ID = Enzyme #Identity of the enzyme that contains this subunit. 
        self.SubUnit_ID = SubUnit #ID of this subunit.
        self.Sequence = Promoter_Sequence #Promoter sequence of this subunit.
        self.Putative_Site = [] #Putative sites that has an impact on this subunit. In form of a list of dictionaries.
        self.Hash = Hash

    def read_site(self, Site_ID, Site_Location, Site_Sequence, Site_RC):
        #Takes input from file to create a list of dictionaries.                 
        self.Putative_Site.append({"Site_ID": Site_ID, "Start_Location": int(Site_Location),
                                        "End_Location": int(Site_Location) + len(Site_Sequence),
                                        "Site_Sequence": Site_Sequence, "Reverse_Complement": Site_RC})
    
    def align(self, Compared_Object, start=0, end=0, OpG_Penalty=-10, ExG_Penalty=-1, EndG_Penalty = False): 
        #Align a given subunit or putative site to a part of this object. 
        end = len(self.Sequence) if end == 0 else end
        Interest_Sq = self.Sequence[start:end].upper() #Sequence that we are interested to check. Is equal to entire sequence if no start and end location is given.
        Compared_Sq = Compared_Object if type(Compared_Object) == str else Compared_Object.Sequence.upper() #Sequence that we comparing.
        if Interest_Sq == "":
            return[["-","-",0]]
        Needleman = pairwise2.align.globalms(Interest_Sq, Compared_Sq, 5, -4, OpG_Penalty, ExG_Penalty, penalize_end_gaps = EndG_Penalty)
        return Needleman

    def find_site(self, Reference_Sq, Known_Coordinates, Index=0, Slash_Count=0, Is_End=0, Aligned=0):       
        Updated_Start = 0
        Updated_End = 0
        
        if Aligned == 0:
            #Count number of slashes before the putative site. Do not count slashes as characters so that the final value of the while loop is 
            # equal to initial location of the site.
            Known_Location = Known_Coordinates[1] if Is_End == 1 else Known_Coordinates[0] + 1
            while Index < Known_Location:
                if Reference_Sq[Index + Slash_Count] == "-":
                    Slash_Count = Slash_Count + 1
                else:
                    Index = Index + 1
            Updated_Start = Known_Location + Slash_Count if Is_End == 1 else Known_Location + Slash_Count - 1
            #After finding the initial location, call the function recursively but this time to find the new end location. 
            if Is_End == 0:
                Updated_End = self.find_site(Reference_Sq, Known_Coordinates, Index, Slash_Count, Is_End=1)
            if Is_End == 1:
                return Updated_Start   
        else:
            #If transferring from aligned to actual, remove slashes from before and within the sequence. 
            Slash_Before_Start = Reference_Sq[0:Known_Coordinates[0]+1].count("-")
            Slash_Before_End = Reference_Sq[0:Known_Coordinates[1]].count("-")
            Shift = Slash_Before_End - Slash_Before_Start #If any gaps where removed from within the corresponding region, pull 1 char from both sides to make up for it. 
            if Shift == Known_Coordinates[1] - Known_Coordinates[0] - 1:  #If there are no matching sites corresponding to the sequence just return an emptry string location.
                return [0, 0]
            Updated_Start = Known_Coordinates[0] - Slash_Before_Start - Shift
            Updated_End = Known_Coordinates[1] - Slash_Before_End + Shift 
        return [Updated_Start, Updated_End]

    def refine(self, ComparisonList):
        RefinedSites = []
        for Comparison in ComparisonList:
            string = Comparison[1]
            slashes = re.finditer("-", string)
            match_pos = [slash.start() for slash in slashes]
            match_pos.insert(0, -1)
            match_pos.append(len(string)+1)

            for ind in range(len(match_pos)):
                if match_pos[ind+1] - match_pos[ind] > 1:
                    start_str, end_str = match_pos[ind], match_pos[ind+1]
                    break;
            RefinedSites.append(Comparison[0][start_str+1:end_str])    
        return RefinedSites

    def compare(self, Compared_Object):
        Alignment = self.align(Compared_Object)
        Alignment = Alignment[0]
        for Item in Compared_Object.Putative_Site:
            Aligned_Site_Locations = self.find_site(Alignment[1] , [Item["Start_Location"], Item["End_Location"]])
            Translated_Site_Locations = self.find_site(Alignment[0], Aligned_Site_Locations, Aligned = 1)
            if Translated_Site_Locations == [0,0]:
                continue;
            Hit_Achieved_Zone = self.Sequence[Translated_Site_Locations[0]:Translated_Site_Locations[1]]
            if len(Hit_Achieved_Zone) == len(Item["Site_Sequence"]):
                Comparison = self.align(Item["Site_Sequence"], Translated_Site_Locations[0], Translated_Site_Locations[1], OpG_Penalty=-100, EndG_Penalty = True)
            else:
                Comparison = self.align(Item["Site_Sequence"], Translated_Site_Locations[0], Translated_Site_Locations[1], OpG_Penalty=-100)
            if Comparison[0][2] >= len(Item["Site_Sequence"])*2 and "-" not in Comparison[0][0]:
                if len(Item["Site_Sequence"]) < len(Hit_Achieved_Zone):
                    Hit_Achieved_Zone = self.refine(Comparison)
                self.Putative_Site.append({"Site_ID": Item["Site_ID"], "Start_Location": Translated_Site_Locations[0],
                                                "End_Location": Translated_Site_Locations[1], "Site_Sequence": Hit_Achieved_Zone,
                                                "Motif": Item["Site_Sequence"], "Reverse_Complement": Item["Reverse_Complement"]})    

    def count_site(self, Database_Directory):
        Site_Counts = {}
        with open(Database_Directory) as Database:
            for Motif in motifs.parse(Database, "jaspar"):
                Count = 0
                for Sites in self.Putative_Site:
                    if Motif.name == Sites["Site_ID"]:
                        Count = Count + 1
                if Count == 0:
                    continue;
                Site_Counts[Motif.name] = Count
        return Site_Counts
    
    def write(self, Database_Directory):
        Site_Counts = self.count_site(Database_Directory)
        Total_Count = sum(Site_Counts.values())
        Record = {"Cell ID": self.Cell_Line, "Pathway": self.Pathway_ID, "Enzyme": self.Enzyme_ID, 
                  "Reaction": self.Reaction_ID, "SubUnit": self.SubUnit_ID, "Promoter": self.Sequence, 
                  "Sites": self.Putative_Site, "Site Counts": Site_Counts, "Site Index": Total_Count/len(self.Sequence), "Hash" : self.Hash}
        return Record

def read_csv(csv_directory, File_Type):
    with open(csv_directory) as csv_file:
        Mode = ""
        if File_Type == "Reference":
            Reference_Dictionary_List = []
            Temporary_Dictionary = {}
            Temporary_Site_List = []
            Reader = csv.reader(csv_file, delimiter = ";")
            for Row in Reader:
                #Read the current program and adjust what program is doing with respect to what you read.
                if Row == []:
                    continue;
                Mode = ("Read Sequence" if Row[0] == "Cell ID" else 
                         ("Read Factor" if Row[0] == "HIT ID" else Mode))
                if (Row[0] == "Enzyme ID") or (Row[0] == "HIT ID"):
                    continue;
                elif Row[0] == "NULL":
                    Temporary_Dictionary["Sites"] = Temporary_Site_List
                    Reference_Dictionary_List.append({Key : Temporary_Dictionary[Key] for Key in Temporary_Dictionary.keys()})
                    Temporary_Dictionary = {}
                    Temporary_Site_List = []
                    continue;
                #If the mode of the program is set, and its not a mode changing location, read the data.
                if Mode == "Read Sequence":
                    Temporary_Dictionary = {"Cell ID": Row[0], "Pathway": Row[1], "Reaction": Row[2],
                                            "Enzyme": Row[3], "SubUnit": Row[4], "Promoter": Row[5], "Hash": Row[7]}
                elif Mode == "Read Factor":
                    factorPosition = int(Row[1]) if int(Row[1]) >= 0 else len(Temporary_Dictionary["Promoter"]) + int(Row[1])
                    reverseComplement = 1 if int(Row[1]) <= 0 else 0
                    Temporary_Site_List.append([Row[0], factorPosition, Row[5], reverseComplement])
            return Reference_Dictionary_List
        elif File_Type == "Comparison":
            Reader = csv.DictReader(csv_file, delimiter = ";")
            Compared_Dictionary_List = []
            Temporary_Dictionary = {"Cell ID": ""}   
            Hash_Tracker = []  
            for Row in Reader:
                if Row["Hash"] in Hash_Tracker:
                    continue;
                Hash_Tracker.append(Row["Hash"])
                if Row["Compared Gene Name"] == "":
                    continue;
                Temporary_Dictionary["Cell ID"] = Row["Compared Cell Line"] if Temporary_Dictionary["Cell ID"] == "" else Temporary_Dictionary["Cell ID"] 
                Temporary_Dictionary["SubUnit"] = Row["Compared Gene Name"]
                Temporary_Dictionary["Promoter"] = Row["Compared Sequence"]
                Temporary_Dictionary["Hash"] = Row["Hash"]
                if Row["Pathway"] != "":
                    Temporary_Dictionary["Pathway"] = Row["Pathway"]
                if Row["Reaction"] != "":
                    Temporary_Dictionary["Reaction"] = Row["Reaction"]
                if Row["Compared Enzyme"] != "":
                    Temporary_Dictionary["Enzyme"] = Row["Compared Enzyme"]
    
                Compared_Dictionary_List.append({Key : Temporary_Dictionary[Key] for Key in Temporary_Dictionary.keys()})
            return Compared_Dictionary_List 
        elif File_Type == "Pathway":
            Pathway_Dict = {}
            Reader = csv.DictReader(csv_file, delimiter = ";")
            for Row in Reader:
                if Row["Pathway"] != "":
                    Pathway_Dict[Row["Pathway"]] = {}
                else:
                    Row["Pathway"] = list(Pathway_Dict.keys())[-1] 
                if Row["Compared Enzyme"] != "":       
                    Pathway_Dict[Row["Pathway"]][Row["Compared Enzyme"]] = {}
                else:
                    Row["Compared Enzyme"] = list(Pathway_Dict[Row["Pathway"]].keys())[-1]
                Pathway_Dict[Row["Pathway"]][Row["Compared Enzyme"]][Row["Compared Gene Name"]] = {}
            return Pathway_Dict        

def measure_key_frequency(dict_of_dicts):
    frequency = {}
    for d in dict_of_dicts.values():
        for key in d.keys():
            frequency[key] = 1 if key not in frequency else frequency[key] + 1
    return frequency

def reverseComplement(sequence):
    complement = {"A": "T", "T": "A", "C": "G", "G":"C"}
    reverseSq = sequence[::-1]
    reverseComplementSq = ""
    for char in reverseSq:
        char = char.upper()
        reverseComplementSq += complement[char]
    return reverseComplementSq

if __name__ == "__main__":    
    Comparison_Directory = input("Please enter Comparison List directory")
    Reference_Directory = input("Please enter Scanner Results' directory")
    Database_Directory = input("Please enter Putative Site database directory")
    start = time.time()

    Reference_Sequence_Dictionary = read_csv(Reference_Directory, "Reference")
    Compared_Sequence_Dictionary = read_csv(Comparison_Directory, "Comparison")
    Pathway_Dictionary = read_csv(Comparison_Directory, "Pathway")
    
    Result_Dictionary = {}
    SubResults = open("SubunitWideAnalysis.csv", "w")
    SubResult_Writer = csv.writer(SubResults, delimiter = ";")
    Enzyme_Results = open("EnzymeWideAnalysis.csv", "w")
    EnzymeResult_Writer = csv.writer(Enzyme_Results, delimiter = ";")
    Pathway_Results = open("PathwayWideAnalysis.csv", "w")
    PathwayResult_Writer = csv.writer(Pathway_Results, delimiter = ";")
 
    for SUBUNIT in Compared_Sequence_Dictionary:
        Compared_Sequence = SubUnit_Object(SUBUNIT["Cell ID"], SUBUNIT["Pathway"], SUBUNIT["Reaction"],
                                           SUBUNIT["Enzyme"], SUBUNIT["SubUnit"], SUBUNIT["Promoter"], SUBUNIT["Hash"])
        for REFERENCE in Reference_Sequence_Dictionary:
            if SUBUNIT["Hash"] == REFERENCE["Hash"]:
                Reference_Sequence = SubUnit_Object(REFERENCE["Cell ID"], REFERENCE["Pathway"], REFERENCE["Reaction"],
                                                    REFERENCE["Enzyme"], REFERENCE["SubUnit"], REFERENCE["Promoter"], REFERENCE["Hash"])
                for Site in REFERENCE["Sites"]:        
                    Reference_Sequence.read_site(Site[0],Site[1],Site[2], Site[3])
                Compared_Sequence.compare(Reference_Sequence)
        Resulting_Object = Compared_Sequence.write(Database_Directory)
        if SUBUNIT["SubUnit"] in Result_Dictionary.keys():
            if type(Result_Dictionary[SUBUNIT["SubUnit"]]["Hash"]) != list:
                print("{}:\t{}".format(Result_Dictionary[SUBUNIT["SubUnit"]]["Hash"], Result_Dictionary[SUBUNIT["SubUnit"]]["Site Counts"]))
                Result_Dictionary[SUBUNIT["SubUnit"]]["Hash"] = [Result_Dictionary[SUBUNIT["SubUnit"]]["Hash"]]
            print("{}:\t{}".format(Resulting_Object["Hash"], Resulting_Object["Site Counts"]))
            Result_Dictionary[SUBUNIT["SubUnit"]]["Sites"] = Result_Dictionary[SUBUNIT["SubUnit"]]["Sites"] + Resulting_Object["Sites"]
            Result_Dictionary[SUBUNIT["SubUnit"]]["Hash"].append(Resulting_Object["Hash"])
            Result_Dictionary[SUBUNIT["SubUnit"]]["Site Counts"] = {k: Result_Dictionary[SUBUNIT["SubUnit"]]["Site Counts"].get(k, 0) + Resulting_Object["Site Counts"].get(k, 0) for k in set(Result_Dictionary[SUBUNIT["SubUnit"]]["Site Counts"]) | set(Resulting_Object["Site Counts"])}
            Result_Dictionary[SUBUNIT["SubUnit"]]["Site Index"] = sum(Result_Dictionary[SUBUNIT["SubUnit"]]["Site Counts"].values())/len(Result_Dictionary[SUBUNIT["SubUnit"]]["Promoter"])
        else:
            Result_Dictionary[SUBUNIT["SubUnit"]] = Resulting_Object

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
                if Enzyme_Key_Frequency[Key] >= len(Enzyme_Sites_Dict)*0.15:
                    Enzyme_Common_Keys.append([Key , Enzyme_Key_Frequency[Key]/len(Enzyme_Sites_Dict)])
            EnzymeResult_Writer.writerow([Pathway_Key, Enzyme_Key, Enzyme_Common_Keys])

        Pathway_Key_Frequency = measure_key_frequency(Pathway_Sites_Dict)
        for Key in Pathway_Key_Frequency:
            if Pathway_Key_Frequency[Key] >= len(Pathway_Sites_Dict)*0.15:
                Pathway_Common_Keys.append([Key, Pathway_Key_Frequency[Key]/len(Pathway_Sites_Dict)])
        PathwayResult_Writer.writerow([Pathway_Key, Pathway_Common_Keys])

    tf_Dict = {}
    for Key in Result_Dictionary.values():
        SubResult_Writer.writerow(list(Key.values()))
        for Site in Key["Sites"]:
            if Site["Site_ID"] in tf_Dict.keys():
                if type(Site["Site_Sequence"]) == str:
                    correctedSq = reverseComplement(Site["Site_Sequence"]) if Site["Reverse_Complement"] == 1 else Site["Site_Sequence"]
                    tf_Dict[Site["Site_ID"]].append(correctedSq)
                else:
                    correctedSq = [reverseComplement(Sq) for Sq in Site["Site_Sequence"]] if Site["Reverse_Complement"] == 1 else Site["Site_Sequence"]
                    tf_Dict[Site["Site_ID"]] = tf_Dict[Site["Site_ID"]] + correctedSq
            else:
                Site_Sequence_List = [Site["Site_Sequence"]] if type(Site["Site_Sequence"]) == str else Site["Site_Sequence"]
                tf_Dict[Site["Site_ID"]] = [reverseComplement(Sq) for Sq in Site_Sequence_List] if Site["Reverse_Complement"] == 1 else Site_Sequence_List
 
    tf_Matrix_Dict = {}
    for Key in tf_Dict.keys():
        seqList = tf_Dict[Key]
        fm = []
        for i in range(len(seqList[0])):
            fm.append({'A':0, 'C':0, 'T':0, 'G':0})
            for site in seqList:
                site = site.upper()
                fm[i][site[i]] = fm[i][site[i]] + 1
        tf_Matrix_Dict[Key] = fm

    with open("site.txt","w") as siteRes:
        for Key in tf_Matrix_Dict.keys():
            siteRes.write("> " +  Key + "\n")
            for base in ["A","C","G","T"]:
                siteRes.write(base + " [")
                for loc in tf_Matrix_Dict[Key]:
                    digit = loc[base]
                    siteRes.write("{:>6}".format(digit))
                siteRes.write(" ]\n")
            siteRes.write("\n")

    with open("targetTFs.txt","w") as targetTFs:
        for Key in ["ADR1", "CAT8", "SIP4", "HAP234", "RDS2", "YBR239C", 
                    "STB5", "MSN2", "MSN4", "MIG1", "TYE7", "GCR1"]:
            if Key in tf_Dict.keys():    
                targetTFs.write("{} \n".format(Key))
                for SiteSequence in tf_Dict[Key]:
                    targetTFs.write("{} \n".format(SiteSequence))
            else:
                continue

    end = time.time()
    print("Run Time: %.3f" %(end-start))
