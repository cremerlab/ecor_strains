import pandas as pd
import os
from Bio import SeqIO
import numpy as np
#############
#define some simple sequencing analysis tools which we repeatedly need
#############
def get_sequences(filenamein):
    if os.path.exists(filenamein) == False:
        print("sequence file not found: "+filenamein)
    #Input: path of fasta file with sequences
    #Output: list of sequences
    fasta_sequences = SeqIO.parse(open(filenamein),'fasta')
    contigname=[]
    contigsequence=[]
    contigseqlength=[]
    for fasta in fasta_sequences:
        name, sequence = fasta.id, str(fasta.seq)
        contigname.append(name)
        contigsequence.append(sequence)
        contigseqlength.append(len(sequence))
    return [contigname,contigsequence,np.array(contigseqlength)]

#get specific subsection of sequence
def get_sequences_region(filenamein,contig,posx,posy):
    #prefent problems in case pos arguments are floats
    posx=int(posx)
    posy=int(posy)
    #load file
    if os.path.exists(filenamein) == False:
        print("sequence file not found: "+filenamein)
        error
    fasta_sequences = SeqIO.parse(open(filenamein),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq
        #print(name)
        #print(sequence)
        if name==contig:
            sequencecur=sequence
            break
    if posx<posy:
        return(sequencecur[posx:posy-1])
    elif posx>posy:
        return(sequencecur[posy:posx-1])
    else:
        error
        
        
#get section of sequence from provided start to the end of the contig:
def get_sequence_to_end(filenamein,contig,posx):
    #prevent problems in case pos arguments are floats
    posx=int(posx)
    #load file
    if os.path.exists(filenamein) == False:
        print("sequence file not found: "+filenamein)
        error
    fasta_sequences = SeqIO.parse(open(filenamein),'fasta')
    sequencecur = ""
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq
        if name==contig:
            sequencecur=sequence
            break
    if sequencecur == "":
        print("contig not found in " + filenamein)
        error
    return(sequencecur[posx:len(sequencecur)])

def get_contig_length(filenamein,contig):
    #load file
    if os.path.exists(filenamein) == False:
        print("sequence file not found: "+filenamein)
        error
    fasta_sequences = SeqIO.parse(open(filenamein),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq
        if name==contig:
            sequencecur=sequence
            break
    return(len(sequencecur))

def get_simple_annotation(row,outputoption="homologuesource"):
    infc=str(row['inference'])
    if "similar to AA sequence:MG1655.gbff" in infc:
        infout="MG1655ref"
        homologinfo=infc.split("similar to AA sequence:MG1655.gbff")[1]
    elif "similar to AA sequence:UniProtKB" in infc:
        infout="UNIPROT"
        homologinfo=infc.split("similar to AA sequence:UniProtKB")[1]
    elif "ab initio prediction:Prodigal" in infc and len(infc) == len("ab initio prediction:Prodigal:002006"):
        infout="no_homologue"
        homologinfo=np.nan
    elif "protein motif:HAMAP:" in infc:
        infout="HAMAP"
        homologinfo=infc.split("protein motif:HAMAP:")[1]
    elif "ISfinder" in infc:
        infout="ISfinder"
        homologinfo=infc.split("ISfinder")[1]
    else:
        infout=infc
        homologinfo=np.nan
    if outputoption=="homologueinfo":
        return homologinfo
    elif outputoption=="homologuesource":
        return infout
    else:
        pass