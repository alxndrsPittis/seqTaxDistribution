import os
from string import strip
from sys import stderr as STDERR

def read_m8(source, cog=True, multiple_hits_per_species=True, multiple_hits_per_sequence=False):
    """ Reads a collection of blast results encoded in m8 format."""

    blast = {}
    # Prepares handle for blast m8 file
    if os.path.isfile(source):
        _source = open(source, "rU")
    else:
        _source = iter(source.split("\n"))

    for line in _source:
        line = line.strip()
        # Blast -m8 fields
        blast_fields = map(strip, line.split())
        query_name = blast_fields[0]
        hit_name = blast_fields[1]
        percent_identity = float(blast_fields[2])
        aln_length = int(blast_fields[3])
        mismatches = int(blast_fields[4])
        gaps = int(blast_fields[5])
        query_start = int(blast_fields[6])
        query_end = int(blast_fields[7])
        hit_start = int(blast_fields[8])
        hit_end = int(blast_fields[9])
        evalue = float(blast_fields[10])
        score = float(blast_fields[11])
        fields = [percent_identity,
                  aln_length,
                  mismatches,
                  gaps,
                  query_start,
                  query_end,
                  hit_start,
                  hit_end,
                  evalue,
                  score]

        if cog:
            #Store to dictionary all m8 fields
            if multiple_hits_per_species and not multiple_hits_per_sequence:                
                blast.setdefault(query_name, {})
                blast[query_name].setdefault(hit_name, fields)
                
        else:
            pass
            
    return blast
        

#blast = read_m8("/home/apittis/Documents/BioProjects/Origins/Fusion/Benchmark/NCBIOrganelles/Tests/6239.NP_006959.blast")