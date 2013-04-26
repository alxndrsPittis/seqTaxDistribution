import sys
import os

from sys import exit
from collections import defaultdict
from string import strip
import sqlite3
import numpy
from pylab import *
from scipy.stats import gaussian_kde

from ete2 import PhyloTree, TreeStyle, NodeStyle, faces, parser
from ncbi_taxonomy.ncbi_query import get_topology, get_taxid_translator
from parser import blastm8
from colors.taxonColors import group2color
from plotting.ggplot import *

module_path = os.path.split(os.path.realpath(__file__))[0]

# DEFAULT THRESHOLDS
THRESHOLD_EVALUE = 1e-3
THRESHOLD_COVERAGE = 0

DATABASE = "uniprot"

#infile_blast = sys.argv[1]
#infile_fasta = sys.argv[2]

# Loading databases

db = sqlite3.connect(os.path.join(module_path, 'mapping/%s/idmapping.sqlite' % DATABASE))

def get_query_lengths(fasta_file):
    fasta = parser.fasta.read_fasta(fasta_file)
    id2length = dict( [ [fields[0], len(fields[1])] for fields in fasta ] )
    return id2length

def get_coverage_simple(start, end, qlength):
    coverage = float(end - start +1) / qlength
    return coverage

def get_ids_taxid(ids):
    filter = ','.join(map(lambda id: '"%s"' %id, [id for id in ids]))
    results = db.execute('select ac,taxon from idmapping where ac IN (%s)' %filter).fetchall()
    return dict(results)

def collapse_subspecies(t, all_taxids):
    species_nodes = [n for n in t.traverse() if n.rank == "species"
                     if int(n.name) in all_taxids]
    for sp_node in species_nodes:
        # Creation of new species node
        lower = sp_node.get_descendants()
        if lower:
            new_node = sp_node.__class__()
            for f in sp_node.features:
                new_node.add_feature(f, getattr(sp_node, f))
            # collapse species and all subcategories to same level    
            for n in lower:
                n.detach()
                n.name = n.name
                sp_node.add_child(n)
            sp_node.add_child(new_node)

def annotate_tree_with_taxonomy(t):
    taxid2name = get_taxid_translator([n.name for n in t.traverse()])
    for n in t.traverse():
        n.add_features(taxid=n.name)
        n.add_features(sci_name= str(taxid2name.get(int(n.name), "Not found")))
        n.name = "%s - %s" % (n.sci_name, n.name)
        if str(n.taxid) in group2color:
            n.add_features(bgcolor=group2color[str(n.taxid)])

def annotate_tree_with_hits(t, query2hits, id2taxid):
    # Create a dictionary with per taxid, hit information
    taxid2hits = defaultdict(list)
    for query in query2hits:
        for hit in query2hits[query]:
            if DATABASE == "uniprot":
                ac = hit.split("|")[1]
                id = hit.split("|")[2]            
                if ac in id2taxid:
                    taxid = id2taxid[ac]
                    entry  = [query, ac]
                    entry.extend(query2hits[query][hit])
                    taxid2hits[taxid].append(entry)    
    # Assigns a list of search results to the tree leaves - species
    for n in t.traverse():
        if int(n.taxid) in taxid2hits:
            results = []
            for hit in taxid2hits[int(n.taxid)]:
                results.append(hit)
            n.add_features(hits=results)

def get_search_values(node, q2length, threshold_evalue=1e-3, threshold_score=0, threshold_coverage=0):
    identities = []
    evalues = []
    scores = []
    coverages = []
    queries = []
    hits = []    
    for n in node.traverse():
        if "hits" in n.features:
            for entry in n.hits:
                identity = entry[2]
                evalue = entry[10]
                score = entry[11]
                query = entry[0]
                hit = entry[1]
                coverage = get_coverage_simple(entry[6], entry[7], q2length[query])
                if coverage >= threshold_coverage and evalue <= threshold_evalue and score >= threshold_score:
                    identities.append(identity)
                    evalues.append(evalue)
                    scores.append(score)
                    coverages.append(coverage)
                    queries.append(query)
                    hits.append(hit)

    return identities, evalues, scores, coverages, queries, hits

def sub_hist(subsets, labels, colors):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    data = []
    for n, sub in enumerate(subsets):
        values = numpy.array(sub)
        minValue = min(values)
        maxValue = max(values)
        data.append(values)
        # kde density estimation
        density = gaussian_kde(values)
        xs = numpy.arange(minValue-(minValue/10.),maxValue+(maxValue/10.),1)
        density.covariance_factor = lambda : .25
        density._compute_covariance()
        ys = density(xs)
        
        ax.plot(xs, ys, antialiased=True, linewidth=1, color="black", alpha=.5)#colors[n])
        ax.fill_between(xs, ys, alpha=.2, zorder=5, antialiased=True, color=colors[n])

    rhist(ax, data, normed=True, label=labels, color=colors, alpha=0.75, linewidth=0.5)

    ax.legend()
    rstyle(ax)
    #ax.invert_xaxis()
    #plt.show()
    plt.savefig("test.png")

def plot_subsets(t, q2length, taxids2plot=["131567"], attr="identities"):
    subsets = []
    names = []
    colors = []
    for taxid in taxids2plot:
        try:
            node = t.search_nodes(taxid=taxid)[0]
        except IndexError:
            print "warn--No hits in %s" % taxid
        identities, evalues, scores, coverages, queries, hits = get_search_values(node, q2length)
        attributes = {
            "identities" : identities,
            "evalues" : evalues,
            "scores" : scores,
            "coverages" : coverages
        }
        subset = attributes[attr]
        subsets.append(subset)
        name = node.name
        names.append(name)
        try:
            color = group2color[str(node.taxid)]
        except KeyError:
            color = group2color["Others"]
        colors.append(color)
    sub_hist(subsets, names, colors)

    #new_taxids = raw_input("new taxids selection (comma separated) ? \n")
    #if new_taxids:
    #    taxids2plot = map(strip, new_taxids.split(","))
    #    plot_subsets(t, q2length, taxids2plot)
        

def get_blast_tree(infile_blast):
    # Loads blast information into a dictionary    
    blast = blastm8.read_m8(infile_blast)

    # Hits and taxids retrieval
    hits = [hit.split("|")[1] for query in blast
            for hit in blast[query]]
    id2taxid = get_ids_taxid(hits)
    all_taxids = set(id2taxid.values())
    print "info--%s hits retrieved in %s taxids" % ( len(hits), len(all_taxids) )

    # Reconstruction of the annotated NCBI tree object
    t = get_topology(all_taxids, intermediate_nodes=False, rank_limit=None)
    collapse_subspecies(t, all_taxids)
    annotate_tree_with_taxonomy(t)
    annotate_tree_with_hits(t, blast, id2taxid)
    return t

    
# t = get_blast_tree(infile_blast) #for getting the annotated NCBI tree
# q2length = get_query_lengths(infile_fasta) #for storing in dictionary the lengths
# plot_subsets(t, q2length, ["taxid1", "taxid2"]) #for plotting test.png