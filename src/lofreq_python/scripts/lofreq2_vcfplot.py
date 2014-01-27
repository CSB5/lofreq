#!/usr/bin/env python
"""Plot characteristics of variants listed in VCF file
"""


__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2014 Genome Institute of Singapore"
__license__ = "GPL2"



# --- standard library imports
#

# imports
import sys
import os
import argparse
import logging
import gzip
from collections import Counter, deque
import itertools

#--- third-party imports
#
import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# only for boxplots
from scipy.stats import gaussian_kde


#--- project specific imports
#
try:
    import lofreq2_local
except ImportError:
    pass

try:
    from lofreq_star import vcf
except ImportError:
    sys.stderr.write("FATAL(%s): Couldn't find LoFreq's vcf module."
                     " Are you sure your PYTHONPATH is set correctly (= %s)?\n" % (
                         (sys.argv[0], os.environ['PYTHONPATH'])))
    sys.exit(1)
from lofreq_star.utils import complement, now


# invocation of ipython on exceptions
#import sys, pdb
#from IPython.core import ultratb
#sys.excepthook = ultratb.FormattedTB(mode='Verbose',
#                                     color_scheme='Linux', call_pdb=1)


# global logger
#
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


COLORS = ["b", "g", "r", "c", "m", "y", "k"]


def subst_type_str(ref, alt, strand_specific=False):
    """
    """
    
    # in case we get a list
    assert len(ref)==1 and len(alt)==1
    ref = ref[0]
    alt = alt[0]

    s = "%s>%s" % (ref, alt)
    if strand_specific:
        return s
    else:
        c = complement(s)
        return '|'.join(sorted([s, c]))
    


def r_ify(axes):
    '''
    source:
    http://stackoverflow.com/questions/14349055/making-matplotlib-graphs-look-like-r-by-default
    ttp://messymind.net/2012/07/making-matplotlib-look-like-ggplot/

    Produce R-style Axes properties
    '''
    xticks = axes.get_xticks()
    yticks = axes.get_yticks()

    #remove right and upper spines
    axes.spines['right'].set_color('none')
    axes.spines['top'].set_color('none')

    #make the background transparent
    axes.set_axis_bgcolor('none')

    #allow space between bottom and left spines and Axes
    axes.spines['bottom'].set_position(('axes', -0.05))
    axes.spines['left'].set_position(('axes', -0.05))

    #allow plot to extend beyond spines
    axes.spines['bottom'].set_bounds(xticks[0], xticks[-1])
    axes.spines['left'].set_bounds(yticks[0], yticks[-1])

    #set tick parameters to be more R-like
    axes.tick_params(direction='out', top=False, right=False, length=10, pad=12, width=1, labelsize='medium')

    #set x and y ticks to include all but the last tick
    axes.set_xticks(xticks[:-1])
    axes.set_yticks(yticks[:-1])

    return axes



def af_hist(ax, af_list, bins=20):
    ax.hist(af_list, bins=bins)
    ax.set_xlim((0, 1.0))
    ax.set_ylabel('#SNVs')
    ax.set_xlabel('AF')
    #r_ify(ax)
    # FIXME would be nice to 'clip' max


def cov_hist(ax, dp_list, bins=20):
    ax.hist(dp_list, bins=bins)
    ax.set_xlim([0, plt.xlim()[1]])
    ax.set_ylabel('#SNVs')
    ax.set_xlabel('DP')
    #r_ify(ax)

    
def subst_perc(ax, subst_type_counts):
    """
    subst_type_counts should be list of array with type as 1st element and count as 2nd
    """

    # FIXME sort by transition/transversion type. Add Ts/Tv ratio to plot
    
    #colors = [cm.jet(1.*i/len(subst_type_counts)) for i in xrange(len(subst_type_counts))]
    colors = [COLORS[i % len(COLORS)] for i in xrange(len(subst_type_counts))]

    count_sum = sum([x[1] for x in subst_type_counts])
    percs = [x[1]/float(count_sum) for x in subst_type_counts]
    ax.bar(xrange(len(subst_type_counts)), percs, color=colors)

    ticks = [x[0] for x in subst_type_counts]
    ax.set_xticks(xrange(len(ticks))) # forced display of all
    ax.set_xticklabels(ticks, rotation=45, ha="left")
    # FIXME rotation=45 doesnt't work
    # FIXME ha="left" doesn't work
    # ax1.set_xticks(ticks)
    # FIXME ticks as string won't work
    ax.set_ylabel('[%]')
    ax.set_xlabel('Type')
    
    # prevent clipping of tick-labels
    #plt.subplots_adjust(bottom=0.15)
    plt.tight_layout()

    
def calc_dist_slow(variants, na=None):
    """FIXME:add-doc
    """

    sys.stderr.write("Don't use me. I'm slow. Use calc_dist() instead\n")

    print "starting at %s" % now()
    dist_to_next = []
    for i in xrange(len(variants)):
        cur_var = variants[i]

        mod = len(variants)/5
        if i % mod == 0:
            print "Alive and munching away at var %d of %d" % (i, len(variants)) # FIXME to log

        prev_var_pos = None
        for j in reversed(range(i)):
            if variants[j].CHROM != cur_var.CHROM:
                break
            assert variants[j].POS <= cur_var.POS, ("vcf file needs to be sorted but isn't!")
            if variants[j].POS != cur_var.POS:
                prev_var_pos = variants[j].POS
                break

        next_var_pos = None
        for j in range(i+1, len(variants)):
            if variants[j].CHROM != cur_var.CHROM:
                break
            assert variants[j].POS >= cur_var.POS, ("vcf file needs to be sorted but isn't!")
            if variants[j].POS != cur_var.POS:
                next_vars_pos = variants[j].POS
                break

        min_dist = na
        if prev_var_pos != None:
            dist_prev = cur_var.POS - prev_var_pos
            min_dist = dist_prev
        if next_var_pos != None:
            dist_next = next_var_pos - cur_var.POS
            if prev_var_pos:
                min_dist = min([dist_prev, dist_next])
        dist_to_next.append(min_dist)

    print "end at %s" % now()
    return dist_to_next


def calc_dist(variants): 
    """Calculated distance to next variant. 

    If a chromosome only contains a single SNV, -1 will be stored as
    dist as we can't use 0 which would mean multi-allelic position.

    Variants need to be sorted (checking via assert here)

    This is several order of magnitudes faster then calc_dist_to_next
    """

    #print "starting at %s" % now()

    dists = []

    # group per chromosome
    processed_chroms = []
    for (chrom, vars_on_chrom) in itertools.groupby(variants, lambda v: v.CHROM):
        assert chrom not in processed_chroms
        processed_chroms.append(chrom)

        # use a queue. fill up with max 3 elements at a time. every
        # time we kick one out report the minimum dist between it and
        # the snv on the left and right (if any)

        deck = deque(itertools.islice(vars_on_chrom, 3))
        if len(deck) == 1:
            dists.append(-1)
            continue

        left_dist = sys.maxint
        for elem in vars_on_chrom:
            right_dist = deck[1].POS - deck[0].POS
            min_dist = min([left_dist, right_dist])
            dists.append(min_dist)
            #print "Popping %s %d with min_dist %d" % (
            #    deck[0].CHROM, deck[0].POS, min_dist)
            deck.popleft()
            deck.append(elem)
            left_dist = right_dist

        # dismantle. same as above without appending and left_dist
        # update
        while len(deck)>1:
            right_dist = deck[1].POS - deck[0].POS
            min_dist = min([left_dist, right_dist])
            dists.append(min_dist)
            #print "Popping %s %d with min_dist %d" % (
            #    deck[0].CHROM, deck[0].POS, min_dist)
            deck.popleft()
            left_dist = right_dist

        dists.append(left_dist)

    assert len(dists) == len(variants)
    #print "end at %s" % now()
    
    return dists


def dist_hist(ax, var_dist, bins=15):
    """FIXME
    """
    x = [np.log10(d) if d>0 else -1 for d in var_dist]
    ax.hist(x, bins=bins)
    #xlim([0, xlim()[1]])
    ax.set_ylabel('#')
    ax.set_xlabel("SNV distance (log10)")
    ax.set_title("(%d bins)" % bins)


def heatmap_dist_vs_x_1(ax, variants, var_dist, info_key, bins=20):
    """FIXME:add-doc

    info_key must be a valid vcf.INFO key, e.g. AF or DP
    """
    #from matplotlib.colors import LogNorm

    assert all([v.INFO.has_key(info_key) for v in variants]), (
        "Not all variants contain INFO key %s", info_key)
    x = [np.log10(d) if d>0 else -1 for d in var_dist]
    y = [v.INFO[info_key] for v in variants]
    ax.dist_vs_dp_plot = plt.hist2d(x, y, bins=bins)

    # FIXME add option to whiten missing values

    #print min(x), max(x)
    #print min(y), max(y)
    ax.set_ylabel(info_key)
    ax.set_xlabel("SNV distance (log10)")
    ax.set_ylim([0, plt.ylim()[1]])

    plt.colorbar()
    

def heatmap_dist_vs_cov_2(ax, variants, var_dists, bins=20):

    # numpy.histogram2d and then heatmap instead of hist2d:
    # http://stackoverflow.com/questions/2369492/generate-a-heatmap-in-matplotlib-using-a-scatter-data-set
    # http://stackoverflow.com/questions/16917836/matplotlib-stretches-histogram2d-vertically
    # smoothing: http://stackoverflow.com/questions/6652671/efficient-method-of-calculating-density-of-irregularly-spaced-points
    # still need to set bins but smoothing should iron out sub-optimal choice
    SMOOTH = False
    if SMOOTH:
        # FIXME
        import scipy.ndimage as ndi

    x = [np.log10(d) if d>0 else -1 for d in var_dists]
    y = [v.INFO['DP'] for v in variants]
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=bins)
    extent = [xedges.min(), xedges.max(), yedges.min(), yedges.max()]
    # http://stackoverflow.com/questions/8762610/heatmap-using-scatter-dataset-python-matplotlib
    data = heatmap.T
    if SMOOTH:
        gk_std = (np.std(x), np.std(y))
        gk_std = (np.std(y), np.std(x))
        #gk_std = (2, 2)
        #gk_std = (min([std(x), std(y)]))
        data = ndi.gaussian_filter(data, (gk_std))
        # small values work better for gaussian kernel stdev
        # swapped(!) std values give good results. not sure why.
        # any value between 1 and 2 worked fine as well

    im = ax.imshow(data, extent=extent,  origin='lower', aspect='auto',
        interpolation=None)#nearest|None|bicubic|bilinear
    ax.set_xlabel("SNV distance (log10)")
    ax.set_ylabel('DP')
    
    return im


def heatmap_freq_vs_cov(ax, variants, bins=50):
    x = [v.INFO['DP'] for v in variants]
    y = [v.INFO['AF'] for v in variants]

    # FIXME add option to whiten missing values
    ax.dist_vs_dp_plot = plt.hist2d(x, y, bins=bins)

    ax.set_xlabel('DP')
    ax.set_ylabel('AF')
    ax.set_xlim([0, ax.get_xlim()[1]])
    #ax1.set_ylim([0, 1])
    plt.colorbar()

    
def violin_plot(ax, data):
    '''
    Create violin plots on an axis

    from http://pyinsci.blogspot.sg/2009/09/violin-plot-with-matplotlib.html
    '''

    # FIXME possible that this needs values between 0 and 1?
    
    w = min(0.15, 0.5)
    k = gaussian_kde(data) # calculates the kernel density
    m = k.dataset.min() #lower bound of violin
    M = k.dataset.max() #upper bound of violin
    x = np.arange(m, M, (M-m)/100.) # support for violin
    v = k.evaluate(x) # violin profile (density curve)
    if v.max():
        v = v/v.max()*w # scaling the violin to the available space
    else:
        LOG.warn("v.max()==0. won't be able to print violin_plot")
        v = 0
    p = 0
    ax.fill_betweenx(x, p, v+p, facecolor='y', alpha=0.3)
    ax.fill_betweenx(x,p, -v+p, facecolor='y', alpha=0.3)
    l = w+w*0.1
    plt.xlim((-l, l))
    #print "DEBUG", w, k, m, M, x, v
    ax.set_xticks([])
    #ax1.set_xticklabels(ticks, rotation=45, ha="left")



def print_overview(ax, text_list):
    """FIXME:add-doc
    """

    # options:
    # - annotate () or text()
    # - tex or text
    
    # See http://jakevdp.github.io/mpl_tutorial/tutorial_pages/tut4.html
    
    #matplotlib.rc('text', usetex=True)
    #table = r'\begin{table} \begin{tabular}{|l|l|l|}  \hline  $\alpha$      & $\beta$        & $\gamma$      \\ \hline   32     & $\alpha$ & 123    \\ \hline   200 & 321    & 50 \\  \hline  \end{tabular} \end{table}'

    ax.axis('off')     
    ax.text(0, 0.8, '\n'.join(text_list), size=14, ha='left', va="top")#, va='center')#, size=50)
    
    # relative to invisible axes
    #ax.annotate('\n'.join(text_list), (0, 1), textcoords='data', size=14)# ha='left', va='center')#, size=50)

    #matplotlib.rc('text', usetex=False)

    
def cmdline_parser():
    """
    creates an OptionParser instance
    """

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument("-v", "--verbose",
                      action="store_true",
                      dest="verbose",
                      help="be verbose")
    parser.add_argument("--debug",
                      action="store_true",
                      dest="debug",
                      help="enable debugging")
    parser.add_argument("-i", "--vcf",
                      dest="vcf",
                      required=False,
                      help="Input vcf file (gzip supported; - for stdin).")
    parser.add_argument("--maxdp",
                      dest="maxdp",
                      type=int,
                      help="Maximum DP")
    parser.add_argument("-o", "--outplot",
                      dest="outplot",
                      required=True,
                      help="Output plot (pdf) filename")
    return parser


def main():
    """main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    for (in_file, descr) in [(args.vcf, "VCF")]:
        if not in_file:
            parser.error("%s input file argument missing." % descr)
            sys.exit(1)
        if not os.path.exists(in_file) and in_file != "-":
            sys.stderr.write(
                "file '%s' does not exist.\n" % in_file)
            sys.exit(1)

    for (out_file, descr) in [(args.outplot, "plot")]:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            sys.stderr.write(
                "Cowardly refusing to overwrite existing"
                " output file '%s'.\n" % out_file)
            sys.exit(1)

            # FIXME example input
    #vcf_fh = open("./Q25/NCOV_30_TCOV500_GC95/lofreq-uniq-relax_bed/normal.vcf")
    #vcf_fh = gzip.open("/Users/wilma/Desktop//mutascope/test_files/dbsnp_135_chr22.vcf.gz")# no LoFreq markup
    #vcf_fh = gzip.open("/Users/wilma/GIS/lofreq//lofreq2-test-data/vcf/WHH021_WHH022_stringent.vcf.gz")# too big
    #vcf_fh = open("/Users/wilma/Desktop/projects/bertrandd/lofreq/CML/lofreq-64420af/WHH019-WHH020/lofreq-64420af_def_somatic_final.vcf")
    #vcf_fh = open("/Users/wilma/scratch/cml/outdat_pass.unique.vcf")
    #vcf_fh = open("/Users/wilma/scratch/cml/lofreq_95727fb_a10-SV_somatic_final.uniq.vcf")
    #vcf_fh = open("/Users/wilma/scratch/cml/lofreq_95727fb_a10_somatic_final.uniq.vcf")

    if args.vcf[-3:] == '.gz':
        vcf_fh = gzip.open(args.vcf)
    else:
        vcf_fh = open(args.vcf)
    vcfreader = vcf.VCFReader(vcf_fh)
    vars = [v for v in vcfreader if v.FILTER in ['PASS', '.']]
    vcf_fh.close()

    summary_txt = []
    summary_txt.append("Reading vars from %s" % args.vcf)
    LOG.info(summary_txt[-1])
    summary_txt.append("Loaded %d (non-filtered) vars" % (len(vars)))
    LOG.info(summary_txt[-1])

    
    filter_list = []
    if args.maxdp:
        filter_list.append((lambda v: v.INFO['DP']<=args.maxdp, "DP<=%d" % args.maxdp))
    #filter_list.append(lambda v: v.CHROM=='chr1')
    filtered_vars = vars
    for (f, n) in filter_list:
        n_in = len(filtered_vars)
        try:
            filtered_vars = [v for v in filtered_vars if f(v)]
        except:
            LOG.fatal("Filter %s failed" % n)
            raise 
        n_out = len(filtered_vars)
        summary_txt.append("Filter '%s' removed %d (more) vars" % (n, n_in-n_out))
        LOG.info(summary_txt[-1])
        
    summary_txt.append("%d vars left after filtering" % (len(filtered_vars)))
    LOG.info(summary_txt[-1])

    vars = filtered_vars
    
    summary_txt.append("#SNVs = %d (%d CONSVARs and %d INDELs)" % (
        len(vars),
        sum([1 for v in vars if v.INFO.has_key('CONSVAR')]),
        sum([1 for v in vars if v.INFO.has_key('INDEL')])))
    LOG.info(summary_txt[-1])

    # np.histogram([v.INFO['DP'] for v in vars if v.INFO['DP']<1000], bins=20)


    pp = PdfPages(args.outplot)
        
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
    #f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False)        
    # multiple columns and rows seem to only make sense with shared x or y
    #ax = fig1.add_subplot(2, 2, 1)


    # FIXME: report quality/p-value distribution

    # create a summary table
    # 
    #matplotlib.rc('text', usetex=False)
    fig = plt.figure()
    ax = plt.subplot(1,1,1)    
    print_overview(ax, summary_txt)
    plt.title('Overview')
    pp.savefig()
    plt.close()


        
    # FIXME: lim wrong
    # FIXME bp overwrites violin or doesn't work at all
        
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    x = [v.INFO['DP'] for v in vars]
    #x = [v.INFO['AF'] for v in vars]
    #x = [log10(d) if d>0 else -1 for d in dist_to_next]
    ax.boxplot(x, notch=1, positions=[0], vert=1)
    violin_plot(ax, x)
    ax.set_xlabel('DP')
    plt.title('DP distribution')
    pp.savefig()
    plt.close()

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    cov_hist(ax, [v.INFO['DP'] for v in vars])
    plt.title('DP Histogram')
    pp.savefig()
    plt.close()
    
            
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    x = [v.INFO['AF'] for v in vars]
    #x = [v.INFO['AF'] for v in vars]
    #x = [log10(d) if d>0 else -1 for d in dist_to_next]
    ax.boxplot(x, notch=1, positions=[0], vert=1)
    violin_plot(ax, x)
    ax.set_xlabel('AF')
    plt.title('AF distribution')
    pp.savefig()
    plt.close()

    
    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    af_hist(ax, [v.INFO['AF'] for v in vars], bins=15)#vars])
    plt.title('AF Histogram')
    pp.savefig()
    plt.close()

        

    # add distance to closest SNV to all vars and keep as separate array (FIXME: could do this for all properties)
    # taken from LoFreq's 0.6.1 lofreq_filter.py
    # determines SNVs on same chrom to left and to the right (ignore multip-allelic, might be none) and record pos
    # FIXME this is slow. Could use groupby for chroms, remove multi-allelic (no need for tests then) etc.
    #
    # # group by chromosome
    # processed_chroms = []
    # for (chrom, group) in itertools.groupby(vars, lambda v: v.CHROM):
    #    assert chrom not in processed_chroms;# to ensure vars where sorted by chrom
    #    processed_chroms.append(chrom)
    #    for var in group:
    #        print v
    #


    dist_to_next = calc_dist(vars)

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    dist_hist(ax, dist_to_next)
    plt.title('Distance Histogram')
    pp.savefig()
    plt.close()

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    heatmap_dist_vs_x_1(ax, vars, dist_to_next, 'DP')
    plt.title('Distance vs DP')
    pp.savefig()
    plt.close()

    if False:
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        im = heatmap_dist_vs_cov_2(ax, vars, dist_to_next)
        fig.colorbar(im)
        plt.title('Distance vs. DF (alternate, smoothed version)')
        pp.savefig()
        plt.close()


    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    heatmap_dist_vs_x_1(ax, vars, dist_to_next, "AF")
    plt.title('Distance vs. AF')
    pp.savefig()
    plt.close()

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    heatmap_freq_vs_cov(ax, [v for v in vars if not v.INFO.has_key('CONSVAR')])
    plt.title('AF vs. Coverage')
    pp.savefig()
    plt.close()


    # FIXME needs percentages
    subst_type_counts = Counter([subst_type_str(v.REF, v.ALT) for v in vars])
    # turn into list of tuples sorted by key
    # subst_type_counts = sorted((k, v/100.0*len(vars)) for (k, v) in subst_type_counts.items())
    subst_type_counts = sorted(subst_type_counts.items())
    # FIXME should go to text report
    #for (k, v) in subst_type_counts:
    #    print "%s %d" % (k, v)
    #print

    fig = plt.figure()
    ax = plt.subplot(1, 1, 1)
    subst_perc(ax, subst_type_counts)
    plt.title('Substitution Types')
    pp.savefig()
    plt.close()

    
    LOG.critical("Put related plots together. See http://blog.marmakoide.org/?p=94")
    """
    # FIXME setup such that they don't share axis
    # FIXME other way to make subplots
    # FIXME create figures separately and add separately?

    NROWS = 8
    NCOLS = 3

    ax = plt.subplot(NROWS, NCOLS, 1)
    cov_hist(ax, [v.INFO['DP'] for v in vars])

    ax = plt.subplot(NROWS, NCOLS, 3)
    af_hist(ax, [v.INFO['AF'] for v in vars])

    ax = plt.subplot(NROWS, NCOLS, 7)
    subst_perc(ax, subst_type_counts)

    ax = plt.subplot(NROWS, NCOLS, 9)
    dist_hist(ax, dist_to_next)

    ax = plt.subplot(NROWS, NCOLS, 13)
    heatmap_dist_vs_cov_1(ax, vars, dist_to_next)

    ax = plt.subplot(NROWS, NCOLS, 15)
    im = heatmap_dist_vs_cov_2(ax, vars, dist_to_next)
    fig.colorbar(im)

    ax = plt.subplot(NROWS, NCOLS, 19)
    heatmap_freq_vs_cov(ax, [v for v in vars if not v.INFO.has_key('CONSVAR')])

    ax = plt.subplot(NROWS, NCOLS, 21)
    x = [v.INFO['DP'] for v in vars]
    #x = [v.INFO['AF'] for v in vars]
    #x = [log10(d) if d>0 else -1 for d in dist_to_next]
    violin_plot(ax, x, True)
    ax.set_xlabel('DP')
    """

    pp.close()

    
if __name__ == "__main__":
    #  FIXME unfinished
    main()
    LOG.info("Successful program exit")
