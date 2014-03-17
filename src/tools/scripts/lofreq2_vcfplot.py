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

import vcf

#--- project specific imports
#
try:
    import lofreq2_local
except ImportError:
    pass

#try:
#    #from lofreq_star import vcf
#except ImportError:
#    sys.stderr.write("FATAL(%s): Couldn't find LoFreq's vcf module."
#                     " Are you sure your PYTHONPATH is set correctly (= %s)?\n" % (
#                         (sys.argv[0], os.environ['PYTHONPATH'])))
#    sys.exit(1)
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



def r_ify(axes):
    '''FIXME:unused

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


def subst_type_str(ref, alt, strand_specific=False):
    """FIXME:add-doc
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
        # FIXME LOG.warn("v.max()==0. won't be able to correctly print violin_plot")
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
    parser.add_argument("--summary-only",
                      action="store_true",
                      help="Don't plot; summarize only")
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

    out_files_and_descr = []
    if not args.summary_only:
        out_files_and_descr = (args.outplot, "plot")
    for (out_file, descr) in []:
        if not out_file:
            parser.error("%s output file argument missing." % descr)
            sys.exit(1)
        if os.path.exists(out_file) and out_file!="-":
            sys.stderr.write(
                "Cowardly refusing to overwrite existing"
                " output file '%s'.\n" % out_file)
            sys.exit(1)


    if args.vcf[-3:] == '.gz':
        vcf_fh = gzip.open(args.vcf)
    else:
        vcf_fh = open(args.vcf)
    vcfreader = vcf.VCFReader(vcf_fh)
    # v.FILTER is empty if not set in pyvcf. LoFreq's vcf.py clone set it to PASS or .
    vars = [v for v in vcfreader if not v.FILTER or v.FILTER in ['PASS', '.']]
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



    # setup props we want to check in all possible combinations
    #
    props = dict()
    for t in ['AF', 'DP']:
        try:
            props[t] = [v.INFO[t] for v in vars]
        except KeyError:
            LOG.critical("Couldn't find %s info tag in all variants"
            " (is %s a LoFreq file?). Won't plot..." % (t, args.vcf))
    props['Distance (log10)'] = [np.log10(d) if d>0 else -1 for d in calc_dist(vars)]


    if args.summary_only:
        for p in [p for p in props.keys()]:
            x = np.array(props[p])
            for (name, val) in [("minimum", np.min(x)),
                                ("1st %ile", np.percentile(x, 1)),
                                ("25th %ile", np.percentile(x, 25)),
                                ("median", np.percentile(x, 50)),
                                ("75th %ile", np.percentile(x, 75)),
                                ("99th %ile", np.percentile(x, 99)),
                                ("maximum", np.max(x))]:
                print "%s\t%s\t%f" % (p, name, val)
            print "%s\trange-min\trange-max\tcount" % (p)
            (hist, bin_edges) = np.histogram(x)
            for (i, val) in enumerate(hist):
                print "%f\t%f\t%d" % (bin_edges[i], bin_edges[i+1], val)
        return
    
    pp = PdfPages(args.outplot)

    # create a summary table
    #
    #matplotlib.rc('text', usetex=False)
    fig = plt.figure()
    ax = plt.subplot(1,1,1)
    print_overview(ax, summary_txt)
    plt.title('Overview')
    pp.savefig()
    plt.close()


    # boxplots and histograms first
    #
    for p in [p for p in props.keys()]:
        # boxplots
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        x = props[p]
        ax.boxplot(x, notch=1, positions=[0], vert=1)
        violin_plot(ax, x)
        ax.set_ylabel('#SNVs')
        ax.set_xlabel(p)
        plt.title('%s Boxplot' % p)
        pp.savefig()
        plt.close()

        # histogram
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        x = props[p]
        ax.hist(x, bins=20)
        ax.set_xlim([0, plt.xlim()[1]])
        ax.set_ylabel('#SNVs')
        ax.set_xlabel(p)
        plt.title('%s Histogram' % p)
        pp.savefig()
        plt.close()

        # scatter plot per positions. assuming snvs are sorted by position!
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)
        y = props[p]
        ax.scatter(range(len(y)), y)
        ax.set_xlim([0, len(y)])
        ax.set_ylabel(p)
        ax.set_xlabel("Neighbourhood")
        #plt.title('%s Histogram' % p)
        pp.savefig()
        plt.close()


    # heatmaps of all combinations
    #
    for (x, y) in itertools.combinations(props.keys(), 2):
        fig = plt.figure()
        ax = plt.subplot(1, 1, 1)

        p = plt.hist2d(props[x], props[y], bins=20)
        ax.set_ylim([0, plt.ylim()[1]])
        ax.set_xlim([0, plt.xlim()[1]])
        plt.colorbar()

        ax.set_xlabel(x)
        ax.set_ylabel(y)
        plt.title('%s vs. %s' % (x, y))
        pp.savefig()
        plt.close()



    # substitution types
    #
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


    # FIXME Put related plots together. See http://blog.marmakoide.org/?p=94"

    pp.close()


if __name__ == "__main__":
    main()
    LOG.info("Successful program exit")
