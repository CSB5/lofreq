#!/usr/bin/env python
"""
NOTE: see also older versions of this file
"""



__author__ = "Andreas Wilm"
__version__ = "0.1"
__email__ = "andreas.wilm@ucd.ie"
__copyright__ = ""
__license__ = ""
__credits__ = [""]
__status__ = ""


# --- standard library imports
#
import os
import sys

# --- third-party imports
#
# /

# --- project specific imports
#
# /


# --- globals
#
FIELDNAME_INFO = "Info"
# support for old field-names
FIELDNAME_STDDEV = "SD"
FIELDNAME_PVALUE = "P-Value"


    
class SNP(object):
    """
    A simple and generic SNP class.

    Original used http://www.hgvs.org/mutnomen/recs.html as template but than got rid of seqtype
    """

    def __eq__(self, other):
        """
        if you define __eq__ you better implement __ne__ as well
        """

        if self.pos == other.pos and \
           self.wildtype == other.wildtype and \
           self.variant == other.variant:
            return True
        else:
            return False


    def __ne__(self, other):
        """
        if you define __ne__ you better implement __eq__ as well
        """

        return not self.__eq__(other)

        
    def __init__(self, pos, wildtype, variant):
        """
        pos should be zero-offset
        """

        assert wildtype != variant        
        self.pos = pos
        self.wildtype = wildtype
        self.variant = variant
    

    def __str__(self):
        """
        http://stackoverflow.com/questions/1436703/difference-between-str-and-repr-in-python:
        The goal of __repr__ is to be unambiguous.
        The goal of __str__ is to be readable.
        """

        outstr = "%d %s>%s" % (self.pos+1, self.wildtype, self.variant)
        return outstr


    def __hash__(self):
        """
        Needed for sets
        """
        return hash("%d %s>%s" % (self.pos+1, self.wildtype, self.variant))




class ExtSNP(SNP):
    """
    Extended SNP representation used for the Dengue SNP caller
    project. Predictions can either come from the 'MFreq' method
    (having a frequency and a standard deviation) or from the 'MQual'
    method (having a frequency and a p-value). Allow arbitrary markup.
    Values will be stored as key=value;[key=value;...]. info should be
    a dictionary of key value pairs.

    NOTE: might be a good idea to associate arbitray data but would
    have to make sure that this is there for all SNPs then
    
    For comparison of two SNPs freq/stddev and pvalue will be ignored,
    ie only the basic SNP class comparison will be used.
    """
 
    def __init__(self, pos, wildtype, variant, freq, info=dict()):
        """
        pos should be zero offset
        """
 
        SNP.__init__(self, pos, wildtype, variant)

        self.freq = freq
        if info:
            assert isinstance(info, dict)
            for v in info:
                # used for constructing info string and therefore not
                # allowed to be in values
                assert ';' not in v
                assert ':' not in v
                

        self.info = info
            
        

    def __str__(self):
        """
        """

        #outstr = SNP.__str__(self)
        #if self.stddev:
        #    outstr = "%s %f" % (outstr, self.freq)
        #elif self.pvalue:
        #   outstr = "%s %f" % (outstr, self.pvalue)

        outstr = '%d %s>%s %g' % (self.pos+1,
                                    self.wildtype, self.variant,
                                    self.freq)
        if self.info:
            info_str = ';'.join(["%s:%s"  % (k, v)
                                 for k, v in sorted(self.info.iteritems())])
            outstr = "%s %s" % (outstr, info_str)
            
        return outstr


    def __hash__(self):
        """
        Needed for sets. Don't just use str which includes extended info.
        SNPs with different freqs should still be the same
        """

        return SNP.__hash__(self)




    def identifier(self):
        """
        Returns an mappable identifier for this SNP, which is it's
        basic SNP class string representation (which does not contain
        freq, pvalue or stddev)
        """

        return SNP.__str__(self)


                
class DengueSNP(ExtSNP):
    """
    For backward compatibility when ExtSNP used to be called ExtSNP
    """
    
    def __init__(self, pos, wildtype, variant, freq,
                 stddev=None, pvalue=None):
        """
        Pos should be zero offset
        """

        ExtSNP.__init__(self, pos, wildtype, variant,
                        stddev=stddev, pvalue=pvalue)
            
    


def write_snp_file(fhandle, snp_list, write_header=True):
    """
    Writes SNP to a filehandle
    """
    
    if write_header:
        header = "Pos SNP Freq Info"
        fhandle.write("%s\n" % header)
    for snp in sorted(snp_list, key=lambda s: s.pos):
        fhandle.write("%s\n" % snp)

    
  
    

def parse_snp_file(filename, extra_fieldname='pvalue', has_header=False):
    """
    Parses a Dengue SNP file and returns parsed SNPs as list of
    ExtSNP instances. Missing wildtypes are allowed and replaced with N's.
    """

    ret_snp_list = []
    
    fhandle = open(filename, 'r')

    for line in fhandle:
        line = line.rstrip(os.linesep)
        if len(line)==0 or line.startswith("#") or line.startswith("Pos"):
            continue
        
        try:
            (pos, snp_str, freq, info_str) = line.split()
            pos = int(pos)-1 # internally zero offset
        except ValueError:
            sys.stderr.write(
                "WARN: Failed to parse line from %s. Line was '%s'" % (filename, line))
        freq = float(freq)
        # old versions don't contain the wildtype. use 'N' instead.
        if ">" in snp_str:
            (wildtype, variant) = snp_str.split(">")
        else:
            wildtype = 'N'
            variant = snp_str

        try:
            info = dict([e.split(':') for e in info_str.split(';')])
        except ValueError:
            info = {"generic-info": info_str}
        new_snp = ExtSNP(pos, wildtype, variant, freq, info)
        ret_snp_list.append(new_snp)
        
    fhandle.close()
    return ret_snp_list



def test(fsnp_in=None, fsnp_out=None):
    """
    a test function
    """

    snp = SNP(111, 'C', 'T')
    print "Got SNP: %s" % snp
    
    snp = SNP(666, 'G', 'A', 'r')
    print "Got SNP: %s" % snp

    info = {'pvalue': 0.025}
    dengue_snp = ExtSNP(12345, 'A', 'U', 2e-4, info)
    print "Got Dengue SNP: %s" % dengue_snp

    if fsnp_in:
        print "Parsing %s" % fsnp_in
        snps = parse_snp_file(fsnp_in)
        for snp in snps:
            print "Parsed SNP: %s" % snp
        if fsnp_out:
            print "Writing to %s" % fsnp_out
            write_snp_file(fsnp_out, snps)


if __name__ == "__main__":
    import sys
    try:
        fsnp_in = sys.argv[1]
    except IndexError:
        fsnp_in = None
    try:
        fsnp_out = sys.argv[2]
    except IndexError:
        fsnp_out = None
    test(fsnp_in, fsnp_out)
    sys.exit(0)
