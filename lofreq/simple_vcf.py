#!/usr/bin/env python
"""Simple module for VCF writing

Had initially planned to use John Dougherty's
https://github.com/jdoughertyii/PyVCF but couldn't figure a way to
create records from scratch and didn't want to introduce another
dependency. Used some concepts from there anyway.

"""


# Copyright (C) 2011, 2012 Genome Institute of Singapore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.



#--- standard library imports
#
from __future__ import division
import collections
import datetime
import sys

#--- third-party imports
#
# /

#--- project specific imports
#
# /

__author__ = "Andreas Wilm"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2011, 2012 Genome Institute of Singapore"
__license__ = "GPL2"

"""
See also VCF 4.0 format: http://www.1000genomes.org/node/101

There are 8 fixed fields per record. All data lines are
tab-delimited. In all cases, missing values are specified with a
dot ("."). Fixed fields are:

1. CHROM chromosome: an identifier from the reference genome. All
entries for a specific CHROM should form a contiguous block within
the VCF file.(Alphanumeric String, Required)

2. POS position: The reference position, with the 1st base having
position 1. Positions are sorted numerically, in increasing order,
within each reference sequence CHROM. (Integer, Required)

3. ID semi-colon separated list of unique identifiers where
available. If this is a dbSNP variant it is encouraged to use the
rs number(s). No identifier should be present in more than one
data record. If there is no identifier available, then the missing
value should be used. (Alphanumeric String)

4. REF reference base(s): Each base must be one of A,C,G,T,N.
Bases should be in uppercase. Multiple bases are permitted. The
value in the POS field refers to the position of the first base in
the String. For InDels, the reference String must include the base
before the event (which must be reflected in the POS field).
(String, Required).

5. ALT comma separated list of alternate non-reference alleles
called on at least one of the samples. Options are base Strings
made up of the bases A,C,G,T,N, or an angle-bracketed ID String
("<ID>"). If there are no alternative alleles, then the missing
value should be used. Bases should be in uppercase. (Alphanumeric
String; no whitespace, commas, or angle-brackets are permitted in
the ID String itself)

6. QUAL phred-scaled quality score for the assertion made in ALT.
i.e. give -10log_10 prob(call in ALT is wrong). If ALT is "." (no
variant) then this is -10log_10 p(variant), and if ALT is not "."
this is -10log_10 p(no variant). High QUAL scores indicate high
confidence calls. Although traditionally people use integer phred
scores, this field is permitted to be a floating point to enable
higher resolution for low confidence calls if desired. (Numeric)

7. FILTER filter: PASS if this position has passed all filters,
i.e. a call is made at this position. Otherwise, if the site has
not passed all filters, a semicolon-separated list of codes for
filters that fail. e.g. "q10;s50" might indicate that at this site
the quality is below 10 and the number of samples with data is
below 50% of the total number of samples. "0" is reserved and
should not be used as a filter String. If filters have not been
applied, then this field should be set to the missing value.
(Alphanumeric String)

8. INFO additional information: (Alphanumeric String) INFO fields
are encoded as a semicolon-separated series of short keys with
optional values in the format: <key>=<data>[,data]. Arbitrary keys
are permitted, although the following sub-fields are reserved
(albeit optional):

  - AA ancestral allele
  - AC allele count in genotypes, for each ALT allele, in the same order as listed
  - AF allele frequency for each ALT allele in the same order as listed: use this when estimated from primary data, not called genotypes
  - AN total number of alleles in called genotypes
  - BQ RMS base quality at this position
  - CIGAR cigar string describing how to align an alternate allele to the reference allele
  - DB dbSNP membership
  - DP combined depth across samples, e.g. DP=154
  - END end position of the variant described in this record (esp. for CNVs)
  - H2 membership in hapmap2
  - MQ RMS mapping quality, e.g. MQ=52
  - MQ0 Number of MAPQ == 0 reads covering this record
  - NS Number of samples with data
  - SB strand bias at this position
  - SOMATIC indicates that the record is a somatic mutation, for cancer genomics
  - VALIDATED validated by follow-up experiment
"""


# Similar to PyVCF
Record = collections.namedtuple('Record', [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
])

Infofield = collections.namedtuple('Infofield', [
    'id', 'number', 'type', 'description'
    ])
REQUIRED_INFO_FIELDS = [
    Infofield('AF', 1, 'Float', 
              "allele frequency for ALT allele"), # 1000genomes says number=".". we require exactly one number
    Infofield('DP', 1, 'Integer', 
              "Raw read depth"), # as in samtools; note that 1000genomes says "read depth" only
    Infofield('DP4', 4, 'Integer', 
              "# high-quality (i.e. after filtering) ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"),
    Infofield('SB', 1, 'Integer',
              "Phred-scaled strand bias at this position") # instead of samtools SP use predefined SB, but explcitely phred-scaled
    ]
# DEV NOTE: changes have to be synced with create_record

# DEV NOTE: check output with vcf-validator

MISSING_VAL = "."



def new_header(myname="LoFreq"):
    """Creates a new vcf header (AKA vcf metainformation); returned as
    string
    """

    isodate = datetime.date.today().isoformat()
    
    header = [] 
    header.append("##fileformat=VCFv4.0")
    header.append("##fileDate=%s" % isodate.replace("-", ""))
    header.append("##source=%s" % myname)

    for infofield in REQUIRED_INFO_FIELDS:
        #header.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Raw read depth">')
        info_line = '##INFO=<ID=%s,Number=%d,Type=%s,Description="%s">' % (
            infofield.id, infofield.number, infofield.type, infofield.description)
        header.append(info_line)
        
    header.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

    return '\n'.join(header)


           
def create_record(rec_chrom, rec_pos, rec_id, 
                      rec_ref, rec_alt, rec_qual,
                      rec_filter, rec_info):
    """Creates a new VCF record. rec_info has to be dict() or None.
    Will perform pedantic checks on input
    """
    

    # 1
    assert len(rec_chrom) != 0

    # 2
    assert rec_pos >= 0 # zero offset internally

    # 3
    if not rec_id:
        rec_id = MISSING_VAL

    # 4
    # we don't allow multiple bases here for the sake of simplicity
    assert rec_ref in ['A', 'C', 'G', 'T', 'N']

    # 5
    # we don't allow multiple bases here and only bases for the sake of simplicity
    assert rec_alt in ['A', 'C', 'G', 'T', 'N']
    assert rec_alt != rec_ref
    
    # 6
    if not rec_qual:
        rec_qual = MISSING_VAL
    else:
        assert isinstance(rec_qual, int) or isinstance(rec_qual, float)
    
    # 7
    if not rec_filter:
        rec_filter = MISSING_VAL
    else:
        assert rec_filter != '0'

    # 8
    if not rec_info:
        rec_info = "."
    else:
        assert isinstance(rec_info, dict)

        for infofield in REQUIRED_INFO_FIELDS:
            assert infofield.id in rec_info.keys(), (
                "Required existance of key %s in INFO dict" % (infofield.id))
                
        rec_info = ';'.join(["%s=%s" % (k, v) 
                             for (k, v) in sorted(rec_info.iteritems())])
   
        
    record = Record._make([rec_chrom, rec_pos, rec_id, 
                                rec_ref, rec_alt, rec_qual,
                                rec_filter, rec_info])

    return record


                     
def write_header(header=None, fh=sys.stdout):
    """Writes a VCF header to stream
    """

    if not header:
        header = new_header()
    fh.write("%s\n" % header)

    
    
def write_record(rec, fh=sys.stdout):
    """Writes a vcf record to stream
    """

    fh.write("%s\n" % '\t'.join([str(x) for x in [
        rec.CHROM, rec.POS+1, rec.ID, 
        rec.REF, rec.ALT, rec.QUAL,
        rec.FILTER, rec.INFO]]))


    
def write(records, fh=sys.stdout):
    """Simple wrapper for write_record, adding a header
    """

    write_header(None, fh)
    for rec in records:
        write_record(rec, fh)
    
    

if __name__ == "__main__":

    print "just testing..."

    print "# testing write_header()"
    write_header()



    print "# testing create_record()"

    records = []

    (rec_chrom, rec_pos, rec_id,
     rec_ref, rec_alt, rec_qual,
     rec_filter, rec_info) = ("chr1", 0, None, "C", "A", None, None, None)
    records.append(create_record(rec_chrom, rec_pos, rec_id,
                            rec_ref, rec_alt, rec_qual,
                            rec_filter, rec_info))

    (rec_chrom, rec_pos, rec_id,
     rec_ref, rec_alt, rec_qual,
     rec_filter, rec_info) = ("chr12.fa", 12323213, None, "G", "A", 666, None,
                              {'AF':0.01, 'DP':100, 'DP4':'1,2,3,4', 'SB':123})
    records.append(create_record(rec_chrom, rec_pos, rec_id,
                            rec_ref, rec_alt, rec_qual,
                            rec_filter, rec_info))

    
    print "# testing write_record()"
    
    for rec in records:
        write_record(rec)
    # FIXME tests missing


    print "# testing write()"

    write(records)
