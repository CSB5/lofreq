#!/usr/bin/env python
'''A VCFv4.0 parser for Python.

The intent of this module is to mimic the ``csv`` module in the Python stdlib,
as opposed to more flexible serialization formats like JSON or YAML.  ``vcf``
will attempt to parse the content of each record based on the data types
specified in the meta-information lines --  specifically the ##INFO and
##FORMAT lines.  If these lines are missing or incomplete, it will check
against the reserved types mentioned in the spec.  Failing that, it will just
return strings.

There is currently one piece of interface: ``VCFReader``.  It takes a file-like
object and acts as a reader::

    >>> import vcf
    >>> vcf_reader = vcf.VCFReader(open('example.vcf', 'rb'))
    >>> for record in vcf_reader:
    ...     print record
    Record(CHROM='20', POS=14370, ID='rs6054257', REF='G', ALT=['A'], QUAL=29,
    FILTER='PASS', INFO={'H2': True, 'NS': 3, 'DB': True, 'DP': 14, 'AF': [0.5]
    }, FORMAT='GT:GQ:DP:HQ', samples=[{'GT': '0', 'HQ': [58, 50], 'DP': 3, 'GQ'
    : 49, 'name': 'NA00001'}, {'GT': '0', 'HQ': [65, 3], 'DP': 5, 'GQ': 3, 'nam
    e' : 'NA00002'}, {'GT': '0', 'DP': 3, 'GQ': 41, 'name': 'NA00003'}])

This produces a great deal of information, but it is conveniently accessed.
The attributes of a Record are the 8 fixed fields from the VCF spec plus two
more.  That is:

    * ``Record.CHROM``
    * ``Record.POS``
    * ``Record.ID``
    * ``Record.REF``
    * ``Record.ALT``
    * ``Record.QUAL``
    * ``Record.FILTER``
    * ``Record.INFO``

plus two more attributes to handle genotype information:

    * ``Record.FORMAT``
    * ``Record.samples``

``samples``, not being the title of any column, is left lowercase.  The format
of the fixed fields is from the spec.  Comma-separated lists in the VCF are
converted to lists.  In particular, one-entry VCF lists are converted to
one-entry Python lists (see, e.g., ``Record.ALT``).  Semicolon-delimited lists
of key=value pairs are converted to Python dictionaries, with flags being given
a ``True`` value. Integers and floats are handled exactly as you'd expect::

    >>> record = vcf_reader.next()
    >>> print record.POS
    17330
    >>> print record.ALT
    ['A']
    >>> print record.INFO['AF']
    [0.017]

``record.FORMAT`` will be a string specifying the format of the genotype
fields.  In case the FORMAT column does not exist, ``record.FORMAT`` is
``None``.  Finally, ``record.samples`` is a list of dictionaries containing the
parsed sample column::

    >>> record = vcf_reader.next()
    >>> for sample in record.samples:
    ...     print sample['GT']
    '1|2'
    '2|1'
    '2/2'

Metadata regarding the VCF file itself can be investigated through the
following attributes:

    * ``VCFReader.metadata``
    * ``VCFReader.infos``
    * ``VCFReader.filters``
    * ``VCFReader.formats``
    * ``VCFReader.samples``

For example::

    >>> vcf_reader.metadata['fileDate']
    20090805
    >>> vcf_reader.samples
    ['NA00001', 'NA00002', 'NA00003']
    >>> vcf_reader.filters
    {'q10': Filter(id='q10', desc='Quality below 10'),
    's50': Filter(id='s50', desc='Less than 50% of samples have data')}
    >>> vcf_reader.infos['AA'].desc
    Ancestral Allele

'''
import collections
import re


# Metadata parsers/constants
RESERVED_INFO = {
    'AA': 'String', 'AC': 'Integer', 'AF': 'Float', 'AN': 'Integer',
    'BQ': 'Float', 'CIGAR': 'String', 'DB': 'Flag', 'DP': 'Integer',
    'END': 'Integer', 'H2': 'Flag', 'MQ': 'Float', 'MQ0': 'Integer',
    'NS': 'Integer', 'SB': 'String', 'SOMATIC': 'Flag', 'VALIDATED': 'Flag'
}

RESERVED_FORMAT = {
    'GT': 'String', 'DP': 'Integer', 'FT': 'String', 'GL': 'Float',
    'GQ': 'Float', 'HQ': 'Float'
}


_Info = collections.namedtuple('Info', ['id', 'num', 'type', 'desc'])
_Filter = collections.namedtuple('Filter', ['id', 'desc'])
_Format = collections.namedtuple('Format', ['id', 'num', 'type', 'desc'])


class _vcf_metadata_parser(object):
    '''Parse the metadat in the header of a VCF file.'''
    def __init__(self, aggressive=False):
        super(_vcf_metadata_parser, self).__init__()
        self.aggro = aggressive
        self.info_pattern = re.compile(r'''\#\#INFO=<
            ID=(?P<id>[^,]+),
            Number=(?P<number>\d+|\.|[AG]),
            Type=(?P<type>Integer|Float|Flag|Character|String),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.filter_pattern = re.compile(r'''\#\#FILTER=<
            ID=(?P<id>[^,]+),
            Description="(?P<desc>[^"]*)"
            >''', re.VERBOSE)
        self.format_pattern = re.compile(r'''\#\#FORMAT=<
            ID=(?P<id>.+),
            Number=(?P<number>\d+|\.|[AG]),
            Type=(?P<type>.+),
            Description="(?P<desc>.*)"
            >''', re.VERBOSE)
        self.meta_pattern = re.compile(r'''##(?P<key>.+)=(?P<val>.+)''')

    def read_info(self, info_string):
        '''Read a meta-information INFO line.'''
        match = self.info_pattern.match(info_string)
        if not match:
            raise SyntaxError(
                "One of the INFO lines is malformed: {}".format(info_string))

        try:
            num = int(match.group('number'))
        except ValueError:
            num = None if self.aggro else '.'

        info = _Info(match.group('id'), num,
                     match.group('type'), match.group('desc'))

        return (match.group('id'), info)

    def read_filter(self, filter_string):
        '''Read a meta-information FILTER line.'''
        match = self.filter_pattern.match(filter_string)
        if not match:
            raise SyntaxError(
                "One of the FILTER lines is malformed: {}".format(
                    filter_string))

        filt = _Filter(match.group('id'), match.group('desc'))

        return (match.group('id'), filt)

    def read_format(self, format_string):
        '''Read a meta-information FORMAT line.'''
        match = self.format_pattern.match(format_string)
        if not match:
            raise SyntaxError(
                "One of the FORMAT lines is malformed: {}".format(
                    format_string))

        try:
            num = int(match.group('number'))
        except ValueError:
            num = None if self.aggro else '.'

        form = _Format(match.group('id'), num,
                       match.group('type'), match.group('desc'))

        return (match.group('id'), form)

    def read_meta(self, meta_string):
        match = self.meta_pattern.match(meta_string)
        return match.group('key'), match.group('val')


# Reader class
class _meta_info(object):
    '''Decorator for a property stored in the header info.'''
    def __init__(self, func):
        self.func = func

    def __call__(self, fself):
        if getattr(fself, "_%s" % self.func.__name__) is None:
            fself._parse_metainfo()

        return self.func(fself)

    def __repr__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

    def __doc__(self):
        '''Return the function's docstring.'''
        return self.func.__doc__

_Record = collections.namedtuple('Record', [
    'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
    'samples'
])


class VCFReader(object):
    '''Read and parse a VCF v 4.0 file'''
    def __init__(self, fsock, aggressive=False):
        super(VCFReader, self).__init__()
        self.aggro = aggressive
        self._metadata = None
        self._infos = None
        self._filters = None
        self._formats = None
        self._samples = None
        self.reader = fsock
        if aggressive:
            self._mapper = self._none_map
        else:
            self._mapper = self._pass_map

    def __iter__(self):
        return self

    @property
    @_meta_info
    def metadata(self):
        '''Return the information from lines starting "##"'''
        return self._metadata

    @property
    @_meta_info
    def infos(self):
        '''Return the information from lines starting "##INFO"'''
        return self._infos

    @property
    @_meta_info
    def filters(self):
        '''Return the information from lines starting "##FILTER"'''
        return self._filters

    @property
    @_meta_info
    def formats(self):
        '''Return the information from lines starting "##FORMAT"'''
        return self._formats

    @property
    @_meta_info
    def samples(self):
        '''Return the names of the genotype fields.'''
        return self._samples

    def _parse_metainfo(self):
        '''Parse the information stored in the metainfo of the VCF.

        The end user shouldn't have to use this.  She can access the metainfo
        directly with ``self.metadata``.'''
        for attr in ('_metadata', '_infos', '_filters', '_formats'):
            setattr(self, attr, {})

        parser = _vcf_metadata_parser()

        line = self.reader.next()
        while line.startswith('##'):
            line = line.strip()
            if line.startswith('##INFO'):
                key, val = parser.read_info(line)
                self._infos[key] = val

            elif line.startswith('##FILTER'):
                key, val = parser.read_filter(line)
                self._filters[key] = val

            elif line.startswith('##FORMAT'):
                key, val = parser.read_format(line)
                self._formats[key] = val

            else:
                key, val = parser.read_meta(line.strip())
                self._metadata[key] = val

            line = self.reader.next()

        fields = line.split()
        self._samples = fields[9:]

    def _none_map(self, func, iterable, bad='.'):
        '''``map``, but make bad values None.'''
        return [func(x) if x != bad else None
                for x in iterable]

    def _pass_map(self, func, iterable, bad='.'):
        '''``map``, but make bad values None.'''
        return [func(x) if x != bad else bad
                for x in iterable]

    def _parse_info(self, info_str):
        '''Parse the INFO field of a VCF entry into a dictionary of Python
        types.

        '''
        entries = info_str.split(';')
        retdict = {}
        for entry in entries:
            entry = entry.split('=')
            ID = entry[0]
            try:
                entry_type = self.infos[ID].type
            except KeyError:
                try:
                    entry_type = RESERVED_INFO[ID]
                except KeyError:
                    if entry[1:]:
                        entry_type = 'String'
                    else:
                        entry_type = 'Flag'

            if entry_type == 'Integer':
                vals = entry[1].split(',')
                val = self._mapper(int, vals)
            elif entry_type == 'Float':
                vals = entry[1].split(',')
                val = self._mapper(float, vals)
            elif entry_type == 'Flag':
                val = True
            elif entry_type == 'String':
                val = entry[1]

            try:
                if self.infos[ID].num == 1:
                    val = val[0]
            except KeyError:
                pass

            retdict[ID] = val

        return retdict

    def _parse_samples(self, samples, samp_fmt):
        '''Parse a sample entry according to the format specified in the FORMAT
        column.'''
        samp_data = []
        samp_fmt = samp_fmt.split(':')
        for sample in samples:
            sampdict = dict(zip(samp_fmt, sample.split(':')))
            for fmt in sampdict:
                vals = sampdict[fmt].split(',')
                try:
                    entry_type = self.formats[fmt].type
                except KeyError:
                    try:
                        entry_type = RESERVED_FORMAT[fmt]
                    except KeyError:
                        entry_type = 'String'

                if entry_type == 'Integer':
                    sampdict[fmt] = self._mapper(int, vals)
                elif entry_type == 'Float' or entry_type == 'Numeric':
                    sampdict[fmt] = self._mapper(float, vals)
                elif sampdict[fmt] == './.' and self.aggro:
                    sampdict[fmt] = None

            samp_data.append(sampdict)

        for name, data in zip(self.samples, samp_data):
            data['name'] = name

        return samp_data

    def next(self):
        '''Return the next record in the file.'''
        if self._samples is None:
            self._parse_metainfo()
        row = self.reader.next().split()
        chrom = row[0]
        pos = int(row[1])

        if row[2] != '.':
            ID = row[2]
        else:
            ID = None if self.aggro else row[2]

        ref = row[3]
        alt = self._mapper(str, row[4].split(','))
        qual = float(row[5]) if '.' in row[5] else int(row[5])
        filt = row[6].split(';') if ';' in row[6] else row[6]
        if filt == 'PASS' and self.aggro:
            filt = None
        info = self._parse_info(row[7])

        try:
            fmt = row[8]
        except IndexError:
            fmt = None
            samples = None
        else:
            samples = self._parse_samples(row[9:], fmt)

        record = _Record(chrom, pos, ID, ref, alt, qual, filt, info, fmt,
                         samples)
        return record


def main():
    '''Parse the example VCF file from the specification and print every
    record.'''
    import contextlib
    import StringIO
    import textwrap
    buff = '''\
        ##fileformat=VCFv4.0
        ##fileDate=20090805
        ##source=myImputationProgramV3.1
        ##reference=1000GenomesPilot-NCBI36
        ##phasing=partial
        ##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
        ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
        ##INFO=<ID=AF,Number=.,Type=Float,Description="Allele Frequency">
        ##INFO=<ID=AA,Number=1,Type=String,Description="Ancestral Allele">
        ##INFO=<ID=DB,Number=0,Type=Flag,Description="dbSNP membership, build 129">
        ##INFO=<ID=H2,Number=0,Type=Flag,Description="HapMap2 membership">
        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
        ##FILTER=<ID=q10,Description="Quality below 10">
        ##FILTER=<ID=s50,Description="Less than 50% of samples have data">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
        ##FORMAT=<ID=HQ,Number=2,Type=Integer,Description="Haplotype Quality">
        #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA00001\tNA00002\tNA00003
        20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51\t1/1:43:5:.,.
        20\t17330\t.\tT\tA\t3\tq10\tNS=3;DP=11;AF=0.017\tGT:GQ:DP:HQ\t0|0:49:3:58,50\t0|1:3:5:65,3\t0/0:41:3
        20\t1110696\trs6040355\tA\tG,T\t67\tPASS\tNS=2;DP=10;AF=0.333,0.667;AA=T;DB\tGT:GQ:DP:HQ\t1|2:21:6:23,27\t2|1:2:0:18,2\t2/2:35:4
        20\t1230237\t.\tT\t.\t47\tPASS\tNS=3;DP=13;AA=T\tGT:GQ:DP:HQ\t0|0:54:7:56,60\t0|0:48:4:51,51\t0/0:61:2
        20\t1234567\tmicrosat1\tGTCT\tG,GTACT\t50\tPASS\tNS=3;DP=9;AA=G\tGT:GQ:DP\t./.:35:4\t0/2:17:2\t1/1:40:3
        '''
    with contextlib.closing(StringIO.StringIO(textwrap.dedent(buff))) as sock:
        vcf_file = VCFReader(sock, aggressive=True)
        for record in vcf_file:
            print record


if __name__ == '__main__':
    main()
