__author__ = 'esteinig'

"""

brigD3
BLAST Ring Image Generator (BRIG) adoption for D3

Alikhan et al. (2012) wrote BRIG:
http://sourceforge.net/projects/brig/

Manual for brigD3 at https://github.com/esteinig/brigD3

Eike J. Steinig
Tong Lab, Menzies School of Health Research
Zenger Lab, James Cook University
eike.steinig@menzies.edu.au, eikejoachim.steinig@my.jcu.edu.au

"""

import json
import os
import csv
import statistics

from subprocess import call
from Bio import SeqIO


class RingGenerator:
    """
    Class: Ring Generator

    Initiate with list of rings and set options for visualization. The generator transforms the ring data into
    a list of dictionaries for writing as JSON. Options can be set via setOptions. The main access is the brigD3 method,
    which initiates the helper class Visualization containing the JS D3 script and sets its parameters in the string.
    The method then writes the JSON and HTML files to working directory.

    Attributes:

    self.rings:     list, ring objects
    self.project:   str, name of output files
    self.radius:    int, radius of the ring center
    self.gap:       int, gap between the rings

    self.data:      list, data dicts for each segment
    """
    def __init__(self, rings, path, project):
        """
        Initialize generator with list of rings.
        """
        self.rings = rings
        self.path = path
        self.project = project
        self.radius = 200
        self.gap = 5
        self.data = []
        self._options = {}
        self._getDataD3()

    def _getDataD3(self):
        """
        Transform data from rings to data appropriate for D3
        """
        d3 = []
        radius = self.radius
        for ring in self.rings:
            for seg in ring.data:
                height = seg.pop('height')
                seg['inner'] = float(radius)
                seg['outer'] = float(radius + height)
                d3.append(seg)
            radius = radius + ring.height + self.gap
        self.data = d3

        return self.data

    def setOptions(self, circle, radius=300, gap=5, project='data', title='brigD3', title_size='300%', title_font='times',
                   ring_opacity=0.8, width=1600, height=800):
        """
        Set options for circle generator and visualization with D3.
        """
        self.radius = radius
        self.gap = gap
        self.project = project
        self._options = {'circle': circle, 'main_title': title, 'title_size': title_size,
                         'title_font': title_font, 'ring_opacity': ring_opacity, 'chart_width': width,
                         'chart_height': height}

    def brigD3(self):
        """
        Write files for brigD3 to working directory.
        """
        print('\nWriting visualization to working directory ...\n')
        viz = Visualization(self._options, self.data, self.path, self.project)
        viz._setScript()
        viz._writeHTML()


class Ring:
    """
    Super Class Ring:

    A ring object reads and transforms data into the data shape accesses by the RingGenerator and JS. The standard Ring
    has base attributes colour, name, height and the tooltip style that can be set by the user via setOptions. The
    standard Ring reads a comma-delimited file via readRing containig raw data on each segment (no header) in
    the columns:

    Start Position, End Position, Color, Height and HTML string for Tooltip

    It writes data in the same format via writeRing. Ring objects can also be merged via mergeRings, which adds multiple
    rings' data (list of ring objects) to the current ring object. Attributes of the ring object which called the
    method are retained.

    Subclasses of the ring object have additional attributes pertaining to their function, as well as different readers
    for data files.

    Attributes:

    self.name:      str, name to be shown in tooltips
    self.color:     str, color of ring, hex or name
    self.height:    int, height of ring
    self.tooltip:   obj, tooltip object to set header and text colors
    """
    def __init__(self):
        """
        Initialize ring object.
        """
        self.color = 'black'
        self.name = 'Genome'
        self.height = 20
        self.path = ''
        self.tooltip = Tooltip()
        self.res_loc = {}
        self._positions = []
        self._popups = []
        self._colors = []
        self._heights = []
        self.data = []

    def mergeRings(self, rings):
        """
        Merge the current ring with a list of ring objects. Add data of ring objects to current ring.
        """
        for ring in rings:
            self.data += ring.data

    def setOptions(self, name='Ring', color='black', height=20, tooltip=None, path='', res_loc={}):
        """
        Set basic attributes for ring object
        """
        self.name = name
        self.color = color
        self.height = height
        self.path = path
        self.res_loc = res_loc
        if tooltip is None:
            self.tooltip = Tooltip()
        else:
            self.tooltip = tooltip

    def _getRing(self):
        """
        Get ring data in dictionary format for Ring Generator and D3.
        """
        resfinder_dictionary = {
            "Fluoroquinolone and aminoglycoside resistance": '#fa0a96',
            "Aminoglycoside resistance": '#6400fa',
            "Beta-lactam resistance": '#d07fff',
            "Phenicol resistance": '#19ff00',
            "Sulphonamide resistance": '#ff6100',
            "Tetracycline resistance": '#00ffe9',
            "Trimethoprim resistance": '#ff0000',
            "Macrolide resistance": '#faff00',
            "Rifampicin resistance": '#007fff',
            "": '#ff7ca8',
            "Inc": '#0000ff'
        }
        print('Generating Ring:', self.name)
        if "_utg" in self.name or "contig_" in self.name:
            for res, loc in self.res_loc.items():
                resgroup = res.split("_")[-1]
                residentity = res.split("_")[-2]
                if "contig_" in self.name:
                    contig_name = res.split("_")[0] + "_" + res.split("_")[1]
                    res_name = res.split("_")[2]
                    contigcheck = str(self.name.strip("\n").split("_")[1]) + "_" + str(self.name.strip("\n").split("_")[2])
                else:
                    contig_name = res.split("_")[0]
                    res_name = res.split("_")[1]
                    contigcheck = str(self.name.strip("\n").split("_")[1])
                if str(contig_name) in contigcheck:
                    try:
                        highest_loc = max(int(x[1]) for x in self._positions)
                        lowest_loc = min(int(x[0]) for x in self._positions)
                        end_loc = int(loc[1])
                        start_loc = int(loc[0])
                        if start_loc >= int(lowest_loc) and end_loc <= int(highest_loc):
                            if resgroup == "Inc":
                                pop = ('<strong><span style="color:#88A2AF">Incompatibility Factor:</span>:</strong><span style="color:white">' + res_name +
                                    '\n</span><br><strong><span style="color:#88A2AF">Location:</span>:</strong><span style="color:white">' + loc[0] + " - " + loc[1] +
                                    '\n</span><br><strong><span style="color:#88A2AF">Identity:</span>:</strong><span style="color:white">' + residentity + "%" +
                                    '\n</span><br>')
                            else:
                                pop = ('<strong><span style="color:#88A2AF">Resistance Gene:</span>:</strong><span style="color:white">' + res_name +
                                    '\n</span><br><strong><span style="color:#88A2AF">Location:</span>:</strong><span style="color:white">' + loc[0] + " - " + loc[1] +
                                    '\n</span><br><strong><span style="color:#88A2AF">ResFinder Identity:</span>:</strong><span style="color:white">' + residentity + "%" +
                                    '\n</span><br>')
                            self._positions.append(loc)
                            self._colors.append(resfinder_dictionary[resgroup])
                            self._popups.append(pop)
                            self._heights.append(20)
                    except ValueError:
                        pass
        n = len(self._positions)
        for i in range(n):
            data_dict = {}
            data_dict['start'] = self._positions[i][0]
            data_dict['end'] = self._positions[i][1]
            data_dict['color'] = self._colors[i]
            data_dict['text'] = self._popups[i]
            data_dict['height'] = self._heights[i]
            self.data.append(data_dict)
        return self.data

    def writeRing(self, file):
        """
        Write raw ring data to comma-delimited file.
        """
        with open(os.path.join(self.path, file), 'w') as outfile:
            w = csv.writer(outfile, delimiter=',')
            d = [[segment['start'], segment['end'], segment['color'], segment['height'], segment['text']]
                 for segment in self.data]
            w.writerows(d)

    def readRing(self, file):
        """
        Read raw ring data from comma-delimited file.
        """
        self._clear()
        with open(os.path.join(self.path, file), 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                data = {}
                self._heights.append(row[3])
                self._colors.append(row[2])
                self._positions.append([row[0], row[1]])
                self._popups.append(row[4])
                data['start'] = row[0]
                data['end'] = row[1]
                data['color'] = row[2]
                data['height'] = row[3]
                data['text'] = row[4]
                self.data.append(data)

    def _clear(self):
        """
        Clear all ring data.
        """
        self._heights = []
        self._colors = []
        self._positions = []
        self._popups = []
        self.data = []


class CoverageRing(Ring):
    """
    Subclass Coverage Ring for depicting coverage matrix across genomes (single or average).
    """
    def __init__(self):
        """
        Initialize super-class ring and attributes for Coverage Ring.
        """
        Ring.__init__(self)
        self.threshold = 0.95
        self.below = '#F1F1F1'

    def readCoverage(self, file, sep='\t', mean=True, n=5):
        """
        Read coverage matrix from file (with header):
        Segment ID, Start Position, End Position, Value Sample1, Value Sample2...
        """
        self._clear()
        with open(os.path.join(self.path, file), 'r') as infile:
            reader = csv.reader(infile, delimiter=sep)
            header = []
            texts = []
            for row in reader:
                if header:
                    start = int(row[1])
                    end = int(row[2])
                    if mean:
                        value = statistics.mean([float(v) for v in row[3:]])
                        cov = 'Mean Coverage'
                    else:
                        value = float(row[n])
                        cov = 'Coverage'
                    self._positions.append((start, end))
                    color = self.below if value < self.threshold else self.color
                    self._colors.append(color)
                    self._heights.append(value*self.height)
                    texts.append([('Genome:', self.name), ('Location:', str(start) + ' - ' + str(end)),
                                 (cov, format(value*100, ".2f") + '%')])
                else:
                    header = row
        self._popups = [self.tooltip.getPopup(text) for text in texts]
        self._getRing()


class AnnotationRing(Ring):
    """
    Sub-class Annotation Ring for depicting genome annotations and SNPs
    """
    def __init__(self):
        """
        Initialize super-class ring and attributes for Annotation Ring.
        """
        Ring.__init__(self)
        self.feature = ['source', 'CDS']
        self.extract = {'organism': 'Organism: ', 'plasmid': 'Plasmid: ', 'gene': 'Gene: ', 'product': 'Product: '}
        self.snp_length = 100
        self.intergenic = 'yellow'
        self.synonymous = 'orange'
        self.non_synonymous = 'red'

    def readSNP(self, file, single=False, n=5):
        """
        Read SNP data from comma_delimited file (without header): SNP ID, Location, Type, Notes
        """
        self._clear()
        with open(os.path.join(self.path, file), 'r') as infile:
            reader = csv.reader(infile)
            for row in reader:
                if single:
                    # If nucleotide is not equal to reference nucleotide or if it is reference nucleotide
                    # Needs a fix, it's a mess...
                    if row[n] != row[4] or n == 4:
                        self._positions.append([int(row[1])-self.snp_length//2, int(row[1])+self.snp_length//2])
                        if row[2] == 'intergenic':
                            self._colors.append(self.intergenic)
                        elif row[2] == 'synonymous':
                            self._colors.append(self.synonymous)
                        elif row[2] == 'non-synonymous':
                            self._colors.append(self.non_synonymous)
                        else:
                            self._colors.append(self.color)
                        self._heights.append(self.height)
                        self._popups.append(self.tooltip.getPopup([('Genome:', self.name), ('SNP:', row[0]),
                                                                   ('Location: ', row[1]), ('Type:', row[2]),
                                                                   ('Note:', row[3])]))
                else:
                    # Single ring for SNPs along Reference, general SNPs across all isolates from SPANDx
                    self._positions.append([int(row[1])-self.snp_length//2, int(row[1])+self.snp_length//2])
                    if row[2] == 'intergenic':
                        self._colors.append(self.intergenic)
                    elif row[2] == 'synonymous':
                        self._colors.append(self.synonymous)
                    elif row[2] == 'non-synonymous':
                        self._colors.append(self.non_synonymous)
                    else:
                        self._colors.append(self.color)
                    self._heights.append(self.height)
                    self._popups.append(self.tooltip.getPopup([('Genome:', self.name), ('SNP:', row[0]), ('Location: ', row[1]),
                                                              ('Type:', row[2]), ('Note:', row[3])]))
        self._getRing()

    def readGenbank(self, file, length):
        """
        Read genbank annotation file and extract relevant features and qualifiers.
        """
        self._clear()
        if file != "":
            genome = SeqIO.read(open(file, "r"), "genbank")
            features = [feature for feature in genome.features if feature.type in self.feature]
            # Only include feature with all qualifiers present.
            clean = []
            for feature in features:
                check = True
                for q in self.extract.keys():
                    if q == 'plasmid':
                        if q in feature.qualifiers:
                            check = True
                        elif q not in feature.qualifiers:
                            check = False
                        else:
                            check = False
                    else:
                        if q not in feature.qualifiers and q != 'organism':
                            check = False
                if check:
                    clean.append(feature)
            # Get tooltips for each extracted feature.
            for feature in features:
                self._positions.append([int(feature.location.start), int(feature.location.end)])
                self._colors.append(self.color)
                self._heights.append(self.height)
                qualifier_texts = []
                for qualifier in self.extract.keys():
                    if qualifier in feature.qualifiers:
                        text_tuple = (self.extract[qualifier], ''.join(feature.qualifiers[qualifier]))
                        qualifier_texts.append(text_tuple)
                qualifier_texts.insert(0, ('Location: ', str(feature.location.start) + '-' + str(feature.location.end)))
                qualifier_texts.insert(0, ('Genome: ', self.name))
                popup = self.tooltip.getPopup(qualifier_texts)
                self._popups.append(popup)
        else:
            self._positions.append([1, int(length)])
            self._colors.append(self.color)
            self._heights.append(self.height)
            qualifier_texts = []
            qualifier_texts.append("plasmid")
            qualifier_texts.insert(0, ('Location: ', str(1) + '-' + str(length)))
            qualifier_texts.insert(0, ('Genome: ', self.name))
            popup = self.tooltip.getPopup(qualifier_texts)
            self._popups.append(popup)
        self._getRing()


class BlastRing(Ring):
    """
    Sub-class Blast Ring, for depicting BLAST comparisons against a reference DB
    """
    def __init__(self):
        """
        Initialize super-class ring and attributes for B<last Ring.
        """
        Ring.__init__(self)
        self.min_identity = 70
        self.min_length = 100
        self.values = []

    def readComparison(self, file):
        """
        Reads comparison files from BLAST output (--outfmt 6)
        """
        self._clear()
        with open(file, 'r') as infile:
            reader = csv.reader(infile, delimiter='\t')
            for row in reader:
                positions = sorted([int(row[8]), int(row[9])])
                if positions[1] - positions[0] >= self.min_length and float(row[2]) >= self.min_identity:
                    self._positions.append(positions)
                    self.values.append(float(row[2]))
        self._colors = [self.color for v in self.values]
        self._heights = [self.height for v in self.values]
        texts = [[('Contig:', self.name), ('BLAST Identity:', str(v) + '%')] for v in self.values]
        self._popups = [self.tooltip.getPopup(text) for text in texts]
        self._getRing()


class Tooltip:
    """
    Tooltip class (under construction), will hold more options to customize Tooltips for D3
    """
    def __init__(self):
        self.text_color = 'white'
        self.head_color = '#88A2AF'

    def getPopup(self, text):
        """
        Converts text - tuple of header and text, i.e. ('Genome:', 'DAR4145') - to HTML string for Tooltip.
        """
        if len(text) > 0:
            popup = ''
            for i in range(len(text)):
                popup += '<strong>' + '<span style="color:' + self.head_color + '">' + text[i][0] +\
                         '</span>' + ':</strong>' + '<span style="color:' + self.text_color + '">' + text[i][1] +\
                         '</span>' + '<br>'
        else:
            popup = '-'
        return popup


class Blaster:
    """
    Class: Blaster

    Convenience module to run BLAST+ (needs to be in $PATH).

    Initialize with a string for the reference genome file (.fasta) and a list of strings for the sequence files
    to be compared (.fasta). Main access is through runBlast. Results attribute hold the string names of the
    output comparison files that can be iterated over to create Blast Rings.

    Attributes:

    self.name_db:       str, name of database to be created from the reference sequence
    self.type:          str, database type, either 'nucl' or 'prot'
    self.mode:          str, type of comparison, either 'blastn' for nucleotide or 'blastp' for protein
    self.path:          str, path of the blast outputs
    """
    def __init__(self, reference, genomes, path):
        """
        Initialize blast object with reference and sequence files
        """
        self.reference = reference
        self.genomes = genomes
        self.path = path
        self.name_db = 'ReferenceDB'
        self.type = 'nucl'
        self.mode = 'blastn'
        self.results = []

    def _getDB(self, path):
        """
        Run makeblastdb to create reference DB.
        """
        call(['makeblastdb', '-in', self.reference, '-dbtype', self.type, '-out', os.path.join(self.path, self.name_db)])
        print('\n')

    def runBLAST(self):
        """
        Blast sequence files against reference DB.
        """
        self._getDB(self.path)
        refname = self.reference.split('.')[0]
        refname_list = refname.split('/')
        referencename = refname_list[-1]
        for genome in self.genomes:
            genname = genome.split('.')[0]
            print('Blasting', genome, 'against Reference DB ...')
            filename = genname + 'vs' + referencename
            call([self.mode, '-query', os.path.join(self.path, genome), '-db', os.path.join(self.path, self.name_db), '-outfmt', '6',
                  '-out', os.path.join(self.path, filename)])
            self.results.append(os.path.join(self.path, filename))
        print('\n')


class Visualization:
    """
    Helper class Visualization, holds script for JS D3. Methods to write replace options from Ring Generator in
    script and write the HTML. Initialize with options dict and data from Ring Generator.
    """
    def __init__(self, options, data, path, project):
        self.options = options
        self.data = data
        self.path = path
        self.project = project
        self.head = '''
                    <!DOCTYPE html>
                    <meta charset="utf-8">
                    <style>

                      div.tooltip {
                        position: absolute;
                        text-align: left;
                        max-width: 300px;
                        padding: 11px;
                        font: 11px sans-serif;
                        background: black;
                        border-radius: 11px;
                        pointer-events: none;
                      }

                    </style>

                    <body>
                    <script src="https://d3js.org/d3.v3.min.js" charset="utf-8"></script>
                    <script type="application/json" id="data">
                    '''
        self.body = '''
                    </script>
                    <script>

                    pi = Math.PI;
                    seqLength = circle

                    var degreeScale = d3.scale.linear()
                                          .domain([0, seqLength])
                                          .range([0,360])
                    ;

                    var data = JSON.parse(document.getElementById('data').innerHTML);

                    var arc = d3.svg.arc()
                        .innerRadius(function(d, i){return d.inner;})
                        .outerRadius(function(d, i){return d.outer;})
                        .startAngle(function(d, i){return degreeScale(d.start) * (pi/180);})
                        .endAngle(function(d, i){return degreeScale(d.end) * (pi/180);})
                    ;

                    var width = chart_width
                    var height = chart_height

                    var chart = d3.select("body").append("svg:svg")
                        .attr("width", width)
                        .attr("height", height)
                        .call(d3.behavior.zoom().on("zoom", function () {
                                chart.attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")")
                                }))
                        .append("svg:g")
                    ;

                    var div = d3.select("body").append("div")
                        .attr("class", "tooltip")
                        .style("opacity", 0)
                    ;

                    var ringShell = chart.append("g")
                                         .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

                    var textShell = chart.append("g")
                                         .attr("transform", "translate(" + width / 2 + "," + height / 2 + ")");

                    ringShell.selectAll("path")
                        .data(data)
                        .enter().append("svg:path")
                        .style("fill", function(d, i){ return d.color; })
                        .style("opacity", ring_opacity)
                        .attr("d", arc)
                        .attr('pointer-events', 'none')
                        .on('mouseover', function(d) {
                                div.transition()
                                    .duration(200)
                                    .style("opacity", .9);
                                div .html(d.text)
                                    .style("left", (d3.event.pageX + 20) + "px")
                                    .style("top", (d3.event.pageY + 10) + "px");
                                })
                        .on('mouseout', function(d) {
                                div.transition()
                                    .duration(200)
                                    .style("opacity", 0)
                            })
                        .attr('pointer-events', 'visible')
                    ;

                    textShell.append("text")
                      .style("opacity", 0)
                      .style("text-anchor", "middle")
                      .style("font-size", "title_size")
                      .style("font-weight", "bold")
                      .style("font-family", "title_font")
                      .attr("class", "inside")
                      .text(function(d) { return 'main_title'; })
                      .transition().duration(5000).style("opacity", 1);
                    ;

                </script>
                </body>

                '''

    def _setScript(self):
        """
        Replace placeholder values in script with given options.
        """
        for placeholder, value in self.options.items():
            self.body = self.body.replace(str(placeholder), str(value))

    def _writeHTML(self):
        """
        Write script to HTML."""
        with open(os.path.join(self.path, self.project + '.html'), 'w') as outfile:
            outfile.write(self.head)
        with open(os.path.join(self.path, self.project + '.html'), 'a') as outfile:
            json.dump(self.data, outfile, indent=4, sort_keys=True)
            outfile.write(self.body)
