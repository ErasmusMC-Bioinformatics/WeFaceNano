import time
import os
import random
import re
import json

from urllib.error import HTTPError
from django.conf import settings
from django.shortcuts import render_to_response, HttpResponseRedirect
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth import authenticate
from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render, redirect
from django.urls import reverse
from django.contrib.auth.models import User
from Bio import Entrez, SeqIO
from PRIMUL import brigD3
from subprocess import call, Popen, PIPE
from json import JSONDecodeError


@csrf_exempt
def signup(request):
    """Let superusers create a new users.

    Arguments:
        request: Request information used to sign up.
    """
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            if os.path.exists(settings.NANOPORE_DRIVE + username):
                pass
            else:
                call(["mkdir", settings.NANOPORE_DRIVE + username])
            if os.path.exists(settings.NANOPORE_DRIVE + username + "/results"):
                pass
            else:
                call(["mkdir", settings.NANOPORE_DRIVE + username + "/results"])
            return redirect('index')
    else:
        superusers = User.objects.filter(
            is_superuser=True).values_list('username')
        su_list = []
        for users in superusers:
            for u in users:
                su_list.append(u)
        if request.session.get("username") in su_list:
            superuser = True
            form = UserCreationForm()
            message = ""
        else:
            superuser = False
            form = None
            message = ("Only logged in superusers can create new users. "
                       "If you are a superuser, please log in.")
    return render(request, 'signup.html', context={
        'form': form,
        'user': request.session.get("username"),
        'super': superuser,
        'message': message})


@csrf_exempt
def index(request):
    """Check if user is superuser and retrieves the users results folder.
    Shows the homepage where you can select a result from the users results folder.

    Arguments:
        request: Request information when user is logged in. 
        Shows information about the users network drive 
        and if user is a superuser.
    """
    superusers = User.objects.filter(is_superuser=True).values_list('username')
    su_list = []
    for users in superusers:
        for u in users:
            su_list.append(u)
    if request.session.get("username") in su_list:
        superuser = True
    else:
        superuser = False
    if request.session.get("username"):
        network_drive = settings.NANOPORE_DRIVE + request.session.get("username")
        folders = os.walk(network_drive).__next__()[1]
        files = os.walk(network_drive).__next__()[2]
        results_folder = (settings.NANOPORE_DRIVE +
                          request.session.get("username") + "/results")
        results = os.walk(results_folder).__next__()[1]
    else:
        folders = []
        files = []
        results = []
    return render_to_response("home.html", context={
        "super": superuser,
        "user": request.session.get("username"),
        "folders": folders,
        "files": files,
        "drive": settings.NANOPORE_DRIVE,
        "results": results})


def draw_plasmid(contigfasta, contigname, genbank, refseq,
                 name, length, path, res_loc):
    """Draw a diagram conataining the plasmid with the aligned contigs.

    Arguments:
        contigfasta: A list of contigsfiles.
        contigname: A list of contig names based on the files.
        genbank: A list of genbank files to retrieve gene info.
        refseq: The reference filename to align the contigs with.
        name: A list of plasmid names.
        length: The length of the plasmid.
        path: The path to store te HTML output.
        res_loc: The names and locations of the 
        antibiotic resistance genes.
    """
    resistance_genes = {}
    for res, loc in res_loc.items():
        if int(loc[1]) <= int(length):
            resistance_genes[res] = [loc[0], loc[1]]
    if path != '':
        ring_gen = brigD3.AnnotationRing()
        ring_gen.setOptions(color='#000000',
                            name=name,
                            height=20,
                            res_loc=resistance_genes)
        ring_gen.readGenbank(genbank, length)
        genomes = contigfasta
        names = contigname
        blaster = brigD3.Blaster(refseq,
                                 genomes=genomes,
                                 path=path)
        blaster.runBLAST()
        blast_rings = []
        for i in range(len(blaster.results)):
            def r(): return random.randint(40, 120)
            color = ('#%02X%02X%02X' % (r(), r(), r()))
            ring_blast = brigD3.BlastRing()
            ring_blast.setOptions(color=color,
                                  name=names[i],
                                  res_loc=resistance_genes)
            ring_blast.min_length = 100
            ring_blast.readComparison(blaster.results[i])
            blast_rings.append(ring_blast)
        rings = [ring_gen] + blast_rings
        generator = brigD3.RingGenerator(rings, path, contigname)
        generator.setOptions(circle=length,
                             project=name,
                             title=name,
                             title_size='100%',
                             radius=200)
        generator.brigD3()


def kmergenie(resultfolder, inputfile):
    """Calculate the kmer-size that will be used when running minimap.
    
    Arguments:
        inputfile: The path to the nanopore reads to calculate the kmer-size
        resultfolder: Path to store the kmergenie output files.
    
    Returns:
        The calculated kmer-size.
    """
    kmer = "15"
    call(["mkdir", resultfolder + "/kmer"])
    kmergenie_out = Popen(
        ["~/kmergenie-1.7051/kmergenie " + inputfile +
            " -l 13 -k 21 -s 1 -o " + resultfolder + "/kmer/kmersizeoutput"],
        stdout=PIPE, shell=True).communicate()[0].decode()
    for line in kmergenie_out.split("\n"):
        if "best k:" in line:
            kmer = line.split(" ")[-1]
    return kmer


def flye(inputfolder, resultfolder, barcode_list, genomesize):
    """Use Flye as the assembler and return the assembly as a FASTA file.
    
    Arguments:
        inputfolder: Name of the folder with the inputfiles.
        resultfolder: Folder where the results will be stored.
        barcode_list: A list of barcodes.
        genomesize: The estimated genome size.
    
    Returns:
        A list of FASTA paths, FASTA filenames and barcodes.
    """
    call(["mkdir", resultfolder + "/assembly/"])
    file_list = []
    unitigs_barcode = []
    nodes = []
    for barcode in barcode_list:
        call(["mkdir", resultfolder + "/assembly/" + barcode])
        assemblyfolder = resultfolder + "/assembly/" + barcode
        call([
            "cat " + inputfolder + "/" + barcode +
            "/trimmed/* > " + inputfolder + "/" + barcode +
            "/trimmed/" + barcode + "_cat.fasta"
        ], shell=True)
        cat_file = inputfolder + "/" + barcode + "/trimmed/" + barcode + "_cat.fasta"
        file_to_use = inputfolder + "/" + barcode + "/trimmed/" + barcode + "_filtered.fasta"
        call(["awk '/^>/{f=!d[$1];d[$1]=1}f' " + cat_file + " > " + file_to_use], shell=True)
        call(["flye --nano-raw " + file_to_use + " -o " + assemblyfolder + " --genome-size " + genomesize], shell=True)
        files = os.listdir(resultfolder + "/assembly/" + barcode)
        for file in files:
            if "assembly.fasta" in file:
                with open(resultfolder + "/assembly/" + barcode + "/" + file) as contigfile:
                    for line in contigfile:
                        if ">" in line:
                            line = line[1:].split(" ")
                            nodes.append(line[0].replace("contig", "edge"))
                if nodes:
                    file_list.append(file)
                    unitigs_barcode.append(barcode)
            if file == "assembly_graph.gfa":
                for node in nodes:
                    call(
                        [
                            "Bandage image " + resultfolder + "/assembly/" + barcode + "/assembly_graph.gfa " +
                            resultfolder + "/assembly/" + barcode + "/" +
                            node.replace("edge", "contig").strip('\n') +
                            ".svg --height 200 --width 200 --iter 4 --colour random --scope aroundnodes --nodes " +
                            node.strip('\n') + " --nodewidth 3"
                        ], shell=True)
    return file_list, unitigs_barcode


def miniasm(inputtype, inputfolder, barcode_list, resultfolder, kmer, mincontig, circularise):
    """Use Miniasm for assembly and return the assembly as a FASTA file. 
    The Kmer-size will be calculated with KmerGenie so the user will not have to 
    manually enter this value. A minimum contig size can be entered to filter out
    all contigs smaller than this value. The default value will be 15000 when 
    no value is entered.

    Arguments:
        inputtype: To see if the inputfiles are FAST5 or FASTQ.
        inputfolder: Name of the folder with the inputfiles.
        barcode_list: A list of barcodes.
        resultfolder: Folder where the results will be stored.
        kmer: Kmer-size calculated with KmerGenie.
        mincontig: The minimum contig length. 
        All contig smaller will be filtered out.

    Returns:
        A list of FASTA paths, FASTA filenames and barcodes.
    """
    file_list = []
    unitigs_barcode = []
    if mincontig == "" or mincontig is None:
        mincontig = str(15000)
    if circularise is None:
        circularise = 0
    call(["mkdir", resultfolder + "/assembly/"])
    for barcode in barcode_list:
        if barcode != "unclassified":
            call(["mkdir", resultfolder + "/assembly/" + barcode])
            if inputtype == "fast5":
                barcode_content = os.listdir(
                    resultfolder + "/workspace/pass/" + barcode)
                if len(barcode_content) > 1:
                    kmer = kmergenie(resultfolder, resultfolder + "/workspace/pass/" + barcode + "/" + barcode + "_cat.fastq")
                    call([
                        "cat " + resultfolder + "/workspace/pass/" + barcode +
                        "/* > " + resultfolder + "/workspace/pass/" + barcode +
                        "/" + barcode + "_cat.fastq"
                    ], shell=True)
                    call([
                        "bash ~/PRIMUL/static/miniasm.bash -v2 -m " + mincontig + " -k " + kmer + " -c " + 
                        str(circularise) + " -i " + resultfolder + "/workspace/pass/" + barcode + "/" +
                        barcode + "_cat.fastq -o " + resultfolder + "/assembly/" + barcode + "/" + "_cat.contigs"
                    ], shell=True)
                else:
                    kmer = kmergenie(resultfolder, resultfolder + "/workspace/pass/" + barcode + "/" + barcode_content[0])
                    call([
                        "bash ~/PRIMUL/static/miniasm.bash -v2 -m " + mincontig + " -k " + kmer + " -c " + 
                        str(circularise) + " -i " + resultfolder + "/workspace/pass/" + barcode + "/" +
                        barcode_content[0] + " -o " + resultfolder + "/assembly/" + barcode + "/" + barcode_content[0] + ".contigs"
                    ], shell=True)
            else:
                fastacount = 0
                barcode_content = os.listdir(inputfolder + "/" + barcode)
                for fasta in barcode_content:
                    if os.path.isfile(inputfolder + "/" + barcode + "/" + fasta):
                        fastacount += 1
                if len(barcode_content) > 1:
                    call([
                        "cat " + inputfolder + "/" + barcode +
                        "/trimmed/* > " + inputfolder + "/" + barcode +
                        "/trimmed/" + barcode + "_cat.fasta"
                    ], shell=True)
                    kmer = kmergenie(resultfolder, str(inputfolder + "/" + barcode + "/trimmed/" + barcode + "_cat.fasta"))
                    call([
                        "bash ~/PRIMUL/static/miniasm.bash -v2 -m " + mincontig + " -k " + kmer + " -c " + 
                        str(circularise) + " -i " + inputfolder + "/" + barcode + "/trimmed/" + barcode +
                        "_cat.fasta -o " + resultfolder + "/assembly/" + barcode + "/" + barcode + "_cat.contigs"
                    ], shell=True)
                    call([
                        "rm", "-r", inputfolder + "/" + barcode +
                        "/trimmed/" + barcode + "_cat.fasta"
                    ])
                else:
                    kmer = kmergenie(resultfolder, inputfolder + "/" + barcode + "/" + barcode_content[0])
                    call([
                        "bash ~/PRIMUL/static/miniasm.bash -v2 -m " + mincontig + " -k " + kmer + " -c " + 
                        str(circularise) + " -i " + inputfolder + "/" + barcode + "/trimmed/" +
                        barcode_content[0] + " -o " + resultfolder + "/assembly/" + barcode + "/" + 
                        barcode_content[0] + ".contigs"
                    ], shell=True)
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file and ".jpg" not in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" +
                              file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    return file_list, unitigs_barcode


def skip_assembly(resultfolder, barcode_list):
    """Skip the assembly if the assembly folder is already available.

    Arguments:
        resultfolder: Folder where the results will be stored.
        barcode_list: A list of barcodes.

    Returns:
        A list of FASTA paths, FASTA filenames and barcodes.
    """
    file_list = []
    unitigs_barcode = []
    for barcode in barcode_list:
        if barcode != "unclassified":
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file and ".jpg" not in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" +
                              file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    return file_list, unitigs_barcode


def plasmidfinder(barcodes, file_list, resultfolder, res_loc):
    """Run plasmidfinder on all FASTA files in the list on the 
    plasmidfinder enterobacteriaceae database.
    
    Arguments:
        barcodes: A list of barcodes
        file_list: A list of files that contains the assembly.
        resultfolder: The ouputfolder to store the results.
    """
    plasmidfinder_dict = {}
    call(["mkdir", resultfolder + "/plasmidfinder"])
    count = 0
    for assemblyfile in file_list:
        call(["mkdir", resultfolder + "/plasmidfinder/" + barcodes[count]])
        inpath = resultfolder + "/assembly/" + \
            barcodes[count] + "/" + assemblyfile
        outpath = str(resultfolder + "/plasmidfinder/" + barcodes[count] + "/")
        plasmidcmd = "python3 ~/plasmidfinder/plasmidfinder.py -i " + inpath + \
            " -o " + outpath + " -p ~/plasmidfinder_db/ -d enterobacteriaceae"
        call([plasmidcmd], shell=True)
        count += 1
        with open(outpath + "data.json", "r") as plasmidfinder:
            try:
                load = json.loads(plasmidfinder.read())
                enterobacteriaceae = load["plasmidfinder"]["results"]["Enterobacteriaceae"]["enterobacteriaceae"]
                for inc in enterobacteriaceae:
                    contig = enterobacteriaceae[inc]["contig_name"].split(" ")[0]
                    if contig in plasmidfinder_dict.keys():
                        res_loc[contig + "_" + enterobacteriaceae[inc]["plasmid"] + "_" + str(enterobacteriaceae[inc]["identity"]) + "_Inc"] = enterobacteriaceae[inc]["positions_in_contig"].split("..")
                    else:
                        res_loc[contig + "_" + enterobacteriaceae[inc]["plasmid"] + "_" + str(enterobacteriaceae[inc]["identity"]) + "_Inc"] = enterobacteriaceae[inc]["positions_in_contig"].split("..")
            except JSONDecodeError:
                res_loc[""] = "No genes found"
    return res_loc


def resfinder(barcodes, file_list, resultfolder,
              resdb, reslength, residentity):
    """Runs resfinder on all FASTA files in the list on the 
    selected database or all databases.

    Arguments:
        barcodes: A list of barcodes.
        file_list: A list of files that contains the assembly.
        resultfolder: The ouputfolder to store the results.
        resdb: The database to search antibiotic resistance genes.
        reslength: The minimum length cutoff when using resfinder.
        residentity: The % of identity cutoff when using resfinder.

    Returns:
        A dictionary with the contigs and the found genes and 
        the gene locations within the contigs.
    """
    call(["mkdir", resultfolder + "/resfinder"])
    db_all = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin",
              "fusidicacid", "macrolide", "nitroimidazole", "oxazolidinone",
              "phenicol", "quinolone", "rifampicin", "sulphonamide",
              "tetracycline", "trimethoprim", "glycopeptide"]
    res_genes = {}
    res_loc = {}
    count = 0
    for file in file_list:
        call(["mkdir", resultfolder + "/resfinder/" + barcodes[count]])
        call([
            "bash ~/PRIMUL/static/splitfasta.sh " + resultfolder +
            "/assembly/" + barcodes[count] + "/" + file + " " + resultfolder +
            "/resfinder/" + barcodes[count]
        ], shell=True)
        unitig_bc = os.listdir(resultfolder + "/resfinder/" + barcodes[count])
        for unitig in unitig_bc:
            fasta_unitig = (resultfolder + "/resfinder/" +
                            barcodes[count] + "/" + unitig)
            if resdb != "all":
                call([
                    "perl ~/PRIMUL/static/resfinder.pl "
                    "-d ~/resfinder_db -i " + fasta_unitig +
                    " -a " + resdb + " -o " + resultfolder + "/resfinder/" +
                    barcodes[count] + " -k " + residentity + " -l " + reslength
                ], shell=True)
            else:
                for db in db_all:
                    call([
                        "perl ~/PRIMUL/static/resfinder.pl -d "
                        "~/resfinder_db -i " + fasta_unitig +
                        " -a " + db + " -o " + resultfolder + "/resfinder/" +
                        barcodes[count] + " -k " + residentity +
                        " -l " + reslength
                    ], shell=True)
            argene = []
            with open(resultfolder + "/resfinder/" + barcodes[count] +
                      "/results_tab.txt") as rf:
                stored_contig = ""
                rcount = 0
                for line in rf:
                    line = line.split("\t")
                    if rcount == 0:
                        stored_contig = line[3]
                    contig = line[3]
                    if contig == stored_contig:
                        argene.append(line[0] + " - " + line[5])
                    else:
                        argene = []
                        argene.append(line[0] + " - " + line[5])
                    res_genes[barcodes[count] + "_" + contig] = argene
                    loc = line[4].split("..")
                    res_loc[contig + "_" + line[0] + "_" + str(rcount) + "_" +
                            line[1] + "_" + line[5]] = [loc[0], loc[1]]
                    stored_contig = line[3]
                    rcount += 1
        count += 1
    return res_genes, res_loc


def run_blast(barcodes, file_list, resultfolder, blastdb, task, res_loc):
    """Runs BLAST on assembled FASTA files.

    FIXME: Function is too long. Split or clean function.

    Arguments:
        barcodes: A list of barcodes that contains contigs.
        file_list: A list of filenames to use with BLAST.
        resultfolder: The folder where the 
        BLAST results will be stored.
        blastdb: The selected BLAST database.
        task: The selected task in BLAST (blastn or megablast)
        res_loc: Locations of resistance genes to draw plasmid.

    Returns:
        Dictionary with BLAST output information (plasmid, description, length).

    Raises:
        HTTPError: No Genbank entry found.
        IndexError, KeyError: No BLAST result.
    """
    blast = {}
    bfiles = []
    call(["mkdir", resultfolder + "/BLAST"])
    blastfolder = resultfolder + "/BLAST"
    bcount = 0
    dict_contigfasta = {}
    dict_draw = {}
    plasmidcount = 0
    ext_plasmidcount = 0
    Entrez.email = "some_email@somedomain.com"
    for fasta in file_list:
        unitig_bc = []
        call(["mkdir", resultfolder + "/BLAST/" + barcodes[bcount]])
        call([
            "bash ~/PRIMUL/static/splitfasta.sh " + resultfolder +
            "/assembly/" + barcodes[bcount] + "/" + fasta + " " +
            blastfolder + "/" + barcodes[bcount]
        ], shell=True)
        raw_unitig_bc = os.listdir(blastfolder + "/" + barcodes[bcount])
        for rub in raw_unitig_bc:
            if ".out" not in rub:
                unitig_bc.append(rub)
        if "unclassified" in os.listdir(blastfolder):
            for utig in os.listdir(blastfolder + "/unclassified/"):
                call([
                    "blastn -db " + blastdb +
                    " -query " + blastfolder + "/unclassified/" + utig +
                    " -out " + blastfolder + "/unclassified/" + utig + ".out "
                    "-max_target_seqs 1 -outfmt 6 -task " + task
                ], shell=True)
            for unc in os.listdir(blastfolder + "/unclassified/"):
                if ".out" in unc:
                    bfiles.append(blastfolder + "/unclassified/" + unc)
        else:
            for utig in unitig_bc:
                call([
                    "blastn -db " + blastdb +
                    " -query " + blastfolder + "/" + barcodes[bcount] + "/" +
                    utig + " -out " + blastfolder + "/" + barcodes[bcount] +
                    "/" + utig + ".out -max_target_seqs 1 -outfmt 6 "
                    "-task " + task
                ], shell=True)
                time.sleep(10)
                with open(blastfolder + "/" + barcodes[bcount] + "/" +
                          utig + ".out") as blastfile:
                    with open(blastfolder + "/" + barcodes[bcount] +
                              "/" + utig) as contigfile:
                        cheader = contigfile.readline()
                        cheader = cheader.split(" ")
                        contigname = barcodes[bcount] + "_" + cheader[0][1:]
                        for record in SeqIO.parse(contigfile.name, "fasta"):
                            bclength = len(record.seq)
                        line = blastfile.readline()
                    if line != "":
                        line = line.split("\t")
                        try:
                            handle = Entrez.efetch(
                                db="nucleotide",
                                id=line[1].split("_")[-1],
                                rettype="gb",
                                retmode="gb")
                            handle_txt = Entrez.efetch(
                                db="nucleotide",
                                id=line[1].split("_")[-1],
                                rettype="gb",
                                retmode="text")
                            gbrecord = SeqIO.read(handle, "genbank")
                            gblength = len(gbrecord.seq)
                            with open(blastfolder + "/" + barcodes[bcount] + "/" +
                                      utig[:-6] + ".refseq.fasta", "w") as ref:
                                ref.write(gbrecord.format("fasta"))
                            with open(blastfolder + "/" + barcodes[bcount] +
                                      "/" + utig[:-6] +
                                      ".genbank.gb", "w") as gb:
                                gb.write(handle_txt.read())
                            genbank = (
                                blastfolder + "/" + barcodes[bcount] + "/" +
                                utig[:-6] + ".genbank.gb"
                            )
                            bplength = gblength
                        except HTTPError:
                            genbank = ""
                        new_path = blastfolder + "/" + barcodes[bcount] + "/"
                        if blastdb == settings.NANOPORE_DRIVE + "plasmidb/plasmidb":
                            for ref_folder in os.listdir(
                                    settings.NANOPORE_DRIVE):
                                if "RB_REF" in ref_folder:
                                   for ref_file in os.listdir(
                                           settings.NANOPORE_DRIVE +
                                           ref_folder):
                                       if line[1].strip("/") in ref_file:
                                           with open(settings.NANOPORE_DRIVE +
                                                     ref_folder + "/" +
                                                     ref_file) as ref_rb_file:
                                               with open(blastfolder + "/" +
                                                         barcodes[bcount] + 
                                                         "/" + utig[:-6] + 
                                                         ".refseq.fasta", 
                                                         "w") as ref:
                                                   ref.write(
                                                       ref_rb_file.read())
                        refseq = (
                            blastfolder + "/" + barcodes[bcount] + "/" +
                            utig[:-6] + ".refseq.fasta"
                        )
                        try:
                            for record in SeqIO.parse(refseq, "fasta"):
                                bplength = len(record)
                            try:
                                for feature in record.features:
                                    if feature.qualifiers.get('plasmid', []):
                                        plasmids = feature.qualifiers.get(
                                            'plasmid')
                                for plasmid in plasmids:
                                    time.sleep(10)
                                    p = plasmid
                                    dict_contigfasta[p + "_" +
                                                    str(plasmidcount)] = utig
                                    dict_draw[p + "_" + str(plasmidcount) + "_" + barcodes[bcount]] = [
                                        contigname,
                                        genbank,
                                        refseq,
                                        bplength,
                                        new_path]
                                    plasmidcount += 1
                                blast[barcodes[bcount] + "_" + line[0]] = [
                                    p, record.description, bplength, bclength]
                            except:
                                ext_plasmidcount += 1
                                time.sleep(10)
                                dict_contigfasta[line[1].split("_")[-1] + "_" +
                                                str(ext_plasmidcount)] = utig
                                dict_draw[line[1].split("_")[-1] + "_" +
                                        str(ext_plasmidcount) + "_" + barcodes[bcount]] = [
                                    contigname,
                                    genbank,
                                    refseq,
                                    bplength,
                                    new_path]
                                try:
                                    blast[barcodes[bcount] + "_" + line[0]] = [
                                        line[1],
                                        record.description,
                                        bplength,
                                        bclength]
                                except:
                                    pass
                        except FileNotFoundError:
                            pass
                    else:
                        blast[contigname] = [
                            "No Accession Number", "No Name", "0", bclength]
        # After running BLAST
        contigkeys = dict_contigfasta.keys()
        usedkeys = []
        for key in contigkeys:
            utgkeys = []
            for dplasmid, dutigs in dict_contigfasta.items():
                if key.split("_")[0] in dplasmid:
                    utgkeys.append(dutigs)
                else:
                    pass
            if key.split("_")[0] not in usedkeys:
                try:
                    contigs = []
                    name = key.split("_")[0]
                    for k, v in dict_draw.items():
                            if name in k:
                                contigs.append(v[0])
                                if barcodes[bcount] in k:
                                    new_path_bc = v[4]
                    draw_plasmid(
                        contigfasta=utgkeys,
                        contigname=contigs,
                        genbank=dict_draw[key + "_" + barcodes[bcount]][1],
                        refseq=dict_draw[key + "_" + barcodes[bcount]][2],
                        name=name,
                        length=dict_draw[key + "_" + barcodes[bcount]][3],
                        path=new_path_bc,
                        res_loc=res_loc)
                except (IndexError, KeyError):
                    pass
            else:
                pass
            usedkeys.append(key.split("_")[0])
        barcodes.pop(0)
        dict_contigfasta = {}
        dict_draw = {}
    return blast


@csrf_exempt
def create_results(request):
    """Run the pipeline with entered data and show the results in an html page.

    Arguments:
        request: Request information to get information to run the pipeline.
    """
    inputfolder = request.POST.get("inputfolder")
    inputtype = request.POST.get("inputtype")
    configuration = request.POST.get("aconfig")
    outfolder = request.POST.get("albacoreout").replace(" ", "_")
    res = request.POST.get("res")
    mlst = request.POST.get("mlst")
    species = request.POST.get("species")
    assembler = request.POST.get("assembler")
    blast = request.POST.get("blast")
    blastdb = request.POST.get("blastdb")
    blasttask = request.POST.get("blasttask")
    resdb = request.POST.get("resdb")
    reslength = request.POST.get("reslength")
    residentity = request.POST.get("residentity")
    kmer = request.POST.get("kmer")
    mincontig = request.POST.get("min-contig")
    genomesize = request.POST.get("gsize")
    circularise = request.POST.get("circularise")
    resultfolder = (settings.NANOPORE_DRIVE + request.session.get("username") +
                    "/results/" + outfolder)
    if res is None:
        res = ''
    if mlst is None:
        mlst = ''
    if species is None:
        species = ''
    request.session['stored_results'] = request.POST
    call(["mkdir", resultfolder])
    if request.POST.get("barcoding") == '1' and inputtype == "fast5":
        call(["mkdir", resultfolder + "/workspace/pass"])
        call(["mkdir", resultfolder + "/qc"])
        call([
            "bash ~/PRIMUL/static/albacore.sh " + res + mlst + species +
            "-b -c " + configuration + " -i" + inputfolder + "/fast5 -o " +
            resultfolder + " -f fastq"
        ], shell=True)
    elif request.POST.get("barcoding") == '0' and inputtype == "fast5":
        call(["mkdir", resultfolder + "/workspace/pass"])
        call(["mkdir", resultfolder + "/qc"])
        call([
            "bash ~/PRIMUL/static/albacore.sh " + res + mlst + species +
            "-c " + configuration + " -i" + inputfolder + "/fast5 -o " +
            resultfolder + " -f fastq"
        ], shell=True)
    if inputtype == "fast5":
        barcode_list = []
        time.sleep(15)
        while True:
            if os.path.isfile(resultfolder + "/finished.dat"):
                break
        albacore_list = os.listdir(resultfolder + "/workspace/pass")
        for al in albacore_list:
            if os.path.isdir(resultfolder + "/workspace/pass/" + al):
                barcode_list.append(al)
        if not barcode_list:
            barcode_list = ["fast5"]
    else:
        if os.path.isdir(resultfolder + "/assembly"):
            barcode_list = os.listdir(resultfolder + "/assembly")
            for bc in os.listdir(inputfolder):
                if not os.path.isdir(inputfolder + "/" + bc + "/trimmed"):
                    trim_reads(inputfolder, barcode_list)
        else:
            barcode_list = os.listdir(inputfolder)
            for bc in barcode_list:
                if not os.path.isdir(inputfolder + "/" + bc + "/trimmed"):
                    trim_reads(inputfolder, barcode_list)
    if os.path.exists(resultfolder + "/assembly/"):
        file_list, barcodes = skip_assembly(
            resultfolder, barcode_list)
    else:
        if assembler == "miniasm":
            file_list, barcodes = miniasm(
                inputtype, inputfolder, barcode_list, resultfolder,
                kmer, mincontig, circularise
            )
        elif assembler == "flye":
            file_list, barcodes = flye(inputfolder, resultfolder, barcode_list, genomesize)
        else:
            return HttpResponseRedirect(reverse("index"))
    # Create QC pages
    for barcode in barcode_list:
        files = os.listdir(inputfolder + "/" + barcode + "/trimmed")
        for file in files:
            file = inputfolder + "/" + barcode + "/trimmed/" + file
            if file.split(".")[1] != "fast5":
                extention = "fasta"
            else:
                extention = "fast5"
            qc_report(extention, barcode, file,
                      resultfolder + "/qc/" + barcode)
    # plasmidfinder(barcodes, file_list, resultfolder)
    if res == '-r ':
        dummyres_genes, res_loc = resfinder(
            barcodes, file_list, resultfolder, resdb, reslength, residentity)
    else:
        res_loc = {}
    res_loc = plasmidfinder(barcodes, file_list, resultfolder, res_loc)
    if blast == "-b ":
        run_blast(
            barcodes, file_list, resultfolder, blastdb, blasttask, res_loc)
    for bc in barcode_list:
        call(["rm", "-rf", inputfolder + "/" + bc + "/trimmed"])
    return HttpResponseRedirect(reverse("index"))


def qc_report(extention, barcode, filename, outfolder):
    """Runs the NanoPot QC tool on the FASTA or FASTQ file. Generates
    an HTML page with QC stats and read distribution plots. When using a
    FASTQ file quality score information will be shown.
    
    Arguments:
        extention: File extention (fasta or fastq).
        barcode: Barcode name.
        filename: Name of the file to run the QC.
        outfolder: Where to store the QC data.
    """
    cmd = ("NanoPlot --" + extention + " " + filename +
           " --N50 --plot kde hex dot -p " + barcode +
           " -o " + outfolder)
    call([cmd], shell=True)


def trim_reads(inputfolder, barcode_list):
    """Trim the adapters from the reads using Porechop.

    Arguments:
        inputfolder: The folder containing the reads.
        barcode_list: A list of barcodes to find the FASTQ files.
    """
    for barcode in barcode_list:
        barcode_folder = os.listdir(inputfolder + "/" + barcode)
        trimmed_folder = inputfolder + "/" + barcode + "/trimmed/"
        if "trimmed" not in barcode_folder:
            for file in barcode_folder:
                file_path = inputfolder + "/" + barcode + "/" + file
                call(["mkdir", trimmed_folder])
                trim_path = trimmed_folder + file
                call(["porechop", "-i", file_path, "-o", trim_path, "--format", "fasta"])
        else:
            pass


def readme(request):
    """Shows the readme page and checks if user is a superuser.

    Arguments:
        request: A request to check if the logged in user is a superuser.
    """
    superusers = User.objects.filter(is_superuser=True).values_list('username')
    su_list = []
    for users in superusers:
        for u in users:
            su_list.append(u)
    if request.session.get("username") in su_list:
        superuser = True
    else:
        superuser = False
    return render_to_response("readme.html", context={
        "user": request.session.get("username"),
        "super": superuser})


@csrf_exempt
def pipeline_start(request):
    """Shows the pipeline page and all the available options.

    Arguments:
        request: Request information to start the plasmid pipeline.

    Raises:
        StopIteration: No Nanopore drive found.
    """
    resdb = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin",
             "fusidicacid", "macrolide", "nitroimidazole", "oxazolidinone",
             "phenicol", "quinolone", "rifampicin", "sulphonamide",
             "tetracycline", "trimethoprim", "glycopeptide"]
    if request.session.get('username') is None:
        return HttpResponseRedirect(reverse("index"))
    else:
        superusers = User.objects.filter(
            is_superuser=True).values_list('username')
        su_list = []
        for users in superusers:
            for u in users:
                su_list.append(u)
        if request.session.get("username") in su_list:
            superuser = True
        else:
            superuser = False
        network_drive = (settings.NANOPORE_DRIVE +
                         request.session.get("username") + "/")
        folders = []
        try:
            blastdb = []
            blastdbfolder = os.listdir(settings.NANOPORE_DRIVE + "blastdb/")
            for db in blastdbfolder:
                if os.path.isdir(settings.NANOPORE_DRIVE + "blastdb/" + db):
                    blastdb.append(db)
            walk = os.walk(network_drive).__next__()[1]
            for folder in walk:
                if "results" not in folder:
                    folders.append(network_drive + folder)
        except StopIteration:
            return HttpResponseRedirect(reverse('index'))
        return render_to_response("pipeline.html", context={
            "folders": folders,
            "user": request.session.get('username'),
            "super": superuser,
            "blastdb": blastdb,
            "drive": settings.NANOPORE_DRIVE,
            "resdb": resdb})


@csrf_exempt
def get_stored_results(request):
    """Gets stored results from the network drive after the user selects a 
    previously run analysis from the homepage.

    Arguments:
        request: Request information to retrieve stored results.
    """
    superusers = User.objects.filter(is_superuser=True).values_list('username')
    su_list = []
    blast_dict = {}
    html_plasmid = {}
    blast_res_dict = {}
    resfinder_dict = {}
    assembly_report = []
    topology = {}
    tools = []
    barcodes = []
    Entrez.email = "some_email@somedomain.com"
    for users in superusers:
        for u in users:
            su_list.append(u)
    if request.POST.get("username") in su_list:
        superuser = True
    else:
        superuser = False
    username = request.POST.get("username")
    results = os.listdir(settings.NANOPORE_DRIVE + username + "/results/")
    res_dict = {}
    selected_result = request.POST.get("resultname")
    for r in results:
        if selected_result != "none":
            if r == selected_result:
                tool_list = []
                tools = os.listdir(settings.NANOPORE_DRIVE +
                                   username + "/results/" + r)
                if tools:
                    for t in tools:
                        if "." not in t:
                            tool_list.append(t)
                            res_dict[r] = tool_list
                        if t == "BLAST":
                            blast_dict = get_stored_blast_results(
                                username, r)[0]
                            blast_res_dict = get_stored_blast_results(
                                username, r)[1]
                            html_plasmid = get_stored_blast_results(
                                username, r)[2]
                            topology = get_stored_blast_results(username, r)[3]
                        if t == "resfinder":
                            resfinder_dict = get_stored_resfinder_results(
                                username, r)
                        if t == "plasmidfinder":
                            plasmidfinder_dict = get_stored_plasmidfinder_results(
                                username, r)
                        if t == "assembly":
                            assembly_report = get_stored_assembly_results(
                                username, r)[0]
                            barcodes = get_stored_assembly_results(
                                username, r)[2]
                            graphs = get_stored_assembly_results(
                                username, r)[3]
                        if t == "qc":
                            qc_html = get_stored_qc(username, r)
                else:
                    err = "Results folder is empty"
                    return render_to_response("home.html", context={
                        "super": superuser,
                        "user": username,
                        'error': err})
    return render_to_response("results.html", context={
        "super": superuser,
        "user": username,
        "results": results,
        "tools": tools,
        "dict": res_dict,
        "blast": blast_dict,
        "plasmid": html_plasmid,
        "resfinder": resfinder_dict,
        "plasmidfinder": plasmidfinder_dict,
        "blastresults": blast_res_dict,
        "assemblyreport": assembly_report,
        "barcodes": barcodes,
        "topology": topology,
        "qc": qc_html,
        "graphs": graphs})


def get_stored_qc(username, r):
    """Get the QC pages stored on the drive.
    
    Arguments:
        username: The user that is logged in.
        r: Name of the selected run.
    
    Returns:
        Dictionary with qc html page.
    """
    qc_dict = {}
    qc_bc = (settings.NANOPORE_DRIVE + username +
             "/results/" + r + "/qc/")
    for bc in os.listdir(qc_bc):
        folder = (settings.NANOPORE_DRIVE + username +
                  "/results/" + r + "/qc/" + bc)
        for qc in os.listdir(folder):
            if ".html" in qc:
                with open(folder + "/" + qc) as qcpage:
                    qc_html = qcpage.read()
                qc_dict[bc] = qc_html
    return qc_dict


def get_stored_blast_results(username, r):
    """Get the stored BLAST result based on the selected run.

    Arguments:
        username: The user that is logged in.
        r: Name of the selected run

    Returns:
        A dictionary with BLAST files, BLAST result information and 
        the HTML code to view the plasmid.

    Raises:
        HTTPError: No Genbank information.
        IndexError: No BLAST data available.
    """
    blast_dict = {}
    html_plasmid = {}
    blast_res_dict = {}
    topology = {}
    coveragelist = []
    total_length = []
    blast_results = os.listdir(
        settings.NANOPORE_DRIVE + username + "/results/" + r + "/BLAST/")
    for br in blast_results:
        bc_blast = os.listdir(settings.NANOPORE_DRIVE +
                              username + "/results/" + r + "/BLAST/" + br)
        for bc in bc_blast:
            if ".html" in bc:
                with open(settings.NANOPORE_DRIVE + username + "/results/" + r +
                          "/BLAST/" + br + "/" + bc) as htmlfile:
                    code = htmlfile.read()
                html_plasmid[br + "_" + bc.split(".")[0]] = code
            if ".out" in bc:
                count = 0
                with open(settings.NANOPORE_DRIVE + username + "/results/" + r +
                          "/BLAST/" + br + "/" + bc) as blastfile:
                    for line in blastfile:
                        line = line.split("\t")
                        if count == 0:
                            first_hit = line[1]
                        if line[1] == first_hit:
                            coveragelist.append(float(line[2]))
                            total_length.append(int(line[3]))
                        count += 1
                    blastfile.seek(0)
                    first_result = blastfile.readline()
                    first_result = first_result.split('\t')
                    try:
                        contig = first_result[0]
                        try:
                            handle = Entrez.efetch(
                                db="nucleotide",
                                id=first_result[1].split("_")[-1],
                                rettype="gb",
                                retmode="gb")
                            record = SeqIO.read(handle, "genbank")
                            topology[record.description] = record.annotations[
                                "topology"]
                            blast_name = record.description
                        except HTTPError:
                            topology[first_result[1]] = "circular"
                            blast_name = first_result[1]
                    except IndexError:
                        blast_name = ""
                        contig = ""
                if contig:
                    blast_res_dict[br + "_" + contig] = [blast_name]
                else:
                    pass
    blast_dict[r] = blast_results
    return blast_dict, blast_res_dict, html_plasmid, topology


def get_stored_resfinder_results(username, r):
    """Get stored resfinder results based on the selected run.

    Arguments:
        username: The user that is logged in.
        r: Name of the selected run.

    Returns:
        A dictionary with resfinder results 
        (Contig name, antibiotic resistance genes and percentag identity).
    
    Raises:
        IndexError: No resfinder information available.
    """
    resfinder_dict = {}
    arglist = []
    resfinder_results = os.listdir(
        settings.NANOPORE_DRIVE + username + "/results/" + r + "/resfinder/")
    for resfolder in resfinder_results:
        contiglist = []
        resbarcode = os.listdir(
            settings.NANOPORE_DRIVE + username + "/results/" + r +
            "/resfinder/" + resfolder)
        if "results_tab.txt" in resbarcode:
            if os.stat(
                settings.NANOPORE_DRIVE + username + "/results/" + r +
                    "/resfinder/" + resfolder + "/results_tab.txt"
            ).st_size > 0:
                with open(settings.NANOPORE_DRIVE + username + "/results/" + r +
                          "/resfinder/" + resfolder +
                          "/results_tab.txt") as resfinderfile:
                    for line in resfinderfile:
                        line = line.split("\t")
                        if line[3] not in contiglist:
                            if contiglist:
                                resfinder_dict[resfolder + "_" +
                                               contiglist[-1]] = arglist
                            arglist = []
                            contiglist.append(line[3])
                        arglist.append(
                            str(
                                line[0] + " - (" + line[1] + "%) - " + line[5])
                        )
                    try:
                        resfinder_dict[resfolder + "_" +
                                       contiglist[-1]] = arglist
                    except IndexError:
                        pass
            else:
                pass
    return resfinder_dict


def get_stored_plasmidfinder_results(username, r):
    """Get stored plasmidfinder results based on the selected run.
    
    Arguments:
        username: The user that is logged in.
        r: Name of the selected run.
    
    Returns:
        A dictionary with contigs and gene names.

    Raises:
        JSONDecodeError: Add no genes found to the plasmidfinder dictionary
        if a JSONDecodeError occurs.
    """
    plasmidfinder_dict = {}
    plasmidfinder_results = (os.listdir(
        settings.NANOPORE_DRIVE + username + "/results/" + r + "/plasmidfinder/"))
    for bc in plasmidfinder_results:
        plasmidfile = (settings.NANOPORE_DRIVE + username +
                       "/results/" + r + "/plasmidfinder/" + bc + "/data.json")
        with open(plasmidfile, "r") as plasmidfinder:
            try:
                load = json.loads(plasmidfinder.read())
                enterobacteriaceae = load["plasmidfinder"]["results"]["Enterobacteriaceae"]["enterobacteriaceae"]
                for inc in enterobacteriaceae:
                    contig = (
                        bc + "_" + enterobacteriaceae[inc]["contig_name"].split(" ")[0])
                    if contig in plasmidfinder_dict.keys():
                        plasmidfinder_dict[contig].append(
                            enterobacteriaceae[inc]["plasmid"] + " - (" +
                            str(enterobacteriaceae[inc]["identity"]) + "%) - " +
                            "<b><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/" +
                            enterobacteriaceae[inc]["accession"] + "\" target=\"_blank\">" +
                            enterobacteriaceae[inc]["accession"] + "</a></b>")
                    else:
                        plasmidfinder_dict[contig] = [
                            enterobacteriaceae[inc]["plasmid"] + " - (" +
                            str(enterobacteriaceae[inc]["identity"]) + "%) - " +
                            "<b><a href=\"https://www.ncbi.nlm.nih.gov/nuccore/" +
                            enterobacteriaceae[inc]["accession"] + "\" target=\"_blank\">" +
                            enterobacteriaceae[inc]["accession"] + "</a></b>"
                        ]
            except JSONDecodeError:
                plasmidfinder_dict[""] = "No genes found"
    return plasmidfinder_dict


def get_stored_assembly_results(username, r):
    """Get stored assembly results based on the selected run.

    Arguments:
        username: The user that is logged in.
        r: Name of the selected run.

    Returns:
        A dictionary with the assembly and contig lengths.
        A list of files containing the assembly results and barcodes.
    """
    graphs = {}
    assembly_report = {}
    assembly_results = os.listdir(
        settings.NANOPORE_DRIVE + username + "/results/" + r + "/assembly/")
    assembly_bc = []
    for assemblyfolder in assembly_results:
        assembly_barcode = os.listdir(
            settings.NANOPORE_DRIVE + username + "/results/" + r +
            "/assembly/" + assemblyfolder)
        for assembly in assembly_barcode:
            if ".contigs.fasta" in assembly or "assembly.fasta" in assembly: 
                if ".jpg" not in assembly:
                    contigcount = 0
                    if os.stat(
                        settings.NANOPORE_DRIVE + username + "/results/" + r +
                        "/assembly/" + assemblyfolder + "/" + assembly
                    ).st_size > 0:
                        assembly_bc.append(assemblyfolder)
                        for record in SeqIO.parse(
                            settings.NANOPORE_DRIVE + username + "/results/" + r +
                                "/assembly/" + assemblyfolder + "/" + assembly,
                                "fasta"):
                            contigcount += 1
                            contiglength = len(record.seq)
                            assembly_report[assemblyfolder +
                                            "_" + record.id] = [record.id, contiglength]
            if ".svg" in assembly:
                with open(settings.NANOPORE_DRIVE + username + "/results/" + r + "/assembly/" + assemblyfolder + "/" + assembly) as svg:
                    graphs[assembly[:-4]] = svg.read()
    return assembly_report, assembly_results, sorted(assembly_bc), graphs


def delete(request):
    """Checks if user is superuser and gets the results created by the user.
    Shows a dropdown menu on the delete page to select a result the user wants 
    to delete.
    
    Arguments:
        request: Request information to visit the delete page.
    """
    superusers = User.objects.filter(is_superuser=True).values_list('username')
    su_list = []
    for users in superusers:
        for u in users:
            su_list.append(u)
    if request.session.get("username") in su_list:
        superuser = True
    else:
        superuser = False
    username = request.session.get('username')
    results = os.listdir(settings.NANOPORE_DRIVE + username + "/results/")
    return render_to_response("delete.html", context={
        'super': superuser,
        'user': username,
        'results': results})


@csrf_exempt
def remove_result(request):
    """Remove the selected result from the users results folder on the
    network drive. This cannot be undone.
    
    Arguments:
        request: Get the username and resultname to remove the 
        selected result from the network drive.
    """
    result = request.POST.get("resultname")
    username = request.POST.get("username")
    cmd = settings.NANOPORE_DRIVE + username + "/results/" + result
    call(["rm", "-rf", cmd])
    return HttpResponseRedirect(reverse("index"))


def logout(request):
    """Flush the session to logout the current user.

    Arguments:
        request: Request the session of the current user to 
        flush the session.
    """
    request.session.flush()
    return HttpResponseRedirect(reverse("index"))


@csrf_exempt
def login(request):
    """Authenticate user. 
    Check if user is in the mysql database and password is correct.
    Store username in a session.

    Arguments:
        request: Request the username and password for login 
        and store information in session.
    """
    if request.session.get('username') is None:
        username = request.POST.get("username")
        password = request.POST.get("password")
        user = authenticate(username=username, password=password)
        if user is None:
            if request.method == "POST":
                error = "USERNAME AND/OR PASSWORD INCORRECT"
            else:
                error = ""
            return render_to_response("login.html", context={"error": error})
        else:
            request.session["username"] = username
            return HttpResponseRedirect(reverse("index"))
    else:
        return HttpResponseRedirect(reverse("index"))
