import time
import os
import random
import re

from urllib.error import HTTPError
from django.shortcuts import render_to_response, HttpResponseRedirect
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth import authenticate
from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render, redirect
from django.urls import reverse
from django.contrib.auth.models import User
from Bio import Entrez, SeqIO
from . import brigD3
from subprocess import call


@csrf_exempt
def signup(request):
    """Let's superusers create a new user to use the pipeline.

    Arguments:
        request -- Request information used to sign up.
    """
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            if os.path.exists("/media/Nanopore/" + username):
                pass
            else:
                call(["mkdir", "/media/Nanopore/" + username])
            if os.path.exists("/media/Nanopore/" + username + "/results"):
                pass
            else:
                call(["mkdir", "/media/Nanopore/" + username + "/results"])
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
    """Information on the homepage and check if user is superuser.

    Arguments:
        request -- Request information when user is logged in. 
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
        network_drive = "/media/Nanopore/" + request.session.get("username")
        folders = os.walk(network_drive).__next__()[1]
        files = os.walk(network_drive).__next__()[2]
        results_folder = ("/media/Nanopore/" +
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
        "results": results})


def draw_plasmid(contigfasta, contigname, genbank, refseq,
                 name, length, path, res_loc):
    """Draw a plasmid and aligned contig based on BLAST results.

    Arguments:
        contigfasta {list} -- A list of contigsfiles.
        contigname {list} -- A list of contig names based on the files.
        genbank {list} -- A list of genbank files to retrieve gene info.
        refseq {str} -- The reference filename to align the contigs with.
        name {list} -- A list of plasmid names.
        length {int} -- The length of the plasmid.
        path {str} -- The path to store te HTML output.
        res_loc {dict} -- The names and locations of the 
        antibiotic resistance genes.
    """
    if path != '':
        ring_gen = brigD3.AnnotationRing()
        ring_gen.setOptions(color='#000000',
                            name=name,
                            height=20,
                            res_loc=res_loc)
        # if genbank != "":
        ring_gen.readGenbank(genbank, length)
        genomes = contigfasta
        names = contigname
        blaster = brigD3.Blaster(refseq,
                                 genomes=genomes,
                                 path=path)
        blaster.runBLAST()
        blast_rings = []
        for i in range(len(blaster.results)):
            def r(): return random.randint(80, 200)
            color = ('#%02X%02X%02X' % (r(), r(), r()))
            ring_blast = brigD3.BlastRing()
            ring_blast.setOptions(color=color,
                                  name=names[i],
                                  res_loc=res_loc)
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


def circularize(assembly, inputfolder, resultfolder, kmer, assembler):
    """Use circulator to check if a contig is circular or linear.

    Arguments:
        assembly {str} -- Filename of the first assembly.
        inputfolder {str} -- The folder containing the reads.
        resultfolder {str} -- Name of the outputfolder to store results.
        kmer {str} -- Kmer-size when using SPAdes assembly.
        assembler {str} -- Selected assembler for use with circlator.
    """
    barcodes = os.listdir(inputfolder)
    out = resultfolder + "/circularize"
    threads = 8
    for fasta in assembly:
        for bc in barcodes:
            call(["mkdir", out + "/" + bc])
            if assembler == "canu":
                if ".fasta" in os.listdir(inputfolder + "/" + bc +
                                          "/trimmed")[0]:
                    cmd = (
                        "circlator all --assembler canu "
                        "--data_type nanopore-raw --bwa_opts '-x ont2d' "
                        "--threads " + str(threads) +
                        " --merge_min_id 85 --merge_breaklen 1000 " + fasta +
                        " " + inputfolder + "/" + bc +
                        "/trimmed/cat_reads.fasta " + out + "/" + bc
                    )
                elif ".fastq" in os.listdir(inputfolder + "/" + bc +
                                            "/trimmed")[0]:
                    cmd = (
                        "circlator all --assembler canu "
                        "--data_type nanopore-raw --bwa_opts '-x ont2d' "
                        "--threads " + str(threads) +
                        " --merge_min_id 85 --merge_breaklen 1000 " + fasta +
                        " " + inputfolder + "/" + bc +
                        "/trimmed/cat_reads.fastq " + out + "/" + bc
                    )
            else:
                if ".fasta" in os.listdir(inputfolder + "/" + bc +
                                          "/trimmed")[0]:
                    cmd = (
                        "circlator all --threads " + str(threads) +
                        " --bwa_opts '-x ont2d' --assemble_spades_k " + kmer +
                        " --merge_min_id 85 --merge_breaklen 1000 " + fasta +
                        " " + inputfolder + "/" + bc +
                        "/trimmed/cat_reads.fasta " + out + "/" + bc
                    )
                elif ".fastq" in os.listdir(inputfolder + "/" + bc +
                                            "/trimmed")[0]:
                    cmd = (
                        "circlator all --threads " + str(threads) +
                        "--bwa_opts '-x ont2d' --assemble_spades_k " + kmer +
                        " --merge_min_id 85 --merge_breaklen 1000 " +
                        fasta + " " + inputfolder + "/" + bc +
                        "/trimmed/cat_reads.fastq " + out + "/" + bc
                    )
            if ".fasta" in os.listdir(inputfolder + "/" + bc + "/trimmed")[0]:
                call([
                    "cat " + inputfolder + "/" + bc + "/trimmed/*.fasta >> " +
                    inputfolder + "/" + bc + "/trimmed/cat_reads.fasta"
                ], shell=True)
            elif ".fastq" in os.listdir(inputfolder + "/" + bc +
                                        "/trimmed")[0]:
                call([
                    "cat " + inputfolder + "/" + bc + "/trimmed/*.fastq >> " +
                    inputfolder + "/" + bc + "/trimmed/cat_reads.fastq"
                ], shell=True)
            call([cmd], shell=True)


def canu(inputtype, inputfolder, barcode_list, resultfolder, gsize):
    """Use Canu for assembly and return the assembly as a FASTA file.

    Arguments:
        inputtype {str} -- To see if the inputfiles are FAST5 or FASTQ.
        inputfolder {str} -- Name of the folder with the inputfiles.
        barcode_list {list} -- A list of barcodes.
        resultfolder {str} -- Folder where the results will be stored.
        gsize {str} -- The genome size for the Canu assembler.

    Returns:
        list -- A list of FASTA paths.
        list -- A list of FASTA filenames
        list -- A list of barcodes.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    call(["mkdir", resultfolder + "/assembly/"])
    for barcode in barcode_list:
        call(["mkdir", resultfolder + "/assembly/" + barcode])
        if barcode != "unclassified":
            if inputtype == "fast5":
                call([
                    "bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " +
                    resultfolder + "/assembly/" + barcode + " -g " + gsize +
                    " -i " + resultfolder + "/workspace/pass/" + barcode
                ], shell=True)
            else:
                call([
                    "bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " +
                    resultfolder + "/assembly/" + barcode + " -g " + gsize +
                    " -i " + inputfolder + "/" + barcode + "/trimmed"
                ], shell=True)
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" +
                              file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        fasta_list.append(
                            resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
        return fasta_list, file_list, unitigs_barcode


def miniasm(inputtype, inputfolder, barcode_list, resultfolder):
    """Use Miniasm for assembly and return the assembly as a FASTA file.

    Arguments:
        inputtype {str} -- To see if the inputfiles are FAST5 or FASTQ.
        inputfolder {str} -- Name of the folder with the inputfiles.
        barcode_list {list} -- A list of barcodes.
        resultfolder {str} -- Folder where the results will be stored.

    Returns:
        list -- A list of FASTA paths.
        list -- A list of FASTA filenames
        list -- A list of barcodes.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    call(["mkdir", resultfolder + "/assembly/"])
    for barcode in barcode_list:
        if barcode != "unclassified":
            call(["mkdir", resultfolder + "/assembly/" + barcode])
            if inputtype == "fast5":
                barcode_content = os.listdir(
                    resultfolder + "/workspace/pass/" + barcode)
                if len(barcode_content) > 1:
                    call([
                        "cat " + resultfolder + "/workspace/pass/" + barcode +
                        "/* > " + resultfolder + "/workspace/pass/" + barcode +
                        "/" + barcode + "_cat.fastq"
                    ], shell=True)
                    call([
                        "bash ~/RISPIC/static/miniasm.bash -v2 -i " +
                        resultfolder + "/workspace/pass/" + barcode + "/" +
                        barcode + "_cat.fastq -o " + resultfolder +
                        "/assembly/" + barcode + "/" + "_cat.contigs"
                    ], shell=True)
                else:
                    call([
                        "bash ~/RISPIC/static/miniasm.bash -v2 -i " +
                        resultfolder + "/workspace/pass/" + barcode + "/" +
                        barcode_content[0] + " -o " + resultfolder +
                        "/assembly/" + barcode + "/" + barcode_content[0] +
                        ".contigs"
                    ], shell=True)
            else:
                barcode_content = os.listdir(inputfolder + "/" + barcode)
                if len(barcode_content) > 1:
                    call([
                        "cat " + inputfolder + "/" + barcode +
                        "/trimmed/* > " + inputfolder + "/" + barcode +
                        "/trimmed/" + barcode + "_cat.fasta"
                    ], shell=True)
                    call([
                        "bash ~/RISPIC/static/miniasm.bash -v2 -i " +
                        inputfolder + "/" + barcode + "/trimmed/" + barcode +
                        "_cat.fasta -o " + resultfolder + "/assembly/" +
                        barcode + "/" + barcode + "_cat.contigs"
                    ], shell=True)
                    call([
                        "rm", "-r", inputfolder + "/" + barcode +
                        "/trimmed/" + barcode + "_cat.fasta"
                    ])
                else:
                    call([
                        "bash ~/RISPIC/static/miniasm.bash -v2 -i " +
                        inputfolder + "/" + barcode + "/trimmed/" +
                        barcode_content[0] + " -o " + resultfolder +
                        "/assembly/" + barcode + "/" + barcode_content[0] +
                        ".contigs"
                    ], shell=True)
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" +
                              file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        fasta_list.append(
                            resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    return fasta_list, file_list, unitigs_barcode


def skip_assembly(resultfolder, barcode_list):
    """Skip the assembly if the assembly folder is already available.

    Arguments:
        resultfolder {str} -- Folder where the results will be stored.
        barcode_list {list} -- A list of barcodes.

    Returns:
        list -- A list of FASTA paths.
        list -- A list of FASTA filenames
        list -- A list of barcodes.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    for barcode in barcode_list:
        if barcode != "unclassified":
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" +
                              file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        fasta_list.append(
                            resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    return fasta_list, file_list, unitigs_barcode


def resfinder(barcodes, file_list, resultfolder,
              resdb, reslength, residentity):
    """Runs resfinder on all FASTA files in the list on the 
    selected database or all databases.

    Arguments:
        barcodes {list} -- A list of barcodes.
        file_list {list} -- A list of files that contains the assembly.
        resultfolder {str} -- The ouputfolder to store the results.
        resdb {str} -- The database to search antibiotic resistance genes.
        reslength {str} -- The minimum length cutoff when using resfinder.
        residentity {str} -- The % of identity cutoff when using resfinder.

    Returns:
        dict -- A dictionary with the contigs and the found genes.
        dict -- A dictionary of the gene locations within the contigs.
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
            "bash ~/RISPIC/static/splitfasta.sh " + resultfolder +
            "/assembly/" + barcodes[count] + "/" + file + " " + resultfolder +
            "/resfinder/" + barcodes[count]
        ], shell=True)
        unitig_bc = os.listdir(resultfolder + "/resfinder/" + barcodes[count])
        for unitig in unitig_bc:
            fasta_unitig = (resultfolder + "/resfinder/" +
                            barcodes[count] + "/" + unitig)
            if resdb != "all":
                call([
                    "perl ~/RISPIC/static/resfinder.pl "
                    "-d ~/resfinder/database -i " + fasta_unitig +
                    " -a " + resdb + " -o " + resultfolder + "/resfinder/" +
                    barcodes[count] + " -k " + residentity + " -l " + reslength
                ], shell=True)
            else:
                for db in db_all:
                    call([
                        "perl ~/RISPIC/static/resfinder.pl -d "
                        "~/resfinder/database -i " + fasta_unitig +
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
                    res_loc[contig + "_" + line[0] + "_" +
                            str(rcount)] = [loc[0], loc[1]]
                    stored_contig = line[3]
                    rcount += 1
        count += 1
    return res_genes, res_loc


def run_blast(barcodes, file_list, resultfolder, blastdb, task, res_loc):
    """Runs BLAST on assembled FASTA files.

    TODO: Use fixed start FASTA if contig is circularized.
    FIXME: Function is too long. Split or clean function.

    Arguments:
        barcodes {list} -- A list of barcodes that contains contigs.
        file_list {list} -- A list of filenames to use with BLAST.
        resultfolder {str} -- The folder where the 
        BLAST results will be stored.
        blastdb {str} -- The selected BLAST database.
        task {str} -- The selected task in BLAST (blastn or megablast)
        res_loc {dict} -- Locations of resistance genes to draw plasmid.

    Returns:
        dict -- BLAST output information (plasmid, description, length).
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
    # plasmidlist = []
    Entrez.email = "some_email@somedomain.com"
    for fasta in file_list:
        unitig_bc = []
        call(["mkdir", resultfolder + "/BLAST/" + barcodes[bcount]])
        call([
            "bash ~/RISPIC/static/splitfasta.sh " + resultfolder +
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
                    "~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb +
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
                    "~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb +
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
                            if int(gblength) < int(bclength):
                                bplength = gblength
                            else:
                                bplength = bclength
                        except HTTPError:
                            bplength = bclength
                            genbank = ""
                        new_path = blastfolder + "/" + barcodes[bcount] + "/"
                        if blastdb == "/media/Nanopore/plasmidb/plasmidb":
                            for ref_folder in os.listdir(
                                    "/media/Nanopore/Rick"):
                                if "RB_REF" in ref_folder:
                                    for ref_file in os.listdir(
                                            "/media/Nanopore/Rick/" +
                                            ref_folder):
                                        if line[1].strip("/") in ref_file:
                                            with open("/media/Nanopore/Rick/" +
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
                        # path=dict_draw[key][4],
                        path=new_path_bc,
                        res_loc=res_loc)
                except (IndexError, KeyError):
                    pass
            else:
                pass
            # usedkeys.append(key[:-2])
            usedkeys.append(key.split("_")[0])
        barcodes.pop(0)
        dict_contigfasta = {}
        dict_draw = {}
    return blast


@csrf_exempt
def create_results(request):
    """Run the pipeline with entered data and show the results in an html page.

    Arguments:
        request -- Request information to get information to run the pipeline.
    """
    inputfolder = request.POST.get("inputfolder")
    inputtype = request.POST.get("inputtype")
    configuration = request.POST.get("aconfig")
    outfolder = request.POST.get("albacoreout").replace(" ", "_")
    res = request.POST.get("res")
    mlst = request.POST.get("mlst")
    species = request.POST.get("species")
    gsize = request.POST.get("gsize")
    assembler = request.POST.get("assembler")
    c_assembler = request.POST.get("c_assembler")
    blast = request.POST.get("blast")
    blastdb = request.POST.get("blastdb")
    blasttask = request.POST.get("blasttask")
    resdb = request.POST.get("resdb")
    reslength = request.POST.get("reslength")
    residentity = request.POST.get("residentity")
    kmer = request.POST.get("kmer")
    circulator = request.POST.get("circularize")
    resultfolder = ("/media/Nanopore/" + request.session.get("username") +
                    "/results/" + outfolder)
    # qscore = "7"
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
            "bash ~/RISPIC/static/albacore.sh " + res + mlst + species +
            "-b -c " + configuration + " -i" + inputfolder + "/fast5 -o " +
            resultfolder + " -f fastq"
        ], shell=True)
        # call([
        #     "Rscript minion_qc/MinionQC.R -i " + resultfolder +
        #     "/sequencing_summary.txt -o " + resultfolder + "/qc -q " + qscore
        # ], shell=True)
    elif request.POST.get("barcoding") == '0' and inputtype == "fast5":
        call(["mkdir", resultfolder + "/workspace/pass"])
        call(["mkdir", resultfolder + "/qc"])
        call([
            "bash ~/RISPIC/static/albacore.sh " + res + mlst + species +
            "-c " + configuration + " -i" + inputfolder + "/fast5 -o " +
            resultfolder + " -f fastq"
        ], shell=True)
        # call([
        #     "Rscript minion_qc/MinionQC.R -i " + resultfolder +
        #     "/sequencing_summary.txt -o " + resultfolder + "/qc -q " + qscore
        # ], shell=True)
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
        fasta_list, file_list, barcodes = skip_assembly(
            resultfolder, barcode_list)
    else:
        if assembler == "canu":
            fasta_list, file_list, barcodes = canu(
                inputtype, inputfolder, barcode_list, resultfolder, gsize)
        elif assembler == "miniasm":
            fasta_list, file_list, barcodes = miniasm(
                inputtype, inputfolder, barcode_list, resultfolder)
        else:
            return HttpResponseRedirect(reverse("index"))
    
    # Create QC pages
    for barcode in barcode_list:
        files = os.listdir(inputfolder + "/" + barcode + "/trimmed")
        for file in files:
            file = inputfolder + "/" + barcode + "/trimmed/" + file
            qc_report(file.split(".")[1], barcode, file,
                      resultfolder + "/qc/" + barcode)

    if circulator == "yes":
        circularize(fasta_list, inputfolder, resultfolder, kmer, c_assembler)
    if res == '-r ':
        dummyres_genes, res_loc = resfinder(
            barcodes, file_list, resultfolder, resdb, reslength, residentity)
    else:
        res_loc = {}
    if blast == "-b ":
        run_blast(
            barcodes, file_list, resultfolder, blastdb, blasttask, res_loc)
    # for bc in barcode_list:
    #     call(["rm", "-rf", inputfolder + "/" + bc + "/trimmed"])
    return HttpResponseRedirect(reverse("index"))


def qc_report(extention, barcode, filename, outfolder):
    """Runs the NanoPot QC tool on the FASTA or FASTQ file. Generates
    an HTML page with QC stats and read distribution plots. When using a
    FASTQ file quality score information will be shown.
    
    Arguments:
        extention {str} -- File extention (fasta or fastq).
        barcode {str} -- Barcode name.
        filename {str} -- Name of the file to run the QC.
        outfolder {str} -- Where to store the QC data.
    """
    cmd = ("NanoPlot --" + extention + " " + filename +
           " --N50 --plot kde hex dot -p " + barcode +
           " -o " + outfolder)
    call([cmd], shell=True)


def trim_reads(inputfolder, barcode_list):
    """Trim the adapters from the reads using Porechop.

    Arguments:
        inputfolder {str} -- The folder containing the reads.
        barcode_list {list} -- A list of barcodes to find the FASTQ files.
    """
    for barcodes in barcode_list:
        barcode_folder = os.listdir(inputfolder + "/" + barcodes)
        if "trimmed" not in barcode_folder:
            for file in barcode_folder:
                file_path = inputfolder + "/" + barcodes + "/" + file
                call(["mkdir", inputfolder + "/" + barcodes + "/trimmed/"])
                trim_path = inputfolder + "/" + barcodes + "/trimmed/" + file
                call(["porechop", "-i", file_path, "-o", trim_path])
        else:
            pass


def readme(request):
    """Shows the readme page and checks if user is a superuser.

    Arguments:
        request -- A request to check if the logged in user is a superuser.
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
        request -- Request information to start the plasmid pipeline.
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
        network_drive = ("/media/Nanopore/" +
                         request.session.get("username") + "/")
        folders = []
        try:
            blastdb = []
            blastdbfolder = os.listdir("/media/Nanopore/blastdb/")
            for db in blastdbfolder:
                if os.path.isdir("/media/Nanopore/blastdb/" + db):
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
            "resdb": resdb})


@csrf_exempt
def get_stored_results(request):
    """Gets stored results from the network drive after the user selects a 
    previously run analysis from the homepage.

    Arguments:
        request -- Request information to retrieve stored results.
    """
    superusers = User.objects.filter(is_superuser=True).values_list('username')
    su_list = []
    blast_dict = {}
    html_plasmid = {}
    blast_res_dict = {}
    resfinder_dict = {}
    assembly_report = []
    topology = {}
    contig_topology = {}
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
    results = os.listdir("/media/Nanopore/" + username + "/results/")
    res_dict = {}
    selected_result = request.POST.get("resultname")
    for r in results:
        if selected_result != "none":
            if r == selected_result:
                tool_list = []
                tools = os.listdir("/media/Nanopore/" +
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
                        if t == "assembly":
                            assembly_report = get_stored_assembly_results(
                                username, r)[0]
                            barcodes = get_stored_assembly_results(
                                username, r)[2]
                        if t == "circularize":
                            contig_topology = get_stored_circularize_results(
                                username, r)
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
        "blastresults": blast_res_dict,
        "assemblyreport": assembly_report,
        "barcodes": barcodes,
        "topology": topology,
        "qc": qc_html,
        "contigtopology": contig_topology})


def get_stored_qc(username, r):
    qc_dict = {}
    qc_bc = ("/media/Nanopore/" + username +
             "/results/" + r + "/qc/")
    for bc in os.listdir(qc_bc):
        folder = ("/media/Nanopore/" + username +
                  "/results/" + r + "/qc/" + bc)
        for qc in os.listdir(folder):
            if ".html" in qc:
                with open(folder + "/" + qc) as qcpage:
                    qc_html = qcpage.read()
                qc_dict[bc] = qc_html
    return qc_dict


def get_stored_circularize_results(username, r):
    """Get the stored circularization results created by circulator.

    Arguments:
        username {str} -- The user that is logged in.
        r {str} -- Name of the selected run

    Returns:
        dict -- Dictionary of the contig topology showing 
        if a contig is circular or linear
    """
    circ_folder = ("/media/Nanopore/" + username +
                   "/results/" + r + "/circularize/")
    circularize_data = os.listdir(circ_folder)
    contig_topology = {}
    for cd in circularize_data:
        if "04.merge.circularise_details.log" in cd:
            with open(circ_folder + cd) as circ_log:
                for line in circ_log:
                    line = line.split("\t")
                    try:
                        if "Circularized: no" in line[2]:
                            contig_topology[line[1]] = "linear"
                        elif "Circularized: yes" in line[2]:
                            contig_topology[line[1]] = "circular"
                    except IndexError:
                        pass
    return contig_topology


def get_stored_blast_results(username, r):
    """Get the stored BLAST result based on the selected run.

    Arguments:
        username {str} -- The user that is logged in.
        r {str} -- Name of the selected run

    Returns:
    blast_dict, blast_res_dict, html_plasmid
        dict -- A dictionary with BLAST files.
        dict -- A dictionary with BLAST result information.
        dict -- A dictionary with the HTML code to view the plasmid.
    """
    blast_dict = {}
    html_plasmid = {}
    blast_res_dict = {}
    topology = {}
    coveragelist = []
    total_length = []
    blast_results = os.listdir(
        "/media/Nanopore/" + username + "/results/" + r + "/BLAST/")
    for br in blast_results:
        bc_blast = os.listdir("/media/Nanopore/" +
                              username + "/results/" + r + "/BLAST/" + br)
        for bc in bc_blast:
            if ".html" in bc:
                with open("/media/Nanopore/" + username + "/results/" + r +
                          "/BLAST/" + br + "/" + bc) as htmlfile:
                    code = htmlfile.read()
                html_plasmid[br + "_" + bc.split(".")[0]] = code
            if ".out" in bc:
                count = 0
                with open("/media/Nanopore/" + username + "/results/" + r +
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
                        # pident = first_result[2]
                        pident = sum(coveragelist)/len(coveragelist)
                        # alignment_length = first_result[3]
                        alignment_length = sum(total_length)
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
                            # return HttpResponseRedirect('index')
                        # record = SeqIO.read(handle, "genbank")
                        # topology[record.description] = record.annotations[
                            # "topology"]
                        # blast_name = record.description
                    except IndexError:
                        blast_name = ""
                        contig = ""
                        pident = ""
                        alignment_length = ""
                if contig:
                    blast_res_dict[br + "_" + contig] = [
                        blast_name, pident, alignment_length]
                else:
                    pass
    blast_dict[r] = blast_results
    return blast_dict, blast_res_dict, html_plasmid, topology


def get_stored_resfinder_results(username, r):
    """Get stored resfinder results based on the selected run

    Arguments:
        username {str} -- The user that is logged in.
        r {str} -- Name of the selected run.

    Returns:
        dict -- A dictionary with resfinder results 
        (Contig name, antibiotic resistance genes and percentag identity).
    """
    resfinder_dict = {}
    arglist = []
    resfinder_results = os.listdir(
        "/media/Nanopore/" + username + "/results/" + r + "/resfinder/")
    for resfolder in resfinder_results:
        contiglist = []
        resbarcode = os.listdir(
            "/media/Nanopore/" + username + "/results/" + r +
            "/resfinder/" + resfolder)
        if "results_tab.txt" in resbarcode:
            if os.stat(
                "/media/Nanopore/" + username + "/results/" + r +
                    "/resfinder/" + resfolder + "/results_tab.txt"
            ).st_size > 0:
                with open("/media/Nanopore/" + username + "/results/" + r +
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


def get_stored_assembly_results(username, r):
    """Get stored assembly results based on the selected run

    Arguments:
        username {str} -- The user that is logged in.
        r {str} -- Name of the selected run.

    Returns:
        dict -- A dictionary with the assembly and contig lengths.
        list -- A list of files containing the assembly results.
        list -- A list of barcodes.
    """
    assembly_report = {}
    assembly_results = os.listdir(
        "/media/Nanopore/" + username + "/results/" + r + "/assembly/")
    assembly_bc = []
    for assemblyfolder in assembly_results:
        assembly_barcode = os.listdir(
            "/media/Nanopore/" + username + "/results/" + r +
            "/assembly/" + assemblyfolder)
        for assembly in assembly_barcode:
            if ".contigs.fasta" in assembly:
                contigcount = 0
                if os.stat(
                    "/media/Nanopore/" + username + "/results/" + r +
                    "/assembly/" + assemblyfolder + "/" + assembly
                ).st_size > 0:
                    assembly_bc.append(assemblyfolder)
                    for record in SeqIO.parse(
                        "/media/Nanopore/" + username + "/results/" + r +
                            "/assembly/" + assemblyfolder + "/" + assembly,
                            "fasta"):
                        contigcount += 1
                        contiglength = len(record.seq)
                        assembly_report[assemblyfolder +
                                        "_" + record.id] = contiglength
    return assembly_report, assembly_results, sorted(assembly_bc)


def delete(request):
    """Checks if user is superuser and gets the results created by the user.
    Shows a dropdown menu on the delete page to select a result the user wants 
    to delete.
    
    Arguments:
        request -- Request information to visit the delete page.
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
    results = os.listdir("/media/Nanopore/" + username + "/results/")
    return render_to_response("delete.html", context={
        'super': superuser,
        'user': username,
        'results': results})


@csrf_exempt
def remove_result(request):
    """Remove the selected result from the users results folder on the
    network drive. This cannot be undone.
    
    Arguments:
        request -- Get the username and resultname to remove the 
        selected result from the network drive.
    """
    result = request.POST.get("resultname")
    username = request.POST.get("username")
    cmd = "/media/Nanopore/" + username + "/results/" + result
    call(["rm", "-rf", cmd])
    return HttpResponseRedirect(reverse("index"))


def logout(request):
    """Flush the session to logout the current user.

    Arguments:
        request -- Request the session of the current user to 
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
        request -- Request the username and password for login 
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
