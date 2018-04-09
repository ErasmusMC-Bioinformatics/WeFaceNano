import time
import os
import random

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
    """
    Let's superusers create a new user to use the pipeline.
    :param request: Form details
    :return:
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
        superusers = User.objects.filter(is_superuser=True).values_list('username')
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
            message = "Only logged in superusers can create new users. If you are a superuser, please log in."
    return render(request, 'signup.html', context={'form': form, 'user': request.session.get("username"), "super": superuser, "message": message})


@csrf_exempt
def index(request):
    """
    Information on the homepage and check if user is superuser.
    :param request:
    :return:
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
        results_folder = "/media/Nanopore/" + request.session.get("username") + "/results"
        results = os.walk(results_folder).__next__()[1]
    else:
        folders = []
        files = []
        results = []
    return render_to_response("home.html", context={"super": superuser, "user": request.session.get("username"), "folders": folders,
                                                    "files": files, "results": results})


def draw_plasmid(contigfasta, contigname, genbank, refseq, name, length, path):
    """
    Draw a plasmid and aligned contig based on BLAST results.
    :param contigfasta: FASTA file containing the contig
    :param contigname: The name of the contig
    :param genbank: Genbank information
    :param refseq: Reference sequence to BLAST the contig against
    :param name: Name of the plasmid
    :param length: Length of the plasmid
    :param path: Path to store the HTML file with the drawn plasmid alignment
    :return:
    """
    if path != '':
        ring_gen = brigD3.AnnotationRing()
        ring_gen.setOptions(color='#000000', name=name, height=20)
        ring_gen.readGenbank(genbank)
        genomes = contigfasta
        names = contigname
        blaster = brigD3.Blaster(refseq, genomes=genomes, path=path)
        blaster.runBLAST()
        blast_rings = []
        for i in range(len(blaster.results)):
            r = lambda: random.randint(80, 200)
            color = ('#%02X%02X%02X' % (r(), r(), r()))
            ring_blast = brigD3.BlastRing()
            ring_blast.setOptions(color=color, name=names[i])
            ring_blast.min_length = 100
            ring_blast.readComparison(blaster.results[i])
            blast_rings.append(ring_blast)
        rings = [ring_gen] + blast_rings
        generator = brigD3.RingGenerator(rings, path, contigname)
        generator.setOptions(circle=length, project=name, title=name, title_size='100%', radius=200)
        generator.brigD3()


def filter_reads(ref, fastq):
    """
    todo: Remove this function if useless.
    Filters reads using Bowtie2.
    :param ref: Reference to filter reads.
    :param fastq: Reads to be filtered.
    :return: FASTQ with unaligned reads.
    """
    # Takes a long time and gives a strange result with Miniasm assembly.
    unaligned = "filtered_" + fastq
    call(["bowtie2", "-x", ref, "-U", fastq, "-un", unaligned])
    return unaligned


def circularize(assembly, inputfolder, resultfolder, kmer, assembler):
    """
    todo: Test circularize function
    Use circulator to check if a contig is circular or linear.
    :param assembly: Previously build assembly
    :param inputfolder: Folder with the used reads in FASTQ format.
    :param resultfolder: The output folder to store the circulator results.
    :param kmer: Selected kmer size for SPAdes assembly
    :param assembler: Selected assembler for circlator
    :return: Lists with circular and linear contigs.
    """
    barcodes = os.listdir(inputfolder)
    out = resultfolder + "/circularize"
    threads = 8
    for fasta in assembly:
        for bc in barcodes:
            cmd = ("circlator all --assembler " + assembler + " --threads " + str(threads) + " --assemble_spades_k " + kmer +
                   " --merge_min_id 85 --merge_breaklen 1000 " + fasta + " " + inputfolder + "/" + bc + "/trimmed/cat_reads.fastq " + out)
            call(["cat " + inputfolder + "/" + bc + "/trimmed/*.fastq >> " + inputfolder + "/" + bc + "/trimmed/cat_reads.fastq"], shell=True)
            call([cmd], shell=True)
            # call(["circlator", "all", "--assemble_spades_k", kmer, "--threads", str(4), "--merge_min_id", str(85), "--merge_breaklen", str(1000),
            #       fasta, inputfolder + "/" + bc + "/trimmed/cat_reads.fastq", out])
            # call(["rm", inputfolder + "/trimmed/cat_reads.fastq"])


def canu(inputtype, inputfolder, barcode_list, resultfolder, gsize):
    """

    :param inputtype: To see if the input type is FAST5 or FASTQ
    :param inputfolder: To search for the FASTQ files in the input folder.
    :param barcode_list: A list of barcodes.
    :param resultfolder: Folder where the reads are stored.
    :param gsize: The genome size for assembly.
    :return: A list of fasta paths and filenames.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    call(["mkdir", resultfolder + "/assembly/"])
    for barcode in barcode_list:
        call(["mkdir", resultfolder + "/assembly/" + barcode])
        if barcode != "unclassified":
            if inputtype == "fast5":
                call(["bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                      " -g " + gsize + " -i " + resultfolder + "/workspace/pass/" + barcode], shell=True)
            else:
                call(["bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                      " -g " + gsize + " -i " + inputfolder + "/fastq/" + barcode + "/trimmed"], shell=True)
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" + file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        fasta_list.append(resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
        return fasta_list, file_list, unitigs_barcode


def miniasm(inputtype, inputfolder, barcode_list, resultfolder):
    """

    :param inputtype: To see if the input type is FAST5 or FASTQ
    :param inputfolder: To search for the FASTQ files in the input folder.
    :param barcode_list: A list of barcodes.
    :param resultfolder: Folder where the reads are stored.
    :return: A list of fasta paths and filenames.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    call(["mkdir", resultfolder + "/assembly/"])
    for barcode in barcode_list:
        if barcode != "unclassified":
            call(["mkdir", resultfolder + "/assembly/" + barcode])
            if inputtype == "fast5":
                barcode_content = os.listdir(resultfolder + "/workspace/pass/" + barcode)
                if len(barcode_content) > 1:
                    call(["cat " + resultfolder + "/workspace/pass/" + barcode + "/* > " +
                          resultfolder + "/workspace/pass/" + barcode + "/" + barcode + "_cat.fastq"],
                         shell=True)
                    call(["bash ~/RISPIC/static/miniasm.bash -v2 -i " + resultfolder + "/workspace/pass/" + barcode + "/" + barcode + "_cat.fastq -o " +
                          resultfolder + "/assembly/" + barcode + "/" + "_cat.contigs"], shell=True)
                else:
                    call(["bash ~/RISPIC/static/miniasm.bash -v2 -i " + resultfolder + "/workspace/pass/" + barcode + "/" + barcode_content[0] +
                          " -o " + resultfolder + "/assembly/" + barcode + "/" + barcode_content[0] + ".contigs"], shell=True)
            else:
                barcode_content = os.listdir(inputfolder + "/" + barcode)
                if len(barcode_content) > 1:
                    call(["cat " + inputfolder + "/" + barcode + "/trimmed/* > " + inputfolder + "/" + barcode + "/trimmed/" + barcode + "_cat.fastq"],
                         shell=True)
                    call(["bash ~/RISPIC/static/miniasm.bash -v2 -i " + inputfolder + "/" + barcode + "/trimmed/" + barcode + "_cat.fastq -o " +
                          resultfolder + "/assembly/" + barcode + "/" + barcode + "_cat.contigs"], shell=True)
                    call(["rm", "-r", inputfolder + "/" + barcode + "/trimmed/" + barcode + "_cat.fastq"])
                else:
                    call(["bash ~/RISPIC/static/miniasm.bash -v2 -i " + inputfolder + "/" + barcode + "/trimmed/" + barcode_content[0] + " -o " +
                          resultfolder + "/assembly/" + barcode + "/" + barcode_content[0] + ".contigs"], shell=True)
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" + file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        fasta_list.append(resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    return fasta_list, file_list, unitigs_barcode


def skip_assembly(resultfolder, barcode_list):
    """
    Skip the assembly if the assembly folder is already available.
    :param resultfolder: Folder to store the results.
    :param barcode_list: List of barcodes.
    :return: A list of fasta paths and filenames.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    for barcode in barcode_list:
        if barcode != "unclassified":
            files = os.listdir(resultfolder + "/assembly/" + barcode)
            for file in files:
                if "contigs.fasta" in file:
                    with open(resultfolder + "/assembly/" + barcode + "/" + file) as contigfile:
                        head = contigfile.readline()
                    if head == '' or head is None:
                        pass
                    else:
                        fasta_list.append(resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    return fasta_list, file_list, unitigs_barcode


def resfinder(barcodes, file_list, resultfolder, resdb, reslength, residentity):
    """
    Runs resfinder on all FASTA files in the list on the selected database or all databases.
    :param barcodes: A list of available barcodes.
    :param file_list: A list of file names from canu.
    :param resultfolder: The output folder to store resfinder results.
    :param resdb: The database to search resistance genes.
    :param reslength: The minimum length of the sequence.
    :param residentity: The % of identity.
    :return: A list of files with resfinder results.
    """
    call(["mkdir", resultfolder + "/resfinder"])
    db_all = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin", "fusidicacid", "macrolide", "nitroimidazole",
              "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulphonamide", "tetracycline", "trimethoprim", "glycopeptide"]
    res_genes = {}
    count = 0
    for file in file_list:
        call(["mkdir", resultfolder + "/resfinder/" + barcodes[count]])
        call(["bash ~/RISPIC/static/splitfasta.sh " + resultfolder + "/assembly/" + barcodes[count] + "/" + file +
              " " + resultfolder + "/resfinder/" + barcodes[count]], shell=True)
        unitig_bc = os.listdir(resultfolder + "/resfinder/" + barcodes[count])
        for unitig in unitig_bc:
            fasta_unitig = resultfolder + "/resfinder/" + barcodes[count] + "/" + unitig
            if resdb != "all":
                call(["perl ~/RISPIC/static/resfinder.pl -d ~/resfinder/database -i " + fasta_unitig + " -a " + resdb + " -o " +
                      resultfolder + "/resfinder/" + barcodes[count] + " -k " + residentity + " -l " + reslength], shell=True)
            else:
                for db in db_all:
                    call(["perl ~/RISPIC/static/resfinder.pl -d ~/resfinder/database -i " + fasta_unitig + " -a " + db + " -o " +
                          resultfolder + "/resfinder/" + barcodes[count] + " -k " + residentity + " -l " + reslength], shell=True)
            argene = []
            with open(resultfolder + "/resfinder/" + barcodes[count] + "/results_tab.txt") as rf:
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
                    stored_contig = line[3]
                    rcount += 1
        count += 1
    return res_genes


def run_blast(barcodes, file_list, resultfolder, blastdb, task):
    """
    todo: Fix error handling.
    Runs BLAST on assembled FASTA files.
    :param barcodes: List of barcodes that contain unitigs.
    :param file_list: List of assembled FASTA files.
    :param resultfolder: Output folder.
    :param blastdb: BLAST database to use.
    :param task: Task for megablast or blastn selection.
    :return: A list of BLAST output files.
    """
    blast = {}
    bfiles = []
    call(["mkdir", resultfolder + "/BLAST"])
    blastfolder = resultfolder + "/BLAST"
    count = 0
    dict_contigfasta = {}
    dict_draw = {}
    plasmidcount = 0
    Entrez.email = "some_email@somedomain.com"
    for fasta in file_list:
        call(["mkdir", resultfolder + "/BLAST/" + barcodes[count]])
        call(["bash ~/RISPIC/static/splitfasta.sh " + resultfolder + "/assembly/" + barcodes[count] + "/" + fasta + " " + blastfolder + "/" + barcodes[count]],
             shell=True)
        unitig_bc = os.listdir(blastfolder + "/" + barcodes[count])
        if "unclassified" in os.listdir(blastfolder):
            call(["~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb + " -query " + blastfolder + "/unclassified/" + utig +
                  " -out " + blastfolder + "/unclassified/" + utig + ".out -max_target_seqs 1 -outfmt 6 -task " + task], shell=True)
            for unc in os.listdir(blastfolder + "/unclassified/"):
                if ".out" in unc:
                    bfiles.append(blastfolder + "/unclassified/" + unc)
        else:
            for utig in unitig_bc:
                call(["~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb + " -query " + blastfolder + "/" + barcodes[count] + "/" + utig +
                      " -out " + blastfolder + "/" + barcodes[count] + "/" + utig + ".out -max_target_seqs 1 -outfmt 6 -task " + task], shell=True)
                with open(blastfolder + "/" + barcodes[count] + "/" + utig + ".out") as blastfile:
                    with open(blastfolder + "/" + barcodes[count] + "/" + utig) as contigfile:
                        cheader = contigfile.readline()
                        cheader = cheader.split(" ")
                        contigname = barcodes[count] + "_" + cheader[0][1:]
                        for record in SeqIO.parse(contigfile.name, "fasta"):
                            bclength = len(record.seq)
                        line = blastfile.readline()
                    if line != "":
                        line = line.split("\t")
                        handle = Entrez.efetch(db="nucleotide", id=line[1].split("_")[-1], rettype="gb", retmode="gb")
                        handle_txt = Entrez.efetch(db="nucleotide", id=line[1].split("_")[-1], rettype="gb", retmode="text")
                        record = SeqIO.read(handle, "genbank")
                        gblength = len(record.seq)
                        if int(gblength) < int(bclength):
                            bplength = gblength
                        else:
                            bplength = bclength
                        with open(blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".refseq.fasta", "w") as ref:
                            ref.write(record.format("fasta"))
                        with open(blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".genbank.gb", "w") as gb:
                            gb.write(handle_txt.read())
                        new_path = blastfolder + "/" + barcodes[count] + "/"
                        genbank = blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".genbank.gb"
                        refseq = blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".refseq.fasta"
                        try:
                            for feature in record.features:
                                if feature.qualifiers.get('plasmid', []):
                                    plasmids = feature.qualifiers.get('plasmid')
                            for plasmid in plasmids:
                                p = plasmid
                                dict_contigfasta[p + "_" + str(plasmidcount)] = utig
                                dict_draw[p + "_" + str(plasmidcount)] = [contigname, genbank, refseq, bplength, new_path]
                                plasmidcount += 1
                            blast[barcodes[count] + "_" + line[0]] = [p, record.description, bplength, bclength]
                        except:
                            dict_contigfasta[line[1].split("_")[-1]] = utig
                            dict_draw[line[1].split("_")[-1]] = [contigname, genbank, refseq, bplength, new_path]
                            try:
                                blast[barcodes[count] + "_" + line[0]] = [line[1], record.description, bplength, bclength]
                            except:
                                pass
                    else:
                        blast[contigname] = ["No Accession Number", "No Name", "0", bclength]
            contigkeys = dict_contigfasta.keys()
            usedkeys = []
            for key in contigkeys:
                utgkeys = []
                for dplasmid, dutigs in dict_contigfasta.items():
                    if key[:-2] in dplasmid:
                        utgkeys.append(dutigs)
                    else:
                        pass
                if key[:-2] not in usedkeys:
                    try:
                        contigs = []
                        for utgfasta in utgkeys:
                            with open(blastfolder + "/" + barcodes[count] + "/" + utgfasta) as unitig_fasta:
                                for line in unitig_fasta:
                                    if ">" in line:
                                        contigs.append(line)
                        if "_" in key:
                            name = key.split("_")[0]
                        else:
                            name = key
                        draw_plasmid(contigfasta=utgkeys, contigname=contigs, genbank=dict_draw[key][1],
                                     refseq=dict_draw[key][2], name=name, length=dict_draw[key][3], path=dict_draw[key][4])
                    except IndexError:
                        pass
                else:
                    pass
                usedkeys.append(key[:-2])
        count += 1
    return blast


@csrf_exempt
def create_results(request):
    """
    Run the pipeline with entered data and show the results in an html page.
    :param request: User inputs given in the pipeline.html page.
    :return: Results page in html.
    """
    inputfolder = request.POST.get("inputfolder")
    inputtype = request.POST.get("inputtype")
    configuration = request.POST.get("aconfig")
    outfolder = request.POST.get("albacoreout")
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
    resultfolder = "/media/Nanopore/" + request.session.get("username") + "/results/" + outfolder
    qscore = "7"
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
        call(["bash ~/RISPIC/static/npanalysis.sh " + res + mlst + species + "-b -c " + configuration +
              " -i" + inputfolder + "/fast5 -o " + resultfolder + " -f fastq"], shell=True)
        call(["Rscript minion_qc/MinionQC.R -i " + resultfolder + "/sequencing_summary.txt -o " + resultfolder + "/qc" + " -q " + qscore], shell=True)
    elif request.POST.get("barcoding") == '0' and inputtype == "fast5":
        call(["mkdir", resultfolder + "/workspace/pass"])
        call(["mkdir", resultfolder + "/qc"])
        call(["bash ~/RISPIC/static/npanalysis.sh " + res + mlst + species + "-c " + configuration +
              " -i" + inputfolder + "/fast5 -o " + resultfolder + " -f fastq"], shell=True)
        call(["Rscript minion_qc/MinionQC.R -i " + resultfolder + "/sequencing_summary.txt -o " + resultfolder + "/qc" + " -q " + qscore], shell=True)
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
        else:
            barcode_list = os.listdir(inputfolder)
            trim_reads(inputfolder, barcode_list)
    if os.path.exists(resultfolder + "/assembly/"):
        skip_assembly(barcode_list)
    else:
        if assembler == "canu":
            fasta_list, file_list, barcodes = canu(inputtype, inputfolder, barcode_list, resultfolder, gsize)
        elif assembler == "miniasm":
            fasta_list, file_list, barcodes = miniasm(inputtype, inputfolder, barcode_list, resultfolder)
        else:
            return HttpResponseRedirect(reverse("index"))
    if circulator == "yes":
        circularize(fasta_list, inputfolder, resultfolder, kmer, c_assembler)
    if blast == "-b ":
        run_blast(barcodes, file_list, resultfolder, blastdb, blasttask)
    if res == '-r ':
        resfinder(barcodes, file_list, resultfolder, resdb, reslength, residentity)
    for bc in barcode_list:
        call(["rm", "-rf", inputfolder + "/" + bc + "/trimmed"])
    return HttpResponseRedirect(reverse("index"))


def trim_reads(inputfolder, barcode_list):
    """
    Trim the adapters from the reads using Porechop.
    :param inputfolder: The reads inputfolder.
    :param barcode_list: The list of barcodes to find the FASTQ files.
    :return:
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
    """
    Shows the readme page and checks if user is a superuser.
    :param request:
    :return:
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
    return render_to_response("readme.html", context={"user": request.session.get("username"), "super": superuser})


@csrf_exempt
def pipeline_start(request):
    """
    Shows the pipeline page and all the available options.
    Checks if user is superuser.
    :param request:
    :return:
    """
    resdb = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin", "fusidicacid", "macrolide", "nitroimidazole",
             "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulphonamide", "tetracycline", "trimethoprim", "glycopeptide"]
    if request.session.get('username') is None:
        return HttpResponseRedirect(reverse("index"))
    else:
        superusers = User.objects.filter(is_superuser=True).values_list('username')
        su_list = []
        for users in superusers:
            for u in users:
                su_list.append(u)
        if request.session.get("username") in su_list:
            superuser = True
        else:
            superuser = False
        network_drive = "/media/Nanopore/" + request.session.get("username") + "/"
        folders = []
        try:
            blastdb = []
            blastdbfolder = os.listdir("/media/Nanopore/blastdb/")
            for db in blastdbfolder:
                if os.path.isdir("/media/Nanopore/blastdb/" + db):
                    blastdb.append(db)
            walk = os.walk(network_drive).__next__()[1]
            for folder in walk:
                folders.append(network_drive + folder)
        except StopIteration:
            return HttpResponseRedirect(reverse('index'))
        return render_to_response("pipeline.html", context={"folders": folders, "user": request.session.get('username'),
                                                            "super": superuser, "blastdb": blastdb, "resdb": resdb})


@csrf_exempt
def get_stored_results(request):
    """
    todo: Add circularize results.
    Gets stored results from the network drive after the user selects a previously run analysis from the homepage.
    Shows the resfinder results and the plasmid contig alignment HTML pages.
    :param request:
    :return: Web page with the stored results
    """
    superusers = User.objects.filter(is_superuser=True).values_list('username')
    su_list = []
    blast_dict = {}
    html_plasmid = {}
    blast_res_dict = {}
    resfinder_dict = {}
    assembly_report = []
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
                tools = os.listdir("/media/Nanopore/" + username + "/results/" + r)
                for t in tools:
                    if "." not in t:
                        tool_list.append(t)
                        res_dict[r] = tool_list
                    if t == "BLAST":
                        blast_dict = get_stored_blast_results(username, r)[0]
                        blast_res_dict = get_stored_blast_results(username, r)[1]
                        html_plasmid = get_stored_blast_results(username, r)[2]
                        topology = get_stored_blast_results(username, r)[3]
                    if t == "resfinder":
                        resfinder_dict = get_stored_resfinder_results(username, r)
                    if t == "assembly":
                        assembly_report = get_stored_assembly_results(username, r)[0]
                        barcodes = get_stored_assembly_results(username, r)[2]
                    if t == "circularize":
                        contig_topology = get_stored_circularize_results(username, r)
    return render_to_response("results.html", context={"super": superuser, "user": username, "results": results, "tools": tools, "dict": res_dict,
                                                       "blast": blast_dict, "plasmid": html_plasmid, "resfinder": resfinder_dict,
                                                       "blastresults": blast_res_dict, "assemblyreport": assembly_report, "barcodes": barcodes,
                                                       "topology": topology, "contigtopology": contig_topology})


def get_stored_circularize_results(username, r):
    """
    Get the stored circularization results created by circulator.
    :param username: The user that is logged in
    :param r: Name of the selected run
    :return: Dictionary of the contig topology showing if a contig is circular or linear
    """
    circ_folder = "/media/Nanopore/" + username + "/results/" + r + "/circularize/"
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
    """
    todo: Fix alignment length information from BLAST output.
    todo: Change output headers in local BLAST output file.
    Get the stored BLAST result based on the selected run.
    :param username: The user that is logged in
    :param r: Name of the selected run
    :return: blast_dict, blast_res_dict, html_plasmid
    """
    blast_dict = {}
    html_plasmid = {}
    blast_res_dict = {}
    topology = {}
    blast_results = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/BLAST/")
    for br in blast_results:
        bc_blast = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/BLAST/" + br)
        for bc in bc_blast:
            if ".html" in bc:
                with open("/media/Nanopore/" + username + "/results/" + r + "/BLAST/" + br + "/" + bc) as htmlfile:
                    code = htmlfile.read()
                html_plasmid[bc.split(".")[0]] = code
            if ".out" in bc:
                with open("/media/Nanopore/" + username + "/results/" + r + "/BLAST/" + br + "/" + bc) as blastfile:
                    first_result = blastfile.readline()
                    first_result = first_result.split('\t')
                    try:
                        contig = first_result[0]
                        pident = first_result[2]
                        alignment_length = first_result[3]
                        handle = Entrez.efetch(db="nucleotide", id=first_result[1].split("_")[-1], rettype="gb", retmode="gb")
                        record = SeqIO.read(handle, "genbank")
                        topology[record.description] = record.annotations["topology"]
                        blast_name = record.description
                    except IndexError:
                        blast_name = ""
                        contig = ""
                        pident = ""
                        alignment_length = ""
                if contig:
                    blast_res_dict[br + "_" + contig] = [blast_name, pident, alignment_length]
                else:
                    pass
    blast_dict[r] = blast_results
    return blast_dict, blast_res_dict, html_plasmid, topology


def get_stored_resfinder_results(username, r):
    """
    Get stored resfinder results based on the selected run
    :param username: The user that is logged in
    :param r: Name of the selected run
    :return: Resfinder results (Contig name, antibiotic resistance genes and percentag identity)
    """
    resfinder_dict = {}
    contiglist = []
    arglist = []
    resfinder_results = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/resfinder/")
    for resfolder in resfinder_results:
        resbarcode = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/resfinder/" + resfolder)
        if "results_tab.txt" in resbarcode:
            if os.stat("/media/Nanopore/" + username + "/results/" + r + "/resfinder/" + resfolder + "/results_tab.txt").st_size > 0:
                with open("/media/Nanopore/" + username + "/results/" + r + "/resfinder/" + resfolder + "/results_tab.txt") as resfinderfile:
                    for line in resfinderfile:
                        line = line.split("\t")
                        if line[3] not in contiglist:
                            if contiglist:
                                resfinder_dict[resfolder + "_" + contiglist[-1]] = arglist
                            arglist = []
                            contiglist.append(line[3])
                        arglist.append(str(line[0] + " - (" + line[1] + "%) - " + line[5]))
                    try:
                        resfinder_dict[resfolder + "_" + contiglist[-1]] = arglist
                    except IndexError:
                        pass
            else:
                pass
    return resfinder_dict


def get_stored_assembly_results(username, r):
    """
    Get stored assembly results based on the selected run
    :param username: The user that is logged in
    :param r: Name of the selected run
    :return: Assembly report with the amount of contigs and the length of the assembly
    """
    assembly_report = {}
    assembly_results = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/assembly/")
    assembly_bc = []
    for assemblyfolder in assembly_results:
        assembly_barcode = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/assembly/" + assemblyfolder)
        for assembly in assembly_barcode:
            if ".contigs.fasta" in assembly:
                contigcount = 0
                if os.stat("/media/Nanopore/" + username + "/results/" + r + "/assembly/" + assemblyfolder + "/" + assembly).st_size > 0:
                    assembly_bc.append(assemblyfolder)
                    for record in SeqIO.parse("/media/Nanopore/" + username + "/results/" + r + "/assembly/" + assemblyfolder + "/" + assembly, "fasta"):
                        contigcount += 1
                        contiglength = len(record.seq)
                        assembly_report[assemblyfolder + "_" + record.id] = contiglength
    return assembly_report, assembly_results, sorted(assembly_bc)


def logout(request):
    """
    Flush the session to logout the current user.
    :param request:
    :return:
    """
    request.session.flush()
    return HttpResponseRedirect(reverse("index"))


@csrf_exempt
def login(request):
    """
    Authenticate user. Check if user is in the mysql database and password is correct.
    Store username in a session.
    :param request:
    :return:
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
