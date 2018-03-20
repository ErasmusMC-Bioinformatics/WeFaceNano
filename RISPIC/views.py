import time
import os

from django.shortcuts import render_to_response, HttpResponseRedirect
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth import authenticate
from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render, redirect
from django.urls import reverse
from django.contrib.auth.models import User
from django.http import HttpRequest
from Bio import Entrez, SeqIO
from brigD3 import *
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
        folders = os.walk(network_drive).next()[1]
        files = os.walk(network_drive).next()[2]
        results_folder = "/media/Nanopore/" + request.session.get("username") + "/results"
        results = os.walk(results_folder).next()[1]
    else:
        folders = []
        files = []
        results = []
    return render_to_response("home.html", context={"super": superuser, "user": request.session.get("username"), "folders": folders,
                                                    "files": files, "results": results})


def draw_plasmid(contigfasta, contigname, genbank, refseq, name, length, path):
    """
    todo: Add multiple rings to the same plasmid if contig aligns.
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
        # CDS ring DAR4145, default setting CDS
        ring_gen = AnnotationRing()
        ring_gen.setOptions(color='#2a3e47', name=name, height=20)
        ring_gen.readGenbank(genbank)

        # Misc Feature ring DAR4145
        ring_misc = AnnotationRing()
        ring_misc.setOptions(color='#34596f', name=name, height=20, path=path)
        ring_misc.feature = 'misc_feature'
        ring_misc.extract = {'note': 'Note: '}
        ring_misc.readGenbank(genbank)

        genomes = [contigfasta]
        names = [contigname]
        blaster = Blaster(refseq, genomes=genomes, path=path)  # Blast the contigs against the sequence fasta
        blaster.runBLAST()

        blast_rings = []
        for i in range(len(blaster.results)):
            ring_blast = BlastRing()
            ring_blast.setOptions(color='#88a2af', name=names[i])
            ring_blast.min_length = 100
            ring_blast.readComparison(blaster.results[i])
            blast_rings.append(ring_blast)
        # Combine rings in preferred order
        rings = [ring_gen] + blast_rings + [ring_misc]

        # Initialize ring generator and set options, write as JSON and HTML
        generator = RingGenerator(rings, path, contigname)
        generator.setOptions(circle=length, project=contigname, title=name, title_size='100%', radius=0)
        generator.brigD3()


def assembly(inputtype, inputfolder, barcode_list, resultfolder, gsize, assembler):
    """
    todo: Fix Miniasm barcode and multiple fastq file options.
    Read barcode list and try to assemble all barcodes.
    Return a list of fasta locations and filenames.
    :param inputtype: To see if the input type is FAST5 or FASTQ
    :param inputfolder: To search for the FASTQ files in the input folder.
    :param barcode_list: A list of barcodes created by albacore.
    :param resultfolder: Folder where the reads are stored.
    :param gsize: The genome size for canu assembly.
    :return: A list of fasta paths and filenames.
    """
    fasta_list = []
    file_list = []
    unitigs_barcode = []
    if os.path.exists(resultfolder + "/assembly/"):
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
    else:
        call(["mkdir", resultfolder + "/assembly/"])
        if assembler == "canu":
            for barcode in barcode_list:
                call(["mkdir", resultfolder + "/assembly/" + barcode])
                if barcode != "unclassified":
                    if inputtype == "fast5" and ".fastq" not in barcode:
                        call(["bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                              " -g " + gsize + " -i " + resultfolder + "/workspace/pass/" + barcode + "/"], shell=True)
                    elif inputtype == "fast5" and ".fastq" in barcode:
                        call(["bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                              " -g " + gsize + " -i " + resultfolder + "/workspace/pass/"], shell=True)
                    else:
                        # if ".fastq" not in barcode or "fastq" not in barcode:
                        #     call(["bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                        #           " -g " + gsize + " -i " + inputfolder + "/fastq/" + barcode + "/"], shell=True)
                        # else:
                            call(["bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                                  " -g " + gsize + " -i " + inputfolder + "/fastq"], shell=True)
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
        elif assembler == "miniasm":
            for barcode in barcode_list:
                if barcode != "unclassified":
                    call(["mkdir", resultfolder + "/assembly/" + barcode])
                    call(["cat " + inputfolder + "/" + barcode + "/* > " + inputfolder + "/" + barcode + "/" + barcode + "_cat.fastq"], shell=True)
                    call(["bash ~/RISPIC/static/miniasm.bash -v1 -i " + inputfolder + "/" + barcode + "/" + barcode + "_cat.fastq -o " + resultfolder +
                          "/assembly/" + barcode + "/" + barcode + "_cat.contigs.fasta"], shell=True)
                    call(["rm", "-r", inputfolder + "/" + barcode + "/" + barcode + "_cat.fastq"])
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
    Runs BLAST on assembled FASTA files.
    :param barcodes: List of barcodes that contain unitigs.
    :param file_list: List of assembled FASTA files.
    :param resultfolder: Output folder.
    :param blastdb: BLAST database to use.
    :param task: Task for megablast or blastn selection.
    :return: A list of BLAST output files.
    """
    blast = {} # BLAST dict
    bfiles = [] # BLAST files
    call(["mkdir", resultfolder + "/BLAST"])
    blastfolder = resultfolder + "/BLAST"
    count = 0
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
                        try:
                            bclength = cheader[1][4:]
                        except IndexError:
                            for record in SeqIO.parse(contigfile, "fasta"):
                                bclength = len(record)
                    line = blastfile.readline()
                    if line != "":
                        line = line.split("\t")
                        handle = Entrez.efetch(db="nucleotide", id=line[1], rettype="gb", retmode="gb")
                        handle_txt = Entrez.efetch(db="nucleotide", id=line[1], rettype="gb", retmode="text")
                        record = SeqIO.read(handle, "genbank")
                        bplength = len(record.seq)
                        with open(blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".refseq.fasta", "w") as ref:
                            ref.write(record.format("fasta"))
                        with open(blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".genbank.gb", "w") as gb:
                            gb.write(handle_txt.read())
                        handle.close()
                        handle_txt.close()
                        new_path = blastfolder + "/" + barcodes[count] + "/"
                        contigfasta = blastfolder + "/" + barcodes[count] + "/" + utig
                        genbank = blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".genbank.gb"
                        refseq = blastfolder + "/" + barcodes[count] + "/" + utig[:-6] + ".refseq.fasta"
                        try:
                            for feature in record.features:
                                if feature.qualifiers.get('plasmid', []):
                                    plasmids = feature.qualifiers.get('plasmid')
                            for plasmid in plasmids:
                                p = plasmid
                                print(p)
                            draw_plasmid(contigfasta=contigfasta, contigname=contigname, genbank=genbank, refseq=refseq, name=p, length=bplength, path=new_path)
                            blast[barcodes[count] + "_" + line[0]] = [p, record.description, bplength, bclength]
                        except:
                            draw_plasmid(contigfasta=contigfasta, contigname=contigname, genbank=genbank, refseq=refseq, name=line[1], length=bplength,
                                         path=new_path)
                            try:
                                blast[barcodes[count] + "_" + line[0]] = [line[1], record.description, bplength, bclength]
                            except:
                                pass
                    else:
                        blast[contigname] = ["No Accession Number", "No Name", "0", bclength]
        count += 1
    return blast


@csrf_exempt
def create_results(request):
    """
    todo: Fix the barcode list so it works with Canu and Miniasm
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
    blast = request.POST.get("blast")
    blastdb = request.POST.get("blastdb")
    blasttask = request.POST.get("blasttask")
    resdb = request.POST.get("resdb")
    reslength = request.POST.get("reslength")
    residentity = request.POST.get("residentity")
    resultfolder = "/media/Nanopore/" + request.session.get("username") + "/results/" + outfolder
    qscore = "7"
    # html_dict = {}
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
    fasta_list, file_list, barcodes = assembly(inputtype, inputfolder, barcode_list, resultfolder, gsize, assembler)
    if blast == "-b ":
        blast_dict = run_blast(barcodes, file_list, resultfolder, blastdb, blasttask)
    if res == '-r ':
        res_genes = resfinder(barcodes, file_list, resultfolder, resdb, reslength, residentity)
    else:
        res_genes = []
    return HttpResponseRedirect(reverse("index"))
    #     for b in barcodes:
    #         bfilelist = os.listdir(resultfolder + "/BLAST/" + b)
    #         for bfile in bfilelist:
    #             if ".html" in bfile:
    #                 with open(resultfolder + "/BLAST/" + b + "/" + bfile) as htmlfile:
    #                     html_code = htmlfile.read()
    #                 html_dict[htmlfile.name] = html_code
    # else:
    #     blast_dict = []
    # return render(request, "viewresults.html", context={"fasta": fasta_list, "barcodes": barcodes, "blast": blast_dict, 'blastrun': blast, 'resrun': res,
    #                                                 "user": request.session.get("username"), "res": res_genes, 'html': html_dict})


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
            blastdb = os.listdir("/media/Nanopore/blastdb/")
            walk = os.walk(network_drive).next()[1]
            for folder in walk:
                folders.append(network_drive + folder)
        except StopIteration:
            return HttpResponseRedirect(reverse('index'))
        return render_to_response("pipeline.html", context={"folders": folders, "user": request.session.get('username'),
                                                            "super": superuser, "blastdb": blastdb, "resdb": resdb})


@csrf_exempt
def get_stored_results(request):
    """
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
        if r == selected_result:
            tool_list = []
            tools = os.listdir("/media/Nanopore/" + username + "/results/" + r)
            for t in tools:
                if "." not in t:
                    tool_list.append(t)
                    res_dict[r] = tool_list
                if t == "BLAST":
                    blast_dict = get_stored_blast(username, r)[0]
                    blast_res_dict = get_stored_blast(username, r)[1]
                    html_plasmid = get_stored_blast(username, r)[2]
                if t == "resfinder":
                    resfinder_dict = get_stored_resfinder(username, r)
                if t == "assembly":
                    assembly_report = get_stored_assembly(username, r)[0]
                    barcodes = get_stored_assembly(username, r)[2]
    return render_to_response("results.html", context={"super": superuser, "user": username, "results": results, "tools": tools, "dict": res_dict,
                                                       "blast": blast_dict, "plasmid": html_plasmid, "resfinder": resfinder_dict,
                                                       "blastresults": blast_res_dict, "assemblyreport": assembly_report, "barcodes": barcodes})


def get_stored_blast(username, r):
    """
    Get the stored BLAST result based on the selected run.
    :param username: The user that is logged in
    :param r: Name of the selected run
    :return: blast_dict, blast_res_dict, html_plasmid
    """
    blast_dict = {}
    html_list = []
    html_plasmid = {}
    blast_res_dict = {}
    blast_results = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/BLAST/")
    for br in blast_results:
        bc_blast = os.listdir("/media/Nanopore/" + username + "/results/" + r + "/BLAST/" + br)
        for bc in bc_blast:
            if ".html" in bc:
                with open("/media/Nanopore/" + username + "/results/" + r + "/BLAST/" + br + "/" + bc) as htmlfile:
                    code = htmlfile.read()
                html_list.append(code)
                html_plasmid[bc] = html_list
            html_list = []
            if ".out" in bc:
                with open("/media/Nanopore/" + username + "/results/" + r + "/BLAST/" + br + "/" + bc) as blastfile:
                    first_result = blastfile.readline()
                    first_result = first_result.split('\t')
                    try:
                        contig = first_result[0]
                        pident = first_result[2]
                        alignment_length = first_result[3]
                        handle = Entrez.efetch(db="nucleotide", id=first_result[1], rettype="gb", retmode="gb")
                        record = SeqIO.read(handle, "genbank")
                        blast_name = record.description
                    except IndexError:
                        blast_name = ""
                        contig = ""
                        pident = ""
                        alignment_length = ""
                blast_res_dict[br + "_" + contig] = [blast_name, pident, alignment_length]
    blast_dict[r] = blast_results
    return blast_dict, blast_res_dict, html_plasmid


def get_stored_resfinder(username, r):
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


def get_stored_assembly(username, r):
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
                # with open("/media/Nanopore/" + username + "/results/" + r + "/assembly/" + assemblyfolder + "/" + assembly) as assembly_report_file:
                    # for line in assembly_report_file:
                    #     if ">" in line:
                    #         contigs += 1
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
