import os
import commands
import time

from django.shortcuts import render_to_response, HttpResponseRedirect
from django.views.decorators.csrf import csrf_exempt
from django.contrib.auth import authenticate
from django.contrib.auth.forms import UserCreationForm
from django.shortcuts import render, redirect
from django.core.urlresolvers import reverse
from django.contrib.auth.models import User
from Bio import Entrez
from Bio import SeqIO
from brigD3 import *


@csrf_exempt
def signup(request):
    if request.method == 'POST':
        form = UserCreationForm(request.POST)
        if form.is_valid():
            form.save()
            username = form.cleaned_data.get('username')
            raw_password = form.cleaned_data.get('password1')
            commands.getoutput("mkdir /media/Nanopore/" + username)
            return redirect('index')
    else:
        form = UserCreationForm()
    return render(request, 'signup.html', context={'form': form, 'user': request.session.get("username")})


@csrf_exempt
def index(request):
    superusers = User.objects.filter(is_superuser=True)
    if str(request.session.get("username")) in str(superusers.get()):
        superuser = True
    else:
        superuser = False
    if request.session.get("username"):
        network_drive = "/media/Nanopore/" + request.session.get("username")
        folders = os.walk(network_drive).next()[1]
        files = os.walk(network_drive).next()[2]
    else:
        folders = []
        files = []
    return render_to_response("home.html", context={"super": superuser, "user": request.session.get("username"), "folders": folders, "files": files})


def draw_plasmid(contigfasta, contigname, genbank, refseq, name, length, path):
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
        generator.setOptions(circle=length, project=contigname, title=name, title_size='100%', radius=200)
        generator.brigD3()


def canu(inputtype, inputfolder, barcode_list, resultfolder, gsize):
    """
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
        for barcode in barcode_list:
            if barcode != "unclassified":
                if inputtype == "fast5":
                    os.system("bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                              " -g " + gsize + " -i " + resultfolder + "/workspace/pass/" + barcode + "/")
                else:
                    os.system("bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                              " -g " + gsize + " -i " + inputfolder + "/" + barcode + "/")
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
    commands.getoutput("mkdir " + resultfolder + "/resfinder")
    db_all = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin", "fusidicacid", "macrolide", "nitroimidazole",
              "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulphonamide", "tetracycline", "trimethoprim", "glycopeptide"]
    res_genes = {}
    count = 0
    for file in file_list:
        commands.getoutput("mkdir " + resultfolder + "/resfinder/" + barcodes[count])
        os.system("bash ~/RISPIC/static/splitfasta.sh " + resultfolder + "/assembly/" + barcodes[count] + "/" + file +
                  " " + resultfolder + "/resfinder/" + barcodes[count])
        unitig_bc = os.listdir(resultfolder + "/resfinder/" + barcodes[count])
        for unitig in unitig_bc:
            fasta_unitig = resultfolder + "/resfinder/" + barcodes[count] + "/" + unitig
            if resdb != "all":
                commands.getoutput("perl ~/RISPIC/static/resfinder.pl -d ~/resfinder/database -i " + fasta_unitig + " -a " + resdb + " -o " +
                                   resultfolder + "/resfinder/" + barcodes[count] + " -k " + residentity + " -l " + reslength)
            else:
                for db in db_all:
                    commands.getoutput("perl ~/RISPIC/static/resfinder.pl -d ~/resfinder/database -i " + fasta_unitig + " -a " + db + " -o " +
                                       resultfolder + "/resfinder/" + barcodes[count] + " -k " + residentity + " -l " + reslength)
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
            # fastafiles = os.listdir(resultfolder + "/resfinder/" + barcodes[count])
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
    commands.getoutput("mkdir " + resultfolder + "/BLAST")
    blastfolder = resultfolder + "/BLAST"
    count = 0
    Entrez.email = "some_email@somedomain.com"
    for fasta in file_list:
        commands.getoutput("mkdir " + resultfolder + "/BLAST/" + barcodes[count])
        os.system("bash ~/RISPIC/static/splitfasta.sh " + resultfolder + "/assembly/" + barcodes[count] + "/" + fasta + " " + blastfolder + "/" + barcodes[count])
        unitig_bc = os.listdir(blastfolder + "/" + barcodes[count])
        if "unclassified" in os.listdir(blastfolder):
            os.system("~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb + " -query " + blastfolder + "/unclassified/" + utig +
                      " -out " + blastfolder + "/unclassified/" + utig + ".out -max_target_seqs 1 -outfmt 6 -task " + task)
            for unc in os.listdir(blastfolder + "/unclassified/"):
                if ".out" in unc:
                    bfiles.append(blastfolder + "/unclassified/" + unc)
        else:
            for utig in unitig_bc:
                os.system("~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb + " -query " + blastfolder + "/" + barcodes[count] + "/" + utig +
                          " -out " + blastfolder + "/" + barcodes[count] + "/" + utig + ".out -max_target_seqs 1 -outfmt 6 -task " + task)
                with open(blastfolder + "/" + barcodes[count] + "/" + utig + ".out") as blastfile:
                    with open(blastfolder + "/" + barcodes[count] + "/" + utig) as contigfile:
                        cheader = contigfile.readline()
                        cheader = cheader.split(" ")
                        contigname = barcodes[count] + "_" + cheader[0][1:]
                        bclength = cheader[1][4:]
                    line = blastfile.readline()
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
                        draw_plasmid(contigfasta=contigfasta, contigname=contigname, genbank=genbank, refseq=refseq, name=p, length=bplength, path=new_path)
                        blast[barcodes[count] + "_" + line[0]] = [p, record.description, bplength, bclength]
                    except:
                        draw_plasmid(contigfasta=contigfasta, contigname=contigname, genbank=genbank, refseq=refseq, name=line[1], length=bplength,
                                     path=new_path)
                        try:
                            blast[barcodes[count] + "_" + line[0]] = [line[1], record.description, bplength, bclength]
                        except:
                            pass

            # for files in os.listdir(blastfolder + "/" + barcodes[count] + "/"):
            #     if ".out" in files:
            #         bfiles.append(blastfolder + "/" + barcodes[count] + "/" + files)
            # blast_files = os.listdir(blastfolder + "/" + barcodes[count])
            # for bl in blast_files:
            #     if ".fasta.out" in bl:
            #         with open(blastfolder + "/" + barcodes[count] + "/" + bl) as blfile:
            #             line = blfile.readline()
            #             line = line.split("\t")
            #             handle = Entrez.efetch(db="nucleotide", id=line[1], rettype="gb", retmode="gb")
            #             record = SeqIO.read(handle, "genbank")
            #             bplength = len(record.seq)
            #             with open(blastfolder + "/" + barcodes[count] + "/" + bl[:-4]) as ffile:
            #                 header = ffile.readline()
            #                 header.split("\t")
            #                 bclength = header[1][4:]
            #                 contigname = header[0][1:]
            #             try:
            #                 for feature in record.features:
            #                     if feature.qualifiers.get('plasmid', []):
            #                         plasmids = feature.qualifiers.get('plasmid')
            #                 for plasmid in plasmids:
            #                     p = plasmid
            #                 handle.close()
                            # new_path = blastfolder + "/" + barcodes[count] + "/"
                            # contigfasta = blastfolder + "/" + barcodes[count] + "/" + bl
                            # genbank = blastfolder + "/" + barcodes[count] + "/" + line[1] + ".genbank.gb"
                            # refseq = blastfolder + "/" + barcodes[count] + "/" + line[1] + ".refseq.fasta"
                            # draw_plasmid(contigfasta=contigfasta, contigname=contigname, genbank=genbank, refseq=refseq, name=p, length=bplength, path=new_path)
                            # blast[barcodes[count] + "_" + line[0]] = [p, record.description, bplength, bclength]
                        # except:
                        #     try:
                        #         draw_plasmid(contigfasta=contigfasta, contigname=contigname, genbank=genbank, refseq=refseq, name=line[1], length=bplength,
                        #                      path=new_path)
                        #         blast[barcodes[count] + "_" + line[0]] = [line[1], record.description, bplength, bclength]
                        #     except:
                        #         pass
        count += 1
    return blast


@csrf_exempt
def show_results(request):
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
    blast = request.POST.get("blast")
    blastdb = request.POST.get("blastdb")
    blasttask = request.POST.get("blasttask")
    resdb = request.POST.get("resdb")
    reslength = request.POST.get("reslength")
    residentity = request.POST.get("residentity")
    resultfolder = "/media/Nanopore/" + request.session.get("username") + "/results/" + outfolder
    qscore = "7"
    if res is None:
        res = ''
    if mlst is None:
        mlst = ''
    if species is None:
        species = ''
    request.session['stored_results'] = request.POST
    commands.getoutput("mkdir " + resultfolder)
    if request.POST.get("barcoding") == '1' and inputtype == "fast5":
        commands.getoutput("mkdir " + resultfolder + "/workspace/pass")
        commands.getoutput("mkdir " + resultfolder + "/qc")
        os.system("bash ~/RISPIC/static/npanalysis.sh " + res + mlst + species + "-b -c " +
                  configuration + " -i" + inputfolder + " -o " + resultfolder + " -f fastq")
        commands.getoutput("Rscript minion_qc/MinionQC.R -i " + resultfolder +
                           "/sequencing_summary.txt -o " + resultfolder + "/qc" + " -q " + qscore)
    elif request.POST.get("barcoding") == '0' and inputtype == "fast5":
        commands.getoutput("mkdir " + resultfolder + "/workspace/pass")
        commands.getoutput("mkdir " + resultfolder + "/qc")
        os.system("bash ~/RISPIC/static/npanalysis.sh " + res + mlst + species + "-c " +
                  configuration + " -i" + inputfolder + " -o " + resultfolder + " -f fastq")
        commands.getoutput("Rscript minion_qc/MinionQC.R -i " + resultfolder +
                           "/sequencing_summary.txt -o " + resultfolder + "/qc" + " -q " + qscore)
    if inputtype == "fast5":
        time.sleep(15)
        while True:
            if os.path.isfile(resultfolder + "/finished.dat"):
                break
        barcode_list = os.listdir(resultfolder + "/workspace/pass")
    else:
        if os.path.isdir(resultfolder + "/assembly"):
            barcode_list = os.listdir(resultfolder + "/assembly")
        else:
            barcode_list = os.listdir(inputfolder)
    fasta_list, file_list, barcodes = canu(inputtype, inputfolder, barcode_list, resultfolder, gsize)
    if blast == "-b ":
        html_dict = {}
        blast_dict = run_blast(barcodes, file_list, resultfolder, blastdb, blasttask)
        for b in barcodes:
            bfilelist = os.listdir(resultfolder + "/BLAST/" + b)
            for bfile in bfilelist:
                if ".html" in bfile:
                    with open(resultfolder + "/BLAST/" + b + "/" + bfile) as htmlfile:
                        html_code = htmlfile.read()
                    html_dict[htmlfile.name] = html_code
    else:
        blast_dict = []
    if res == '-r ':
        res_genes = resfinder(barcodes, file_list, resultfolder, resdb, reslength, residentity)
    else:
        res_genes = []
    return render(request, "results.html", context={"fasta": fasta_list, "barcodes": barcodes, "blast": blast_dict, 'blastrun': blast, 'resrun': res,
                                                    "user": request.session.get("username"), "res": res_genes, 'html': html_dict})


def readme(request):
    superusers = User.objects.filter(is_superuser=True)
    if str(request.session.get("username")) in str(superusers.get()):
        superuser = True
    else:
        superuser = False
    return render_to_response("readme.html", context={"user": request.session.get("username"), "super": superuser})


@csrf_exempt
def pipeline_start(request):
    resdb = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin", "fusidicacid", "macrolide", "nitroimidazole",
             "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulphonamide", "tetracycline", "trimethoprim", "glycopeptide"]
    if request.session.get('username') is None:
        return HttpResponseRedirect(reverse("index"))
    else:
        superusers = User.objects.filter(is_superuser=True)
        if str(request.session.get("username")) in str(superusers.get()):
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


def logout(request):
    request.session.flush()
    return HttpResponseRedirect(reverse("index"))


@csrf_exempt
def login(request):
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
            request.session["password"] = password
            return HttpResponseRedirect(reverse("index"))
    else:
        return HttpResponseRedirect(reverse("index"))
