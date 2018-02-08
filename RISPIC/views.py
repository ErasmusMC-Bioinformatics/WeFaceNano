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
    return render(request, 'signup.html', context={'form': form})


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
    # commands.getoutput("mkdir " + resultfolder + "/assembly")
    if os.path.exists(resultfolder + "/assembly/"):
        for barcode in barcode_list:
            if barcode != "unclassified":
                files = os.listdir(resultfolder + "/assembly/" + barcode)
                for file in files:
                    if "contigs.fasta" in file:
                        fasta_list.append(resultfolder + "/assembly/" + barcode + "/" + file)
                        file_list.append(file)
                        unitigs_barcode.append(barcode)
    else:
        for barcode in barcode_list:
            if barcode != "unclassified":
                if inputtype == "fastq":
                    os.system("bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                              " -g " + gsize + " -i " + resultfolder + "/workspace/pass/" + barcode + "/")
                else:
                    os.system("bash ~/RISPIC/static/canu.sh -p " + barcode + " -d " + resultfolder + "/assembly/" + barcode +
                              " -g " + gsize + " -i " + inputfolder + "/" + barcode + "/")
                files = os.listdir(resultfolder + "/assembly/" + barcode)
                for file in files:
                    if "contigs.fasta" in file:
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
    db_all = ["aminoglycoside", "beta-lactam", "colistin", "fosfomycin", "fusidicacid", "macrolide", "nitroimidazole",
              "oxazolidinone", "phenicol", "quinolone", "rifampicin", "sulphonamide", "tetracycline", "trimethoprim", "glycopeptide"]
    res_genes = {}
    count = 0
    # port = 4000
    for file in file_list:
        # time.sleep(10)
        # os.system("fuser -k " + str(port) + "/tcp")
        resfile_list = []
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
                    commands.getoutput("perl ~/resfinder/resfinder.pl -d ~/resfinder/database -i " + fasta_unitig + " -a " + db + " -o " +
                                       resultfolder + "/resfinder/" + barcodes[count] + " -k " + residentity + " -l " + reslength)

            # os.system("jsa.util.streamServer --port=" + str(port) +
            #           " | bwa mem -t12 -x ont2d -K10000 /home/myfair/japsa/ResGene/resFinder/DB.fasta - 2> /dev/null | jsa.np.rtResistGenes --output=" +
            #           fasta_unitig + ".out --bamfile=- --resDB=/home/myfair/japsa/ResGene/resFinder --time=60 --tmp=" +
            #           resultfolder + "/tmp/resTemp 2> /dev/null | sleep 5 | jsa.util.streamClient -i " + fasta_unitig + " --server=0.0.0.0:" + str(port))
            # resfile_list.append(fasta_unitig + ".out")
            # port += 1
            # time.sleep(60)
            # os.system("fuser -k " + str(port) + "/tcp")
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
                        argene.append(line[0])
                    else:
                        argene = []
                        argene.append(line[0])
                    res_genes[barcodes[count] + "_" + contig] = argene
                    stored_contig = line[3]
                    rcount += 1
            # with open(fasta_unitig + ".out") as rf:
            #     for line in rf:
            #         if "##" not in line:
            #             line = line.split("\t")
            #             argene = line[5]
            #             with open(fasta_unitig) as ufile:
            #                 header = ufile.readline()
            #                 header = header.split(" ")
            #             res_genes[barcodes[count] + "_" + header[0][1:]] = argene
        count += 1
    return res_genes


def run_blast(barcodes, file_list, resultfolder, blastdb, task):
    """
    Runs BLAST on assembled FASTA files.
    :param barcodes: List of barcodes that contain unitigs.
    :param file_list: List of assembled FASTA files.
    :param resultfolder: Output folder.
    :return: A list of BLAST output files.
    """
    blast = {} # BLAST dict
    bfiles = [] # BLAST files
    commands.getoutput("mkdir " + resultfolder + "/BLAST")
    blastfolder = resultfolder + "/BLAST"
    count = 0
    for fasta in file_list:
        commands.getoutput("mkdir " + resultfolder + "/BLAST/" + barcodes[count])
        os.system("bash ~/RISPIC/static/splitfasta.sh " + resultfolder + "/assembly/" + barcodes[count] + "/" + fasta + " " + blastfolder + "/" + barcodes[count])
        unitig_bc = os.listdir(blastfolder + "/" + barcodes[count])
        if "unclassified" in os.listdir(blastfolder):
            time.sleep(10)
            os.system("~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb + " -query " + blastfolder + "/unclassified/" + utig +
                      " -out " + blastfolder + "/unclassified/" + utig + ".out -max_target_seqs 1 -outfmt 6 -task " + task)
            for unc in os.listdir(blastfolder + "/unclassified/"):
                time.sleep(10)
                if ".out" in files:
                    bfiles.append(blastfolder + "/unclassified/" + files)
        else:
            for utig in unitig_bc:
                time.sleep(10)
                os.system("~/ncbi-blast-2.7.1+/bin/blastn -db " + blastdb + " -query " + blastfolder + "/" + barcodes[count] + "/" + utig +
                          " -out " + blastfolder + "/" + barcodes[count] + "/" + utig + ".out -max_target_seqs 1 -outfmt 6 -task " + task)
            for files in os.listdir(blastfolder + "/" + barcodes[count] + "/"):
                time.sleep(10)
                if ".out" in files:
                    bfiles.append(blastfolder + "/" + barcodes[count] + "/" + files)
            blast_files = os.listdir(blastfolder + "/" + barcodes[count])
            for bl in blast_files:
                time.sleep(10)
                if ".fasta.out" in bl:
                    with open(blastfolder + "/" + barcodes[count] + "/" + bl) as blfile:
                        line = blfile.readline()
                        line = line.split("\t")
                        try:
                            Entrez.email = "some_email@somedomain.com"
                            handle = Entrez.efetch(db="nucleotide", id=line[1], rettype="gb", retmode="gb")
                            record = SeqIO.read(handle, "genbank")
                            for feature in record.features:
                                if feature.qualifiers.get('plasmid', []):
                                    plasmids = feature.qualifiers.get('plasmid')
                            for plasmid in plasmids:
                                p = plasmid
                            handle.close()
                            blast[barcodes[count] + "_" + line[0]] = [p, record.description]
                        except:
                            try:
                                print record.description
                                blast[barcodes[count] + "_" + line[0]] = [line[1], record.description]
                            except:
                                pass
                else:
                    commands.getoutput("rm " + blastfolder + "/" + barcodes[count] + "/" + bl)
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
    gsa = request.POST.get("gsa")
    qscore = request.POST.get("qscore")
    gsize = request.POST.get("gsize")
    blast = request.POST.get("blast")
    blastdb = request.POST.get("blastdb")
    blasttask = request.POST.get("blasttask")
    resdb = request.POST.get("resdb")
    reslength = request.POST.get("reslength")
    residentity = request.POST.get("residentity")
    resultfolder = "/media/Nanopore/" + request.session.get("username") + "/results/" + outfolder

    if res is None:
        res = ''
    if mlst is None:
        mlst = ''
    if species is None:
        species = ''
    request.session['stored_results'] = request.POST
    commands.getoutput("mkdir " + resultfolder)
    commands.getoutput("mkdir " + resultfolder + "/workspace/pass")
    commands.getoutput("mkdir " + resultfolder + "/resfinder")
    if request.POST.get("barcoding") == '1' and inputtype == "fast5":
        commands.getoutput("mkdir " + resultfolder + "/qc")
        os.system("bash ~/RISPIC/static/npanalysis.sh " + res + mlst + species + "-b -c " +
                  configuration + " -i" + inputfolder + " -o " + resultfolder + " -f fastq")
        if qscore is None or qscore == "":
            qscore = "7"
        commands.getoutput("Rscript minion_qc/MinionQC.R -i " + resultfolder +
                           "/sequencing_summary.txt -o " + resultfolder + "/qc" + " -q " + qscore)
    elif request.POST.get("barcoding") == '0' and inputtype == "fast5":
        commands.getoutput("mkdir " + resultfolder + "/qc")
        os.system("bash ~/RISPIC/static/npanalysis.sh " + res + mlst + species + "-c " +
                  configuration + " -i" + inputfolder + " -o " + resultfolder + " -f fastq")
        if qscore is None or qscore == "":
            qscore = "7"
        commands.getoutput("Rscript minion_qc/MinionQC.R -i " + resultfolder +
                           "/sequencing_summary.txt -o " + resultfolder + "/qc" + " -q " + qscore)
    if inputtype != "fast5":
        time.sleep(15)
        while True:
            if os.path.isfile(resultfolder + "/finished.dat"):
                break
        barcode_list = os.listdir(resultfolder + "/workspace/pass")
        print barcode_list
    else:
        barcode_list = os.listdir(inputfolder)
    fasta_list, file_list, barcodes = canu(inputtype, inputfolder, barcode_list, resultfolder, gsize)
    if blast == "-b ":
        blast_dict = run_blast(barcodes, file_list, resultfolder, blastdb, blasttask)
    else:
        blast_dict = {}
    if res == '-r ':
        res_genes = resfinder(barcodes, file_list, resultfolder, resdb, reslength, residentity)
    else:
        res_genes = []
    return render(request, "results.html", context={"fasta": fasta_list, "barcodes": barcodes, "blast": blast_dict,
                                                    "user": request.session.get("username"), "res": res_genes})


def readme(request):
    superusers = User.objects.filter(is_superuser=True)
    if str(request.session.get("username")) in str(superusers.get()):
        superuser = True
    else:
        superuser = False
    return render_to_response("readme.html", context={"user": request.session.get("username"), "super": superuser})


@csrf_exempt
def pipeline_start(request):
    resdb = ["aminoglycoside", "beta - lactam", "colistin", "fosfomycin", "fusidicacid", "macrolide", "nitroimidazole",
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
            # error = "There is no network folder called " + request.session.get("username") + ". Please contact an administrator."
            return HttpResponseRedirect(reverse('index'))
            # return render_to_response("home.html", context={"error": error})
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
