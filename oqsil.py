#!/usr/bin/env python3

import os
import pandas as pd
import subprocess
import matplotlib.pyplot as plt
from glob import glob

#got  my AP from NCBI, should now work
os.environ['NCBI_API_KEY'] = 'cb6b860872dfccb5f04396d6f056d4b36c08'

###Greeting the user and explaining the names of the program 
while True:
    print("Welcome to OQSIL-Online protein Quality aSsesIng tooL")
    #and taking values from the user
    oq_prot = input("Please enter the name of the protein you would like to analyze: ")
    takson = input("Please enter the Taxon ID(txid): ")
    ##Winning some time
    print(f">Connecting to NCBI\n\n>Searching for {oq_prot} fastas for txid{takson}")
    # sending request to ncbi for the with protein name and takson
    zapros = f'{oq_prot}[Protein Name] AND txid{takson}[organism]'
    qidiruv_komandasi = f'esearch -db protein -query "{zapros}" | efetch -format fasta'
    ##calculating the number of found fastas
    kolichestvo = os.popen(f'{qidiruv_komandasi} | grep ">" | wc -l').read().strip()
    # clarifying whether user they wants to download all fastas
    if kolichestvo == "0":
        print("Error: No results found for the specified organizm or txid. \nEither there is no fastas presented for this specific protein and/or txid. \nPlease check https://www.ncbi.nlm.nih.gov/, specifying protein {oq_prot} and txid{txid}\n \nIf these data are presented in NCBI, then Please check your input and try again. \nLets start from the beginning!!!")
        continue
    else:
        print(f">Generating results \nOverall:{kolichestvo} {oq_prot} fastas found. \nPLEASE KEEP IN MIND. If the number of proteins is too high, further generated results, might be not that much beautiful and informative as You expect")
    
    #Downloading all what user wants if he/she wants
    skachay_vsyo = input(f'Do you want to download all of the {oq_prot} fastas for for txid {takson} ? (y/n): ')

    if skachay_vsyo.lower() == 'y':
        print("\n>Downloading all the fastas")
        # downloading all fastas
        yuklab_olamiz = f'{qidiruv_komandasi} > {oq_prot}.fa'
        os.system(yuklab_olamiz)
        ##Saving the  output
        print(f'\n\n>All {kolichestvo} fastas for txid {takson} have been saved to {oq_prot}.fa')           
        oq_protlar = f"{oq_prot}.fa"                              
        with open(oq_protlar, "r") as f:
            lines = f.readlines()  
        # Creating a dictionary to hold the protein names and their sequences
        proteins = {}
        for line in lines:
            if line.startswith(">"):
                protein_name = line.strip()[1:]
                proteins[protein_name] = ""
            else:
                proteins[protein_name] += line.strip()
        
	    #MIN/MAX/MEAN
                
                
        # Creating a pandas DataFrame from the dictionary
        oq_protcha = pd.DataFrame.from_dict(proteins, orient="index", columns=["Sequence"])
        # Counting the number of nucleotides for each protein
        oq_protcha["Nucleotide Count"] = oq_protcha["Sequence"].apply(len)
        smallest_one = oq_protcha["Nucleotide Count"].min()
        longest_one = oq_protcha["Nucleotide Count"].max()
        oq_prot_mean_length = oq_protcha["Nucleotide Count"].mean()
        print(f"The min. length of these fastas is {smallest_one} bp and the longest sequence comprises of {longest_one} bp. \n Mean length of {oq_prot} sequences is {oq_prot_mean_length} bp")
                           
        with open(f'{oq_prot}.fa', 'r') as oq_protlar:
            fasta_tarkibi = oq_protlar.read()     

        while True:
            keyingi_etap = input('\n\n>>>What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
            if keyingi_etap.lower() == 'a':
                continue
            elif keyingi_etap.lower() == 'b':
            # Again specifying the input and output files for clustalo
            ###So this code will be taking the downloaded file and creating two files
            ##one multiple sequence aligned and one distance matrix. 
                fasta_fayli = f"{oq_prot}.fa"
                clust_fayli = fasta_fayli[:-3] + ".msf"
                d_file = f"{oq_prot}.dist"
                # Generating Clustalo matrix with the similarity scores for sequence
                clustalo_path = "/usr/bin/clustalo"
                os.environ["PATH"] += os.pathsep + clustalo_path
                clustalo_run = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force --full --distmat-out={d_file}"
                try:
                    subprocess.run(clustalo_run, shell=True, check=True, capture_output=True, text=True)
                    print("\n\n>>>Clustal Omega alignment is finished. Generating the results. Please wait")
                    # Taking out  similarity scores from the generated results
                    similarities_between = []
                    with open(clust_fayli, "r") as f:
                        for line in f:
                            if line.startswith('#'):
                                score = float(line.split(':')[-1].strip())
                                similarities_between.append(f"{score:.2%}")
                    print(f"Similarity scores between fastas: {similarities_between}")
                    # Distance matrix going to screen and might be saved as file at a later stage
                    with open(d_file, "r") as f:
                        print(f.read())
                except subprocess.CalledProcessError as e:
                    print(f"Unfortunately, there was an error running Clustal Omega: {e}")
                # Visualising the aligned sequences to the screen
                with open(clust_fayli, "r") as f:
                    print(f.read())
                #Clarifying whether the user wants the result to be saved?
                result_desicision = input("Would you like to save the output? (y/n)")
                #if the user wants to save then downloading results into two different files
                #the aligned pairs and the distance matrix
                if result_desicision.lower() == "y":
                    # At this point we are saving the output files with the protein {oq_prot} name
                    ##clustal one
                    yangi_clustalo = f"{oq_prot}.msf"
                    os.rename(clust_fayli, yangi_clustalo)
                    print(f"Alignment result saved as {yangi_clustalo}.")
                    ###and the dist matrix
                    new_dist_msf = f"{oq_prot}.matrix"
                    os.rename(d_file, new_dist_msf)
                    print(f"Similarity scores saved as {new_dist_msf}.")
                    # Prompt the user to plot the level of conservation between protein sequences
                    plot_con = input("Do you want to plot the level of conservation between protein sequences? (y/n): ")

                    if plot_con.lower() == "y":
                        # Run the command to plot conservation level
                        conservation_plotting = f"plotcon -sformat msf {yangi_clustalo} -winsize 16 -graph pdf"
                        
                        os.system(conservation_plotting)
                        # view the generated PDF file on screen
                        view_pdf = "gs plotcon.pdf"
                        os.system(view_pdf)
                        # Display the plot
                        plt.show()
                        print("Conservation plot saved as plotcon.pdf")
                       
                        ### PATMATMOTIFS 
                       
                        
                        # asking the  user if they want to compare protein motifs
                        motif_comparison = input("Do you want to identify protein motifs in your set? (y/n): ")
                        if motif_comparison.lower() == "y":
                            # Setting the path to the patmatmotifs command
                            
                            which_patmatmotifs = "/usr/bin/patmatmotifs"
                            #The output will look like this
                            oq_prot_patmat=f"{oq_prot}.patmatmotifs"
                            #Tho patmat command itself
                            patmatmotifs_komandasi = f"patmatmotifs -sequence {yangi_clustalo} -outfile {oq_prot_patmat} -full"
                            print(patmatmotifs_komandasi)  # Check the command before running it
                            #Running the command, fingers crossed
                            try:
                                subprocess.run(patmatmotifs_komandasi, shell=True, check=True)
                                print("patmatmotifs finished successfully.")
                                print(f"The results saved into {oq_prot}.patmatmotifs")
                            except subprocess.CalledProcessError as e:
                                print(f"Error running patmatmotifs: {e}")
                                #To analyse the chemical properties of aminoacid content in the proteins i decided to utilise pepinfo
                            pepinfo = input("Would You like to analyse the chemical properties of your proteins as well? (y/n):")
                            
                            if pepinfo.lower()=="y":
                                with open(f"{oq_prot}.fa", "r") as f:
                                    content = f.read().split(">")[1:]
                                    for fasta in content:
                                        header = fasta.split("\n", 1)[0].strip()
                                        seq_id = header.split()[0]
                                        with open(f"{seq_id}.fasta", "w") as f:
                                            f.write(">" + fasta)
                            #setting the pepinfor directory
                            pepinfo_dir = "/usr/bin/pepinfo"
                            #fasta output be like this
                            fasta_files=glob("*.fasta")
                            #fasta file names going to screen                       
                            print(fasta_files)
                            # sys.exit()
                             ##So i decided it to take the input from the user, as the code above generated the fastas and sends them to screen

                            user_fasta=input("Please enter the file name of any protein sequence \njust created in current dir, which you would like to analyse\nJust copy paste any sequence from the list above!\n>>>")
                            pepinfo_dir = "/usr/bin/pepinfo"
                            #the output will be like this 
                            pep_file=f"{oq_prot}.pepinfo"
                            #and the command itself as this
                            pepinfo_cmd=f"pepinfo -sequence {user_fasta} -outfile {pep_file} -graph pdf"
                            print(pepinfo_cmd)
                            try:
                                subprocess.run(pepinfo_cmd, shell=True, check=True)
                                print("Pepinfo finished successfully.")
                                print("The results saved into pepinfo.pdf")
                                print("Thanks for using this tool! Hope it will be helpful!")
                            except subprocess.CalledProcessError as e:
                                print(f"Error running pepinfo: {e}")
                        #    break 
                
                    
                    
                    
                    # Ask the user if they want to analyze another protein/the same protein for another txid or exit
                    another_search = input("Do you want to analyze another protein/the same protein for another txid or exit? (a/n/e): ")

                    if another_search.lower() == "a":
                        continue
                    elif another_search.lower() == "n":
                        break
                    elif another_search.lower() == "e":
                        print("Thank you for using oq_prot.")
                        break
                    else:
                        print("Invalid input. Please enter 'a' to analyze another protein/the same protein for another txid, 'n' to exit, or 'e' to exit.")


            elif keyingi_etap.lower() == 'c':
                break


    else:
        # identifying whether the user wants to filter the data for specific organizm
          organizm_kirit = input('Do you want to specify the organizm? (y/n): ')
	       ##asking the user whether or not they want spefici fastas for specific organizm within this taxon
          if organizm_kirit.lower() == 'y':
             ##if they say yes, then asking them to enter the name o f that organizm
             organizm = input('Enter organizm name: ')
	         ###Sending query for this specific organims
             zapros += f' AND {organizm}[organism]'
             qidiruv_komandasi = f'esearch -db protein -query "{zapros}" | efetch -format fasta'
             #counting the number of found results
             kolichestvo = os.popen(f'{qidiruv_komandasi} | grep ">" | wc -l').read().strip()
             if kolichestvo == "0":
                 print("Error: No results found for the specified organizm or txid. \nEither there are no fastas presented for this specific protein and/or txid. \nPlease check https://www.ncbi.nlm.nih.gov/ specifying protein {oq_prot}, species {organizm}\nIf these data are presented in NCBI, then Please check your input and try again. \nLets start again!")
                 continue
             else:
                 print(f'Number of results found for {organizm}: {kolichestvo}\n PLEASE KEEP IN MIND. If the number of proteins is too high, further generated results, might be not that much beautiful and informative as You expect')
             #again giving a word to user()
             turlab_yuklash = input(f'Do you want to download all {kolichestvo} {oq_prot} fastas for {organizm}? (y/n): ')
             if turlab_yuklash.lower() == 'y': # downloading filteres fastas
                 print(f"Downloading {oq_prot} fasta files for {organizm} ")
                 yuklab_olamiz = f'{qidiruv_komandasi} > {oq_prot}.fa'
                 os.system(yuklab_olamiz)
                 print(f'All {kolichestvo} {oq_prot} fastas for {organizm} have been downloaded to {oq_prot}.fa')
                 
                 
                 oq_protlar = f"{oq_prot}.fa"                              
                 with open(oq_protlar, "r") as f:
                     lines = f.readlines()  
                 # Creating a dictionary to hold the protein names and their sequences
                 proteins = {}
                 for line in lines:
                     if line.startswith(">"):
                         protein_name = line.strip()[1:]
                         proteins[protein_name] = ""
                     else:
                         proteins[protein_name] += line.strip()

                 #MIN/MAX/MEAN
                         
                         
                 # Creating a pandas DataFrame from the dictionary
                 oq_protcha = pd.DataFrame.from_dict(proteins, orient="index", columns=["Sequence"])
                 # Counting the number of nucleotides for each protein
                 oq_protcha["Nucleotide Count"] = oq_protcha["Sequence"].apply(len)
                 smallest_one = oq_protcha["Nucleotide Count"].min()
                 longest_one = oq_protcha["Nucleotide Count"].max()
                 oq_prot_mean_length = oq_protcha["Nucleotide Count"].mean()
                 print(f"The min. length of these fastas is {smallest_one} bp and the longest sequence comprises of {longest_one} bp. \n Mean length of {oq_prot} sequences is {oq_prot_mean_length} bp")
                                    
                 with open(f'{oq_prot}.fa', 'r') as oq_protlar:
                     fasta_tarkibi = oq_protlar.read()

                     # assessing the the length and completeness of each sequence
                 while True:
                     keyingi_etap = input('What would you like to do next? (a) Look for another protein, (b) Continue for similarity analysis, (c) Exit: ')
                     if keyingi_etap.lower() == 'a':
                         continue
                     elif keyingi_etap.lower() == 'b':
                         fasta_fayli = f"{oq_prot}.fa"
                         clust_fayli = fasta_fayli[:-3] + ".msf"
                         d_file = f"{oq_prot}.dist"
                         # Generating dist matrix with the similarity scores for sequence
                         clustalo_run = f"clustalo -v -i {fasta_fayli} -o {clust_fayli} --outfmt=msf --threads=20 --force --full --distmat-out={d_file}"

                         try:
                             subprocess.run(clustalo_run, shell=True, check=True, capture_output=True, text=True)
                             print("Clustal Omega alignment is finished. Generating the results. Please wait")
                             # Taking out  similarity scores from the generated results
                             similarities_between = []
                             with open(clust_fayli, "r") as f:
                                 for line in f:
                                     if line.startswith('#'):
                                         score = float(line.split(':')[-1].strip())
                                         similarities_between.append(f"{score:.2%}")
                             print(f"Similarity scores between fastas: {similarities_between}")
                             # Distance matrix going to screen and might be saved as file at a later stage
                             with open(d_file, "r") as f:
                                 print(f.read())
                         except subprocess.CalledProcessError as e:
                             print(f"Unfortunately, there was an error running Clustal Omega: {e}")

                         # Visualising the aligned sequences to the screen
                         with open(clust_fayli, "r") as f:
                             print(f.read())

                         #Clarifying whether the user wants the result to be saved?
                         result_desicision = input("Would you like to save the output? (y/n)")
                         #if the user wants to save then downloading results into two different files
                         #the aligned pairs and the distance matrix
                         if result_desicision.lower() == "y":
                             # At this point we are saving the output files with the protein {oq_prot} name
                             ##clustal one
                             yangi_clustalo = f"{oq_prot}.msf"
                             os.rename(clust_fayli, yangi_clustalo)
                             print(f"\nAlignment result saved as {yangi_clustalo}.")
                             ###and the dist matrix
                             new_dist_msf = f"{oq_prot}.matrix"
                             os.rename(d_file, new_dist_msf)
                             print(f"Similarity scores saved as {new_dist_msf}.")
                             # Prompt the user to plot the level of conservation between protein sequences
                             
                             plot_con = input("\nDo you want to plot the level of conservation between protein sequences? (y/n): ")
                        
                             if plot_con.lower() == "y":
                                # Visualising the plot of conservation level on screen
                                 conservation_plotting = f"plotcon -sformat msf {yangi_clustalo} -winsize 16 -graph pdf" 
                                 os.system(conservation_plotting)
                                 # view the generated PDF file on screen
                                 view_pdf = "gs plotcon.pdf"
                                 os.system(view_pdf)
                                 # Display the plot
                                 plt.show()
                                 print("Conservation plot saved as plotcon.pdf")
                               
                                 
                              
                                
                              #PATMATMOTIFS
                                
                               
                                
                               
                                # asking the  user if they want to compare protein motifs
                                 motif_comparison = input("Do you want to identify protein motifs in your set? (y/n): ")
                                 if motif_comparison.lower() == "y":
                                     # In our case the patmatmotifs file is located here
                                     which_patmatmotifs = "/usr/bin/patmatmotifs"
                                     #if the user patmatmotif is not here that better to make sure its installed
                                     oq_prot_patmat=f"{oq_prot}.patmatmotifs"
                                     patmatmotifs_komandasi = f"patmatmotifs -sequence {yangi_clustalo} -outfile {oq_prot_patmat}"
                                     print(patmatmotifs_komandasi)  # Check the command before running it
                                     try:
                                         subprocess.run(patmatmotifs_komandasi, shell=True, check=True)
                                         print("patmatmotifs finished successfully.")
                                         print(f"The results saved into {oq_prot}.patmatmotifs")
                                     except subprocess.CalledProcessError as e:
                                         print(f"Error running patmatmotifs: {e}")
                                         ##########      
                                         
                                    #To analyse the chemical properties of aminoacid content in the proteins i decided to utilise pepinfo
                                     pepinfo = input("Would You like to analyse the chemical properties of your proteins as well? (y/n): ")
                                    
                                     if pepinfo.lower()=="y":
                                         with open(f"{oq_prot}.fa", "r") as f:
                                             content = f.read().split(">")[1:]
                                             for fasta in content:
                                                 header = fasta.split("\n", 1)[0].strip()
                                                 seq_id = header.split()[0]
                                                 with open(f"{seq_id}.fasta", "w") as f:
                                                     f.write(">" + fasta)
                                     #pepinfo is here
                                     pepinfo_dir = "/usr/bin/pepinfo"
                                     #the output will be like this
                                     fasta_files=glob("*.fasta")
                                    
                                     print(fasta_files)
                                     # sys.exit()
                                     ##So i decided it to take the input from the user, as the code above generated the fastas and sends them to screen
                                     user_fasta=input("Please enter the file name of any protein sequence \njust created in current dir, which you would like to analyse\nJust copy paste any sequence from the list above!\n>>>")
                                     print("Just copy paste any sequence from the list above!)")
                                     pepinfo_dir = "/usr/bin/pepinfo"
                                     #the generated outputs be like this and saved as protein name.pepinfo 
                                     pep_file=f"{oq_prot}.pepinfo"
                                     #while being processed with this command
                                     pepinfo_cmd=f"pepinfo -sequence {user_fasta} -outfile {pep_file} -graph pdf"
                                     print(pepinfo_cmd)
                                     try:
                                         subprocess.run(pepinfo_cmd, shell=True, check=True)
                                         print("Pepinfo finished successfully.")
                                         print("The results saved into pepinfo.pdf")
                                         print("Thanks for using this tool! Hope it will be helpful!")
                                     except subprocess.CalledProcessError as e:
                                         print(f"Error running pepinfo: {e}")
                                    #    break 
                                    


                                         
                     ################                 
                           # Asking the user, if they want to analyze another protein/the same protein for another txid or exit
                             another_search = input("Do you want to analyze another protein/the same protein for another txid or exit? (a/n/e): ")
                        
                             if another_search.lower() == "a":
                                 continue
                             elif another_search.lower() == "n":
                                 break
                             elif another_search.lower() == "e":
                                 print("Thank you for using OQSIL! BYE!")
                                 break
                             else:
                                 print("Invalid input. Please enter 'a' to analyze another protein/the same protein for another txid, 'n' to exit, or 'e' to exit.")
                     elif keyingi_etap.lower() == 'c':
                         break

