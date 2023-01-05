# -*- coding: utf-8 -*-
# This code is used to introduce enzyme concentration constraint in GEMs
# by COBRApy and to calculate the parameters that need to be entered
# during the construction of the enzyme-constrained model.
#from warnings import warn

import cobra
import math
import random
import statistics
import sys
import pandas as pd
import json
import re
import copy
from copy import deepcopy
import numpy as np
from typing import Any, Dict, List
from cobra.core import Reaction
from cobra.util.solver import set_objective
from cobra.io.dict import model_to_dict
from urllib.parse import urlencode
from urllib.request import urlopen, Request
from optlang.symbolics import Zero, add
import plotly.graph_objects as go

def standardize_folder(folder: str) -> str:
    """Returns for the given folder path is returned in a more standardized way.

    I.e., folder paths with potential \\ are replaced with /. In addition, if
    a path does not end with / will get an added /.

    Argument
    ----------
    * folder: str ~ The folder path that shall be standardized.
    """
    # Standardize for \ or / as path separator character.
    folder = folder.replace("\\", "/")

    # If the last character is not a path separator, it is
    # added so that all standardized folder path strings
    # contain it.
    if folder[-1] != "/":
        folder += "/"

    return folder

def convert_to_irreversible(model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    Arguments
    ----------
    * model: cobra.Model ~ A Model object which will be modified in place.

    """
    #warn("deprecated, not applicable for optlang solvers", DeprecationWarning)
    reactions_to_add = []
    coefficients = {}
    for reaction in model.reactions:
        # If a reaction is reverse only, the forward reaction (which
        # will be constrained to 0) will be left in the model.
        if reaction.lower_bound < 0 and reaction.upper_bound > 0:
            reverse_reaction = Reaction(reaction.id + "_reverse")
            reverse_reaction.lower_bound = max(0, -reaction.upper_bound)
            reverse_reaction.upper_bound = -reaction.lower_bound
            coefficients[
                reverse_reaction] = reaction.objective_coefficient * -1
            reaction.lower_bound = max(0, reaction.lower_bound)
            reaction.upper_bound = max(0, reaction.upper_bound)
            # Make the directions aware of each other
            reaction.notes["reflection"] = reverse_reaction.id
            reverse_reaction.notes["reflection"] = reaction.id
            reaction_dict = {k: v * -1
                             for k, v in reaction.metabolites.items()}
            reverse_reaction.add_metabolites(reaction_dict)
            reverse_reaction._model = reaction._model
            reverse_reaction._genes = reaction._genes
            reverse_reaction.name = reaction.name+' reverse'
            reverse_reaction.annotation = reaction.annotation
            for gene in reaction._genes:
                gene._reaction.add(reverse_reaction)
            reverse_reaction.subsystem = reaction.subsystem
            reverse_reaction.gene_reaction_rule = reaction.gene_reaction_rule
            reactions_to_add.append(reverse_reaction)
    model.add_reactions(reactions_to_add)
    set_objective(model, coefficients, additive=True)
    
def get_genes_and_gpr(model,gene_outfile,gpr_outfile):
    """Retrieving genes and gene_reaction_rule from GEM.

    Arguments
    ----------
    * model: cobra.Model ~ A genome scale metabolic network model for
        constructing the enzyme-constrained model.

    :return: all genes and gpr in model.
    """
    model_dict = model_to_dict(model, sort=False)
    genes = pd.DataFrame(model_dict['genes']).set_index(['id'])
    genes.to_csv(gene_outfile)
    all_gpr = pd.DataFrame(model_dict['reactions']).set_index(['id'])
    all_gpr.to_csv(gpr_outfile)
    return [genes, all_gpr]

def get_reaction_mw(sbml_path,project_folder,project_name,json_output_file):
        model = cobra.io.read_sbml_model(sbml_path)
        basepath: str = project_folder + project_name
        # READ REACTIONS<->KEGG ID XLSX
        protein_id_mass_mapping: Dict[str, float] = json_load(
            basepath + "_protein_id_mass_mapping.json")
        #subunit_num	1 and 1 and 1 
        reaction_mw={}
        for r in model.reactions:
                if re.search('_num',r.id):
                    r_id=r.id.split('_num')[0]
                else:
                    r_id=r.id            
            #if r_id in reactions_kcat_mapping_database.keys():
                #print(r.id,r.gene_reaction_rule)
                mass_sum = .0
                if re.search(' and ',r.gene_reaction_rule):
                    genelist=r.gene_reaction_rule.split(' and ')
                    for eachgene in genelist:
                        enzyme_unit_number=1
                        if eachgene in protein_id_mass_mapping.keys():
                            mass_sum += protein_id_mass_mapping[eachgene] * enzyme_unit_number
                    #print(mass_sum)
                    reaction_mw[r.id]=mass_sum
                else:  # Single enzyme
                    eachgene=r.gene_reaction_rule
                    enzyme_unit_number = 1
                    if eachgene in protein_id_mass_mapping.keys():
                        #print(protein_id_mass_mapping[eachgene] * enzyme_unit_number)
                        reaction_mw[r.id]=protein_id_mass_mapping[eachgene] * enzyme_unit_number
        json_write(json_output_file, reaction_mw) 
                        
# PUBLIC FUNCTIONS
def get_reaction_kcat_mw(model_file: str,project_folder: str, project_name: str,enzyme_unit_number_file: str,
                                       type_of_default_kcat_selection: str = "median") -> None:
    """Adds proteomic constraints according to sMOMENT to the given stoichiometric model and stores it as SBML.

    Arguments
    ----------

    * model: cobra.Model ~ A cobra Model representation of the metabolic network. This model will
      be changed using cobrapy functions in order to add the proteomic constraints.
    * project_folder: str ~ The folder in which the spreadsheets and JSONs with the model's supplemental
      data can be found.
    * project_name: str ~ The sMOMENTed model creation's name, which will be added at the beginning
      of the created SBML's name.
    * type_of_default_kcat_selection: str ~ The type of selection of default kcat values. Can be "mean",
      "median" or "random". Is "median" by default.

    Output
    ----------
    An SBML in the given folder with the given name, which describes the given stoichiometric model
    enhanced by the protein constraint introduction with this function.
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # Set folder path for newly created SBML and name for the reaction ID addition (added at the end,
    # and used in order to have a programatically convinient way to separate additions such as 'reverse'
    # from the 'actual' reaction ID).
    basepath: str = project_folder + project_name
    id_addition: str = "_num"
    # READ REACTIONS<->KEGG ID XLSX
    protein_id_mass_mapping: Dict[str, float] = json_load(
        basepath + "_protein_id_mass_mapping.json")
    if enzyme_unit_number_file != 'none':
        enzyme_unit_number=json_load(enzyme_unit_number_file) 
    # Make model irreversible, separating all reversible reactions to which a gene rule is given
    # in order to save some reactions.
    model = cobra.io.read_sbml_model(model_file)
    convert_to_irreversible(model)
    #split isoenzyme
    model = isoenzyme_split(model)
    
    # Read reaction <-> kcat mapping :-)
    reactions_kcat_mapping_database = json_load(
        basepath + "_reactions_kcat_mapping_combined.json")

    # sMOMENT :D
    # Get all kcats which are not math.nan and calculate the median of them, which will be used as default kcat
    all_kcats = [x["forward"] for x in reactions_kcat_mapping_database.values()] + \
                [x["reverse"] for x in reactions_kcat_mapping_database.values()]
    all_kcats = [x for x in all_kcats if not math.isnan(x)]

    if type_of_default_kcat_selection == "median":
        default_kcat = statistics.median(all_kcats)
    elif type_of_default_kcat_selection == "mean":
        default_kcat = statistics.mean(all_kcats)
    elif type_of_default_kcat_selection == "max":
        default_kcat = np.max(all_kcats)
    elif type_of_default_kcat_selection == "random":
        default_kcat = random.choice(all_kcats)
    else:
        default_kcat = 'Null'

    print(f"Default kcat is: {default_kcat}")

    # Get all reaction IDs of the given model
    model_reaction_ids = [x.id for x in model.reactions]

    # Main loop :D, add enzyme constraints to reactions \o/
    reaction_kcat_mw={}
    for model_reaction_id in model_reaction_ids:
        # Get the reaction and split the ID at the ID addition
        reaction = model.reactions.get_by_id(model_reaction_id)
        splitted_id = reaction.id.split(id_addition)

        # If the reaction has no name, ignore it
        if splitted_id[0] == "":
            continue
        # Take the reaction ID from the first part of the split
        reaction_id = splitted_id[0]
        # Remove GPRSPLIT name addition from reactions with measured protein concentrations
        if "_GPRSPLIT_" in reaction_id:
            reaction_id = reaction_id.split("_GPRSPLIT_")[0]

        # Retrieve the reaction's forward and reverse kcats from the given reaction<->kcat database
        if re.search('_reverse',reaction_id):
            reaction_id=reaction_id.split('_reverse')[0]
        if reaction_id not in reactions_kcat_mapping_database.keys():
            continue
        forward_kcat = reactions_kcat_mapping_database[reaction_id]["forward"]
        reverse_kcat = reactions_kcat_mapping_database[reaction_id]["reverse"]

        # If the given reaction<->kcat database contains math.nan as the reaction's kcat,
        # set the default kcat as math.nan means that no kcat could be found.
        if math.isnan(forward_kcat):
            forward_kcat = default_kcat
        if math.isnan(reverse_kcat):
            reverse_kcat = default_kcat

        # Add the given forward or reverse kcat is the reaction was
        # splitted due to its reversibility.
        # If the reaction is not splitted, add the forward kcat (this
        # is the only possible direction for non-splitted=non-reversible
        # reactions)
        if model_reaction_id.endswith(id_addition + "forward"):
            reaction_kcat = forward_kcat
        elif model_reaction_id.endswith(id_addition + "reverse"):
            reaction_kcat = reverse_kcat
        else:
            reaction_kcat = forward_kcat

        reaction_kcat_mw[model_reaction_id]={}
        if reaction_kcat=='Null':
            continue
        reaction_kcat_mw[model_reaction_id]['kcat']=reaction_kcat
        
        #MW
        #subunit_num	1 and 1 and 1 
        reaction_mw={}
        for r in model.reactions:
                #print(r.id)
                if re.search('_num',r.id):
                    r_id=r.id.split('_num')[0]
                else:
                    r_id=r.id            
            #if r_id in reactions_kcat_mapping_database.keys():
                #print(r.id,r.gene_reaction_rule)
                mass_sum = .0
                if re.search(' and ',r.gene_reaction_rule):
                    genelist=r.gene_reaction_rule.split(' and ')
                    for eachgene in genelist:
                        if eachgene in protein_id_mass_mapping.keys():
                            if enzyme_unit_number_file != 'none':
                                if r.id in enzyme_unit_number.keys():
                                    if eachgene in enzyme_unit_number[r.id].keys():
                                        mass_sum += protein_id_mass_mapping[eachgene] * int(enzyme_unit_number[r.id][eachgene])
                            else:
                                #print(eachgene)
                                mass_sum += protein_id_mass_mapping[eachgene]
                    #print(mass_sum)
                    if mass_sum>0:
                        reaction_mw[r.id]=mass_sum
                else:  # Single enzyme
                    eachgene=r.gene_reaction_rule
                    #enzyme_unit_number = 1
                    if eachgene in protein_id_mass_mapping.keys():
                        #print(protein_id_mass_mapping[eachgene] * enzyme_unit_number)
                        #reaction_mw[r.id]=protein_id_mass_mapping[eachgene]['mw'] * enzyme_unit_number
                        if enzyme_unit_number_file != 'none':
                            if r.id in enzyme_unit_number.keys():
                                if eachgene in enzyme_unit_number[r.id].keys(): 
                                    reaction_mw[r.id]=protein_id_mass_mapping[eachgene] * int(enzyme_unit_number[r.id][eachgene])
                        else:
                            reaction_mw[r.id]=protein_id_mass_mapping[eachgene]
                # isoform enzyme's mw should be modified
                if re.search('ASPK',r.id):
                    reaction_mw[r.id] = 126584
                if re.search('ACGS',r.id):
                    reaction_mw[r.id] = 79428
                if re.search('ASP1DC',r.id):
                    reaction_mw[r.id] = 56588   
        #print(model_reaction_id,reaction_mw.keys())
        
        if model_reaction_id in reaction_mw.keys():
            # print(reaction_mw.keys())
            #print(model_reaction_id,reaction_mw[model_reaction_id])
            reaction_kcat_mw[model_reaction_id]['MW']=reaction_mw[model_reaction_id]
            
            reaction_kcat_mw[model_reaction_id]['kcat_MW']=reaction_kcat_mw[model_reaction_id]['kcat']*3600000/reaction_mw[model_reaction_id]
        
    reaction_kcat_mw_df = pd.DataFrame(reaction_kcat_mw)
    reaction_kcat_mw_df_T=reaction_kcat_mw_df.T
    reaction_kcat_mw_df_T_select=reaction_kcat_mw_df_T[abs(reaction_kcat_mw_df_T['kcat_MW'])>0]
    reaction_kcat_mw_df_T_select.to_csv(project_folder + 'reaction_kcat_MW.csv')
    
def isoenzyme_split(model):
    """Split isoenzyme reaction to mutiple reaction

    Arguments
    ----------
    * model: cobra.Model.
    
    :return: new cobra.Model.
    """  
    for r in model.reactions:
        if re.search(" or ", r.gene_reaction_rule):
            rea = r.copy()
            gene = r.gene_reaction_rule.split(" or ")
            for index, value in enumerate(gene):
                if index == 0:
                    r.id = r.id + "_num1"
                    r.gene_reaction_rule = value
                else:
                    r_add = rea.copy()
                    r_add.id = rea.id + "_num" + str(index+1)
                    r_add.gene_reaction_rule = value
                    model.add_reaction(r_add)
    for r in model.reactions:
        r.gene_reaction_rule = r.gene_reaction_rule.strip("( )")
    return model

def trans_model2enz_json_model_split_isoenzyme(model_file, reaction_kcat_mw_file, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json mode with  
    enzyme concentration constraintat.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 

    """
    model = cobra.io.read_sbml_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    reaction_kcat_mw = pd.read_csv(reaction_kcat_mw_file, index_col=0)
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        #if re.search('_num',reaction_id):
        #    reaction_id=reaction_id.split('_num')[0]
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    json_write(json_output_file, dictionary_model)
    
def get_enzyme_constraint_model(json_model_file):
    """using enzyme concentration constraint
    json model to create a COBRApy model.

    Arguments
    ----------
    * json_model_file: json Model file.

    :return: Construct an enzyme-constrained model.
    """

    dictionary_model = json_load(json_model_file)
    model = cobra.io.json.load_json_model(json_model_file)

    coefficients = dict()
    for rxn in model.reactions:
        for eachr in dictionary_model['reactions']:
            if rxn.id == eachr['id']:
                if eachr['kcat_MW']:
                    coefficients[rxn.forward_variable] = 1 / float(eachr['kcat_MW'])
                break

    lowerbound = dictionary_model['enzyme_constraint']['lowerbound']
    upperbound = dictionary_model['enzyme_constraint']['upperbound']
    #print(upperbound)
    constraint = model.problem.Constraint(0, lb=lowerbound, ub=upperbound)
    model.add_cons_vars(constraint)
    model.solver.update()
    constraint.set_linear_coefficients(coefficients=coefficients)
    return model

def get_fluxes_detail_in_model(model,model_pfba_solution,fluxes_outfile,json_model_file):
    """Get the detailed information of each reaction

    Arguments
    ----------
    * model: cobra.Model.
    * fluxes_outfile: reaction flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.

    :return: fluxes, kcat, MW and kcat_MW in dataframe.
    """
    dictionary_model = json_load(json_model_file)
    model_pfba_solution = model_pfba_solution.to_frame()
    model_pfba_solution_detail = pd.DataFrame()
    for index, row in model_pfba_solution.iterrows():
        reaction_detail = model.reactions.get_by_id(index)
        model_pfba_solution_detail.loc[index, 'fluxes'] = row['fluxes']
        for eachreaction in dictionary_model['reactions']:
            if index ==eachreaction['id']:
                if 'annotation' in eachreaction.keys():
                    if 'ec-code' in eachreaction['annotation'].keys():
                        if isinstance (eachreaction['annotation']['ec-code'],list):
                            model_pfba_solution_detail.loc[index, 'ec-code'] = (',').join(eachreaction['annotation']['ec-code'])
                        else:
                            model_pfba_solution_detail.loc[index, 'ec-code'] = eachreaction['annotation']['ec-code']    
                if 'kcat_MW' in eachreaction.keys():
                    if eachreaction['kcat_MW']:
                        model_pfba_solution_detail.loc[index, 'kcat_MW'] = eachreaction['kcat_MW']
                        model_pfba_solution_detail.loc[index, 'E'] = float(row['fluxes'])/float(eachreaction['kcat_MW'])
                break
        model_pfba_solution_detail.loc[index, 'equ'] = reaction_detail.reaction
    print('Enzyme cost total is:'+str(np.sum(model_pfba_solution_detail['E'])))
    model_pfba_solution_detail.to_csv(fluxes_outfile)
    return model_pfba_solution_detail

def json_load(path):
    """Loads the given JSON file and returns it as dictionary.

    Arguments
    ----------
    * path: The path of the JSON file
    """
    with open(path) as f:
        dictionary = json.load(f)
    return dictionary

def json_write(path, dictionary):
    """Writes a JSON file at the given path with the given dictionary as content.

    Arguments
    ----------
    * path:   The path of the JSON file that shall be written
    * dictionary: The dictionary which shalll be the content of
      the created JSON file
    """
    json_output = json.dumps(dictionary, indent=4)
    with open(path, "w", encoding="utf-8") as f:
        f.write(json_output)
        
def GENENAME_2_ACC_from_uniprot(query,outfile):
    #print(' '.join(query).replace('511145.',''))
    url = 'https://legacy.uniprot.org/uploadlists/'
    params = {
        'from': 'GENENAME',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,database(PDB)'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    outFile = open(outfile,'w') 
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()
    
def GENENAME_2_ACC_from_uniprot_byID(query,outfile):
    #print(' '.join(query).replace('511145.',''))
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
        'from': 'ID',
        'to': 'ACC',
        'format': 'tab',
        'query': ' '.join(query),
        'columns':'id,entry name,protein names,genes,organism,ec,mass,database(PDB)'
    }
    data = urlencode(params).encode()
    request = Request(url, data)
    # Please set your email address here to help us debug in case of problems.
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urlopen(request)
    page = response.read()
    outFile = open(outfile,'w') 
    namesRegex = re.compile(r'yourlist:(.*)\n')
    outFile.write(namesRegex.sub('Gene ID\n',page.decode('utf-8')))
    #print(namesRegex.sub('Protein AC\t',page.decode('utf-8')))
    outFile.close()    
    
def calculate_f(uni_model_gene_list, gene_abundance_file, gene_mw_file, gene_mw_colname,gene_abundance_colname):
    """Calculating f (the mass fraction of enzymes that are accounted
    in the model out of all proteins) based on the protein abundance
    which can be obtained from PAXdb database.

    Arguments
    ----------
    * genes: All the genes in the model.
    * gene_abundance_file: The protein abundance of each gene
     in the E. coli genome.
    * gene_mw_file: The molecular weight of the protein expressed by each gene.

    :return: The enzyme mass fraction f.
    """
    gene_abundance = pd.read_csv(gene_abundance_file, index_col=0)
    gene_mw = pd.read_csv(gene_mw_file, sep='\t', index_col=gene_mw_colname)
    enzy_abundance = 0
    pro_abundance = 0
    for gene_i in gene_abundance.index:
        if gene_i in gene_mw.index:
            if isinstance(gene_mw.loc[gene_i,'Mass'],str):
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(gene_mw.loc[gene_i,'Mass'].replace(',',''))
            else:
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * int(gene_mw.loc[gene_i,'Mass'][0].replace(',',''))
            pro_abundance += abundance
            if gene_i in uni_model_gene_list:
                enzy_abundance += abundance
    f = enzy_abundance/pro_abundance
    return f

def calculate_f_special(gene_abundance_file,modelgene2uniprot_file,paxdbgene2uniprot_file,gene_abundance_colname):
    """Calculating f (the mass fraction of enzymes that are accounted
    in the model out of all proteins) based on the protein abundance
    which can be obtained from PAXdb database.

    Arguments
    ----------
    * genes: All the genes in the model.
    * gene_abundance_file: The protein abundance of each gene
     in the E. coli genome.
    * gene_mw_file: The molecular weight of the protein expressed by each gene.

    :return: The enzyme mass fraction f.
    """
    gene_abundance = pd.read_csv(gene_abundance_file, index_col='gene_id',sep='\t')
    model_gene = pd.read_csv(modelgene2uniprot_file, sep='\t', index_col='Gene ID')
    paxdb_gene = pd.read_csv(paxdbgene2uniprot_file, sep='\t', index_col='Gene ID')
    paxdb_gene=paxdb_gene[paxdb_gene['Organism']=='Corynebacterium glutamicum (strain ATCC 13032 / DSM 20300 / BCRC 11384 / JCM 1318 / LMG 3730 / NCIMB 10025)']
    model_gene=model_gene[model_gene['Organism']=='Corynebacterium glutamicum (strain ATCC 13032 / DSM 20300 / BCRC 11384 / JCM 1318 / LMG 3730 / NCIMB 10025)']
    enzy_abundance = 0
    pro_abundance = 0
    for gene_i in gene_abundance.index:
        #print(gene_i)
        if gene_i in list(paxdb_gene.index):
            #print(gene_i)
            if str(paxdb_gene.loc[gene_i,'Mass']) !='nan':
                abundance = gene_abundance.loc[gene_i, gene_abundance_colname] * float(paxdb_gene.loc[gene_i,'Mass'].replace(',',''))
                #print(abundance)
                pro_abundance += abundance
                if gene_i in list(model_gene.index):
                    #print(gene_i)
                    enzy_abundance += abundance
    #print(enzy_abundance,pro_abundance)
    f = enzy_abundance/pro_abundance
    return f

def get_model_substrate_obj(use_model):
    ATPM='No' 
    substrate_list=[]
    concentration_list=[]
    EX_exclude_reaction_list=['EX_pi_e','EX_h_e','EX_fe3_e','EX_mn2_e','EX_co2_e','EX_fe2_e','EX_h2_e','EX_zn2_e',\
                             'EX_mg2_e','EX_ca2_e','EX_so3_e','EX_ni2_e','EX_no_e','EX_cu2_e','EX_hg2_e','EX_cd2_e',\
                             'EX_h2o2_e','EX_h2o_e','EX_no2_e','EX_nh4_e','EX_so4_e','EX_k_e','EX_na1_e','EX_o2_e',\
                             'EX_o2s_e','EX_ag_e','EX_cu_e','EX_so2_e','EX_cl_e','EX_n2o_e','EX_cs1_e','EX_cobalt2_e']
    EX_exclude_reaction_list=EX_exclude_reaction_list+[i+'_reverse' for i in EX_exclude_reaction_list]
    for r in use_model.reactions:
        if r.objective_coefficient == 1:
            obj=r.id #Product name
        #elif not r.lower_bound==0 and not r.lower_bound==-1000 and not r.lower_bound==-999999 and abs(r.lower_bound)>0.1:#排除很小的值
        elif not r.upper_bound==0 and not r.upper_bound==1000 and not r.upper_bound==999999 and abs(r.upper_bound)>0.1:#排除很小的值
            #print(r.id,r.upper_bound,r.lower_bound)
            if r.id=='ATPM':
                if r.upper_bound>0:
                    ATPM='Yes' #ATP maintenance requirement
            elif r.id not in EX_exclude_reaction_list:
                #print(r.id,r.upper_bound,r.lower_bound)
                #substrate=r.id #Substrate name
                substrate_list.append(r.id)
                #concentration=r.upper_bound #Substrate uptake rate  
                concentration_list.append(r.upper_bound)
    return(obj,substrate_list,concentration_list,ATPM)

def parse_sabio_rk_for_eclist(ec_numbers_list: List[str], json_output_path: str, bigg_id_name_mapping_path: str) -> None:
    """Retrieves kcats from SABIO-RK for the given model and stores it in a JSON for the given model in the given path.

    Algorithm
    ----------
    Using the SABIO-RK REST API (as of 2019/30/04, it is explained under
    http://sabiork.h-its.org/layouts/content/docuRESTfulWeb/RESTWebserviceIntro.gsp),


    Arguments
    ----------
    * eclist: List[str] ~ eclist.
    * json_output_path: str ~ The path of the JSON that shall be created

    Output
    ----------
    * A JSON in the given project folder with the following structure:
    <pre>
        {
            "$EC_NUMBER_OR_KEGG_REACTION_ID": {
                "$SUBSTRATE_WITH_BIGG_ID_1": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                },
                (...),
                "REST": {
                    "$ORGANISM_1": [
                        $kcat_1,
                        (...)
                        $kcat_n,
                    ]
                }
            }
            (...),
        }
    </pre>
    'REST' stands for a substrate without found BIGG ID.
    """
    # GET KCATS FOR EC NUMBERS
    ec_number_kcat_mapping = get_ec_number_kcats_wildcard_search(
        ec_numbers_list, bigg_id_name_mapping_path)

    json_write(json_output_path, ec_number_kcat_mapping)
      
def get_protein_mass_mapping_from_local(sbml_path: str, project_folder: str, project_name: str,uniprot_data_file: str) -> None:
    """Returns a JSON with a mapping of protein IDs as keys, and as values the protein mass in kDa.

    The protein masses are calculated using the amino acid sequence from UniProt (retrieved using
    UniProt's REST API).

    Arguments
    ----------
    * model: cobra.Model ~ The model in the cobrapy format
    * project_folder: str ~ The folder in which the JSON shall be created
    * project_name: str ~ The beginning of the JSON's file name
    * uniprot_data_file: str ~ The gene information obtained from uniprot
    Output
    ----------
    A JSON file with the path project_folder+project_name+'_protein_id_mass_mapping.json'
    and the following structure:
    <pre>
    {
        "$PROTEIN_ID": $PROTEIN_MASS_IN_KDA,
        (...),
    }
    </pre>
    """
    # Standardize project folder
    project_folder = standardize_folder(project_folder)

    # The beginning of the created JSON's path :D
    basepath: str = project_folder + project_name
    model= cobra.io.read_sbml_model(sbml_path)
    # GET UNIPROT ID - PROTEIN MAPPING
    uniprot_id_protein_id_mapping: Dict[str, List[str]] = {}
    for gene in model.genes:
        # Without a UniProt ID, no mass mapping can be found
        if "uniprot" not in gene.annotation:
            continue
        uniprot_id = gene.annotation["uniprot"]
        if uniprot_id in uniprot_id_protein_id_mapping.keys():
            uniprot_id_protein_id_mapping[uniprot_id].append(gene.id)
        else:
            uniprot_id_protein_id_mapping[uniprot_id] = [gene.id]

    # GET UNIPROT ID<->PROTEIN MASS MAPPING
    uniprot_id_protein_mass_mapping = json_load(uniprot_data_file)
    
    # Create the final protein ID <-> mass mapping
    protein_id_mass_mapping: Dict[str, float] = {}
    for uniprot_id in list(uniprot_id_protein_mass_mapping.keys()):
        try:
            protein_ids = uniprot_id_protein_id_mapping[uniprot_id]
        except Exception:
            #print(f"No mass found for {uniprot_id}!")
            continue
        for protein_id in protein_ids:
            protein_id_mass_mapping[protein_id] = uniprot_id_protein_mass_mapping[uniprot_id]

    # Write protein mass list JSON :D
    #print("Protein ID<->Mass mapping done!")
    json_write(basepath+"_protein_id_mass_mapping.json", protein_id_mass_mapping)
    
def change_enz_model_by_enz_usage(enz_ratio,json_model_path, reaction_flux_file, EC_max_file, reaction_kcat_mw, need_change_reaction_list, changed_reaction_list,f, ptot, sigma, lowerbound, upperbound, json_output_file):

    """Get new enzyme model using enzyme usage to calibration

    Arguments
    ----------
    * enz_ratio: enzyme ratio which needed change.
    * json_model_path: The file storing json model.
    * reaction_flux_file: reaction-flux file.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_enz_usage_file： enzyme usage of each reaction.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * model_file: cobra model.
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint.  
    * json_output_file: json file store json model
    * reaction_mw_outfile: changed file stored reaction kcat/mw.

    :return: new enzyme model
    """ 
    reaction_fluxes = pd.read_csv(reaction_flux_file, index_col=0)
    reaction_fluxes['enz ratio'] = reaction_fluxes['E']/np.sum(reaction_fluxes['E'])
    reaction_fluxes=reaction_fluxes.sort_values(by="enz ratio", axis=0, ascending=False)
    i=0
    select_reaction = reaction_fluxes.index[0]
    # if select_reaction == 'PPNDH':
    #     select_reaction = reaction_fluxes.index[1]
    while (select_reaction in need_change_reaction_list):
        i=i+1
        #print(i)
        select_reaction = reaction_fluxes.index[i+1]
        
    print('Need changing reaction: ')
    print(select_reaction)
    [need_change_reaction_list, changed_reaction_list,reaction_kcat_mw] = adj_reaction_kcat_by_database(json_model_path,[select_reaction], need_change_reaction_list, changed_reaction_list, EC_max_file, reaction_kcat_mw)
    print('Changed reaction: ')
    print(changed_reaction_list)

    adj_trans_model2enz_model(json_model_path, reaction_kcat_mw, f, ptot, sigma, lowerbound, upperbound, json_output_file)

    enz_model = get_enzyme_constraint_model(json_output_file)
    print('Enzyme cost total is:'+str(np.sum(reaction_fluxes['E'])))
    return (enz_model,reaction_kcat_mw,need_change_reaction_list, changed_reaction_list)

def adj_reaction_kcat_by_database(json_model_path,select_reactionlist, need_change_reaction_list, changed_reaction_list,EC_max_file, reaction_kcat_mw):
    """Use the kcat in database to change reaction kcat in model

    Arguments
    ----------
    * json_model_path: The file storing json model.
    * select_reactionlist: reaction list need to change.
    * kcat_database_combined_file: combined kcat file got from autoPACMEN.
    * reaction_kcat_mw_file: reaction kcat/mw file.
    * reaction_kapp_change_file: changed file stored reaction kcat/mw.

    :return: a dataframe stored new reaction kcat/mw .
    """
    Brenda_sabio_combined_select = json_load(EC_max_file)
    
    json_model=cobra.io.load_json_model(json_model_path)
    for eachreaction in select_reactionlist:
        need_change_reaction_list.append(eachreaction)
        select_reaction = json_model.reactions.get_by_id(eachreaction)
        if "ec-code" in select_reaction.annotation.keys():
            ec_number = select_reaction.annotation["ec-code"]
            kcat_max_list = []
            if isinstance(ec_number, str):
                if ec_number in Brenda_sabio_combined_select.keys():
                    reaction_kcat_max = Brenda_sabio_combined_select[ec_number]['kcat_max']
                    if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max:
                        reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max#h_1
                        reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                        changed_reaction_list.append(eachreaction) 
            else:
                for eachec in ec_number:
                    if eachec in Brenda_sabio_combined_select.keys():
                        kcat_max_list.append(Brenda_sabio_combined_select[eachec]['kcat_max'])
                reaction_kcat_max = np.max(kcat_max_list)     
                if reaction_kcat_mw.loc[eachreaction, 'kcat'] < reaction_kcat_max:
                    reaction_kcat_mw.loc[eachreaction,'kcat'] = reaction_kcat_max
                    reaction_kcat_mw.loc[eachreaction, 'kcat_MW'] = reaction_kcat_max * 3600*1000/reaction_kcat_mw.loc[eachreaction, 'MW']
                    changed_reaction_list.append(eachreaction)    
                
    return(need_change_reaction_list,changed_reaction_list,reaction_kcat_mw)

def adj_trans_model2enz_model(model_file, reaction_kcat_mw, f, ptot, sigma, lowerbound, upperbound, json_output_file):
    """Tansform cobra model to json mode with  
    enzyme concentration constraintat.

    Arguments
    ----------
    * model_file:   The path of sbml model
    * reaction_kcat_mw_file: The path of storing kcat/MW value of the enzyme catalyzing each
     reaction in the GEM model
    * f: The enzyme mass fraction 
    * ptot: The total protein fraction in cell.  
    * sigma: The approximated average saturation of enzyme. 
    * lowerbound:  Lowerbound  of enzyme concentration constraint. 
    * upperbound:  Upperbound  of enzyme concentration constraint. 

    """
    if re.search('\.xml',model_file):
        model = cobra.io.read_sbml_model(model_file)
    elif re.search('\.json',model_file):
        model = cobra.io.json.load_json_model(model_file)
    convert_to_irreversible(model)
    model = isoenzyme_split(model)
    model_name = model_file.split('/')[-1].split('.')[0]
    json_path = "./model/%s_irreversible.json" % model_name
    cobra.io.save_json_model(model, json_path)
    dictionary_model = json_load(json_path)
    dictionary_model['enzyme_constraint'] = {'enzyme_mass_fraction': f, 'total_protein_fraction': ptot,
                                             'average_saturation': sigma, 'lowerbound': lowerbound, 'upperbound': upperbound}
    # Reaction-kcat_mw file.
    # eg. AADDGT,49389.2889,40.6396,1215.299582180927
    for eachreaction in range(len(dictionary_model['reactions'])):
        reaction_id = dictionary_model['reactions'][eachreaction]['id']
        #if re.search('_num',reaction_id):
        #    reaction_id=reaction_id.split('_num')[0]
        if reaction_id in reaction_kcat_mw.index:
            dictionary_model['reactions'][eachreaction]['kcat'] = reaction_kcat_mw.loc[reaction_id, 'kcat']
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = reaction_kcat_mw.loc[reaction_id, 'kcat_MW']
        else:
            dictionary_model['reactions'][eachreaction]['kcat'] = ''
            dictionary_model['reactions'][eachreaction]['kcat_MW'] = ''
    json_write(json_output_file, dictionary_model)
    
def draw_cdf_fig(data_cdf_data,output_file,x_name,y_name,y_index,nticks):
    trace0 = go.Scatter(x=data_cdf_data,y=y_index,mode='lines',marker={'color': 'blue'},xaxis='x2',yaxis="y2")
    trace1 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines')
    data1 = [trace0, trace1]
    layout = go.Layout(plot_bgcolor='lightgrey',
            xaxis=dict(title=dict(text=x_name,font=dict(size=20, family='Times New Roman')),
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',ticks='inside',tickcolor='black',nticks=nticks,zeroline=False,
                   showexponent = 'all',exponentformat =  "power"),
            xaxis2=dict(linecolor='black',showticklabels=False,type="log",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   rangemode="tozero",overlaying='x', side='top', zeroline=False,
                   ),
            yaxis=dict(title=dict(text=y_name,font=dict(size=20, family='Times New Roman')),range=[0, 1],
                   showgrid=False,zeroline=False,rangemode="tozero", tickfont=dict(color='black', size=20, family='Times New Roman'),
                   ticks='inside',tickcolor='black',linecolor='black'),
            yaxis2=dict(range=[0, 1],linecolor='black',showgrid=False,zeroline=False,tickfont=dict(color='black', size=20, family='Times New Roman'),
                    showticklabels=False,overlaying='y',side='right'),
            showlegend=False,height=450,width=750,margin=go.layout.Margin(l=10, r=10, b=10, t=10))

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5,line_width=2,line_color="orange")
    fig.write_image(output_file)
    return fig

def draw_cdf_fig_kcat(data_cdf_data,output_file,x_name,y_name,y_index,nticks):
    trace0 = go.Scatter(x=data_cdf_data,y=y_index,mode='lines',marker={'color': 'blue'})
    trace1 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',line={'color': 'blue', 'width': 3},xaxis='x2',yaxis="y2")
    data1 = [trace0, trace1]
    layout = go.Layout(plot_bgcolor='lightgrey',xaxis=dict(title=dict(text=x_name,font=dict(size=20, family='Times New Roman')),
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',ticks='inside',tickcolor='black',zeroline=False,
                   showexponent = 'all',exponentformat =  "power", gridcolor="yellow"),
            xaxis2=dict(linecolor='black',showticklabels=False,type="log",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   rangemode="tozero",overlaying='x', side='top',nticks=nticks, zeroline=False,
                   showexponent = 'all', exponentformat =  "power",gridcolor="white"),
            yaxis=dict(title=dict(text=y_name,font=dict(size=20, family='Times New Roman')),range=[0, 1],
                   showgrid=False,zeroline=False,rangemode="tozero", tickfont=dict(color='black', size=20, family='Times New Roman'),
                   ticks='inside',tickcolor='black',linecolor='black'),
            yaxis2=dict(range=[0, 1],linecolor='black',showgrid=False,zeroline=False,tickfont=dict(color='black', size=20, family='Times New Roman'),
                    showticklabels=False,overlaying='y',side='right'),
            showlegend=False,height=450,width=750,margin=go.layout.Margin(l=10, r=10, b=10, t=10))

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5,line_width=2,line_color="orange")
    fig.write_image(output_file)
    return fig

def get_PhPP_data(model_file,model_type,obj,number,outputfile):
    if model_type=='GEM':
        cg_model = cobra.io.json.load_json_model(model_file)
    else:
        cg_model=get_enzyme_constraint_model(model_file)
    objlist = []
    glclist = []
    o2list = []
    exlisto2 = []
    exlist = list(np.linspace(0,10,number))
    exlisto2 = list(np.linspace(0,10,number))
    print(exlist)

    exlistn = []
    exlistno2 = []
    for i in exlist:
        i = format(i,'.4f')
        exlistn.append(float(i))
    exlistm=deepcopy(exlistn)
    for i in exlisto2:
        i = format(i,'.4f')
        exlistno2.append(float(i))
    exlistmo2=deepcopy(exlistno2)    
    
    exlistn.insert(0,0.0)
    exlistno2.insert(0,0.0)

    ectest = pd.DataFrame(exlistn)
    o2df = pd.DataFrame(exlistno2).T
    df = pd.concat([o2df,ectest],axis=0)
    df = df.reset_index(drop = True)
    df1 = df.drop(1,axis=0)
    df1 = df1.reset_index(drop = True)
    df1
    k=1
    v=1
    for i in exlistmo2:
        condi = i
        for j in exlistm:
            condj = -j

            cg_model.reactions.get_by_id('EX_o2_e').bounds=(0,0)
            cg_model.reactions.get_by_id('EX_o2_e_reverse').bounds=(0,condi)
            cg_model.reactions.get_by_id('EX_glc_e').bounds=(condj,0)
            cg_model.reactions.get_by_id('EX_glc_e_reverse').bounds=(0,0)

            cg_model.objective=obj
            enz_model_pfba_solution = cobra.flux_analysis.pfba(cg_model)
            # print(k,v)
            df1.iloc[k,v] = enz_model_pfba_solution.fluxes[obj]
            k = k+1
            if k == number+1:
                k=1
                v = v+1
        if v == number+1:
            break
#     print(df1)
    df1.to_csv(outputfile)
    return(df1)

def draw_cdf_fig_mw(data_cdf_data,output_file,x_name,y_name,y_index,nticks):
    trace0 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',line={'color': 'blue', 'width': 3},xaxis='x2',yaxis="y2")
    trace1 = go.Scatter(x=data_cdf_data,y=y_index,marker={'color': 'blue'},mode='lines',line={'color': 'blue', 'width': 3},xaxis='x2',yaxis="y2")
    data1 = [trace0,trace1]
    layout = go.Layout(plot_bgcolor='lightgrey',
            xaxis=dict(title=dict(text=x_name,font=dict(size=20, family='Times New Roman')),range=[0.89, 3.3],
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),
                   linecolor='black',ticks='inside',tickcolor='black',zeroline=False,showgrid=False,
                   showexponent = 'all',exponentformat =  "power"),
            xaxis2=dict(showticklabels=False,
                   type="log",rangemode="tozero",tickfont=dict(color='black', size=20, family='Times New Roman'),range=[0.89, 3.3],
                   linecolor='black',overlaying='x', side='top',tickcolor='black', zeroline=False,nticks = nticks,
                   showexponent = 'all', exponentformat =  "power"),
            yaxis=dict(title=dict(text=y_name,font=dict(size=20, family='Times New Roman')),range=[0, 1],
                   showgrid=False,zeroline=False,rangemode="tozero", tickfont=dict(color='black', size=20, family='Times New Roman'),
                   ticks='inside',tickcolor='black',linecolor='black'),
            yaxis2=dict(range=[0, 1],linecolor='black',showgrid=False,zeroline=False,tickfont=dict(color='black', size=20, family='Times New Roman'),
                    showticklabels=False,overlaying='y',side='right'),
            showlegend=False,height=450,width=750,margin=go.layout.Margin(l=10, r=10, b=10, t=10))

    fig = go.Figure(data1, layout=layout)
    fig.add_hline(y=0.5,line_width=2,line_color="orange")
    fig.write_image(output_file)
    return fig 

def get_min_enzyme_cost(model, dict_coeff):
    """Get model flux using Minimum enzyme cost algorithm

    Arguments
    ----------
    * model: cobra model.
    * dict_coeff: {reaction ID: coeffient}.
    
    :return: cobra solution.
    """
    with model:
        bounds = (model.slim_optimize(), model.slim_optimize())
        cons_obj = model.problem.Constraint(
            model.objective.expression,
            lb=min(bounds), ub=max(bounds))
        model.add_cons_vars(cons_obj)

        dict_obj = dict()
        for r in model.reactions:
            if r.id in list(dict_coeff.index):
                #print(dict_coeff.loc[r.id,'kcat_MW'])
                dict_obj[r.forward_variable] = 1 / dict_coeff.loc[r.id,'kcat_MW']

        model_obj = model.problem.Objective(Zero, direction="min", sloppy=True)
        model.objective = model_obj
        model.objective.set_linear_coefficients(dict_obj)

        solution = model.optimize()
    return solution

def draw_3d_rbas(z_data,out_fig_file):
    layout = go.Layout(template="none",plot_bgcolor='lightgrey')
    exlist = list(np.linspace(0,10,11))
    exlisto2 = list(np.linspace(0,10,11))
    fig = go.Figure(data=[go.Surface(y=exlist,x=exlisto2,z=z_data.values)],layout=layout)
    fig.update_layout(scene = dict(
        xaxis = dict(range=[0,10],tickfont=dict(size=13, family='Times New Roman'),backgroundcolor = "lightgrey",title=dict(text="<b>O2 uptake rates<br>(mmol/gDW/h)</b>",font=dict(size=18, family='Times New Roman'))),
        yaxis = dict(range=[0,10],tickfont=dict(size=13, family='Times New Roman'), backgroundcolor = "lightgrey",title=dict(text="<b>Glucose uptake rates<br>(mmol/gDW/h)</b>",font=dict(size=18, family='Times New Roman'))),
        zaxis = dict(range=[0,0.65],tickfont=dict(size=13, family='Times New Roman'), backgroundcolor = "grey", gridcolor = "white", title=dict(text="<b>Growth rates<br>(1/h)</b>",font=dict(size=18, family='Times New Roman')))))

    fig.update_traces(contours_z=dict(usecolormap=True, highlightcolor="mistyrose", project_z=True))
    # fig.update_xaxes(title_standoff = 1)
    # fig.update_traces(hovertemplate="none")
    fig.update_layout(autosize=False,scene_camera_eye=dict(x=-0.8, y=-2.1, z=0.3),
        width=850, height=850,margin=dict(l=20, r=20, b=20, t=20))
    fig.update_scenes(yaxis_tickangle=0)
    fig.update_scenes(xaxis_tickangle=0)

    fig.write_image(out_fig_file) 
    return fig
