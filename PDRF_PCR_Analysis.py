#!/bin/sh

import os
import re
import pickle
import time
import argparse


nucleotide_symbols = {'A': ['A'], 'C': ['C'], 'G': ['G'], 'T': ['T'], 'R': ['G', 'A'], 'Y': ['C', 'T'], 'K': ['G', 'T'], 'M': ['A', 'C'], 'S': ['G', 'C'], 'W': ['A', 'T'],
                      'B': ['G', 'T', 'C'], 'D': ['G', 'A', 'T'], 'H': ['A', 'C', 'T'], 'V': ['G', 'C', 'A'], 'N': ['A', 'G', 'C', 'T']}

default_dic_enzymes = {'AarI': 'CACCTGC', 'AatII': 'GACGTC', 'Aba13301I': 'GCAAAC', 'Aba6411II': 'CRRTAAG', 'AbaB8342IV': 'CATTAG', 'AbaCIII': 'CTATCAV',
                       'AbaPBA3II': 'CAYGAC', 'AbaUMB2I': 'YCCGSS', 'Acc65V': 'GACGCA', 'AccI': 'GTMKAC', 'AccIX': 'GACRAC', 'AccX': 'GGARCA', 'AceIII': 'CAGCTC',
                       'AchA6III': 'AGCCAG', 'AciI': 'CCGC', 'AclI': 'AACGTT', 'Aco12261II': 'CCRGAG', 'AcoY31II': 'TAGCRAB', 'AcyI': 'GRCGYC', 'AflII': 'CTTAAG',
                       'AflIII': 'ACRYGT', 'AgeI': 'ACCGGT', 'AgsI': 'TTSAA', 'AhaIII': 'TTTAAA', 'AhyRBAHI': 'GCYYGAC', 'AhyYL17I': 'YAAMGAG', 'AluI': 'AGCT',
                       'AmaCSI': 'GCTCCA', 'AoxI': 'GGCC', 'ApaI': 'GGGCCC', 'ApaLI': 'GTGCAC', 'ApoI': 'RAATTY', 'ApyPI': 'ATCGAC', 'AquIV': 'GRGGAAG',
                       'Asl11923II': 'GGGABCC', 'Asp103I': 'CGRAGGC', 'Asp114pII': 'AGCABCC', 'Asp337I': 'CARABGG', 'AspAMDIV': 'ACCCAC', 'AspJHL3II': 'CGCCCAG',
                       'AspNIH4III': 'AAGAACB', 'AspSLV7III': 'GTCTCA', 'Asu14238IV': 'CGTRAC', 'AsuII': 'TTCGAA', 'AteTI': 'GGGRAG', 'AvaI': 'CYCGRG',
                       'AvaII': 'GGWCC', 'AvaIII': 'ATGCAT', 'Avi249I': 'CTGCA', 'AvrII': 'CCTAGG', 'Awo1030IV': 'GCCRAG', 'Bag18758I': 'CCCGAG',
                       'BalI': 'TGGCCA', 'BamHI': 'GGATCC', 'BanLI': 'RTCAGG', 'Bau1417V': 'GTTCAG', 'Bbr52II': 'GGCGAG', 'Bbr57III': 'GTRAAYG',
                       'Bbr7017II': 'CGGGAG', 'Bbr7017III': 'GGRCAG', 'BbuB31II': 'CGRKA', 'BbvCI': 'CCTCAGC', 'BbvI': 'GCAGC', 'BbvII': 'GAAGAC',
                       'BccI': 'CCATC', 'Bce3081I': 'TAGGAG', 'Bce83I': 'CTTGAG', 'BcefI': 'ACGGC', 'BciVI': 'GTATCC', 'BclI': 'TGATCA',
                       'Bco11035III': 'GAAGCY', 'BetI': 'WCCGGW', 'BfiI': 'ACTGGG', 'Bga514I': 'GTRAAG', 'BglII': 'AGATCT', 'BinI': 'GGATC',
                       'Ble402II': 'GRAGCAG', 'BloAII': 'GAGGAC', 'BmgI': 'GKGCCC', 'BsaAI': 'YACGTR', 'BsbI': 'CAACAC', 'BscGI': 'CCCGT',
                       'BscXI': 'GCAGGC', 'BseMII': 'CTCAG', 'BsePI': 'GCGCGC', 'BseRI': 'GAGGAG', 'BseSI': 'GKGCMC', 'BseYI': 'CCCAGC', 'BsgI': 'GTGCAG',
                       'BsiI': 'CACGAG', 'BsmAI': 'GTCTC', 'BsmI': 'GAATGC', 'Bsp1407I': 'TGTACA', 'Bsp3004IV': 'CCGCAT', 'Bsp460III': 'CGCGCAG', 'BspGI': 'CTGGAC',
                       'BspHI': 'TCATGA', 'BspHII': 'GTAGAT', 'BspLU11I': 'ACATGT', 'BspMI': 'ACCTGC', 'BspMII': 'TCCGGA', 'BspNCI': 'CCAGA', 'BsrBI': 'CCGCTC',
                       'BsrDI': 'GCAATG', 'BsrI': 'ACTGG', 'BtgZI': 'GCGATG', 'BtrI': 'CACGTC', 'BtsI': 'GCAGTG', 'BtsIMutI': 'CAGTG', 'Cal14237I': 'GGTTAG',
                       'CalB3II': 'GRTTRAG', 'Cau10061II': 'GTTAAT', 'CauII': 'CCSGG', 'Cba13II': 'AGGAAT', 'Cbo67071IV': 'GCRGAAG', 'CchII': 'GGARGA',
                       'CchIII': 'CCCAAG', 'Cco11366VI': 'GAAGAA', 'Cco14983V': 'GGGTDA', 'Cco14983VI': 'GCYGA', 'CcrNAIII': 'CGACCAG', 'Cdi11397I': 'GCGCAG',
                       'Cdi13746V': 'RGAAAGR', 'Cdi13750III': 'CCGATCC', 'CdiI': 'CATCG', 'CdpI': 'GCGGAG', 'Cdu23823II': 'GTGAAG', 'Cfa8380I': 'GRGGAY',
                       'Cfr10I': 'RCCGGY', 'CfrI': 'YGGCCR', 'CfrMH16VI': 'CTAAAG', 'Cfupf3II': 'GARCAG', 'Cgl13032I': 'GGCGCA', 'Cgl13032II': 'ACGABGG',
                       'Cin11811I': 'TGKMCA', 'Cje265V': 'GKAAGC', 'Cje54107III': 'GKAAYC', 'CjeFIII': 'GCAAGG', 'CjeFV': 'GGRCA', 'CjeNIII': 'GKAAYG',
                       'CjeNV': 'CCYGA', 'Cko11077IV': 'TGACAG', 'Cla11845III': 'GCGAA', 'ClaI': 'ATCGAT', 'Cly7489II': 'AAAAGRG', 'Cre7908I': 'GCGGGA',
                       'Csp2014I': 'GGAGGC', 'CspL61I': 'TYGAYCT', 'CspX1II': 'ACCCCA', 'CstMI': 'AAGGAG', 'CviJI': 'RGCY', 'CviRI': 'TGCA', 'Dde51507I': 'CCWGG',
                       'Dpi3069I': 'GACAG', 'Dpi3084I': 'CGRAG', 'Dpi3090II': 'AAGRAG', 'DpnI': 'GATC', 'DrdII': 'GAACCA', 'DrdIV': 'TACGAC', 'DrdVI': 'GCAGCC',
                       'DrdVIII': 'ARGAGC', 'DsaI': 'CCRYGG', 'DspS02II': 'TGCCGAC', 'EciI': 'GGCGGA', 'Ecl35734I': 'GAAAYTC', 'Eco31I': 'GGTCTC',
                       'Eco4174I': 'GCACAG', 'Eco43896II': 'CRARCAG', 'Eco4465II': 'GAAABCC', 'Eco47III': 'AGCGCT', 'Eco57I': 'CTGAAG', 'Eco57MI': 'CTYCAG',
                       'Eco8164I': 'GCCKAG', 'Eco9020I': 'CGAABTT', 'Eco9699II': 'TAGARC', 'EcoBLMcrX': 'RCSRC', 'EcoE1140I': 'ACCYAC', 'EcoHSI': 'GGTAAG',
                       'EcoNIH6II': 'ATGAAG', 'EcoP15I': 'CAGCAG', 'EcoRI': 'GAATTC', 'EcoRII': 'CCWGG', 'EcoRV': 'GATATC', 'Eli8509II': 'CCGGAG', 'EsaSSI': 'GACCAC',
                       'Esp3007I': 'CAGAAG', 'Esp3I': 'CGTCTC', 'FaiI': 'YATR', 'FauI': 'CCCGC', 'Fba202Z8II': 'AGAAGG', 'Fco1691IV': 'GCVGAG', 'FinI': 'GGGAC',
                       'Fna13121I': 'TTGAYC', 'Fnu11326IV': 'CTTAATT', 'FnuDII': 'CGCG', 'FokI': 'GGATG', 'FspPK15I': 'GARGAAG', 'FtnUV': 'GAAACA',
                       'Gba708II': 'ATGCAC', 'GdiII': 'CGGCCR', 'GlaI': 'GCGC', 'Gru56503II': 'CARABGC', 'GsuI': 'CTGGAG', 'GsuPI': 'GTACAG', 'HaeI': 'WGGCCW',
                       'HaeII': 'RGCGCY', 'HaeIII': 'GGCC', 'HbaII': 'GCCCAG', 'HgaI': 'GACGC', 'HgiAI': 'GWGCWC', 'HgiCI': 'GGYRCC', 'HgiJII': 'GRGCYC',
                       'HhaI': 'GCGC', 'Hin4II': 'GAAGG', 'HindII': 'GTYRAC', 'HindIII': 'AAGCTT', 'HpaI': 'GTTAAC', 'HpaII': 'CCGG', 'HphI': 'GGTGA',
                       'Hpy99I': 'CGWCG', 'Hpy99XIII': 'GCCTA', 'Hpy99XIV': 'GGWTAA', 'HpyAXIV': 'GCGTA', 'HpyAXVI-mut1': 'CRTTAA', 'HpyG272XV': 'GAAAAG',
                       'HpyUM032XIV': 'GAAAG', 'HspMHR1II': 'GAGCAGC', 'Jma19592II': 'GRGCRAC', 'Kas9737III': 'CCCRAG', 'Kor51II': 'RTCGAG', 'Kpn156V': 'CRTGATT',
                       'Kpn327I': 'GACATC', 'Kpn9644II': 'GRACRAC', 'KpnI': 'GGTACC', 'KpnNH25III': 'CTRGAG', 'KpnNIH50I': 'GCYAAG', 'Kro7512II': 'ARCAGKC',
                       'KroI': 'GCCGGC', 'Ksp632I': 'CTCTTC', 'Lde4408II': 'ACAAAG', 'Lla047I': 'CTCCA', 'Lla047II': 'AGAAG', 'LlaG50I': 'CCGTKA',
                       'Lme32I': 'CTYCAA', 'LmnI': 'GGAGC', 'Lmo370I': 'AGCGCCG', 'Lmo911II': 'TAGRAG', 'Lpl1004II': 'AGGRAG', 'Lpn11417II': 'ACGAAT', 'LpnPI': 'CCDG',
                       'LsaDS4I': 'TGGAAT', 'Lsp48III': 'AGCACC', 'Lsp6406VI': 'CRAGCAC', 'NaeI': 'GCCGGC', 'Nal45188II': 'ACCAGC', 'NarI': 'GGCGCC',
                       'Nbr128II': 'ACCGAC', 'NcoI': 'CCATGG', 'NdeI': 'CATATG', 'NgoAVII': 'GCCGC', 'NhaXI': 'CAAGRAG', 'NheI': 'GCTAGC', 'NhoI': 'GCWGC',
                       'NlaCI': 'GTGATG', 'NlaIII': 'CATG', 'NmeA6CIII': 'GCCGAC', 'NmeAIII': 'GCCGAG', 'NpeUS61II': 'GATCGAC', 'NruI': 'TCGCGA', 'NspBII': 'CMGCKG',
                       'NspES21II': 'CRTTCAG', 'NspI': 'RCATGY', 'ObaBS10I': 'ACGAG', 'OspHL35III': 'YAGGAG', 'Pac19842II': 'CCTTGA', 'PacIII': 'GTAATC',
                       'Pae10662III': 'TGACGAG', 'Pae8506I': 'CATCGAR', 'PaeCFIORFAP': 'ACGACC', 'PaePA99III': 'AAGAYC', 'Pal408I': 'CCRTGAG', 'PasI': 'CCCWGGG',
                       'Pba2294I': 'GTAAG', 'Pbu13063II': 'GTATYC', 'PcaII': 'GACGAG', 'Pcr308II': 'CCAAAG', 'Pdu1735I': 'CACCAC', 'PenI': 'GCAGT',
                       'Pfl10783II': 'GCGTCAG', 'Pfl1108I': 'TCGTAG', 'PflPt14I': 'RGCCCAC', 'PfrJS12V': 'GGCGGAG', 'PgaP73III': 'TTCGAG', 'Pin17FIII': 'GGYGAB',
                       'PinP23II': 'CTRKCAG', 'PleI': 'GAGTC', 'PliMI': 'CGCCGAC', 'PmaCI': 'CACGTG', 'Pme10899I': 'GACAGG', 'PpiP13II': 'CGCRGAC', 'PpuMI': 'RGGWCCY',
                       'Pse18267I': 'RCCGAAG', 'PsiI': 'TTATAA', 'Psp0357II': 'GCGAAG', 'PspD7DII': 'CCGCGAG', 'PspMR102II': 'CAAGAAC', 'PspOMII': 'CGCCCAR',
                       'PspPRI': 'CCYCAG', 'PspR84I': 'TACYCAC', 'Pst145I': 'CTAMRAG', 'Pst273I': 'GATCGAG', 'PstI': 'CTGCAG', 'PstII': 'CTGATG', 'PsuGI': 'BBCGD',
                       'PvuI': 'CGATCG', 'PvuII': 'CAGCTG', 'Ran11014IV': 'GAAAGAG', 'Rba2021I': 'CACGAGH', 'RceI': 'CATCGAC', 'RdeGBI': 'CCGCAG', 'RdeGBII': 'ACCCAG',
                       'RdeGBIII': 'TGRYCA', 'Rer8036II': 'CCGAKGG', 'RflFIII': 'CGCCAG', 'Rgo13296IV': 'GRAAGCG', 'Rho5650I': 'AACGAG', 'RlaII': 'ACACAG',
                       'RleAI': 'CCCACA', 'Rmu369III': 'GGCYAC', 'RpaB5I': 'CGRGGAC', 'RpaBI': 'CCCGCAG', 'RpaI': 'GTYGGAG', 'RpaTI': 'GRTGGAG', 'RsaI': 'GTAC',
                       'Rsp008IV': 'ACGCAG', 'Rsp008V': 'GCCCAT', 'Rsp531II': 'CACACG', 'RspPBTS2III': 'CTTCGAG', 'RsrII': 'CGGWCCG', 'SacI': 'GAGCTC',
                       'SacII': 'CCGCGG', 'Sag901I': 'GCAAAT', 'SalI': 'GTCGAC', 'SanDI': 'GGGWCCC', 'SapI': 'GCTCTTC', 'Sau5656II': 'GTTGCA', 'Sbo46I': 'TGAAC',
                       'ScaI': 'AGTACT', 'ScoDS2II': 'GCTAAT', 'SdeAI': 'CAGRAG', 'SduI': 'GDGCHC', 'Sen17963III': 'CCAAAC', 'Sen5794III': 'ACGAACB',
                       'Sen6480IV': 'GTTCAT', 'SenSARA26III': 'ACRCAG', 'SenTFIV': 'GATCAG', 'Sep11964I': 'CGYCAT', 'SetI': 'ASST', 'SexAI': 'ACCWGGT',
                       'SfaNI': 'GCATC', 'SfeI': 'CTRYAG', 'Sgr7807I': 'GCCGAGG', 'SgrAII': 'CGAGATC', 'SgrTI': 'CCDS', 'SimI': 'GGGTC', 'Sma10259II': 'CAAAGA',
                       'Sma325I': 'ARCCCT', 'SmaI': 'CCCGGG', 'SmaUMH5I': 'CTTGAC', 'SmaUMH8I': 'GCGAACB', 'SmlI': 'CTYRAG', 'Sna507VIII': 'CRTTGAG',
                       'SnaBI': 'TACGTA', 'SnaI': 'GTATAC', 'Sno506I': 'GGCCGAG', 'Spe19205IV': 'GGACY', 'SpeI': 'ACTAGT', 'SphI': 'GCATGC', 'SplI': 'CGTACG',
                       'SpnRII': 'TCGAG', 'SpoDI': 'GCGGRAG', 'Sse8647I': 'AGGWCCT', 'Ssp6803IV': 'GAAGGC', 'Ssp714II': 'CGCAGCG', 'SspI': 'AATATT',
                       'SstE37I': 'CGAAGAC', 'Sth132I': 'CCCG', 'Sth20745III': 'GGACGAC', 'SthSt3II': 'GAAGT', 'StuI': 'AGGCCT', 'StyI': 'CCWWGG',
                       'SurP32aII': 'ACRGAG', 'TagI': 'ACGT', 'TaqI': 'TCGA', 'TaqII': 'GACCGA', 'TaqIII': 'CACCCA', 'TatI': 'WGTACW', 'TauI': 'GCSGC',
                       'TfiI': 'GAWTC', 'TkoII': 'TTCAAG', 'TpyTP2I': 'ACCAAG', 'TseI': 'GCWGC', 'TsoI': 'TARCCA', 'Tsp45I': 'GTSAC', 'TspARh3I': 'GRACGAC',
                       'TspDTI': 'ATGAA', 'TspEI': 'AATT', 'TspGWI': 'ACGGA', 'TspRI': 'CASTG', 'TsuI': 'GCGAC', 'Tth111II': 'CAARCA', 'UbaF11I': 'TACGA',
                       'UbaPI': 'CGAACG', 'Van9116I': 'CCKAAG', 'VchE4II': 'RTAAAYG', 'VpaSKIII': 'CGTCAG', 'VspI': 'ATTAAT', 'Vtu19109I': 'CACRAYC',
                       'WviI': 'CACRAG', 'XbaI': 'TCTAGA', 'Xca85IV': 'TACGAG', 'XhoI': 'CTCGAG', 'XhoII': 'RGATCY', 'XmaIII': 'CGGCCG', 'Yps3606I': 'CGGAAG',
                       'Yru12986I': 'AGGAAG'}


def create_regular_seq(seq: str):
    "Функция принимает на вход ДНК последовательность и возращает в виде регулярного выражения для использования всех вхождений подстроки. ПРИМЕР: 'ACYGT' -> 'AC[CT]GT' "
    result = ''
    for nucleotide in seq:
        s = nucleotide_symbols[nucleotide]
        if len(s) == 1:
            result += str(*s)
        else:
            result += f"[{''.join(s)}]"
    return result


def create_dict_restriction_enzymes(list_restriction_enzymes, filter=[]):
    result = {}
    if type(list_restriction_enzymes) == dict:
        if len(filter) != 0:
            return {k: v for k, v in list_restriction_enzymes.items() if k in filter}
        else:
            return list_restriction_enzymes
    elif all([type(list_restriction_enzymes) != dict, type(list_restriction_enzymes) == list, len(list_restriction_enzymes) > 1]):
        result = dict(zip([list_restriction_enzymes[i] for i in range(len(list_restriction_enzymes)) if i % 2 == 0], [
                      list_restriction_enzymes[i] for i in range(len(list_restriction_enzymes)) if i % 2 != 0]))
        return result
    elif all([type(list_restriction_enzymes) != dict, len(list_restriction_enzymes) > 1]):
        path = list_restriction_enzymes
        with open(path, 'rb') as file_dict_enzymes:
            loaded_dict = pickle.load(file_dict_enzymes)
            result = dict(loaded_dict)
        result = {k: v for k, v in result.items() if set(
            v).issubset(set(nucleotide_symbols.keys()))}
        return result


def create_matrix_count_restriction_sites_in_seq(dir_with_sequences, dict_enzymes, filter=[]):
    result = {}
    path_dir_with_files_sequences = dir_with_sequences
    files_sequences = os.listdir(path_dir_with_files_sequences)
    for file in files_sequences:
        if any(tag in file for tag in filter):
            sequence = ''
            with open(path_dir_with_files_sequences+'\\'+file, 'r') as file_gene:
                f_line = file_gene.readline()
                sequence = file_gene.read()
            for k, v in sorted(dict_enzymes.items()):
                match = re.finditer(create_regular_seq(v), sequence)
                result.setdefault(k, [])
                result[k] += [len([m.start() for m in match])]
    return result


def create_table_of_restriction_sites_in_sequences(dir_with_sequences, dict_enzymes, filter=[]):
    result = []
    path_dir_with_files_sequences = dir_with_sequences
    files_sequences = os.listdir(path_dir_with_files_sequences)
    for file in files_sequences:
        if any(tag in file for tag in filter):
            sequence = ''
            with open(path_dir_with_files_sequences+'\\'+file, 'r') as file_gene:
                f_line = file_gene.readline()
                sequence = file_gene.read()
            count = ''
            for k, v in sorted(dict_enzymes.items()):
                match = re.finditer(create_regular_seq(v), sequence)
                finded_matches = [m.start()+1 for m in match]
                count += str(len(finded_matches)) + str(finded_matches) + ';'
            result.append(file + ';' + str(len(sequence)) + ';' + count + ';')
    return result


def list_resrictase(dic_enzymes):
    print(dic_enzymes)


parser = argparse.ArgumentParser(
    description="Create table site restrictions for DNA sequences")
subparser = parser.add_subparsers(dest='command')
analyze_parser = subparser.add_parser('analyze', help="Analyze files")
analyze_parser.add_argument(
    'indir', type=str, help='Input path to files (dir) or file')
analyze_parser.add_argument('-o', '--outdir', type=str, default=os.getcwd(
), help='Output path file for save data table. For save results default used directory that contains script file')
analyze_parser.add_argument('-rs', '--restrictase', type=str, nargs='+', default=default_dic_enzymes,
                            help='By default, it contains an internal dictionary of restriction sites. The list of restriction sites can be viewed using the command list sites')
analyze_parser.add_argument('-ff', '--filter_files', type=str, default=['.fna'], nargs='+',
                            help="Filter for name or type file, default .fna. Example: Seq .fna. You will analysis files that contains in the name Seq and type .fna: Seq234.fna, Seq234.fna. Files name's 123.fna, Seq.gbk not alnalysis, or others files not contains in the name Seq and .fna")
analyze_parser.add_argument('-fcs', '--filter_count_sites', default=['100', '100'], type=str, nargs='+',
                            help='Filter count sites of enzyme in each sequence. Exmaple 0 0. First number is count sites enzymes in analyzed sequences, second number is max count sites positions restriction enzymes in sequence.')
analyze_parser.add_argument(
    '-fr', '--filter_restrictase', type=str, default=[], nargs='+')
restrictase_parser = subparser.add_parser(
    'enzymes', help="List of default restrictase")
args = parser.parse_args()


if __name__ == '__main__':
    if args.command == 'enzymes':
        list_resrictase(default_dic_enzymes)
    elif args.command == 'analyze':
        path_to_file_dict_enzymes = args.restrictase
        path_to_files_seqences = args.indir
        output_data_dir = args.outdir
        filter_files = args.filter_files
        filter_sites = list(map(int, args.filter_count_sites))
        filter_restrictase = args.filter_restrictase
        # Загрузка словаря ферментов с сайтами рестрикциями через файл формата .pkl
        dict_enzymes = create_dict_restriction_enzymes(
            path_to_file_dict_enzymes, filter_restrictase)
        print(f'Фильтр файлов: {filter_files}')
        print(f'Список рестриктаз: {dict_enzymes}')
        # Создание матрицы количества встречаемости сайтов рестрикции для каждого фермента в каждом файле с последовательностью.
        matrix_sites_restriction = create_matrix_count_restriction_sites_in_seq(
            path_to_files_seqences, dict_enzymes, filter_files)
        print(
            f'Количество сайтов рестрикатаз в анализе: {len(matrix_sites_restriction)}')
        # Фильтрация ферментов по значениям в матрице, по максимальному количеству сайтов для фермента и количеству всех сайтов рестрикций ферментов в словаре для последовательности
        # k = имя фермента рестрикции, v = список с количеством сайтов рестрикции для каждой последовательности
        # v.count() = можно задать фильтр по количеству встречаемости сайтов рестрикции в последовательностях
        # max(v) = фильтр по максимальной встречаемости сайта рестрикции
        filtered_matrix_sites_restriction = {k: v for k, v in matrix_sites_restriction.items(
        ) if v.count(0) <= filter_sites[0] and max(v) <= filter_sites[1]}
        print(
            f'Количесвто отобранных сайтов рестрикции: {len((filtered_matrix_sites_restriction))}')
        # Создание нового словаря ферментов после фильтрации
        filtered_dict_enzymes = {k: v for k, v in dict_enzymes.items(
        ) if k in filtered_matrix_sites_restriction.keys()}
        # Создание матрицы всех отобранных ферментов по количеству сайтов рестриктаз их начальной позиции в каждой входящей последовательности
        output_data = create_table_of_restriction_sites_in_sequences(
            path_to_files_seqences, filtered_dict_enzymes, filter_files)
        # Запись результата в файл формата .csv
        output_data_name_enzymes_seq = ';' + 'Lenght sequence' + ';' + \
            ';'.join(
                sorted(list({k+'-'+v for k, v in filtered_dict_enzymes.items()})))
        get_time = time.gmtime()
        date = f'{get_time.tm_mday}{get_time.tm_mon}{get_time.tm_hour}{get_time.tm_min}{get_time.tm_sec}{get_time.tm_year}'

        if args.outdir == os.getcwd():
            output_path = output_data_dir+'\\' + 'Result_' + \
                ''.join(filter_files) + str(date) + '.csv'
        else:
            output_path = output_data_dir
        with open(output_path, '+a') as output_file:
            output_file.write(output_data_name_enzymes_seq+'\n')
            for string in output_data:
                output_file.write(string+'\n')
        print(f'Данные записаны в {output_path}')
