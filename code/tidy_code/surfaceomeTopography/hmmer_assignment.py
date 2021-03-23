#assignment of domains based on selection of best from overlapping

def hmmer_import(filepath, filter = True, Evalue = 0.0001, i_Evalue = 2):
    hmmer_dom = pd.read_csv(filepath, skiprows = 3, header = None,
                    delim_whitespace= True, skipfooter = 10, engine='python',
                   names = ['target_name','accession', 'tlen', 'query_name','Pfam', 'qlen', 'E-value_full',
                            'score_full', 'bias_full', 'dom_num', 'type_count', 'c-Evalue_dom','i-Evalue_dom',
                            'score_dom', 'bias_dom', 'hmm_from', 'hmm_to', 'align_from', 'align_to',
                            'env_from', 'env_to', 'accuracy','description', 'nothing'])
    hmmer_dom.drop(columns = ['accession', 'description', 'nothing'], inplace = True)

    if filter:
        #filtering for specific domain assignments and then lest specific individual domains seems to give best agreement with notch1
        filtered_dom = hmmer_dom[(hmmer_dom['E-value_full']  < Evalue) & (hmmer_dom['i-Evalue_dom']  < i_Evalue)]
        filtered_dom = filtered_dom.assign(Pfam = filtered_dom['Pfam'].str.split('.').str[0])

        #length of domain assignement must cover 90% of the domain hmm profile
        filtered = filtered_dom[filtered_dom['tlen'] > 0.9*filtered_dom['qlen']]

        return filtered

    else:
        return hmmer_dom
