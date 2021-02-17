#domain overlap and selection
#handling of overlaping domains and selecting the best domain assigment
#also working with combining annotation/modeling methods
#creation of full height estimates


## Domain assigment for Human Extracellular domains

#Using the domain level output from the hmmscan of the Pfam-A hmm profiles, assign each extracellular domain, non overlapping domains.

#Steps:
#1. Find best domain assignment for each ECD
#2. pick that domain
#3. find next best, check for overlap
#4. if no overlap keep, if overlap detected discard

def parse_hmm_scan(scan_path):

    hmmer_dom = pd.read_csv(scan_path, skiprows = 3, header = None,
                    delim_whitespace= True, skipfooter = 10, engine='python',
                   names = ['target_name','accession', 'tlen', 'query_name','Pfam', 'qlen', 'E-value_full',
                            'score_full', 'bias_full', 'dom_num', 'type_count', 'c-Evalue_dom','i-Evalue_dom',
                            'score_dom', 'bias_dom', 'hmm_from', 'hmm_to', 'align_from', 'align_to',
                            'env_from', 'env_to', 'accuracy','description', 'nothing'])
`   hmmer_dom.drop(columns = ['accession', 'description', 'nothing'], inplace = True)

    return hmmer_dom

def threshold_domains(hmmer_dom, full_sig == 0.0001, dom_sig = 2, coverage = 0.9):
    #filtering for specific domain assignments and then lest specific individual domains seems to give best agreement with notch1
    filtered_dom = hmmer_dom[(hmmer_dom['E-value_full']  < full_sig) & (hmmer_dom['i-Evalue_dom']  < dom_sig)]

    #assign pfam names
    filtered_dom = filtered_dom.assign(Pfam = filtered_dom['Pfam'].str.split('.').str[0])

    #length of domain assignement must cover 90% of the domain hmm profile
    filtered_tm = filtered_dom[filtered_dom['tlen'] > coverage*filtered_dom['qlen']]

    return filtered_tm

def connections(graph, id):
    def dict_to_df(d):
        df = pd.DataFrame(data=[d.keys(), d.values()], index=['ID', 'Subgraph']).T
        df['id'] = id
        return df[['id', 'Subgraph', 'ID']]

    def dfs(node, num):
        visited[node] = num
        for _node in graph.loc[node].iloc[0]:
            if _node not in visited:
                dfs(_node, num)

    visited = {}
    graph = graph.loc[id]
    for (num, node) in enumerate(graph.index):
        if node not in visited:
            dfs(node, num)

    return dict_to_df(visited)
