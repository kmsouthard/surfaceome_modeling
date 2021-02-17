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

#Finding the best domains with a little graph therory:
#https://stackoverflow.com/questions/54969074/python-3-remove-overlaps-in-table

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

def find_overlaps(data, allowed_overlap = 4):
    #make copy of data so not mutating it
    df = data.copy()
    #set an id for each row
    df.loc[:,'ID'] = range(df.shape[0])
    #use pd.interval for intervals (add/subtract allowed overlap(total allowed overlap = 8))
    df.loc[:,'Interval'] = df.apply(lambda x: pd.Interval(x['align_from']+allowed_overlap, x['align_to']-allowed_overlap, closed='neither'), axis=1)

    return df

def overlap_graph(overlap_df):

    columns = ['target_name', 'Interval', 'ID']
    connected = overlap_df[columns].merge(overlap_df[columns], on='target_name')
    connected['Overlap'] = connected.apply(lambda x: x['Interval_x'].overlaps(x['Interval_y']), axis=1)
    connected = connected.loc[connected['Overlap'] == True, ['target_name', 'ID_x', 'ID_y']]

    graph = connected.groupby(['target_name', 'ID_x']).agg(list)

    return graph

def find_connections():
    dfs = []
    for id in graph.index.get_level_values(0).unique():
        dfs.append(connections(graph, id))

    conns = pd.concat(dfs)

    return conns

def annotate_subgraphs(connections, overlaps):
    #merge the connections and overlaps dataframes so each subgraph is annotated
    data = overlaps.merge(connections[['Subgraph', 'ID']], on=['ID'])

def select_most_sig(x):
    m = x['i-Evalue_dom'].min()
    if len(x) > 1 and (x['i-Evalue_dom'] == m).all():
        return -1
    else:
        return x['i-Evalue_dom'].idxmin()

def select_hightest_score(x):
    m = x['score_dom'].min()
    if len(x) > 1 and (x['score_dom'] == m).all():
        return -1
    else:
        return x['score_dom'].idxmax()

def select_something(x, something):
    m = x[something].min()
    if len(x) > 1 and (x['score_dom'] == m).all():
        return -1
    else:
        return x['score_dom'].idxmax()


#selected = data.groupby(['target_name', 'Subgraph'])['score_dom', 'ID'].apply(select_max)
#selected = selected[selected >= 0]
