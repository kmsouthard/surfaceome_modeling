#glycosite identification

def parse_glycosites(gff, surfaceome_accessions):
    glycosylation = gff[(gff['td'] == 'Glycosylation') & (gff['ID link'].isin(surfaceome_accessions))]

    glycosylation = glycosylation.assign(glycosite = glycosylation['start'])

    glycosylation.drop(columns = ['start', 'end'], inplace = True)

    return glycosylation

def surface_glycosites(glycosylation, surfaceome):

    surfaceome_gly = surfaceome[['ID link', 'start', 'end']].merge(glycosylation, on ='ID link', how ='left')

    ecd_glycosites = surfaceome_gly[
                    (surfaceome_gly['glycosite'] <= surfaceome_gly['end'])
                    &(surfaceome_gly['glycosite'] >= surfaceome_gly['start'])]

    return ecd_glycosites
