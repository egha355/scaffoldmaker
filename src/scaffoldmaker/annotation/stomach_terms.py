"""
Common resource for stomach annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
stomach_terms = [
    ( "antrum on serosa split margin", "None"),
    ( "body of stomach", "UBERON:0001161", " FMA:14560", "ILX:0724929"),
    ( "body on serosa split margin", "None"),
    ( "cardia of stomach", "UBERON:0001162", " FMA:14561", "ILX:0729096"),
    ( "dorsal stomach", "None"),
    ( "duodenum", "UBERON:0002114", " FMA:7206", "ILX:0726125"),
    ( "duodenum on greater curvature", "None"),
    ( "esophagus", "UBERON:0001043", "FMA: 7131", "ILX:0735017"),
    ( "esophagus on serosa split margin", "None"),
    ( "esophagogastric junction", "UBERON:0007650", "FMA: 9434", "ILX:0733910"),
    ( "forestomach-glandular stomach junction", "UBERON:0012270", "ILX:0729974"),
    ( "forestomach-glandular stomach junction on inner wall", "None"),
    ( "forestomach-glandular stomach junction on outer wall", "None"),
    ( "fundus of stomach", "UBERON:0001160", " FMA:14559", "ILX:0724443"),
    ( "fundus on serosa split margin", "None"),
    ( "gastro-esophagal junction on lesser curvature", "None"),
    #( "greater curvature of stomach", "UBERON:0001164", "FMA: 14574", "ILX:0724395"),
    ( "junction between fundus and body on greater curvature", "None"),
    # ( "lesser curvature of stomach", "UBERON:0001163", "FMA: 14572", "ILX:0733753"),
    ( "limiting ridge on greater curvature", "None"),
    ( "mucosa of stomach", "UBERON:0001199", "FMA:14907", "ILX:0736669"),
    ( "pylorus", "UBERON:0001166", " FMA:14581", "ILX:0734150"),
    ( "pyloric antrum", "UBERON:0001165", " FMA:14579", "ILX:0728672"),
    ( "pylorus on greater curvature", "None"),
    ( "pylorus on serosa split margin", "None"),
    ( "serosa of dorsal stomach", "None"),
    ( "serosa of stomach", "UBERON:0001201", "FMA:14914", "ILX:0735818"),
    ( "serosa split margin", "None"),
    ( "serosa of ventral stomach", "None"),
    ( "stomach", "UBERON:0000945", "FMA:7148", "ILX:0736697"),
    ( "ventral stomach", "None"),
    ]

def get_stomach_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in stomach_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Stomach annotation term '" + name + "' not found.")
