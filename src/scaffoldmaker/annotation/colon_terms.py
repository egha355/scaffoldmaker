"""
Common resource for colon annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
colon_terms = [
    ( "ascending colon", "UBERON:0001156", "FMA:14545", "ILX:0734393"),
    ( "descending colon", "UBERON:0001158", "FMA:14547", "ILX:0724444"),
    ( "caecum", "UBERON:0001153", "FMA:14541", "ILX:0732270"),
    ( "colon", "UBERON:0001155", "FMA:14543", "ILX:0736005"),
    ( "colonic mucosa", "UBERON:0000317", "FMA:14984", "ILX:0731046"),
    ( "distal colon", "UBERON:0008971", "ILX:0727523"),
    ( "mesenteric zone", None),
    ( "non-mesenteric zone", None),
    ( "proximal colon", "UBERON:0008972", "ILX:0733240"),
    ( "serosa of colon", "UBERON:0003335", "FMA:14990", "ILX:0736932"),
    ( "spiral colon", "UBERON:0010239", "ILX:0735018"),
    ( "submucosa of colon", "UBERON:0003331", "FMA:14985", "ILX:0728041"),
    ( "taenia coli", "UBERON:0012419", "FMA:15041", "ILX:0731555"),
    ( "taenia libera", "ILX:0739285"),
    ( "taenia mesocolica", "ILX:0739284"),
    ( "taenia omentalis", "ILX:0739286"),
    ( "transverse colon", "UBERON:0001157", "FMA:14546", "ILX:0728767")
    ]

def get_colon_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in colon_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Colon annotation term '" + name + "' not found.")
