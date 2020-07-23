"""
Common resource for small intestine annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
smallintestine_terms = [
    ( "duodenum", "UBERON:0002114", "FMA:7206", "ILX:0726125"),
    ( "ileum", "UBERON:0002116", "FMA:7208", "ILX:0728151"),
    ( "jejunum", "UBERON:0002115", "FMA:7207", "ILX:0724224"),
    ( "small intestine", "UBERON:0002108", "FMA:7200", "ILX:0726770")
    ]

def get_smallintestine_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in smallintestine_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Small intestine annotation term '" + name + "' not found.")
