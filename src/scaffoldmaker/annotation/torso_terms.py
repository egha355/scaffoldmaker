"""
Common resource for torso annotation terms.
"""

# convention: preferred name, preferred id, followed by any other ids and alternative names
torso_terms = [
    ("anterior of torso", "FMA:15900"),
    ("posterior of torso", "FMA:15901"),
   ]

def get_torso_term(name : str):
    """
    Find term by matching name to any identifier held for a term.
    Raise exception if name not found.
    :return ( preferred name, preferred id )
    """
    for term in torso_terms:
        if name in term:
            return ( term[0], term[1] )
    raise NameError("Torso annotation term '" + name + "' not found.")