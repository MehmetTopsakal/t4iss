from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.chemenv.coordination_environments\
     .coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments\
     .chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments\
     .structure_environments import LightStructureEnvironments


def get_cesym(lgf, structure, site):

    # doc: http://pymatgen.org/_modules/pymatgen/analysis/chemenv/
    #      coordination_environments/coordination_geometry_finder.html
    lgf.setup_structure(structure)

    # doc: http://pymatgen.org/_modules/pymatgen/analysis/chemenv/
    #      coordination_environments/
    #      chemenv_strategies.html#MultiWeightsChemenvStrategy.
    #      stats_article_weights_parameters
    strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
    
    # returns all information about the structure; se is a structure object
    se = lgf.compute_structure_environments(maximum_distance_factor=1.2,
                                            only_cations=False,
                                            only_indices=[site])

    lse = LightStructureEnvironments.\
          from_structure_environments(strategy=strategy,
                                      structure_environments=se)
    coor = lse.coordination_environments
    
    # ce = chemical environment
    # csm = continuous symmetry measure
    # from Waroquiers et al (verbatim)
    # DOI: 10.1021/acs.chemmater.7b02766
    # "The environment of the atom is then the model polyhedron for which 
    # the similarity is the highest, that is, for which the CSM is the lowest."
    # in this case, it looks like O:6 (octahedral?)
    return [coor[site][0]['ce_symbol'], coor[site][0]['csm']]


def main():
    pass

if __name__ == '__main__':
    main()