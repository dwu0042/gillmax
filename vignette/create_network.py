import networkx as nx
import numpy as np

def create_twolayer_bipartite_network(n_indvs, n_hh, n_comm, p_comm, seed=None) -> nx.Graph:
    
    rng = np.random.default_rng(seed=seed)

    # group 1 layer represents households
    household_size = np.round(rng.dirichlet([1]*n_hh) * (n_indvs-n_hh)) + 1

    ns_hh = np.sum(household_size)

    while ns_hh != n_indvs:
        print("fixing household size sample")
        household_size[-1] = household_size[-1] + n_indvs - ns_hh
        if household_size[-1] < 1:
            household_size = household_size[:-1]
        ns_hh = np.sum(household_size)

    household_size = household_size.astype(int)

    household = nx.bipartite.generators.configuration_model(
        aseq=[1]*n_indvs,
        bseq=household_size,
        create_using=nx.Graph,
        seed=rng,
    )

    # group 2 layer is the community layer
    community = nx.bipartite.generators.preferential_attachment_graph(
        aseq=[n_comm]*n_indvs,
        p=p_comm,
        create_using=nx.Graph,
        seed=rng
    )

    # relabel nodes
    household = nx.relabel_nodes(household, {i: (i if i < n_indvs else f"HH{i}") for i in household.nodes})
    community = nx.relabel_nodes(community, {i: (i if i < n_indvs else f"CC{i}") for i in community.nodes})

    combined = nx.compose(household, community)

    return combined

def paper_network() -> nx.Graph:
    return create_twolayer_bipartite_network(
        n_indvs=100,
        n_hh=25,
        n_comm=2,
        p_comm=0.25,
        seed=0x12b342b71
    )